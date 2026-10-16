"""Input sky.

v0.1.1: Isotropic backgrounds observed with the exact Boost (ADR 0002),

    I_obs(ν, n̂) = (ν / ν′)³ I_rest(ν′),   ν′ = ν γ (1 − β·n̂),

with β = (v_sun + v_orbit(t)) / c and n̂ the observed direction. The beam
acts on the anisotropic part only, through b₁ = exp(−σ²).

v0.1.2: the Static sky, following the design of litebird_sim's
``SkyGenerator``: pysm3 maps per frequency -> alm, Gaussian beam, rotation to
the J2000 ecliptic, cached on disk, sampled per IFG with
``ducc0.sht.synthesis_general``. It is defined in the barycentric frame and
receives the first-order orbital Doppler term (ADR 0003).
"""

import hashlib
import json
from collections.abc import Callable
from pathlib import Path

import ducc0
import h5py
import healpy as hp
import numpy as np
from astropy import units as u
from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord

from .config import EXTRAGALACTIC_STATIC_MODELS, SkyConfig, StaticSkyConfig
from .constants import C_LIGHT_M_S
from .spectra import (
    blackbody_mjy_sr,
    cib_mjy_sr,
    extragalactic_co_mjy_sr,
    mu_distortion_mjy_sr,
    relativistic_sz_correction_mjy_sr,
    y_distortion_mjy_sr,
)

C_LIGHT_KM_S = C_LIGHT_M_S / 1e3


def isotropic_background_mjy_sr(freq_ghz: np.ndarray, sky: SkyConfig) -> np.ndarray:
    """Sum of the Isotropic backgrounds in the CMB rest frame, in MJy/sr."""
    t_cmb_k = sky.t_cmb_k
    total = np.zeros(np.shape(freq_ghz))
    if sky.include_cmb_blackbody:
        total += blackbody_mjy_sr(freq_ghz, t_cmb_k)
    if sky.y:
        total += y_distortion_mjy_sr(freq_ghz, sky.y, t_cmb_k)
        total += relativistic_sz_correction_mjy_sr(
            freq_ghz, sky.y, sky.y_electron_temperature_kev, t_cmb_k
        )
    if sky.mu:
        total += mu_distortion_mjy_sr(freq_ghz, sky.mu, t_cmb_k)
    if sky.cib_amplitude:
        total += cib_mjy_sr(
            freq_ghz,
            sky.cib_amplitude,
            sky.cib_spectral_index,
            sky.cib_temperature_k,
            sky.cib_reference_freq_ghz,
        )
    if sky.co_amplitude:
        total += extragalactic_co_mjy_sr(freq_ghz, sky.co_amplitude)
    return total


def solar_velocity_km_s(sky: SkyConfig) -> np.ndarray:
    """Velocity of the Sun with respect to the CMB, J2000 ecliptic, in km/s."""
    direction = SkyCoord(
        l=sky.solar_velocity_gal_lon_deg * u.deg,
        b=sky.solar_velocity_gal_lat_deg * u.deg,
        frame="galactic",
    ).transform_to(BarycentricMeanEcliptic(equinox="J2000"))
    return sky.solar_velocity_km_s * direction.cartesian.xyz.value


def dipole_beam_factor(beam_fwhm_deg: float) -> float:
    """Gaussian beam window at ℓ = 1, b₁ = exp(−σ²).

    >>> round(1 - dipole_beam_factor(1.6), 6)
    0.000141
    """
    sigma_rad = np.deg2rad(beam_fwhm_deg) / np.sqrt(8 * np.log(2))
    return float(np.exp(-(sigma_rad**2)))


def boosted_isotropic_background_mjy_sr(
    freq_ghz: np.ndarray,
    direction: np.ndarray,
    velocity_km_s: np.ndarray,
    sky: SkyConfig,
    beam_factor: float = 1.0,
) -> np.ndarray:
    """Observed Isotropic backgrounds along each direction.

    Parameters
    ----------
    freq_ghz : ndarray, shape (n_freq,)
    direction : ndarray, shape (n, 3)
        Observed directions (unit vectors).
    velocity_km_s : ndarray, shape (n, 3) or (3,)
        Observer velocity with respect to the CMB.
    sky : SkyConfig
    beam_factor : float
        Factor applied to the anisotropic part (b₁ of the beam, 1 = no beam).

    Returns
    -------
    ndarray, shape (n, n_freq), in MJy/sr
    """
    beta = np.broadcast_to(velocity_km_s, direction.shape) / C_LIGHT_KM_S
    gamma = 1.0 / np.sqrt(1.0 - np.sum(beta**2, axis=1))
    factor = gamma * (1.0 - np.sum(beta * direction, axis=1))  # ν′ / ν
    observed = (
        isotropic_background_mjy_sr(freq_ghz[None, :] * factor[:, None], sky)
        / factor[:, None] ** 3
    )
    if beam_factor != 1.0:
        rest = isotropic_background_mjy_sr(freq_ghz, sky)[None, :]
        observed = rest + beam_factor * (observed - rest)
    return observed


# Relative accuracy of ducc0's synthesis at arbitrary positions
SYNTHESIS_EPSILON = 1e-9
# ponytail: all cores for ducc0; set explicitly if several MPI ranks share a node
DUCC0_NTHREADS = 0

# (freq_ghz) -> (Galactic map, extragalactic map), MJy/sr, Galactic coordinates
EmissionFunction = Callable[[float], tuple[np.ndarray, np.ndarray]]


def pysm3_emission(models: tuple[str, ...], nside: int) -> EmissionFunction:
    """pysm3 total intensity of the Galactic and of the extragalactic models."""
    import pysm3
    import pysm3.units as pysm_units

    def sky(names):
        if not names:
            return None
        return pysm3.Sky(nside=nside, preset_strings=list(names), output_unit="MJy/sr")

    galactic = sky([m for m in models if not m.startswith(EXTRAGALACTIC_STATIC_MODELS)])
    extragalactic = sky(
        [m for m in models if m.startswith(EXTRAGALACTIC_STATIC_MODELS)]
    )
    zero = np.zeros(hp.nside2npix(nside))

    def emission(freq_ghz: float) -> tuple[np.ndarray, np.ndarray]:
        freq = freq_ghz * pysm_units.GHz
        return tuple(
            zero if s is None else np.asarray(s.get_emission(freq)[0].value)
            for s in (galactic, extragalactic)
        )

    return emission


def static_sky_cache_key(
    config: StaticSkyConfig, freq_ghz: np.ndarray, fwhm_deg: float
) -> str:
    """Hash of everything that determines the Static-sky cube."""
    import pysm3

    content = {
        "format": 1,
        "models": list(config.models),
        "nside_eval": 2 * config.nside,
        "lmax": config.lmax,
        "fwhm_deg": fwhm_deg,
        "freq_ghz": [float(f) for f in freq_ghz],
        "pysm3": pysm3.__version__,
    }
    return hashlib.sha256(json.dumps(content, sort_keys=True).encode()).hexdigest()


def build_static_sky_alm(
    config: StaticSkyConfig,
    freq_ghz: np.ndarray,
    fwhm_deg: float,
    emission: EmissionFunction | None = None,
) -> np.ndarray:
    """Static-sky alm cube, shape (n_freq, 1, n_alm), MJy/sr, J2000 ecliptic.

    For each frequency: map -> alm (3 iterations), monopole of the
    extragalactic part removed, Gaussian beam, rotation Galactic -> ecliptic.
    """
    if emission is None:
        emission = pysm3_emission(config.models, 2 * config.nside)
    lmax = config.lmax
    beam = hp.gauss_beam(np.deg2rad(fwhm_deg), lmax=lmax)
    rotator = hp.Rotator(coord=["G", "E"])
    cube = np.empty((freq_ghz.size, 1, hp.Alm.getsize(lmax)), dtype=np.complex128)
    for i, freq in enumerate(freq_ghz):
        galactic_map, extragalactic_map = emission(float(freq))
        alm = hp.map2alm(galactic_map, lmax=lmax, iter=3)
        extragalactic_alm = hp.map2alm(extragalactic_map, lmax=lmax, iter=3)
        extragalactic_alm[0] = (
            0.0  # a_00: the monopole belongs to the Isotropic backgrounds
        )
        alm = hp.almxfl(alm + extragalactic_alm, beam)
        cube[i, 0] = rotator.rotate_alm(alm, lmax=lmax)
    return cube


class StaticSky:
    """Static sky cube, loaded from (or written to) the cache, sampled per direction.

    Parameters
    ----------
    config : StaticSkyConfig
    freq_ghz : ndarray
        Frequency grid.
    fwhm_deg : float
        Gaussian beam FWHM.
    emission : callable, optional
        Map generator (default: pysm3 with ``config.models``).
    """

    def __init__(
        self,
        config: StaticSkyConfig,
        freq_ghz: np.ndarray,
        fwhm_deg: float,
        emission: EmissionFunction | None = None,
    ) -> None:
        self.config = config
        self.freq_ghz = np.asarray(freq_ghz)
        key = static_sky_cache_key(config, self.freq_ghz, fwhm_deg)
        self.path = Path(config.cache_dir) / f"static_sky_{key[:16]}.h5"
        if self.path.exists():
            with h5py.File(self.path, "r") as f:
                if f.attrs["cache_key"] != key:
                    raise RuntimeError(f"Cache key mismatch in {self.path}")
                self.alm = f["alm"][()]
        else:
            self.alm = build_static_sky_alm(config, self.freq_ghz, fwhm_deg, emission)
            self.path.parent.mkdir(parents=True, exist_ok=True)
            tmp = self.path.with_suffix(".tmp")
            with h5py.File(tmp, "w") as f:
                f.create_dataset("alm", data=self.alm)
                f.create_dataset("freq_ghz", data=self.freq_ghz)
                f.attrs["cache_key"] = key
                f.attrs["unit"] = "MJy/sr"
                f.attrs["coord"] = "J2000 ecliptic"
            tmp.replace(self.path)
        self.maps = None
        if config.sampling == "nearest":
            self.maps = np.stack(
                [hp.alm2map(a[0], config.nside, lmax=config.lmax) for a in self.alm]
            )

    def sample(self, direction: np.ndarray) -> np.ndarray:
        """Static sky along each direction, shape (n, n_freq), MJy/sr."""
        theta, phi = hp.vec2ang(direction)
        if self.maps is not None:
            pixels = hp.ang2pix(self.config.nside, theta, phi)
            return self.maps[:, pixels].T.copy()
        loc = np.column_stack([theta, phi])
        result = np.empty((direction.shape[0], self.freq_ghz.size))
        for i in range(self.freq_ghz.size):
            result[:, i] = ducc0.sht.synthesis_general(
                alm=self.alm[i],
                spin=0,
                lmax=self.config.lmax,
                loc=loc,
                epsilon=SYNTHESIS_EPSILON,
                nthreads=DUCC0_NTHREADS,
            )[0]
        return result


def orbital_doppler_mjy_sr(
    freq_ghz: np.ndarray,
    intensity_mjy_sr: np.ndarray,
    direction: np.ndarray,
    orbital_velocity_km_s: np.ndarray,
) -> np.ndarray:
    """First-order Doppler term of the Static sky, ΔI = (β·n̂)(3 I − dI/d ln ν) (ADR 0003)."""
    beta_dot_n = np.sum(orbital_velocity_km_s * direction, axis=1) / C_LIGHT_KM_S
    d_intensity = np.gradient(intensity_mjy_sr, np.log(freq_ghz), axis=1)
    return beta_dot_n[:, None] * (3 * intensity_mjy_sr - d_intensity)
