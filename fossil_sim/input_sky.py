"""Input sky.

v0.1.1: Isotropic backgrounds observed with the exact Boost (ADR 0002),

    I_obs(ν, n̂) = (ν / ν′)³ I_rest(ν′),   ν′ = ν γ (1 − β·n̂),

with β = (v_sun + v_orbit(t)) / c and n̂ the observed direction. The beam
acts on the anisotropic part only, through b₁ = exp(−σ²).

v0.1.3: zodiacal light from zodipy, seen from the spacecraft position.

v0.1.2: the Static sky, following the design of litebird_sim's
``SkyGenerator``: pysm3 maps per frequency -> alm, Gaussian beam, rotation to
the J2000 ecliptic, cached on disk, sampled per IFG with
``ducc0.sht.synthesis_general``. Galactic models are defined in the
barycentric frame and receive the first-order orbital Doppler term; CMB
anisotropies and extragalactic models are defined in the CMB rest frame and
receive the solar + orbital one (ADR 0003).
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

from .config import EXTRAGALACTIC_STATIC_MODELS, Config, SkyConfig, StaticSkyConfig
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
        It multiplies all of ``observed − rest``, β² terms included (quadrupole
        should take b₂, the boost monopole no factor): error ~1e-10 relative.

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
# Half step in ln ν of the central difference giving dI/d ln ν from pysm3:
# for I ∝ ν^α the slope is off by α³ε²/6 · I (1.3e-4 I for α = 2)
DOPPLER_LOG_STEP = 1e-2

# (freq_ghz) -> (barycentric map, cosmological map), MJy/sr, Galactic coordinates
EmissionFunction = Callable[[float], tuple[np.ndarray, np.ndarray]]


def is_cosmological_model(model: str) -> bool:
    """CMB anisotropies (``c1``…``c7``) and extragalactic models: CMB rest frame.

    >>> [is_cosmological_model(m) for m in ("c1", "cib1", "tsz1", "d2", "s1")]
    [True, True, True, False, False]
    """
    cmb = model[0] == "c" and model[1:].isdigit()
    return cmb or model.startswith(EXTRAGALACTIC_STATIC_MODELS)


def pysm3_emission(models: tuple[str, ...], nside: int) -> EmissionFunction:
    """pysm3 total intensity of the barycentric (Galactic) and cosmological models."""
    import pysm3
    import pysm3.units as pysm_units

    def sky(names):
        if not names:
            return None
        return pysm3.Sky(nside=nside, preset_strings=list(names), output_unit="MJy/sr")

    barycentric = sky([m for m in models if not is_cosmological_model(m)])
    cosmological = sky([m for m in models if is_cosmological_model(m)])
    zero = np.zeros(hp.nside2npix(nside))

    def emission(freq_ghz: float) -> tuple[np.ndarray, np.ndarray]:
        freq = freq_ghz * pysm_units.GHz
        return tuple(
            zero if s is None else np.asarray(s.get_emission(freq)[0].value)
            for s in (barycentric, cosmological)
        )

    return emission


def static_sky_cache_key(
    config: StaticSkyConfig, freq_ghz: np.ndarray, fwhm_deg: float
) -> str:
    """Hash of everything that determines the Static-sky cubes."""
    import pysm3

    content = {
        "format": 2,
        "models": list(config.models),
        "nside_eval": 2 * config.nside,
        "lmax": config.lmax,
        "fwhm_deg": fwhm_deg,
        "freq_ghz": [float(f) for f in freq_ghz],
        "doppler_log_step": DOPPLER_LOG_STEP if config.doppler else None,
        "pysm3": pysm3.__version__,
    }
    return hashlib.sha256(json.dumps(content, sort_keys=True).encode()).hexdigest()


def build_static_sky_alm(
    config: StaticSkyConfig,
    freq_ghz: np.ndarray,
    fwhm_deg: float,
    emission: EmissionFunction | None = None,
) -> tuple[np.ndarray, np.ndarray | None]:
    """Static-sky alm cubes, MJy/sr, J2000 ecliptic.

    For each frequency: maps -> alm (3 iterations), monopole of the
    cosmological part removed, Gaussian beam, rotation Galactic -> ecliptic.

    Returns
    -------
    alm : ndarray, shape (n_freq, 1, n_alm)
        Total intensity I.
    doppler_alm : ndarray, shape (2, n_freq, 1, n_alm), or None
        K = 3 I − dI/d ln ν of the barycentric and of the cosmological part
        (ADR 0003), with dI/d ln ν from pysm3 at ν e^{±ε}; None when
        ``config.doppler`` is false.
    """
    if emission is None:
        emission = pysm3_emission(config.models, 2 * config.nside)
    lmax = config.lmax
    beam = hp.gauss_beam(np.deg2rad(fwhm_deg), lmax=lmax)
    rotator = hp.Rotator(coord=["G", "E"])

    def to_alm(barycentric_map, cosmological_map):
        alm = hp.map2alm(barycentric_map, lmax=lmax, iter=3)
        cosmological = hp.map2alm(cosmological_map, lmax=lmax, iter=3)
        cosmological[0] = 0.0  # a_00: the monopole belongs to the Isotropic backgrounds
        return [
            rotator.rotate_alm(hp.almxfl(a, beam), lmax=lmax)
            for a in (alm, cosmological)
        ]

    n_alm = hp.Alm.getsize(lmax)
    cube = np.empty((freq_ghz.size, 1, n_alm), dtype=np.complex128)
    doppler = (
        np.empty((2, freq_ghz.size, 1, n_alm), dtype=np.complex128)
        if config.doppler
        else None
    )
    eps = DOPPLER_LOG_STEP
    for i, freq in enumerate(freq_ghz):
        maps = emission(float(freq))
        barycentric, cosmological = to_alm(*maps)
        cube[i, 0] = barycentric + cosmological
        if doppler is not None:
            up = emission(float(freq) * np.exp(eps))
            down = emission(float(freq) * np.exp(-eps))
            k_maps = [3 * m - (u - d) / (2 * eps) for m, u, d in zip(maps, up, down)]
            doppler[:, i, 0] = to_alm(*k_maps)
    return cube, doppler


class StaticSky:
    """Static sky cubes, loaded from (or written to) the cache, sampled per direction.

    Parameters
    ----------
    config : StaticSkyConfig
    freq_ghz : ndarray
        Frequency grid.
    fwhm_deg : float
        Gaussian beam FWHM.
    emission : callable, optional
        Map generator, (barycentric, cosmological) per frequency (default:
        pysm3 with ``config.models``).
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
                self.doppler_alm = f["doppler_alm"][()] if "doppler_alm" in f else None
        else:
            self.alm, self.doppler_alm = build_static_sky_alm(
                config, self.freq_ghz, fwhm_deg, emission
            )
            self.path.parent.mkdir(parents=True, exist_ok=True)
            tmp = self.path.with_suffix(".tmp")
            with h5py.File(tmp, "w") as f:
                f.create_dataset("alm", data=self.alm)
                if self.doppler_alm is not None:
                    f.create_dataset("doppler_alm", data=self.doppler_alm)
                f.create_dataset("freq_ghz", data=self.freq_ghz)
                f.attrs["cache_key"] = key
                f.attrs["unit"] = "MJy/sr"
                f.attrs["coord"] = "J2000 ecliptic"
            tmp.replace(self.path)
        self.maps = self.doppler_maps = None
        if config.sampling == "nearest":
            self.maps = self._alm2map(self.alm)
            if self.doppler_alm is not None:
                self.doppler_maps = np.stack(
                    [self._alm2map(c) for c in self.doppler_alm]
                )

    def _alm2map(self, cube: np.ndarray) -> np.ndarray:
        nside, lmax = self.config.nside, self.config.lmax
        return np.stack([hp.alm2map(a[0], nside, lmax=lmax) for a in cube])

    def _synthesize(self, cube, maps, direction: np.ndarray) -> np.ndarray:
        """One cube along each direction, shape (n, n_freq)."""
        theta, phi = hp.vec2ang(direction)
        if maps is not None:
            pixels = hp.ang2pix(self.config.nside, theta, phi)
            return maps[:, pixels].T.copy()
        loc = np.column_stack([theta, phi])
        result = np.empty((direction.shape[0], self.freq_ghz.size))
        for i in range(self.freq_ghz.size):
            result[:, i] = ducc0.sht.synthesis_general(
                alm=cube[i],
                spin=0,
                lmax=self.config.lmax,
                loc=loc,
                epsilon=SYNTHESIS_EPSILON,
                nthreads=DUCC0_NTHREADS,
            )[0]
        return result

    def sample(self, direction: np.ndarray) -> np.ndarray:
        """Static sky along each direction, shape (n, n_freq), MJy/sr."""
        return self._synthesize(self.alm, self.maps, direction)

    def sample_doppler(self, direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """K = 3 I − dI/d ln ν of the barycentric and cosmological parts, each (n, n_freq)."""
        if self.doppler_alm is None:
            raise ValueError("Static sky built with doppler = false")
        maps = [None, None] if self.doppler_maps is None else self.doppler_maps
        return tuple(
            self._synthesize(cube, m, direction)
            for cube, m in zip(self.doppler_alm, maps)
        )


def static_doppler_mjy_sr(
    k_mjy_sr: np.ndarray, direction: np.ndarray, velocity_km_s: np.ndarray
) -> np.ndarray:
    """First-order Doppler of a Static-sky part, ΔI = (β·n̂) K (ADR 0003).

    ``k_mjy_sr`` is K = 3 I − dI/d ln ν of that part (:meth:`StaticSky.sample_doppler`),
    ``velocity_km_s`` the observer velocity in its frame, (n, 3) or (3,).
    """
    velocity_km_s = np.broadcast_to(velocity_km_s, direction.shape)
    beta_dot_n = np.sum(velocity_km_s * direction, axis=1) / C_LIGHT_KM_S
    return beta_dot_n[:, None] * k_mjy_sr


class ZodiacalLight:
    """Zodiacal light per IFG, from zodipy, evaluated Ring by Ring.

    The Scans of a Ring retrace the same circle and only the observer moves
    (about 0.4° in 9.2 h). The emission is evaluated on the ``ifgs_per_scan``
    directions of each Ring at its start and end (for the last Ring, the end
    of the Recycling window) and linearly interpolated in time for each IFG.
    The observer is the spacecraft position of :class:`~fossil_sim.orbit.L2Orbit`.
    No beam smoothing and no Doppler effect from the dust motion are applied.
    """

    def __init__(self, config: Config, orbit, freq_ghz: np.ndarray) -> None:
        import zodipy

        self.config = config
        self.orbit = orbit
        self.freq_ghz = np.asarray(freq_ghz)
        zc = config.zodi
        self.models = [
            zodipy.Model(
                f * u.GHz,
                name=zc.model,
                extrapolate=zc.extrapolate,
                gauss_quad_degree=zc.gauss_quad_degree,
            )
            for f in self.freq_ghz
        ]

    def evaluate(self, direction: np.ndarray, time_s: np.ndarray) -> np.ndarray:
        """Zodiacal light along directions at the given times, (n, n_freq), MJy/sr."""
        from astropy.time import TimeDelta

        obstime = self.config.simulation.mission_start_time + TimeDelta(
            time_s, format="sec"
        )
        lon_rad = np.arctan2(direction[:, 1], direction[:, 0])
        lat_rad = np.arcsin(np.clip(direction[:, 2], -1.0, 1.0))
        coords = SkyCoord(
            lon_rad * u.rad,
            lat_rad * u.rad,
            frame=BarycentricMeanEcliptic(equinox="J2000"),
            obstime=obstime,
        )
        obspos = (self.orbit.heliocentric_position_km(time_s).T * u.km).to(u.AU)
        result = np.empty((direction.shape[0], self.freq_ghz.size))
        for i, model in enumerate(self.models):
            emission = model.evaluate(
                coords, obspos=obspos, nprocesses=self.config.zodi.nprocesses
            )
            result[:, i] = emission.to_value(u.MJy / u.sr)
        return result

    def for_cycle(
        self, cycle_id: int, ring_axes: np.ndarray, table: np.ndarray
    ) -> np.ndarray:
        """Zodiacal light of the IFGs of one cycle (``table`` rows), (n, n_freq)."""
        from .scanning import ring_circles

        tl = self.config.timeline
        n_rings, n_phase = ring_axes.shape[0], tl.ifgs_per_scan
        ring_duration_s = tl.scans_per_ring * tl.scan_duration_s
        start_s = cycle_id * tl.cycle_duration_s + np.arange(n_rings) * ring_duration_s
        end_s = start_s + ring_duration_s
        end_s[-1] = (cycle_id + 1) * tl.cycle_duration_s  # last Ring covers Recycling
        circles = ring_circles(self.config, ring_axes)  # (n_rings, n_phase, 3)
        directions = np.concatenate([circles, circles], axis=0).reshape(-1, 3)
        times_s = np.repeat(np.concatenate([start_s, end_s]), n_phase)
        values = self.evaluate(directions, times_s).reshape(2, n_rings, n_phase, -1)

        ring = table["ring_in_cycle"].astype(np.int64)
        phase = table["ifg_in_scan"].astype(np.int64)
        weight = (table["time_mid_s"] - start_s[ring]) / (end_s[ring] - start_s[ring])
        first, last = values[0, ring, phase], values[1, ring, phase]
        return first + weight[:, None] * (last - first)
