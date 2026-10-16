"""Input sky.

v0.1.1: Isotropic backgrounds observed with the exact Boost (ADR 0002),

    I_obs(ν, n̂) = (ν / ν′)³ I_rest(ν′),   ν′ = ν γ (1 − β·n̂),

with β = (v_sun + v_orbit(t)) / c and n̂ the observed direction. The beam
acts on the anisotropic part only, through b₁ = exp(−σ²).

The Static sky (pysm3 anisotropies as alm, sampled with ducc0) follows the
design of litebird_sim's ``SkyGenerator`` and comes in v0.1.2.
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import BarycentricMeanEcliptic, SkyCoord

from .config import SkyConfig
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
