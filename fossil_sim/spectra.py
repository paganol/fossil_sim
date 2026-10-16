"""Analytic spectra in MJy/sr."""

import numpy as np

from .constants import C_LIGHT_M_S, H_PLANCK_J_S, K_BOLTZMANN_J_K, MJY_SR_IN_SI


def _x(freq_ghz, temperature_k):
    return H_PLANCK_J_S * np.asarray(freq_ghz) * 1e9 / (K_BOLTZMANN_J_K * temperature_k)


def blackbody_mjy_sr(freq_ghz: np.ndarray, temperature_k: float) -> np.ndarray:
    """Planck brightness B_ν(T) in MJy/sr.

    >>> round(float(blackbody_mjy_sr(160.0, 2.72548)), 1)  # CMB peak region
    383.9
    """
    freq_hz = np.asarray(freq_ghz) * 1e9
    prefactor = 2 * H_PLANCK_J_S * freq_hz**3 / C_LIGHT_M_S**2
    return prefactor / np.expm1(_x(freq_ghz, temperature_k)) / MJY_SR_IN_SI


def y_distortion_mjy_sr(
    freq_ghz: np.ndarray, y: float, temperature_k: float
) -> np.ndarray:
    """Non-relativistic Compton-y distortion of a blackbody of temperature T, in MJy/sr.

    ΔI = y B_ν(T) x eˣ / (eˣ − 1) [x coth(x/2) − 4], with x = hν / kT
    (Zeldovich & Sunyaev 1969).
    """
    x = _x(freq_ghz, temperature_k)
    shape = x * np.exp(x) / np.expm1(x) * (x / np.tanh(x / 2) - 4)
    return y * blackbody_mjy_sr(freq_ghz, temperature_k) * shape
