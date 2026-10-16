"""Rest-frame spectra of the Isotropic backgrounds, in MJy/sr.

Formulas follow bisou_sky (``generate_sky.py``) with explicit sources; all
functions accept frequency arrays of any shape.
"""

from importlib.resources import files

import numpy as np
from numba import vectorize

from .constants import C_LIGHT_M_S, H_PLANCK_J_S, K_BOLTZMANN_J_K, MJY_SR_IN_SI

ELECTRON_REST_ENERGY_KEV = 510.998_950  # CODATA 2018


def _x(freq_ghz, temperature_k):
    return H_PLANCK_J_S * np.asarray(freq_ghz) * 1e9 / (K_BOLTZMANN_J_K * temperature_k)


def _i0_mjy_sr(temperature_k: float) -> float:
    """2 (kT)³ / (hc)², the intensity scale of x-based spectra."""
    kt_j = K_BOLTZMANN_J_K * temperature_k
    return 2 * kt_j**3 / (H_PLANCK_J_S * C_LIGHT_M_S) ** 2 / MJY_SR_IN_SI


def _dbdt_shape(x):
    """x⁴ eˣ / (eˣ − 1)², the shape of ∂B/∂T (in units of I0 / T)."""
    return x**4 * np.exp(x) / np.expm1(x) ** 2


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
    """Non-relativistic Compton-y distortion, in MJy/sr.

    ΔI = y I0 x⁴ eˣ / (eˣ − 1)² [x coth(x/2) − 4], x = hν / kT
    (Zeldovich & Sunyaev 1969).
    """
    x = _x(freq_ghz, temperature_k)
    return y * _i0_mjy_sr(temperature_k) * _dbdt_shape(x) * (x / np.tanh(x / 2) - 4)


@vectorize(["float64(float64, float64)"], cache=True)
def _itoh_series(x, theta_e):
    """Σ_k Y_k θ_e^k of Itoh, Kohyama & Nozawa 1998, k = 1..4."""
    xt = x / np.tanh(0.5 * x)
    st = x / np.sinh(0.5 * x)
    y1 = -10.0 + 23.5 * xt - 8.4 * xt**2 + 0.7 * xt**3 + st**2 * (-4.2 + 1.4 * xt)
    y2 = (
        -7.5
        + 127.875 * xt
        - 173.6 * xt**2
        + 65.8 * xt**3
        - 8.8 * xt**4
        + 11.0 / 30.0 * xt**5
        + st**2 * (-86.8 + 131.6 * xt - 48.4 * xt**2 + 143.0 / 30.0 * xt**3)
        + st**4 * (-8.8 + 187.0 / 60.0 * xt)
    )
    y3 = (
        7.5
        + 313.125 * xt
        - 1419.6 * xt**2
        + 1425.3 * xt**3
        - 18594.0 / 35.0 * xt**4
        + 12059.0 / 140.0 * xt**5
        - 128.0 / 21.0 * xt**6
        + 16.0 / 105.0 * xt**7
        + st**2
        * (
            -709.8
            + 2850.6 * xt
            - 102267.0 / 35.0 * xt**2
            + 156767.0 / 140.0 * xt**3
            - 1216.0 / 7.0 * xt**4
            + 64.0 / 7.0 * xt**5
        )
        + st**4
        * (
            -18594.0 / 35.0
            + 205003.0 / 280.0 * xt
            - 1920.0 / 7.0 * xt**2
            + 1024.0 / 35.0 * xt**3
        )
        + st**6 * (-544.0 / 21.0 + 992.0 / 105.0 * xt)
    )
    y4 = (
        -4.21875
        + 237.3046875 * xt
        - 6239.1 * xt**2
        + 15368.175 * xt**3
        - 12438.9 * xt**4
        + 4446.2875 * xt**5
        - 16568.0 / 21.0 * xt**6
        + 7516.0 / 105.0 * xt**7
        - 22.0 / 7.0 * xt**8
        + 11.0 / 210.0 * xt**9
        + st**2
        * (
            -3119.55
            + 30736.35 * xt
            - 68413.95 * xt**2
            + 57801.7375 * xt**3
            - 157396.0 / 7.0 * xt**4
            + 30064.0 / 7.0 * xt**5
            - 2717.0 / 7.0 * xt**6
            + 2761.0 / 210.0 * xt**7
        )
        + st**4
        * (
            -12438.9
            + 37793.44375 * xt
            - 248520.0 / 7.0 * xt**2
            + 481024.0 / 35.0 * xt**3
            - 15972.0 / 7.0 * xt**4
            + 18689.0 / 140.0 * xt**5
        )
        + st**6
        * (
            -70414.0 / 21.0
            + 465992.0 / 105.0 * xt
            - 11792.0 / 7.0 * xt**2
            + 19778.0 / 105.0 * xt**3
        )
        + st**8 * (-682.0 / 7.0 + 7601.0 / 210.0 * xt)
    )
    return y1 * theta_e + y2 * theta_e**2 + y3 * theta_e**3 + y4 * theta_e**4


def relativistic_sz_correction_mjy_sr(
    freq_ghz: np.ndarray,
    y: float,
    electron_temperature_kev: float,
    temperature_k: float,
) -> np.ndarray:
    """Relativistic corrections to the y distortion, up to θ_e⁴, in MJy/sr.

    Itoh, Kohyama & Nozawa 1998 (ApJ 502, 7) expansion, with
    θ_e = kT_e / m_e c². Add to :func:`y_distortion_mjy_sr` for the total.
    """
    x = _x(freq_ghz, temperature_k)
    theta_e = electron_temperature_kev / ELECTRON_REST_ENERGY_KEV
    series = _itoh_series(x, theta_e)
    return y * _i0_mjy_sr(temperature_k) * _dbdt_shape(x) * series


# Zero of the μ-distortion shape, β_μ ≈ 2.1923 (Sunyaev & Zeldovich 1970)
MU_BETA = 2.1923


def mu_distortion_mjy_sr(
    freq_ghz: np.ndarray, mu: float, temperature_k: float
) -> np.ndarray:
    """Chemical-potential (μ) distortion, ΔI = μ I0 x⁴ eˣ / (eˣ − 1)² (1/β_μ − 1/x), in MJy/sr."""
    x = _x(freq_ghz, temperature_k)
    return mu * _i0_mjy_sr(temperature_k) * _dbdt_shape(x) * (1 / MU_BETA - 1 / x)


def cib_mjy_sr(
    freq_ghz: np.ndarray,
    amplitude: float,
    spectral_index: float,
    temperature_k: float,
    reference_freq_ghz: float,
) -> np.ndarray:
    """CIB monopole as a modified blackbody A (ν/ν₀)^β B_ν(T), in MJy/sr (Fixsen et al. 1998)."""
    freq_ghz = np.asarray(freq_ghz)
    shape = (freq_ghz / reference_freq_ghz) ** spectral_index
    return amplitude * shape * blackbody_mjy_sr(freq_ghz, temperature_k)


def _load_co_template() -> tuple[np.ndarray, np.ndarray]:
    path = files("fossil_sim").joinpath("data/extragalactic_co_template.txt")
    freq_ghz, intensity_mjy_sr = np.loadtxt(str(path), unpack=True)
    return np.log(freq_ghz), np.log(intensity_mjy_sr)


_CO_LOG_FREQ, _CO_LOG_INTENSITY = _load_co_template()


def extragalactic_co_mjy_sr(freq_ghz: np.ndarray, amplitude: float) -> np.ndarray:
    """Sky-averaged extragalactic CO (Abitbol et al. 2017), in MJy/sr.

    Log-log interpolation of the tabulated template; zero outside its range
    (10–615 GHz) instead of extrapolating.
    """
    log_freq = np.log(np.asarray(freq_ghz, dtype=np.float64))
    inside = (log_freq >= _CO_LOG_FREQ[0]) & (log_freq <= _CO_LOG_FREQ[-1])
    values = np.exp(np.interp(log_freq, _CO_LOG_FREQ, _CO_LOG_INTENSITY))
    return amplitude * np.where(inside, values, 0.0)
