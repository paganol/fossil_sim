"""Ideal-FTS interferograms (illustration level; systematics from v0.4).

For a DU with output sign s and band mask w(ν), the detected signal at
mirror position δ is (Step2)

    P(δ) = Σ_ν w(ν) [I₁(ν) + I₂(ν) + s (I₁(ν) − I₂(ν)) cos(4π ν δ / c)] Δν,

in MJy/sr × GHz (optical efficiency and étendue not included).
"""

import numpy as np
from numba import njit

from .config import IFGConfig
from .constants import C_LIGHT_M_S
from .instrument import DetectorUnit


def mirror_positions_m(config: IFGConfig) -> np.ndarray:
    """Symmetric mirror scan from −δ_max to +δ_max."""
    return np.linspace(
        -config.max_mirror_position_m,
        config.max_mirror_position_m,
        config.samples_per_ifg,
    )


@njit(cache=True)
def _ifg_kernel(freq_ghz, total, difference, mirror_position_m, freq_step_ghz, out):
    for j in range(mirror_position_m.size):
        acc = 0.0
        for k in range(freq_ghz.size):
            phase = 4 * np.pi * freq_ghz[k] * 1e9 * mirror_position_m[j] / C_LIGHT_M_S
            acc += total[k] + difference[k] * np.cos(phase)
        out[j] = acc * freq_step_ghz


def interferogram(
    freq_ghz: np.ndarray,
    port1_mjy_sr: np.ndarray,
    port2_mjy_sr: np.ndarray,
    du: DetectorUnit,
    mirror_position_m: np.ndarray,
) -> np.ndarray:
    """Interferogram of one DU, in MJy/sr × GHz, for the spectra entering Ports 1 and 2."""
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = float(freq_ghz[1] - freq_ghz[0])
    total = np.where(mask, port1_mjy_sr + port2_mjy_sr, 0.0)
    difference = np.where(mask, du.output_sign * (port1_mjy_sr - port2_mjy_sr), 0.0)
    out = np.empty(mirror_position_m.size)
    _ifg_kernel(freq_ghz, total, difference, mirror_position_m, freq_step_ghz, out)
    return out


def modulation_matrix(
    freq_ghz: np.ndarray, mirror_position_m: np.ndarray
) -> np.ndarray:
    """cos(4π ν δ / c), shape (n_samples, n_freq)."""
    return np.cos(4 * np.pi * np.outer(mirror_position_m, freq_ghz * 1e9) / C_LIGHT_M_S)


def differential_spectrum(
    freq_ghz: np.ndarray,
    ifg: np.ndarray,
    du: DetectorUnit,
    mirror_position_m: np.ndarray,
) -> np.ndarray:
    """Recover s (I₁ − I₂) in the DU band from an interferogram, by least squares.

    Bins outside the band are returned as NaN.
    """
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = float(freq_ghz[1] - freq_ghz[0])
    design = np.column_stack(
        [
            np.ones_like(mirror_position_m),
            modulation_matrix(freq_ghz[mask], mirror_position_m),
        ]
    )
    solution = np.linalg.lstsq(design, ifg / freq_step_ghz, rcond=None)[0]
    result = np.full(freq_ghz.size, np.nan)
    result[mask] = solution[1:]
    return result
