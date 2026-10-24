"""Ideal-FTS interferograms (illustration level; systematics from v0.4).

For a DU with output sign s and band mask w(ν), the detected signal at
optical path difference x is

    P(x) = m Σ_ν w(ν) [I₁(ν) + I₂(ν) + s (I₁(ν) − I₂(ν)) cos(2π ν x / c)] Δν,

in MJy/sr × GHz (optical efficiency and étendue not included), with m the
``modulation_efficiency`` (provisional, see ``fossil_sim.noise``). The mirror
position is δ = x / ``opd_per_mirror_position``; the Step2 formula
cos(4π ν δ / c) corresponds to a factor 2 (TBD with the FTS team).
"""

import numpy as np
from numba import njit

from .config import IFGConfig
from .constants import C_LIGHT_M_S
from .instrument import DetectorUnit


def opd_samples_m(config: IFGConfig) -> np.ndarray:
    """Symmetric scan of the optical path difference, from −x_max to +x_max."""
    return np.linspace(-config.max_opd_m, config.max_opd_m, config.samples_per_ifg)


def mirror_positions_m(config: IFGConfig) -> np.ndarray:
    """Mirror positions δ corresponding to :func:`opd_samples_m`."""
    return opd_samples_m(config) / config.opd_per_mirror_position


@njit(cache=True)
def _ifg_kernel(freq_ghz, total, difference, opd_m, freq_step_ghz, out):
    for j in range(opd_m.size):
        acc = 0.0
        for k in range(freq_ghz.size):
            phase = 2 * np.pi * freq_ghz[k] * 1e9 * opd_m[j] / C_LIGHT_M_S
            acc += total[k] + difference[k] * np.cos(phase)
        out[j] = acc * freq_step_ghz


def interferogram(
    freq_ghz: np.ndarray,
    port1_mjy_sr: np.ndarray,
    port2_mjy_sr: np.ndarray,
    du: DetectorUnit,
    opd_m: np.ndarray,
    modulation_efficiency: float = 1.0,
) -> np.ndarray:
    """Interferogram of one DU, in MJy/sr × GHz, for the spectra entering Ports 1 and 2."""
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = float(freq_ghz[1] - freq_ghz[0])
    total = np.where(mask, port1_mjy_sr + port2_mjy_sr, 0.0)
    difference = np.where(mask, du.output_sign * (port1_mjy_sr - port2_mjy_sr), 0.0)
    out = np.empty(opd_m.size)
    _ifg_kernel(freq_ghz, total, difference, opd_m, freq_step_ghz, out)
    return modulation_efficiency * out


def modulation_matrix(freq_ghz: np.ndarray, opd_m: np.ndarray) -> np.ndarray:
    """cos(2π ν x / c), shape (n_samples, n_freq)."""
    return np.cos(2 * np.pi * np.outer(opd_m, freq_ghz * 1e9) / C_LIGHT_M_S)


def differential_spectrum(
    freq_ghz: np.ndarray,
    ifg: np.ndarray,
    du: DetectorUnit,
    opd_m: np.ndarray,
    modulation_efficiency: float = 1.0,
) -> np.ndarray:
    """Recover s (I₁ − I₂) in the DU band from an interferogram, by least squares.

    Bins outside the band are returned as NaN.
    """
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = float(freq_ghz[1] - freq_ghz[0])
    design = np.column_stack(
        [
            np.ones_like(opd_m),
            modulation_matrix(freq_ghz[mask], opd_m),
        ]
    )
    scale = freq_step_ghz * modulation_efficiency
    solution = np.linalg.lstsq(design, ifg / scale, rcond=None)[0]
    result = np.full(freq_ghz.size, np.nan)
    result[mask] = solution[1:]
    return result
