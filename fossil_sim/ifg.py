"""Ideal-FTS interferograms (systematics from v0.4).

For a DU with output sign s and band mask w(ν), the detected signal at
optical path difference x is

    P(x) = m Σ_ν w(ν) [I₁(ν) + I₂(ν) + s (I₁(ν) − I₂(ν)) cos(2π ν x / c)] Δν,

in MJy/sr × GHz (optical efficiency and étendue not included), with m the
``modulation_efficiency`` (provisional, see ``fossil_sim.noise``). The mirror
position is δ = x / ``opd_per_mirror_position``; the Step2 formula
cos(4π ν δ / c) corresponds to a factor 2 (TBD with the FTS team).

The batched functions (:func:`du_interferograms`, :func:`differential_spectra`)
serve the IFG-domain simulation (``ifg.simulate``): the sky of each IFG is
given on a short window of a fine grid of spin phase and interpolated to the
sample times (:func:`sample_weights`), so that it can change during the IFG.
Samples are equispaced in time, sample j at (j + ½) δt, with the OPD linear in
time from −x_max to +x_max (TBD — placeholder mirror trajectory).
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


def sample_weights(samples_per_ifg: int, grid_per_ifg: int) -> np.ndarray:
    """Cubic Lagrange weights from a sky-window grid to the IFG samples.

    The window of an IFG holds ``grid_per_ifg + 3`` grid points: the
    ``grid_per_ifg`` steps of the IFG (grid point g at phase g / grid_per_ifg
    of the IFG, in IFG units from its start) plus one point before and two
    after. Sample j, at phase (j + ½) / samples_per_ifg, is interpolated from
    its four neighbours. Shape (samples_per_ifg, grid_per_ifg + 3); each row
    sums to 1, so a constant sky passes through exactly.
    """
    u = (np.arange(samples_per_ifg) + 0.5) * grid_per_ifg / samples_per_ifg
    base = np.floor(u).astype(int)
    t = u - base
    weights = np.zeros((samples_per_ifg, grid_per_ifg + 3))
    rows = np.arange(samples_per_ifg)
    # window index of grid point g is g + 1; neighbours g = base - 1 .. base + 2
    for offset, w in enumerate(
        (
            -t * (t - 1) * (t - 2) / 6,
            (t + 1) * (t - 1) * (t - 2) / 2,
            -(t + 1) * t * (t - 2) / 2,
            (t + 1) * t * (t - 1) / 6,
        )
    ):
        weights[rows, base + offset] = w
    return weights


def du_interferograms(
    freq_ghz: np.ndarray,
    port1_mjy_sr: np.ndarray,
    port2_mjy_sr: np.ndarray,
    du: DetectorUnit,
    opd_m: np.ndarray,
    weights: np.ndarray,
    modulation_efficiency: float = 1.0,
) -> np.ndarray:
    """Interferograms of one DU for a block of IFGs, (n_ifg, n_samples).

    ``port1_mjy_sr`` and ``port2_mjy_sr`` are the Port spectra on the sky
    window of each IFG, (n_ifg, n_window, n_freq); ``weights`` maps the window
    to the samples (:func:`sample_weights`). Same units as :func:`interferogram`.
    """
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = float(freq_ghz[1] - freq_ghz[0])
    total = (port1_mjy_sr[..., mask] + port2_mjy_sr[..., mask]).sum(axis=-1)
    difference = port1_mjy_sr[..., mask] - port2_mjy_sr[..., mask]
    # modulated part on every (window point, sample), then interpolated
    modulated = difference @ modulation_matrix(freq_ghz[mask], opd_m).T
    ifgs = total @ weights.T + du.output_sign * np.einsum(
        "jw,iwj->ij", weights, modulated
    )
    return modulation_efficiency * freq_step_ghz * ifgs


def differential_spectra(
    freq_ghz: np.ndarray,
    ifgs: np.ndarray,
    du: DetectorUnit,
    opd_m: np.ndarray,
    modulation_efficiency: float = 1.0,
) -> np.ndarray:
    """Recover s (I₁ − I₂) on the DU band bins from interferograms (n_ifg, n_samples).

    Least squares on a constant plus the cosines of the band bins (no
    apodization), shape (n_ifg, n_band).
    """
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = float(freq_ghz[1] - freq_ghz[0])
    design = np.column_stack(
        [np.ones_like(opd_m), modulation_matrix(freq_ghz[mask], opd_m)]
    )
    solver = np.linalg.pinv(design)[1:] / (freq_step_ghz * modulation_efficiency)
    return ifgs @ solver.T


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
    result = np.full(freq_ghz.size, np.nan)
    result[du.band_mask(freq_ghz)] = differential_spectra(
        freq_ghz, ifg[None], du, opd_m, modulation_efficiency
    )[0]
    return result
