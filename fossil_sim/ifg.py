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
time between the first and last sample, at −x_max and +x_max
(:func:`opd_samples_m`): the mirror turns around in the half sample at each end
of the IFG (TBD — placeholder mirror trajectory).
"""

import numpy as np
from numba import njit, prange

from .config import IFGConfig
from .constants import C_LIGHT_M_S
from .instrument import DetectorUnit


def opd_samples_m(config: IFGConfig) -> np.ndarray:
    """Symmetric scan of the optical path difference, from −x_max to +x_max."""
    return np.linspace(-config.max_opd_m, config.max_opd_m, config.samples_per_ifg)


def mirror_positions_m(config: IFGConfig) -> np.ndarray:
    """Mirror positions δ corresponding to :func:`opd_samples_m`."""
    return opd_samples_m(config) / config.opd_per_mirror_position


def _frequency_step_ghz(freq_ghz: np.ndarray, freq_step_ghz: float | None) -> float:
    """Use the supplied bin width, or infer it from a grid with at least two bins.

    With two or more bins the supplied width must equal the grid spacing.
    """
    spacing = np.diff(freq_ghz)
    if freq_step_ghz is None:
        if freq_ghz.size < 2:
            raise ValueError(
                "freq_step_ghz is required for a single-bin frequency grid"
            )
        freq_step_ghz = float(spacing[0])
    if not freq_step_ghz > 0:  # also NaN
        raise ValueError("freq_step_ghz must be > 0")
    if spacing.size and not np.allclose(spacing, freq_step_ghz, rtol=1e-9, atol=0):
        raise ValueError("freq_step_ghz must equal the spacing of a uniform grid")
    return freq_step_ghz


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
    *,
    freq_step_ghz: float | None = None,
) -> np.ndarray:
    """Interferogram of one DU, in MJy/sr × GHz, for the spectra entering Ports 1 and 2.

    ``freq_step_ghz`` defaults to the grid spacing; supply it for a single bin.
    """
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = _frequency_step_ghz(freq_ghz, freq_step_ghz)
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
    *,
    freq_step_ghz: float | None = None,
    gain: np.ndarray | None = None,
    offset_mjy_sr: np.ndarray | None = None,
    transmission: np.ndarray | None = None,
    constant_window: np.ndarray | None = None,
    modulation: np.ndarray | None = None,
) -> np.ndarray:
    """Interferograms of one DU for a block of IFGs, (n_ifg, n_samples).

    ``port1_mjy_sr`` and ``port2_mjy_sr`` are the Port spectra on the sky
    window of each IFG, (n_ifg, n_window, n_freq); ``weights`` maps the window
    to the samples (:func:`sample_weights`). Same units as :func:`interferogram`.
    ``freq_step_ghz`` defaults to the grid spacing; supply it for a single bin.
    ``transmission`` (on the whole grid, default: the top-hat
    :meth:`~fossil_sim.instrument.DetectorUnit.band_mask`) weights both terms;
    the band bins are those where it is positive. ``gain`` and
    ``offset_mjy_sr`` (on the band bins) give the DU response of
    :mod:`fossil_sim.response`: modulated amplitude g (s (I₁ − I₂) + o), and
    g also on the unmodulated term. ``constant_window`` flags the IFGs whose
    Port spectra do not change over the window (detected when not given).
    ``modulation`` is :func:`modulation_matrix` on the band bins (computed
    when not given): it depends only on the configuration, pass it to reuse it.
    """
    mask, factor = _band_factor(freq_ghz, du, gain, transmission)
    freq_step_ghz = _frequency_step_ghz(freq_ghz, freq_step_ghz)
    band = _band_index(mask)
    port1, port2 = port1_mjy_sr[..., band], port2_mjy_sr[..., band]
    total = port1 @ factor + port2 @ factor
    difference = port1 - port2  # (n_ifg, n_window, n_band)
    if offset_mjy_sr is not None:
        difference += du.output_sign * offset_mjy_sr
    difference *= factor
    modulated = _interpolated_modulation(
        difference,
        _band_modulation(modulation, freq_ghz, mask, opd_m),
        weights,
        constant_window,
    )
    ifgs = total @ weights.T + du.output_sign * modulated
    return modulation_efficiency * freq_step_ghz * ifgs


def static_basis_interferograms(
    freq_ghz: np.ndarray,
    spectra: np.ndarray,
    du: DetectorUnit,
    opd_m: np.ndarray,
    *,
    gain: np.ndarray | None = None,
    transmission: np.ndarray | None = None,
    modulation: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """IFG terms of each basis spectrum U_k(ν) of the Static sky, for one DU.

    ``spectra`` is (n_freq, K). Returns the modulated term Σ_ν cos(2πνx/c)
    f(ν) U_k(ν), (n_samples, K), and the unmodulated one Σ_ν f(ν) U_k(ν),
    (K,), with f the transmission times the gain (as in
    :func:`du_interferograms`, without m Δν). Computed once per cycle, they
    replace the sum over the frequencies in :func:`du_interferograms_low_rank`.
    """
    mask, factor = _band_factor(freq_ghz, du, gain, transmission)
    band = spectra[mask] * factor[:, None]
    modulation = _band_modulation(modulation, freq_ghz, mask, opd_m)
    return modulation @ band, band.sum(axis=0)


def du_interferograms_low_rank(
    freq_ghz: np.ndarray,
    port1_mjy_sr: np.ndarray,
    port2_mjy_sr: np.ndarray,
    static_window: np.ndarray,
    sky_sign: np.ndarray,
    sky_ports: np.ndarray,
    basis: tuple[np.ndarray, np.ndarray],
    du: DetectorUnit,
    opd_m: np.ndarray,
    weights: np.ndarray,
    modulation_efficiency: float = 1.0,
    *,
    freq_step_ghz: float | None = None,
    gain: np.ndarray | None = None,
    offset_mjy_sr: np.ndarray | None = None,
    transmission: np.ndarray | None = None,
    modulation: np.ndarray | None = None,
) -> np.ndarray:
    """Interferograms of one DU, (n_ifg, n_samples), Static sky in a low-rank basis.

    The Port spectra are ``port1_mjy_sr`` and ``port2_mjy_sr`` (n_ifg,
    n_freq), constant over the IFG (sky without the Static sky, or BBIR),
    plus on each Port that sees the sky the Static sky Σ_k a_k U_k(ν), with
    the amplitudes a_k on the window of each IFG, ``static_window`` (n_ifg,
    n_window, K). ``sky_sign`` is (Port 1 on sky) − (Port 2 on
    sky) and ``sky_ports`` their sum; ``basis`` comes from
    :func:`static_basis_interferograms` with the same gain and transmission.
    Equal to :func:`du_interferograms` on the full Port spectra, at a cost per
    sample proportional to K instead of the number of band frequencies.
    ``modulation`` as in :func:`du_interferograms`.
    """
    mask, factor = _band_factor(freq_ghz, du, gain, transmission)
    freq_step_ghz = _frequency_step_ghz(freq_ghz, freq_step_ghz)
    band = _band_index(mask)
    port1, port2 = port1_mjy_sr[:, band], port2_mjy_sr[:, band]
    total = port1 @ factor + port2 @ factor
    difference = port1 - port2
    if offset_mjy_sr is not None:
        difference += du.output_sign * offset_mjy_sr
    difference *= factor
    modulation = _band_modulation(modulation, freq_ghz, mask, opd_m)
    modulated = difference @ modulation.T
    # Static sky, only on the IFGs with a Port on the sky
    rows = np.flatnonzero(sky_ports)
    sky_row = np.full(len(total), -1)
    sky_row[rows] = np.arange(rows.size)
    window = static_window if rows.size == len(total) else static_window[rows]
    basis_modulated, basis_unmodulated = basis
    varying = _interpolated_modulation(
        window, basis_modulated, weights, np.zeros(rows.size, dtype=bool)
    )
    unmodulated = window @ basis_unmodulated  # (n_rows, n_window)
    ifgs = np.empty((len(total), weights.shape[0]))
    _low_rank_kernel(
        total,
        modulated,
        weights,
        sky_row,
        sky_ports,
        sky_sign,
        unmodulated,
        varying,
        float(du.output_sign),
        modulation_efficiency * freq_step_ghz,
        ifgs,
    )
    return ifgs


@njit(parallel=True, cache=True)
def _low_rank_kernel(
    total,
    modulated,
    weights,
    sky_row,
    sky_ports,
    sky_sign,
    unmodulated,
    varying,
    sign,
    scale,
    out,
):
    """Interferograms of :func:`du_interferograms_low_rank`, one IFG per thread.

    out[i, j] = scale (W_j (total_i + s M_ij) + p_i Σ_w W_jw u_wr + s σ_i V_rj),
    with W_j the window sum, r = sky_row[i] (−1: no Port on the sky) and u,
    V the unmodulated and modulated Static sky of row r.
    """
    n_ifg, n_samples = out.shape
    n_window = weights.shape[1]
    window_sum = np.zeros(n_samples)
    first = np.zeros(n_samples, dtype=np.int64)  # non-zero weights: first:stop
    stop = np.zeros(n_samples, dtype=np.int64)
    for j in range(n_samples):
        first[j] = n_window
        for w in range(n_window):
            window_sum[j] += weights[j, w]
            if weights[j, w] != 0.0:
                first[j] = min(first[j], w)
                stop[j] = w + 1
    for i in prange(n_ifg):
        r = sky_row[i]
        for j in range(n_samples):
            value = window_sum[j] * (total[i] + sign * modulated[i, j])
            if r >= 0:
                static = 0.0
                for w in range(first[j], stop[j]):
                    static += weights[j, w] * unmodulated[r, w]
                value += sky_ports[i] * static + sign * sky_sign[i] * varying[r, j]
            out[i, j] = scale * value


def _band_modulation(
    modulation: np.ndarray | None,
    freq_ghz: np.ndarray,
    mask: np.ndarray,
    opd_m: np.ndarray,
) -> np.ndarray:
    """``modulation`` if given, else :func:`modulation_matrix` on the band bins."""
    if modulation is None:
        return modulation_matrix(freq_ghz[mask], opd_m)
    return modulation


def _band_factor(
    freq_ghz: np.ndarray,
    du: DetectorUnit,
    gain: np.ndarray | None,
    transmission: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Band bins of a DU and the factor on them: transmission (or 1) times gain."""
    if transmission is None:
        mask = du.band_mask(freq_ghz)
        factor = np.ones(np.count_nonzero(mask))
    else:
        mask = transmission > 0
        factor = transmission[mask]
    if gain is not None:
        factor = factor * gain
    return mask, factor


def _band_index(mask: np.ndarray) -> slice | np.ndarray:
    """``mask`` as a slice when its True entries are contiguous (a view, not a copy)."""
    index = np.flatnonzero(mask)
    if index.size and index[-1] - index[0] + 1 == index.size:
        return slice(index[0], index[-1] + 1)
    return mask


def _window_groups(weights: np.ndarray) -> list[tuple[slice | np.ndarray, slice]]:
    """Samples grouped by the window points they use: (samples, window columns).

    Each row of ``weights`` is non-zero on a few consecutive columns (four
    for :func:`sample_weights`); samples starting on the same column form a
    group, with the columns spanned by its rows.
    """
    nonzero = weights != 0
    first = nonzero.argmax(axis=1)
    last = weights.shape[1] - 1 - nonzero[:, ::-1].argmax(axis=1)
    groups = []
    for start in np.unique(first):
        samples = np.flatnonzero(first == start)
        if np.array_equal(samples, np.arange(samples[0], samples[-1] + 1)):
            samples = slice(samples[0], samples[-1] + 1)
        groups.append((samples, slice(start, last[first == start].max() + 1)))
    return groups


def _interpolated_modulation(
    difference: np.ndarray,
    modulation: np.ndarray,
    weights: np.ndarray,
    constant_window: np.ndarray | None = None,
) -> np.ndarray:
    """Σ_w weights[j, w] Σ_k difference[i, w, k] modulation[j, k], (n_ifg, n_samples).

    ``difference`` is (n_ifg, n_window, n_freq). Equal, to rounding, to
    ``einsum("jw,iwk,jk->ij", weights, difference, modulation)``, without the
    products on window points a sample does not use: IFGs whose window is
    constant (``constant_window``, detected when not given) take one product,
    the others one per group of samples (:func:`_window_groups`), on its few
    window points with the weights folded into the modulation.
    """
    n_ifg = difference.shape[0]
    n_samples = modulation.shape[0]
    if constant_window is None:
        constant_window = np.all(difference == difference[:, :1], axis=(1, 2))
    out = np.empty((n_ifg, n_samples))
    if constant_window.any():
        out[constant_window] = (
            difference[constant_window, 0] @ modulation.T
        ) * weights.sum(axis=1)
    varying = np.flatnonzero(~constant_window)
    if not varying.size:
        return out
    rows = difference if varying.size == n_ifg else difference[varying]
    # one product per group, (IFGs, window points × frequencies) @ B with
    # the weights folded into the modulation, B[(w, k), j] = W_jw C_jk
    result = np.empty((varying.size, n_samples))
    for samples, columns in _window_groups(weights):
        folded = (
            weights[samples, columns].T[:, :, None] * modulation[samples][None]
        ).transpose(0, 2, 1)  # (window points, frequencies, samples)
        result[:, samples] = rows[:, columns].reshape(varying.size, -1) @ (
            folded.reshape(-1, folded.shape[2])
        )
    out[varying] = result
    return out


# Norton & Beer (1976) apodizations, A(u) = Σ_i c_i (1 − u²)^i with Σ c_i = 1
NORTON_BEER_COEFFICIENTS = {
    "norton_beer_weak": (0.384093, -0.087577, 0.703484),
    "norton_beer_medium": (0.152442, -0.136176, 0.983734),
    "norton_beer_strong": (0.045335, 0.0, 0.554883, 0.0, 0.399782),
}


# Kaiser β when not given (FWHM ≈ 1.96 c / (2 x_max), sidelobes 0.65 %)
DEFAULT_KAISER_BETA = 6.0


def apodization_window(
    name: str, opd_m: np.ndarray, kaiser_beta: float = DEFAULT_KAISER_BETA
) -> np.ndarray:
    """Apodization A(x) on the OPD samples, A(0) = 1, u = |x| / max|x|.

    ``"none"`` (A = 1), the Norton–Beer weak, medium and strong functions,
    ``"hann"``, cos²(π u / 2), ``"triangular"``, 1 − u (line shape sinc²,
    FWHM 1.772 c / (2 x_max), sidelobes positive), or ``"kaiser"``,
    I₀(β √(1 − u²)) / I₀(β): β = 0 is no apodization, a larger β a wider
    line shape with lower sidelobes. A(0) = 1 keeps the level of a smooth
    spectrum: the instrumental line shape integrates to 1.
    """
    u = np.abs(opd_m) / np.max(np.abs(opd_m))
    if name == "none":
        return np.ones_like(u)
    if name == "hann":
        return np.cos(np.pi * u / 2) ** 2
    if name == "triangular":
        return 1 - u
    if name == "kaiser":
        return np.i0(kaiser_beta * np.sqrt(1 - u**2)) / np.i0(kaiser_beta)
    if name in NORTON_BEER_COEFFICIENTS:
        return sum(
            c * (1 - u**2) ** i for i, c in enumerate(NORTON_BEER_COEFFICIENTS[name])
        )
    raise ValueError(f"unknown apodization {name!r}")


def differential_spectrum_solver(
    freq_ghz: np.ndarray,
    du: DetectorUnit,
    opd_m: np.ndarray,
    modulation_efficiency: float = 1.0,
    *,
    freq_step_ghz: float | None = None,
    apodization: str = "none",
    kaiser_beta: float = DEFAULT_KAISER_BETA,
) -> np.ndarray:
    """Linear map from IFG samples to the DU band spectra, (n_band, n_samples).

    ``ifgs @ solver.T`` is :func:`differential_spectra`; the same matrix
    propagates sample noise (:func:`fossil_sim.noise.sky_estimate_noise_variance`).
    It depends only on the configuration: compute it once per DU and reuse it.

    Without apodization, the least squares on a constant plus the cosines of
    the band bins, D⁺. With apodization A (:func:`apodization_window`), the
    constant a₀ of that fit is removed first, then the apodized samples are
    projected on the cosines: D⁺[1:] diag(A) (1 − 𝟙 D⁺[0]). Each bin is then
    the spectrum seen through the line shape of A (wider, lower sidelobes)
    and neighbouring bins are correlated; the rows still annihilate a
    constant.
    """
    mask = du.band_mask(freq_ghz)
    freq_step_ghz = _frequency_step_ghz(freq_ghz, freq_step_ghz)
    design = np.column_stack(
        [np.ones_like(opd_m), modulation_matrix(freq_ghz[mask], opd_m)]
    )
    pinv = np.linalg.pinv(design)
    fit = pinv[1:]
    if apodization != "none":
        window = apodization_window(apodization, opd_m, kaiser_beta)
        fit = fit * window - np.outer(fit @ window, pinv[0])
    return fit / (freq_step_ghz * modulation_efficiency)


def differential_spectra(
    freq_ghz: np.ndarray,
    ifgs: np.ndarray,
    du: DetectorUnit,
    opd_m: np.ndarray,
    modulation_efficiency: float = 1.0,
    *,
    freq_step_ghz: float | None = None,
    apodization: str = "none",
    kaiser_beta: float = DEFAULT_KAISER_BETA,
) -> np.ndarray:
    """Recover s (I₁ − I₂) on the DU band bins from interferograms (n_ifg, n_samples).

    Least squares on a constant plus the cosines of the band bins, with
    optional ``apodization`` (:func:`differential_spectrum_solver`; β of
    ``"kaiser"`` in ``kaiser_beta``), shape
    (n_ifg, n_band). ``freq_step_ghz`` defaults to the grid spacing; supply
    it for a single bin.
    """
    solver = differential_spectrum_solver(
        freq_ghz,
        du,
        opd_m,
        modulation_efficiency,
        freq_step_ghz=freq_step_ghz,
        apodization=apodization,
        kaiser_beta=kaiser_beta,
    )
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
    *,
    freq_step_ghz: float | None = None,
    apodization: str = "none",
    kaiser_beta: float = DEFAULT_KAISER_BETA,
) -> np.ndarray:
    """Recover s (I₁ − I₂) in the DU band from an interferogram, by least squares.

    Bins outside the band are returned as NaN.
    ``freq_step_ghz`` defaults to the grid spacing; supply it for a single bin.
    ``apodization`` and ``kaiser_beta`` as in :func:`apodization_window`.
    """
    result = np.full(freq_ghz.size, np.nan)
    result[du.band_mask(freq_ghz)] = differential_spectra(
        freq_ghz,
        ifg[None],
        du,
        opd_m,
        modulation_efficiency,
        freq_step_ghz=freq_step_ghz,
        apodization=apodization,
        kaiser_beta=kaiser_beta,
    )[0]
    return result
