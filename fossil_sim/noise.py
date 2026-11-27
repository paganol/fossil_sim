"""Spectral-domain white noise, per IFG, frequency bin and DU.

Detector noise white in the IFG maps to a noise flat in power across the
synthesized bins. An IFG of duration t, N samples of δt = t / N with noise
NEP / √(2 δt) each, and a DU power

    P(x) = AΩ α f_opt m Σ_ν [(I₁ + I₂) + s (I₁ − I₂) cos(2π ν x / c)] Δν

(``fossil_sim.ifg``) give, for the least-squares amplitude of each cosine
(variance 2 σ² / N = NEP² / t), per bin

    δS_ν = NEP / (√t · AΩ · Δν · α · f_opt · m)

(no apodization; checked by Monte Carlo in ``test/test_noise.py``). PIXIE
(Kogut et al. 2011; arXiv:2405.20403v2, eq. 4–6) writes √(2t): its δP_ν is the complex
Fourier coefficient, a/2 for a cosine of amplitude a.

``modulation_efficiency`` (m) is a provisional knob for a convention still to
be fixed with the instrument team: 1 when one DU sees the full difference,
1/4 for a single-polarization bolometer of a Martin-Puplett FTS (PIXIE 2025,
eq. 1), 1/2 for a DU summing two. Once fixed, fold it into the signal model
and remove the key. IFGs in the Recycling window get ``recycling_noise_factor``
times this level. Realizations are keyed by (component, cycle, DU), never by
MPI rank (ADR 0001).

The NEP is either given per band (``nep_model = "fixed"``) or computed for
each DU from its optical loading (``"photon"``, :func:`photon_nep_w_sqrt_hz`,
the model of Coulon, Aghanim & Maffei 2026, arXiv:2609.17829, Eq. 15–16).

With ``ifg.simulate`` the noise is drawn instead on the IFG samples,
:func:`sample_noise_sigma`, and reaches the spectrum through the least-squares
fit; it is drawn lazily per block (:class:`SampleNoise`), keyed by (component,
cycle, DU, chunk of IFGs). There it can have a 1/f component
(``fknee_lf/hf_mhz``, generated with ``ducc0.misc.OofaNoise``).
"""

import ducc0
import numpy as np

from .config import Config
from .constants import C_LIGHT_M_S, H_PLANCK_J_S
from .ifg import differential_spectrum_solver, opd_samples_m
from .input_sky import isotropic_background_mjy_sr
from .instrument import DetectorUnit
from .spectra import blackbody_mjy_sr

# First entry of the SeedSequence key, separating noise from future components
NOISE_COMPONENT_ID = 1
IFG_NOISE_COMPONENT_ID = 2
# IFGs per seed of the sample noise (fixes the realization: do not change)
SAMPLE_NOISE_CHUNK_IFGS = 256
# Zero padding of the fit rows for the 1/f variance (converged to 1e-8 at 16)
ONE_OVER_F_PADDING = 16

# 1 MJy/sr = 1e-20 W m⁻² Hz⁻¹ sr⁻¹
SI_TO_MJY_SR = 1e20
# Frequency step of the photon-NEP integral (smooth integrand: 0.1 % at 1 GHz)
PHOTON_NEP_STEP_GHZ = 0.5


def du_nep_w_sqrt_hz(config: Config, du: DetectorUnit) -> float:
    """NEP of one DU: given per band, or the photon NEP of its loading."""
    nc = config.noise
    if nc.nep_model == "photon":
        return photon_nep_w_sqrt_hz(config, du)
    return nc.nep_lf_w_sqrt_hz if du.name.startswith("LF") else nc.nep_hf_w_sqrt_hz


def optical_loading_w_hz(
    config: Config, du: DetectorUnit, freq_ghz: np.ndarray
) -> np.ndarray:
    """Spectral power on one DU, P_ν in W/Hz (Coulon et al. 2026, Eq. 15).

    P_ν = m α AΩ w(ν) [f_opt (I_sky + I_ref) + Σ_k ε_k B_ν(T_k)]: the
    unmodulated term of the DU power (``fossil_sim.ifg``, all frequencies
    adding incoherently: an upper bound on the photon noise), with w the DU
    transmission, I_sky the Isotropic backgrounds (CMB and distortions, CIB,
    extragalactic CO; no Galactic or zodiacal emission), I_ref the BBIR at the
    mean set point, and the emission of the optical elements reaching the
    detectors unattenuated. m = 1/2 is the Martin-Puplett value for a DU
    absorbing both polarisations (Coulon et al. 2026, Eq. 11–12).
    """
    nc = config.noise
    sky = isotropic_background_mjy_sr(freq_ghz, config.sky)
    reference = blackbody_mjy_sr(freq_ghz, float(np.mean(config.bbir.setpoints_k)))
    radiance = nc.optical_efficiency * (sky + reference)
    for emissivity, temperature in zip(
        nc.optics_emissivity, nc.optics_temperature_k, strict=True
    ):
        radiance = radiance + emissivity * blackbody_mjy_sr(freq_ghz, temperature)
    transmission = du.transmission(freq_ghz, config.instrument.band_edge_fraction)
    return (
        nc.modulation_efficiency
        * nc.absorptivity
        * nc.etendue_cm2_sr
        * 1e-4
        * transmission
        * radiance
        / SI_TO_MJY_SR
    )


def photon_nep_w_sqrt_hz(config: Config, du: DetectorUnit) -> float:
    """Photon NEP of one DU from its loading (Coulon et al. 2026, Eq. 16).

    NEP² = ∫ [2 h ν P_ν + 2 P_ν² / N(ν)] dν over the DU band, shot and wave
    noise, with N = AΩ ν² / c² modes and P_ν from :func:`optical_loading_w_hz`.
    """
    edge = 1 + config.instrument.band_edge_fraction / 2
    freq_ghz = np.arange(
        du.band_min_ghz / edge, du.band_max_ghz * edge, PHOTON_NEP_STEP_GHZ
    )
    power = optical_loading_w_hz(config, du, freq_ghz)
    freq_hz = freq_ghz * 1e9
    etendue_m2_sr = config.noise.etendue_cm2_sr * 1e-4
    integrand = 2 * H_PLANCK_J_S * freq_hz * power + 2 * power**2 * (
        C_LIGHT_M_S**2 / (freq_hz**2 * etendue_m2_sr)
    )
    return float(np.sqrt(np.trapezoid(integrand, freq_hz)))


def noise_sigma_mjy_sr(config: Config, du: DetectorUnit) -> float:
    """Nominal white-noise standard deviation per IFG and bin of one DU, MJy/sr.

    Exact for the spectral-domain noise model. For simulated IFGs the fit
    changes the bin variances; see :func:`sky_estimate_noise_variance`.
    """
    nc = config.noise
    sigma_si = du_nep_w_sqrt_hz(config, du) / (
        np.sqrt(config.timeline.ifg_duration_s)
        * nc.etendue_cm2_sr
        * 1e-4
        * config.instrument.freq_step_ghz
        * 1e9
        * nc.absorptivity
        * nc.optical_efficiency
        * nc.modulation_efficiency
    )
    return sigma_si * SI_TO_MJY_SR


def sample_noise_sigma(config: Config, du: DetectorUnit) -> float:
    """White-noise standard deviation per IFG sample of one DU, MJy/sr × GHz.

    NEP / √(2 δt) in the units of :func:`fossil_sim.ifg.du_interferograms`,
    giving :func:`noise_sigma_mjy_sr` for orthogonal cosine amplitudes.
    """
    return (
        noise_sigma_mjy_sr(config, du)
        * config.instrument.freq_step_ghz
        * config.noise.modulation_efficiency
        * np.sqrt(config.ifg.samples_per_ifg / 2)
    )


def du_fknee_hz(config: Config, du: DetectorUnit) -> float:
    """Knee frequency of the 1/f noise of one DU, Hz (0: white noise)."""
    nc = config.noise
    fknee_mhz = nc.fknee_lf_mhz if du.name.startswith("LF") else nc.fknee_hf_mhz
    return fknee_mhz * 1e-3


class SampleNoise:
    """Sample noise of one DU over a Cooler cycle, (n_ifg, samples_per_ifg).

    Sliced like the array it stands for (``noise[start:stop]``), without
    holding it: rows are drawn in chunks of ``SAMPLE_NOISE_CHUNK_IFGS``, each
    from its own seed, so an IFG gets the same realization whatever the slicing
    and the order of the calls. As for :func:`du_noise_mjy_sr`, the row index
    is the position in the cycle table.

    With a knee frequency (:func:`du_fknee_hz`) the white stream of the cycle,
    the IFGs back to back, goes through ``ducc0.misc.OofaNoise`` (as in
    litebird_sim), PSD σ² [(f² + f_knee²) / (f² + f_min²)]^(α/2): σ is the
    white level, unchanged. The filter keeps its state between slices, so
    contiguous slices in order cost one pass; a slice before the last one
    restarts the filter from the first IFG of the cycle. The stream starts at
    the cycle start (one realization per cycle, ADR 0001), and the Recycling
    factor scales the filtered noise of its IFGs.
    """

    def __init__(
        self,
        config: Config,
        du: DetectorUnit,
        du_index: int,
        cycle_id: int,
        in_recycling: np.ndarray,
    ):
        self.key = [
            config.simulation.master_seed,
            IFG_NOISE_COMPONENT_ID,
            cycle_id,
            du_index,
        ]
        self.samples_per_ifg = config.ifg.samples_per_ifg
        self.sigma = sample_noise_sigma(config, du) * np.where(
            in_recycling, config.noise.recycling_noise_factor, 1.0
        )
        self.fknee_hz = du_fknee_hz(config, du)
        self.alpha = config.noise.one_over_f_alpha
        self.fmin_hz = config.noise.one_over_f_fmin_hz
        self.sample_rate_hz = self.samples_per_ifg / config.timeline.ifg_duration_s
        self._filter = None
        self._next_row = 0

    def __len__(self) -> int:
        return self.sigma.size

    def __getitem__(self, rows: slice) -> np.ndarray:
        start, stop, step = rows.indices(len(self))
        if step != 1:
            raise ValueError("SampleNoise supports contiguous slices only")
        if stop <= start:
            return np.empty((0, self.samples_per_ifg))
        if self.fknee_hz > 0:
            noise = self._one_over_f(start, stop)
        else:
            noise = self._white(start, stop)
        return noise * self.sigma[start:stop, None]

    def _white(self, start: int, stop: int) -> np.ndarray:
        """Unit white noise of rows start:stop, from their chunk seeds."""
        size = SAMPLE_NOISE_CHUNK_IFGS
        first, last = start // size, (stop - 1) // size + 1
        return np.concatenate(
            [
                np.random.default_rng(
                    np.random.SeedSequence(self.key + [chunk])
                ).standard_normal((size, self.samples_per_ifg))
                for chunk in range(first, last)
            ]
        )[start - first * size : stop - first * size]

    def _one_over_f(self, start: int, stop: int) -> np.ndarray:
        """Unit-white-level 1/f noise of rows start:stop, continuing the stream."""
        if self._filter is None or start < self._next_row:
            self._filter = ducc0.misc.OofaNoise(
                sigmawhite=1.0,
                f_knee=self.fknee_hz,
                f_min=self.fmin_hz,
                f_samp=self.sample_rate_hz,
                slope=-self.alpha,
            )
            self._next_row = 0
        # advance over skipped rows, a chunk at a time
        for row in range(self._next_row, start, SAMPLE_NOISE_CHUNK_IFGS):
            end = min(row + SAMPLE_NOISE_CHUNK_IFGS, start)
            self._filter.filterGaussian(self._white(row, end).ravel())
        noise = self._filter.filterGaussian(self._white(start, stop).ravel())
        self._next_row = stop
        return noise.reshape(stop - start, self.samples_per_ifg)


def du_noise_mjy_sr(
    config: Config,
    du: DetectorUnit,
    du_index: int,
    cycle_id: int,
    in_recycling: np.ndarray,
    num_bins: int,
) -> np.ndarray:
    """Noise realization of one DU over a whole Cooler cycle, (n_ifg, num_bins).

    ``in_recycling`` marks the IFGs of the cycle in the Recycling window.
    The realization depends on the number of IFGs, so ``recycling_data =
    "off"`` and ``"flagged"`` give different noise on the same IFG.
    """
    seed = np.random.SeedSequence(
        [config.simulation.master_seed, NOISE_COMPONENT_ID, cycle_id, du_index]
    )
    noise = np.random.default_rng(seed).standard_normal((in_recycling.size, num_bins))
    scale = np.where(in_recycling, config.noise.recycling_noise_factor, 1.0)
    return noise * (noise_sigma_mjy_sr(config, du) * scale)[:, None]


def one_over_f_variance_factor(
    config: Config, du: DetectorUnit, solver: np.ndarray
) -> np.ndarray:
    """Variance of each fitted bin with the 1/f noise, relative to white noise.

    For stationary sample noise of PSD σ² S(f) the fit row L_k has variance
    σ² ∫ S(f) |L̂_k(f)|² df; by Parseval this is σ² Σ_n L_kn² for S = 1, so the
    factor is the S-weighted mean of |L̂_k|² over the frequencies (rows zero
    padded). The fitted constant makes L̂_k vanish at f = 0, so the noise below
    the IFG rate barely reaches the spectrum: only bins whose electrical
    frequency ν v_OPD / c is near f_knee are raised. Correlations between IFGs
    (from 1/f below the IFG rate) are neglected; they are below 1 % for
    f_knee up to 2 Hz (Monte Carlo, ``test/test_noise.py``). Ones without 1/f.
    """
    fknee_hz = du_fknee_hz(config, du)
    if fknee_hz == 0:
        return np.ones(solver.shape[0])
    nc = config.noise
    size = solver.shape[1] * ONE_OVER_F_PADDING
    rate_hz = config.ifg.samples_per_ifg / config.timeline.ifg_duration_s
    freq_hz = np.abs(np.fft.fftfreq(size, 1 / rate_hz))
    psd = ((freq_hz**2 + fknee_hz**2) / (freq_hz**2 + nc.one_over_f_fmin_hz**2)) ** (
        nc.one_over_f_alpha / 2
    )
    power = np.abs(np.fft.fft(solver, n=size, axis=1)) ** 2
    return (power @ psd) / size / np.sum(solver**2, axis=1)


def ifg_solvers(
    config: Config, freq_ghz: np.ndarray, dus: list[DetectorUnit]
) -> dict[str, np.ndarray]:
    """:func:`fossil_sim.ifg.differential_spectrum_solver` of each DU, by name."""
    opd_m = opd_samples_m(config.ifg)
    return {
        du.name: differential_spectrum_solver(
            freq_ghz,
            du,
            opd_m,
            config.noise.modulation_efficiency,
            freq_step_ghz=config.instrument.freq_step_ghz,
            apodization=config.ifg.apodization,
            kaiser_beta=config.ifg.kaiser_beta,
        )
        for du in dus
    }


def sky_estimate_noise_variance(
    config: Config,
    freq_ghz: np.ndarray,
    dus: list[DetectorUnit],
    solvers: dict[str, np.ndarray] | None = None,
    gain: dict[str, np.ndarray] | None = None,
) -> np.ndarray:
    """Noise variance per IFG of the sky estimate (mean over the DUs of each bin).

    ``gain`` maps DU names to the calibration gain ĝ on the band bins
    (:mod:`fossil_sim.response`): the detector noise is outside the gain, so
    the calibrated variance is divided by ĝ².

    With simulated IFGs, propagate the sample noise through the actual
    least-squares fit: diag(C) = σ_sample² sum_j solver[k, j]² for white noise,
    times :func:`one_over_f_variance_factor` with 1/f noise.
    ``solvers`` maps DU names to :func:`fossil_sim.ifg.differential_spectrum_solver`
    (computed here if not given).
    """
    variance = np.zeros(freq_ghz.size)
    count = np.zeros(freq_ghz.size)
    if config.ifg.simulate and solvers is None:
        solvers = ifg_solvers(config, freq_ghz, dus)
    for du in dus:
        mask = du.band_mask(freq_ghz)
        if config.ifg.simulate:
            solver = solvers[du.name]
            du_variance = (
                sample_noise_sigma(config, du) ** 2
                * np.sum(solver**2, axis=1)
                * one_over_f_variance_factor(config, du, solver)
            )
        else:
            du_variance = noise_sigma_mjy_sr(config, du) ** 2
        if gain is not None:
            du_variance = du_variance / gain[du.name] ** 2
        variance[mask] += du_variance
        count[mask] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return variance / count**2
