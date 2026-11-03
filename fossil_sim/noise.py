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

With ``ifg.simulate`` the noise is drawn instead on the IFG samples,
:func:`sample_noise_sigma`, and reaches the spectrum through the least-squares
fit; it is drawn lazily per block (:class:`SampleNoise`), keyed by (component,
cycle, DU, chunk of IFGs).
"""

import numpy as np

from .config import Config
from .ifg import _differential_spectrum_solver, opd_samples_m
from .instrument import DetectorUnit

# First entry of the SeedSequence key, separating noise from future components
NOISE_COMPONENT_ID = 1
IFG_NOISE_COMPONENT_ID = 2
# IFGs per seed of the sample noise (fixes the realization: do not change)
SAMPLE_NOISE_CHUNK_IFGS = 256

# 1 MJy/sr = 1e-20 W m⁻² Hz⁻¹ sr⁻¹
SI_TO_MJY_SR = 1e20


def noise_sigma_mjy_sr(config: Config, du: DetectorUnit) -> float:
    """Nominal white-noise standard deviation per IFG and bin of one DU, MJy/sr.

    Exact for the spectral-domain noise model. For simulated IFGs the fit
    changes the bin variances; see :func:`sky_estimate_noise_variance`.
    """
    nc = config.noise
    nep = nc.nep_lf_w_sqrt_hz if du.name.startswith("LF") else nc.nep_hf_w_sqrt_hz
    sigma_si = nep / (
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


class SampleNoise:
    """Sample noise of one DU over a Cooler cycle, (n_ifg, samples_per_ifg).

    Sliced like the array it stands for (``noise[start:stop]``), without
    holding it: rows are drawn in chunks of ``SAMPLE_NOISE_CHUNK_IFGS``, each
    from its own seed, so an IFG gets the same realization whatever the slicing
    and the order of the calls. As for :func:`du_noise_mjy_sr`, the row index
    is the position in the cycle table.
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

    def __len__(self) -> int:
        return self.sigma.size

    def __getitem__(self, rows: slice) -> np.ndarray:
        start, stop, step = rows.indices(len(self))
        if step != 1:
            raise ValueError("SampleNoise supports contiguous slices only")
        if stop <= start:
            return np.empty((0, self.samples_per_ifg))
        size = SAMPLE_NOISE_CHUNK_IFGS
        first, last = start // size, (stop - 1) // size + 1
        noise = np.concatenate(
            [
                np.random.default_rng(
                    np.random.SeedSequence(self.key + [chunk])
                ).standard_normal((size, self.samples_per_ifg))
                for chunk in range(first, last)
            ]
        )[start - first * size : stop - first * size]
        return noise * self.sigma[start:stop, None]


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


def sky_estimate_noise_variance(
    config: Config, freq_ghz: np.ndarray, dus: list[DetectorUnit]
) -> np.ndarray:
    """Noise variance per IFG of the sky estimate (mean over the DUs of each bin).

    With simulated IFGs, propagate the independent sample noise through the
    actual least-squares fit: diag(C) = σ_sample² sum_j solver[k, j]².
    """
    variance = np.zeros(freq_ghz.size)
    count = np.zeros(freq_ghz.size)
    opd_m = opd_samples_m(config.ifg) if config.ifg.simulate else None
    for du in dus:
        mask = du.band_mask(freq_ghz)
        if config.ifg.simulate:
            solver = _differential_spectrum_solver(
                freq_ghz,
                du,
                opd_m,
                config.noise.modulation_efficiency,
                freq_step_ghz=config.instrument.freq_step_ghz,
            )
            variance[mask] += sample_noise_sigma(config, du) ** 2 * np.sum(
                solver**2, axis=1
            )
        else:
            variance[mask] += noise_sigma_mjy_sr(config, du) ** 2
        count[mask] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return variance / count**2
