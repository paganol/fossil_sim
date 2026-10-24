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
"""

import numpy as np

from .config import Config
from .instrument import DetectorUnit

# First entry of the SeedSequence key, separating noise from future components
NOISE_COMPONENT_ID = 1

# 1 MJy/sr = 1e-20 W m⁻² Hz⁻¹ sr⁻¹
SI_TO_MJY_SR = 1e20


def noise_sigma_mjy_sr(config: Config, du: DetectorUnit) -> float:
    """White-noise standard deviation per IFG and per bin of one DU, MJy/sr."""
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
    """Noise variance per IFG of the sky estimate (mean over the DUs of each bin)."""
    variance = np.zeros(freq_ghz.size)
    count = np.zeros(freq_ghz.size)
    for du in dus:
        mask = du.band_mask(freq_ghz)
        variance[mask] += noise_sigma_mjy_sr(config, du) ** 2
        count[mask] += 1
    with np.errstate(invalid="ignore", divide="ignore"):
        return variance / count**2
