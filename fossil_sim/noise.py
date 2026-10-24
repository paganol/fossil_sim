"""Spectral-domain white noise, per IFG, frequency bin and DU.

Detector noise white in the IFG maps to a noise flat in power across the
synthesized bins. Following PIXIE (Kogut et al. 2011, §2), an IFG of duration
t with NEP gives per bin

    δS_ν = NEP / (√(2 t) · AΩ · Δν · α · f_opt)

(no apodization). IFGs in the Recycling window get ``recycling_noise_factor``
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
        np.sqrt(2.0 * config.timeline.ifg_duration_s)
        * nc.etendue_cm2_sr
        * 1e-4
        * config.instrument.freq_step_ghz
        * 1e9
        * nc.absorptivity
        * nc.optical_efficiency
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
