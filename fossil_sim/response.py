"""Response of each detector unit (DU): gain, offset and calibration solution.

A DU with output sign s measures, in each bin of its band,

    d(ν) = g(ν) [s (I₁ − I₂)(ν) + o(ν)] + n(ν),

with g the gain relative to the nominal response (dimensionless, 1 =
nominal), o an offset referred to the sky (MJy/sr) and n the detector noise,
outside the gain. The offset is the emission imbalance between the arms,
o(ν) = ε B_ν(T_inst), with ε signed per DU. The gain is a per-DU scalar times
1 + a Gaussian fluctuation per bin of rms ``gain_rms``, constant over the
mission and keyed by (seed, component, DU) (ADR 0001).

The sky is reconstructed from calibrated spectra d / ĝ − ô, with (ĝ, ô) a
calibration solution: ``nominal`` (ĝ = 1, ô = 0) or ``true`` (the simulated
response, an oracle).
"""

from dataclasses import dataclass

import numpy as np

from .config import Config
from .instrument import DetectorUnit
from .spectra import blackbody_mjy_sr

# SeedSequence component of the gain fluctuations (noise uses 1 and 2)
RESPONSE_COMPONENT_ID = 3


@dataclass(frozen=True)
class DUResponse:
    """Gain and offset (MJy/sr) of each DU on its band bins, by DU name."""

    gain: dict[str, np.ndarray]
    offset_mjy_sr: dict[str, np.ndarray]

    def apply(self, du_spectra: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """Measured spectra g (x + o) of the differential spectra x = s (I₁ − I₂)."""
        return {
            name: self.gain[name] * (x + self.offset_mjy_sr[name])
            for name, x in du_spectra.items()
        }

    def calibrate(self, du_spectra: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
        """Calibrated spectra d / g − o, the inverse of :meth:`apply`."""
        return {
            name: d / self.gain[name] - self.offset_mjy_sr[name]
            for name, d in du_spectra.items()
        }


def nominal_response(freq_ghz: np.ndarray, dus: list[DetectorUnit]) -> DUResponse:
    """Unit gain and zero offset."""
    sizes = {du.name: int(du.band_mask(freq_ghz).sum()) for du in dus}
    return DUResponse(
        {name: np.ones(n) for name, n in sizes.items()},
        {name: np.zeros(n) for name, n in sizes.items()},
    )


def true_response(
    config: Config, freq_ghz: np.ndarray, dus: list[DetectorUnit]
) -> DUResponse:
    """Simulated response of each DU from ``[response]``."""
    rc = config.response
    gain, offset = {}, {}
    for i, du in enumerate(dus):
        band_ghz = freq_ghz[du.band_mask(freq_ghz)]
        fluctuation = np.zeros(band_ghz.size)
        if rc.gain_rms > 0:
            seed = np.random.SeedSequence(
                [config.simulation.master_seed, RESPONSE_COMPONENT_ID, i]
            )
            fluctuation = rc.gain_rms * np.random.default_rng(seed).standard_normal(
                band_ghz.size
            )
        gain[du.name] = rc.gain[i] * (1.0 + fluctuation)
        if np.any(gain[du.name] <= 0):
            raise ValueError(
                f"response: non-positive gain in {du.name}, lower gain_rms"
            )
        offset[du.name] = rc.offset_emissivity[i] * blackbody_mjy_sr(
            band_ghz, rc.instrument_temperature_k
        )
    return DUResponse(gain, offset)


def generation_response(
    config: Config,
    freq_ghz: np.ndarray,
    response: DUResponse,
    generation_freq_ghz: np.ndarray,
    dus: list[DetectorUnit],
) -> DUResponse:
    """``response`` (on the band bins of ``freq_ghz``) on the IFG generation grid.

    The bins of each DU are the generation frequencies it transmits
    (:meth:`~fossil_sim.instrument.DetectorUnit.transmission`), each mapped
    to its nearest band bin. The gain is that of the bin, so its fluctuation
    is constant within a bin; the offset, ε B_ν(T_inst), is that of the bin
    scaled by B_ν(T_inst) / B_bin(T_inst), exact at every frequency.
    """
    rc, fraction = config.response, config.instrument.band_edge_fraction
    gain, offset = {}, {}
    for i, du in enumerate(dus):
        bins_ghz = freq_ghz[du.band_mask(freq_ghz)]
        fine_ghz = generation_freq_ghz[
            du.transmission(generation_freq_ghz, fraction) > 0
        ]
        if not bins_ghz.size:  # a DU without bins on a narrow grid
            gain[du.name] = np.full(fine_ghz.size, rc.gain[i])
            offset[du.name] = np.zeros(fine_ghz.size)
            continue
        nearest = np.abs(fine_ghz[:, None] - bins_ghz[None, :]).argmin(axis=1)
        gain[du.name] = response.gain[du.name][nearest]
        t_k = rc.instrument_temperature_k
        offset[du.name] = (
            response.offset_mjy_sr[du.name][nearest]
            * blackbody_mjy_sr(fine_ghz, t_k)
            / blackbody_mjy_sr(bins_ghz[nearest], t_k)
        )
    return DUResponse(gain, offset)


def calibration_solution(
    config: Config, freq_ghz: np.ndarray, dus: list[DetectorUnit]
) -> DUResponse:
    """Calibration used to reconstruct the sky, from ``calibration.solution``."""
    if config.calibration.solution == "true":
        return true_response(config, freq_ghz, dus)
    return nominal_response(freq_ghz, dus)
