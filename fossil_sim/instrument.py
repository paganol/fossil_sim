"""Frequency grid and detector units (DUs).

Each FTS output (L, R) is split by a dichroic into a low-frequency and a
high-frequency DU: LF-L, HF-L, LF-R, HF-R. The modulated signal of the L
output is +(I₁ − I₂), the one of the R output −(I₁ − I₂) (CONTEXT.md, "Mode").
"""

from dataclasses import dataclass

import numpy as np

from .config import InstrumentConfig


@dataclass(frozen=True)
class DetectorUnit:
    name: str
    output_sign: int  # +1 for the L output, -1 for the R output
    band_min_ghz: float
    band_max_ghz: float

    def band_mask(self, freq_ghz: np.ndarray) -> np.ndarray:
        """Ideal top-hat band: True where ``band_min_ghz <= ν < band_max_ghz``."""
        return (freq_ghz >= self.band_min_ghz) & (freq_ghz < self.band_max_ghz)


def frequency_grid_ghz(config: InstrumentConfig) -> np.ndarray:
    """Frequency bin centres, as ``np.arange(min, max, step)``.

    >>> from fossil_sim.config import InstrumentConfig
    >>> grid = frequency_grid_ghz(InstrumentConfig(60.0, 2000.0, 15.0, 50.0, 2000.0, 300.0, 1.6))
    >>> grid[0], grid[-1], grid.size
    (np.float64(60.0), np.float64(1995.0), 130)
    """
    return np.arange(config.freq_min_ghz, config.freq_max_ghz, config.freq_step_ghz)


def detector_units(config: InstrumentConfig) -> list[DetectorUnit]:
    """The four DUs, with ideal dichroic split."""
    low = (config.band_min_ghz, config.dichroic_split_ghz)
    high = (config.dichroic_split_ghz, config.band_max_ghz)
    return [
        DetectorUnit("LF-L", +1, *low),
        DetectorUnit("HF-L", +1, *high),
        DetectorUnit("LF-R", -1, *low),
        DetectorUnit("HF-R", -1, *high),
    ]
