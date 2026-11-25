"""Frequency grid and detector units (DUs).

Each FTS output (L, R) is split by a dichroic into a low-frequency and a
high-frequency DU: LF-L, HF-L, LF-R, HF-R. The modulated signal of the L
output is +(I₁ − I₂), the one of the R output −(I₁ − I₂) (CONTEXT.md, "Mode").
"""

from dataclasses import dataclass

import numpy as np

from .config import Config, InstrumentConfig

# Tolerance, in bins, on (max - min) / step when counting the grid bins
GRID_TOLERANCE = 1e-9


@dataclass(frozen=True)
class DetectorUnit:
    name: str
    output_sign: int  # +1 for the L output, -1 for the R output
    band_min_ghz: float
    band_max_ghz: float

    def band_mask(self, freq_ghz: np.ndarray) -> np.ndarray:
        """Ideal top-hat band: True where ``band_min_ghz <= ν < band_max_ghz``."""
        return (freq_ghz >= self.band_min_ghz) & (freq_ghz < self.band_max_ghz)

    def transmission(self, freq_ghz: np.ndarray, edge_fraction: float) -> np.ndarray:
        """Band transmission in [0, 1], with cosine edges.

        Each edge at frequency e rises (or falls) as ½ (1 + sin(π t)),
        t = (ν − e) / (``edge_fraction`` × e), over |t| <= ½, so the
        transmission is ½ at e. The two DUs of an output share the dichroic
        edge and sum to 1. ``edge_fraction = 0`` gives :meth:`band_mask`.
        """
        freq_ghz = np.asarray(freq_ghz, dtype=float)
        if edge_fraction == 0:
            return self.band_mask(freq_ghz).astype(float)
        return _edge(freq_ghz, self.band_min_ghz, edge_fraction) * _edge(
            -freq_ghz, -self.band_max_ghz, edge_fraction
        )


def frequency_grid_ghz(config: InstrumentConfig) -> np.ndarray:
    """Frequency bin centres ``min + k step`` below ``max`` (excluded).

    The bin count is an integer computed with a tolerance, so rounding in
    ``(max - min) / step`` adds no bin (``np.arange(1, 1.3, 0.1)`` has four).

    >>> from fossil_sim.config import InstrumentConfig
    >>> grid = frequency_grid_ghz(InstrumentConfig(60.0, 2000.0, 15.0, 50.0, 2000.0, 300.0, 1.6))
    >>> grid[0], grid[-1], grid.size
    (np.float64(60.0), np.float64(1995.0), 130)
    """
    span = (config.freq_max_ghz - config.freq_min_ghz) / config.freq_step_ghz
    n_bins = int(np.ceil(span - GRID_TOLERANCE))
    return config.freq_min_ghz + config.freq_step_ghz * np.arange(n_bins)


def _edge(freq_ghz: np.ndarray, edge_ghz: float, fraction: float) -> np.ndarray:
    """Rising cosine edge at ``edge_ghz`` (pass negated values for a falling one)."""
    t = (freq_ghz - edge_ghz) / (fraction * abs(edge_ghz))
    return 0.5 * (1.0 + np.sin(np.pi * np.clip(t, -0.5, 0.5)))


def generation_frequency_grid_ghz(config: Config) -> np.ndarray:
    """Frequency grid the interferograms are generated on (docs: Interferogram).

    Without ``ifg.generation_freq_step_ghz`` it is the instrument grid
    (:func:`frequency_grid_ghz`). Otherwise the multiples k × step up to
    ``ifg.generation_freq_max_ghz`` where some DU transmits; each point
    stands for a band ``step`` wide, and the instrument bin centres are
    among the points.
    """
    step = config.ifg.generation_freq_step_ghz
    if step is None:
        return frequency_grid_ghz(config.instrument)
    n = int(np.floor(config.ifg.generation_freq_max_ghz / step + GRID_TOLERANCE))
    freq_ghz = step * np.arange(1, n + 1)
    fraction = config.instrument.band_edge_fraction
    transmits = np.zeros(freq_ghz.size, dtype=bool)
    for du in detector_units(config.instrument):
        transmits |= du.transmission(freq_ghz, fraction) > 0
    return freq_ghz[transmits]


def sky_frequency_grid_ghz(config: Config) -> np.ndarray:
    """Frequencies the sky is evaluated on: generation grid with ``ifg.simulate``."""
    if config.ifg.simulate:
        return generation_frequency_grid_ghz(config)
    return frequency_grid_ghz(config.instrument)


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
