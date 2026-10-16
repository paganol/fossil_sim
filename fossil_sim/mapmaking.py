"""Hit maps, accumulated in streaming over Cooler cycles."""

from pathlib import Path

import healpy as hp
import numpy as np
from numba import njit

from .config import MODES


@njit(cache=True)
def _accumulate_hits(pixels, mode, time_s, hits, first_hit_time_s):
    for i in range(pixels.size):
        hits[mode[i], pixels[i]] += 1
        first_hit_time_s[pixels[i]] = min(first_hit_time_s[pixels[i]], time_s[i])


class HitMaps:
    """Mode-resolved hit maps and first-hit times at several resolutions.

    Only valid IFGs (``flags == 0``) are counted. Maps are in J2000 ecliptic
    coordinates, RING ordering.

    Attributes
    ----------
    hits : dict[int, ndarray]
        ``hits[nside]`` has shape ``(len(MODES), npix)``; row ``i`` counts
        IFGs taken in Mode ``MODES[i]``.
    first_hit_time_s : dict[int, ndarray]
        Time of the first valid hit of each pixel (``inf`` if never hit).
    """

    def __init__(self, nsides: tuple[int, ...]) -> None:
        self.hits = {
            n: np.zeros((len(MODES), hp.nside2npix(n)), dtype=np.int64) for n in nsides
        }
        self.first_hit_time_s = {n: np.full(hp.nside2npix(n), np.inf) for n in nsides}

    def add(self, boresight: np.ndarray, mode: np.ndarray, time_s: np.ndarray) -> None:
        """Accumulate valid IFGs (already selected by the caller)."""
        for nside, hits in self.hits.items():
            pixels = hp.vec2pix(
                nside, boresight[:, 0], boresight[:, 1], boresight[:, 2]
            )
            _accumulate_hits(pixels, mode, time_s, hits, self.first_hit_time_s[nside])

    def total(self, nside: int) -> np.ndarray:
        """Hits summed over Modes."""
        return self.hits[nside].sum(axis=0)

    def coverage_fraction(self, nside: int, time_s: np.ndarray) -> np.ndarray:
        """Fraction of the sky with at least one valid hit before each time."""
        first_s = np.sort(self.first_hit_time_s[nside])
        return np.searchsorted(first_s, time_s, side="right") / first_s.size

    def save(self, output_dir: str | Path) -> list[Path]:
        """Write one FITS file per nside: total hits first, then one map per Mode."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        paths = []
        for nside, hits in self.hits.items():
            path = output_dir / f"hits_nside{nside:04d}.fits"
            hp.write_map(
                path,
                np.vstack([hits.sum(axis=0), hits]),
                column_names=["TOTAL", *MODES],
                coord="E",
                dtype=np.int64,
                overwrite=True,
            )
            paths.append(path)
        return paths
