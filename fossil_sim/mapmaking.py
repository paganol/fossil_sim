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


@njit(cache=True)
def _accumulate_spectra(pixels, values, sums, counts):
    for i in range(pixels.size):
        counts[pixels[i]] += 1
        sums[pixels[i]] += values[i]


DIAGNOSTIC_OUTPUTS = ("L", "R")
DIAGNOSTIC_MODES = ("B10", "B01")


class SpectralMaps:
    """Binned spectral maps, accumulated in streaming (J2000 ecliptic, RING).

    - ``sky``: the reconstructed sky spectrum (perfect calibration), averaged
      over all valid science IFGs and DUs;
    - ``raw``: the differential spectrum s (I₁ − I₂) per (output, Mode), the
      diagnostic product (B10 and B01 have opposite sky sign).

    Sums are stored as (npix, n_freq) for contiguous accumulation.
    ``noise_variance_per_ifg`` (per bin, MJy²/sr²), when given, yields the
    white-noise variance maps of the sky, variance / hits.
    """

    def __init__(
        self,
        nsides: tuple[int, ...],
        freq_ghz: np.ndarray,
        noise_variance_per_ifg: np.ndarray | None = None,
    ) -> None:
        self.freq_ghz = np.asarray(freq_ghz)
        self.noise_variance_per_ifg = noise_variance_per_ifg
        n_freq = self.freq_ghz.size
        self.sky_sum = {n: np.zeros((hp.nside2npix(n), n_freq)) for n in nsides}
        self.sky_count = {n: np.zeros(hp.nside2npix(n), dtype=np.int64) for n in nsides}
        self.raw_sum = {
            n: np.zeros(
                (
                    len(DIAGNOSTIC_OUTPUTS),
                    len(DIAGNOSTIC_MODES),
                    hp.nside2npix(n),
                    n_freq,
                )
            )
            for n in nsides
        }
        self.raw_count = {
            n: np.zeros(
                (len(DIAGNOSTIC_OUTPUTS), len(DIAGNOSTIC_MODES), hp.nside2npix(n)),
                dtype=np.int64,
            )
            for n in nsides
        }

    def add(
        self,
        boresight: np.ndarray,
        sky_estimate: np.ndarray,
        raw_by_output: dict[str, np.ndarray],
        mode_names: np.ndarray,
    ) -> None:
        """Accumulate valid IFGs.

        ``sky_estimate`` is (n, n_freq), NaN for non-science IFGs;
        ``raw_by_output[o]`` is the full-band differential spectrum of output o.
        """
        science = ~np.isnan(sky_estimate[:, 0])
        for nside, sums in self.sky_sum.items():
            pixels = hp.vec2pix(
                nside, boresight[:, 0], boresight[:, 1], boresight[:, 2]
            )
            _accumulate_spectra(
                pixels[science], sky_estimate[science], sums, self.sky_count[nside]
            )
            for i_out, output in enumerate(DIAGNOSTIC_OUTPUTS):
                for i_mode, mode in enumerate(DIAGNOSTIC_MODES):
                    sel = mode_names == mode
                    _accumulate_spectra(
                        pixels[sel],
                        raw_by_output[output][sel],
                        self.raw_sum[nside][i_out, i_mode],
                        self.raw_count[nside][i_out, i_mode],
                    )

    @staticmethod
    def _average(sums: np.ndarray, counts: np.ndarray) -> np.ndarray:
        with np.errstate(invalid="ignore", divide="ignore"):
            maps = (sums / counts[..., None]).swapaxes(-1, -2)
        return np.where(np.isfinite(maps), maps, hp.UNSEEN)

    def sky_map(self, nside: int) -> np.ndarray:
        """Reconstructed sky, shape (n_freq, npix), MJy/sr; UNSEEN where not observed."""
        return self._average(self.sky_sum[nside], self.sky_count[nside])

    def sky_variance_map(self, nside: int) -> np.ndarray:
        """White-noise variance of the sky map, (n_freq, npix), (MJy/sr)²; UNSEEN where not observed."""
        counts = self.sky_count[nside]
        with np.errstate(invalid="ignore", divide="ignore"):
            maps = self.noise_variance_per_ifg[:, None] / counts[None, :]
        return np.where(np.isfinite(maps), maps, hp.UNSEEN)

    def raw_maps(self, nside: int) -> np.ndarray:
        """Differential spectra, shape (2 outputs, 2 Modes, n_freq, npix), MJy/sr."""
        return self._average(self.raw_sum[nside], self.raw_count[nside])

    def save(self, output_dir: str | Path) -> list[Path]:
        """Write the sky and diagnostic maps, one FITS file each per nside."""
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        freq_names = [f"{f:07.2f}GHz" for f in self.freq_ghz]
        paths = []
        for nside in self.sky_sum:
            path = output_dir / f"sky_nside{nside:04d}.fits"
            hp.write_map(
                path,
                self.sky_map(nside),
                column_names=freq_names,
                column_units="MJy/sr",
                coord="E",
                overwrite=True,
            )
            paths.append(path)
            if self.noise_variance_per_ifg is not None:
                path = output_dir / f"sky_variance_nside{nside:04d}.fits"
                hp.write_map(
                    path,
                    self.sky_variance_map(nside),
                    column_names=freq_names,
                    column_units="(MJy/sr)^2",
                    coord="E",
                    overwrite=True,
                )
                paths.append(path)
            raw = self.raw_maps(nside)
            names = [
                f"{o}_{m}_{f}"
                for o in DIAGNOSTIC_OUTPUTS
                for m in DIAGNOSTIC_MODES
                for f in freq_names
            ]
            path = output_dir / f"raw_nside{nside:04d}.fits"
            hp.write_map(
                path,
                raw.reshape(-1, raw.shape[-1]),
                column_names=names,
                column_units="MJy/sr",
                coord="E",
                overwrite=True,
            )
            paths.append(path)
        return paths
