"""Spectral-domain signal: per-IFG Port spectra -> differential spectra per DU.

Port k sees the sky when its BSM points to the sky (Mode ``Bxy``, CONTEXT.md),
otherwise the BBIR spectrum of that Port (ideal blackbody at the Port
temperature, emissivity 1). The L output measures +(I₁ − I₂), the R output
−(I₁ − I₂), each DU restricted to its band.
"""

import numpy as np

from .instrument import DetectorUnit
from .operations import mode_sky_ports
from .spectra import blackbody_mjy_sr


def bbir_spectra_mjy_sr(
    freq_ghz: np.ndarray, port_temperature_k: np.ndarray
) -> np.ndarray:
    """BBIR spectrum seen by each Port, shape (n_ifg, 2, n_freq).

    Temperatures repeat a few set-points, so blackbodies are computed once per
    distinct value.
    """
    unique_k, inverse = np.unique(port_temperature_k, return_inverse=True)
    spectra = np.stack([blackbody_mjy_sr(freq_ghz, t) for t in unique_k])
    return spectra[inverse.reshape(port_temperature_k.shape)]


def port_spectra_mjy_sr(
    sky_mjy_sr: np.ndarray, bbir_mjy_sr: np.ndarray, mode: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Spectra entering Port 1 and Port 2, each of shape (n_ifg, n_freq).

    ``sky_mjy_sr`` may carry extra axes, (n_ifg, ..., n_freq), with
    ``bbir_mjy_sr`` (n_ifg, 2, ..., n_freq) broadcasting to them.
    """
    port1_on_sky, port2_on_sky = mode_sky_ports(mode)
    shape = (-1,) + (1,) * (sky_mjy_sr.ndim - 1)
    port1 = np.where(port1_on_sky.reshape(shape), sky_mjy_sr, bbir_mjy_sr[:, 0])
    port2 = np.where(port2_on_sky.reshape(shape), sky_mjy_sr, bbir_mjy_sr[:, 1])
    return port1, port2


def du_differential_spectra(
    freq_ghz: np.ndarray,
    port1_mjy_sr: np.ndarray,
    port2_mjy_sr: np.ndarray,
    dus: list[DetectorUnit],
) -> dict[str, np.ndarray]:
    """Differential spectrum s (I₁ − I₂) of each DU on its band bins, (n_ifg, n_band)."""
    difference = port1_mjy_sr - port2_mjy_sr
    return {
        du.name: du.output_sign * difference[:, du.band_mask(freq_ghz)] for du in dus
    }


def sky_estimate_mjy_sr(
    freq_ghz: np.ndarray,
    du_spectra: dict[str, np.ndarray],
    dus: list[DetectorUnit],
    bbir_mjy_sr: np.ndarray,
    mode: np.ndarray,
) -> np.ndarray:
    """Sky spectrum reconstructed from the DUs, assuming perfect calibration.

    Only science Modes (one Port on the sky) give an estimate; other IFGs are
    NaN. For each DU, Ŝ = B_k + σ s d with B_k the BBIR of the Port on the
    BBIR and σ = +1 (B10) or −1 (B01); L and R estimates are averaged.
    """
    port1_on_sky, port2_on_sky = mode_sky_ports(mode)
    science = port1_on_sky ^ port2_on_sky
    sigma = np.where(port1_on_sky, 1.0, -1.0)[:, None]
    bbir_other = np.where(port1_on_sky[:, None], bbir_mjy_sr[:, 1], bbir_mjy_sr[:, 0])
    total = np.zeros((mode.size, freq_ghz.size))
    count = np.zeros(freq_ghz.size)
    for du in dus:
        mask = du.band_mask(freq_ghz)
        total[:, mask] += (
            bbir_other[:, mask] + sigma * du.output_sign * du_spectra[du.name]
        )
        count[mask] += 1
    estimate = total / count
    estimate[~science] = np.nan
    return estimate
