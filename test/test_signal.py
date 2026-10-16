import numpy as np

from fossil_sim.instrument import detector_units, frequency_grid_ghz
from fossil_sim.operations import MODE_CODES
from fossil_sim.signal import (
    bbir_spectra_mjy_sr,
    du_differential_spectra,
    port_spectra_mjy_sr,
    sky_estimate_mjy_sr,
)
from fossil_sim.spectra import blackbody_mjy_sr


def test_modes_signs_and_sky_estimate(config):
    freq_ghz = frequency_grid_ghz(config.instrument)
    dus = detector_units(config.instrument)
    modes = np.array([MODE_CODES[m] for m in ("B10", "B01", "B00", "B11")])
    sky = np.stack([blackbody_mjy_sr(freq_ghz, 2.7 + 0.01 * i) for i in range(4)])
    port_k = np.array([[2.6, 2.8]] * 4)
    bbir = bbir_spectra_mjy_sr(freq_ghz, port_k)
    port1, port2 = port_spectra_mjy_sr(sky, bbir, modes)
    spectra = du_differential_spectra(freq_ghz, port1, port2, dus)

    lf = dus[0].band_mask(freq_ghz)
    b1, b2 = bbir[0, 0, lf], bbir[0, 1, lf]
    assert np.allclose(spectra["LF-L"][0], sky[0, lf] - b2)  # B10: +(S − B₂)
    assert np.allclose(spectra["LF-L"][1], -(sky[1, lf] - b1))  # B01: −(S − B₁)
    assert np.allclose(spectra["LF-L"][2], b1 - b2)  # B00
    assert np.allclose(spectra["LF-L"][3], 0.0)  # B11
    assert np.allclose(spectra["LF-R"], -spectra["LF-L"])

    estimate = sky_estimate_mjy_sr(freq_ghz, spectra, dus, bbir, modes)
    assert np.allclose(estimate[:2], sky[:2], rtol=1e-12)
    assert np.all(np.isnan(estimate[2:]))
