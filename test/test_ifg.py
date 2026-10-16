import numpy as np

from fossil_sim.ifg import (
    differential_spectrum,
    interferogram,
    mirror_positions_m,
    opd_samples_m,
)
from fossil_sim.instrument import detector_units, frequency_grid_ghz
from fossil_sim.spectra import blackbody_mjy_sr, y_distortion_mjy_sr


def test_y_distortion_null_at_217_ghz(config):
    t_cmb_k = config.sky.t_cmb_k
    freq_ghz = np.array([150.0, 217.5, 300.0])
    delta = y_distortion_mjy_sr(freq_ghz, 1e-6, t_cmb_k)
    assert delta[0] < 0 < delta[2]
    assert abs(delta[1]) < 1e-3 * abs(delta[0])


def test_ifg_round_trip_and_output_signs(config):
    t_cmb_k = config.sky.t_cmb_k
    freq_ghz = frequency_grid_ghz(config.instrument)
    delta_m = opd_samples_m(config.ifg)
    sky = blackbody_mjy_sr(freq_ghz, t_cmb_k) + y_distortion_mjy_sr(
        freq_ghz, 1e-5, t_cmb_k
    )
    bbir = blackbody_mjy_sr(freq_ghz, 2.725)
    dus = {du.name: du for du in detector_units(config.instrument)}
    for name in ("LF-L", "HF-R"):
        du = dus[name]
        ifg = interferogram(freq_ghz, sky, bbir, du, delta_m)
        recovered = differential_spectrum(freq_ghz, ifg, du, delta_m)
        mask = du.band_mask(freq_ghz)
        expected = du.output_sign * (sky - bbir)
        assert np.allclose(recovered[mask], expected[mask], rtol=0, atol=1e-6)
        assert np.all(np.isnan(recovered[~mask]))
    # L and R modulated parts are opposite
    ifg_l = interferogram(freq_ghz, sky, bbir, dus["LF-L"], delta_m)
    ifg_r = interferogram(freq_ghz, sky, bbir, dus["LF-R"], delta_m)
    total = ifg_l + ifg_r
    assert np.allclose(total, total[0], rtol=1e-12)


def test_mirror_positions_follow_opd_factor(config):
    # Step2 formula cos(4π ν δ / c): δ_max = c / (4 Δν) = 5 mm for 15 GHz bins
    assert np.max(mirror_positions_m(config.ifg)) == 0.005
    assert np.allclose(mirror_positions_m(config.ifg) * 2, opd_samples_m(config.ifg))
