import numpy as np
import pytest

from fossil_sim.ifg import (
    _interpolated_modulation,
    differential_spectrum,
    interferogram,
    mirror_positions_m,
    opd_samples_m,
    sample_weights,
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


@pytest.mark.parametrize("frequency_ghz", [100.0, 300.0])
def test_single_bin_interferogram_uses_explicit_width(config, frequency_ghz):
    freq_ghz = np.array([frequency_ghz])
    du = next(
        du for du in detector_units(config.instrument) if du.band_mask(freq_ghz)[0]
    )
    opd_m = opd_samples_m(config.ifg)
    sky, bbir = np.array([10.0]), np.array([3.0])
    width_ghz, efficiency = 15.0, 0.25
    ifg = interferogram(
        freq_ghz, sky, bbir, du, opd_m, efficiency, freq_step_ghz=width_ghz
    )
    from fossil_sim.constants import C_LIGHT_M_S

    expected = (
        efficiency
        * width_ghz
        * (13.0 + 7.0 * np.cos(2 * np.pi * frequency_ghz * 1e9 * opd_m / C_LIGHT_M_S))
    )
    np.testing.assert_allclose(ifg, expected, atol=1e-12)
    recovered = differential_spectrum(
        freq_ghz, ifg, du, opd_m, efficiency, freq_step_ghz=width_ghz
    )
    np.testing.assert_allclose(recovered, [7.0], atol=1e-12)
    with pytest.raises(ValueError, match="freq_step_ghz is required"):
        interferogram(freq_ghz, sky, bbir, du, opd_m)


@pytest.mark.parametrize("width_ghz", [np.nan, 10.0])
def test_explicit_width_must_match_grid(config, width_ghz):
    freq_ghz = frequency_grid_ghz(config.instrument)
    du = detector_units(config.instrument)[0]
    opd_m = opd_samples_m(config.ifg)
    port = np.ones(freq_ghz.size)
    with pytest.raises(ValueError, match="freq_step_ghz"):
        interferogram(freq_ghz, port, port, du, opd_m, freq_step_ghz=width_ghz)


def test_single_bin_rejects_nan_width(config):
    freq_ghz = np.array([300.0])
    du = detector_units(config.instrument)[1]
    with pytest.raises(ValueError, match="freq_step_ghz"):
        interferogram(
            freq_ghz,
            np.ones(1),
            np.ones(1),
            du,
            opd_samples_m(config.ifg),
            freq_step_ghz=np.nan,
        )


@pytest.mark.parametrize(
    "fmin, fmax, step, size", [(1.0, 1.3, 0.1, 3), (30.0, 2000.0, 15.0, 132)]
)
def test_frequency_grid_excludes_maximum(config, fmin, fmax, step, size):
    from dataclasses import replace

    inst = replace(
        config.instrument, freq_min_ghz=fmin, freq_max_ghz=fmax, freq_step_ghz=step
    )
    grid = frequency_grid_ghz(inst)
    assert grid.size == size
    assert grid[-1] < fmax
    np.testing.assert_allclose(grid, fmin + step * np.arange(size), rtol=0, atol=0)


@pytest.mark.parametrize("zeros", [False, True])
def test_interpolated_modulation_equals_the_direct_product(zeros):
    """Grouped and constant-window shortcuts equal the product on every window point."""
    rng = np.random.default_rng(1)
    n_ifg, n_freq, n_samples, grid = 12, 40, 64, 8
    weights = sample_weights(n_samples, grid)
    if zeros:  # rows starting or ending with an exact zero form their own groups
        weights[3, np.flatnonzero(weights[3])[0]] = 0.0
        weights[10, np.flatnonzero(weights[10])[-1]] = 0.0
    difference = rng.standard_normal((grid + 3, n_ifg, n_freq))  # window-major
    difference[:, ::3] = difference[:1, ::3]  # constant windows
    modulation = rng.standard_normal((n_samples, n_freq))
    direct = np.einsum("jw,wij->ij", weights, difference @ modulation.T)
    detected = _interpolated_modulation(difference, modulation, weights)
    np.testing.assert_allclose(detected, direct, rtol=1e-12, atol=1e-12)
    # flags from the caller; a varying window flagged as varying is always exact
    flags = np.zeros(n_ifg, dtype=bool)
    flags[::6] = True
    given = _interpolated_modulation(difference, modulation, weights, flags)
    np.testing.assert_allclose(given, direct, rtol=1e-12, atol=1e-12)
