from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.instrument import detector_units, frequency_grid_ghz
from fossil_sim.noise import sky_estimate_noise_variance
from fossil_sim.operations import mode_sky_ports
from fossil_sim.response import (
    calibration_solution,
    nominal_response,
    true_response,
)
from fossil_sim.simulation import SignalModel, build_ifg_table
from fossil_sim.spectra import blackbody_mjy_sr

GAIN = (1.02, 0.97, 1.05, 0.99)
EMISSIVITY = (1e-3, -2e-3, -1e-3, 5e-4)


def response_config(config, solution="nominal", gain_rms=0.01, simulate=False):
    return replace(
        config,
        static_sky=replace(config.static_sky, models=()),
        zodi=replace(config.zodi, enabled=False),
        noise=replace(config.noise, enabled=False),
        ifg=replace(config.ifg, simulate=simulate, sky_variation=False),
        response=replace(
            config.response,
            gain=GAIN,
            gain_rms=gain_rms,
            offset_emissivity=EMISSIVITY,
        ),
        calibration=replace(config.calibration, solution=solution),
    )


@pytest.fixture(scope="module")
def table(config, orbit):
    table = build_ifg_table(config, orbit, 0)
    # science IFGs in both Modes and Recycling (B00) IFGs
    return np.concatenate([table[90:110], table[8580:8600], table[-20:]])


def test_baseline_is_the_ideal_instrument(config):
    freq_ghz = frequency_grid_ghz(config.instrument)
    dus = detector_units(config.instrument)
    response = true_response(config, freq_ghz, dus)
    nominal = nominal_response(freq_ghz, dus)
    for du in dus:
        assert np.all(response.gain[du.name] == 1.0)
        assert np.all(response.offset_mjy_sr[du.name] == 0.0)
        assert np.all(nominal.gain[du.name] == 1.0)


def test_true_response_gain_and_offset(config):
    cfg = response_config(config, gain_rms=0.05)
    freq_ghz = frequency_grid_ghz(cfg.instrument)
    dus = detector_units(cfg.instrument)
    response = true_response(cfg, freq_ghz, dus)
    again = true_response(cfg, freq_ghz, dus)
    other_seed = true_response(
        replace(cfg, simulation=replace(cfg.simulation, master_seed=12345)),
        freq_ghz,
        dus,
    )
    for i, du in enumerate(dus):
        band_ghz = freq_ghz[du.band_mask(freq_ghz)]
        gain = response.gain[du.name]
        assert gain.shape == band_ghz.shape
        np.testing.assert_array_equal(gain, again.gain[du.name])  # fixed per seed
        assert not np.array_equal(gain, other_seed.gain[du.name])
        assert 0.02 < np.std(gain / GAIN[i]) < 0.08 or band_ghz.size < 20
        np.testing.assert_allclose(
            response.offset_mjy_sr[du.name],
            EMISSIVITY[i] * blackbody_mjy_sr(band_ghz, 2.725),
            rtol=1e-14,
        )
    # the DUs draw independent fluctuations
    assert not np.allclose(
        response.gain["HF-L"] / GAIN[1], response.gain["HF-R"] / GAIN[3]
    )


def test_apply_and_calibrate_are_inverse(config):
    cfg = response_config(config)
    freq_ghz = frequency_grid_ghz(cfg.instrument)
    dus = detector_units(cfg.instrument)
    response = true_response(cfg, freq_ghz, dus)
    rng = np.random.default_rng(0)
    x = {
        du.name: rng.standard_normal((3, int(du.band_mask(freq_ghz).sum())))
        for du in dus
    }
    measured = response.apply(x)
    for du in dus:
        np.testing.assert_allclose(
            measured[du.name],
            response.gain[du.name] * (x[du.name] + response.offset_mjy_sr[du.name]),
        )
    back = response.calibrate(measured)
    for name, value in x.items():
        np.testing.assert_allclose(back[name], value, atol=1e-12)


def test_calibration_solution_choice(config):
    freq_ghz = frequency_grid_ghz(config.instrument)
    dus = detector_units(config.instrument)
    cfg = response_config(config, "true")
    true = calibration_solution(cfg, freq_ghz, dus)
    expected = true_response(cfg, freq_ghz, dus)
    nominal = calibration_solution(response_config(config, "nominal"), freq_ghz, dus)
    for du in dus:
        np.testing.assert_array_equal(true.gain[du.name], expected.gain[du.name])
        np.testing.assert_array_equal(
            true.offset_mjy_sr[du.name], expected.offset_mjy_sr[du.name]
        )
        assert np.all(nominal.gain[du.name] == 1.0)
        assert np.all(nominal.offset_mjy_sr[du.name] == 0.0)


@pytest.mark.parametrize("simulate", [False, True])
def test_true_calibration_recovers_the_ideal_sky(config, orbit, table, simulate):
    ideal = SignalModel(response_config(config, gain_rms=0.0, simulate=simulate), orbit)
    ideal.response = nominal_response(ideal.freq_ghz, ideal.dus)
    model = SignalModel(response_config(config, "true", simulate=simulate), orbit)
    reference, sky_reference, _ = ideal.compute(table)
    measured, sky, _ = model.compute(table)
    for name, x in reference.items():
        expected = model.response.gain[name] * (x + model.response.offset_mjy_sr[name])
        np.testing.assert_allclose(
            measured[name], expected, rtol=0, atol=1e-9 * np.abs(expected).max()
        )
    science = ~np.isnan(sky_reference[:, 0])
    np.testing.assert_allclose(sky[science], sky_reference[science], rtol=1e-9)
    assert np.all(np.isnan(sky[~science]))


def test_nominal_calibration_bias(config, orbit, table):
    """Unit gain, offset only: Ŝ − S = σ (o_L − o_R) / 2 per bin covered by L and R."""
    cfg = response_config(config, "nominal", gain_rms=0.0)
    cfg = replace(cfg, response=replace(cfg.response, gain=(1.0, 1.0, 1.0, 1.0)))
    model = SignalModel(cfg, orbit)
    ideal = SignalModel(replace(cfg, response=config.response), orbit)
    sky = model.compute(table)[1]
    sky_reference = ideal.compute(table)[1]
    port1_on_sky, port2_on_sky = mode_sky_ports(table["mode"])
    science = port1_on_sky ^ port2_on_sky
    sigma = np.where(port1_on_sky, 1.0, -1.0)[science, None]
    o = model.response.offset_mjy_sr
    lf = model.dus[0].band_mask(model.freq_ghz)
    expected = np.empty(model.freq_ghz.size)
    expected[lf] = (o["LF-L"] - o["LF-R"]) / 2
    expected[~lf] = (o["HF-L"] - o["HF-R"]) / 2
    np.testing.assert_allclose(
        sky[science] - sky_reference[science], sigma * expected, atol=1e-12
    )


def test_noise_variance_is_divided_by_the_calibration_gain(config):
    freq_ghz = frequency_grid_ghz(config.instrument)
    dus = detector_units(config.instrument)
    gain = {du.name: np.full(int(du.band_mask(freq_ghz).sum()), 2.0) for du in dus}
    nominal = sky_estimate_noise_variance(config, freq_ghz, dus)
    scaled = sky_estimate_noise_variance(config, freq_ghz, dus, gain=gain)
    np.testing.assert_allclose(scaled, nominal / 4)


@pytest.mark.parametrize(
    "key, value",
    [
        ("gain", [1.0, 1.0, 1.0]),
        ("gain", [1.0, 0.0, 1.0, 1.0]),
        ("offset_emissivity", [0.0]),
        ("gain_rms", -0.1),
        ("instrument_temperature_k", 0.0),
    ],
)
def test_invalid_response_is_rejected(config, key, value):
    with pytest.raises(ValueError):
        replace(
            config,
            response=replace(
                config.response,
                **{key: tuple(value) if isinstance(value, list) else value},
            ),
        )


def test_invalid_calibration_solution_is_rejected(config):
    with pytest.raises(ValueError):
        replace(config, calibration=replace(config.calibration, solution="bbir"))
