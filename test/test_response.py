from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.instrument import detector_units, frequency_grid_ghz
from fossil_sim.noise import sky_estimate_noise_variance
from fossil_sim.operations import mode_sky_ports
from fossil_sim.response import (
    calibration_solution,
    gain_drift,
    generation_response,
    nominal_response,
    true_response,
)
from fossil_sim.simulation import SignalModel, build_ifg_table, run_simulation
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
    if simulate:
        ideal.ifg_response = generation_response(
            ideal.config, ideal.freq_ghz, ideal.response, ideal.sky_freq_ghz, ideal.dus
        )
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
        ("gain_drift_sigma", -1e-3),
        ("gain_drift_fknee_mhz", -1.0),
        ("gain_drift_alpha", 0.0),
        ("gain_drift_fmin_hz", 0.0),
        ("gain_drift_fknee_mhz", 1e-3),  # below gain_drift_fmin_hz (1e-5 Hz)
        ("gain_drift_fknee_mhz", 250.0),  # at the IFG Nyquist frequency (0.25 Hz)
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


def drift_config(config, sigma=1e-3, fknee_mhz=0.0, alpha=1.0):
    return replace(
        config,
        response=replace(
            config.response,
            gain_drift_sigma=sigma,
            gain_drift_fknee_mhz=fknee_mhz,
            gain_drift_alpha=alpha,
        ),
    )


def test_gain_drift_keyed_and_scaled(config):
    assert gain_drift(config, 0, 0, 100) is None
    cfg = drift_config(config, sigma=2e-3)
    drift = gain_drift(cfg, 1, 7, 50000)
    np.testing.assert_array_equal(drift, gain_drift(cfg, 1, 7, 50000))
    assert not np.array_equal(drift, gain_drift(cfg, 1, 8, 50000))
    assert not np.array_equal(drift, gain_drift(cfg, 2, 7, 50000))
    assert drift.mean() == pytest.approx(1.0, abs=1e-4)
    assert drift.std() == pytest.approx(2e-3, rel=0.02)


@pytest.mark.parametrize("alpha", [1.0, 2.0])
def test_gain_drift_spectrum(config, alpha):
    from scipy.signal import welch

    cfg = drift_config(config, sigma=1e-3, fknee_mhz=10.0, alpha=alpha)
    rate_hz = 1 / cfg.timeline.ifg_duration_s
    delta = (gain_drift(cfg, 0, 0, 2**20) - 1) / 1e-3
    freq_hz, psd = welch(delta, fs=rate_hz, nperseg=2**14)
    psd = psd * rate_hz / 2  # one-sided PSD of unit white noise = 2 / rate
    for f in (1e-3, 3e-3, 1e-2, 3e-2, 0.1):
        i = np.argmin(np.abs(freq_hz - f))
        expected = ((f**2 + 1e-4) / f**2) ** (alpha / 2)
        assert psd[i] == pytest.approx(expected, rel=0.1)


@pytest.mark.parametrize("simulate", [False, True])
def test_gain_drift_scales_signal_not_noise(config, orbit, table, simulate):
    cfg = drift_config(response_config(config, simulate=simulate), sigma=0.01)
    cfg = replace(cfg, noise=replace(cfg.noise, enabled=True))
    model = SignalModel(cfg, orbit)
    drift = {
        name: g[: len(table)]
        for name, g in model.gain_drift_for_cycle(
            0, build_ifg_table(cfg, orbit, 0)
        ).items()
    }
    noise = model.noise_for_cycle(0, table)
    noise = {name: n[0 : len(table)] for name, n in noise.items()}
    clean = model.compute(table)[0]
    noisy = model.compute(table, noise=noise)[0]
    drifting = model.compute(table, noise=noise, drift=drift)[0]
    for name, x in clean.items():
        expected = noisy[name] + (drift[name][:, None] - 1) * x
        np.testing.assert_allclose(
            drifting[name], expected, rtol=0, atol=1e-9 * np.abs(x).max()
        )


def test_gain_drift_in_tod(config, tmp_path):
    import h5py

    cfg = replace(
        response_config(config, gain_rms=0.0),
        output=replace(config.output, save_tod=True),
    )
    run_simulation(cfg, num_of_cycles=1, output_dir=tmp_path / "steady")
    run_simulation(
        drift_config(cfg, sigma=1e-3, fknee_mhz=10.0),
        num_of_cycles=1,
        output_dir=tmp_path / "drift",
    )
    name = "cycle_00000.h5"
    with (
        h5py.File(tmp_path / "steady" / "tod" / name) as steady,
        h5py.File(tmp_path / "drift" / "tod" / name) as drifting,
    ):
        assert "gain_drift" not in steady["response"]
        for du in detector_units(cfg.instrument):
            drift = drifting[f"response/gain_drift/{du.name}"][:]
            assert drift.size == steady[f"spectra/{du.name}"].shape[0]
            np.testing.assert_allclose(
                drifting[f"spectra/{du.name}"][:],
                steady[f"spectra/{du.name}"][:] * drift[:, None],
                rtol=1e-12,
            )
