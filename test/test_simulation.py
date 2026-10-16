from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.orbit import L2Orbit
from fossil_sim.simulation import build_ifg_table, run_simulation


def test_ifg_table(config, orbit):
    table = build_ifg_table(config, orbit, 0)
    assert table.size == 87 * 720
    assert np.all(np.diff(table["time_mid_s"]) == 2.0)
    assert np.all(table["ring_id"][: 69 * 720] >= 0)


def test_short_run_statistics(config):
    result = run_simulation(config, num_of_cycles=3, spectral=False)
    lost = result.lost_fraction()
    assert result.num_ifgs == 3 * 87 * 720
    assert lost["recycling"] == pytest.approx(18 / 87)
    assert lost["mode_switch"] == pytest.approx((3 * 69 + 3) / (3 * 87 * 720))
    assert lost["saa_exceeded"] == 0.0
    assert (
        result.mode_counts["B10"] + result.mode_counts["B01"] == result.num_valid_ifgs
    )
    assert abs(result.mode_counts["B10"] - result.mode_counts["B01"]) <= 720
    for nside in config.maps.nsides:
        assert result.hit_maps.total(nside).sum() == result.num_valid_ifgs


def test_full_sky_in_six_months(config):
    num_of_cycles = int(np.ceil(183 * 86400 / config.timeline.cycle_duration_s))
    result = run_simulation(config, num_of_cycles=num_of_cycles, spectral=False)
    maps = result.hit_maps
    assert maps.coverage_fraction(64, np.array([183 * 86400.0]))[0] >= 0.999
    assert maps.coverage_fraction(64, np.array([150 * 86400.0]))[0] < 0.9
    for mode_row in (1, 2):  # B01, B10 separately
        assert np.mean(maps.hits[64][mode_row] > 0) >= 0.999
    assert result.max_saa_deg < 0.36


def test_spectral_maps_recover_kinematic_dipole(config, tmp_path):
    import h5py
    import healpy as hp

    from fossil_sim.input_sky import isotropic_background_mjy_sr, solar_velocity_km_s
    from fossil_sim.spectra import blackbody_mjy_sr

    cfg = replace(config, output=replace(config.output, save_tod=True))
    result = run_simulation(cfg, num_of_cycles=2, output_dir=tmp_path)
    maps = result.spectral_maps
    freq_ghz = maps.freq_ghz
    sky_map = maps.sky_map(64)
    seen = sky_map[0] != hp.UNSEEN
    assert np.array_equal(seen, result.hit_maps.total(64) > 0)

    # ΔT = ΔI / (∂B/∂T) at 150 GHz; fit monopole + dipole on the observed pixels
    i = np.argmin(np.abs(freq_ghz - 150.0))
    t_k = config.sky.t_cmb_k
    dbdt = (
        blackbody_mjy_sr(freq_ghz[i], t_k + 1e-4)
        - blackbody_mjy_sr(freq_ghz[i], t_k - 1e-4)
    ) / 2e-4
    delta_t = np.full(sky_map.shape[1], hp.UNSEEN)
    rest = isotropic_background_mjy_sr(freq_ghz[i], config.sky)
    delta_t[seen] = (sky_map[i, seen] - rest) / dbdt
    _, dipole_k = hp.fit_dipole(delta_t)
    orbit_vel = L2Orbit(
        config.orbit, config.simulation.mission_start_time, 2 * 87 * 1440.0
    )
    v_km_s = (
        solar_velocity_km_s(config.sky) + orbit_vel.position_velocity(87 * 1440.0)[1][0]
    )
    expected_k = t_k * np.linalg.norm(v_km_s) / 299792.458
    # Two cycles cover two thin strips: the fit is degenerate-free only to ~1%
    assert np.linalg.norm(dipole_k) == pytest.approx(expected_k, rel=0.02)

    # Diagnostic maps: B10 and B01 carry opposite sky sign on the same output
    raw = maps.raw_maps(64)  # (output, mode, freq, pix)
    both = (raw[0, 0, i] != hp.UNSEEN) & (raw[0, 1, i] != hp.UNSEEN)
    assert np.all(np.sign(raw[0, 0, i, both]) == -np.sign(raw[0, 1, i, both]))
    assert np.allclose(raw[1, :, i, both], -raw[0, :, i, both])

    with h5py.File(tmp_path / "tod" / "cycle_00001.h5") as f:
        assert f["ifg_table"].shape == (87 * 720,)
        assert f["spectra/HF-R"].shape[0] == 87 * 720
        assert "config_hash" in f.attrs
    assert len(maps.save(tmp_path / "maps")) == 4
