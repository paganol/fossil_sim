import numpy as np
import pytest

from fossil_sim.simulation import build_ifg_table, run_simulation


def test_ifg_table(config, orbit):
    table = build_ifg_table(config, orbit, 0)
    assert table.size == 87 * 720
    assert np.all(np.diff(table["time_mid_s"]) == 2.0)
    assert np.all(table["ring_id"][: 69 * 720] >= 0)


def test_short_run_statistics(config):
    result = run_simulation(config, num_of_cycles=3)
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
    result = run_simulation(config, num_of_cycles=num_of_cycles)
    maps = result.hit_maps
    assert maps.coverage_fraction(64, np.array([183 * 86400.0]))[0] >= 0.999
    assert maps.coverage_fraction(64, np.array([150 * 86400.0]))[0] < 0.9
    for mode_row in (1, 2):  # B01, B10 separately
        assert np.mean(maps.hits[64][mode_row] > 0) >= 0.999
    assert result.max_saa_deg < 0.36
