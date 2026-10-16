from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.input_sky import ZodiacalLight
from fossil_sim.scanning import ring_spin_axes
from fossil_sim.simulation import build_ifg_table

FREQ_GHZ = np.array([100.0, 545.0, 1500.0])


@pytest.fixture(scope="module")
def zodi_and_table(config, orbit):
    zodi = ZodiacalLight(config, orbit, FREQ_GHZ)
    table = build_ifg_table(config, orbit, 5)
    return zodi, table


def test_ring_interpolation_matches_direct_evaluation(config, orbit, zodi_and_table):
    zodi, table = zodi_and_table
    interpolated = zodi.for_cycle(5, ring_spin_axes(config, orbit, 5), table)
    # Mid-Ring IFGs (largest interpolation error) and Recycling IFGs, every 7th phase
    mid_ring = (table["scan_id"] - 5 * 87) == 11
    recycling = (table["scan_id"] - 5 * 87) == 80
    sel = np.flatnonzero((mid_ring | recycling) & (table["ifg_in_scan"] % 7 == 0))
    direct = zodi.evaluate(table["boresight"][sel], table["time_mid_s"][sel])
    assert np.all(direct > 0)
    relative = np.abs(interpolated[sel] - direct) / direct
    assert relative.max() < 1e-4  # measured ~2e-5


def test_ecliptic_brighter_than_pole(config, orbit):
    zodi = ZodiacalLight(config, orbit, FREQ_GHZ)
    anti_sun = orbit.anti_sun_direction(np.array([0.0]))[0]
    side = np.cross([0.0, 0.0, 1.0], anti_sun)  # 90° elongation on the ecliptic
    values = zodi.evaluate(np.array([side, [0.0, 0.0, 1.0]]), np.zeros(2))
    assert np.all(values[0] > 2 * values[1])
    assert np.all(np.diff(values[0]) > 0)  # rising with frequency up to 1.5 THz


def test_disabled_zodi(config, orbit):
    from fossil_sim.simulation import SignalModel

    cfg = replace(config, zodi=replace(config.zodi, enabled=False))
    assert SignalModel(cfg, orbit).zodi_for_cycle(0, None) is None
