from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.operations import cycle_timeline, num_cycles
from fossil_sim.orbit import L2Orbit
from fossil_sim.scanning import cycle_pointing, ring_spin_axes

BEAM_FWHM_DEG = 1.6  # Step2


def angle_deg(a, b):
    return np.rad2deg(np.arccos(np.clip(np.sum(a * b, axis=-1), -1, 1)))


@pytest.fixture(scope="module")
def cycle0(config, orbit):
    tl = cycle_timeline(config, 0)
    return tl, cycle_pointing(config, orbit, tl)


def test_boresight_geometry(config, cycle0):
    tl, pointing = cycle0
    assert np.allclose(np.linalg.norm(pointing.boresight, axis=1), 1.0)
    axis = pointing.spin_axis[tl.ring_in_cycle]
    assert np.allclose(angle_deg(pointing.boresight, axis), 90.0)
    # Spin phase 0 is the circle point closest to the north ecliptic pole
    first_scan = pointing.boresight[tl.scan_in_cycle == 0]
    assert np.argmax(first_scan[:, 2]) in (0, 719)


def test_scan_sampling(config, cycle0):
    tl, pointing = cycle0
    scan = pointing.boresight[tl.scan_in_cycle == 0]
    step_deg = angle_deg(scan[1:], scan[:-1])
    assert np.allclose(step_deg, 0.5)  # 0.25 deg/s × 2 s
    assert BEAM_FWHM_DEG / step_deg.mean() == pytest.approx(3.2)
    # The 23 Scans of a Ring repeat the same circle
    ring0 = pointing.boresight[tl.ring_id == 0].reshape(23, 720, 3)
    assert np.allclose(ring0, ring0[0])


def test_ring_steps(config, orbit):
    cycle_s = config.timeline.cycle_duration_s
    for cycle_id in (0, 120):  # near perihelion and near aphelion
        axes = ring_spin_axes(config, orbit, cycle_id)
        next_axes = ring_spin_axes(config, orbit, cycle_id + 1)
        steps_deg = angle_deg(np.vstack([axes[1:], next_axes[:1]]), axes)
        sun_lon_rad = orbit.anti_sun_longitude_rad(
            cycle_id * cycle_s + np.array([0, cycle_s])
        )
        solar_motion_deg = np.rad2deg(np.diff(sun_lon_rad)[0])
        # ~1.44 deg per 34.8 h, in 3 equal steps of ~0.48 deg ("~0.5 deg", Step2)
        assert 1.38 < solar_motion_deg < 1.48
        assert np.allclose(steps_deg, solar_motion_deg / 3, atol=1e-3)


@pytest.mark.parametrize(
    "rule, max_saa_deg", [("equal_per_cycle", 0.36), ("anchor_ring_start", 0.72)]
)
def test_saa_stays_small(config, orbit, rule, max_saa_deg):
    cfg = replace(config, scanning=replace(config.scanning, ring_stepping=rule))
    for cycle_id in range(0, num_cycles(cfg), 97):
        saa_deg = np.rad2deg(
            cycle_pointing(cfg, orbit, cycle_timeline(cfg, cycle_id)).saa_rad
        )
        assert saa_deg.max() < max_saa_deg


def test_fixed_stepping_drifts(config, orbit):
    cfg = replace(config, scanning=replace(config.scanning, ring_stepping="fixed"))
    saa_deg = np.rad2deg(cycle_pointing(cfg, orbit, cycle_timeline(cfg, 250)).saa_rad)
    assert saa_deg.min() > 10.0  # ~1 yr of 0.06 deg/cycle drift


def test_precession_cone(config):
    sc = replace(
        config.scanning, precession_amplitude_deg=5.0, precession_period_days=30.0
    )
    cfg = replace(config, scanning=sc)
    orbit = L2Orbit(cfg.orbit, cfg.simulation.mission_start_time, 90 * 86400.0)
    for cycle_id in range(0, 60, 7):
        axes = ring_spin_axes(cfg, orbit, cycle_id)
        nominal = ring_spin_axes(config, orbit, cycle_id)
        assert np.allclose(angle_deg(axes, nominal), 5.0)
