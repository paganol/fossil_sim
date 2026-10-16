"""Cross-check of the pointing against litebird_sim (optional oracle).

litebird_sim composes R_z(λ) R_x(ψ) R_y(π/2 − α) R_z(ω t): for a frozen
Ring we pass λ = Spin-axis longitude, α = Precession amplitude and
ψ = Precession phase − π/2. Its spin phase 0 points to the south side of the
circle, ours to the north, hence ω t = φ + π.
"""

from dataclasses import replace

import numpy as np
import pytest

lbs = pytest.importorskip("litebird_sim")

from fossil_sim.operations import cycle_timeline
from fossil_sim.scanning import cycle_pointing, ring_axis_longitudes_rad, ring_spin_axes

pytestmark = pytest.mark.heavy


def litebird_vector(longitude_rad, alpha_rad, psi_rad, spin_angle_rad, vector):
    quat = np.array(lbs.quat_rotation_z(spin_angle_rad), dtype=np.float64)
    lbs.quat_left_multiply(quat, *lbs.quat_rotation_y(np.pi / 2 - alpha_rad))
    lbs.quat_left_multiply(quat, *lbs.quat_rotation_x(psi_rad))
    lbs.quat_left_multiply(quat, *lbs.quat_rotation_z(longitude_rad))
    result = np.empty(3)
    lbs.rotate_vector(result, *quat, np.asarray(vector, dtype=np.float64))
    return result


def test_boresight_matches_litebird_sim(config, orbit):
    timeline = cycle_timeline(config, 3)
    pointing = cycle_pointing(config, orbit, timeline)
    longitude_rad = ring_axis_longitudes_rad(config, orbit, 3)
    beta_rad = np.deg2rad(config.scanning.spin_boresight_angle_deg)
    detector = [np.sin(beta_rad), 0.0, np.cos(beta_rad)]
    for i in np.linspace(0, 69 * 720 - 1, 300).astype(int):
        phase_rad = (
            2 * np.pi * (timeline.ifg_in_scan[i] + 0.5) / config.timeline.ifgs_per_scan
        )
        expected = litebird_vector(
            longitude_rad[timeline.ring_in_cycle[i]],
            0.0,
            0.0,
            phase_rad + np.pi,
            detector,
        )
        assert np.allclose(pointing.boresight[i], expected, atol=1e-9)  # < 1 µrad


def test_precessing_axis_matches_litebird_sim(config, orbit):
    sc = replace(
        config.scanning, precession_amplitude_deg=5.0, precession_period_days=20.0
    )
    cfg = replace(config, scanning=sc)
    tl = cfg.timeline
    alpha_rad = np.deg2rad(sc.precession_amplitude_deg)
    for cycle_id in (0, 5, 11):
        axes = ring_spin_axes(cfg, orbit, cycle_id)
        longitude_rad = ring_axis_longitudes_rad(cfg, orbit, cycle_id)
        ring_mid_s = (
            cycle_id * tl.cycle_duration_s
            + (np.arange(3) + 0.5) * tl.scans_per_ring * tl.scan_duration_s
        )
        phase_rad = 2 * np.pi * ring_mid_s / (sc.precession_period_days * 86400)
        for ring in range(3):
            expected = litebird_vector(
                longitude_rad[ring],
                alpha_rad,
                phase_rad[ring] - np.pi / 2,
                0.0,
                [0, 0, 1.0],
            )
            assert np.allclose(axes[ring], expected, atol=1e-9)
