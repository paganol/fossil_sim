from dataclasses import replace
from pathlib import Path

import healpy as hp
import numpy as np
import pytest

from fossil_sim.config import load_config
from fossil_sim.orbit import L2Orbit, icrs_to_ecliptic_matrix

CONFIG = load_config(Path(__file__).parents[1] / "configs" / "baseline.toml")
DURATION_S = 400 * 86400.0


@pytest.fixture(scope="module")
def orbit():
    return L2Orbit(CONFIG.orbit, CONFIG.simulation.mission_start_time, DURATION_S)


def test_icrs_to_ecliptic_matches_healpy():
    # healpy "C" -> "E" rotation; agreement to ~0.1 arcsec (frame-bias level)
    healpy_matrix = hp.Rotator(coord=["C", "E"]).mat
    assert np.allclose(icrs_to_ecliptic_matrix(), healpy_matrix, atol=1e-6)


def test_velocity_is_derivative_of_position(orbit):
    time_s = np.linspace(1e5, DURATION_S - 1e5, 50)
    step_s = 60.0
    pos_plus, _ = orbit.position_velocity(time_s + step_s)
    pos_minus, _ = orbit.position_velocity(time_s - step_s)
    _, vel_km_s = orbit.position_velocity(time_s)
    numerical_km_s = (pos_plus - pos_minus) / (2 * step_s)
    assert np.allclose(numerical_km_s, vel_km_s, atol=1e-6)  # 1 mm/s


def test_l2_distance_and_speed(orbit):
    no_lissajous = replace(
        CONFIG.orbit, lissajous_radius1_km=0.0, lissajous_radius2_km=0.0
    )
    plain = L2Orbit(no_lissajous, CONFIG.simulation.mission_start_time, DURATION_S)
    time_s = np.linspace(0, DURATION_S, 20)
    pos_km, vel_km_s = plain.position_velocity(time_s)
    dist_au = np.linalg.norm(pos_km, axis=1) / 1.495978707e8
    assert np.all((dist_au > 0.98) & (dist_au < 1.03))
    speed_km_s = np.linalg.norm(vel_km_s, axis=1)
    assert np.all((speed_km_s > 29.5) & (speed_km_s < 31.0))
    # The Lissajous term stays within its radii
    full_pos_km, _ = orbit.position_velocity(time_s)
    offset_km = np.linalg.norm(full_pos_km - pos_km, axis=1)
    assert offset_km.max() <= np.hypot(244450.0, 137388.0) + 1.0


def test_anti_sun_longitude_rate(orbit):
    time_s = np.linspace(0.0, 365.25 * 86400, 1000)
    longitude_rad = orbit.anti_sun_longitude_rad(time_s)
    # One Julian year: ~one full turn (Earth sidereal period is 365.256 d)
    assert np.rad2deg(longitude_rad[-1] - longitude_rad[0]) == pytest.approx(
        360.0, abs=0.2
    )
