from dataclasses import replace

import healpy as hp
import numpy as np
import pytest

from fossil_sim.input_sky import (
    C_LIGHT_KM_S,
    boosted_isotropic_background_mjy_sr,
    dipole_beam_factor,
    isotropic_background_mjy_sr,
    solar_velocity_km_s,
)
from fossil_sim.spectra import blackbody_mjy_sr

FREQ_GHZ = np.arange(60.0, 2000.0, 15.0)


def blackbody_only(sky):
    return replace(sky, y=0.0, mu=0.0, cib_amplitude=0.0, co_amplitude=0.0)


def test_boost_of_blackbody_is_exact(config):
    sky = blackbody_only(config.sky)
    rng = np.random.default_rng(1)
    direction = rng.normal(size=(50, 3))
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    velocity_km_s = np.array([300.0, -150.0, 80.0])
    observed = boosted_isotropic_background_mjy_sr(
        FREQ_GHZ, direction, velocity_km_s, sky
    )
    beta = velocity_km_s / C_LIGHT_KM_S
    gamma = 1 / np.sqrt(1 - beta @ beta)
    t_obs_k = sky.t_cmb_k / (gamma * (1 - direction @ beta))
    expected = np.stack([blackbody_mjy_sr(FREQ_GHZ, t) for t in t_obs_k])
    assert np.allclose(observed, expected, rtol=1e-12, atol=0)


def test_no_velocity_gives_rest_frame(config):
    direction = np.array([[0.0, 0.0, 1.0], [1.0, 0.0, 0.0]])
    observed = boosted_isotropic_background_mjy_sr(
        FREQ_GHZ, direction, np.zeros(3), config.sky
    )
    assert np.allclose(
        observed, isotropic_background_mjy_sr(FREQ_GHZ, config.sky), rtol=1e-14
    )


def test_beam_factor_scales_dipole(config):
    sky = blackbody_only(config.sky)
    direction = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    velocity_km_s = np.array([370.0, 0.0, 0.0])
    b1 = dipole_beam_factor(1.6)
    plain = boosted_isotropic_background_mjy_sr(FREQ_GHZ, direction, velocity_km_s, sky)
    smoothed = boosted_isotropic_background_mjy_sr(
        FREQ_GHZ, direction, velocity_km_s, sky, b1
    )
    assert np.allclose(
        (smoothed[0] - smoothed[1]) / (plain[0] - plain[1]), b1, rtol=1e-10
    )


def test_solar_velocity_in_ecliptic(config):
    velocity = solar_velocity_km_s(config.sky)
    assert np.linalg.norm(velocity) == pytest.approx(369.816)
    # Independent check: healpy Galactic -> ecliptic rotation
    galactic = hp.ang2vec(np.deg2rad(90 - 48.253), np.deg2rad(264.021))
    expected = hp.Rotator(coord=["G", "E"])(galactic)
    assert np.allclose(velocity / np.linalg.norm(velocity), expected, atol=1e-5)
    lat_deg = np.rad2deg(np.arcsin(velocity[2] / np.linalg.norm(velocity)))
    assert lat_deg == pytest.approx(-11.15, abs=0.02)  # (λ, β) ≈ (171.66°, −11.15°)
