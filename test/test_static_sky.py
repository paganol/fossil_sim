from dataclasses import replace

import healpy as hp
import numpy as np
import pytest

from fossil_sim.config import config_from_dict
from fossil_sim.input_sky import StaticSky, orbital_doppler_mjy_sr

FREQ_GHZ = np.array([100.0, 200.0, 400.0])
NSIDE = 16


def small_config(static_sky_config, tmp_path, **changes):
    return replace(
        static_sky_config,
        nside=NSIDE,
        lmax=3 * NSIDE - 1,
        cache_dir=str(tmp_path),
        **changes,
    )


def synthetic_emission(galactic_monopole=2.0, extragalactic_monopole=5.0):
    """Galactic: monopole + cos(θ_gal) pattern scaling as ν; extragalactic: monopole."""
    nside_eval = 2 * NSIDE
    z_gal = hp.pix2vec(nside_eval, np.arange(hp.nside2npix(nside_eval)))[2]

    def emission(freq_ghz):
        galactic = (galactic_monopole + z_gal) * freq_ghz / 100.0
        return galactic, np.full_like(z_gal, extragalactic_monopole)

    return emission


def random_directions(n=200, seed=3):
    v = np.random.default_rng(seed).normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def test_monopoles_rotation_and_beam(config, tmp_path):
    cfg = small_config(config.static_sky, tmp_path)
    sky = StaticSky(cfg, FREQ_GHZ, 1.6, emission=synthetic_emission())
    direction = random_directions()
    sampled = sky.sample(direction)
    # Expected: Galactic monopole kept, extragalactic monopole removed, dipole × b₁,
    # Galactic north pole rotated to ecliptic coordinates
    z_gal_in_ecl = hp.Rotator(coord=["G", "E"])(np.array([0.0, 0.0, 1.0]))
    b1 = hp.gauss_beam(np.deg2rad(1.6), lmax=1)[1]
    expected = (2.0 + b1 * direction @ z_gal_in_ecl)[:, None] * FREQ_GHZ / 100.0
    assert np.allclose(sampled, expected, atol=1e-6)


def test_cache_is_reused_and_keyed(config, tmp_path):
    cfg = small_config(config.static_sky, tmp_path)
    first = StaticSky(cfg, FREQ_GHZ, 1.6, emission=synthetic_emission())
    # A different emission is ignored: the cached cube is read back
    second = StaticSky(
        cfg, FREQ_GHZ, 1.6, emission=synthetic_emission(galactic_monopole=9.0)
    )
    assert np.array_equal(first.alm, second.alm)
    other = StaticSky(
        replace(cfg, models=("d1",)), FREQ_GHZ, 1.6, emission=synthetic_emission(9.0)
    )
    assert other.path != first.path


def test_nearest_sampling_close_to_alm(config, tmp_path):
    cfg = small_config(config.static_sky, tmp_path)
    exact = StaticSky(cfg, FREQ_GHZ, 1.6, emission=synthetic_emission())
    nearest = StaticSky(
        replace(cfg, sampling="nearest"), FREQ_GHZ, 1.6, emission=synthetic_emission()
    )
    direction = random_directions()
    pixel_rad = hp.nside2resol(NSIDE)
    assert np.allclose(
        nearest.sample(direction), exact.sample(direction), atol=4 * pixel_rad
    )


def test_orbital_doppler_power_law():
    freq_ghz = np.geomspace(100.0, 1000.0, 200)
    alpha = 1.7
    intensity = (freq_ghz / 100.0) ** alpha * np.ones((2, 1))
    direction = np.array([[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]])
    velocity_km_s = np.array([[30.0, 0.0, 0.0]] * 2)
    delta = orbital_doppler_mjy_sr(freq_ghz, intensity, direction, velocity_km_s)
    beta = 30.0 / 299792.458
    inner = slice(1, -1)  # np.gradient is first order at the edges
    assert np.allclose(
        delta[0, inner], beta * (3 - alpha) * intensity[0, inner], rtol=1e-3
    )
    assert np.allclose(delta[1], -delta[0])


def test_forbidden_models_rejected(config):
    values = config.to_dict()
    for model in ("dip1", "co2"):
        values["static_sky"]["models"] = [model]
        with pytest.raises(ValueError):
            config_from_dict(values)


@pytest.mark.heavy
def test_pysm3_static_sky(config, tmp_path):
    """Real pysm3 models (downloads templates): dust dominates at high frequency."""
    cfg = small_config(config.static_sky, tmp_path, models=("d1", "s1"))
    sky = StaticSky(cfg, np.array([100.0, 1000.0]), 1.6)
    plane = hp.Rotator(coord=["G", "E"])(np.array([[1.0, 0.0, 0.0]]).T).T
    values = sky.sample(plane)
    assert values[0, 1] > 10 * values[0, 0] > 0
