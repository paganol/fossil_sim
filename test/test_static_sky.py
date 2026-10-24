from dataclasses import replace

import healpy as hp
import numpy as np
import pytest

from fossil_sim.config import config_from_dict
from fossil_sim.input_sky import StaticSky, static_doppler_mjy_sr
from fossil_sim.simulation import SignalModel

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


def synthetic_emission(
    galactic_monopole=2.0, extragalactic_monopole=5.0, alpha=1.0, cosmological=0.0
):
    """Barycentric: (monopole + cos θ_gal) ν^alpha; cosmological: monopole + cosmological × cos θ_gal, ν^-1."""
    nside_eval = 2 * NSIDE
    z_gal = hp.pix2vec(nside_eval, np.arange(hp.nside2npix(nside_eval)))[2]

    def emission(freq_ghz):
        galactic = (galactic_monopole + z_gal) * (freq_ghz / 100.0) ** alpha
        cosmic = (extragalactic_monopole + cosmological * z_gal) * 100.0 / freq_ghz
        return galactic, cosmic

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


def test_doppler_slope_on_coarse_grid_with_edges(config, tmp_path):
    # Baseline-like 15 GHz grid from 30 GHz: K = 3I − dI/d ln ν = (3 − α) I
    # exactly for power laws, edges included (np.gradient gave −0.08 instead of
    # +1 times I at 30 GHz for α = 2)
    freq_ghz = np.array([30.0, 45.0, 60.0])
    cfg = small_config(config.static_sky, tmp_path)
    emission = synthetic_emission(alpha=2.0, cosmological=3.0)
    sky = StaticSky(cfg, freq_ghz, 1.6, emission=emission)
    direction = random_directions()
    k_barycentric, k_cosmological = sky.sample_doppler(direction)
    z = direction @ hp.Rotator(coord=["G", "E"])(np.array([0.0, 0.0, 1.0]))
    b1 = hp.gauss_beam(np.deg2rad(1.6), lmax=1)[1]
    galactic = (2.0 + b1 * z)[:, None] * (freq_ghz / 100.0) ** 2
    cosmic = (3.0 * b1 * z)[:, None] * 100.0 / freq_ghz  # monopole removed
    # central difference in ln ν: dI/d ln ν = α (1 + α²ε²/6) I, i.e. 1.3e-4 of K here
    assert np.allclose(k_barycentric, (3 - 2) * galactic, rtol=3e-4, atol=1e-9)
    assert np.allclose(k_cosmological, (3 + 1) * cosmic, rtol=3e-4, atol=1e-9)
    nearest = StaticSky(replace(cfg, sampling="nearest"), freq_ghz, 1.6, emission)
    pixel_rad = hp.nside2resol(NSIDE)
    for exact, near in zip(
        (k_barycentric, k_cosmological), nearest.sample_doppler(direction)
    ):
        assert np.allclose(near, exact, atol=40 * pixel_rad)


def test_doppler_velocities_per_frame(config, orbit, tmp_path):
    # Barycentric part moves with the orbital velocity only, cosmological part
    # with solar + orbital
    freq_ghz = np.array([100.0, 200.0])
    cfg = replace(
        config,
        static_sky=small_config(config.static_sky, tmp_path),
        zodi=replace(config.zodi, enabled=False),
        sky=replace(
            config.sky,
            include_cmb_blackbody=False,
            y=0.0,
            mu=0.0,
            cib_amplitude=0.0,
            co_amplitude=0.0,
        ),
    )

    class OnlyDoppler:
        def __init__(self, part):
            self.part = part

        def sample(self, direction):
            return np.zeros((len(direction), freq_ghz.size))

        def sample_doppler(self, direction):
            k = np.ones((len(direction), freq_ghz.size))
            return (k, 0 * k) if self.part == 0 else (0 * k, k)

    table = np.zeros(4, dtype=[("boresight", float, 3), ("time_mid_s", float)])
    table["boresight"] = random_directions(4)
    table["time_mid_s"] = 86400.0 * np.arange(4)
    orbital = orbit.position_velocity(table["time_mid_s"])[1]
    model = SignalModel(cfg, orbit, OnlyDoppler(0))
    model.freq_ghz = freq_ghz
    for part, velocity in ((0, orbital), (1, model.solar_velocity_km_s + orbital)):
        model.static_sky = OnlyDoppler(part)
        expected = static_doppler_mjy_sr(np.ones((4, 2)), table["boresight"], velocity)
        assert np.allclose(model.sky_mjy_sr(table), expected, rtol=1e-12, atol=0)


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
