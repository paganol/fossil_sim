from dataclasses import replace

import healpy as hp
import numpy as np
import pytest

from fossil_sim.instrument import detector_units
from fossil_sim.noise import du_noise_mjy_sr, noise_sigma_mjy_sr
from fossil_sim.simulation import run_simulation


def test_sigma_matches_pixie_formula(config):
    lf = detector_units(config.instrument)[0]
    # NEP / (√(2 t) AΩ Δν α f), in W m⁻² Hz⁻¹ sr⁻¹, then MJy/sr
    expected = 4.8e-16 / (np.sqrt(4.0) * 4e-4 * 15e9 * 0.57 * 0.9) * 1e20
    assert noise_sigma_mjy_sr(config, lf) == pytest.approx(expected)


def test_realization_keyed_and_scaled(config):
    du = detector_units(config.instrument)[1]
    recycling = np.zeros(20000, dtype=bool)
    recycling[10000:] = True
    a = du_noise_mjy_sr(config, du, 1, 7, recycling, 5)
    assert np.array_equal(a, du_noise_mjy_sr(config, du, 1, 7, recycling, 5))
    assert not np.array_equal(a, du_noise_mjy_sr(config, du, 1, 8, recycling, 5))
    sigma = noise_sigma_mjy_sr(config, du)
    assert a[:10000].std() == pytest.approx(sigma, rel=0.02)
    assert a[10000:].std() == pytest.approx(3 * sigma, rel=0.02)


def test_noise_map_matches_variance_map(config):
    quiet = replace(
        config,
        static_sky=replace(config.static_sky, models=()),
        zodi=replace(config.zodi, enabled=False),
    )
    noisy = run_simulation(quiet, num_of_cycles=1).spectral_maps
    clean = run_simulation(
        replace(quiet, noise=replace(quiet.noise, enabled=False)), num_of_cycles=1
    ).spectral_maps
    seen = clean.sky_count[64] > 0
    residual = noisy.sky_map(64)[:, seen] - clean.sky_map(64)[:, seen]
    chi = residual / np.sqrt(noisy.sky_variance_map(64)[:, seen])
    assert chi.std() == pytest.approx(1.0, rel=0.02)
    assert np.all(noisy.sky_variance_map(64)[:, ~seen] == hp.UNSEEN)
