from dataclasses import replace

import healpy as hp
import numpy as np
import pytest

from fossil_sim.ifg import differential_spectrum, interferogram, opd_samples_m
from fossil_sim.instrument import detector_units, frequency_grid_ghz
from fossil_sim.noise import du_noise_mjy_sr, noise_sigma_mjy_sr
from fossil_sim.simulation import run_simulation


def test_sigma_formula(config):
    lf = detector_units(config.instrument)[0]
    # NEP / (√t AΩ Δν α f m), in W m⁻² Hz⁻¹ sr⁻¹, then MJy/sr
    expected = 4.8e-16 / (np.sqrt(2.0) * 4e-4 * 15e9 * 0.57 * 0.9 * 1.0) * 1e20
    assert noise_sigma_mjy_sr(config, lf) == pytest.approx(expected)


@pytest.mark.parametrize("modulation_efficiency", [1.0, 0.25])
def test_sigma_from_interferogram_samples(config, modulation_efficiency):
    # End to end: NEP / √(2 δt) on each IFG sample -> least-squares spectrum;
    # its scatter must equal noise_sigma_mjy_sr (PIXIE's √(2t) gives 1/√2 of it)
    cfg = replace(
        config, noise=replace(config.noise, modulation_efficiency=modulation_efficiency)
    )
    nc, t_s = cfg.noise, cfg.timeline.ifg_duration_s
    freq_ghz = frequency_grid_ghz(cfg.instrument)
    du = detector_units(cfg.instrument)[0]
    opd_m = opd_samples_m(cfg.ifg)
    # W per (MJy/sr × GHz) of interferogram
    w_per_unit = (
        nc.etendue_cm2_sr * 1e-4 * nc.absorptivity * nc.optical_efficiency * 1e-11
    )
    sigma_sample = nc.nep_lf_w_sqrt_hz / np.sqrt(2 * t_s / opd_m.size) / w_per_unit
    clean = interferogram(
        freq_ghz,
        np.full(freq_ghz.size, 300.0),
        np.full(freq_ghz.size, 299.0),
        du,
        opd_m,
        modulation_efficiency,
    )
    rng = np.random.default_rng(5)
    recovered = np.array(
        [
            differential_spectrum(
                freq_ghz,
                clean + rng.normal(0.0, sigma_sample, opd_m.size),
                du,
                opd_m,
                modulation_efficiency,
            )[du.band_mask(freq_ghz)]
            for _ in range(1000)
        ]
    )
    assert recovered.mean() == pytest.approx(1.0, abs=0.05 * recovered.std())
    ratio = recovered.std(axis=0) / noise_sigma_mjy_sr(cfg, du)
    assert np.mean(ratio) == pytest.approx(1.0, abs=0.03)


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
