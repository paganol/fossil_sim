from dataclasses import replace

import healpy as hp
import numpy as np
import pytest

from fossil_sim.ifg import (
    differential_spectra,
    differential_spectrum,
    interferogram,
    opd_samples_m,
)
from fossil_sim.instrument import detector_units, frequency_grid_ghz
from fossil_sim.noise import (
    du_noise_mjy_sr,
    noise_sigma_mjy_sr,
    sky_estimate_noise_variance,
)
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


@pytest.mark.parametrize("samples_per_ifg", [268, 512])
@pytest.mark.parametrize("modulation_efficiency", [0.25, 1.0])
def test_ifg_variance_propagates_sample_noise(
    config, samples_per_ifg, modulation_efficiency
):
    cfg = replace(
        config,
        ifg=replace(config.ifg, simulate=True, samples_per_ifg=samples_per_ifg),
        noise=replace(config.noise, modulation_efficiency=modulation_efficiency),
    )
    freq_ghz = frequency_grid_ghz(cfg.instrument)
    dus = detector_units(cfg.instrument)
    opd_m = opd_samples_m(cfg.ifg)
    variance = sky_estimate_noise_variance(cfg, freq_ghz, dus)
    # Responses to unit impulses give the exact covariance of the real fit.
    # Compute sample noise directly from NEP, independently of the map formula.
    nc = cfg.noise
    w_per_unit = (
        nc.etendue_cm2_sr * 1e-4 * nc.absorptivity * nc.optical_efficiency * 1e-11
    )
    for du in dus[:2]:  # LF-L and HF-L: the R DUs have the same noise
        nep = nc.nep_lf_w_sqrt_hz if du.name.startswith("LF") else nc.nep_hf_w_sqrt_hz
        sigma_sample = (
            nep
            / np.sqrt(2 * cfg.timeline.ifg_duration_s / samples_per_ifg)
            / w_per_unit
        )
        responses = differential_spectra(
            freq_ghz, np.eye(samples_per_ifg), du, opd_m, modulation_efficiency
        )
        expected = sigma_sample**2 * np.sum(responses**2, axis=0) / 2
        np.testing.assert_allclose(
            variance[du.band_mask(freq_ghz)], expected, rtol=1e-12
        )


def coulon_reference_config(config):
    """FOSSIL reference of Coulon, Aghanim & Maffei 2026 (arXiv:2609.17829), Table 2.

    η = 0.5, AΩ = 3 cm² sr (8 modes at 50 GHz), one optical element (1 %, 4.5 K),
    Martin-Puplett m = 1/2 (Eq. 11–12), FOSSIL bands 50–300–2000 GHz; sky and
    reference at the CMB temperature.
    """
    return replace(
        config,
        instrument=replace(
            config.instrument,
            freq_min_ghz=50.0,
            band_min_ghz=50.0,
            dichroic_split_ghz=300.0,
            band_edge_fraction=0.0,
        ),
        sky=replace(
            config.sky,
            t_cmb_k=2.725,
            y=0.0,
            mu=0.0,
            cib_amplitude=0.0,
            co_amplitude=0.0,
        ),
        bbir=replace(config.bbir, setpoints_k=(2.725,)),
        noise=replace(
            config.noise,
            nep_model="photon",
            etendue_cm2_sr=3.0,
            absorptivity=0.5,
            optical_efficiency=1.0,
            modulation_efficiency=0.5,
            optics_emissivity=(0.01,),
            optics_temperature_k=(4.5,),
        ),
    )


def test_photon_nep_matches_direct_integral(config):
    """Eq. 15–16 of Coulon et al. 2026 evaluated independently of the module."""
    from fossil_sim.noise import photon_nep_w_sqrt_hz

    h, k, c = 6.62607015e-34, 1.380649e-23, 299792458.0
    cfg = coulon_reference_config(config)
    for du in detector_units(cfg.instrument)[:2]:
        nu = np.linspace(du.band_min_ghz, du.band_max_ghz, 200001) * 1e9

        def planck(temperature, nu=nu):
            return 2 * h * nu**3 / c**2 / np.expm1(h * nu / (k * temperature))

        p = 0.5 * 0.5 * 3e-4 * (2 * planck(2.725) + 0.01 * planck(4.5))
        nep2 = np.trapezoid(2 * h * nu * p + 2 * c**2 / (nu**2 * 3e-4) * p**2, nu)
        assert photon_nep_w_sqrt_hz(cfg, du) == pytest.approx(np.sqrt(nep2), rel=2e-3)


def test_photon_noise_gives_coulon_sky_sensitivity(config):
    """Mission sensitivity per channel, Eq. 17–18: √2 NEP / (√t AΩ Δν η T_tot).

    Our per-IFG noise of one DU (m = 1/2), averaged over the L and R outputs
    and over t_obs / t_IFG IFGs, is the same quantity: 1.14 Jy/sr (LF) and
    0.91 Jy/sr (HF) for t_obs = 2.8 yr.
    """
    cfg = coulon_reference_config(config)
    t_obs_s = 2.8 * 365.25 * 86400
    dus = detector_units(cfg.instrument)
    for du, expected_jy_sr in zip(dus[:2], (1.14, 0.91), strict=True):
        per_ifg = noise_sigma_mjy_sr(cfg, du) / np.sqrt(2)  # mean of L and R
        mission = per_ifg * np.sqrt(cfg.timeline.ifg_duration_s / t_obs_s) * 1e6
        assert mission == pytest.approx(expected_jy_sr, rel=0.01)


def test_photon_nep_grows_with_optics_temperature(config):
    from fossil_sim.noise import photon_nep_w_sqrt_hz

    cfg = coulon_reference_config(config)
    hf = detector_units(cfg.instrument)[1]
    neps = [
        photon_nep_w_sqrt_hz(
            replace(cfg, noise=replace(cfg.noise, optics_temperature_k=(t,))), hf
        )
        for t in (3.0, 4.5, 10.0)
    ]
    assert neps[0] < neps[1] < neps[2]


def test_noise_model_validated(config):
    from fossil_sim.config import config_from_dict

    for changes in (
        {"nep_model": "other"},
        {"optics_emissivity": [0.01], "optics_temperature_k": []},
        {"optics_emissivity": [1.5], "optics_temperature_k": [4.5]},
        {"optics_emissivity": [0.01], "optics_temperature_k": [0.0]},
    ):
        values = config.to_dict()
        values["noise"].update(changes)
        with pytest.raises(ValueError):
            config_from_dict(values)


def test_step2_config_uses_photon_noise():
    from pathlib import Path

    from fossil_sim.config import load_config
    from fossil_sim.noise import du_nep_w_sqrt_hz

    cfg = load_config(Path(__file__).parents[1] / "configs" / "step2.toml")
    assert cfg.noise.nep_model == "photon"
    lf, hf = detector_units(cfg.instrument)[:2]
    # CMB-dominated loading, plus CIB and CO in the HF band
    assert 1.6e-16 < du_nep_w_sqrt_hz(cfg, lf) < 1.8e-16
    assert 1.3e-16 < du_nep_w_sqrt_hz(cfg, hf) < 1.5e-16
