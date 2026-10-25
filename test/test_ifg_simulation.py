from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.ifg import sample_weights
from fossil_sim.noise import noise_sigma_mjy_sr
from fossil_sim.scanning import ring_circles, ring_spin_axes
from fossil_sim.simulation import SignalModel, build_ifg_table


class ProfileSky:
    """Static-sky stub: spectrum A(ν) times a profile g(n̂), no Doppler."""

    def __init__(self, freq_ghz, profile):
        self.amplitude = 1e-2 * (freq_ghz / 100.0) ** 1.5
        self.profile = profile

    def sample(self, direction):
        return self.profile(direction)[:, None] * self.amplitude

    def sample_doppler(self, direction):
        zero = np.zeros((len(direction), self.amplitude.size))
        return zero, zero


def blob(center, fwhm_deg=1.6):
    """Gaussian blob of the beam size."""
    sigma = np.deg2rad(fwhm_deg) / np.sqrt(8 * np.log(2))
    return lambda n: np.exp(-np.sum((n - center) ** 2, axis=1) / (2 * sigma**2))


def ifg_config(config, sky_variation=True, noise=False):
    return replace(
        config,
        ifg=replace(config.ifg, simulate=True, sky_variation=sky_variation),
        zodi=replace(config.zodi, enabled=False),
        noise=replace(config.noise, enabled=noise),
    )


@pytest.fixture(scope="module")
def table(config, orbit):
    table = build_ifg_table(config, orbit, 0)
    # science IFGs in both Modes and Recycling (B00) IFGs
    return np.concatenate([table[90:110], table[8580:8600], table[-20:]])


def test_sample_weights_are_exact_for_cubics():
    q, n = 8, 512
    weights = sample_weights(n, q)
    g = np.arange(-1, q + 2)  # window grid, IFG units × q
    u = (np.arange(n) + 0.5) * q / n
    for degree in range(4):
        assert np.allclose(weights @ g**degree, u**degree, rtol=1e-12, atol=1e-9)


def test_ifg_domain_without_sky_variation_matches_spectral(config, orbit, table):
    center = table["boresight"][10]
    spectral = SignalModel(replace(ifg_config(config), ifg=config.ifg), orbit, None)
    spectral.static_sky = ProfileSky(spectral.freq_ghz, blob(center))
    model = SignalModel(ifg_config(config, sky_variation=False), orbit, None)
    model.static_sky = spectral.static_sky
    expected, estimate_expected, _ = spectral.compute(table)
    du_spectra, estimate, _ = model.compute(table)
    for name, s in expected.items():
        assert np.allclose(du_spectra[name], s, rtol=0, atol=1e-9 * np.abs(s).max())
    science = ~np.isnan(estimate_expected[:, 0])
    assert np.allclose(estimate[science], estimate_expected[science], rtol=1e-12)


def test_sky_on_the_ring_is_interpolated_to_the_samples(config, orbit, table):
    cfg = ifg_config(config)
    model = SignalModel(cfg, orbit, None)
    model.static_sky = ProfileSky(model.freq_ghz, blob(table["boresight"][10]))
    circles = model.static_circles(0)
    static = (
        model.sky_windows(table, None, circles)
        - model.sky_mjy_sr(table, include_static=False)[:, None]
    )
    interpolated = np.einsum("jw,iwf->ijf", model.sample_weights, static)
    # exact sky at the sample directions: phase (k + (j + ½) / N_s) / N_I
    n_samples = cfg.ifg.samples_per_ifg
    directions = ring_circles(cfg, ring_spin_axes(cfg, orbit, 0), n_samples)
    k = table["ifg_in_scan"].astype(int)[:, None] * n_samples + np.arange(n_samples)
    exact = model.static_sky.sample(
        directions[table["ring_in_cycle"].astype(int)[:, None], k].reshape(-1, 3)
    ).reshape(interpolated.shape)
    peak = np.abs(exact).max()
    assert np.abs(interpolated - exact).max() < 1e-4 * peak
    # the blob really changes during the IFG it is centred on
    assert np.ptp(exact[10, :, 0]) > 0.03 * np.abs(exact[10, :, 0]).max()


def test_linear_sky_change_during_the_ifg_cancels(config, orbit, table):
    """A symmetric scan rejects the part of the sky odd in time."""
    i = 10
    tangent = np.cross(
        ring_spin_axes(config, orbit, 0)[table["ring_in_cycle"][i]],
        table["boresight"][i],
    )
    sky = lambda freq: ProfileSky(freq, lambda n: n @ tangent)
    model = SignalModel(ifg_config(config), orbit, None)
    model.static_sky = sky(model.freq_ghz)
    fixed = SignalModel(ifg_config(config, sky_variation=False), orbit, None)
    fixed.static_sky = model.static_sky
    varying = model.compute(table[i : i + 1], circles=model.static_circles(0))[0]
    centre = fixed.compute(table[i : i + 1])[0]
    change = np.deg2rad(0.5) * np.abs(model.static_sky.amplitude).max()
    for name, s in varying.items():
        assert np.abs(s - centre[name]).max() < 1e-3 * change


def test_sample_noise_gives_the_spectral_noise_level(config, orbit, table):
    cfg = ifg_config(config, sky_variation=False, noise=True)
    model = SignalModel(cfg, orbit, None)
    rows = np.resize(table[:40], 400)  # science IFGs (B10, B01)
    noiseless = SignalModel(ifg_config(config, sky_variation=False), orbit, None)
    clean = noiseless.compute(rows)[0]
    noisy = model.compute(rows, noise=model.noise_for_cycle(0, rows))[0]
    for du in model.dus:
        ratio = (noisy[du.name] - clean[du.name]) / noise_sigma_mjy_sr(cfg, du)
        assert np.std(ratio) == pytest.approx(1.0, abs=0.03)


def test_ifg_matches_sample_by_sample_evaluation(config, orbit, table):
    """Interpolated sky windows give the IFG of the exact sky at each sample."""
    from fossil_sim.ifg import differential_spectra, interferogram
    from fossil_sim.signal import bbir_spectra_mjy_sr, port_spectra_mjy_sr

    cfg = ifg_config(config)
    model = SignalModel(cfg, orbit, None)
    model.static_sky = ProfileSky(model.freq_ghz, blob(table["boresight"][10]))
    fixed = SignalModel(ifg_config(config, sky_variation=False), orbit, None)
    fixed.static_sky = model.static_sky
    row = table[10:11]
    n_samples = cfg.ifg.samples_per_ifg
    start = int(row["ifg_in_scan"][0]) * n_samples
    directions = ring_circles(cfg, ring_spin_axes(cfg, orbit, 0), n_samples)[
        int(row["ring_in_cycle"][0]), start : start + n_samples
    ]
    sky = model.sky_mjy_sr(row, include_static=False) + model.static_sky.sample(
        directions
    )
    bbir = bbir_spectra_mjy_sr(model.freq_ghz, row["bbir_port_temperature_k"])
    port1, port2 = port_spectra_mjy_sr(
        sky, np.repeat(bbir, n_samples, axis=0), np.repeat(row["mode"], n_samples)
    )
    varying = model.compute(row, circles=model.static_circles(0))[0]
    centre = fixed.compute(row)[0]
    for du in model.dus:
        exact_ifg = np.array(
            [
                interferogram(
                    model.freq_ghz, port1[j], port2[j], du, model.opd_m[j : j + 1]
                )[0]
                for j in range(n_samples)
            ]
        )
        exact = differential_spectra(model.freq_ghz, exact_ifg[None], du, model.opd_m)
        effect = np.abs(exact - centre[du.name]).max()
        assert effect > 1e-4  # the sky change during the IFG is visible...
        # ...and reproduced by the interpolated windows (cubic interpolation,
        # 8 grid points per IFG: ~1e-4 of the blob on the sky)
        assert np.abs(varying[du.name] - exact).max() < 3e-3 * effect
