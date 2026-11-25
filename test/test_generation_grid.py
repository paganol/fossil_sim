"""IFG generation on a frequency grid finer and wider than the instrument bins."""

from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.config import config_from_dict
from fossil_sim.instrument import (
    detector_units,
    frequency_grid_ghz,
    generation_frequency_grid_ghz,
    sky_frequency_grid_ghz,
)
from fossil_sim.response import generation_response, true_response
from fossil_sim.simulation import SignalModel, build_ifg_table
from fossil_sim.spectra import blackbody_mjy_sr

STEP_GHZ, MAX_GHZ = 3.0, 3825.0


def fine_config(config, edge_fraction=0.0, **ifg):
    return replace(
        config,
        static_sky=replace(config.static_sky, models=()),
        zodi=replace(config.zodi, enabled=False),
        noise=replace(config.noise, enabled=False),
        instrument=replace(config.instrument, band_edge_fraction=edge_fraction),
        ifg=replace(
            config.ifg,
            simulate=True,
            sky_variation=False,
            generation_freq_step_ghz=STEP_GHZ,
            generation_freq_max_ghz=MAX_GHZ,
            **ifg,
        ),
    )


@pytest.fixture(scope="module")
def table(config, orbit):
    table = build_ifg_table(config, orbit, 0)
    return np.concatenate([table[90:110], table[8580:8600]])  # B10 and B01


def test_top_hat_transmission_is_the_band_mask(config):
    freq_ghz = np.arange(1.0, 2500.0, 0.5)
    for du in detector_units(config.instrument):
        np.testing.assert_array_equal(
            du.transmission(freq_ghz, 0.0), du.band_mask(freq_ghz)
        )


def test_cosine_edges(config):
    inst, fraction = config.instrument, 0.1
    lf, hf = detector_units(inst)[:2]
    freq_ghz = np.arange(1.0, 2500.0, 0.25)
    t_lf, t_hf = (
        lf.transmission(freq_ghz, fraction),
        hf.transmission(freq_ghz, fraction),
    )
    assert np.all((t_lf >= 0) & (t_lf <= 1) & (t_hf >= 0) & (t_hf <= 1))
    # the dichroic splits the power: the DUs of an output sum to 1 between the edges
    inner = (freq_ghz > inst.band_min_ghz * 1.05) & (
        freq_ghz < inst.band_max_ghz * 0.95
    )
    np.testing.assert_allclose(t_lf[inner] + t_hf[inner], 1.0, atol=1e-12)
    edges = (inst.band_min_ghz, inst.dichroic_split_ghz, inst.band_max_ghz)
    for du in (lf, hf):
        for edge in (du.band_min_ghz, du.band_max_ghz):
            t = du.transmission(np.array([edge]), fraction)[0]
            assert t == pytest.approx(0.5)
        outside = np.ones(freq_ghz.size, dtype=bool)
        for edge in edges:
            outside &= np.abs(freq_ghz - edge) >= fraction * edge / 2
        far = np.where(du.band_mask(freq_ghz), 1.0, 0.0)
        np.testing.assert_allclose(
            du.transmission(freq_ghz, fraction)[outside], far[outside], atol=1e-12
        )


def test_generation_grid(config):
    # default: the instrument grid, in both domains
    np.testing.assert_array_equal(
        generation_frequency_grid_ghz(config), frequency_grid_ghz(config.instrument)
    )
    cfg = fine_config(config, edge_fraction=0.1)
    grid = generation_frequency_grid_ghz(cfg)
    np.testing.assert_allclose(np.diff(grid), STEP_GHZ)
    assert grid[-1] <= MAX_GHZ
    assert np.all(np.isin(frequency_grid_ghz(cfg.instrument), grid))
    # only frequencies some DU transmits: the band, widened by half an edge
    assert grid[0] > cfg.instrument.band_min_ghz * 0.95
    assert cfg.instrument.band_max_ghz < grid[-1] < cfg.instrument.band_max_ghz * 1.05
    np.testing.assert_array_equal(sky_frequency_grid_ghz(cfg), grid)
    # IFG-domain options in the spectral domain: ignored, with a warning
    with pytest.warns(UserWarning, match="only in the IFG domain"):
        spectral = replace(cfg, ifg=replace(cfg.ifg, simulate=False))
    np.testing.assert_array_equal(
        sky_frequency_grid_ghz(spectral), frequency_grid_ghz(cfg.instrument)
    )


@pytest.mark.parametrize(
    "ifg, instrument",
    [
        ({"generation_freq_step_ghz": 3.0}, {}),  # max missing
        ({"generation_freq_step_ghz": 4.0, "generation_freq_max_ghz": 3000.0}, {}),
        ({"generation_freq_step_ghz": 3.0, "generation_freq_max_ghz": 1500.0}, {}),
        ({"generation_freq_step_ghz": 3.0, "generation_freq_max_ghz": 3900.0}, {}),
        ({}, {"band_edge_fraction": -0.1}),
        ({}, {"band_edge_fraction": 2.0}),  # overlapping edges
    ],
)
def test_invalid_generation_settings(config, ifg, instrument):
    values = config.to_dict()
    values["ifg"].update(simulate=True, **ifg)
    values["instrument"].update(instrument)
    with pytest.raises(ValueError):
        config_from_dict(values)


def test_generation_response(config):
    cfg = replace(
        config,
        response=replace(
            config.response,
            gain=(1.02, 0.97, 1.05, 0.99),
            gain_rms=0.05,
            offset_emissivity=(1e-3, -2e-3, -1e-3, 5e-4),
        ),
    )
    cfg = fine_config(cfg, edge_fraction=0.1)
    freq_ghz = frequency_grid_ghz(cfg.instrument)
    fine_ghz = generation_frequency_grid_ghz(cfg)
    dus = detector_units(cfg.instrument)
    response = true_response(cfg, freq_ghz, dus)
    fine = generation_response(cfg, freq_ghz, response, fine_ghz, dus)
    rc = cfg.response
    for i, du in enumerate(dus):
        transmits = fine_ghz[du.transmission(fine_ghz, 0.1) > 0]
        bins = freq_ghz[du.band_mask(freq_ghz)]
        at_bins = np.isin(transmits, bins)
        # the gain at a bin centre is the bin gain; constant within a bin
        np.testing.assert_array_equal(
            fine.gain[du.name][at_bins], response.gain[du.name]
        )
        np.testing.assert_allclose(
            fine.offset_mjy_sr[du.name],
            rc.offset_emissivity[i]
            * blackbody_mjy_sr(transmits, rc.instrument_temperature_k),
            rtol=1e-12,
        )


def test_fine_grid_matches_spectral_domain_for_smooth_skies(config, orbit, table):
    """Smooth sky, top-hat bands: the fine-grid IFG recovers the bin values.

    The difference is the average of the spectrum over each 15 GHz bin; the
    edge bin of the band, cut by the top-hat at band_max_ghz, is excluded.
    """
    spectral = SignalModel(replace(fine_config(config), ifg=config.ifg), orbit)
    fine = SignalModel(fine_config(config), orbit)
    expected = spectral.compute(table)[1]
    estimate = fine.compute(table)[1]
    freq_ghz = fine.freq_ghz
    inner = (freq_ghz < 1500.0) & (freq_ghz > 30.0)
    scale = np.abs(expected).max(axis=0)
    error = np.abs(estimate - expected).max(axis=0) / scale
    assert np.all(error[inner] < 3e-3)
    # the effect is real: the fine grid is not the instrument grid
    assert error[inner].max() > 1e-6


def test_band_edges_act_most_on_the_edge_bins(config, orbit, table):
    """Cosine edges change the edge bins; the rest only through fit leakage.

    The least squares (no apodization) spreads the change of the spectrum at
    the edges over the whole band, decreasing away from them.
    """
    top_hat = SignalModel(fine_config(config), orbit).compute(table)[1]
    tapered = SignalModel(fine_config(config, edge_fraction=0.1), orbit).compute(table)[
        1
    ]
    freq_ghz = frequency_grid_ghz(config.instrument)
    inst = config.instrument
    near_edge = np.zeros(freq_ghz.size, dtype=bool)
    for edge in (inst.band_min_ghz, inst.dichroic_split_ghz, inst.band_max_ghz):
        near_edge |= np.abs(freq_ghz - edge) < 0.05 * edge + inst.freq_step_ghz
    change = np.abs(tapered - top_hat).max(axis=0) / np.abs(top_hat).max(axis=0)
    assert change[near_edge].max() > 1e-2
    assert change[~near_edge].max() < 0.2 * change[near_edge].max()
    # mid HF band, far from every edge
    mid = (freq_ghz > 800.0) & (freq_ghz < 1500.0)
    assert np.all(change[mid] < 1e-3)
