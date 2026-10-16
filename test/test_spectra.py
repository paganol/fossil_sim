from importlib.resources import files

import numpy as np
import pytest

from fossil_sim.spectra import (
    MU_BETA,
    blackbody_mjy_sr,
    cib_mjy_sr,
    extragalactic_co_mjy_sr,
    mu_distortion_mjy_sr,
    relativistic_sz_correction_mjy_sr,
    y_distortion_mjy_sr,
)

T_K = 2.72548
CO_TEMPLATE = files("fossil_sim").joinpath("data/extragalactic_co_template.txt")
FREQ_GHZ = np.arange(60.0, 2000.0, 15.0)


def x_to_freq_ghz(x):
    return x * 1.380649e-23 * T_K / 6.62607015e-34 / 1e9


def test_mu_zero_crossing():
    freq_ghz = x_to_freq_ghz(np.array([MU_BETA * 0.99, MU_BETA * 1.01]))
    values = mu_distortion_mjy_sr(freq_ghz, 1e-8, T_K)
    assert values[0] < 0 < values[1]


def test_rsz_is_small_correction_of_y():
    y = y_distortion_mjy_sr(FREQ_GHZ, 1e-6, T_K)
    rsz = relativistic_sz_correction_mjy_sr(FREQ_GHZ, 1e-6, 1.24, T_K)
    peak = np.argmax(np.abs(y))
    assert 0 < np.abs(rsz[peak]) / np.abs(y[peak]) < 0.05
    assert np.all(relativistic_sz_correction_mjy_sr(FREQ_GHZ, 1e-6, 0.0, T_K) == 0.0)


def test_cib_and_co():
    cib = cib_mjy_sr(FREQ_GHZ, 1.3e-5, 0.64, 18.5, 3000.0)
    assert np.all(cib > 0)
    assert (
        1000 < FREQ_GHZ[np.argmax(cib)] < 1600
    )  # modified blackbody peaks at ~1.3 THz
    co = extragalactic_co_mjy_sr(np.array([5.0, 100.0, 615.0, 700.0]), 1.0)
    assert co[0] == 0.0 and co[3] == 0.0 and co[1] > 0 and co[2] > 0
    assert np.all(blackbody_mjy_sr(FREQ_GHZ, T_K) > 0)


def test_against_bisou_sky():
    pytest.importorskip("pysm3")
    gs = pytest.importorskip("bisou_sky.generate_sky")
    assert np.allclose(
        blackbody_mjy_sr(FREQ_GHZ, T_K), gs.planck(FREQ_GHZ, T_K), rtol=1e-6
    )
    assert np.allclose(
        y_distortion_mjy_sr(FREQ_GHZ, 1e-6, T_K),
        gs.deltaI_y_distortions(FREQ_GHZ, 1e-6, T_K),
        rtol=1e-6,
        atol=1e-12,
    )
    assert np.allclose(
        relativistic_sz_correction_mjy_sr(FREQ_GHZ, 1e-6, 1.24, T_K),
        gs.deltaI_y_distortions_relativistic_corrections(FREQ_GHZ, 1e-6, 1.24, T_K),
        rtol=1e-6,
        atol=1e-14,
    )
    assert np.allclose(
        mu_distortion_mjy_sr(FREQ_GHZ, 1e-8, T_K),
        gs.deltaI_mu_distortions(FREQ_GHZ, 1e-8, T_K),
        rtol=1e-6,
        atol=1e-14,
    )
    assert np.allclose(
        cib_mjy_sr(FREQ_GHZ, 1.3e-5, 0.64, 18.5, 3000.0), gs.CIB(FREQ_GHZ), rtol=1e-6
    )
    inside = FREQ_GHZ <= 615
    assert np.allclose(
        extragalactic_co_mjy_sr(FREQ_GHZ, 1.0)[inside],
        # bisou_sky wheels miss their data file: pass the (identical) template
        gs.extragalactic_CO(FREQ_GHZ, 1.0, template=str(CO_TEMPLATE))[inside],
        rtol=1e-6,
    )


test_against_bisou_sky = pytest.mark.heavy(test_against_bisou_sky)
