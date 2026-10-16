from dataclasses import replace

import numpy as np

from fossil_sim.operations import (
    FLAG_MODE_SWITCH,
    FLAG_RECYCLING,
    MODE_CODES,
    bbir_temperatures_k,
    cycle_timeline,
    mode_sky_ports,
    num_cycles,
    operation_flags,
    scan_modes,
)

B10, B01, B00 = MODE_CODES["B10"], MODE_CODES["B01"], MODE_CODES["B00"]


def test_timeline_structure(config):
    tl = cycle_timeline(config, 2)
    assert tl.ifg_index.size == 87 * 720
    assert tl.time_start_s[0] == 2 * 87 * 1440.0
    scan_duration_s = np.ptp(tl.time_start_s[tl.scan_in_cycle == 5]) + 2.0
    assert scan_duration_s == 1440.0
    assert np.all(np.bincount(tl.scan_in_cycle) == 720)
    observing = ~tl.in_recycling
    assert np.all(np.bincount(tl.ring_id[observing] - 6) == 23 * 720)
    assert np.all(tl.ring_id[tl.in_recycling] == -1)
    assert tl.in_recycling.sum() == 18 * 720


def test_last_cycle_is_truncated(config):
    last = num_cycles(config) - 1
    tl = cycle_timeline(config, last)
    assert 0 < tl.ifg_index.size <= 87 * 720
    assert tl.time_start_s[-1] < config.simulation.duration_days * 86400.0


def test_mode_alternation_is_continuous(config):
    observing_scans = np.concatenate([np.arange(69) + 87 * c for c in range(3)])
    modes = scan_modes(config, observing_scans)
    assert np.all(modes[::2] == B10) and np.all(modes[1::2] == B01)
    # Rings alternate their first Mode; cycle 1 starts with B01 (69 is odd)
    assert scan_modes(config, np.array([23]))[0] == B01
    assert scan_modes(config, np.array([87]))[0] == B01
    assert np.all(scan_modes(config, np.arange(69, 87)) == B00)


def test_flags(config):
    tl = cycle_timeline(config, 1)
    mode, flags = operation_flags(config, tl)
    switch = (flags & FLAG_MODE_SWITCH) > 0
    # One dead IFG at the start of every observing Scan, and entering Recycling
    assert switch.sum() == 69 + 1
    assert np.all(tl.ifg_in_scan[switch] == 0)
    assert np.all(((flags & FLAG_RECYCLING) > 0) == tl.in_recycling)
    assert np.all(mode[tl.in_recycling] == B00)


def test_bbir_schedule_changes_only_in_recycling(config):
    cfg = replace(
        config,
        bbir=replace(config.bbir, setpoints_k=(2.6, 2.8), port_offsets_k=(0.0, 1e-3)),
    )
    for cycle_id in (0, 1):
        tl = cycle_timeline(cfg, cycle_id)
        port_k, hf_k = bbir_temperatures_k(cfg, tl)
        assert np.all(port_k[~tl.in_recycling, 0] == cfg.bbir.setpoints_k[cycle_id])
        assert np.all(port_k[tl.in_recycling, 0] == cfg.bbir.setpoints_k[1 - cycle_id])
        assert np.allclose(port_k[:, 1] - port_k[:, 0], 1e-3)
        assert np.all(hf_k == 0.0)


def test_mode_sky_ports():
    codes = np.array([MODE_CODES[m] for m in ("B00", "B01", "B10", "B11")])
    port1, port2 = mode_sky_ports(codes)
    assert port1.tolist() == [False, False, True, True]
    assert port2.tolist() == [False, True, False, True]
