from dataclasses import replace

import numpy as np
import pytest

from fossil_sim.operations import (
    FLAG_MODE_SWITCH,
    FLAG_RECYCLING,
    MODE_CODES,
    bbir_temperatures_k,
    cycle_timeline,
    ifg_modes,
    mode_sky_ports,
    num_cycles,
    operation_flags,
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


SCAN_S, RING_S = 1440.0, 23 * 1440.0


def with_modes(config, interval_s, pattern, dead_time_s=300.0):
    ops = replace(
        config.operations,
        mode_switch_interval_s=interval_s,
        mode_pattern=tuple(pattern),
        mode_switch_dead_time_s=dead_time_s,
    )
    return replace(config, operations=ops)


def switch_ifgs(config, cycle_id):
    """Local index of the first dead IFG of each switch, and the dead count."""
    tl = cycle_timeline(config, cycle_id)
    _, flags = operation_flags(config, tl)
    dead = np.flatnonzero(flags & FLAG_MODE_SWITCH)
    starts = dead[np.diff(dead, prepend=-2) > 1]
    return starts, dead.size


def test_baseline_switches_mid_ring_at_south_pole(config):
    # half Rings (B10, B01, B01, B10): every Ring in both Modes, switches only
    # at spin phase π of the middle Scan (11.5 Scans into each Ring)
    tl = cycle_timeline(config, 0)
    mode, flags = operation_flags(config, tl)
    halves = np.split(mode[~tl.in_recycling], 6)
    assert [np.unique(h).tolist() for h in halves] == [
        [m] for m in (B10, B01, B01, B10, B10, B01)
    ]
    assert np.all(mode[tl.in_recycling] == B00)
    starts, n_dead = switch_ifgs(config, 0)
    assert starts.tolist() == [(23 * r + 11) * 720 + 360 for r in range(3)]
    assert n_dead == 3 * 150
    assert np.all(((flags & FLAG_RECYCLING) > 0) == tl.in_recycling)


@pytest.mark.parametrize(
    "interval_s, pattern, cycle_modes, switches",
    [
        # one cycle in B01, the next in B10: switches only across Recycling
        (3 * RING_S, ["B01", "B10"], [[B01], [B10], [B01]], []),
        (RING_S, ["B10"], [[B10], [B10], [B10]], []),
        (RING_S, ["B10", "B01"], [[B10, B01, B10], [B01, B10, B01]], [23, 46]),
        (
            RING_S / 3,
            ["B10", "B01", "B00"],
            [[B10, B01, B00] * 3],
            [23 / 3 * k for k in range(1, 9)],
        ),
        (1.5 * RING_S, ["B10", "B01"], [[B10, B01], [B10, B01]], [34.5]),
    ],
)
def test_mode_blocks(config, interval_s, pattern, cycle_modes, switches):
    cfg = with_modes(config, interval_s, pattern)
    for cycle_id, expected in enumerate(cycle_modes):
        tl = cycle_timeline(cfg, cycle_id)
        observing = ifg_modes(cfg, tl)[~tl.in_recycling]
        segments = np.array_split(observing, len(expected))
        assert [np.unique(seg).tolist() for seg in segments] == [[m] for m in expected]
        starts, n_dead = switch_ifgs(cfg, cycle_id)
        assert starts.tolist() == [round(s * 720) for s in switches]
        assert n_dead == 150 * len(switches)


def test_scan_alternation_and_dead_time(config):
    # Scan by Scan; dead time rounded up to whole IFGs; 0 = instantaneous
    cfg = with_modes(config, SCAN_S, ["B10", "B01"], dead_time_s=3.0)
    starts, n_dead = switch_ifgs(cfg, 0)
    assert starts.tolist() == [s * 720 for s in range(1, 69)] and n_dead == 2 * 68
    assert switch_ifgs(with_modes(cfg, SCAN_S, ["B10", "B01"], 0.0), 0)[1] == 0
    tl = cycle_timeline(cfg, 1)
    modes = ifg_modes(cfg, tl)[~tl.in_recycling][::720]
    assert modes[0] == B01  # 69 observing Scans per cycle: odd, continuous alternation


def test_block_resumes_after_recycling(config):
    # Two Rings per block: the clock stops in Recycling, so the second block
    # spans Ring 2 of cycle 0 and Ring 0 of cycle 1, with no switch (no dead
    # time) at the start of cycle 1
    cfg = with_modes(config, 2 * RING_S, ["B10", "B01"])
    rings = []
    for cycle_id in (0, 1):
        tl = cycle_timeline(cfg, cycle_id)
        rings += np.split(ifg_modes(cfg, tl)[~tl.in_recycling], 3)
    assert [np.unique(r).tolist() for r in rings] == [
        [m] for m in (B10, B10, B01, B01, B10, B10)
    ]
    assert switch_ifgs(cfg, 0)[0].tolist() == [46 * 720]
    assert switch_ifgs(cfg, 1)[0].tolist() == [23 * 720]


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
