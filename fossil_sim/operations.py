"""Observation timeline: IFG/Scan/Ring/Cooler-cycle indices, Modes, flags, BBIR.

Every quantity is a pure function of global indices (ADR 0001), so any
Cooler cycle can be built in isolation.

Timeline of one Cooler cycle (``n = rings_per_cycle`` Rings, then the
Recycling window)::

    | Ring 0 | Ring 1 | ... | Ring n-1 | Recycling window |
      each Ring = scans_per_ring Scans, each Scan = ifgs_per_scan IFGs

Mode convention (see CONTEXT.md, "Mode"): ``Bxy`` with x = BSM1, y = BSM2,
1 = sky, 0 = BBIR.
"""

from dataclasses import dataclass

import numpy as np

from .config import MODES, Config

# Flag bits (uint8)
FLAG_RECYCLING = np.uint8(1)
FLAG_MODE_SWITCH = np.uint8(2)
FLAG_SAA_EXCEEDED = np.uint8(4)
FLAG_SLEW = np.uint8(8)  # reserved: slews are instantaneous in v0.1
FLAG_NAMES = {
    "recycling": FLAG_RECYCLING,
    "mode_switch": FLAG_MODE_SWITCH,
    "saa_exceeded": FLAG_SAA_EXCEEDED,
    "slew": FLAG_SLEW,
}

MODE_CODES = {name: np.int8(code) for code, name in enumerate(MODES)}


def mode_sky_ports(mode_code: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Whether Port 1 and Port 2 look at the sky, for integer Mode codes.

    >>> mode_sky_ports(np.array([MODE_CODES["B10"], MODE_CODES["B01"]]))
    (array([ True, False]), array([False,  True]))
    """
    port1_on_sky = np.array([name[1] == "1" for name in MODES])
    port2_on_sky = np.array([name[2] == "1" for name in MODES])
    return port1_on_sky[mode_code], port2_on_sky[mode_code]


def num_cycles(config: Config) -> int:
    """Number of (possibly truncated) Cooler cycles in the mission."""
    duration_s = config.simulation.duration_days * 86400.0
    return int(np.ceil(duration_s / config.timeline.cycle_duration_s))


@dataclass
class CycleTimeline:
    """Per-IFG indices and times of one Cooler cycle (all arrays of length n_ifg)."""

    cycle_id: int
    ifg_index: np.ndarray  # global IFG index
    time_start_s: np.ndarray
    time_mid_s: np.ndarray
    ifg_in_scan: np.ndarray
    scan_in_cycle: np.ndarray
    ring_in_cycle: np.ndarray  # rings_per_cycle - 1 during the Recycling window
    in_recycling: np.ndarray
    scan_id: np.ndarray  # global Scan index (Recycling Scans included)
    ring_id: np.ndarray  # global Ring index; -1 during the Recycling window


def cycle_timeline(config: Config, cycle_id: int) -> CycleTimeline:
    """Build the timeline of one Cooler cycle, truncated at the mission end."""
    tl = config.timeline
    first_ifg = cycle_id * tl.ifgs_per_cycle
    ifg_index = first_ifg + np.arange(tl.ifgs_per_cycle, dtype=np.int64)
    time_start_s = ifg_index * tl.ifg_duration_s
    keep = time_start_s < config.simulation.duration_days * 86400.0
    ifg_index, time_start_s = ifg_index[keep], time_start_s[keep]

    local_ifg = ifg_index - first_ifg
    scan_in_cycle = local_ifg // tl.ifgs_per_scan
    ring_raw = scan_in_cycle // tl.scans_per_ring
    in_recycling = ring_raw >= tl.rings_per_cycle
    ring_in_cycle = np.minimum(ring_raw, tl.rings_per_cycle - 1)
    return CycleTimeline(
        cycle_id=cycle_id,
        ifg_index=ifg_index,
        time_start_s=time_start_s,
        time_mid_s=time_start_s + 0.5 * tl.ifg_duration_s,
        ifg_in_scan=local_ifg % tl.ifgs_per_scan,
        scan_in_cycle=scan_in_cycle,
        ring_in_cycle=ring_in_cycle,
        in_recycling=in_recycling,
        scan_id=cycle_id * tl.scans_per_cycle + scan_in_cycle,
        ring_id=np.where(
            in_recycling, -1, cycle_id * tl.rings_per_cycle + ring_in_cycle
        ),
    )


def mode_switch_dead_ifgs(config: Config) -> int:
    """IFGs discarded after a Mode switch: the dead time rounded up to whole IFGs."""
    tl, ops = config.timeline, config.operations
    return int(np.ceil(ops.mode_switch_dead_time_s / tl.ifg_duration_s - 1e-9))


def _mode_blocks(config: Config, timeline: CycleTimeline):
    """Per IFG: global Mode-block index, and IFG offset from the block start.

    Blocks of ``mode_switch_interval_s`` tile the observing time only: the
    clock stops during Recycling windows, so a block may resume after one.
    """
    tl = config.timeline
    block_ifgs = round(config.operations.mode_switch_interval_s / tl.ifg_duration_s)
    local_ifg = timeline.ifg_index - timeline.cycle_id * tl.ifgs_per_cycle
    observing_ifg = timeline.cycle_id * tl.observing_ifgs_per_cycle + local_ifg
    block, offset = np.divmod(observing_ifg, block_ifgs)
    return block, offset


def _mode_pattern(config: Config) -> np.ndarray:
    return np.array([MODE_CODES[m] for m in config.operations.mode_pattern], np.int8)


def ifg_modes(config: Config, timeline: CycleTimeline) -> np.ndarray:
    """Mode code of each IFG.

    The observing time is split into Mode blocks of ``mode_switch_interval_s``
    and ``mode_pattern`` runs over the block index. Recycling IFGs get
    ``recycling_mode``.
    """
    block, _ = _mode_blocks(config, timeline)
    pattern = _mode_pattern(config)
    return np.where(
        timeline.in_recycling,
        MODE_CODES[config.operations.recycling_mode],
        pattern[block % pattern.size],
    ).astype(np.int8)


def bbir_temperatures_k(
    config: Config, timeline: CycleTimeline
) -> tuple[np.ndarray, np.ndarray]:
    """BBIR main-absorber temperature per Port and HF-cavity temperature.

    The set-point of cycle ``c`` is ``setpoints_k[c % len]``; it changes only
    during Recycling windows, so the window of cycle ``c`` already carries
    the set-point of cycle ``c + 1`` (the transient is not modelled).

    Returns
    -------
    port_temperature_k : ndarray, shape (n_ifg, 2)
    hf_cavity_temperature_k : ndarray, shape (n_ifg,)
    """
    bbir = config.bbir
    setpoints_k = np.asarray(bbir.setpoints_k)
    cycle = timeline.cycle_id + timeline.in_recycling.astype(np.int64)
    main_k = setpoints_k[cycle % len(setpoints_k)]
    port_temperature_k = main_k[:, None] + np.asarray(bbir.port_offsets_k)[None, :]
    return port_temperature_k, np.full(main_k.shape, bbir.hf_cavity_k)


def operation_flags(
    config: Config, timeline: CycleTimeline
) -> tuple[np.ndarray, np.ndarray]:
    """Mode code and flags (without the SAA bit) of each IFG.

    After a Mode switch within an observing period the first
    :func:`mode_switch_dead_ifgs` IFGs are flagged. The start of an observing
    period (mission start, end of Recycling) needs no switch: the BSMs are
    set during the Recycling window.
    """
    mode = ifg_modes(config, timeline)
    block, offset = _mode_blocks(config, timeline)
    pattern = _mode_pattern(config)
    local_ifg = timeline.ifg_index - timeline.cycle_id * config.timeline.ifgs_per_cycle
    # the block started inside this observing period, not at its start
    started_here = local_ifg - offset > 0
    switched = started_here & (
        pattern[block % pattern.size] != pattern[(block - 1) % pattern.size]
    )
    dead = switched & ~timeline.in_recycling & (offset < mode_switch_dead_ifgs(config))
    flags = np.zeros(mode.shape, dtype=np.uint8)
    flags[timeline.in_recycling] |= FLAG_RECYCLING
    flags[dead] |= FLAG_MODE_SWITCH
    return mode, flags
