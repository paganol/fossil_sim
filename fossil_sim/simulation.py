"""Streaming driver: one Cooler cycle at a time, maps accumulated on the fly.

Cycles are independent (ADR 0001): with an MPI communicator, rank ``r``
processes cycles ``r, r + size, ...`` and maps and statistics are summed at
the end, giving results identical to a serial run.
"""

from dataclasses import dataclass, field

import numpy as np

from .config import MODES, Config
from .mapmaking import HitMaps
from .operations import (
    FLAG_NAMES,
    FLAG_SAA_EXCEEDED,
    bbir_temperatures_k,
    cycle_timeline,
    num_cycles,
    operation_flags,
)
from .orbit import L2Orbit
from .scanning import cycle_pointing

IFG_DTYPE = np.dtype(
    [
        ("time_start_s", np.float64),
        ("time_mid_s", np.float64),
        ("cycle_id", np.int32),
        ("ring_id", np.int32),  # -1 during the Recycling window
        ("scan_id", np.int32),
        ("boresight", np.float64, 3),  # unit vector, J2000 ecliptic
        ("mode", np.int8),  # index into config.MODES
        ("flags", np.uint8),
        ("saa_rad", np.float32),
        ("bbir_port_temperature_k", np.float64, 2),
        ("bbir_hf_temperature_k", np.float64),
    ]
)


def build_ifg_table(config: Config, orbit: L2Orbit, cycle_id: int) -> np.ndarray:
    """Per-IFG record of one Cooler cycle (structured array, ``IFG_DTYPE``).

    The pointing is shared by all detector units; per-DU data (spectra,
    interferograms) are stored separately from v0.2 on.
    """
    timeline = cycle_timeline(config, cycle_id)
    pointing = cycle_pointing(config, orbit, timeline)
    mode, flags = operation_flags(config, timeline)
    flags[pointing.saa_rad > np.deg2rad(config.scanning.max_saa_deg)] |= (
        FLAG_SAA_EXCEEDED
    )
    port_k, hf_k = bbir_temperatures_k(config, timeline)

    table = np.empty(timeline.ifg_index.size, dtype=IFG_DTYPE)
    table["time_start_s"] = timeline.time_start_s
    table["time_mid_s"] = timeline.time_mid_s
    table["cycle_id"] = cycle_id
    table["ring_id"] = timeline.ring_id
    table["scan_id"] = timeline.scan_id
    table["boresight"] = pointing.boresight
    table["mode"] = mode
    table["flags"] = flags
    table["saa_rad"] = pointing.saa_rad
    table["bbir_port_temperature_k"] = port_k
    table["bbir_hf_temperature_k"] = hf_k
    if config.operations.recycling_data == "off":
        table = table[~timeline.in_recycling]
    return table


@dataclass
class SimulationResult:
    """Accumulated products of a run."""

    hit_maps: HitMaps
    num_ifgs: int = 0
    num_valid_ifgs: int = 0
    flag_counts: dict[str, int] = field(default_factory=dict)  # IFGs with each bit set
    mode_counts: dict[str, int] = field(default_factory=dict)  # valid IFGs per Mode
    max_saa_deg: float = 0.0

    def lost_fraction(self) -> dict[str, float]:
        """Fraction of IFGs carrying each flag, plus the total of flagged IFGs."""
        fractions = {
            name: count / self.num_ifgs for name, count in self.flag_counts.items()
        }
        fractions["any"] = 1.0 - self.num_valid_ifgs / self.num_ifgs
        return fractions


def run_simulation(
    config: Config, num_of_cycles: int | None = None, comm=None
) -> SimulationResult:
    """Run the mission (or its first ``num_of_cycles`` cycles) and accumulate hit maps.

    Parameters
    ----------
    config : Config
    num_of_cycles : int, optional
        Number of Cooler cycles to simulate; default: the whole mission.
    comm : mpi4py communicator, optional
        If given, cycles are distributed round-robin over ranks and the
        results are reduced on every rank.
    """
    total_cycles = num_cycles(config) if num_of_cycles is None else num_of_cycles
    rank, size = (0, 1) if comm is None else (comm.rank, comm.size)
    orbit = L2Orbit(
        config.orbit,
        config.simulation.mission_start_time,
        total_cycles * config.timeline.cycle_duration_s,
    )
    result = SimulationResult(hit_maps=HitMaps(config.maps.nsides))
    flag_counts = np.zeros(len(FLAG_NAMES), dtype=np.int64)
    mode_counts = np.zeros(len(MODES), dtype=np.int64)
    num_ifgs = num_valid = 0
    max_saa_rad = 0.0

    for cycle_id in range(rank, total_cycles, size):
        table = build_ifg_table(config, orbit, cycle_id)
        valid = table[table["flags"] == 0]
        result.hit_maps.add(valid["boresight"], valid["mode"], valid["time_mid_s"])
        num_ifgs += table.size
        num_valid += valid.size
        flag_counts += [
            np.count_nonzero(table["flags"] & bit) for bit in FLAG_NAMES.values()
        ]
        mode_counts += np.bincount(valid["mode"], minlength=len(MODES))
        max_saa_rad = max(max_saa_rad, float(table["saa_rad"].max()))

    if comm is not None:
        from mpi4py import MPI

        for array in (*result.hit_maps.hits.values(), flag_counts, mode_counts):
            comm.Allreduce(MPI.IN_PLACE, array, op=MPI.SUM)
        for array in result.hit_maps.first_hit_time_s.values():
            comm.Allreduce(MPI.IN_PLACE, array, op=MPI.MIN)
        num_ifgs = comm.allreduce(num_ifgs, op=MPI.SUM)
        num_valid = comm.allreduce(num_valid, op=MPI.SUM)
        max_saa_rad = comm.allreduce(max_saa_rad, op=MPI.MAX)

    result.num_ifgs, result.num_valid_ifgs = num_ifgs, num_valid
    result.flag_counts = dict(zip(FLAG_NAMES, flag_counts.tolist(), strict=True))
    result.mode_counts = dict(zip(MODES, mode_counts.tolist(), strict=True))
    result.max_saa_deg = float(np.rad2deg(max_saa_rad))
    return result
