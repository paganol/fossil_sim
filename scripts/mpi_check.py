"""MPI equivalence check: ``mpirun -n 3 python scripts/mpi_check.py``.

The reduced hit maps must equal those of a serial run (ADR 0001).
"""

from dataclasses import replace
from pathlib import Path

import numpy as np
from mpi4py import MPI

from fossil_sim.config import load_config
from fossil_sim.simulation import run_simulation

config = load_config(Path(__file__).parents[1] / "configs" / "baseline.toml")
# ponytail: no Static sky here (pysm3 downloads); its MPI path is rank-0 build + barrier
config = replace(config, static_sky=replace(config.static_sky, models=()))
parallel = run_simulation(config, num_of_cycles=7, comm=MPI.COMM_WORLD)
if MPI.COMM_WORLD.rank == 0:
    serial = run_simulation(config, num_of_cycles=7)
    for nside in config.maps.nsides:
        assert np.array_equal(
            parallel.hit_maps.hits[nside], serial.hit_maps.hits[nside]
        )
        assert np.array_equal(
            parallel.hit_maps.first_hit_time_s[nside],
            serial.hit_maps.first_hit_time_s[nside],
        )
    assert parallel.lost_fraction() == serial.lost_fraction()
    assert parallel.max_saa_deg == serial.max_saa_deg
    print("MPI run matches serial run")
