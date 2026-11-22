"""MPI equivalence check: ``mpirun -n 3 python scripts/mpi_check.py``.

The reduced hit maps must equal those of a serial run (ADR 0001); the
spectral maps, noise included, must agree to rounding (sums in another order).
"""

from dataclasses import replace
from pathlib import Path

import healpy as hp
import numpy as np
from mpi4py import MPI

from fossil_sim.config import load_config
from fossil_sim.simulation import run_simulation

config = load_config(Path(__file__).parents[1] / "configs" / "baseline.toml")
# ponytail: no Static sky (pysm3 downloads) nor zodi (slow) here; both are per-cycle
# functions of the configuration, so they do not change the MPI logic
config = replace(
    config,
    static_sky=replace(config.static_sky, models=()),
    zodi=replace(config.zodi, enabled=False),
)


def maps_agree(a: np.ndarray, b: np.ndarray, rtol: float = 1e-10) -> bool:
    """Same unseen pixels; seen ones equal to rtol of each map's largest value."""
    seen = a != hp.UNSEEN
    if not np.array_equal(seen, b != hp.UNSEEN):
        return False
    scale = np.abs(np.where(seen, a, 0.0)).max(axis=-1, keepdims=True)
    return bool(np.all(np.abs(np.where(seen, a - b, 0.0)) <= rtol * scale))


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
    assert config.noise.enabled  # noise realizations must match too
    for nside in parallel.spectral_maps.sky_count:
        p, s = parallel.spectral_maps, serial.spectral_maps
        assert np.array_equal(p.sky_count[nside], s.sky_count[nside])
        assert maps_agree(p.sky_map(nside), s.sky_map(nside))
        assert maps_agree(p.raw_maps(nside), s.raw_maps(nside))
    assert parallel.lost_fraction() == serial.lost_fraction()
    assert parallel.max_saa_deg == serial.max_saa_deg
    print("MPI run matches serial run")
