"""Build (or check) the cached Static-sky cube for a configuration.

Usage::

    uv run python scripts/build_static_sky.py [config.toml]

Run it once before parallel simulations, so that MPI ranks only read the cache.
"""

import sys
import time
from pathlib import Path

from fossil_sim.config import load_config
from fossil_sim.input_sky import StaticSky
from fossil_sim.instrument import frequency_grid_ghz

path = (
    Path(sys.argv[1])
    if len(sys.argv) > 1
    else Path(__file__).parents[1] / "configs" / "baseline.toml"
)
config = load_config(path)
start = time.perf_counter()
sky = StaticSky(
    config.static_sky,
    frequency_grid_ghz(config.instrument),
    config.instrument.beam_fwhm_deg,
)
print(f"{sky.path}: {sky.alm.shape}, {time.perf_counter() - start:.0f} s")
