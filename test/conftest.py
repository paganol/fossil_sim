from pathlib import Path

import pytest

from fossil_sim.config import load_config
from fossil_sim.orbit import L2Orbit

BASELINE_PATH = Path(__file__).parents[1] / "configs" / "baseline.toml"


@pytest.fixture(scope="session")
def config():
    return load_config(BASELINE_PATH)


@pytest.fixture(scope="session")
def orbit(config):
    return L2Orbit(
        config.orbit,
        config.simulation.mission_start_time,
        config.simulation.duration_days * 86400.0,
    )
