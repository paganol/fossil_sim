import tomllib
from pathlib import Path

import pytest

from fossil_sim.config import config_from_dict, load_config

BASELINE = Path(__file__).parents[1] / "configs" / "baseline.toml"


def baseline_dict() -> dict:
    with open(BASELINE, "rb") as f:
        return tomllib.load(f)


def test_baseline_timeline():
    tl = load_config(BASELINE).timeline
    assert tl.scan_duration_s == 1440.0
    assert tl.scans_per_cycle == 23 * 3 + 18
    assert tl.cycle_duration_s / 3600 == pytest.approx(34.8)


def test_hash_is_stable_and_content_based():
    a = load_config(BASELINE)
    b = config_from_dict(baseline_dict())
    assert a.hash() == b.hash()
    values = baseline_dict()
    values["simulation"]["master_seed"] += 1
    assert config_from_dict(values).hash() != a.hash()
    assert set(a.provenance()) >= {"config_hash", "git_sha", "fossil_sim_version"}


@pytest.mark.parametrize(
    "section, key, value",
    [
        ("scanning", "precession_amplitude_deg", 9.5),  # SAA limit
        ("scanning", "ring_stepping", "continuous"),
        ("operations", "ring_mode_pattern", ["B10", "B2"]),
        ("operations", "recycling_data", "degraded"),
        ("bbir", "port_offsets_k", [0.0]),
        ("timeline", "scans_per_ring", 0),
        ("simulation", "unknown_key_s", 1.0),
    ],
)
def test_invalid_values_are_rejected(section, key, value):
    values = baseline_dict()
    values[section][key] = value
    with pytest.raises(ValueError):
        config_from_dict(values)
