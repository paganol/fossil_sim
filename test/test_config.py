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
        ("operations", "mode_pattern", ["B10", "B2"]),
        ("operations", "mode_switch_interval_s", 1441.0),  # not whole IFGs
        ("operations", "mode_switch_interval_s", 0.0),
        ("operations", "mode_switch_dead_time_s", -1.0),
        ("operations", "recycling_data", "degraded"),
        ("bbir", "port_offsets_k", [0.0]),
        ("timeline", "scans_per_ring", 0),
        ("simulation", "unknown_key_s", 1.0),
        ("simulation", "duration_days", 0.0),
        ("operations", "mode_switch_dead_time_s", 16560.0),  # whole block
        ("bbir", "hf_cavity_k", 20.0),
        ("instrument", "freq_min_ghz", 15.0),  # grid outside the band
        ("instrument", "freq_step_ghz", 0.0),
        ("instrument", "dichroic_split_ghz", 2500.0),
        ("maps", "nsides", [64, 100]),
        ("ifg", "max_opd_m", 0.005),
        ("ifg", "samples_per_ifg", 200),  # OPD step above Nyquist at 2 THz
        ("noise", "modulation_efficiency", 1.5),
    ],
)
def test_invalid_values_are_rejected(section, key, value):
    values = baseline_dict()
    values[section][key] = value
    with pytest.raises(ValueError):
        config_from_dict(values)


def test_uncovered_poles_warn():
    values = baseline_dict()
    values["scanning"]["spin_boresight_angle_deg"] = 85.0
    with pytest.warns(UserWarning, match="poles"):
        config_from_dict(values)
