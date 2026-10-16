"""Simulation configuration: TOML file -> frozen dataclasses.

Every numerical parameter lives in the TOML file (see
``configs/baseline.toml``), each with a comment giving its source. Unit
suffixes in the key names follow the litebird_sim convention.
"""

import hashlib
import json
import subprocess
import tomllib
import warnings
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path

from astropy.time import Time

from .version import __version__

MODES = ("B00", "B01", "B10", "B11")
RING_STEPPING_RULES = ("equal_per_cycle", "anchor_ring_start", "fixed")

# Largest Spin-axis/anti-Sun offset produced by the "equal_per_cycle" Ring
# stepping (half a cycle of solar motion, ~1.44 deg / 2), plus margin
SAA_STEPPING_MARGIN_DEG = 0.75


@dataclass(frozen=True)
class SimulationConfig:
    master_seed: int
    mission_start: str
    duration_days: float

    @property
    def mission_start_time(self) -> Time:
        """``mission_start`` as an astropy Time (TDB scale)."""
        return Time(self.mission_start, scale="tdb")


@dataclass(frozen=True)
class TimelineConfig:
    ifg_duration_s: float
    ifgs_per_scan: int
    scans_per_ring: int
    rings_per_cycle: int
    recycling_scans: int

    @property
    def scan_duration_s(self) -> float:
        return self.ifg_duration_s * self.ifgs_per_scan

    @property
    def scans_per_cycle(self) -> int:
        return self.scans_per_ring * self.rings_per_cycle + self.recycling_scans

    @property
    def ifgs_per_cycle(self) -> int:
        return self.scans_per_cycle * self.ifgs_per_scan

    @property
    def cycle_duration_s(self) -> float:
        return self.scans_per_cycle * self.scan_duration_s


@dataclass(frozen=True)
class ScanningConfig:
    spin_boresight_angle_deg: float
    ring_stepping: str
    fixed_ring_step_deg: float
    precession_amplitude_deg: float
    precession_period_days: float
    precession_phase_deg: float
    max_saa_deg: float


@dataclass(frozen=True)
class OrbitConfig:
    earth_l2_distance_km: float
    lissajous_radius1_km: float
    lissajous_radius2_km: float
    lissajous_freq1_cycles_yr: float
    lissajous_freq2_cycles_yr: float
    lissajous_phase_deg: float


@dataclass(frozen=True)
class OperationsConfig:
    ring_mode_pattern: tuple[str, ...]
    mode_switch_dead_ifgs: int
    recycling_data: str
    recycling_mode: str


@dataclass(frozen=True)
class BBIRConfig:
    setpoints_k: tuple[float, ...]
    port_offsets_k: tuple[float, float]
    hf_cavity_k: float


@dataclass(frozen=True)
class MapsConfig:
    nsides: tuple[int, ...]


@dataclass(frozen=True)
class InstrumentConfig:
    freq_min_ghz: float
    freq_max_ghz: float
    freq_step_ghz: float
    band_min_ghz: float
    band_max_ghz: float
    dichroic_split_ghz: float
    beam_fwhm_deg: float


@dataclass(frozen=True)
class IFGConfig:
    max_opd_m: float
    opd_per_mirror_position: float
    samples_per_ifg: int


@dataclass(frozen=True)
class Config:
    simulation: SimulationConfig
    timeline: TimelineConfig
    scanning: ScanningConfig
    orbit: OrbitConfig
    operations: OperationsConfig
    bbir: BBIRConfig
    maps: MapsConfig
    instrument: InstrumentConfig
    ifg: IFGConfig
    source_path: str | None = field(default=None, compare=False)

    def __post_init__(self) -> None:
        validate(self)

    def to_dict(self) -> dict:
        """Configuration as a plain dict (without ``source_path``)."""
        result = asdict(self)
        result.pop("source_path")
        return result

    def hash(self) -> str:
        """SHA-256 of the configuration content, independent of key order."""
        text = json.dumps(self.to_dict(), sort_keys=True)
        return hashlib.sha256(text.encode()).hexdigest()

    def provenance(self) -> dict:
        """Metadata identifying a run: config hash, code version, git commit."""
        try:
            git_sha = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                check=True,
                cwd=Path(__file__).parent,
            ).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            git_sha = "unknown"
        return {
            "config_hash": self.hash(),
            "config_path": self.source_path,
            "fossil_sim_version": __version__,
            "git_sha": git_sha,
        }


def _build(cls, values: dict):
    """Build a section dataclass, turning lists into tuples."""
    known = {f.name for f in fields(cls)}
    unknown = set(values) - known
    if unknown:
        raise ValueError(f"Unknown keys in [{cls.__name__}]: {sorted(unknown)}")
    return cls(**{k: tuple(v) if isinstance(v, list) else v for k, v in values.items()})


def config_from_dict(values: dict, source_path: str | None = None) -> Config:
    """Build a validated :class:`Config` from a nested dict (as read from TOML)."""
    sections = {
        f.name: _build(f.type, values.get(f.name, {}))
        for f in fields(Config)
        if f.name != "source_path"
    }
    return Config(**sections, source_path=source_path)


def load_config(path: str | Path) -> Config:
    """Read and validate a TOML configuration file."""
    with open(path, "rb") as f:
        return config_from_dict(tomllib.load(f), source_path=str(path))


def validate(config: Config) -> None:
    """Raise ``ValueError`` on inconsistent or out-of-range parameters."""
    tl, sc, ops, bbir = config.timeline, config.scanning, config.operations, config.bbir
    for name in ("ifgs_per_scan", "scans_per_ring", "rings_per_cycle"):
        if getattr(tl, name) < 1:
            raise ValueError(f"timeline.{name} must be >= 1")
    if tl.recycling_scans < 0 or tl.ifg_duration_s <= 0:
        raise ValueError("timeline: recycling_scans >= 0 and ifg_duration_s > 0")
    if sc.ring_stepping not in RING_STEPPING_RULES:
        raise ValueError(f"scanning.ring_stepping must be one of {RING_STEPPING_RULES}")
    if sc.ring_stepping != "fixed":
        worst_saa_deg = sc.precession_amplitude_deg + SAA_STEPPING_MARGIN_DEG
        if worst_saa_deg > sc.max_saa_deg:
            raise ValueError(
                f"Precession amplitude + stepping offset ({worst_saa_deg} deg) "
                f"exceeds max_saa_deg ({sc.max_saa_deg} deg)"
            )
    pole_gap_deg = abs(90.0 - sc.spin_boresight_angle_deg) - sc.precession_amplitude_deg
    if pole_gap_deg > 0:
        warnings.warn(
            f"Ecliptic poles never observed: |90° − spin_boresight_angle_deg| exceeds "
            f"precession_amplitude_deg by {pole_gap_deg} deg",
            stacklevel=2,
        )
    if sc.precession_amplitude_deg > 0 and sc.precession_period_days <= 0:
        raise ValueError("scanning.precession_period_days must be > 0")
    bad_modes = set(ops.ring_mode_pattern) | {ops.recycling_mode}
    if not ops.ring_mode_pattern or not bad_modes <= set(MODES):
        raise ValueError(f"Modes must be in {MODES}, got {sorted(bad_modes)}")
    if ops.mode_switch_dead_ifgs < 0:
        raise ValueError("operations.mode_switch_dead_ifgs must be >= 0")
    if ops.recycling_data not in ("flagged", "off"):
        raise ValueError("operations.recycling_data must be 'flagged' or 'off'")
    if config.ifg.opd_per_mirror_position <= 0 or config.ifg.max_opd_m <= 0:
        raise ValueError("ifg: max_opd_m and opd_per_mirror_position must be > 0")
    if not bbir.setpoints_k or len(bbir.port_offsets_k) != 2:
        raise ValueError("bbir: setpoints_k non-empty, port_offsets_k of length 2")
    Time(config.simulation.mission_start, scale="tdb")  # raises if malformed
