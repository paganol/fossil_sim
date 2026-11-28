"""Simulation configuration: TOML file -> frozen dataclasses.

Every numerical parameter lives in the TOML file (see
``configs/baseline.toml``), each with a comment giving its source. Unit
suffixes in the key names follow the litebird_sim convention.
"""

import hashlib
import json
import os
import subprocess
import tomllib
import warnings
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path

import healpy as hp
from astropy.time import Time

from .constants import C_LIGHT_M_S
from .version import __version__

MODES = ("B00", "B01", "B10", "B11")
RING_STEPPING_RULES = ("equal_per_cycle", "anchor_ring_start", "fixed")
# apodization of the IFG samples before the spectrum fit (fossil_sim.ifg)
APODIZATIONS = (
    "none",
    "norton_beer_weak",
    "norton_beer_medium",
    "norton_beer_strong",
    "hann",
    "triangular",
    "kaiser",
)

# Largest Spin-axis/anti-Sun offset produced by the "equal_per_cycle" Ring
# stepping (half a cycle of solar motion, ~1.44 deg / 2), plus margin
SAA_STEPPING_MARGIN_DEG = 0.75


@dataclass(frozen=True)
class SimulationConfig:
    master_seed: int
    mission_start: str
    duration_days: float
    # threads of the threaded libraries (ducc0, numba, BLAS), performance
    # only; None (key omitted): OMP_NUM_THREADS if set, else 1
    nthreads: int | None = None

    @property
    def resolved_nthreads(self) -> int:
        """Number of threads to use, see :func:`resolve_nthreads`."""
        return resolve_nthreads(self.nthreads)

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

    @property
    def observing_ifgs_per_cycle(self) -> int:
        """IFGs of a cycle outside the Recycling window."""
        return self.scans_per_ring * self.rings_per_cycle * self.ifgs_per_scan


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
    mode_switch_interval_s: float
    mode_pattern: tuple[str, ...]
    mode_switch_dead_time_s: float
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
    # width of the cosine band edges / edge frequency, 0 = top-hat (IFG domain only)
    band_edge_fraction: float = 0.0


@dataclass(frozen=True)
class SkyConfig:
    t_cmb_k: float
    include_cmb_blackbody: bool
    y: float
    y_electron_temperature_kev: float
    mu: float
    cib_amplitude: float
    cib_spectral_index: float
    cib_temperature_k: float
    cib_reference_freq_ghz: float
    co_amplitude: float
    solar_velocity_km_s: float
    solar_velocity_gal_lat_deg: float
    solar_velocity_gal_lon_deg: float
    include_orbital_velocity: bool
    apply_beam_to_dipole: bool


@dataclass(frozen=True)
class StaticSkyConfig:
    models: tuple[str, ...]
    nside: int
    lmax: int
    doppler: bool
    cache_dir: str
    sampling: str
    # relative accuracy of the low-rank factorization in frequency (sampling
    # "alm"); 0 = no factorization, one synthesis per frequency
    rank_epsilon: float = 1e-9


# pysm3 models that would double count the Isotropic backgrounds or are not
# suited to a discrete frequency grid (Galactic CO lines)
FORBIDDEN_STATIC_MODELS = ("dip", "co")
# Extragalactic pysm3 models: their monopole is removed (ADR 0002)
EXTRAGALACTIC_STATIC_MODELS = ("cib", "tsz", "ksz", "rg")


@dataclass(frozen=True)
class ZodiConfig:
    enabled: bool
    model: str
    extrapolate: bool
    gauss_quad_degree: int


def resolve_nthreads(nthreads: int | None = None) -> int:
    """``nthreads``, or ``OMP_NUM_THREADS`` when ``None`` (1 when unset).

    Of a nested ``OMP_NUM_THREADS`` list (e.g. ``"4,2"``) the first level is used.
    """
    if nthreads is not None:
        return nthreads
    env = os.environ.get("OMP_NUM_THREADS", "").split(",")[0].strip()
    if not env:
        return 1
    if not env.isdigit() or int(env) < 1:
        raise ValueError(f"OMP_NUM_THREADS must be a positive integer, got {env!r}")
    return int(env)


@dataclass(frozen=True)
class OutputConfig:
    save_tod: bool


@dataclass(frozen=True)
class NoiseConfig:
    enabled: bool
    nep_lf_w_sqrt_hz: float
    nep_hf_w_sqrt_hz: float
    etendue_cm2_sr: float
    absorptivity: float
    optical_efficiency: float
    modulation_efficiency: float
    recycling_noise_factor: float
    # "fixed": nep_lf/hf_w_sqrt_hz; "photon": photon NEP of each DU from its
    # optical loading (fossil_sim.noise.photon_nep_w_sqrt_hz, Coulon et al. 2026)
    nep_model: str = "fixed"
    # emitting optical elements (grey bodies) seen by the detectors, "photon" only
    optics_emissivity: tuple[float, ...] = ()
    optics_temperature_k: tuple[float, ...] = ()
    # 1/f detector noise on the IFG samples (ifg.simulate only), PSD
    # σ² [(f² + f_knee²) / (f² + f_min²)]^(α/2) (ducc0.misc.OofaNoise);
    # f_knee = 0: white noise
    fknee_lf_mhz: float = 0.0
    fknee_hf_mhz: float = 0.0
    one_over_f_alpha: float = 1.0
    one_over_f_fmin_hz: float = 1e-5


@dataclass(frozen=True)
class IFGConfig:
    max_opd_m: float
    opd_per_mirror_position: float
    samples_per_ifg: int
    simulate: bool
    sky_variation: bool
    # frequency grid the IFGs are generated on (None: the instrument grid)
    generation_freq_step_ghz: float | None = None
    generation_freq_max_ghz: float | None = None
    # apodization of the samples before the spectrum fit, one of APODIZATIONS
    apodization: str = "none"
    # β of the Kaiser window (apodization = "kaiser" only)
    kaiser_beta: float = 6.0


@dataclass(frozen=True)
class ResponseConfig:
    """True DU response (``fossil_sim.response``); DU order LF-L, HF-L, LF-R, HF-R."""

    gain: tuple[float, ...] = (1.0, 1.0, 1.0, 1.0)
    gain_rms: float = 0.0
    offset_emissivity: tuple[float, ...] = (0.0, 0.0, 0.0, 0.0)
    instrument_temperature_k: float = 2.725
    # gain drift g → g (1 + δg(t)), one value per IFG and DU, δg with PSD
    # σ² [(f² + f_knee²) / (f² + f_min²)]^(α/2) (ducc0.misc.OofaNoise);
    # σ = 0: no drift, f_knee = 0: white
    gain_drift_sigma: float = 0.0
    gain_drift_fknee_mhz: float = 0.0
    gain_drift_alpha: float = 1.0
    gain_drift_fmin_hz: float = 1e-5


CALIBRATION_SOLUTIONS = ("nominal", "true")
NEP_MODELS = ("fixed", "photon")


@dataclass(frozen=True)
class CalibrationConfig:
    solution: str = "nominal"


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
    sky: SkyConfig
    static_sky: StaticSkyConfig
    zodi: ZodiConfig
    output: OutputConfig
    noise: NoiseConfig
    response: ResponseConfig = field(default_factory=ResponseConfig)
    calibration: CalibrationConfig = field(default_factory=CalibrationConfig)
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
        """Metadata identifying a run: config hash, code version, git commit.

        ``git_sha`` ends with ``-dirty`` when the checkout has uncommitted changes.
        """
        repo = Path(__file__).parents[1]
        try:
            if not (repo / ".git").exists():  # installed copy: no checkout to ask
                raise OSError
            git_sha = subprocess.run(
                ["git", "rev-parse", "HEAD"],
                capture_output=True,
                text=True,
                check=True,
                cwd=repo,
            ).stdout.strip()
            # tracked files only, as `git describe --dirty`
            dirty = subprocess.run(
                ["git", "diff", "--quiet", "HEAD"], check=False, cwd=repo
            )
            if dirty.returncode == 1:
                git_sha += "-dirty"
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
    if config.simulation.duration_days <= 0:
        raise ValueError("simulation.duration_days must be > 0")
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
    bad_modes = set(ops.mode_pattern) | {ops.recycling_mode}
    if not ops.mode_pattern or not bad_modes <= set(MODES):
        raise ValueError(f"Modes must be in {MODES}, got {sorted(bad_modes)}")
    block_ifgs = ops.mode_switch_interval_s / tl.ifg_duration_s
    if block_ifgs < 1 or abs(block_ifgs - round(block_ifgs)) > 1e-9:
        raise ValueError(
            "operations.mode_switch_interval_s must be a positive multiple of "
            "timeline.ifg_duration_s"
        )
    if not 0 <= ops.mode_switch_dead_time_s < ops.mode_switch_interval_s:
        raise ValueError(
            "operations.mode_switch_dead_time_s must be in [0, mode_switch_interval_s)"
        )
    if ops.recycling_data not in ("flagged", "off"):
        raise ValueError("operations.recycling_data must be 'flagged' or 'off'")
    if config.ifg.opd_per_mirror_position <= 0 or config.ifg.max_opd_m <= 0:
        raise ValueError("ifg: max_opd_m and opd_per_mirror_position must be > 0")
    if config.ifg.samples_per_ifg < 2:
        raise ValueError("ifg.samples_per_ifg must be >= 2")
    opd_step_m = 2 * config.ifg.max_opd_m / (config.ifg.samples_per_ifg - 1)
    if (
        config.ifg.simulate
        and opd_step_m * config.instrument.freq_max_ghz * 1e9 > 0.5 * C_LIGHT_M_S
    ):
        raise ValueError(
            "ifg.samples_per_ifg too few to sample instrument.freq_max_ghz"
        )
    _validate_generation_grid(config, opd_step_m)
    if config.ifg.apodization not in APODIZATIONS:
        raise ValueError(f"ifg.apodization must be one of {APODIZATIONS}")
    if not 0 < config.ifg.kaiser_beta <= 50:
        raise ValueError("ifg.kaiser_beta must be in (0, 50]")
    if not config.ifg.simulate and config.ifg.apodization != "none":
        warnings.warn(
            "ifg.apodization acts only in the IFG domain (ifg.simulate = true): "
            "ignored in the spectral domain",
            stacklevel=2,
        )
    st = config.static_sky
    for model in st.models:
        if model.startswith(FORBIDDEN_STATIC_MODELS):
            raise ValueError(f"static_sky.models: '{model}' is not allowed (ADR 0002)")
    if st.sampling not in ("alm", "nearest"):
        raise ValueError("static_sky.sampling must be 'alm' or 'nearest'")
    if st.lmax < 1 or st.nside < 1:
        raise ValueError("static_sky: nside and lmax must be >= 1")
    if not 0 <= st.rank_epsilon < 1:  # also NaN
        raise ValueError("static_sky.rank_epsilon must be in [0, 1) (0 = off)")
    if config.zodi.gauss_quad_degree < 1:
        raise ValueError("zodi.gauss_quad_degree must be >= 1")
    if config.simulation.nthreads is not None and config.simulation.nthreads < 1:
        raise ValueError(
            "simulation.nthreads must be >= 1 (omit it: OMP_NUM_THREADS or 1)"
        )
    if not bbir.setpoints_k or len(bbir.port_offsets_k) != 2:
        raise ValueError("bbir: setpoints_k non-empty, port_offsets_k of length 2")
    if config.timeline.recycling_scans == 0 and len(set(bbir.setpoints_k)) > 1:
        raise ValueError(
            "bbir.setpoints_k: set-point changes need Recycling windows "
            "(timeline.recycling_scans > 0), the transient is not modelled"
        )
    if bbir.hf_cavity_k != 0:
        raise ValueError("bbir.hf_cavity_k: the HF cavity is not modelled yet, use 0")
    inst = config.instrument
    if inst.freq_step_ghz <= 0 or inst.freq_min_ghz >= inst.freq_max_ghz:
        raise ValueError(
            "instrument: freq_step_ghz > 0 and freq_min_ghz < freq_max_ghz"
        )
    if inst.freq_min_ghz < inst.band_min_ghz or inst.freq_max_ghz > inst.band_max_ghz:
        raise ValueError("instrument: the frequency grid must lie inside the band")
    if not inst.band_min_ghz < inst.dichroic_split_ghz < inst.band_max_ghz:
        raise ValueError("instrument: dichroic_split_ghz must lie inside the band")
    for nside in config.maps.nsides:
        if not hp.isnsideok(nside, nest=True):  # powers of 2 only
            raise ValueError(f"maps.nsides: {nside} is not a power of 2")
    # spectral resolution c / (2 x_max) must resolve the grid step
    if config.ifg.max_opd_m * inst.freq_step_ghz * 1e9 < 0.5 * C_LIGHT_M_S * (1 - 1e-9):
        raise ValueError("ifg.max_opd_m too short for instrument.freq_step_ghz")
    nc = config.noise
    if (
        min(
            nc.nep_lf_w_sqrt_hz,
            nc.nep_hf_w_sqrt_hz,
            nc.etendue_cm2_sr,
            nc.absorptivity,
            nc.optical_efficiency,
            nc.modulation_efficiency,
            nc.recycling_noise_factor,
        )
        <= 0
    ):
        raise ValueError("noise: all parameters must be > 0")
    if nc.modulation_efficiency > 1:
        raise ValueError("noise.modulation_efficiency must be in (0, 1]")
    if nc.nep_model not in NEP_MODELS:
        raise ValueError(f"noise.nep_model must be one of {NEP_MODELS}")
    if len(nc.optics_emissivity) != len(nc.optics_temperature_k):
        raise ValueError(
            "noise: optics_emissivity and optics_temperature_k need equal lengths"
        )
    if any(not 0 <= e <= 1 for e in nc.optics_emissivity) or any(
        not t > 0 for t in nc.optics_temperature_k
    ):
        raise ValueError("noise: optics emissivities in [0, 1] and temperatures > 0")
    fknees_hz = [nc.fknee_lf_mhz * 1e-3, nc.fknee_hf_mhz * 1e-3]
    if min(fknees_hz) < 0 or nc.one_over_f_alpha <= 0 or nc.one_over_f_fmin_hz <= 0:
        raise ValueError(
            "noise: fknee_lf/hf_mhz >= 0, one_over_f_alpha and _fmin_hz > 0"
        )
    if max(fknees_hz) > 0:
        if not config.ifg.simulate:
            raise ValueError("noise: 1/f noise (fknee_*_mhz > 0) needs ifg.simulate")
        sample_rate_hz = config.ifg.samples_per_ifg / config.timeline.ifg_duration_s
        if any(
            0 < f <= nc.one_over_f_fmin_hz or f >= sample_rate_hz / 2 for f in fknees_hz
        ):
            raise ValueError(
                "noise: fknee_*_mhz must lie between one_over_f_fmin_hz and the "
                "Nyquist frequency of the IFG samples"
            )
    rc = config.response
    if len(rc.gain) != 4 or len(rc.offset_emissivity) != 4:
        raise ValueError(
            "response: gain and offset_emissivity need one value per DU (4)"
        )
    if min(rc.gain) <= 0 or rc.gain_rms < 0 or rc.instrument_temperature_k <= 0:
        raise ValueError(
            "response: gain > 0, gain_rms >= 0 and instrument_temperature_k > 0"
        )
    if (
        min(rc.gain_drift_sigma, rc.gain_drift_fknee_mhz) < 0
        or rc.gain_drift_alpha <= 0
        or rc.gain_drift_fmin_hz <= 0
    ):
        raise ValueError(
            "response: gain_drift_sigma, _fknee_mhz >= 0; _alpha, _fmin_hz > 0"
        )
    fknee_hz = rc.gain_drift_fknee_mhz * 1e-3
    if fknee_hz > 0 and not (
        rc.gain_drift_fmin_hz < fknee_hz < 0.5 / config.timeline.ifg_duration_s
    ):
        raise ValueError(
            "response.gain_drift_fknee_mhz must lie between gain_drift_fmin_hz "
            "and the Nyquist frequency of the IFGs"
        )
    if config.calibration.solution not in CALIBRATION_SOLUTIONS:
        raise ValueError(f"calibration.solution must be one of {CALIBRATION_SOLUTIONS}")
    Time(config.simulation.mission_start, scale="tdb")  # raises if malformed


def _is_multiple(value: float, step: float) -> bool:
    """``value`` is an integer multiple of ``step``, to rounding."""
    ratio = value / step
    return abs(ratio - round(ratio)) < 1e-9 * max(1.0, abs(ratio))


def _validate_generation_grid(config: Config, opd_step_m: float) -> None:
    """IFG generation grid and band edges (docs: Interferogram)."""
    ifg, inst = config.ifg, config.instrument
    step, max_ghz = ifg.generation_freq_step_ghz, ifg.generation_freq_max_ghz
    if (step is None) != (max_ghz is None):
        raise ValueError(
            "ifg: set both generation_freq_step_ghz and generation_freq_max_ghz, "
            "or neither"
        )
    if step is not None:
        if step <= 0:
            raise ValueError("ifg.generation_freq_step_ghz must be > 0")
        # every instrument bin centre is a generation frequency
        if not (
            _is_multiple(inst.freq_step_ghz, step)
            and _is_multiple(inst.freq_min_ghz, step)
        ):
            raise ValueError(
                "ifg.generation_freq_step_ghz must divide instrument.freq_step_ghz "
                "and instrument.freq_min_ghz"
            )
        if max_ghz < inst.freq_max_ghz:
            raise ValueError(
                "ifg.generation_freq_max_ghz must be >= instrument.freq_max_ghz"
            )
        if ifg.simulate and opd_step_m * max_ghz * 1e9 > 0.5 * C_LIGHT_M_S:
            raise ValueError(
                "ifg.generation_freq_max_ghz above the Nyquist frequency of the "
                "OPD sampling, c / (2 OPD step)"
            )
    if not ifg.simulate and (step is not None or inst.band_edge_fraction > 0):
        warnings.warn(
            "ifg.generation_freq_* and instrument.band_edge_fraction act only in "
            "the IFG domain (ifg.simulate = true): ignored in the spectral domain",
            stacklevel=3,
        )
    f = inst.band_edge_fraction
    lo, split, hi = inst.band_min_ghz, inst.dichroic_split_ghz, inst.band_max_ghz
    if f < 0 or f * (lo + split) / 2 > split - lo or f * (split + hi) / 2 > hi - split:
        raise ValueError(
            "instrument.band_edge_fraction must be >= 0, with non-overlapping edges"
        )
