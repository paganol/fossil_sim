"""Streaming driver: one Cooler cycle at a time, maps accumulated on the fly.

Cycles are independent (ADR 0001): with an MPI communicator, rank ``r``
processes cycles ``r, r + size, ...`` and maps and statistics are summed at
the end, giving results identical to a serial run.
"""

from dataclasses import dataclass, field
from pathlib import Path

import h5py
import numpy as np

from .config import MODES, Config
from .input_sky import (
    StaticSky,
    ZodiacalLight,
    boosted_isotropic_background_mjy_sr,
    dipole_beam_factor,
    orbital_doppler_mjy_sr,
    solar_velocity_km_s,
)
from .instrument import detector_units, frequency_grid_ghz
from .mapmaking import DIAGNOSTIC_OUTPUTS, HitMaps, SpectralMaps
from .noise import du_noise_mjy_sr, sky_estimate_noise_variance
from .operations import (
    FLAG_NAMES,
    FLAG_SAA_EXCEEDED,
    bbir_temperatures_k,
    cycle_timeline,
    num_cycles,
    operation_flags,
)
from .orbit import L2Orbit
from .scanning import cycle_pointing, ring_spin_axes
from .signal import (
    bbir_spectra_mjy_sr,
    du_differential_spectra,
    port_spectra_mjy_sr,
    sky_estimate_mjy_sr,
)

# IFGs processed together when computing spectra (bounds the memory per step)
SIGNAL_BLOCK_IFGS = 8192

IFG_DTYPE = np.dtype(
    [
        ("time_start_s", np.float64),
        ("time_mid_s", np.float64),
        ("cycle_id", np.int32),
        ("ring_id", np.int32),  # -1 during the Recycling window
        ("scan_id", np.int32),
        ("ring_in_cycle", np.int8),  # last Ring during the Recycling window
        ("ifg_in_scan", np.int16),
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
    table["ring_in_cycle"] = timeline.ring_in_cycle
    table["ifg_in_scan"] = timeline.ifg_in_scan
    table["boresight"] = pointing.boresight
    table["mode"] = mode
    table["flags"] = flags
    table["saa_rad"] = pointing.saa_rad
    table["bbir_port_temperature_k"] = port_k
    table["bbir_hf_temperature_k"] = hf_k
    if config.operations.recycling_data == "off":
        table = table[~timeline.in_recycling]
    return table


class SignalModel:
    """Spectral-domain signal of each IFG, optionally with white noise."""

    def __init__(
        self, config: Config, orbit: L2Orbit, static_sky: StaticSky | None = None
    ) -> None:
        self.config = config
        self.orbit = orbit
        self.static_sky = static_sky
        self.freq_ghz = frequency_grid_ghz(config.instrument)
        self.zodi = (
            ZodiacalLight(config, orbit, self.freq_ghz) if config.zodi.enabled else None
        )
        self.dus = detector_units(config.instrument)
        self.solar_velocity_km_s = solar_velocity_km_s(config.sky)
        self.beam_factor = (
            dipole_beam_factor(config.instrument.beam_fwhm_deg)
            if config.sky.apply_beam_to_dipole
            else 1.0
        )

    def velocity_km_s(self, orbital_km_s: np.ndarray) -> np.ndarray:
        """Spacecraft velocity with respect to the CMB, (n, 3)."""
        velocity = np.broadcast_to(self.solar_velocity_km_s, orbital_km_s.shape).copy()
        if self.config.sky.include_orbital_velocity:
            velocity += orbital_km_s
        return velocity

    def zodi_for_cycle(self, cycle_id: int, table: np.ndarray) -> np.ndarray | None:
        """Zodiacal light of a whole cycle (evaluated Ring by Ring), or None if disabled."""
        if self.zodi is None:
            return None
        axes = ring_spin_axes(self.config, self.orbit, cycle_id)
        return self.zodi.for_cycle(cycle_id, axes, table)

    def sky_mjy_sr(
        self, table: np.ndarray, zodi: np.ndarray | None = None
    ) -> np.ndarray:
        """Total sky along the boresight, (n, n_freq).

        Isotropic backgrounds (boosted) + Static sky (with orbital Doppler) +
        zodiacal light (precomputed for the block, see :meth:`zodi_for_cycle`).
        """
        direction = table["boresight"]
        orbital_km_s = self.orbit.position_velocity(table["time_mid_s"])[1]
        sky = boosted_isotropic_background_mjy_sr(
            self.freq_ghz,
            direction,
            self.velocity_km_s(orbital_km_s),
            self.config.sky,
            self.beam_factor,
        )
        if self.static_sky is not None:
            static = self.static_sky.sample(direction)
            if self.config.static_sky.orbital_doppler:
                static += orbital_doppler_mjy_sr(
                    self.freq_ghz, static, direction, orbital_km_s
                )
            sky += static
        if zodi is not None:
            sky += zodi
        return sky

    def noise_for_cycle(
        self, cycle_id: int, table: np.ndarray
    ) -> dict[str, np.ndarray] | None:
        """White noise of each DU over a whole cycle, or None if disabled."""
        if not self.config.noise.enabled:
            return None
        in_recycling = table["ring_id"] == -1
        return {
            du.name: du_noise_mjy_sr(
                self.config,
                du,
                i,
                cycle_id,
                in_recycling,
                int(du.band_mask(self.freq_ghz).sum()),
            )
            for i, du in enumerate(self.dus)
        }

    def compute(
        self,
        table: np.ndarray,
        zodi: np.ndarray | None = None,
        noise: dict[str, np.ndarray] | None = None,
    ) -> tuple[dict[str, np.ndarray], np.ndarray, np.ndarray]:
        """DU differential spectra, sky estimate and BBIR spectra for a block of IFGs.

        ``noise`` (per DU, same shapes as the spectra) is added to the DU
        spectra, hence propagated to the sky estimate.
        """
        sky = self.sky_mjy_sr(table, zodi)
        bbir = bbir_spectra_mjy_sr(self.freq_ghz, table["bbir_port_temperature_k"])
        port1, port2 = port_spectra_mjy_sr(sky, bbir, table["mode"])
        du_spectra = du_differential_spectra(self.freq_ghz, port1, port2, self.dus)
        if noise is not None:
            du_spectra = {name: s + noise[name] for name, s in du_spectra.items()}
        estimate = sky_estimate_mjy_sr(
            self.freq_ghz, du_spectra, self.dus, bbir, table["mode"]
        )
        return du_spectra, estimate, bbir

    def output_spectra(
        self, du_spectra: dict[str, np.ndarray]
    ) -> dict[str, np.ndarray]:
        """Full-band differential spectrum of each output (L, R), assembled from its DUs."""
        result = {}
        for output in DIAGNOSTIC_OUTPUTS:
            full = np.full(
                (next(iter(du_spectra.values())).shape[0], self.freq_ghz.size), np.nan
            )
            for du in self.dus:
                if du.name.endswith(output):
                    full[:, du.band_mask(self.freq_ghz)] = du_spectra[du.name]
            result[output] = full
        return result


def write_cycle_tod(
    path: Path, config: Config, model: SignalModel, table, du_spectra
) -> None:
    """Write one Cooler cycle (IFG table + DU spectra) to HDF5."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        f.create_dataset("ifg_table", data=table)
        f.create_dataset("freq_ghz", data=model.freq_ghz)
        for du in model.dus:
            dataset = f.create_dataset(f"spectra/{du.name}", data=du_spectra[du.name])
            dataset.attrs["unit"] = "MJy/sr"
            dataset.attrs["freq_ghz"] = model.freq_ghz[du.band_mask(model.freq_ghz)]
        for key, value in config.provenance().items():
            f.attrs[key] = str(value)


@dataclass
class SimulationResult:
    """Accumulated products of a run."""

    hit_maps: HitMaps
    spectral_maps: SpectralMaps | None = None
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
    config: Config,
    num_of_cycles: int | None = None,
    comm=None,
    spectral: bool = True,
    output_dir: str | Path | None = None,
    static_sky: StaticSky | None = None,
) -> SimulationResult:
    """Run the mission (or its first ``num_of_cycles`` cycles) and accumulate maps.

    Parameters
    ----------
    config : Config
    num_of_cycles : int, optional
        Number of Cooler cycles to simulate; default: the whole mission.
    comm : mpi4py communicator, optional
        If given, cycles are distributed round-robin over ranks and the
        results are reduced on every rank.
    spectral : bool
        Compute the signal and the spectral maps (otherwise hit maps only).
    output_dir : path, optional
        Where per-cycle TOD files go when ``config.output.save_tod`` is true.
    static_sky : StaticSky, optional
        Pre-built Static sky; by default it is loaded from (or built into) the
        cache when ``config.static_sky.models`` is not empty. With MPI, rank 0
        builds it and the other ranks read the cache.
    """
    if config.output.save_tod and (output_dir is None or not spectral):
        raise ValueError("save_tod needs spectral=True and an output_dir")
    total_cycles = num_cycles(config) if num_of_cycles is None else num_of_cycles
    rank, size = (0, 1) if comm is None else (comm.rank, comm.size)
    orbit = L2Orbit(
        config.orbit,
        config.simulation.mission_start_time,
        total_cycles * config.timeline.cycle_duration_s,
    )
    if spectral and static_sky is None and config.static_sky.models:
        static_sky = _load_static_sky(config, comm)
    model = SignalModel(config, orbit, static_sky) if spectral else None
    result = SimulationResult(
        hit_maps=HitMaps(config.maps.nsides),
        spectral_maps=SpectralMaps(
            config.maps.nsides,
            model.freq_ghz,
            sky_estimate_noise_variance(config, model.freq_ghz, model.dus)
            if config.noise.enabled
            else None,
        )
        if spectral
        else None,
    )
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
        if model is not None:
            _process_signal(
                config, model, table, result.spectral_maps, cycle_id, output_dir
            )

    if comm is not None:
        from mpi4py import MPI

        summed = [*result.hit_maps.hits.values(), flag_counts, mode_counts]
        if result.spectral_maps is not None:
            sm = result.spectral_maps
            for maps in (sm.sky_sum, sm.sky_count, sm.raw_sum, sm.raw_count):
                summed += list(maps.values())
        for array in summed:
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


def _process_signal(config, model, table, spectral_maps, cycle_id, output_dir) -> None:
    """Signal and spectral maps of one cycle, in blocks of IFGs."""
    mode_names = np.asarray(MODES)
    tod_blocks = []
    zodi = model.zodi_for_cycle(cycle_id, table)
    noise = model.noise_for_cycle(cycle_id, table)
    for start in range(0, table.size, SIGNAL_BLOCK_IFGS):
        chunk = slice(start, start + SIGNAL_BLOCK_IFGS)
        block = table[chunk]
        block_zodi = None if zodi is None else zodi[chunk]
        block_noise = (
            None if noise is None else {name: n[chunk] for name, n in noise.items()}
        )
        du_spectra, estimate, _ = model.compute(block, block_zodi, block_noise)
        if config.output.save_tod:
            tod_blocks.append(du_spectra)
        valid = block["flags"] == 0
        outputs = model.output_spectra(du_spectra)
        spectral_maps.add(
            block["boresight"][valid],
            estimate[valid],
            {name: spectrum[valid] for name, spectrum in outputs.items()},
            mode_names[block["mode"][valid]],
        )
    if config.output.save_tod:
        du_spectra = {
            du.name: np.concatenate([b[du.name] for b in tod_blocks])
            for du in model.dus
        }
        path = Path(output_dir) / "tod" / f"cycle_{cycle_id:05d}.h5"
        write_cycle_tod(path, config, model, table, du_spectra)


def _load_static_sky(config: Config, comm) -> StaticSky:
    """Build the Static sky on rank 0 (if not cached), then load it everywhere."""
    args = (
        config.static_sky,
        frequency_grid_ghz(config.instrument),
        config.instrument.beam_fwhm_deg,
    )
    if comm is None:
        return StaticSky(*args)
    if comm.rank == 0:
        StaticSky(*args)
    comm.Barrier()
    return StaticSky(*args)
