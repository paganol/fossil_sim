"""Generate and execute the example notebooks.

Usage::

    uv run python scripts/build_notebooks.py [output_dir]   # default: notebooks/

The notebooks are committed with their outputs (the documentation shows them
without re-executing); rerun this script when results change.
"""

import sys
from pathlib import Path

import nbformat
from nbclient import NotebookClient

out_dir = (
    Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parents[1] / "notebooks"
)
md, code = nbformat.v4.new_markdown_cell, nbformat.v4.new_code_cell

scanning = [
    md("""# FOSSIL scanning strategy

This notebook shows how the FOSSIL pointing is built: one Cooler cycle of
3 Rings × 23 Scans, followed by the Recycling window. See `CONTEXT.md` for the
vocabulary and `configs/baseline.toml` for the parameters and their sources."""),
    code("""from pathlib import Path

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

from fossil_sim.config import load_config
from fossil_sim.operations import MODES, cycle_timeline, operation_flags
from fossil_sim.orbit import L2Orbit
from fossil_sim.scanning import cycle_pointing
from fossil_sim.simulation import run_simulation

config = load_config(Path.cwd().parent / "configs" / "baseline.toml")
tl = config.timeline
print(f"Scan {tl.scan_duration_s:.0f} s, Ring {tl.scans_per_ring * tl.scan_duration_s / 3600:.1f} h, "
      f"Cooler cycle {tl.cycle_duration_s / 3600:.1f} h")"""),
    md("## Orbit and the anti-Sun direction"),
    code("""orbit = L2Orbit(config.orbit, config.simulation.mission_start_time, 365.25 * 86400)
days = np.linspace(0, 365, 2000)
pos_km, vel_km_s = orbit.position_velocity(days * 86400)
speed_km_s = np.linalg.norm(vel_km_s, axis=1)

fig, axes = plt.subplots(1, 2, figsize=(11, 3.5))
axes[0].plot(days, speed_km_s)
axes[0].set(xlabel="days since mission start", ylabel="|v| [km/s]", title="Barycentric speed")
lon_rate_deg_day = np.rad2deg(np.gradient(orbit.anti_sun_longitude_rad(days * 86400), days))
axes[1].plot(days, lon_rate_deg_day)
axes[1].set(xlabel="days since mission start", ylabel="deg / day", title="Anti-Sun longitude rate")
fig.tight_layout()"""),
    md("""## One Cooler cycle

Spin phase 0 (start of each Scan) is the point of the circle closest to the
north ecliptic pole. The 23 Scans of a Ring retrace the same circle; Rings are
stepped by one third of the solar motion over the cycle."""),
    code("""timeline = cycle_timeline(config, 0)
pointing = cycle_pointing(config, orbit, timeline)
mode, flags = operation_flags(config, timeline)
theta_rad, phi_rad = hp.vec2ang(pointing.boresight)
lat_deg, lon_deg = 90 - np.rad2deg(theta_rad), np.rad2deg(phi_rad)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))
first_scan = timeline.scan_id == 0
axes[0].plot(timeline.ifg_in_scan[first_scan], lat_deg[first_scan], ".", ms=2)
axes[0].set(xlabel="IFG index in Scan", ylabel="ecliptic latitude [deg]", title="First Scan")
equator = np.abs(lat_deg) < 2
for ring in range(tl.rings_per_cycle):
    sel = equator & (timeline.ring_in_cycle == ring) & ~timeline.in_recycling
    axes[1].plot(lon_deg[sel], lat_deg[sel], ".", ms=3, label=f"Ring {ring}")
axes[1].set(xlabel="ecliptic longitude [deg]", ylabel="latitude [deg]",
            title="Ring crossings of the ecliptic", xlim=(lon_deg[equator].min() - 0.3, lon_deg[equator].min() + 2))
axes[1].legend()
fig.tight_layout()
print("Ring steps [deg]:", np.rad2deg(np.arccos(np.sum(pointing.spin_axis[1:] * pointing.spin_axis[:-1], axis=1))))
print(f"Max SAA in this cycle: {np.rad2deg(pointing.saa_rad.max()):.3f} deg")"""),
    md("## Modes and flags along the cycle"),
    code("""hours = timeline.time_mid_s / 3600
fig, ax = plt.subplots(figsize=(11, 2.8))
ax.plot(hours, mode, drawstyle="steps-post", lw=0.5)
ax.set_yticks(range(len(MODES)), MODES)
ax.fill_between(hours, -0.5, 3.5, where=flags > 0, color="C3", alpha=0.3, step="post", label="flagged")
ax.set(xlabel="hours", title="Mode sequence (B10/B01 alternate Scan by Scan; B00 in Recycling)")
ax.legend(loc="upper left")
fig.tight_layout()"""),
    md("## Hit map after 30 days"),
    code("""result = run_simulation(config, num_of_cycles=21, spectral=False)
hp.mollview(result.hit_maps.total(64), title="Valid hits, 30 days, nside 64", unit="hits")
print({name: f"{100 * frac:.2f}%" for name, frac in result.lost_fraction().items()})"""),
]

ifg = [
    md("""# Ideal FTS interferograms

Illustration of the FOSSIL measurement principle with an ideal FTS (no
systematics): each output sees

$P(x) = \\sum_\\nu w(\\nu)\\,[I_1 + I_2 + s\\,(I_1 - I_2)\\cos(2\\pi\\nu x/c)]\\,\\Delta\\nu$,\n\nwith $x$ the optical path difference ($x = 2\\delta$ for mirror position $\\delta$ in the\nStep2 formula; the factor is a TBD config parameter),
$s = +1$ for the L output and $-1$ for the R output, and $w(\\nu)$ the
band of the detector unit (LF below 300 GHz, HF above). OPD range and
sampling are TBD placeholders from `configs/baseline.toml`."""),
    code("""from pathlib import Path

import matplotlib.pyplot as plt

from fossil_sim.config import load_config
from fossil_sim.ifg import differential_spectrum, interferogram, mirror_positions_m, opd_samples_m
from fossil_sim.instrument import detector_units, frequency_grid_ghz
from fossil_sim.spectra import blackbody_mjy_sr, y_distortion_mjy_sr

config = load_config(Path.cwd().parent / "configs" / "baseline.toml")
freq_ghz = frequency_grid_ghz(config.instrument)
opd_m = opd_samples_m(config.ifg)
delta_m = mirror_positions_m(config.ifg)
dus = {du.name: du for du in detector_units(config.instrument)}
print(f"{freq_ghz.size} frequency bins, {delta_m.size} samples, "
      f"OPD_max = {1e3 * opd_m.max():.1f} mm, mirror δ_max = {1e3 * delta_m.max():.1f} mm "
      f"(OPD / δ = {config.ifg.opd_per_mirror_position}, TBD)")"""),
    md("""## Mode B10: sky on Port 1, BBIR on Port 2

Sky = CMB blackbody plus a Compton-y distortion with $y = 10^{-5}$ (exaggerated
for visibility); BBIR at the 2.725 K set-point."""),
    code("""y = 1e-5
sky_mjy_sr = blackbody_mjy_sr(freq_ghz, config.sky.t_cmb_k) + y_distortion_mjy_sr(freq_ghz, y, config.sky.t_cmb_k)
bbir_mjy_sr = blackbody_mjy_sr(freq_ghz, config.bbir.setpoints_k[0])

fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharex=True)
for ax, name in zip(axes.ravel(), dus):
    ifg = interferogram(freq_ghz, sky_mjy_sr, bbir_mjy_sr, dus[name], opd_m)
    ax.plot(1e3 * delta_m, ifg - ifg.mean())
    ax.set(title=name, ylabel="P − ⟨P⟩ [MJy/sr GHz]")
for ax in axes[1]:
    ax.set_xlabel("mirror position δ [mm]")
fig.suptitle("Modulated part of the interferograms (L and R have opposite sign)")
fig.tight_layout()"""),
    md("""## Back to the spectrum

The recovered differential spectrum is dominated by the residual monopole
(the BBIR set-point differs from $T_\\mathrm{CMB}$ by 0.5 mK). Subtracting the
known blackbody difference leaves the y distortion."""),
    code("""blackbody_difference_mjy_sr = blackbody_mjy_sr(freq_ghz, config.sky.t_cmb_k) - bbir_mjy_sr
fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for name in ("LF-L", "HF-L"):
    du = dus[name]
    ifg = interferogram(freq_ghz, sky_mjy_sr, bbir_mjy_sr, du, opd_m)
    recovered = differential_spectrum(freq_ghz, ifg, du, opd_m)
    axes[0].plot(freq_ghz, recovered, "o", ms=3, label=name)
    axes[1].plot(freq_ghz, recovered - blackbody_difference_mjy_sr, "o", ms=3, label=name)
axes[0].plot(freq_ghz, sky_mjy_sr - bbir_mjy_sr, "k-", lw=1, label="sky − BBIR (input)")
axes[0].set(title="Recovered s (I₁ − I₂), Mode B10")
axes[1].plot(freq_ghz, y_distortion_mjy_sr(freq_ghz, y, config.sky.t_cmb_k), "k-", lw=1, label=f"y = {y:g} (input)")
axes[1].set(title="After removing the blackbody difference")
for ax in axes:
    ax.axhline(0, color="gray", lw=0.5)
    ax.set(xscale="log", xlabel="ν [GHz]", ylabel="ΔI [MJy/sr]")
    ax.legend()
fig.tight_layout()"""),
]

spectral = [
    md("""# Spectral maps: Isotropic backgrounds, Static sky and the kinematic dipole

The sky is the sum of the Isotropic backgrounds (CMB blackbody, y with
relativistic corrections, μ, CIB, extragalactic CO), observed with the exact
Boost from the solar plus orbital velocity (ADR 0002), and of the Static sky
(pysm3 foregrounds and CMB anisotropies, with the orbital Doppler term,
ADR 0003), plus the zodiacal light seen from the spacecraft (zodipy, planck18).
No noise yet. The sky estimate assumes perfect
calibration. The first run builds the Static-sky cache (about 30 s)."""),
    code("""from dataclasses import replace
from pathlib import Path

import healpy as hp
import matplotlib.pyplot as plt
import numpy as np

from fossil_sim.config import load_config
from fossil_sim.input_sky import isotropic_background_mjy_sr, solar_velocity_km_s
from fossil_sim.orbit import L2Orbit
from fossil_sim.simulation import run_simulation
from fossil_sim.spectra import (
    blackbody_mjy_sr,
    cib_mjy_sr,
    extragalactic_co_mjy_sr,
    mu_distortion_mjy_sr,
    relativistic_sz_correction_mjy_sr,
    y_distortion_mjy_sr,
)

config = load_config(Path.cwd().parent / "configs" / "baseline.toml")
config = replace(config, zodi=replace(config.zodi, nprocesses=4))
sky = config.sky"""),
    md("## Rest-frame spectra of the Isotropic backgrounds"),
    code("""freq_ghz = np.geomspace(30, 2000, 400)
t_k = sky.t_cmb_k
components = {
    "y (non-relativistic)": y_distortion_mjy_sr(freq_ghz, sky.y, t_k),
    "rSZ correction": relativistic_sz_correction_mjy_sr(freq_ghz, sky.y, sky.y_electron_temperature_kev, t_k),
    "μ": mu_distortion_mjy_sr(freq_ghz, sky.mu, t_k),
    "CIB": cib_mjy_sr(freq_ghz, sky.cib_amplitude, sky.cib_spectral_index, sky.cib_temperature_k, sky.cib_reference_freq_ghz),
    "extragalactic CO": extragalactic_co_mjy_sr(freq_ghz, sky.co_amplitude),
}
fig, ax = plt.subplots(figsize=(9, 4.5))
ax.plot(freq_ghz, blackbody_mjy_sr(freq_ghz, t_k), "k", label=f"CMB blackbody, {t_k} K")
for name, values in components.items():
    ax.plot(freq_ghz, np.abs(values), label=f"|{name}|")
ax.set(xscale="log", yscale="log", ylim=(1e-6, 1e3), xlabel="ν [GHz]", ylabel="|I| [MJy/sr]")
ax.legend(fontsize=8)
fig.tight_layout()"""),
    md("""## Dipole amplitude along the year

The observed dipole follows |v_sun + v_orbit(t)|: the orbital term modulates it
by about ±30 km/s over the year."""),
    code("""orbit = L2Orbit(config.orbit, config.simulation.mission_start_time, 365.25 * 86400)
days = np.linspace(0, 365, 500)
v_km_s = solar_velocity_km_s(sky) + orbit.position_velocity(days * 86400)[1]
dipole_mk = 1e3 * t_k * np.linalg.norm(v_km_s, axis=1) / 299792.458
fig, ax = plt.subplots(figsize=(7, 3))
ax.plot(days, dipole_mk)
ax.set(xlabel="days since mission start", ylabel="dipole amplitude [mK]")
fig.tight_layout()"""),
    md("""## Static sky

The cached alm cube sampled on a HEALPix grid (ecliptic coordinates), and its
orbital Doppler term: on the Galactic plane at THz frequencies it is far
larger than the y distortion."""),
    code("""from fossil_sim.input_sky import StaticSky, orbital_doppler_mjy_sr
from fossil_sim.instrument import frequency_grid_ghz

grid_ghz = frequency_grid_ghz(config.instrument)
static = StaticSky(config.static_sky, grid_ghz, config.instrument.beam_fwhm_deg)
nside = 64
vectors = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))).T
static_maps = static.sample(vectors)
v_orbit = np.broadcast_to(orbit.position_velocity(0.0)[1][0], vectors.shape)
doppler = orbital_doppler_mjy_sr(grid_ghz, static_maps, vectors, v_orbit)
fig = plt.figure(figsize=(12, 7))
for k, freq in enumerate((150, 1005)):
    i = np.argmin(np.abs(grid_ghz - freq))
    hp.mollview(static_maps[:, i], fig=fig, sub=(2, 2, k + 1), norm="hist",
                title=f"Static sky, {freq} GHz", unit="MJy/sr")
    hp.mollview(doppler[:, i], fig=fig, sub=(2, 2, k + 3),
                title=f"Orbital Doppler, {freq} GHz", unit="MJy/sr")"""),
    md("""## Zodiacal light

zodipy (planck18 model) seen from the spacecraft at mission start: brightest
along the ecliptic and towards the Sun, rising steeply with frequency. FOSSIL
only observes at about 90° elongation (the great circle perpendicular to the
anti-Sun Spin axis)."""),
    code("""from fossil_sim.input_sky import ZodiacalLight

zodi = ZodiacalLight(config, orbit, np.array([545.0, 1500.0]))
zodi_maps = zodi.evaluate(vectors, np.zeros(len(vectors)))
fig = plt.figure(figsize=(12, 3.5))
for k, freq in enumerate((545, 1500)):
    hp.mollview(zodi_maps[:, k], fig=fig, sub=(1, 2, k + 1),
                title=f"Zodiacal light, {freq} GHz (ecliptic)", unit="MJy/sr", norm="log")"""),
    md("## 10 days of spectral maps"),
    code("""result = run_simulation(config, num_of_cycles=7)
maps = result.spectral_maps
sky_map = maps.sky_map(64)
seen = sky_map[0] != hp.UNSEEN
i150 = np.argmin(np.abs(maps.freq_ghz - 150))
f = maps.freq_ghz[i150]
dbdt = (blackbody_mjy_sr(f, t_k + 1e-4) - blackbody_mjy_sr(f, t_k - 1e-4)) / 2e-4
delta_t_mk = np.full(sky_map.shape[1], hp.UNSEEN)
delta_t_mk[seen] = 1e3 * (sky_map[i150, seen] - isotropic_background_mjy_sr(f, sky)) / dbdt
hp.mollview(delta_t_mk, title=f"ΔT at {f:.0f} GHz, 10 days", unit="mK")
# Fit the dipole away from the Galactic plane (|b| > 30°)
gal_vectors = hp.Rotator(coord=["E", "G"])(np.array(hp.pix2vec(64, np.arange(hp.nside2npix(64)))))
high_latitude = np.abs(gal_vectors[2]) > np.sin(np.deg2rad(30))
_, dipole = hp.fit_dipole(np.where(seen & high_latitude, delta_t_mk, hp.UNSEEN))
print(f"Fitted dipole (|b| > 30°): {np.linalg.norm(dipole):.4f} mK")"""),
    md("""## Diagnostic maps

The same output (L) in B10 and B01 sees the sky with opposite sign; the
difference doubles the sky and cancels any common instrumental offset."""),
    code("""raw = maps.raw_maps(64)
fig = plt.figure(figsize=(12, 3.5))
for k, mode in enumerate(("B10", "B01")):
    hp.mollview(raw[0, k, i150], fig=fig, sub=(1, 2, k + 1), title=f"L output, {mode}, {f:.0f} GHz",
                unit="MJy/sr", min=-2, max=2)"""),
]

for name, cells in (
    ("fossil_sim_scanning.ipynb", scanning),
    ("fossil_sim_interferogram.ipynb", ifg),
    ("fossil_sim_spectral_maps.ipynb", spectral),
):
    nb = nbformat.v4.new_notebook(cells=cells)
    nb.metadata["kernelspec"] = {
        "name": "python3",
        "display_name": "Python 3",
        "language": "python",
    }
    NotebookClient(
        nb, timeout=600, resources={"metadata": {"path": str(out_dir)}}
    ).execute()
    nbformat.write(nb, out_dir / name)
    print("wrote", name)
