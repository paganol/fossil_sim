"""Generate and execute the example notebooks.

Usage::

    uv run python scripts/build_notebooks.py [output_dir]   # default: notebooks/

On every push to main, the ``notebooks`` workflow runs this script and
commits the notebooks with their outputs; GitHub and the documentation show
them executed. Pull requests only edit this script (CI checks it runs).
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
ax.set(xlabel="hours", title="Mode sequence (switch in the middle of each Ring; B00 in Recycling)")
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
band of the detector unit (LF below 200 GHz, HF above). OPD range and
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
    md("""## Detector noise: white and 1/f

With `ifg.simulate = true` the detector noise is drawn on the IFG samples,
optionally with a 1/f component (`noise.fknee_lf_mhz`, `fknee_hf_mhz`),
generated as in litebird_sim with `ducc0.misc.OofaNoise`. The mirror
modulation puts a bin at optical frequency $\\nu$ at the electrical frequency
$f = \\nu\\,v_\\mathrm{OPD}/c$, so the 1/f noise only reaches the bins whose
$f$ is near or below $f_\\mathrm{knee}$; the fit removes the constant of each
IFG, and with it the noise below the IFG rate. Left: power spectrum of the
sample noise of LF-L (unit white level), with the electrical frequencies of
the LF band shaded. Right: variance of each fitted bin relative to white
noise (`one_over_f_variance_factor`), the factor carried by the variance
maps."""),
    code("""from dataclasses import replace

import numpy as np
from scipy.signal import welch

from fossil_sim.noise import SampleNoise, ifg_solvers, one_over_f_variance_factor, sample_noise_sigma

simulated = replace(config, ifg=replace(config.ifg, simulate=True))
rate_hz = config.ifg.samples_per_ifg / config.timeline.ifg_duration_s
hz_per_ghz = 1e9 * 2 * opd_m.max() / config.timeline.ifg_duration_s / 299792458.0
lf = dus["LF-L"]
solvers = ifg_solvers(simulated, freq_ghz, list(dus.values()))
n_ifg = 2048

fig, axes = plt.subplots(1, 2, figsize=(12, 4))
for fknee_mhz, color in ((0.0, "k"), (500.0, "C0"), (2000.0, "C1")):
    cfg = replace(simulated, noise=replace(simulated.noise, fknee_lf_mhz=fknee_mhz, fknee_hf_mhz=fknee_mhz))
    samples = SampleNoise(cfg, lf, 0, 0, np.zeros(n_ifg, dtype=bool))[0:n_ifg]
    f_hz, psd = welch(samples.ravel() / sample_noise_sigma(cfg, lf), fs=rate_hz, nperseg=2**15)
    label = "white" if fknee_mhz == 0 else f"f_knee = {fknee_mhz / 1e3:g} Hz"
    axes[0].loglog(f_hz[1:], psd[1:] * rate_hz / 2, color=color, lw=1, label=label)
    if fknee_mhz > 0:
        for name in ("LF-L", "HF-L"):
            band = dus[name].band_mask(freq_ghz)
            axes[1].plot(freq_ghz[band], one_over_f_variance_factor(cfg, dus[name], solvers[name]),
                         "o-", ms=3, color=color, label=label if name == "LF-L" else None)
axes[0].axvspan(lf.band_min_ghz * hz_per_ghz, lf.band_max_ghz * hz_per_ghz, color="gray", alpha=0.2, label="LF band")
axes[0].axvline(1 / config.timeline.ifg_duration_s, color="gray", ls=":", label="IFG rate")
axes[0].set(xlabel="electrical frequency [Hz]", ylabel="PSD / white", title="Sample noise, LF-L")
axes[1].axhline(1, color="k", lw=0.5)
axes[1].set(xscale="log", xlabel="ν [GHz]", ylabel="Var / white", title="Variance of the fitted bins")
for ax in axes:
    ax.legend(fontsize=8)
fig.tight_layout()"""),
    md("""## Noise correlation between bins

The fit can correlate the frequency bins:
`sky_estimate_noise_covariance` gives the covariance per IFG of the sky
estimate, whose diagonal is the variance map; a pixel seen $N$ times has
$C/N$. Without apodization the bins are nearly independent; apodization
correlates neighbouring bins (Norton–Beer medium: 0.57; Kaiser with a large
β, more bins). 1/f noise adds at most a few percent between the lowest LF
bins. LF and HF bins come from different DUs and are uncorrelated: the
blocks break at the 200 GHz dichroic split (zoom on 105–300 GHz)."""),
    code("""from fossil_sim.noise import sky_estimate_noise_covariance

fig, axes = plt.subplots(1, 3, figsize=(13, 4))
cases = (("none", 8.0), ("norton_beer_medium", 8.0), ("kaiser", 12.0))
for ax, (apodization, beta) in zip(axes, cases):
    cfg = replace(simulated, ifg=replace(simulated.ifg, apodization=apodization, kaiser_beta=beta))
    covariance = sky_estimate_noise_covariance(cfg, freq_ghz, list(dus.values()))
    sigma = np.sqrt(np.diag(covariance))
    zoom = (freq_ghz >= 105) & (freq_ghz <= 300)  # around the 200 GHz dichroic split
    half = config.instrument.freq_step_ghz / 2
    extent = (freq_ghz[zoom][0] - half, freq_ghz[zoom][-1] + half) * 2
    image = ax.imshow((covariance / np.outer(sigma, sigma))[np.ix_(zoom, zoom)], cmap="RdBu_r",
                      vmin=-1, vmax=1, origin="lower", extent=extent)
    title = f"kaiser, β = {beta:g}" if apodization == "kaiser" else apodization
    ax.set(title=title, xlabel="ν [GHz]", ylabel="ν [GHz]")
fig.colorbar(image, ax=axes, label="correlation", shrink=0.8)"""),
]

spectral = [
    md("""# Spectral maps: Isotropic backgrounds, Static sky and the kinematic dipole

The sky is the sum of the Isotropic backgrounds (CMB blackbody, y with
relativistic corrections, μ, CIB, extragalactic CO), observed with the exact
Boost from the solar plus orbital velocity (ADR 0002), and of the Static sky
(pysm3 foregrounds and CMB anisotropies, with the orbital Doppler term,
ADR 0003), plus the zodiacal light seen from the spacecraft (zodipy, planck18).
White noise is added per IFG, bin and DU, at the
photon-noise level of the optical loading of each DU. The sky estimate assumes perfect calibration. The first run builds the Static-sky cache (about 30 s)."""),
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
config = replace(config, simulation=replace(config.simulation, nthreads=4))
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
first-order Doppler term: orbital velocity on the Galactic foregrounds (on the
plane at THz frequencies far larger than the y distortion), solar + orbital on
the CMB anisotropies (their Doppler modulation, dominant near the CMB peak)."""),
    code("""from fossil_sim.input_sky import StaticSky, static_doppler_mjy_sr
from fossil_sim.instrument import frequency_grid_ghz

grid_ghz = frequency_grid_ghz(config.instrument)
static = StaticSky(config.static_sky, grid_ghz, config.instrument.beam_fwhm_deg)
nside = 64
vectors = np.array(hp.pix2vec(nside, np.arange(hp.nside2npix(nside)))).T
static_maps = static.sample(vectors)
v_orbit = orbit.position_velocity(0.0)[1][0]
k_galactic, k_cmb = static.sample_doppler(vectors)
doppler = static_doppler_mjy_sr(k_galactic, vectors, v_orbit) + static_doppler_mjy_sr(
    k_cmb, vectors, solar_velocity_km_s(sky) + v_orbit)
fig = plt.figure(figsize=(12, 7))
for k, freq in enumerate((150, 1005)):
    i = np.argmin(np.abs(grid_ghz - freq))
    hp.mollview(static_maps[:, i], fig=fig, sub=(2, 2, k + 1), norm="hist",
                title=f"Static sky, {freq} GHz", unit="MJy/sr")
    hp.mollview(doppler[:, i], fig=fig, sub=(2, 2, k + 3),
                title=f"Doppler, {freq} GHz", unit="MJy/sr")"""),
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
    md("""## White noise

Noise per IFG and DU, flat in each band (photon noise of the loading of each
DU, `noise.nep_model = "photon"`), the
corresponding mission monopole sensitivity, and the noise per pixel of the
10-day map (square root of the variance map)."""),
    code("""from fossil_sim.instrument import detector_units
from fossil_sim.noise import noise_sigma_mjy_sr

tl = config.timeline
observing_scans = tl.scans_per_ring * tl.rings_per_cycle  # Recycling flagged
valid_fraction = observing_scans / (observing_scans + tl.recycling_scans)
valid_science = config.simulation.duration_days * 86400 / tl.ifg_duration_s * valid_fraction
sigma_ifg = np.sqrt(maps.noise_variance_per_ifg)
fig, ax = plt.subplots(figsize=(8, 4))
for du in detector_units(config.instrument)[:2]:
    band = du.band_mask(maps.freq_ghz)
    ax.plot(maps.freq_ghz[band], np.full(band.sum(), 1e6 * noise_sigma_mjy_sr(config, du)), label=f"{du.name}, per IFG")
ax.plot(maps.freq_ghz, 1e6 * sigma_ifg / np.sqrt(valid_science), "k", label="sky monopole, full mission")
ax.set(xscale="log", yscale="log", xlabel="ν [GHz]", ylabel="δS [Jy/sr]")
ax.legend()
fig.tight_layout()
noise_jy_sr = np.full(seen.size, hp.UNSEEN)
noise_jy_sr[seen] = 1e6 * np.sqrt(maps.sky_variance_map(64)[i150, seen])
hp.mollview(noise_jy_sr, title=f"Noise per pixel, {f:.0f} GHz, 10 days", unit="Jy/sr")"""),
    md("""## Gain drift

The gain of each DU can drift in time, $g \\to g\\,[1 + \\gamma(t)]$, one
$\\gamma$ per IFG and DU, drawn with `ducc0.misc.OofaNoise` at the IFG rate
(`response.gain_drift_sigma`, `gain_drift_fknee_mhz`, ...). The calibration
does not follow it, so the sky estimate carries
$\\gamma(t)\\,[s(I_1 - I_2) + o]$: in the science Modes the sky minus the
BBIR, dominated by the dipole, scaled by a drift correlated over many Scans.
Here $\\sigma_\\gamma = 10^{-4}$ per IFG and $f_\\mathrm{knee}$ = 10 mHz, on the
same 10 days without noise, Static sky and zodiacal light, so the residual
map is the drift alone."""),
    code("""from fossil_sim.response import gain_drift

drift_response = replace(config.response, gain_drift_sigma=1e-4, gain_drift_fknee_mhz=10.0)
gamma = gain_drift(replace(config, response=drift_response), 0, 0, 40000) - 1
hours = np.arange(gamma.size) * tl.ifg_duration_s / 3600
fig, ax = plt.subplots(figsize=(9, 3))
ax.plot(hours, 1e4 * gamma, lw=0.3, color="0.7", label="per IFG")
scan = tl.ifgs_per_scan
ax.plot(hours[: gamma.size // scan * scan].reshape(-1, scan).mean(axis=1),
        1e4 * gamma[: gamma.size // scan * scan].reshape(-1, scan).mean(axis=1), "C1o-", ms=3, label="mean per Scan")
ax.legend(fontsize=8)
ax.set(xlabel="hours since cycle start", ylabel="γ [10⁻⁴]", title="Gain drift of LF-L, first Cooler cycle")
fig.tight_layout()

quiet = replace(config, static_sky=replace(config.static_sky, models=()),
                zodi=replace(config.zodi, enabled=False), noise=replace(config.noise, enabled=False))
steady_maps = run_simulation(quiet, num_of_cycles=7).spectral_maps
drift_maps = run_simulation(replace(quiet, response=drift_response), num_of_cycles=7).spectral_maps
residual_jy_sr = np.full(seen.size, hp.UNSEEN)
residual_jy_sr[seen] = 1e6 * (drift_maps.sky_map(64)[i150, seen] - steady_maps.sky_map(64)[i150, seen])
hp.mollview(residual_jy_sr, title=f"Gain-drift residual, {f:.0f} GHz, 10 days", unit="Jy/sr")
print(f"rms of the residual: {np.std(residual_jy_sr[seen]):.1f} Jy/sr; "
      f"median noise per pixel: {np.median(noise_jy_sr[seen]):.0f} Jy/sr")"""),
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
