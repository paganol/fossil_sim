# Changelog

Changes are grouped by version, newest first, with links to the pull requests
that introduced them. Versions before 0.1.5 were development versions, without
tagged releases.

## HEAD

Changes on `main` that are not in a release yet.

## 0.1.11 (2026-11-28)

- Apodization: the plain least-squares fit is now followed by the line shape
  of the window at the bin spacing, ℓ(m) = ∫₀¹ A(u) cos(π m u) du (new
  `fossil_sim.ifg.apodization_line_shape`), instead of weighting the samples
  before the projection. On the sampled OPD grid the cosines are not exactly
  orthogonal (both endpoints sampled, c / (2 x_max) = 14.99 GHz), and the old
  solver leaked the first bin of each DU band into every bin with sign
  (−1)^k: a smooth spectrum came back with an alternating residual of about
  2 × 10⁻³ (10³–10⁴ Jy/sr on region-averaged step2 maps). Line shapes, noise
  per bin and correlations are unchanged to 0.1 GHz / two digits; results
  without apodization are unchanged.
  Reference PR: [#43](https://github.com/paganol/fossil_sim/pull/43).

## 0.1.10 (2026-11-28)

- Faster IFG-domain simulation with the low-rank Static sky (about 2.7× on
  the IFG synthesis, 6.6 → 4.4 s per baseline cycle with 4 threads, same
  results to rounding): the window interpolation takes one matrix product per
  group of samples, with the weights folded into the modulation, and the
  terms of each IFG are combined in a numba kernel parallel over the IFGs
  (`simulation.nthreads`). The sky windows are IFG-major,
  (n_ifg, n_window, …), in `SignalModel.static_window`,
  `du_interferograms_low_rank` and `_interpolated_modulation`.
  Reference PR: [#40](https://github.com/paganol/fossil_sim/pull/40).

- One thread setting for the threaded libraries: `simulation.nthreads` sets the
  concurrent ducc0 syntheses of the Static sky, the numba and BLAS threads
  (numpy matrix products, through the new dependency `threadpoolctl`; not
  Accelerate on macOS) and the chunks of IFG white noise, drawn concurrently
  with the same realization for any number of threads. Omitted (default), it
  is `OMP_NUM_THREADS`, or 1 when unset. zodipy always runs serially
  (`nprocesses = 1`): its process pool, rebuilt at every call, costs more than
  it saves, and the default single-pass `planck18` evaluation does not use it.
  **Config change:** `static_sky.nthreads` and `zodi.nprocesses` are removed
  (unknown keys now raise); the Static sky no longer uses all cores by
  default. `StaticSky` takes `nthreads`.
  Reference PR: [#40](https://github.com/paganol/fossil_sim/pull/40).

- Noise covariance between frequency bins of the sky map:
  `sky_estimate_noise_covariance` gives the covariance per IFG (diagonal in
  the spectral domain; in the IFG domain propagated through the fit, with
  apodization and 1/f noise), the variance maps being its diagonal. A pixel
  seen N times has C / N: `SpectralMaps.sky_covariance` and
  `sky_correlation`, and `sky_noise_covariance.h5` among the outputs
  (covariance per IFG, frequencies, hits per nside). `SpectralMaps` takes
  `noise_covariance_per_ifg`; `noise_variance_per_ifg` still works alone.
  Reference PR: [#39](https://github.com/paganol/fossil_sim/pull/39).

- Optional apodization of the IFG samples before the spectrum fit
  (`ifg.simulate = true`): Norton–Beer weak, medium, strong, Hann,
  triangular (sinc² line shape, FWHM 1.772 c / (2 x_max)) or Kaiser (β in
  `ifg.kaiser_beta`, default 6), with
  A(0) = 1. The constant of the plain fit is removed, then the apodized
  samples are projected on the cosines; the noise variance and the 1/f
  factor follow the new solver. It suppresses the sinc sidelobes (−22 % →
  1.4 % with Norton–Beer medium) at the cost of a wider line shape (18 →
  25 GHz), correlated neighbouring bins and low edge bins of each DU; docs
  (Interferogram) give its effect on the band-edge systematics.
  **Config change:** new keys `ifg.apodization` (default `"none"`, results
  unchanged) and `ifg.kaiser_beta` (default 6, used only by `"kaiser"`);
  the configuration hash changes for every configuration.
  Reference PR: [#35](https://github.com/paganol/fossil_sim/pull/35).

- Gain drift in time, g → g (1 + γ(t)): one γ per IFG and DU, the same in
  every bin, drawn as in litebird_sim with `ducc0.misc.OofaNoise` at the IFG
  rate (white rms per IFG, knee, slope, f_min), independent DUs, one stream
  per Cooler cycle. It multiplies signal and offset, not the noise, in the
  spectral and IFG domains; the calibration solution does not follow it, so
  it leaves stripes in the maps. The TOD files store 1 + γ under
  `response/gain_drift/`.
  **Config change:** new keys `response.gain_drift_sigma` (default 0, no
  drift), `response.gain_drift_fknee_mhz` (0), `response.gain_drift_alpha`
  (1), `response.gain_drift_fmin_hz` (1e-5); the configuration hash changes
  for every configuration.
  Reference PR: [#34](https://github.com/paganol/fossil_sim/pull/34).

- 1/f detector noise on the IFG samples (`ifg.simulate = true`), generated
  as in litebird_sim with `ducc0.misc.OofaNoise` on the white stream of each
  Cooler cycle: PSD σ² [(f² + f_knee²) / (f² + f_min²)]^(α/2), knee per band.
  The variance maps include its effect on each fitted bin. The FTS
  modulation keeps it out of the maps' large scales: the spectra of
  successive IFGs stay uncorrelated, only the bins at electrical frequencies
  ν v_OPD / c near f_knee get noisier (LF 30 GHz ≈ 1 Hz in the baseline).
  **Config change:** new keys `noise.fknee_lf_mhz`, `noise.fknee_hf_mhz`
  (default 0, white noise), `noise.one_over_f_alpha` (1),
  `noise.one_over_f_fmin_hz` (1e-5); the configuration hash changes for
  every configuration, the noise realizations do not.
  Reference PR: [#33](https://github.com/paganol/fossil_sim/pull/33).

- Static sky sampled through a low-rank factorization in frequency of each
  alm cube (SVD, every frequency map kept within `static_sky.rank_epsilon`
  of its rms, default 10⁻⁹ = the synthesis accuracy; 0 = off, one synthesis
  per frequency): 87 syntheses per cycle instead of 2070 on the
  690-frequency generation grid, and the cubes are no longer loaded by each
  rank (about 2.4 GB less per rank). The factorization is cached next to the
  cubes, one file per threshold (`*_low_rank_1e-09.h5`), built by
  `scripts/build_static_sky.py` or by the first run. Same outputs to the
  synthesis accuracy.
  **Config change:** new key `static_sky.rank_epsilon` (optional, default
  1e-9); the configuration hash changes for every configuration (the
  Static-sky cache key does not).
  Reference PR: [#29](https://github.com/paganol/fossil_sim/pull/29).
- IFG domain with sky variation: the Static sky enters the interferograms
  through the IFGs of its low-rank basis spectra, computed once per cycle and
  DU, instead of sky windows on the generation grid. Same outputs to
  10⁻¹⁰; the full-mission run on the 3 GHz generation grid (690 frequencies,
  cosine edges, 48 ranks on G100) goes from 32 to 10 min.
  Reference PR: [#29](https://github.com/paganol/fossil_sim/pull/29).
- IFG domain: the matrix cos(2πνx/c) on the band bins of each DU is built
  once by `SignalModel` instead of at every block (248 times per cycle).
  Identical outputs; a fine-grid cycle on one rank goes from 16.6 to 15.8 s.
  Reference PR: [#30](https://github.com/paganol/fossil_sim/pull/30).
- Photon noise from the optical loading: `noise.nep_model = "photon"`
  computes the NEP of each DU from the power it receives (Isotropic
  backgrounds, BBIR, grey-body optics `optics_emissivity` /
  `optics_temperature_k`), shot and wave noise, as in Coulon, Aghanim &
  Maffei 2026 (arXiv:2609.17829, Eq. 15–16). Their FOSSIL reference
  sensitivity (1.14 / 0.91 Jy/sr per 15 GHz channel) is reproduced to 1 %.
  New `configs/step2.toml`: FOSSIL bands 50–300–2000 GHz, 4-yr mission,
  photon noise with Coulon's reference values and `modulation_efficiency =
  0.5`. `baseline.toml` is unchanged (`nep_model = "fixed"` is the default).
  **Config change:** new optional `[noise]` keys `nep_model`,
  `optics_emissivity`, `optics_temperature_k`; the configuration hash
  changes for every configuration.
  Reference PR: [#31](https://github.com/paganol/fossil_sim/pull/31).

- The baseline uses the photon noise from the optical loading
  (`nep_model = "photon"`, Coulon et al. reference values: η = 0.5,
  AΩ = 3 cm² sr, optics 1 % at 4.5 K, `modulation_efficiency = 0.5`) instead
  of the fixed NEPs of the M8 proposal, whose HF value (1.3e-14 W/√Hz) is not
  a photon NEP. Bands and mission are unchanged (30–200–2000 GHz, 3.5 yr).
  **Output change:** the baseline noise changes (per IFG and DU, LF 1.1e4 →
  7.6e3 Jy/sr, HF 3.0e5 → 1.2e4 Jy/sr).
  Reference PR: [#32](https://github.com/paganol/fossil_sim/pull/32).

## 0.1.9 (2026-11-25)

- Add a response for each DU, d = g (s (I₁ − I₂) + o) + n: gain per DU with a
  per-bin Gaussian fluctuation fixed by the seed, and an offset from the
  emission imbalance between the arms, ε B_ν(T_inst) (new `[response]`
  section; the baseline is the ideal instrument, so outputs are unchanged).
  The sky estimate and its noise variance use a calibration solution,
  `calibration.solution = "nominal"` (g = 1, o = 0) or `"true"` (oracle).
  Both domains: in the IFG domain the gain also scales the unmodulated term.
  TOD files store the uncalibrated spectra and the true response under
  `response/`.
  **Config change:** new `[response]` and `[calibration]` sections (optional,
  defaults as the baseline); the configuration hash changes for every
  configuration.
  Reference PR: [#24](https://github.com/paganol/fossil_sim/pull/24).
- IFG domain without sky variation (`ifg.simulate = true`,
  `ifg.sky_variation = false`): the Static sky is read from the Ring circles
  at the IFG centres, as in the spectral domain, instead of being synthesized
  from the alm at every IFG. Same outputs (to the synthesis accuracy), a
  full-mission baseline run no longer exceeds 1 h on 48 ranks.
  Reference PR: [#25](https://github.com/paganol/fossil_sim/pull/25).
- IFG generation on a finer, wider frequency grid than the fit:
  `ifg.generation_freq_step_ghz` and `ifg.generation_freq_max_ghz` (e.g. 3 GHz
  up to 3825 GHz, below the Nyquist frequency of the OPD sampling), and band
  transmissions with cosine edges, `instrument.band_edge_fraction`. The sky
  entering the instrument (and the Static-sky cache) is evaluated on the
  generation grid; the fit, the sky estimate and the maps stay on the
  instrument bins. Shows bin averaging, leakage of the unapodized fit and
  out-of-band signal (docs: Interferogram). IFG domain only; off in the
  baseline, so outputs are unchanged.
  **Config change:** new optional keys; the configuration hash changes for
  every configuration.
  Reference PR: [#26](https://github.com/paganol/fossil_sim/pull/26).
- Faster IFG-domain simulation: the modulated term of each sample is computed
  only on the four sky-window points it interpolates from, and once per IFG
  when the window is constant (sky of the IFG centre, both Ports on the BBIR);
  fewer array copies in `ifg.du_interferograms`. Same results to rounding; a
  baseline cycle on the 690-frequency generation grid goes from 135 s to 60 s
  on CINECA Galileo100.
  Reference PR: [#27](https://github.com/paganol/fossil_sim/pull/27).

## 0.1.8 (2026-11-22)

- Harden configuration and numerical edge cases found in the v0.1.7 review;
  the baseline outputs are unchanged.
  Reference PR: [#21](https://github.com/paganol/fossil_sim/pull/21).
  - Build the frequency grid from an integer bin count, so rounding in
    `(max - min) / step` no longer adds a bin (e.g. 1–1.3 GHz in 0.1 GHz steps
    gave four bins).
  - Reject a NaN `freq_step_ghz` in the interferogram helpers, and an explicit
    width that differs from the spacing of a multi-bin grid.
  - Evaluate boosted backgrounds directly when the range of ν′/ν exceeds the
    Chebyshev accuracy range (far above physical velocities).
  - **Config change:** reject several `bbir.setpoints_k` with
    `timeline.recycling_scans = 0`: set-point changes need a Recycling window.
  - Write the Static-sky cache through a per-process temporary file, so
    concurrent jobs sharing a cache directory do not collide.
  - Pin `zodipy<1.2` (the single-pass zodiacal evaluation uses its private API)
    and compare its shared parameters element-wise.
  - Reject `ifg.samples_per_ifg < 2`.
  - Record the full commit hash in `git_sha` (`git rev-parse HEAD`, plus
    `-dirty` for uncommitted changes to tracked files) instead of the
    `git describe` string, so reports show the commit rather than the last tag.
- Extend test coverage: an end-to-end IFG-domain test that noisy maps follow
  the propagated variance, Mode-switch flags with and without Recycling, and
  an MPI check that also compares sky and raw spectral maps with noise. The
  MPI job is now required in CI.
  Reference PR: [#21](https://github.com/paganol/fossil_sim/pull/21).
- Compute the least-squares solver of each DU once per run in IFG-domain
  simulations and share it with the noise variance; it is public as
  `fossil_sim.ifg.differential_spectrum_solver` (results unchanged).
  Reference PR: [#21](https://github.com/paganol/fossil_sim/pull/21).
- Reuse one thread pool for the Static-sky synthesis, sized by the new
  optional `static_sky.nthreads` (0 = all cores, the previous behaviour).
  **Config change:** the configuration hash changes for every configuration.
  Reference PR: [#21](https://github.com/paganol/fossil_sim/pull/21).
- **TOD format change:** store `ring_in_cycle` and `ifg_in_scan` as `int32`,
  like the other indices of the IFG table.
  Reference PR: [#21](https://github.com/paganol/fossil_sim/pull/21).
- Derive the valid-data fraction in the spectral-maps notebook from the
  timeline instead of a hard-coded 69/87, and select hit maps by Mode name in
  the hit-map example.
  Reference PR: [#21](https://github.com/paganol/fossil_sim/pull/21).

## 0.1.7 (2026-11-03)

- Add IFG-domain simulation with `ifg.simulate = true`: generate interferogram
  samples, add sample noise and recover differential spectra by least squares.
  Optional `ifg.sky_variation` follows the Static sky during each IFG using
  Ring-circle interpolation. The spectral-domain path remains the default.
  Reference PR: [#15](https://github.com/paganol/fossil_sim/pull/15).
- Make IFG sample noise reproducible regardless of block slicing, processing
  order or repeated calls. Apply the OPD Nyquist check only in IFG-domain
  simulations and clarify the mirror trajectory convention.
  Reference PR: [#16](https://github.com/paganol/fossil_sim/pull/16).
- Speed up spectral-domain simulations by sampling the Static sky once per Ring
  circle and evaluating zodiacal-light quadrature once for all frequencies.
  Reference PR: [#17](https://github.com/paganol/fossil_sim/pull/17).
- Further reduce spectral-domain cost with Chebyshev interpolation of boosted
  smooth backgrounds, matrix products for zodiacal emission and threaded
  Static-sky synthesis. Bound zodiacal-light evaluation memory by processing
  directions in blocks, including unsorted times and single instants.
  Reference PR: [#18](https://github.com/paganol/fossil_sim/pull/18).
- Fix IFG-domain variance maps to propagate sample noise through the actual
  least-squares solver, including finite-sampling effects; the spectral-domain
  noise model is unchanged. Store local timeline indices as `int64`, preserve
  Mode-switch dead time across Cooler cycles without Recycling, and support
  single-bin interferograms through `freq_step_ghz`.
  Reference PR: [#19](https://github.com/paganol/fossil_sim/pull/19).
- Organize the documentation into six top-level sections with separate concept
  pages for the sky components, noise and instrument model.
  Reference PR: [#14](https://github.com/paganol/fossil_sim/pull/14).

## 0.1.6 (2026-10-24)

- Correct the Static-sky Doppler slope using nearby PySM frequencies rather
  than finite differences on the coarse simulation grid. Apply orbital velocity
  to Galactic emission and solar plus orbital velocity to CMB anisotropies and
  extragalactic emission. **Config change:** `static_sky.orbital_doppler` becomes
  `static_sky.doppler`; the Static-sky cache is rebuilt automatically.
  Reference PR: [#11](https://github.com/paganol/fossil_sim/pull/11).
- Correct spectral white-noise normalization to match the interferogram model:
  the noise standard deviation increases by a factor of sqrt(2). Add the
  provisional `noise.modulation_efficiency` parameter (default 1), with a
  Monte Carlo validation from sample NEP to recovered spectra.
  Reference PR: [#12](https://github.com/paganol/fossil_sim/pull/12).
- Bump the package version for the 0.1.6 release.
  Reference PR: [#13](https://github.com/paganol/fossil_sim/pull/13).

## 0.1.5 (2026-10-24)

First tagged release.

- Define Mode switches by observing-time blocks, a repeating Mode pattern and
  dead time that can cross Scan boundaries; the clock pauses during Recycling.
  **Config change:** `operations.mode_switch_interval_s`, `mode_pattern` and
  `mode_switch_dead_time_s` replace `ring_mode_pattern` and
  `mode_switch_dead_ifgs`. The baseline observes every Ring in both B10 and B01.
  Reference PR: [#10](https://github.com/paganol/fossil_sim/pull/10).
- Select science IFGs by Mode when making spectral maps, validate unsupported
  configurations and mission lengths, and record dirty source checkouts in
  output provenance.
  Reference PR: [#9](https://github.com/paganol/fossil_sim/pull/9).
- Add a spectral-maps example comparing noisy and noiseless six-month runs,
  with FITS outputs and a report validating sky coverage and map variance.
  Reference PR: [#8](https://github.com/paganol/fossil_sim/pull/8).
- Execute and commit example-notebook outputs automatically on `main`, and
  validate notebooks using their generator in CI.
  Reference PR: [#7](https://github.com/paganol/fossil_sim/pull/7).

## 0.1.4.dev0

- Add spectral-domain white noise per IFG, frequency bin and detector unit,
  reproducible component seeds, Recycling noise scaling and sky-variance maps.
  Adopt the M8 30–2000 GHz band and 200 GHz dichroic split. The original noise
  normalization is corrected in 0.1.6.
  Reference PR: [#6](https://github.com/paganol/fossil_sim/pull/6).

## 0.1.3.dev0

- Add zodiacal light through zodipy using the spacecraft heliocentric position.
  Evaluate each Ring at its start and end and interpolate per IFG; add model
  configuration, validation and diagnostic maps.
  Reference PR: [#5](https://github.com/paganol/fossil_sim/pull/5).

## 0.1.2.dev0

- Add a cached Static-sky harmonic cube from PySM, beam smoothing and coordinate
  rotation, exact ducc0 or nearest-pixel sampling, and first-order orbital
  Doppler. Build the cache once on rank 0 for MPI simulations. The Doppler
  treatment is revised in 0.1.6.
  Reference PR: [#4](https://github.com/paganol/fossil_sim/pull/4).

## 0.1.1.dev0

- Add Isotropic backgrounds with the exact solar and orbital Boost,
  Mode-dependent port signals, differential spectra, streaming spectral maps
  and optional HDF5 TOD with provenance.
  Reference PR: [#2](https://github.com/paganol/fossil_sim/pull/2).
- Introduce Sphinx/MyST documentation, Read the Docs configuration, concept and
  mathematics pages, tutorials, API reference and a strict documentation CI job.
  Reference PR: [#3](https://github.com/paganol/fossil_sim/pull/3).

## 0.1.0.dev0

- Introduce the simulation package, baseline configuration, approximate L2 and
  Lissajous orbit, observation timeline, scanning strategy, streaming hit maps,
  ideal-FTS interferograms and optional MPI driver. Add examples, notebooks,
  tests, CI, domain terminology and the first design decision.
  Reference PR: [#1](https://github.com/paganol/fossil_sim/pull/1).
