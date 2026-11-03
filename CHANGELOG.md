# Changelog

Changes are grouped by version, newest first, with links to the pull requests
that introduced them. Versions before 0.1.5 were development versions, without
tagged releases.

## HEAD

Changes on `main` that are not in a release yet.

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
