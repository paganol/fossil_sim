# Validation

Every module has unit tests (`test/`) and doctests; CI runs them on Python
3.12 together with the notebooks. The checks below are those that validate the
physics and the mission design.

## Timeline and scanning

- Scan of 1440 s and 720 IFGs, 0.5° per IFG, 3.2 samples per 1.6° beam; the
  23 Scans of a Ring retrace the same circle.
- Ring steps equal to one third of the actual solar motion over each cycle
  (0.46–0.49°, about 1.44° per cycle), near perihelion and aphelion.
- Maximum Solar Aspect Angle below 0.36° (`equal_per_cycle`) and 0.72°
  (`anchor_ring_start`) over the mission; the `fixed` rule drifts above 10°
  within a year.
- Full sky at nside 64 within 183 days, for B10 and B01 separately
  (half-Ring Mode blocks).
- Time lost: 20.69 % in Recycling windows, 0.72 % in mode switches
  (half-Ring Mode blocks, 300 s dead time).

## Orbit and frames

- Analytic velocity equal to the numerical derivative of the position to
  1 mm s⁻¹; L2 distance and speed within physical ranges.
- ICRS → ecliptic rotation equal to healpy's `C → E` rotation; solar dipole
  direction equal to healpy's `G → E` rotation.

## Sky and signal

- Boosted blackbody equal to $B_\nu(T_0 / \gamma(1 - \boldsymbol\beta\cdot\hat{\mathbf n}))$ to $10^{-12}$.
- Rest-frame spectra equal to bisou_sky to $10^{-6}$ (optional test).
- Mode table and output signs; exact sky reconstruction in B10 and B01.
- Kinematic dipole recovered from the spectral maps within 2 %.
- Static sky (synthetic emission): Galactic monopole kept, extragalactic
  monopole removed, Galactic → ecliptic rotation and beam $b_1$ exact to
  $10^{-6}$; cache reused and keyed by content; nearest-pixel sampling
  consistent with the exact alm sampling; first-order Doppler
  $(3 - \alpha)\,\beta\cdot\hat{\mathbf n}\,I$ for power laws.
- Zodiacal light: Ring-by-Ring interpolation equal to the direct zodipy
  evaluation to $2\times10^{-5}$ (mid-Ring and Recycling IFGs); brighter on
  the ecliptic than at the pole, rising with frequency.
- Static sky with real pysm3 models: dust dominates at THz on the Galactic
  plane (optional test).
- Ideal-FTS interferogram round trip exact; L + R constant.
- Noise: per-IFG level equal to the PIXIE formula; realizations reproducible
  and keyed by cycle; Recycling IFGs × `recycling_noise_factor`; (noisy −
  noiseless) sky map normalized by the variance map has unit standard
  deviation to 2 %.

## Spectral maps over six months

`examples/spectral_maps.py` with the baseline (127 Cooler cycles, 184 days,
Static sky `c1 d2 s1 a1 f1` + zodi + noise; about 25 min per run, serial):

- 7.96 M IFGs, 6.30 M valid (20.80 % lost); maximum Solar Aspect Angle 0.35°.
- Sky coverage 100 % at nside 64 (median 92 hits per pixel), 86.5 % at nside
  128 (median 23): the Ring step (~0.48°) is close to the nside 128 pixel size,
  so full coverage at nside 128 needs the whole year.
- (noisy − noiseless) / √variance has standard deviation between 0.991 and
  1.005 in every frequency bin at nside 64 (0.995–1.004 at nside 128): the
  noise integrates down as 1/√hits over the mission.

## Cross-checks

- Pointing equal to litebird_sim's quaternion composition to better than
  1 nrad, with and without Precession (optional test).
- MPI run with 3 ranks identical to the serial run (`scripts/mpi_check.py`).
