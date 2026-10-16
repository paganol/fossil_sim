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
- Full sky at nside 64 within 183 days, for B10 and B01 separately.
- Time lost: 20.69 % in Recycling windows, 0.11 % in mode switches.

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
- Static sky with real pysm3 models: dust dominates at THz on the Galactic
  plane (optional test).
- Ideal-FTS interferogram round trip exact; L + R constant.

## Cross-checks

- Pointing equal to litebird_sim's quaternion composition to better than
  1 nrad, with and without Precession (optional test).
- MPI run with 3 ranks identical to the serial run (`scripts/mpi_check.py`).
