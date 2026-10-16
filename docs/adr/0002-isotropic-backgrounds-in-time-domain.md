# Isotropic backgrounds are boosted per IFG, never stored in the static sky

All backgrounds that are isotropic in the CMB rest frame (CMB blackbody, y with relativistic corrections, μ, CIB monopole, extragalactic CO) are evaluated per IFG with the exact Lorentz boost I_obs(ν, n̂) = (ν/ν′)³ I_rest(ν′), ν′ = ν γ (1 − β·n̂), β = β_sun + β_orbit(t). The static sky (pysm3 anisotropies, stored as beam-smoothed alm) contains no monopole or dipole. One formula thus replaces separate CMB/CIB dipole maps, gives the kinematic dipoles of every distortion for free, and keeps the orbital dipole exact in time — needed if it becomes a primary calibrator. The beam acts only on the anisotropic part of the boosted signal, through b₁ = exp(−σ²) (higher multipoles neglected).

## Considered Options

- Monopole + solar dipole baked into the static cube (as bisou_sky `get_sky`): rejected, the orbital dipole is time dependent and the cube would duplicate it.
