# The Static sky is defined in the barycentric frame

pysm3 templates derive from Planck/WMAP maps observed from L2 and averaged over years, i.e. effectively in the Solar-System barycentric frame: the CMB solar dipole is removed, but foreground anisotropies are not corrected for the solar motion. fossil_sim therefore treats the Static sky as given in the barycentric frame and applies only the orbital velocity of the spacecraft (β_orb ≈ 10⁻⁴) to it, as a first-order Doppler term ΔI = (β_orb·n̂)(3 − d ln I/d ln ν) I, with the spectral index from finite differences on the frequency cube. This is the orbital dipole of the foregrounds (up to ~0.1 MJy/sr on the Galactic plane at THz frequencies, far above y), needed for orbital-dipole calibration studies. Aberration (~20″) and O(β²) terms are neglected.

## Considered Options

- No boost of the Static sky: rejected, it misses the foreground orbital dipole.
- Exact boost: requires spectral interpolation of the cube off the frequency grid; its gain, O(β²) ~ 10⁻⁸ I, does not justify the interpolation error.
