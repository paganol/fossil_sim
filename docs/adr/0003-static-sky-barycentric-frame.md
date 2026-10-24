# The Static sky is split by rest frame: barycentric (Galactic) and CMB (cosmological)

pysm3 Galactic templates derive from Planck/WMAP maps observed from L2 and averaged over years, i.e. effectively in the Solar-System barycentric frame: the CMB solar dipole is removed, but foreground anisotropies are not corrected for the solar motion. CMB anisotropies (`c*`) and extragalactic models (`cib*`, `tsz*`, `ksz*`, `rg*`) are instead defined in the CMB rest frame. fossil_sim therefore keeps two parts of the Static sky and applies to each the first-order Doppler term ΔI = (β·n̂) K, K = 3I − dI/d ln ν, with the velocity of the observer in its frame:

- barycentric part: β = β_orb (≈ 10⁻⁴), the orbital dipole of the foregrounds (up to ~0.1 MJy/sr on the Galactic plane at THz frequencies, far above y), needed for orbital-dipole calibration studies;
- cosmological part: β = β_sun + β_orb (≈ 1.2 × 10⁻³), the Doppler modulation of the CMB anisotropies (tens of Jy/sr per pixel around the CMB peak).

K is computed when the cache is built, from pysm3 evaluated at ν e^{±ε} (ε = 10⁻², relative error ~10⁻⁴), and stored as two extra alm cubes. Aberration (~20″ orbital, ~4′ solar, small against the 1.6° beam) and O(β²) terms are neglected.

## Considered Options

- No boost of the Static sky: rejected, it misses the foreground orbital dipole.
- Exact boost: requires spectral interpolation of the cube off the frequency grid; its gain, O(β²) ~ 10⁻⁸ I (barycentric) or ~10⁻⁶ I (cosmological), does not justify the interpolation error.
- dI/d ln ν by finite differences on the frequency grid (v0.1.2–0.1.5): rejected, the 15 GHz grid is coarse in ln ν at low frequency (0.41 between 30 and 45 GHz) and one-sided at the edges: −8 instead of +100 Jy/sr at 30 GHz for I ∝ ν², errors up to ~30 Jy/sr on the Galactic plane for the baseline sky.
- A single Static sky in the barycentric frame (v0.1.2–0.1.5): rejected, it omits the solar Doppler modulation of the CMB anisotropies.
