# Interferogram

An ideal-FTS interferogram model ({mod}`fossil_sim.ifg`). By default the
simulation works in the spectral domain; with `ifg.simulate = true` it builds
the samples of every interferogram and recovers the spectrum from them
(below). Instrumental systematics (mirror trajectory, time constants, phase
errors, arm imbalance) are planned on top of the IFG-domain simulation.

## Model

On a frequency grid $\nu_k$ with spacing $\Delta\nu$, the signal of a DU with
output sign $s$ and band mask $w_k$ at optical path difference $x$ is

$$
P(x) = \sum_k w_k\left[I_1(\nu_k) + I_2(\nu_k) + s\,\big(I_1(\nu_k) - I_2(\nu_k)\big)\cos\frac{2\pi\nu_k x}{c}\right]\Delta\nu,
$$

in MJy sr⁻¹ GHz (optical efficiency and étendue not included). The mirror
position is $\delta = x / \kappa$, with $\kappa$ = `opd_per_mirror_position`;
the form $\cos(4\pi\nu\delta/c)$ of the mission document corresponds to
$\kappa = 2$ (a placeholder until the FTS design fixes it).

## Resolution and sampling

A symmetric scan $-x_\mathrm{max} \le x \le x_\mathrm{max}$ resolves
frequencies separated by

$$
\Delta\nu = \frac{c}{2\,x_\mathrm{max}},
$$

so 15 GHz bins need $x_\mathrm{max} = 10$ mm. Sampling the highest frequency
$\nu_\mathrm{max}$ without aliasing requires an OPD step

$$
\Delta x \le \frac{c}{2\,\nu_\mathrm{max}} = 75\ \mu\mathrm{m}\ \text{at 2 THz};
$$

the baseline uses 512 samples, a step of 39 µm.

## Spectrum recovery

The differential spectrum is recovered by linear least squares on the model

$$
P(x_j) = a_0 + \Delta\nu \sum_{k \in \mathrm{band}} a_k \cos\frac{2\pi\nu_k x_j}{c},
$$

where $a_0$ absorbs the unmodulated term and $a_k = s\,(I_1 - I_2)(\nu_k)$.
With the grid matched to the scan ($x_\mathrm{max} = c / 2\Delta\nu$) the
cosines are close to orthogonal and the recovery is exact for a noiseless,
band-limited input.

## IFG-domain simulation

With `ifg.simulate = true` each IFG of each DU is built sample by sample, white
noise $\mathrm{NEP}/\sqrt{2\,\delta t}$ is added to the samples
({doc}`noise`) and the differential spectrum is recovered by the least squares
above; the rest of the chain (sky estimate, maps, flags) is unchanged. The
samples are equispaced in time, sample $j$ at $(j + \tfrac12)\,\delta t$, with
the OPD linear in time from $-x_\mathrm{max}$ to $+x_\mathrm{max}$ (a
placeholder mirror trajectory).

### Sky during the IFG

In 2 s the boresight moves by 0.5° along the Ring circle, a third of the
beam. With `ifg.sky_variation = true` the Static sky follows the boresight
during the IFG. The Scans of a Ring retrace the same circle, so along a Ring
the Static sky is a function of the spin phase only, band-limited by the beam.
It is synthesized once per Ring on a grid of 8 points per IFG (0.0625°, 5760
per circle) for $I$ and the two Doppler $K$ cubes, and each sample takes the
cubic (Lagrange) interpolation of its four neighbours. The interpolation error
is below $10^{-4}$ of a beam-size source; the synthesis is cheaper than the
per-IFG sampling of the spectral domain, and a cycle costs about 17 s instead
of 15 s. The Isotropic backgrounds, the zodiacal light and $\boldsymbol\beta\cdot\hat{\mathbf n}$
in the Doppler term, smooth on 0.5°, keep their value at the IFG centre. With
`sky_variation = false` the whole sky is that of the IFG centre: without noise
the result then equals the spectral domain to numerical precision.

A precession within the Ring would break the circle symmetry; the general
alternative is the interpolator of ducc0 (`totalconvolve.Interpolator`), as in
litebird_sim.

### What the sky change does

A symmetric scan rejects any part of the signal odd in time about the IFG
centre: a sky changing linearly during the IFG gives the same spectrum as the
sky at the centre. The even part does not cancel, and it enters mostly through
the unmodulated term: $\sum_k w_k (I_1 + I_2)\,\Delta\nu$, the sky integrated
over the whole band, varies during the IFG and leaks into the cosine
amplitudes of every bin. A beam-size source crossed at the IFG centre changes
the recovered spectrum by a few per cent of its amplitude, and the effect is
largest on the bins where the source is faint compared with the band total.
