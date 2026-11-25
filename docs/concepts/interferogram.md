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

The interferogram and recovery functions accept an optional `freq_step_ghz`
bin width. It defaults to the spacing between the first two grid frequencies;
for a single bin it must be supplied explicitly. The simulation always passes
`instrument.freq_step_ghz`, so single-bin configurations are supported in both
the spectral and IFG domains.

## Generation grid and band edges

The fit above is exact only because, by default, the interferogram is built
on the same frequencies as the fit, one value per bin: a discrete input on the
fit basis. The sky is continuous, and an instrument with resolution
$c / 2x_\mathrm{max}$ sees each bin through its line shape, with leakage
between bins and from outside the band. With
`ifg.generation_freq_step_ghz` and `ifg.generation_freq_max_ghz` the
interferogram is generated on a finer, wider grid
({func}`fossil_sim.instrument.generation_frequency_grid_ghz`): the multiples
$k\,\delta\nu$ up to the maximum where some DU transmits, each point standing
for a band $\delta\nu$ wide, so $\Delta\nu$ in $P(x)$ becomes $\delta\nu$.
The fit stays on the instrument bins, whose centres are among the generation
frequencies. Both keys are off in the baseline (the generation grid is the
instrument grid); a typical choice is $\delta\nu$ = 3 GHz up to 3825 GHz,
below the Nyquist frequency of the OPD sampling ($c / 2\Delta x$ = 3830 GHz
at baseline), checked by the configuration.

The sky entering the instrument (Isotropic backgrounds, Static sky, zodiacal
light, BBIR) is evaluated on the generation grid, so the Static-sky cache is
built on it. The cache key includes the frequencies: a new generation grid
(or switching `ifg.simulate` on with one set) needs a new cache, about
$3 \times 16$ bytes per frequency and $a_{\ell m}$ ($I$ and the two Doppler
$K$ cubes), 2.4 GB for 690 frequencies at $\ell_\mathrm{max}$ = 383 against
0.5 GB for the 132 bins. Build it once before parallel runs with
`scripts/build_static_sky.py config.toml`, where pysm3 can download its
templates (on clusters, a node with network access).

Cost, baseline with `sky_variation = true` on CINECA Galileo100, one thread
per MPI rank: the 690-frequency cache takes 13 min to build (3 min for the 132
bins); a Cooler cycle takes about 60 s, with a peak of 5.2 GB per rank
(2.6–4.2 GB on the instrument grid). A full-mission run (882 cycles) on 48
ranks then takes about 20 min. Each sample interpolates the sky from four
window points, so the modulated term is computed only on those, by groups of
samples, and once per IFG when the window is constant (sky of the IFG centre,
or both Ports on the BBIR); the rest of a cycle is mostly the Static-sky
synthesis on the Ring circles and array operations on the sky windows. The DU response maps each generation frequency to its nearest
bin: the gain of the bin, and the offset $\varepsilon B_\nu(T_\mathrm{inst})$
at the generation frequency.

The band mask $w$ becomes a transmission with cosine edges
({meth}`fossil_sim.instrument.DetectorUnit.transmission`): at each edge $e$
(band limits and dichroic split) it varies as $\tfrac12(1 + \sin\pi t)$,
$t = (\nu - e) / (f e)$, over $|t| \le \tfrac12$, with $f$ =
`instrument.band_edge_fraction` (0 = top-hat, the baseline; 0.1 is a
placeholder until the filter specifications). The two DUs of an output share
the dichroic edge and sum to 1. Band edges and the generation grid act only in
the IFG domain; the spectral domain keeps one value per bin and top-hat bands.

The size of these effects, measured on the Isotropic backgrounds alone (no
Static sky, no zodiacal light, no noise), in science IFGs of the first cycle,
generation grid of 3 GHz up to 3825 GHz (`test/test_generation_grid.py`),
relative to the largest value of each bin:

- top-hat bands: the recovered bins differ from the bin values of the spectral
  domain by the average of the spectrum over each bin, $10^{-5}$–$10^{-3}$ up
  to 1.5 THz and about $10^{-2}$ on the Wien tail of the CMB near 1.9 THz. The
  bin at 1995 GHz, cut by the top-hat at 2000 GHz, loses about 20 % of its
  power;
- cosine edges ($f$ = 0.1) instead of top-hat: up to 3 % on the bins next to
  the 2000 GHz edge, and on every other bin, through the unapodized least
  squares, $2\times10^{-5}$–$4\times10^{-3}$, decreasing away from the edges
  ($5\times10^{-4}$ MJy/sr at 1 THz).

These are effects of the instrument model, not of the sky: the Static sky
adds structure (dust, synchrotron) but no new mechanism.

## IFG-domain simulation

With `ifg.simulate = true` each IFG of each DU is built sample by sample, white
noise $\mathrm{NEP}/\sqrt{2\,\delta t}$ is added to the samples
({doc}`noise`) and the differential spectrum is recovered by the least squares
above; the rest of the chain (sky estimate, maps, flags) is unchanged. The
samples are equispaced in time, sample $j$ at $(j + \tfrac12)\,\delta t$, with
the OPD linear in time between the first and last sample, at $-x_\mathrm{max}$
and $+x_\mathrm{max}$; the mirror turns around in the half sample at each end
of the IFG (a placeholder mirror trajectory).

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
