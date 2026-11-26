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

### Apodization

Without apodization each bin sees the spectrum through the line shape of the
finite scan, a sinc of FWHM $1.21\,c/2x_\mathrm{max}$ (18 GHz at baseline)
with a first sidelobe of −22 %. On the fit basis this is invisible (the
recovery is exact), but structure off the basis rings over the whole band:
a sharp band edge, an unresolved line, the sky changing during the IFG.
`ifg.apodization` weights the samples by $A(x)$, $A(0) = 1$, before the fit
({func}`fossil_sim.ifg.differential_spectrum_solver`): the constant $a_0$ of
the plain fit is removed first, then the apodized samples are projected on
the cosines,

$$
R_A = D^+_{1:}\,\mathrm{diag}(A)\,\bigl(\mathbb 1 - \mathbf 1\,D^+_{0}\bigr),
$$

with $D$ the design matrix of the plain fit. The rows still annihilate a
constant, so the 1/f factor ({doc}`noise`) applies unchanged. The choices
are the Norton–Beer weak, medium and strong functions {cite}`norton1976`,
$A = \sum_i c_i (1 - u^2)^i$, and Hann, $\cos^2(\pi u/2)$, with
$u = |x|/x_\mathrm{max}$. On the baseline (bin at 1005 GHz):

| `apodization` | FWHM | first sidelobe | σ per bin | ρ(k, k+1) | ρ(k, k+2) |
|---|---|---|---|---|---|
| `none` | 18.0 GHz | −22 % | 1 | 0.00 | 0.00 |
| `norton_beer_weak` | 21.7 GHz | 5.9 % | 0.74 | 0.40 | 0.03 |
| `norton_beer_medium` | 25.3 GHz | 1.4 % | 0.66 | 0.57 | 0.08 |
| `norton_beer_strong` | 28.9 GHz | 0.4 % | 0.61 | 0.67 | 0.18 |
| `hann` | 29.9 GHz | 2.7 % | 0.61 | 0.67 | 0.17 |

The noise per bin falls only because the bins become correlated: at equal
resolution apodization costs sensitivity. The variance maps carry the
diagonal of the covariance only. $A(0) = 1$ keeps the level of a smooth
spectrum (within $10^{-3}$ away from the edges), but the line shape of the
last bins of each DU reaches past the fitted band: for a flat spectrum the
outermost bin is recovered 15–25 % low, and the nominal calibration does
not correct it.

Effect on the systematics of the generation grid (3 GHz grid, 600 science
IFGs at $|b| < 2°$ in the first Cooler cycle, no noise; sky estimate, the
mean of L and R):

- **Sky change during the IFG, top-hat bands.** The modulated part falls
  10–40 times inside the HF band with Norton–Beer medium or Hann (Hann:
  $7\times10^{-4} \to 2\times10^{-5}$ MJy/sr at 1500 GHz) and 3–5 times in
  the 1995 GHz bin (0.15 → 0.03–0.05 MJy/sr); the 1980 GHz bin stays at
  0.05 MJy/sr, its wider line shape now catching the 2000 GHz step.
- **Sky change, cosine edges** (`band_edge_fraction = 0.1`). Already small
  ($\le 4\times10^{-3}$ MJy/sr), and apodization leaves it unchanged except
  next to the dichroic (195–210 GHz, about 3 times lower).
- **Unmodulated leakage in one DU.** The even part of the sky change grows
  as $t^2$, largest at the ends of the scan, where $A$ is small: the leakage
  falls 3.5 times (Norton–Beer medium) to 6 times (Hann) with either edge;
  it cancels in the sky estimate anyway.
- **Static sky, top-hat bands.** Without apodization the 2000 GHz step rings
  hundreds of GHz into the band (0.14 MJy/sr at 1500 GHz, 0.6 at 1800, 3.2 at
  1950, rms of the error on the bin-centre sky); Norton–Beer medium brings it
  to 0.08, 0.14 and 0.24 MJy/sr. In mid-band, away from any edge, the wider
  line shape instead raises the error from 0.01–0.05 to about 0.08 MJy/sr,
  a local smoothing of the dust curvature that a model can include.
- **Static sky, cosine edges.** Apodization changes little; the error is
  dominated by an offset of about +0.1 MJy/sr, nearly constant over the HF
  bins and the same with any apodization (not understood yet).

Apodization is therefore useful against sharp spectral structure (a sharp
filter edge, lines); with smooth filter edges its gain is small, and it
costs resolution, bin correlations and the edge bins. The default is
`"none"`.

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
or both Ports on the BBIR).

When the Static sky has its low-rank factorization ({doc}`static_sky`,
`sampling = "alm"`), the IFGs avoid the sky windows on the generation grid
altogether. The model is linear in the Port spectra, so each Port that sees
the sky is split into the rest of the sky (Isotropic backgrounds, zodiacal
light, at the IFG centre) and the Static sky $\sum_k a_k U_k(\nu)$, with the
amplitudes $a_k$ on the window (Doppler parts scaled by
$\boldsymbol\beta\cdot\hat{\mathbf n}$). The IFG of each basis spectrum
$U_k$ is computed once per cycle and DU
({func}`fossil_sim.ifg.static_basis_interferograms`), and the Static-sky term
of every IFG is their combination with the interpolated amplitudes
({func}`fossil_sim.ifg.du_interferograms_low_rank`): a cost per sample
proportional to $K \approx 87$ instead of the about 640 band frequencies, with
no approximation beyond the factorization. Measured with cosine edges
(`band_edge_fraction = 0.1`) and noise on the 690-frequency grid: a cycle on a
single rank went from 63 s to 44 s with the low-rank synthesis, then to 17 s
with the basis IFGs; the full-mission run on 48 ranks from 32 min to 24 min,
then to 10 min (peak 7.9 → 5.8 GB per rank), maps equal to
$10^{-5}\sigma$ in every pixel. Without the factorization (`rank_epsilon = 0`,
or nearest-pixel sampling) the sky windows are built on the generation grid
as before.

The DU response maps each generation frequency to its nearest
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
sky at the centre. The even part does not cancel, and it enters through both
terms of the interferogram.

**Unmodulated term.** $\sum_k w_k (I_1 + I_2)\,\Delta\nu$, the sky integrated
over the whole band, varies during the IFG and leaks into the cosine
amplitudes of every bin. In the spectrum of one DU this is the larger part (in
the Galactic plane, about 0.1 MJy/sr in every HF bin of the baseline). A
beam-size source crossed at the IFG centre changes the recovered spectrum of a
DU by a few per cent of its amplitude, most on the bins where the source is
faint compared with the band total. The leakage has the same sign in the L
and R outputs, while the sky estimate multiplies the spectrum by $s$, so it
cancels in the mean of the two outputs and does not reach the maps.

**Modulated term.** A sky changing during the scan changes the amplitude of
the fringes along the OPD; in frequency this is a convolution of the
spectrum with a narrow kernel. Inside a band where the spectrum is smooth it
moves little power, but a sharp band edge is a step in the measured spectrum,
which the convolution smears asymmetrically onto the bins next to it. This
part survives in the sky estimate. It is largest at the top of the HF band,
where the band ends on the bright dust spectrum. In the Galactic plane, for
the baseline with top-hat bands, it grows by about 4 times per bin towards
the 2000 GHz edge: 0.0014 MJy/sr at 1500 GHz, 0.015 at 1950 GHz, 0.31 at
1995 GHz. In the full-mission maps with photon noise (nside 64) it reaches
650σ per pixel in the 1995 GHz bin. It is not numerical: generating each
sample from the exact sky in its direction gives the same values as the
interpolated windows.

Smooth band edges remove most of it. On the 3 GHz generation grid, the same
IFGs give 0.20 MJy/sr at 1995 GHz with top-hat bands and 0.0003 MJy/sr with
cosine edges (`band_edge_fraction = 0.1`); with cosine edges only about
0.005 MJy/sr remains, spread over the transition (1900–2100 GHz), and the
bins inside the band improve 10–60 times. The full-mission maps confirm it
(baseline with photon noise, nside 64, rms over the sky of the difference
with and without sky change, in units of the noise per pixel): with cosine
edges the HF bins away from the edges go from 3σ to 0.05–0.1σ and the
1995 GHz bin from 650σ to 1.3σ, while the bins inside the transitions keep
a smaller effect (up to 8σ at 1905–1950 GHz, 1.9σ at 195 GHz, dominated by
the Galactic plane); the noise is unchanged. Studies of the sky change during
the IFG therefore need the generation grid with realistic band edges; the
baseline (top-hat on the instrument grid) is the worst case. For the
instrument, the shape of the filter at the top of the band matters: a sharp
cut-off turns the Galactic sky change during the IFG into a systematic error
in the highest channels.
