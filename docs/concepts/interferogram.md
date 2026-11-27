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
$A = \sum_i c_i (1 - u^2)^i$, Hann, $\cos^2(\pi u/2)$, triangular,
$1 - u$, and Kaiser, $I_0\bigl(\beta\sqrt{1 - u^2}\bigr)/I_0(\beta)$ with
$\beta$ = `ifg.kaiser_beta` (default 6), with $u = |x|/x_\mathrm{max}$. On the baseline (bin at 1005 GHz):

| `apodization` | FWHM | first sidelobe | σ per bin | ρ(k, k+1) | ρ(k, k+2) |
|---|---|---|---|---|---|
| `none` | 18.0 GHz | −22 % | 1 | 0.00 | 0.00 |
| `norton_beer_weak` | 21.7 GHz | 5.9 % | 0.74 | 0.40 | 0.03 |
| `norton_beer_medium` | 25.3 GHz | 1.4 % | 0.66 | 0.57 | 0.08 |
| `norton_beer_strong` | 28.9 GHz | 0.4 % | 0.61 | 0.67 | 0.18 |
| `hann` | 29.9 GHz | 2.7 % | 0.61 | 0.67 | 0.17 |
| `triangular` | 26.4 GHz | +4.7 % | 0.58 | 0.61 | 0.15 |
| `kaiser`, β = 4 | 24.9 GHz | 3.2 % | 0.68 | 0.54 | 0.06 |
| `kaiser`, β = 6 | 29.2 GHz | 0.65 % | 0.61 | 0.67 | 0.18 |
| `kaiser`, β = 8 | 33.2 GHz | 0.12 % | 0.56 | 0.74 | 0.29 |

The triangular line shape is $\mathrm{sinc}^2$, the classical FWHM of
$1.772\,c/(2x_\mathrm{max})$ (26.6 GHz; 26.4 GHz on the 512 samples): its
sidelobes are never negative, but the first is higher than with Norton–Beer
medium at a similar width, and they decay only as $1/\nu^2$. Their area
beyond the edges of the fitted band is missing, so a flat spectrum comes out
about 1 % low in mid-band and 5 % low a few bins from the edges.

The Kaiser window is a one-parameter family: $\beta = 0$ is no apodization,
and a larger $\beta$ widens the line shape, lowers the sidelobes and the
noise per bin, and correlates more bins. It makes the trade-off continuous,
for scanning it; at a given width Norton–Beer, optimized for this line
shape, has the lower sidelobes ($\beta = 4$ against medium, $\beta = 6$
against strong).

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

**Unmodulated term.** $\sum_k f_k (I_1 + I_2)\,\Delta\nu$, the sky integrated
over the whole band, varies during the IFG and leaks into the cosine
amplitudes of every bin. In the spectrum of one DU this is the larger part (in
the Galactic plane, about 0.1 MJy/sr in every HF bin of the baseline). A
beam-size source crossed at the IFG centre changes the recovered spectrum of a
DU by a few per cent of its amplitude, most on the bins where the source is
faint compared with the band total. The leakage has the same sign in the L
and R outputs, while the sky estimate takes $s\,\tilde d$ (step 7), so it
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

## From samples to maps

This section follows one block of IFGs through the IFG-domain chain
({meth}`fossil_sim.simulation.SignalModel.compute`), from the sky to the
binned maps. The spectral domain (`ifg.simulate = false`) skips steps 3–5:
it writes $d = g\,[s(I_1 - I_2) + o] + n$ directly, with Gaussian noise per
bin, and shares steps 6–8.

### 0. Set-up, once per run

Everything that depends only on the configuration is built once in the
constructor of {class}`~fossil_sim.simulation.SignalModel`:

- the OPD samples $x_j = -x_\mathrm{max} + 2x_\mathrm{max}\,j/(N-1)$,
  $j = 0, \dots, N-1$ ({func}`fossil_sim.ifg.opd_samples_m`);
- two frequency grids: the generation grid $\nu_k$, on which the sky and the
  IFGs are computed, and the instrument grid $\nu_b$ (bins of width
  $\Delta\nu$), on which the spectrum is fitted. They coincide unless the
  generation grid is set (above);
- for each DU, the modulation matrix on the generation frequencies it
  transmits ({func}`fossil_sim.ifg.modulation_matrix`),

  $$
  M_{jk} = \cos\frac{2\pi\nu_k x_j}{c};
  $$

- the interpolation weights $W_{jw}$ from the sky window to the samples
  ({func}`fossil_sim.ifg.sample_weights`), below;
- the least-squares solver $R$ of each DU (step 5);
- the true DU response on the generation grid
  ({func}`fossil_sim.response.generation_response`): gain $g_k$ and offset
  $o_k$, and the band factor $f_k = T_k\,g_k$, with $T_k$ the transmission.

### 1. Sky on the window of each IFG

Each IFG $i$ of the table ({func}`fossil_sim.simulation.build_ifg_table`)
carries the boresight at its centre, its Mode, the BBIR Port temperatures and
its flags. The Static sky is taken on a window of $G + 3 = 11$ points of the
Ring grid ($G = 8$ points per IFG, plus one before and two after),
$S_{iw}(\nu_k)$, $w = 0, \dots, G + 2$
({meth}`~fossil_sim.simulation.SignalModel.sky_windows`); the other
components are those of the IFG centre. Sample $j$ sits at the window
coordinate $u_j = (j + \tfrac12)\,G/N$ and sees

$$
S_i(t_j, \nu_k) = \sum_w W_{jw}\, S_{iw}(\nu_k), \qquad \sum_w W_{jw} = 1,
$$

with $W_{jw}$ the cubic Lagrange weights of the four grid points around
$u_j$ (four non-zero entries per row, so a constant sky passes exactly).

### 2. Port spectra

Each Port sees the sky when its BSM points to the sky, otherwise the BBIR
blackbody at the Port temperature
({func}`fossil_sim.signal.port_spectra_mjy_sr`):

| Mode | $I_1$ | $I_2$ |
|------|-------|-------|
| B10  | $S$ | $B_\nu(T_2)$ |
| B01  | $B_\nu(T_1)$ | $S$ |
| B11  | $S$ | $S$ |
| B00  | $B_\nu(T_1)$ | $B_\nu(T_2)$ |

### 3. Interferogram

For a DU with output sign $s$ ($+1$ for L, $-1$ for R) the samples are
({func}`fossil_sim.ifg.du_interferograms`)

$$
P_{ij} = m\,\Delta\nu \sum_w W_{jw}\left[\sum_k f_k\,(I_1 + I_2)_{iwk}
+ s \sum_k M_{jk}\, f_k \big((I_1 - I_2)_{iwk} + s\,o_k\big)\right],
$$

with $\Delta\nu$ the spacing of the generation grid and $m$ the modulation
efficiency ({doc}`noise`). The unmodulated term $(I_1 + I_2)$ gives a level
that changes only with the sky; the modulated term is the band-weighted
difference projected on the cosines. The gain multiplies both terms, the
offset $o = \varepsilon B_\nu(T_\mathrm{inst})$ only the modulated one, so the
cosine amplitudes are $g\,[s(I_1 - I_2) + o]$.

The model is linear in the Port spectra, which the code exploits twice. When
the window is constant (sky of the IFG centre, or both Ports on the BBIR),
$\sum_w W_{jw} D_w = D \sum_w W_{jw}$ and one product per IFG suffices. With
the low-rank Static sky $S(\hat{\mathbf n}, \nu) = \sum_\kappa
a_\kappa(\hat{\mathbf n})\,U_\kappa(\nu)$, the basis IFGs

$$
B_{j\kappa} = \sum_k M_{jk}\, f_k\, U_\kappa(\nu_k), \qquad
b_\kappa = \sum_k f_k\, U_\kappa(\nu_k)
$$

are computed once per cycle ({func}`fossil_sim.ifg.static_basis_interferograms`)
and only the amplitudes $a_\kappa$ are interpolated on the window
({func}`fossil_sim.ifg.du_interferograms_low_rank`).

### 4. Sample noise

Independent Gaussian noise $n_{ij} \sim \mathcal N(0, \sigma_s^2)$ is added to
every sample ({func}`fossil_sim.noise.sample_noise_sigma`), with

$$
\sigma_s = \sigma_\nu\,\Delta\nu\,m\,\sqrt{N/2},
$$

$\sigma_\nu$ the per-bin noise of the spectral domain ({doc}`noise`): the
least-squares amplitude of a cosine sampled $N$ times has variance
$2\sigma_s^2/N$, so both domains give the same noise per bin.

### 5. Least-squares spectrum

On the band bins $\nu_b$ of the DU (instrument grid), the design matrix, its
pseudo-inverse and the solver are
({func}`fossil_sim.ifg.differential_spectrum_solver`)

$$
D = \Big[\,\mathbf 1 \;\Big|\; \cos\frac{2\pi\nu_b x_j}{c}\,\Big]_{N \times (1 + n_\mathrm{band})},
\qquad
D^+ = (D^\mathsf{T}D)^{-1}D^\mathsf{T},
\qquad
R = \frac{D^+_{1:}}{m\,\Delta\nu_\mathrm{inst}},
\qquad
d_i = R\,(P_i + n_i),
$$

where the first row of the pseudo-inverse, $D^+_0$ (the constant, which
absorbs the unmodulated term), is dropped. Only cosines enter: the scan is
symmetric and the FTS ideal, with no phase. With the generation grid
equal to the instrument grid the input lies on the fit basis and, noise
aside, $d = g\,[s(I_1 - I_2) + o]$ exactly; on a finer grid $d$ is a weighted
bin average, with leakage and out-of-band signal (above).

With `ifg.apodization` set (above), the constant of the plain fit is removed
first, then the samples are weighted by the window $A(x_j)$ and projected on
the cosines:

$$
R = \frac{D^+_{1:}\,\mathrm{diag}\big(A(x_j)\big)\,
\big(\mathbb 1_N - \mathbf 1\,D^+_0\big)}{m\,\Delta\nu_\mathrm{inst}}.
$$

The rows of $R$ still annihilate a constant, but each bin is now the spectrum
seen through the line shape of $A$, so $d = g\,[s(I_1 - I_2) + o]$ holds
only for a spectrum smooth on the scale of that line shape, and neighbouring
bins are correlated. In both cases $R$ depends only on the configuration and
is shared with the noise propagation (step 8).

### 6. Calibration

The calibration solution $(\hat g, \hat o)$ inverts the response
({meth}`fossil_sim.response.DUResponse.calibrate`):

$$
\tilde d = d/\hat g - \hat o \;\approx\; s\,(I_1 - I_2).
$$

### 7. Sky estimate per IFG

Only science Modes (one Port on the sky) give an estimate, with the BBIR
temperature taken as known ({func}`fossil_sim.signal.sky_estimate_mjy_sr`).
In B10 $\tilde d = s(S - B_2)$, in B01 $\tilde d = s(B_1 - S)$, so

$$
\hat S = B_\mathrm{other} + \sigma\, s\, \tilde d, \qquad
\sigma = \begin{cases} +1 & \text{B10} \\ -1 & \text{B01} \end{cases},
$$

with $B_\mathrm{other}$ the BBIR of the Port not on the sky. Each bin is the
mean over the DUs that cover it (L and R):

$$
\hat S_i(\nu_b) = \frac{1}{N_\mathrm{DU}} \sum_\mathrm{DU}
\big[B_\mathrm{other} + \sigma\, s\, \tilde d_\mathrm{DU}\big].
$$

### 8. Maps

Valid (`flags == 0`) science IFGs are binned with uniform weights in the
HEALPix pixel of the boresight at the IFG centre, $p_i$
({class}`fossil_sim.mapmaking.SpectralMaps`, {doc}`mapmaking`):

$$
m_p(\nu) = \frac{1}{N_p}\sum_{i \in p} \hat S_i(\nu).
$$

The white-noise variance map propagates the sample noise through the actual
solver ({func}`fossil_sim.noise.sky_estimate_noise_variance`):

$$
\sigma_{\hat S}^2(\nu_b) = \frac{1}{N_\mathrm{DU}^2}\sum_\mathrm{DU}
\frac{\sigma_s^2 \sum_j R_{bj}^2}{\hat g_b^2},
\qquad
\mathrm{Var}[m_p(\nu_b)] = \frac{\sigma_{\hat S}^2(\nu_b)}{N_p}.
$$

The raw maps average the uncalibrated $d$ per output (L, R) and Mode (B10,
B01), a diagnostic product with opposite sky sign in the two Modes.

Two consequences of this chain: the pixel is that of the IFG centre, while the
samples see the sky along 0.5° of scan, so the smearing along the scan is in
the data but not in the pointing; and there is no destriping or GLS map-maker,
the pixel mean being the maximum-likelihood map only for white noise without
drifts.
