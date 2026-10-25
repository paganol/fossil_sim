# Noise

Detector noise white in the IFG becomes, after the Fourier transform, a noise
flat in power across the synthesized bins. An IFG of duration $t$ has $N$
samples of $\delta t = t / N$, each with noise $\mathrm{NEP}/\sqrt{2\,\delta t}$,
and the DU power is ({mod}`fossil_sim.ifg`)

$$
P(x) = A\Omega\,\alpha\,f_\mathrm{opt}\,m \sum_\nu \left[I_1 + I_2 + s\,(I_1 - I_2)\cos(2\pi\nu x / c)\right]\Delta\nu .
$$

The least-squares amplitude of each cosine has variance
$2\sigma^2/N = \mathrm{NEP}^2/t$, so each bin of width $\Delta\nu$ gets
independent Gaussian noise of standard deviation

$$
\delta S_\nu = \frac{\mathrm{NEP}}{\sqrt{t}\, A\Omega\, \Delta\nu\, \alpha\, f_\mathrm{opt}\, m},
$$

with $A\Omega$ the étendue, $\alpha$ the detector absorptivity and
$f_\mathrm{opt}$ the optical efficiency (no apodization). A Monte Carlo test
(NEP → IFG samples → least-squares spectrum) checks this to 3 %. PIXIE
{cite}`kogut2011` writes $\sqrt{2t}$ instead of $\sqrt{t}$: its noise is the
modulus of the complex Fourier coefficient, which for a cosine of amplitude
$a$ is $a/2$.

$m$ (`modulation_efficiency`) is a provisional parameter, to be removed once
the convention is fixed with the instrument team: the amplitude of the
modulated difference in the DU power. It is 1 when one DU sees the full
difference (the convention used so far), $1/4$ for a single-polarization
bolometer of a Martin-Puplett FTS (PIXIE, Kogut et al., arXiv:2405.20403v2, eq. 1), $1/2$ for a DU summing
two of them; the NEP must refer to the same detector.

The noise is added to the differential spectrum $d(\nu)$ of each DU, so that
it propagates to the sky estimate (with `ifg.simulate = true` it is drawn on the
IFG samples instead, $\mathrm{NEP}/\sqrt{2\,\delta t}$, and reaches $d(\nu)$
through the fit, {doc}`interferogram`); IFGs of the Recycling window get
`recycling_noise_factor` times this level. The NEPs are the totals of the M8
proposal (Table 2.3, photon-noise dominated), read as one value per DU (four
pixels combined); AΩ, $\alpha$ and $f_\mathrm{opt}$ are PIXIE placeholders.

With the baseline values ($m = 1$) $\delta S_\nu \simeq 1.1\times10^{4}$ Jy/sr
(LF) and $3.0\times10^{5}$ Jy/sr (HF) per IFG and DU, i.e. a mission monopole
sensitivity of $\simeq 1.2$ Jy/sr below 200 GHz. This is above the
0.64–0.8 Jy/sr requirement of the proposal, which the earlier $\sqrt{2t}$
formula matched (0.83 Jy/sr): the proposal numbers likely follow the PIXIE
convention, to be checked together with $m$ and the NEP definition.

Each realization is drawn from a `numpy.random.SeedSequence` keyed by
(master seed, component, cycle, DU), never by MPI rank (ADR 0001): a run gives
the same noise for any number of ranks.
