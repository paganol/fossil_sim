# Instrument and signal

## Ports, outputs and detector units

The FTS has two inputs, Port 1 and Port 2, fed by BSM1 and BSM2. Port $p$ sees
either the sky through Telescope $p$ or the BBIR. Writing $I_1(\nu)$ and
$I_2(\nu)$ for the spectra entering the Ports, the two outputs measure

$$
P_L(x) \propto \int \left[I_1 + I_2 + (I_1 - I_2)\cos\frac{2\pi\nu x}{c}\right] w(\nu)\,d\nu,
$$

$$
P_R(x) \propto \int \left[I_1 + I_2 - (I_1 - I_2)\cos\frac{2\pi\nu x}{c}\right] w(\nu)\,d\nu,
$$

with $x$ the optical path difference and $w(\nu)$ the band. A dichroic at each
output splits the band at 200 GHz (30–2000 GHz band,
M8 proposal) into a low-frequency and a high-frequency
detector unit: LF-L, HF-L, LF-R, HF-R. In the spectral domain each DU measures
the differential spectrum

$$
d(\nu) = s\,\left[I_1(\nu) - I_2(\nu)\right], \qquad s = +1\ (\mathrm{L}),\ -1\ (\mathrm{R}),
$$

on the frequency bins of its band. Spectra are evaluated at the bin centres
(no channel response or instrumental line shape in the spectral domain).

## Modes

A Mode is written $\mathrm{B}xy$, with $x$ the state of BSM1 and $y$ the
state of BSM2 (1 = sky, 0 = BBIR). With $S$ the sky and $B_p$ the BBIR
spectrum seen by Port $p$, the L output measures

| Mode | Port 1 | Port 2 | $d_L = I_1 - I_2$ |
|---|---|---|---|
| B00 | BBIR | BBIR | $B_1 - B_2$ |
| B10 | sky | BBIR | $+(S - B_2)$ |
| B01 | BBIR | sky | $-(S - B_1)$ |
| B11 | sky | sky | $S_1 - S_2$ |

and the R output the opposite. B10 and B01 are the science Modes: they see
the sky with opposite signs, so that their difference doubles the sky signal
and cancels any common instrumental offset. B00 cross-calibrates the two
Ports; B11 is a null test (zero for identical, co-aligned telescopes).

## BBIR

The BBIR spectrum seen by Port $p$ is an ideal blackbody,
$B_p(\nu) = B_\nu(T + \delta T_p)$, with $T$ the main-absorber set-point and
$\delta T_p$ a per-Port offset; the data model already allows a different
effective spectrum per Port. The high-frequency cavity is off in the baseline.

## Sky reconstruction

Assuming perfect calibration (unit gain, BBIR temperatures known exactly), the
sky spectrum is recovered from each DU in a science Mode as

$$
\hat S(\nu) = B_q(\nu) + \sigma\, s\, d(\nu), \qquad
\sigma = \begin{cases} +1 & \mathrm{B10},\ q = 2 \\ -1 & \mathrm{B01},\ q = 1 \end{cases}
$$

with $q$ the Port looking at the BBIR. In each frequency bin the estimates of
the L and R DUs covering it are averaged. B00 and B11 IFGs give no sky
estimate.

## White noise

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
it propagates to the sky estimate; IFGs of the Recycling window get
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
