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

### Cross term in B11 (neglected)

In B11 the FTS coherently combines the fields of two separate apertures, so on
top of $I_1 - I_2$ the output carries a cross-correlation term between the two
telescopes. In a Martin-Puplett FTS it interferes orthogonal polarizations of
the two Ports, $\operatorname{Im}\langle E_{x1}E_{y2}^*\rangle\,\sin(2\pi\nu x/c)$:

- it needs sky polarization or differential cross-polarization between the two
  arms, and vanishes for an unpolarized sky and ideal grids;
- it is odd in $x$, so it lands in the quadrature (sine) spectrum and reaches
  the cosine spectrum only through phase errors;
- the monopole and dipole do not contribute for non-overlapping apertures:
  only sky structure on angular scales near $\lambda/b$, with $b$ the baseline,
  couples, through the edges of the beam.

With the baseline sky and placeholder optics it is $\sim 10^{-5}$ MJy/sr per
IFG, 3–4 orders of magnitude below the noise per IFG, with zero mean over the
sky. The simulation therefore sets $d = 0$ in B11 (identical telescopes). The
term is centred at $x = \pm\Delta$, with $\Delta$ the path imbalance between
the two telescope arms, and becomes relevant only with arm imbalance and
phase errors in the IFG domain.

## BBIR

The BBIR spectrum seen by Port $p$ is an ideal blackbody,
$B_p(\nu) = B_\nu(T + \delta T_p)$, with $T$ the main-absorber set-point and
$\delta T_p$ a per-Port offset; the data model already allows a different
effective spectrum per Port. The high-frequency cavity is off in the baseline.

## Detector-unit response

Each DU measures, in each bin of its band,

$$
d(\nu) = g(\nu)\,\bigl[s\,(I_1 - I_2)(\nu) + o(\nu)\bigr] + n(\nu),
$$

with $g$ the gain relative to the nominal response (dimensionless), $o$ an
offset referred to the sky and $n$ the detector noise ({doc}`noise`), which is
outside the gain. The gain is a scalar per DU times $1 + \delta g(\nu)$, with
$\delta g$ Gaussian of rms `response.gain_rms`, independent between bins and
DUs and fixed by the seed for the whole mission. The offset is the emission
imbalance between the two arms, $o(\nu) = \varepsilon\, B_\nu(T_\mathrm{inst})$,
with a signed emissivity $\varepsilon$ per DU: an imbalance that enters like
$I_1 - I_2$ has $\varepsilon_R = -\varepsilon_L$. The baseline is the ideal
instrument, $g = 1$ and $o = 0$; gain drifts in time are not modelled.

In the IFG domain the gain multiplies both the modulated and the unmodulated
term of the interferogram, and the offset adds $g\,o$ to the modulated
amplitude, so the fit recovers the same $d$ as in the spectral domain.

## Sky reconstruction

The sky is reconstructed from the calibrated spectra
$\tilde d = d / \hat g - \hat o$, with $(\hat g, \hat o)$ the calibration
solution chosen by `calibration.solution`: `nominal` ($\hat g = 1$,
$\hat o = 0$) or `true` (the simulated response, an oracle). With the BBIR
temperatures known exactly, the sky spectrum is recovered from each DU in a
science Mode as

$$
\hat S(\nu) = B_q(\nu) + \sigma\, s\, \tilde d(\nu), \qquad
\sigma = \begin{cases} +1 & \mathrm{B10},\ q = 2 \\ -1 & \mathrm{B01},\ q = 1 \end{cases}
$$

with $q$ the Port looking at the BBIR. In each frequency bin the estimates of
the L and R DUs covering it are averaged. B00 and B11 IFGs give no sky
estimate.

A calibration error leaves, per DU,

$$
\hat S - S = \left(\frac{g}{\hat g} - 1\right)(S - B_q)
+ \sigma\, s \left(\frac{g}{\hat g}\, o - \hat o\right):
$$

the gain error scales the difference between sky and BBIR, so it vanishes
when the BBIR matches the sky (the null of an absolute FTS), while the offset
error flips sign between B10 and B01 and cancels in maps with equal time in the
two Modes. The noise variance of the estimate is divided by $\hat g^2$.
