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
output splits the band at 300 GHz into a low-frequency and a high-frequency
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
