# Interferogram

An ideal-FTS interferogram model ({mod}`fossil_sim.ifg`) illustrates the
measurement; instrumental systematics (mirror trajectory, time constants,
phase errors, arm imbalance) are planned for the IFG-domain simulation.

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
