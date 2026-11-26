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
two of them (Coulon et al. 2026, Eq. 11–12); the NEP must refer to the same
detector.

The noise is added to the differential spectrum $d(\nu)$ of each DU, so that
it propagates to the sky estimate (with `ifg.simulate = true` it is drawn on the
IFG samples instead, $\mathrm{NEP}/\sqrt{2\,\delta t}$, and reaches $d(\nu)$
through the fit, {doc}`interferogram`); IFGs of the Recycling window get
`recycling_noise_factor` times this level. The baseline computes the NEP of
each DU from its loading (`nep_model = "photon"`, below); with
`nep_model = "fixed"` it is given per band. The NEP_total of the M8 proposal
(Table 2.3: 4.8×10⁻¹⁶ LF, 1.3×10⁻¹⁴ W/√Hz HF) is not a photon NEP of the
FOSSIL DUs: the HF value is about 70 times the photon noise of the CMB
loading, and the FOSSIL sensitivity curve (arXiv:2609.17618, Fig. 1) excludes
it. The Step2 values (Table 4.1: 4.6 and 3.8×10⁻¹⁶) are of the right order.

The formula for $\delta S_\nu$ assumes orthogonal cosines. With
`ifg.simulate = true`, the variance maps instead propagate the sample noise
through the actual least-squares solver $L$:
$\operatorname{Var}(d_k) = \sigma_\mathrm{sample}^2 \sum_j L_{kj}^2$.
This includes the fitted constant and the finite OPD sampling; near Nyquist
the bin variances can differ substantially from the nominal formula. The fit
can also correlate frequency bins; the variance maps store only the diagonal
of that covariance. The spectral-domain noise model retains independent bins
with the nominal variance.

With the baseline (photon noise, M8 bands 30–200–2000 GHz, $m = 1/2$) the
NEP is 1.2×10⁻¹⁶ (LF) and 1.9×10⁻¹⁶ W/√Hz (HF), and
$\delta S_\nu \simeq 7.6\times10^{3}$ and $1.2\times10^{4}$ Jy/sr per IFG
and DU. The HF DU is noisier than the LF one here because the 200 GHz split
puts the CMB peak in its band (with the 300 GHz split of `configs/step2.toml`
the order is reversed).

## Photon noise from the optical loading

With `nep_model = "photon"` the NEP of each DU is not given but computed from
its optical loading, following Coulon, Aghanim & Maffei (2026,
arXiv:2609.17829, Sect. 4), the model behind the FOSSIL sensitivity curves
(FOSSIL Collaboration 2026, arXiv:2609.17618, Fig. 1)
({func}`fossil_sim.noise.photon_nep_w_sqrt_hz`). The spectral power on the
DU is the unmodulated term of its interferogram, all frequencies adding
incoherently (an upper bound on the photon noise),

$$
P_\nu = m\,\alpha\,A\Omega\,w(\nu)\left[f_\mathrm{opt}\,(I_\mathrm{sky} + I_\mathrm{BBIR}) + \sum_k \varepsilon_k B_\nu(T_k)\right],
$$

with $w$ the DU transmission, $I_\mathrm{sky}$ the Isotropic backgrounds (no
Galactic or zodiacal emission), $I_\mathrm{BBIR}$ the BBIR at the mean set
point and the optical elements as grey bodies (`optics_emissivity`,
`optics_temperature_k`) whose emission reaches the detectors unattenuated.
Shot and wave noise give

$$
\mathrm{NEP}^2 = \int \left[2h\nu P_\nu + \frac{2 c^2}{\nu^2 A\Omega} P_\nu^2\right] d\nu
$$

over the band of the DU. For a Martin-Puplett FTS with DUs absorbing both
polarisations $m = 1/2$ (Coulon et al., Eq. 11–12):
$P = \tfrac12\left[(I_1 + I_2) \pm (I_1 - I_2)\cos\right]$. With this $m$,
$\delta S_\nu$ above, averaged over the L and R outputs and over the mission,
is Coulon's sensitivity per channel,
$\sqrt2\,\mathrm{NEP}/(\sqrt{t_\mathrm{obs}}\,A\Omega\,\Delta\nu\,\eta\,T_\mathrm{tot})$
(Eq. 17–18). Their reference (Table 2: $\eta = 0.5$, $A\Omega = 3$ cm² sr,
one optical element at 1 % and 4.5 K, $t_\mathrm{obs}$ = 2.8 yr), with sky
and reference at the CMB temperature, gives NEP = 1.7 (1.35)×10⁻¹⁶ W/√Hz and
1.14 (0.91) Jy/sr in the 50–300 (300–2000) GHz band, which the tests
reproduce to 1 %; `configs/step2.toml` uses these values. The step of the
FOSSIL curve just above 300 GHz is not reproduced: it needs frequency-dependent
efficiencies, not yet modelled.

Each realization is drawn from a `numpy.random.SeedSequence` keyed by
(master seed, component, cycle, DU), never by MPI rank (ADR 0001): a run gives
the same noise for any number of ranks.
