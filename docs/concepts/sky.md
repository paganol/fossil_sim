# Sky model

The sky is split into two parts ({doc}`../adr/0002-isotropic-backgrounds-in-time-domain`):

- the **Isotropic backgrounds**, isotropic in the CMB rest frame, evaluated
  for every IFG with the exact Boost to the spacecraft frame;
- the **Static sky**, the time-independent anisotropies (Galactic
  foregrounds, CMB anisotropies), defined in the barycentric frame.

## Boost

### Transformation of the specific intensity

The photon occupation number, $I_\nu / \nu^3$ up to constants, is a Lorentz
invariant. A photon observed at frequency $\nu$ arriving from the direction
$\hat{\mathbf n}$ (the observed direction, pointing to the sky) by an observer
moving with velocity $\boldsymbol\beta c$ with respect to the CMB has, in the
CMB rest frame, frequency

$$
\nu' = \nu\,\gamma\,(1 - \boldsymbol\beta\cdot\hat{\mathbf n}), \qquad
\gamma = (1 - \beta^2)^{-1/2}.
$$

The observed intensity of any isotropic rest-frame spectrum
$I_\mathrm{rest}(\nu)$ is therefore

$$
I_\mathrm{obs}(\nu, \hat{\mathbf n}) = \left(\frac{\nu}{\nu'}\right)^3 I_\mathrm{rest}(\nu').
$$

fossil_sim applies this to the *sum* of all Isotropic backgrounds, which
yields their monopole, dipole and all higher kinematic multipoles at once.
Aberration does not affect an isotropic background.

### Blackbody case

For a blackbody at $T_0$ the occupation number $1/(e^{h\nu'/kT_0} - 1)$ equals
that of a blackbody at

$$
T(\hat{\mathbf n}) = \frac{T_0}{\gamma\,(1 - \boldsymbol\beta\cdot\hat{\mathbf n})}
= T_0\left[1 + \beta\mu + \beta^2\left(\mu^2 - \tfrac12\right) + O(\beta^3)\right],
\qquad \mu = \hat{\boldsymbol\beta}\cdot\hat{\mathbf n},
$$

the familiar dipole $\Delta T = T_0\,\beta\cos\theta$ (about 3.36 mK for the
solar velocity) plus a kinematic quadrupole and a monopole shift at order
$\beta^2$. This exact result is used to test the implementation.

### Velocity

$$
\boldsymbol\beta(t) = \frac{\mathbf v_\odot + \mathbf v_\mathrm{sc}(t)}{c},
$$

with $\mathbf v_\odot$ the solar velocity with respect to the CMB (369.816
km s⁻¹ towards Galactic $(l, b) = (264.021°, 48.253°)$, {cite}`planck2018i`;
ecliptic $(\lambda, \beta) \simeq (171.66°, -11.15°)$) and
$\mathbf v_\mathrm{sc}$ the barycentric velocity of the spacecraft
({doc}`orbit`). The orbital term modulates the dipole amplitude by about
±30 km s⁻¹ over the year.

### Beam

The beam does not change a monopole but smooths the dipole. For a Gaussian
beam of width $\sigma = \mathrm{FWHM}/\sqrt{8\ln 2}$ the window function is
$b_\ell = e^{-\ell(\ell+1)\sigma^2/2}$, so the dipole is multiplied by

$$
b_1 = e^{-\sigma^2},
$$

$1 - b_1 = 1.4\times 10^{-4}$ for FWHM = 1.6°. fossil_sim applies $b_1$ to the
anisotropic part of the boosted signal,

$$
I = I_\mathrm{rest}(\nu) + b_1\left[I_\mathrm{obs}(\nu, \hat{\mathbf n}) - I_\mathrm{rest}(\nu)\right],
$$

which is exact for the dipole and neglects the beam on the $O(\beta^2)$ terms.

## Rest-frame spectra

With $x = h\nu / kT_0$ and $I_0 = 2(kT_0)^3/(hc)^2$
(270.06 MJy sr⁻¹ for $T_0 = 2.72548$ K {cite}`fixsen2009`):

CMB blackbody
: $B_\nu(T_0) = I_0\, x^3 / (e^x - 1)$.

Compton-$y$ distortion {cite}`zeldovich1969`
: $\Delta I_y = y\, I_0\, \dfrac{x^4 e^x}{(e^x - 1)^2}\left[x \coth\frac{x}{2} - 4\right]$,
  with a null at 217.5 GHz.

Relativistic corrections {cite}`itoh1998`
: $\Delta I_\mathrm{rSZ} = y\, I_0\, \dfrac{x^4 e^x}{(e^x - 1)^2} \sum_{k=1}^{4} Y_k(x)\,\theta_e^k$,
  $\theta_e = kT_e / m_e c^2$, where the $Y_k$ are polynomials in
  $\tilde X = x\coth(x/2)$ and $\tilde S = x / \sinh(x/2)$; the first is
  $Y_1 = -10 + \tfrac{47}{2}\tilde X - \tfrac{42}{5}\tilde X^2 + \tfrac{7}{10}\tilde X^3 + \tilde S^2\left(-\tfrac{21}{5} + \tfrac{7}{5}\tilde X\right)$
  and $Y_2$–$Y_4$ are given in the reference and in
  {func}`fossil_sim.spectra.relativistic_sz_correction_mjy_sr`.

Chemical-potential $\mu$ distortion {cite}`sunyaev1970`
: $\Delta I_\mu = \mu\, I_0\, \dfrac{x^4 e^x}{(e^x - 1)^2}\left[\frac{1}{\beta_\mu} - \frac{1}{x}\right]$,
  $\beta_\mu = 2.1923$.

Cosmic infrared background {cite}`fixsen1998`
: $I_\mathrm{CIB} = A\,(\nu / \nu_0)^{\beta}\,B_\nu(T_\mathrm{CIB})$,
  with $A = 1.3\times 10^{-5}$, $\beta = 0.64$, $T_\mathrm{CIB} = 18.5$ K,
  $\nu_0 = 3000$ GHz.

Extragalactic CO {cite}`abitbol2017,mashian2016`
: a tabulated sky-averaged spectrum (10–615 GHz), interpolated in
  $\log I$–$\log\nu$ and set to zero outside the tabulated range.

The baseline amplitudes are the ΛCDM fiducial values $y = 1.77\times10^{-6}$,
$kT_e = 1.24$ keV {cite}`hill2015` and $\mu = 2\times 10^{-8}$
{cite}`chluba2016`. Any component is switched off by setting its amplitude to
zero.

## Static sky

The Static sky contains the anisotropic, time-independent emission: Galactic
foregrounds and CMB anisotropies from pysm3, following the design of
litebird_sim's `SkyGenerator`.

### Generation

For every frequency bin $\nu_k$ (a delta bandpass at the bin centre,
consistent with the Isotropic backgrounds):

1. pysm3 evaluates the total intensity of the chosen models in MJy sr⁻¹ at
   $N_\mathrm{side}^\mathrm{eval} = 2\,N_\mathrm{side}$, Galactic coordinates;
2. the map is expanded in spherical harmonics up to $\ell_\mathrm{max}$
   (`healpy.map2alm`, 3 iterations);
3. for extragalactic models (`cib*`, `tsz*`, `ksz*`, `rg*`) the monopole is
   removed, $a_{00} = 0$, since it belongs to the Isotropic backgrounds;
   Galactic models keep their physical monopole and dipole; `dip*` and `co*`
   are rejected by the configuration;
4. the Gaussian beam is applied, $a_{\ell m} \to b_\ell\, a_{\ell m}$;
5. the alm are rotated from Galactic to J2000 ecliptic coordinates.

The cube $a_{\ell m}(\nu_k)$ is cached in HDF5 under a name derived from a
hash of everything that determines it (models, resolutions, beam, frequency
grid, pysm3 version), so it is built once and reused.

### Sampling

Each IFG samples the cube at the exact boresight direction,

$$
I_\mathrm{static}(\nu_k, \hat{\mathbf n}) = \sum_{\ell \le \ell_\mathrm{max}} \sum_{m} a_{\ell m}(\nu_k)\, Y_{\ell m}(\hat{\mathbf n}),
$$

with `ducc0.sht.synthesis_general` (relative accuracy $10^{-9}$). No pixel
window is involved. A faster nearest-pixel sampling of maps at
$N_\mathrm{side}$ is available for quick runs.

### Orbital Doppler

The pysm3 templates derive from maps averaged over years of observations from
L2, i.e. in the barycentric frame ({doc}`../adr/0003-static-sky-barycentric-frame`).
Relative to it the spacecraft moves with the orbital velocity
$\boldsymbol\beta_\mathrm{orb}$ only. Expanding the Boost to first order,
$\nu' \simeq \nu\,(1 - \boldsymbol\beta_\mathrm{orb}\cdot\hat{\mathbf n})$ and

$$
I_\mathrm{obs}(\nu) \simeq \left(1 + 3\,\boldsymbol\beta_\mathrm{orb}\cdot\hat{\mathbf n}\right)\left[I(\nu) - \boldsymbol\beta_\mathrm{orb}\cdot\hat{\mathbf n}\, \frac{dI}{d\ln\nu}\right]
\simeq I(\nu) + \boldsymbol\beta_\mathrm{orb}\cdot\hat{\mathbf n}\left[3\,I(\nu) - \frac{dI}{d\ln\nu}\right],
$$

with $dI/d\ln\nu$ computed by finite differences on the frequency grid. For
a power law $I \propto \nu^\alpha$ the correction is
$(3 - \alpha)\,\boldsymbol\beta_\mathrm{orb}\cdot\hat{\mathbf n}\, I$. Aberration (about 20″)
and $O(\beta^2)$ terms are neglected.
