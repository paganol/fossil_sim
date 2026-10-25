# Boost

The Isotropic backgrounds are evaluated in the CMB rest frame and boosted
exactly to the spacecraft frame for every IFG.

## Transformation of the specific intensity

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

## Blackbody case

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

## Velocity

$$
\boldsymbol\beta(t) = \frac{\mathbf v_\odot + \mathbf v_\mathrm{sc}(t)}{c},
$$

with $\mathbf v_\odot$ the solar velocity with respect to the CMB (369.816
km s⁻¹ towards Galactic $(l, b) = (264.021°, 48.253°)$, {cite}`planck2018i`;
ecliptic $(\lambda, \beta) \simeq (171.66°, -11.15°)$) and
$\mathbf v_\mathrm{sc}$ the barycentric velocity of the spacecraft
({doc}`orbit`). The orbital term modulates the dipole amplitude by about
±30 km s⁻¹ over the year.

## Beam

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

