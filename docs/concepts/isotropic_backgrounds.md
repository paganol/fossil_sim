# Isotropic backgrounds

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

