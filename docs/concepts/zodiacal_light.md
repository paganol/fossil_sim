# Zodiacal light

The thermal emission of the interplanetary dust is computed with zodipy, which
integrates each component of the dust cloud along the line of sight,

$$
I_\mathrm{zodi}(\nu, \hat{\mathbf n}, t) = \sum_c \int_0^\infty n_c(\mathbf r)\, E_c(\nu)\, B_\nu\big(T(r)\big)\, ds,
\qquad \mathbf r = \mathbf r_\mathrm{sc}(t) - \mathbf r_\odot(t) + s\,\hat{\mathbf n},
$$

with $n_c$ the density of component $c$, $E_c$ its emissivity and $T(r)$ the
dust temperature at heliocentric distance $r$. The observer is the
spacecraft position of {doc}`orbit`, passed explicitly to zodipy, so that
orbit, pointing and zodiacal light are consistent.

The default model is `planck18`, calibrated on Planck HFI between 100 and
857 GHz and extrapolated outside that range (an uncertainty of about 20 %
above 857 GHz). The `dirbe` model, calibrated above 1249 GHz, becomes negative
when extrapolated below about 500 GHz, and blending the two would introduce an
artificial spectral feature between 857 and 1249 GHz.

## Evaluation Ring by Ring

A direct evaluation for every IFG and frequency costs about 3 minutes per
Cooler cycle. Since the Scans of a Ring retrace the same circle, only the
observer moves during a Ring (about 0.4° around the Sun in 9.2 h). For each
Ring $r$ the emission is therefore evaluated on the $N_I$ directions of the
circle at the start $t_r^a$ and end $t_r^b$ of the Ring (for the last Ring,
the end of the Recycling window, during which the axis stays fixed), and
interpolated linearly for IFG $i$ at phase $k$:

$$
I_\mathrm{zodi}(\nu, i) = I_\mathrm{zodi}(\nu, \hat{\mathbf b}_k, t_r^a)
+ \frac{t_i - t_r^a}{t_r^b - t_r^a}\left[I_\mathrm{zodi}(\nu, \hat{\mathbf b}_k, t_r^b) - I_\mathrm{zodi}(\nu, \hat{\mathbf b}_k, t_r^a)\right].
$$

This reduces the cost to about 11 s per cycle; the interpolation error,
measured against the direct evaluation in the middle of a Ring, is below
$2\times10^{-5}$ of the signal. The zodiacal light is smooth on the scale of
the beam, which is not applied; the Doppler effect of the dust motion is not
modelled.
