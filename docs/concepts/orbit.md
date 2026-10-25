# Orbit

The spacecraft position is needed for the anti-Sun direction (scanning
strategy, Solar Aspect Angle) and its velocity for the kinematic dipole
({doc}`boost`). The model, implemented in {class}`fossil_sim.orbit.L2Orbit`,
follows litebird_sim: an approximate L2 point plus a Lissajous orbit.

## Ephemerides

The barycentric position $\mathbf r_\oplus(t)$ and velocity
$\mathbf v_\oplus(t)$ of the Earth and the position of the Sun
$\mathbf r_\odot(t)$ are computed with astropy (built-in ERFA ephemerides) on
an hourly grid, rotated from ICRS to the J2000 mean ecliptic with a fixed
matrix, and interpolated with cubic splines.

## Approximate L2

L2 is placed on the line from the barycentre through the Earth, at distance
$d_{L2}$ beyond the Earth:

$$
\mathbf r_{L2} = s\,\mathbf r_\oplus, \qquad s = 1 + \frac{d_{L2}}{|\mathbf r_\oplus|}.
$$

Differentiating, with $\hat{\mathbf r} = \mathbf r_\oplus / |\mathbf r_\oplus|$
and $v_r = \mathbf v_\oplus \cdot \hat{\mathbf r}$,

$$
\mathbf v_{L2} = s\,\mathbf v_\oplus - (s - 1)\, v_r\, \hat{\mathbf r}:
$$

the radial velocity is unchanged and the tangential one is scaled by $s$.

## Lissajous orbit

Around L2 the spacecraft follows

$$
\Delta\mathbf r = R_1 \cos(\omega_1 t)\,\hat{\mathbf t} + R_2 \sin(\omega_2 t + \phi)\,\hat{\mathbf z},
\qquad \hat{\mathbf t} = \hat{\mathbf z} \times \hat{\mathbf r},
$$

an along-track and an out-of-ecliptic oscillation. Its exact velocity is

$$
\Delta\mathbf v = -R_1 \omega_1 \sin(\omega_1 t)\,\hat{\mathbf t}
+ R_1 \cos(\omega_1 t)\,\dot{\hat{\mathbf t}}
+ R_2 \omega_2 \cos(\omega_2 t + \phi)\,\hat{\mathbf z},
$$

$$
\dot{\hat{\mathbf t}} = \hat{\mathbf z} \times \dot{\hat{\mathbf r}}, \qquad
\dot{\hat{\mathbf r}} = \frac{\mathbf v_\oplus - v_r\,\hat{\mathbf r}}{|\mathbf r_\oplus|}.
$$

The default amplitudes and frequencies are those of the WMAP orbit, used as
placeholders. The model has no radial component, whereas real L2 orbits have
one of comparable amplitude (about 0.1 km s⁻¹ in velocity); a reference orbit
can replace {class}`~fossil_sim.orbit.L2Orbit` by providing the same methods.

## Anti-Sun direction

The anti-Sun direction seen from the spacecraft uses the true position of the
Sun, which differs from the barycentre by up to about 0.01 AU:

$$
\hat{\mathbf n}_{\bar\odot}(t) = \frac{\mathbf r_\mathrm{sc}(t) - \mathbf r_\odot(t)}{|\mathbf r_\mathrm{sc}(t) - \mathbf r_\odot(t)|},
\qquad \mathbf r_\mathrm{sc} = \mathbf r_{L2} + \Delta\mathbf r.
$$

Its ecliptic longitude $\lambda_{\bar\odot}(t)$ drives the scanning strategy;
it advances by about 0.986° per day, between 0.953° (July) and 1.019°
(January) because of the eccentricity of the Earth orbit.
