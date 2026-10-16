# Scanning strategy

The geometry follows litebird_sim's spinning strategy (spin, cone around the
anti-Sun direction, yearly motion), with one difference: the Spin-axis
longitude and the Precession phase are frozen within each Ring, so that its
Scans repeat the same circle.

## Spin axis

Ring $r$ of cycle $c$ has a nominal axis in the ecliptic plane,
$\hat{\mathbf a}_0 = (\cos\lambda_r, \sin\lambda_r, 0)$, with longitude
$\lambda_r$ given by the Ring-stepping rule below. With Precession of opening
$A$ the axis is

$$
\hat{\mathbf a} = \cos A\,\hat{\mathbf a}_0 + \sin A\,(\cos\varphi_p\,\hat{\mathbf e} + \sin\varphi_p\,\hat{\mathbf z}),
\qquad \hat{\mathbf e} = (-\sin\lambda_r, \cos\lambda_r, 0),
$$

$$
\varphi_p = \frac{2\pi\,t_\mathrm{ring}}{P} + \varphi_0,
$$

with $t_\mathrm{ring}$ the mid-time of the Ring and $P$ the Precession
period: the axis lies on a cone of opening $A$ around the nominal axis.

## Ring stepping

Let $\Delta = \lambda_{\bar\odot}(t_{c+1}) - \lambda_{\bar\odot}(t_c)$ be the
motion of the anti-Sun longitude over cycle $c$ (about 1.44°) and
$n = N_R$.

`equal_per_cycle` (default)
: $\lambda_r = \lambda_{\bar\odot}(t_c) + (r + o)\,\Delta / n$: equal steps of
  $\Delta / n \approx 0.48°$, also across the Recycling window, which gives a
  uniform longitude coverage.

`anchor_ring_start`
: $\lambda_r = \lambda_{\bar\odot}$ at the start of each Ring; steps of
  about 0.38° within the cycle and 0.67° across the Recycling window.

`fixed`
: $\lambda_r = \lambda_{\bar\odot}(0) + (c\,n + r)\,\delta$ with a fixed step
  $\delta$; without re-anchoring the axis drifts away from the anti-Sun
  direction by $n\delta - \Delta$ per cycle (about 0.06° for
  $\delta = 0.5°$, i.e. 10° in about 7 months). For diagnostics only.

During the Recycling window the axis stays at the last Ring.

### Offset of the equal-step rule

Assume a uniform solar motion over the cycle and let
$f = N_S / N_\mathrm{cyc}$ be the fraction of the cycle spent in one Ring. In
units of $\Delta$, the anti-Sun longitude relative to the axis of Ring $r$
spans

$$
\left[\, r f - \frac{r + o}{n},\; (r+1) f - \frac{r + o}{n} \,\right],
$$

and during the Recycling window (axis of Ring $n - 1$) it reaches
$1 - (n - 1 + o)/n$. Since $f < 1/n$, the minimum is at the start of the last
Ring, $(n-1)(f - 1/n) - o/n$, and the maximum at the end of recycling,
$(1 - o)/n$. Centring the excursion gives

$$
o = \frac{(n-1)(n f - 1) + 1}{2},
\qquad
\max |\lambda_{\bar\odot} - \lambda_\mathrm{axis}| = \frac{1 - (n-1)(n f - 1)}{2n}\,\Delta .
$$

For the baseline ($n = 3$, $f = 23/87$): $o = 0.293$ and a maximum offset of
$0.236\,\Delta \approx 0.34°$. With `anchor_ring_start` the maximum offset is
the solar motion over one Ring plus the Recycling window,
$(N_S + N_\mathrm{rec})/N_\mathrm{cyc}\,\Delta \approx 0.68°$.

## Boresight

The boresight is at angle $\theta_b$ from the Spin axis
(`spin_boresight_angle_deg`, 90° for great circles). With the reference
vectors

$$
\hat{\mathbf u} = \frac{\hat{\mathbf z} - (\hat{\mathbf z}\cdot\hat{\mathbf a})\,\hat{\mathbf a}}{|\hat{\mathbf z} - (\hat{\mathbf z}\cdot\hat{\mathbf a})\,\hat{\mathbf a}|},
\qquad \hat{\mathbf w} = \hat{\mathbf a} \times \hat{\mathbf u},
$$

the boresight of IFG $k$ in its Scan is

$$
\hat{\mathbf b} = \cos\theta_b\,\hat{\mathbf a} + \sin\theta_b\,(\cos\phi_k\,\hat{\mathbf u} + \sin\phi_k\,\hat{\mathbf w}),
\qquad \phi_k = \frac{2\pi\,(k + \frac12)}{N_I}.
$$

Spin phase $\phi = 0$ is the point of the circle closest to the north
ecliptic pole and the spin is right-handed about $\hat{\mathbf a}$. The angular
step between consecutive IFGs is $\simeq 2\pi \sin\theta_b / N_I$ (exact
for $\theta_b = 90°$), 0.5° for the baseline, i.e. 3.2 samples per 1.6° beam.

The same pointing is obtained with litebird_sim's quaternion sequence
$R_z(\lambda_r)\,R_x(\varphi_p - \pi/2)\,R_y(\pi/2 - A)\,R_z(\phi + \pi)$ applied
to $(\sin\theta_b, 0, \cos\theta_b)$; the optional oracle test checks this to
better than 1 nrad.

## Solar Aspect Angle

$$
\mathrm{SAA}(t) = \arccos\left(\hat{\mathbf a} \cdot \hat{\mathbf n}_{\bar\odot}(t)\right),
$$

evaluated for every IFG and flagged above `max_saa_deg`. The configuration
rejects Precession amplitudes with $A + 0.75° > \mathrm{SAA}_\mathrm{max}$, and
warns when $|90° - \theta_b| > A$, since the ecliptic poles are then never
observed.
