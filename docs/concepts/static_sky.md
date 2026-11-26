# Static sky

The Static sky contains the anisotropic, time-independent emission: Galactic
foregrounds and CMB anisotropies from pysm3, following the design of
litebird_sim's `SkyGenerator`.

## Generation

For every frequency bin $\nu_k$ (a delta bandpass at the bin centre,
consistent with the Isotropic backgrounds):

1. pysm3 evaluates the total intensity of the chosen models in MJy sr⁻¹ at
   $N_\mathrm{side}^\mathrm{eval} = 2\,N_\mathrm{side}$, Galactic coordinates;
2. the map is expanded in spherical harmonics up to $\ell_\mathrm{max}$
   (`healpy.map2alm`, 3 iterations);
3. the models are split by rest frame
   ({doc}`../adr/0003-static-sky-barycentric-frame`): Galactic models in the
   barycentric frame, CMB anisotropies (`c*`) and extragalactic models
   (`cib*`, `tsz*`, `ksz*`, `rg*`) in the CMB frame; for the latter the
   monopole is removed, $a_{00} = 0$, since it belongs to the Isotropic
   backgrounds; Galactic models keep their physical monopole and dipole;
   `dip*` and `co*` are rejected by the configuration;
4. the Gaussian beam is applied, $a_{\ell m} \to b_\ell\, a_{\ell m}$;
5. the alm are rotated from Galactic to J2000 ecliptic coordinates.

The cube $a_{\ell m}(\nu_k)$ is cached in HDF5 under a name derived from a
hash of everything that determines it (models, resolutions, beam, frequency
grid, pysm3 version, Doppler step), so it is built once and reused. With
`doppler = true` it also holds the two cubes of $K$ (below).

## Sampling

Each IFG samples the cube at the exact boresight direction,

$$
I_\mathrm{static}(\nu_k, \hat{\mathbf n}) = \sum_{\ell \le \ell_\mathrm{max}} \sum_{m} a_{\ell m}(\nu_k)\, Y_{\ell m}(\hat{\mathbf n}),
$$

with `ducc0.sht.synthesis_general` (relative accuracy $10^{-9}$). No pixel
window is involved. A faster nearest-pixel sampling of maps at
$N_\mathrm{side}$ is available for quick runs.

The synthesis does not run once per frequency: the spectra of the pysm3
models span a space of few dimensions, so each cube is factorized in
frequency,

$$
a_{\ell m}(\nu_k) \simeq \sum_{j=1}^{K} c_j(\nu_k)\, b^{(j)}_{\ell m},
$$

from the SVD of the cube with each frequency normalized to its own norm
({func}`fossil_sim.input_sky.compress_alm_cube`). $K$ is the smallest rank
that keeps the map of every frequency within `rank_epsilon` of its rms
(default $10^{-9}$, the accuracy of the synthesis): with real $c_j$ the $b^{(j)}$ are the alm of real
maps, and only the $K$ basis maps are synthesized. At baseline (models c1,
d2, s1, a1, f1; 690 frequencies up to 2.1 THz) $K$ = 37 for $I$, 49 and 1
for the barycentric and cosmological $K$ cubes: 87 syntheses instead of 2070.
The factorization is cached next to the cubes, one file per threshold
(`*_low_rank_1e-09.h5`, about two minutes of SVD to build), and once cached
the cubes are not loaded. `rank_epsilon = 0` turns it off: one synthesis per
frequency, and in the IFG domain the sky windows on the full generation grid
(the reference, about three times slower on the 690-frequency grid).

Since the Scans of a Ring retrace the same circle, the boresight of an IFG
depends only on its Ring and its phase in the Scan. In the spectral domain the
cubes are therefore sampled once per cycle on the $N_I$ IFG centres of each
Ring circle and read back for every IFG: the result equals sampling each IFG
within the synthesis accuracy, at $1/30$ of the cost at baseline.

## Doppler

Each part of the Static sky is seen from a frame moving with velocity
$\boldsymbol\beta$ relative to its own rest frame
({doc}`../adr/0003-static-sky-barycentric-frame`): the orbital velocity
$\boldsymbol\beta_\mathrm{orb}$ for the Galactic (barycentric) part,
$\boldsymbol\beta_\mathrm{sun} + \boldsymbol\beta_\mathrm{orb}$ for the
cosmological part. Expanding the Boost to first order,
$\nu' \simeq \nu\,(1 - \boldsymbol\beta\cdot\hat{\mathbf n})$ and

$$
I_\mathrm{obs}(\nu) \simeq \left(1 + 3\,\boldsymbol\beta\cdot\hat{\mathbf n}\right)\left[I(\nu) - \boldsymbol\beta\cdot\hat{\mathbf n}\, \frac{dI}{d\ln\nu}\right]
\simeq I(\nu) + \boldsymbol\beta\cdot\hat{\mathbf n}\; K(\nu),
\qquad K = 3\,I - \frac{dI}{d\ln\nu}.
$$

$K$ is built with the cube: pysm3 is evaluated at $\nu\,e^{\pm\varepsilon}$
($\varepsilon = 10^{-2}$) and
$dI/d\ln\nu \simeq [I(\nu e^{\varepsilon}) - I(\nu e^{-\varepsilon})]/(2\varepsilon)$,
exact to $\alpha^3\varepsilon^2/6 \sim 10^{-4}$ of $I$ for a power law
$I \propto \nu^\alpha$, at every bin including the band edges; the maps of $K$
go through the same beam and rotation as $I$. For a power law
$K = (3 - \alpha)\,I$. Each IFG samples $I$ and the two $K$ cubes (three
syntheses). With the baseline sky the orbital term reaches $10^{-2}$ MJy sr⁻¹
on the Galactic plane at 1 THz; the solar term, the Doppler modulation of the
CMB anisotropies, is about 40–80 Jy sr⁻¹ rms at 150–300 GHz.
Aberration (about 20″ orbital, 4′ solar) and $O(\beta^2)$ terms are
neglected.

