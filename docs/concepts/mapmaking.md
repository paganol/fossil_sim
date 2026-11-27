# Map-making

All maps are accumulated in streaming, one Cooler cycle at a time, from Valid
IFGs only, in J2000 ecliptic coordinates (HEALPix RING). With MPI each rank
accumulates its own cycles and the sums are reduced at the end.

## Hit maps

For each resolution and Mode $m$, $N_{p,m}$ counts the Valid IFGs whose
boresight falls in pixel $p$. The time of the first valid hit of each pixel,
$t_p$, gives the coverage curve

$$
f(t) = \frac{1}{N_\mathrm{pix}}\,\#\{p : t_p \le t\}.
$$

The boresight-only hit maps sample the sky at the pixel containing the
boresight; at nside 128 the pixel size (0.46°) is smaller than the Ring step
(0.48°), so about 13 % of the pixels stay empty even though the 1.6° beam
samples the sky completely.

## Spectral maps

The reconstructed-sky map is the binned average

$$
m_p(\nu) = \frac{1}{N_p} \sum_{i \in p} \hat S_i(\nu),
$$

over the Valid science IFGs $i$ in pixel $p$. Diagnostic maps bin the
differential spectrum $d(\nu)$ of each output (L, R) separately for B10 and
B01. Pixels never observed are set to `healpy.UNSEEN`.

When noise is on, the white-noise variance map of the sky is
$\sigma^2(\nu) / N_p$, with $\sigma^2(\nu)$ the variance per IFG of
$\hat S(\nu)$ (the mean of the L and R estimates: half the per-DU variance).
It is exact because the noise level is the same for all Valid IFGs.
For simulated interferograms the per-bin variance is propagated through the
actual least-squares solver, including its finite-sampling effects
({doc}`noise`).

The noise of different bins can be correlated (the fit, apodization, 1/f
noise; {doc}`noise`). For the same reason the covariance between bins of a
pixel is $C/N_p$, with $C$ the covariance per IFG of $\hat S$
({func}`fossil_sim.noise.sky_estimate_noise_covariance`): one
$n_\mathrm{freq}\times n_\mathrm{freq}$ matrix for the whole map
({meth}`~fossil_sim.mapmaking.SpectralMaps.sky_covariance`,
{meth}`~fossil_sim.mapmaking.SpectralMaps.sky_correlation`). The variance
maps are its diagonal. With noise on, the outputs include
`sky_noise_covariance.h5`: $C$ (`covariance_per_ifg`), the bin frequencies
and the hits $N_p$ of the sky map per nside (`hits/nside_NNNN`).

This is the optimal estimator for white noise of equal variance and no
time-correlated systematics; more general map-makers will be needed once
noise and calibration errors are introduced.
