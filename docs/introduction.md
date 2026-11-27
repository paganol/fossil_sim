# Introduction

## The FOSSIL measurement

FOSSIL observes the sky with an absolute FTS in the spirit of PIXIE
{cite}`kogut2011`. The spectrometer has two inputs (Ports) and two outputs;
at each output a dichroic splits the band at 200 GHz, giving four detector
units (DUs). The modulated part of each output is proportional to the
*difference* of the two input spectra, with opposite signs at the two outputs.
Two beam-switching mechanisms (BSMs) send each Port either to the sky, through
its own telescope, or to an internal blackbody reference (BBIR). The
combination of the two BSM states is the observing Mode.

The spacecraft orbits the Sun–Earth L2 point and spins about an axis kept near
the anti-Sun direction, so that each spin traces a (nearly) great circle on the
sky; the axis is stepped to follow the Sun and covers the full sky in about six
months.

## What fossil_sim simulates

fossil_sim generates the mission one Cooler cycle at a time
({doc}`adr/0001-stateless-cooler-cycles`), and for every interferogram (IFG)
computes:

- the observation timeline: indices, Mode, flags, BBIR temperatures;
- the pointing: spacecraft orbit, Spin axis, boresight, Solar Aspect Angle;
- the sky signal entering each Port and the differential spectrum measured by
  each DU, through its response (gain, offset, gain drift);
- the detector noise: photon noise from the optical loading, white or with a
  1/f component;
- the sky estimate from the calibrated spectra, streaming hit maps and binned
  spectral maps with their variance.

It works at two levels. In the **spectral domain** each IFG is represented
by its spectrum, with Gaussian noise per bin, for mission-length studies. In
the **IFG domain** (`ifg.simulate = true`) the interferogram samples are
generated, with the sky changing during the IFG and noise on the samples,
and the spectrum of each DU is recovered by least squares
({doc}`concepts/interferogram`).

## Status

| Version | Content |
|---|---|
| 0.1 | Configuration, orbit, scanning strategy, timeline, hit maps |
| 0.1.1 | Isotropic backgrounds with exact Boost, signal for all Modes, spectral maps |
| 0.1.2 | Static sky: pysm3 foregrounds as alm, sampled with ducc0, orbital Doppler |
| 0.1.3 | Zodiacal light (zodipy) seen from the spacecraft, evaluated Ring by Ring |
| 0.1.4 | Spectral-domain white noise per DU, variance maps |
| 0.1.5 | Mode switches by observing-time blocks with a repeating pattern, science IFGs selected by Mode; first tagged release |
| 0.1.6 | Static-sky Doppler slope and spectral white-noise normalization corrected |
| 0.1.7 | IFG domain: interferograms, sample noise and least-squares spectra; faster spectral domain |
| 0.1.8 | Configuration and numerical edge cases hardened, end-to-end IFG-domain tests |
| 0.1.9 | DU response (gain, offset) and calibration solution; IFG generation on a finer, wider grid |
| `main` (unreleased) | Photon noise from the optical loading (baseline), Static sky through a low-rank factorization, 1/f detector noise, gain drift |

The {doc}`changelog` lists every change.
Several design ideas come from litebird_sim, which is used only as an optional
validation oracle, and from bisou_sky for the sky model; fossil_sim depends on
neither.
