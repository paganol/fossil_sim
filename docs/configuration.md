# Configuration

A simulation is fully defined by a TOML file, read into frozen dataclasses by
{func}`fossil_sim.config.load_config` and validated on load. Every numerical
parameter lives there, with a comment giving its source; values still to be
fixed are marked "TBD — placeholder". The configuration hash and the git commit
identify a run ({meth}`fossil_sim.config.Config.provenance`).

| Section | Content |
|---|---|
| `[simulation]` | master seed, start of science operations, mission duration |
| `[timeline]` | IFG duration, IFGs per Scan, Scans per Ring, Rings per cycle, Recycling Scans |
| `[scanning]` | spin–boresight angle, Ring-stepping rule, Precession, maximum SAA |
| `[orbit]` | L2 distance and Lissajous parameters |
| `[operations]` | Mode pattern, mode-switch dead time, Recycling data and Mode |
| `[bbir]` | set-point schedule, per-Port offsets, HF cavity |
| `[maps]` | HEALPix resolutions |
| `[instrument]` | frequency grid, band, dichroic split, beam |
| `[sky]` | Isotropic backgrounds, solar velocity, beam on the dipole |
| `[static_sky]` | pysm3 models, resolution, cache, sampling, Doppler |
| `[zodi]` | zodipy model, extrapolation, quadrature, parallelism |
| `[noise]` | NEPs, étendue, efficiencies, Recycling noise factor |
| `[output]` | per-cycle TOD files |
| `[ifg]` | OPD range, samples per IFG, OPD per mirror displacement |

## Baseline

```{literalinclude} ../configs/baseline.toml
:language: toml
```
