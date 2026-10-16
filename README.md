# fossil_sim

Simulation framework for FOSSIL, a proposed ESA M-class mission measuring CMB
spectral distortions with an absolute Fourier Transform Spectrometer.

```sh
uv sync --extra dev
uv run pytest
```

- Domain language: [CONTEXT.md](CONTEXT.md)
- Design decisions: [docs/adr](docs/adr)

## Conventions

Same rules as [litebird_sim](https://github.com/litebird/litebird_sim), from
which several design ideas are taken (code is rewritten, it is not a
dependency):

- measure units are always explicit in variable, parameter and config-key
  names, lowercase suffix: `freq_ghz`, `time_s`, `temperature_k`,
  `theta_rad`, `v_km_s`, `intensity_mjy_sr`;
- every function has a numpy-style docstring; every module has tests in `test/`;
- `ruff format` and `ruff check` must pass;
- numba is used for per-sample kernels.

Internal units: GHz, MJy/sr, s (relative to `mission_start`), rad, K, km/s.
