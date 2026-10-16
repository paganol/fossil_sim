# fossil_sim

Simulation framework for FOSSIL, a proposed ESA M-class mission measuring CMB
spectral distortions with an absolute Fourier Transform Spectrometer.

```sh
uv sync --extra dev
uv run pytest                                   # unit tests + doctests
uv run pytest --nbmake notebooks                # example notebooks
uv run --with litebird_sim pytest -m heavy      # optional litebird_sim oracle
uv run --with mpi4py mpirun -n 3 python scripts/mpi_check.py
uv run python examples/hitmaps.py               # hit maps + report in output/hitmaps
```

Notebooks: [scanning](notebooks/fossil_sim_scanning.ipynb),
[interferogram](notebooks/fossil_sim_interferogram.ipynb).

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
