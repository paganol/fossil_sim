# Tutorials

The notebooks are stored executed in `notebooks/`; they are regenerated with
`scripts/build_notebooks.py` and executed in CI.

```{toctree}
:maxdepth: 1

notebooks/fossil_sim_scanning
notebooks/fossil_sim_interferogram
notebooks/fossil_sim_spectral_maps
```

## Hit maps and report

`examples/hitmaps.py` runs 30 days and 6 months of observations and writes
FITS hit maps (total and per Mode), figures and a markdown report:

```sh
uv run python examples/hitmaps.py [config.toml] [output_dir]
```
