# Installation

fossil_sim requires Python ≥ 3.12. With [uv](https://docs.astral.sh/uv/):

```sh
git clone https://github.com/paganol/fossil_sim.git
cd fossil_sim
uv sync --extra dev
```

or with pip: `pip install -e ".[dev]"`.

Optional extras:

`mpi`
: `mpi4py`, to distribute Cooler cycles over MPI ranks.

`docs`
: Sphinx and extensions, to build this documentation.

`oracle`
: litebird_sim, for the optional pointing cross-check.

## Running the tests

```sh
uv run pytest                                  # unit tests and doctests
uv run pytest --nbmake notebooks               # example notebooks
uv run --with litebird_sim --with "git+https://github.com/paganol/BISOU-sky" \
    pytest -m heavy                            # cross-checks (litebird_sim, bisou_sky)
uv run --with mpi4py mpirun -n 3 python scripts/mpi_check.py
```

## Building the documentation

```sh
uv run --extra docs sphinx-build -W docs docs/_build/html
```
