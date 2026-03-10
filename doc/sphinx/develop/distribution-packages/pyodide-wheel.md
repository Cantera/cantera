# Building a Pyodide Wheel for Cantera

Cantera can be built as a Pyodide-compatible wheel so the Python module can run in
browser-based environments (for example, JupyterLite).

This page documents the local and CI workflow used in this repository for generating and
smoke-testing that wheel.

## Overview

The Pyodide build follows the same two-stage packaging model as the standard Python
wheel flow:

1. Build an sdist from the Cantera source tree using SCons.
2. Build the wheel from that sdist using the Pyodide toolchain.

The SCons entry point for this is:

```sh
scons pyodide-wheel
```

Output wheels are written to:

```text
build/pyodide_dist/
```

## Requirements

The Pyodide build uses the dedicated `pyodide` Pixi environment in `pixi.toml`.
At minimum, this environment needs:

- Python 3.13
- SCons
- `pyodide-build`
- `build`
- `scikit-build-core`
- CMake and Ninja
- Boost headers (`boost-cpp`)
- Cantera's Python build dependencies (`cython`, `numpy`, `ruamel.yaml`, `pip`)

```{note}
When using `pyodide venv`, the generated interpreter requires `node` to be available on
`PATH`. This can be achieved by running Pyodide venv commands through Pixi as
`pixi run -e pyodide ...` so the right runtime is present.
```

## Local Build

From the repository root:

```sh
pixi run -e pyodide scons pyodide-wheel
```

This creates a wheel named something like:

```text
build/pyodide_dist/cantera-<version>-cp313-cp313-pyodide_<abi>_wasm32.whl
```

## Local Testing

To test the Pyodide wheel locally, you can set up a Pyodide venv and install the wheel
and run one of the Cantera examples:

```sh
pixi run -e pyodide pyodide venv .pyodide-venv
WHEEL=$(ls build/pyodide_dist/cantera-*.whl)
pixi run -e pyodide .pyodide-venv/bin/pip install "$WHEEL"
pixi run -e pyodide .pyodide-venv/bin/python samples/python/thermo/rankine.py
```

To run other examples, you may need to install additional packages like Matplotlib
and Pandas in this venv.

## CI Job and Artifacts

The `pyodide-wheel` CI job in `.github/workflows/main.yml` does the following:

1. Builds the Pyodide wheel via `scons pyodide-wheel`.
2. Creates a Pyodide virtual environment.
3. Installs the built wheel.
4. Uploads the wheel as a GitHub artifact (`cantera-wheel-pyodide`).
5. Runs the Python unit test suite using Pyodide.

## Common Pitfalls

- **Wrong import target during testing**:
  If `pytest` imports `build/python/cantera` instead of the wheel, clear path overrides
  such as `PYTHONPATH` and verify `import cantera; print(cantera.__file__)`.
- **`No node executable found on the path`**:
  Run `.pyodide-venv/bin/python` through `pixi run -e pyodide ...` so the Node runtime
  from that environment is available.
- **JupyterLite `--piplite-wheels` path confusion**:
  `--piplite-wheels` paths are resolved relative to `--lite-dir`, not the current
  working directory.
