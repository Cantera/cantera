# Cantera.jl

A Julia interface to [Cantera](https://cantera.org) for chemical kinetics,
thermodynamics, and transport. It links Cantera's native `libcantera` directly
through the generated CLib API — it is **not** a reimplementation of Cantera and
has **no** Python dependency.

!!! warning "Experimental"
    This interface targets Cantera's experimental CLib backend. The API may
    change between releases.

## Installation

Requires Julia 1.9+ and a compiled `libcantera`. Point the package at the
library and the mechanism data:

```bash
export CANTERA_LIBRARY_PATH=/path/to/cantera/lib   # or use a conda env
export CANTERA_DATA=/path/to/cantera/data          # for gri30.yaml, etc.
```

Generate the low-level CLib bindings from the built `cantera_clib` headers
(again whenever those headers change), then instantiate the environment:

```bash
julia interfaces/julia/generate/generate_bindings.jl
julia --project=interfaces/julia -e 'using Pkg; Pkg.instantiate()'
```

## Quickstart

```julia
using Cantera

gas = Solution("gri30.yaml")
set_TPX!(gas, 1000.0, one_atm, "H2:2, O2:1, N2:4")
net = ReactorNet(IdealGasReactor(gas))
advance!(net, 1e-3)
temperature(gas)          # ignited temperature [K]
```

See the [`examples/`](https://github.com/Cantera/cantera/tree/main/interfaces/julia/examples)
directory for more runnable usage, and the [API Reference](reference.md) for the
full list of documented functions.
