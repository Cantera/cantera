# Cantera.jl

A Julia interface to [Cantera](https://cantera.org) for chemical
kinetics, thermodynamics, and transport. It links Cantera's native `libcantera`
directly through the generated CLib API — it is **not** a reimplementation of
Cantera and has **no** Python dependency.

## Install

Requires Julia 1.9+ and a compiled `libcantera`. Point the package at the
library (a directory or a full path) and the mechanism data:

```bash
export CANTERA_LIBRARY_PATH=/path/to/cantera/lib   # or use a conda env
export CANTERA_DATA=/path/to/cantera/data          # for gri30.yaml, etc.
```

```julia
using Pkg
Pkg.activate("interfaces/julia")
Pkg.instantiate()
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

## Documentation

Build the API reference with [Documenter](https://documenter.juliadocs.org):

```julia
julia --project=docs docs/make.jl
```

## Status

Experimental. This interface targets Cantera's experimental CLib backend, and
APIs may change.
