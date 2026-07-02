# Cantera.jl

A Julia interface to [Cantera](https://cantera.org), calling `libcantera`
directly through Cantera's generated CLib API.

```julia
using Cantera

gas = Solution("gri30.yaml")
set_TPY!(gas, 1200.0, one_atm, "CH4:1.0, O2:2.0, N2:7.52")
temperature(gas), pressure(gas), net_production_rates(gas)
```

## API reference

```@autodocs
Modules = [Cantera]
```
