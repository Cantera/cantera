# Regenerating test reference values

The numeric reference values hard-coded in `runtests.jl` are produced with the
**Python** Cantera interface (the reference implementation), so that the Julia
interface is validated against Cantera's own results.

They were generated with Cantera 3.2.0 and `gri30.yaml`. To regenerate:

```python
import cantera as ct

gas = ct.Solution("gri30.yaml")
gas.TPX = 1200.0, ct.one_atm, "CH4:1.0, O2:2.0, N2:7.52"
print("density        ", repr(gas.density))
print("cp_mass        ", repr(gas.cp_mass))
print("cv_mass        ", repr(gas.cv_mass))
print("enthalpy_mass  ", repr(gas.enthalpy_mass))
print("entropy_mass   ", repr(gas.entropy_mass))
print("mean_mw        ", repr(gas.mean_molecular_weight))
print("viscosity      ", repr(gas.viscosity))
print("conductivity   ", repr(gas.thermal_conductivity))
print("wdot(CH4)      ", repr(gas.net_production_rates[gas.species_index("CH4")]))
print("kf[0]          ", repr(gas.forward_rate_constants[0]))
print("eq[0]          ", gas.reaction(0).equation)

# constant-volume H2/air ignition reference
g2 = ct.Solution("gri30.yaml")
g2.TPX = 1000.0, ct.one_atm, "H2:2, O2:1, N2:4"
r = ct.IdealGasReactor(g2, clone=False)
net = ct.ReactorNet([r])
net.advance(1e-3)
print("reactor T@1ms  ", repr(r.T))
```

Notes:

* Python species indices are 0-based; the Julia interface is 1-based, so
  `gas.species_index("OH") == 4` in Python corresponds to `species_index(gas,
  "OH") == 5` in Julia.
* Tolerances in `runtests.jl` are chosen to accommodate benign floating-point
  differences between the C++ core called directly (Julia) and through Python.
