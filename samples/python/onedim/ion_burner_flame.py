"""
A burner-stabilized premixed methane-air flame with charged species.

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, burner-stabilized flame, plasma, premixed flame
"""

from pathlib import Path
import cantera as ct


p = ct.one_atm
tburner = 600.0
reactants = 'CH4:1.0, O2:2.0, N2:7.52'  # premixed gas composition
width = 0.05  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('gri30_ion.yaml')
gas.TPX = tburner, p, reactants
mdot = 0.15 * gas.density

f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show()

f.solve(loglevel, auto=True)
f.solve(loglevel=loglevel, stage=2)

if "native" in ct.hdf_support():
    output = Path() / "ion_burner_flame.h5"
else:
    output = Path() / "ion_burner_flame.yaml"
output.unlink(missing_ok=True)

f.save(output, name="mix", description="solution with mixture-averaged transport")

f.save('ion_burner_flame.csv', basis="mole", overwrite=True)
