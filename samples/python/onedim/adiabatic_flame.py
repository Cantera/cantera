"""
A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, premixed flame, multicomponent transport,
          saving output
"""

from pathlib import Path
import cantera as ct


# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'H2:1.1, O2:1, AR:5'  # premixed gas composition
width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# Solution object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('h2o2.yaml')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.show()

# Solve with mixture-averaged transport model
f.transport_model = 'mixture-averaged'
f.solve(loglevel=loglevel, auto=True)

if "native" in ct.hdf_support():
    output = Path() / "adiabatic_flame.h5"
else:
    output = Path() / "adiabatic_flame.yaml"
output.unlink(missing_ok=True)

# Solve with the energy equation enabled
f.save(output, name="mix", description="solution with mixture-averaged transport")

f.show()
print(f"mixture-averaged flamespeed = {f.velocity[0]:7f} m/s")

# Solve with multi-component transport properties
f.transport_model = 'multicomponent'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show()
print(f"multicomponent flamespeed = {f.velocity[0]:7f} m/s")
f.save(output, name="multi", description="solution with multicomponent transport")

# write the velocity, temperature, density, and mole fractions to a CSV file
f.save('adiabatic_flame.csv', basis="mole", overwrite=True)
