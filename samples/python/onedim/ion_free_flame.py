"""
A freely-propagating, premixed methane-air flat flame with charged species.

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, burner-stabilized flame, plasma, premixed flame
"""

from pathlib import Path
import cantera as ct


if "native" in ct.hdf_support():
    output = Path() / "ion_free_flame.h5"
else:
    output = Path() / "ion_free_flame.yaml"
output.unlink(missing_ok=True)

# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:1, O2:2, N2:7.52'  # premixed gas composition
width = 0.05  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# Solution object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('gri30_ion.yaml')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.IonFreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.05, curve=0.1)
f.show_solution()

# stage one
f.solve(loglevel=loglevel, auto=True)

# stage two
f.solve(loglevel=loglevel, stage=2, enable_energy=True)
f.save(output, name="ion", description="solution with ionized gas transport")

f.show_solution()
print(f"mixture-averaged flamespeed = {f.velocity[0]:7f} m/s")

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('ion_free_flame.csv', quiet=False)
