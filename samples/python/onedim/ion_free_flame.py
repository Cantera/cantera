"""
Freely-propagating flame with charged species
=============================================

A freely-propagating, premixed methane-air flat flame with charged species.

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, burner-stabilized flame, plasma, premixed flame
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct


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
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.05, curve=0.1)
f.show()

# stage one
f.solve(loglevel=loglevel, auto=True)

# stage two
f.solve(loglevel=loglevel, stage=2)

if "native" in ct.hdf_support():
    output = Path() / "ion_free_flame.h5"
else:
    output = Path() / "ion_free_flame.yaml"
output.unlink(missing_ok=True)

f.save(output, name="ion", description="solution with ionized gas transport")

f.show()
print(f"mixture-averaged flamespeed = {f.velocity[0]:7f} m/s")

# write the velocity, temperature, density, and mole fractions to a CSV file
f.save('ion_free_flame.csv', basis="mole", overwrite=True)

# %%
# Temperature and Heat Release Rate
# ---------------------------------

# Find the region that covers most of the temperature rise
z = 1000 * f.grid  # convert to mm
i_left = np.where(f.T > f.T[0] + 0.01 * (f.T[-1] - f.T[0]))[0][0]
i_right = np.where(f.T > f.T[0] + 0.95 * (f.T[-1] - f.T[0]))[0][0]
z_left = z[i_left]
z_right = z[i_right]
dz = z_right - z_left
z_left -= 0.3 * dz
z_right += 0.3 * dz

fig, ax1 = plt.subplots()
ax1.plot(z, f.heat_release_rate / 1e6, color='C4')
ax1.set_ylabel('heat release rate [MW/m³]', color='C4')
ax1.set(xlabel='flame coordinate [mm]', xlim=[z_left, z_right])

ax2 = ax1.twinx()
ax2.plot(z, f.T, color='C3')
ax2.set_ylabel('temperature [K]', color='C3')
plt.show()

# %%
# Major Species Profiles
# ----------------------
fig, ax = plt.subplots()
major = ('O2', 'CH4', 'H2O', 'CO2')
states = f.to_array()
ax.plot(z, states(*major).X, label=major)
ax.set(xlabel='flame coordinate [mm]', ylabel='mole fractions')
ax.set_xlim(z_left, z_right)
ax.legend()
plt.show()

# %%
# Minor / Ionized Species Profiles
# --------------------------------
fig, ax = plt.subplots()
minor = ('OH', 'H', 'O', 'E', 'HCO+', 'H3O+')

ax.semilogy(z, states(*minor).X, label=minor, linestyle='--')
ax.set(xlabel='flame coordinate [mm]', ylabel='mole fractions')
ax.set_xlim(z_left, z_right)
ax.legend()
plt.show()

#sphinx_gallery_thumbnail_number = -1
