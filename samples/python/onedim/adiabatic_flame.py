"""
Laminar flame speed calculation
===============================

A freely-propagating, premixed hydrogen flat flame with multicomponent
transport properties.

Requires: cantera >= 3.2, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, premixed flame, multicomponent transport,
          saving output
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
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
# Compute diffusive fluxes using a mass fraction-based gradient ("mass")
# or mole fraction-based gradient ("mole", default)
f.flux_gradient_basis = "mass" # only relevant for mixture-averaged model
f.solve(loglevel=loglevel, auto=True)

f.show()
print(f"mixture-averaged flamespeed = {f.velocity[0]:7f} m/s")

# Solve with mixture-averaged transport model and Soret diffusion
f.soret_enabled = True
f.solve(loglevel=loglevel) # don't use 'auto' on subsequent solves

f.show()
print("mixture-averaged flamespeed with Soret diffusion"
      f" = {f.velocity[0]:7f} m/s")

if "native" in ct.hdf_support():
    output = Path() / "adiabatic_flame.h5"
else:
    output = Path() / "adiabatic_flame.yaml"
output.unlink(missing_ok=True)

f.save(output, name="mix", description="solution with mixture-averaged "
                                       "transport and Soret diffusion")

# Solve with multi-component transport properties
# but without Soret diffusion
f.transport_model = 'multicomponent'
f.soret_enabled = False
f.solve(loglevel)
f.show()
print("multicomponent flamespeed without Soret diffusion"
      f" = {f.velocity[0]:7f} m/s")

# Solve with multi-component transport properties and Soret diffusion
f.soret_enabled = True
f.solve(loglevel)
f.show()
print("multicomponent flamespeed with Soret diffusion"
      f" = {f.velocity[0]:7f} m/s")

f.save(output, name="multi", description="solution with multicomponent transport "
                                         "and Soret diffusion")

# write the velocity, temperature, density, and mole fractions to a CSV file
f.save('adiabatic_flame.csv', basis="mole", overwrite=True)


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
major = ('O2', 'H2', 'H2O')
states = f.to_array()
ax.plot(z, states(*major).X, label=major)
ax.set(xlabel='flame coordinate [mm]', ylabel='mole fractions')
ax.set_xlim(z_left, z_right)
ax.legend()
plt.show()

# %%
# Minor Species Profiles
# ----------------------
fig, ax = plt.subplots()
minor = ('OH', 'H', 'O')

ax.plot(z, states(*minor).X, label=minor, linestyle='--')
ax.set(xlabel='flame coordinate [mm]', ylabel='mole fractions')
ax.set_xlim(z_left, z_right)
ax.legend()
plt.show()
