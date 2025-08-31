"""
Burner-stabilized flame including ionized species
=================================================

A burner-stabilized premixed methane-air flame with charged species.

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, burner-stabilized flame, plasma, premixed flame
"""

from pathlib import Path
import matplotlib.pyplot as plt
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

f.electric_field_enabled = False
f.solve(loglevel, auto=True)
f.electric_field_enabled = True
f.solve(loglevel=loglevel)

if "native" in ct.hdf_support():
    output = Path() / "ion_burner_flame.h5"
else:
    output = Path() / "ion_burner_flame.yaml"
output.unlink(missing_ok=True)

f.save(output, name="mix", description="solution with mixture-averaged transport")

f.save('ion_burner_flame.csv', basis="mole", overwrite=True)

# %%
# Temperature and Heat Release Rate
# ---------------------------------
fig, ax1 = plt.subplots()

ax1.plot(f.grid * 1000, f.heat_release_rate / 1e6, color='C4')
ax1.set_ylabel('heat release rate [MW/m³]', color='C4')
ax1.set_xlim(0, 3.0)
ax1.set(xlabel='distance from burner [mm]')

ax2 = ax1.twinx()
ax2.plot(f.grid * 1000, f.T, color='C3')
ax2.set_ylabel('temperature [K]', color='C3')
plt.show()

# %%
# Major Species Profiles
# ----------------------
fig, ax = plt.subplots()
major = ('O2', 'CH4', 'H2O', 'CO2')
states = f.to_array()
ax.plot(states.grid * 1000, states(*major).X, label=major)
ax.set(xlabel='distance from burner [mm]', ylabel='mole fractions')
ax.set_xlim(0, 3.0)
ax.legend()
plt.show()

#sphinx_gallery_thumbnail_number = 2

# %%
# Ionized Species Profiles
# ------------------------
fig, ax = plt.subplots()
minor = ('E', 'H3O+', 'HCO+')

ax.semilogy(states.grid * 1000, states(*minor).X, label=minor, linestyle='--')
ax.set(xlabel='distance from burner [mm]', ylabel='mole fractions', )
ax.set_xlim(0, 3.0)
ax.legend()
plt.show()
