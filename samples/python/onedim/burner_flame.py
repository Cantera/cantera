"""
Burner-stabilized flame
=======================

A burner-stabilized lean premixed hydrogen-oxygen flame at low pressure.

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, premixed flame, saving output,
          multicomponent transport
"""

from pathlib import Path
import matplotlib.pyplot as plt
import cantera as ct

p = 0.05 * ct.one_atm
tburner = 373.0
mdot = 0.06
reactants = 'H2:1.5, O2:1, AR:7'  # premixed gas composition
width = 0.5  # m
loglevel = 1  # amount of diagnostic output (0 to 5)

gas = ct.Solution('h2o2.yaml')
gas.TPX = tburner, p, reactants

f = ct.BurnerFlame(gas, width=width)
f.burner.mdot = mdot
f.set_refine_criteria(ratio=3.0, slope=0.05, curve=0.1)
f.show()

f.transport_model = 'mixture-averaged'
f.solve(loglevel, auto=True)

if "native" in ct.hdf_support():
    output = Path() / "burner_flame.h5"
else:
    output = Path() / "burner_flame.yaml"
output.unlink(missing_ok=True)

f.save(output, name="mix", description="solution with mixture-averaged transport")

f.transport_model = 'multicomponent'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show()
f.save(output, name="multi", description="solution with multicomponent transport")

f.save('burner_flame.csv', basis="mole", overwrite=True)

# %%
# Temperature and Heat Release Rate
# ---------------------------------
fig, ax1 = plt.subplots()

ax1.plot(f.grid, f.heat_release_rate / 1e6, color='C4')
ax1.set_ylabel('heat release rate [MW/m³]', color='C4')
ax1.set_xlim(0, 0.2)
ax1.set(xlabel='distance from burner [m]')

ax2 = ax1.twinx()
ax2.plot(f.grid, f.T, color='C3')
ax2.set_ylabel('temperature [K]', color='C3')
plt.show()

# %%
# Major Species Profiles
# ----------------------
fig, ax = plt.subplots()
major = ('O2', 'H2', 'H2O')
states = f.to_array()
ax.plot(states.grid, states(*major).X, label=major)
ax.set(xlabel='distance from burner [m]', ylabel='mole fractions')
ax.set_xlim(0, 0.2)
ax.legend()
plt.show()

# %%
# Minor Species Profiles
# ----------------------
fig, ax = plt.subplots()
minor = ('OH', 'H', 'O')

ax.plot(states.grid, states(*minor).X, label=minor, linestyle='--')
ax.set(xlabel='distance from burner [m]', ylabel='mole fractions', )
ax.set_xlim(0, 0.2)
ax.legend()
plt.show()
