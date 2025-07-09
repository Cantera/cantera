"""
Critical state properties
=========================

Print the critical state properties for the fluids for which Cantera has
built-in liquid/vapor equations of state.

.. tags:: Python, thermodynamics, multiphase, non-ideal fluid
"""

import cantera as ct
import matplotlib.pyplot as plt

# %%
# Create `PureFluid` objects:
fluids = {'water': ct.Water(),
          'nitrogen': ct.Nitrogen(),
          'methane': ct.Methane(),
          'hydrogen': ct.Hydrogen(),
          'oxygen': ct.Oxygen(),
          'carbon dioxide': ct.CarbonDioxide(),
          'heptane': ct.Heptane(),
          'HFC-134a': ct.Hfc134a()
          }

# %%
# Plot critical properties and print tabulated values:
fig, ax = plt.subplots()

print('Critical State Properties')
print(f"{'Fluid':^16s}   {'Tc [K]':^7s}   {'Pc [Pa]':^10s}   {'Zc':^7s}")
print(f"{'-'*16}   {'-'*7}   {'-'*10}   {'-'*7}")
for name in fluids:
    f = fluids[name]
    tc = f.critical_temperature
    pc = f.critical_pressure
    rc = f.critical_density
    mw = f.mean_molecular_weight
    zc = pc * mw / (rc * ct.gas_constant * tc)
    ax.plot(tc, pc, 'o')
    ax.annotate(name, (tc, pc), (4, 2), textcoords='offset points', size=9)
    print(f'{name:16s}   {tc:7.2f}   {pc:10.4g}   {zc:7.4f}')

ax.set(xlabel='Critical Temperature [K]', ylabel='Critical Pressure [Pa]')
ax.set(xlim=(0, 750))
ax.grid(True)
plt.show()
