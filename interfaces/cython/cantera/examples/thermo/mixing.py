"""
Mixing two streams using `Quantity` objects.

In this example, air and methane are mixed in stoichiometric proportions. This
is a simpler, steady-state version of the example ``reactors/mix1.py``.

Since the goal is to simulate a continuous flow system, the mixing takes place
at constant enthalpy and pressure.

Requires: cantera >= 2.5.0
"""

import cantera as ct

gas = ct.Solution('gri30.yaml')

# Stream A (air)
A = ct.Quantity(gas, constant='HP')
A.TPX = 300.0, ct.one_atm, 'O2:0.21, N2:0.78, AR:0.01'

# Stream B (methane)
B = ct.Quantity(gas, constant='HP')
B.TPX = 300.0, ct.one_atm, 'CH4:1'

# Set the molar flow rates corresponding to stoichiometric reaction,
# CH4 + 2 O2 -> CO2 + 2 H2O
A.moles = 1
nO2 = A.X[A.species_index('O2')]
B.moles = nO2 * 0.5

# Compute the mixed state
M = A + B
print(M.report())

# Show that this state corresponds to stoichiometric combustion
M.equilibrate('TP')
print(M.report())
