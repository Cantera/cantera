# This script exists to show to use fictive species in Cantera.
# The implementation makes it's use similar to AVBP/
# Note: only constant Schmidt number per species is possible


import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
import os
import pytest
os.chdir(os.path.dirname(os.path.abspath(__file__)))


# Fictive species and soot cannot be activated at the same time 
# Test error message for activation of fictive and soot at the same time
gas_test = ct.Solution('./inputs/BISETTI.yaml') 
with pytest.raises(ct.CanteraError, match="Cannot use both fictive species and soot sections at the same time"):
    f_test = ct.FreeFlame(gas_test, width=0.02, sections=2, fictives=1)


loglevel = 0 # For clarity

# Creating gas object
P, T = 1e5, 300
gas = ct.Solution('./inputs/BISETTI.yaml')

phi = 2.0                                       # equivalence ratio
fuel = {'C2H4': 1}                              # Ethylene composition
oxidizer = {'O2': 1, 'N2': 3.76}                # Oxygen composition

print("First flame: a rich premixed ethylene/air flame !")

gas.TP = T, P
gas.set_equivalence_ratio(phi, fuel, oxidizer)

# -------------- First, a flame without fictive species ------------------ #

f = ct.FreeFlame(gas, width=0.02)   # Create the free laminar premixed flame specifying the width of the grid
f.inlet.X = gas.X                   # Inlet condition on mass fraction
f.inlet.T = gas.T                   # Inlet condition on temperature


f.energy_enabled = True
tol_ss = [1.0e-5, 1.0e-9]  # [rtol atol] for steady-state problem
tol_ts = [1.0e-5, 1.0e-9]  # [rtol atol] for time stepping
f.flame.set_steady_tolerances(default=tol_ss)
f.flame.set_transient_tolerances(default=tol_ts)
f.set_max_jac_age(50, 50)
f.set_time_step(1.0e-5, [2, 5, 10, 20, 80])  # s

f.set_refine_criteria(ratio=2.0, slope=0.05, curve=0.05)

print("Solving flame without fictive")
f.solve(loglevel=1, refine_grid='refine')
f.save('./RESULTS/BISETTI.yaml','without_fic',overwrite=True)

# For the purpose of showing the fictive species, we store the source term of the soot precursor.
Omega_fic = np.zeros(f.flame.n_points,'d') 
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    Omega_fic[n] = gas.net_production_rates[gas.species_index('A2')] * gas.molecular_weights[gas.species_index('A2')]
grid = f.flame.grid

# ------------- Second, fictive species are declared and the previosu flame is restored -------------------- #
f = ct.FreeFlame(gas,  width=0.02, fictives=1)

f.restore('./RESULTS/BISETTI.yaml', 'without_fic')

# Fictive species properties are set like this (remember the order !)
# Note: the Schmidt of A2 is not the real one but sounds close..
f.add_fic(
    fictive_schmidt          = [2.0], 
    fictive_fuel_inlet_Y     = [0.0],
    fictive_oxidizer_inlet_Y = [0.0]
)


# Add source term to the ficitve species (here, this is the "real" source term of A2)
f.flame.set_fictive_source_term_profile('Yfic_0', grid, Omega_fic)
print("Solving flame with fictive")
f.solve(loglevel=1, refine_grid='disabled') # Grid refinement must be deactivated !
f.save('./RESULTS/BISETTI.yaml','with_fic',overwrite=True)

# Retrieve data
fictive_data  = f.fic_Y
Y_A2_fic     = fictive_data[0,:]

Y_A2 = np.zeros(f.flame.n_points,'d') 
for n in range(f.flame.n_points):
    f.set_gas_state(n)
    Y_A2[n]=gas.Y[gas.species_index('A2')]
plt.figure()
plt.plot(f.grid, Y_A2_fic, marker="*", label='A2-Fictive')
plt.plot(f.grid, Y_A2, label='A2-Real')
plt.legend()
# plt.savefig('plot_premixed_with_fictive.png', dpi=500)
plt.show()
