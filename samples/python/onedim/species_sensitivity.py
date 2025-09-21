r"""
Sensitivity analysis
========================================
In this example we simulate an impinging jet flame, premixed ammonia/hydrogen-air flame,
calculating the sensitivity of N2O to the A-factor reaction rate constants.
Requires: cantera >= 3.0.0, pandas
.. tags:: Python, combustion, 1D flow, species reaction, premixed flame,
          sensitivity analysis, plotting
"""
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define the gas mixture
gas = ct.Solution("example_data/nakamura.yaml")  # Use the Nakamura mechanism for ammonia combustion blends
gas.TP = 295, ct.one_atm
air = "O2:0.21,N2:0.79"
gas.set_equivalence_ratio(phi=0.6, fuel="NH3:0.7, H2:0.3", oxidizer=air)
flame = ct.ImpingingJet(gas=gas, width=0.02)
flame.inlet.mdot = 0.255 * gas.density
flame.surface.T = 493.5
flame.set_initial_guess("equil")


# Refine grid to improve accuracy
flame.set_refine_criteria(ratio=3, slope=0.025, curve=0.05)

# Solve the flame
flame.solve(loglevel=1, auto=True)  # Increase loglevel for more output

# Plot temperature profile
plt.figure(figsize=(8, 6))
plt.plot(flame.grid * 1e3, flame.T, label="Flame Temperature", color="red")
plt.xlabel("Distance (mm)")
plt.ylabel("Temperature (K)")
plt.title("Temperature Profile of a Flame")
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend()
plt.tight_layout()
plt.show()

# Create a DataFrame to store sensitivity-analysis data
sens = pd.DataFrame(index=gas.reaction_equations(), columns=["sensitivity"])

# Use the adjoint method to calculate species sensitivities at a set distance in the flame domain
distance = 0.02
species = "N2O"
sens.sensitivity = flame.get_species_reaction_sensitivities(species, distance)

sens = sens.iloc[(-sens['sensitivity'].abs()).argsort()]
fig, ax = plt.subplots()

sens.head(15).plot.barh(ax=ax, legend=None)
ax.invert_yaxis()  # put the largest sensitivity on top
ax.set_title(f"Sensitivities for {species} Using the Nakamura mechanism")
ax.set_xlabel(r"Sensitivity: $\frac{\partial\:\ln X_{N2O}}{\partial\:\ln k}$")
ax.grid(axis='x')
plt.tight_layout()
plt.show()


# Forward sensitivities
dk = 3e-4

# get index in the grid at distance
ix = np.argmin(np.abs(flame.grid - distance))

Su0 = flame.X[gas.species_index(species), ix]
fwd = []
for m in range(flame.gas.n_reactions):
    flame.gas.set_multiplier(1.0)  # reset all multipliers
    flame.gas.set_multiplier(1 + dk, m)  # perturb reaction m
    flame.solve(loglevel=0, refine_grid=False)
    Suplus = flame.X[gas.species_index(species), ix]
    flame.gas.set_multiplier(1 - dk, m)  # perturb reaction m
    flame.solve(loglevel=0, refine_grid=False)
    Suminus = flame.X[gas.species_index(species), ix]
    fwd.append((Suplus - Suminus) / (2 * Su0 * dk))
sens = pd.DataFrame(index=gas.reaction_equations(), columns=["sensitivity with forward"], data=fwd)
sens = sens.iloc[(-sens['sensitivity with forward'].abs()).argsort()]
fig, ax = plt.subplots()

sens.head(15).plot.barh(ax=ax, legend=None)
ax.invert_yaxis()  # put the largest sensitivity on top
ax.set_title(f"Sensitivities for {species} Using Nakamura Mech")
ax.set_xlabel(r"Sensitivity: $\frac{\partial\:\ln X_{i}}{\partial\:\ln k}$")
ax.grid(axis='x')
plt.tight_layout()
plt.show()