"""
A freely-propagating, premixed methane-air flame.
Examples of saving and loading a flame and restarting
with different initial guesses.

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, flame speed, premixed flame, saving output
"""
import os
import cantera as ct
import pandas as pd
import numpy as np

data_directory = "flame_initial_guess_data"
os.makedirs(data_directory, exist_ok=True)


# Simulation parameters
p = ct.one_atm  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = "CH4:0.45, O2:1.0, N2:3.76"

width = 0.03  # m

# Solution object used to compute mixture properties
gas = ct.Solution("gri30.yaml")
gas.TPX = Tin, p, reactants

# Flame object
f = ct.FreeFlame(gas, width=width)
# These refine criteria are for a relatively quick convergence for the example
# and may not be suitable for a publication-quality flame speed calculation
refine_criteria = {"ratio": 3, "slope": 0.1, "curve": 0.2}
f.set_refine_criteria(**refine_criteria)
f.solve(loglevel=1, auto=True)

def describe(flame):
    "Print a short description of the flame, with a few properties."
    print(f"Flame speed                    = {flame.velocity[0] * 100 :.2f} cm/s")
    print(f"Maximum temperature            = {flame.T.max() :.0f} K")
    # Find the location of the peak OH mole fraction
    oh_index = flame.gas.species_index("OH")
    grid_index = flame.X[oh_index].argmax()
    oh_peak = flame.grid[grid_index]
    print(f"Peak OH mole fraction location = {oh_peak * 100 :.2f} cm")
    print(f"Solved with {flame.grid.size} grid points")


describe(f)

# Save the flame in a few different formats
print("Save YAML")
yaml_filepath = os.path.join(data_directory, "flame.yaml")
f.save(yaml_filepath, name="solution", description="Initial methane flame")

print("Save CSV")
csv_filepath = os.path.join(data_directory, "flame.csv")
f.write_csv(csv_filepath)

print("Save HDF")
hdf_filepath = os.path.join(data_directory, "flame.h5")
f.write_hdf(
    hdf_filepath,
    group="flame",
    mode="w",
    quiet=False,
    description=("Initial methane flame"),
)


print("\nRestore solution from YAML")
gas.TPX = Tin, p, reactants
f2 = ct.FreeFlame(gas, width=width)
f2.restore(yaml_filepath, name="solution", loglevel=0)
describe(f2)

print("\nRestore solution from HDF")
gas.TPX = Tin, p, reactants
f2 = ct.FreeFlame(gas, width=width)
f2.read_hdf(hdf_filepath, group="flame")
describe(f2)

print("\nLoad initial guess from CSV file directly")
csv_filepath = os.path.join(data_directory, "flame.csv")
f.write_csv(csv_filepath)
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=csv_filepath)
describe(f2)

print("\nLoad initial guess from CSV file via Pandas")
csv_filepath = os.path.join(data_directory, "flame.csv")
f.write_csv(csv_filepath)
df = pd.read_csv(csv_filepath)
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=df)
describe(f2)


print("\nLoad initial guess from CSV file via Pandas and SolutionArray")
# In Cantera v2.6.0 passing set_initial_guess a Pandas DataFrame crashed so you
# must use this work-around to create a SolutionArray from your DataFrame
# and then pass the SolutionArray to the set_initial_guess method.
df = pd.read_csv(csv_filepath)
arr2 = ct.SolutionArray(gas)
arr2.from_pandas(df)
gas.TPX = Tin, p, reactants
# set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=arr2)
describe(f2)


print("\nLoad initial guess from CSV file via Pandas, with modifications.")
df = pd.read_csv(csv_filepath)
print("\nModify the Pandas dataframe, removing half the grid points")
df_pruned = df[::2]  # remove half of the grid points
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We wouldn't expect the flame solutions to be exactly the same
describe(f2)

print(
    "\nModify the Pandas dataframe, removing half the grid points "
    "and all but the first 20 species"
)
df_pruned = df.iloc[::2, :24]
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We wouldn't expect the flame solutions to be exactly the same
describe(f2)

print(
    "\nModify the Pandas dataframe, removing half the grid points, "
    "and raise the T by 50 K"
)
df_pruned = df.iloc[::2]
# set the gas T back to the (new) inlet before making new flame
gas.TPX = (Tin + 50, p, reactants)
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We expect thees flames to be different because we raised the temperature.
describe(f2)

print("\nAll done")
