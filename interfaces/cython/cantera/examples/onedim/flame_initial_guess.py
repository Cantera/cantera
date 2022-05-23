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

# Try to speed things up by loading the solution
yaml_filepath = os.path.join(data_directory, "flame.yaml")
try:
    f.restore(yaml_filepath, name="solution", loglevel=0)
except:
    print("Couldn't load yaml solution, so starting from scratch")

f.solve(loglevel=1, auto=True)
print("mixture-averaged flame speed = {:7f} m/s\n".format(f.velocity[0]))


def compare_flames(f1, f2):
    """
    Compare two flame solutions f1 and f2 and check that they're the same.

    Expects the grids to be exactly the same, as if you have restored
    the solution.
    """
    assert np.allclose(f1.grid, f2.grid), "Grid is different"
    assert np.allclose(f1.T, f2.T), "Temperature profile is different"
    assert np.allclose(f1.X[f.gas.species_index("CH4")], f2.X[f2.gas.species_index("CH4")] ), "Methane profile is different"

print("Save YAML and restore solution")
yaml_filepath = os.path.join(data_directory, "flame.yaml")
f.save(yaml_filepath, name="solution", description="Initial methane flame")
gas.TPX = Tin, p, reactants
f2 = ct.FreeFlame(gas, width=width)
f2.restore(yaml_filepath, name="solution", loglevel=0)
compare_flames(f, f2)

print("Save HDF and restore solution")
hdf_filepath = os.path.join(data_directory, "flame.h5")
f.write_hdf(hdf_filepath, group="flame", mode="w", quiet=False,
            description=("Initial methane flame"))
gas.TPX = Tin, p, reactants
f2 = ct.FreeFlame(gas, width=width)
f2.read_hdf(hdf_filepath, group="flame")
compare_flames(f, f2)


print("Save CSV then load initial guess via Pandas and SolutionArray")
# In Cantera v2.6.0 passing set_initial_guess a Pandas DataFrame crashes so you
# must use this work-around to create a SolutionArray from your DataFrame
# and then pass the SolutionArray to the set_initial_guess method.
csv_filepath = os.path.join(data_directory, "flame.csv")
f.write_csv(csv_filepath)
df = pd.read_csv(csv_filepath)
arr2 = ct.SolutionArray(gas)
arr2.from_pandas(df)
gas.TPX = Tin, p, reactants # must set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=arr2)
compare_flames(f, f2)

print("Save CSV then load initial guess via Pandas")
csv_filepath = os.path.join(data_directory, "flame.csv")
f.write_csv(csv_filepath)
df = pd.read_csv(csv_filepath)
gas.TPX = Tin, p, reactants # must set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=df)
compare_flames(f, f2)

print("Save CSV and load initial guess directly")
# This one fails
csv_filepath = os.path.join(data_directory, "flame.csv")
f.write_csv(csv_filepath)
gas.TPX = Tin, p, reactants # must set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=csv_filepath)
compare_flames(f, f2)

print("Save HDF and load initial guess directly")
# This one fails
hdf_filepath = os.path.join(data_directory, "flame.h5")
f.write_hdf(hdf_filepath, group="flame", mode="w", quiet=False,
            description=("Initial methane flame"))
gas.TPX = Tin, p, reactants # must set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
try:
    f2.set_initial_guess(data=hdf_filepath)
except KeyError:
    print("As expected, this doesn't work, though the documentation suggests it should")
else:
    print("If this now works, something has been fixed and the example file should be updated")
    compare_flames(f, f2)


print("Save CSV then load initial guess via Pandas, with modifications")
csv_filepath = os.path.join(data_directory, "flame.csv")
f.write_csv(csv_filepath)
df = pd.read_csv(csv_filepath)
print("Modify the Pandas dataframe, removing half the grid points")
df_pruned = df[::2] # remove half of the grid points
gas.TPX = Tin, p, reactants # must set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We wouldn't expect the flame solutions to be exactly the same, so we just check a few other properties
assert np.isclose(f.velocity[0], f2.velocity[0], 1e-2 ), "Flame speeds not within 1%"
assert np.allclose(f.X.max(axis=1), f2.X.max(axis=1), rtol=0.2, atol=1e-20), "Max mole fractions not within 20%"

print("Modify the Pandas dataframe, removing half the grid points and all but the first 20 species")
df_pruned = df.iloc[::2,:24]
gas.TPX = Tin, p, reactants # must set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We wouldn't expect the flame solutions to be exactly the same, so we just check a few other properties
assert np.isclose(f.velocity[0], f2.velocity[0], 1e-2 ), "Flame speeds not within 1%"
assert np.allclose(f.X.max(axis=1), f2.X.max(axis=1), rtol=0.2, atol=1e-20), "Max mole fractions not within 20%"

print("Modify the Pandas dataframe, removing half the grid points, and raise the T by 50 K")
df_pruned = df.iloc[::2]
gas.TPX = Tin + 50, p, reactants # must set the gas T back to the (new) inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We expect thees flames to be different because we raised the temperature.
assert  f2.velocity[0] / f.velocity[0] > 1.1, "Flame speed hasn't increased by more than 10%"

print("All done")
