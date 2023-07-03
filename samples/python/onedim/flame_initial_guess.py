"""
A freely-propagating, premixed methane-air flame.
Examples of saving and loading a flame and restarting
with different initial guesses.

Requires: cantera >= 3.0
Keywords: combustion, 1D flow, flame speed, premixed flame, saving output
"""
import sys
from pathlib import Path
import cantera as ct
try:
    import pandas as pd
except ImportError:
    pd = None


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
    """Print a short description of the flame, with a few properties."""
    print(f"\nFlame speed                    = {flame.velocity[0] * 100 :.2f} cm/s")
    print(f"Maximum temperature            = {flame.T.max() :.0f} K")
    # Find the location of the peak OH mole fraction
    oh_index = flame.gas.species_index("OH")
    grid_index = flame.X[oh_index].argmax()
    oh_peak = flame.grid[grid_index]
    print(f"Peak OH mole fraction location = {oh_peak * 100 :.2f} cm")
    print(f"Solved with {flame.grid.size} grid points\n")


describe(f)

# Save the flame in a few different formats

output_path = Path() / "flame_initial_guess_data"
output_path.mkdir(parents=True, exist_ok=True)

print("Save YAML")
yaml_filepath = output_path / "flame.yaml"
f.save(yaml_filepath, name="solution", description="Initial methane flame")

print("Save CSV")
csv_filepath = output_path / "flame.csv"
f.save(csv_filepath, basis="mole", overwrite=True)

if "native" in ct.hdf_support():
    # HDF is not a required dependency
    hdf_filepath = output_path / "flame.h5"
    hdf_filepath.unlink(missing_ok=True)
    f.save(hdf_filepath, name="freeflame", description="Initial methane flame")
    print("Save HDF\n")
else:
    print(f"Skipping HDF: Cantera compiled without HDF support\n")
    hdf_filepath = None

# Restore the flame from different formats

print("Restore solution from YAML")
gas.TPX = Tin, p, reactants
f2 = ct.FreeFlame(gas, width=width)
f2.restore(yaml_filepath, name="solution")
describe(f2)

if hdf_filepath:
    print("Restore solution from HDF")
    gas.TPX = Tin, p, reactants
    f2 = ct.FreeFlame(gas, width=width)
    f2.restore(hdf_filepath, name="freeflame")
    describe(f2)

# Restore the flame via initial guess

print("Load initial guess from CSV file directly")
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=csv_filepath)
describe(f2)

if hdf_filepath:
    print("Load initial guess from HDF file directly")
    gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
    f2 = ct.FreeFlame(gas, width=width)
    f2.set_initial_guess(data=hdf_filepath, group="freeflame")
    describe(f2)

    print("Load initial guess from HDF file via SolutionArray")
    arr2 = ct.SolutionArray(gas)
    # the flame domain needs to be specified as subgroup
    arr2.restore(hdf_filepath, name="freeflame", sub="flame")
    gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
    f2 = ct.FreeFlame(gas, width=width)
    f2.set_initial_guess(data=arr2)
    describe(f2)

if pd is None:
    # skip remaining examples, as optional dependency 'pandas' is not installed
    print("All done")
    sys.exit()

print("Load initial guess from CSV file via Pandas")
df = pd.read_csv(csv_filepath)
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=df)
describe(f2)

print("Load initial guess from CSV file via Pandas and SolutionArray")
df = pd.read_csv(csv_filepath)
arr2 = ct.SolutionArray(gas)
arr2.from_pandas(df)
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_initial_guess(data=arr2)
describe(f2)

# Restart flame simulations with modified initial guesses

print("Load initial guess from CSV file via Pandas, with modifications.\n")
df = pd.read_csv(csv_filepath)
print("Modify the Pandas dataframe, removing half the grid points")
df_pruned = df[::2]  # remove half of the grid points
gas.TPX = Tin, p, reactants  # set the gas T back to the inlet before making new flame
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We wouldn't expect the flame solutions to be exactly the same
describe(f2)

print(
    "Modify the Pandas dataframe, removing half the grid points "
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
    "Modify the Pandas dataframe, removing half the grid points, "
    "and raise the T by 50 K"
)
df_pruned = df.iloc[::2]
# set the gas T back to the (new) inlet before making new flame
gas.TPX = (Tin + 50, p, reactants)
f2 = ct.FreeFlame(gas, width=width)
f2.set_refine_criteria(**refine_criteria)
f2.set_initial_guess(data=df_pruned)
f2.solve()
# We expect these flames to be different because we raised the temperature.
describe(f2)

print("All done")
