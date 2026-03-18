"""
Plasma Reactor
==============
"""

# %%
from pathlib import Path
import cantera as ct

mechanism = Path(__file__).parent / "simple.yaml"

plasma = ct.Solution(
    mechanism,
    "isotropic-electron-energy-plasma",
    transport_model=None,
)

# %%
P = ct.one_atm
T = 300.0
Te = T  # K
X = [0.8, 0.1, 0.1]  # O, O+, e
plasma.Te = Te
plasma.TPX = T, P, f"O: {X[0]}, O+: {X[1]}, e-: {X[2]}"

rho = plasma.density

print(f"Tg: {plasma.T:.2f} K")
print(f"Te: {plasma.Te:.2f} K")
print(f"Total pressure: {plasma.P:.5e} Pa")
print(f"Total mass density: {plasma.density:.5e} kg/m^3")
print(f"Mass fractions: {plasma.Y}")
print("-" * 40)

# %%
P = ct.one_atm
T = 300.0
Te = 30000.0  # K
X = [0.8, 0.1, 0.1]  # O, O+, e
plasma.Te = Te
plasma.TPX = T, P, f"O: {X[0]}, O+: {X[1]}, e-: {X[2]}"

print(f"Tg: {plasma.T:.2f} K")
print(f"Te: {plasma.Te:.2f} K")
print(f"Total pressure: {plasma.P:.5e} Pa")
print(f"Total mass density: {plasma.density:.5e} kg/m^3")
print(f"Mass fractions: {plasma.Y}")
print("-" * 40)

# %%
plasma.Te = Te
plasma.TDX = T, rho, f"O: {X[0]}, O+: {X[1]}, e-: {X[2]}"

print(f"Tg: {plasma.T:.2f} K")
print(f"Te: {plasma.Te:.2f} K")
print(f"Total pressure: {plasma.P:.5e} Pa")
print(f"Total mass density: {plasma.density:.5e} kg/m^3")
print(f"Mass fractions: {plasma.Y}")
print("-" * 40)

# %%
c_p_O = 3 * ct.gas_constant  # J/kmol/K
c_p_O_plus = 8 * ct.gas_constant  # J/kmol/K
c_p_e = 2.5 * ct.gas_constant  # J/kmol/K

cp_mole = X[0] * c_p_O + X[1] * c_p_O_plus + X[2] * c_p_e
print(f"cp_mole: {cp_mole:.2f} J/kmol/K")

print(f"cp_mole (from plasma): {plasma.cp_mole:.2f} J/kmol/K")

print(plasma.cp_mole_e)

# %%
h_O = 3 * ct.gas_constant * T + 3e4 * ct.gas_constant  # J/kmol
h_O_plus = 8 * ct.gas_constant * T + 2e5 * ct.gas_constant  # J/kmol
h_e = 2.5 * ct.gas_constant * Te - 750 * ct.gas_constant  # J/kmol

h_mole = X[0] * h_O + X[1] * h_O_plus + X[2] * h_e
print(f"h_mole: {h_mole:.2f} J/kmol")
print(f"h_mole (from plasma): {plasma.enthalpy_mole:.2f} J/kmol")


# %%

print("Testing Reactor")
reactor = ct.Reactor(plasma)
print("Reactor created successfully")

plasma()
