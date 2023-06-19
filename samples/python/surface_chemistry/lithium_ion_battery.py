"""
This example calculates the cell voltage of a lithium-ion battery at
given temperature, pressure, current, and range of state of charge (SOC).

The thermodynamics are based on a graphite anode and a LiCoO2 cathode,
modeled using the 'BinarySolutionTabulatedThermo' class.
Further required cell parameters are the electrolyte ionic resistance, the
stoichiometry ranges of the active materials (electrode balancing), and the
surface area of the active materials.

The functionality of this example is presented in greater detail in a jupyter
notebook as well as the reference (which also describes the derivation of the
'BinarySolutionTabulatedThermo' class):

Reference:
M. Mayur, S. C. DeCaluwe, B. L. Kee, W. G. Bessler, “Modeling and simulation
of the thermodynamics of lithium-ion battery intercalation materials in the
open-source software Cantera,” Electrochim. Acta 323, 134797 (2019),
https://doi.org/10.1016/j.electacta.2019.134797

Requires: cantera >= 2.6.0, matplotlib >= 2.0
Keywords: surface chemistry, kinetics, electrochemistry, battery, plotting
"""

import cantera as ct
import numpy as np

# Parameters
samples = 101
soc = np.linspace(0., 1., samples)  # [-] Input state of charge (0...1)
current = -1  # [A] Externally-applied current, negative for discharge
T = 293  # [K] Temperature
P = ct.one_atm  # [Pa] Pressure

# Cell properties
input_file = "lithium_ion_battery.yaml"  # Cantera input file name
R_electrolyte = 0.0384  # [Ohm] Electrolyte resistance
area_cathode = 1.1167  # [m^2] Cathode total active material surface area
area_anode = 0.7824  # [m^2] Anode total active material surface area

# Electrode balancing: The "balancing" of the electrodes relates the chemical
# composition (lithium mole fraction in the active materials) to the macroscopic
# cell-level state of charge.
X_Li_anode_0 = 0.01  # [-] anode Li mole fraction at SOC = 0
X_Li_anode_1 = 0.75  # [-] anode Li mole fraction at SOC = 100
X_Li_cathode_0 = 0.99  # [-] cathode Li mole fraction at SOC = 0
X_Li_cathode_1 = 0.49  # [-] cathode Li mole fraction at SOC = 100

# Calculate mole fractions from SOC
X_Li_anode = (X_Li_anode_1 - X_Li_anode_0) * soc + X_Li_anode_0
X_Li_cathode = (X_Li_cathode_0 - X_Li_cathode_1) * (1 - soc) + X_Li_cathode_1

# Import all Cantera phases
anode = ct.Solution(input_file, "anode")
cathode = ct.Solution(input_file, "cathode")
metal = ct.Solution(input_file, "electron")
electrolyte = ct.Solution(input_file, "electrolyte")
anode_int = ct.Interface(
    input_file, "edge_anode_electrolyte", adjacent=[anode, metal, electrolyte])
cathode_int = ct.Interface(
    input_file, "edge_cathode_electrolyte", adjacent=[cathode, metal, electrolyte])

# Set the temperatures and pressures of all phases
for phase in [anode, cathode, metal, electrolyte, anode_int, cathode_int]:
    phase.TP = T, P


# Root finding function
def newton_solve(f, xstart, C=0.0):
    """
    Solve f(x) = C by Newton iteration using initial guess *xstart*
    """
    f0 = f(xstart) - C
    x0 = xstart
    dx = 1.0e-6
    n = 0
    while n < 200:
        ff = f(x0 + dx) - C
        dfdx = (ff - f0) / dx
        step = - f0 / dfdx

        # avoid taking steps too large
        if abs(step) > 0.1:
            step = 0.1 * step / abs(step)

        x0 += step
        emax = 0.00001  # 0.01 mV tolerance
        if abs(f0) < emax and n > 8:
            return x0
        f0 = f(x0) - C
        n += 1
    raise Exception("no root!")


# This function returns the Cantera calculated anode current (in A)
def anode_current(phi_s, phi_l, X_Li_anode):
    """
    Current from the anode as a function of anode potential relative to
    electrolyte.
    """
    # Set the active material mole fraction
    anode.X = {"Li[anode]": X_Li_anode, "V[anode]": 1 - X_Li_anode}

    # Set the electrode and electrolyte potential
    metal.electric_potential = phi_s
    electrolyte.electric_potential = phi_l

    # Get the net reaction rate at the anode-side interface
    # Reaction according to input file:
    # Li+[electrolyte] + V[anode] + electron <=> Li[anode]
    r = anode_int.net_rates_of_progress[0]  # [kmol/m2/s]

    # Calculate the current. Should be negative for cell discharge.
    return r * ct.faraday * area_anode


# This function returns the Cantera calculated cathode current (in A)
def cathode_current(phi_s, phi_l, X_Li_cathode):
    """
    Current to the cathode as a function of cathode potential relative to electrolyte
    """
    # Set the active material mole fractions
    cathode.X = {"Li[cathode]": X_Li_cathode, "V[cathode]": 1 - X_Li_cathode}

    # Set the electrode and electrolyte potential
    metal.electric_potential = phi_s
    electrolyte.electric_potential = phi_l

    # Get the net reaction rate at the cathode-side interface
    # Reaction according to input file:
    # Li+[electrolyte] + V[cathode] + electron <=> Li[cathode]
    r = cathode_int.net_rates_of_progress[0]  # [kmol/m2/s]

    # Calculate the current. Should be negative for cell discharge.
    return - r * ct.faraday * area_cathode


# Calculate cell voltage, separately for each entry of the input vectors
V_cell = np.zeros_like(soc)
phi_l_anode = 0
phi_s_cathode = 0
for i in range(samples):
    # Set anode electrode potential to 0
    phi_s_anode = 0

    # Calculate anode electrolyte potential
    phi_l_anode = newton_solve(
        lambda E: anode_current(phi_s_anode, E, X_Li_anode[i]),
        phi_l_anode, C=current)

    # Calculate cathode electrolyte potential
    phi_l_cathode = phi_l_anode + current * R_electrolyte

    # Calculate cathode electrode potential
    phi_s_cathode = newton_solve(
        lambda E: cathode_current(E, phi_l_cathode, X_Li_cathode[i]),
        phi_s_cathode, C=current)

    # Calculate cell voltage
    V_cell[i] = phi_s_cathode - phi_s_anode

try:
    import matplotlib.pyplot as plt

    # Plot the cell voltage, as a function of the state of charge
    plt.plot(soc * 100, V_cell)
    plt.xlabel("State of charge / %")
    plt.ylabel("Cell voltage / V")
    plt.show()

except ImportError:
    print("Install matplotlib to plot the outputs")
