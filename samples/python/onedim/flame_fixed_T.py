"""
Burner-stabilized flame with imposed temperature profile
========================================================

A burner-stabilized, premixed methane/air flat flame with multicomponent
transport properties and a specified temperature profile.

Requires: cantera >= 3.0, matplotlib >= 2.0

.. tags:: Python, combustion, 1D flow, burner-stabilized flame, premixed flame, plotting,
          saving output
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import cantera as ct


# %%
# parameter values
p = ct.one_atm  # pressure
tburner = 373.7  # burner temperature
mdot = 0.04  # kg/m^2/s
comp = 'CH4:0.65, O2:1, N2:3.76'  # premixed gas composition

# The solution domain is chosen to be 1 cm
width = 0.01  # m

loglevel = 1  # amount of diagnostic output (0 to 5)
refine_grid = True  # 'True' to enable refinement

# %%
# Create the gas object. This object will be used to evaluate all thermodynamic,
# kinetic, and transport properties
gas = ct.Solution('gri30.yaml')

# set its state to that of the unburned gas at the burner
gas.TPX = tburner, p, comp

# create the BurnerFlame object.
f = ct.BurnerFlame(gas=gas, width=width)

# set the mass flow rate at the burner
f.burner.mdot = mdot

# %%
# Specify temperature vs. position data. For the purposes of this example, the data
# is embedded directly. For more general use, the temperature profile should be stored
# as a separate data file and read in, for example using `np.genfromtxt`.

zloc = np.array([
    0.        , 0.00015625, 0.00023437, 0.00039063, 0.00046875, 0.00050781,
    0.00054688, 0.000625  , 0.00066406, 0.00070312, 0.00074219, 0.00078125,
    0.00082031, 0.00085938, 0.00089844, 0.0009375 , 0.00101563, 0.00105469,
    0.00109375, 0.00113281, 0.00117187, 0.00121094, 0.00125   , 0.00128906,
    0.00132813, 0.00136719, 0.00140625, 0.00144531, 0.00148438, 0.00152344,
    0.0015625 , 0.00160156, 0.00164062, 0.00171875, 0.00175781, 0.00179688,
    0.00183594, 0.001875  , 0.00191406, 0.00195312, 0.00199219, 0.00203125,
    0.00207031, 0.00210938, 0.00214844, 0.0021875 , 0.00222656, 0.00226562,
    0.00230469, 0.00234375, 0.00238281, 0.00242187, 0.00246094, 0.0025    ,
    0.00257813, 0.00265625, 0.00273437, 0.0028125 , 0.00289062, 0.00296875,
    0.00304688, 0.003125  , 0.00328125, 0.0034375 , 0.00359375, 0.00375   ,
    0.00390625, 0.0087    , 0.01
])

tvalues = np.array([
     373.7      ,  465.4070428,  510.4311676,  599.5552837,  643.8342938,
     665.9335545,  688.0122338,  732.1284327,  754.1744755,  776.2170662,
     798.2588757,  820.3020011,  842.348001 ,  864.3979228,  886.4523159,
     908.5112198,  952.6396629,  974.7018199,  996.7515831, 1018.777651 ,
    1040.765863 , 1062.69948  , 1084.558639 , 1106.320078 , 1127.956918 ,
    1149.438472 , 1170.730129 , 1191.793309 , 1212.585506 , 1233.060477 ,
    1253.168589 , 1272.857384 , 1292.072391 , 1328.859767 , 1346.323998 ,
    1363.101361 , 1379.147594 , 1394.425274 , 1408.905834 , 1422.569115 ,
    1435.40408  , 1447.410648 , 1458.597668 , 1468.982722 , 1478.590978 ,
    1487.453914 , 1495.607879 , 1503.092709 , 1509.950449 , 1516.224147 ,
    1521.956853 , 1527.19079  , 1531.966722 , 1536.32348  , 1543.891739 ,
    1550.203579 , 1555.480771 , 1559.908135 , 1563.637879 , 1566.794144 ,
    1569.477867 , 1571.77099  , 1575.385829 , 1578.108169 , 1580.194856 ,
    1581.820666 , 1583.106578 , 1589.51315  , 1589.578955
])

zloc /= max(zloc)

# set the temperature profile to the values read in
f.flame.set_fixed_temp_profile(zloc, tvalues)

# %%
# show the initial estimate for the solution
f.show()

# don't solve the energy equation
f.energy_enabled = False

# first solve the flame with mixture-averaged transport properties
f.transport_model = 'mixture-averaged'
f.set_refine_criteria(ratio=3.0, slope=0.3, curve=1)

f.solve(loglevel, refine_grid)

if "native" in ct.hdf_support():
    output = Path() / "flame_fixed_T.h5"
else:
    output = Path() / "flame_fixed_T.yaml"
output.unlink(missing_ok=True)

f.save(output, name="mix", description="solution with mixture-averaged transport")

print('\n\n switching to multicomponent transport...\n\n')
f.transport_model = 'multicomponent'

f.set_refine_criteria(ratio=3.0, slope=0.1, curve=0.2)
f.solve(loglevel, refine_grid)
f.save(output, name="multi", description="solution with multicomponent transport")

# write the velocity, temperature, density, and mole fractions to a CSV file
f.save('flame_fixed_T.csv', basis="mole", overwrite=True)
f.show_stats()

# %%
# Temperature and Heat Release Rate
# ---------------------------------
fig, ax1 = plt.subplots()

ax1.plot(f.grid, f.heat_release_rate / 1e6, color='C4')
ax1.set_ylabel('heat release rate [MW/m³]', color='C4')
ax1.set_xlim(0, 0.01)
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
ax.set_xlim(0, 0.01)
ax.legend()
plt.show()

# %%
# Minor Species Profiles
# ----------------------
fig, ax = plt.subplots()
minor = ('OH', 'H', 'O')

ax.plot(states.grid, states(*minor).X, label=minor, linestyle='--')
ax.set(xlabel='distance from burner [m]', ylabel='mole fractions', )
ax.set_xlim(0, 0.01)
ax.legend()
plt.show()
