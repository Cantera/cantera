"""
Sound speeds
============

Compute the "equilibrium" and "frozen" sound speeds for a gas

Requires: cantera >= 3.0.0, matplotlib

.. tags:: Python, thermodynamics, equilibrium
"""

import cantera as ct
import math
import numpy as np
import matplotlib.pyplot as plt


def equilSoundSpeeds(gas, rtol=1.0e-8, max_iter=5000):
    """
    Returns a tuple containing the equilibrium and frozen sound speeds for a
    gas with an equilibrium composition.  The gas is first set to an
    equilibrium state at the temperature and pressure of the gas, since
    otherwise the equilibrium sound speed is not defined.
    """

    # set the gas to equilibrium at its current T and P
    gas.equilibrate('TP', rtol=rtol, max_iter=max_iter)

    # save properties
    s0 = gas.s
    p0 = gas.P
    r0 = gas.density

    # perturb the pressure
    p1 = p0*1.0001

    # set the gas to a state with the same entropy and composition but
    # the perturbed pressure
    gas.SP = s0, p1

    # frozen sound speed
    afrozen = math.sqrt((p1 - p0)/(gas.density - r0))

    # now equilibrate the gas holding S and P constant
    gas.equilibrate('SP', rtol=rtol, max_iter=max_iter)

    # equilibrium sound speed
    aequil = math.sqrt((p1 - p0)/(gas.density - r0))

    # check against the built-in sound speed function
    afrozen2 = gas.sound_speed

    return aequil, afrozen, afrozen2

# %%
# Calculate sound speed at a range of temperatures
gas = ct.Solution('gri30_highT.yaml')
gas.X = 'CH4:1.00, O2:2.0, N2:7.52'
T_range = np.arange(300, 5901, 200)
data = []
print('   T     Equilibrium  Frozen manual  Frozen built-in')
print('  [K]       [m/s]         [m/s]           [m/s]')
print('-------  -----------  -------------  ---------------')
for T in T_range:
    gas.TP = T, ct.one_atm
    aequil, afrozen, afrozen2 = equilSoundSpeeds(gas)
    data.append((T, aequil, afrozen, afrozen2))
    print(f'{T:6.1f}  {aequil:12.2f}  {afrozen:13.2f}  {afrozen2:15.2f}')

# %%
# Plot results
data = np.array(data)
fig, ax = plt.subplots()
ax.plot(data[:,0], data[:,1], label='Equilibrium')
ax.plot(data[:,0], data[:,2], label='Frozen (manual)')
ax.plot(data[:,0], data[:,3], linestyle='none', marker='.', label='Frozen (built-in)')
ax.set(xlabel='T [K]', ylabel='Sound Speed [m/s]')
ax.legend()
plt.show()
