"""
Isentropic, adiabatic flow example - calculate area ratio vs. Mach number curve

Requires: cantera >= 2.5.0, matplotlib >= 2.0
"""

import cantera as ct
import math
import numpy as np


def soundspeed(gas):
    """The speed of sound. Assumes an ideal gas."""

    gamma = gas.cp / gas.cv
    return math.sqrt(gamma * ct.gas_constant
                     * gas.T / gas.mean_molecular_weight)


def isentropic(gas=None):
    """
    In this example, the area ratio vs. Mach number curve is computed. If a gas
    object is supplied, it will be used for the calculations, with the
    stagnation state given by the input gas state. Otherwise, the calculations
    will be done for a 10:1 hydrogen/nitrogen mixture with stagnation T0 = 1200
    K, P0 = 10 atm.
    """
    if gas is None:
        gas = ct.Solution('gri30.yaml')
        gas.TPX = 1200.0, 10.0*ct.one_atm, 'H2:1,N2:0.1'

    # get the stagnation state parameters
    s0 = gas.s
    h0 = gas.h
    p0 = gas.P

    mdot = 1  # arbitrary
    amin = 1.e14

    data = np.zeros((200, 4))

    # compute values for a range of pressure ratios
    for r in range(200):

        p = p0*(r+1)/201.0
        # set the state using (p,s0)
        gas.SP = s0, p

        v2 = 2.0*(h0 - gas.h)      # h + V^2/2 = h0
        v = math.sqrt(v2)
        area = mdot/(gas.density*v)    # rho*v*A = constant
        amin = min(amin, area)
        data[r, :] = [area, v/soundspeed(gas), gas.T, p/p0]

    data[:, 0] /= amin

    return data


if __name__ == "__main__":
    print(__doc__)
    data = isentropic()
    try:
        import matplotlib.pyplot as plt
        plt.plot(data[:, 1], data[:, 0])
        plt.ylabel('Area Ratio')
        plt.xlabel('Mach Number')
        plt.title('Isentropic Flow: Area Ratio vs. Mach Number')
        plt.show()

    except ImportError:
        print('area ratio,   Mach number,   temperature,   pressure ratio')
        print(data)
