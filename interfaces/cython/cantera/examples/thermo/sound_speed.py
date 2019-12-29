"""
Compute the "equilibrium" and "frozen" sound speeds for a gas

Requires: cantera >= 2.5.0
"""

import cantera as ct
import math


def equilSoundSpeeds(gas, rtol=1.0e-6, maxiter=5000):
    """
    Returns a tuple containing the equilibrium and frozen sound speeds for a
    gas with an equilibrium composition.  The gas is first set to an
    equilibrium state at the temperature and pressure of the gas, since
    otherwise the equilibrium sound speed is not defined.
    """

    # set the gas to equilibrium at its current T and P
    gas.equilibrate('TP', rtol=rtol, maxiter=maxiter)

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
    gas.equilibrate('SP', rtol=rtol, maxiter=maxiter)

    # equilibrium sound speed
    aequil = math.sqrt((p1 - p0)/(gas.density - r0))

    # compute the frozen sound speed using the ideal gas expression as a check
    gamma = gas.cp/gas.cv
    afrozen2 = math.sqrt(gamma * ct.gas_constant * gas.T /
                         gas.mean_molecular_weight)

    return aequil, afrozen, afrozen2


# test program
if __name__ == "__main__":
    gas = ct.Solution('gri30.yaml')
    gas.X = 'CH4:1.00, O2:2.0, N2:7.52'
    for n in range(27):
        T = 300.0 + 100.0 * n
        gas.TP = T, ct.one_atm
        print(T, equilSoundSpeeds(gas))
