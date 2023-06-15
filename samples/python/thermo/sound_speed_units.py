"""
Compute the "equilibrium" and "frozen" sound speeds for a gas. Uses the pint library to
include customized units in the calculation.

Requires: Cantera >= 3.0.0, pint
Keywords: thermodynamics, equilibrium, units
"""

import cantera.with_units as ctu
import numpy as np

# This sets the default output format of the units to have 2 significant digits
# and the units are printed with a Unicode font. See:
# https://pint.readthedocs.io/en/stable/formatting.html#unit-format-types
ctu.units.default_format = ".2F~P"

def equilibrium_sound_speeds(gas, rtol=1.0e-6, max_iter=5000):
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
    afrozen = np.sqrt((p1 - p0)/(gas.density - r0)).to("ft/s")

    # now equilibrate the gas holding S and P constant
    gas.equilibrate('SP', rtol=rtol, max_iter=max_iter)

    # equilibrium sound speed
    aequil = np.sqrt((p1 - p0)/(gas.density - r0)).to("ft/s")

    # check against the built-in sound speed function
    afrozen2 = gas.sound_speed.to("ft/s")

    return aequil, afrozen, afrozen2

# test program
if __name__ == "__main__":
    gas = ctu.Solution('gri30.yaml')
    gas.X = 'CH4:1.00, O2:2.0, N2:7.52'
    T_range = np.linspace(80, 4880, 25) * ctu.units.degF
    print("Temperature      Equilibrium Sound Speed     Frozen Sound Speed      Frozen Sound Speed Check")
    for T in T_range:
        gas.TP = T, 1.0 * ctu.units.atm
        print(T.to("degF"), *equilibrium_sound_speeds(gas), sep = "               ")
