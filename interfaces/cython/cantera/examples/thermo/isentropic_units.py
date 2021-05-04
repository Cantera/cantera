"""
Isentropic, adiabatic flow example - calculate area ratio vs. Mach number curve

Requires: cantera >= 2.6.0
"""

import cantera.units as ct
import numpy as np
ct.units.default_format = ".2F~P"
label_string = "area ratio\tMach number\ttemperature\tpressure ratio"
output_string = "{0:.2E~P}\t{1}            {2}\t{3:.2E~P}"


def soundspeed(gas):
    """The speed of sound. Assumes an ideal gas."""

    gamma = gas.cp / gas.cv
    return np.sqrt(gamma * ct.units.molar_gas_constant
                     * gas.T / gas.mean_molecular_weight).to("m/s")


def isentropic(gas=None):
    """
    In this example, the area ratio vs. Mach number curve is computed. If a gas
    object is supplied, it will be used for the calculations, with the
    stagnation state given by the input gas state. Otherwise, the calculations
    will be done for a 10:1 hydrogen/nitrogen mixture with stagnation T0 = 1700.33
    degrees Fahrenheit, P0 = 10 atm.
    """
    if gas is None:
        gas = ct.Solution('gri30.yaml')
        gas.TPX = 2160 * ct.units.degR, 10.0 * ct.units.atm, 'H2:1,N2:0.1'

    # get the stagnation state parameters
    s0 = gas.s
    h0 = gas.h
    p0 = gas.P

    mdot = 1 * ct.units.kg / ct.units.s  # arbitrary
    amin = 1.e14 * ct.units.m**2

    data = []

    # compute values for a range of pressure ratios
    p_range = np.logspace(-3, 0, 200) * p0
    for p in p_range:

        # set the state using (p,s0)
        gas.SP = s0, p

        v = np.sqrt(2.0*(h0 - gas.h)).to("m/s")     # h + V^2/2 = h0
        area = mdot/(gas.density*v)    # rho*v*A = constant
        amin = min(amin, area)
        data.append([area, v/soundspeed(gas), gas.T, p/p0])

    return data, amin


if __name__ == "__main__":
    print(__doc__)
    data, amin = isentropic()
    print(label_string)
    for row in data:
        print(output_string.format(row[0] / amin, row[1], row[2], row[3]))
