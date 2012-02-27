from Cantera import *
from Cantera.num import zeros
import math


def soundspeed(gas):
    """The speed of sound. Assumes an ideal gas."""

    #    if gas.isIdealGas():
    gamma = gas.cp_mass()/gas.cv_mass()
    return math.sqrt(gamma * GasConstant
                     * gas.temperature() / gas.meanMolecularWeight() )
    #else:
    #    raise "non-ideal not implemented."



def isentropic(g = None):
    """
    ISENTROPIC  isentropic, adiabatic flow example

    In this example, the area ratio vs. Mach number curve is
    computed. If a gas object is supplied, it will be used for the
    calculations, with the stagnation state given by the input gas
    state. Otherwise, the calculations will be done for a 10:1
    hydrogen/nitrogen mixture with stagnation T0 = 1200 K, P0 = 10
    atm.

    """

    if g == None:
        gas = GRI30()
        gas.set(T = 1200.0,P = 10.0*OneAtm,X = 'H2:1,N2:0.1')
    else:
        gas = g


    # get the stagnation state parameters
    s0 = gas.entropy_mass()
    h0 = gas.enthalpy_mass()
    p0 = gas.pressure()

    mdot = 1  # arbitrary
    amin = 1.e14

    data = zeros((200,4),'d')

    # compute values for a range of pressure ratios
    for r in range(200):
        p = p0*(r+1)/201.0

        # set the state using (p,s0)
        gas.set(S = s0, P = p)

        h = gas.enthalpy_mass()
        rho = gas.density()

        v2 = 2.0*(h0 - h)      #   h + V^2/2 = h0
        v = math.sqrt(v2)
        area = mdot/(rho*v);   #   rho*v*A = constant
        if area < amin: amin = area
        data[r,:] = [area, v/soundspeed(gas), gas.temperature(), p/p0]

    data[:,0] /= amin

    return data



if __name__ == "__main__":
    print isentropic.__doc__
    data = isentropic()
    try:
        from pylab import *
        clf
        plot(data[:,1], data[:,0])
        ylabel('Area Ratio')
        xlabel('Mach Number')
        title('Isentropic Flow: Area Ratio vs. Mach Number')
        show()

    except:
        print 'area ratio,   Mach number,   temperature,   pressure ratio'
        print data
