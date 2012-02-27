from Cantera import *
import math

def equilSoundSpeeds(gas, rtol = 1.0e-6, maxiter = 5000):

    """Returns a tuple containing the equilibrium and frozen sound
    speeds for a gas with an equilibrium composition.  The gas is
    first set to an equilibrium state at the temperature and pressure
    of the gas, since otherwise the equilibrium sound speed is not
    defined.

    """

    # set the gas to equilibrium at its current T and P
    gas.equilibrate('TP', rtol = rtol, maxiter = maxiter)

    # save properties
    s0 = gas.entropy_mass()
    p0 = gas.pressure()
    r0 = gas.density()

    # perturb the pressure
    p1 = p0*1.0001

    # set the gas to a state with the same entropy and composition but
    # the perturbed pressure
    gas.set(S = s0, P = p1)

    # save the density for this case for the frozen sound speed
    rho_frozen = gas.density()

    # now equilibrate the gas holding S and P constant
    gas.equilibrate("SP", loglevel=0, rtol = rtol, maxiter = maxiter) # , rtol = 1.0e-3, maxsteps=10000)

    r1 = gas.density()

    # equilibrium sound speed
    aequil = math.sqrt((p1 - p0)/(r1 - r0));

    # frozen sound speed
    afrozen = math.sqrt((p1 - p0)/(rho_frozen - r0));

    # compute the frozen sound speed using the ideal gas
    # expression as a check
    cp = gas.cp_mass()
    cv = gas.cv_mass()
    gamma = cp/cv
    afrozen2 = math.sqrt(gamma*GasConstant*gas.temperature()
                     /gas.meanMolecularWeight())
    return (aequil, afrozen, afrozen2)


# test program
if __name__ == "__main__":

    gas = GRI30()
    gas.set(X = 'CH4:1.00, O2:2.0, N2:7.52')
    for n in range(27):
        temp = 300.0 + n*100.0
        gas.set(T = temp, P = OneAtm)
        print temp, equilSoundSpeeds(gas)
