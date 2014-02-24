#
# A Rankine vapor power cycle
#

from Cantera import *
from Cantera.liquidvapor import Water


########################################################
#
# parameters
#

eta_pump = 0.6     # pump isentropic efficiency
eta_turbine = 0.8  # turbine isentropic efficiency
pmax = 8.0e5       # maximum pressure


########################################################
#
# some useful functions
#

def pump(fluid, pfinal, eta):
    """Adiabatically pump a fluid to pressure pfinal, using
    a pump with isentropic efficiency eta."""
    h0 = fluid.enthalpy_mass()
    s0 = fluid.entropy_mass()
    fluid.set(S = s0, P = pfinal)
    h1s = fluid.enthalpy_mass()
    isentropic_work = h1s - h0
    actual_work = isentropic_work / eta
    h1 = h0 + actual_work
    fluid.set(H = h1, P = pfinal)
    return actual_work


def expand(fluid, pfinal, eta):
    """Adiabatically expand a fluid to pressure pfinal, using
    a turbine with isentropic efficiency eta."""
    h0 = fluid.enthalpy_mass()
    s0 = fluid.entropy_mass()
    fluid.set(S = s0, P = pfinal)
    h1s = fluid.enthalpy_mass()
    isentropic_work = h0 - h1s
    actual_work = isentropic_work * eta
    h1 = h0 - actual_work
    fluid.set(H = h1, P = pfinal)
    return actual_work


def printState(n, fluid):
    print '\n\n***************** State '+`n`+' ******************\n', fluid


###############################################################


# create an object representing water
w = Water()

# start with saturated liquid water at 300 K
w.set(T = 300.0, Vapor = 0.0)
h1 = w.enthalpy_mass()
p1 = w.pressure()
printState(1,w)

# pump it adiabatically to pmax
pump_work = pump(w, pmax, eta_pump)
h2 = w.enthalpy_mass()
printState(2,w)

# heat it at constant pressure until it reaches the
# saturated vapor state at this pressure
w.set(P = pmax, Vapor = 1.0)
h3 = w.enthalpy_mass()
heat_added = h3 - h2
printState(3,w)

# expand back to p1
turbine_work = expand(w, p1, eta_turbine)
printState(4,w)

# efficiency
eff = (turbine_work - pump_work)/heat_added

print 'efficiency = ',eff
