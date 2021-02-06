"""
A Rankine vapor power cycle

Requires: Cantera >= 2.5.0
"""

import cantera as ct

# parameters
eta_pump = 0.6     # pump isentropic efficiency
eta_turbine = 0.8  # turbine isentropic efficiency
p_max = 8.0e5       # maximum pressure


def pump(fluid, p_final, eta):
    """Adiabatically pump a fluid to pressure p_final, using
    a pump with isentropic efficiency eta."""
    h0 = fluid.h
    s0 = fluid.s
    fluid.SP = s0, p_final
    h1s = fluid.h
    isentropic_work = h1s - h0
    actual_work = isentropic_work / eta
    h1 = h0 + actual_work
    fluid.HP = h1, p_final
    return actual_work


def expand(fluid, p_final, eta):
    """Adiabatically expand a fluid to pressure p_final, using
    a turbine with isentropic efficiency eta."""
    h0 = fluid.h
    s0 = fluid.s
    fluid.SP =s0, p_final
    h1s = fluid.h
    isentropic_work = h0 - h1s
    actual_work = isentropic_work * eta
    h1 = h0 - actual_work
    fluid.HP = h1, p_final
    return actual_work


def printState(n, fluid):
    print('\n***************** State {0} ******************'.format(n))
    print(fluid.report())


if __name__ == '__main__':
    # create an object representing water
    w = ct.Water()

    # start with saturated liquid water at 300 K
    w.TQ = 300.0, 0.0
    h1 = w.h
    p1 = w.P
    printState(1, w)

    # pump it adiabatically to p_max
    pump_work = pump(w, p_max, eta_pump)
    h2 = w.h
    printState(2, w)

    # heat it at constant pressure until it reaches the saturated vapor state
    # at this pressure
    w.PQ = p_max, 1.0
    h3 = w.h
    heat_added = h3 - h2
    printState(3, w)

    # expand back to p1
    turbine_work = expand(w, p1, eta_turbine)
    printState(4, w)

    # efficiency
    eff = (turbine_work - pump_work)/heat_added

    print('efficiency = ', eff)
