"""
Rankine cycle (with units)
==========================

Calculate the efficiency of a Rankine vapor power cycle using a pure fluid model
for water. Includes the units of quantities in the calculations.

Requires: Cantera >= 3.0.0, pint

.. tags:: Python, thermodynamics, thermodynamic cycle, non-ideal fluid, units
"""

import cantera.with_units as ctu

# %%
# Parameters
# ----------
eta_pump = 0.6 * ctu.units.dimensionless  # pump isentropic efficiency
eta_turbine = 0.8 * ctu.units.dimensionless  # turbine isentropic efficiency
p_max = 116.03 * ctu.units.psi  # maximum pressure

# %%
# Helper Functions
# ----------------
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


def print_state(n, fluid):
    print('\n***************** State {0} ******************'.format(n))
    print(fluid.report())

# %%
# Evaluate Rankine Cycle
# ----------------------

# Create an object representing water:
w = ctu.Water()

# %%
# Start with saturated liquid water at 80.33 degrees Fahrenheit:
w.TQ = ctu.Q_(80.33, "degF"), 0.0 * ctu.units.dimensionless
h1 = w.h
p1 = w.P
print_state(1, w)

# %%
# Pump it adiabatically to ``p_max``:
pump_work = pump(w, p_max, eta_pump)
h2 = w.h
print_state(2, w)

# %%
# Heat it at constant pressure until it reaches the saturated vapor state
# at this pressure:
w.PQ = p_max, 1.0 * ctu.units.dimensionless
h3 = w.h
heat_added = h3 - h2
print_state(3, w)

# %%
# expand back to ``p1``:
turbine_work = expand(w, p1, eta_turbine)
print_state(4, w)

# %%
# Calculate the efficiency:
eff = (turbine_work - pump_work)/heat_added
print('efficiency = ', eff)
