#
# an Rankine cycle
#
from Cantera import *
from Cantera.liquidvapor import Water

# parameters
eta_pump = 0.6     # pump isentropic efficiency
et_turbine = 0.8   # turbine isentropic efficiency
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
    fluid.setState_SP(s0, pfinal)
    h1s = fluid.enthalpy_mass()
    isentropic_work = h1s - h0
    actual_work = isentropic_work / eta
    h1 = h0 + actual_work
    fluid.setState_HP(h1, pfinal)
    return actual_work

def expand(fluid, pfinal, eta):
    """Adiabatically expand a fluid to pressure pfinal, using
    a turbine with isentropic efficiency eta."""
    h0 = fluid.enthalpy_mass()
    s0 = fluid.entropy_mass()
    fluid.setState_SP(s0, pfinal)
    h1s = fluid.enthalpy_mass()
    isentropic_work = h0 - h1s
    actual_work = isentropic_work * eta
    h1 = h0 - actual_work
    fluid.setState_HP(h1, pfinal)
    return actual_work

###############################################################


# create an object representing water
w = Water()

# start with saturated liquid water at 300 K
w.setTemperature(300.0)
w.setState_Tsat(0.0)
hf = w.enthalpy_mass()
print w
w.setState_Tsat(1.0)
hv = w.enthalpy_mass()
print hv - hf

print w

# pump it adiabatically to pmax
pump_work = pump(w, pmax, eta_pump)
print pump_work

# heat it at constant pressure until it reaches the
# saturated vapor state at this pressure
#w.setState_Psat(1.0)
#print w

w.setTemperature(273.16)
w.setState_Tsat(0.0)
h0 = w.enthalpy_mass()
for t in [300.0, 350.0, 400.0, 450.0, 500.0]:
    w.setTemperature(t)
    w.setState_Tsat(0.0)
    hf = w.enthalpy_mass()
    w.setState_Tsat(1.0)
    hv = w.enthalpy_mass()
    print t, 0.001*(hf - h0), 0.001*(hv - h0), 0.001*(hv - hf)

for t in [750.0, 800.0, 850.0, 1150.0]:
    w.setState_TP(t, 2.0e4)
    print t, w.enthalpy_mass() - h0
    
