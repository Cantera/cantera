#
# an ideal Rankine cycle
#
from Cantera import *
from Cantera.pureFluids import Water

w = Water()

# start with saturated liquid water at 300 K
w.setTemperature(300.0)
w.setState_satLiquid()
h1 = w.enthalpy_mass()
p1 = w.pressure()

# pump it isentropically to 10 MPa
w.setState_SP(w.entropy_mass(), 1.0e7)
h2 = w.enthalpy_mass()

pump_work = h2 - h1

# heat at constant pressure to 1500 K
w.setState_TP(1500.0, w.pressure())
h3 = w.enthalpy_mass()

heat_in = h3 - h2

# expand isentropically back to 300 K
w.setState_SP(w.entropy_mass(), p1)
h4 = w.enthalpy_mass()

work_out = h3 - h4
heat_out = h4 - h1

efficiency = (work_out - pump_work)/heat_in

print 'efficiency = ',efficiency

