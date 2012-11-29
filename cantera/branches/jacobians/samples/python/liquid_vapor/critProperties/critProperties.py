"""
Print the critical state properties for the fluids for which Cantera has
built-in liquid/vapor equations of state.
"""

from Cantera import *
from Cantera.liquidvapor import *

fluids = {'water':Water(),
          'nitrogen':Nitrogen(),
          'methane':Methane(),
          'hydrogen':Hydrogen(),
          'oxygen':Oxygen(),
          'carbondioxide':CarbonDioxide(),
          'heptane':Heptane()
          }

print 'Critical State Properties'
print '%20s  %10s  %10s  %10s' % ('Fluid','Tc [K]', 'Pc [Pa]', 'Zc')
for name in fluids.keys():
    f = fluids[name]
    tc = f.critTemperature()
    pc = f.critPressure()
    rc = f.critDensity()
    mw = f.meanMolecularWeight()
    zc = pc*mw/(rc*GasConstant*tc)
    print '%20s   %10.4g   %10.4G  %10.4G' % (name, tc, pc, zc)
