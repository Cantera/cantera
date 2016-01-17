"""
Print the critical state properties for the fluids for which Cantera has
built-in liquid/vapor equations of state.
"""

import cantera as ct

fluids = {'water': ct.Water(),
          'nitrogen': ct.Nitrogen(),
          'methane': ct.Methane(),
          'hydrogen': ct.Hydrogen(),
          'oxygen': ct.Oxygen(),
          'carbon dioxide': ct.CarbonDioxide(),
          'heptane': ct.Heptane(),
          'hfc134a': ct.Hfc134a()
          }

print('Critical State Properties')
print('%20s  %10s  %10s  %10s' % ('Fluid','Tc [K]', 'Pc [Pa]', 'Zc'))
for name in fluids:
    f = fluids[name]
    tc = f.critical_temperature
    pc = f.critical_pressure
    rc = f.critical_density
    mw = f.mean_molecular_weight
    zc = pc * mw / (rc * ct.gas_constant * tc)
    print('%20s   %10.4g   %10.4G  %10.4G' % (name, tc, pc, zc))
