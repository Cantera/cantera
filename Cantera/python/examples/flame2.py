#
# rich methane/air flame
#
import os
from Cantera.flame import *
from Cantera import units
#from Cantera.gases import H_O_AR


gas = GRI30(transport = 'Mix')

flame = BurnerFlame(
    domain = (0, 0.01),
    fuel = 'CH4:1',
    oxidizer = 'O2:1, N2:3.76',
    gas = gas,
    grid = [0.0, 0.0025, 0.005, 0.0075, 0.0099, 0.01]
    )

flame.set(mdot            = 0.04,
          equiv_ratio     = 1.3,
          T_burner        = 373.7,
          pressure        = 1.0 * units.atm,
          tol             = (1.e-4, 1.e-9),
          rtol            = (1.e-5, 1.e-5),
          timesteps       = ([1,2,5,10,20], 1.e-5),
          refine          = (10.0, 1, 1),
          jac_age         = (50, 50),
          nsteps          = [1,2,5,10,20]
          )

flame.set(energy = 'off')
flame.solve(1)
flame.save('no_energy','solution with the energy equation disabled',
           'ch4_flame1.xml')

flame.set(energy = 'on', refine = (3.0, 0.1, 0.2))
flame.solve(1)

flame.save('energy','solution with the energy equation enabled',
           'ch4_flame1.xml')

# write plot files
flame.plot(plotfile = 'flame2.dat', title = 'methane/air flame',
           fmt = 'TECPLOT')
flame.plot(plotfile = 'flame2.csv', fmt = 'EXCEL')
print '  Solution written to TECPLOT file flame2.dat and Excel CSV file flame2.csv'
print '  Directory: '+os.getcwd()
flame.showStatistics()
