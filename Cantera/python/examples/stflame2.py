"""
 A stagnation-point flame using GRI-Mech 3.0.

 In this script, a hydrogen / oxygen / argon flame is first solved,
 and then used as the starting estimate for a methane / oxygen / argon
 flame.

"""

from Cantera.flame import *
from Cantera import units

# start with only a hydrogen/oxygen mechanism
gas = IdealGasMix('h2o2.xml', transport='Mix')

flame = StagnationFlame(
    domain = (0, 0.02),
    fuel = 'H2:1',
    oxidizer = 'O2:1, AR:4',
    gas = gas,
    grid = [0.0, 0.0025, 0.005, 0.0075, 0.01, 0.016, 0.02]    
    )

flame.set(mdot            = 0.4,
          equiv_ratio     = 0.9,
          T_burner        = 373.7,
          T_surface       = 600.0,
          pressure        = 1.0 * units.atm,
          tol             = (1.e-7, 1.e-11),
          timesteps       = ([1,2,5,10,20], 1.e-5),
          refine          = (10.0, 1.0, 1.0),
          jac_age         = (50, 50),
          )

# first solve the fixed-temperature problem
flame.set(energy = 'off')
flame.solve(1)

# now enable the energy equation, and specify that rough
# grid refinement should be done
flame.set(energy = 'on', refine = (3.0, 0.9, 0.9))
flame.solve(1)

flame.save('energy','h2/o2 soln with energy equation', 'h2.xml')


#-----------------------------------------------------------

# Now construct the methane flame using GRI-Mech 3.0
gas2 = GRI30(transport = 'Mix')

flame2 = StagnationFlame(
    domain = (0, 0.02),
    fuel = 'CH4:1',
    oxidizer = 'O2:1, AR:4',
    gas = gas2,
    grid = [0.0, 0.0025, 0.005, 0.0075, 0.01, 0.016, 0.02]    
    )

flame2.set(mdot            = 0.8,
          equiv_ratio     = 0.9,
          T_burner        = 373.7,
          T_surface       = 600.0,
          pressure        = 1.0 * units.atm,
          tol             = (1.e-7, 1.e-11),
          timesteps       = ([1,2,5], 1.e-5),
          refine          = (2.0, 0.1, 0.2),
          jac_age         = (50, 50),
          )

# Use the hydrogen results as the starting guess. 
flame2.restore(src = 'h2.xml', solution = 'energy')

flame2.set(energy = 'on')
flame2.solve(1)
flame2.show()

# write plot files
flame2.plot(plotfile = 'stflame2.dat', title = 'methane/air flame',
           fmt = 'TECPLOT')
flame2.plot(plotfile = 'stflame2.csv', fmt = 'EXCEL')
print 'Solution written to TECPLOT file stflame2.dat and Excel CSV file stflame2.csv'
flame2.showStatistics()

