########################################################
#
#   A burner-stabilized hydrogen/oxygen flame
#
########################################################

# note that SI units (m, kg, J, kmol) are used, not cgs units.
import os
from Cantera import units
from Cantera.flame import *

gas = IdealGasMix(src = 'h2o2.xml', transport='Mix')

# create a burner-stabilized flame in the domain z = 0 to z = 20 cm,
# define the fuel to be pure hydrogen, and the oxidizer to be
# oxygen diluted in argon.

flame = BurnerFlame(
    domain = (0, 0.4),   
    fuel = 'H2:1',
    oxidizer = 'O2:1, AR:7',
    gas = gas,
    grid = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2, 0.49, 0.5]
    )

# Set some parameters.
#     mdot     --    mass flow rate in kg/m^2/s
#     T0       --    burner temperature
#     pressure --    P in pascals
#     tol      --    (relative, absolute)
#     timesteps --   ( [sequence of number of steps], initial step size )
#     refine   --    (max size ratio between adj cells, slope parameter,
#                     curvature parameter)
#     jac_age  --    (steady age, transient age)

flame.set(mdot            = 0.04,
          equiv_ratio     = 0.9,
          T_burner        = 373.0,
          pressure        = 0.05 * units.atm,
          tol             = (1.e-5, 1.e-12),
          timesteps       = ([1,2,5,10,20], 1.e-5),
          refine          = (2.0, 0.8, 0.9),
          jac_age         = (20, 10),
          )


# if you want to start from a previously saved solution, uncomment
# this line and modify as necessary
#flame.restore(src = 'h2o2_flame1.xml', solution = 'energy_1')


# turn the energy equation off (default)
flame.set(energy = 'off')

# solve the flame, with output level 1
flame.solve(1)

# save the solution
flame.save('no_energy','solution with the energy equation disabled',
           'h2o2_flame1.xml')

# turn the energy equation on, and change the grid refinement parameters
flame.set(energy = 'on', refine = (2.0, 0.05, 0.1))

# solve it again
flame.solve(1)

# save it to the same file, but with a different solution id.
flame.save('energy','solution with the energy equation enabled',
           'h2o2_flame1.xml')

# write plot files
flame.plot(plotfile = 'flame1.dat', title = 'H2/O2 flame', fmt = 'TECPLOT')
flame.plot(plotfile = 'flame1.csv', title = 'H2/O2 flame', fmt = 'EXCEL')
print '  TECPLOT file flame1.dat and Excel CSV file flame1.csv written'
print '  Directory: '+os.getcwd()

# show statistics -- number of Jacobians, etc.
flame.showStatistics()
