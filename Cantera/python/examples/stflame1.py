"""

A hydrogen/oxygen flame stabilized in an axisymmetric stagnation
flow.

"""

from Cantera import units
from Cantera.flame import *

# Import the hydrogen/oxygen reaction mechanism
# The input file is in directory 'data/inputs'.

gas = IdealGasMix('h2o2.cti')


# Create a stagnation-point flame in the domain z = 0 (the inlet) to z
# = 20 cm (the surface).  The fuel stream will be pure hydrogen,
# and the oxidizer stream oxygen diluted in argon.

flame = StagnationFlame(
    domain      = (0, 0.2),   
    fuel        = 'H2:1',
    oxidizer    = 'O2:1, AR:7',
    gas         = gas,
    grid        = [0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.15, 0.2] # initial grid
    )


# Set some parameters.
#     mdot      --    mass flow rate in kg/m^2/s
#     T_burner  --    burner temperature
#     T_surface --    surface temperature
#     pressure  --    P in pascals
#     tol       --    (relative, absolute)
#     timesteps --    ( [sequence of number of steps], initial step size )
#     refine    --    (max size ratio between adj cells, slope parameter,
#                       curvature parameter)
#     jac_age   --    (steady age, transient age)

flame.set(mdot            = 0.1,
          equiv_ratio     = 1.2,
          T_burner        = 373.0,
          T_surface       = 600.0,          
          pressure        = 0.05 * units.atm,
          tol             = (1.e-7, 1.e-9),
          timesteps       = ([1,2,5,10], 1.e-5),
          refine          = (2.0, 0.5, 0.5),
          jac_age         = (20, 10),
          )


# if you want to start from a previously saved solution, uncomment
# this line and modify as necessary
# flame.restore(src = 'h2o2_flame1.xml', solution = 'energy_1')

# turn the energy equation off (default)
flame.set(energy = 'off')
flame.show()

# solve the flame, with output level 1
flame.solve(1)
flame.show()

# save the solution
flame.save('no_energy','solution with the energy equation disabled',
           'h2o2_stflame1.xml')

# turn the energy equation on, and change the grid refinement parameters
flame.set(energy = 'on', refine = (2.0, 0.1, 0.2))

# solve it again
flame.solve(1)

# save it to the same file, but with a different solution id.
flame.save('energy','solution with the energy equation enabled',
           'h2o2_stflame1.xml')

# write a TECPLOT plot file
flame.plot(plotfile = 'stflame1.dat', title = 'H2/O2 flame', fmt = 'TECPLOT')

# write an Excel CSV file
flame.plot(plotfile = 'stflame1.csv', title = 'H2/O2 flame', fmt = 'EXCEL')

print 'TECPLOT file stflame1.dat and Excel CSV file stflame1.csv written'

# show statistics -- number of Jacobians, etc.
flame.showStatistics()

