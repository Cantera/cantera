# NPFLAME1 - A nonpremixed counterflow flame.
#
#    This script computes an atmospheric-pressure ethane/air
#    counterflow flame using GRI-Mech 3.0.
#    Run time on a Mac G4: ~ 5 minutes
#
from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.CounterFlame import CounterFlame
from Cantera.num import array

##################################################################
# parameter values
#
# These are grouped here to simplify changing flame conditions

p          =   OneAtm               # pressure
tin_f      =   300.0                # fuel inlet temperature
tin_o      =   300.0                # oxidizer inlet temperature
mdot_o     =   0.72                 # kg/m^2/s
mdot_f     =   0.24                 # kg/m^2/s

comp_o       =  'O2:0.21, N2:0.78, AR:0.01';   # air composition
comp_f       =  'C2H6:1';                      # fuel composition

# distance between inlets is 2 cm; start with an evenly-spaced 6-point
# grid
initial_grid = 0.02*array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0],'d')


tol_ss    = [1.0e-5, 1.0e-9]        # [rtol, atol] for steady-state
                                    # problem
tol_ts    = [1.0e-3, 1.0e-9]        # [rtol, atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0
                                    # to 5)

refine_grid = 1                     # 1 to enable refinement, 0 to
                                    # disable


################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic,
# and transport properties
#

# Here we use GRI-Mech 3.0 with mixture-averaged transport
# properties. To use your own mechanism, use function
# IdealGasMix('mech.cti') to read a mechanism in Cantera format.  If
# you need to convert from Chemkin format, use the ck2cti utility
# program first.
gas = GRI30('Mix')
gas.setPressure(p)

# create an object representing the counterflow flame configuration,
# which consists of a fuel inlet on the left, the flow in the middle,
# and the oxidizer inlet on the right. Class CounterFlame creates this
# configuration.

f = CounterFlame(gas = gas, grid = initial_grid)

# Set the state of the two inlets

f.fuel_inlet.set(massflux = mdot_f,
                 mole_fractions = comp_f,
                 temperature = tin_f)

f.oxidizer_inlet.set(massflux = mdot_o,
                     mole_fractions = comp_o,
                     temperature = tin_o)

# set the error tolerances
f.set(tol = tol_ss, tol_time = tol_ts)

# construct the initial solution estimate. To do so, it is necessary
# to specify the fuel species. If a fuel mixture is being used,
# specify a representative species here for the purpose of
# constructing an initial guess.
f.init(fuel = 'C2H6')

# show the starting estimate
f.showSolution()

# First disable the energy equation and solve the problem without
# refining the grid
f.set(energy = 'off')
f.solve(loglevel, 0)

# Now specify grid refinement criteria, turn on the energy equation,
# and solve the problem again. The ratio parameter controls the
# maximum size ratio between adjacent cells; slope and curve should be
# between 0 and 1 and control adding points in regions of high
# gradients and high curvature, respectively. If prune > 0, points
# will be removed if the relative slope and curvature for all
# components fall below the prune level. Set prune < min(slope,
# curve), or to zero to disable removing grid points.
f.setRefineCriteria(ratio = 200.0, slope = 0.1, curve = 0.2, prune = 0.02)
f.set(energy = 'on')
f.solve(1)

# Save the solution
f.save('npflame1.xml')

# write the velocity, temperature, and mole fractions to a CSV file
z = f.flame.grid()
T = f.T()
u = f.u()
V = f.V()
fcsv = open('npflame1.csv','w')
writeCSV(fcsv, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)']
         + list(gas.speciesNames()))
for n in range(f.flame.nPoints()):
    f.setGasState(n)
    writeCSV(fcsv, [z[n], u[n], V[n], T[n]]+list(gas.moleFractions()))
fcsv.close()

print 'solution saved to npflame1.csv'

f.showSolution()
f.showStats(0)
