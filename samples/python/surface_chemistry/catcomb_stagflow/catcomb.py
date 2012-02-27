# CATCOMB  -- Catalytic combustion of methane on platinum.
#
# This script solves a catalytic combustion problem. A stagnation flow
# is set up, with a gas inlet 10 cm from a platinum surface at 900
# K. The lean, premixed methane/air mixture enters at ~ 6 cm/s (0.06
# kg/m2/s), and burns catalytically on the platinum surface. Gas-phase
# chemistry is included too, and has some effect very near the
# surface.
#
# The catalytic combustion mechanism is from Deutschman et al., 26th
# Symp. (Intl.) on Combustion,1996 pp. 1747-1754
#
#

from Cantera import *
from Cantera.OneD import *
#from Cantera.OneD.StagnationFlow import StagnationFlow
import math

###############################################################
#
#  Parameter values are collected here to make it easier to modify
#  them

p          =   OneAtm              # pressure
tinlet     =   300.0               # inlet temperature
tsurf      =   900.0               # surface temperature
mdot       =   0.06                # kg/m^2/s

transport  =  'Mix'                # transport model


# We will solve first for a hydrogen/air case to
# use as the initial estimate for the methane/air case

# composition of the inlet premixed gas for the hydrogen/air case
comp1       =  'H2:0.05, O2:0.21, N2:0.78, AR:0.01'

# composition of the inlet premixed gas for the methane/air case
comp2       =  'CH4:0.095, O2:0.21, N2:0.78, AR:0.01'

# the initial grid, in meters. The inlet/surface separation is 10 cm.
initial_grid = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1]  # m


# numerical parameters
tol_ss    = [1.0e-5, 1.0e-9]       # [rtol, atol] for steady-state problem
tol_ts    = [1.0e-4, 1.0e-9]       # [rtol, atol] for time stepping

loglevel  = 1                      # amount of diagnostic output
                                   # (0 to 5)

refine_grid = 1                    # 1 to enable refinement, 0 to
                                   # disable

################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic,
# and transport properties
#
# The gas phase will be taken from the definition of phase 'gas' in
# input file 'ptcombust.cti,' which is a stripped-down version of
# GRI-Mech 3.0.
gas = importPhase('ptcombust.cti','gas')
gas.set(T = tinlet, P = p, X = comp1)


################ create the interface object ##################
#
# This object will be used to evaluate all surface chemical production
# rates. It will be created from the interface definition 'Pt_surf'
# in input file 'ptcombust.cti,' which implements the reaction
# mechanism of Deutschmann et al., 1995 for catalytic combustion on
# platinum.
#
surf_phase = importInterface('ptcombust.cti','Pt_surf', [gas])
surf_phase.setTemperature(tsurf)


# integrate the coverage equations in time for 1 s, holding the gas
# composition fixed to generate a good starting estimate for the
# coverages.
surf_phase.advanceCoverages(1.0)

# create the object that simulates the stagnation flow, and specify an
# initial grid
sim = StagnationFlow(gas = gas, surfchem = surf_phase,
                     grid = initial_grid)

# Objects of class StagnationFlow have members that represent the gas inlet ('inlet') and the surface ('surface'). Set some parameters of these objects.
sim.inlet.set(mdot = mdot, T = tinlet, X = comp1)
sim.surface.set(T = tsurf)

# Set error tolerances
sim.set(tol = tol_ss, tol_time = tol_ts)

# Method 'init' must be called before beginning a simulation
sim.init()

# Show the initial solution estimate
sim.showSolution()



# Solving problems with stiff chemistry coulpled to flow can require
# a sequential approach where solutions are first obtained for
# simpler problems and used as the initial guess for more difficult
# problems.

# start with the energy equation on (default is 'off')
sim.set(energy = 'on')

# disable the surface coverage equations, and turn off all gas and
# surface chemistry.
sim.surface.setCoverageEqs('off')
surf_phase.setMultiplier(0.0);
gas.setMultiplier(0.0);

# solve the problem, refining the grid if needed, to determine the
# non-reacting velocity and temperature distributions
sim.solve(loglevel, refine_grid)

# now turn on the surface coverage equations, and turn the
# chemistry on slowly
sim.surface.setCoverageEqs('on')
for iter in range(6):
    mult = math.pow(10.0,(iter - 5));
    surf_phase.setMultiplier(mult);
    gas.setMultiplier(mult);
    print 'Multiplier = ',mult
    sim.solve(loglevel, refine_grid);

# At this point, we should have the solution for the hydrogen/air
# problem.
sim.showSolution()

# Now switch the inlet to the methane/air composition.
sim.inlet.set(X = comp2)

# set more stringent grid refinement criteria
sim.setRefineCriteria(100.0, 0.15, 0.2, 0.0)

# solve the problem for the final time
sim.solve(loglevel, refine_grid)

# show the solution
sim.showSolution()

# save the solution in XML format. The 'restore' method can be used to restart
# a simulation from a solution stored in this form.
sim.save("catcomb.xml", "soln1")

# save selected solution components in a CSV file for plotting in
# Excel or MATLAB.

# These methods return arrays containing the values at all grid points
z = sim.flow.grid()
T = sim.T()
u = sim.u()
V = sim.V()

f = open('catcomb.csv','w')
writeCSV(f, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)', 'rho (kg/m3']
         + list(gas.speciesNames()))
for n in range(sim.flow.nPoints()):
    sim.setGasState(n)
    writeCSV(f, [z[n], u[n], V[n], T[n], gas.density()]
             +list(gas.moleFractions()))

# write the surface coverages to the CSV file
cov = sim.coverages()
names = surf_phase.speciesNames()
for n in range(len(names)):
    writeCSV(f, [names[n], cov[n]])

f.close()

print 'solution saved to catcomb.csv'

# show some statistics
sim.showStats()
