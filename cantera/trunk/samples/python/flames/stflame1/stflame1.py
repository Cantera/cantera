#
# STFLAME1 - A detached flat flame stabilized at a stagnation point
#

#    This script simulates a lean hydrogen-oxygen flame stabilized in
#    a strained flowfield at an axisymmetric stagnation point on a
#    non-reacting surface. The solution begins with a flame attached
#    to the inlet (burner), and the mass flow rate is progressively
#    increased, causing the flame to detach and move closer to the
#    surface. This example illustrates use of the new 'prune' grid
#    refinement parameter, which allows grid points to be removed if
#    they are no longer required to resolve the solution. This is
#    important here, since the flamefront moves as the mass flowrate
#    is increased. Without using 'prune', a large number of grid
#    points would be concentrated upsteam of the flame, where the
#    flamefront had been previously. (To see this, try setting prune
#    to zero.)

from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.StagnationFlow import StagnationFlow

################################################################
#
# parameter values
#
p          =   0.05*OneAtm          # pressure
tburner    =   373.0                # burner temperature
tsurf      =   600.0

# each mdot value will be solved to convergence, with grid refinement,
# and then that solution will be used for the next mdot
mdot       =   [0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12]                  # kg/m^2/s

rxnmech    =  'h2o2.cti'            # reaction mechanism file
comp       =  'H2:1.8, O2:1, AR:7'  # premixed gas composition

# The solution domain is chosen to be 50 cm, and a point very near the
# downstream boundary is added to help with the zero-gradient boundary
# condition at this boundary.
initial_grid = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1,
                0.15, 0.2]  # m

tol_ss    = [1.0e-5, 1.0e-13]        # [rtol atol] for steady-state
                                    # problem
tol_ts    = [1.0e-4, 1.0e-9]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0
                                    # to 5)

refine_grid = 1                     # 1 to enable refinement, 0 to
                                    # disable
ratio = 5.0
slope = 0.1
curve = 0.2
prune = 0.05



################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic,
# and transport properties
#
gas = IdealGasMix(rxnmech)

# set its state to that of the unburned gas at the burner
gas.setState_TPX(tburner, p, comp)

# Create the stagnation flow object with a non-reactive surface.  (To
# make the surface reactive, supply a surface reaction mechanism.  see
# example catcomb.py for how to do this.)
f = StagnationFlow(gas = gas, grid = initial_grid)

# set the properties at the inlet
f.inlet.set(massflux = mdot[0], mole_fractions = comp, temperature = tburner)

# set the surface state
f.surface.setTemperature(tsurf)

f.set(tol = tol_ss, tol_time = tol_ts)
f.setMaxJacAge(5, 10)
f.set(energy = 'off')
f.init(products = 'equil') # assume adiabatic equilibrium products
f.showSolution()

f.solve(loglevel, refine_grid)

f.setRefineCriteria(ratio = ratio, slope = slope,
                    curve = curve, prune = prune)
f.set(energy = 'on')

m = 0
for md in mdot:
    f.inlet.set(mdot = md)
    f.solve(loglevel,refine_grid)
    m = m + 1
    f.save('stflame1.xml','mdot'+`m`,'mdot = '+`md`+' kg/m2/s')


    # write the velocity, temperature, and mole fractions to a CSV file
    z = f.flow.grid()
    T = f.T()
    u = f.u()
    V = f.V()
    fcsv = open('stflame1_'+`m`+'.csv','w')
    writeCSV(fcsv, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)']
             + list(gas.speciesNames()))
    for n in range(f.flow.nPoints()):
        f.setGasState(n)
        writeCSV(fcsv, [z[n], u[n], V[n], T[n]]+list(gas.moleFractions()))
    fcsv.close()

    print 'solution saved to flame1.csv'

f.showStats()
