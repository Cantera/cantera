#
# A freely-propagating premixed hydrogen/air flame
#
#
from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.FreeFlame import FreeFlame

################################################################
#
# parameter values
#
p          =   OneAtm          # pressure
tin        =   300.0                # unburned gas temperature

rxnmech    =  'ohn.cti'            # reaction mechanism file
mix        =  'gas'              # gas mixture model
comp       =  'H2:2, O2:1, N2:3.76'  # premixed gas composition

# The solution domain is chosen to be 50 cm, and a point very near the
# downstream boundary is added to help with the zero-gradient boundary
# condition at this boundary.
initial_grid = [0.0, 0.001, 0.02, 0.04, 0.07, 0.099, 0.1]  # m

tol_ss    = [1.0e-5, 1.0e-13]        # [rtol atol] for steady-state
                                    # problem
tol_ts    = [1.0e-4, 1.0e-9]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0
                                    # to 5)

refine_grid = 1                     # 1 to enable refinement, 0 to
                                    # disable


################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic,
# and transport properties
#
gas = IdealGasMix(rxnmech, mix)

# set its state to that of the unburned gas at the burner
gas.set(T = tin, P = p, X = comp)

f = FreeFlame(gas = gas, grid = initial_grid)

# set the properties at the inlet
f.inlet.set(mole_fractions = comp, temperature = tin)

f.set(tol = tol_ss, tol_time = tol_ts)
f.setMaxJacAge(5, 10)
f.set(energy = 'off')
#f.init()
f.showSolution()

f.solve(loglevel, refine_grid)

f.setRefineCriteria(ratio = 5.0, slope = 0.05, curve = 0.005, prune = 0.0)
f.set(energy = 'on')
f.solve(loglevel,refine_grid)

f.save('freeflame1.xml')
f.showSolution()


# write the velocity, temperature, and mole fractions to a CSV file
z = f.flame.grid()
T = f.T()
u = f.u()
V = f.V()
fcsv = open('freeflame1.csv','w')
writeCSV(fcsv, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)', 'rho (kg/m3)']
         + list(gas.speciesNames()))
for n in range(f.flame.nPoints()):
    f.setGasState(n)
    writeCSV(fcsv, [z[n], u[n], V[n], T[n], gas.density()]
             +list(gas.moleFractions()))
fcsv.close()

print 'solution saved to freeflame1.csv'
print 'flamespeed = ',u[0],'m/s'
f.showStats()
