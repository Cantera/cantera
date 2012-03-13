#
# FLAME1 - A burner-stabilized flat flame
#
#    This script simulates a burner-stablized lean hydrogen-oxygen flame
#    at low pressure.
#
from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.BurnerFlame import BurnerFlame

################################################################
#
# parameter values
#
p          =   0.05*OneAtm          # pressure
tburner    =   373.0                # burner temperature
mdot       =   0.06                 # kg/m^2/s

rxnmech    =  'h2o2.cti'            # reaction mechanism file
mix        =  'ohmech'              # gas mixture model
comp       =  'H2:1.8, O2:1, AR:7'  # premixed gas composition

# The solution domain is chosen to be 50 cm, and a point very near the
# downstream boundary is added to help with the zero-gradient boundary
# condition at this boundary.
initial_grid = [0.0, 0.02, 0.04, 0.06, 0.08, 0.1,
                0.15, 0.2, 0.4, 0.49, 0.5]  # m

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
gas.set(T = tburner, P = p, X = comp)

f = BurnerFlame(gas = gas, grid = initial_grid)

# set the properties at the burner
f.burner.set(massflux = mdot, mole_fractions = comp, temperature = tburner)

f.set(tol = tol_ss, tol_time = tol_ts)
f.setMaxJacAge(5, 10)
f.set(energy = 'off')
f.init()
f.showSolution()

f.solve(loglevel, refine_grid)

f.setRefineCriteria(ratio = 200.0, slope = 0.05, curve = 0.1)
f.set(energy = 'on')
f.solve(loglevel,refine_grid)

f.save('flame1.xml')
f.showSolution()


# write the velocity, temperature, and mole fractions to a CSV file
z = f.flame.grid()
T = f.T()
u = f.u()
V = f.V()
fcsv = open('flame1.csv','w')
writeCSV(fcsv, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)', 'rho (kg/m3)']
         + list(gas.speciesNames()))
for n in range(f.flame.nPoints()):
    f.setGasState(n)
    writeCSV(fcsv, [z[n], u[n], V[n], T[n], gas.density()]
             +list(gas.moleFractions()))
fcsv.close()

print 'solution saved to flame1.csv'

f.showStats()
