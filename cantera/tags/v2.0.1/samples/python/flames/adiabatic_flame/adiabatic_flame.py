#
# ADIABATIC_FLAME - A freely-propagating, premixed methane/air flat
#          flame with multicomponent transport properties
#
from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.FreeFlame import FreeFlame

################################################################
#
# parameter values
#
p          =   OneAtm               # pressure
tin        =   300.0                # unburned gas temperature
mdot       =   0.04                 # kg/m^2/s

comp       =  'CH4:0.45, O2:1, N2:3.76'  # premixed gas composition

initial_grid = [0.0, 0.001, 0.01, 0.02, 0.029, 0.03] # m

tol_ss    = [1.0e-5, 1.0e-9]        # [rtol atol] for steady-state
                                    # problem
tol_ts    = [1.0e-5, 1.0e-9]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0
                                    # to 5)

refine_grid = 1                     # 1 to enable refinement, 0 to
                                    # disable


gas = GRI30('Mix')
gas.addTransportModel('Multi')

# set its state to that of the unburned gas
gas.setState_TPX(tin, p, comp)

f = FreeFlame(gas = gas, grid = initial_grid, tfix = 600.0)

# set the upstream properties
f.inlet.set(mole_fractions = comp, temperature = tin)

f.set(tol = tol_ss, tol_time = tol_ts)
f.showSolution()

f.set(energy = 'off')
f.setRefineCriteria(ratio = 10.0, slope = 1, curve = 1)
f.setMaxJacAge(50, 50)
f.setTimeStep(1.0e-5, [2, 5, 10, 20, 50])

f.solve(loglevel, refine_grid)
f.save('ch4_adiabatic.xml','no_energy',
       'solution with the energy equation disabled')

f.set(energy = 'on')
f.setRefineCriteria(ratio = 3.0, slope = 0.1, curve = 0.2)
f.solve(loglevel, refine_grid)
f.save('ch4_adiabatic.xml','energy',
       'solution with the energy equation enabled')
print 'mixture-averaged flamespeed = ',f.u()[0]

gas.switchTransportModel('Multi')
f.flame.setTransportModel(gas)
f.solve(loglevel, refine_grid)
f.save('ch4_adiabatic.xml','energy_multi',
       'solution with the energy equation enabled and multicomponent transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
z = f.flame.grid()
T = f.T()
u = f.u()
V = f.V()
fcsv = open('adiabatic_flame.csv','w')
writeCSV(fcsv, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)', 'rho (kg/m3)']
         + list(gas.speciesNames()))
for n in range(f.flame.nPoints()):
    f.setGasState(n)
    writeCSV(fcsv, [z[n], u[n], V[n], T[n], gas.density()]
             +list(gas.moleFractions()))
fcsv.close()

print 'solution saved to adiabatic_flame.csv'
print 'multicomponent flamespeed = ',u[0]
f.showStats()
