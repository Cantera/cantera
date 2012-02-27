#
# FIXED_T_FLAME - A burner-stabilized, premixed methane/air flat flame
#                 with multicomponent transport properties and a specified
#                 temperature profile
#
from Cantera import *
from Cantera.OneD import *
from Cantera.OneD.BurnerFlame import BurnerFlame
from string import atof


# read temperature vs. position data from a file.
# The file is assumed to have one z, T pair per line, separated by a comma.

def getTempData(filename):
    # open the file containing the temperature data for reading
    f = open(filename)
    z = []
    t = []
    lines = f.readlines()

    # check for unix/Windows/Mac line ending problems
    if len(lines) == 1:
        print 'Warning: only one line found. Possible text file line-ending'
        print 'problem?'
        print 'The one line found is: ',lines[0]


    for line in lines:
        if line[0] == '#':   # use '#' as the comment character
            pass
        else:
            try:
                zval, tval = line.split(',')
                z.append(atof(zval))
                t.append(atof(tval))
            except:
                pass
    print 'read',len(z),'temperature values.'
    f.close()

    # convert z values into non-dimensional relative positions.
    n = len(z)
    zmax = z[n-1]
    for i in range(n):
        z[i] = z[i]/zmax

    return [z,t]



################################################################
#
# parameter values
#
p          =   OneAtm               # pressure
tburner    =   373.7                # burner temperature
mdot       =   0.04                 # kg/m^2/s

comp       =  'CH4:0.65, O2:1, N2:3.76'  # premixed gas composition

# The solution domain is chosen to be 1 cm, and a point very near the
# downstream boundary is added to help with the zero-gradient boundary
# condition at this boundary.
initial_grid = [0.0, 0.0025, 0.005, 0.0075, 0.0099, 0.01] # m

tol_ss    = [1.0e-5, 1.0e-9]        # [rtol atol] for steady-state
                                    # problem
tol_ts    = [1.0e-5, 1.0e-4]        # [rtol atol] for time stepping

loglevel  = 1                       # amount of diagnostic output (0
                                    # to 5)

refine_grid = 1                     # 1 to enable refinement, 0 to
                                    # disable


################ create the gas object ########################
#
# This object will be used to evaluate all thermodynamic, kinetic, and
# transport properties. It is created with two transport managers, to
# enable switching from mixture-averaged to multicomponent transport
# on the last solution.
gas = GRI30('Mix')
gas.addTransportModel('Multi')

# set its state to that of the unburned gas at the burner
gas.setState_TPX(tburner, p, comp)

# create the BurnerFlame object.
f = BurnerFlame(gas = gas, grid = initial_grid)

# set the properties at the burner
f.burner.set(massflux = mdot, mole_fractions = comp, temperature = tburner)

# read in the fixed temperature profile
[zloc, tvalues] = getTempData('tdata.dat')

# set the temperature profile to the values read in
f.flame.setFixedTempProfile(zloc, tvalues)

f.set(tol = tol_ss, tol_time = tol_ts)

# show the initial estimate for the solution
f.showSolution()

# don't solve the energy equation
f.set(energy = 'off')

# first solve the flame with mixture-averaged transport properties
f.setRefineCriteria(ratio = 3.0, slope = 0.3, curve = 1)
f.setMaxJacAge(50, 50)
f.setTimeStep(1.0e-5, [1, 2, 5, 10, 20])

f.solve(loglevel, refine_grid)
f.save('ch4_flame_fixed_T.xml','mixav',
       'solution with mixture-averaged transport')

print '\n\n switching to multicomponent transport...\n\n'
gas.switchTransportModel('Multi')
f.flame.setTransportModel(gas)

f.setRefineCriteria(ratio = 3.0, slope = 0.1, curve = 0.2)
f.solve(loglevel, refine_grid)
f.save('ch4_flame_fixed_T.xml','multi',
       'solution with  multicomponent transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
z = f.flame.grid()
T = f.T()
u = f.u()
V = f.V()
fcsv = open('flame_fixed_T.csv','w')
writeCSV(fcsv, ['z (m)', 'u (m/s)', 'V (1/s)', 'T (K)', 'rho (kg/m3)']
         + list(gas.speciesNames()))
for n in range(f.flame.nPoints()):
    f.setGasState(n)
    writeCSV(fcsv, [z[n], u[n], V[n], T[n], gas.density()]
             +list(gas.moleFractions()))
fcsv.close()

print 'solution saved to flame_fixed_T.csv'

f.showStats()
