"""
Adiabatic flame temperature and equilibrium composition for a
fuel/air mixture as a function of equivalence ratio,
including formation of solid carbon.
"""

from Cantera import *
import sys

##############################################################
#
# Edit these parameters to change the initial temperature, the
# pressure, and the phases in the mixture
#
###############################################################

temp = 300.0
pres = 101325.0

# phases
gas  = importPhase('gri30.cti')
carbon = importPhase('graphite.cti')

# the phases that will be included in the calculation, and their
# initial moles
mix_phases = [ (gas, 1.0), (carbon, 0.0) ]

# gaseous fuel species
fuel_species = 'CH4'

# air composition
air_N2_O2_ratio = 3.76

# equivalence ratio range
phi_min = 0.3
phi_max = 3.5
npoints = 50

##################################################


mix = Mixture(mix_phases)
nsp = mix.nSpecies()

# create some arrays to hold the data
phi = zeros(npoints,'d')
tad = zeros(npoints,'d')
xeq = zeros([nsp,npoints],'d')

# find fuel, nitrogen, and oxygen indices
ifuel, io2, in2 = gas.speciesIndex([fuel_species, 'O2', 'N2'])
if ifuel < 0:
    raise "fuel species "+fuel_species+" not present!"

if gas.nAtoms(fuel_species,'O') > 0 or  gas.nAtoms(fuel_species,'N') > 0:
    raise "Error: only hydrocarbon fuels are supported."

stoich_O2 = gas.nAtoms(fuel_species,'C') + 0.25*gas.nAtoms(fuel_species,'H')


for i in range(npoints):
    phi[i] = phi_min + (phi_max - phi_min)*i/(npoints - 1)
    x = zeros(nsp,'d')
    x[ifuel] = phi[i]
    x[io2] = stoich_O2
    x[in2] = stoich_O2*air_N2_O2_ratio

    # set the gas state
    gas.set(T = temp, P = pres, X = x)

    # create a mixture of 1 mole of gas, and 0 moles of solid carbon.
    mix = Mixture(mix_phases)
    mix.setTemperature(temp)
    mix.setPressure(pres)

    # equilibrate the mixture adiabatically at constant P
    #
    #    mix.equilibrate('HP', maxsteps = 1000,
    #                err = 1.0e-6, maxiter = 200, loglevel=0)
    mix.vcs_equilibrate('HP', maxsteps = 1000,
                        rtol = 1.0e-6, maxiter = 200, loglevel=0)

    tad[i] = mix.temperature();
    print 'At phi = %12.4g, Tad = %12.4g'  % (phi[i],tad[i])
    xeq[:,i] = mix.speciesMoles()


# write output CSV file for importing into Excel

csvfile = 'adiabatic.csv'
f = open(csvfile,'w')
writeCSV(f,['phi','T (K)']+mix.speciesNames())
for n in range(npoints):
    writeCSV(f,[phi[n], tad[n]]+list(xeq[:,n]))
f.close()
print 'output written to '+csvfile


# make plots
if '--plot' in sys.argv:
    import plotting
    plotting.plotEquilData(mix, phi, tad, xeq)
