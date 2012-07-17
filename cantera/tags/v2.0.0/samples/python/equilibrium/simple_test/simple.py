# homogeneous equilibrium of a gas.

from Cantera import *

# create an object representing the gas phase
gas = importPhase("h2o2.cti")

# set the initial state
gas.set(T = 2000.0, P = 0.1*OneAtm, X = "H2:1.0, O2:0.6")

# equilibrate the gas holding T and P fixed
gas.equilibrate("TP")

# print a summary of the results
print gas

# Individual properties can also be retrieved...
x = gas.moleFractions()
y = gas.massFractions()
mu = gas.chemPotentials()
names = gas.speciesNames()
for n in range(gas.nSpecies()):
    print "%20s  %10.4g   %10.4g   %10.4g " % (names[n], x[n], y[n], mu[n])
