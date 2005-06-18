# This script is used to test handling of non-integral product
# stoichiometric coefficients. See file frac.cti for more information.

from Cantera import *
gas = importPhase('frac.cti')
gas.set(T = 2000, P = OneAtm, X = 'H2O:1.0, OH:0.1')
fwd_rop = gas.fwdRatesOfProgress()

cdot = gas.creationRates()

nsp = gas.nSpecies()
nr = gas.nReactions()

# print the reaction equations
for i in range(nr):
    print gas.reactionEqn(i)

# print the creation rates, and check that the creation rates have the
# correct relationship to the reaction rates of progress
for k in range(nsp):
    print '%12s  %10.4e  %10.4e ' % (gas.speciesName(k),
                                     cdot[k], cdot[k]/fwd_rop[0])

# print the arrays of reactant and product stoichiometric coefficients
print gas.reactantStoichCoeffs()
print gas.productStoichCoeffs()

