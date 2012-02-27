# This script is used to test handling of non-integral product
# stoichiometric coefficients. See file frac.cti for more information.

print 'testing handling of non-integral stoichiometric coefficients...'

from Cantera import *
gas = importPhase('frac.cti')
gas.set(T = 2000, P = OneAtm, X = 'H2O:1.0, OH:0.1, H:0.2, O2:0.3, H2:0.4')

ih2, ih, io, io2, ioh, ih2o = gas.speciesIndex(['H2','H','O','O2','OH','H2O'])

# forward rates of progress
fwd_rop = gas.fwdRatesOfProgress()

# species creation and destruction rates
cdot = gas.creationRates()
ddot = gas.destructionRates()

nsp = gas.nSpecies()
nr = gas.nReactions()

# print the reaction equations
print 'Reaction Equations:'
for i in range(nr):
    print gas.reactionEqn(i)
print

# print the creation rates, and check that the creation rates have the
# correct relationship to the reaction rates of progress
print 'Creation Rates: '
for k in range(nsp-1):
    print '%12s  %10.4e  %10.4e ' % (gas.speciesName(k),
                                     cdot[k], cdot[k]/fwd_rop[0])
print '%12s  %10.4e  %10.4e ' % (gas.speciesName(ih2o),
                                 cdot[ih2o], cdot[ih2o]/fwd_rop[1])

# print the destruction rates, and check that the destruction rates have the
# correct relationship to the reaction rates of progress
print '\nDestruction Rates:'
for k in range(nsp-1):
    print '%12s  %10.4e  %10.4e ' % (gas.speciesName(k),
                                     ddot[k], ddot[k]/fwd_rop[1])
print '%12s  %10.4e  %10.4e ' % (gas.speciesName(ih2o),
                                 ddot[ih2o], ddot[ih2o]/fwd_rop[0])
print

# print the arrays of reactant and product stoichiometric coefficients

x = gas.moleFractions()
c = gas.molarDensity() * x

# rxn 2 orders from frac.cti
order_H2 = 0.8
order_OH = 2.0
order_O2 = 1.0

kf = gas.fwdRateConstants()
print '\nForward rate constants:'
print kf

cproduct = pow(c[ih2],order_H2) * pow(c[ioh], order_OH) * pow(c[io2], order_O2)
print '\nFwd rate of progress, kf*concentration product, difference:'
r1 = fwd_rop[1]
r2 = cproduct*kf[1]
diff12 = (r1 - r2)/(r1 + r2)
if (abs(diff12) < 1.0E-10) :
    diff12 = 0.0
print r1, r2, diff12
print
print 'Reactant stoichiometric coefficients:'
print gas.reactantStoichCoeffs()
print 'Product stoichiometric coefficients:'
print gas.productStoichCoeffs()
