########################################################################
####  NOTE: with the changes made to the ChemEquil solver in version 1.7,
####  it now converges, and the Multiphase solver is no longer invoked
#### in this demo
########################################################################

# Equilibrium of a (nearly) stoichiometric hydrogen/oxygen mixture at
# fixed temperature.

# Cantera has 3 different equilibrium solvers. The 'ChemEquil' solver
# uses the element potential method for homogeneous equilibrium in gas
# mixtures. It is fast, but sometimes doesn't converge. The
# 'MultiPhaseEquil' solver uses the VCS algorithm (Gibbs
# minimization), which is slower but more robust. As the name
# suggests, it can also handle multiple phases. Here we'll solve a
# problem for which the ChemEquil solver fails, but the
# MultiPhaseEquil solver has no problem.

from Cantera import *

# create an object representing the gas phase
gas = importPhase("h2o2.cti")

temp = 400.0

# make the composition very close to stoichiometric
comp = "H2:1.00000001, O2:0.5"

# set the initial state
gas.set(T = temp, P = OneAtm, X = comp)

# equilibrate the gas holding T and P fixed. First try the default
# (ChemEquil) solver...  (This will fail, throwing an exception that
# will be caught in the 'except' block, where we will try the other
# solver.)

####################################################################
# Note: We are setting solver = 0 here to demonstrate the difference
# between the two solvers. If you do not set 'solver', or set it to a
# negative value, then ChemEquil will be tried first, and if it fails
# the MultiPhaseEquil solver will be tried. In most cases this will
# give the best results.
####################################################################

try:
    gas.equilibrate("TP", solver = 0)       # use the ChemEquil (0) solver
except:
    print "ChemEquil solver failed! Try the MultiPhaseEquil solver..."

    # Try again. Reset the gas to the initial state
    gas.set(T = temp, P = OneAtm, X = comp)

    # The MultiPhaseEquil solver is used to equilibrate 'Mixture'
    # objects, since these may have more than one phase. Here we'll
    # create a Mixture object containing only the gas. Some other
    # useful parameters are rtol (relative error tolerance, default =
    # 1.0e-9), max_steps (default = 1000), loglevel (default = 0).
    mix = Mixture([(gas,1.0)])
    mix.equilibrate("TP", loglevel=4)

    # Note: another way to do this is:
    # gas.equilibrate("TP", solver = 1, loglevel = 4)

# print a summary of the results
print gas



# To check that this is an equilibrium state, verify that the chemical
# potentials may be computed by summing the element potentials for each atom.
# (The element potentials are the chemical potentials of the atomic vapors.)

mu_H2, mu_OH, mu_H2O, mu_O2, lambda_H, lambda_O = gas.chemPotentials(
    ["H2", "OH", "H2O", "O2", "H", "O"])


print
print "    Comparison between Chem potentials and element potentials:"
print
s_mu_H2       = "%11.4e" % mu_H2
s_lam_mu_H2   = "%11.4e" % (2.0*lambda_H)
print "mu_H2   : ", s_mu_H2, ",    ", s_lam_mu_H2


s_mu_O2       = "%11.4e" % mu_O2
s_lam_mu_O2   = "%11.4e" % (2.0*lambda_O)
print "mu_O2   : ", s_mu_O2, ",    ", s_lam_mu_O2


s_mu_OH       = "%11.4e" % mu_OH
s_lam_mu_OH   = "%11.4e" % (lambda_H + lambda_O)
print "mu_OH   : ", s_mu_OH, ",    ", s_lam_mu_OH


s_mu_H2O      = "%11.4e" % mu_H2O
s_lam_mu_H2O  = "%11.4e" % (2.0 * lambda_H + lambda_O)
print "mu_H2O  : ", s_mu_H2O, ",    ", s_lam_mu_H2O
