# Equilibrium of a (nearly) stoichiometric hydrogen/oxygen mixture at
# fixed temperature.

# Cantera has 2 different equilibrium solvers. The 'ChemEquil' solver
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
try:
    gas.equilibrate("TP")   
except:
    print "ChemEquil solver failed! Try the MultiPhaseEquil solver..."
    
    # Try again. Reset the gas to the initial state
    gas.set(T = temp, P = OneAtm, X = comp)        

    # setting parameter 'solver' to 1 requests that the
    # MultiPhaseEquil solver be used (specifying 0 would cause
    # ChemEquil to be used). Some other useful parameters are rtol
    # (relative error tolerance, default = 1.0e-9), max_steps (default = 1000),
    # loglevel (default = 0).
    gas.equilibrate("TP", solver = 1, rtol = 1.0e-10, loglevel = 4)
                    
# print a summary of the results
print gas



# To check that this is an equilibrium state, verify that the chemical
# potentials may be computed by summing the element potentials for each atom.
# (The element potentials are the chemical potentials of the atomic vapors.)

mu_H2, mu_OH, mu_H2O, mu_O2, lambda_H, lambda_O = gas.chemPotentials(
    ["H2", "OH", "H2O", "O2", "H", "O"])

print mu_H2, 2.0*lambda_H
print mu_O2, 2.0*lambda_O
print mu_OH, lambda_H + lambda_O
print mu_H2O, 2.0*lambda_H + lambda_O

