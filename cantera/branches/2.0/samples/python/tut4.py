#################################################################
print """

   Tutorial 4:   Chemical Equilibrium

"""
#################################################################

# To set a gas mixture to a state of chemical equilibrium, use the
# equilibrate method.
#
from Cantera import *
g = GRI30()
g.set(T = 300.0, P = OneAtm, X = 'CH4:0.95,O2:2,N2:7.52')
g.equilibrate('TP')

# The above statement sets the state of object 'g' to the state of
# chemical equilibrium holding temperature and pressure
# fixed. Alternatively, the specific enthalpy and pressure can be held
# fixed:

g.set(T = 300.0, P = OneAtm, X = 'CH4:0.95,O2:2,N2:7.52')
g.equilibrate('HP')

# Other options are
#     'UV'   fixed specific internal energy and specific volume
#     'SV'   fixed specific entropy and specific volume
#     'SP'   fixed specific entropy and pressure

g.set(T = 300.0, P = OneAtm, X = 'CH4:0.95,O2:2,N2:7.52')
g.equilibrate('UV')
print g

g.set(T = 300.0, P = OneAtm, X = 'CH4:0.95,O2:2,N2:7.52')
g.equilibrate('SV')
print g

g.set(T = 300.0, P = OneAtm, X = 'CH4:0.95,O2:2,N2:7.52')
g.equilibrate('SP')
print g

# How can you tell if 'equilibrate' has correctly found the
# chemical equilibrium state? One way is verify that the net rates of
# progress of all reversible reactions are zero.

# Here is the code to do this:
g.set(T = 300.0, P = OneAtm, X = 'CH4:0.95,O2:2,N2:7.52')
g.equilibrate('HP')

rf = g.fwdRatesOfProgress()
rr = g.revRatesOfProgress()
for i in range(g.nReactions()):
    if g.isReversible(i) and rf[i] <> 0.0:
        print ' %4i  %10.4g  ' % (i, (rf[i] - rr[i])/rf[i])
# If the magnitudes of the numbers in this list are all very small,
# then each reversible reaction is very nearly equilibrated, which
# only occurs if the gas is in chemical equilibrium.

# You might be wondering how 'equilibrate' works. (Then again, you might
# not, in which case you can go on to the next tutorial now.)  Method
# 'equilibrate' invokes Cantera's chemical equilibrium solver, which
# uses an element potential method. The element potential method is
# one of a class of equivalent 'nonstoichiometric' methods that all
# have the characteristic that the problem reduces to solving a set of
# M nonlinear algebraic equations, where M is the number of elements
# (not species). The so-called 'stoichiometric' methods, on the other
# hand, (including Gibbs minimization), require solving K nonlinear
# equations, where K is the number of species (usually K >> M). See
# Smith and Missen, "Chemical Reaction Equilibrium Analysis" for more
# information on the various algorithms and their characteristics.
#
# Cantera uses a damped Newton method to solve these equations, and
# does a few other things to generate a good starting guess and to
# produce a reasonably robust algorithm. If you want to know more
# about the details, look at the on-line documented source code of
# Cantera C++ class 'ChemEquil.h' at http://www.cantera.org.
