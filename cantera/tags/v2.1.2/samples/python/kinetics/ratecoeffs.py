# print forward and reverse rate coefficients for all reactions

from Cantera import *


def show_rate_coefficients(mech, doIrrev = 0):

    """ Print the rate coefficients for a reaction mechanism.  The
    rate coefficients include all factors that multiply the reactant
    concentrations in the law of mass action. The rate coefficients
    may depend on temperature and/or pressure, and are evaluated at
    the temperature and pressure of the object 'mech.'

    If doIrrev has a non-zero value, then the reverse rate coefficients
    will be computed for all reactions, including irreversible
    ones. Otherwise, the rate coefficients for irreversible reactions
    will be set to zero.  """

    # get the array of forward rate coefficients
    kf = mech.fwdRateConstants()

    # get the array of reverse rate coefficients
    kr = mech.revRateConstants(doIrreversible = doIrrev)

    nr = mech.nReactions()

    # get the array of reaction equations
    eqs = mech.reactionEqn(range(nr))

    for i in range(nr):
        print '%40s  %12.5g  %12.5g ' % (eqs[i], kf[i], kr[i])

    print 'units: kmol, m, s'



if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        mech = importPhase(sys.argv[1])
    elif len(sys.argv) > 2:
        mech = importPhase(sys.argv[1], sys.argv[2])
    else:
        mech = GRI30()

    show_rate_coefficients(mech)
