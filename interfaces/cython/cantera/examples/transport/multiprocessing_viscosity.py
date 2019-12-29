"""
This example demonstrates how Cantera can be used with the 'multiprocessing'
module.

Because Cantera Python objects are built on top of C++ objects which cannot be
passed between Python processes, it is necessary to set up the computation so
that each process has its own copy of the relevant Cantera objects. One way to
do this is by storing the objects in (module) global variables, which are
initialized once per worker process.

Requires: cantera >= 2.5.0
"""

import multiprocessing
import numpy as np
import cantera as ct
import itertools
from time import time

# Global storage for Cantera Solution objects
gases = {}


def init_process(mech):
    """
    This function is called once for each process in the Pool. We use it to
    initialize any Cantera objects we need to use.
    """
    gases[mech] = ct.Solution(mech)
    gases[mech].transport_model = 'Multi'


def get_thermal_conductivity(args):
    # Pool.imap only permits a single argument, so we pack all of the needed
    # arguments into the tuple 'args'
    mech, T, P, X = args
    gas = gases[mech]
    gas.TPX = T, P, X
    return gas.thermal_conductivity


def get_viscosity(args):
    # Pool.imap only permits a single argument, so we pack all of the needed
    # arguments into the tuple 'args'
    mech, T, P, X = args
    gas = gases[mech]
    gas.TPX = T, P, X
    return gas.viscosity


def parallel(mech, predicate, nProcs, nTemps):
    """
    Call the function ``predicate`` on ``nProcs`` processors for ``nTemps``
    different temperatures.
    """
    P = ct.one_atm
    X = 'CH4:1.0, O2:1.0, N2:3.76'
    pool = multiprocessing.Pool(processes=nProcs,
                                initializer=init_process,
                                initargs=(mech,))

    y = pool.map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(300, 900, nTemps),
                     itertools.repeat(P),
                     itertools.repeat(X)))
    return y


def serial(mech, predicate, nTemps):
    P = ct.one_atm
    X = 'CH4:1.0, O2:1.0, N2:3.76'
    init_process(mech)
    y = list(map(predicate,
                 zip(itertools.repeat(mech),
                     np.linspace(300, 900, nTemps),
                     itertools.repeat(P),
                     itertools.repeat(X))))
    return y


if __name__ == '__main__':
    nPoints = 5000
    nProcs = 4

    # For functions where the work done in each subprocess is substantial,
    # significant speedup can be obtained using the multiprocessing module.
    print('Thermal conductivity')
    t1 = time()
    parallel('gri30.yaml', get_thermal_conductivity, nProcs, nPoints)
    t2 = time()
    print('Parallel: {0:.3f} seconds'.format(t2-t1))

    t1 = time()
    serial('gri30.yaml', get_thermal_conductivity, nPoints)
    t2 = time()
    print('Serial: {0:.3f} seconds'.format(t2-t1))

    # On the other hand, if the work done per call to the predicate function is
    # small, there may be no advantage to using multiprocessing.
    print('\nViscosity')
    t1 = time()
    parallel('gri30.yaml', get_viscosity, nProcs, nPoints)
    t2 = time()
    print('Parallel: {0:.3f} seconds'.format(t2-t1))

    t1 = time()
    serial('gri30.yaml', get_viscosity, nPoints)
    t2 = time()
    print('Serial: {0:.3f} seconds'.format(t2-t1))
