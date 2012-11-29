""" Solve a steady-state problem by combined damped Newton iteration
 and time integration. Function solve is no longer used, now that the
 functional equivalent has been added to the Cantera C++ kernel.  """

from Cantera import CanteraError
from Cantera.num import array
import math, types

print
"""
module solve is deprecated, and may be removed in a future release. If you
use it and do not want it removed, send an e-mail to cantera-help@caltech.edu.
"""

def solve(sim, loglevel = 0, refine_grid = 1, plotfile = '', savefile = ''):
    """
    Solve a steady-state problem by combined damped Newton iteration
    and time integration.
    """

    new_points = 1

    # get options
    dt = sim.option('timestep')
    ft = sim.option('ftime')

    # sequence of timesteps
    _steps = sim.option('nsteps')
    if type(_steps) == types.IntType: _steps = [_steps]

    len_nsteps = len(_steps)
    dt = sim.option('timestep')

    ll = loglevel
    soln_number = -1
    max_timestep = sim.option('max_timestep')

    sim.collect()

    # loop until refine adds no more points
    while new_points > 0:

        istep = 0
        nsteps = _steps[istep]

        # loop until Newton iteration succeeds
        ok = 0
        while ok == 0:

            # Try to solve the steady-state problem by damped
            # Newton iteration.
            try:
                if loglevel > 0:
                    print 'Attempt Newton solution of ',\
                          'steady-state problem...',
                sim.newton_solve(loglevel-1)

                if loglevel > 0:
                    print 'success.\n\n'
                    print '%'*79+'\n'
                    print 'Problem solved on ',sim.npts,' point grid(s).\n'
                    print '%'*79+'\n'
                ok = 1
                soln_number += 1
                sim.finish()


            except CanteraError:

                # Newton iteration failed.
                if loglevel > 0: print '\n'

                # Take nsteps time steps, starting with step size
                # dt. The final dt may be smaller than the initial
                # value if one or more steps fail.

                if loglevel == 1:
                    print 'Take',nsteps,' timesteps',

                dt = sim.py_timeStep(nsteps,dt,loglevel=ll-1)
                if loglevel == 1: print dt, math.log10(sim.ssnorm())
                istep += 1
                if istep >=  len_nsteps:
                    nsteps = _steps[-1]
                    dt *= 2.0
                else:
                    nsteps = _steps[istep]
                if dt > max_timestep: dt = max_timestep



        # A converged solution was found. Save and/or plot it, then
        # check whether the grid should be refined.

        # Add the solution to the plot file
        if plotfile:
            sim.outputTEC(plotfile,"flame","p"+`sim.npts`,append=soln_number)

        # If a filename has been specified for a save file, add
        # the solution to this file
        if savefile:
            sim.save(savefile, soln_name+'_'+`sim.npts`+'_points')

        if loglevel > 2: sim.show()

        if refine_grid:

            # Call refine to add new points, if needed
            new_points = sim.refine(loglevel = loglevel - 1)

        else:
            new_points = 0
