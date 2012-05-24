/**
 * @file equilibrate.cpp
 * Driver routines for the chemical equilibrium solvers.
 *
 */
#include "cantera/equil/equil.h"
#include "cantera/equil/ChemEquil.h"
#include "cantera/equil/MultiPhaseEquil.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/base/global.h"

namespace Cantera
{

/*
 * Set a multiphase mixture to a state of chemical equilibrium.
 * This is the top-level driver for multiphase equilibrium. It
 * doesn't do much more than call the equilibrate method of class
 * MultiPhase, except that it adds some messages to the logfile,
 * if loglevel is set > 0.
 *
 * @ingroup equil
 */
doublereal equilibrate(MultiPhase& s, const char* XY,
                       doublereal tol, int maxsteps, int maxiter,
                       int loglevel)
{
    if (loglevel > 0) {
        beginLogGroup("equilibrate",loglevel);
        addLogEntry("multiphase equilibrate function");
        beginLogGroup("arguments");
        addLogEntry("XY",XY);
        addLogEntry("tol",tol);
        addLogEntry("maxsteps",maxsteps);
        addLogEntry("maxiter",maxiter);
        addLogEntry("loglevel",loglevel);
        endLogGroup("arguments");
    }
    s.init();
    int ixy = _equilflag(XY);
    if (ixy == TP || ixy == HP || ixy == SP || ixy == TV) {
        try {
            double err = s.equilibrate(ixy, tol, maxsteps, maxiter, loglevel);
            if (loglevel > 0) {
                addLogEntry("Success. Error",err);
                endLogGroup("equilibrate");

            }
            return err;
        } catch (CanteraError& err) {
            err.save();
            if (loglevel > 0) {
                addLogEntry("Failure.",lastErrorMessage());
                endLogGroup("equilibrate");
            }
            throw err;
        }
    } else {
        if (loglevel > 0) {
            addLogEntry("multiphase equilibrium can be done only for TP, HP, SP, or TV");
            endLogGroup("equilibrate");
        }
        throw CanteraError("equilibrate","unsupported option");
        return -1.0;
    }
    return 0.0;
}

/*
 *  Set a single-phase chemical solution to chemical equilibrium.
 *  This is a convenience function that uses one or the other of
 *  the two chemical equilibrium solvers.
 *
 *  @param s The object to set to an equilibrium state
 *
 *  @param XY An integer specifying the two properties to be held
 *  constant.
 *
 *  @param solver The equilibrium solver to use. If solver = 0,
 *  the ChemEquil solver will be used, and if solver = 1, the
 *  MultiPhaseEquil solver will be used (slower than ChemEquil,
 *  but more stable). If solver < 0 (default, then ChemEquil will
 *  be tried first, and if it fails MultiPhaseEquil will be tried.
 *
 *  @param maxsteps The maximum number of steps to take to find
 *  the solution.
 *
 *  @param maxiter For the MultiPhaseEquil solver only, this is
 *  the maximum number of outer temperature or pressure iterations
 *  to take when T and/or P is not held fixed.
 *
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *  = 0 suppresses diagnostics, and increasingly-verbose messages
 *  are written as loglevel increases. The messages are written to
 *  a file in HTML format for viewing in a web browser.
 *  @see HTML_logs
 *
 *  @ingroup equil
 */
int equilibrate(thermo_t& s, const char* XY, int solver,
                doublereal rtol, int maxsteps, int maxiter, int loglevel)
{
    MultiPhase* m = 0;
    ChemEquil* e = 0;
    bool redo = true;
    int retn = -1;
    int nAttempts = 0;
    int retnSub = 0;


    if (loglevel > 0) {
        beginLogGroup("equilibrate", loglevel);
        addLogEntry("Single-phase equilibrate function");
        {
            beginLogGroup("arguments");
            addLogEntry("phase",s.id());
            addLogEntry("XY",XY);
            addLogEntry("solver",solver);
            addLogEntry("rtol",rtol);
            addLogEntry("maxsteps",maxsteps);
            addLogEntry("maxiter",maxiter);
            addLogEntry("loglevel",loglevel);
            endLogGroup("arguments");
        }
    }
    while (redo) {
        if (solver >= 2) {
#ifdef WITH_VCSNONIDEAL
            int printLvlSub = 0;
            int estimateEquil = 0;
            try {
                MultiPhase m;
                m.addPhase(&s, 1.0);
                m.init();
                nAttempts++;
                vcs_equilibrate(m, XY, estimateEquil, printLvlSub, solver,
                                rtol, maxsteps, maxiter, loglevel-1);
                redo = false;
                if (loglevel > 0) {
                    addLogEntry("VCSnonideal solver succeeded.");
                }
                retn = nAttempts;
            } catch (CanteraError& err) {
                err.save();
                if (loglevel > 0) {
                    addLogEntry("VCSnonideal solver failed.");
                }
                if (nAttempts < 2) {
                    if (loglevel > 0) {
                        addLogEntry("Trying single phase ChemEquil solver.");
                    }
                    solver = -1;
                } else {
                    if (loglevel > 0) {
                        endLogGroup("equilibrate");
                    }
                    throw err;
                }
            }
#else
            throw CanteraError("equilibrate",
                               "VCSNonIdeal solver called, but not compiled");
#endif
        } else if (solver == 1) {
            try {
                MultiPhase m;
                m.addPhase(&s, 1.0);
                m.init();
                nAttempts++;
                equilibrate(m, XY, rtol, maxsteps, maxiter, loglevel-1);
                redo = false;
                if (loglevel > 0) {
                    addLogEntry("MultiPhaseEquil solver succeeded.");
                }
                retn = nAttempts;
            } catch (CanteraError& err) {
                err.save();
                if (loglevel > 0) {
                    addLogEntry("MultiPhaseEquil solver failed.");
                }
                if (nAttempts < 2) {
                    if (loglevel > 0) {
                        addLogEntry("Trying single phase ChemEquil solver.");
                    }
                    solver = -1;
                } else {
                    if (loglevel > 0) {
                        endLogGroup("equilibrate");
                    }
                    throw err;
                }
            }
        } else {      // solver <= 0
            /*
             * Call the element potential solver
             */
            try {
                ChemEquil e;
                e.options.maxIterations = maxsteps;
                e.options.relTolerance = rtol;
                nAttempts++;
                bool useThermoPhaseElementPotentials = true;
                retnSub = e.equilibrate(s, XY, useThermoPhaseElementPotentials,
                                        loglevel-1);
                if (retnSub < 0) {
                    if (loglevel > 0) {
                        addLogEntry("ChemEquil solver failed.");
                    }
                    if (nAttempts < 2) {
                        if (loglevel > 0) {
                            addLogEntry("Trying MultiPhaseEquil solver.");
                        }
                        solver = 1;
                    } else {
                        throw CanteraError("equilibrate",
                                           "Both equilibrium solvers failed");
                    }
                }
                retn = nAttempts;
                s.setElementPotentials(e.elementPotentials());
                redo = false;
                if (loglevel > 0) {
                    addLogEntry("ChemEquil solver succeeded.");
                }
            }

            catch (CanteraError& err) {
                err.save();
                if (loglevel > 0) {
                    addLogEntry("ChemEquil solver failed.");
                }
                // If ChemEquil fails, try the MultiPhase solver
                if (solver < 0) {
                    if (loglevel > 0) {
                        addLogEntry("Trying MultiPhaseEquil solver.");
                    }
                    solver = 1;
                } else {
                    redo = false;
                    if (loglevel > 0) {
                        endLogGroup("equilibrate");
                    }
                    throw err;
                }
            }
        }
    } // while (redo)
    /*
     * We are here only for a success
     */
    if (loglevel > 0) {
        endLogGroup("equilibrate");
    }
    return retn;
}
}
