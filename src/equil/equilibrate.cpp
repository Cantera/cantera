/**
 * @file equilibrate.cpp Driver routines for the chemical equilibrium solvers.
 */

#include "cantera/equil/equil.h"
#include "cantera/equil/ChemEquil.h"
#include "cantera/equil/MultiPhaseEquil.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/base/global.h"

namespace Cantera
{

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

int equilibrate(thermo_t& s, const char* XY, int solver,
                doublereal rtol, int maxsteps, int maxiter, int loglevel)
{
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
