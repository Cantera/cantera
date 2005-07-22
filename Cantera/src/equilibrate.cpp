/**
 * @file equilibrate.cpp
 *
 * Driver routines for the chemical equilibrium solvers.
 *
 */

#include "ChemEquil.h"
#include "MultiPhaseEquil.h"

namespace Cantera {


 /**
     * Set a mixture to a state of chemical equilibrium. The flag 'XY'
     * determines the two properties that will be held fixed in the
     * calculation.
     */
    doublereal equilibrate(MultiPhase& s, const char* XY, 
        doublereal tol = 1.0e-9, int maxsteps = 1000, int maxiter = 100, 
        int loglevel = -99) {

        beginLogGroup("equilibrate",loglevel);
        addLogEntry("multiphase equilibrate function");
        beginLogGroup("arguments");
        addLogEntry("XY",XY);
        addLogEntry("tol",tol);
        addLogEntry("maxsteps",maxsteps);
        addLogEntry("maxiter",maxiter);
        addLogEntry("loglevel",loglevel);
        endLogGroup("arguments");
        
        s.init();
        int ixy = _equilflag(XY);
        if (ixy == TP || ixy == HP || ixy == SP || ixy == TV) {
            try {
                double err = s.equilibrate(ixy, tol, maxsteps, maxiter);
                addLogEntry("Success. Error",err);
                endLogGroup("equilibrate");
                return err;
            }
            catch (CanteraError e) {
                addLogEntry("Failure.",lastErrorMessage());
                endLogGroup("equilibrate");
                throw e;
            }
        }
        else {
            addLogEntry("multiphase equilibrium can be done only for TP, HP, SP, or TV");
            endLogGroup("equilibrate");
            throw CanteraError("equilibrate","unsupported option");
            return -1.0;
        }
    }

    /// Set a single-phase chemical solution to chemical equilibrium.
    /// This is a convenience function that uses one or the other of
    /// the two chemical equilibrium solvers.  @param The object to
    /// set to an equilibrium state @param XY An integer specifying
    /// the two properties to be held constant.  @param solver The
    /// equilibrium solver to use. If solver = 0, the ChemEquil solver
    /// will be used, and if solver = 1, the MultiPhaseEquil solver
    /// will be used (slower than ChemEquil, but more stable). If
    /// solver < 0 (default, then ChemEquil will be tried first, and
    /// if it fails MultiPhaseEquil will be tried.  @param maxsteps
    /// The maximum number of steps to take to find the solution.
    /// @param maxiter For the MultiPhaseEquil solver only, this is
    /// the maximum number of outer temperature or pressure iterations
    /// to take when T and/or P is not held fixed.  @param loglevel
    /// Controls amount of diagnostic output. loglevel = 0 suppresses
    /// diagnostics, and increasingly-verbose messages are written as
    /// loglevel increases. The messages are written to a file in HTML
    /// format for viewing in a web browser.

    void equilibrate(thermo_t& s, const char* XY, int solver,
        doublereal rtol, int maxsteps, int maxiter, int loglevel) {
        MultiPhase* m = 0;
        ChemEquil* e = 0;
        bool redo = true;

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
        while (redo) {
            if (solver > 0) {
                m = new MultiPhase;
                try { 
                    m->addPhase(&s, 1.0);
                    m->init();
                    equilibrate(*m, XY, rtol, maxsteps, maxiter, loglevel); 
                    redo = false;
                    addLogEntry("MultiPhaseEquil solver succeeded.");
                    delete m;
                }
                catch (CanteraError err) {
                    addLogEntry("MultiPhaseEquil solver failed.");
                    endLogGroup("equilibrate");
                    delete m;
                    throw err;
                }
            }
            else {        // solver <= 0
                e = new ChemEquil;
                try {
                    e->options.maxIterations = maxsteps;
                    e->options.relTolerance = rtol;
                    e->equilibrate(s,XY);
                    s.setElementPotentials(e->elementPotentials());
                    redo = false;
                    delete e;
                    addLogEntry("ChemEquil solver succeeded.");
                }

                catch (CanteraError err) {
                    delete e;
                    addLogEntry("ChemEquil solver failed.");
                    // If ChemEquil fails, try the MultiPhase solver
                    if (solver < 0) {
                        addLogEntry("Trying MultiPhaseEquil solver.");
                        solver = 1;
                    }
                    else {
                        redo = false;
                        endLogGroup("equilibrate");
                        throw err;
                    }
                }
            } 
        } // while (redo)
        endLogGroup("equilibrate");
    }
}
