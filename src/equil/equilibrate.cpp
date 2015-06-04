/**
 * @file equilibrate.cpp Driver routines for the chemical equilibrium solvers.
 */

#include "cantera/equil/ChemEquil.h"
#include "cantera/equil/vcs_MultiPhaseEquil.h"

namespace Cantera
{

doublereal equilibrate(MultiPhase& s, const char* XY,
                       doublereal tol, int maxsteps, int maxiter,
                       int loglevel)
{
    warn_deprecated("equilibrate(MultiPhase&, ...)",
        "Use MultiPhase::equilibrate instead. To be removed after Cantera 2.2.");
    s.init();
    int ixy = _equilflag(XY);
    if (ixy == TP || ixy == HP || ixy == SP || ixy == TV) {
        try {
            double err = s.equilibrate(ixy, tol, maxsteps, maxiter, loglevel);
            return err;
        } catch (CanteraError& err) {
            err.save();
            throw err;
        }
    } else {
        throw CanteraError("equilibrate","unsupported option");
        return -1.0;
    }
    return 0.0;
}

int equilibrate(thermo_t& s, const char* XY, int solver,
                doublereal rtol, int maxsteps, int maxiter, int loglevel)
{
    warn_deprecated("equilibrate(ThermoPhase&, ...)",
        "Use ThermoPhase::equilibrate instead. To be removed after Cantera 2.2.");
    bool redo = true;
    int retn = -1;
    int nAttempts = 0;
    int retnSub = 0;

    while (redo) {
        if (solver >= 2) {
            int printLvlSub = loglevel;
            int estimateEquil = 0;
            try {
                MultiPhase m;
                m.addPhase(&s, 1.0);
                m.init();
                nAttempts++;
                vcs_equilibrate(m, XY, estimateEquil, printLvlSub, solver,
                                rtol, maxsteps, maxiter, loglevel-1);
                redo = false;
                retn = nAttempts;
            } catch (CanteraError& err) {
                err.save();
                if (nAttempts < 2) {
                    solver = -1;
                } else {
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
                retn = nAttempts;
            } catch (CanteraError& err) {
                err.save();
                if (nAttempts < 2) {
                    solver = -1;
                } else {
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
                    if (nAttempts < 2) {
                        solver = 1;
                    } else {
                        throw CanteraError("equilibrate",
                                           "Both equilibrium solvers failed");
                    }
                }
                retn = nAttempts;
                s.setElementPotentials(e.elementPotentials());
                redo = false;
            }

            catch (CanteraError& err) {
                err.save();
                // If ChemEquil fails, try the MultiPhase solver
                if (solver < 0) {
                    solver = 1;
                } else {
                    redo = false;
                    throw err;
                }
            }
        }
    } // while (redo)
    /*
     * We are here only for a success
     */
    return retn;
}
}
