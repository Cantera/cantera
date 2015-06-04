/**
 *  @file vcs_equilibrate.cpp
 *    Driver routines for equilibrium solvers
 */
/*
 * Copyright (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include "cantera/equil/vcs_MultiPhaseEquil.h"
#include "cantera/equil/equil.h"

#include "cantera/base/stringUtils.h"
#include "cantera/equil/ChemEquil.h"

using namespace std;

namespace Cantera
{
int vcs_equilibrate(thermo_t& s, const char* XY,
                    int estimateEquil,  int printLvl,
                    int solver,
                    doublereal rtol, int maxsteps, int maxiter,
                    int loglevel)
{
    warn_deprecated("vcs_equilibrate", "Use ThermoPhase::equilibrate instead. "
        "To be removed after Cantera 2.2.");
    MultiPhase* m = 0;
    int retn = 1;

    if (solver == 2) {
        m = new MultiPhase;
        try {
            /*
             *  Set the kmoles of the phase to 1.0, arbitrarily.
             *  It actually doesn't matter.
             */
            m->addPhase(&s, 1.0);
            m->init();

            retn = vcs_equilibrate(*m, XY, estimateEquil, printLvl, solver,
                                   rtol, maxsteps, maxiter, loglevel);
            delete m;
        } catch (CanteraError& err) {
            err.save();
            delete m;
            throw err;
        }
    } else if (solver == 1) {
        m = new MultiPhase;
        try {
            m->addPhase(&s, 1.0);
            m->init();
            (void) equilibrate(*m, XY, rtol, maxsteps, maxiter, loglevel-1);
            delete m;
            retn = 1;
        } catch (CanteraError& err) {
            err.save();
            delete m;
            throw err;
        }
    } else if (solver == 0) {
        ChemEquil* e = new ChemEquil;
        try {
            e->options.maxIterations = maxsteps;
            e->options.relTolerance = rtol;
            bool useThermoPhaseElementPotentials = false;
            if (estimateEquil == 0) {
                useThermoPhaseElementPotentials = true;
            }
            int retnSub = e->equilibrate(s, XY,
                                     useThermoPhaseElementPotentials, loglevel-1);
            if (retnSub < 0) {
                delete e;
                throw CanteraError("equilibrate",
                                   "ChemEquil equilibrium solver failed");
            }
            retn = 1;
            s.setElementPotentials(e->elementPotentials());
            delete e;
        } catch (CanteraError& err) {
            err.save();
            delete e;
            throw err;
        }
    } else {
        throw CanteraError("vcs_equilibrate",
                           "unknown solver");
    }

    /*
     * We are here only for a success
     */
    return retn;
}

int vcs_equilibrate(MultiPhase& s, const char* XY,
                    int estimateEquil, int printLvl, int solver,
                    doublereal tol, int maxsteps, int maxiter,
                    int loglevel)
{
    int ixy = _equilflag(XY);
    int retn = vcs_equilibrate_1(s, ixy, estimateEquil, printLvl, solver,
                                 tol, maxsteps, maxiter, loglevel);
    return retn;
}

int vcs_equilibrate_1(MultiPhase& s, int ixy,
                      int estimateEquil, int printLvl, int solver,
                      doublereal tol, int maxsteps, int maxiter, int loglevel)
{
    warn_deprecated("vcs_equilibrate_1", "Use MultiPhase::equilibrate instead. "
        "To be removed after Cantera 2.2.");
    static int counter = 0;
    int retn = 1;

    int printLvlSub = std::max(0, printLvl-1);

    s.init();

    if (solver == 2) {
        try {
            vcs_MultiPhaseEquil* eqsolve =  new vcs_MultiPhaseEquil(&s, printLvlSub);
            int err = eqsolve->equilibrate(ixy, estimateEquil, printLvlSub, tol, maxsteps, loglevel);
            if (err != 0) {
                retn = -1;
            }
            // hard code a csv output file.
            if (printLvl > 0) {
                string reportFile = "vcs_equilibrate_res.csv";
                if (counter > 0) {
                    reportFile = "vcs_equilibrate_res_" + int2str(counter) + ".csv";
                }
                eqsolve->reportCSV(reportFile);
                counter++;
            }
            delete eqsolve;
        } catch (CanteraError& e) {
            e.save();
            retn = -1;
            throw e;
        }
    } else if (solver == 1) {
        if (ixy == TP || ixy == HP || ixy == SP || ixy == TV) {
            try {
                s.equilibrate(ixy, tol, maxsteps, maxiter, loglevel);
                return 0;
            } catch (CanteraError& e) {
                e.save();
                throw e;
            }
        } else {
            throw CanteraError("equilibrate","unsupported option");
        }
    } else {
        throw CanteraError("vcs_equilibrate_1", "unknown solver");
    }
    return retn;
}

int vcs_determine_PhaseStability(MultiPhase& s, int iphase,
                                 double& funcStab, int printLvl, int loglevel)
{
    int iStab = 0;
    static int counter = 0;
    int printLvlSub = std::max(0, printLvl-1);

    s.init();
    try {
        vcs_MultiPhaseEquil* eqsolve = new vcs_MultiPhaseEquil(&s, printLvlSub);
        iStab = eqsolve->determine_PhaseStability(iphase, funcStab, printLvlSub, loglevel);
        // hard code a csv output file.
        if (printLvl > 0) {
            string reportFile = "vcs_phaseStability.csv";
            if (counter > 0) {
                reportFile = "vcs_phaseStability_" + int2str(counter) + ".csv";
            }
            eqsolve->reportCSV(reportFile);
            counter++;
        }
        delete eqsolve;
    } catch (CanteraError& e) {
        throw e;
    }
    return iStab;
}

}
