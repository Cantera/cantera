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
#include "cantera/equil/vcs_prob.h"
#include "cantera/equil/vcs_internal.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_species_thermo.h"
#include "cantera/equil/vcs_SpeciesProperties.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/equil.h"

#include "cantera/base/ct_defs.h"
#include "cantera/base/stringUtils.h"
#include "cantera/thermo/mix_defs.h"
#include "cantera/thermo/speciesThermoTypes.h"
#include "cantera/thermo/IdealSolidSolnPhase.h"
#include "cantera/thermo/IdealMolalSoln.h"
#include "cantera/equil/ChemEquil.h"

#include <string>
#include <vector>

using namespace std;

namespace Cantera
{

/*
 *  Set a single-phase chemical solution to chemical equilibrium.
 *  This is a convenience function that uses one or the other of
 *  the two chemical equilibrium solvers.
 *
 *  @param s The object to set to an equilibrium state
 *
 *  @param XY An integer specifying the two properties to be held
 *            constant.
 *
 *  @param estimateEquil integer indicating whether the solver
 *                   should estimate its own initial condition.
 *                   If 0, the initial mole fraction vector
 *                   in the %ThermoPhase object is used as the
 *                   initial condition.
 *                   If 1, the initial mole fraction vector
 *                   is used if the element abundances are
 *                   satisfied.
 *                   if -1, the initial mole fraction vector
 *                   is thrown out, and an estimate is
 *                   formulated.
 *
 *  @param printLvl Determines the amount of printing that
 *                  gets sent to stdout from the vcs package
 *                  (Note, you may have to compile with debug
 *                   flags to get some printing).
 *
 *  @param solver The equilibrium solver to use. If solver = 0,
 *                the ChemEquil solver will be used, and if
 *                solver = 1, the vcs_MultiPhaseEquil solver will
 *                be used (slower than ChemEquil,
 *                but more stable). If solver < 0 (default, then
 *                ChemEquil will be tried first, and if it fails
 *                vcs_MultiPhaseEquil will be tried.
 *
 *  @param maxsteps The maximum number of steps to take to find
 *                  the solution.
 *
 *  @param maxiter For the MultiPhaseEquil solver only, this is
 *                 the maximum number of outer temperature or
 *                 pressure iterations to take when T and/or P is
 *                 not held fixed.
 *
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 */
int vcs_equilibrate(thermo_t& s, const char* XY,
                    int estimateEquil,  int printLvl,
                    int solver,
                    doublereal rtol, int maxsteps, int maxiter,
                    int loglevel)
{
    MultiPhase* m = 0;
    int retn = 1;
    int retnSub = 0;

    beginLogGroup("equilibrate", loglevel);
    // retry:
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
            if (retn == 1) {
                addLogEntry("MultiPhaseEquil solver succeeded.");
            } else {
                addLogEntry("MultiPhaseEquil solver returned an error code: ", retn);
            }
            delete m;
        } catch (CanteraError& err) {
            err.save();
            addLogEntry("MultiPhaseEquil solver failed.");
            delete m;
            throw err;
        }
    } else if (solver == 1) {
        m = new MultiPhase;
        try {
            m->addPhase(&s, 1.0);
            m->init();
            (void) equilibrate(*m, XY, rtol, maxsteps, maxiter, loglevel-1);
            if (loglevel > 0) {
                addLogEntry("MultiPhaseEquil solver succeeded.");
            }
            delete m;
            retn = 1;
        } catch (CanteraError& err) {
            err.save();
            if (loglevel > 0) {
                addLogEntry("MultiPhaseEquil solver failed.");
            }
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
            retnSub = e->equilibrate(s, XY,
                                     useThermoPhaseElementPotentials, loglevel-1);
            if (retnSub < 0) {
                if (loglevel > 0) {
                    addLogEntry("ChemEquil solver failed.");
                }
                delete e;
                throw CanteraError("equilibrate",
                                   "ChemEquil equilibrium solver failed");
            }
            retn = 1;
            s.setElementPotentials(e->elementPotentials());
            delete e;
            if (loglevel > 0) {
                addLogEntry("ChemEquil solver succeeded.");
            }
        } catch (CanteraError& err) {
            err.save();
            if (loglevel > 0) {
                addLogEntry("ChemEquil solver failed.");
            }
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
    endLogGroup("equilibrate");
    return retn;
}

//  Set a multi-phase chemical solution to chemical equilibrium.
/*
 *  This function uses the vcs_MultiPhaseEquil interface to the
 *  vcs solver.
 *  The function uses the element abundance vector that is
 *  currently consistent with the composition within the phases
 *  themselves. Two other thermodynamic quantities, determined by the
 *  XY string,  are held constant during the equilibration.
 *
 *  @param s The object to set to an equilibrium state
 *
 *  @param XY A character string specifying the two properties to
 *            be held constant
 *
 *  @param estimateEquil integer indicating whether the solver
 *                   should estimate its own initial condition.
 *                   If 0, the initial mole fraction vector
 *                   in the %ThermoPhase object is used as the
 *                   initial condition.
 *                   If 1, the initial mole fraction vector
 *                   is used if the element abundances are
 *                   satisfied.
 *                   if -1, the initial mole fraction vector
 *                   is thrown out, and an estimate is
 *                   formulated.
 *
 *  @param printLvl Determines the amount of printing that
 *                  gets sent to stdout from the vcs package
 *                  (Note, you may have to compile with debug
 *                   flags to get some printing).
 *
 *  @param maxsteps The maximum number of steps to take to find
 *                  the solution.
 *
 *  @param maxiter For the MultiPhaseEquil solver only, this is
 *                 the maximum number of outer temperature or
 *                 pressure iterations to take when T and/or P is
 *                 not held fixed.
 *
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 *
 *  @ingroup equilfunctions
 */
int vcs_equilibrate(MultiPhase& s, const char* XY,
                    int estimateEquil, int printLvl, int solver,
                    doublereal tol, int maxsteps, int maxiter,
                    int loglevel)
{
    int ixy = _equilflag(XY);
    int retn = vcs_equilibrate_1(s, ixy, estimateEquil, printLvl, solver,
                                 tol, maxsteps, maxiter, loglevel);
    return retn;
};


//  Set a multi-phase chemical solution to chemical equilibrium.
/*
 *  This function uses the vcs_MultiPhaseEquil interface to the
 *  vcs solver.
 *  The function uses the element abundance vector that is
 *  currently consistent with the composition within the phases
 *  themselves. Two other thermodynamic quantities, determined by the
 *  XY string,  are held constant during the equilibration.
 *
 *  @param s The object to set to an equilibrium state
 *
 *  @param XY An integer specifying the two properties to be held
 *            constant.
 *
 *  @param estimateEquil integer indicating whether the solver
 *                   should estimate its own initial condition.
 *                   If 0, the initial mole fraction vector
 *                   in the %ThermoPhase object is used as the
 *                   initial condition.
 *                   If 1, the initial mole fraction vector
 *                   is used if the element abundances are
 *                   satisfied.
 *                   if -1, the initial mole fraction vector
 *                   is thrown out, and an estimate is
 *                   formulated.
 *
 *  @param printLvl Determines the amount of printing that
 *                  gets sent to stdout from the vcs package
 *                  (Note, you may have to compile with debug
 *                   flags to get some printing).
 *
 *  @param maxsteps The maximum number of steps to take to find
 *                  the solution.
 *
 *  @param maxiter For the MultiPhaseEquil solver only, this is
 *                 the maximum number of outer temperature or
 *                 pressure iterations to take when T and/or P is
 *                 not held fixed.
 *
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 *
 *  @ingroup equilfunctions
 */
int vcs_equilibrate_1(MultiPhase& s, int ixy,
                      int estimateEquil, int printLvl, int solver,
                      doublereal tol, int maxsteps, int maxiter, int loglevel)
{
    static int counter = 0;
    int retn = 1;

    beginLogGroup("equilibrate",loglevel);
    addLogEntry("multiphase equilibrate function");
    beginLogGroup("arguments");
    addLogEntry("XY",ixy);
    addLogEntry("tol",tol);
    addLogEntry("maxsteps",maxsteps);
    addLogEntry("maxiter",maxiter);
    addLogEntry("loglevel",loglevel);
    endLogGroup("arguments");

    int printLvlSub = std::max(0, printLvl-1);

    s.init();

    if (solver == 2) {
        try {
            VCSnonideal::vcs_MultiPhaseEquil* eqsolve =  new VCSnonideal::vcs_MultiPhaseEquil(&s, printLvlSub);
            int err = eqsolve->equilibrate(ixy, estimateEquil, printLvlSub, tol, maxsteps, loglevel);
            if (err != 0) {
                retn = -1;
                addLogEntry("vcs_equilibrate Error   - ", err);
            } else {
                addLogEntry("vcs_equilibrate Success - ", err);
            }
            endLogGroup("equilibrate");
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
            addLogEntry("Failure.", lastErrorMessage());
            endLogGroup("equilibrate");
            throw e;
        }
    } else if (solver == 1) {
        if (ixy == TP || ixy == HP || ixy == SP || ixy == TV) {
            try {
                double err = s.equilibrate(ixy, tol, maxsteps, maxiter, loglevel);

                addLogEntry("Success. Error",err);
                endLogGroup("equilibrate");

                return 0;
            } catch (CanteraError& e) {
                e.save();
                addLogEntry("Failure.",lastErrorMessage());
                endLogGroup("equilibrate");

                throw e;
            }
        } else {

            addLogEntry("multiphase equilibrium can be done only for TP, HP, SP, or TV");
            endLogGroup("equilibrate");

            throw CanteraError("equilibrate","unsupported option");
            //return -1.0;
        }
    } else {
        throw CanteraError("vcs_equilibrate_1", "unknown solver");
    }
    return retn;
}

//====================================================================================================================
// Determine the phase stability of a single phase given the current conditions
// in a MultiPhase object
/*
 *
 *  @param s         The MultiPhase object to be set to an equilibrium state
 *  @param iphase    Phase index within the multiphase object to be
 *                   tested for stability.
 *  @param funcStab  Function value that tests equilibrium. > 0 indicates stable
 *                   < 0 indicates unstable
 *
 *  @param printLvl   Determines the amount of printing that
 *                  gets sent to stdout from the vcs package
 *                  (Note, you may have to compile with debug
 *                   flags to get some printing).
 *
 *  @param loglevel Controls amount of diagnostic output. loglevel
 *                  = 0 suppresses diagnostics, and increasingly-verbose
 *                  messages are written as loglevel increases. The
 *                  messages are written to a file in HTML format for viewing
 *                  in a web browser. @see HTML_logs
 */
int vcs_determine_PhaseStability(MultiPhase& s, int iphase,
                                 double& funcStab, int printLvl, int loglevel)
{
    int iStab = 0;
    static int counter = 0;
    beginLogGroup("PhaseStability",loglevel);
    addLogEntry("multiphase phase stability function");
    beginLogGroup("arguments");
    addLogEntry("iphase",iphase);
    addLogEntry("loglevel",loglevel);
    endLogGroup("arguments");

    int printLvlSub = std::max(0, printLvl-1);

    s.init();
    try {
        VCSnonideal::vcs_MultiPhaseEquil* eqsolve = new VCSnonideal::vcs_MultiPhaseEquil(&s, printLvlSub);
        iStab = eqsolve->determine_PhaseStability(iphase, funcStab, printLvlSub, loglevel);
        if (iStab != 0) {
            addLogEntry("Phase is stable  - ", iphase);
        } else {
            addLogEntry("Phase is not stable - ", iphase);
        }
        endLogGroup("PhaseStability");
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
        addLogEntry("Failure.", lastErrorMessage());
        endLogGroup("equilibrate");
        throw e;
    }
    return iStab;
}
//====================================================================================================================
}
