/**
 *
 *  @file DAE_Solver.h
 *
 *  Header file for class DAE_Solver
 */

/*  
 *  $Date: 2009/03/27 23:39:26 $
 *  $Revision: 1.4 $
 *
 *  Copyright 2006 California Institute of Technology
 *
 */

#undef DAE_DEVEL 

#ifndef CT_DAE_Solver_H
#define CT_DAE_Solver_H

#include <vector>

#include "ct_defs.h"
#include "ResidEval.h"
#include "global.h"

namespace Cantera {

#ifdef DAE_DEVEL
    /**
     * @defgroup numerics  Numerical Utilities within Cantera
     *
     *
     */

    class Jacobian {
    public:
        Jacobian(){}
        virtual ~Jacobian(){}
        virtual bool supplied() { return false; }
        virtual bool isBanded() { return false; }
        virtual int lowerBandWidth() { return 0; }
        virtual int upperBandWidth() { return 0; }
    };

    class BandedJacobian : public Jacobian {
    public:
        BandedJacobian(int ml, int mu) {
            m_ml = ml; m_mu = mu;
        }
        virtual bool supplied() { return false; }
        virtual bool isBanded() { return true; }
        virtual int lowerBandWidth() { return m_ml; }
        virtual int upperBandWidth() { return m_mu; }
    protected:
        int m_ml, m_mu;
    };

    const int cDirect = 0;
    const int cKrylov = 1;



    /**
     * Wrapper for DAE solvers
     */
    class DAE_Solver {
    public:

        DAE_Solver(ResidEval& f) : m_resid(f),
                                   m_neq(f.nEquations()), 
                                   m_time(0.0) {}

        virtual ~DAE_Solver(){}

        /** 
         * Set error tolerances. This version specifies a scalar
         * relative tolerance, and a vector absolute tolerance.
         */
        virtual void setTolerances(doublereal reltol, 
            doublereal* abstol) {
            warn("setTolerances");
        }

        /** 
         * Set error tolerances. This version specifies a scalar
         * relative tolerance, and a scalar absolute tolerance.
         */
        virtual void setTolerances(doublereal reltol, doublereal abstol) {
            warn("setTolerances");
        }

        /** 
         * Specify a Jacobian evaluator. If this method is not called,
         * the Jacobian will be computed by finite difference.
         */
        void setJacobian(Jacobian& jac) {
            warn("setJacobian");
        }

        virtual void setLinearSolverType(int solverType) {
            warn("setLinearSolverType");
        }

        virtual void setDenseLinearSolver() {
            warn("setDenseLinearSolver");
        }

        virtual void setBandedLinearSolver(int m_upper, int m_lower) {
            warn("setBandedLinearSolver");
        }
        virtual void setMaxTime(doublereal tmax) {
            warn("setMaxTime");
        }
        virtual void setMaxStepSize(doublereal dtmax) {
            warn("setMaxStepSize");
        }
        virtual void setMaxOrder(int n) {
            warn("setMaxOrder");
        }
        virtual void setMaxNumSteps(int n) {
            warn("setMaxNumSteps");
        }
        virtual void setInitialStepSize(doublereal h0) {
            warn("setInitialStepSize");
        }
        virtual void setStopTime(doublereal tstop) {
            warn("setStopTime");
        }
        virtual void setMaxErrTestFailures(int n) {
            warn("setMaxErrTestFailures");
        }
        virtual void setMaxNonlinIterations(int n) {
            warn("setMaxNonlinIterations");
        }
        virtual void setMaxNonlinConvFailures(int n) {
            warn("setMaxNonlinConvFailures");
        }
        virtual void inclAlgebraicInErrorTest(bool yesno) {
            warn("inclAlgebraicInErrorTest");
        }

        virtual void correctInitial_Y_given_Yp() {
            warn("correctInitial_Y_given_Yp");
        }

        virtual void correctInitial_YaYp_given_Yd() {
            warn("correctInitial_YaYp_given_Yd");
        }

        /**
         * Solve the system of equations up to time tout.
         */
        virtual int solve(doublereal tout) {
            warn("solve"); return 0;
        }

        /**
         * Take one internal step.
         */
        virtual doublereal step(doublereal tout) {
            warn("step"); return 0;
        }

        /// Number of equations.
        int nEquations() const { return m_resid.nEquations(); }

        /**
         * initialize. Base class method does nothing.
         */
        virtual void init(doublereal t0) {}

        /**
         * Set a solver-specific input parameter.
         */
        virtual void setInputParameter(int flag, doublereal value) {
            warn("setInputParameter");
        }

        /**
         * Get the value of a solver-specific output parameter.
         */
        virtual doublereal getOutputParameter(int flag) const {
            warn("getOutputParameter"); return 0.0;
        }

        /// the current value of solution component k.
        virtual doublereal solution(int k) const {
            warn("solution"); return 0.0;
        }

        virtual const doublereal* solutionVector() const {
            warn("solutionVector"); return &m_dummy;
        }

        /// the current value of the derivative of solution component k.
        virtual doublereal derivative(int k) const {
            warn("derivative"); return 0.0;
        }

        virtual const doublereal* derivativeVector() const {
            warn("derivativeVector"); return &m_dummy;
        }

    protected:

        doublereal m_dummy;

        ResidEval& m_resid;

        integer    m_neq;
        doublereal m_time;

        
    private:
        void warn(std::string msg) const {
            writelog(">>>> Warning: method "+msg+" of base class "
                +"DAE_Solver called. Nothing done.\n");
        }
    };

#endif

}

#endif
