/**
 *
 *  @file IDA_Solver.h
 *
 *  Header file for class IDA_Solver
 */

/*  $Author: hkmoffa $
 *  $Date: 2008/04/08 20:19:27 $
 *  $Revision: 1.3 $
 *
 *  Copyright 2006 California Institute of Technology
 *
 */

#ifndef CT_IDA_Solver_H
#define CT_IDA_Solver_H

#include <vector>

#include "DAE_Solver.h"
#include "ctexceptions.h"

namespace Cantera {

    /**
     * Exception thrown when a IDA error is encountered.
     */
    class IDA_Err : public CanteraError {
    public:
        IDA_Err(std::string msg) : CanteraError("IDA_Solver", msg){}
    };


    class ResidData; // forward reference

    class IDA_Solver : public DAE_Solver {
    public:

        IDA_Solver(ResidEval& f);

        virtual ~IDA_Solver();

        /** 
         * Set error tolerances. This version specifies a scalar
         * relative tolerance, and a vector absolute tolerance.
         */
        virtual void setTolerances(doublereal reltol, 
            doublereal* abstol);

        /** 
         * Set error tolerances. This version specifies a scalar
         * relative tolerance, and a scalar absolute tolerance.
         */
        virtual void setTolerances(doublereal reltol, doublereal abstol);

        virtual void setLinearSolverType(int solverType);

        virtual void setDenseLinearSolver();
        virtual void setBandedLinearSolver(int m_upper, int m_lower);

        virtual void setMaxTime(doublereal tmax);

        virtual void setMaxOrder(int n);

        virtual void setMaxNumSteps(int n);
        virtual void setInitialStepSize(doublereal h0);
        virtual void setStopTime(doublereal tstop);
        virtual void setMaxErrTestFailures(int n);
        virtual void setMaxNonlinIterations(int n);
        virtual void setMaxNonlinConvFailures(int n);
        virtual void inclAlgebraicInErrorTest(bool yesno);

        virtual void setInputParameter(int flag, doublereal value);
        virtual doublereal getOutputParameter(int flag);


        /**
         * This method may be called if the initial conditions do not
         * satisfy the residual equation F = 0. Given the derivatives
         * of all variables, this method computes the initial y
         * values.
         */
        virtual void correctInitial_Y_given_Yp(doublereal* y, doublereal* yp, 
            doublereal tout);

        /**
         * This method may be called if the initial conditions do not
         * satisfy the residual equation F = 0. Given the initial
         * values of all differential variables, it computes the
         * initial values of all algebraic variables and the initial
         * derivatives of all differential variables.
         */
        virtual void correctInitial_YaYp_given_Yd(doublereal* y, doublereal* yp,
            doublereal tout);


        virtual int solve(doublereal tout);

        virtual doublereal step(doublereal tout);

        virtual void init(doublereal t0);

        /// the current value of solution component k.
        virtual doublereal solution(int k) const;

        virtual const doublereal* solutionVector() const;

        /// the current value of the derivative of solution component k.
        virtual doublereal derivative(int k) const;

        virtual const doublereal* derivativeVector() const;

    protected:

	int m_neq;
        void* m_ida_mem;
        doublereal m_t0;
        void *m_y, *m_ydot, *m_id, *m_constraints, *m_abstol;
        int m_type;
        int m_itol;
        int m_iter;
        doublereal m_reltol;
        doublereal m_abstols;
        int m_nabs;
        doublereal m_hmax, m_hmin;
        int m_maxsteps, m_maxord;
        ResidData* m_fdata;
        int m_mupper, m_mlower;        
    };

}

#endif

