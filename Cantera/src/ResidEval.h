/**
 *  @file FuncEval.h
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology

#ifndef CT_RESIDEVAL_H
#define CT_RESIDEVAL_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

namespace Cantera {


    /**
     *
     *  Virtual base class for DAE residual function evaluators.
     *
     */
    class ResidEval {

    public:

        /**
         * Evaluate the residual function.  
         * @param t time (input, do not modify) 
         * @param y solution vector (input, do not modify)
         * @param ydot rate of change of solution vector. (input, do
         * not modify)
         */
	virtual void evalResid(double t, const double deltaT, 
			       const double* y,
			       const double* ydot,
			       double* resid)=0;

        /** Number of equations. */
        virtual int neq()=0;

        //        virtual bool init(double t0)=0;
        virtual double* solution()=0;
        virtual double* solution_dot()=0;

        /**
         *       Fill the solution vector with the initial conditions
         *       at initial time t0.
         */
        virtual void getInitialConditions(double t0, size_t leny, double* y)=0;

        /**
         *      Report the analytical result if available and describe
         *      how close the numerical solution is to the analytical
         *      solution
         */
        virtual void analytical_solution(double, double *)=0;
        virtual void analytical_solution_dot(double, double *)=0;

        /**
         *      Write out to a file or to standard output the current solution
         *      ievent  is a description of the event that caused this
         *      function to be called.
         */
        virtual void writeSolution(int ievent, double t,
                                   const double *y, const double *ydot) = 0;

        /**
         *     Apply a filter
         */
        virtual void applyFilter(const double time, const double deltaT,
                                 const double *ydot, const double *y,
                                 double *delta_y) = 0;

        /**
         *      Write out to a file or to standard output the current solution
         *      ievent  is a description of the event that caused this
         *      function to be called.
         */
        virtual void writeSolutionFilter(int ievent, double t,
                                         const double * const y,
                                         const double * const ydot,
                                         const double * const delta_y,
                                         const char * const fname) = 0;

	/**
         * This function returns a value of the delta t constraint
         * that may exist in the problem.
         *
         */
        virtual double delta_t_constraint(const double t, const double *y,
                                          const double *ydot) = 0;

        virtual void adjustAtol(double *) = 0;

	/**
	 * This function may be used to create output at various points in the
	 * execution of an application.
	 *
	 */
        virtual void user_out(const int ifunc, const double t, const double *y,
                              const double *ydot) = 0;

    protected:

    private:

    };
}
#endif
