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
	virtual void evalResid(double t, const double* y, 
            const double* ydot, double* resid)=0;

        /** Number of equations. */
        virtual int neq()=0;

        //        virtual bool init(double t0)=0;
        virtual double* solution()=0;
        virtual double* solution_dot()=0;

    protected:

    private:

    };

}

#endif
