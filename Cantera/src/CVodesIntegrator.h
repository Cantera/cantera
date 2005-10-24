/**
 *  @file CVodesWrapper.h
 */

/* $Author$
 * $Date$
 * $Revision$
 */

// Copyright 2005  California Institute of Technology


#ifndef CT_CVODESWRAPPER_H
#define CT_CVODESWRAPPER_H

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "Integrator.h"
#include "FuncEval.h"
#include "ctexceptions.h"
#include "ct_defs.h"

namespace Cantera {

    /**
     * Exception thrown when a CVODES error is encountered.
     */
    class CVodesErr : public CanteraError {
    public:
        CVodesErr(string msg) : CanteraError("CVodesIntegrator", msg){}
    };


    /**
     *  Wrapper class for 'cvodes' integrator from LLNL.
     *
     * @see FuncEval.h. Classes that use CVodeInt:
     * ImplicitChem, ImplicitSurfChem, Reactor
     *
     */
    class CVodesIntegrator : public Integrator {

    public:

        CVodesIntegrator();
        virtual ~CVodesIntegrator();
        virtual void setTolerances(double reltol, int n, double* abstol);
        virtual void setTolerances(double reltol, double abstol);
        virtual void setProblemType(int probtype);
        virtual void initialize(double t0, FuncEval& func);
        virtual void reinitialize(double t0, FuncEval& func);
        virtual void integrate(double tout);
        virtual doublereal step(double tout);
	virtual double& solution(int k);
	virtual double* solution();
	virtual int nEquations() const { return m_neq;}
        virtual int nEvals() const;
        virtual void setMaxOrder(int n) { m_maxord = n; }
        virtual void setMethod(MethodType t);
        virtual void setIterator(IterType t);
        virtual void setMaxStepSize(double hmax);
        virtual void setMinStepSize(double hmin);
        virtual void setMaxSteps(int nmax);

    private:

	int m_neq;
        void* m_cvode_mem;
        double m_t0;
        void *m_y, *m_abstol;
        int m_type;
        int m_itol;
        int m_method;
        int m_iter;
        int m_maxord;
        double m_reltol;
        double m_abstols;
        int m_nabs;
        double m_hmax, m_hmin;
        int m_maxsteps;
        void* m_data;
    };

}    // namespace

#endif 
