/**
 *  @file ImplicitSurfChem.cpp
 *
 * Implicit integration of surface site density equations
 *
 * $Author$
 * $Revision$
 * $Date$
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ImplicitSurfChem.h"
#include "CVode.h"

namespace Cantera {

    ImplicitSurfChem::ImplicitSurfChem(SurfKinetics& kin) 
        : FuncEval(), m_kin(&kin), m_integ(0),
          m_atol(1.e-15), m_rtol(1.e-7), m_maxstep(0.0)
    {
        m_integ = new CVodeInt;
        m_surf = &kin.sphase();

        // use backward differencing, with a full Jacobian computed
        // numerically, and use a Newton linear iterator

        m_integ->setMethod(BDF_Method);
        m_integ->setProblemType(DENSE + NOJAC);
        m_integ->setIterator(Newton_Iter);
        m_nsp = m_surf->nSpecies();
        m_work.resize(m_kin->nTotal());
    }

    // overloaded method of FuncEval. Called by the integrator to
    // get the initial conditions.
    void ImplicitSurfChem::getInitialConditions(double t0, size_t lenc, 
        double* c) 
    {
        m_surf->getCoverages(c);
    }


    /**
     *  Must be called before calling method 'advance'
     */
    void ImplicitSurfChem::initialize(doublereal t0) {
        m_integ->setTolerances(m_rtol, m_atol);
        m_integ->initialize(t0, *this);
    }

    void ImplicitSurfChem::updateState(doublereal* c) {
        m_surf->setCoverages(c);
    }

    /**
     * Called by the integrator to evaluate ydot given y at time 'time'.
     */
    void ImplicitSurfChem::eval(doublereal time, doublereal* y, 
        doublereal* ydot) 
    {
        updateState(y);   // synchronize the surface state with y
        doublereal rs0 = 1.0/m_surf->siteDensity();
        m_kin->getNetProductionRates(m_work.begin());
        int k;
        int kbulk = m_kin->nTotal() - m_nsp;
        doublereal sum = 0.0;
        for (k = 0; k < m_nsp; k++) {
            ydot[k] = m_work[kbulk + k] * rs0 * m_surf->size(k);
        }
    }

}
