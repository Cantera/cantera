/**
 *  @file ImplicitSurfChem.cpp
 *
 * Implicit integration of surface site density equations
 *
 * $Author: dggoodwin $
 * $Revision: 1.10 $
 * $Date: 2006/04/28 17:22:23 $
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ImplicitSurfChem.h"

#include "Integrator.h"

namespace Cantera {

    ImplicitSurfChem::ImplicitSurfChem(vector<InterfaceKinetics*> k) 
        : FuncEval(),  m_nv(0), m_integ(0), 
          m_atol(1.e-14), m_rtol(1.e-7), m_maxstep(0.0)
    {
        m_nsurf = static_cast<int>(k.size());
        int ns;
        int nt, ntmax = 0;
        for (int n = 0; n < m_nsurf; n++) {
            m_kin.push_back(k[n]);
            ns = k[n]->surfacePhaseIndex();
            if (ns < 0) 
                throw CanteraError("ImplicitSurfChem",
                    "kinetics manager contains no surface phase");
            m_surfindex.push_back(ns);
            m_surf.push_back((SurfPhase*)&k[n]->thermo(ns));
            m_nsp.push_back(m_surf.back()->nSpecies());
            m_nv += m_nsp.back();
            nt = k[n]->nTotalSpecies();
            if (nt > ntmax) ntmax = nt;
        }
        m_integ = newIntegrator("CVODE");// CVodeInt;

        // use backward differencing, with a full Jacobian computed
        // numerically, and use a Newton linear iterator

        m_integ->setMethod(BDF_Method);
        m_integ->setProblemType(DENSE + NOJAC);
        m_integ->setIterator(Newton_Iter);
        m_work.resize(ntmax);
    }


    // overloaded method of FuncEval. Called by the integrator to
    // get the initial conditions.
    void ImplicitSurfChem::getInitialConditions(double t0, size_t lenc, 
        double* c) 
    {
        int loc = 0;
        for (int n = 0; n < m_nsurf; n++) {
            m_surf[n]->getCoverages(c + loc);
            loc += m_nsp[n];
        }
    }


    /**
     *  Must be called before calling method 'advance'
     */
    void ImplicitSurfChem::initialize(doublereal t0) {
        m_integ->setTolerances(m_rtol, m_atol);
        m_integ->initialize(t0, *this);
    }


    void ImplicitSurfChem::updateState(doublereal* c) {
        int loc = 0;
        for (int n = 0; n < m_nsurf; n++) {
            m_surf[n]->setCoverages(c + loc);
            loc += m_nsp[n];
        }
    }


    /**
     * Called by the integrator to evaluate ydot given y at time 'time'.
     */
    void ImplicitSurfChem::eval(doublereal time, doublereal* y, 
        doublereal* ydot, doublereal* p) 
    {
        int n;
        updateState(y);   // synchronize the surface state(s) with y
        doublereal rs0, sum;
        int loc, k, kstart;
        for (n = 0; n < m_nsurf; n++) {
            rs0 = 1.0/m_surf[n]->siteDensity();
            m_kin[n]->getNetProductionRates(DATA_PTR(m_work));
            kstart = m_kin[n]->kineticsSpeciesIndex(0,m_surfindex[n]);
            sum = 0.0;
            loc = 0;
            for (k = 1; k < m_nsp[n]; k++) {
                ydot[k + loc] = m_work[kstart + k] * rs0 * m_surf[n]->size(k);
                sum -= ydot[k];
            }
            ydot[loc] = sum;
            loc += m_nsp[n];
        }
    }

}
