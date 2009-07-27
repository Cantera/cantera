/**
 *  @file ImplicitChem.cpp
 */

/* $Author: dggoodwin $
 * $Revision: 1.1 $
 * $Date: 2007/05/04 14:27:23 $
 */

// Copyright 2001  California Institute of Technology


#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ImplicitChem.h"
#include "Integrator.h"

namespace Cantera {

    ImplicitChem::ImplicitChem(Kinetics& kin, ThermoPhase& therm) 
        : FuncEval(), m_kin(&kin), m_thermo(&therm), m_integ(0),
          m_atol(1.e-15), m_rtol(1.e-7), m_maxstep(0.0), m_energy(false)
    {
        m_integ = newIntegrator("CVODE"); //CVodeInt;
        //m_mix = &kin.phase();
        m_wt = m_thermo->molecularWeights();

        // use backward differencing, with a full Jacobian computed
        // numerically, and use a Newton linear iterator
        m_integ->setMethod(BDF_Method);
        m_integ->setProblemType(DENSE + NOJAC);
        m_integ->setIterator(Newton_Iter);
        m_nsp = m_thermo->nSpecies();
    }

    // overloaded method of FuncEval. Called by the integrator to
    // get the initial conditions.
    void ImplicitChem::getInitialConditions(double t0, size_t leny, double* y) 
    {
        m_thermo->getMassFractions(y);
        m_h0 = m_thermo->enthalpy_mass();
        m_rho = m_thermo->density();
        m_press = m_thermo->pressure();
    }


    /**
     *  Must be called before calling method 'advance'
     */
    void ImplicitChem::initialize(doublereal t0) {
        m_integ->setTolerances(m_rtol, m_atol);
        //        m_integ->setMaxStep(m_maxstep);
        m_integ->initialize(t0, *this);
    }

    
    void ImplicitChem::updateState(doublereal* y) {
        m_thermo->setMassFractions(y);
        if (m_energy) {
            doublereal delta, temp = m_thermo->temperature();
            do {
                delta = -(m_thermo->enthalpy_mass() - m_h0)/m_thermo->cp_mass();
                temp += delta;
                m_thermo->setTemperature(temp);
            }
            while (fabs(delta) > 1.e-7);
        }
        m_thermo->setPressure(m_press);
    }

    /**
     * Called by the integrator to evaluate ydot given y at time 'time'.
     */
    void ImplicitChem::eval(doublereal time, doublereal* y, 
        doublereal* ydot, doublereal* p) 
    {
        updateState(y);   // synchronize the mixture state with y
        m_thermo->setPressure(m_press);
        m_kin->getNetProductionRates(ydot);   // "omega dot"
        int k;
        for (k = 0; k < m_nsp; k++) {
            ydot[k] *= m_wt[k]/m_rho;
        }
    }

}
