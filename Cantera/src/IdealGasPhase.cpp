/**
 *
 *  @file IdealGasPhase.cpp
 *
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "mix_defs.h"
#include "IdealGasPhase.h"
#include "SpeciesThermo.h"

namespace Cantera {

    void IdealGasPhase::getChemPotentials(doublereal* mu) const {
        doublereal logp = log(pressure()/m_spthermo->refPressure());
        doublereal xx;
        doublereal rt = temperature() * GasConstant;
        const array_fp& g_RT = gibbs_RT();
        for (int k = 0; k < m_kk; k++) {
            xx = fmaxx(SmallNumber, moleFraction(k));
            mu[k] = rt*(g_RT[k] + log(xx) + logp);
        }
    }

    // new methods defined here


    void IdealGasPhase::initThermo() {
        m_kk = nSpecies();
        m_mm = nElements();
        doublereal tmin = m_spthermo->minTemp();
        doublereal tmax = m_spthermo->maxTemp();
        if (tmin > 0.0) m_tmin = tmin;
        if (tmax > 0.0) m_tmax = tmax;
        m_p0 = refPressure();

        int leng = m_kk;
        m_h0_RT.resize(leng);
        m_g0_RT.resize(leng);
        m_expg0_RT.resize(leng);
        m_cp0_R.resize(leng);
        m_s0_R.resize(leng);
        m_pe.resize(leng, 0.0);
        m_pp.resize(leng);
    }


    /** 
     * Set mixture to an equilibrium state consistent with specified 
     * element potentials and temperature.
     *
     * @param lambda_RT vector of non-dimensional element potentials
     * \f[ \lambda_m/RT \f].
     * @param t temperature in K.
     * @param work. Temporary work space. Must be dimensioned at least
     * as large as the number of species. 
     *
     */
    void IdealGasPhase::setToEquilState(const doublereal* lambda_RT) 
    {
        const array_fp& grt = gibbs_RT();

        // set the pressure and composition to be consistent with
        // the temperature, 
        doublereal pres = 0.0;
        for (int k = 0; k < m_kk; k++) {
            m_pp[k] = -grt[k];
            for (int m = 0; m < m_mm; m++) {
                m_pp[k] += nAtoms(k,m)*lambda_RT[m];
                //cout << "m = " << m << "  k = " << k << "  " << 
                //    nAtoms(k,m) << "  " << lambda_RT[m] << endl;
            } 
            //cout << "m_pp = " << m_pp[k] << endl;
            m_pp[k] = m_p0 * exp(m_pp[k]);
            //cout << "after exp, m_pp[k] = " << m_pp[k] << endl;
            pres += m_pp[k];
        }
        //cout << "pres = " << pres << endl;
        // set state
        setState_PX(pres, m_pp.begin());
    }

    void IdealGasPhase::_updateThermo() const {
        doublereal tnow = temperature();
        if (m_tlast != tnow) {
            m_spthermo->update(tnow, m_cp0_R.begin(), m_h0_RT.begin(), 
                m_s0_R.begin());
            m_tlast = tnow;
            //            doublereal rrt = 1.0 / (GasConstant * tnow);
            int k;
            //doublereal deltaE;
            for (k = 0; k < m_kk; k++) {
                //deltaE = rrt * m_pe[k];
                //m_h0_RT[k] += deltaE;
                m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
            }
            m_logc0 = log(m_p0/(GasConstant * tnow));
            m_tlast = tnow;
        }
    }
}




