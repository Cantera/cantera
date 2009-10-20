/**
 *  @file ConstDensityThermo.cpp
 * Declarations for a Thermo manager for incompressible ThermoPhases 
 * (see \ref thermoprops and \link Cantera::ConstDensityThermo ConstDensityThermo
\endlink).
 */
/*
 * $Id: ConstDensityThermo.cpp,v 1.5 2008/12/17 17:04:47 hkmoffa Exp $
 *
 *  Copyright 2002 California Institute of Technology
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "mix_defs.h"
#include "ConstDensityThermo.h"
#include "SpeciesThermo.h"

#include <cmath>

namespace Cantera {

  ConstDensityThermo::ConstDensityThermo() : m_tlast(0.0) {
  }


  ConstDensityThermo::ConstDensityThermo(const ConstDensityThermo &right)
    : m_tlast(0.0) {
    *this = operator=(right);
  }

  ConstDensityThermo& ConstDensityThermo::operator=(const ConstDensityThermo &right) {
    if (&right == this) return *this;

    m_mm            = right.m_mm;
    m_tmin          = right.m_tmin;
    m_tmax          = right.m_tmax;
    m_p0            = right.m_p0;
    m_tlast         = right.m_tlast;
    m_h0_RT         = right.m_h0_RT;
    m_cp0_R         = right.m_cp0_R;
    m_g0_RT         = right.m_g0_RT;
    m_s0_R          = right.m_s0_R;
    m_expg0_RT      = right.m_expg0_RT;
    m_pe            = right.m_pe;
    m_pp            = right.m_pp;
    
    return *this;

  }

  SpeciesThermo *ConstDensityThermo::duplMyselfAsSpeciesThermo() const {
     ConstDensityThermo *cdt = new ConstDensityThermo(*this);
     return (SpeciesThermo *) cdt;
  }
    int ConstDensityThermo::
    eosType() const { return cIncompressible; }

    doublereal ConstDensityThermo::enthalpy_mole() const {
        doublereal p0 = m_spthermo->refPressure();
        return GasConstant * temperature() * 
            mean_X(&enthalpy_RT()[0]) 
            + (pressure() - p0)/molarDensity();
    }

    doublereal ConstDensityThermo::intEnergy_mole() const {
        doublereal p0 = m_spthermo->refPressure();
        return GasConstant * temperature() * 
            mean_X(&enthalpy_RT()[0]) 
            - p0/molarDensity();
    }

    doublereal ConstDensityThermo::entropy_mole() const {
        return GasConstant * (mean_X(&entropy_R()[0]) -
            sum_xlogx());
    }

    doublereal ConstDensityThermo::gibbs_mole() const {
        return enthalpy_mole() - temperature() * entropy_mole();
    }

    doublereal ConstDensityThermo::cp_mole() const {
        return GasConstant * mean_X(&cp_R()[0]);
    }

    doublereal ConstDensityThermo::cv_mole() const {
        return cp_mole();
    }

    doublereal ConstDensityThermo::pressure() const {
        return m_press;
    }

    void ConstDensityThermo::setPressure(doublereal p) {
        m_press = p;
    }

    void ConstDensityThermo::getActivityConcentrations(doublereal* c) const {
        getConcentrations(c);
    }

    void ConstDensityThermo::getActivityCoefficients(doublereal* ac) const {
        for (int k = 0; k < m_kk; k++) {
 	  ac[k] = 1.0;
	}
    }

    doublereal ConstDensityThermo::standardConcentration(int k) const {
        return molarDensity();
    }

    doublereal ConstDensityThermo::logStandardConc(int k) const {
        return log(molarDensity());
    }

    void ConstDensityThermo::getChemPotentials(doublereal* mu) const {
        doublereal vdp = (pressure() - m_spthermo->refPressure())/
                         molarDensity();
        doublereal xx;
        doublereal rt = temperature() * GasConstant;
        const array_fp& g_RT = gibbs_RT();
        for (int k = 0; k < m_kk; k++) {
            xx = fmaxx(SmallNumber, moleFraction(k));
            mu[k] = rt*(g_RT[k] + log(xx)) + vdp;
        }
    }


    void ConstDensityThermo::getStandardChemPotentials(doublereal* mu0) const {
        getPureGibbs(mu0);
    }

    void ConstDensityThermo::initThermo() {
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


    void ConstDensityThermo::setToEquilState(const doublereal* lambda_RT) {
        throw CanteraError("setToEquilState","not yet impl.");
    }

    void ConstDensityThermo::_updateThermo() const {
        doublereal tnow = temperature();
        if (m_tlast != tnow) {
            m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], 
                &m_s0_R[0]);
            m_tlast = tnow;
            int k;
            for (k = 0; k < m_kk; k++) {
                m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
            }
            m_tlast = tnow;
        }
    }

    void ConstDensityThermo::setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","Incompressible");
        doublereal rho = getFloat(eosdata, "density", "toSI");
        setDensity(rho);
    }

}




