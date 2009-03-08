/**
 *
 *  @file StoichSubstance.cpp
 *
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ct_defs.h"
#include "mix_defs.h"
#include "StoichSubstance.h"
#include "SpeciesThermo.h"

namespace Cantera {

    void StoichSubstance::initThermo() {
        m_kk = nSpecies();
        if (m_kk > 1) {
            throw CanteraError("initThermo",
                "stoichiometric substances may only contain one species.");
        } 
        doublereal tmin = m_spthermo->minTemp();
        doublereal tmax = m_spthermo->maxTemp();
        if (tmin > 0.0) m_tmin = tmin;
        if (tmax > 0.0) m_tmax = tmax;
        m_p0 = refPressure();

        int leng = m_kk;
        m_h0_RT.resize(leng);
        m_cp0_R.resize(leng);
        m_s0_R.resize(leng);
    }


    void StoichSubstance::_updateThermo() const {
        doublereal tnow = temperature();
        if (m_tlast != tnow) {
            m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0], 
                &m_s0_R[0]);
            m_tlast = tnow;
        }
    }

    void StoichSubstance::
    getUnitsStandardConc(double *uA, int k, int sizeUA) {
	for (int i = 0; i < sizeUA; i++) {
	  uA[i] = 0.0;
	}
    }

    void StoichSubstance::setParameters(int n, double * c) {
        double rho = c[0];
        setDensity(rho);
    }

    void StoichSubstance::getParameters(int &n, double * const c) {
        double rho = density();
        c[0] = rho;
    }

    void StoichSubstance::setParametersFromXML(const XML_Node& eosdata) {
        eosdata._require("model","StoichSubstance");
        doublereal rho = getFloat(eosdata, "density", "-");
        setDensity(rho);
    }

}




