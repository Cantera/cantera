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
#include "SolidCompound.h"
#include "SpeciesThermo.h"

namespace Cantera {

    void SolidCompound::initThermo() {
        m_kk = nSpecies();
        if (m_kk > 1) {
            throw CanteraError("initThermo",
                "solid compounds may only contain one species.");
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


    void SolidCompound::_updateThermo() const {
        doublereal tnow = temperature();
        if (m_tlast != tnow) {
            m_spthermo->update(tnow, m_cp0_R.begin(), m_h0_RT.begin(), 
                m_s0_R.begin());
            m_tlast = tnow;
        }
    }
}




