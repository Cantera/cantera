/**
 *
 *  @file SolidTransport.cpp
 */

/* $Author$
 * $Revision$
 * $Date$
 */

// copyright 2003 California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "../ThermoPhase.h"
#include "SolidTransport.h"

#include "utilities.h"
#include <iostream>


namespace Cantera {

    //////////////////// class SolidTransport methods //////////////

    SolidTransport::SolidTransport() {}

    void SolidTransport::setParameters(int n, int k, double* p) {
        switch (n) {
        case 0:
            m_sp.push_back(k);
            m_Adiff.push_back(p[0]);
            m_Ndiff.push_back(p[1]);
            m_Ediff.push_back(p[2]);
            m_nmobile = m_sp.size();
            break;
        case 1:
            m_lam = p[0];
            break;
        default:
            ;
        }
    }


    /*********************************************************
     *
     *                Public methods
     *
     *********************************************************/


    void SolidTransport::getMobilities(doublereal* mobil) {
        int k;
        getMixDiffCoeffs(mobil);
        doublereal t = m_thermo->temperature();
        int nsp = m_thermo->nSpecies();
        doublereal c1 = ElectronCharge / (Boltzmann * t);
        for (k = 0; k < nsp; k++) {
            mobil[k] *= c1 * m_thermo->charge(k);
        }
    } 
        

    doublereal SolidTransport::thermalConductivity() {
        return m_lam;
    }


    void SolidTransport::getMixDiffCoeffs(doublereal* d) {
        doublereal temp = m_thermo->temperature();
        int nsp = m_thermo->nSpecies();
        int k;
        for (k = 0; k < nsp; k++) d[k] = 0.0;
        for (k = 0; k < m_nmobile; k++) {
            d[m_sp[k]] = 
                m_Adiff[k] * pow(temp, m_Ndiff[k]) * exp(-m_Ediff[k]/temp);
        }
    }
}
