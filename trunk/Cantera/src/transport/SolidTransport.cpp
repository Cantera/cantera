/**
 *
 *  @file SolidTransport.cpp
 */

/* $Author: dggoodwin $
 * $Revision: 1.5 $
 * $Date: 2005/11/10 15:06:33 $
 */

// copyright 2003 California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "../ThermoPhase.h"
#include "SolidTransport.h"

#include "../utilities.h"
#include <iostream>


namespace Cantera {

    //////////////////// class SolidTransport methods //////////////

    SolidTransport::SolidTransport() {}

    void SolidTransport::setParameters(int n, int k, double* p) {
        switch (n) {

            // set the Arrhenius parameters for the diffusion coefficient
            // of species k.
        case 0:
            m_sp.push_back(k);
            m_Adiff.push_back(p[0]);
            m_Ndiff.push_back(p[1]);
            m_Ediff.push_back(p[2]);
            m_nmobile = m_sp.size();
            break;

            // set the thermal conductivity.
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


    /**
     * Compute the mobilities of the species from the diffusion coefficients, 
     * using the Einstein relation.
     */
    void SolidTransport::getMobilities(doublereal* mobil) {
        int k;
        getMixDiffCoeffs(mobil);
        doublereal t = m_thermo->temperature();
        int nsp = m_thermo->nSpecies();
        doublereal c1 = ElectronCharge / (Boltzmann * t);
        for (k = 0; k < nsp; k++) {
            mobil[k] *= c1 * fabs(m_thermo->charge(k));
        }
    } 
        

    doublereal SolidTransport::thermalConductivity() {
        return m_lam;
    }


    /**
     * The diffusion coefficients are computed from 
     *
     * \f[
     * D_k = A_k T^{n_k} \exp(-E_k/RT).
     * \f]
     *
     * The diffusion coefficients are only non-zero for species for
     * which parameters have been specified using method
     * setParameters.
     */
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

//     void SolidTransport::electricalConductivity() {
//         getMobilities(m_work.begin());
//         int nsp = m_thermo->nSpecies();
//         int k;
//         doublereal sum = 0.0;
//         for (k = 0; k < nsp; n++) {
//             sum += m_thermo->charge(k)*m_thermo->moleFraction(k)*m_work[k];
//         }
//         return sum * m_thermo->molarDensity();
//     }

}
