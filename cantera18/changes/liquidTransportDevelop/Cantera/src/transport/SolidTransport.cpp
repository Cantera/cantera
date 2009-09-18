/**
 *
 *  @file SolidTransport.cpp
 */

/* $Author: hkmoffa $
 * $Revision: 1.10 $
 * $Date: 2009/03/27 18:24:39 $
 */

// copyright 2008 California Institute of Technology


// turn off warnings under Windows
#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "ThermoPhase.h"
#include "SolidTransport.h"

#include "utilities.h"
#include <iostream>

using namespace std;

namespace Cantera {

    SolidTransport::SolidTransport() {}

    void SolidTransport::setParameters(const int n, const int k, const double* const p) {
        switch (n) {
            
        case 0:
            // set the Arrhenius parameters for the diffusion coefficient
            // of species k.
            m_sp.push_back(k);
            m_Adiff.push_back(p[0]);
            m_Ndiff.push_back(p[1]);
            m_Ediff.push_back(p[2]);
            m_nmobile = m_sp.size();
            break;

        case 1:
            // set the thermal conductivity Arrhenius parameters.
            m_Alam = p[0];
            m_Nlam = p[2];            
            m_Elam = p[2];            
            break;

        default:
            ;
        }
    }


    /**
     * Compute the mobilities of the species from the diffusion coefficients, 
     * using the Einstein relation.
     */
    void SolidTransport::getMobilities(doublereal* const mobil) {
        int k;
        getMixDiffCoeffs(mobil);
        doublereal t = m_thermo->temperature();
        int nsp = m_thermo->nSpecies();
        doublereal c1 = ElectronCharge / (Boltzmann * t);
        for (k = 0; k < nsp; k++) {
            mobil[k] *= c1 * fabs(m_thermo->charge(k));
        }
    } 
        
    /**
     * Thermal Conductivity.
     * \f[
     * \lambda = A T^n \exp(-E/RT)
     */
    doublereal SolidTransport::thermalConductivity() {
        doublereal t = m_thermo->temperature();
        return m_Alam *pow(t, m_Nlam) * exp(-m_Elam/t);
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
    void SolidTransport::getMixDiffCoeffs(doublereal* const d) {
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
