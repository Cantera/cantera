/**
 *  @file UnityLewisNumberTransport.cpp
 *  Mixture-averaged transport properties for ideal gas mixtures with
 *  diffusion coefficients set based on the assumption that Lewis #  = 1.
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/transport/UnityLewisNumberTransport.h"
#include "cantera/base/stringUtils.h"

using namespace std;

namespace Cantera
{

UnityLewisTransport::UnityLewisTransport()
{
}

void UnityLewisTransport::init(ThermoPhase* thermo, int mode, int log_level)
{
    MixTransport::init(thermo, mode, log_level);
}


void UnityLewisTransport::getMixDiffCoeffs(doublereal* const d)
{
    MixTransport::update_T();
    MixTransport::update_C();
    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    //Unity Lewis number Le_k = alpha/Dm_k. If Le_k = 1, then
    //Dm_k = alpha = lambda / (rho * cp)
    double rho = m_thermo->density();
    double cp = m_thermo->cp_mass();
    double lambda = MixTransport::thermalConductivity();
    double thermal_diffusivity = lambda / (rho * cp);
    for (size_t k = 0; k < m_nsp; k++) {
        d[k] = thermal_diffusivity;
    }
}


void UnityLewisTransport::getMixDiffCoeffsMole(doublereal* const d)
{
    update_T();
    update_C();
    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    //Unity Lewis number Le_k = alpha/Dm_k. If Le_k = 1, then
    //Dm_k = alpha = lambda / (rho * cp)
    double rho = m_thermo->density();
    double cp = m_thermo->cp_mass();
    double lambda = MixTransport::thermalConductivity();
    double thermal_diffusivity = lambda / (rho * cp);
    for (size_t k = 0; k < m_nsp; k++) {
        d[k] = thermal_diffusivity; 
    }

}


void UnityLewisTransport::getMixDiffCoeffsMass(doublereal* const d)
{
    update_T();
    update_C();
    // update the binary diffusion coefficients if necessary
    if (!m_bindiff_ok) {
        updateDiff_T();
    }

    //Unity Lewis number Le_k = alpha/Dm_k. If Le_k = 1, then
    //Dm_k = alpha = lambda / (rho * cp)
    double rho = m_thermo->density();
    double cp = m_thermo->cp_mass();
    double lambda = MixTransport::thermalConductivity();
    double thermal_diffusivity = lambda / (rho * cp);
    for (size_t k = 0; k < m_nsp; k++) {
        d[k] = thermal_diffusivity; 
    }


}




}
