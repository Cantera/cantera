/**
 *  @file vcs_nondim.cpp
 *     Nondimensionalization routines within VCSnonideal
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at http://www.cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctexceptions.h"

namespace Cantera
{

void VCS_SOLVE::vcs_nondim_TP()
{
    if (m_unitsState == VCS_DIMENSIONAL_G) {
        m_unitsState = VCS_NONDIMENSIONAL_G;
        double tf = 1.0 / (GasConstant * m_temperature);
        for (size_t i = 0; i < m_nsp; ++i) {
            // Modify the standard state and total chemical potential data,
            // FF(I), to make it dimensionless, i.e., mu / RT. Thus, we may
            // divide it by the temperature.
            m_SSfeSpecies[i] *= tf;
            m_deltaGRxn_new[i] *= tf;
            m_deltaGRxn_old[i] *= tf;
            m_feSpecies_old[i] *= tf;
        }

        m_Faraday_dim = ElectronCharge * Avogadro / (m_temperature * GasConstant);
    }
}

void VCS_SOLVE::vcs_redim_TP()
{
    if (m_unitsState != VCS_DIMENSIONAL_G) {
        m_unitsState = VCS_DIMENSIONAL_G;
        double tf = m_temperature * GasConstant;
        for (size_t i = 0; i < m_nsp; ++i) {

            // Modify the standard state and total chemical potential data,
            // FF(I), to make it have units, i.e. mu = RT * mu_star
            m_SSfeSpecies[i] *= tf;
            m_deltaGRxn_new[i] *= tf;
            m_deltaGRxn_old[i] *= tf;
            m_feSpecies_old[i] *= tf;
        }
        m_Faraday_dim *= tf;
    }
}

}
