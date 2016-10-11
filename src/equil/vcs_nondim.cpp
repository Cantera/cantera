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
        for (size_t i = 0; i < m_numSpeciesTot; ++i) {
            // Modify the standard state and total chemical potential data,
            // FF(I), to make it dimensionless, i.e., mu / RT. Thus, we may
            // divide it by the temperature.
            m_SSfeSpecies[i] *= tf;
            m_deltaGRxn_new[i] *= tf;
            m_deltaGRxn_old[i] *= tf;
            m_feSpecies_old[i] *= tf;
        }

        m_Faraday_dim = ElectronCharge * Avogadro / (m_temperature * GasConstant);

        // Scale the total moles if necessary: First find out the total moles
        double tmole_orig = vcs_tmoles();

        // Then add in the total moles of elements that are goals. Either one or
        // the other is specified here.
        double esum = 0.0;
        for (size_t i = 0; i < m_numElemConstraints; ++i) {
            if (m_elType[i] == VCS_ELEM_TYPE_ABSPOS) {
                esum += fabs(m_elemAbundancesGoal[i]);
            }
        }
        tmole_orig += esum;

        // Ok now test out the bounds on the total moles that this program can
        // handle. These are a bit arbitrary. However, it would seem that any
        // reasonable input would be between these two numbers below.
        if (tmole_orig < 1.0E-200 || tmole_orig > 1.0E200) {
            throw CanteraError("VCS_SOLVE::vcs_nondim_TP",
                "Total input moles, {} is outside the range handled by vcs.\n",
                tmole_orig);
        }

        // Determine the scale of the problem
        if (tmole_orig > 1.0E4) {
            m_totalMoleScale = tmole_orig / 1.0E4;
        } else if (tmole_orig < 1.0E-4) {
            m_totalMoleScale = tmole_orig / 1.0E-4;
        } else {
            m_totalMoleScale = 1.0;
        }

        if (m_totalMoleScale != 1.0) {
            if (m_debug_print_lvl >= 2) {
                plogf("  --- vcs_nondim_TP() called: USING A MOLE SCALE OF %g until further notice", m_totalMoleScale);
                plogendl();
            }
            for (size_t i = 0; i < m_numSpeciesTot; ++i) {
                if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                    m_molNumSpecies_old[i] *= (1.0 / m_totalMoleScale);
                }
            }
            for (size_t i = 0; i < m_numElemConstraints; ++i) {
                m_elemAbundancesGoal[i] *= (1.0 / m_totalMoleScale);
            }

            for (size_t iph = 0; iph < m_numPhases; iph++) {
                TPhInertMoles[iph] *= (1.0 / m_totalMoleScale);
                if (TPhInertMoles[iph] != 0.0) {
                    vcs_VolPhase* vphase = m_VolPhaseList[iph];
                    vphase->setTotalMolesInert(TPhInertMoles[iph]);
                }
            }
            vcs_tmoles();
        }
    }
}

void VCS_SOLVE::vcs_redim_TP()
{
    if (m_unitsState != VCS_DIMENSIONAL_G) {
        m_unitsState = VCS_DIMENSIONAL_G;
        double tf = m_temperature * GasConstant;
        for (size_t i = 0; i < m_numSpeciesTot; ++i) {

            // Modify the standard state and total chemical potential data,
            // FF(I), to make it have units, i.e. mu = RT * mu_star
            m_SSfeSpecies[i] *= tf;
            m_deltaGRxn_new[i] *= tf;
            m_deltaGRxn_old[i] *= tf;
            m_feSpecies_old[i] *= tf;
        }
        m_Faraday_dim *= tf;
    }
    if (m_totalMoleScale != 1.0) {
        if (m_debug_print_lvl >= 2) {
            plogf("  --- vcs_redim_TP() called: getting rid of mole scale of %g", m_totalMoleScale);
            plogendl();
        }
        for (size_t i = 0; i < m_numSpeciesTot; ++i) {
            if (m_speciesUnknownType[i] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
                m_molNumSpecies_old[i] *= m_totalMoleScale;
            }
        }
        for (size_t i = 0; i < m_numElemConstraints; ++i) {
            m_elemAbundancesGoal[i] *= m_totalMoleScale;
        }

        for (size_t iph = 0; iph < m_numPhases; iph++) {
            TPhInertMoles[iph] *= m_totalMoleScale;
            if (TPhInertMoles[iph] != 0.0) {
                vcs_VolPhase* vphase = m_VolPhaseList[iph];
                vphase->setTotalMolesInert(TPhInertMoles[iph]);
            }
        }
        vcs_tmoles();
    }
}

}
