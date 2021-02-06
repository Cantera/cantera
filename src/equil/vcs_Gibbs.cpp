/**
 * @file vcs_Gibbs.cpp
 *   Functions which calculate the extrinsic Gibbs Free energies
 */

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"

namespace Cantera
{
double VCS_SOLVE::vcs_Total_Gibbs(double* molesSp, double* chemPot,
                                  double* tPhMoles)
{
    double g = 0.0;

    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* Vphase = m_VolPhaseList[iph].get();
        if ((TPhInertMoles[iph] > 0.0) && (tPhMoles[iph] > 0.0)) {
            g += TPhInertMoles[iph] *
                 log(TPhInertMoles[iph] / tPhMoles[iph]);
            if (Vphase->m_gasPhase) {
                g += TPhInertMoles[iph] * log(m_pressurePA/(1.01325E5));
            }
        }
    }

    for (size_t kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
        if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            g += molesSp[kspec] * chemPot[kspec];
        }
    }

    return g;
}

double VCS_SOLVE::vcs_GibbsPhase(size_t iphase, const double* const w,
                                 const double* const fe)
{
    double g = 0.0;
    double phaseMols = 0.0;
    for (size_t kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
        if (m_phaseID[kspec] == iphase && m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
            g += w[kspec] * fe[kspec];
            phaseMols += w[kspec];
        }
    }

    if (TPhInertMoles[iphase] > 0.0) {
        phaseMols += TPhInertMoles[iphase];
        g += TPhInertMoles[iphase] * log(TPhInertMoles[iphase] / phaseMols);
        vcs_VolPhase* Vphase = m_VolPhaseList[iphase].get();
        if (Vphase->m_gasPhase) {
            g += TPhInertMoles[iphase] * log(m_pressurePA/1.01325E5);
        }
    }

    return g;
}

}
