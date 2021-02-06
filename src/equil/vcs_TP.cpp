//! @file vcs_TP.cpp

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/equil/vcs_solve.h"
#include "cantera/equil/vcs_VolPhase.h"

namespace Cantera
{
int VCS_SOLVE::vcs_TP(int ipr, int ip1, int maxit, double T_arg, double pres_arg)
{
    // Store the temperature and pressure in the private global variables
    m_temperature = T_arg;
    m_pressurePA = pres_arg;
    m_Faraday_dim = Faraday / (m_temperature * GasConstant);

    // Evaluate the standard state free energies
    // at the current temperatures and pressures.
    int iconv = vcs_evalSS_TP(ipr, ip1, m_temperature, pres_arg);

    // Prep the fe field
    vcs_fePrep_TP();

    // Decide whether we need an initial estimate of the solution If so, go get
    // one. If not, then
    if (m_doEstimateEquil) {
        int retn = vcs_inest_TP();
        if (retn != VCS_SUCCESS) {
            plogf("vcs_inest_TP returned a failure flag\n");
        }
    }

    // Solve the problem at a fixed Temperature and Pressure (all information
    // concerning Temperature and Pressure has already been derived. The free
    // energies are now in dimensionless form.)
    iconv = vcs_solve_TP(ipr, ip1, maxit);

    // Return the convergence success flag.
    return iconv;
}

int VCS_SOLVE::vcs_evalSS_TP(int ipr, int ip1, double Temp, double pres)
{
    for (size_t iph = 0; iph < m_numPhases; iph++) {
        vcs_VolPhase* vph = m_VolPhaseList[iph].get();
        vph->setState_TP(m_temperature, m_pressurePA);
        vph->sendToVCS_GStar(&m_SSfeSpecies[0]);
    }
    for (size_t k = 0; k < m_nsp; k++) {
        m_SSfeSpecies[k] /= GasConstant * m_temperature;
    }

    return VCS_SUCCESS;
}

void VCS_SOLVE::vcs_fePrep_TP()
{
    for (size_t i = 0; i < m_nsp; ++i) {
        // For single species phases, initialize the chemical potential with the
        // value of the standard state chemical potential. This value doesn't
        // change during the calculation
        if (m_SSPhase[i]) {
            m_feSpecies_old[i] = m_SSfeSpecies[i];
            m_feSpecies_new[i] = m_SSfeSpecies[i];
        }
    }
}

}
