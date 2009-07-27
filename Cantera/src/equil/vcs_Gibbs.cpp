/**
 * @file vcs_Gibbs.cpp
 *   Functions which calculate the extrinsic Gibbs Free energies
 */
/*
 * $Id: vcs_Gibbs.cpp,v 1.12 2009/03/27 21:32:33 hkmoffa Exp $
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "vcs_solve.h"
#include "vcs_internal.h" 
#include "vcs_VolPhase.h"
#include "math.h"

namespace VCSnonideal {

  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/

  double VCS_SOLVE::vcs_Total_Gibbs(double *molesSp, double *chemPot,
				    double *tPhMoles)
   
    /*************************************************************************
     *
     * vcs_Total_Gibbs:
     *
     *    Calculate the total dimensionless Gibbs free energy
     *     -> Inert species are handled as if they had a standard free
     *        energy of zero.
     *    Note, for this algorithm this function should be MONOTONICALLY
     *    DECREASING.
     *************************************************************************/
  {
    double g = 0.0;
  
	for (int iph = 0; iph < m_numPhases; iph++) {
		vcs_VolPhase *Vphase = m_VolPhaseList[iph];
		if ((TPhInertMoles[iph] > 0.0) && (tPhMoles[iph] > 0.0)) {
			g += TPhInertMoles[iph] *
				log(TPhInertMoles[iph] / tPhMoles[iph]);
			if (Vphase->m_gasPhase) {
				g += TPhInertMoles[iph] * log(m_pressurePA/(1.01325E5));
			}
		}
	}

	for (int kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
		if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
			g += molesSp[kspec] * chemPot[kspec];
		}
	}
  
    return g;
  }

  // Calculate the total dimensionless Gibbs free energy of a single phase
  /*
   *     -> Inert species are handled as if they had a standard free
   *        energy of zero and if they obeyed ideal solution/gas theory
   *
   * @param iphase   ID of the phase
   * @param w        Species mole number vector
   * @param fe       vector of partial molar free energies of the species.
   */
  double VCS_SOLVE::vcs_GibbsPhase(int iphase, const double * const w,
				   const double * const fe) {
	double g = 0.0;
    double phaseMols = 0.0;
    for (int kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
      if (m_phaseID[kspec] == iphase) {
	    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	      g += w[kspec] * fe[kspec];
	      phaseMols += w[kspec];
	    }
      }
    }

    if (TPhInertMoles[iphase] > 0.0) {
      phaseMols += TPhInertMoles[iphase];
      g += TPhInertMoles[iphase] * log(TPhInertMoles[iphase] / phaseMols);
      vcs_VolPhase *Vphase = m_VolPhaseList[iphase];
      if (Vphase->m_gasPhase) {
        g += TPhInertMoles[iphase] * log(m_pressurePA/1.01325E5);
      }
    }

    return g;
  }

}
