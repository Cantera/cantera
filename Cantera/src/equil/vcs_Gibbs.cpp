/**
 * @file vcs_Gibbs.cpp
 *   Functions which calculate the extrinsic Gibbs Free energies
 */
/*
 * $Id$
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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

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
  
    for (int iph = 0; iph < NPhase; iph++) {
      vcs_VolPhase *Vphase = VPhaseList[iph];
      if ((TPhInertMoles[iph] > 0.0) && (tPhMoles[iph] > 0.0)) {
	g += TPhInertMoles[iph] *
	    log(TPhInertMoles[iph] / tPhMoles[iph]);
	if (Vphase->GasPhase) {
	  g += TPhInertMoles[iph] * log(Pres);
	}
      }
    }

    for (int kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
      g += molesSp[kspec] * chemPot[kspec];
    }
  
    return g;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

double VCS_SOLVE::vcs_GibbsPhase(int iphase, double *w, double *fe)
   
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
 
    vcs_VolPhase *Vphase = VPhaseList[iphase];
    if ((TPhInertMoles[iphase] > 0.0) && (TPhMoles[iphase] > 0.0)) {
      g += TPhInertMoles[iphase] *
	  log(TPhInertMoles[iphase] / TPhMoles[iphase]);
      if (Vphase->GasPhase == iphase) {
	g += TPhInertMoles[iphase] * log(Pres);
      }
    }
    
    for (int kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
      if (PhaseID[kspec] == iphase) {
	g += w[kspec] * fe[kspec];
      }
    }
  
    return g;
}

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
}


