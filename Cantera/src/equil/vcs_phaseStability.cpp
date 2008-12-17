/* ======================================================================= */
/* -------------------------------------------------- */
/* | RCS Head Information on zuzax.pchem.sandia.gov | */
/* -------------------------------------------------- */
/* $RCSfile$ */
/* $Author$ */
/* $Date$ */
/* $Revision$ */
/* ======================================================================= */


#include "vcs_solve.h"
#include "vcs_internal.h" 
#include "vcs_species_thermo.h"
#include "vcs_VolPhase.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>

using namespace std;

namespace VCSnonideal {
  
  int VCS_SOLVE::vcs_phaseStabilityTest(const int iph) {

    /*
     * We will use the _new state calc here
     */
    int kspec, irxn, k;
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    double deltaGRxn;
    // We will do a full newton calculation later, but for now, ...
    bool doSuccessiveSubstitution = true;
    int res =  0;

    vector<double> X_est(Vphase->nSpecies(), 0.0);
    vector<double> X_est_old(Vphase->nSpecies(), 0.0);
    vector<double> delX(Vphase->nSpecies(), 0.0);
    vector<double> E_phi(Vphase->nSpecies(), 0.0);
    double damp = 1.0;
    double normUpdate = 1.0;
    double normUpdateOld = 1.0;

    // get the activity coefficients
    Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(m_actCoeffSpecies_new));
    
    
    if (doSuccessiveSubstitution) {

      for (int its = 0; its < 20; its++) {


	normUpdateOld = normUpdate;
	for (k = 0; k < Vphase->nSpecies(); k++) {
	  X_est_old[k] = X_est[k];
	}

	double poly = -1.0;
	for (k = 0; k < Vphase->nSpecies(); k++) {
	  kspec = Vphase->spGlobalIndexVCS(k);
	  irxn = kspec - m_numComponents;
	  deltaGRxn = m_deltaGRxn_old[irxn];
	  // We may need to look at deltaGRxn for components!
	  if (irxn >= 0) {
	    if (deltaGRxn >  50.0) deltaGRxn =  50.0;
	    if (deltaGRxn < -50.0) deltaGRxn = -50.0;
	    E_phi[k] = exp(-deltaGRxn)/m_actCoeffSpecies_new[kspec];
	    poly += E_phi[k];
	  }
	}
	double sum = poly + 1.0;

	for (k = 0; k <  Vphase->nSpecies(); k++) {
	  delX[k] = E_phi[k]/sum - X_est_old[k];
	}
	normUpdate = vcs_l2norm(delX);

	// Figure out the damping coefficient
	double ratio = normUpdate / normUpdateOld;
	if (ratio < 0.4) {
	  damp = 1.0;
	} else if (ratio > 1.0) {
	  damp = 0.03;
	} else {
	  damp = 0.1;
	}

	
	for (k = 0; k <  Vphase->nSpecies(); k++) {
	  X_est[k] = X_est_old[k] + damp * delX[k];
	}

	for (k = 0; k < Vphase->nSpecies(); k++) {
	  kspec = Vphase->spGlobalIndexVCS(k);
	  m_molNumSpecies_new[kspec] = X_est[k];
	}
	Vphase->setMolesFromVCS(VCS_STATECALC_NEW);
	

      }
    } else {
      printf("not done yet\n");
      exit(-1);
    }
    
    return res;
  }

}

