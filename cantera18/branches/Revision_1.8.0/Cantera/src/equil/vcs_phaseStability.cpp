/* ======================================================================= */
/* $RCSfile: vcs_phaseStability.cpp,v $ */
/* $Date: 2009/04/18 00:42:51 $ */
/* $Revision: 1.8 $ */
/* ======================================================================= */


#include "vcs_solve.h"
#include "vcs_internal.h" 
#include "vcs_species_thermo.h"
#include "vcs_VolPhase.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstring>

using namespace std;

namespace VCSnonideal {


  // Utility function that evaluates whether a phase can be popped
  // into existence
  /*
   * A phase can be popped iff the stoichiometric coefficients for the
   * component species, whose concentrations will be lowered during the
   * process, are positive by at least a small degree.
   * 
   * @param iphasePop  id of the phase, which is currently zeroed,
   *        
   * @return Returns true if the phase can come into existence
   *         and false otherwise.
   */
  bool VCS_SOLVE::vcs_popPhasePossible(const int iphasePop) const {
    
    vcs_VolPhase *Vphase = m_VolPhaseList[iphasePop];

#ifdef DEBUG_MODE
    int existence = Vphase->exists();
    if (existence > 0) {
      printf("ERROR vcs_popPhasePossible called for a phase that exists!");
      std::exit(-1);
    }
#endif

    /*
     * Loop through all of the species in the phase. We say the phase
     * can be popped, if there is one species in the phase that can be
     * popped.
     */
    for (int k = 0; k < Vphase->nSpecies(); k++) {
      int kspec = Vphase->spGlobalIndexVCS(k);
      int irxn = kspec - m_numComponents;
      if (irxn >= 0) {
	int iPopPossible = true;
	for (int j = 0; j < m_numComponents; ++j) {
	  if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
	    double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
	    if (stoicC != 0.0) {    
	      double negChangeComp = - stoicC * 1.0;
	      if (negChangeComp > 0.0) {
		// TODO: We may have to come up with a tolerance here
		if (m_molNumSpecies_old[j] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.1) {
		  iPopPossible = false;
		}
	      }
	    }
	  }
	}
	if (iPopPossible == true) {
	  return true;
	}
      }
    }
    return false;
  }

  // Decision as to whether a phase pops back into existence
  /*
   * @return returns the phase id of the phase that pops back into 
   *         existence. Returns -1 if there are no phases
   */
  int VCS_SOLVE::vcs_popPhaseID() {
    int iphasePop = -1;
    int iph;
    int irxn, kspec;
    doublereal FephaseMax = -1.0E30;
    doublereal Fephase = -1.0E30;
    vcs_VolPhase *Vphase = 0;

  
#ifdef DEBUG_MODE
    char anote[128];
    if (m_debug_print_lvl >= 2) {
      plogf("   --- vcs_popPhaseID() called\n");
      plogf("   ---   Phase               Status      F_e         MoleNum\n");
      plogf("   --------------------------------------------------------------\n");
    }
#endif
    for (iph = 0; iph < m_numPhases; iph++) {
      Vphase = m_VolPhaseList[iph];
      int existence = Vphase->exists();
#ifdef DEBUG_MODE
      strcpy(anote, "");
#endif
      if (existence > 0) {
	
#ifdef DEBUG_MODE
	if (m_debug_print_lvl >= 2) {
	  plogf("  ---    %18s %5d         NA     %11.3e\n", 
		Vphase->PhaseName.c_str(),
		existence,
		m_tPhaseMoles_old[iph]); 
	}
#endif
      } else {
	if (Vphase->m_singleSpecies) {
	  /***********************************************************************
	   *
	   *  Single Phase Stability Resolution
	   *
	   ***********************************************************************/
	  kspec = Vphase->spGlobalIndexVCS(0);
	  irxn = kspec - m_numComponents;
	  doublereal deltaGRxn = m_deltaGRxn_old[irxn];
	  Fephase = exp(-deltaGRxn) - 1.0;
	  if (Fephase > 0.0) {
#ifdef DEBUG_MODE
	    strcpy(anote," (ready to be birthed)");
#endif
	    if (Fephase > FephaseMax) {
	      iphasePop = iph;
	      FephaseMax = Fephase;
#ifdef DEBUG_MODE
	      strcpy(anote," (chosen to be birthed)");
#endif
	    }
	  }

#ifdef DEBUG_MODE
	  if (m_debug_print_lvl >= 2) {
	    plogf("  ---    %18s %5d         NA     %11.3g %s\n", 
		  Vphase->PhaseName.c_str(),
		  existence, Fephase,
		  m_tPhaseMoles_old[iph], anote); 
	  }
#endif

	} else {
	  /***********************************************************************
	   *
	   * MultiSpecies Phase Stability Resolution
	   *
	   ***********************************************************************/
	  if (vcs_popPhasePossible(iph)) {
	    Fephase = vcs_phaseStabilityTest(iph);
	    if (Fephase > 0.0) {
	      if (Fephase > FephaseMax) {
		iphasePop = iph;
		FephaseMax = Fephase;
	      }
	    } else {
	      if (Fephase > FephaseMax) {
		FephaseMax = Fephase;
	      }
	    }
#ifdef DEBUG_MODE
	    if (m_debug_print_lvl >= 2) {
	      plogf("  ---    %18s %5d  %11.3g %11.3g\n",
		    Vphase->PhaseName.c_str(),
		    existence, Fephase,
		  m_tPhaseMoles_old[iph]); 
	    }
#endif
	  } else {
#ifdef DEBUG_MODE
	    if (m_debug_print_lvl >= 2) {
	      plogf("  ---    %18s %5d   blocked  %11.3g\n",
		    Vphase->PhaseName.c_str(),
		    existence, m_tPhaseMoles_old[iph]); 
	    }
#endif
	  }
	}
      }
    }


#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
      plogf("   --------------------------------------------------------------\n");
    }
#endif
    return iphasePop;
  }

  // Calculates the deltas of the reactions due to phases popping
  // into existence
  /*
   * @param iphasePop  Phase id of the phase that will come into existence
   *
   * Output 
   * ------- 
   * m_deltaMolNumSpecies(irxn) : reaction adjustments, where irxn refers 
   *                              to the irxn'th species
   *                              formation reaction. This  adjustment 
   *                              is for species
   *                               irxn + M, where M is the number
   *                              of components.
   *
   * @return  Returns an int representing the status of the step
   *            -  0 : normal return
   *            -  1 : A single species phase species has been zeroed out
   *                   in this routine. The species is a noncomponent 
   *            -  2 : Same as one but, the zeroed species is a component. 
   *            -  3 : Nothing was done because the phase couldn't be birthed
   *                   because a needed component is zero.
   */
  int VCS_SOLVE::vcs_popPhaseRxnStepSizes(const int iphasePop) {
    vcs_VolPhase *Vphase = m_VolPhaseList[iphasePop];
    // Identify the first species in the phase
    int kspec = Vphase->spGlobalIndexVCS(0);
    // Identify the formation reaction for that species
    int irxn = kspec - m_numComponents;
  
    doublereal s;
    int j, k;
    // Calculate the initial moles of the phase being born.
    //   Here we set it to 10x of the value which would cause the phase to be
    //   zeroed out within the algorithm.  We may later adjust the value.
    doublereal tPhaseMoles = 10. * m_totalMolNum * VCS_DELETE_PHASE_CUTOFF;
   
    
#ifdef DEBUG_MODE
    int existence = Vphase->exists();
    if (existence > 0) {
      printf("ERROR vcs_popPhaseRxnStepSizes called for a phase that exists!");
      exit(-1);
    }
    char anote[256];
    if (m_debug_print_lvl >= 2) {
      plogf("  ---  vcs_popPhaseRxnStepSizes() called to pop phase %s %d into existence\n",
	    Vphase->PhaseName.c_str(), iphasePop);   
    }
#endif
   // Section for a single-species phase
   //
   if (Vphase->m_singleSpecies) {
     s = 0.0;
     double *dnPhase_irxn = m_deltaMolNumPhase[irxn];
     for (j = 0; j < m_numComponents; ++j) {
       if (!m_SSPhase[j]) {
	 if (m_molNumSpecies_old[j] > 0.0) {
	   s += SQUARE(m_stoichCoeffRxnMatrix[irxn][j]) / m_molNumSpecies_old[j];
	 }
       }
     }
     for (j = 0; j < m_numPhases; j++) {
       Vphase = m_VolPhaseList[j];
       if (! Vphase->m_singleSpecies) {
	 if (m_tPhaseMoles_old[j] > 0.0) 
	   s -= SQUARE(dnPhase_irxn[j]) / m_tPhaseMoles_old[j];
       }
     }
     if (s != 0.0) {
       double s_old = s;
       s = vcs_Hessian_diag_adj(irxn, s_old);
#ifdef DEBUG_MODE
       if (s_old != s) {
	 sprintf(anote, "Normal calc: diag adjusted from %g "
		 "to %g due to act coeff",  s_old, s);
       }
#endif
       m_deltaMolNumSpecies[kspec] = -m_deltaGRxn_new[irxn] / s;
     } else {
       // Ok, s is equal to zero. We can not apply a sophisticated theory
       // to birth the phase. Just pick a small delta and go with it.
       m_deltaMolNumSpecies[kspec] = tPhaseMoles;
     }

     /*
      * section to do damping of the m_deltaMolNumSpecies[] 
      */
     for (j = 0; j < m_numComponents; ++j) {
       double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
       if (stoicC != 0.0) {
	 if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
	   double negChangeComp = - stoicC * m_deltaMolNumSpecies[kspec];
	   if (negChangeComp > m_molNumSpecies_old[j]) {
	     if (m_molNumSpecies_old[j] > 0.0) {
#ifdef DEBUG_MODE
	       sprintf(anote, "Delta damped from %g "
		       "to %g due to component %d (%10s) going neg", m_deltaMolNumSpecies[kspec],
		       -m_molNumSpecies_old[j]/stoicC, j,  m_speciesName[j].c_str());
#endif
	       m_deltaMolNumSpecies[kspec] = - 0.5 * m_molNumSpecies_old[j] / stoicC; 
	     } else {
#ifdef DEBUG_MODE
	       sprintf(anote, "Delta damped from %g "
		       "to %g due to component %d (%10s) zero", m_deltaMolNumSpecies[kspec],
		       -m_molNumSpecies_old[j]/stoicC, j,  m_speciesName[j].c_str());
#endif
	       m_deltaMolNumSpecies[kspec] = 0.0;
	     }
	   } 
	 }
       }
     }
     // Implement a damping term that limits m_deltaMolNumSpecies to the size of the mole number
     if (-m_deltaMolNumSpecies[kspec] > m_molNumSpecies_old[kspec]) {
#ifdef DEBUG_MODE
       sprintf(anote, "Delta damped from %g "
	       "to %g due to %s going negative", m_deltaMolNumSpecies[kspec],
	       -m_molNumSpecies_old[kspec],  m_speciesName[kspec].c_str());
#endif
       m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[kspec];
     }


   } else {
     vector<doublereal> fracDelta(Vphase->nSpecies());
     vector<doublereal> X_est(Vphase->nSpecies());
     fracDelta = Vphase->fractionCreationDeltas();
 
     double sumFrac = 0.0;
     for (k = 0; k < Vphase->nSpecies(); k++) {
       sumFrac += fracDelta[k];
     }
     for (k = 0; k < Vphase->nSpecies(); k++) {
       X_est[k] = fracDelta[k] / sumFrac;
     }

     doublereal deltaMolNumPhase = tPhaseMoles;
     doublereal damp = 1.0;
     m_deltaGRxn_tmp = m_molNumSpecies_old;
     double * molNumSpecies_tmp = DATA_PTR(m_deltaGRxn_tmp);


     for (k = 0; k < Vphase->nSpecies(); k++) {
       kspec = Vphase->spGlobalIndexVCS(k);
       double delmol =  deltaMolNumPhase * X_est[k];
       irxn = kspec - m_numComponents;
       if (kspec >= m_numComponents) {
	 for (j = 0; j < m_numComponents; ++j) {
	   double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
	   if (stoicC != 0.0) {
	     if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
	       molNumSpecies_tmp[j] +=  stoicC * delmol;
	     }
	   }
	 }
       }
     }

     doublereal ratioComp = 0.0;
     for (j = 0; j < m_numComponents; ++j) {
       double deltaJ = m_molNumSpecies_old[j] - molNumSpecies_tmp[j];
       if (molNumSpecies_tmp[j] < 0.0) {
	 ratioComp = 1.0;
	 if (deltaJ > 0.0) {
	   double delta0 = m_molNumSpecies_old[j];
	   double dampj = delta0 / deltaJ * 0.9;
	   if (dampj < damp) {
	     damp = dampj;
	   }
	 }
       } else {
	 if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
	   int jph =  m_phaseID[j];
	   if ((jph != iphasePop) && (!m_SSPhase[j])) {
	     double fdeltaJ = fabs(deltaJ);
	     if ( m_molNumSpecies_old[j] > 0.0) {
	       ratioComp = MAX(ratioComp, fdeltaJ/ m_molNumSpecies_old[j]);
	     }
	   }
	 }
       }
     }
     
     // We may have greatly underestimated the deltaMoles for the phase pop
     // Here we create a damp > 1 to account for this possibility.
     // We adjust upwards to make sure that a component in an existing multispecies
     // phase is modified by a factor of 1/1000.
     if (ratioComp > 1.0E-30) {
       if (ratioComp < 0.001) {
	 damp = 0.001 / ratioComp; 
       }
     }


     if (damp <= 1.0E-6) {
       return 3;
     }

     for (k = 0; k < Vphase->nSpecies(); k++) {
       kspec = Vphase->spGlobalIndexVCS(k);
       if (kspec < m_numComponents) {
	 m_speciesStatus[kspec] = VCS_SPECIES_COMPONENT;
       } else {
	 m_deltaMolNumSpecies[kspec] = deltaMolNumPhase * X_est[k] * damp;
	 if (X_est[k] > 1.0E-3) {
	   m_speciesStatus[kspec] = VCS_SPECIES_MAJOR;
	 } else {
	   m_speciesStatus[kspec] = VCS_SPECIES_MINOR;
	 }
       }
     }
    
   }

    return 0;
  }
  

  //

  double VCS_SOLVE::vcs_phaseStabilityTest(const int iph) {

    /*
     * We will use the _new state calc here
     */
    int kspec, irxn, k, i, kc, kc_spec;
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    doublereal deltaGRxn;

    // We will do a full newton calculation later, but for now, ...
    bool doSuccessiveSubstitution = true;
    double funcPhaseStability;
    vector<doublereal> X_est(Vphase->nSpecies(), 0.0);
    vector<doublereal> delFrac(Vphase->nSpecies(), 0.0);
    vector<doublereal> E_phi(Vphase->nSpecies(), 0.0);
    vector<doublereal> fracDelta_new(Vphase->nSpecies(), 0.0);
    vector<doublereal> fracDelta_old(Vphase->nSpecies(), 0.0);
    vector<doublereal> fracDelta_raw(Vphase->nSpecies(), 0.0);

    vector<doublereal> m_feSpecies_Deficient(m_numComponents, 0.0); 
    doublereal damp = 1.0;
    doublereal dampOld = 1.0;
    doublereal normUpdate = 1.0;
    doublereal normUpdateOld = 1.0;
    doublereal sum = 0.0;
    doublereal dirProd = 0.0;
    doublereal dirProdOld = 0.0;

    // get the activity coefficients
    Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(m_actCoeffSpecies_new));

    // Get the storred estimate for the composition of the phase if 
    // it gets created
    fracDelta_new = Vphase->fractionCreationDeltas();
  

    bool oneIsComponent = false;
    std::vector<int> componentList;

    for (k = 0; k < Vphase->nSpecies(); k++) {
      kspec = Vphase->spGlobalIndexVCS(k);
      if (kspec < m_numComponents) {
	oneIsComponent = true;
        componentList.push_back(k);
      }
    }

    for (k = 0; k < m_numComponents; k++) {
      m_feSpecies_Deficient[k]  = m_feSpecies_old[k];
    }
    normUpdate = 0.1 * vcs_l2norm(fracDelta_new);
    damp = 1.0E-2;

    if (doSuccessiveSubstitution) {

#ifdef DEBUG_MODE
      int KP = 0;
      if (m_debug_print_lvl >= 2) {
	plogf("   --- vcs_phaseStabilityTest() called\n");
	plogf("   ---  Its   X_old[%2d]  FracDel_old[%2d]  deltaF[%2d] FracDel_new[%2d]"
	      "  normUpdate     damp     FuncPhaseStability\n", KP, KP, KP, KP);
	plogf("   --------------------------------------------------------------"
	      "--------------------------------------------------------\n");
      } else {
	plogf("   --- vcs_phaseStabilityTest() called for phase %d\n", iph);
      }
#endif
      bool converged = false;
      for (int its = 0; its < 200  && (!converged); its++) {
	
	dampOld = damp;
	normUpdateOld = normUpdate;
	fracDelta_old = fracDelta_new;
	dirProdOld = dirProd;

	// Given a set of fracDelta's, we calculate the fracDelta's
	// for the component species, if any
	for (i = 0; i < (int) componentList.size(); i++) {
	  kc = componentList[i];
	  kc_spec = Vphase->spGlobalIndexVCS(kc);
	  fracDelta_old[kc] = 0.0;
	  for (k = 0; k <  Vphase->nSpecies(); k++) {
	    kspec = Vphase->spGlobalIndexVCS(k);
	    irxn = kspec - m_numComponents;
	    if (irxn >= 0) {
	      fracDelta_old[kc] += m_stoichCoeffRxnMatrix[irxn][kc_spec] *  fracDelta_old[k];
	    }
	  } 
	}

	// Now, calculate the predicted mole fractions, X_est[k]
	double sumFrac = 0.0;
	for (k = 0; k < Vphase->nSpecies(); k++) {
	  sumFrac += fracDelta_old[k];
	}
	double sum_Xcomp = 0.0;
	for (k = 0; k < Vphase->nSpecies(); k++) {
	  X_est[k] = fracDelta_old[k] / sumFrac;
	  kc_spec = Vphase->spGlobalIndexVCS(k);
	  if (kc_spec < m_numComponents) {
	    sum_Xcomp += X_est[k];
	  }
	}
	

	/*
	 * Feed the newly formed estimate of the mole fractions back into the
	 * ThermoPhase object
	 */
	Vphase->setMoleFractionsState(0.0, VCS_DATA_PTR(X_est), VCS_STATECALC_PHASESTABILITY); 
	
	/*
	 *   get the activity coefficients
	 */
	Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(m_actCoeffSpecies_new));

	/*
	 * first Calculate altered chemical potentials for component species
	 * belonging to this phase.
	 */
	for (i = 0; i < (int) componentList.size(); i++) {
	  kc = componentList[i];
	  kc_spec = Vphase->spGlobalIndexVCS(kc);
	  if ( X_est[kc] > VCS_DELETE_MINORSPECIES_CUTOFF) {
	    m_feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec] 
	      + log(m_actCoeffSpecies_new[kc_spec] * X_est[kc]);
	  } else {
	    m_feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec] 
	      + log(m_actCoeffSpecies_new[kc_spec] * VCS_DELETE_MINORSPECIES_CUTOFF);
	  }
	}

	for (i = 0; i < (int) componentList.size(); i++) {
	  kc = componentList[i];
	  kc_spec = Vphase->spGlobalIndexVCS(kc);
	  
	  for (k = 0; k <  Vphase->nSpecies(); k++) {
	    kspec = Vphase->spGlobalIndexVCS(k);
	    irxn = kspec - m_numComponents;
	    if (irxn >= 0) {
	      if (i == 0) {
		m_deltaGRxn_Deficient[irxn] = m_deltaGRxn_old[irxn];
	      }
	      double *dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
	      if (dtmp_ptr[kc_spec] != 0.0) {
		m_deltaGRxn_Deficient[irxn] += 
		  dtmp_ptr[kc_spec] * (m_feSpecies_Deficient[kc_spec]- m_feSpecies_old[kc_spec]);
	      }
	    }

	  }
	} 
	
	/*
	 *  Calculate the E_phi's
	 */
	sum = 0.0;
	funcPhaseStability = sum_Xcomp - 1.0;
	for (k = 0; k <  Vphase->nSpecies(); k++) {
	  kspec = Vphase->spGlobalIndexVCS(k);
	  irxn = kspec - m_numComponents;
	  if (irxn >= 0) {
	    deltaGRxn = m_deltaGRxn_Deficient[irxn];
	    if (deltaGRxn >  50.0) deltaGRxn =  50.0;
	    if (deltaGRxn < -50.0) deltaGRxn = -50.0;
	    E_phi[k] = std::exp(-deltaGRxn) / m_actCoeffSpecies_new[kspec];
	    sum +=  E_phi[k];
	    funcPhaseStability += E_phi[k];
	  } else {
	    E_phi[k] = 0.0;
	  }
	}

	/*
	 * Calculate the raw estimate of the new fracs
	 */
	for (k = 0; k <  Vphase->nSpecies(); k++) {
	  kspec = Vphase->spGlobalIndexVCS(k);
	  irxn = kspec - m_numComponents;
	  double b =  E_phi[k] / sum * (1.0 - sum_Xcomp);
	  if (irxn >= 0) {
	    fracDelta_raw[k] = b;
	  } 
	}


	// Given a set of fracDelta's, we calculate the fracDelta's
	// for the component species, if any
	for (i = 0; i < (int) componentList.size(); i++) {
	  kc = componentList[i];
	  kc_spec = Vphase->spGlobalIndexVCS(kc);
	  fracDelta_raw[kc] = 0.0;
	  for (k = 0; k <  Vphase->nSpecies(); k++) {
	    kspec = Vphase->spGlobalIndexVCS(k);
	    irxn = kspec - m_numComponents;
	    if (irxn >= 0) {
	      fracDelta_raw[kc] += m_stoichCoeffRxnMatrix[irxn][kc_spec] * fracDelta_raw[k];
	    }
	  } 
	}



	/*
	 * Now possibly dampen the estimate.
	 */
	doublereal sumADel = 0.0;
	for (k = 0; k <  Vphase->nSpecies(); k++) {
	  delFrac[k] = fracDelta_raw[k] - fracDelta_old[k];
	  sumADel += fabs(delFrac[k]);
	}
	normUpdate = vcs_l2norm(delFrac);

	dirProd = 0.0;
     	for (k = 0; k <  Vphase->nSpecies(); k++) {
	  dirProd += fracDelta_old[k] * delFrac[k];
	}
	bool crossedSign = false;
	if (dirProd * dirProdOld < 0.0) {
	  crossedSign = true;
	}


	damp = 0.5;
	if (dampOld < 0.25) {
	  damp = 2.0 * dampOld;
	}
	if (crossedSign) {
	  if (normUpdate *1.5 > normUpdateOld) {
	    damp = 0.5 * dampOld;
	  } else if (normUpdate *2.0 > normUpdateOld) {
	    damp = 0.8 * dampOld;
	  }
	} else {
	  if (normUpdate > normUpdateOld * 2.0) {
	    damp = 0.6 * dampOld;
	  } else if (normUpdate > normUpdateOld * 1.2) {
	    damp = 0.9 * dampOld;
	  }
	}

	for (k = 0; k < Vphase->nSpecies(); k++) {
	  if (fabs(damp * delFrac[k]) > 0.3*fabs(fracDelta_old[k])) {
	    damp = MAX(0.3*fabs(fracDelta_old[k]) / fabs( delFrac[k]), 
		       1.0E-8/fabs( delFrac[k]));
	  }
	  if (delFrac[k] < 0.0) {
	    if (2.0 * damp * (-delFrac[k]) > fracDelta_old[k]) {
	      damp = fracDelta_old[k] / (2.0 * (-delFrac[k]));
	    }
	  }
	  if (delFrac[k] > 0.0) {
	    if (2.0 * damp * delFrac[k] > fracDelta_old[k]) {
	      damp = fracDelta_old[k] / (2.0 * delFrac[k]);
	    }
	  }
	}
	if (damp < 0.000001) {
	  damp = 0.000001;
	}

	for (k = 0; k <  Vphase->nSpecies(); k++) {
	  fracDelta_new[k] = fracDelta_old[k] + damp * (delFrac[k]);
	}

#ifdef DEBUG_MODE
	if (m_debug_print_lvl >= 2) {
	  plogf("  --- %3d %12g %12g %12g %12g %12g %12g %12g\n", its, X_est[KP], fracDelta_old[KP],
		delFrac[KP], fracDelta_new[KP], normUpdate, damp, funcPhaseStability);
	}
#endif
	
	if (normUpdate < 1.0E-5) {
	  converged = true;
	}

      }

      if (converged) {
	Vphase->setMoleFractionsState(0.0, VCS_DATA_PTR(X_est), 
				      VCS_STATECALC_PHASESTABILITY);
	Vphase->setFractionCreationDeltas( VCS_DATA_PTR(fracDelta_new));
      }


    } else {
      printf("not done yet\n");
      exit(-1);
    }
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
      plogf("  ------------------------------------------------------------"
	    "-------------------------------------------------------------\n");
    } else if (m_debug_print_lvl == 1) {
      if (funcPhaseStability > 0.0) {
	plogf("  --- phase %d with func = %g is to be born\n", iph, funcPhaseStability);
      } else {
	plogf("  --- phase %d with func = %g stays dead\n", iph, funcPhaseStability);
      }
    }
#endif
    return funcPhaseStability;
  }

}

