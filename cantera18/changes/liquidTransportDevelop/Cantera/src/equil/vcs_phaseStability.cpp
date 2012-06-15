/**
 * @file vcs_phaseStability.cpp
 *  Implementation class for functions associated with determining the stability of a phase 
 *   (see Class \link Cantera::VCS_SOLVE VCS_SOLVE\endlink and \ref equilfunctions ).
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
#include "vcs_species_thermo.h"
#include "vcs_VolPhase.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <cstring>

using namespace std;

namespace VCSnonideal {

  //====================================================================================================================
  // Utility function that evaluates whether a phase can be popped into existence
  /*
   * A phase can be popped iff the stoichiometric coefficients for the
   * component species, whose concentrations will be lowered during the
   * process, are positive by at least a small degree.
   *  
   * If one of the phase species is a zeroed component, then the phase can
   * be popped if the component increases in mole number as the phase moles
   * are increased.
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
#ifdef DEBUG_MODE
      if (m_molNumSpecies_old[kspec] > 0.0) {
	printf("ERROR vcs_popPhasePossible we shouldn't be here %d %g > 0.0",
	       kspec, m_molNumSpecies_old[kspec]);
	exit(-1);
      }
#endif
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
		if (m_molNumSpecies_old[j] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
		  iPopPossible = false;
		}
	      }
	    }
	  }
	}
	if (iPopPossible == true) {
	  return true;
	}
      } else {
	/*
	 * We are here when the species in the phase is a component. Its mole number is zero.
	 * We loop through the regular reaction looking for a reaction that can pop the
	 * component.
	 */
        //printf("WE are here at new logic - CHECK\n");
	for (int jrxn = 0; jrxn < m_numRxnRdc; jrxn++) {
	  bool foundJrxn = false;
	  // First, if the component is a product of the reaction
	  if (m_stoichCoeffRxnMatrix[jrxn][kspec] > 0.0) {
	    foundJrxn = true;
	    for (int kcomp = 0; kcomp < m_numComponents; kcomp++) {
	      if (m_stoichCoeffRxnMatrix[jrxn][kcomp] < 0.0) {
		if (m_molNumSpecies_old[kcomp] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
		  foundJrxn = false;
		}
	      }
	    }
            if (foundJrxn) {
	      //printf("We have found a component phase pop! CHECK1 \n");
	      return true;
	    }
	  }
	  // Second we are here if the component is a reactant in the reaction, and the reaction goes backwards.
	  else if (m_stoichCoeffRxnMatrix[jrxn][kspec] < 0.0) {
	    foundJrxn = true;
	    int jspec = jrxn + m_numComponents;
	    if (m_molNumSpecies_old[jspec] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
	      foundJrxn = false;
	      continue;
	    }
	    for (int kcomp = 0; kcomp < m_numComponents; kcomp++) {
	      if (m_stoichCoeffRxnMatrix[jrxn][kcomp] > 0.0) {
		if (m_molNumSpecies_old[kcomp] <= VCS_DELETE_ELEMENTABS_CUTOFF*0.5) {
		  foundJrxn = false;
		}
	      }
	    }
            if (foundJrxn) {
	      //printf("We have found a component phase pop! CHECK2 \n");
	      return true;
	    }
	  }
	}
      }
    }
    return false;
  }
  //====================================================================================================================
 
  int inList(const std::vector<int> &list, int val)
  {
    for (int i = 0; i < (int) list.size(); i++) {
      if (val == list[i]) {
	return i;
      }
    }
    return -1;
  }
  
  //====================================================================================================================
  // Determine the list of problems that need to be checked to see if there are any phases pops
  /*
   *  This routine evaluates and fills in the following quantities
   *              phasePopProblemLists_
   *
   *  Need to work in species that are zeroed by element constraints
   *
   *  @return    Returns the number of problems that must be checked.
   */
  int  VCS_SOLVE::vcs_phasePopDeterminePossibleList() {

    int nfound = 0;
    int irxn, kspec;
    vcs_VolPhase *Vphase = 0;
    int iph, j, k;
    int nsp;
    double stoicC;
    double molComp;
    std::vector<int> linkedPhases;
    phasePopProblemLists_.clear();

    /*
     *  This is a vector over each component.
     *  For zeroed components it lists the phases, which are currently zeroed,
     *     which have a species with a positive stoichiometric value wrt the component.
     *     Therefore, we could pop the component species and pop that phase at the same time
     *     if we considered no other factors than keeping the component mole number positve.
     *
     *     It does not count species with positive stoichiometric values if that species
     *     already has a positive mole number. The phase is already popped. 
     */
    std::vector< std::vector<int> > zeroedComponentLinkedPhasePops(m_numComponents);
    /*
     *  The logic below calculates zeroedComponentLinkedPhasePops
     */
    for (j = 0; j < m_numComponents; j++) {
      if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
	molComp = m_molNumSpecies_old[j];
	if (molComp <= 0.0) {	
	  std::vector<int> &jList = zeroedComponentLinkedPhasePops[j];
	  iph = m_phaseID[j];
	  jList.push_back(iph);
	  for (irxn = 0; irxn < m_numRxnTot; irxn++) {
	    kspec = irxn +  m_numComponents;
	    iph = m_phaseID[kspec];
	    Vphase = m_VolPhaseList[iph];
	    int existence = Vphase->exists();
	    if (existence < 0) {
	      stoicC = m_stoichCoeffRxnMatrix[irxn][j];
	      if (stoicC > 0.0) {
		if (inList(jList, iph) != -1) {
		  jList.push_back(iph);
		}
	      }
	    }
	  }
	}
      }
    }
    /*
     *   This is a vector over each zeroed phase
     *   For zeroed phases, it lists the components, which are currently zereoed,
     *     which have a species with a negative stoichiometric value wrt one or more species in the phase.
     *     Cut out components which have a pos stoichiometric value with another species in the phase.
     */
    std::vector< std::vector<int> > zeroedPhaseLinkedZeroComponents(m_numPhases);
    /*
     *   The logic below calculates  zeroedPhaseLinkedZeroComponents
     */
    for (iph = 0; iph < m_numPhases; iph++) {
      std::vector<int> &iphList = zeroedPhaseLinkedZeroComponents[iph];
      iphList.clear();
      Vphase = m_VolPhaseList[iph];
      int existence = Vphase->exists();
      if (existence < 0) {
       
	linkedPhases.clear();
	nsp = Vphase->nSpecies();
	for (k = 0; k < nsp; k++) {

	  kspec = Vphase->spGlobalIndexVCS(k);
	  irxn = kspec - m_numComponents;
	  
	  for (j = 0; j < m_numComponents; j++) {
	    if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
	      molComp = m_molNumSpecies_old[j];
	      if (molComp <= 0.0) {
		stoicC = m_stoichCoeffRxnMatrix[irxn][j];
		if (stoicC < 0.0) {
		  bool foundPos = false;
		  for (int kk = 0; kk < nsp; kk++) {
		    int kkspec  = Vphase->spGlobalIndexVCS(kk);
		    int iirxn = kkspec - m_numComponents;
		    if (iirxn >= 0) {
		      if (m_stoichCoeffRxnMatrix[iirxn][j] > 0.0) {
			foundPos = true;
		      }
		    }
		  }
		  if (!foundPos) {
		    if (inList(iphList, j) != -1) {
		      iphList.push_back(j);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  
    /*
     *  Now fill in the   phasePopProblemLists_  list.
     * 
     */
    for (iph = 0; iph < m_numPhases; iph++) {
      Vphase = m_VolPhaseList[iph];
      int existence = Vphase->exists();
      if (existence < 0) {
	std::vector<int> &iphList = zeroedPhaseLinkedZeroComponents[iph];
	std::vector<int> popProblem(0);
	popProblem.push_back(iph);
	for (int i = 0; i < (int) iphList.size(); i++) {
	  j = iphList[i];
	  std::vector<int> &jList = zeroedComponentLinkedPhasePops[j];
	  for (int jjl = 0; jjl < (int) jList.size(); jjl++) {
	    int jph = jList[jjl];
	    if (inList(popProblem, jph) != -1) {
	      popProblem.push_back(jph);
	    }
	  }
	}
        phasePopProblemLists_.push_back(popProblem);
      }
    }

    return nfound;
  }


  //====================================================================================================================
  // Decision as to whether a phase pops back into existence
  /*
   * @return returns the phase id of the phases that pops back into 
   *         existence. Returns -1 if there are no phases
   */
  int VCS_SOLVE::vcs_popPhaseID(std::vector<int> & phasePopPhaseIDs) {
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
      plogf("   ---   Phase                 Status       F_e        MoleNum\n");
      plogf("   --------------------------------------------------------------------------\n");
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
	  plogf("  ---    %18s %5d           NA       %11.3e\n", 
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
	  if (Fephase < 0.0) {
	    strcpy(anote," (not stable)");
	    if (m_tPhaseMoles_old[iph] > 0.0) {
	      printf("shouldn't be here\n");
	      exit(-1);
	    }
	  }
#endif

#ifdef DEBUG_MODE
	  if (m_debug_print_lvl >= 2) {
	    plogf("  ---    %18s %5d %10.3g %10.3g %s\n", 
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
    phasePopPhaseIDs.resize(0);
    if (iphasePop >= 0) {
      phasePopPhaseIDs.push_back(iphasePop);
    }

    /*
     *   Insert logic here to figure out if phase pops are linked together. Only do one linked
     *   pop at a time.
     */

#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
      plogf("   ---------------------------------------------------------------------\n");
    }
#endif
    return iphasePop;
  }
  //====================================================================================================================
  // Calculates the deltas of the reactions due to phases popping
  // into existence
  /*
   * @param iphasePop  Phase id of the phase that will come into existence
   *
   * Output 
   * ------- 
   * m_deltaMolNumSpecies[kspec] : Adjustments to the species vector
   *                               due to a phase or group of phases popping
   *                               back into existence.
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
    std::vector<int> creationGlobalRxnNumbers;
  
    doublereal s;
    int j, k;
    //  Calculate the initial moles of the phase being born.
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
     fracDelta = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);
 
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
 
  //====================================================================================================================
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
    vector<int>        creationGlobalRxnNumbers(Vphase->nSpecies(), -1);
    vcs_dcopy(VCS_DATA_PTR(m_deltaGRxn_Deficient), VCS_DATA_PTR(m_deltaGRxn_old), m_numRxnRdc);

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
    fracDelta_new = Vphase->creationMoleNumbers(creationGlobalRxnNumbers);
  

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
      } else if (m_debug_print_lvl == 1) {
	plogf("   --- vcs_phaseStabilityTest() called for phase %d\n", iph);
      }
#endif

      for (k = 0; k < Vphase->nSpecies(); k++) {
	if (fracDelta_new[k] < 1.0E-13) {
	  fracDelta_new[k] = 1.0E-13;
	}
      }
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
        // Necessary because this can be identically zero. -> we need to fix this algorithm!
	if (sumFrac <= 0.0) {
	  sumFrac = 1.0;
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
	 * First calculate altered chemical potentials for component species
	 * belonging to the phase considered for deletion.
	 */
	for (i = 0; i < (int) componentList.size(); i++) {
	  kc = componentList[i];
	  kc_spec = Vphase->spGlobalIndexVCS(kc);
	  if (X_est[kc] > VCS_DELETE_MINORSPECIES_CUTOFF) {
	    m_feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec] + log(m_actCoeffSpecies_new[kc_spec] * X_est[kc]);
	  } else {
	    m_feSpecies_Deficient[kc_spec] = m_feSpecies_old[kc_spec] + log(m_actCoeffSpecies_new[kc_spec] * VCS_DELETE_MINORSPECIES_CUTOFF);
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
		m_deltaGRxn_Deficient[irxn] += dtmp_ptr[kc_spec] * (m_feSpecies_Deficient[kc_spec]- m_feSpecies_old[kc_spec]);
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
	Vphase->setCreationMoleNumbers(VCS_DATA_PTR(fracDelta_new), creationGlobalRxnNumbers);
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
  //====================================================================================================================
  //   Check to see if an existing phase passes the phase non-stability test
  /*
   *   If the number returned is less than zero, then the phase can probably be eliminated.
   *
   *  Complicated example: test how this routine does when trying to get rid of phase 1.
   *
   *    components         PhaseID
   *          CH            1
   *          O             1
   *
   *     non-components
   *          CHO2          1
   *          CHO(S)        2
   * 
   *  @param iph   PhaseID (will be replaced later by a phasePop structure)
   *
   *  @return    returns the stability function value.
   *             If phase pop is not possible this function returns 0.0 exactly.
   */
  double VCS_SOLVE::vcs_phaseStabilityFE_value(int iph) {
    int k, kspec, kc_spec, irxn;
    double funcPhaseStability = 0.0;

    
#ifdef DEBUG_MODE_NOT
   
    if (m_VCount->Its > 4685) {
      printf("we are here\n");
    }
#endif
   
    /*
     *  right now we assume one phase. This will be relaxed shortly.
     */
    // iph = popProb->phaseList_[0];
    /*
     *  Get the volume phase object
     */
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    int nsp = Vphase->nSpecies();
    /*
     *  Get the integer representing whether the phase exists or not. If it exists existence is greater than zero.
     *  This function will only be used for phases which exist 
     */
    int existence = Vphase->exists();
    if (existence <= 0) {
      /*
       *  This is the current way to test for phase pops for phases which don't exist. We'll put them in here
       *  even though this approach is currently not used here.
       */
      if (vcs_popPhasePossible(iph)) {
	funcPhaseStability = vcs_phaseStabilityTest(iph);
	return funcPhaseStability;
      } else {
	return 0.0;
      }
    }

    //vector<doublereal> m_feSpecies_Deficient(m_numComponents, 0.0); 
    /*
     * Pull in the mole fractions for the phase to be tested
     */
    const std::vector<double> &X_est =  Vphase->moleFractions();
    /*
     *  Store an initial value of deltaG deficient from the storred value of deltaGRxn_old , which is current
     */
    vcs_dcopy(VCS_DATA_PTR(m_deltaGRxn_Deficient), VCS_DATA_PTR(m_deltaGRxn_old), m_numRxnRdc);

    /*
     *  Get a list of species in the current phase which are components
     */
    double sum_Xcomp = 0.0;
    bool oneIsComponent = false;
    std::vector<int> componentNoRxnList(0);
    std::vector<int> componentRxnList(0);
    // List of reaction ids to use for each component.
    std::vector<int> componentIRxnList(0);
    std::vector<double> componentdeltaGRxnList(0);
    std::vector<double> componentStoichC(0);
    
    /*
     *  Try to find reactions for each component that have non-zero moles
     */
    for (k = 0; k < Vphase->nSpecies(); k++) {
      kspec = Vphase->spGlobalIndexVCS(k);
      if (kspec < m_numComponents) {
	oneIsComponent = true;
	/*
	 *  Find a reaction that eliminates the component without causing trouble
	 */
	int iirxnRxnList = -1;
	int iirxnRxnAlt = -1;
	double gmax = -1.0E300;
	double stoichCRxnList = 0.0;;
	if (X_est[k] > 0.0) {
	  for (irxn = 0; irxn < m_numRxnRdc; irxn++) {
	    kc_spec = irxn + m_numComponents;
	    double *dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
	    // Choose a reaction that has a non-zero stoich coefficient for the component
	    if (dtmp_ptr[kspec] != 0.0) {
	      //  Store the id of the candidate reaction
	      int iirxn = irxn;
	      //  If the reaction creates a non-stoich species in the deletion phase, eliminate it from consideration
	      //  because it doesn't represent a different degree of freedom.
	      if (m_phaseID[kc_spec] == iph) {
		iirxn = -1;
	      }
	      //  If the reaction involves a component that has a zero mole number, eliminate it from consideration
	      //  as the component is from a phase that hasn't popped yet, or can't pop.
	      for (int kk = 0; kk < m_numComponents; kk++) {
		if (kk != kspec) {
		  if (dtmp_ptr[kk] != 0.0) {
		    if (m_molNumSpecies_old[kk] <= 0.0) {
		      iirxn = -1;
		    }
		  }
		}
	      }
	      //  If the reaction creates a species in a phase that doesn't exist, eliminate it from consideration.
	      int jph = m_phaseID[kc_spec];
	      vcs_VolPhase *VJphase = m_VolPhaseList[jph];
	      int exists = VJphase->exists();
	      if (exists <= 0) {
		iirxn = -1;
	      }
	      //  If the reaction has alreay been used before for another component in the phaseStability problem, eliminate it from consideration.
	      //  Note, we didn't do any optimization here about what reaction to use for which component.
	      //  So, if we can't find a full range of component reactions, we are in a more complicated situation.
	      //  The issue is that the composition of the exiting phase has to belong to a diminished range space supplied by the
	      //  exiting component reactions. An example of this situation is supplied in the notes.
	      //  We can only really check on this condition in the accompanying routine which goes through the math to actually eliminate the phase.
	      //  
	      for (int ii = 0; ii < (int) componentIRxnList.size(); ii++) {
		if (irxn == componentIRxnList[ii]) {
		  iirxn = -1; 
		  iirxnRxnAlt = irxn;
		}
	      }
	      if (iirxn < 0) {
		continue;
	      }
	      /*
	       *  we're good if we last till here.
	       */
	      double stoichC = dtmp_ptr[kspec];
	      double gtmp =  m_deltaGRxn_old[irxn] / stoichC;
	      if (gtmp > gmax) {
		iirxnRxnList = iirxn;
		gmax = gtmp;
		stoichCRxnList = stoichC;
	      }
	    }
	  }
	}
	/*
	 *  We now separate the components into two groups
	 *      1) First group contains components that we found a reaction for
	 *      2) Second group contains components that we didn't find a reaction for, or for some components involving
	 *         reduced range space exit cases.
	 */
	if (iirxnRxnList >= 0) {
	  componentRxnList.push_back(k);
	  componentIRxnList.push_back(iirxnRxnList);
	  componentdeltaGRxnList.push_back(gmax);
	  componentStoichC.push_back(stoichCRxnList);
	} else {
	  componentNoRxnList.push_back(k);
	  // We only allow non-zero mole fractions for the particular case of reduced range spaces of components reactions.
	  sum_Xcomp += X_est[k];
	  // If the component is nonzero and not attached to any valid exit reaction, not even one involving a reduced range space,
	  // then it is a terminal condition for popping the phase. Return with a zero.
	  if (X_est[k] > 0.0  && iirxnRxnAlt < 0) {
	    return 0.0;
	  }
	}
      }
    }
    // If all species in the phase are nonreacting components, then we can't eliminate the phase.
    if ((int) componentNoRxnList.size() == nsp) {
      return 0.0;
    }

    /*
     * Get the activity coefficients
     */
    //  Vphase->sendToVCS_ActCoeff(VCS_STATECALC_OLD, VCS_DATA_PTR(m_actCoeffSpecies_old));
    /*
     * First calculate altered chemical potentials for component species belonging to this phase.
     */
    // for (i = 0; i < (int) componentNoRxnList.size(); i++) {
    //    kc = componentNoRxnList[i];
    //    kc_spec = Vphase->spGlobalIndexVCS(kc);
    //    if (X_est[kc] > VCS_DELETE_MINORSPECIES_CUTOFF) {
    //      m_feSpecies_Deficient[kc_spec] = m_SSfeSpecies[kc_spec] + log(m_actCoeffSpecies_old[kc_spec] * X_est[kc]);
    //    } else {
    // 	    m_feSpecies_Deficient[kc_spec] = m_SSfeSpecies[kc_spec] + log(m_actCoeffSpecies_old[kc_spec] * VCS_DELETE_MINORSPECIES_CUTOFF);
    //    }
    //  }

    /*
     *  Calculate the E_phi's
     */
    double E_phi_k;
    funcPhaseStability = sum_Xcomp - 1.0;
    for (k = 0; k <  Vphase->nSpecies(); k++) {
      kspec = Vphase->spGlobalIndexVCS(k);
      irxn = kspec - m_numComponents;
      if (irxn >= 0) {
	double deltaGRxn_deficient = m_deltaGRxn_old[irxn] - log(X_est[k] * m_actCoeffSpecies_old[kspec]);
	if (deltaGRxn_deficient >  50.0) deltaGRxn_deficient =  50.0;
	if (deltaGRxn_deficient < -50.0) deltaGRxn_deficient = -50.0;
	E_phi_k = std::exp(-deltaGRxn_deficient) / m_actCoeffSpecies_old[kspec];
	funcPhaseStability += E_phi_k;
      }
    }
    // Now do the components that have exit channels.
    for (int kk = 0; kk < (int) componentRxnList.size(); kk++) {
      k = componentRxnList[kk];
      kspec = Vphase->spGlobalIndexVCS(k);
    
      double deltaGcomp = componentdeltaGRxnList[kk];
      double deltaGcomp_deficient =  deltaGcomp - log(X_est[k] * m_actCoeffSpecies_old[kspec]);
      if (deltaGcomp_deficient >  50.0) deltaGcomp_deficient =  50.0;
      if (deltaGcomp_deficient < -50.0) deltaGcomp_deficient = -50.0;
      E_phi_k = std::exp(-deltaGcomp_deficient) / m_actCoeffSpecies_old[kspec];
   
      funcPhaseStability += E_phi_k;
    }

    return funcPhaseStability;
  }
 
  //=========================================================================================================================
  // Calculate a set of reaction delta and species deltas that eliminates a phase or group of phases
  /*
   *  This has repeatedly been shown to be an important and somewhat hard undertaking. 
   *  This is needed to thwart large number of iterations at stability boundaries for phases.
   *
   *  @return  1 successful
   *          -1 unsuccessful can't be done under the current conditions without triggering an unwanted condition
   *           0 unsuccessful because we can't find a case where the total gibbs free energy decreases. We did however
   *             find reaction coordinates to get rid of the phase
   */
  int  VCS_SOLVE::vcs_popPhase_calcDeleteDelta(const int iph)
  {
    int k, kspec, irxn, j, jph, exists, ifound;
    double dx;
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    int success = -1;
    double total_new, deltaG;
    std::vector<double>	RxnDeltaBest;
    std::vector<int>	componentListBest;
    std::vector<int>	componentRxnListBest ;
#ifdef DEBUG_MODE
    bool ttt;
    if (m_VCount->Its >= 326) {
      // printf("we are here\n");
    }
#endif
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
      plogf("  ---  vcs_popPhase_calcDeleteDelta called to find a finite method to delete phase %d (%-21.21s)\n",
	    iph, Vphase->PhaseName.c_str());
    }
#endif
     
    /*
     *  Get a list of species in the current phase which are components
     */

    std::vector<int> componentNoRxnList;
    std::vector<int> componentRxnList;
    std::vector<int> componentIRxnList;
    std::vector<double> componentdeltaGRxnList;
    std::vector<double> componentStoichC;
    std::vector<int> componentMatrix;
    std::vector<int> componentRxnMatrix;  
    std::vector<int> componentListTrial; 
    std::vector<int>  componentListGlobal;
    std::vector< std::vector<int> > componentListGlobal_PossibleRxns;
    double *sc_irxn;
    std::vector<double> RxnDelta(m_numSpeciesTot, 0.0);
    double total_old = vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_feSpecies_old), 
				       VCS_DATA_PTR(m_tPhaseMoles_old));
    /*
     * This is the final extent of reaction delta
     */
    std::vector<double> extentRxnVector(m_numSpeciesTot, 0.0);
    

    vcs_dcopy(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_molNumSpecies_old), m_numSpeciesTot);
    vcs_dzero(VCS_DATA_PTR(m_deltaMolNumSpecies), m_numSpeciesTot);
    vcs_dzero(VCS_DATA_PTR(RxnDelta), m_numSpeciesTot);

    /*
     *  Calculate the mole numbers of all species after we have eliminated all of the noncomponents of phases
     *  that we want to eliminate
     */
    for (k = 0; k < Vphase->nSpecies(); k++) {
      kspec = Vphase->spGlobalIndexVCS(k);
      if (kspec >= m_numComponents) {
	irxn = kspec - m_numComponents;
	sc_irxn = m_stoichCoeffRxnMatrix[irxn];
	dx = m_molNumSpecies_old[kspec];
	m_molNumSpecies_new[kspec] = 0.0;
	m_deltaMolNumSpecies[kspec] = -dx;
	RxnDelta[irxn] = -dx;
	if (dx != 0.0) {
	  for (j = 0; j < m_numComponents; ++j) {
	    m_molNumSpecies_new[j] -= sc_irxn[j] * dx;
	    m_deltaMolNumSpecies[j] -= sc_irxn[j] * dx;      
	  }
	}
      }
    }

    /*
     * At the end of this step, we have found out what components need to be zeroed.
     *     We store this in the list componentListGlobal. If we can't zero the mole numbers
     *     in this list, we can't delete the phase.
     */
    for (k = 0; k < Vphase->nSpecies(); k++) {
      kspec = Vphase->spGlobalIndexVCS(k);
      /*
       *  If a species is zeroed, we don't need to find a way to eliminate its mole number
       */
      if (m_molNumSpecies_old[kspec] > 0.0) {
	if (kspec < m_numComponents) {
	  /*
	   * Will add in a check for uniqueness when this gets to be a multiphase process.
	   */
	  componentListGlobal.push_back(kspec);
	}
      }
    }

    /*
     *  Make a list of the components which are good, i.e., are in phases which currently exist.
     *  This is called componentListGood.
     */
    std::vector<int> componentListGood(m_numComponents, 1);
    for (int jj = 0; jj < m_numComponents; jj++) {
      jph = m_phaseID[jj];
      int exist = Vphase->exists();
      if (exist <= 0) {
	componentListGood[jj] = 0;
      }
    }

    /*
     *  Make a matrix of the good reactions for all components that need to be zeroed.
     *   This is called  componentListGlobal_PossibleRxns.
     *   The solution to the problem will consist of picking reactions for each component in componentListGlobal
     *   
     */
    for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
      int kcomp = componentListGlobal[jj];
      int rgood = 1;
      std::vector<int> goodRxns(0);
      for (k = m_numComponents; k < m_numSpeciesTot; k++) {
	irxn = k - m_numComponents;
	sc_irxn = m_stoichCoeffRxnMatrix[irxn];
	rgood = 0;
	if (sc_irxn[kcomp] != 0.0) {
	  rgood = 1;
	  jph = m_phaseID[k];
	  vcs_VolPhase *VJphase = m_VolPhaseList[jph];
	  exists = VJphase->exists();
	  if (exists <= 0) {
	    rgood = 0;
	  }
	  for (int j = 0; j < m_numComponents; j++) {
	    if (sc_irxn[j] != 0.0) {
	      if (!componentListGood[j])
		rgood = 0;
	    }
	  }
	}
	if (rgood) {
	  goodRxns.push_back(irxn);
	}
      }
      /*
       *  If we can't find any reaction for a particular nonzero component in the deletion
       *  phase, we have failed. Return!
       */
      if (goodRxns.size() == 0) {
	return success;
      }
      componentListGlobal_PossibleRxns.push_back(goodRxns);
    }

    /*
     *   componentRxnListTrial[] is the current trial 
     */
    std::vector<int> componentRxnListTrial(componentListGlobal.size(), 0);
    std::vector<int> componentListLocal(componentListGlobal);

    /*
     *  Find the total number of combinations to check
     */
    int nCombos = 1;
    for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
      const std::vector<int> & grs = componentListGlobal_PossibleRxns[jj];
      nCombos *= grs.size();
    }
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 3) {
      plogf("  ----    List of possible exit reactions for component species in the current phasePop phases:\n");
      for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	int icomp = componentListGlobal[jj];
	const std::vector<int> & grs = componentListGlobal_PossibleRxns[jj];
	plogf("  ---            %-21.21s:", m_speciesName[icomp].c_str());
	for (int jjj = 0; jjj <  (int) grs.size(); jjj++) {
	  int kk = grs[jjj] + m_numComponents;
	  plogf(" %3d(%-15.15s),  ", grs[jjj], m_speciesName[kk].c_str());
	}
	plogf("\n");
      }
    }
#endif
    /*
     *  Determine the maximum size of the range space and store it in rsc
     *  This is the maximum number of different reactions in the reaction vector.
     *  We will only consider combinations that have this property, as the rest will
     *  be degenerate.
     */
    int rsc = 0;
    for (int ii = 1; ii <= nCombos; ii++) {
      int iNum = ii;
      for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	const std::vector<int> & grs = componentListGlobal_PossibleRxns[jj];
	int sz = grs.size();
	int jNum = (iNum -1) % sz;
	componentRxnListTrial[jj] = grs[jNum];
	iNum = 1 + (iNum - 1)/sz;
      }
      int nFound = 0;
      ifound = -1;
      for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	if (componentRxnListTrial[jj] >= 0) {
	  ifound = componentRxnListTrial[jj];
	  componentRxnListTrial[jj] = -1;
	  nFound++;
	}
	for (int kk = jj+1; kk < (int) componentListGlobal.size(); kk++) {
	  if (componentRxnListTrial[kk] == ifound) {
	    componentRxnListTrial[kk] = -1;
	  }
	}
      }
      rsc = MAX(rsc, nFound);
      if (rsc == (int) componentListGlobal.size()) {
	break;
      }
    }


    int bestCombo = -1;
    double lowestDelta = 1.0E300;
    int cs = (int) componentListGlobal.size();
    double sm[300];
    double rhs[15];

    for (int iCombo = 1; iCombo <= nCombos; iCombo++) {
      int iNum = iCombo;
      for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	const std::vector<int> & grs = componentListGlobal_PossibleRxns[jj];
	int sz = grs.size();
	int jNum = (iNum -1) % sz;
	componentRxnListTrial[jj] = grs[jNum];
	iNum = 1 + (iNum - 1)/sz;
      }

#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 3) {
      plogf("  ----    Trial # %d:\n", iCombo);
      for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	int icomp = componentListGlobal[jj];
	int irxn = componentRxnListTrial[jj];
	int kk = irxn + m_numComponents;
        printf("  ---                      Component            Reaction\n");
	printf("  ---            %3d(%-21.21s)      %3d(%-15.15s)  \n",
	       icomp, m_speciesName[icomp].c_str(), irxn,  m_speciesName[kk].c_str());
      }
    }
#endif

      /*
       * Check to see that the current combination has a maximum range space
       */
      std::vector<int> ccc(componentRxnListTrial);
      int nFound = 0;
      for (int jj = 0; jj < (int) ccc.size(); jj++) {
	if (ccc[jj] >= 0) {
	  ifound = ccc[jj];
	  ccc[jj] = -1;
	  nFound++;
	}
	for (int kk = jj+1; kk < (int) ccc.size(); kk++) {
	  if (ccc[kk] == ifound) {
	    ccc[kk] = -1;
	  }
	}
      }
      // if less than max, then go to the next trial
      if (nFound < rsc) {
#ifdef DEBUG_MODE
	    if (m_debug_print_lvl >= 3) {
	      plogf("          combo %d disallowed because couldn't find %d distinct reactions. Found %d instead\n",  iCombo, rsc,
		    nFound);
	    }
#endif
	goto NEXTTRIAL;
      }

      // If we have the same reaction used for multiple components, we will get a singular matrix to invert
      // If we have different reactions used for multiple components, we may get a singular matrix
      // anyway, because the exit channels don't cover the full range space.
      // This needs to be discovered below.
      componentListTrial = componentListGlobal;
      for (int jj = 0; jj < (int) componentListGlobal.size() - 1; jj++) {
	int i1 = componentRxnListTrial[jj];
	for (int iii = jj+1; iii < (int) componentListGlobal.size(); iii++) {
	  if (i1 ==  componentRxnListTrial[iii]) {
	    componentRxnListTrial[iii] = -1;
	    componentListTrial[iii] = -1;
#ifdef DEBUG_MODE
	    if (m_debug_print_lvl >= 3) {
	      plogf("          combo %d disallowed because of duplicate exit reactions\n", iCombo);
	    }
#endif
	    goto NEXTTRIAL;
	  }
	}
      }

      /*
       *  componentMatrix (vector of length number of rows in the matrix) containing component indeces in the matrix problem
       */
      componentMatrix.clear();
      for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	if (componentListTrial[jj] >= 0) {
	  componentMatrix.push_back(componentListTrial[jj]);
	}
      }  
      /*
       *  componentRxnMatrix (vector of length number of cols in the matrix) containing component Rxn indeces in the matrix problem
       */
   
      componentRxnMatrix.clear();
      for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	if (componentRxnListTrial[jj] >= 0) {
	  componentRxnMatrix.push_back(componentRxnListTrial[jj]);
	}
      }

      
     
      vcs_dcopy(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_molNumSpecies_old), m_numSpeciesTot);
      vcs_dzero(VCS_DATA_PTR(m_deltaMolNumSpecies), m_numSpeciesTot);
      vcs_dzero(VCS_DATA_PTR(RxnDelta), m_numSpeciesTot);
      vcs_dcopy(VCS_DATA_PTR(m_tPhaseMoles_new), VCS_DATA_PTR(m_tPhaseMoles_old), m_numPhases);
      /*
       *  Calculate the mole numbers of all species after we have eliminated all of the noncomponents of phases
       *  that we want to eliminate
       */
      for (k = 0; k < Vphase->nSpecies(); k++) {
	kspec = Vphase->spGlobalIndexVCS(k);
	if (kspec >= m_numComponents) {
	  irxn = kspec - m_numComponents;
	  sc_irxn = m_stoichCoeffRxnMatrix[irxn];
	  dx = m_molNumSpecies_old[kspec];
	  m_molNumSpecies_new[kspec] = 0.0;
	  m_deltaMolNumSpecies[kspec] = -dx;
	  RxnDelta[irxn] = -dx;
	  int kph = m_phaseID[kspec];
	  m_tPhaseMoles_new[kph] -= dx;     
	  if (dx != 0.0) {
	    for (j = 0; j < m_numComponents; ++j) {
	      m_molNumSpecies_new[j] -= sc_irxn[j] * dx;
	      m_deltaMolNumSpecies[j] -= sc_irxn[j] * dx;   
	      int jph = m_phaseID[j];
	      m_tPhaseMoles_new[jph] -= sc_irxn[j] * dx;        
	    }
	  }
	}
      }
      if (rsc > 0) {

	/*
	 *  Formulate the matrix problem
	 */
	int j = 0;
	for (int j1 = 0; j1 < cs; j1++) {  
	  int jrxn = componentRxnListTrial[j1];
	  if (jrxn >= 0) {
	    double *sc_jrxn = m_stoichCoeffRxnMatrix[jrxn];
	    int i = 0;
	    for (int i1 = 0; i1 < cs; ++i1) {
	      irxn = componentRxnListTrial[i1];
	      if (irxn >= 0) {
		int icomp = componentListGlobal[i1];
		sm[i + j *rsc] = sc_jrxn[icomp];
		i++;
	      }   
	    }
	    int jcomp = componentListGlobal[j1];
	    rhs[j] = m_molNumSpecies_new[jcomp];
	    j++;
	  }	  
	}
	/*
	 *  Now calculate the rank
	 */
	vector<int> compRes;
	vector<int> elemComp;
	int usedZeroedSpecies;
	int rank = vcs_rank(rhs, rsc, DATA_PTR(sm), rsc, compRes, elemComp, &usedZeroedSpecies);
	vcs_heapsort(compRes);
	vcs_heapsort(elemComp);

#ifdef DEBUG_MODE
	if (m_debug_print_lvl >= 3) {
	  plogf("          combo %d has rank %d\n", iCombo, rank);
	}
#endif
	// reformulate the problem if the rank is less than rsc
	cs = rsc;
	if (rank < rsc) {

	  /*
	   *  Get rid of rows that lead to the deficient rank.
	   */
	  int j = 0;
	  for (int iii = 0; iii <(int) componentMatrix.size(); iii++) {
	    bool ifound = false;
	    for (int i = 0; i < (int) compRes.size(); i++) {
	      if (compRes[i] == iii) {
		ifound = true;
	      }
	    }
	    if (ifound) {
	      componentMatrix[j] = componentMatrix[iii];
	      j++;
	    }
	  }
	  componentMatrix.resize(compRes.size());
	  /*
	   *  Get rid of columns that lead to the deficient rank., because it was culled.
	   */
	  j = 0;
	  for (int iii = 0; iii < (int)componentRxnMatrix.size(); iii++) {
	    bool ifound = false;
	    for (int i = 0; i < (int) elemComp.size(); i++) {
	      if (elemComp[i] == iii) {
		ifound = true;
	      }
	    }
	    if (ifound) {
	      componentRxnMatrix[j] = componentRxnMatrix[iii];
	      j++;
	    }
	  }
	  componentRxnMatrix.resize(elemComp.size());

	  /*
	   *  Fix the componentListTrial vector, because it was culled.
	   */ 
	  for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	    int icomp = componentListTrial[jj];
	    if (icomp >= 0) {
	      int ifound = false;
	      for (int kk = 0; kk <(int)  componentMatrix.size(); kk++) {
		int kcomp = componentMatrix[kk];
		if (kcomp == icomp) {
		  ifound = true;
		}
	      }
	      if (!ifound) {
		componentListTrial[jj] = -1;
	      }
	    }
	  }
	  
	  /*
	   * Fix the componentRxnListTrial vector
	   */
	  for (int jj = 0; jj < (int) componentListGlobal.size(); jj++) {
	    int jrxn = componentRxnListTrial[jj];
	    if (jrxn >= 0) {
	      int ifound = false;
	      for (int kk = 0; kk <(int)  componentRxnMatrix.size(); kk++) {
		irxn = componentRxnMatrix[kk];
		if (jrxn == irxn) {
		  ifound = true;
		}
	      }
	      if (!ifound) {
		componentRxnListTrial[jj] = -1;
	      }
	    }
	  }

	  /*
	   * Reformulate the matrix
	   */
	  cs = componentRxnMatrix.size();
	  for (int j1 = 0; j1 < cs; j1++) {  
	    int jrxn = componentRxnMatrix[j1];
	    double *sc_jrxn = m_stoichCoeffRxnMatrix[jrxn];
	  
	    for (int i1 = 0; i1 < cs; ++i1) {
	      int icomp = componentMatrix[i1];
	      sm[i1 + j1 *cs] = sc_jrxn[icomp];
	    }
	  }
	 
	}

	int res = vcsUtil_mlequ(DATA_PTR(sm), cs, cs, rhs, 1);
	if (res == 1) {
	  printf ("res = %d\n", res);
	  exit(-1);
	}

	int i = 0;
	for (int i1 = 0; i1 < rank; ++i1) {
	  irxn = componentRxnMatrix[i1];
	  //	  irxn = componentRxnListTrial[i1];
	  if (irxn >= 0) {
	    RxnDelta[irxn] = rhs[i1];
	    i++;
	  }
	}

	for (int kph = 0; kph < m_numPhases; kph++) {
	  m_tPhaseMoles_new[kph] =  m_tPhaseMoles_old[kph];
	}
	for (kspec = 0; kspec < m_numSpeciesTot; kspec++) {
	  m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];

	}
	for (kspec = m_numComponents; kspec < m_numSpeciesTot; kspec++) {
	  irxn = kspec - m_numComponents;	 
	  sc_irxn = m_stoichCoeffRxnMatrix[irxn];
	  dx = RxnDelta[irxn];
	  m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
	  m_deltaMolNumSpecies[kspec] += dx;
	  int kph = m_phaseID[kspec];
	  m_tPhaseMoles_new[kph] += dx;     
	  if (dx != 0.0) {
	    for (j = 0; j < m_numComponents; ++j) {
	      m_molNumSpecies_new[j] += sc_irxn[j] * dx;
	      m_deltaMolNumSpecies[j] += sc_irxn[j] * dx;
	      int kph = m_phaseID[j];
	      m_tPhaseMoles_new[kph] += sc_irxn[j] * dx;     
	    }
	  }
	  
	}
      }
      /*
       *  Check for nonnegativity of proposed delta
       *        (HKM again cutoff values are important. We want values to be comfortably inside zero. here.
       *         Phases are only accurate to VCS_DELETE_PHASE_CUTOFF, which is based on being comfortably within the  roundoff dp roundoff
       *         limit. However, species are much more accurate. The 10-7 is really saying that we must have 7 digits
       *         of accuracy on the mole fraction for this routine to work.
       *         If we are not within these limits, it's fine, because this is an acceleration routine.) 
       *         (this algorithm failed once because of too loose tolerances here).
       */
      for (kspec = 0; kspec < m_numSpeciesTot; kspec++) {
	if (m_molNumSpecies_new[kspec] < -1.0E-7 * VCS_DELETE_PHASE_CUTOFF * m_totalMolNum) {
#ifdef DEBUG_MODE
	    if (m_debug_print_lvl >= 3) {
	      plogf("          Combo %d disallowed because species %-21.21s went negative %g \n",  iCombo, m_speciesName[kspec].c_str(),
		    m_molNumSpecies_new[kspec]);
	    }
#endif
	  goto NEXTTRIAL; 
	}
      }

      /*
       *  Check for zeroing of the phase. If it isn't zeroed, then give up on trial.
       *  This is part of the normal running of the routine. The range space may be deficient and the code below
       *  will pick this up.
       */
      for (k = 0; k < Vphase->nSpecies(); k++) {
	kspec = Vphase->spGlobalIndexVCS(k);
	if (fabs(m_molNumSpecies_new[kspec]) > 1.0E-7 * VCS_DELETE_PHASE_CUTOFF * m_totalMolNum) {
#ifdef DEBUG_MODE
	    if (m_debug_print_lvl >= 3) {
	      plogf("          Combo %d disallowed because species %-21.21s not zeroed %g \n",  iCombo, m_speciesName[kspec].c_str(),
		    m_molNumSpecies_new[kspec]);
	    }
#endif
	  goto NEXTTRIAL; 
	}
      }

      /*
       *  Now calculate the total Gibbs free energy of the new state.
       */
      vcs_setFlagsVolPhases(false, VCS_STATECALC_NEW);
      vcs_dfe(VCS_STATECALC_NEW, 0, 0, m_numSpeciesTot);
      total_new = vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_feSpecies_new), 
					 VCS_DATA_PTR(m_tPhaseMoles_new));
      deltaG = total_new - total_old;
#ifdef DEBUG_MODE
	    if (m_debug_print_lvl >= 3) {
	      plogf("          Combo %d produced DeltaG = %g total_new = %g total_old = %g\n", iCombo, deltaG, total_new, total_old);

	    }
#endif


      if (deltaG < lowestDelta) {
	bestCombo = iCombo;
	lowestDelta = deltaG;
	RxnDeltaBest = RxnDelta;
	componentListBest = 	componentListTrial;
	componentRxnListBest = 	componentRxnListTrial;
      }
      
    NEXTTRIAL:
      ;

    }
    /*
     * Get the best combo
     */
    if (bestCombo < 0) {
      goto RESTORE_OLD;
    }
 

   
    for (int kph = 0; kph < m_numPhases; kph++) {
      m_tPhaseMoles_new[kph] =  m_tPhaseMoles_old[kph];
    }
    for (kspec = 0; kspec < m_numSpeciesTot; kspec++) {
      m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
      m_deltaMolNumSpecies[kspec] = 0.0;

    }
    for (kspec = m_numComponents; kspec < m_numSpeciesTot; kspec++) {
      irxn = kspec - m_numComponents;	 
      sc_irxn = m_stoichCoeffRxnMatrix[irxn];
      dx = RxnDeltaBest[irxn];
      m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
      m_deltaMolNumSpecies[kspec] += dx;
      int kph = m_phaseID[kspec];
      m_tPhaseMoles_new[kph] += dx;     
      if (dx != 0.0) {
	for (j = 0; j < m_numComponents; ++j) {
	  m_molNumSpecies_new[j] += sc_irxn[j] * dx;
	  m_deltaMolNumSpecies[j] += sc_irxn[j] * dx;
	  int kph = m_phaseID[j];
	  m_tPhaseMoles_new[kph] += sc_irxn[j] * dx;     
	}
      }
	  
    }

#ifdef DEBUG_MODE
    for (kspec = 0; kspec < m_numSpeciesTot; kspec++) {
	if (m_molNumSpecies_new[kspec] < -1.0E-7 * VCS_DELETE_PHASE_CUTOFF * m_totalMolNum) {
	    if (m_debug_print_lvl >= 3) {
	      plogf("          Combo %d disallowed because species %-21.21s went negative %g \n",  bestCombo, m_speciesName[kspec].c_str(),
		    m_molNumSpecies_new[kspec]);
	    }
	    printf("Shouldn't be here\n");
	    exit(-1);
	}
      }

      /*
       *  Check for zeroing of the phase. If it isn't zeroed, then give up on trial.
       *  This is part of the normal running of the routine. The range space may be deficient and the code below
       *  will pick this up.
       */
      for (k = 0; k < Vphase->nSpecies(); k++) {
	kspec = Vphase->spGlobalIndexVCS(k);
	if (fabs(m_molNumSpecies_new[kspec]) > 1.0E-7 * VCS_DELETE_PHASE_CUTOFF * m_totalMolNum) {
	  if (m_debug_print_lvl >= 3) {
	    plogf("          Combo %d disallowed because species %-21.21s not zeroed %g \n",  bestCombo, m_speciesName[kspec].c_str(),
		  m_molNumSpecies_new[kspec]);
	  }	  
	  printf("Shouldn't be here\n");
	    exit(-1);
	    
	}
      }

    /*
     *  Now calculate the total Gibbs free energy of the new state.
     */
    vcs_setFlagsVolPhases(false, VCS_STATECALC_NEW);
    vcs_dfe(VCS_STATECALC_NEW, 0, 0, m_numSpeciesTot);
    total_new = vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_feSpecies_new), 
				VCS_DATA_PTR(m_tPhaseMoles_new));
    deltaG = total_new - total_old;
    if (deltaG < 0.0) {
      success = 1;
    } else {
      success = -1;
    }
    if (m_debug_print_lvl >= 5) {

      printf("   --- Subroutine vcs_popPhase_calcDeleteDelta Results: Getting rid of phase %d: %-15.15s\n", iph, (Vphase->PhaseName).c_str());
      printf("   ---    SPECIES          OLD_MOLES      RXN_DELTA        NEW_MOLES       DELTA_MOLES     Rxn_ID    \n");
      for (k = 0; k < m_numComponents; k++) {
	printf("   ---    %-15.15s % -11.4E                     % -11.4E  ", m_speciesName[k].c_str(),   m_molNumSpecies_old[k], m_molNumSpecies_new[k]);
	int kph = m_phaseID[k];
	int irxn = -1;
	printf("  % -11.4E  ",  m_deltaMolNumSpecies[k]);
	if (kph == iph) {
	  for (int i = 0; i < cs; ++i) {
	    int icomp = componentListGlobal[i];
	    if (icomp == k) {
	      irxn = componentRxnListBest[i];
	    }
	  }
	}
	if (irxn >= 0) {
	  printf("  %4d", irxn);
	} else {
	  printf("        ");
	}
	printf("\n");
      }
	
      for (kspec = m_numComponents; kspec < m_numSpeciesTot; kspec++) {
	printf("   ---    %-15.15s % -11.4E ",m_speciesName[kspec].c_str(), m_molNumSpecies_old[kspec]);
	int kph = m_phaseID[kspec];
	irxn = kspec -  m_numComponents;
	if (kph == iph) {
	  printf("   % -11.4E   ",   RxnDeltaBest[irxn]);
	} else {
	  bool found = false;
	  for (int i = 0; i < cs; ++i) {
	    int icomp = componentListGlobal[i];
	   
	    int iirxn = componentRxnListBest[i];
	    if (iirxn == irxn) {
	      printf("%3d:% -11.4E  ", icomp, RxnDelta[iirxn]);
	      found = true;
	    }
	      
	  }
	  if (!found) {
	    printf("                 ");
	  }

	}
	printf("   % 11.4E  ",  m_molNumSpecies_new[kspec]);
	printf("  % -11.4E  ",  m_deltaMolNumSpecies[kspec]);
	kph = m_phaseID[kspec];
	if (kph == iph) {
	  printf("  %4d", irxn);
	} else {
	  printf("        ");
	}
	printf("\n");
      }
      printf("   ---    TotalG          %-13.7E                  %-13.7E   (%13.7E)\n", total_old, total_new, deltaG);
    }
    if (m_debug_print_lvl >= 2) {
      if (deltaG >= 0.0) {
	printf("   --- vcs_popPhase_calcDeleteDelta() Overall result for phase %d is to RETAIN PHASE\n", iph);
      } else {
	printf("   --- vcs_popPhase_calcDeleteDelta() Overall result for phase %d is to DELETE PHASE\n", iph);
      }
    }
    ttt = false;
    for (kspec = 0; kspec < m_numSpeciesTot; kspec++) {
      if (m_molNumSpecies_new[kspec] < -1.0E-12) {
	printf(" Error in update vector m_molNumSpecies_new[%d] = %g %s\n", kspec, m_molNumSpecies_new[kspec], m_speciesName[kspec].c_str() );
	ttt = true;
      }
      int kph = m_phaseID[kspec];
      if (kph == iph) {
	if (fabs(m_molNumSpecies_new[kspec]) > 1.0E-14) {
	  printf(" Error in update vector - should have been zereod m_molNumSpecies_new[%d] = %g %s\n", kspec, m_molNumSpecies_new[kspec], m_speciesName[kspec].c_str() );
	  ttt = true;
	}
      }
    }
    if (ttt) {
      exit(-1);
    }


#endif
  RESTORE_OLD:
    /*
     *  needed to set the phase existence flags back again >!>
     */
    vcs_setFlagsVolPhases(false, VCS_STATECALC_OLD);
    vcs_dfe(VCS_STATECALC_OLD, 0, 0, m_numSpeciesTot);
    return success;
  }

}
//======================================================================================================================
