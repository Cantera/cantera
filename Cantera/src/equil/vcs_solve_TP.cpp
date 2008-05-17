/**
 * @file vcs_solve_TP.cpp Implementation file that contains the
 *     main algorithm for finding an equilibrium
 */
/*
 * $Id$
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "vcs_solve.h"
#include "vcs_internal.h" 
#include "vcs_VolPhase.h"
#include "vcs_species_thermo.h"

#include "clockWC.h"

#ifdef WIN32
#pragma warning(disable:4996)
#endif

using namespace std;

namespace VCSnonideal {


  /***************************************************************************/
  /************ Prototypes for static functions ******************************/


  static void print_space(int num);


#ifdef DEBUG_MODE 
  //static double minor_alt_calc(int, int, int *, char *); 
#else
  //static double minor_alt_calc(int, int, int *);
#endif
#ifdef DEBUG_MODE
#  ifdef DEBUG_NOT
  static void prneav(void);
  static int  prnfm(void);
#  endif 
#endif
  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

#ifdef DEBUG_MODE
  void VCS_SOLVE::checkDelta1(double * const dsLocal, 
			      double * const delTPhMoles, int kspec) {
    std::vector<double> dchange(m_numPhases, 0.0);
    for (int k = 0; k < kspec; k++) {
      if (m_speciesUnknownType[k] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	int iph = m_phaseID[k];
	dchange[iph] += dsLocal[k];
      }
    }
    for (int iphase = 0; iphase < m_numPhases; iphase++) {
      double denom = MAX(m_totalMolNum, 1.0E-4);
      if (!vcs_doubleEqual(dchange[iphase]/denom, delTPhMoles[iphase]/denom)) {
	plogf("checkDelta1: we have found a problem\n");
	exit(-1);
      }
    }
  }
#endif
  /*****************************************************************************/

  // Main routine that solves for equilibrium at constant T and P
  // using a variant of the VCS method
  /*
   * This is the main routine  taht solves for equilibrium at constant T and P
   * using a variant of the VCS method. Nonideal phases can be accommodated
   * as well.
   *
   * Any number of single-species phases and  multi-species phases 
   * can be handled by the present version.
   *
   *     Input
   * ------------
   *   @param print_lvl     1 -> Print results to standard output 
   *                        0 -> don't report on anything 
   *
   *   @param printDetails  1 -> Print intermediate results. 
   *
   *   @param maxit         Maximum number of iterations for the algorithm 
   *
   *   @return     0 = Equilibrium Achieved
   *               1 = Range space error encountered. The element abundance criteria are
   *                   only partially satisfied. Specifically, the first NC= (number of
   *                   components) conditions are satisfied. However, the full NE 
   *                   (number of elements) conditions are not satisfied. The equilibrirum
   *                   condition is returned.
   *              -1 = Maximum number of iterations is exceeded. Convergence was not
   *                   found.
   */
  int VCS_SOLVE::vcs_solve_TP(int print_lvl, int printDetails, int maxit) {
    int conv = FALSE, retn = VCS_SUCCESS;
    double test, RT;
    int j, k, l, solveFail, l1, kspec, irxn, im, forced, iph;
    double dx, xx, par;
    int liqphase = FALSE, numSpecliquid = 0;
    int dofast, soldel, ll, it1;
    int lec, npb, iti, i, lnospec;
    int rangeErrorFound = 0;
    bool giveUpOnElemAbund = false;
    int finalElemAbundAttempts = 0;
    bool MajorSpeciesHaveConverged = false;
    int uptodate_minors = TRUE;
    bool justDeletedMultiPhase = FALSE;
    int usedZeroedSpecies; /* return flag from basopt indicating that 
			      one of the components had a zero concentration */

    vcs_VolPhase *Vphase;
    double     *sc_irxn = NULL;  /* Stoichiometric coefficients for cur rxn  */
    double *dnPhase_irxn;
#ifdef DEBUG_MODE
    char ANOTE[128];
    /*
     * Set the debug print lvl to the same as the print lvl.
     */
    vcs_debug_print_lvl = printDetails;
#endif
    if (printDetails > 0 && print_lvl == 0) { 
      print_lvl = 1;
    }
    /*
     *    Initialize and set up all counters
     */
    vcs_counters_init(0);
    Cantera::clockWC ticktock;
   
    /*
     *  Malloc temporary space for usage in this routine and in
     *  subroutines
     *        sm[ne*ne]
     *        ss[ne]
     *        sa[ne]
     *        aw[m]
     *        wx[ne]
     *        xy[m]
     */
 

    std::vector<double> sm(m_numElemConstraints*m_numElemConstraints, 0.0);
    std::vector<double> ss(m_numElemConstraints, 0.0);
    std::vector<double> sa(m_numElemConstraints, 0.0);

    std::vector<double> aw(m_numSpeciesTot, 0.0);
    std::vector<double> wx(m_numElemConstraints, 0.0);
   
    solveFail = FALSE;
    im = FALSE;
   
    /* ****************************************************** */
    /* **** Evaluate the elemental composition         ****** */
    /* ****************************************************** */
    vcs_elab();
   
    /* ******************************************************* */
    /* **** Printout the initial conditions for problem ****** */
    /* ******************************************************* */
    if (m_numPhases > 1) {
      if (! m_VolPhaseList[1]->SingleSpecies) {
	liqphase = TRUE;
	numSpecliquid = m_VolPhaseList[1]->NVolSpecies;
      }
    }
    if (print_lvl != 0) {
      plogf("VCS CALCULATION METHOD\n\n ");
      plogf("%s\n", m_title.c_str());
      plogf("\n\n%5d SPECIES%8d ELEMENTS", m_numSpeciesTot, m_numElemConstraints);
      plogf("%16d COMPONENTS\n%5d PHASE1 SPECIES", m_numComponents,
	    ((m_VolPhaseList[0])->NVolSpecies));
      plogf("%10d PHASE2 SPECIES%8d SINGLE SPECIES PHASES\n\n", 
	    numSpecliquid, 
	    m_numSpeciesTot - (m_VolPhaseList[0])->NVolSpecies - numSpecliquid);  
      //string punits = "atm";
      //if (m_VCS_UnitsFormat == 3) {
      //	punits = "Pa ";
      //}
      plogf(" PRESSURE%22.8g %3s\n", m_pressurePA, "Pa "); 
      plogf(" TEMPERATURE%19.3f K\n", m_temperature);
      Vphase = m_VolPhaseList[0];
      if (Vphase->NVolSpecies > 0) {
	plogf(" PHASE1 INERTS%17.3f\n", TPhInertMoles[0]);
      }
      if (liqphase) {
	plogf(" PHASE2 INERTS%17.3f\n", TPhInertMoles[1]);
      }
      plogf("\n ELEMENTAL ABUNDANCES             CORRECT");
      plogf("          FROM ESTIMATE           Type\n\n");
      for (i = 0; i < m_numElemConstraints; ++i) { 
	print_space(26); plogf("%-2.2s", (m_elementName[i]).c_str());
	plogf("%20.12E%20.12E     %3d\n", m_elemAbundancesGoal[i], m_elemAbundances[i], 
	      m_elType[i]);
      }
      if (m_doEstimateEquil < 0) {
	plogf("\n MODIFIED LINEAR PROGRAMMING ESTIMATE OF EQUILIBRIUM - forced\n");
      } else if (m_doEstimateEquil > 0) {
	plogf("\n MODIFIED LINEAR PROGRAMMING ESTIMATE OF EQUILIBRIUM - where necessary\n");
      }
      if (m_doEstimateEquil == 0) {
	plogf("\n USER ESTIMATE OF EQUILIBRIUM\n");
      }
      if (m_VCS_UnitsFormat == VCS_UNITS_KCALMOL) {
	plogf(" Stan. Chem. Pot. in kcal/mole\n");
      }
      if (m_VCS_UnitsFormat == VCS_UNITS_UNITLESS) {
	plogf(" Stan. Chem. Pot. is MU/RT\n");
      }
      if (m_VCS_UnitsFormat == VCS_UNITS_KJMOL) {
	plogf(" Stan. Chem. Pot. in KJ/mole\n");
      } 
      if (m_VCS_UnitsFormat == VCS_UNITS_KELVIN) {
	plogf(" Stan. Chem. Pot. in Kelvin\n");
      } 
      if (m_VCS_UnitsFormat == VCS_UNITS_MKS) {
	plogf(" Stan. Chem. Pot. in J/kmol\n");
      }
      plogf("\n SPECIES       FORMULA VECTOR");
      print_space(29);
      plogf("   STAN_CHEM_POT   EQUILIBRIUM_EST.  Species_Type\n\n");
      print_space(14);
      for (i = 0; i < m_numElemConstraints; ++i) plogf(" %-2.2s", m_elementName[i].c_str());
      plogf(" SI(I)\n");
      RT = vcs_nondimMult_TP(m_VCS_UnitsFormat, m_temperature);
      for (i = 0; i < m_numSpeciesTot; ++i) {
	plogf(" %-12s", m_speciesName[i].c_str());
	for (j = 0; j < m_numElemConstraints; ++j) {
	  plogf("%3g", m_formulaMatrix[j][i]);
	}
	if  (m_phaseID[i] == 0) {
	  plogf("  1");
	} else if (m_phaseID[i] == 1) {
	  if (liqphase) plogf("  2");
	  else          plogf("  0");
	} else { 
	  plogf("  0");
	}
	print_space(47-m_numElemConstraints*3);
	plogf("%12.5E  %12.5E", RT * m_SSfeSpecies[i], m_molNumSpecies_old[i]);
	if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
	  plogf("       Mol_Num");
	} else if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  plogf("       Voltage");
	} else {
	  plogf("       Unknown");
	}
	plogendl();
      }
    }

    for (i = 0; i < m_numSpeciesTot; ++i) {
      if (m_molNumSpecies_old[i] < 0.0) {
	plogf("On Input species %-12s has a "
	      "negative MF, setting it small",
	      m_speciesName[i].c_str());	
	plogendl();
	m_molNumSpecies_old[i] = VCS_DELETE_SPECIES_CUTOFF;
      }
    }

   
    /* 
     *  Evaluate the total moles of species in the problem
     */
    vcs_tmoles();

    /* ***************************************************************************** */
    /* **** EVALUATE ALL CHEMICAL POTENTIALS AT THE OLD (CURRENT) MOLE NUMBERS ***** */
    /* ***************************************************************************** */
    vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);

    /*
     *  HKM -> If there was a machine estimate, we used to branch
     *         to the code segment which determined whether we needed a
     *         new component basis. If we did, we would go to L429.
     *         If we didn't, we would go to a point below basopt() below.
     *         I have taken this section out of the code for simplicity's
     *         sake. It's not need for speed, since in any recursive
     *         call to this subroutine we would have an initial estimate
     *         of the solution. And, we don't need to optimize the
     *         startup of nonrecursive calls to this subroutine.
     */
    /* *********************************************************** */
    /* **** DETERMINE BASIS SPECIES, EVALUATE STOICHIOMETRY ****** */
    /* *********************************************************** */
    /*
     *   This is an entry point for later in the calculation 
     */
  L_COMPONENT_CALC: ;
    test = -1.0e-10;
    retn = vcs_basopt(FALSE, VCS_DATA_PTR(aw), VCS_DATA_PTR(sa),
		      VCS_DATA_PTR(sm), VCS_DATA_PTR(ss), 
		      test, &usedZeroedSpecies);
    if (retn != VCS_SUCCESS) return retn;

    if (conv) {
      goto L_RETURN_BLOCK;
    }
    it1 = 1;
    MajorSpeciesHaveConverged = false;
    /*************************************************************************/
    /************** EVALUATE INITIAL MAJOR-MINOR VECTOR **********************/
    /*************************************************************************/
    m_numRxnMinorZeroed = 0;
    for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
      kspec = m_indexRxnToSpecies[irxn];
      m_rxnStatus[irxn] = vcs_species_type(kspec);
      if (m_rxnStatus[irxn] == VCS_SPECIES_MINOR) {
	m_rxnStatus[irxn] = VCS_SPECIES_MAJOR;
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  plogf("   --- Minor species changed to major: ");
	  plogf("%-12s", m_speciesName[kspec].c_str());
	  plogendl();
	}
#endif
      }
      if (m_rxnStatus[irxn] != VCS_SPECIES_MAJOR) {
	++m_numRxnMinorZeroed;
      }
    }
    im = (m_numRxnMinorZeroed == m_numRxnRdc);
    lec = FALSE;
    if (! vcs_elabcheck(0)) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- Element Abundance check failed");
 	plogendl();
      }
#endif
      vcs_elcorr(VCS_DATA_PTR(sm), VCS_DATA_PTR(wx));
      vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
    }
#ifdef DEBUG_MODE	
    else {
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- Element Abundance check passed");	
	plogendl();
      }
    }
#endif
    // Update the phase objects with the contents of the soln vector
    vcs_updateVP(0);
    vcs_deltag(0, false);
    iti = 0;
    goto L_MAINLOOP_ALL_SPECIES;

    /* ********************************************************* */
    /* **** SET INITIAL VALUES FOR ITERATION ******************* */
    /* **** EVALUATE REACTION ADJUSTMENTS    ******************* */
    /* ********************************************************* */
    /*
     *  This is the top of the loop ----------------------------------------
     *  Every 4th iteration ITI = 0. Else, It's equal to a negative number
     */
  L_MAINLOOP_MM4_SPECIES: ;
    iti = ((it1/4) *4) - it1;
    /*
     *      Entry point when the code wants to force an ITI=0 calculation 
     */
  L_MAINLOOP_ALL_SPECIES: ;
    if (iti == 0) {
      /* 
       *          Evaluate the minor non-componenent species chemical 
       *          potentials and delta G for their formation reactions 
       *          We have already evaluated the major non-components 
       */
      if (uptodate_minors == FALSE) {
	vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
	vcs_deltag(1, false);
      }
      uptodate_minors = TRUE;
    } else {
      uptodate_minors = FALSE;
    }

    if (printDetails) {
      plogf("\n"); vcs_print_line("=", 110); 
      plogf(" Iteration = %3d, Iterations since last evaluation of "
	    "optimal basis = %3d",
	    m_VCount->Its, it1 - 1);
      if (iti == 0) {
	plogf(" (all species)\n");
      } else {
	plogf(" (only major species)\n");
      }
    }

    vcs_dcopy(VCS_DATA_PTR(m_feSpecies_old), VCS_DATA_PTR(m_feSpecies_curr), m_numSpeciesRdc);
    vcs_dcopy(VCS_DATA_PTR(m_feSpecies_new), VCS_DATA_PTR(m_feSpecies_curr), m_numSpeciesRdc);
    vcs_dcopy(VCS_DATA_PTR(m_actCoeffSpecies_old), VCS_DATA_PTR(m_actCoeffSpecies_new), m_numSpeciesRdc);
    vcs_dcopy(VCS_DATA_PTR(m_deltaGRxn_old), VCS_DATA_PTR(m_deltaGRxn_new), m_numRxnRdc);
 
    /*        Go find a new reaction adjustment -> 
     *         i.e., change in extent of reaction for each reaction. 
     *
     *     Zero out the entire vector of updates. We sometimes would
     *     query these values below, and we want to be sure that no
     *     information is left from previous iterations.
     */
    vcs_dzero(VCS_DATA_PTR(m_deltaMolNumSpecies), m_numSpeciesTot);
    /*
     * Figure out whether we will calculate new reaction step sizes
     * for the major species.
     *    -> We won't if all species are minors (im), OR
     *       all major species have already converged 
     */
    if (!(MajorSpeciesHaveConverged) && ! im) {
      soldel = vcs_RxnStepSizes();
      /* -         If SOLDEL is true then we encountered a reaction between */
      /* -         single-species-phase species, only, and have adjusted */
      /* -         the mole number vector, W(), directly. In this case, */
      /* -         we should immediately go back and recompute a new */
      /* -         component basis, if the species that was zeroed was */
      /* -         a component. SOLDEL is true when this is so. */
      if (soldel > 0) {
	/* -           We have changed the base mole number amongst single- */
	/* -           species-phase species. However, we don't need to */
	/* -           recaculate their chemical potentials because they */
	/* -           are constant, anyway! */
	if (soldel == 2) {
	  goto L_COMPONENT_CALC;
	}
	/* -           We have not changed the actual DG values for */
	/* -           any species, even the one we deleted. Thus, */
	/* -           we don't need to start over. */
      }
    } else {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	if (im) {
	  plogf("   --- vcs_RxnStepSizes not called because all"
		"species are minors\n");
	} else {
	  plogf("   --- vcs_RxnStepSizes not called because "
		"all majors have converged\n");
	}
      }
#endif
    }
   
    lec = FALSE;
    /*
     *    Zero out the net change in moles of multispecies phases 
     */
    vcs_dzero(VCS_DATA_PTR(m_deltaPhaseMoles), m_numPhases);

    /*
     * Check on too many iterations. 
     *         If we have too many iterations, Clean up and exit code even though we haven't 
     *         converged. -> we have run out of iterations! 
     */
    if (m_VCount->Its > maxit) {
      solveFail = -1;
      goto L_RETURN_BLOCK;
    }

    /* ********************************************************************** */
    /* ***************** MAIN LOOP IN CALCULATION *************************** */
    /* ***************** LOOP OVER IRXN TO DETERMINE STEP SIZE ************** */
    /* ********************************************************************** */
    /*
     *   Loop through all of the reactions, irxn, pertaining to the
     *   formation reaction for species kspec in canonical form.
     *
     *   At the end of this loop, we will have a new estimate for the
     *   mole numbers for all species consistent with an extent
     *   of reaction for all noncomponent species formation
     *   reactions. We will have also ensured that all predicted
     *   non-component mole numbers are greater than zero.
     *
     *         Old_Solution               New_Solution             Description
     * -----------------------------------------------------------------------------
     *  m_molNumSpecies_old[kspec]   m_molNumSpecies_new[kspec]   Species Mole Numbers
     *                               m_deltaMolNumSpecies[kspec]  Delta in the Species Mole Numbers
     *
     *
     *
     */  
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Main Loop Treatment of each non-component species ");
      if (iti == 0) plogf("- Full Calculation:\n");
      else          plogf("- Major Components Calculation:\n");
      plogf("   --- Species     IC    ");
      plogf(" KMoles  Tent_KMoles Rxn_Adj   |    Comment \n");
    }
#endif
   
    for (irxn = 0; irxn < m_numRxnRdc; irxn++) {
      kspec = m_indexRxnToSpecies[irxn];
      sc_irxn = m_stoichCoeffRxnMatrix[irxn];
      iph = m_phaseID[kspec];
      Vphase = m_VolPhaseList[iph];
#ifdef DEBUG_MODE
      ANOTE[0] = '\0';	 
#endif
    
      if (m_rxnStatus[irxn] == VCS_SPECIES_INTERFACIALVOLTAGE) {
	/********************************************************************/
	/************************ VOLTAGE SPECIES ***************************/
	/********************************************************************/
#ifdef DEBUG_MODE	 
	dx = minor_alt_calc(kspec, irxn, &soldel, ANOTE); 
#else
	dx = minor_alt_calc(kspec, irxn, &soldel);
#endif
	m_deltaMolNumSpecies[kspec] = dx;
      }
      else if (m_rxnStatus[irxn] < VCS_SPECIES_MINOR) {
	/********************************************************************/
	/********************** ZEROED OUT SPECIES **************************/
	/********************************************************************/
	bool resurrect = true;
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 3) {
	  plogf("   --- %s currently zeroed (SpStatus=%-2d):", 
		m_speciesName[kspec].c_str(), m_rxnStatus[irxn]);
	  plogf("%3d DG = %11.4E WT = %11.4E W = %11.4E DS = %11.4E\n",
		irxn, m_deltaGRxn_new[irxn], m_molNumSpecies_new[kspec],
		m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec]);
	}
#endif
	// HKM Alternative is to not allow ds[] = 0.0 phases
	// to pop back into existence. For esthetics, I'm allowing this.
	// so that dg < 0.0 phases with zero mole numbers become components.
	// This is also better, because that component will be the first
	// one to pop into existence if there is a minute quantity of the element.
	// This could change in the future.
	//if (dg[irxn] >= 0.0 || ds[kspec] <= 0.0) {
	if (m_deltaGRxn_new[irxn] >= 0.0 ) {
	  m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
	  m_deltaMolNumSpecies[kspec] = 0.0;
	  resurrect = false;
#ifdef DEBUG_MODE
	  sprintf(ANOTE, "Species stays zeroed: DG = %11.4E",
		  m_deltaGRxn_new[irxn]);
	  if (m_deltaGRxn_new[irxn] < 0.0) {
	    sprintf(ANOTE, "Species stays zeroed even though dg neg:DG = %11.4E, ds zeroed ",
		    m_deltaGRxn_new[irxn]);
	  }
	  //if (vcs_debug_print_lvl >= 2) {
	  //plogf("   --- "); plogf("%-12s", m_speciesName[kspec]);
	  //plogf("%3d%11.4E%11.4E%11.4E | %s\n", 
	  // m_rxnStatus[irxn], w[kspec], wt[kspec],
	  // ds[kspec], ANOTE);
	  //}
#endif
	} else {
	  for (int j = 0; j < m_numElemConstraints; ++j) {
	    int elType = m_elType[j];
	    if (elType == VCS_ELEM_TYPE_ABSPOS) {
	      double atomComp = m_formulaMatrix[j][kspec];
	      if (atomComp > 0.0) {
		double maxPermissible = m_elemAbundancesGoal[j] / atomComp;
		if (maxPermissible < VCS_DELETE_MINORSPECIES_CUTOFF) {
#ifdef DEBUG_MODE
		  sprintf(ANOTE, "Species stays zeroed even though dG neg, because of %s elemAbund",
			  m_elementName[j].c_str());
#endif
		  resurrect = false;
		  break;
		}
	      }
	    }
	  }
	}
	/*
	 * Resurrect the species
	 */ 
	if (resurrect) {
	  if (Vphase->Existence == 0) Vphase->Existence = 1;
	  --m_numRxnMinorZeroed;
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {
	    plogf("   --- Zeroed species changed to major: ");
	    plogf("%-12s\n", m_speciesName[kspec].c_str());
	  }
#endif
	  m_rxnStatus[irxn] = VCS_SPECIES_MAJOR;
	  im = FALSE;
	  MajorSpeciesHaveConverged = false;
	  if (m_deltaMolNumSpecies[kspec] > 0.0) {
	    dx = m_deltaMolNumSpecies[kspec] * 0.01;
	  
	    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
	  } else {
	    m_molNumSpecies_new[kspec] = m_totalMolNum * VCS_DELETE_PHASE_CUTOFF * 10.;
	    dx = m_molNumSpecies_new[kspec] - m_molNumSpecies_old[kspec];
	  }
	  m_deltaMolNumSpecies[kspec] = dx;
#ifdef DEBUG_MODE
	  sprintf(ANOTE, "Born:IC=-1 to IC=1:DG=%11.4E", m_deltaGRxn_new[irxn]);
#endif
	} else {
	  m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
	  m_deltaMolNumSpecies[kspec] = 0.0;
	  dx = 0.0;
	}
      } else if (m_rxnStatus[irxn] == VCS_SPECIES_MINOR) {
	/********************************************************************/
	/***************************** MINOR SPECIES ************************/
	/********************************************************************/
	/* 
	 *    Unless ITI isn't equal to zero we zero out changes 
	 *    to minor species. 
	 */
	if (iti != 0) {
	  m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
	  m_deltaMolNumSpecies[kspec] = 0.0;
	  dx = 0.0;
#ifdef DEBUG_MODE
	  sprintf(ANOTE,"minor species not considered");
	  if (vcs_debug_print_lvl >= 2) {
	    plogf("   --- "); plogf("%-12s", m_speciesName[kspec].c_str());
	    plogf("%3d%11.4E%11.4E%11.4E | %s", 
		  m_rxnStatus[irxn], m_molNumSpecies_old[kspec], m_molNumSpecies_new[kspec],
		  m_deltaMolNumSpecies[kspec], ANOTE);
	    plogendl();
	  }
#endif
	  continue;
	}
	/*
	 *        Minor species alternative calculation 
	 *       --------------------------------------- 
	 *    This is based upon the following approximation: 
	 *    The mole fraction changes due to these reactions don't affect 
	 *    the mole numbers of the component species. Therefore the 
	 *    following approximation is valid for an ideal solution 
	 *       0 = DG(I) + log(WT(I)/W(I))
	 *       (DG contains the contribution from FF(I) + log(W(I)/TL) ) 
	 *    Thus, 
	 *        WT(I) = W(I) EXP(-DG(I)) 
	 *    If soldel is true on return, then we branch to the section
	 *    that deletes a species from the current set of active species.
	 */
#ifdef DEBUG_MODE	 
	dx = minor_alt_calc(kspec, irxn, &soldel, ANOTE); 
#else
	dx = minor_alt_calc(kspec, irxn, &soldel);
#endif
	m_deltaMolNumSpecies[kspec] = dx;
	if (soldel) {
	  /*******************************************************************/
	  /*****  DELETE MINOR SPECIES LESS THAN  VCS_DELETE_SPECIES_CUTOFF  */
	  /*****  MOLE NUMBER                                                */
	  /*******************************************************************/
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {
	    plogf("   --- Delete minor species in multispec phase: %-12s",
		  m_speciesName[kspec].c_str());
	    plogendl();
	  }
#endif
	  m_deltaMolNumSpecies[kspec] = 0.0;
	  /*
	   *       Delete species, kspec. The alternate return is for the case
	   *       where all species become deleted. Then, we need to 
	   *       branch to the code where we reevaluate the deletion 
	   *       of all species.
	   */
	  lnospec = delete_species(kspec);
	  if (lnospec) goto L_RECHECK_DELETED;
	  /*
	   *       Go back to consider the next species in the list.
	   *       Note, however, that the next species in the list is now 
	   *       in slot l. In deleting the previous species L, We have 
	   *       exchanged slot MR with slot l, and then have 
	   *       decremented MR. 
	   *       Therefore, we will decrement the species counter, here.
	   */
	  --irxn;
#ifdef DEBUG_MODE
	  goto L_MAIN_LOOP_END_NO_PRINT;
#else
	  goto L_MAIN_LOOP_END;
#endif	    
	}
      } else {
	/********************************************************************/
	/*********************** MAJOR SPECIES ******************************/
	/********************************************************************/
#ifdef DEBUG_MODE
	sprintf(ANOTE, "Normal Major Calc");
#endif
	/*
	 * Check for superconvergence of the formation reaction. Do 
	 * nothing if it is superconverged. Skip to the end of the
	 * irxn loop if it is superconverged.
	 */
	if (fabs(m_deltaGRxn_new[irxn]) <= m_tolmaj2) {
	  m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
	  m_deltaMolNumSpecies[kspec] = 0.0;
	  dx = 0.0;
#ifdef DEBUG_MODE
	  sprintf(ANOTE, "major species is converged");
	  if (vcs_debug_print_lvl >= 2) {
	    plogf("   --- "); plogf("%-12s", m_speciesName[kspec].c_str());
	    plogf("%3d%11.4E%11.4E%11.4E | %s", 
		  m_rxnStatus[irxn], m_molNumSpecies_old[kspec], m_molNumSpecies_new[kspec],
		  m_deltaMolNumSpecies[kspec], ANOTE);
	    plogendl();
	  }
#endif
	  continue;
	}
	/*
	 *      Set the initial step size, dx, equal to the value produced
	 *      by the routine, vcs_RxnStepSize().
	 *
	 *          Note the multiplition logic is to make sure that
	 *          dg[] didn't change sign due to w[] changing in the
	 *          middle of the iteration. (it can if a single species
	 *          phase goes out of existence).
	 */
	if ((m_deltaGRxn_new[irxn] * m_deltaMolNumSpecies[kspec]) <= 0.0) {
	  dx = m_deltaMolNumSpecies[kspec];
	} else {
	  dx = 0.0;
	  m_deltaMolNumSpecies[kspec] = 0.0;
#ifdef DEBUG_MODE
	  sprintf(ANOTE, "dx set to 0, DG flipped sign due to "
		  "changed initial point");
#endif
	}
	/*
	 *      Form a tentative value of the new species moles 
	 */
	m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;

	/*
	 *      Check for non-positive mole fraction of major species.
	 *      If we find one, we branch to a section below. Then,
	 *      depending upon the outcome, we branch to sections below,
	 *      or we restart the entire iteration.
	 */
	if (m_molNumSpecies_new[kspec] <= 0.0) {
#ifdef DEBUG_MODE
	  sprintf(ANOTE, "initial nonpos kmoles= %11.3E",
		  m_molNumSpecies_new[kspec]);
#endif
	  /* ************************************************* */
	  /* *** NON-POSITIVE MOLES OF MAJOR SPECIES ********* */
	  /* ************************************************* */
	  /*
	   *          We are here when a tentative value of a mole fraction 
	   *          created by a tentative value of M_DELTAMOLNUMSPECIES(*) is negative. 
	   *          We branch from here depending upon whether this
	   *          species is in a single species phase or in 
	   *          a multispecies phase.
	   */
	  if (! (m_SSPhase[kspec])) {
	    /* 
	     *   Section for multispecies phases:
	     *     - Cut reaction adjustment for positive kmoles of 
	     *       major species in multispecies phases.
	     *       Decrease its concentration by a factor of 10.
	     */
	    dx = -0.9 * m_molNumSpecies_old[kspec];
	    m_deltaMolNumSpecies[kspec] = dx;
	    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
	    /*
	     *        Change major to minor if the current species
	     *        has a mole number that is less than 1/100 of the
	     *        total kmoles in the problem.
	     *        However, it also has to be a small species within its
	     *        own phase as well.
	     *        we can't call vcs_species_type() because the phase moles
	     *        would be wrong.
	     */
	    if (m_molNumSpecies_new[kspec] < 0.005 * m_totalMolNum) {
	      iph = m_phaseID[kspec];
	      if (m_molNumSpecies_new[kspec] < (m_tPhaseMoles_old[iph] * 0.01)) {
#ifdef DEBUG_MODE
		if (vcs_debug_print_lvl >= 2) {
		  plogf("   --- Major species changed to minor: ");
		  plogf("%-12s", m_speciesName[kspec].c_str());
		  plogendl();
		}
#endif
		m_rxnStatus[irxn] = VCS_SPECIES_MINOR;
		++m_numRxnMinorZeroed;
		im = (m_numRxnMinorZeroed == m_numRxnRdc);
	      }
	    }
	  } else {
	    /* 
	     *   Section for single species phases:
	     *       Calculate a dx that will wipe out the 
	     *       moles in the phase.
	     */
	    dx = -m_molNumSpecies_old[kspec];
	    /*
	     *       Calculate an update that doesn't create a negative mole
	     *       number for a component species. Actually, restrict this 
	     *       a little more so that the component values can only be
	     *       reduced by two 99%,
	     */
	    for (j = 0; j < m_numComponents; ++j) {
	      if (sc_irxn[j] != 0.0) {
		wx[j] = m_molNumSpecies_old[j] + sc_irxn[j] * dx;
		if (wx[j] <= m_molNumSpecies_old[j] * 0.01 - 1.0E-150) {
		  dx = MAX(dx,  m_molNumSpecies_old[j] * -0.99 / sc_irxn[j]);
		}
	      } else {
		wx[j] = m_molNumSpecies_old[j];
	      }
	    }
	    m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + dx;
	    if (m_molNumSpecies_new[kspec] > 0.0) {
	      m_deltaMolNumSpecies[kspec] = dx;
#ifdef DEBUG_MODE
	      sprintf(ANOTE, 
		      "zeroing SS phase created a neg component species "
		      "-> reducing step size instead");
#endif 
	    } else {
	      /*
	       *     We are going to zero the single species phase.
	       *     Set the existence flag
	       */
	      iph = m_phaseID[kspec];
	      Vphase = m_VolPhaseList[iph];
	      Vphase->Existence = 0;
#ifdef DEBUG_MODE
	      sprintf(ANOTE, "zero SS phase: moles went neg");
#endif
	      /*
	       *     Change the base mole numbers for the iteration.
	       *     We need to do this here, because we have decided 
	       *     to eliminate  the phase in this special section 
	       *     outside the main loop.
	       */
	      m_molNumSpecies_old[kspec] = 0.0;
	      for (j = 0; j < m_numComponents; ++j) {
		m_molNumSpecies_old[j] = wx[j];
	      }
	      /*
	       *     Change the total number of moles in all phases due to
	       *     the reaction that wil be zeroing out the pure species
	       *     phase. Make sure the moles in the current ss phase is
	       *     identically zero.
	       */
	      dnPhase_irxn = m_deltaMolNumPhase[irxn];
	      for (int iphase = 0; iphase < m_numPhases; iphase++) {
		m_tPhaseMoles_old[iphase] += dnPhase_irxn[iphase] * dx;
	      }
	      m_tPhaseMoles_old[iph] = 0.0;
	      vcs_updateVP(0);
	      /*
	       *       Recalcuate the chemical potentials, FE(), and the 
	       *       reaction free energy changes, DG(), for the current
	       *       set of reactions being considered. The set of reactions
	       *       is determined by the value of iti.
	       */
	      vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, iti, 0, m_numSpeciesRdc);
	      vcs_deltag(iti, false);
	      /*
	       *       Redefine the starting conditions for noncomponents
	       *       which have yet to be processed in the main loop
	       */
	      for (ll = kspec+1; ll < m_numSpeciesRdc; ++ll) {
		m_feSpecies_old[ll] = m_feSpecies_curr[ll];
	      }
	      for (ll = irxn+1; ll < m_numRxnRdc; ++ll) {
		m_deltaGRxn_old[ll] = m_deltaGRxn_new[ll];
	      }
#ifdef DEBUG_MODE
	      if (vcs_debug_print_lvl >= 2) {
		if (m_rxnStatus[irxn] >= 0) {
		  plogf("   --- SS species changed to zeroedss: ");
		  plogf("%-12s", m_speciesName[kspec].c_str());
		  plogendl();
		}
	      }
#endif
	      m_rxnStatus[irxn] = VCS_SPECIES_ZEROEDSS;
	      ++m_numRxnMinorZeroed;
	      im = (m_numRxnMinorZeroed == m_numRxnRdc);
	      if (im && iti != 0) {
		goto L_EQUILIB_CHECK;
	      }
	      m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec];
	      m_deltaMolNumSpecies[kspec] = 0.0;
	      dx = 0.0;
	    }
	  }
	}
	/*********************************************************************/
	/*** LINE SEARCH ALGORITHM FOR MAJOR SPECIES IN NON-IDEAL PHASES *****/
	/*********************************************************************/
	/*
	 * Skip the line search if we are birthing a species
	 */
	if (dx != 0.0 && (m_molNumSpecies_old[kspec] > 0.0) &&
	    (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE)) {
	  double dx_old = dx;
#ifdef DEBUG_MODE
	  dx = vcs_line_search(irxn, dx_old, ANOTE);
#else
	  dx = vcs_line_search(irxn, dx_old);
#endif
	}
	m_deltaMolNumSpecies[kspec] = dx;

      } /* End of Loop on ic[irxn] -> the type of species */

      /***********************************************************************/
      /****** CALCULATE KMOLE NUMBER CHANGE FOR THE COMPONENT BASIS **********/
      /***********************************************************************/
      if (dx != 0.0 && (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE)) {
	/*
	 *         Change the amount of the component compounds according 
	 *         to the reaction delta that we just computed. 
	 *         This should keep the amount of material constant. 
	 */
#ifdef DEBUG_MODE
	if (fabs(m_deltaMolNumSpecies[kspec] -dx) > 1.0E-14*(fabs(m_deltaMolNumSpecies[kspec]) + fabs(dx) + 1.0E-32)) {
	  plogf(" ds[kspec] = %20.16g dx = %20.16g , kspec = %d\n", m_deltaMolNumSpecies[kspec], dx, kspec);
	  plogf("we have a problem!");	
	  plogendl();
	  exit(-1);
	}
#endif
	for (k = 0; k < m_numComponents; ++k) {
	  m_deltaMolNumSpecies[k] += sc_irxn[k] * dx;
	}
	/*
	 *         Calculate the tentative change in the total number of 
	 *         moles in all of the phases 
	 */
	dnPhase_irxn = m_deltaMolNumPhase[irxn];
	for (iph = 0; iph < m_numPhases; iph++) {
	  m_deltaPhaseMoles[iph] += dx * dnPhase_irxn[iph];
	}
      }

#ifdef DEBUG_MODE
      checkDelta1(VCS_DATA_PTR(m_deltaMolNumSpecies), VCS_DATA_PTR(m_deltaPhaseMoles), kspec+1);
#endif
      /*
       *          Branch point for returning -
       */
#ifndef DEBUG_MODE
    L_MAIN_LOOP_END: ;
#endif
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + m_deltaMolNumSpecies[kspec];
	plogf("   --- "); plogf("%-12.12s", m_speciesName[kspec].c_str());
	plogf("%3d%11.4E%11.4E%11.4E | %s", 
	      m_rxnStatus[irxn], m_molNumSpecies_old[kspec], m_molNumSpecies_new[kspec],
	      m_deltaMolNumSpecies[kspec], ANOTE);
	plogendl();
      }
    L_MAIN_LOOP_END_NO_PRINT: ;
#endif
      
    }  /**************** END OF MAIN LOOP OVER FORMATION REACTIONS ************/

#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      for (k = 0; k < m_numComponents; k++) {
	plogf("   --- ");  plogf("%-12.12s", m_speciesName[k].c_str());
	plogf("  c%11.4E%11.4E%11.4E |\n",
	      m_molNumSpecies_old[k], m_molNumSpecies_old[k]+m_deltaMolNumSpecies[k], m_deltaMolNumSpecies[k]);
      }
      plogf("   "); vcs_print_line("-", 80);
      plogf("   --- Finished Main Loop");
      plogendl();
    }
#endif

    /*************************************************************************/
    /*********** LIMIT REDUCTION OF BASIS SPECIES TO 99% *********************/
    /*************************************************************************/
    /*
     *        We have a tentative M_DELTAMOLNUMSPECIES(L=1,MR). Now apply other criteria 
     *        to limit it's magnitude. 
     */
    par = 0.5;
    for (k = 0; k < m_numComponents; ++k) {
      if (m_molNumSpecies_old[k] > 0.0) {
	xx = -m_deltaMolNumSpecies[k] / m_molNumSpecies_old[k];
	if (par < xx) {
	  par = xx;
#ifdef DEBUG_MODE
	  ll = k;	    
#endif
	}
      } else {
	if (m_deltaMolNumSpecies[k] < 0.0) {
	  /*
	   * If we are here, we then do a step which violates element
	   * conservation.
	   */
	  iph = m_phaseID[k];
	  m_deltaPhaseMoles[iph] -= m_deltaMolNumSpecies[k];
	  m_deltaMolNumSpecies[k] = 0.0;
	}
      }
    }
    par = 1.0 / par;
    if (par <= 1.01 && par > 0.0) {
      /* Reduce the size of the step by the multiplicative factor, par */
      par *= 0.99;
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- Reduction in step size due to component ");
	plogf("%s", m_speciesName[ll].c_str());
	plogf(" going negative = %11.3E", par); 
	plogendl();
      }
#endif
      for (i = 0; i < m_numSpeciesTot; ++i) {
	m_deltaMolNumSpecies[i] *= par;
      }
      for (iph = 0; iph < m_numPhases; iph++) {
	m_deltaPhaseMoles[iph] *= par;	 
      }
    } else {
      par = 1.0;
    }
#ifdef DEBUG_MODE
    checkDelta1(VCS_DATA_PTR(m_deltaMolNumSpecies), 
		VCS_DATA_PTR(m_deltaPhaseMoles), m_numSpeciesTot);
#endif
   
    /*
     *      Now adjust the wt[kspec]'s so that the reflect the decrease in 
     *      the overall length of m_deltaMolNumSpecies[kspec] just calculated. At the end 
     *      of this section wt[], m_deltaMolNumSpecies[], tPhMoles, and tPhMoles1 should all be 
     *      consistent with a new estimate of the state of the system. 
     */
    for (kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
      m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + m_deltaMolNumSpecies[kspec];
      if (m_molNumSpecies_new[kspec] < 0.0 && (m_speciesUnknownType[kspec] 
					       != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE)) {
	plogf("vcs_solve_TP: ERROR on step change wt[%d:%s]: %g < 0.0",
	      kspec, m_speciesName[kspec].c_str(), m_molNumSpecies_new[kspec]);
	plogendl();
	exit(-1);
      }
    }
   
    /*
     *        Calculate the tentative total mole numbers for each phase
     */
    for (iph = 0; iph < m_numPhases; iph++) {
      m_tPhaseMoles_new[iph] = m_tPhaseMoles_old[iph] + m_deltaPhaseMoles[iph];
    }
    /*
     *         Calculate the new chemical potentials using the tentative 
     *         solution values. We only calculate a subset of these, because 
     *         we have only updated a subset of the W(). 
     */
    vcs_updateVP(1);
    //vcs_dfe(VCS_DATA_PTR(wt), VCS_STATECALC_NEW, iti, 0, m_numSpeciesTot);
    vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_new), VCS_STATECALC_NEW, 0, 0, m_numSpeciesTot);
    /*
     *         Evaluate DeltaG for all components if ITI=0, and for 
     *         major components only if ITI NE 0 
     */
    // if (iti == 0) vcs_deltag(0, false);
    //else          vcs_deltag(-1, false);
    vcs_deltag(0, false);
 
    // Actually always need to calculate this 
    // or else nonprintouts get different results and sometimes
    // fail in the line search algorithm -> Why is this?
    //vcs_dfe(VCS_DATA_PTR(wt), VCS_STATECALC_NEW, 1, 0, m_numSpeciesRdc);
    //if (iti != 0) {
    // vcs_deltag(1, false);
    //}


    /* *************************************************************** */
    /* **** CONVERGENCE FORCER SECTION ******************************* */
    /* *************************************************************** */
    if (printDetails) {
      plogf("   --- Total Old       Dimensionless Gibbs Free Energy = %20.13E\n", 
	    vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_feSpecies_old), 
			    VCS_DATA_PTR(m_tPhaseMoles_old)));
      plogf("   --- Total tentative Dimensionless Gibbs Free Energy = %20.13E", 
	    vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_feSpecies_curr), 
			    VCS_DATA_PTR(m_tPhaseMoles_new)));
      plogendl();
    }

    forced = globStepDamp(iti);

    /*
     *       Print out the changes to the solution that FORCER produced 
     */
    if (printDetails && forced) {
      
      plogf(" -----------------------------------------------------\n");
      plogf("   --- FORCER SUBROUTINE changed the solution:\n");
      plogf("   --- SPECIES Status INIT MOLES TENT_MOLES");
      plogf("  FINAL KMOLES  INIT_DEL_G/RT  TENT_DEL_G/RT  FINAL_DELTA_G/RT\n");
      for (i = 0; i < m_numComponents; ++i) {
	plogf("  --- %-12.12s", m_speciesName[i].c_str());
	plogf("    %14.6E %14.6E %14.6E\n",  m_molNumSpecies_old[i], m_molNumSpecies_old[i] + m_deltaMolNumSpecies[i], m_molNumSpecies_new[i]);
      }
      for (kspec = m_numComponents; kspec < m_numSpeciesRdc; ++kspec) {
	irxn = kspec - m_numComponents;
	plogf("  --- %-12.12s", m_speciesName[kspec].c_str());
	plogf(" %2d %14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n", m_rxnStatus[irxn],
	      m_molNumSpecies_old[kspec], m_molNumSpecies_old[kspec]+m_deltaMolNumSpecies[kspec], 
	      m_molNumSpecies_new[kspec], m_deltaGRxn_old[irxn], 
	      m_deltaGRxn_tmp[irxn], m_deltaGRxn_new[irxn]);
      }
      print_space(26); 
      plogf("Norms of Delta G():%14.6E%14.6E\n",
	    l2normdg(VCS_DATA_PTR(m_deltaGRxn_old)),
	    l2normdg(VCS_DATA_PTR(m_deltaGRxn_new)));
      plogf("   Total kmoles of gas    = %15.7E\n", m_tPhaseMoles_old[0]);
      if ((m_numPhases > 1) && (! (m_VolPhaseList[1])->SingleSpecies)) { 
	plogf("   Total kmoles of liquid = %15.7E\n", m_tPhaseMoles_old[1]); 
      } else {
	plogf("   Total kmoles of liquid = %15.7E\n", 0.0);
      }
      plogf("   Total New Dimensionless Gibbs Free Energy = %20.13E\n", 
	    vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_feSpecies_curr),
			    VCS_DATA_PTR(m_tPhaseMoles_new)));
      plogf(" -----------------------------------------------------");
      plogendl();
    }

    /* *************************************************************** */
    /* **** ITERATION SUMMARY PRINTOUT SECTION *********************** */
    /* *************************************************************** */

    if (printDetails) {
      plogf("   "); vcs_print_line("-", 103); 
      plogf("   --- Summary of the Update ");
      if (iti == 0) {
	plogf(" (all species):\n");
      } else {
	plogf(" (only major species):\n");
      }
      plogf("   ---      Species Status Initial_KMoles Final_KMoles Initial_Mu/RT");
      plogf("     Mu/RT     Init_Del_G/RT   Delta_G/RT\n");
      for (i = 0; i < m_numComponents; ++i) {
	plogf("   ---   %-12.12s", m_speciesName[i].c_str()); plogf("    "); 
	plogf("%14.6E%14.6E%14.6E%14.6E\n", m_molNumSpecies_old[i],
	      m_molNumSpecies_new[i], m_feSpecies_old[i], m_feSpecies_curr[i]);
      }
      for (i = m_numComponents; i < m_numSpeciesRdc; ++i) {
	l1 = i - m_numComponents;
	plogf("   ---   %-12.12s", m_speciesName[i].c_str());
	plogf(" %2d %14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
	      m_rxnStatus[l1], m_molNumSpecies_old[i],
	      m_molNumSpecies_new[i], m_feSpecies_old[i], m_feSpecies_curr[i],
	      m_deltaGRxn_old[l1], m_deltaGRxn_new[l1]);
      }
      for (kspec = m_numSpeciesRdc; kspec < m_numSpeciesTot; ++kspec) {
	l1 = kspec - m_numComponents;
	plogf("   ---   %-12.12s", m_speciesName[kspec].c_str());
	plogf(" %2d %14.6E%14.6E%14.6E%14.6E%14.6E%14.6E\n",
	      m_rxnStatus[l1], m_molNumSpecies_old[kspec],
	      m_molNumSpecies_new[kspec], m_feSpecies_old[kspec], m_feSpecies_curr[kspec],
	      m_deltaGRxn_old[l1], m_deltaGRxn_new[l1]);
      }
      plogf("   ---"); print_space(56);
      plogf("Norms of Delta G():%14.6E%14.6E",
	    l2normdg(VCS_DATA_PTR(m_deltaGRxn_old)),
	    l2normdg(VCS_DATA_PTR(m_deltaGRxn_new)));
      plogendl();
      
      plogf("   ---           Phase_Name    KMoles(after update)\n");
      plogf("   ---   "); vcs_print_line("-", 50);
      for (iph = 0; iph < m_numPhases; iph++) {
	Vphase = m_VolPhaseList[iph];
	plogf("   ---   %18s = %15.7E\n", Vphase->PhaseName.c_str(), m_tPhaseMoles_new[iph]);
      }
      plogf("   "); vcs_print_line("-", 103);
      plogf("   --- Total Old Dimensionless Gibbs Free Energy = %20.13E\n", 
	    vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_feSpecies_old), 
			    VCS_DATA_PTR(m_tPhaseMoles_old)));
      plogf("   --- Total New Dimensionless Gibbs Free Energy = %20.13E", 
	    vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_new), VCS_DATA_PTR(m_feSpecies_curr), 
			    VCS_DATA_PTR(m_tPhaseMoles_new)));
      plogendl();
      if (m_VCount->Its > 550) {
	plogf("   --- Troublesome solve"); 
	plogendl();
      }

    }

    /*************************************************************************/
    /******************* RESET VALUES AT END OF ITERATION ********************/
    /******************* UPDATE MOLE NUMBERS *********************************/
    /*************************************************************************/
    /*
     *       If the solution wasn't changed in the forcer routine,
     *       then copy the tentative mole numbers and Phase moles
     *       into the actual mole numbers and phase moles.
     *       We will consider this current step to be completed.
     *
     *   Accept the step. -> the tentative solution now becomes 
     *                the real solution. If FORCED is true, then 
     *                we have already done this inside the FORCED 
     *                loop. 
     */
    vcs_dcopy(VCS_DATA_PTR(m_tPhaseMoles_old), VCS_DATA_PTR(m_tPhaseMoles_new), m_numPhases);
    vcs_dcopy(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_molNumSpecies_new), m_numSpeciesRdc);
    vcs_dcopy(VCS_DATA_PTR(m_deltaGRxn_old), VCS_DATA_PTR(m_deltaGRxn_new), m_numRxnRdc);
    vcs_dcopy(VCS_DATA_PTR(m_feSpecies_old), VCS_DATA_PTR(m_feSpecies_curr), m_numSpeciesRdc);
      
    vcs_updateVP(0);
    /*
     *       Increment the iteration counters
     */
    ++(m_VCount->Its);
    ++it1;
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Increment counter increased, step is accepted: %4d",
	    m_VCount->Its);
      plogendl();
    }
#endif
    /*************************************************************************/
    /******************* HANDLE DELETION OF MULTISPECIES PHASES **************/
    /*************************************************************************/
    /*
     *   We delete multiphases, when the total moles in the multiphase
     *   is reduced below a relative threshold.
     *   Set microscopic multispecies phases with total relative 
     *   number of moles less than  VCS_DELETE_PHASE_CUTOFF to
     *   absolute zero.
     */
    justDeletedMultiPhase = FALSE;
    for (iph = 0; iph < m_numPhases; iph++) {
      Vphase = m_VolPhaseList[iph];
      if (!(Vphase->SingleSpecies)) {
	if (m_tPhaseMoles_old[iph] != 0.0 &&
	    m_tPhaseMoles_old[iph]/m_totalMolNum <= VCS_DELETE_PHASE_CUTOFF) {
	  soldel = 1;
	  for (kspec = 0; kspec < m_numSpeciesRdc; kspec++) {
	    if (m_phaseID[kspec] == iph && m_molNumSpecies_old[kspec] > 0.0) {
	      irxn = kspec - m_numComponents;
	      if (kspec < m_numComponents) {
		if (m_molNumSpecies_old[kspec] > VCS_DELETE_SPECIES_CUTOFF) {
		  soldel = 0;
		  break;
		}
	      } else {
		for (k = 0; k < m_numComponents; k++) {
		  if (m_stoichCoeffRxnMatrix[irxn][k] != 0.0) {
		    if (m_molNumSpecies_old[kspec]/m_molNumSpecies_old[k] > VCS_DELETE_PHASE_CUTOFF) {
		      soldel = 0;
		      break;
		    }
		  }
		}
	      }
	    }
	  }
	  if (soldel) {
#ifdef DEBUG_MODE
	    if (vcs_debug_print_lvl >= 1) {
	      plogf("   --- Setting microscopic phase %d to zero", iph);
	      plogendl();
	    }
#endif
	    justDeletedMultiPhase = TRUE;
	    delete_multiphase(iph);
	  }
	} 
      }
    }
    /*
     *       If we have deleted a multispecies phase because the
     *       equilibrium moles decreased, then we will update all
     *       the component basis calculation, and therefore all
     *       of the thermo functions just to be safe.
     */
    if (justDeletedMultiPhase) {
      justDeletedMultiPhase = FALSE;
      retn = vcs_basopt(FALSE, VCS_DATA_PTR(aw), VCS_DATA_PTR(sa),
			VCS_DATA_PTR(sm), VCS_DATA_PTR(ss), test, 
			&usedZeroedSpecies);
      if (retn != VCS_SUCCESS) return retn;
      vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
      vcs_deltag(0, true);
      uptodate_minors = TRUE;
      if (conv) {
	/*
	 *      HKM -> I don't understand why the code would just give
	 *             up here in some cases.
	 *             This should probably be taken out
	 */	    
	plogf(" DELETION OF MULTISPECIES PHASE. ");
	plogf("Convergence to number of positive n(i) less than C.\n");
	plogf("Check results to follow carefully.   \n");
	plogendl();
	goto L_RETURN_BLOCK;
      }
    }
    /*************************************************************************/
    /***************** CHECK FOR ELEMENT ABUNDANCE****************************/
    /*************************************************************************/
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Normal element abundance check");
    }
#endif
    vcs_elab();
    if (! vcs_elabcheck(0)) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf(" - failed -> redoing element abundances.");
	plogendl();
      }
#endif
      vcs_elcorr(VCS_DATA_PTR(sm), VCS_DATA_PTR(wx));
      vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
      vcs_deltag(0, true);
      uptodate_minors = TRUE;
    }
#ifdef DEBUG_MODE
    else {
      if (vcs_debug_print_lvl >= 2) {
	plogf(" - passed");
	plogendl();
      }
    }
#endif
    /*************************************************************************/
    /***************** CHECK FOR OPTIMUM BASIS *******************************/
    /*************************************************************************/
    /*
     *   HKM -> We first evaluate whether the components species are
     *          ordered according to their mole numbers. If they are,
     *          then we can essential do an order(NR) operation instead
     *          of an order(NR*NC) operation to determine whether
     *          a new basis is needed.
     *
     *   HKM -> This section used to be branched to initially if
     *          there was a machine estimate. I took it out to simplify
     *          the code logic.
     */
    dofast = (m_numComponents != 1);
    for (i = 1; i < m_numComponents; ++i) {
      if ((m_molNumSpecies_old[i - 1] * m_spSize[i-1]) < (m_molNumSpecies_old[i] * m_spSize[i])) {
	dofast = FALSE;
	break;
      }
    }
    dofast = false;
    if (dofast) {
      for (i = 0; i < m_numRxnRdc; ++i) {
	l = m_indexRxnToSpecies[i];
	for (j = m_numComponents - 1; j >= 0; j--) {
	  bool doSwap = false;
	  if (m_SSPhase[j]) {
	    doSwap = (m_molNumSpecies_old[l] * m_spSize[l]) > (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
	    if (!m_SSPhase[i]) {
	      if (doSwap) {
		doSwap = (m_molNumSpecies_old[l]) > (m_molNumSpecies_old[j] * 1.01);
	      }
	    }
	  } else {
	    if (m_SSPhase[i]) {
	      doSwap = (m_molNumSpecies_old[l] * m_spSize[l]) > (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
	      if (!doSwap) {
		doSwap = (m_molNumSpecies_old[l]) > (m_molNumSpecies_old[j] * 1.01);
	      }
	    } else {
	      doSwap = (m_molNumSpecies_old[l] * m_spSize[l]) > (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
	    }
	  }
	  if (doSwap) {
	    if (m_stoichCoeffRxnMatrix[i][j] != 0.0) {
#ifdef DEBUG_MODE
	      if (vcs_debug_print_lvl >= 2) {
		plogf("   --- Get a new basis because %s", m_speciesName[l].c_str());
		plogf(" is better than comp %s", m_speciesName[j].c_str());
		plogf(" and share nonzero stoic: %-9.1f", 
		      m_stoichCoeffRxnMatrix[i][j]);
		plogendl();
	      }
#endif		     
	      goto L_COMPONENT_CALC;
	    }
	  } else {
	    break;
	  }
#ifdef DEBUG_NOT
	  if (m_rxnStatus[i] == VCS_SPECIES_ZEROEDMS) {
	    if (m_molNumSpecies_old[j] == 0.0) {
	      if (m_stoichCoeffRxnMatrix[i][j] != 0.0) {
		if (dg[i] < 0.0) {
#ifdef DEBUG_MODE
		  if (vcs_debug_print_lvl >= 2) {
		    plogf("   --- Get a new basis because %s", m_speciesName[l].c_str());
		    plogf(" has dg < 0.0 and comp %s has zero mole num", m_speciesName[j].c_str());
		    plogf(" and share nonzero stoic: %-9.1f", 
			  m_stoichCoeffRxnMatrix[i][j]);
		    plogendl();
		  }
#endif		     
		  goto L_COMPONENT_CALC;
		}
	      }
	    }
	  }
#endif
	}
      }
    } else {
      for (i = 0; i < m_numRxnRdc; ++i) {
	l = m_indexRxnToSpecies[i];
	for (j = 0; j < m_numComponents; ++j) {
	  bool doSwap = false;
	  if (m_SSPhase[j]) {
	    doSwap = (m_molNumSpecies_old[l] * m_spSize[l]) > (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
	    if (!m_SSPhase[l]) {
	      if (doSwap) {
		doSwap = (m_molNumSpecies_old[l]) > (m_molNumSpecies_old[j] * 1.01);
	      }
	    }
	  } else {
	    if (m_SSPhase[l]) {
	      doSwap = (m_molNumSpecies_old[l] * m_spSize[l]) > (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
	      if (!doSwap) {
		doSwap = (m_molNumSpecies_old[l]) > (m_molNumSpecies_old[j] * 1.01);
	      }
	    } else {
	      doSwap = (m_molNumSpecies_old[l] * m_spSize[l]) > (m_molNumSpecies_old[j] * m_spSize[j] * 1.01);
	    }
	  }
	  if (doSwap) {
	    if (m_stoichCoeffRxnMatrix[i][j] != 0.0) {
#ifdef DEBUG_MODE
	      if (vcs_debug_print_lvl >= 2) {
		plogf("   --- Get a new basis because ");
		plogf("%s", m_speciesName[l].c_str());
		plogf(" is better than comp ");
		plogf("%s", m_speciesName[j].c_str());
		plogf(" and share nonzero stoic: %-9.1f", 
		      m_stoichCoeffRxnMatrix[i][j]);
		plogendl();
	      }
#endif			     
	      goto L_COMPONENT_CALC;
	    }
	  }
#ifdef DEBUG_NOT
	  if (m_rxnStatus[i] == VCS_SPECIES_ZEROEDMS) {
	    if (m_molNumSpecies_old[j] == 0.0) {
	      if (m_stoichCoeffRxnMatrix[i][j] != 0.0) {
		if (dg[i] < 0.0) {
#ifdef DEBUG_MODE
		  if (vcs_debug_print_lvl >= 2) {
		    plogf("   --- Get a new basis because %s", m_speciesName[l].c_str());
		    plogf(" has dg < 0.0 and comp %s has zero mole num", m_speciesName[j].c_str());
		    plogf(" and share nonzero stoic: %-9.1f", 
			  m_stoichCoeffRxnMatrix[i][j]);
		    plogendl();
		  }
#endif		     
		  goto L_COMPONENT_CALC;
		}
	      }
	    }
	  }
#endif
	}
      }
    }
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Check for an optimum basis passed");
      plogendl();
    }
#endif
    /*************************************************************************/
    /********************** RE-EVALUATE MAJOR-MINOR VECTOR IF NECESSARY ******/
    /*************************************************************************/
    /*
     *     Skip this section if we haven't done a full calculation. 
     *     Go right to the check equilibrium section 
     */
    if (iti == 0) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- Reevaluate major-minor status of noncomponents:\n");
      }
#endif
      m_numRxnMinorZeroed = 0;
      for (irxn = 0; irxn < m_numRxnRdc; irxn++) {
	kspec = m_indexRxnToSpecies[irxn];
	 
	int speciesType = vcs_species_type(kspec);
	if (speciesType < VCS_SPECIES_MINOR) {
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {
	    if (m_rxnStatus[irxn]  >= VCS_SPECIES_MINOR) {
	      plogf("   ---    major/minor species is now zeroed out: %s\n", 
		    m_speciesName[kspec].c_str());
	    }
	  }
#endif	    
	  ++m_numRxnMinorZeroed;	    
	} else if (speciesType == VCS_SPECIES_MINOR) {
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {
	    if (m_rxnStatus[irxn] != VCS_SPECIES_MINOR) {
	      if (m_rxnStatus[irxn] == VCS_SPECIES_MAJOR) {
		plogf("   ---   Noncomponent turned from major to minor: ");
	      } else if (m_rxnStatus[irxn] == VCS_SPECIES_COMPONENT) {
		plogf("   ---   Component turned into a minor species: ");
	      } else {
		plogf("   ---   Zeroed Species turned into a "
		      "minor species: ");       
	      }
	      plogf("%s\n", m_speciesName[kspec].c_str());
	    }
	  }
#endif
	  ++m_numRxnMinorZeroed;
	} else if (speciesType == VCS_SPECIES_MAJOR) {
	  if (m_rxnStatus[irxn] != VCS_SPECIES_MAJOR) {
#ifdef DEBUG_MODE
	    if (vcs_debug_print_lvl >= 2) {
	      if (m_rxnStatus[irxn] ==  VCS_SPECIES_MINOR) {
		plogf("   ---   Noncomponent turned from minor to major: ");
	      } else if (m_rxnStatus[irxn] == VCS_SPECIES_COMPONENT) {
		plogf("   ---   Component turned into a major: ");	       
	      } else {
		plogf("   ---   Noncomponent turned from zeroed to major: ");
	      }
	      plogf("%s\n", m_speciesName[kspec].c_str());
	    }
#endif
	    m_rxnStatus[irxn] = VCS_SPECIES_MAJOR;
	    /*
	     *   For this special case, we must reevaluate thermo functions
	     */
	    if (iti != 0) {
	      vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 0, kspec, kspec+1);
	      vcs_deltag(0, false);
	    }
	  }
	}
	m_rxnStatus[irxn] = speciesType;
      }
      /*
       *         This logical variable indicates whether all current 
       *         non-component species are minor or nonexistent 
       */
      im = (m_numRxnMinorZeroed == m_numRxnRdc);
    }
    /*************************************************************************/
    /***************** EQUILIBRIUM CHECK FOR MAJOR SPECIES *******************/
    /*************************************************************************/
  L_EQUILIB_CHECK: ;
    if (! im) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- Equilibrium check for major species: ");
      }
#endif
      for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
	if (m_rxnStatus[irxn] == VCS_SPECIES_MAJOR && (fabs(m_deltaGRxn_new[irxn]) > m_tolmaj)) {
	  if (m_VCount->Its >= maxit) {
	    solveFail = -1;
	    /* 
	     *         Clean up and exit code even though we haven't 
	     *         converged. -> we have run out of iterations! 
	     */
	    goto L_RETURN_BLOCK;
	  } else {
#ifdef DEBUG_MODE
	    if (vcs_debug_print_lvl >= 2) {
	      plogf("%s failed\n", m_speciesName[m_indexRxnToSpecies[irxn]].c_str());
	    }
#endif
	    /*
	     *  Set MajorSpeciesHaveConverged to false to indicate that
	     *  convergence amongst 
	     *  major species has not been achieved
	     */
	    MajorSpeciesHaveConverged = false;
	    /*
	     *   Go back and do another iteration with variable ITI 
	     */
	    goto L_MAINLOOP_MM4_SPECIES;
	  }
	}
      }
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf(" MAJOR SPECIES CONVERGENCE achieved");
	plogendl();
      }
#endif
    }
#ifdef DEBUG_MODE
    else {
      if (vcs_debug_print_lvl >= 2) {
	plogf(" MAJOR SPECIES CONVERGENCE achieved "
	      "(because there are no major species)");
	plogendl();
      }
    }
#endif
    /*
     *  Set MajorSpeciesHaveConverged to true to indicate
     * that convergence amongst major species has been achieved
     */
    MajorSpeciesHaveConverged = true;
    /*************************************************************************/
    /*************** EQUILIBRIUM CHECK FOR MINOR SPECIES *********************/
    /*************************************************************************/
    if (m_numRxnMinorZeroed != 0) {
      /*
       *       Calculate the chemical potential and reaction DeltaG 
       *       for minor species, if needed.
       */
      if (iti != 0) {
	vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
	vcs_deltag(1, false);
	uptodate_minors = TRUE;
      }
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- Equilibrium check for minor species: ");
      }
#endif
      for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
	if (m_rxnStatus[irxn] == VCS_SPECIES_MINOR && (fabs(m_deltaGRxn_new[irxn]) > m_tolmin)) {
	  if (m_VCount->Its >= maxit) {
	    solveFail = -1;
	    /*
	     *       Clean up and exit code. -> Even though we have not
	     *        converged, we have run out of iterations !        
	     */
	    goto L_RETURN_BLOCK;
	  }
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {
	    plogf("%s failed\n", m_speciesName[m_indexRxnToSpecies[irxn]].c_str());
	  }
#endif
	  /*
	   *  Set iti to zero to force a full calculation, and go back
	   *  to the main loop to do another iteration.
	   */
	  iti = 0;
	  goto L_MAINLOOP_ALL_SPECIES;
	}
      }
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf(" CONVERGENCE achieved\n");
      }
#endif
    }
    /*************************************************************************/
    /*********************** FINAL ELEMENTAL ABUNDANCE CHECK *****************/
    /*************************************************************************/
    /*
     *    Recalculate the element abundance vector again
     */
    vcs_updateVP(0);
    vcs_elab();
   
    /*        LEC is only true when we are near the end game */
    if (lec) {
      if (!giveUpOnElemAbund) {
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  plogf("   --- Check the Full Element Abundances: ");
	}
#endif
	/*
	 *  Final element abundance check:
	 *        If we fail then we need to go back and correct
	 *        the element abundances, and then go do a major step
	 */
	if (! vcs_elabcheck(1) ) {
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {
	    if (! vcs_elabcheck(0)) {
	      plogf(" failed\n");
	    } else {
	      plogf(" passed for NC but failed for NE: RANGE ERROR\n");
	    }
	  }
#endif
	  // delete?
	  goto L_ELEM_ABUND_CHECK;
	}
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  plogf(" passed\n");
	}
#endif
      }
      /*
       *   If we have deleted a species then we need to recheck the
       *   the deleted species, before exiting
       */
      if (m_numSpeciesRdc != m_numSpeciesTot) {
	goto L_RECHECK_DELETED;
      }
      /* - Final checks are passed -> go check out */
      goto L_RETURN_BLOCK;
    }
    lec = TRUE;
    /* *************************************************** */
    /* **** CORRECT ELEMENTAL ABUNDANCES ***************** */
    /* *************************************************** */
  L_ELEM_ABUND_CHECK: ;
    /*
     *  HKM - Put in an element abundance check. The element abundances  
     *        were being corrected even if they were perfectly OK to 
     *        start with. This is actually an expensive operation, so 
     *        I took it out. Also vcs_dfe() doesn't need to be called if 
     *        no changes were made. 
     */
    rangeErrorFound = 0; 
    if (! vcs_elabcheck(1)) {
      int ncBefore = vcs_elabcheck(0);
      vcs_elcorr(VCS_DATA_PTR(sm), VCS_DATA_PTR(wx));
      int ncAfter = vcs_elabcheck(0);
      int neAfter = vcs_elabcheck(1);
      /*
       *      Go back to evaluate the total moles of gas and liquid. 
       */
      vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 0, 0, m_numSpeciesRdc);
      vcs_deltag(0, false);
      /*
       * 
       */
      if (!ncBefore) {
	if (ncAfter) {
	  /*
	   * We have breathed new life into the old problem. Now the
	   * element abundances up to NC agree. Go back and
	   * restart the main loop calculation, resetting the
	   * end conditions.
	   */
	  lec = FALSE;
	  iti = 0;
	  goto L_MAINLOOP_ALL_SPECIES;
	} else {
	  /*
	   * We are still hosed 
	   */
	  if (finalElemAbundAttempts >= 3) {
	    giveUpOnElemAbund = true;
	    goto L_EQUILIB_CHECK;
	  } else {
	    finalElemAbundAttempts++;
	    lec = FALSE;
	    iti = 0;
	    goto L_MAINLOOP_ALL_SPECIES;
	  }
	}
      } else {
	if (ncAfter) {
	  if (neAfter) {
	    /*
	     * Recovery of end element abundances
	     * -> go do equilibrium check again and then
	     *    check out.
	     */
	    goto L_EQUILIB_CHECK;
	  } else {
	    /*
	     * Probably an unrecoverable range error
	     */
#ifdef DEBUG_MODE
	    if (vcs_debug_print_lvl >= 2) {
	      plogf(" ---  vcs_solve_tp: RANGE SPACE ERROR ENCOUNTERED\n");
	      plogf(" ---  vcs_solve_tp: - Giving up on NE Element Abundance satisfaction \n");
	      plogf(" ---  vcs_solve_tp: - However, NC Element Abundance criteria is satisfied \n");
	      plogf(" ---  vcs_solve_tp: - Returning the calculated equilibrium condition ");
	      plogendl();
	    }
#endif
	    rangeErrorFound = 1;
	    giveUpOnElemAbund = true;
	    goto L_EQUILIB_CHECK;
	  }
	}
      }    
    }
    // Calculate delta g's
    vcs_deltag(0, false);
    // Go back to equilibrium check as a prep to eventually checking out
    goto L_EQUILIB_CHECK;

    /* *************************************************** */
    /* **** RECHECK DELETED SPECIES ********************** */
    /* *************************************************** */
    /*
     *         We are here for two reasons. One is if we have
     *         achieved convergence, but some species have been eliminated
     *         from the problem because they were in multispecies phases
     *         and their mole fractions drifted less than 
     *         VCS_DELETE_SPECIES_CUTOFF .
     *         The other reason why we are here is because all of the
     *         non-component species in the problem have been eliminated
     *         for one reason or another.
     */
  L_RECHECK_DELETED: ;
    npb = recheck_deleted();
    /*
     *        If we haven't found any species that needed adding we are done.
     */
    if (npb <= 0) {
      goto L_RETURN_BLOCK_B;
    }
    /*
     *        If we have found something to add, recalculate everything 
     *        for minor species and go back to do a full iteration
     */
    MajorSpeciesHaveConverged = true;
    vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
    vcs_deltag(0, false);
    iti = 0;
    goto L_MAINLOOP_ALL_SPECIES;
    /*************************************************************************/
    /******************** CLEANUP AND RETURN BLOCK ***************************/
    /*************************************************************************/
  L_RETURN_BLOCK: ;
    
    npb = recheck_deleted();
    /*
     *        If we haven't found any species that needed adding we are done.
     */
    if (npb > 0) {
      /*
       *        If we have found something to add, recalculate everything 
       *        for minor species and go back to do a full iteration
       */
      MajorSpeciesHaveConverged = true;
      vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 1, 0, m_numSpeciesRdc);
      vcs_deltag(0, false);
      iti = 0;
      goto L_MAINLOOP_ALL_SPECIES;
    }

  L_RETURN_BLOCK_B: ;

    /*
     *  Add back deleted species in non-zeroed phases. Estimate their
     *  mole numbers.
     */
    npb = vcs_add_all_deleted();
    if (npb > 0) {
      MajorSpeciesHaveConverged = true;
      iti = 0;
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 1) {
	plogf("  --- add_all_deleted(): some rxns not converged. RETURNING TO LOOP!");
	plogendl();
      }
#endif
      goto L_MAINLOOP_ALL_SPECIES;
    }

    /*
     * Make sure the volume phase objects hold the same state and 
     * information as the vcs object. This also update the Cantera objects
     * with this information.
     */
    vcs_updateVP(0);
    /*
     *    Store the final Delta G values for each non-component species
     *    in the species slot rather than the reaction slot
     */
    kspec = m_numSpeciesTot;
    i = m_numRxnTot;
    for (irxn = 0; irxn < m_numRxnTot; ++irxn) {
      --kspec;
      --i;
      m_deltaGRxn_new[kspec] = m_deltaGRxn_new[i];
    }
    vcs_dzero(VCS_DATA_PTR(m_deltaGRxn_new), m_numComponents);
    /* 
     *       Evaluate the final mole fractions
     *        storring them in wt[]
     */
    vcs_vdzero(m_molNumSpecies_new, m_numSpeciesTot);
    for (kspec = 0; kspec < m_numSpeciesTot; ++kspec) {
      if (m_SSPhase[kspec]) {
	m_molNumSpecies_new[kspec] = 1.0;
      } else {
	iph = m_phaseID[kspec];
	if (m_tPhaseMoles_old[iph] != 0.0) {
	  m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] / m_tPhaseMoles_old[iph];
	} else {
	  /*
	   * For MultiSpecies phases that are zeroed out,
	   * return the mole fraction vector from the VolPhase object.
	   * This contains the mole fraction that would be true if
	   * the phase just pops into existence.
	   */
	  i = m_speciesLocalPhaseIndex[kspec];
	  Vphase = m_VolPhaseList[iph];
	  m_molNumSpecies_new[kspec] = Vphase->molefraction(i);
	}
      }
    }
    // Return an error code if a Range Space Error is thought to have occurred.
    if (rangeErrorFound) {
      solveFail = 1;
    }
    /*
     *           Free temporary storage used in this routine
     *           and increment counters
     */
    /*
     *           Calculate counters
     */
    double tsecond = ticktock.secondsWC();
    m_VCount->Time_vcs_TP = tsecond;
    m_VCount->T_Time_vcs_TP += m_VCount->Time_vcs_TP;
    (m_VCount->T_Calls_vcs_TP)++;
    m_VCount->T_Its += m_VCount->Its;
    m_VCount->T_Basis_Opts += m_VCount->Basis_Opts;
    m_VCount->T_Time_basopt += m_VCount->Time_basopt;
    /*
     *          Return a Flag indicating whether convergence occurred
     */
    return solveFail;
  } /* vcs_solve_TP() **********************************************************/

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  double VCS_SOLVE::minor_alt_calc(int kspec, int irxn, int *do_delete
#ifdef DEBUG_MODE
				   , char *ANOTE  
#endif
				   )
   
    /**************************************************************************
     *
     * minor_alt_calc:
     *
     *        Minor species alternative calculation 
     *       --------------------------------------- 
     * 
     *  This is based upon the following approximation: 
     *    The mole fraction changes due to these reactions don't affect 
     *    the mole numbers of the component species. Therefore the following 
     *    approximation is valid for an ideal solution phase:
     *       0 = M_DELTAGRXN_NEW(I) + log(WT(I)/W(I))
     *
     *          W(i) = Old mole number of species i in the phase
     *          WT(i) = Trial new mole number of species i in the pahse
     * 
     *       (M_DELTAGRXN_NEW contains the contribution from
     *         FF(I) + log(ActCoeff[i] * W(I)/Total_Moles) ) 
     *    Thus, 
     *        WT(I) = W(I) EXP(-M_DELTAGRXN_NEW(I)) 
     * 
     *    Most of this section is mainly restricting the update to reasonable 
     *    values.
     *
     *
     *    Note: This routine was generalized to incorporate
     *          nonideal phases.
     *
     *    Input:
     *     ------
     *     kspec, irxn = the current species and corresponding formation
     *                   reaction number.
     *    Output:
     *    ---------
     *     return value: dx = the change in mole number
     *     do_delete:  BOOLEAN which if true on return, then we branch 
     *                      to the section that deletes a species from the
     *                      current set of active species.
     *************************************************************************/
  {
    double dx;
    double  w_kspec  = m_molNumSpecies_old[kspec];
    double *wt_kspec = VCS_DATA_PTR(m_molNumSpecies_new) + kspec;
    double  wTrial;
    double *ds_kspec = VCS_DATA_PTR(m_deltaMolNumSpecies) + kspec;
    double  dg_irxn  = m_deltaGRxn_new[irxn];
    int iphase = m_phaseID[kspec];
    vcs_VolPhase *Vphase = m_VolPhaseList[iphase];
    *do_delete = FALSE;
    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
      if (w_kspec <= 0.0) {
	w_kspec = VCS_DELETE_MINORSPECIES_CUTOFF;
      }
      if (dg_irxn < -20.) {
	dg_irxn = -20.;
      }
#ifdef DEBUG_MODE
      sprintf(ANOTE,"minor species alternative calc");
#endif
      if (dg_irxn >= 82.0) {
	(*wt_kspec) = w_kspec * 1.0e-6;
	if (w_kspec < VCS_DELETE_MINORSPECIES_CUTOFF) {
	  goto L_ZERO_SPECIES;
	}
      } else {
	if (fabs(dg_irxn) <= m_tolmin2) {
	  (*wt_kspec) = w_kspec;
	  (*ds_kspec) = 0.0;
	  return 0.0;
	}
	//  c = log(ActCoeff[kspec] * w_kspec) - dg_irxn;
    


      }

      if (dg_irxn > 10.0) {
	(*wt_kspec) = w_kspec * 1.0e-5;
	if (w_kspec <  VCS_DELETE_MINORSPECIES_CUTOFF) {
	  goto L_ZERO_SPECIES;
	}
      } else {
	double ac0 = m_actCoeffSpecies_new[kspec];
	double ac  = ac0;
	double w0 = w_kspec;
	double dd = exp(-dg_irxn);

	wTrial = w0 * ac0 / ac * dd;
	*wt_kspec = wTrial;
	Vphase->setMolesFromVCS(VCS_DATA_PTR(m_molNumSpecies_new));
	Vphase->sendToVCSActCoeff(VCS_DATA_PTR(m_actCoeffSpecies_new));
	double ac1 = m_actCoeffSpecies_new[kspec];
	double acprime = 0.0;
	if (fabs(wTrial - w0) > 1.0E-8 * w0) {
	  acprime = (ac1 - ac0) / (wTrial - w0);
	}
	double jac = acprime * wTrial + ac1;
	double fTrial = ac1 * wTrial - ac0*w0*dd;
	double w2 = wTrial - fTrial / jac;
	if (w2 > 100.*w0) {
	  *wt_kspec = 100.0 * w0;
	} else if (100. * w2 < w0) {
	  *wt_kspec = 0.01 * w0;
	} else {
	  *wt_kspec = w2;
	}
      }

      if ((*wt_kspec) <  VCS_DELETE_MINORSPECIES_CUTOFF) {
	goto L_ZERO_SPECIES;
      }
      dx = (*wt_kspec) - w_kspec;
      (*ds_kspec) = dx;
      return dx;
      /*
       *
       *  Alternate return based for cases where we need to delete the species 
       *  from the current list of active species, because its concentration
       *  has gotten too small.
       */
    L_ZERO_SPECIES: ;
      *do_delete = TRUE;
      dx = - w_kspec;
      (*ds_kspec) = dx;
      return dx;
    } 
    else {
      /*
       * Voltage calculation 
       *      HKM -> Need to check the sign
       */
      dx = m_deltaGRxn_new[irxn]/ Faraday_dim;
#ifdef DEBUG_MODE
      sprintf(ANOTE,"voltage species alternative calc");
#endif
    }
    return dx;
  }
  /*****************************************************************************/

  int VCS_SOLVE::delta_species(int kspec, double *delta_ptr)

    /************************************************************************
     *
     * delta_species():
     *
     *  Change the concentration of a species by delta moles. 
     *  Make sure to conserve
     *  elements and keep track of the total kmoles in all phases.
     *
     *  return:
     *      1: succeeded
     *      0: failed.
     ************************************************************************/
  {
    int irxn = kspec - m_numComponents;
    int retn = 1;
    int j;
    double tmp;
    double delta = *delta_ptr;
    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
      /*
       * Attempt the given dx. If it doesn't work, try to see if a smaller
       * one would work,
       */
      double dx = delta;
      double *sc_irxn = m_stoichCoeffRxnMatrix[irxn];
      for (j = 0; j < m_numComponents; ++j) {
	if (m_molNumSpecies_old[j] > 0.0) {
	  tmp = sc_irxn[j] * dx;
	  if (-tmp > m_molNumSpecies_old[j]) {
	    retn = 0;
	    dx = MIN(dx, - m_molNumSpecies_old[j] / sc_irxn[j]); 
	  }
	}
	/*
	 * If the component has a zero concentration and is a reactant
	 * in the formation reaction, then dx == 0.0, and we just return.
	 */
	if (m_molNumSpecies_old[j] <= 0.0) {
	  if (sc_irxn[j] < 0.0) {
	    *delta_ptr = 0.0;
	    return 0;
	  }
	}
      }
      /*
       * ok, we found a positive dx. implement it.
       */
      *delta_ptr = dx;
      m_molNumSpecies_old[kspec] += dx;
      int iph = m_phaseID[kspec];
      m_tPhaseMoles_old[iph] += dx;
      for (j = 0; j < m_numComponents; ++j) {
	iph = m_phaseID[j];
	tmp = sc_irxn[j] * dx;
	m_molNumSpecies_old[j] += tmp;
	m_tPhaseMoles_old[iph] += tmp;
	if (m_molNumSpecies_old[j] < 0.0) {
	  m_molNumSpecies_old[j] = 0.0;
	}
      }
    }
    return retn;
  }

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  int VCS_SOLVE::zero_species(int kspec) 

    /************************************************************************
     *
     * zero_species:
     *
     *  Zero out the concentration of a species. Make sure to conserve
     *  elements and keep track of the total moles in all phases.
     *       w[]
     *       m_tPhaseMoles_old[]
     *
     *  return:
     *      1: succeeded
     *      0: failed.
     ************************************************************************/
  {
    int retn = 1;
    /*
     * Calculate a delta that will eliminate the species.
     */
    if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
      double dx = -(m_molNumSpecies_old[kspec]);
      if (dx != 0.0) {
	retn = delta_species(kspec, &dx);
	if (!retn) {
	  plogf("zero_species: Couldn't zero the species %d, "
		"did delta of %g. orig conc of %g\n",
		kspec, dx, m_molNumSpecies_old[kspec] + dx);
	}
      }
    }
    return retn;
  }

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/
   
  int VCS_SOLVE::delete_species(int kspec)
   
    /************************************************************************
     *
     * delete_species:
     *
     * Rearrange data when species is added or removed. The Lth species is 
     * moved to the back of the species vector. The back of the species 
     * vector is indicated by the value of MR, the current number of 
     * active species in the mechanism. 
     *
     * Input 
     *     kspec  = species number 
     * Return value 
     *     The return is true when the current number of 
     * noncomponent species is equal to zero. A recheck of deleted species 
     * is carried out in the main code. 
     *************************************************************************/
  {
    int klast = m_numSpeciesRdc - 1;
    int iph = m_phaseID[kspec];
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    int irxn = kspec - m_numComponents;     /* This is the noncomponent rxn index */
    /*
     * Zero the concentration of the species.
     *     -> This zeroes w[kspec] and modifies m_tPhaseMoles_old[]
     */
    int retn = zero_species(kspec);
    if (! retn) {
      plogf("Failed to delete a species!\n");
      exit(-1);
    }
    /*
     *    Decrement the minor species counter if the current species is
     *    a minor species
     */
    if (m_rxnStatus[irxn] != VCS_SPECIES_MAJOR) --(m_numRxnMinorZeroed);
    m_rxnStatus[irxn] = VCS_SPECIES_DELETED;
    m_deltaGRxn_new[irxn] = 0.0;
    m_deltaGRxn_old[irxn] = 0.0;
    m_feSpecies_curr[kspec] = 0.0;
    m_feSpecies_old[kspec] = 0.0;
    m_molNumSpecies_new[kspec] = 0.0;
    /*
     *    Rearrange the data if the current species isn't the last active
     *    species.
     */
    if (kspec != klast) {
      vcs_switch_pos(TRUE, klast, kspec);
    }
    /*  
     *       Adjust the total moles in a phase downwards. 
     */
    Vphase->setMolesFromVCSCheck(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_tPhaseMoles_old));
  
    /* 
     *    Adjust the current number of active species and reactions counters 
     */
    --(m_numRxnRdc);
    --(m_numSpeciesRdc);
   
    /*
     *    Check to see whether we have just annihilated a multispecies phase.
     *    If it is extinct, call the delete_multiphase() function.
     */
    if (! m_SSPhase[klast]) {
      if (Vphase->Existence != 2) {
	Vphase->Existence = 0;
	for (kspec = 0; kspec < m_numSpeciesRdc; kspec++) {
	  if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	    if (m_phaseID[kspec] == iph) {
	      if (m_molNumSpecies_old[kspec] > 0.0) {
		Vphase->Existence = 1;
		break;
	      }
	    }
	  }
	}
	if (Vphase->Existence == 0) {
	  delete_multiphase(iph);
	}
      }
    }
    /* 
     *    When the total number of noncomponent species is zero, we 
     *    have to signal the calling code 
     */
    return (m_numRxnRdc == 0);
  } 
  /***************************************************************************/

  /*
   *  
   *  reinsert_deleted():
   *
   * irxn = id of the noncomponent species formation reaction for the
   *        species to be added in.
   *
   * We make decisions on the initial mole number, and major-minor status 
   * here. We also fix up the total moles in a phase.
   *
   * The algorithm proceeds to implement these decisions in the previous 
   * position of the species. Then, vcs_switch_pos is called to move the
   * species into the last active species slot, incrementing the number 
   * of active species at the same time.
   *
   * This routine is responsible for the global data manipulation only.
   */
  void VCS_SOLVE::vcs_reinsert_deleted(int kspec) { 
    int i, k, irxn = kspec - m_numComponents;
    int *phaseID = VCS_DATA_PTR(m_phaseID);
    double dx;
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Add back a deleted species: %-12s\n", m_speciesName[kspec].c_str());
    }
#endif
    /*
     * Set the species back to minor species status
     *  this adjusts m_molNumSpecies_old[] and m_tPhaseMoles_old[]
     * HKM -> make this a relative mole number!
     */
    dx = VCS_DELETE_SPECIES_CUTOFF * 10.;
    delta_species(kspec, &dx);
    m_rxnStatus[irxn] = VCS_SPECIES_MINOR;

    if (m_SSPhase[kspec]) {
      m_rxnStatus[irxn] = VCS_SPECIES_MAJOR;
      --(m_numRxnMinorZeroed);
    }
    int iph = m_phaseID[kspec];
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    Vphase->setMolesFromVCSCheck(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_tPhaseMoles_old));
    /*
     *   We may have popped a multispecies phase back 
     *   into existence. If we did, we have to check 
     *   the other species in that phase.
     *       Take care of the m_rxnStatus[] flag.
     *       The value of m_rxnStatus[] must change from 
     *       VCS_SPECIES_ZEROEDPHASE to VCS_SPECIES_ZEROEDMS
     *       for those other species.
     */
    if (! m_SSPhase[kspec]) {
      if (Vphase->Existence == 0) {
	Vphase->Existence = 1;
	for (k = 0; k < m_numSpeciesTot; k++) {
	  if (phaseID[k] == iph) {
	    i = k - m_numComponents;
	    if (m_rxnStatus[i] == VCS_SPECIES_ZEROEDPHASE) 
	      m_rxnStatus[i] = VCS_SPECIES_ZEROEDMS;
	  }
	}
      }
    } else {
      Vphase->Existence = 1;
    }
   
    ++(m_numRxnRdc);
    ++(m_numSpeciesRdc);
    ++(m_numRxnMinorZeroed);
   
    if (kspec != (m_numSpeciesRdc - 1)) {
      /*
       *  Rearrange both the species and the non-component global data
       */
      vcs_switch_pos(TRUE, (m_numSpeciesRdc - 1), kspec);
    }
  } /* vcs_reinsert_deleted() */

  /****************************************************************************
   *
   * delete_multiphase():
   *
   *  This routine handles the bookkeepking involved with the
   *  deletion of multiphase phases from the 
   *  problem. When they are deleted, all of their species become active
   *  species, even though their mole numbers are set to zero.
   *  The routine does not make the decision to eliminate multiphases.
   *
   *   Note, species in phases with zero mole numbers are still
   *   considered active. Whether the phase pops back into
   *   existence or not is checked as part of the main iteration
   *   loop.
   */
  void VCS_SOLVE::delete_multiphase(int iph) {
    int kspec, j, irxn;
    double dx;
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    /*
     * set the phase existence flag to dead
     */
    Vphase->Existence = 0;
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- delete_multiphase %d, %s\n", iph, Vphase->PhaseName.c_str());
    }
#endif
    /*
     * Zero out the total moles counters for the phase
     */
    m_tPhaseMoles_old[iph] = 0.0;
    m_tPhaseMoles_new[iph] = 0.0;
    m_deltaPhaseMoles[iph] = 0.0;
   
    /*
     * Loop over all of the active species in the phase.
     */
    for (kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
      if (m_phaseID[kspec] == iph) {
	if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  irxn = kspec - m_numComponents;
	  /*
	   * calculate an extent of rxn, dx, that zeroes out the species.
	   */
	  dx = - (m_molNumSpecies_old[kspec]);
	  /*
	   * Set the mole numbers of that species to zero.
	   */
	  m_molNumSpecies_old[kspec]  = 0.0;
	  m_molNumSpecies_new[kspec] = 0.0;
	  m_deltaMolNumSpecies[kspec] = 0.0;
	  /*
	   * Change the status flag of the species to that of an
	   * zeroed phase
	   */
	  m_rxnStatus[irxn] = VCS_SPECIES_ZEROEDPHASE;
	  /*
	   *  changed the component mole numbers to account for the
	   *  final extent of reaction. Make sure to keep component
	   *  mole numbers constant.
	   *  HKM -> note, this will cause a loss of moles!
	   */
	  for (j = 0; j < m_numComponents; ++j) {
	    m_molNumSpecies_old[j] += m_stoichCoeffRxnMatrix[irxn][j] * dx;
	    if (m_molNumSpecies_old[j] < 0.0) {
	      m_molNumSpecies_old[j] = 0.0;
	    }
	  }
	}
      }
    }
    /*
     *   Loop over all of the inactive species in the phase:
     *   Right now we reinstate all species in a deleted multiphase.
     *   We may only want to reinstate the "major ones" in the future.
     *   Note, species in phases with zero mole numbers are still
     *   considered active. Whether the phase pops back into
     *   existence or not is checked as part of the main iteration
     *   loop.
     */
    for (kspec = m_numSpeciesRdc; kspec < m_numSpeciesTot; ++kspec) {
      if (m_phaseID[kspec] == iph) {
	irxn = kspec - m_numComponents;
	m_molNumSpecies_old[kspec]  = 0.0;
	m_molNumSpecies_new[kspec] = 0.0;
	m_deltaMolNumSpecies[kspec] = 0.0;
	m_rxnStatus[irxn] = VCS_SPECIES_ZEROEDPHASE;
	 
	++(m_numRxnRdc);
	++(m_numSpeciesRdc);
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  plogf("   ---    Make %s", m_speciesName[kspec].c_str()); 
	  plogf(" an active but zeroed species because its phase "
		"was zeroed\n");
	}
#endif
	if (kspec != (m_numSpeciesRdc - 1)) {
	  /*
	   *  Rearrange both the species and the non-component global data
	   */
	  vcs_switch_pos(TRUE, (m_numSpeciesRdc - 1), kspec);
	}
      }
    }
    /*
     * Upload the state to the VP object
     */
    Vphase->setMolesFromVCSCheck(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_tPhaseMoles_old), iph);

  } /* delete_multiphase() *****************************************************/
   
  /*****************************************************************************
   *
   * recheck_deleted:
   *
   * Recheck deleted species in multispecies phases.
   *
   * HKM -> This algorithm needs to be updated for activity coefficients
   */
  int VCS_SOLVE::recheck_deleted(void)
  {
    int iph, kspec, irxn, npb;
    double *xtcutoff = VCS_DATA_PTR(TmpPhase);
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Start rechecking deleted species in multispec phases\n");
    }
#endif
    if (m_numSpeciesRdc == m_numSpeciesTot) return 0;
    /*
     * Use the standard chemical potentials for the chemical potentials
     * of deleted species. Then, calculate Delta G for 
     * for formation reactions
     */
    for (kspec = m_numSpeciesRdc; kspec < m_numSpeciesTot; ++kspec) {
      m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec];
    }
    /*
     *      Recalculate the DeltaG's of the formation reactions for the
     *      deleted species in the mechanism
     */
    vcs_deltag(0, true);
   
    for (iph = 0; iph < m_numPhases; iph++) {
      if (m_tPhaseMoles_old[iph] > 0.0)  
	xtcutoff[iph] = log (m_tPhaseMoles_old[iph] / VCS_DELETE_SPECIES_CUTOFF);
      else
	xtcutoff[iph] = 0.0;
    }
    /*
     *  
     *   We are checking the equation:
     *
     *         sum_u = sum_j_comp [ sigma_i_j * u_j ] 
     *               = u_i_O + log((AC_i * W_i)/m_tPhaseMoles_old) 
     *
     *   by first evaluating: 
     *
     *          DG_i_O = u_i_O - sum_u. 
     *
     *   Then, if TL is zero, the phase pops into existence if DG_i_O < 0.
     *   Also, if the phase exists, then we check to see if the species
     *   can have a mole number larger than  VCS_DELETE_SPECIES_CUTOFF 
     *   (default value = 1.0E-32).
     *
     *   HKM:
     *       This seems to be an inconsistency in the algorithm here that needs
     *   correcting. The requirement above may bypass some multiphases which
     *   should exist. The real requirement for the phase to exist is:
     *
     *           sum_i_in_phase [ exp(-DG_i_O) ] >= 1.0
     *   
     *   Thus, we need to amend th code. Also nonideal solutions will tend to 
     *   complicate matters severely also.
     */
    npb = 0;
    for (irxn = m_numRxnRdc; irxn < m_numRxnTot; ++irxn) {
      kspec = m_indexRxnToSpecies[irxn];
      iph = m_phaseID[kspec];
      if (m_tPhaseMoles_old[iph] == 0.0) {
	if (m_deltaGRxn_new[irxn] < 0.0) {
	  vcs_reinsert_deleted(kspec);
	  npb++;
	} else {
	  m_molNumSpecies_old[kspec] = 0.0;
	}
      } else if (m_tPhaseMoles_old[iph] > 0.0) {
	if (m_deltaGRxn_new[irxn] < xtcutoff[iph]) {
	  vcs_reinsert_deleted(kspec);
	  npb++;
	}
      }
    }
    return npb;
  }
  /***********************************************************************************/

  //  Provide an estimate for the deleted species in phases that
  //  are not zeroed out
  /*
   *  Try to add back in all deleted species. An estimate of the kmol numbers
   *  are obtained and the species is added back into the equation system,
   *  into the old state vector.
   */
  int VCS_SOLVE::vcs_add_all_deleted() {
    int iph, kspec, retn;
    if (m_numSpeciesRdc == m_numSpeciesTot) return 0;
    /*
     * Use the standard chemical potentials for the chemical potentials
     * of deleted species. Then, calculate Delta G for 
     * for formation reactions.
     *     We are relying here on a old saved value of m_actCoeffSpecies_old[kspec]
     *  being sufficiently good. Note, we will recalculate everything at the
     *  end of the routine.
     */
    for (kspec = m_numSpeciesRdc; kspec < m_numSpeciesTot; ++kspec) {
      iph = m_phaseID[kspec];
      m_feSpecies_curr[kspec] = (m_SSfeSpecies[kspec] + log(m_actCoeffSpecies_old[kspec])
				 - SpecLnMnaught[kspec] 
				 + m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iph]);
    }
    /*
     *      Recalculate the DeltaG's of the formation reactions for the
     *      deleted species in the mechanism
     */
    vcs_deltag(0, true);
  
    for (int irxn = m_numRxnRdc; irxn < m_numRxnTot; ++irxn) {
      kspec = m_indexRxnToSpecies[irxn];
      iph = m_phaseID[kspec];
      if (m_tPhaseMoles_old[iph] > 0.0) {
	double maxDG = MIN(m_deltaGRxn_new[irxn], 300.0);

	double dx = m_tPhaseMoles_old[iph] * exp(- maxDG);
	retn = delta_species(kspec, &dx);
	if (retn == 0) {
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl) {
	    plogf("  --- add_deleted(): delta_species() failed for species %s (%d) with mol number %g\n",
		 m_speciesName[kspec].c_str(), kspec, dx);
	  }
#endif
	  if (dx > 1.0E-50) {
	    dx = 1.0E-50;
	    retn = delta_species(kspec, &dx);
#ifdef DEBUG_MODE
	    if (retn == 0) {
	      if (vcs_debug_print_lvl) {
		plogf("  --- add_deleted(): delta_species() failed for species %s (%d) with mol number %g\n",
		      m_speciesName[kspec].c_str(), kspec, dx);
	      }
	    }
#endif
	  }
	}
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  if (retn != 0) {
	    plogf("  --- add_deleted():  species %s added back in with mol number %g",
		  m_speciesName[kspec].c_str(), dx);
	    plogendl();
	  } else {
	    plogf("  --- add_deleted():  species %s failed to be added  back in");
	    plogendl();
	  }
	}
#endif
      }
    }

    vcs_dfe(VCS_DATA_PTR(m_molNumSpecies_old), VCS_STATECALC_OLD, 0, 0, m_numSpeciesTot);
    vcs_deltag(0, true);

    retn = 0;
    for (int irxn = m_numRxnRdc; irxn < m_numRxnTot; ++irxn) {
      kspec = m_indexRxnToSpecies[irxn];
      iph = m_phaseID[kspec];
      if (m_tPhaseMoles_old[iph] > 0.0) {
	if (fabs(m_deltaGRxn_old[irxn]) > m_tolmin) {
	  retn++;
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {	
	    plogf("  --- add_deleted():  species %s with mol number %g not converged: DG = %g",
		  m_speciesName[kspec].c_str(), m_molNumSpecies_old[kspec],
		  m_deltaGRxn_old[irxn]);
	    plogendl();
	  }
#endif
	}
      }
    }
    return retn;
  }
  /***********************************************************************************/

  /* globalStepDamp
   *
   *  Convergence Forcer:
   *
   *  This routine optimizes the minimization of the total gibbs free
   *  energy by making sure the slope of the following functional stays
   *  negative:
   *
   *                 d_Gibbs/ds = sum_k( m_deltaGRxn * m_deltaMolNumSpecies[k] )
   *
   *  along the current direction m_deltaMolNumSpecies[], by choosing a value, al: (0<al<1)
   *  such that the a parabola approximation to Gibbs(al) fit to the 
   *  end points al = 0 and al = 1 is minimizied.
   *      s1 = slope of Gibbs function at al = 0, which is the previous
   *           solution = d(Gibbs)/d(al).
   *      s2 = slope of Gibbs function at al = 1, which is the current
   *           solution = d(Gibbs)/d(al).
   *  Only if there has been an inflection point (i.e., s1 < 0 and s2 > 0),
   *  does this code section kick in. It finds the point on the parabola
   *  where the slope is equal to zero.
   *
   */
  int VCS_SOLVE::globStepDamp(int iti) {
    double s1, s2, al;
    int irxn, kspec, iph;
    double *dptr = VCS_DATA_PTR(m_deltaGRxn_new);
      
    /* *************************************************** */
    /* **** CALCULATE SLOPE AT END OF THE STEP  ********** */
    /* *************************************************** */
    s2 = 0.0;
    for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
      kspec = irxn + m_numComponents;
      s2 += dptr[irxn] * m_deltaMolNumSpecies[kspec];
    }


    /* *************************************************** */
    /* **** CALCULATE ORIGINAL SLOPE ********************* */
    /* ************************************************** */
    s1 = 0.0;
    dptr = VCS_DATA_PTR(m_deltaGRxn_old);
    for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
      kspec = irxn + m_numComponents;
      s1 += dptr[irxn] * m_deltaMolNumSpecies[kspec];
    }

#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- subroutine FORCE: Beginning Slope = %g\n", s1);
      plogf("   --- subroutine FORCE: End Slope       = %g\n", s2);
    }
#endif

    if (s1 > 0.0) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- subroutine FORCE produced no adjustments,");
	if (s1 < 1.0E-40) {
	  plogf(" s1 positive but really small");
	} else {
	  plogf(" failed s1 test");
	}
	plogendl();
      }
#endif
      return FALSE;
    }

    if (s2 <= 0.0) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- subroutine FORCE produced no adjustments, s2 < 0");
	plogendl();
      }
#endif     
      return FALSE;
    }

    /* *************************************************** */
    /* **** FIT PCJ2822ARABOLA ********************************* */
    /* *************************************************** */
    al = 1.0;
    if (fabs(s1 -s2) > 1.0E-200) {
      al = s1 / (s1 - s2);
    }
    if (al >= 0.95 || al < 0.0) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("   --- subroutine FORCE produced no adjustments (al = %g)\n", al);
      }
#endif
      return FALSE;
    }
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- subroutine FORCE produced a damping factor = %g\n", al);
    }
#endif

    /* *************************************************** */
    /* **** ADJUST MOLE NUMBERS, CHEM. POT *************** */
    /* *************************************************** */
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      vcs_dcopy(VCS_DATA_PTR(m_deltaGRxn_tmp), VCS_DATA_PTR(m_deltaGRxn_new),  
		m_numRxnRdc);
    }
#endif
   
    dptr = VCS_DATA_PTR(m_molNumSpecies_new);
    for (kspec = 0; kspec < m_numSpeciesRdc; ++kspec) {
      m_molNumSpecies_new[kspec] = m_molNumSpecies_old[kspec] + al * m_deltaMolNumSpecies[kspec];
    }
    for (iph = 0; iph < m_numPhases; iph++) {
      m_tPhaseMoles_new[iph] = m_tPhaseMoles_old[iph] + al * m_deltaPhaseMoles[iph];
    }
    vcs_updateVP(1);
   
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- subroutine FORCE adjusted the mole "
	    "numbers, AL = %10.3f\n", al);
    }
#endif
    /*
     *           Because we changed the mole numbers, we need to 
     *           calculate the chemical potentials again. If a major-
     *           only step is being carried out, then we don't need to
     *           update the minor noncomponents. 
     */
    // vcs_dfe(dptr, VCS_STATECALC_NEW, iti, 0, m_numSpeciesRdc);
    vcs_dfe(dptr, VCS_STATECALC_NEW, 0, 0, m_numSpeciesRdc);
    /*
     *           Evaluate DeltaG for all components if ITI=0, and for 
     *           major components only if ITI NE 0 
     */
    // vcs_deltag(iti, false);
    vcs_deltag(0, false);

    dptr = VCS_DATA_PTR(m_deltaGRxn_new);
    s2 = 0.0;
    for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
      kspec = irxn + m_numComponents;
      s2 += dptr[irxn] * m_deltaMolNumSpecies[kspec];
    }


#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- subroutine FORCE: Adj End Slope   = %g", s2);
      plogendl();
    }
#endif
    return TRUE;
  }
  


  /*
   * vcs_RxnStepSizes():
   *
   * Calculates formation reaction step sizes.
   * This is equation 6.4-16, p. 143 in Smith and Missen. 
   *
   * Output 
   * ------- 
   * m_deltaMolNumSpecies(I) : reaction adjustments, where I refers to the Ith species
   *         formation reaction. This is adjustment is for species
   *         i + M, where M is the number of components.
   * Special branching occurs sometimes. This causes the component basis 
   * to be reevaluated 
   *     return = 0 : normal return
   *              1 : A single species phase species has been zeroed out
   *                  in this routine. The species is a noncomponent 
   *              2 : Same as one but, the zeroed species is a component. 
   */
  int VCS_SOLVE::vcs_RxnStepSizes() {  
    int  j, irxn, kspec, soldel = 0, iph;
    double s, xx, dss;
    int k = 0;
    vcs_VolPhase *Vphase = 0;
    double *dnPhase_irxn;
#ifdef DEBUG_MODE
    char ANOTE[128];
    if (vcs_debug_print_lvl >= 2) {
      plogf("   "); for (j = 0; j < 82; j++) plogf("-"); plogf("\n");
      plogf("   --- Subroutine vcs_RxnStepSizes called - Details:\n");
      plogf("   "); for (j = 0; j < 82; j++) plogf("-"); plogf("\n");
      plogf("   --- Species        KMoles     Rxn_Adjustment    DeltaG"
	    "   | Comment\n");
    }
#endif
    /*
     * We update the matrix dlnActCoeffdmolNumber[][] at the
     * top of the loop, when necessary
     */
    if (m_useActCoeffJac) {
      vcs_CalcLnActCoeffJac(VCS_DATA_PTR(m_molNumSpecies_old));
    }
    /************************************************************************
     ******** LOOP OVER THE FORMATION REACTIONS *****************************
     ************************************************************************/

    for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
#ifdef DEBUG_MODE
      sprintf(ANOTE,"Normal Calc");
#endif

      kspec = m_indexRxnToSpecies[irxn];

      if (m_speciesUnknownType[kspec] !=  VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {

	dnPhase_irxn = m_deltaMolNumPhase[irxn];
      
	if (m_molNumSpecies_old[kspec] == 0.0 && (! m_SSPhase[kspec])) {
	  /********************************************************************/
	  /******* MULTISPECIES PHASE WITH total moles equal to zero *********/
	  /*******************************************************************/
	  /* 
	   *   If dg[irxn] is negative, then the multispecies phase should
	   *   come alive again. Add a small positive step size to 
	   *   make it come alive. 
	   */
	  if (m_deltaGRxn_new[irxn] < -1.0e-4) {
	    /*
	     * First decide if this species is part of a multiphase that
	     * is nontrivial in size.
	     */
	    iph = m_phaseID[kspec];
	    double tphmoles = m_tPhaseMoles_old[iph];
	    double trphmoles = tphmoles / m_totalMolNum;
	    if (trphmoles > VCS_DELETE_PHASE_CUTOFF) {
	      m_deltaMolNumSpecies[kspec] = m_totalMolNum * VCS_SMALL_MULTIPHASE_SPECIES;
#ifdef DEBUG_MODE
	      sprintf(ANOTE,
		      "MultSpec: small species born again DG = %11.3E", 
		      m_deltaGRxn_new[irxn]);
#endif
	    } else {
#ifdef DEBUG_MODE
	      sprintf(ANOTE, "MultSpec: phase come alive DG = %11.3E", 
		      m_deltaGRxn_new[irxn]);   
#endif
	      Vphase = m_VolPhaseList[iph];
	      int numSpPhase = Vphase->NVolSpecies;
	      m_deltaMolNumSpecies[kspec] = 
		m_totalMolNum * 10.0 * VCS_DELETE_PHASE_CUTOFF / numSpPhase;
	    }
	    --(m_numRxnMinorZeroed);
	  } else {
#ifdef DEBUG_MODE
	    sprintf(ANOTE, "MultSpec: still dead DG = %11.3E", m_deltaGRxn_new[irxn]);       
#endif
	    m_deltaMolNumSpecies[kspec] = 0.0;
	  }
	} else {
	  /********************************************************************/
	  /************************* REGULAR PROCESSING            ************/
	  /********************************************************************/
	  /*
	   *     First take care of cases where we want to bail out
	   *
	   *
	   *     Don't bother if superconvergence has already been achieved 
	   *     in this mode.
	   */
	  if (fabs(m_deltaGRxn_new[irxn]) <= m_tolmaj2) {
#ifdef DEBUG_MODE
	    sprintf(ANOTE,"Skipped: superconverged DG = %11.3E", m_deltaGRxn_new[irxn]);
	    if (vcs_debug_print_lvl >= 2) {
	      plogf("   --- %-12.12s", m_speciesName[kspec].c_str()); 
	      plogf("  %12.4E %12.4E %12.4E | %s\n",  
		    m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
		    m_deltaGRxn_new[irxn], ANOTE);
	    }
#endif		    
	    continue;
	  }
	  /*
	   *     Don't calculate for minor or nonexistent species if      
	   *     their values are to be decreasing anyway.                
	   */
	  if ((m_rxnStatus[irxn] != VCS_SPECIES_MAJOR) && (m_deltaGRxn_new[irxn] >= 0.0)) {
#ifdef DEBUG_MODE
	    sprintf(ANOTE,"Skipped: IC = %3d and DG >0: %11.3E", 
		    m_rxnStatus[irxn], m_deltaGRxn_new[irxn]);
	    if (vcs_debug_print_lvl >= 2) {
	      plogf("   --- %-12.12s", m_speciesName[kspec].c_str());
	      plogf("  %12.4E %12.4E %12.4E | %s\n", 
		    m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec], 
		    m_deltaGRxn_new[irxn], ANOTE);
	    }
#endif		    
	    continue;
	  }
	  /*
	   *     Start of the regular processing
	   */
	  if (m_SSPhase[kspec]) {
	    s = 0.0; 
	  } else {
	    s = 1.0 / m_molNumSpecies_old[kspec] ;
	  }
	  for (j = 0; j < m_numComponents; ++j) {
	    if (!m_SSPhase[j]) {
	      if (m_molNumSpecies_old[j] > 0.0) {
		s += SQUARE(m_stoichCoeffRxnMatrix[irxn][j]) / m_molNumSpecies_old[j];
	      }
	    }
	  }
	  for (j = 0; j < m_numPhases; j++) {
	    Vphase = m_VolPhaseList[j];
	    if (! Vphase->SingleSpecies) {
	      if (m_tPhaseMoles_old[j] > 0.0) 
		s -= SQUARE(dnPhase_irxn[j]) / m_tPhaseMoles_old[j];
	    }
	  }
	  if (s != 0.0) {
	    /*
	     *  Take into account of the
	     *  derivatives of the activity coefficients with respect to the
	     *  mole numbers, even in our diagonal approximation.
	     */
	    if (m_useActCoeffJac) {
	      double s_old = s;
	      s = vcs_Hessian_diag_adj(irxn, s_old);
#ifdef DEBUG_MODE
	      if (s_old != s) {
		sprintf(ANOTE, "Normal calc: diag adjusted from %g "
			"to %g due to act coeff",  s_old, s);
	      }
#endif
	    }
	  
	    m_deltaMolNumSpecies[kspec] = -m_deltaGRxn_new[irxn] / s; 
	    // New section to do damping of the m_deltaMolNumSpecies[] 
	    /*
	     * 
	     */
	    for (j = 0; j < m_numComponents; ++j) {
	      double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
	      if (stoicC != 0.0) {
		double negChangeComp = - stoicC * m_deltaMolNumSpecies[kspec];
		if (negChangeComp > m_molNumSpecies_old[j]) {
		  if (m_molNumSpecies_old[j] > 0.0) {
#ifdef DEBUG_MODE
		    sprintf(ANOTE, "Delta damped from %g "
			    "to %g due to component %d (%10s) going neg", m_deltaMolNumSpecies[kspec],
			    -m_molNumSpecies_old[j]/stoicC, j,  m_speciesName[j].c_str());
#endif
		    m_deltaMolNumSpecies[kspec] = - m_molNumSpecies_old[j] / stoicC; 
		  } else {
#ifdef DEBUG_MODE
		    sprintf(ANOTE, "Delta damped from %g "
			    "to %g due to component %d (%10s) zero", m_deltaMolNumSpecies[kspec],
			    -m_molNumSpecies_old[j]/stoicC, j,  m_speciesName[j].c_str());
#endif
		    m_deltaMolNumSpecies[kspec] = 0.0;
		  }
		}
	      }
	    }
	    // Implement a damping term that limits m_deltaMolNumSpecies to the size of the mole number
	    if (-m_deltaMolNumSpecies[kspec] > m_molNumSpecies_old[kspec]) {
#ifdef DEBUG_MODE
	      sprintf(ANOTE, "Delta damped from %g "
		      "to %g due to %s going negative", m_deltaMolNumSpecies[kspec],
		      -m_molNumSpecies_old[kspec],  m_speciesName[kspec].c_str());
#endif
	      m_deltaMolNumSpecies[kspec] = -m_molNumSpecies_old[kspec];
	    }

	  } else {
	    /* ************************************************************ */
	    /* **** REACTION IS ENTIRELY AMONGST SINGLE SPECIES PHASES **** */
	    /* **** DELETE ONE OF THE PHASES AND RECOMPUTE BASIS  ********* */
	    /* ************************************************************ */
	    /* 
	     *     Either the species L will disappear or one of the 
	     *     component single species phases will disappear. The sign 
	     *     of DG(I) will indicate which way the reaction will go. 
	     *     Then, we need to follow the reaction to see which species 
	     *     will zero out first. 
	     *      -> The species to be zeroed out will be "k".
	     */
	    if (m_deltaGRxn_new[irxn] > 0.0) {
	      dss = m_molNumSpecies_old[kspec];
	      k = kspec;
	      for (j = 0; j < m_numComponents; ++j) {
		if (m_stoichCoeffRxnMatrix[irxn][j] > 0.0) {
		  xx = m_molNumSpecies_old[j] / m_stoichCoeffRxnMatrix[irxn][j];
		  if (xx < dss) {
		    dss = xx;
		    k = j;
		  }
		}
	      }
	      dss = -dss;
	    } else {
	      dss = 1.0e10;
	      for (j = 0; j < m_numComponents; ++j) {
		if (m_stoichCoeffRxnMatrix[irxn][j] < 0.0) {
		  xx = -m_molNumSpecies_old[j] / m_stoichCoeffRxnMatrix[irxn][j];
		  if (xx < dss) {
		    dss = xx;
		    k = j;
		  }
		}
	      }
	    }
	    /*
	     *          Here we adjust the mole fractions 
	     *          according to DSS and the stoichiometric array 
	     *          to take into account that we are eliminating 
	     *          the kth species. DSS contains the amount 
	     *          of moles of the kth species that needs to be 
	     *          added back into the component species. 
	     */
	    if (dss != 0.0) {
	      m_molNumSpecies_old[kspec] += dss;
	      m_tPhaseMoles_old[m_phaseID[kspec]] +=  dss;
	      for (j = 0; j < m_numComponents; ++j) {
		m_molNumSpecies_old[j] += dss * m_stoichCoeffRxnMatrix[irxn][j];
		m_tPhaseMoles_old[m_phaseID[j]] +=  dss * m_stoichCoeffRxnMatrix[irxn][j];
	      }
	      m_molNumSpecies_old[k] = 0.0;
	      iph = m_phaseID[k];
	      Vphase = m_VolPhaseList[iph];
	      Vphase->Existence = 0;
	      m_tPhaseMoles_old[iph] = 0.0;
#ifdef DEBUG_MODE
	      if (vcs_debug_print_lvl >= 2) {
		plogf("   --- vcs_RxnStepSizes Special section to delete %s\n",
		      m_speciesName[k].c_str());
		plogf("   ---   Immediate return - Restart iteration\n");
	      }
#endif
	      /*
	       *            We need to immediately recompute the 
	       *            component basis, because we just zeroed 
	       *            it out. 
	       */
	      if (k != kspec) soldel = 2;
	      else            soldel = 1;
	      return soldel;
	    }
	  }
	} /* End of regular processing */
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  plogf("   --- %-12.12s", m_speciesName[kspec].c_str());
	  plogf("  %12.4E %12.4E %12.4E | %s\n", 
		m_molNumSpecies_old[kspec], m_deltaMolNumSpecies[kspec],
		m_deltaGRxn_new[irxn], ANOTE);
	}
#endif	
      } /* End of loop over m_speciesUnknownType */
    } /* End of loop over non-component stoichiometric formation reactions */
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   "); vcs_print_line("-", 82);
    }
#endif
    return soldel;
  }
  /*****************************************************************************/

  /*
   *
   * vcs_deltag:
   *
   * This subroutine calculates reaction free energy changes for 
   * all noncomponent formation reactions. Formation reactions are 
   * reactions which create each noncomponent species from the component 
   * species. M_STOICHCOEFFRXNMATRIX(J,I) are the stoichiometric coefficients for these 
   * reactions. A stoichiometric coefficient of one is assumed for 
   * species I in this reaction. 
   *
   *  INPUT 
   *    L < 0   :  Calculate reactions corresponding to 
   *               major noncomponent and zeroed species only 
   *    L = 0   :  Do all noncomponent reactions, i, between 
   *               0 <= i < irxnl 
   *    L > 0   :  Calculate reactions corresponding to 
   *               minor noncomponent and zeroed species only 
   *    irxnl : used with L = 0 to indicate upper limit.
   *
   * Note we special case one important issue.
   * If the component has zero moles, then we do not
   * allow deltaG < 0.0 for formation reactions which
   * would lead to the loss of more of the component.
   * This dG < 0.0 feeds back into the algorithm in several
   * places, and leads to a infinite loop in at least one case. 
   */
  void VCS_SOLVE::vcs_deltag(int l, bool doDeleted) {
    int iph;
    int   lneed, irxn, kspec;
    double *dtmp_ptr;
    int icase = 0;
    int irxnl = m_numRxnRdc;
    if (doDeleted) {
      irxnl = m_numRxnTot;
    }

#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Subroutine vcs_deltag called for ");
      if (l < 0) {
	plogf("major noncomponents\n");
      } else if (l == 0) {
	plogf("all noncomponents\n");     
      } else {
	plogf("minor noncomponents\n");
      }
    }
#endif
    /* ************************************************* */
    /* **** MAJORS and ZEREOD SPECIES ONLY ************* */
    /* ************************************************* */
    if (l < 0) {
      for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
	if (m_rxnStatus[irxn] != VCS_SPECIES_MINOR) {
	  icase = 0;
	  m_deltaGRxn_new[irxn] = m_feSpecies_curr[m_indexRxnToSpecies[irxn]];
	  dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
	  for (kspec = 0; kspec < m_numComponents; ++kspec) {
	    m_deltaGRxn_new[irxn] += dtmp_ptr[kspec] * m_feSpecies_curr[kspec];
	    if (m_molNumSpecies_old[kspec] < VCS_DELETE_MINORSPECIES_CUTOFF && dtmp_ptr[kspec] < 0.0) {
	      icase = 1;
	    } 
	  }
	  if (icase) {
	    m_deltaGRxn_new[irxn] = MAX(0.0, m_deltaGRxn_new[irxn]);
	  }
	}
      }
    } else if (l == 0) {
      /* ************************************************* */
      /* **** ALL REACTIONS ****************************** */
      /* ************************************************* */
      for (irxn = 0; irxn < irxnl; ++irxn) {
	icase = 0;
	m_deltaGRxn_new[irxn] = m_feSpecies_curr[m_indexRxnToSpecies[irxn]];
	dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
	for (kspec = 0; kspec < m_numComponents; ++kspec) {
	  m_deltaGRxn_new[irxn] += dtmp_ptr[kspec] * m_feSpecies_curr[kspec];
	  if (m_molNumSpecies_old[kspec] < VCS_DELETE_MINORSPECIES_CUTOFF && 
	      dtmp_ptr[kspec] < 0.0) {
	    icase = 1;
	  }
	}
	if (icase) {
	  m_deltaGRxn_new[irxn] = MAX(0.0, m_deltaGRxn_new[irxn]);
	}
      }
    } else {
      /* ************************************************* */
      /* **** MINORS AND ZEROED SPECIES ****************** */
      /* ************************************************* */
      for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
	if (m_rxnStatus[irxn] <= VCS_SPECIES_MINOR) {
	  icase = 0;
	  m_deltaGRxn_new[irxn] = m_feSpecies_curr[m_indexRxnToSpecies[irxn]];
	  dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
	  for (kspec = 0; kspec < m_numComponents; ++kspec) {
	    m_deltaGRxn_new[irxn] += dtmp_ptr[kspec] * m_feSpecies_curr[kspec];
	    if (m_molNumSpecies_old[kspec] < VCS_DELETE_MINORSPECIES_CUTOFF && 
		dtmp_ptr[kspec] < 0.0) {
	      icase = 1;
	    } 
	  }
	  if (icase) {
	    m_deltaGRxn_new[irxn] = MAX(0.0, m_deltaGRxn_new[irxn]);
	  }
	}
      }
    }
    /* ************************************************* */
    /* **** MULTISPECIES PHASES WITH ZERO MOLES ******** */
    /* ************************************************* */
    /*
     *    Massage the free energies for species with zero mole fractions 
     *  in multispecies phases.  This section implements the 
     *  Equation 3.8-5 in Smith and Missen, p.59. 
     *  A multispecies phase will exist iff 
     *           1 < sum_i(exp(-dg_i)/AC_i) 
     *  If DG is negative then that species wants to be reintroduced into 
     *  the calculation. 
     *  For small dg_i, the expression below becomes: 
     *      1 - sum_i(exp(-dg_i)/AC_i) ~ sum_i((dg_i-1)/AC_i)  + 1
     * 
     *  So, what we are doing here is equalizing all DG's in a multispecies
     *  phase whose total mole number has already been zeroed out. 
     *  It must have to do with the case where a complete multispecies 
     *  phase is currently zeroed out. In that case, when one species 
     *  in that phase has a negative DG, then the phase should kick in. 
     *  This code section will cause that to happen, because a negative 
     *  DG will dominate the calculation of SDEL. Then, DG(I) for all 
     *  species in that phase will be forced to be equal and negative. 
     *  Thus, all species in that phase will come into being at the 
     *  same time. 
     *
     *  HKM -> The ratio of mole fractions at the reinstatement
     *         time should be equal to the normalized weighting
     *         of exp(-dg_i) / AC_i. This should be implemented.
     *
     *  HKM -> There is circular logic here. ActCoeff depends on the
     *         mole fractions of a phase that does not exist. In actuality
     *         the proto-mole fractions should be selected from the
     *         solution of a nonlinear problem with NsPhase unknowns
     *
     *              X_i = exp(-dg[irxn]) / ActCoeff_i / denom
     *
     *              where 
     *               denom = sum_i[  exp(-dg[irxn]) / ActCoeff_i  ]
     *      
     *         This can probably be solved by successive iteration.
     *         This should be implemented.
     */
    int k;
    for (iph = 0; iph < m_numPhases; iph++) {
      lneed = FALSE;
      vcs_VolPhase *Vphase = m_VolPhaseList[iph];
      if (! Vphase->SingleSpecies) {
	double sum = 0.0;
	for (k = 0; k < Vphase->NVolSpecies; k++) {
	  kspec = Vphase->IndSpecies[k];
	  if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	    sum += m_molNumSpecies_old[kspec];
	  }
	  if (sum > 0.0) break;
	}
	if (sum == 0.0) {
	  lneed = TRUE;
	}
      }
      
      if (lneed) {
	double poly = 0.0;
	for (k = 0; k < Vphase->NVolSpecies; k++) {
	  kspec = Vphase->IndSpecies[k];
	  irxn = kspec - m_numComponents;
	  if (m_deltaGRxn_new[irxn] >  50.0) m_deltaGRxn_new[irxn] =  50.0;
	  if (m_deltaGRxn_new[irxn] < -50.0) m_deltaGRxn_new[irxn] = -50.0;
	  poly += exp(-m_deltaGRxn_new[irxn])/m_actCoeffSpecies_new[kspec];
	}
	/*
	 *      Calculate m_deltaGRxn_new[] for each species in a zeroed multispecies phase.
	 *      All of the m_deltaGRxn_new[]'s will be equal. If m_deltaGRxn_new[] is negative, then
	 *      the phase will come back into existence.
	 */
	for (k = 0; k < Vphase->NVolSpecies; k++) {
	  kspec = Vphase->IndSpecies[k];
	  irxn = kspec - m_numComponents;
	  m_deltaGRxn_new[irxn] = 1.0 - poly;
	}

      }
    }


#ifdef DEBUG_NOT
    for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
      checkFinite(m_deltaGRxn_new[irxn]);
    }
#endif
  }
  /*****************************************************************************/

  int VCS_SOLVE::vcs_basopt(int ifirst, double aw[], double sa[], double sm[], 
			    double ss[], double test, int *usedZeroedSpecies)
   
    /**************************************************************************
     * Choose the optimum basis for the calculations. This is done by 
     * choosing the species with the largest mole fraction 
     * not currently a linear combination of the previous components. 
     * Then, calculate the stoichiometric coefficient matrix for that 
     * basis. 
     *
     * Calculates the identity of the component species in the mechanism. 
     * Rearranges the solution data to put the component data at the 
     * front of the species list. 
     *
     * Then, calculates M_STOICHCOEFFRXNMATRIX(J,I) the formation reactions for all noncomponent 
     *
     * species in the mechanism. 
     * Also calculates DNG(I) and DNL(I), the net mole change for each 
     * formation reaction. 
     * Also, initializes IR(I) to the default state. 
     *
     * Input 
     * --------- 
     * IFIRST = If true, the M_STOICHCOEFFRXNMATRIX, DNG, and DNL are not calculated. 
     * TEST   = This is a small negative number dependent upon whether 
     *          an estimate is supplied or not. 
     * W(I)   = Mole fractions which will be used to construct an 
     *          optimal basis from. 
     * 
     * Output 
     * --------- 
     * usedZeroedSpecies = If true, then a species with a zero concentration
     *                     was used as a component. The problem may be
     *                     converged.
     *
     * Other Variables 
     *  aw[i] = Mole fraction work space      (# species in length)
     *  sa[j] = Gramm-Schmidt orthog work space (nc in length)
     *  ss[j] = Gramm-Schmidt orthog work space (nc in length)
     *  sm[i+j*ne] = QR matrix work space (nc*ne in length)
     *
     *************************************************************************/
  {
    int  j, k, l, i, jl, ml, jr, lindep, irxn, kspec;
    int ncTrial;
    int juse  = -1;
    int jlose = -1;
    double *dptr, *scrxn_ptr;
    Cantera::clockWC tickTock;
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   "); for(i=0; i<77; i++) plogf("-"); plogf("\n");
      plogf("   --- Subroutine BASOPT called to ");
      if (ifirst)   plogf("calculate the number of components\n");
      else          plogf("reevaluate the components\n");
      if (vcs_debug_print_lvl >= 2) {
	plogf("\n");
	plogf("   ---       Formula Matrix used in BASOPT calculation\n");
	plogf("   ---      Active |   ");
	for (j = 0; j < m_numElemConstraints; j++) {
	  plogf("  %1d  ", m_elementActive[j]);
	}
	plogf("\n");
	plogf("   ---     Species | ");
	for (j = 0; j < m_numElemConstraints; j++) {
	  plogf(" ");
	  vcs_print_stringTrunc(m_elementName[j].c_str(), 8, 1);
	}
	plogf("\n");
	for (k = 0; k < m_numSpeciesTot; k++) {
	  plogf("   --- ");
	  vcs_print_stringTrunc(m_speciesName[k].c_str(), 11, 1);
	  plogf(" | ");
	  for (j = 0; j < m_numElemConstraints; j++) {
	    plogf(" %8.2g", m_formulaMatrix[j][k]);
	  }
	  plogf("\n");
	}
	plogendl();
      }
    }
#endif
   
    /*
     *  Calculate the maximum value of the number of components possible
     *     It's equal to the minimum of the number of elements and the
     *     number of total species.
     */
    ncTrial = MIN(m_numElemConstraints, m_numSpeciesTot);
    m_numComponents = ncTrial;
    *usedZeroedSpecies = FALSE;

    /* 
     *     Use a temporary work array for the mole numbers, aw[] 
     */
    vcs_dcopy(aw, VCS_DATA_PTR(m_molNumSpecies_old), m_numSpeciesTot);
    /*
     * Take out the Voltage unknowns from consideration
     */
    for (k = 0; k < m_numSpeciesTot; k++) {
      if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	aw[k] = test;
      }
    }

    jr = -1;
    /*
     *   Top of a loop of some sort based on the index JR. JR is the 
     *   current number of component species found. 
     */
    do {
      ++jr;
      /* - Top of another loop point based on finding a linearly */
      /* - independent species */
      do {
	/*
	 *    Search the remaining part of the mole fraction vector, AW, 
	 *    for the largest remaining species. Return its identity in K. 
	 *    The first search criteria is always the largest positive
	 *    magnitude of the mole number.
	 */
	//	k = vcs_optMax(aw, VCS_DATA_PTR(m_spSize), jr, m_numSpeciesTot);
	k = vcs_basisOptMax(aw, jr, m_numSpeciesTot);

	/*
	 * The fun really starts when you have run out of species that have a significant
	 * concentration. It becomes extremely important to make a good choice of which
	 * species you want to pick to fill out the basis. Basically, you don't want to
	 * use species with elements abundances which aren't pegged to zero. This means
	 * that those modes will never be allowed to grow. You want to have the
	 * best chance that the component will grow positively.
	 *
	 * Suppose  you start with CH4, N2, as the only species with nonzero compositions.
	 * You have the following abundances:
	 *
	 *  Abundances:
	 * ----------------
	 *     C 2.0
	 *     N 2.0
	 *     H 4.0
	 *     O 0.0
	 *
	 * For example, Make the following choice:
	 *
	 *   CH4   N2  O choose ->  OH
	 * or
	 *   CH4   N2  O choose ->  H2
	 *  
	 *  OH and H2 both fill out the basis. They will pass the algorithm. However,
	 *  choosing OH as the next species will create a situation where H2 can not
	 *  grow in concentration. This happened in practice, btw. The reason is that
	 *  the formation reaction for H2 will cause one of the component species
	 *  to go negative. 
	 *  
	 *  The basic idea here is to pick a simple species whose mole number 
	 *  can grow according to the element compositions. Candidates are still
	 *  filtered according to their linear independence.
	 *  
	 *   Note, if there is electronic charge and the electron species,
	 *   you should probably pick the electron as a component, if it
	 *   linearly independent. The algorithm below will do this automagically.
	 *
	 */
	if ((aw[k] != test) && aw[k] < VCS_DELETE_MINORSPECIES_CUTOFF) {
	  *usedZeroedSpecies = TRUE;

	  double maxConcPossKspec = 0.0;
	  double maxConcPoss = 0.0;
	  int kfound = -1;
	  int minNonZeroes = 100000;
	  int nonZeroesKspec = 0;
	  for (kspec = ncTrial; kspec < m_numSpeciesTot; kspec++) {
	    if (aw[kspec] >= 0.0) {
	      if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
		maxConcPossKspec = 1.0E10;
		nonZeroesKspec = 0;
		for (int j = 0; j < m_numElemConstraints; ++j) {
		  if (m_elementActive[j]) {
		    if (m_elType[j] == VCS_ELEM_TYPE_ABSPOS) {
		      double nu = m_formulaMatrix[j][kspec];
		      if (nu != 0.0) {
			nonZeroesKspec++;
			maxConcPossKspec = MIN(m_elemAbundancesGoal[j] / nu, maxConcPossKspec);
		      }
		    }
		  }
		}
		if ((maxConcPossKspec >= maxConcPoss) || (maxConcPossKspec > 1.0E-5)) {
		  if (nonZeroesKspec <= minNonZeroes) {
		    if (kfound < 0 || nonZeroesKspec < minNonZeroes) {
		      kfound = kspec;
		    } else {
		      // ok we are sitting pretty equal here decide on the raw ss Gibbs energy
		      if (m_SSfeSpecies[kspec] <= m_SSfeSpecies[kfound]) {
			kfound = kspec;
		      }
		    }
		  }
		  if (nonZeroesKspec < minNonZeroes) {
		    minNonZeroes = nonZeroesKspec;
		  }
		  if (maxConcPossKspec > maxConcPoss) {
		    maxConcPoss = maxConcPossKspec;
		  }
		}
	      }
	    }
	  }
	  if (kfound == -1) {
	    double gmin = 0.0;
	    kfound = k;
	    for (kspec = ncTrial; kspec < m_numSpeciesTot; kspec++) {
	      if (aw[kspec] >= 0.0) {
		irxn = kspec - ncTrial;
		if (m_deltaGRxn_new[irxn] < gmin) {
		  gmin = m_deltaGRxn_new[irxn];
		  kfound = kspec;
		}
	      }
	    }
	  }
	  k = kfound;
	}

    
	if (aw[k] == test) {
	  m_numComponents = jr;
	  ncTrial = m_numComponents;
	  int numPreDeleted = m_numRxnTot - m_numRxnRdc;
	  if (numPreDeleted != (m_numSpeciesTot - m_numSpeciesRdc)) {
	    plogf("we shouldn't be here\n");
	    exit(-1);
	  }
	  m_numRxnTot = m_numSpeciesTot - ncTrial;
	  m_numRxnRdc = m_numRxnTot - numPreDeleted;
	  m_numSpeciesRdc = m_numSpeciesTot - numPreDeleted;
	  for (i = 0; i < m_numSpeciesTot; ++i) {
	    m_indexRxnToSpecies[i] = ncTrial + i;
	  }
#ifdef DEBUG_MODE
	  if (vcs_debug_print_lvl >= 2) {
	    plogf("   ---   Total number of components found = %3d (ne = %d)\n ", 
		  ncTrial, m_numElemConstraints);
	  }
#endif
	  goto L_END_LOOP;
	}
	/*
	 *  Assign a small negative number to the component that we have
	 *  just found, in order to take it out of further consideration.
	 */
	aw[k] = test;
	/* *********************************************************** */
	/* **** CHECK LINEAR INDEPENDENCE WITH PREVIOUS SPECIES ****** */
	/* *********************************************************** */
	/*    
	 *          Modified Gram-Schmidt Method, p. 202 Dalquist 
	 *          QR factorization of a matrix without row pivoting. 
	 */
	jl = jr;
	for (j = 0; j < m_numElemConstraints; ++j) {
	  sm[j + jr*m_numElemConstraints] = m_formulaMatrix[j][k];
	}
	if (jl > 0) {
	  /*
	   *         Compute the coefficients of JA column of the 
	   *         the upper triangular R matrix, SS(J) = R_J_JR 
	   *         (this is slightly different than Dalquist) 
	   *         R_JA_JA = 1 
	   */
	  for (j = 0; j < jl; ++j) {
	    ss[j] = 0.0;
	    for (i = 0; i < m_numElemConstraints; ++i) {
	      ss[j] += sm[i + jr*m_numElemConstraints] * sm[i + j*m_numElemConstraints];
	    }
	    ss[j] /= sa[j];
	  }
	  /* 
	   *     Now make the new column, (*,JR), orthogonal to the 
	   *     previous columns
	   */
	  for (j = 0; j < jl; ++j) {
	    for (l = 0; l < m_numElemConstraints; ++l) {
	      sm[l + jr*m_numElemConstraints] -= ss[j] * sm[l + j*m_numElemConstraints];
	    }
	  }
	}
	/*
	 *        Find the new length of the new column in Q. 
	 *        It will be used in the denominator in future row calcs. 
	 */
	sa[jr] = 0.0;
	for (ml = 0; ml < m_numElemConstraints; ++ml) {
	  sa[jr] += SQUARE(sm[ml + jr*m_numElemConstraints]);
	}
	/* **************************************************** */
	/* **** IF NORM OF NEW ROW  .LT. 1E-3 REJECT ********** */
	/* **************************************************** */
	if (sa[jr] < 1.0e-6)  lindep = TRUE;
	else                  lindep = FALSE;
      } while(lindep);
      /* ****************************************** */
      /* **** REARRANGE THE DATA ****************** */
      /* ****************************************** */
      if (jr != k) {
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  plogf("   ---   %-12.12s", (m_speciesName[k]).c_str());
	  plogf("(%9.2g) replaces %-12.12s", m_molNumSpecies_old[k], m_speciesName[jr].c_str());
	  plogf("(%9.2g) as component %3d\n", m_molNumSpecies_old[jr], jr);
	}
#endif
	vcs_switch_pos(FALSE, jr, k);
	vcsUtil_dsw(aw, jr, k);
      }
#ifdef DEBUG_MODE
      else {
	if (vcs_debug_print_lvl >= 2) {
	  plogf("   ---   %-12.12s", m_speciesName[k].c_str());
	  plogf("(%9.2g) remains            ", m_molNumSpecies_old[k]);
	  plogf("              as component %3d\n", jr);
	}
      }
#endif
      /* - entry point from up above */
    L_END_LOOP: ;
      /*
       *      If we haven't found enough components, go back 
       *      and find some more. (nc -1 is used below, because
       *      jr is counted from 0, via the C convention.
       */
    } while (jr < (ncTrial-1));
   
    if (ifirst) goto L_CLEANUP;
    /* ****************************************************** */
    /* **** EVALUATE THE STOICHIOMETRY ********************** */
    /* ****************************************************** */
    /*
     *  Formulate the matrix problem for the stoichiometric
     *  coefficients. CX + B = 0
     *      C will be an nc x nc matrix made up of the formula 
     * vectors for the components.
     *      n rhs's will be solved for. Thus, B is an nc x n
     * matrix. 
     *
     * BIG PROBLEM 1/21/99:
     *
     *    This algorithm makes the assumption that the
     * first nc rows of the formula matrix aren't rank deficient.
     * However, this might not be the case. For example, assume
     * that the first element in m_formulaMatrix[] is argon. Assume that
     * no species in the matrix problem actually includes argon.
     * Then, the first row in sm[], below will be indentically
     * zero. bleh. 
     *    What needs to be done is to perform a rearrangement
     * of the ELEMENTS -> i.e. rearrange, m_formulaMatrix, sp, and m_elemAbundancesGoal, such
     * that the first nc elements form in combination with the
     * nc components create an invertible sm[]. not a small
     * project, but very doable.
     *    An alternative would be to turn the matrix problem
     * below into an ne x nc problem, and do QR elimination instead
     * of Gauss-Jordon elimination.
     *    Note the rearrangement of elements need only be done once
     * in the problem. It's actually very similar to the top of 
     * this program with ne being the species and nc being the
     * elements!!
     */
    for (j = 0; j < ncTrial; ++j) {
      for (i = 0; i < ncTrial; ++i) {
	sm[i + j*m_numElemConstraints] = m_formulaMatrix[i][j];
      }
    }
    for (i = 0; i < m_numRxnTot; ++i) {
      k = m_indexRxnToSpecies[i];
      for (j = 0; j < ncTrial; ++j) {
	m_stoichCoeffRxnMatrix[i][j] = m_formulaMatrix[j][k];
      }
    }
    /*
     *     Use Gauss-Jordon block elimination to calculate
     *     the reaction matrix, m_stoichCoeffRxnMatrix[][].
     */
    j = vcsUtil_mlequ(sm, m_numElemConstraints, ncTrial, m_stoichCoeffRxnMatrix[0], m_numRxnTot);
    if (j == 1) {
      plogf("vcs_solve_TP ERROR: mlequ returned an error condition\n");
      return  VCS_FAILED_CONVERGENCE;
    }

    /*
     * NOW, if we have interfacial voltage unknowns, what we did
     * was just wrong -> hopefully it didn't blow up. Redo the problem.
     * Search for inactive E
     */
    juse  = -1;
    jlose = -1;
    for (j = 0; j < m_numElemConstraints; j++) {
      if (! (m_elementActive[j])) {
	if (!strcmp((m_elementName[j]).c_str(), "E")) {
	  juse = j;
	}
      }
    }
    for (j = 0; j < m_numElemConstraints; j++) {
      if (m_elementActive[j]) {
	if (!strncmp((m_elementName[j]).c_str(), "cn_", 3)) {
	  jlose = j;
	}
      }
    }
    for (k = 0; k < m_numSpeciesTot; k++) {
      if (m_speciesUnknownType[k] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
      
	for (j = 0; j < ncTrial; ++j) {
	  for (i = 0; i < ncTrial; ++i) {
	    if (i == jlose) {
	      sm[i + j*m_numElemConstraints] = m_formulaMatrix[juse][j];
	    } else {
	      sm[i + j*m_numElemConstraints] = m_formulaMatrix[i][j];
	    }
	  }
	}
	for (i = 0; i < m_numRxnTot; ++i) {
	  k = m_indexRxnToSpecies[i];
	  for (j = 0; j < ncTrial; ++j) {
	    if (j == jlose) {
	      aw[j] = m_formulaMatrix[juse][k];
	    } else {
	      aw[j] = m_formulaMatrix[j][k];
	    }
	  }
	}  
	j = vcsUtil_mlequ(sm, m_numElemConstraints, ncTrial, aw, 1);
	if (j == 1) {
	  plogf("vcs_solve_TP ERROR: mlequ returned an error condition\n");
	  return  VCS_FAILED_CONVERGENCE;
	} 
	i = k - ncTrial;
	for (j = 0; j < ncTrial; j++) {
	  m_stoichCoeffRxnMatrix[i][j] = aw[j];
	}
      }
    }


    /*
     * Calculate the szTmp array for each formation reaction
     */
    for (i = 0; i < m_numRxnTot; i++) {
      double szTmp = 0.0;
      for (j = 0; j < ncTrial; j++) {
	szTmp += fabs(m_stoichCoeffRxnMatrix[i][j]);
      }
      m_scSize[i] = szTmp;
    }


#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   ---                Components:");
      for (j = 0; j < ncTrial; j++) {
	plogf("     %3d  ", j);
      }
      plogf("\n   ---          Components Moles:");
      for (j = 0; j < ncTrial; j++) {
	plogf("%10.3g", m_molNumSpecies_old[j]);
      }
      plogf("\n   ---   NonComponent|   Moles  |       ");
      for (j = 0; j < ncTrial; j++) {
	plogf("%-10.10s", m_speciesName[j].c_str());
      }
      //plogf("|    m_scSize");
      plogf("\n");
      for (i = 0; i < m_numRxnTot; i++) {
	plogf("   --- %3d ", m_indexRxnToSpecies[i]);
	plogf("%-10.10s", m_speciesName[m_indexRxnToSpecies[i]].c_str());
	plogf("|%10.3g|", m_molNumSpecies_old[m_indexRxnToSpecies[i]]);
	for (j = 0; j < ncTrial; j++) {
	  plogf("     %6.2f", m_stoichCoeffRxnMatrix[i][j]);
	}
	//plogf(" |  %6.2f", m_scSize[i]);
	plogf("\n");
      }
      plogf("   "); for(i=0; i<77; i++) plogf("-"); plogf("\n");
    }
#endif
    /* **************************************************** */
    /* **** EVALUATE DELTA N VALUES *********************** */
    /* **************************************************** */
    /*
     *       Evaluate the change in gas and liquid total moles 
     *       due to reaction vectors, DNG and DNL.
     */

    /*
     *  Zero out the change of Phase Moles array
     */
    vcs_dzero(m_deltaMolNumPhase[0], (NSPECIES0)*(NPHASE0));
    vcs_izero(m_phaseParticipation[0], (NSPECIES0)*(NPHASE0));
    /*
     *  Loop over each reaction, creating the change in Phase Moles
     *  array, m_deltaMolNumPhase[irxn][iphase],
     *  and the phase participation array, PhaseParticipation[irxn][iphase]
     */
    for (irxn = 0; irxn < m_numRxnTot; ++irxn) {
      scrxn_ptr = m_stoichCoeffRxnMatrix[irxn];
      dptr = m_deltaMolNumPhase[irxn];
      kspec = m_indexRxnToSpecies[irxn];
      int iph = m_phaseID[kspec];
      int *pp_ptr = m_phaseParticipation[irxn];
      dptr[iph] = 1.0;
      pp_ptr[iph]++;
      for (j = 0; j < ncTrial; ++j) {
	iph = m_phaseID[j];
	if (fabs(scrxn_ptr[j]) <= 1.0e-6) {  
	  scrxn_ptr[j] = 0.0;
	} else {  
	  dptr[iph] += scrxn_ptr[j];
	  pp_ptr[iph]++;
	}
      }
    }
    
  L_CLEANUP: ;
    double tsecond = tickTock.secondsWC();
    m_VCount->Time_basopt += tsecond;
    (m_VCount->Basis_Opts)++;
    return VCS_SUCCESS;
  } /* vcs_basopt() ************************************************************/


  int 
  VCS_SOLVE::vcs_basisOptMax(const double * const x, const int j, const int n) {
    int i;
    int largest = j;
    double big = x[j];

    assert(m_spSize[j] > 0.0);
    big *= m_spSize[j];
 
    for (i = j + 1; i < n; ++i) {
      assert(m_spSize[i] > 0.0);

      bool doSwap = false;
      if (m_SSPhase[j]) {
	doSwap = (x[i] * m_spSize[i]) > (big);
	if (!m_SSPhase[i]) {
	  if (doSwap) {
	    doSwap = (x[i]) > (x[largest]);
	  }
	}
      } else {
	if (m_SSPhase[i]) {
	  doSwap = (x[i] * m_spSize[i]) > (big);
	  if (!doSwap) {
	    doSwap = (x[i]) > (x[largest]);
	  }
	} else {
	  doSwap = (x[i] * m_spSize[i]) > (big);
	}
      }

      if (doSwap) {
	largest = i;
	big = x[i] * m_spSize[i];
      }
    } 
    
    return largest;

  }



  int VCS_SOLVE::vcs_species_type(int kspec)
   
    /*************************************************************************
     *
     * vcs_species_type:
     *
     *     Evaluate the species category for the input species
     *     return the type in the return variable
     *************************************************************************/
  {
    int irxn = kspec - m_numComponents;
    int iph, k;
   
    if (kspec < m_numComponents) return VCS_SPECIES_COMPONENT;
    if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
      return  VCS_SPECIES_INTERFACIALVOLTAGE;
    }
    iph = m_phaseID[kspec];
    if (m_molNumSpecies_old[kspec] <= 0.0) {
      if (m_deltaGRxn_new[irxn] >= 0.0) {
	/*
	 *  We are here when the species is or should be zeroed out
	 */
	if (m_SSPhase[kspec]) {
	  return VCS_SPECIES_ZEROEDSS;
	} else {
	  if (m_tPhaseMoles_old[iph] == 0.0) return VCS_SPECIES_ZEROEDPHASE;
	  else                         return VCS_SPECIES_ZEROEDMS;
	}
      }
      /*
       *     The Gibbs free energy for this species is such that
       *     it will pop back into existence.
       *     -> Set it to a major species in anticipation.
       *     -> One exception to this is if a needed component
       *        is also zeroed out. Then, don't pop the phase back into
       *        existence.
       *     -> Another exception to this is if a needed regular element
       *        is also zeroed out. Then, don't pop the phase or the species back into
       *        existence.
       */
      for (int j = 0; j < m_numComponents; ++j) {
	double stoicC = m_stoichCoeffRxnMatrix[irxn][j];
	if (stoicC != 0.0) {
	  double negChangeComp = - stoicC;
	  if (negChangeComp > 0.0) {
	    if (m_molNumSpecies_old[j] < 1.0E-60) {
#ifdef DEBUG_MODE
	      if (vcs_debug_print_lvl >= 2) {
		plogf("   ---   %s would have popped back into existance but"
		      " needed component %s is zero\n",
		      m_speciesName[kspec].c_str(),  m_speciesName[j].c_str());
	      }
#endif
	      if (m_SSPhase[kspec]) {
		return VCS_SPECIES_ZEROEDSS;
	      } else {
		return VCS_SPECIES_ZEROEDMS;
	      }
	    }
	  }
	}
      }

      for (int j = 0; j < m_numElemConstraints; ++j) {
	int elType = m_elType[j];
	if (elType == VCS_ELEM_TYPE_ABSPOS) {
	  double atomComp = m_formulaMatrix[j][kspec];
	  if (atomComp > 0.0) {
	    double maxPermissible = m_elemAbundancesGoal[j] / atomComp;
	    if (maxPermissible < VCS_DELETE_MINORSPECIES_CUTOFF) {
#ifdef DEBUG_MODE
	      if (vcs_debug_print_lvl >= 2) {
		plogf("   ---   %s would have popped back into existance but"
		      " needed element %s is zero\n",
		      m_speciesName[kspec].c_str(),  (m_elementName[j]).c_str());
	      }
#endif
	      if (m_SSPhase[kspec]) {
		return VCS_SPECIES_ZEROEDSS;
	      } else {
		return VCS_SPECIES_ZEROEDMS;
	      }
	    }
	  }
	}
      }

      return VCS_SPECIES_MAJOR;
    } 
    /*
     *   Always treat species in single species phases as majors
     */
    if (m_SSPhase[kspec]) return VCS_SPECIES_MAJOR;
    /*
     *   Check to see whether the current species is a major component
     *   of its phase. If it is, it is a major component
     */
    if (m_molNumSpecies_old[kspec] > (m_tPhaseMoles_old[iph] * 0.1)) return VCS_SPECIES_MAJOR;
    /*
     *   Main check in the loop:
     *      Check to see if there is a component with a mole number that is
     *      within a factor of 100 of the current species.
     *      If there is and that component is not part of a single species
     *      phase and shares a non-zero stoichiometric coefficient, then
     *      the current species is a major species.
     */
    double szAdj = m_scSize[irxn] * std::sqrt((double)m_numRxnTot);
    for (k = 0; k < m_numComponents; ++k) {
      if (!(m_SSPhase[k])) {
	if (m_stoichCoeffRxnMatrix[irxn][k] != 0.0) {
	  if (m_molNumSpecies_old[kspec] * szAdj >= m_molNumSpecies_old[k] * 0.01) {
	    return VCS_SPECIES_MAJOR;
	  }
	}
      }
    }  
    return VCS_SPECIES_MINOR;
  }

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  void VCS_SOLVE::vcs_chemPotPhase(int iph, const double *const molNum, 
				   double * const ac, double * const mu_i,
				   bool do_deleted)

    /**************************************************************************
     *
     * vcs_chemPotPhase:
     *
     * We calculate the dimensionless chemical potentials of all species 
     * in a single phase.
     *
     * Formula: 
     * --------------- 
     *
     *     Ideal Mixtures:
     *
     *          fe(I) = ff(I) + ln(z(I)) - ln(tPhMoles_ptr[iph])
     *
     *              ( This is equivalent to the adding the log of the 
     *                mole fraction onto the standard chemical 
     *                potential. ) 
     *
     *     Non-Ideal Mixtures:
     *        ActivityConvention = 0:
     *           fe(I) = ff(I) + ln(ActCoeff[i]z(I)) - ln(tPhMoles_ptr[iph])
     *  
     *              ( This is equivalent to the adding the log of the 
     *                mole fraction multiplied by the activity coefficient
     *                onto the standard chemical potential. ) 
     *
     *        ActivityConvention = 1: -> molality activity formulation
     *           fe(I) = ff(I) + ln(ActCoeff[i]z(I)) - ln(tPhMoles_ptr[iph])
     *                     - ln(Mnaught * m_units)
     *
     *     note:   z(I)/tPhMoles_ptr[iph] = Xmol[i] is the mole fraction
     *                                     of i in the phase.
     *
     *  NOTE:
     *  As per the discussion in vcs_dfe(), for small species where the mole
     *  fraction 
     *           z(i) < VCS_DELETE_MINORSPECIES_CUTOFF
     *   The chemical potential is calculated as:
     *         fe(I) = ff(I) + ln(ActCoeff[i](VCS_DELETE_MINORSPECIES_CUTOFF))
     *
     * Input 
     * -------- 
     *     iph    : Phase to be calculated
     *     molNum(i)   : Number of moles of species i 
     *               (VCS species order)
     *     ff     : standard state chemical potentials. These are the
     *              chemical potentials of the standard states at
     *              the same T and P as the solution.
     *                (VCS species order)
     * Output
     * -------
     *     ac[]   : Activity coefficients for species in phase
     *               (VCS species order)
     *    mu_i[]  : Dimensionless chemical potentials for phase species
     *              (VCS species order)
     * 
     *************************************************************************/
  {
    vcs_VolPhase *Vphase = m_VolPhaseList[iph];
    int nkk = Vphase->NVolSpecies;
    int k, kspec;

#ifdef DEBUG_MODE
    //if (vcs_debug_print_lvl >= 2) {
    // plogf("   --- Subroutine vcs_chemPotPhase called for phase %d\n",
    //     iph);
    //}
#endif
    double tMoles = TPhInertMoles[iph]; 
    for (k = 0; k < nkk; k++) {
      kspec = Vphase->IndSpecies[k];
      tMoles += molNum[kspec];
    }
    double tlogMoles = 0.0;
    if (tMoles > 0.0) {
      tlogMoles = log(tMoles);
    }

    Vphase->setMolesFromVCS(molNum);
    Vphase->sendToVCSActCoeff(ac);

    double phi = Vphase->electricPotential();
    double Faraday_phi = Faraday_dim * phi;

    for (k = 0; k < nkk; k++) {
      kspec = Vphase->IndSpecies[k];
      if (kspec >= m_numComponents) {
	int irxn =  kspec - m_numComponents;
	if (!do_deleted &&
	    (m_rxnStatus[irxn] == VCS_SPECIES_DELETED)) {
	  continue;
	}
      }
      if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
#ifdef DEBUG_MODE
	if (molNum[kspec] != phi) {
	  plogf("We have an inconsistency!\n");
	  exit(-1);
	}
	if (m_chargeSpecies[kspec] != -1.0) {
	  plogf("We have an unexpected situation!\n");
	  exit(-1);
	}
#endif
	mu_i[kspec] = m_SSfeSpecies[kspec] + m_chargeSpecies[kspec] * Faraday_phi;
      } else {
	if (m_SSPhase[kspec]) { 
	  mu_i[kspec] = m_SSfeSpecies[kspec] + m_chargeSpecies[kspec] * Faraday_phi;
	} else if (molNum[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
	  mu_i[kspec] = m_SSfeSpecies[kspec] + log(ac[kspec] * VCS_DELETE_MINORSPECIES_CUTOFF)
	    - tlogMoles - SpecLnMnaught[kspec] + m_chargeSpecies[kspec] * Faraday_phi;
	} else {
	  mu_i[kspec] = m_SSfeSpecies[kspec] + log(ac[kspec] * molNum[kspec])
	    - tlogMoles - SpecLnMnaught[kspec] + m_chargeSpecies[kspec] * Faraday_phi;
	}
      }
    }
  }
  /*****************************************************************************/
 
  // Calculalte the dimensionless chemical potentials of all species or
  // of certain groups of species, at a fixed temperature and pressure.
  /*
   * We calculate the dimensionless chemical potentials of all species 
   * or certain groups of species here, at a fixed temperature and pressure,
   * for the input mole vector z[] in the parameter list.
   * Nondimensionalization is achieved by division by RT.
   *
   * Note, for multispecies phases which are currently zeroed out,
   * the chemical potential is filled out with the standard chemical
   * potential.
   *
   * For species in multispecies phases whose concentration is zero,
   * we need to set the mole fraction to a very low value.
   * It's chemical potential
   * is then calculated using the VCS_DELETE_MINORSPECIES_CUTOFF concentration
   * to keep numbers positive.
   *
   *
   * Formula: 
   * --------------- 
   *
   *     Ideal Mixtures:
   *
   *          m_feSpecies(I) = m_SSfeSpecies(I) + ln(z(I)) - ln(m_tPhaseMoles[iph])
   *                            + m_chargeSpecies[I] * Faraday_dim * m_phasePhi[iphase]; 
   *
   *              ( This is equivalent to the adding the log of the 
   *                mole fraction onto the standard chemical 
   *                potential. ) 
   *
   *     Non-Ideal Mixtures:  -> molar activity formulation
   *        ActivityConvention = 0:
   *
   *          m_feSpecies(I) = m_SSfeSpecies(I)
   *                           + ln(ActCoeff[I] * z(I)) - ln(m_tPhaseMoles[iph])
   *                            + m_chargeSpecies[I] * Faraday_dim * m_phasePhi[iphase]; 
   *  
   *              ( This is equivalent to the adding the log of the 
   *                mole fraction multiplied by the activity coefficient
   *                onto the standard chemical potential. ) 
   *
   *                 note:   z(I)/tPhMoles_ptr[iph] = Xmol[i] is the mole fraction
   *                                                  of i in the phase.
   *
   *        ActivityConvention = 1: -> molality activity formulation
   *
   *          m_feSpecies(I) = m_SSfeSpecies(I)
   *                           + ln(ActCoeff[I] * z(I)) - ln(m_tPhaseMoles[iph])
   *                           - ln(Mnaught * m_units)
   *                            + m_chargeSpecies[I] * Faraday_dim * m_phasePhi[iphase]; 
   *
   *                  note:   m_SSfeSpecies(I) is the molality based standard state.
   *                          However, ActCoeff[I] is the molar based activity coefficient
   *                          We have used the formulas;
   *
   *                                ActCoeff_M[I] =  ActCoeff[I] / Xmol[N]
   *                                      where Xmol[N] is the mole fraction of the solvent 
   *                                            ActCoeff_M[I] is the molality based act coeff.
   *
   *                                m_feSpecies(I) = m_SSfeSpecies(I)
   *                                                   + ln(ActCoeff_M[I] * m(I))
   *                                                   + m_chargeSpecies[I] * Faraday_dim * m_phasePhi[iphase]
   *                                       where m[I] is the molality of the ith solute
   * 
   *                                m[I] = Xmol[I] / ( Xmol[N] * Mnaught * m_units)
   * 
   *
   *  Handling of Small Species:
   * ------------------------------
   *  As per the discussion above, for small species where the mole
   *  fraction 
   *
   *           z(i) < VCS_DELETE_MINORSPECIES_CUTOFF
   *
   *   The chemical potential is calculated as:
   *
   *         m_feSpecies(I) = m_SSfeSpecies(I) + ln(ActCoeff[i](VCS_DELETE_MINORSPECIES_CUTOFF))
   *
   *   Handling of "Species" Representing Interfacial Voltages
   *  ---------------------------------------------------------
   *
   *      These species have species types of  VCS_SPECIES_TYPE_INTERFACIALVOLTAGE
   *      The chemical potentials for these "species" refer to electrons in
   *      metal electrodes. They have the following formula
   *
   *          m_feSpecies(I) = m_SSfeSpecies(I) - F z[I] / RT
   *
   *      F is Faraday's constant.
   *      R = gas constant
   *      T = temperature
   *      V = potential of the interface = phi_electrode - phi_solution
   *
   *      For these species, the solution vector unknown, z[I], is V, the phase voltage,
   *      in volts.    
   *
   * Input 
   * -------- 
   * @param ll     Determine which group of species gets updated
   *     ll =  0: Calculate for all species 
   *        <  0: calculate for components and for major non-components 
   *           1: calculate for components and for minor non-components 
   *
   * @param lbot    Restricts the calculation of the chemical potential 
   *                to the species between LBOT <= i < LTOP. Usually 
   *                 LBOT and LTOP will be equal to 0 and MR, respectively.
   * @param ltop    Top value of the loops
   *  
   * @param z   z[i]   : Number of moles of species i 
   *                   -> This can either be the current solution vector WT() 
   *                      or the actual solution vector W() 
   *
   * @param kk   Determines whether z is old or new or tmp:
   *            VCS_STATECALC_NEW: Use the tentative values for the total number of 
   *                               moles in the phases, i.e., use TG1 instead of TG etc. 
   *            VCS_STATECALC_OLD: Use the base values of the total number of 
   *                               moles in each system. 
   *
   *  Also needed:
   *     ff     : standard state chemical potentials. These are the
   *              chemical potentials of the standard states at
   *              the same T and P as the solution.
   *     tg     : Total Number of moles in the phase.
   */
  void VCS_SOLVE::vcs_dfe(double const * const z, int kk, int ll, int lbot, int ltop) {
    int l1, l2, iph, kspec, irxn;
    int iphase;
    double *tPhMoles_ptr;
    double *actCoeff_ptr;
    double *tlogMoles;
    vcs_VolPhase *Vphase;
    VCS_SPECIES_THERMO *st_ptr;

#ifdef DEBUG_MODE
    if (kk != VCS_STATECALC_OLD && kk != VCS_STATECALC_NEW) {
      plogf("   --- Subroutine vcs_dfe called with bad kk value: %d", kk);
      plogendl();
      exit(-1);
    }
#endif

#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      if (ll == 0) {
	if (lbot != 0) {
	  plogf("   --- Subroutine vcs_dfe called for one species: ");
	  plogf("%-12.12s", m_speciesName[lbot].c_str());
	} else {
	  plogf("   --- Subroutine vcs_dfe called for all species");
	}
      } else if (ll > 0) {
	plogf("   --- Subroutine vcs_dfe called for components and minors");
      } else {
	plogf("   --- Subroutine vcs_dfe called for components and majors");
      }
      if (kk == VCS_STATECALC_NEW)  plogf(" using tentative solution");
      plogendl();
    }
#endif
    if (kk <= VCS_STATECALC_OLD) {
      tPhMoles_ptr = VCS_DATA_PTR(m_tPhaseMoles_old);
      actCoeff_ptr = VCS_DATA_PTR(m_actCoeffSpecies_old);
    } else {
      tPhMoles_ptr = VCS_DATA_PTR(m_tPhaseMoles_new);
      actCoeff_ptr = VCS_DATA_PTR(m_actCoeffSpecies_new);
    }
    tlogMoles = VCS_DATA_PTR(TmpPhase);
    /*
     * Might as well recalculate the phase mole vector
     * and compare to the storred one. They should be correct.
     */
    double *tPhInertMoles = VCS_DATA_PTR(TPhInertMoles);
    for (iph = 0; iph < m_numPhases; iph++) {
      tlogMoles[iph] = tPhInertMoles[iph];
 
    }
    for (kspec = 0; kspec < m_numSpeciesTot; kspec++) {
      if(m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	iph = m_phaseID[kspec];
	tlogMoles[iph] += z[kspec];
      }
    }
#ifdef DEBUG_MODE
    for  (iph = 0; iph < m_numPhases; iph++) {
      if (! vcs_doubleEqual(tlogMoles[iph], tPhMoles_ptr[iph])) {
	plogf("phase Moles may be off, iph = %d, %20.14g %20.14g \n",
	      iph, tlogMoles[iph], tPhMoles_ptr[iph]);
	exit(0);
      }
    }
#endif
    vcs_dzero(tlogMoles, m_numPhases);
    for (iph = 0; iph < m_numPhases; iph++) {
      if (tPhMoles_ptr[iph] > 0.0) {
	tlogMoles[iph] = log(tPhMoles_ptr[iph]);
      }
    }
    /*
     *    Zero the indicator that that tells us the activity coefficients
     *    are current
     */
    vcs_izero(VCS_DATA_PTR(CurrPhAC), m_numPhases);
    
    if (ll != 0) {
      l1 = lbot;
      l2 = m_numComponents;
    } else {
      l1 = lbot;
      l2 = ltop;
    }

    /*
     *  Calculate activity coefficients for all phases that are
     *  not current
     */
    for (iphase = 0; iphase < m_numPhases; iphase++) {
      if (!CurrPhAC[iphase]) {
	Vphase = m_VolPhaseList[iphase];
	if (!Vphase->SingleSpecies) {
	  Vphase->setMolesFromVCS(z);
	  Vphase->sendToVCSActCoeff(VCS_DATA_PTR(actCoeff_ptr));
	}
	m_phasePhi[iphase] = Vphase->electricPotential();
	CurrPhAC[iphase] = 1;
      }
    }
    /* ************************************************************** */
    /* **** ALL SPECIES, OR COMPONENTS ****************************** */
    /* ************************************************************** */
    /*
     *     Do all of the species when LL = 0. Then we are done for the routine 
     *     When LL ne 0., just do the initial components. We will then 
     *     finish up below with loops over either the major noncomponent 
     *     species or the minor noncomponent species.
     */
    for (kspec = l1; kspec < l2; ++kspec) {
      iphase = m_phaseID[kspec];
      if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
#ifdef DEBUG_MODE
	if (z[kspec] != m_phasePhi[iphase]) {
	  plogf("We have an inconsistency!\n");
	  exit(-1);
	}
	if (m_chargeSpecies[kspec] != -1.0) {
	  plogf("We have an unexpected situation!\n");
	  exit(-1);
	}
#endif
	m_feSpecies_curr[kspec] =
	  m_SSfeSpecies[kspec] + m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iphase];
      } else {
	if (m_SSPhase[kspec]) {
	  m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec];
	} else {
	  if (z[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
	    iph = m_phaseID[kspec];
	    if (tPhMoles_ptr[iph] > 0.0) { 
	      m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec] 
		+ log(actCoeff_ptr[kspec] * VCS_DELETE_MINORSPECIES_CUTOFF)
		- tlogMoles[m_phaseID[kspec]] - SpecLnMnaught[kspec] 
		+ m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iphase];
	    } else {
	      m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec];
	    }
	  } else {
	    m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec] + log(actCoeff_ptr[kspec] * z[kspec])
	      - tlogMoles[m_phaseID[kspec]] - SpecLnMnaught[kspec] 
	      + m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iphase]; 
	  }
	}
      }
    }
    /* ************************************************ */
    /* **** MAJORS ONLY ******************************* */
    /* ************************************************ */
    if (ll < 0) {
      for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
	if (m_rxnStatus[irxn] != VCS_SPECIES_MINOR) {
	  kspec = m_indexRxnToSpecies[irxn];
	  iphase = m_phaseID[kspec];
	  if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
#ifdef DEBUG_MODE
	    if (z[kspec] != m_phasePhi[iphase]) {
	      plogf("We have an inconsistency!\n");
	      exit(-1);
	    }
	    if (m_chargeSpecies[kspec] != -1.0) {
	      plogf("We have an unexpected situation!\n");
	      exit(-1);
	    }
#endif
	    m_feSpecies_curr[kspec] = 
	      m_SSfeSpecies[kspec] + m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iphase];   
	  } else {
	    if (m_SSPhase[kspec]) {
	      m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec];
	    } else {
	      if (z[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
		iph = m_phaseID[kspec];
		if (tPhMoles_ptr[iph] > 0.0) { 
		  m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec] 
		    + log(actCoeff_ptr[kspec] * VCS_DELETE_MINORSPECIES_CUTOFF)
		    - tlogMoles[m_phaseID[kspec]] - SpecLnMnaught[kspec]
		    + m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iphase]; ;
		} else {
		  m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec];
		}
	      } else {
		m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec] + log(actCoeff_ptr[kspec] * z[kspec]) 
		  - tlogMoles[m_phaseID[kspec]] - SpecLnMnaught[kspec] 
		  + m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iphase]; 
	      }
	    }
	  }
	}
      }
      /* ************************************************ */
      /* **** MINORS ONLY ******************************* */
      /* ************************************************ */
    } else if (ll > 0) {
      for (irxn = 0; irxn < m_numRxnRdc; ++irxn) {
	if (m_rxnStatus[irxn] == VCS_SPECIES_MINOR) {
	  kspec = m_indexRxnToSpecies[irxn];
	  iphase = m_phaseID[kspec];
	  if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
#ifdef DEBUG_MODE
	    if (z[kspec] != m_phasePhi[iphase]) {
	      plogf("We have an inconsistency!\n");
	      exit(-1);
	    }
	    if (m_chargeSpecies[kspec] != -1.0) {
	      plogf("We have an unexpected situation!\n");
	      exit(-1);
	    }
#endif
	    m_feSpecies_curr[kspec] = 
	      m_SSfeSpecies[kspec] + m_chargeSpecies[kspec] * Faraday_dim * m_phasePhi[iphase]; ;
	  } else {
	    if (m_SSPhase[kspec]) {
	      m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec];
	    } else {
	      if (z[kspec] <= VCS_DELETE_MINORSPECIES_CUTOFF) {
		iph = m_phaseID[kspec];
		if (tPhMoles_ptr[iph] > 0.0) { 
		  m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec]
		    + log(actCoeff_ptr[kspec] * VCS_DELETE_MINORSPECIES_CUTOFF)
		    - tlogMoles[m_phaseID[kspec]] - SpecLnMnaught[kspec];
		} else {
		  m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec];
		}
	      } else {
		st_ptr = m_speciesThermoList[kspec];
		m_feSpecies_curr[kspec] = m_SSfeSpecies[kspec] + log(actCoeff_ptr[kspec] * z[kspec]) 
		  - tlogMoles[m_phaseID[kspec]] - SpecLnMnaught[kspec]; 
	      }
	    }
	  }
	}
      }
    }
#ifdef DEBUG_NOT
    for (kspec = 0; kspec < m_numSpeciesRdc; kspec++) {
      checkFinite(fe[kspec]);
    } 
#endif
  }
  /*****************************************************************************/

#ifdef DEBUG_MODE
  //! Print out and check the elemental abundance vector
  void VCS_SOLVE::prneav() const {
    int kerr, j;
    std::vector<double> eav(m_numElemConstraints, 0.0);

    for (j = 0; j < m_numElemConstraints; ++j) {
      for (int i = 0; i < m_numSpeciesTot; ++i) {
	if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  eav[j] += m_formulaMatrix[j][i] * m_molNumSpecies_old[i];
	}
      }
    }
    kerr = FALSE;
    plogf( "--------------------------------------------------");
    plogf("ELEMENT ABUNDANCE VECTOR:\n");
    plogf(" Element   Now      Orignal      Deviation          Type\n");
    for (j = 0; j < m_numElemConstraints; ++j) {
      plogf("  ");  plogf("%-2.2s", (m_elementName[j]).c_str());
      plogf(" = %15.6E     %15.6E     %15.6E  %3d\n",
	    eav[j], m_elemAbundancesGoal[j], eav[j] - m_elemAbundancesGoal[j], m_elType[j]);
      if (m_elemAbundancesGoal[j] != 0.) {
	if (fabs(eav[j] - m_elemAbundancesGoal[j]) > m_elemAbundancesGoal[j] * 5.0e-9)
	  kerr = TRUE;
      } else {
	if (fabs(eav[j]) > 1.0e-10) kerr = TRUE;
      }
    }
    if (kerr) {
      plogf("Element abundance check failure\n");
    }
    plogf("--------------------------------------------------");
    plogendl();
  } 
#endif
  /*****************************************************************************/

  // Calculate the norm of a deltaGibbs free energy vector
  /*
   *   Positive DG for species which don't exist are ignored. 
   *
   * @param dgLocal  Vector of local delta G's.
   */
  double VCS_SOLVE::l2normdg(double dgLocal[]) const {
    double tmp;
    int irxn;
    if (m_numRxnRdc <= 0) return 0.0;
    for (irxn = 0, tmp = 0.0; irxn < m_numRxnRdc; ++irxn) {
      if (m_rxnStatus[irxn] == VCS_SPECIES_MAJOR || m_rxnStatus[irxn] == VCS_SPECIES_MINOR ||
	  dgLocal[irxn] < 0.0) {
	if (m_rxnStatus[irxn] != VCS_SPECIES_ZEROEDMS) {
	  tmp += dgLocal[irxn] * dgLocal[irxn];
	}
      }
    }
    return (std::sqrt(tmp / m_numRxnRdc));
  }
  /*****************************************************************************/

  //  Calculates the total number of moles of species in all phases.
  /*
   *  Calculates the total number of moles in all phases.
   *  Reconciles Phase existence flags with total moles in each phase.
   */
  void VCS_SOLVE::vcs_tmoles() {
    int i;
    double sum;
    vcs_VolPhase *Vphase;
    for (i = 0; i < m_numPhases; i++) {
      m_tPhaseMoles_old[i] = TPhInertMoles[i];
    }
    for (i = 0; i < m_numSpeciesTot; i++) {
      if (m_speciesUnknownType[i] == VCS_SPECIES_TYPE_MOLNUM) {
	m_tPhaseMoles_old[m_phaseID[i]] += m_molNumSpecies_old[i];
      }
    }
    sum = 0.0;
    for (i = 0; i < m_numPhases; i++) {
      sum += m_tPhaseMoles_old[i];
      Vphase = m_VolPhaseList[i];
      // Took out because we aren't updating mole fractions in Vphase
      // Vphase->TMoles = m_tPhaseMoles_old[i];
      if (m_tPhaseMoles_old[i] == 0.0) {
	Vphase->Existence = 0;
      } else {
	if (TPhInertMoles[i] > 0.0) {
	  Vphase->Existence = 2;
	} else {
	  Vphase->Existence = 1;
	}
      }
    }  
    m_totalMolNum = sum;
  }
  /*****************************************************************************/

  //  This routine uploads the state of the system into all of the 
  //  vcs_VolPhase objects in the current problem.
  /*
   *  @param vcsState Determines where to get the mole numbers from.
   *                -  VCS_STATECALC_OLD -> from m_molNumSpecies_old
   *                -  VCS_STATECALC_NEW -> from m_molNumSpecies_new
   *
   */
  void VCS_SOLVE::vcs_updateVP(const int vcsState) {
    vcs_VolPhase *Vphase;
    for (int i = 0; i < m_numPhases; i++) {
      Vphase = m_VolPhaseList[i];
      if (vcsState == VCS_STATECALC_OLD) {
	Vphase->setMolesFromVCSCheck(VCS_DATA_PTR(m_molNumSpecies_old),
				     VCS_DATA_PTR(m_tPhaseMoles_old), i);
      } else if (vcsState == VCS_STATECALC_NEW) {
	Vphase->setMolesFromVCSCheck(VCS_DATA_PTR(m_molNumSpecies_new),
				     VCS_DATA_PTR(m_tPhaseMoles_new), i);
      }
#ifdef DEBUG_MODE
      else {
	plogf("we shouldn't be here");
	plogendl();
	std::exit(-1);
      }
#endif
    }
  }

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  void VCS_SOLVE::vcs_switch2D(double * const * const Jac, int k1, int k2)
   
    /**************************************************************************
     *  vcs_switch2D:
     *
     *        Switch rows and columns of a square matrix
     *************************************************************************/
  {
    int i;
    register double dtmp;
    for (i = 0; i < m_numSpeciesTot; i++) {
      SWAP(Jac[k1][i], Jac[k2][i], dtmp);  
    }
    for (i = 0; i < m_numSpeciesTot; i++) {
      SWAP(Jac[i][k1], Jac[i][k2], dtmp);  
    }
  }
  /*****************************************************************************/

  //  Swaps the indecises for all of the global data for two species, k1
  //  and k2.
  /*
   *  @param  ifunc:  If true, switch the species data and the noncomponent reaction
   *                  data. This must be called for a non-component species only.
   *                  If false, switch the species data only. Typically, we use this
   *                  option when determining the component species and at the
   *                  end of the calculation, when we want to return unscrambled 
   *                  results. All rxn data will be out-of-date.
   *
   *  @param k1        First species index
   *
   *  @param k2        Second species index
   */
  void VCS_SOLVE::vcs_switch_pos(const int ifunc, const int k1, const int k2) {
    register int j;
    register double t1 = 0.0;
    int i1, i2, iph, kp1, kp2;
    vcs_VolPhase *pv1, *pv2;
    VCS_SPECIES_THERMO *st_tmp;
    if (k1 == k2) return;
#ifdef DEBUG_MODE
    if (k1 < 0 || k1 > (m_numSpeciesTot - 1) ||
	k2 < 0 || k2 > (m_numSpeciesTot - 1)    ) {
      plogf("vcs_switch_pos: ifunc = 0: inappropriate args: %d %d\n",
	    k1, k2);
    }
#endif
    /*
     *     Handle the index pointer in the phase structures first
     */
    pv1 = m_VolPhaseList[m_phaseID[k1]];
    pv2 = m_VolPhaseList[m_phaseID[k2]];

    kp1 = m_speciesLocalPhaseIndex[k1];
    kp2 = m_speciesLocalPhaseIndex[k2];
#ifdef DEBUG_MODE
    if (pv1->IndSpecies[kp1] != k1) {
      plogf("Indexing error in program\n");
      exit(-1);
    }
    if (pv2->IndSpecies[kp2] != k2) {
      plogf("Indexing error in program\n");
      exit(-1);
    }
#endif
    pv1->IndSpecies[kp1] = k2;
    pv2->IndSpecies[kp2] = k1;
   
    vcsUtil_stsw(m_speciesName,  k1, k2);
    SWAP(m_molNumSpecies_old[k1], m_molNumSpecies_old[k2], t1);
    SWAP(m_speciesUnknownType[k1], m_speciesUnknownType[k2], j);
    SWAP(m_molNumSpecies_new[k1], m_molNumSpecies_new[k2], t1);
    SWAP(m_SSfeSpecies[k1], m_SSfeSpecies[k2], t1);
    SWAP(m_spSize[k1], m_spSize[k2], t1);
    SWAP(m_feSpecies_curr[k1], m_feSpecies_curr[k2], t1);
    SWAP(m_deltaMolNumSpecies[k1], m_deltaMolNumSpecies[k2], t1);
    SWAP(m_feSpecies_old[k1], m_feSpecies_old[k2], t1);
    SWAP(m_feSpecies_new[k1], m_feSpecies_new[k2], t1);
    SWAP(m_SSPhase[k1], m_SSPhase[k2], j);
    SWAP(m_phaseID[k1], m_phaseID[k2], j);
    SWAP(m_speciesMapIndex[k1], m_speciesMapIndex[k2], j); 
    SWAP(m_speciesLocalPhaseIndex[k1], m_speciesLocalPhaseIndex[k2], j);
    SWAP(SpecActConvention[k1], SpecActConvention[k2], j);
    SWAP(SpecLnMnaught[k1], SpecLnMnaught[k2], t1);
    SWAP(m_actCoeffSpecies_new[k1], m_actCoeffSpecies_new[k2], t1);
    SWAP(m_actCoeffSpecies_old[k1], m_actCoeffSpecies_old[k2], t1);
    SWAP(m_wtSpecies[k1], m_wtSpecies[k2], t1);
    SWAP(m_chargeSpecies[k1], m_chargeSpecies[k2], t1);
    SWAP(m_speciesThermoList[k1], m_speciesThermoList[k2], st_tmp);
    SWAP(m_PMVolumeSpecies[k1], m_PMVolumeSpecies[k2], t1);

    for (j = 0; j < m_numElemConstraints; ++j) {
      SWAP(m_formulaMatrix[j][k1], m_formulaMatrix[j][k2], t1);
    }   
    if (m_useActCoeffJac) {
      vcs_switch2D(m_dLnActCoeffdMolNum.baseDataAddr(), k1, k2);
    }
   
    /*
     *     Handle the index pointer in the phase structures
     */

   
    if (ifunc) {
      /*
       * Find the Rxn indecises corresponding to the two species
       */
      i1 = k1 - m_numComponents;
      i2 = k2 - m_numComponents;
#ifdef DEBUG_MODE
      if (i1 < 0 || i1 > (m_numRxnTot - 1) ||
	  i2 < 0 || i2 > (m_numRxnTot - 1)    ) {
	plogf("switch_pos: ifunc = 1: inappropriate noncomp values: %d %d\n",
	      i1 , i2);
      }
#endif
      for (j = 0; j < m_numComponents; ++j) {
	SWAP(m_stoichCoeffRxnMatrix[i1][j], m_stoichCoeffRxnMatrix[i2][j], t1);
      }
      SWAP(m_scSize[i1], m_scSize[i2], t1);
      for (iph = 0; iph < m_numPhases; iph++) {
	SWAP(m_deltaMolNumPhase[i1][iph], m_deltaMolNumPhase[i2][iph], t1);
	SWAP(m_phaseParticipation[i1][iph], 
	     m_phaseParticipation[i2][iph], j);
      }
      SWAP(m_deltaGRxn_new[i1],  m_deltaGRxn_new[i2],  t1);
      SWAP(m_deltaGRxn_old[i1], m_deltaGRxn_old[i2], t1);
      SWAP(m_deltaGRxn_tmp[i1], m_deltaGRxn_tmp[i2], t1);
      SWAP(m_rxnStatus[i1],  m_rxnStatus[i2],  j);

      /*
       *   We don't want to swap ir[], because the values of ir should 
       *   stay the same after the swap
       *
       * vcs_isw(ir, i1, i2); 
       */
    }
  } /* vcs_switch_pos() ********************************************************/

  static void print_space(int num)
  {
    int j;
    for (j = 0; j < num; j++) plogf(" ");
  }

  /****************************************************************************
   *
   *  vcs_deltag_Phase():
   *
   *     Calculate deltag of formation for all species in a single
   *     phase. It is assumed that the fe[] is up to date for all species.
   *     Howevever, if the phase is currently zereoed out, a subproblem
   *     is calculated to solve for AC[i] and pseudo-X[i] for that 
   *     phase. 
   */
  void VCS_SOLVE::vcs_deltag_Phase(int iphase, bool doDeleted) {
    int iph;
    int  irxn, kspec, kcomp;
    double *dtmp_ptr;
    int irxnl = m_numRxnRdc;
    if (doDeleted) irxnl = m_numRxnTot;
    vcs_VolPhase *vPhase = m_VolPhaseList[iphase];

#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("   --- Subroutine vcs_deltag_Phase called for phase %d\n",
	    iphase);
    }
#endif

    /*
     * Single species Phase
     */
    if (vPhase->SingleSpecies) {
      kspec = vPhase->IndSpecies[0];
#ifdef DEBUG_MODE
      if (iphase != m_phaseID[kspec]) {
	plogf("vcs_deltag_Phase index error\n");
	exit(-1);
      }
#endif
      if (kspec >= m_numComponents) {
	irxn = kspec - m_numComponents;
	m_deltaGRxn_new[irxn] = m_feSpecies_curr[kspec];
	dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
	for (kcomp = 0; kcomp < m_numComponents; ++kcomp) {
	  m_deltaGRxn_new[irxn] += dtmp_ptr[kcomp] * m_feSpecies_curr[kcomp];
	}
      }
    } 
    /*
     * Multispecies Phase
     */
    else {
      bool zeroedPhase = TRUE;

      for (irxn = 0; irxn < irxnl; ++irxn) {
	kspec = m_indexRxnToSpecies[irxn];
	if (m_speciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  iph = m_phaseID[kspec];
	  if (iph == iphase ) {
	    if (m_molNumSpecies_old[kspec] > 0.0) zeroedPhase = FALSE;
	    m_deltaGRxn_new[irxn] = m_feSpecies_curr[kspec];
	    dtmp_ptr = m_stoichCoeffRxnMatrix[irxn];
	    for (kcomp = 0; kcomp < m_numComponents; ++kcomp) {
	      m_deltaGRxn_new[irxn] += dtmp_ptr[kcomp] * m_feSpecies_curr[kcomp];
	    }
	  }
	}
      }
  
      /*
       * special section for zeroed phases
       */
      /* ************************************************* */
      /* **** MULTISPECIES PHASES WITH ZERO MOLES************ */
      /* ************************************************* */
      /*
       *    Massage the free energies for species with zero mole fractions 
       *  in multispecies phases.  This section implements the 
       *  Equation 3.8-5 in Smith and Missen, p.59. 
       *  A multispecies phase will exist iff 
       *           1 < sum_i(exp(-dg_i)/AC_i) 
       *  If DG is negative then that species wants to be reintroduced into 
       *  the calculation. 
       *  For small dg_i, the expression below becomes: 
       *      1 - sum_i(exp(-dg_i)/AC_i) ~ sum_i((dg_i-1)/AC_i)  + 1
       * 
       *
       *  HKM -> The ratio of mole fractions at the reinstatement
       *         time should be equal to the normalized weighting
       *         of exp(-dg_i) / AC_i. This should be implemented.
       *
       *  HKM -> There is circular logic here. ActCoeff depends on the
       *         mole fractions of a phase that does not exist. In actuality
       *         the proto-mole fractions should be selected from the
       *         solution of a nonlinear problem with NsPhase unknowns
       *
       *              X_i = exp(-dg[irxn]) / ActCoeff_i / denom
       *
       *              where 
       *               denom = sum_i[  exp(-dg[irxn]) / ActCoeff_i  ]
       *      
       *         This can probably be solved by successive iteration.
       *         This should be implemented.
       */
      /*
       *      Calculate dg[] for each species in a zeroed multispecies phase.
       *      All of the dg[]'s will be equal. If dg[] is negative, then
       *      the phase will come back into existence.
       */
      if (zeroedPhase) {
	double phaseDG = 1.0;
	for (irxn = 0; irxn < irxnl; ++irxn) {
	  kspec = m_indexRxnToSpecies[irxn];
	  iph = m_phaseID[kspec];
	  if (iph == iphase) {
	    if (m_deltaGRxn_new[irxn] >  50.0) m_deltaGRxn_new[irxn] =  50.0;
	    if (m_deltaGRxn_new[irxn] < -50.0) m_deltaGRxn_new[irxn] = -50.0;
	    phaseDG -= exp(-m_deltaGRxn_new[irxn])/m_actCoeffSpecies_new[kspec];
	  }
	}
	/*
	 * Overwrite the individual dg's with the phase DG.
	 */
	for (irxn = 0; irxn < irxnl; ++irxn) {
	  kspec = m_indexRxnToSpecies[irxn];
	  iph = m_phaseID[kspec];
	  if (iph == iphase) {
	    m_deltaGRxn_new[irxn] = 1.0 - phaseDG;
	  }
	}
      }
    }

  }

  /****************************************************************************
   *
   *  vcs_birthGuess
   *
   *     Birth guess returns the number of moles of a species 
   *     that is coming back to life. or -> whose concentration has
   *     been forced to zero by a constraint for some reason, and needs
   *     to be reinitialized.
   */
  double VCS_SOLVE::vcs_birthGuess(int kspec) {
    int irxn = kspec - m_numComponents;
    int soldel = false;
    double dx = 0.0;
    if (m_speciesUnknownType[kspec] == VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
      return dx;
    }
    double w_kspec = VCS_DELETE_SPECIES_CUTOFF;
    // Check to make sure that species is zero in the solution vector
    // If it isn't, we don't know what's happening
    if (m_molNumSpecies_old[kspec] != 0.0) {
      w_kspec = 0.0;
      plogf("we shouldn't be here\n");
      exit(-1);
    }
    int ss = m_SSPhase[kspec];
    if (!ss) {
      /*
       * Logic to handle species in multiple species phases
       */
#ifdef DEBUG_MODE	 
      char ANOTE[32];
      double dxm = minor_alt_calc(kspec, irxn, &soldel, ANOTE); 
#else
      double dxm = minor_alt_calc(kspec, irxn, &soldel);
#endif
      dx = w_kspec + dxm;
      if (dx > 1.0E-15) {
	dx = 1.0E-15;
      }
    } else {
      /*
       * Logic to handle single species phases
       */
      dx = VCS_DELETE_SPECIES_CUTOFF * 100.;
    }
    
    /*
     * Check to see if the current value of the components 
     * allow the dx.
     * If we are in danger of zeroing a component,
     * only go 1/3 the way to zeroing the component with
     * this dx. Note, this may mean that dx= 0 coming
     * back from this routine. This evaluation should
     * be respected.
     */ 
    double *sc_irxn = m_stoichCoeffRxnMatrix[irxn];
    for (int j = 0; j < m_numComponents; ++j) {
      // Only loop over element contraints that involve positive def. constraints
      if (m_speciesUnknownType[j] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	if (m_molNumSpecies_old[j] > 0.0) {
	  double tmp = sc_irxn[j] * dx;
	  if (3.0*(-tmp) > m_molNumSpecies_old[j]) {
	    dx = MIN(dx, - 0.3333* m_molNumSpecies_old[j] / sc_irxn[j]); 
	  }
	}
	if (m_molNumSpecies_old[j] <= 0.0) {
	  if (sc_irxn[j] < 0.0) {
	    dx = 0.0;
	  }
	}
      }
    }
    return dx;
  }
  /*****************************************************************/
}

