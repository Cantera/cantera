/**
 * @file vcs_inest.cpp
 *   Implementation methods for obtaining a good initial guess
 */
/*  $Author$
 *  $Date$
 *  $Revision$
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
  
#include "vcs_solve.h"
#include "vcs_internal.h" 
#include "vcs_VolPhase.h"

#include "clockWC.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

namespace VCSnonideal {

  static char pprefix[20] = "   --- vcs_inest: ";
 
  // Estimate equilibrium compositions
  /*
   *  Estimates equilibrium compositions.
   *  Algorithm covered in a section of Smith and Missen's Book.
   *
   *  Linear programming module is based on using dbolm.
   *
   *  @param aw   aw[i[  Mole fraction work space        (ne in length)
   *  @param sa   sa[j] = Gramm-Schmidt orthog work space (ne in length)
   *  @param sm   sm[i+j*ne] = QR matrix work space (ne*ne in length)
   *  @param ss   ss[j] = Gramm-Schmidt orthog work space (ne in length)
   *  @param test This is a small negative number.
   */
  void VCS_SOLVE::inest(double * const aw, double * const sa, double * const sm, 
			double * const ss, double test) {
    int conv, k,  lt, ikl, kspec, iph, irxn;
    double s;
    double s1 = 0.0;
    double xl, par;
    int finished;
    int     nspecies   = m_numSpeciesTot;
    int     nrxn       = m_numRxnTot;
    vcs_VolPhase *Vphase = 0;
 
    double *molNum   = VCS_DATA_PTR(m_molNumSpecies_old);
    double TMolesMultiphase;
    double *xtphMax = VCS_DATA_PTR(TmpPhase);
    double *xtphMin = VCS_DATA_PTR(TmpPhase2);
  
    ikl = 0;
    lt = 0;
   

    /*
     *       CALL ROUTINE TO SOLVE MAX(CC*molNum) SUCH THAT AX*molNum = BB 
     *           AND molNum(I) .GE. 0.0
     *
     *   Note, both of these programs do this.
     */
#ifdef ALTLINPROG
    vcs_setMolesLinProg();
#else
    int j, jj;
    std::vector<double> ax(m_numElemConstraints*nspecies, 0.0);
    std::vector<double> bb(m_numElemConstraints, 0.0);
    std::vector<double> cc(nspecies, 0.0);

    int neActive = 0;
    jj = 0;
    for (j = 0; j < m_numElemConstraints; j++) {
      if (ElActive[j]) {
	neActive++;
	bb[jj] = gai[j];
	jj++;
      }
    }
    for (kspec = 0; kspec < nspecies; ++kspec) {
      cc[kspec] = -ff[kspec];
      jj = 0;
      for (j = 0; j < m_numElemConstraints; ++j) {
	if (ElActive[j]) {
	  ax[jj + kspec * neActive] = FormulaMatrix[j][kspec];
	  jj++;
	}
      }
    }
    linprogmax(molNum, VCS_DATA_PTR(cc), VCS_DATA_PTR(ax),
	       VCS_DATA_PTR(bb), neActive, nspecies, neActive);
#endif

#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("%s Mole Numbers returned from linear programming (vcs_inest initial guess):\n",
	     pprefix);
      plogf("%s     SPECIES          MOLE_NUMBER      -SS_ChemPotential\n", pprefix);
      for (kspec = 0; kspec < nspecies; ++kspec) {
	plogf("%s     ", pprefix); plogf("%-12.12s", SpName[kspec].c_str());
	plogf(" %15.5g  %12.3g\n", molNum[kspec], -m_SSfeSpecies[kspec]);
      }
      plogf("%s Element Abundance Agreement returned from linear "
	     "programming (vcs_inest initial guess):",
	     pprefix);
      plogendl();
      plogf("%s     Element           Goal         Actual\n", pprefix);
      int jj = 0;
      for (int j = 0; j < m_numElemConstraints; j++) {
	if (ElActive[j]) {
	  double tmp = 0.0;
	  for (kspec = 0; kspec < nspecies; ++kspec) {
	    tmp +=  FormulaMatrix[j][kspec] * molNum[kspec];
	  }
	  plogf("%s     ", pprefix); plogf("   %-9.9s", (ElName[j]).c_str());
	  plogf(" %12.3g %12.3g\n", m_elemAbundancesGoal[j], tmp);
	  jj++;
	}
      }
      plogendl();
    }
#endif
   

    /*
     *     Make sure all species have positive definite mole numbers
     *     Set voltages to zero for now, until we figure out what to do
     */
    vcs_dzero(VCS_DATA_PTR(m_deltaMolNumSpecies), nspecies);
    for (kspec = 0; kspec < nspecies; ++kspec) {
      iph = PhaseID[kspec];
      Vphase = VPhaseList[iph];
      if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	if (molNum[kspec] <= 0.0) {
	  /*
	   * HKM Should eventually include logic here for non SS phases
	   */
	  if (!SSPhase[kspec]) {
	    molNum[kspec] = 1.0e-30;
	  }
	}
      } else {
	molNum[kspec] = 0.0;
      }
      if (molNum[kspec] > 0.0) {
	if (Vphase->Existence == 0) {
	  Vphase->Existence = 1;
	}
      } else if (SSPhase[kspec]) {
	Vphase->Existence = 0;
      }
    }
   
    /*
     *      Now find the optimized basis that spans the stoichiometric
     *      coefficient matrix
     */
    (void) vcs_basopt(FALSE, aw, sa, sm, ss, test, &conv);

    /* ***************************************************************** */
    /* **** CALCULATE TOTAL GASEOUS AND LIQUID MOLES, ****************** */
    /* **** CHEMICAL POTENTIALS OF BASIS              ****************** */
    /* ***************************************************************** */
    /*
     * Calculate TMoles and m_tPhaseMoles_old[]
     */
    vcs_tmoles();
    /*
     * m_tPhaseMoles_new[] will consist of just the component moles
     */
    for (iph = 0; iph < NPhase; iph++) {
      m_tPhaseMoles_new[iph] = TPhInertMoles[iph] + 1.0E-20;
    }
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
      if (SpeciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
	m_tPhaseMoles_new[PhaseID[kspec]] += molNum[kspec];
      }
    }
    TMolesMultiphase = 0.0;
    for (iph = 0; iph < NPhase; iph++) {
      if (! VPhaseList[iph]->SingleSpecies) {
	TMolesMultiphase += m_tPhaseMoles_new[iph];
      }
    }     
    vcs_dcopy(VCS_DATA_PTR(m_molNumSpecies_new), molNum,  nspecies);
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
      if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_MOLNUM) {
	m_molNumSpecies_new[kspec] = 0.0;
      }
    }
    vcs_dcopy(VCS_DATA_PTR(m_feSpecies_curr), VCS_DATA_PTR(m_SSfeSpecies),
	      nspecies);
 
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
      if (SpeciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
	if (! SSPhase[kspec]) {
	  iph = PhaseID[kspec];
	  m_feSpecies_curr[kspec] += log(m_molNumSpecies_new[kspec] / m_tPhaseMoles_old[iph]);
	}
      } else {
	m_molNumSpecies_new[kspec] = 0.0;
      }
    }
    vcs_deltag(0, true);
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      for (kspec = 0; kspec < nspecies; ++kspec) {
	plogf("%s", pprefix); plogf("%-12.12s", SpName[kspec].c_str());
	if (kspec < m_numComponents)
	  plogf("fe* = %15.5g ff = %15.5g\n", m_feSpecies_curr[kspec], 
		m_SSfeSpecies[kspec]);
	else
	  plogf("fe* = %15.5g ff = %15.5g dg* = %15.5g\n", 
		 m_feSpecies_curr[kspec], m_SSfeSpecies[kspec], m_deltaGRxn_new[kspec-m_numComponents]);
      }
    }
#endif
    /* ********************************************************** */
    /* **** ESTIMATE REACTION ADJUSTMENTS *********************** */
    /* ********************************************************** */
    vcs_dzero(VCS_DATA_PTR(m_deltaPhaseMoles), NPhase);
    for (iph = 0; iph < NPhase; iph++) {
      xtphMax[iph] = log(m_tPhaseMoles_new[iph] * 1.0E32);
      xtphMin[iph] = log(m_tPhaseMoles_new[iph] * 1.0E-32);
    }
    for (irxn = 0; irxn < nrxn; ++irxn) {
      kspec = ir[irxn];
      /*
       * For single species phases, we will not estimate the
       * mole numbers. If the phase exists, it stays. If it
       * doesn't exist in the estimate, it doesn't come into
       * existence here.
       */
      if (! SSPhase[kspec]) {
	iph = PhaseID[kspec];
	if (m_deltaGRxn_new[irxn] > xtphMax[iph]) m_deltaGRxn_new[irxn] = 0.8 * xtphMax[iph];
	if (m_deltaGRxn_new[irxn] < xtphMin[iph]) m_deltaGRxn_new[irxn] = 0.8 * xtphMin[iph];
	/*
	 *   HKM -> The TMolesMultiphase is a change of mine.
	 *          It more evenly distributes the initial moles amongst
	 *          multiple multispecies phases according to the
	 *          relative values of the standard state free energies.
	 *          There is no change for problems with one multispecies
	 *          phase.
	 *            It cut diamond4.vin iterations down from 62 to 14.
	 */
	m_deltaMolNumSpecies[kspec] = 0.5 * (m_tPhaseMoles_new[iph] + TMolesMultiphase) 
	  * exp(-m_deltaGRxn_new[irxn]);
	 
	for (k = 0; k < m_numComponents; ++k) {
	  m_deltaMolNumSpecies[k] += m_stoichCoeffRxnMatrix[irxn][k] * m_deltaMolNumSpecies[kspec];
	}
	
	for (iph = 0; iph < NPhase; iph++) {
	  m_deltaPhaseMoles[iph] += DnPhase[irxn][iph] * m_deltaMolNumSpecies[kspec];
	}
      }
    }
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      for (kspec = 0; kspec < nspecies; ++kspec) {
	if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  plogf("%sdirection (", pprefix); plogf("%-12.12s", SpName[kspec].c_str());
	  plogf(") = %g", m_deltaMolNumSpecies[kspec]);
	  if (SSPhase[kspec]) {
	    if (molNum[kspec] > 0.0) {
	      plogf(" (ssPhase exists at w = %g moles)", molNum[kspec]);
	    } else {
	      plogf(" (ssPhase doesn't exist -> stability not checked)");
	    }
	  }
          plogendl();
	}
      }
    }
#endif
    /* *********************************************************** */
    /* **** KEEP COMPONENT SPECIES POSITIVE ********************** */
    /* *********************************************************** */
    par = 0.5;
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
      if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	if (par < -m_deltaMolNumSpecies[kspec] / m_molNumSpecies_new[kspec]) {
	  par = -m_deltaMolNumSpecies[kspec] / m_molNumSpecies_new[kspec];
	}
      }
    }
    par = 1. / par;
    if (par <= 1.0 && par > 0.0) {
      par *= 0.8;
    } else {
      par = 1.0;
    }
    /* ******************************************** */
    /* **** CALCULATE NEW MOLE NUMBERS ************ */
    /* ******************************************** */
    finished = FALSE;
    do {
      for (kspec = 0; kspec < m_numComponents; ++kspec) {
	if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  molNum[kspec] = m_molNumSpecies_new[kspec] + par * m_deltaMolNumSpecies[kspec];
	} else {
	  m_deltaMolNumSpecies[kspec] = 0.0;
	}
      }
      for (kspec = m_numComponents; kspec < nspecies; ++kspec) {
	if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  if (m_deltaMolNumSpecies[kspec] != 0.0) molNum[kspec] = m_deltaMolNumSpecies[kspec] * par;
	}
      }
      /*
       * We have a new w[] estimate, go get the 
       * TMoles and m_tPhaseMoles_old[] values
       */
      vcs_tmoles();
      if (lt > 0)  goto finished;
      /* ******************************************* */
      /* **** CONVERGENCE FORCING SECTION ********** */
      /* ******************************************* */
      vcs_dfe(molNum, 0, 0, 0, nspecies);
      for (kspec = 0, s = 0.0; kspec < nspecies; ++kspec) {
	s += m_deltaMolNumSpecies[kspec] * m_feSpecies_curr[kspec];
      }
      if (s == 0.0) {
	finished = TRUE; continue;
      }
      if (s < 0.0) {
	if (ikl <= 0) {
	  finished = TRUE; continue;
	}
      }
      /* ***************************************** */
      /* *** TRY HALF STEP SIZE ****************** */
      /* ***************************************** */
      if (ikl <= 0) {
	s1 = s;
	par *= 0.5;
	ikl = 1;
	continue;
      }
      /* **************************************************** */
      /* **** FIT PARABOLA THROUGH HALF AND FULL STEPS ****** */
      /* **************************************************** */
      xl = (1.0 - s / (s1 - s)) * 0.5;
      if (xl < 0.0) {
	/* *************************************************** */
	/* *** POOR DIRECTION, REDUCE STEP SIZE TO 0.2 ******* */
	/* *************************************************** */
	par *= 0.2;
      } else {
	if (xl > 1.0) {
	  /* *************************************************** */
	  /* **** TOO BIG A STEP, TAKE ORIGINAL FULL STEP ****** */
	  /* *************************************************** */
	  par *= 2.0;
	} else {
	  /* *************************************************** */
	  /* **** ACCEPT RESULTS OF FORCER ********************* */
	  /* *************************************************** */
	  par = par * 2.0 * xl;
	}
      }
      lt = 1;
    } while (!finished);
  finished:
    ;
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("%s     Final Mole Numbers produced by inest:\n",
	     pprefix);
      plogf("%s     SPECIES      MOLE_NUMBER\n", pprefix);
      for (kspec = 0; kspec < nspecies; ++kspec) {
	plogf("%s     ", pprefix); plogf("%-12.12s", SpName[kspec].c_str());
	plogf(" %g", molNum[kspec]);			     
        plogendl();
      }
    }
#endif
  } /* inest() *****************************************************************/

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  int VCS_SOLVE::vcs_inest_TP(void)
   
    /**************************************************************************
     *
     * vcs_inest_TP:
     *
     *   Create an initial estimate of the solution to the thermodynamic 
     *   equilibrium problem. 
     *
     *    Return value:
     *
     *      0: successful initial guess
     *     -1: Unsuccessful initial guess, the elemental abundances aren't 
     *         satisfied.
     ***************************************************************************/
  {
    int retn = 0;
    double    test;
    Cantera::clockWC tickTock;
    test = -1.0E20;
    /*
     *  Malloc temporary space for usage in this routine and in
     *  subroutines
     *        sm[ne*ne]
     *        ss[ne]
     *        sa[ne]
     *        aw[m]
     */

  
    std::vector<double> sm(m_numElemConstraints*m_numElemConstraints, 0.0);
    std::vector<double> ss(m_numElemConstraints, 0.0);
    std::vector<double> sa(m_numElemConstraints, 0.0);
    std::vector<double> aw(m_numSpeciesTot+ m_numElemConstraints, 0.0);
    /*
     *  Go get the estimate of the solution
     */
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("%sGo find an initial estimate for the equilibrium problem",
	     pprefix);
      plogendl();
    }
#endif
    inest(VCS_DATA_PTR(aw), VCS_DATA_PTR(sa), VCS_DATA_PTR(sm),
	  VCS_DATA_PTR(ss), test);
    /*
     *      Calculate the elemental abundances
     */
    vcs_elab();
    /*
     *      If we still fail to achieve the correct elemental abundances, 
     *      try to fix the problem again by calling the main elemental abundances
     *      fixer routine, used in the main program. What this does, is that it
     *      attempts to tweak the mole numbers of the component species to
     *      satisfy the element abundance constraints. 
     *
     *       Note: We won't do this unless we have to since it involves inverting
     *             a matrix.
     */
    int rangeCheck  = vcs_elabcheck(1);
    if (!vcs_elabcheck(0)) {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	plogf("%sInitial guess failed element abundances\n", pprefix);  
	plogf("%sCall vcs_elcorr to attempt fix", pprefix);
        plogendl();
      }
#endif
      vcs_elcorr(VCS_DATA_PTR(sm), VCS_DATA_PTR(aw));
      rangeCheck  = vcs_elabcheck(1);
      if (!vcs_elabcheck(0)) {
	plogf("%sInitial guess still fails element abundance equations\n",
	       pprefix);
	plogf("%s - Inability to ever satisfy element abundance "
	       "constraints is probable", pprefix);
        plogendl();
	retn = -1;
      } else {
#ifdef DEBUG_MODE
	if (vcs_debug_print_lvl >= 2) {
	  if (rangeCheck) {
	    plogf("%sInitial guess now satisfies element abundances", pprefix);
            plogendl();
	  } else {
	    plogf("%sElement Abundances RANGE ERROR\n", pprefix);
	    plogf("%s - Initial guess satisfies NC=%d element abundances, "
		   "BUT not NE=%d element abundances", pprefix,
		   m_numComponents, m_numElemConstraints);
            plogendl();
	  }
	}
#endif 
      }
    }
    else {
#ifdef DEBUG_MODE
      if (vcs_debug_print_lvl >= 2) {
	if (rangeCheck) {
	  plogf("%sInitial guess satisfies element abundances", pprefix);
          plogendl();
	} else {
	  plogf("%sElement Abundances RANGE ERROR\n", pprefix);
	  plogf("%s - Initial guess satisfies NC=%d element abundances, "
		 "BUT not NE=%d element abundances", pprefix, 
		 m_numComponents, m_numElemConstraints);
          plogendl();
	}
      }
#endif
    } 
      
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      plogf("%sTotal Dimensionless Gibbs Free Energy = %15.7E", pprefix,
	     vcs_Total_Gibbs(VCS_DATA_PTR(m_molNumSpecies_old), VCS_DATA_PTR(m_feSpecies_curr), 
			     VCS_DATA_PTR(m_tPhaseMoles_old)));   
      plogendl();
    }
#endif

    /*
     *      Free malloced memory
     */
    double tsecond = tickTock.secondsWC();
    m_VCount->T_Time_inest += tsecond;
    (m_VCount->T_Calls_Inest)++;
    return retn;
  }

}

