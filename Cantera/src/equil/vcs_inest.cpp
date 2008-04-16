/**
 * @file vcs_inest.cpp
 *   Methods for obtaining a good initial guess
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
  
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vcs_solve.h"
#include "vcs_internal.h" 
#include "vcs_VolPhase.h"

#include "clockWC.h"

namespace VCSnonideal {

  static char pprefix[20] = "   --- vcs_inest: ";
  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  void VCS_SOLVE::inest(double *aw, double *sa, double *sm, 
			double *ss, double test)
   
    /**************************************************************************
     *
     * inest:
     *
     *    Estimates equilibrium compositions.
     * Algorithm covered in a section of Smith and Missen's Book.
     *
     * Linear programming module is based on using dbolm.
     ***************************************************************************/
  {
    int conv, k,  lt, ikl, kspec, iph, irxn;
    double s;
    double s1 = 0.0;
    double xl, par;
    int finished;
    int     nspecies   = m_numSpeciesTot;
    int     nrxn       = m_numRxnTot;
    vcs_VolPhase *Vphase = 0;
 
    double *molNum   = VCS_DATA_PTR(soln);
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
	plogf(" %15.5g  %12.3g\n", molNum[kspec], -ff[kspec]);			     
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
	  plogf(" %12.3g %12.3g\n", gai[j], tmp);
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
    vcs_dzero(VCS_DATA_PTR(ds), nspecies);
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
     * Calculate TMoles and TPhMoles[]
     */
    vcs_tmoles();
    /*
     * TPhMoles1[] will consist of just the component moles
     */
    for (iph = 0; iph < NPhase; iph++) {
      TPhMoles1[iph] = TPhInertMoles[iph] + 1.0E-20;
    }
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
      if (SpeciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
	TPhMoles1[PhaseID[kspec]] += molNum[kspec];
      }
    }
    TMolesMultiphase = 0.0;
    for (iph = 0; iph < NPhase; iph++) {
      if (! VPhaseList[iph]->SingleSpecies) {
	TMolesMultiphase += TPhMoles1[iph];
      }
    }     
    vcs_dcopy(VCS_DATA_PTR(wt), molNum,  nspecies);
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
      if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_MOLNUM) {
	wt[kspec] = 0.0;
      }
    }
    vcs_dcopy(VCS_DATA_PTR(m_gibbsSpecies), VCS_DATA_PTR(ff), nspecies);
 
    for (kspec = 0; kspec < m_numComponents; ++kspec) {
      if (SpeciesUnknownType[kspec] == VCS_SPECIES_TYPE_MOLNUM) {
	if (! SSPhase[kspec]) {
	  iph = PhaseID[kspec];
	  m_gibbsSpecies[kspec] += log(wt[kspec] / TPhMoles[iph]);
	}
      } else {
	wt[kspec] = 0.0;
      }
    }
    vcs_deltag(0, true);
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      for (kspec = 0; kspec < nspecies; ++kspec) {
	plogf("%s", pprefix); plogf("%-12.12s", SpName[kspec].c_str());
	if (kspec < m_numComponents)
	  plogf("fe* = %15.5g ff = %15.5g\n", m_gibbsSpecies[kspec], ff[kspec]);
	else
	  plogf("fe* = %15.5g ff = %15.5g dg* = %15.5g\n", 
		 m_gibbsSpecies[kspec], ff[kspec], dg[kspec-m_numComponents]);
      }
    }
#endif
    /* ********************************************************** */
    /* **** ESTIMATE REACTION ADJUSTMENTS *********************** */
    /* ********************************************************** */
    vcs_dzero(VCS_DATA_PTR(DelTPhMoles), NPhase);
    for (iph = 0; iph < NPhase; iph++) {
      xtphMax[iph] = log(TPhMoles1[iph] * 1.0E32);
      xtphMin[iph] = log(TPhMoles1[iph] * 1.0E-32);
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
	if (dg[irxn] > xtphMax[iph]) dg[irxn] = 0.8 * xtphMax[iph];
	if (dg[irxn] < xtphMin[iph]) dg[irxn] = 0.8 * xtphMin[iph];
	/*
	 *   HKM -> The TMolesMultiphase is a change of mine.
	 *          It more evenly distributes the initial moles amongst
	 *          multiple multispecies phases according to the
	 *          relative values of the standard state free energies.
	 *          There is no change for problems with one multispecies
	 *          phase.
	 *            It cut diamond4.vin iterations down from 62 to 14.
	 */
	ds[kspec] = 0.5 * (TPhMoles1[iph] + TMolesMultiphase) 
	  * exp(-dg[irxn]);
	 
	for (k = 0; k < m_numComponents; ++k) {
	  ds[k] += sc[irxn][k] * ds[kspec];
	}
	
	for (iph = 0; iph < NPhase; iph++) {
	  DelTPhMoles[iph] += DnPhase[irxn][iph] * ds[kspec];
	}
      }
    }
#ifdef DEBUG_MODE
    if (vcs_debug_print_lvl >= 2) {
      for (kspec = 0; kspec < nspecies; ++kspec) {
	if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  plogf("%sdirection (", pprefix); plogf("%-12.12s", SpName[kspec].c_str());
	  plogf(") = %g", ds[kspec]);
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
	if (par < -ds[kspec] / wt[kspec]) par = -ds[kspec] / wt[kspec];
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
	  molNum[kspec] = wt[kspec] + par * ds[kspec];
	} else {
	  ds[kspec] = 0.0;
	}
      }
      for (kspec = m_numComponents; kspec < nspecies; ++kspec) {
	if (SpeciesUnknownType[kspec] != VCS_SPECIES_TYPE_INTERFACIALVOLTAGE) {
	  if (ds[kspec] != 0.0) molNum[kspec] = ds[kspec] * par;
	}
      }
      /*
       * We have a new w[] estimate, go get the 
       * TMoles and TPhMoles[] values
       */
      vcs_tmoles();
      if (lt > 0)  goto finished;
      /* ******************************************* */
      /* **** CONVERGENCE FORCING SECTION ********** */
      /* ******************************************* */
      vcs_dfe(molNum, 0, 0, 0, nspecies);
      for (kspec = 0, s = 0.0; kspec < nspecies; ++kspec) {
	s += ds[kspec] * m_gibbsSpecies[kspec];
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
	     vcs_Total_Gibbs(VCS_DATA_PTR(soln), VCS_DATA_PTR(m_gibbsSpecies), 
			     VCS_DATA_PTR(TPhMoles)));   
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
  }/**** vcs_inest() ***********************************************************/

}

