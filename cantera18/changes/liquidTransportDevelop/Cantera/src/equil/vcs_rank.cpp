
/*!
 * @file vcs_rank.cpp
 *    Header file for the internal class that holds the problem.
 */
/*
 * $Id: vcs_solve.cpp 735 2011-07-25 14:44:41Z hkmoffa $
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */


#include "vcs_solve.h"
#include "vcs_Exception.h"
#include "vcs_internal.h"
#include "vcs_prob.h"

#include "vcs_VolPhase.h"
#include "vcs_SpeciesProperties.h"
#include "vcs_species_thermo.h"

#include "clockWC.h"

#include <string>
#include "math.h"
using namespace std;

namespace VCSnonideal {
  //====================================================================================================================
  static int basisOptMax1(const double * const molNum,
		  const int n) {
    // int largest = 0;
    
    for (int i = 0; i < n; ++i) {
 
      if (molNum[i] > -1.0E200 && fabs(molNum[i]) > 1.0E-13) {
	return i;
      }
    }
    for (int i = 0; i < n; ++i) {
 
      if (molNum[i] > -1.0E200) {
	return i;
      }
    }
    return n-1;
  }
  //====================================================================================================================
  //  Calculate the rank of a matrix and return the rows and columns that will generate an independent basis
  //  for that rank
  /*
   * Choose the optimum component species basis for the calculations, finding the rank and 
   * set of linearly independent rows for that calculation.
   * Then find the set of linearly indepedent element columns that can support that rank.
   * This is done by taking the transpose of the matrix and redoing the same calculation.
   * (there may be a better way to do this. I don't know.)
   *
   *
   * Input 
   * --------- 
   *
   * @param awtmp      Vector of mole numbers which will be used to construct a 
   *                   ranking for how to pick the basis species. This is largely ignored
   *                   here.
   *
   * @param numSpecies Number of species. This is the number of rows in the matrix.
   *
   * @param matrix     Matrix. This is the formula matrix. Nominally, the rows are species, while
   *                   the columns are element compositions. However, this routine
   *                   is totally general, so that the rows and columns can be anything.
   *
   * @param numElemConstraints Number of element constraints
   *
   * Output 
   * --------- 
   * @param usedZeroedSpecies = If true, then a species with a zero concentration
   *                            was used as a component.
   *
   *
   * @param compRes    Vector of rows which are linearly independent. (these are the components)
   *
   * @param elemComp   Vector of columns which are linearly independent (These are the actionable element
   *                   constraints).
   *
   * @return        Returns number of components. This is the rank of the matrix
   */
  int VCS_SOLVE::vcs_rank(const double * awtmp, int numSpecies,  const double matrix[], int numElemConstraints,
			  std::vector<int> &compRes, std::vector<int>& elemComp, int * const usedZeroedSpecies) const 
  {

    int  j, k, l, i, jl, ml, lindep;
    int numComponents = 0;

    compRes.clear();
    elemComp.clear();
    vector<double> sm(numElemConstraints*numSpecies);
    vector<double> sa(numSpecies);
    vector<double> ss(numSpecies);

    double test = -0.2512345E298;
#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
      plogf("   "); for(i=0; i<77; i++) plogf("-"); plogf("\n");
      plogf("   --- Subroutine vcs_rank called to ");
      plogf("calculate the rank and independent rows /colums of the following matrix\n");     
      if (m_debug_print_lvl >= 5) {
	plogf("   ---     Species |  ");
	for (j = 0; j < numElemConstraints; j++) {
	  plogf(" "); 
	  plogf("   %3d  ", j);
	}
	plogf("\n");
	plogf("   ---     -----------");
	for (j = 0; j < numElemConstraints; j++) {
	  plogf("---------");
	}
	plogf("\n");
	for (k = 0; k < numSpecies; k++) {
	  plogf("   --- ");
	  plogf("  %3d  ", k);
	  plogf("     |");
	  for (j = 0; j < numElemConstraints; j++) {
	    plogf(" %8.2g", matrix[j*numSpecies + k]);
	  }
	  plogf("\n");
	}
	plogf("   ---");
	plogendl();
      }
    }
#endif
   
    /*
     *  Calculate the maximum value of the number of components possible
     *     It's equal to the minimum of the number of elements and the
     *     number of total species.
     */
    int ncTrial = MIN(numElemConstraints, numSpecies);
    numComponents = ncTrial;
    *usedZeroedSpecies = FALSE;

    /* 
     *     Use a temporary work array for the mole numbers, aw[] 
     */
    std::vector<double> aw(numSpecies);
    for (j = 0; j < numSpecies; j++) {
      aw[j] = awtmp[j];
    }

    int jr = -1;
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
	 *    Search the remaining part of the mole number vector, AW, 
	 *    for the largest remaining species. Return its identity in K. 
	 *    The first search criteria is always the largest positive
	 *    magnitude of the mole number.
	 */
	k = basisOptMax1(VCS_DATA_PTR(aw), numSpecies);

	if ((aw[k] != test) && fabs(aw[k]) == 0.0) {
	  *usedZeroedSpecies = TRUE;
	}

    
	if (aw[k] == test) {
	  numComponents = jr;

	  goto L_CLEANUP;
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
	for (j = 0; j < numElemConstraints; ++j) {
	  sm[j + jr*numElemConstraints] = matrix[j*numSpecies + k];
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
	    for (i = 0; i < numElemConstraints; ++i) {
	      ss[j] += sm[i + jr* numElemConstraints] * sm[i + j* numElemConstraints];
	    }
	    ss[j] /= sa[j];
	  }
	  /* 
	   *     Now make the new column, (*,JR), orthogonal to the 
	   *     previous columns
	   */
	  for (j = 0; j < jl; ++j) {
	    for (l = 0; l < numElemConstraints; ++l) {
	      sm[l + jr*numElemConstraints] -= ss[j] * sm[l + j*numElemConstraints];
	    }
	  }
	}
	/*
	 *        Find the new length of the new column in Q. 
	 *        It will be used in the denominator in future row calcs. 
	 */
	sa[jr] = 0.0;
	for (ml = 0; ml < numElemConstraints; ++ml) {
	  sa[jr] += SQUARE(sm[ml + jr * numElemConstraints]);
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
      compRes.push_back(k);
      elemComp.push_back(jr);
 
    } while (jr < (ncTrial-1));

  L_CLEANUP: ;
 
    if (numComponents  == ncTrial && numElemConstraints == numSpecies) {
      return numComponents;
    }
  

    int  numComponentsR = numComponents;

    ss.resize(numElemConstraints);
    sa.resize(numElemConstraints);
 
 
    elemComp.clear();
  
    aw.resize(numElemConstraints);
    for (j = 0; j < numSpecies; j++) {
      aw[j] = 1.0;
    }

    jr = -1;

    do {
      ++jr;

      do {

	k = basisOptMax1(VCS_DATA_PTR(aw), numElemConstraints);
    
	if (aw[k] == test) {
	  numComponents = jr;
	  goto LE_CLEANUP;
	}
	aw[k] = test;


	jl = jr;
	for (j = 0; j < numSpecies; ++j) {
	  sm[j + jr*numSpecies] = matrix[k*numSpecies + j];
	}
	if (jl > 0) {

	  for (j = 0; j < jl; ++j) {
	    ss[j] = 0.0;
	    for (i = 0; i < numSpecies; ++i) {
	      ss[j] += sm[i + jr* numSpecies] * sm[i + j* numSpecies];
	    }
	    ss[j] /= sa[j];
	  }

	  for (j = 0; j < jl; ++j) {
	    for (l = 0; l < numSpecies; ++l) {
	      sm[l + jr*numSpecies] -= ss[j] * sm[l + j*numSpecies];
	    }
	  }
	}

	sa[jr] = 0.0;
	for (ml = 0; ml < numSpecies; ++ml) {
	  sa[jr] += SQUARE(sm[ml + jr * numSpecies]);
	}

	if (sa[jr] < 1.0e-6)  lindep = TRUE;
	else                  lindep = FALSE;
      } while(lindep);
  
      elemComp.push_back(k);

    } while (jr < (ncTrial-1));
    numComponents = jr;
  LE_CLEANUP: ;

#ifdef DEBUG_MODE
    if (m_debug_print_lvl >= 2) {
      plogf("   --- vcs_rank found rank %d\n", numComponents);
      if (m_debug_print_lvl >= 5) {
	if (compRes.size() == elemComp.size()) {
	  printf("   ---       compRes    elemComp\n");
	  for (int i = 0; i < (int) compRes.size(); i++) {
	    printf("   ---          %d          %d \n", compRes[i], elemComp[i]);
	  }
	} else {
	  for (int i = 0; i < (int) compRes.size(); i++) {
	    printf("   ---   compRes[%d] =   %d \n", i, compRes[i]);
	  }
	  for (int i = 0; i < (int) elemComp.size(); i++) {
	    printf("   ---   elemComp[%d] =   %d \n", i, elemComp[i]);
	  }
	} 
      }
    }
#endif

    if (numComponentsR != numComponents) {
      printf("vcs_rank ERROR: number of components are different: %d %d\n", numComponentsR,  numComponents);
      throw vcsError("vcs_rank ERROR:",
			 " logical inconsistency");
      exit(-1);
    }
    return numComponents;
  }
  

}
