/**
 * @file basopt.cpp
 *
 *  $Author$
 *  $Date$
 *  $Revision$
 */

#include "ct_defs.h"
#include "ThermoPhase.h"
#include "MultiPhase.h"

using namespace Cantera;
using namespace std;
#ifdef DEBUG_HKM
extern int debug_print_lvl;
static void print_stringTrunc(const char *str, int space, int alignment);
#endif
static int amax(double *x, int j, int n);
static void switch_pos(vector_int &orderVector, int jr, int kspec);
static int mlequ(double *c, int idem, int n, double *b, int m);

#ifndef MIN
#define MIN(x,y) (( (x) < (y) ) ? (x) : (y))
#endif

/**
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
 * Then, calculates SC(J,I) the formation reactions for all noncomponent 
 * species in the mechanism. 
 *
 * Input 
 * --------- 
 * mphase
 * orderVectorElement
 * orderVectorSpecies
 * 
 * Output 
 * --------- 
 * usedZeroedSpecies = If true, then a species with a zero concentration
 *                     was used as a component. The problem may be
 *                     converged.
 * formRxnMatrix 
 *
 * Return
 * --------------
 * returns the number of components.
 *
 *
 */
int Cantera::BasisOptimize(int *usedZeroedSpecies, bool doFormRxn,
		    MultiPhase *mphase, vector_int & orderVectorSpecies,
		    vector_int & orderVectorElements, 
                    vector_fp & formRxnMatrix) {
  int  j, jj, k, kk, l, i, jl, ml;
  bool lindep;
  std::string ename;
  std::string sname;
  /*
   * Get the total number of elements defined in the multiphase object
   */
  int ne = mphase->nElements();
  /*
   * Get the total number of species in the multiphase object
   */
  int nspecies = mphase->nSpecies();
  doublereal tmp;
  doublereal const USEDBEFORE = -1;
 
  /*
   * Perhaps, initialize the element ordering
   */
  if ((int) orderVectorElements.size() < ne) {
    orderVectorElements.resize(ne);
    for (j = 0; j < ne; j++) {
      orderVectorElements[j] = j;
    }
  }

  /*
   * Perhaps, initialize the species ordering
   */
  if ((int) orderVectorSpecies.size() != nspecies) {
    orderVectorSpecies.resize(nspecies);
    for (k = 0; k < nspecies; k++) {
      orderVectorSpecies[k] = k;
    }
  }
 
#ifdef DEBUG_HKM
  if (debug_print_lvl >= 1) {
    printf("   "); for(i=0; i<77; i++) printf("-"); printf("\n");
    printf("   --- Subroutine BASOPT called to ");
    printf("calculate the number of components and ");
    printf("evaluate the formation matrix\n");
    if (debug_print_lvl > 0) {
      printf("   ---\n");
      
      printf("   ---      Formula Matrix used in BASOPT calculation\n");
      printf("   ---      Species | Order | ");
      for (j = 0; j < ne; j++) {
	jj = orderVectorElements[j];
	printf(" ");
	ename = mphase->elementName(jj);
	print_stringTrunc(ename.c_str(), 4, 1);
	printf("(%1d)", j);
      }
      printf("\n");
      for (k = 0; k < nspecies; k++) {
	kk = orderVectorSpecies[k];
	printf("   --- ");
	sname = mphase->speciesName(kk);
	print_stringTrunc(sname.c_str(), 11, 1);
	printf(" |   %4d |", k);
	for (j = 0; j < ne; j++) {
	  jj = orderVectorElements[j]; 
	  double num = mphase->nAtoms(kk,jj);
	  printf("%6.1g  ", num);
	}
	printf("\n");
      }
      printf("   --- \n");
    }
  }
#endif
   
  /*
   *  Calculate the maximum value of the number of components possible
   *     It's equal to the minimum of the number of elements and the
   *     number of total species.
   */
  int nComponents = MIN(ne, nspecies);
  int nNonComponents = nspecies - nComponents;
  /*
   * Set this return variable to false
   */
  *usedZeroedSpecies = false;

  /*
   * Create an array of mole numbers
   */
  vector_fp molNum(nspecies,0.0);
  mphase->getMoles(DATA_PTR(molNum));

  /*
   * Other workspace
   */
  vector_fp sm(ne*ne, 0.0);
  vector_fp ss(ne, 0.0);
  vector_fp sa(ne, 0.0);
  if ((int) formRxnMatrix.size() < nspecies*ne) {
    formRxnMatrix.resize(nspecies*ne, 0.0);
  }

#ifdef DEBUG_HKM
  /*
   * For debugging purposes keep an unmodified copy of the array.
   */
  vector_fp molNumBase(molNum);
#endif
  

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
       *    Search the remaining part of the mole number vector, molNum
       *    for the largest remaining species. Return its identity. 
       */
      k = amax(DATA_PTR(molNum), jr, nspecies);
      if (molNum[k] == 0.0) *usedZeroedSpecies = true;
      /*
       * If the largest molNum is negative, then we are done.
       */
      if (molNum[k] == USEDBEFORE) {
	nComponents = jr;
	nNonComponents = nspecies - nComponents;
	goto L_END_LOOP;
      }
      /*
       *  Assign a small negative number to the component that we have
       *  just found, in order to take it out of further consideration.
       */
      molNum[k] = USEDBEFORE;
      /* *********************************************************** */
      /* **** CHECK LINEAR INDEPENDENCE WITH PREVIOUS SPECIES ****** */
      /* *********************************************************** */
      /*    
       *          Modified Gram-Schmidt Method, p. 202 Dalquist 
       *          QR factorization of a matrix without row pivoting. 
       */
      jl = jr;
      for (j = 0; j < ne; ++j) {
	jj = orderVectorElements[j];
	sm[j + jr*ne] = mphase->nAtoms(k,jj);
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
	  for (i = 0; i < ne; ++i) {
	    ss[j] += sm[i + jr*ne] * sm[i + j*ne];
	  }
	  ss[j] /= sa[j];
	}
	/* 
	 *     Now make the new column, (*,JR), orthogonal to the 
	 *     previous columns
	 */
	for (j = 0; j < jl; ++j) {
	  for (l = 0; l < ne; ++l) {
	    sm[l + jr*ne] -= ss[j] * sm[l + j*ne];
	  }
	}
      }
      /*
       *        Find the new length of the new column in Q. 
       *        It will be used in the denominator in future row calcs. 
       */
      sa[jr] = 0.0;
      for (ml = 0; ml < ne; ++ml) {
	tmp = sm[ml + jr*ne];
	sa[jr] += tmp * tmp;
      }
      /* **************************************************** */
      /* **** IF NORM OF NEW ROW  .LT. 1E-3 REJECT ********** */
      /* **************************************************** */
      if (sa[jr] < 1.0e-6)  lindep = true;
      else                  lindep = false;
    } while(lindep);
    /* ****************************************** */
    /* **** REARRANGE THE DATA ****************** */
    /* ****************************************** */
    if (jr != k) {
#ifdef DEBUG_HKM
      if (debug_print_lvl >= 1) {
	kk = orderVectorSpecies[k];
	sname = mphase->speciesName(kk);
	printf("   ---   %-12.12s", sname.c_str());
	jj = orderVectorSpecies[jr];
	ename = mphase->speciesName(jj);
	printf("(%9.2g) replaces %-12.12s", molNum[kk], ename.c_str());
	printf("(%9.2g) as component %3d\n", molNum[jj], jr);
      }
#endif
      switch_pos(orderVectorSpecies, jr, k);
    }
    /* - entry point from up above */
  L_END_LOOP: ;
    /*
     *      If we haven't found enough components, go back 
     *      and find some more. (nc -1 is used below, because
     *      jr is counted from 0, via the C convention.
     */
  } while (jr < (nComponents-1));
   

 if (! doFormRxn) return nComponents;

  /* ****************************************************** */
  /* **** EVALUATE THE STOICHIOMETRY ********************** */
  /* ****************************************************** */
  /*
   *  Formulate the matrix problem for the stoichiometric
   *  coefficients. CX + B = 0
   *      C will be an nc x nc matrix made up of the formula 
   * vectors for the components. Each component's formular
   * vector is a column. The rows are the elements.
   *      n rhs's will be solved for. Thus, B is an nc x n
   * matrix. 
   *
   * BIG PROBLEM 1/21/99:
   *
   *    This algorithm makes the assumption that the
   * first nc rows of the formula matrix aren't rank deficient.
   * However, this might not be the case. For example, assume
   * that the first element in FormulaMatrix[] is argon. Assume that
   * no species in the matrix problem actually includes argon.
   * Then, the first row in sm[], below will be indentically
   * zero. bleh. 
   *    What needs to be done is to perform a rearrangement
   * of the ELEMENTS -> i.e. rearrange, FormulaMatrix, sp, and gai, such
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
  for (k = 0; k < nComponents; ++k) {
    kk = orderVectorSpecies[k];
    for (j = 0; j < nComponents; ++j) {
      jj = orderVectorElements[j];
      sm[j + k*ne] = mphase->nAtoms(kk, jj);
    }
  }

  for (i = 0; i < nNonComponents; ++i) {
    k = nComponents + i;
    kk = orderVectorSpecies[k];
    for (j = 0; j < nComponents; ++j) {
      jj = orderVectorElements[j];
      formRxnMatrix[j + i * ne] = mphase->nAtoms(kk, jj);
    }
  }
  /*
   *     Use Gauss-Jordon block elimination to calculate
   *     the reaction matrix 
   */
  j = mlequ(DATA_PTR(sm), ne, nComponents, DATA_PTR(formRxnMatrix), nNonComponents);
  if (j == 1) {
    printf("ERROR: mlequ returned an error condition\n");
    throw CanteraError("basopt", "mlequ returned an error condition");
  }
    
#ifdef DEBUG_HKM
  if (debug_print_lvl >= 1) {
    printf("   ---\n");
    printf("   ---  Number of Components = %d\n", nComponents);
    printf("   ---  Formula Matrix:\n");
    printf("   ---                      Components:    ");
    for (k = 0; k < nComponents; k++) {
      kk = orderVectorSpecies[k];
      printf(" %3d (%3d) ", k, kk);
    }
    printf("\n   ---                Components Moles:       ");
    for (k = 0; k < nComponents; k++) {
      kk = orderVectorSpecies[k];
      printf("%-11.3g", molNumBase[kk]);
    }
    printf("\n   ---        NonComponent |   Moles  |       ");
    for (i = 0; i < nComponents; i++) {
      kk = orderVectorSpecies[i];
      sname = mphase->speciesName(kk);
      printf("%-11.10s", sname.c_str());
    }
    printf("\n");
  
    for (i = 0; i < nNonComponents; i++) {
      k = i + nComponents;
      kk = orderVectorSpecies[k];
      printf("   --- %3d (%3d) ", k, kk);
      sname = mphase->speciesName(kk);
      printf("%-10.10s", sname.c_str());
      printf("|%10.3g|", molNumBase[kk]);
      /*
       * Print the negative of formRxnMatrix[]; it's easier to interpret.
       */
      for (j = 0; j < nComponents; j++) {
	printf("     %6.2f", - formRxnMatrix[j + i * ne]);
      }
      printf("\n");
    }
    printf("   "); for (i=0; i<77; i++) printf("-"); printf("\n");
  }
#endif

  return nComponents;
} /* basopt() ************************************************************/



#ifdef DEBUG_HKM
static void print_stringTrunc(const char *str, int space, int alignment)

   /***********************************************************************
    *  vcs_print_stringTrunc():
    *
    *     Print a string within a given space limit. This routine
    *     limits the amount of the string that will be printed to a
    *     maximum of "space" characters.
    *
    *     str = String -> must be null terminated.
    *     space = space limit for the printing.
    *     alignment = 0 centered
    *           1 right aligned
    *           2 left aligned
    ***********************************************************************/
{
  int i, ls=0, rs=0;
  int len = strlen(str);
  if ((len) >= space) {
    for (i = 0; i < space; i++) {
      printf("%c", str[i]);
    }
  } else {
    if (alignment == 1) {
      ls = space - len;
    } else if (alignment == 2) {
      rs = space - len;
    } else {
      ls = (space - len) / 2;
      rs = space - len - ls;
    }
    if (ls != 0) {
      for (i = 0; i < ls; i++) printf(" ");
    }
    printf("%s", str);
    if (rs != 0) {
      for (i = 0; i < rs; i++) printf(" ");
    }
  }
}
#endif

/*
 * Finds the location of the maximum component in a double vector
 * INPUT
 *    x(*) - Vector to search
 *    j <= i < n     : i is the range of indecises to search in X(*)
 *
 * RETURN
 *    return index of the greatest value on X(*) searched
 */
static int amax(double *x, int j, int n) {
  int i;
  int largest = j;
  double big = x[j];
  for (i = j + 1; i < n; ++i) {
    if (x[i] > big) {
      largest = i;
      big = x[i];
    }
  }
  return largest;
}


 static void switch_pos(vector_int &orderVector, int jr, int kspec) {
   int kcurr = orderVector[jr];
   orderVector[jr] = orderVector[kspec];
   orderVector[kspec] = kcurr;
 }

 /*
  * vcs_mlequ:
  *
  *  Invert an nxn matrix and solve m rhs's
  *
  *    Solve         C X + B = 0;
  *
  * This routine uses Gauss elimination and is optimized for the solution
  * of lots of rhs's.
  * A crude form of row pivoting is used here.
  *
  *
  * c[i+j*idem] = c_i_j = Matrix to be inverted: i = row number
  *                                              j = column number
  * b[i+j*idem] = b_i_j = vectors of rhs's:      i = row number
  *                                              j = column number
  *            (each column is a new rhs)
  * n = number of rows and columns in the matrix
  * m = number of rhs to be solved for
  * idem = first dimension in the calling routine
  *        idem >= n must be true
  *
  * Return Value
  *      1 : Matrix is singluar
  *      0 : solution is OK
  *
  *      The solution is returned in the matrix b.
  */
 static int mlequ(double *c, int idem, int n, double *b, int m) {
   int i, j, k, l;
   double R;

   /*
    * Loop over the rows
    *    -> At the end of each loop, the only nonzero entry in the column
    *       will be on the diagonal. We can therfore just invert the
    *       diagonal at the end of the program to solve the equation system.
    */
   for (i = 0; i < n; ++i) {
     if (c[i + i * idem] == 0.0) {
       /*
	*   Do a simple form of row pivoting to find a non-zero pivot
	*/
       for (k = i + 1; k < n; ++k) {
	 if (c[k + i * idem] != 0.0) goto FOUND_PIVOT;
       }
       printf("vcs_mlequ ERROR: Encountered a zero column: %d\n", i);
       return 1;
     FOUND_PIVOT: ;
       for (j = 0; j < n; ++j) c[i + j * idem] += c[k + j * idem];
       for (j = 0; j < m; ++j) b[i + j * idem] += b[k + j * idem];
     }

     for (l = 0; l < n; ++l) {
       if (l != i && c[l + i * idem] != 0.0) {
	 R = c[l + i * idem] / c[i + i * idem];
	 c[l + i * idem] = 0.0;
	 for (j = i+1; j < n; ++j) c[l + j * idem] -= c[i + j * idem] * R;
	 for (j = 0; j < m; ++j)   b[l + j * idem] -= b[i + j * idem] * R;
       }
     }
   }
   /*
    *  The negative in the last expression is due to the form of B upon
    *  input
    */
   for (i = 0; i < n; ++i) {
     for (j = 0; j < m; ++j)
       b[i + j * idem] = -b[i + j * idem] / c[i + i*idem];
   }
   return 0;
 } /* mlequ() *************************************************************/
