/**
 *  @file vcs_util.cpp
 *  Internal definitions for utility functions for the VCSnonideal package
 */
/*
 * $Id: vcs_util.cpp,v 1.13 2009/07/06 22:30:58 hkmoffa Exp $
 */
/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#include <cstdlib>
#include <cmath>
#include <cassert>

#include "vcs_internal.h" 
#include <string.h>
#include <cstdlib>

using namespace std;

namespace VCSnonideal {

  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/
#ifndef USE_MEMSET
  void vcs_dzero(double *vector, int length)
   
    /**************************************************************************
     *
     *  vcs_dzero:
     *
     *     Zeroes a double vector
     *************************************************************************/
  {
    int i;
    for (i = 0; i < length; i++) vector[i] = 0.0; 
  } /* vcs_dzero() ***********************************************************/
#endif
  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/
#ifndef USE_MEMSET
  void vcs_izero(int *vector, int length)
   
    /**************************************************************************
     *
     *  vcs_izero:
     *
     *     Zeroes an int vector
     *************************************************************************/
  {
    int i;
    for (i = 0; i < length; i++) vector[i] = 0; 
  } /* vcs_izero() ***********************************************************/
#endif
  /***************************************************************************/
  /***************************************************************************/
  /***************************************************************************/

#ifndef USE_MEMSET
  void vcs_dcopy(double *vec_to, double *vec_from, int length)
   
    /**************************************************************************
     *
     *  vcs_dcopy:
     *
     *     Copies a double vector
     ***************************************************************************/
  {
    int i;
    for (i = 0; i < length; i++) vec_to[i] = vec_from[i]; 
  } /* vcs_dzero() *************************************************************/
#endif
  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

#ifndef USE_MEMSET
  void vcs_icopy(int *vec_to, int *vec_from, int length)
   
    /**************************************************************************
     *
     *  vcs_icopy:
     *
     *     copies an int vector
     ***************************************************************************/
  {
    int i;
    for (i = 0; i < length; i++) vec_to[i] = vec_from[i]; 
  } /* vcs_dzero() *************************************************************/
#endif

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

#ifndef USE_MEMSET 
  /*
   *  vcs_vdzero
   *
   *    zeroes a double vector
   */
  void vcs_vdzero(std::vector<double> &vvv, int len) {
    if (len < 0) {
      std::fill(vvv.begin(), vvv.end(), 0.0);
    } else {
      std::fill_n(vvv.begin(), len, 0.0);
    }
  }
#endif

  double vcs_l2norm(const std::vector<double> vec) {
     int len = vec.size();
     if (len == 0) {
       return 0.0;
     }
     double sum = 0.0;
     std::vector<double>::const_iterator pos;
     for (pos = vec.begin(); pos != vec.end(); ++pos) {
       sum += (*pos) * (*pos);
     }
    return std::sqrt(sum/len);
  }

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

#ifndef USE_MEMSET 
  /*
   *  vcs_vizero
   *
   *    zeroes a double vector
   */
  void vcs_vizero(std::vector<int> &vvv, int len) {
    if (len < 0) {
      std::fill(vvv.begin(), vvv.end(), 0.0);
    } else {
      std::fill_n(vvv.begin(), len, 0.0);
    }
  }
#endif



#ifndef USE_MEMSET 
  /*
   *  vcs_vdcopy
   *
   *    copies a vector of doubles to another vector of doubles
   *
   * @param vec_to      Vector to be copied to
   * @param vec_from    Vector to be copied from
   * @param length      Length of the copy
   */
  void vcs_vdcopy(std::vector<double> &vec_to, 
		  const std::vector<double> & vec_from, int length) {
    std::copy(vec_from.begin(), vec_from.begin() + length, vec_to.begin());
  }
#endif


#ifndef USE_MEMSET
  /*
   *  vcs_vicopy
   *
   *    copies a vector to another vector
   *
   * @param vec_to      Vector to be copied to
   * @param vec_from    Vector to be copied from
   * @param length      Length of the copy
   */
  void vcs_vicopy(std::vector<int> &vec_to, 
		  const std::vector<int> & vec_from, int length) {
    std::copy(vec_from.begin(), vec_from.begin() + length, vec_to.begin());
  }
#endif

  /*
   *
   * Finds the location of the maximum component in a double vector 
   * INPUT 
   *    x(*) - Vector to search
   *    xSize(*) if nonnull, this is the multiplier vector to be
   *             multiplied into x(*) before making the decision. 
   *    j <= i < n     : i is the range of indecises to search in X(*) 
   * 
   * RETURN
   *    return index of the greatest value on X(*) searched 
   */
  int vcs_optMax(const double *x, const double * xSize, int j, int n) {
    int i;
    int largest = j;
    double big = x[j];
    if (xSize) {
      assert(xSize[j] > 0.0);
      big *= xSize[j];
      for (i = j + 1; i < n; ++i) {
	assert(xSize[i] > 0.0);
	if ((x[i]*xSize[i]) > big) {
	  largest = i;
	  big = x[i]*xSize[i];
	}
      } 
    } else {
      for (i = j + 1; i < n; ++i) {
	if (x[i] > big) {
	  largest = i;
	  big = x[i];
	}
      }
    }
    return largest;
  } 
  
  int vcs_max_int(const int *vector, int length)
   
    /**************************************************************************
     *
     *  vcs_max_int:
     *
     *     returns the maximum integer in a list.
     ***************************************************************************/
  {
    int i, retn;
    if (vector == NULL || length <= 0) return 0;
    retn = vector[0];
    for (i = 1; i < length; i++) {
      retn = MAX( retn, vector[i]);
    }
    return retn;
  }
  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  // Swap values in a std vector string
  /*
   * Switches the value of vecStrings[i1] with vecStrings[i2]
   * 
   * @param vecStrings  Vector of integers
   * @param i1 first index
   * @param i2 second index
   */
  void vcsUtil_stsw(std::vector<std::string> & vstr, int i1, int i2) {
    std::string tmp(vstr[i2]);
    vstr[i2] = vstr[i1];
    vstr[i1] = tmp;
  }

  // Swap values in vector of doubles
  /*
   * Switches the value of x[i1] with x[i2]
   * 
   * @param x  Vector of doubles
   * @param i1 first index
   * @param i2 second index
   */
  void vcsUtil_dsw(double x[], int i1, int i2) {
    double t = x[i1];
    x[i1] = x[i2];
    x[i2] = t;
  }

  // Swap values in an integer array
  /*
   * Switches the value of x[i1] with x[i2]
   * 
   * @param x  Vector of integers
   * @param i1 first index
   * @param i2 second index
   */
  void vcsUtil_isw(int x[], int i1, int i2) {
    int t = x[i1];
    x[i1] = x[i2];
    x[i2] = t;
  }

  // Invert an n x n matrix and solve m rhs's
  /*
   * Solve a square matrix with multiple right hand sides
   *
   * \f[
   *     C X + B = 0;
   * \f]
   *
   * This routine uses Gauss elimination and is optimized for the solution
   * of lots of rhs's. A crude form of row pivoting is used here. 
   * The matrix C is destroyed.
   *
   * @return Routine returns an integer representing success:
   *     -   1 : Matrix is singluar
   *     -   0 : solution is OK
   *    The solution x[] is returned in the matrix b.
   *
   *  @param c  Matrix to be inverted. c is in fortran format, i.e., rows
   *            are the inner loop. Row  numbers equal to idem.
   *            c[i+j*idem] = c_i_j = Matrix to be inverted: i = row number
   *                                                         j = column number
   *  @param idem number of row dimensions in c
   *  @param n  Number of rows and columns in c
   *  @param b  Multiple RHS. Note, b is actually the negative of 
   *            most formulations.  Row  numbers equal to idem.
   *             b[i+j*idem] = b_i_j = vectors of rhs's:      i = row number
   *                                                          j = column number
   *            (each column is a new rhs)
   *  @param m  number of rhs's
   */
  int vcsUtil_mlequ(double *c, int idem, int n, double *b, int m) {
    int i, j, k, l;
    double R;
    if (n > idem || n <= 0) {
      plogf("vcsUtil_mlequ ERROR: badly dimensioned matrix: %d %d\n", n, idem);
      return 1;
    }
   
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
	plogf("vcsUtil_mlequ ERROR: Encountered a zero column: %d\n", i);
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
      for (j = 0; j < m; ++j) {
	b[i + j * idem] = -b[i + j * idem] / c[i + i*idem];
      }
    }
    return 0;
  }

  //  Returns the value of the gas constant in the units specified by a parameter
  /*
   *  @param mu_units Specifies the units. 
   *           -  VCS_UNITS_KCALMOL: kcal gmol-1 K-1
   *           -  VCS_UNITS_UNITLESS:  1.0 K-1
   *           -  VCS_UNITS_KJMOL:   kJ gmol-1 K-1
   *           -  VCS_UNITS_KELVIN:    1.0 K-1
   *           -  VCS_UNITS_MKS:   joules  kmol-1 K-1 =  kg m2 s-2 kmol-1 K-1 
   */
  double vcsUtil_gasConstant(int mu_units) {
    double r;
    switch (mu_units) {
    case VCS_UNITS_KCALMOL:
      r =  0.008314472/4.184;
      break;
    case VCS_UNITS_UNITLESS: 
      r = 1.0;
      break;
    case VCS_UNITS_KJMOL: 
      r = 0.008314472;
      break;
    case VCS_UNITS_KELVIN:
      r = 1.0;
      break;
    case VCS_UNITS_MKS:
      /* joules / kg-mol K = kg m2 / s2 kg-mol K */
      r = 8.314472E3;
      break;
    default:
      plogf("vcs_gasConstant error: uknown units: %d\n", 
	    mu_units);
      exit(EXIT_FAILURE);
    }
    return r;
  }



  void vcs_print_line(const char *string, int num)
   
    /**************************************************************************
     *
     * vcs_print_char:
     *
     *      Print a line consisting of a multiple of the same string
     *
     ***************************************************************************/
  {
    if (string) {
      for (int j = 0; j < num; j++) plogf("%s", string);
    }
    plogendl();
  }


  const char *vcs_speciesType_string(int speciesStatus, int length) {
    const char *sss;
    switch (speciesStatus) {
    case VCS_SPECIES_COMPONENT:
      sss = "Component Species";
      break;
    case VCS_SPECIES_MAJOR:
      sss ="Major Species";
      break;
    case VCS_SPECIES_MINOR:
      sss ="Minor Species";
      break;
    case VCS_SPECIES_ZEROEDPHASE:
      if (length < 48) {
	sss = "Set Zeroed-Phase";
      } else {
	sss ="Purposely Zeroed-Phase Species (not in problem)";
      }
      break;
    case VCS_SPECIES_ZEROEDMS:
      if (length < 23) {
	sss = "Zeroed-MS Phase";
      } else {
	sss ="Zeroed-MS Phase Species";
      }
      break;
    case VCS_SPECIES_ZEROEDSS:
      if (length < 23) {
	sss = "Zeroed-SS Phase";
      } else {
	sss ="Zeroed-SS Phase Species";
      }
      break;
    case VCS_SPECIES_DELETED:
      if (length < 22) {
	sss = "Deleted Species";
      } else if (length < 40) {
	sss = "Deleted-Small Species";
      } else {
	sss ="Deleted-Small Species in a MS phase";
      }
      break;
    case VCS_SPECIES_ACTIVEBUTZERO:
      if (length < 47) {
	sss = "Tmp Zeroed in MS";
      } else {
	sss ="Zeroed Species in an active MS phase (tmp)";
      }
      break;
    case VCS_SPECIES_STOICHZERO:
      if (length < 56) {
	sss = "Stoich Zeroed in MS";
      } else {
	sss ="Zeroed Species in an active MS phase (Stoich Constraint)";
      }
      break;
    case VCS_SPECIES_INTERFACIALVOLTAGE:
      if (length < 29) {
	sss = "InterfaceVoltage";
      } else {
	sss ="InterfaceVoltage Species";
      }
      break;
    default:
      sss = "unknown species type";
    }
    return sss;
  }


  /************************************************************************ **/

  void vcs_print_stringTrunc(const char *str, int space, int alignment) 

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
	plogf("%c", str[i]);
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
	for (i = 0; i < ls; i++) plogf(" ");
      }
      plogf("%s", str);
      if (rs != 0) {
	for (i = 0; i < rs; i++) plogf(" ");
      }
    } 
  }

  /*****************************************************************************/
  /*****************************************************************************/
  /*****************************************************************************/

  bool vcs_doubleEqual(double d1, double d2) 

    /*************************************************************************
     * vcs_doubleEqual()
     *
     *  Simple routine to check whether two doubles are equal up to
     *  roundoff error. Currently it's set to check for 10 digits of 
     *  accuracy.
     *************************************************************************/
  {
    double denom = fabs(d1) + fabs(d2) + 1.0;
    double fac = fabs(d1 - d2) / denom;
    if (fac > 1.0E-10) {
      return false;
    }
    return true;
  }

}
