/*
 * @file sortAlgorithms.h
 *
 * $Author: hkmoffa $
 * $Revision: 1.1 $
 * $Date: 2006/07/06 21:13:39 $
 */
/*
 * Copywrite 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#include "sortAlgorithms.h"

/**************************************************************/

void sort_dbl_1(double * const x, const int n) {
    double rra;
    int ll = n/2;
    int iret = n - 1;
    while (1 > 0) {
      if (ll > 0) {
        ll--;
        rra = x[ll];
      } else {
        rra = x[iret];
        x[iret] = x[0];
        iret--;
        if (iret == 0) {
          x[0] = rra;
          return;
        }
      }     
      int i = ll;
      int j = ll + ll + 1;
      while (j <= iret) {
        if (j < iret) {
          if (x[j] < x[j+1])
	      j++;
        }
        if (rra < x[j]) {
          x[i] = x[j];
          i = j;
          j = j + j + 1;
        } else {
          j = iret + 1;
        }
      }
      x[i] = rra;
    }
}
/*****************************************************/
