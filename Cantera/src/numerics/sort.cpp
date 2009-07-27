/**
 *  @file sort.cpp
 *
 *  $Id: sort.cpp,v 1.1 2007/05/04 14:40:27 dggoodwin Exp $
 */

#ifdef WIN32
#pragma warning(disable:4786)
#endif

#include "sort.h"

namespace Cantera {

  // sort (x,y) pairs by x

  void heapsort(vector_fp& x, vector_int& y) {
    int n = x.size();
    if (n < 2) return;
    doublereal rra;
    integer rrb;
    int ll = n/2;
    int iret = n-1;
    
    while (1 > 0) {
      if (ll > 0) {
	ll--;
	rra = x[ll];
	rrb = y[ll];
      }
      else {
	rra = x[iret];
	rrb = y[iret];
	x[iret] = x[0];
	y[iret] = y[0];
	iret--;
	if (iret == 0) {
	  x[0] = rra;
	  y[0] = rrb;
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
	  y[i] = y[j];
	  i = j;
	  j = j + j + 1;
	}
	else {
	  j = iret + 1;
	}
      }
      x[i] = rra;
      y[i] = rrb;
    }
  }

  void heapsort(vector_fp& x, vector_fp& y) {
    int n = x.size();
    if (n < 2) return;
    doublereal rra;
    doublereal rrb;
    int ll = n/2;
    int iret = n-1;
    
    while (1 > 0) {
      if (ll > 0) {
	ll--;
	rra = x[ll];
	rrb = y[ll];
      }
      else {
	rra = x[iret];
	rrb = y[iret];
	x[iret] = x[0];
	y[iret] = y[0];
	iret--;
	if (iret == 0) {
	  x[0] = rra;
	  y[0] = rrb;
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
	  y[i] = y[j];
	  i = j;
	  j = j + j + 1;
	}
	else {
	  j = iret + 1;
	}
      }
      x[i] = rra;
      y[i] = rrb;
    }
  }

}

