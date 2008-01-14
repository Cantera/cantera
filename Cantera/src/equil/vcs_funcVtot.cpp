/**
 * @file vcs_funcVtot.cpp
 *  Routine to calculate tht total volume of an extrinsic system.
 */
/*
 *  $Id$
 */
/*
 * Copywrite (2006) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vcs_solve.h"
#include "vcs_internal.h" 

#define TOL_CONV 1.0E-5

namespace VCSnonideal {
   
/**************************************************************************
 *
 *  vcs_funcVtot:
 *
 *       This is the rootfinder function call for the function vcs_TV().
 *
 */
double vcs_funcVtot(double xval, double Vtarget, int varID, void *fptrPassthrough, int *err)
{
  VCS_SOLVE *vptr = (VCS_SOLVE *) fptrPassthrough;
  static int first_time = TRUE;
  int retn;
  double vol;
  if (varID == 0) {
    vptr->T = xval;  
  } else if (varID == 1) {
    vptr->Pres = xval;
  }
#ifdef DEBUG_MODE
  retn = vptr->vcs_TP(1, 1, 10000, vptr->T, vptr->Pres);
#else
  retn = vptr->vcs_TP(0, 0, 10000, vptr->T, vptr->Pres);
#endif
  if (retn != VCS_SUCCESS) {
    plogf("vcs_funcVtot ERROR: vcs_TP returned error condition, %d\n",
	   retn);
    *err = retn;
  }
  vol = vptr->vcs_VolTotal(vptr->T, vptr->Pres, VCS_DATA_PTR(vptr->soln), 
			   VCS_DATA_PTR(vptr->VolPM));
#ifdef DEBUG_MODE
  vptr->vcs_report(retn);
#endif

   
  if (first_time) {
    first_time = FALSE;
    vptr->iest = FALSE;
  }
  return (vol - Vtarget);
}
/*****************************************************************************/
}

