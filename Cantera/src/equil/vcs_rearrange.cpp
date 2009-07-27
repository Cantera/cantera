/**
 * @file vcs_rearrange.c
 *    implementation file for rearranging species.
 */
/*
 * $Id: vcs_rearrange.cpp,v 1.5 2008/12/17 16:34:18 hkmoffa Exp $
 */
/*
 * Copywrite (2007) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "vcs_solve.h"
#include "vcs_internal.h" 

#include <cstdio>
#include <cstdlib>
#include <cmath>


namespace VCSnonideal {
 

  //  Switch all species data back to the original order.
  /*
   *  This destroys the data based on reaction ordering.
   */
  int VCS_SOLVE::vcs_rearrange() {
    int i, l, j;
    int k1 = 0;
 
    /* - Loop over all of the species */
    for (i = 0; i < m_numSpeciesTot; ++i) {
      /*
       *         Find the index of I in the index vector m_speciesIndexVector[]. 
       *         Call it k1 and continue.
       */
      for (j = 0; j < m_numSpeciesTot; ++j) {
	l = m_speciesMapIndex[j];
	k1 = j;
	if (l == i) break;
      }
      /* 
       * - Switch the species data back from k1 into i
       *     -> because we loop over all species, reaction data
       *        are now permanently hosed.
       */
      vcs_switch_pos(FALSE, i, k1);
    }
    return 0;
  }

}

