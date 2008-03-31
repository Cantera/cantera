/* ======================================================================= */
/* -------------------------------------------------- */
/* | RCS Head Information on zuzax.pchem.sandia.gov | */
/* -------------------------------------------------- */
/* $RCSfile$ */
/* $Author$ */
/* $Date$ */
/* $Revision$ */
/* ======================================================================= */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "vcs_solve.h"
#include "vcs_internal.h" 

namespace VCSnonideal {
 
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int VCS_SOLVE::vcs_rearrange(void)
   /**************************************************************************
   *
   *  vcs_rearrange:
   *
   *     Switch all species data back to the original order. This destroys
   *  the data based on reaction ordering.
   **************************************************************************/
{
   int i, l, j;
   int k1 = 0;
   /* ********************************************************* */
   /* **** RE-ARRANGE INPUT DATA ****************************** */
   /* ********************************************************* */
   /* - Loop over all of the species */
   for (i = 0; i < m_numSpeciesTot; ++i) {
      /*
      *         Find the index of I in the index vector IND. 
      *         Call it K1 and continue. 
      */
      for (j = 0; j < m_numSpeciesTot; ++j) {
	 l = ind[j];
	 k1 = j;
	 if (l == i) break;  
      }
      /* 
      * - Switch the species data back from K1 into I
      * -> because we loop over all species, reaction data
      *    are now permanently hosed.
      */
      vcs_switch_pos(FALSE, i, k1);
   }
   return 0;
} /* vcs_rearrange() *********************************************************/

}

