/* ======================================================================= */
/* -------------------------------------------------- */
/* | RCS Head Information on zuzax.pchem.sandia.gov | */
/* -------------------------------------------------- */
/* $RCSfile$ */
/* $Author$ */
/* $Date$ */
/* $Revision$ */
/* ======================================================================= */

#include <cstdlib>
#include <cmath>

#include "vcs_internal.h" 

namespace VCSnonideal {

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifndef USE_MEMSET
void vcs_dzero(double *vector, int length)
   
   /**************************************************************************
   *
   *  vcs_dzero:
   *
   *     Zeroes a double vector
   ***************************************************************************/
{
   int i;
   for (i = 0; i < length; i++) vector[i] = 0.0; 
} /* vcs_dzero() *************************************************************/
#endif
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/
#ifndef USE_MEMSET
void vcs_izero(int *vector, int length)
   
   /**************************************************************************
   *
   *  vcs_izero:
   *
   *     Zeroes an int vector
   ***************************************************************************/
{
   int i;
   for (i = 0; i < length; i++) vector[i] = 0; 
} /* vcs_izero() *************************************************************/
#endif
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

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
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int vcs_amax(const double *x, int j, int n)
   
   /**************************************************************************
   *
   *  vcs_amax:
   *
   * Finds the location of the maximum component in a double vector 
   * INPUT 
   *    x(*) - Vector to search 
   *    j <= i < n     : i is the range of indecises to search in X(*) 
   * 
   * RETURN
   *    return index of the greatest value on X(*) searched 
   ***************************************************************************/
{
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
} /* vcs_amax() **************************************************************/
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

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

static void printSev(int severity)
{
   switch (severity) {
      case 1:
	 plogf("   VCS NOTE : ");        break;
      case 2:
	 plogf("   VCS WARNING : ");     break;
      case 3:
	 plogf("   VCS ERROR : ");       break;
      case 4:
	 plogf("   VCS FATAL ERROR : "); break;
   }   
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

int vcsUtil_err_check(VCS_ERR_STRUCT &vcsE, char *string1, int ival)
   
   /**************************************************************************
   *
   *  vcs_err_check:
   *
   *     This routine checks the error flag, prints a message, and then
   *     resets the error counters.
   *
   * Input
   *   severity:  Affects what happens in the routine.
   *           0  = Ignore warning, don't print anything
   *           1  = Note
   *           2  = Warning
   *           3  = Error -> Probably fatal
   *           4  = Immediate Fatal Error
   *   ival       = Another integer input that means different things
   *                depending on the error condition. Frequently, it will 
   *                point to a species or reaction which is the culprit.
   *          
   * Return
   *   Returns the error flag value
   ***************************************************************************/
{
   int flag = vcsE.Flag, doPrint, severity;
   /*
   *            Do a fast exit if there is no error encountered
   */
   if (flag == VCS_SUCCESS) return flag;
   switch (flag) {
      case VCS_FAILED_CONVERGENCE:
	 severity = 1;
	 break;
      case VCS_THERMO_OUTOFRANGE:
	 severity = 2;
	 break;
      case VCS_NOMEMORY:
      case VCS_PUB_BAD:
	 severity = 3;
	 break;
      case VCS_SHOULDNT_BE_HERE:
	 severity = 4;
	 break;
   }
   doPrint = (vcsE.PrintLevel <= severity);
   if (doPrint) {
      printSev(severity);
      if (string1 != NULL) plogf("%s ", string1);
   }
   switch (flag) {
      case VCS_NOMEMORY:
	 if (doPrint) {
	    plogf("Out of Memory");
	 }
	 break;
      case VCS_FAILED_CONVERGENCE:
	 if (doPrint) {
	    plogf("Failed Convergence");
	 }
	 break;
      case VCS_SHOULDNT_BE_HERE:
	 if (doPrint) {
	    plogf("Shouldn't be here, internal vcsc error");
	 }
	 break;
      case VCS_PUB_BAD:
	 if (doPrint) {
	    plogf("Public data structure is corrupt");
	 }
	 break;
      case VCS_THERMO_OUTOFRANGE:
	 if (doPrint) {
	    plogf("Thermo Data out of Range");
	    if (vcsE.Species1 >= 0) plogf(", Species = %d", vcsE.Species1);
	    if (ival >= 0)          plogf(", Species = %d", ival);
	    if (vcsE.Value1 != VCS_ERR_NOVALUE)
	       plogf(", TKelvin = %g", vcsE.Value1);
	    if (vcsE.Value1 != VCS_ERR_NOVALUE)
	       plogf(", Tlimit = %g", vcsE.Value2);
	 }
	 break;
      default:
	 if (doPrint) {
	    plogf("Unknown Error Condition = %d ??", flag);
	 }
	 break;
   }
   
   if (doPrint) {
      plogf("\n");
      if (vcsE.Mess[0] != '\0') {
        printSev(severity);
	plogf(" %s\n", vcsE.Mess);
      }
      fflush(stdout);
   }
   if (severity > 3) exit(-1);
   vcsUtil_err_reset(vcsE);
   return flag;
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

void vcsUtil_err_reset (VCS_ERR_STRUCT &vcsE)
   
   /**************************************************************************
   *
   *  vcs_err_reset:
   *
   *     Sets the error handler back to default conditions. 
   ***************************************************************************/
{
   vcsE.Flag    = VCS_SUCCESS;
   vcsE.Species1 = -1;
   vcsE.Species2 = -1;
   vcsE.Value1   = VCS_ERR_NOVALUE;
   vcsE.Value2   = VCS_ERR_NOVALUE;
   vcsE.Mess[0] = '\0';
   vcsE.Mess[119] = '\0';
}
/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

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
  plogf("\n");
}

/***************************************************************************/
/************************************************************************ **/
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

/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

}
