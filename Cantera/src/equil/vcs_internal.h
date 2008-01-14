/*
 * $Id$
 */

/*
 * Copywrite (2005) Sandia Corporation. Under the terms of 
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#ifndef _VCS_INTERNAL_H
#define _VCS_INTERNAL_H

#include "vcs_defs.h"
#include "vcs_DoubleStarStar.h"
#include "vcs_Exception.h"

#include "global.h"

namespace VCSnonideal {

//! Points to the data in a std::vector<> object
#define VCS_DATA_PTR(vvv) (&(vvv[0]))

//! define this Cantera function to replace printf
/*!
 * We can replace this with printf easily
 */
#define plogf  Cantera::writelogf


/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/

/************ ERROR HANDLING *****************/

struct VCS_ERR {
  int Flag;
  int Species1;
  int Species2;
  int PrintLevel;
  double Value1;
  double Value2;
  char Mess[120];
};
typedef struct VCS_ERR VCS_ERR_STRUCT;


extern int vcsUtil_err_check(VCS_ERR_STRUCT &vcsE, char *string1, int ival);
extern void vcsUtil_err_reset(VCS_ERR_STRUCT &vcsE);


/*********************************/
/* Function Pointer Typedefs */
/*********************************/

typedef double (*VCS_FUNC_PTR)(double, double, int, void *, int *);

/*
 * Forward references
 */
class VCS_SPECIES_THERMO;
class VCS_PROB;


/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

//!  Amount of extra printing that is done while in debug mode.
/*!
 *                     0 -> none
 *                     1 -> some    
 *                     2 -> alot      (default)
 *                     3 -> everything
 */

class VCS_COUNTERS {
public:
  int T_Its;         /* Total number of iterations in the main loop 
			of vcs_TP() to solve for thermo equilibrium */
  int Its;               /* Current number of iterations in the main loop 
			    of vcs_TP() to solve for thermo equilibrium */
  int T_Basis_Opts;  /* Total number of optimizations of the 
			components basis set done */
  int Basis_Opts;        /* Total number of optimizations of the 
			    components basis set done */
  int T_Calls_Inest;       /* Current number of times the initial thermo
			      equilibrium estimator has been called */
  int T_Calls_vcs_TP;      /* Current number of calls to vcs_TP */
  double T_Time_vcs_TP;  /* Total time spent in vcs_TP */
  double Time_vcs_TP;    /* Current time spent in vcs_TP */
  double T_Time_basopt;  /* Total Time spent in basopt () */
  double Time_basopt;    /* Current Time spent in basopt () */
  double T_Time_inest;     /* Time spent in initial estimator */
  double T_Time_vcs;       /* Time spent in the vcs suite of programs */
};

/*****************************************************************************/
/**************** Prototypes *************************************************/
/*****************************************************************************/

/* Externals for vcs_funcVtot.c */

extern double vcs_funcVtot(double, double, int, void *, int *);

/* Externals defined in vcs_nondim.c */

extern double vcsUtil_gasConstant(int);

/* Externals in vcs_solve_TP.c */

extern int    vcsUtil_mlequ(double *, int, int, double *, int);
extern void   vcsUtil_dsw(double *, int ,int);
extern void   vcsUtil_isw(int [], int, int);
extern void   vcsUtil_ssw(char **, int, int); 
extern void   vcsUtil_stsw(std::vector<std::string> & vecStrings, int, int); 

/* Externals for vcs_root1d.c */

extern int vcsUtil_root1d(double, double, int, VCS_FUNC_PTR , void *, 
			  double , int, double *, int printLvl = 0);

/* Externals defined in vcs_timer_generic.c */

extern double vcs_second(void);


/* Externals defined in vcs_util.c */

#define USE_MEMSET
#ifdef USE_MEMSET
#include <string.h>
#  define vcs_dzero(vector, length) (void) memset((void *) (vector), 0, \
                                    (length) * sizeof(double))
#  define vcs_izero(vector, length) (void)  memset((void *) (vector), 0, \
                                    (length) * sizeof(int))
#  define vcs_dcopy(vec_to, vec_from, length) \
                     (void) memcpy((void *) (vec_to), (const void *) (vec_from), \
                                   (length) * sizeof(double))
#  define vcs_icopy(vec_to, vec_from, length) \
                     (void) memcpy((void *) (vec_to), (const void *) (vec_from), \
                                   (length) * sizeof(int))
#  define vcs_vdzero(vector, length) (void) memset(VCS_DATA_PTR(vector), 0, \
                                    (length) * sizeof(double))
#  define vcs_vizero(vector, length) (void)  memset(VCS_DATA_PTR(vector), 0, \
                                    (length) * sizeof(int))

inline void  vcs_vdcopy(std::vector<double> & vec_to, 
                        const std::vector<double> & vec_from, int length) {
 (void) memcpy((void *)&(vec_to[0]), (const void *) &(vec_from[0]), 
               (length) * sizeof(double));
}
inline void  vcs_vicopy(std::vector<int> & vec_to, 
                        const std::vector<int> & vec_from, int length) {
 (void) memcpy((void *)&(vec_to[0]), (const void *) &(vec_from[0]), 
               (length) * sizeof(int));
}
#else
extern void vcs_dzero(double *, int);
extern void vcs_izero(int *, int);
extern void vcs_dcopy(double *, double *, int);
extern void vcs_icopy(int *, int *, int);
extern void vcs_vdzero(std::vector<double> &vvv, int len = -1);
extern void vcs_vizero(std::vector<double> &vvv, int len = -1);
void vcs_vdcopy(std::vector<double> &vec_to, 
                const std::vector<double> vec_from, int len = -1);
void vcs_vicopy(std::vector<int> &vec_to, 
                const std::vector<int> vec_from, int len = -1);
#endif
extern int  vcs_amax(double *, int, int);
extern int  vcs_max_int(int *, int);
extern void vcs_print_line(const char *, int);
extern void vcs_print_stringTrunc(const char *, int, int);
extern bool vcs_doubleEqual(double, double);

/****************************************************************************/
}

#endif
/****************************************************************************/
