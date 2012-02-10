#ifndef MDP_ALLO_H
#define MDP_ALLO_H

#include <string.h>
/*****************************************************************************/
/*
*  If we have array_alloc() from another Sandia program, we will not use
*  the one from this mdp_array_alloc. Instead we will redefine the names
*/

#ifdef HAVE_ARRAY_ALLOC
#  define mdp_array_alloc array_alloc
#  define mdp_safe_free   safe_free
#endif

/*
 *  These are a poor man's way of specifying whether a value should be
 *  initialized. These are seldom used numbers whic
 * $Name:  $
 *====================h can be used in place
 *  of real ints and dbls to indicate that initialization shouldn't take
 *  place.
 */
#define MDP_INT_NOINIT          -68361
#define MDP_DBL_NOINIT          -1.241E11


#ifndef _C16_NAME_DEF
#  define _C16_NAME_DEF
typedef char    C16_NAME[16]; /* Character array used for fortran names */
typedef char    C16_NAME_STR[17];
#endif
/*****************************************************************************/
/*
 * Externals that should be set by the calling program.
 * These are only used for debugging purposes.
 */
#ifdef MDP_MPDEBUGIO
extern int MDP_MP_Nprocs;
extern int MDP_MP_myproc;
#endif
/*****************************************************************************/

#define mdp_alloc_struct(x, num) (x *) mdp_array_alloc(1, (num), sizeof(x))


/* function declarations for dynamic array allocation */

extern double* mdp_array_alloc(int numdim, ...);
extern void    mdp_safe_free(void**);

extern int*    mdp_alloc_int_1(int, const int);
extern void    mdp_safe_alloc_int_1(int**, int, const int);
extern void    mdp_realloc_int_1(int**, int, int,
                                 const int defval = MDP_INT_NOINIT);
extern int**   mdp_alloc_int_2(int, int, const int);

extern double* mdp_alloc_dbl_1(int, const double);
extern void    mdp_safe_alloc_dbl_1(double**, int, const double);
extern void    mdp_realloc_dbl_1(double**, int, int, const double);
extern void    mdp_realloc_dbl_2(double***, int, int, int, int,
                                 const double);

extern char*   mdp_alloc_char_1(int, const char);
extern void    mdp_safe_alloc_char_1(char**, int, const char);
extern char**  mdp_alloc_VecFixedStrings(int, int);
extern void    mdp_safe_alloc_VecFixedStrings(char***, int, int);
extern void    mdp_realloc_VecFixedStrings(char***, int,  int, int);

extern double** mdp_alloc_dbl_2(int, int, const double);
extern void    mdp_safe_alloc_dbl_2(double***, int, int,
                                    const double = MDP_DBL_NOINIT);

extern C16_NAME* mdp_alloc_C16_NAME_1(int, const int);
extern void    mdp_safe_alloc_C16_NAME_1(C16_NAME**, int, const int);

extern void**  mdp_alloc_ptr_1(int);
extern void    mdp_safe_alloc_ptr_1(void***, int);
extern void    mdp_realloc_ptr_1(void***, int, int);
extern char*   mdp_copy_C16_NAME_to_string(const C16_NAME);
extern char*   mdp_copy_string(const char*);
extern void    mdp_safe_copy_string(char**, const char*);

extern void*** mdp_alloc_ptr_2(int, int);

/*****************************************************************************/
#endif
/*****************************************************************************/

