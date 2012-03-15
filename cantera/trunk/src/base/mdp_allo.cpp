/**
 * @file mdp_allo.cpp
 *   Definitions for dynamic allocation of multidimensional pointer arrays
 */
/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */
#include <new>

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstdarg>

#include "cantera/base/mdp_allo.h"

using namespace std;

namespace mdp
{
/*
 *  Allocate global storage for 2 debugging ints that are used in IO of
 *  error information.
 */
#ifdef MDP_MPDEBUGIO
int MDP_MP_Nprocs = 1;
int MDP_MP_myproc = 0;
#endif
/*
 *      Error Handling
 *         7 print and exit
 *         6 exit
 *         5 print and create a divide by zero for stack trace analysis.
 *         4 create a divide by zero for stack analysis trace
 *         3 print a message and throw the bad_alloc exception.
 *         2 throw the bad_alloc exception and be quite
 *         1 print a message and return from package with the NULL pointer
 *         0 Keep completely silent about the matter and return with
 *           a null pointer.
 *
 *   -> Right now, the only way to change this option is to right here
 */
int MDP_ALLO_errorOption = 3;

inline int MinI(const int& x, const int& y)
{
    return ((x < y) ? x : y);
}
inline int MaxI(const int& x, const int& y)
{
    return ((x > y) ? x : y);
}


const int MDP_ALLOC_INTERFACE_ERROR = -230346;

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void mdp_alloc_eh(const char* const rname, const int bytes)
/*************************************************************************
 *
 * mdp_alloc_eh:
 *
 *      Error Handling
 *         7 print and exit
 *         6 exit
 *         5 print and create a divide by zero for stack trace analysis.
 *         4 create a divide by zero for stack trace analysis.
 *         3 print a message and throw the bad_alloc exception.
 *         2 throw the bad_alloc exception and be quite
 *         1 print a message and return from package with the NULL pointer
 *         0 Keep completely silent about the matter and return with
 *           a null pointer.
 **************************************************************************/
{
    double cd = 0.0;
    static char mesg[64];
    if (bytes == MDP_ALLOC_INTERFACE_ERROR) {
#ifdef MDP_MPDEBUGIO
        sprintf(mesg,"MDP_ALLOC Interface ERROR P_%d: %s", MDP_MP_my_proc,
                rname);
#else
        sprintf(mesg,"MDP_ALLOC Interface ERROR: %s", rname);
#endif
    } else {
        sprintf(mesg,"%s ERROR: out of memory while mallocing %d bytes",
                rname, bytes);
    }
    if (MDP_ALLO_errorOption % 2 == 1) {
        fprintf(stderr, "\n%s", mesg);
#ifdef MDP_MPDEBUGIO
        if (MDP_MP_Nprocs > 1) {
            fprintf(stderr,": proc = %d\n", MDP_MP_myproc);
        } else {
            fprintf(stderr,"\n");
        }
#else
        fprintf(stderr,"\n");
#endif
    }
    fflush(stderr);
    if (MDP_ALLO_errorOption == 2 || MDP_ALLO_errorOption == 3) {
        throw std::bad_alloc();
    }
    if (MDP_ALLO_errorOption == 4 || MDP_ALLO_errorOption == 5) {
        cd = 1.0 / cd;
    }
    if (MDP_ALLO_errorOption > 5) {
        std::exit(-1);
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static void mdp_alloc_eh2(const char* const rname)

/*************************************************************************
 *
 *  mdp_alloc_eh2:
 *
 *      Second Level Error Handling
 *         This routine is used at the second level.
 *         It will be called within the routine whenever another routine
 *         returns a NULL pointer.
 *
 **************************************************************************/
{
    if (MDP_ALLO_errorOption == 1 ||
            MDP_ALLO_errorOption == 3 ||
            MDP_ALLO_errorOption == 5 ||
            MDP_ALLO_errorOption == 7) {
#ifdef MDP_MPDEBUGIO
        if (MDP_MP_Nprocs > 1) {
            fprintf(stderr,": proc = %d, %s ERROR: returning with null pointer\n",
                    MDP_MP_myproc, rname);
        } else {
            fprintf(stderr,"%s ERROR: returning with null pointer", rname);
        }
#else
        fprintf(stderr,"%s ERROR: returning with null pointer", rname);
#endif
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/****************************************************************************/
#ifndef HAVE_ARRAY_ALLOC
/****************************************************************************/
static double* smalloc(size_t n);

/*****************************************************************************
 *
 *                    Dynamic Allocation of Multidimensional Arrays
 *----------------------------------------------------------------------------
 *
 * Example Usage:
 *
 *     typedef  struct
 *       {      int     bus1;
 *              int     bus2;
 *              int     dest;
 *      }       POINT;
 *
 *      POINT    **points, corner;
 *
 *      points = (POINT **) array_alloc (2, x, y, sizeof(POINT));
 *                                       ^  ^  ^
 *                                       |  |  |
 *         number of dimensions----------+  |  |
 *                                          |  |
 *          first dimension max-------------+  |
 *                                             |
 *         second dimension max----------------+
 *
 *         (points may be now be used as if it were declared
 *          POINT points[x][y])
 *
 *          This particular version is limited to dimensions of 3 or less.
 *
 *      corner = points[2][3]; (refer to the structure as you would any array)
 *
 *      free (points); (frees the entire structure in one fell swoop)
 *
 ****************************************************************************/
/*****************************************************************************
 *       The following section is a commented section containing
 *       an example main code:
 ******************************************************************************
 *double *array_alloc();
 *main()
 *{
 *  int ***temp;
 *   int *temp2;
 *   int i, j, k;
 *   int il, jl, kl;
 *
 *   malloc_debug(2);
 *   il = 2;
 *   jl = 3;
 *   kl = 3;
 *   temp = (int ***) array_alloc(3,il,jl,kl,sizeof(int));
 *   for (i=0; i<il; i++) {
 *      for (j=0; j<jl; j++) {
 *         for (k=0; k<kl; k++) temp[i][j][k] = 1;
 *      }
 *   }
 *
 *   temp2 = (int *) malloc(10*sizeof(int));
 *   for (i=0; i<10; i++) temp2[i] = 0;
 *
 *   for (i=0; i<il; i++) {
 *      for (j=0; j<jl; j++) {
 *         for (k=0; k<kl; k++) (void) printf(" %d\n", temp[i][j][k]);
 *      }
 *   }
 *   malloc_verify();
 *}
 *****************************************************************************/

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

double* mdp_array_alloc(int numdim, ...)

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
{
    int i, j;
    struct dim {
        long index;  /* Number of elements in the dimension  */
        long total;  /* Total number of elements             */
        long size;   /* Size of a single element in bytes    */
        long off;    /* offset from beginning of array       */
    } dim[4];      /* Info about each dimension            */


    long    total;        /* Total size of the array              */
    double* dfield;       /* ptr to avoid lint complaints         */
    char*   field;        /* The multi-dimensional array          */
    char**  ptr;          /* Pointer offset                       */
    char*   data;         /* Data offset                          */
    va_list va;           /* Current pointer in the argument list */

    va_start(va, numdim);

    if (numdim <= 0) {
        (void) fprintf(stderr,
                       "mdp_array_alloc ERROR: number of dimensions, %d, is <=0\n",
                       numdim);
        return NULL;
    } else if (numdim > 4) {
        (void) fprintf(stderr,
                       "mdp_array_alloc ERROR: number of dimensions, %d, is > 4\n",
                       numdim);
        return NULL;
    }

    dim[0].index = va_arg(va, int);

    if (dim[0].index <= 0) {
#ifdef MDP_MEMDEBUG
        (void) fprintf(stderr, "WARNING: mdp_array_alloc called with first "
                       "dimension <= 0, %d\n\twill return the nil pointer\n",
                       (int)(dim[0].index));
#endif
        return((double*) NULL);
    }

    dim[0].total = dim[0].index;
    dim[0].size  = sizeof(void*);
    dim[0].off   = 0;
    for (i = 1; i < numdim; i++) {
        dim[i].index = va_arg(va, int);
        if (dim[i].index <= 0) {
            (void) fprintf(stderr,
                           "WARNING: mdp_array_alloc called with dimension %d <= 0, "
                           "%d\n", i+1, (int)(dim[i].index));
            fprintf(stderr, "\twill return the nil pointer\n");
            return((double*) NULL);
        }
        dim[i].total = dim[i-1].total * dim[i].index;
        dim[i].size  = sizeof(void*);
        dim[i].off   = dim[i-1].off + dim[i-1].total * dim[i-1].size;
    }

    dim[numdim-1].size = va_arg(va, int);
    va_end(va);

    /*
     * Round up the last offset value so data is properly aligned.
     */

    dim[numdim-1].off = dim[numdim-1].size *
                        ((dim[numdim-1].off+dim[numdim-1].size-1)/dim[numdim-1].size);

    total = dim[numdim-1].off + dim[numdim-1].total * dim[numdim-1].size;

    dfield = (double*) smalloc((size_t) total);
    field  = (char*) dfield;

    for (i = 0; i < numdim - 1; i++) {
        ptr  = (char**)(field + dim[i].off);
        data = (char*)(field + dim[i+1].off);
        for (j = 0; j < dim[i].total; j++) {
            ptr[j] = data + j * dim[i+1].size * dim[i+1].index;
        }
    }

    return dfield;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

static double* smalloc(size_t n)

/**************************************************************************
 *  smalloc: safe version of malloc
 *
 *  This version of smalloc assigns space in even chunks of 8 bytes only
 *
 **************************************************************************/
{
#ifdef MDP_MEMDEBUG
    static int firsttime = 1;
    FILE* file;
#endif
    double* pntr;
    if (n == 0) {
        pntr = NULL;
    } else {
        n = ((n - 1) / 8);
        n = (n + 1) * 8;
        pntr = (double*) malloc((size_t) n);
    }
    if (pntr == NULL && n != 0) {
        if (MDP_ALLO_errorOption == 7 ||
                MDP_ALLO_errorOption == 5 || MDP_ALLO_errorOption == 3  ||
                MDP_ALLO_errorOption == 1) {
            fprintf(stderr, "smalloc : Out of space - number of bytes "
                    "requested = %d\n", (int) n);
        }
    }
#ifdef MDP_MEMDEBUG
    if (firsttime) {
        firsttime = 0;
        file = fopen("memops.txt", "w");
    } else {
        file = fopen("memops.txt", "a");
    }
    fprintf(file, "%x  %d   malloc\n", pntr, n);
    fclose(file);
#endif
    return pntr;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_safe_free(void** ptr)

/*************************************************************************
 *
 * mdp_safe_free():
 *
 *  This version of free calls the system's free function
 *  with maximum error checking. It also doesn't call free if ptr is
 *  the NULL pointer already.
 *  It will then set the freed pointer to NULL. Thus, a convention may
 *  be established wherein all pointers that can be malloced can be
 *  set to NULL if they are not malloced.
 **************************************************************************/
{
#ifdef MDP_MEMDEBUG
    FILE* file;
#endif
    if (ptr == NULL) {
        mdp_alloc_eh("mdp_safe_free: handle is NULL", MDP_ALLOC_INTERFACE_ERROR);
    }
    if (*ptr != NULL) {
#ifdef MDP_MEMDEBUG
        file = fopen("memops.txt", "a");
        Fprintf(file, "%x   free\n", *ptr);
        fflush(file);
        if ((int) *ptr == 0x00000001) {
            Fprintf(stderr, "FOUND IT!\n");
            exit(-1);
        }
        fclose(file);
#endif
        free(*ptr);
        /*
         *  Set the value of ptr to NULL, so that further references
         *  to it will be flagged.
         */
        *ptr = NULL;
    }
}
/****************************************************************************/
#endif
/*****************************************************************************
 *
 * Wrapper Functions
 * --------------------
 *
 *    The function definitions below are wrappers around array_alloc for
 *  common operations. The following principles are followed:
 *
 *  Argument dimensions le 0 are increased to 1 before calling array_alloc.
 *  Thus, something is always malloced during a call. The reason for this is
 *  that it minimizes the number of special cases in the calling program.
 *
 *  A pointer to something else other than NULL indicates that that pointer
 *  has been previously malloced. Thus, it can be freed. Note, after a free
 *  operation, this package always sets the pointer to NULL before returning.
 *
 *  "safe_alloc" routines try to free the pointer if nonNULL, before calling
 *  the base alloc_int_#() routines.
 *
 *  The regular routines  initialize the malloced space, unless told to not
 *  do so.
 *
 *  The function memcpy and memset are used where possible to increase speed.
 *
 *  Naming Convention:
 *
 *   The names of functions have the following format:
 *
 *      mdp_[safe]_[operation]_[underlyingDataType]_[Dimension]()
 *
 *         safe -> optional name added to the function name. These
 *                 functions differ in that a handle to the pointer is passed
 *                 to the function instead of just the pointer being
 *                 returned. A check to see whether the pointer has been
 *                 previously malloced is made, before a new allocation
 *                 is made. If previously malloced, the pointer is freed
 *                 first.
 *
 *         operation -> name of the operation, such as alloc, realloc,
 *                      init, or copy.
 *
 *         underlyingDataType -> dble, int, ptr, or char
 *
 *         dimension -> 1 or 2. For dimensions greater than 1, the data
 *                      is laid out in F77 or blas-compatible format
 *                      where data is contiguous wrt rows being the
 *                      inner loop.
 *
 *****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int* mdp_alloc_int_1(int nvalues, const int val)

/**************************************************************************
 *
 *  mdp_alloc_int_1:
 *
 *    Allocate and initialize a one dimensional array of integers.
 *
 *    Input
 *    -------
 *        nvalues = Length of the array
 *        val     = initialization value
 *    Return
 *    ------
 *        Pointer to the initialized integer array
 *        Failures are indicated by returning the NULL pointer.
 **************************************************************************/
{
    int* array;
    if (nvalues <= 0) {
        nvalues = 1;
    }
    array= (int*) mdp_array_alloc(1, nvalues, sizeof(int));
    if (array != NULL) {
        if (val != MDP_INT_NOINIT) {
            if (val == 0) {
                (void) memset(array, 0, sizeof(int)*nvalues);
            } else {
                for (int i = 0; i < nvalues; i++) {
                    array[i] = val;
                }
            }
        }
    } else {
        mdp_alloc_eh("mdp_alloc_int_1", nvalues * sizeof(int));
    }
    return array;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_safe_alloc_int_1(int** array_hdl, int nvalues, const int val)

/*************************************************************************
 *
 *  mdp_safe_alloc_int_1:
 *
 *    Allocates and/or initializse a one dimensional array of integers.
 *
 *    Input
 *    -------
 *        *array_hdl = Previous value of pointer. If non-NULL will try
 *                     to free the memory at this address.
 *        nvalues = Length of the array
 *        val     = intialization value
 *    Output
 *    ------
 *        *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 **************************************************************************/
{
    if (array_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_alloc_int_1: handle is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (*array_hdl != NULL) {
        mdp_safe_free((void**) array_hdl);
    }
    *array_hdl = mdp_alloc_int_1(nvalues, val);
    if (*array_hdl == NULL) {
        mdp_alloc_eh2("mdp_safe_alloc_int_1");
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void
mdp_realloc_int_1(int** array_hdl, int new_length, int old_length,
                  const int defval)

/*************************************************************************
 *
 *  mdp_realloc_int_1_(array_hdl, new_num_ptrs, old_num_ptrs);
 *
 *    Reallocates a one dimensional array of ints.
 *    This routine always allocates space for at least one int.
 *    Calls the smalloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *    Input
 *    -------
 *        array_hdl    = Pointer to the global variable that
 *                       holds the old and (eventually new)
 *                       address of the array of integers to be reallocated
 *        new_length   = Length of the array
 *        old_length   = Length of the old array
 **************************************************************************/
{
    if (new_length == old_length) {
        return;
    }
    if (new_length <= 0) {
#ifdef MDP_MPDEBUGIO
        fprintf(stderr,
                "Warning: mdp_realloc_int_1 P_%d: called with n = %d\n",
                MDP_MP_myproc, new_length);
#else
        fprintf(stderr,
                "Warning: mdp_realloc_int_1: called with n = %d\n",
                new_length);
#endif
        new_length = 1;
    }
    if (old_length < 0) {
        old_length = 0;
    }
    if (new_length == old_length) {
        return;
    }
    size_t bytenum = new_length * sizeof(int);
    int* array = (int*) smalloc(bytenum);
    if (array != NULL) {
        if (*array_hdl) {
            if (old_length > 0) {
                bytenum = sizeof(int) * old_length;
            } else {
                bytenum = 0;
            }
            if (new_length < old_length) {
                bytenum = sizeof(int) * new_length;
            }
            if (bytenum > 0) {
                (void) memcpy((void*) array, (const void*) *array_hdl, bytenum);
            }
            mdp_safe_free((void**) array_hdl);
        } else {
            old_length = 0;
        }
        *array_hdl = array;
        if ((defval != MDP_INT_NOINIT) && (new_length > old_length)) {
            if (defval == 0) {
                bytenum = sizeof(int) * (new_length - old_length);
                (void) memset((void*)(array+old_length), 0, bytenum);
            } else {
                for (int i = old_length; i < new_length; i++) {
                    array[i] = defval;
                }
            }
        }
    } else {
        mdp_alloc_eh("mdp_realloc_int_1", static_cast<int>(bytenum));
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

int** mdp_alloc_int_2(int ndim1, int ndim2, const int val)

/*************************************************************************
 *
 *  mdp_alloc_int_2:
 *
 *    Allocate and initialize a two dimensional array of ints.
 *
 *    Input
 *    -------
 *        ndim1 = Length of the first dimension of the array
 *        ndim2 = Length of the second dimension of the array
 *        val   = initialization value
 *    Return
 *    ------
 *        Pointer to the initialized integer array
 *        Failures are indicated by returning the NULL pointer.
 **************************************************************************/
{
    int i;
    int** array, *dptr;
    if (ndim1 <= 0) {
        ndim1 = 1;
    }
    if (ndim2 <= 0) {
        ndim2 = 1;
    }
    array = (int**) mdp_array_alloc(2, ndim1, ndim2, sizeof(int));
    if (array != NULL) {
        if (val != MDP_INT_NOINIT) {
            if (val == 0) {
                (void) memset((void*) array[0], 0, ndim1 * ndim2 * sizeof(int));
            } else {
                dptr = &(array[0][0]);
                for (i = 0; i < ndim1 * ndim2; i++) {
                    dptr[i] = val;
                }
            }
        }
    } else {
        mdp_alloc_eh("mdp_alloc_int_2",
                     sizeof(int) * ndim1 * ndim2 +
                     ndim1 * sizeof(void*));
    }
    return array;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

double* mdp_alloc_dbl_1(int nvalues, const double val)

/*************************************************************************
 *
 *  mdp_alloc_dbl_1:
 *
 *    Allocate and initialize a one dimensional array of doubles.
 *
 *    Input
 *    -------
 *        nvalues = Length of the array
 *        val     = initialization value
 *    Return
 *    ------
 *        Pointer to the initialized array of doubles
 *        Failures are indicated by returning the NULL pointer.
 **************************************************************************/
{
    int i;
    double* array;
    if (nvalues <= 0) {
        nvalues = 1;
    }
    array = (double*) mdp_array_alloc(1, nvalues, sizeof(double));
    if (array != NULL) {
        if (val != MDP_DBL_NOINIT) {
            if (val == 0.0) {
                (void) memset((void*) array, 0, nvalues * sizeof(double));
            } else {
                for (i = 0; i < nvalues; i++) {
                    array[i] = val;
                }
            }
        }
    } else {
        mdp_alloc_eh("mdp_alloc_dbl_1", nvalues * sizeof(double));
    }
    return array;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_safe_alloc_dbl_1(double** array_hdl, int nvalues, const double val)

/*************************************************************************
 *
 *  mdp_safe_alloc_dbl_1:
 *
 *    Allocates and/or initializes a one dimensional array of doubles.
 *
 *    Input
 *    -------
 *        *array_hdl = Previous value of pointer. If non-NULL will try
 *                     to free the memory at this address.
 *        nvalues = Length of the array
 *        val     = initialization value
 *    Output
 *    ------
 *        *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 **************************************************************************/
{
    if (array_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_alloc_dbl_1: handle is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (*array_hdl != NULL) {
        mdp_safe_free((void**) array_hdl);
    }
    *array_hdl = mdp_alloc_dbl_1(nvalues, val);
    if (*array_hdl == NULL) {
        mdp_alloc_eh2("mdp_safe_alloc_dbl_1");
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_realloc_dbl_1(double** array_hdl, int new_length,
                       int old_length, const double defval)

/*************************************************************************
 *
 *  mdp_realloc_dbl_1_(array_hdl, new_num_ptrs, old_num_ptrs);
 *
 *    Reallocates a one dimensional array of doubles.
 *    This routine always allocates space for at least one dbl.
 *    Calls the smalloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *    Input
 *    -------
 *        array_hdl    = Pointer to the global variable that
 *                       holds the old and (eventually new)
 *                       address of the array of doubles to be reallocated
 *        new_length   = Length of the array
 *        old_length   = Length of the old array
 **************************************************************************/
{
    if (new_length == old_length) {
        return;
    }
    if (new_length <= 0) {
#ifdef MDP_MPDEBUGIO
        fprintf(stderr, "Warning: mdp_realloc_dbl_1 P_%d: called with n = %d\n",
                MDP_MP_myproc, new_length);
#else
        fprintf(stderr, "Warning: mdp_realloc_dbl_1: called with n = %d\n",
                new_length);
#endif
        new_length = 1;
    }
    if (old_length < 0) {
        old_length = 0;
    }
    if (new_length == old_length) {
        return;
    }
    size_t bytenum = new_length * sizeof(double);
    double* array = (double*) smalloc(bytenum);
    if (array != NULL) {
        if (*array_hdl) {
            if (old_length > 0) {
                bytenum = sizeof(double) * old_length;
            } else {
                bytenum = 0;
            }
            if (new_length < old_length) {
                bytenum = sizeof(double) * new_length;
            }
            (void) memcpy((void*) array, (const void*) *array_hdl, bytenum);
            mdp_safe_free((void**) array_hdl);
        } else {
            old_length = 0;
        }
        *array_hdl = array;
        if ((defval != MDP_DBL_NOINIT) && (new_length > old_length)) {
            if (defval == 0) {
                bytenum = sizeof(double) * (new_length - old_length);
                (void) memset((void*)(array+old_length), 0, bytenum);
            } else {
                for (int i = old_length; i < new_length; i++) {
                    array[i] = defval;
                }
            }
        }
    } else {
        mdp_alloc_eh("mdp_realloc_dbl_1", static_cast<int>(bytenum));
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

double** mdp_alloc_dbl_2(int ndim1, int ndim2, const double val)

/*************************************************************************
 *
 *  mdp_alloc_dbl_2:
 *
 *    Allocate a two dimensional array of doubles. The array is in
 *    fortran order and can be accessed via the following form:
 *
 *          dblArray[ndim1][ndim2]
 *
 *   Note, ndim2 is the inner dimension. i.e., the array is
 *   in column ordering.
 *
 *    Input
 *    -------
 *        ndim1 = Length of the first dimension of the array
 *        ndim2 = Length of the second dimension of the array
 *        val   = initialization value
 *    Return
 *    ------
 *        Pointer to the initialized integer array
 *        Failures are indicated by returning the NULL pointer.
 **************************************************************************/
{
    int i;
    double** array, *dptr;
    if (ndim1 <= 0) {
        ndim1 = 1;
    }
    if (ndim2 <= 0) {
        ndim2 = 1;
    }
    array = (double**) mdp_array_alloc(2, ndim1, ndim2, sizeof(double));
    if (array != NULL) {
        if (val != MDP_DBL_NOINIT) {
            if (val == 0.0) {
                (void) memset((void*) array[0], 0, ndim1*ndim2 * sizeof(double));
            } else {
                dptr = &(array[0][0]);
                for (i = 0; i < ndim1*ndim2; i++) {
                    dptr[i] = val;
                }
            }
        }
    } else {
        mdp_alloc_eh("mdp_alloc_dbl_2",
                     sizeof(double) * ndim1 * ndim2 +
                     ndim1 * sizeof(void*));
    }
    return array;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

char* mdp_alloc_char_1(int nvalues, const char val)

/*************************************************************************
 *
 *  mdp_alloc_char_1:
 *
 *    Allocate and initialize a one dimensional array of characters.
 *
 *    Input
 *    -------
 *        nvalues = Length of the array
 *        val     = initialization value
 *    Return
 *    ------
 *        Pointer to the initialized character array
 *        Failures are indicated by returning the NULL pointer.
 **************************************************************************/
{
    int i;
    char* array;
    if (nvalues <= 0) {
        nvalues = 1;
    }
    array = (char*) mdp_array_alloc(1, nvalues, sizeof(char));
    if (array != NULL) {
        if (val == '\0') {
            (void) memset((void*)array, 0, sizeof(char)*nvalues);
        } else {
            for (i = 0; i < nvalues; i++) {
                array[i] = val;
            }
        }
    } else {
        mdp_alloc_eh("mdp_alloc_char_1", nvalues * sizeof(char));
    }
    return array;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_safe_alloc_char_1(char** array_hdl, int nvalues, const char val)

/*************************************************************************
 *
 *  mdp_safe_alloc_char_1:
 *
 *   Allocate and initialize a one dimensional array of characters,
 *    deallocating the space before hand.
 *
 *    This routine will free any old memory that was located at that
 *    position, before it will allocate a new vector.
 *    Allocates and/or initializes a one dimensional array of characters.
 *
 *    Input
 *    -------
 *        array_hdl =  Previous value of pointer. If non-NULL will try
 *                     to free the memory at this address.
 *        nvalues = Length of the array
 *        val     = initialization value
 *    Output
 *    ------
 *        *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 **************************************************************************/
{
    if (array_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_alloc_char_1: handle is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (*array_hdl != NULL) {
        mdp_safe_free((void**) array_hdl);
    }
    *array_hdl = mdp_alloc_char_1(nvalues, val);
    if (*array_hdl == NULL) {
        mdp_alloc_eh2("mdp_safe_alloc_char_1");
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/



void mdp_safe_alloc_dbl_2(double** *array_hdl, int ndim1, int ndim2,
                          const double val)

/*************************************************************************
 *
 *  mdp_safe_alloc_dbl_2:
 *
 *    Allocate and initialize a two dimensional array of doubles.
 *
 *    Input
 *    -------
 *        *array_hdl = Previous value of pointer. If non-NULL will try
 *                     to free the memory at this address.
 *        ndim1 = Length of the array
 *        ndim2 = Length of inner loop of the array
 *        val     = initialization value
 *    Return
 *    ------
 *        *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 **************************************************************************/
{
    if (array_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_alloc_dbl_2: handle is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (*array_hdl != NULL) {
        mdp_safe_free((void**) array_hdl);
    }
    *array_hdl = mdp_alloc_dbl_2(ndim1, ndim2, val);
    if (*array_hdl == NULL) {
        mdp_alloc_eh2("mdp_safe_alloc_dbl_2");
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_realloc_dbl_2(double** *array_hdl,  int ndim1, int ndim2,
                       int ndim1Old, int ndim2Old, const double val)

/*************************************************************************
 *
 *  mdp_realloc_dbl_2:
 *
 *  mdp_realloc_dbl_2(array_hdl, int ndim1, int ndim2,
 *                    int ndim1Old, int ndim2Old, const double val)
 *
 *    Reallocates a two dimensional array of doubles.
 *    This routine will then copy the pertinent information from
 *    the old array to the new array.
 *
 *    If both old dimensions are set to zero or less, then this routine
 *    will free the old memory before mallocing the new memory. This may
 *    be a benefit for extremely large mallocs.
 *    In all other cases, the new and the old malloced arrays will
 *    exist for a short time together.
 *
 *    Input
 *    -------
 *        array_hdl    = Pointer to the global variable that
 *                       holds the old and (eventually new)
 *                       address of the array of doubles to be reallocated
 *        ndim1        = First dimension of the new array
 *        ndim2        = Second dimension of the new array
 *        ndim1Old     = First dimension of the old array
 *        ndim2Old     = Second dimension of the old array
 *        val          = Default fill value.
 **************************************************************************/
{
    if (ndim1 <= 0) {
        ndim1 = 1;
    }
    if (ndim2 <= 0) {
        ndim2 = 1;
    }
    ndim1Old = MaxI(ndim1Old, 0);
    ndim2Old = MaxI(ndim2Old, 0);
    /*
     * One way to do it, if old information isn't needed. In this algorithm
     * the arrays are never malloced at the same time.
     */
    if ((*array_hdl == NULL) || (ndim1Old <= 0 && ndim2Old <= 0)) {
        mdp_safe_free((void**) array_hdl);
        *array_hdl = mdp_alloc_dbl_2(ndim1, ndim2, val);
        if (*array_hdl == NULL) {
            mdp_alloc_eh2("mdp_realloc_dbl_2");
        }
    }
    /*
     * Other way to do when old information is available and needed
     */
    else {
        double** array_old = *array_hdl;
        *array_hdl = (double**) mdp_array_alloc(2, ndim1, ndim2, sizeof(double));
        if (*array_hdl == NULL) {
            mdp_alloc_eh2("mdp_realloc_dbl_2");
        } else {
            /*
             * Now, let's initialize the arrays
             */
            int ndim1Min = MinI(ndim1, ndim1Old);
            int ndim2Min = MinI(ndim2, ndim2Old);
            double** array_new = *array_hdl;
            /*
             * When the second dimensions are equal, we can copy blocks
             * using the very efficient bit moving kernels.
             */
            if (ndim2 == ndim2Old) {
                size_t sz = ndim1Min * ndim2 * sizeof(double);
                (void) memcpy((void*) array_new[0], (const void*) array_old[0], sz);
            }
            /*
             * If the second dimensions aren't equal, then we have to
             * break up the bit operations even more
             */
            else {
                size_t sz = ndim2Min * sizeof(double);
                size_t sz2 = (ndim2 - ndim2Min) * sizeof(double);
                for (int i = 0; i < ndim1Min; i++) {
                    (void) memcpy((void*) array_new[i],
                                  (const void*) array_old[i], sz);
                    if (ndim2 > ndim2Min && val != MDP_DBL_NOINIT) {
                        if (val == 0.0) {
                            (void) memset((void*)(array_new[i] + ndim2Min), 0, sz2);
                        } else {
                            double* dptr = array_new[i];
                            for (int j = ndim2Min; j < ndim2; j++) {
                                dptr[j] = val;
                            }
                        }
                    }
                }
            }
            /*
             * finish up initializing the rest of the array
             */
            if (ndim1 > ndim1Min && val != MDP_DBL_NOINIT) {
                if (val == 0.0) {
                    size_t sz = (ndim1 - ndim1Min) * ndim2 * sizeof(double);
                    (void) memset((void*) array_new[ndim1Min], 0, sz);
                } else {
                    double* dptr = array_new[ndim1Min];
                    int num = (ndim1 - ndim1Min) * ndim2;
                    for (int i = 0; i < num; i++) {
                        dptr[i] = val;
                    }
                }
            }
            /*
             * Free the old array
             */
            mdp_safe_free((void**) &array_old);
        }
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

char** mdp_alloc_VecFixedStrings(int numStrings, int lenString)

/*************************************************************************
 *
 *  mdp_alloc_VecFixedStrings:
 *
 *    Allocate and initialize a vector of fixed-length
 *    strings. Each string is initialized to the NULL string.
 *
 *    Input
 *    -------
 *        numStrings = Number of strings
 *        lenString =  Length of each string including the trailing null
 *                     character
 *    Return
 *    ------
 *        This value is initialized to the correct address of the array.
 *        A NULL value in the position indicates an error.
 **************************************************************************/
{
    int i;
    char** array;
    if (numStrings <= 0) {
        numStrings = 1;
    }
    if (lenString <= 0) {
        lenString  = 1;
    }
    array = (char**) mdp_array_alloc(2, numStrings, lenString, sizeof(char));
    if (array != NULL) {
        for (i = 0; i < numStrings; i++) {
            array[i][0] = '\0';
            array[i][lenString - 1] = '\0';
        }
    } else {
        mdp_alloc_eh("mdp_alloc_VecFixedStrings",
                     sizeof(char) * numStrings * lenString +
                     numStrings * sizeof(void*));
    }
    return array;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_realloc_VecFixedStrings(char** *array_hdl, int numStrings,
                                 int numOldStrings, int lenString)

/*************************************************************************
 *
 *  mdp_realloc_VecFixedStrings:
 *
 *    Reallocate and initialize a vector of fixed-length
 *    strings. Each new string is initialized to the NULL string.
 *    old strings are copied.
 *
 *    Input
 *    -------
 *        ***array_hdl = The pointer to the char ** location holding
 *                       the data to be reallocated.
 *        numStrings = Number of strings
 *        numOldStrings = Number of old strings
 *        lenString =  Length of each string including the trailing null
 *                     character
 **************************************************************************/
{
    int i;
    char** array, **ao;
    if (numStrings <= 0) {
        numStrings = 1;
    }
    if (numStrings == numOldStrings) {
        return;
    }
    if (lenString <= 0) {
        lenString  = 1;
    }
    array = (char**) mdp_array_alloc(2, numStrings, lenString, sizeof(char));
    if (array != NULL) {
        int len = MinI(numStrings, numOldStrings);
        ao = *array_hdl;
        if (ao) {
            for (i = 0; i < len; i++) {
                strncpy(array[i], ao[i], lenString);
            }
        }
        if (numStrings > numOldStrings) {
            for (i = numOldStrings; i < numStrings; i++) {
                array[i][0] = '\0';
                array[i][lenString - 1] = '\0';
            }
        }
        mdp_safe_free((void**) array_hdl);
        *array_hdl = array;

    } else {
        mdp_alloc_eh("mdp_realloc_VecFixedStrings",
                     sizeof(char) * numStrings * lenString +
                     numStrings * sizeof(void*));
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_safe_alloc_VecFixedStrings(char** *array_hdl,
                                    int numStrings, int lenString)

/*************************************************************************
 *
 *   mdp_safe_alloc_VecFixedStrings
 *
 *    Allocate and initialize an array of strings of fixed length
 *
 *    Input
 *    -------
 *        *array_hdl = Previous value of pointer. If non-NULL will try
 *                     to free the memory at this address.
 *        numStrings = Number of strings
 *        lenString =  Length of each string including the trailing null
 *                     character
 *    Output
 *    ------
 *        *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 **************************************************************************/
{
    if (array_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_alloc_VecFixedStrings: handle is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (*array_hdl != NULL) {
        mdp_safe_free((void**) array_hdl);
    }
    *array_hdl = mdp_alloc_VecFixedStrings(numStrings, lenString);
    if (*array_hdl == NULL) {
        mdp_alloc_eh2("mdp_safe_alloc_VecFixedStrings");
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void** mdp_alloc_ptr_1(int numPointers)

/*************************************************************************
 *
 *  mdp_alloc_ptr_1:
 *
 *    Allocate and initialize a vector of pointers
 *    of type pointer to void. All pointers are initialized to the NULL
 *    value.
 *
 *    Input
 *    -------
 *        numPointers = Number of pointers
 *    Return
 *    ------
 *        This value is initialized to the correct address of the vector.
 *        A NULL value in the position indicates an error.
 **************************************************************************/
{
    int i;
    void** array;
    if (numPointers <= 0) {
        numPointers = 1;
    }
    array = (void**) mdp_array_alloc(1, numPointers, sizeof(void*));
    if (array != NULL) {
        for (i = 0; i < numPointers; i++) {
            array[i] = NULL;
        }
    } else {
        mdp_alloc_eh("mdp_alloc_ptr_1",
                     sizeof(void*) * numPointers);
    }
    return array;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_safe_alloc_ptr_1(void** *array_hdl, int numPointers)

/**************************************************************************
 *
 *  mdp_safe_alloc_ptr_1:
 *
 *    Allocate and initialize a vector of pointers
 *    of type pointer to void. All pointers are initialized to the NULL
 *    value.
 *
 *    Input
 *    -------
 *        *array_hdl = Previous value of pointer. If non-NULL will try
 *                     to free the memory at this address.
 *        numPointers = Number of pointers
 *    Output
 *    ------
 *        *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 **************************************************************************/
{
    if (array_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_alloc_ptr_1: handle is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (*array_hdl != NULL) {
        mdp_safe_free((void**) array_hdl);
    }
    *array_hdl = mdp_alloc_ptr_1(numPointers);
    if (*array_hdl == NULL) {
        mdp_alloc_eh2("mdp_safe_alloc_ptr_1");
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_realloc_ptr_1(void** *array_hdl, int numLen, int numOldLen)

/*************************************************************************
 *
 *  mdp_realloc__ptr_1:
 *
 *    Reallocate and initialize a vector of pointers
 *    Each new pointer is initialized to NULL.
 *    old Pointers are copied.
 *
 *    Input
 *    -------
 *        ***array_hdl = The pointer to the char ** location holding
 *                       the data to be reallocated.
 *        numLen = Number of strings
 *        numOldLen = Number of old strings
 **************************************************************************/
{
    if (array_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_alloc_ptr_1: handle is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (numLen <= 0) {
        numLen = 1;
    }
    if (numOldLen < 0) {
        numOldLen = 0;
    }
    if (numLen == numOldLen) {
        return;
    }
    size_t bytenum = sizeof(void*) * numLen;
    void** array = (void**) smalloc(bytenum);
    if (array != NULL) {
        int len = MinI(numLen, numOldLen);
        if (*array_hdl) {
            void** ao = *array_hdl;
            for (int i = 0; i < len; i++) {
                array[i] = ao[i];
            }
        } else {
            numOldLen = 0;
        }
        if (numLen > numOldLen) {
            bytenum = sizeof(void*) * (numLen - numOldLen);
            (void) memset((void*)(array + numOldLen), 0, bytenum);
        }
        mdp_safe_free((void**) array_hdl);
        *array_hdl = array;
    } else {
        mdp_alloc_eh("mdp_realloc_ptr_1", sizeof(void*) * numLen);
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_copy_VecFixedStrings(char** const copyTo,
                              const char** const copyFrom,
                              int numStrings, size_t maxLenString)

/*************************************************************************
 *
 *     mdp_copy_VecFixedStrings
 *
 *   Copies an array of string vectors.
 *
 *    Input
 *    -------
 *        copyFrom = vector of C strings. It should be null terminated
 *        numStrings = number of strings
 *        maxLenString = maximum of the size of the string arrays,
 *                       copyTo and copyFrom. This is used as the
 *                       argument to strncpy() function.
 *
 *    Output
 *         copyTo   = vector of strings
 *
 **************************************************************************/
{
    if (maxLenString > 0) {
        char* dTo;
        const char* dFrom;
        for (int i = 0; i < numStrings; i++) {
            dTo = copyTo[i];
            dFrom = copyFrom[i];
            if (!dTo) {
                mdp_alloc_eh2("mdp_copy_VecFixedStrings");
            }
            if (!dFrom) {
                mdp_alloc_eh2("mdp_copy_VecFixedStrings");
            }
            (void) strncpy(dTo, dFrom, maxLenString);
            dTo[maxLenString-1] = '\0';
        }
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

char* mdp_copy_string(const char* const copyFrom)

/*************************************************************************
 *
 *  mdp_copy_string:
 *
 *    Allocates space for and copies a string
 *
 *    Input
 *    -------
 *        copyFrom = null terminated string. If NULL is supplied, then
 *                   nothing is malloced and a NULL value is returned.
 *    Return
 *    ------
 *        This value is initialized to the correct address of the array.
 *        A NULL value in the position either indicates an error, or
 *        that the original pointer to the string was NULL.
 **************************************************************************/
{
    char* cptr;
    if (copyFrom == NULL) {
        return NULL;
    }
    cptr = (char*) mdp_array_alloc(1, strlen(copyFrom) + 1, sizeof(char));
    if (cptr != NULL) {
        (void) strcpy(cptr, copyFrom);
    } else {
        mdp_alloc_eh("mdp_copy_string",
                     static_cast<int>(sizeof(char) * (strlen(copyFrom) + 1)));
    }
    return cptr;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_safe_copy_string(char** string_hdl, const char* copyFrom)

/*************************************************************************
 *
 *  mdp_safe_copy_string:
 *
 *    Allocates space for and copies a string
 *
 *    Input
 *    -------
 *        *string_hdl = Previous value of pointer. If non-NULL will try
 *                     to free the memory at this address.
 *        *copyFrom = String to be copied
 *    Output
 *    ------
 *        *string_hdl = Pointer to the copied string
 *                     A NULL value in the position indicates an error.
 **************************************************************************/
{
    if (string_hdl == NULL) {
        mdp_alloc_eh("mdp_safe_copy_string: string_hdl is NULL",
                     MDP_ALLOC_INTERFACE_ERROR);
        return;
    }
    if (*string_hdl != NULL) {
        mdp_safe_free((void**) string_hdl);
    }
    if (copyFrom == NULL) {
        *string_hdl = NULL;
        return;
    }
    *string_hdl = mdp_copy_string(copyFrom);
    if (*string_hdl == NULL) {
        mdp_alloc_eh2("mdp_safe_copy_string");
    }
    return;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_copy_dbl_1(double* const copyTo,
                    const double* const copyFrom, const int len)

/**************************************************************************
 *
 * mdp_copy_dbl_1:
 *
 * Copies one contiguous double vector into another double vector
 *
 * Input
 * -------------
 *  *copyFrom = Vector of values to be copied
 *  len       = Length of the vector
 *
 * Output
 * ------------
 * *copyTo   = Vector of values to receive the copy
 *
 **************************************************************************/
{
    if (len > 0) {
        (void) memcpy((void*)copyTo, (const void*)copyFrom, len * sizeof(double));
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_copy_dbl_2(double** const copyTo,
                    const double** const copyFrom, const int len1,
                    const int len2)

/**************************************************************************
 *
 * mdp_copy_dbl_2:
 *
 * Copies one double array into another double array
 *
 * Input
 * -------------
 *  *copyFrom = Vector of values to be copied
 *  len1       = Length of the first array
 *  len2       = length of the second array
 *
 * Output
 * ------------
 * *copyTo   = array of values to receive the copy
 *
 **************************************************************************/
{
    if (len1 > 0 && len2 > 0) {
        size_t bytelen = len1 * len2 * sizeof(double);
        const double* const dFrom = copyFrom[0];
        double* const dTo = copyTo[0];
        (void) memcpy((void*)dTo, (const void*)dFrom, bytelen);
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_copy_int_1(int* const copyTo,
                    const int* const copyFrom, int len)

/**************************************************************************
 *
 * mdp_copy_int_1:
 *
 * Copies one int vector into int double vector
 *
 * Input
 * -------------
 *  *copyFrom = Vector of values to be copied
 *  len       = Length of the vector
 *
 * Output
 * ------------
 * *copyTo   = Vector of values to receive the copy
 *
 **************************************************************************/
{
    if (len > 0) {
        size_t bytelen = len * sizeof(int);
        (void) memcpy((void*)copyTo, (const void*)copyFrom, bytelen);
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_copy_int_2(int** const copyTo,
                    const int** const copyFrom, const int len1,
                    const int len2)

/**************************************************************************
 *
 * mdp_copy_int_2:
 *
 * Copies one 2D int array into another 2D int array
 *
 * Input
 * -------------
 *  *copyFrom = Vector of values to be copied
 *  len1       = Length of the first array
 *  len2       = Length of the second array
 *
 * Output
 * ------------
 * *copyTo   = array of values to receive the copy
 *
 **************************************************************************/
{
    if (len1 > 0 && len2 > 0) {
        size_t bytelen = len1 * len2 * sizeof(int);
        const int* const dFrom = copyFrom[0];
        int* const dTo = copyTo[0];
        (void) memcpy((void*)dTo, (const void*)dFrom, bytelen);
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_copy_ptr_1(void* const copyTo,
                    const void* const copyFrom, int len)

/**************************************************************************
 *
 * mdp_copy_ptr_1:
 *
 * Copies one ptr vector into another ptr vector
 *
 * Input
 * -------------
 *  *copyFrom = Vector of ptr values to be copied
 *  len       = Length of the vector
 *
 * Output
 * ------------
 * *copyTo   = Vector of values to receive the copy
 *
 **************************************************************************/
{
    if (len > 0) {
        size_t bytelen = len * sizeof(void*);
        (void) memcpy((void*)copyTo, (const void*)copyFrom, bytelen);
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void** mdp_dupl_ptr_1(const void* const copyFrom, int len)

/**************************************************************************
 *
 * mdp_dupl_ptr_1:
 *
 * duplicates one ptr vector into another ptr vector
 *
 * Input
 * -------------
 *  *copyFrom = Vector of ptr values to be copied
 *  len       = Length of the vector
 *
 * Output
 * ------------
 * *copyTo   = Vector of values to receive the copy
 *
 **************************************************************************/
{
    if (len > 0) {
        void** array = mdp_alloc_ptr_1(len);
        if (copyFrom) {
            size_t bytelen = len * sizeof(void*);
            (void) memcpy((void*)array, (const void*)copyFrom, bytelen);
        }
        return array;
    }
    return 0;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_init_dbl_1(double* const v, const double value, const int len)

/**************************************************************************
 *
 * mdp_init_dbl_1:
 *
 * Assigns a single value to a double vector
 *
 * Input
 * -------------
 *  v = Vector of values to be assigned
 *  value = value to assign with
 *  len       = Length of the vector
 *
 **************************************************************************/
{
    if (len > 0) {
        if (value == 0.0) {
            (void) memset((void*)v, 0, len * sizeof(double));
        } else {
            int m = len % 7;
            if (m != 0) {
                for (int i = 0; i < m; i++) {
                    v[i] = value;
                }
                if (len < 7) {
                    return;
                }
            }
            for (int i = m; i < len; i += 7) {
                v[i  ] = value;
                v[i+1] = value;
                v[i+2] = value;
                v[i+3] = value;
                v[i+4] = value;
                v[i+5] = value;
                v[i+6] = value;
            }
        }
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_zero_dbl_1(double* const v, const int len)

/**************************************************************************
 *
 * mdp_zero_dbl_1:
 *
 * Zeroes out a double vector (special form of mdp_allo_dbl_1())
 *
 * Input
 * -------------
 *  v = Vector of values to be set to zero
 *  len       = Length of the vector
 **************************************************************************/
{
    if (len > 0) {
        (void) memset((void*)v, 0, len * sizeof(double));
    }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_zero_int_1(int* const v, const int len)

/**************************************************************************
 *
 * mdp_zero_int_1:
 *
 * Zeroes out an int vector
 *
 * Input
 * -------------
 *  v = Vector of values to be set to zero
 *  len       = Length of the vector
 **************************************************************************/
{
    if (len > 0) {
        (void) memset((void*)v, 0, len * sizeof(int));
    }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_init_dbl_2(double** const v, double value, int len1, int len2)

/**************************************************************************
 *
 * mdp_init_dbl_2:
 *
 * Assigns a single value to a double matrix. Contiguous data for the
 * matrix is assumed.
 *
 * Input
 * -------------
 *  v = matrix of values to be assigned
 *  value = value to assign with
 *  len       = Length of the vector
 *
 **************************************************************************/
{
    int len = len1 * len2;
    if (len > 0 && len1 > 0 && v) {
        double* const dstart = v[0];
        if (value == 0.0) {
            size_t bytelen = len * sizeof(double);
            (void) memset((void*)dstart, 0, bytelen);
        } else {
            int m = len % 7;
            if (m != 0) {
                for (int i = 0; i < m; i++) {
                    dstart[i] = value;
                }
                if (len < 7) {
                    return;
                }
            }
            for (int i = m; i < len; i += 7) {
                dstart[i  ] = value;
                dstart[i+1] = value;
                dstart[i+2] = value;
                dstart[i+3] = value;
                dstart[i+4] = value;
                dstart[i+5] = value;
                dstart[i+6] = value;
            }
        }
    }
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

void mdp_init_int_1(int* const v, const int value, const int len)

/**************************************************************************
 *
 * mdp_init_int_1:
 *
 * Assigns a single value to an int vector
 *
 * Input
 * -------------
 *  v = Vector of values to be assigned
 *  value = value to assign with
 *  len       = Length of the vector
 *
 **************************************************************************/
{
    if (len > 0) {
        if (value == 0) {
            size_t bytelen = len * sizeof(int);
            (void) memset((void*)v, 0, bytelen);
        } else {
            int m = len % 7;
            if (m != 0) {
                for (int i = 0; i < m; i++) {
                    v[i] = value;
                }
                if (len < 7) {
                    return;
                }
            }
            for (int i = m; i < len; i += 7) {
                v[i  ] = value;
                v[i+1] = value;
                v[i+2] = value;
                v[i+3] = value;
                v[i+4] = value;
                v[i+5] = value;
                v[i+6] = value;
            }
        }
    }
}
}
/****************************************************************************/
/*                      END of mdp_allo.cpp                                 */
/****************************************************************************/
