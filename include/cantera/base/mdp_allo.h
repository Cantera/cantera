/**
 * @file mdp_allo.h
 *   Declarations for Multi Dimensional Pointer (mdp) malloc routines, which
 *   allow for dimensioning of arbitrarily dimensioned pointer arrays using
 *   one call.
 */

/*
 * Copyright 2004 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * See file License.txt for licensing information.
 */

#ifndef MDP_ALLO_H
#define MDP_ALLO_H

#include <stdexcept>

/*
 * Include the header here in order to pick up size_t definition
 */
#include <cstring>

/**
 *  The mdp routines are extremely lightweight and fast fortran compatibile
 *  malloc routines for allocating multiple dimensioned arrays of doubles
 *  ints, char, and pointers using a single call. These routines don't
 *  use the std C+ lib routines.
 *
 *  All calls are essentially wrappers around the routine mdp_alloc_array()
 *  which allocates multidimensioned arrays. The arrays contain room for
 *  the data and the pointer information that is used to access data
 *  in the object.
 *
 *  One convention that is always used is that a pointer that is
 *  not malloced always has a value of zero. If the pointer is nonnull,
 *  then it may be freed, always (and vica-versa).
 *
 *  Where possible, the low leve routines
 *   memcpy and memset are used to copy or zero memory.
 *
 *  No array bounds checking is ever done within these routines. buyer beware.
 *  The bounds of arrays are not carried with the array object, ever.
 *  Bounds of arrays are input via the parameter list. Normally,
 *  for applications this is not a problem, since the application
 *  knows what the array dimensions are.
 *
 *  There are several other general principles.
 *
 *     If an allocation size is set to 0, then the actual
 *     allocation size is set to 1 within the program.
 *     Something is always allocated whenever a call is made
 *     to an mdp routine.
 *
 *     All checks for allocations are always checked for success.
 *     If a failure is found, thene mdp_alloc_eh() is called
 *     for disposition of the error condition.
 *
 *   Error handling behavior is set by the MDP_ALLO_errorOption external
 *   int. The default error behavior is to print an error message to stderr, and
 *   then throw an exception that inherited from std::exception. Usually
 *   std::bad_alloc() is thrown whenever there is a problem.
 *
 */
namespace mdp
{

/**
 *  If we have array_alloc() from another Sandia program, we will not use
 *  the one from this mdp_array_alloc. Instead we will redefine the names
 */
#ifdef HAVE_ARRAY_ALLOC
#  define mdp_array_alloc array_alloc
#  define mdp_safe_free   safe_free
#endif

/*!
 *  MDP_INT_NOINIT is a poor man's way of specifying whether a value should be
 *  initialized. These are seldom used numbers which can be used in place
 *  of real ints and dbls to indicate that initialization shouldn't take
 *  place.
 */
#define MDP_INT_NOINIT          -68361

/*!
 *  MDP_DBL_NOINIT is a poor man's way of specifying whether a value should be
 *  initialized. These are seldom used numbers which can be used in place
 *  of real ints and dbls to indicate that initialization shouldn't take
 *  place.
 */
#define MDP_DBL_NOINIT          -1.241E11

/*!
 *      Error Handling
 *         7 print and exit
 *         6 exit
 *         5 print and create a divide by zero for stack trace analysis.
 *         4 create a divide by zero for stack analysis trace
 *         3 print a message and throw the std::bad_alloc() exception.
 *         2 throw the std::bad_alloc() exception and be quiet.
 *         1 print a message and return from package with the NULL pointer
 *         0 Keep completely silent about the matter and return with
 *           a null pointer.
 *
 *   -> Right now, the only way to change this option is to right here
 *
 *   The default is to set it to 3.
 */
extern int MDP_ALLO_errorOption;


/****************************************************************************/
/*
 * Externals that should be set by the calling program.
 * These are only used for debugging purposes.
 */
#ifdef MDP_MPDEBUGIO
extern int MDP_MP_Nprocs;
extern int MDP_MP_myproc;
#endif

/****************************************************************************/

#define mdp_alloc_struct(x, num) (x *) mdp_array_alloc(1, (num), sizeof(x))


/* function declarations for dynamic array allocation */

//! allocates multidimensional pointer arrays of arbitrary length
//! via a single malloc call
/*!
 *  The first dimension is the number of dimensions in the allocation
 */
extern double* mdp_array_alloc(int numdim, ...);

//! Free a vector and set its value to 0
/*!
 *  This function carries out the following operation
 *  @code
 *    free(*hndVec);
 *    *hndVec = 0;
 *  @endcode
 *
 * @param hndVec This is the address of the pointer, expressed as
 *               a void **
 *
 * Note, a key idea behind the mdp suite of routines is that
 * a pointer that can be malloced is either malloced or its value
 * is 0. This routine enforces this convention.
 */
extern void mdp_safe_free(void** hndVec);

//! Allocate a vector of integers
/*!
 *  The vector is initialized, unless the default int value is set
 * to MDP_INT_NOINIT
 *
 * @param len  Length of the vector
 * @param defval Default value for the int, defaults to MDP_INT_NOINIT
 *
 * @return returns a pointer to the vector
 */
extern int* mdp_alloc_int_1(int len, const int defval = MDP_INT_NOINIT);

//! Allocate a vector of integers, potentially freeing memory first
/*!
 *  The vector is initialized, unless the default int value is set
 * to MDP_INT_NOINIT
 *
 * Input
 * --------------
 * @param array_hdl Previous value of pointer. If non-NULL will try
 *                to free the memory at this address before doing
 *                a new alloc
 * @param len  Length of the vector
 * @param defval Default value for the int, defaults to MDP_INT_NOINIT
 *
 * Output
 * ---------
 * @return  *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 */
extern void mdp_safe_alloc_int_1(int** array_hdl, int len,
                                 const int defval = MDP_INT_NOINIT);

//!    Reallocates a one dimensional array of ints, copying old
//!    information to the new array
/*!
 *    Reallocates a one dimensional array of ints.
 *    This routine always allocates space for at least one int.
 *    Calls the smalloc() routine to ensure that all malloc
 *    calls go through one location. This routine will then copy
 *    the pertinent information from the old array to the
 *    new array.
 *
 *     NewArray[0:old_len-1] = OldArray[0:old_len-1];
 *     NewArray[old_len:new_len-1] = defval;
 *
 * Input
 * --------------
 * @param array_hdl Previous value of pointer. If non-NULL will try
 *                to free the memory at this address before doing
 *                a new alloc
 * @param new_len  New Length of the vector
 * @param old_len  New Length of the vector
 * @param defval Default value for the int, defaults to MDP_INT_NOINIT
 *
 * Output
 * ---------
 * @return  *array_hdl = This value is initialized to the correct address
 *                     of the array.
 *                     A NULL value in the position indicates an error.
 */
extern void mdp_realloc_int_1(int** array_hdl, int new_len, int old_len,
                              const int defval = MDP_INT_NOINIT);


//! Allocate a 2D matrix of integers
/*!
 *  The matrix is initialized, unless the default int value is set
 * to MDP_INT_NOINIT, which is the default.
 *
 *     matrix[len1][len2]
 *
 *  All int data entries are contiguous. Therefore, it may
 *  be used as input into BLAS matrix function calls.
 *  This can be considered to be in fortran order format with
 *  len2 as the number of rows, and len1 as the number of columns.
 *
 *  matrix[jcol] refers to the jcol column of the matrix.
 *  Therefore, matrix[0] is a pointer to the beginning of the
 *  data portion of the structure.
 *  The structure will have len1 pointers at the beginning
 *  that holds pointers into the top of the columns of the
 *  contiguous data.
 *
 *  The entire structure may be deallocated via one free call.
 *
 * @param len1  Outer Length of the vector
 * @param len2  Inner length of the matrix
 * @param defval Default value for the int, defaults to MDP_INT_NOINIT
 *
 * @return returns a pointer to the matrix
 */
extern int** mdp_alloc_int_2(int len1, int len2,
                             const int defval = MDP_INT_NOINIT);


//! Allocate and initialize a one dimensional array of doubles.
/*!
 *   As per the convention in mdp, this routine always initializes
 *   at least one slot.
 *
 *      @param nvalues    Length of the array. If this number is
 *                        less than one, it is set to one. Therefore,
 *                        This routine always initializes at least
 *                        one double.
 *      @param val        initialization value. Set it to the
 *                        constant MDP_DBL_NOINIT if you don't
 *                        want any initialization. memset() is
 *                        used for zero initialization for fast
 *                        execution speed.
 *
 *     @return    Pointer to the initialized array of doubles
 *                Failures are indicated by returning the NULL pointer.
 */
extern double* mdp_alloc_dbl_1(int nvalues, const double val=MDP_DBL_NOINIT);

//!    Allocates and/or initializes a one dimensional array of doubles.
/*!
 *  This routine will free any old memory that was located at that
 *  position, before it will allocate a new vector.
 *
 *     @param hndVec     Previous value of pointer. If non-NULL will try
 *                       to free the memory at this address. On output,
 *                       this value is initialized to the correct address
 *                       of the array.  A NULL value in the position
 *                       indicates an error.
 *     @param nvalues    Length of the array
 *     @param val        initialization value
 *
 */
extern void mdp_safe_alloc_dbl_1(double** hndVec, int nvalues,
                                 const double val=MDP_DBL_NOINIT);

//! Reallocate a vector of doubles possibly retaining a subset of values
/*!
 *  Reallocates the array and sets:
 *
 *  (*hndVec)[0:oldLen-1]    = oldVec[0:oldLen-1]
 *  (*hndVec)[oldLen:newLen] = defVal
 *
 * Input
 * ------
 *   @param newLen  New Length of the vector
 *
 * Output
 * -------
 *   @param oldLen  Old  Length of the vector
 *
 */
extern void mdp_realloc_dbl_1(double** hndVec, int newLen, int oldLen,
                              const double defVal=MDP_DBL_NOINIT);

//!    Allocate and initialize a two dimensional array of doubles.
/*!
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
 *      @param  ndim1     Length of the first dimension of the array
 *      @param  ndim2     Length of the second dimension of the array
 *      @param  val       Initialization value
 *
 *
 *     @return       Pointer to the initialized array of doubles.
 *                   Failures are indicated by returning the NULL pointer.
 */
extern double** mdp_alloc_dbl_2(int ndim1, int ndim2, const double val);

//!    Allocate and initialize a two dimensional array of doubles.
/*!
 *    Allocate a two dimensional array of doubles. The array is in
 *    fortran order and can be accessed via the following form:
 *
 *          (*arrayHndl)[ndim1][ndim2]
 *
 *   Note, ndim2 is the inner dimension. i.e., the array is
 *   in column ordering.
 *
 *    Input
 *    -------
 *      @param  ndim1     Length of the first dimension of the array
 *      @param  ndim2     Length of the second dimension of the array
 *      @param  val       Initialization value
 *      @param  arrayHndl Handle to the array. If nonnull, the array
 *                        is first freed. Failures are indicated
 *                        by returning the NULL pointer.
 */
extern void mdp_safe_alloc_dbl_2(double** *arrayHndl, int ndim1, int ndim2,
                                 const double val = MDP_DBL_NOINIT);

//!  Reallocates a two dimensional array of doubles to a new set of
//!  dimensions, copying the old results into the new array.
/*!
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
 *   @param hndArray     = Pointer to the global variable that
 *                         holds the old and (eventually new)
 *                         address of the array of doubles to be reallocated
 *   @param ndim1        = First dimension of the new array
 *   @param ndim2        = Second dimension of the new array
 *   @param ndim1Old     = First dimension of the old array
 *   @param ndim2Old     = Second dimension of the old array
 *   @param defVal       = Default fill value.
 *
 *   Output
 *   --------------
 *   The resulting vector looks like this:
 *
 *  (*hndArray)[0:ndim1Old-1][0:ndim2Old-1] = oldVec[0:ndim1Old-1][0:ndim2Old-1]
 *  (*hndArray)[ndim1Old:ndim1][ndim2Old:ndim2] = defVal
 *
 */
extern void mdp_realloc_dbl_2(double** * hndArray, int ndim1, int ndim2,
                              int ndim1Old, int ndim2Old,
                              const double defVal=MDP_DBL_NOINIT);

//!    Allocate and initialize a one dimensional array of characters.
/*!
 *  The array is always initialized.
 *
 *    Input
 *    -------
 *      @param  nvalues   Length of the array
 *      @param  val       initialization value. defaults to the NULL char
 *
 *    @return    Pointer to the initialized character array
 *               Failures are indicated by returning the NULL pointer.
 */
extern char*   mdp_alloc_char_1(int nvalues, const char val = '\0');

//!    Allocate and initialize a one dimensional array of characters,
//!    deallocating the space before hand.
/*!
 *  This routine will free any old memory that was located at that
 *  position, before it will allocate a new vector.
 *
 *  The array is always initialized.
 *
 *      @param arrayHnd   Pointer to the global variable that
 *                         holds the old and (eventually new)
 *                         address of the array of char to be reallocated
 *      @param  nvalues   Length of the array
 *      @param  val       initialization value. defaults to the NULL char
 *
 *    @return    Pointer to the initialized character array
 *               Failures are indicated by returning the NULL pointer.
 */
extern void mdp_safe_alloc_char_1(char** arrayHnd, int nvalues,
                                  const char val = '\0');

//!    Allocate and initialize a vector of fixed-length
//!    strings. Each string is initialized to the NULL string.
/*!
 *      @param  numStrings   Number of strings
 *      @param  lenString    Length of each string including the trailing null
 *                           character
 *
 *     @return      Value is initialized to the correct address
 *                           of the new array on exit.
 *                           A NULL value in the position indicates an error.
 */
extern char**  mdp_alloc_VecFixedStrings(int numStrings, int lenString);

//!  Allocate and initialize an array of strings of fixed length
/*!
 *      @param  numStrings   Number of strings
 *      @param  lenString    Length of each string including the trailing null
 *                           character
 *
 *     @param array_hdl      Value is initialized to the correct address
 *                           of the new array on exit.
 *                           A NULL value in the position indicates an error.
 *                           If non-NULL on entry the routine will first
 *                           free the memory at the address.
 */
extern void mdp_safe_alloc_VecFixedStrings(char** *arrayHnd, int numStrings,
        int lenString);

//!   Reallocate and initialize a vector of fixed-length strings.
/*!
 *  Each new string is initialized to the NULL string.
 *  Old strings are copied.
 *
 *  @param array_hdl        The pointer to the char ** location holding
 *                          the data to be reallocated.
 *  @param numStrings       Number of strings
 *  @param numOldStrings    Number of old strings
 *  @param lenString        Length of each string including the trailing null
 *                          character
 */
extern void mdp_realloc_VecFixedStrings(char** *array_hdl, int numStrings,
                                        int numOldStrings, int lenString);

//! Allocate and initialize a vector of pointers of type pointer to void.
/*!
 *   All pointers are initialized to the NULL  value.
 *
 *   @param  numPointers  Number of pointers
 *
 *   @return       This value is initialized to the correct address of the vector.
 *                 A NULL value in the position indicates an error.
 */
extern void**  mdp_alloc_ptr_1(int numPointers);

//! Allocate and initialize a vector of pointers of type pointer to void.
/*!
 *    All pointers are initialized to the NULL value.
 *
 *    @param numPointers  Number of pointers
 *    @param array_hdl    This value is initialized to the correct address
 *                        of the array.
 *                        A NULL value in the position indicates an error.
 *                        Previous value of pointer. If non-NULL will try
 *                        to free the memory at this address.
 */
extern void mdp_safe_alloc_ptr_1(void** *array_hnd, int numPointers);

//!    Reallocate and initialize a vector of pointers
/*!
 *  All old pointers are copied
 *  Each new pointer not associated with an old pointer is
 *  initialized to NULL.
 *
 *   @param array_hdl    The pointer to the char ** location holding
 *                       the data to be reallocated.
 *   @param numLen       Number of new pointers
 *   @param numOldLen    Number of old pointers
 */
extern void  mdp_realloc_ptr_1(void** *array_hdl, int numLen, int numOldLen);

//! Copies one ptr vector into another ptr vector
/*!
 *
 *  @param copyFrom      Vector of ptr values to be copied
 *  @param len           Length of the vector
 *
 *  @param copyTo        Vector of values to receive the copy
 */
extern void mdp_copy_ptr_1(void* const copyTo,
                           const void* const copyFrom, const int len);

//!  Duplicates one ptr vector into another ptr vector
/*!
 *  Mallocs a copy of one vector of pointers and returns the pointer
 *  to the copy.
 *
 * Input
 * -------------
 * @param  *copyFrom   Vector of ptr values to be copied
 * @param   len        Length of the vector
 *
 * Output
 * ------------
 * @return  Vector of values to receive the copy
 */
extern void**  mdp_dupl_ptr_1(const void* const copyFrom, int len);

//!   Copies an array of string vectors from one char ** vector to
//!   another
/*!
 *  The space must have already been allocated within copyFrom and copyTo
 *  arrays. Overwrites are prevented by the proper application of the
 *  variable maxLenString. Strings are forced to be null-terminated
 *  Therefore copyTo[maxLenString-1] = '/0'
 *
 *    Input
 *    -------
 *     @param   copyFrom       vector of C strings. It should be null terminated
 *     @param   numStrings     number of strings
 *     @param   maxLenString   maximum of the size of the string arrays,
 *                             copyTo and copyFrom. This is used as the
 *                             argument to strncpy() function.
 *
 *    Output
 *    ------
 *     @param  copyTo          vector of strings
 *
 */
extern void mdp_copy_VecFixedStrings(char** const copyTo,
                                     const char** const copyFrom,
                                     int numStrings, size_t maxLenString);

//!    Allocates space for and copies a string
/*!
 *
 *    @param  copyFrom   null terminated string. If NULL is supplied, then
 *                       nothing is malloced and a NULL value is returned.
 *
 *    @return    This value is initialized to the correct address of the array.
 *               A NULL value in the position either indicates an error, or
 *               that the original pointer to the string was NULL.
 */
extern char* mdp_copy_string(const char* const copyFrom);

//!  Allocates space for and copies a string
/*!
 *    @param stringHnd  Previous value of pointer. If non-NULL will try
 *                      to free the memory at this address.
 *
 *    @param  copyFrom   null terminated string. If NULL is supplied, then
 *                       nothing is malloced and a NULL value is returned.
 *
 *    @return    This value is initialized to the correct address of the array.
 *               A NULL value in the position either indicates an error, or
 *               that the original pointer to the string was NULL.
 */
extern void mdp_safe_copy_string(char** stringHnd, const char* const copyFrom);

//! Copy a double vector to a double vector
/*!
 *  copyTo[len] = copyFrom[len]
 *
 * Input
 * -------
 * @param copyFrom Vector to  copy ( length >= len)
 * @param len      Length of the copy
 *
 * Output
 * -------
 * @param copyTo   Vector to receive the copy ( length >= len)
 */
extern void mdp_copy_dbl_1(double* const copyTo,
                           const double* const copyFrom,
                           const int len);

//! Copy a double array to a double array
/*!
 *   This routine carries out a straight copy on the effective 1D description
 *   of each of the arrays. It copies
 *   the first len1*len2 doubless stored within copyFrom into the
 *   the first len1*len2 double slots in copyTo. It does not account
 *   for the actual dimensions of the two arrays.
 *
 * Input
 * --------
 * @param copyFrom Vector to  copy ( length >= len1 * len2)
 * @param len1      Length of the first dimension
 * @param len2      Length of the second dimension
 *
 * Output
 * ----------
 * @param copyTo   Array to receive the copy ( length >= len1 * len2)
 */
extern void  mdp_copy_dbl_2(double** const copyTo, const double** const copyFrom,
                            const int len1, const int len2);

//! Copies one int vector into one int vector
/*!
 * Input
 * -------------
 *  @param  copyFrom    Vector of values to be copied
 *  @param  len         Length of the vector
 *
 * Output
 * ------------
 * *copyTo   = Vector of values to receive the copy
 */
extern void mdp_copy_int_1(int* const copyTo, const int* const copyFrom,
                           const int len);

//! Copies one 2D int array into another 2D int array
/*!
 *   This routine carries out a straight copy. Actually it copies
 *   the first len1*len2 ints stored within copyFrom into the
 *   the first len1*len2 int slots in copyTo. It does not account
 *   for the actual dimensions of the two arrays.
 *
 *  @param  copyFrom     Vector of values to be copied
 *  @param  len1         Length of the first array
 *  @param  len2         Length of the second array
 * Output
 * ------------
 * @param  copyTo        Vector of values to receive the copy
 */
extern void mdp_copy_int_2(int** const copyTo,  const int** const copyFrom,
                           const int len1, const int len2);

//! Assigns a single value to a double vector
/*!
 *  @param v        Vector of values to be assigned
 *  @param value    Value to assign with
 *  @param len      Length of the vector
 */
extern void mdp_init_dbl_1(double* const v, const double value, const int len);

//! Zeroes a double vector
/*!
 *  @param v = Vector of values to be assigned
 *  @param len       = Length of the vector
 */
extern void mdp_zero_dbl_1(double* const v , const int len);

//! Zeroes an int vector
/*!
 *  @param v = Vector of values to be assigned
 *  @param len       = Length of the vector
 */
extern void mdp_zero_int_1(int* const v , const int len);

//!  Assigns a single value to a double matrix. Contiguous data for the
//!  matrix is assumed.
/*!
 * Input
 * -------------
 *  @param v          matrix of values to be assigned
 *  @param value      value to assign with
 *  @param len1       Length one of the vector
 *  @param len2       length two of the vector
 */
extern void mdp_init_dbl_2(double** const v, const double value,
                           const int len1, const int len2);

//! Assigns a single value to an int vector
/*!
 *  @param v          Vector of values to be assigned
 *  @param value      Value to assign with
 *  @param len        Length of the vector
 */
extern void mdp_init_int_1(int* const v, const int value, const int len);


/*
 *        Utility routines to check that a number is finite
 */

//! Utility routine to check to see that a number is neither zero
//! nor indefinite.
/*!
 *  This check can be used before using the number in a denominator.
 *
 *  @param tmp number to be checked
 */
extern void checkZeroFinite(const double tmp);

//! Utility routine to check to see that a number is finite.
/*!
 *  @param tmp number to be checked
 */
extern void checkFinite(const double tmp);

//! Utility routine to link checkFinte() to fortran program
/*!
 *  This routine is accessible from fortran, usually
 *
 * @param tmp Pointer to the number to check
 *
 * @todo link it into the usual way Cantera handles Fortran calls
 */
extern "C" void checkfinite_(double* tmp);

//! utility routine to check that a double stays bounded
/*!
 *   This routine checks to see if a number stays bounded. The absolute
 *   value of the number is required to stay below the trigger.
 *
 * @param tmp     Number to be checked
 * @param trigger bounds on the number. Defaults to 1.0E20
 */
extern void checkMagnitude(const double tmp, const double trigger = 1.0E20);

}  /* end of mdp namespace */
/****************************************************************************/
#endif
/****************************************************************************/

