/* ../config.h.  Generated by configure.  */
//
//  Run the 'configure' script to generate 'config.h' from this input file.
//
#ifndef CT_CONFIG_H
#define CT_CONFIG_H


//------------------------ Development flags ------------------//
//
// These flags turn on or off features that are still in
// development and are not yet stable.

#undef DEV_EQUIL


//------------------------ Fortran settings -------------------//


// define types doublereal, integer, and ftnlen to match the
// corresponding Fortran data types on your system. The defaults
// are OK for most systems

typedef  double       doublereal;   // Fortran double precision
typedef  int          integer;      // Fortran integer
typedef  int          ftnlen;       // Fortran hidden string length type



// Fortran compilers pass character strings in argument lists by
// adding a hidden argement with the length of the string. Some
// compilers add the hidden length argument immediately after the
// CHARACTER variable being passed, while others put all of the hidden
// length arguments at the end of the argument list. Define this if
// the lengths are at the end of the argument list. This is usually the
// case for most unix Fortran compilers, but is (by default) false for
// Visual Fortran under Windows.
#define STRING_LEN_AT_END


// Define this if Fortran adds a trailing underscore to names in object files.
// For linux and most unix systems, this is the case.
#define FTN_TRAILING_UNDERSCORE


//-------- LAPACK / BLAS ---------

// Define if you are using LAPACK and BLAS from the Intel Math Kernel
// Library
/* #undef HAVE_INTEL_MKL */

#define LAPACK_FTN_STRING_LEN_AT_END
#define LAPACK_NAMES_LOWERCASE
#define LAPACK_FTN_TRAILING_UNDERSCORE


//--------- operating system --------------------------------------

// The configure script defines this if the operatiing system is Mac
// OS X, This used to add some Mac-specific directories to the default
// data file search path.
#define DARWIN 0
/* #undef HAS_SSTREAM */

// Identify whether the operating system is cygwin's overlay of
// windows, with gcc being used as the compiler.
/* #undef CYGWIN */

// Identify whether the operating system is windows based, with
// microsoft vc++ being used as the compiler
#define WINMSVC

//--------- Fonts for reaction path diagrams ----------------------
#define RXNPATH_FONT "Helvetica"

//--------------------- Python ------------------------------------
// If this is defined, the Cantera Python interface will use the
// Numeric package; otherwise, it will use numarray.
/* #define HAS_NUMERIC 1 */

#define INCL_PURE_FLUIDS 1

//--------------------- compile options ----------------------------
/* #define USE_PCH 1 */

#endif
