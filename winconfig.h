/* config.h.  Generated from config.h.in by configure.  */
//
//  Run the 'configure' script to generate 'config.h' from this input file.
//
#ifndef CT_CONFIG_H
#define CT_CONFIG_H


//------------------------ Development flags ------------------//
//
// These flags turn on or off features that are still in 
// development and are not yet stable. 

#define DEV_EQUIL  

// Compile in additional debug printing where available.
// Note, the printing may need to be turned on via a switch.
// This just compiles in the code.
/* #undef DEBUG_MODE */

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


#define HAS_SUNDIALS 1


//-------- LAPACK / BLAS ---------

// Define if you are using LAPACK and BLAS from the Intel Math Kernel
// Library 
/* #undef HAVE_INTEL_MKL */

#define LAPACK_FTN_STRING_LEN_AT_END
#define LAPACK_NAMES_LOWERCASE
#define LAPACK_FTN_TRAILING_UNDERSCORE


//--------- operating system --------------------------------------

/* #undef DARWIN */
#define HAS_SSTREAM 1

// Cantera version
#define CANTERA_VERSION "1.7.0"

// Identify whether the operating system is cygwin's overlay of
// windows, with gcc being used as the compiler.
/* #undef CYGWIN */

// Identify whether the operating system is windows based, with
// microsoft vc++ being used as the compiler
#define WINMSVC 1

#define HAS_SUNDIALS 1
#define SUNDIALS_VERSION_23

//--------- Fonts for reaction path diagrams ----------------------
#define RXNPATH_FONT "Helvetica"

//--------------------- Python ------------------------------------
// This path to the python executable is created during
// Cantera's setup. It identifies the python executable 
// used to run Python to process .cti files. Note that this is only
// used if environment variable PYTHON_CMD is not set.
#define PYTHON_EXE "/usr/bin/python"

// If this is defined, the Cantera Python interface will use the
// Numeric package; otherwise, it will use numarray.
/* #undef HAS_NUMERIC */

// If this is defined, then python will not be assumed to be
// present to support conversions
/* #undef HAS_NO_PYTHON */

#define WITH_HTML_LOGS

// define STORE_MOLE_FRACTIONS if you want Cantera to internally
// represent the composition of a mixture as mole fractions. Usually
// the best choice.
#define STORE_MOLE_FRACTIONS

// define STORE_MASS_FRACTIONS if you want Cantera to internally
// represent the composition of a mixture as mass fractions. Usually
// results in slightly worse performance, but may not in all cases.
//#define STORE_MASS_FRACTIONS
/* #undef STORE_MASS_FRACTIONS */

#define WITH_METAL
#define WITH_STOICH_SUBSTANCE
#define WITH_PURE_FLUIDS
#define WITH_LATTICE_SOLID

#endif
