/* config.h.  Generated from config.h.in by configure.  */
//
//  Run the 'preconfig' script to generate 'config.h' from this input file.
//
#ifndef CT_CONFIG_H
#define CT_CONFIG_H


//------------------------ Development flags ------------------//
//
// Compile in additional debug printing where available.
// Note, the printing may need to be turned on via a switch.
// This just compiles in the code.
/* #undef DEBUG_MODE */

// Compiling with PURIFY instrumentation
/* #undef PURIFY_MODE */

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
#define FTN_TRAILING_UNDERSCORE 1


#define HAS_SUNDIALS 1
/* #undef SUNDIALS_VERSION_22 */
#define SUNDIALS_VERSION_23 1

//-------- LAPACK / BLAS ---------

#define LAPACK_FTN_STRING_LEN_AT_END 1
#define LAPACK_NAMES_LOWERCASE 1
#define LAPACK_FTN_TRAILING_UNDERSCORE 1


//--------- operating system --------------------------------------

// The configure script defines this if the operatiing system is Mac
// OS X, This used to add some Mac-specific directories to the default
// data file search path.
/* #undef DARWIN */
#define HAS_SSTREAM 1

// Cantera version
#define CANTERA_VERSION "1.8.0"

// Identify whether the operating system is cygwin's overlay of
// windows, with gcc being used as the compiler.
/* #undef CYGWIN */

// Identify whether the operating system is windows based, with
// microsoft vc++ being used as the compiler
#define WINMSVC 1

// Identify whether the operating system is solaris 
// with a native compiler 
/* #undef SOLARIS */

//--------- Fonts for reaction path diagrams ----------------------
#define RXNPATH_FONT "Helvetica"

//---------- C++ Compiler Variations ------------------------------

// This define is needed to account for the variability for how
// static variables in templated classes are defined. Right now
// this is only turned on for the SunPro compiler on solaris.
// in that system , you need to declare the static storage variable.
// with the following line in the include file
//
//    template<class M> Cabinet<M>* Cabinet<M>::__storage;
//
// Note, on other systems that declaration is treated as a definition
// and this leads to multiple defines at link time
/* #undef NEEDS_GENERIC_TEMPL_STATIC_DECL */

//--------------------- Python ------------------------------------
// This path to the python executable is created during
// Cantera's setup. It identifies the python executable 
// used to run Python to process .cti files. Note that this is only
// used if environment variable PYTHON_CMD is not set.
#define PYTHON_EXE "c:/python26/python.exe"

// If this is defined, the Cantera Python interface will use the
// Numeric package
/* #undef HAS_NUMERIC */

// If this is defined, the Cantera Python interface will use the
// numarray package
/* #undef HAS_NUMARRAY */

// If this is defined, the Cantera Python interface will use the
// numpy package
#define HAS_NUMPY 1

// If this is defined, then python will not be assumed to be
// present to support conversions
/* #undef HAS_NO_PYTHON */

//--------------------- Cantera ----------------------------------- 
// This is the data pathway used to locate the top of the 
// build directory.
/* #undef CANTERA_ROOT */

// This data pathway is used to locate a directory where datafiles
// are to be found. Note, the local directory is always searched
// as well. 
#define CANTERA_DATA "C:/cygwin/usr/local/cantera/data"


#define WITH_HTML_LOGS 1

// define STORE_MOLE_FRACTIONS if you want Cantera to internally
// represent the composition of a mixture as mole fractions. Usually
// the best choice.
#define STORE_MOLE_FRACTIONS

// define STORE_MASS_FRACTIONS if you want Cantera to internally
// represent the composition of a mixture as mass fractions. Usually
// results in slightly worse performance, but may not in all cases.
//#define STORE_MASS_FRACTIONS
/* #undef STORE_MASS_FRACTIONS */

//--------------------- compile options ----------------------------
/* #undef USE_PCH */

/* #undef THREAD_SAFE_CANTERA */

//--------------------- optional phase models ----------------------
//    This define indicates the enabling of the inclusion of
//    accurate liquid/vapor equations
//    of state for several fluids, including water, nitrogen, hydrogen,
//    oxygen, methane, andd HFC-134a.
#define INCL_PURE_FLUIDS 1
#define WITH_PURE_FLUIDS 1

#define WITH_LATTICE_SOLID 1
#define WITH_METAL 1
#define WITH_STOICH_SUBSTANCE 1
//    Enable expanded thermodynamic capabilities, adding
//    ideal solid solutions
#define WITH_IDEAL_SOLUTIONS 1
//    Enable expanded electrochemistry capabilities, include thermo
//    models for electrolyte solutions.
#define WITH_ELECTROLYTES 1

/* #undef WITH_PRIME */

//    Enable the VCS NonIdeal equilibrium solver. This is
//    accessed by specifying the solver=2 option
#define WITH_VCSNONIDEAL 1

//-------------- Optional Cantera Capabilities ----------------------

//    Enable sensitivity analysis via changing H298 directly
//    for species
/* #undef H298MODIFY_CAPABILITY */

#endif
