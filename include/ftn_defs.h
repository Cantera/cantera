
#ifndef CT_FTNDEFS_H
#define CT_FTNDEFS_H

// These definitions are required for Fortran/C mixed-language programming.


// Fortran compilers pass character strings in argument lists by
// adding a hidden argement with the length of the string. Some
// compilers add the hidden length argument immediately after the
// CHARACTER variable being passed, while others put all of the hidden
// length arguments at the end of the argument list. Define this if 
// the lengths are at the end of the argument list.
#define STRING_LEN_AT_END

// Define this if Fortran adds a trailing underscore to names in object files.
// For linux and most unix systems, this is the case.
#define FTN_TRAILING_UNDERSCORE

#endif
