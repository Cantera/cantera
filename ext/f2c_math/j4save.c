/* j4save.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#ifdef _cpluscplus
extern "C" {
#endif
#include "f2c.h"

/* DECK J4SAVE */
integer j4save_(integer *iwhich, integer *ivalue, logical *iset)
{
    /* Initialized data */

    static integer iparam[9] = { 0,2,0,10,1,0,0,0,0 };

    /* System generated locals */
    integer ret_val;

/* ***BEGIN PROLOGUE  J4SAVE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Save or recall global variables needed by error */
/*            handling routines. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***TYPE      INTEGER (J4SAVE-I) */
/* ***KEYWORDS  ERROR MESSAGES, ERROR NUMBER, RECALL, SAVE, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/*     Abstract */
/*        J4SAVE saves and recalls several global variables needed */
/*        by the library error handling routines. */

/*     Description of Parameters */
/*      --Input-- */
/*        IWHICH - Index of item desired. */
/*                = 1 Refers to current error number. */
/*                = 2 Refers to current error control flag. */
/*                = 3 Refers to current unit number to which error */
/*                    messages are to be sent.  (0 means use standard.) */
/*                = 4 Refers to the maximum number of times any */
/*                     message is to be printed (as set by XERMAX). */
/*                = 5 Refers to the total number of units to which */
/*                     each error message is to be written. */
/*                = 6 Refers to the 2nd unit for error messages */
/*                = 7 Refers to the 3rd unit for error messages */
/*                = 8 Refers to the 4th unit for error messages */
/*                = 9 Refers to the 5th unit for error messages */
/*        IVALUE - The value to be set for the IWHICH-th parameter, */
/*                 if ISET is .TRUE. . */
/*        ISET   - If ISET=.TRUE., the IWHICH-th parameter will BE */
/*                 given the value, IVALUE.  If ISET=.FALSE., the */
/*                 IWHICH-th parameter will be unchanged, and IVALUE */
/*                 is a dummy parameter. */
/*      --Output-- */
/*        The (old) value of the IWHICH-th parameter will be returned */
/*        in the function value, J4SAVE. */

/* ***SEE ALSO  XERMSG */
/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   790801  DATE WRITTEN */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900205  Minor modifications to prologue.  (WRB) */
/*   900402  Added TYPE section.  (WRB) */
/*   910411  Added KEYWORDS section.  (WRB) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  J4SAVE */
/* ***FIRST EXECUTABLE STATEMENT  J4SAVE */
    ret_val = iparam[(0 + (0 + (*iwhich - 1 << 2))) / 4];
    if (*iset) {
	iparam[*iwhich - 1] = *ivalue;
    }
    return ret_val;
} /* j4save_ */

#ifdef _cpluscplus
}
#endif
