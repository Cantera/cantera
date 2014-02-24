/* xersve.f -- translated by f2c (version 20030320).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;

/* DECK XERSVE */
/* Subroutine */ int xersve_(char *librar, char *subrou, char *messg, integer 
	*kflag, integer *nerr, integer *level, integer *icount, ftnlen 
	librar_len, ftnlen subrou_len, ftnlen messg_len)
{
    /* Initialized data */

    static integer kountx = 0;
    static integer nmsg = 0;

    /* Format strings */
    static char fmt_9000[] = "(\0020          ERROR MESSAGE SUMMARY\002/\002"
	    " LIBRARY    SUBROUTINE MESSAGE START             NERR\002,\002  "
	    "   LEVEL     COUNT\002)";
    static char fmt_9010[] = "(1x,a,3x,a,3x,a,3i10)";
    static char fmt_9020[] = "(\0020OTHER ERRORS NOT INDIVIDUALLY TABULATED "
	    "= \002,i10)";
    static char fmt_9030[] = "(1x)";

    /* System generated locals */
    integer i__1, i__2;

    /* Builtin functions */
    integer s_wsfe(cilist *), e_wsfe(void), do_fio(integer *, char *, ftnlen);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);
    integer s_cmp(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__;
    char lib[8], mes[20], sub[8];
    integer lun[5], iunit, kunit, nunit;
    static integer kount[10];
    extern integer i1mach_(integer *);
    static char libtab[8*10], mestab[20*10];
    static integer nertab[10], levtab[10];
    static char subtab[8*10];
    extern /* Subroutine */ int xgetua_(integer *, integer *);

    /* Fortran I/O blocks */
    static cilist io___7 = { 0, 0, 0, fmt_9000, 0 };
    static cilist io___9 = { 0, 0, 0, fmt_9010, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_9020, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_9030, 0 };


/* ***BEGIN PROLOGUE  XERSVE */
/* ***SUBSIDIARY */
/* ***PURPOSE  Record that an error has occurred. */
/* ***LIBRARY   SLATEC (XERROR) */
/* ***CATEGORY  R3 */
/* ***TYPE      ALL (XERSVE-A) */
/* ***KEYWORDS  ERROR, XERROR */
/* ***AUTHOR  Jones, R. E., (SNLA) */
/* ***DESCRIPTION */

/* *Usage: */

/*        INTEGER  KFLAG, NERR, LEVEL, ICOUNT */
/*        CHARACTER * (len) LIBRAR, SUBROU, MESSG */

/*        CALL XERSVE (LIBRAR, SUBROU, MESSG, KFLAG, NERR, LEVEL, ICOUNT) */

/* *Arguments: */

/*        LIBRAR :IN    is the library that the message is from. */
/*        SUBROU :IN    is the subroutine that the message is from. */
/*        MESSG  :IN    is the message to be saved. */
/*        KFLAG  :IN    indicates the action to be performed. */
/*                      when KFLAG > 0, the message in MESSG is saved. */
/*                      when KFLAG=0 the tables will be dumped and */
/*                      cleared. */
/*                      when KFLAG < 0, the tables will be dumped and */
/*                      not cleared. */
/*        NERR   :IN    is the error number. */
/*        LEVEL  :IN    is the error severity. */
/*        ICOUNT :OUT   the number of times this message has been seen, */
/*                      or zero if the table has overflowed and does not */
/*                      contain this message specifically.  When KFLAG=0, */
/*                      ICOUNT will not be altered. */

/* *Description: */

/*   Record that this error occurred and possibly dump and clear the */
/*   tables. */

/* ***REFERENCES  R. E. Jones and D. K. Kahaner, XERROR, the SLATEC */
/*                 Error-handling Package, SAND82-0800, Sandia */
/*                 Laboratories, 1982. */
/* ***ROUTINES CALLED  I1MACH, XGETUA */
/* ***REVISION HISTORY  (YYMMDD) */
/*   800319  DATE WRITTEN */
/*   861211  REVISION DATE from Version 3.2 */
/*   891214  Prologue converted to Version 4.0 format.  (BAB) */
/*   900413  Routine modified to remove reference to KFLAG.  (WRB) */
/*   900510  Changed to add LIBRARY NAME and SUBROUTINE to calling */
/*           sequence, use IF-THEN-ELSE, make number of saved entries */
/*           easily changeable, changed routine name from XERSAV to */
/*           XERSVE.  (RWC) */
/*   910626  Added LIBTAB and SUBTAB to SAVE statement.  (BKS) */
/*   920501  Reformatted the REFERENCES section.  (WRB) */
/* ***END PROLOGUE  XERSVE */
/* ***FIRST EXECUTABLE STATEMENT  XERSVE */

    if (*kflag <= 0) {

/*        Dump the table. */

	if (nmsg == 0) {
	    return 0;
	}

/*        Print to each unit. */

	xgetua_(lun, &nunit);
	i__1 = nunit;
	for (kunit = 1; kunit <= i__1; ++kunit) {
	    iunit = lun[kunit - 1];
	    if (iunit == 0) {
		iunit = i1mach_(&c__4);
	    }

/*           Print the table header. */

	    io___7.ciunit = iunit;
	    s_wsfe(&io___7);
	    e_wsfe();

/*           Print body of table. */

	    i__2 = nmsg;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		io___9.ciunit = iunit;
		s_wsfe(&io___9);
		do_fio(&c__1, libtab + ((i__ - 1) << 3), (ftnlen)8);
		do_fio(&c__1, subtab + ((i__ - 1) << 3), (ftnlen)8);
		do_fio(&c__1, mestab + (i__ - 1) * 20, (ftnlen)20);
		do_fio(&c__1, (char *)&nertab[i__ - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&levtab[i__ - 1], (ftnlen)sizeof(
			integer));
		do_fio(&c__1, (char *)&kount[i__ - 1], (ftnlen)sizeof(integer)
			);
		e_wsfe();
/* L10: */
	    }

/*           Print number of other errors. */

	    if (kountx != 0) {
		io___16.ciunit = iunit;
		s_wsfe(&io___16);
		do_fio(&c__1, (char *)&kountx, (ftnlen)sizeof(integer));
		e_wsfe();
	    }
	    io___17.ciunit = iunit;
	    s_wsfe(&io___17);
	    e_wsfe();
/* L20: */
	}

/*        Clear the error tables. */

	if (*kflag == 0) {
	    nmsg = 0;
	    kountx = 0;
	}
    } else {

/*        PROCESS A MESSAGE... */
/*        SEARCH FOR THIS MESSG, OR ELSE AN EMPTY SLOT FOR THIS MESSG, */
/*        OR ELSE DETERMINE THAT THE ERROR TABLE IS FULL. */

	s_copy(lib, librar, (ftnlen)8, librar_len);
	s_copy(sub, subrou, (ftnlen)8, subrou_len);
	s_copy(mes, messg, (ftnlen)20, messg_len);
	i__1 = nmsg;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (s_cmp(lib, libtab + ((i__ - 1) << 3), (ftnlen)8, (ftnlen)8) ==
		    0 && s_cmp(sub, subtab + ((i__ - 1) << 3), (ftnlen)8, (
		    ftnlen)8) == 0 && s_cmp(mes, mestab + (i__ - 1) * 20, (
		    ftnlen)20, (ftnlen)20) == 0 && *nerr == nertab[i__ - 1] &&
		     *level == levtab[i__ - 1]) {
		++kount[i__ - 1];
		*icount = kount[i__ - 1];
		return 0;
	    }
/* L30: */
	}

	if (nmsg < 10) {

/*           Empty slot found for new message. */

	    ++nmsg;
	    s_copy(libtab + ((i__ - 1) << 3), lib, (ftnlen)8, (ftnlen)8);
	    s_copy(subtab + ((i__ - 1) << 3), sub, (ftnlen)8, (ftnlen)8);
	    s_copy(mestab + (i__ - 1) * 20, mes, (ftnlen)20, (ftnlen)20);
	    nertab[i__ - 1] = *nerr;
	    levtab[i__ - 1] = *level;
	    kount[i__ - 1] = 1;
	    *icount = 1;
	} else {

/*           Table is full. */

	    ++kountx;
	    *icount = 0;
	}
    }
    return 0;

/*     Formats. */

} /* xersve_ */

