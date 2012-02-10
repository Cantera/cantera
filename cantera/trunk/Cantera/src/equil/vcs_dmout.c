/* dmout.f -- translated by f2c (version 20031025).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Table of constant values */

static integer c__2 = 2;
static integer c__1 = 1;

/* DECK DMOUT */
/* Subroutine */
int dmout_(integer* m, integer* n, integer* lda, doublereal *
           a, char* ifmt, integer* idigit, ftnlen ifmt_len)
{
    /* Initialized data */

    static char icol[3] = "COL";

    /* Format strings */
    static char fmt_1010[] = "(10x,10(4x,a,i4,1x))";
    static char fmt_1009[] = "(1x,\002ROW\002,i4,2x,1p10d12.3)";
    static char fmt_1000[] = "(10x,8(5x,a,i4,2x))";
    static char fmt_1004[] = "(1x,\002ROW\002,i4,2x,1p8d14.5)";
    static char fmt_1001[] = "(10x,5(9x,a,i4,6x))";
    static char fmt_1005[] = "(1x,\002ROW\002,i4,2x,1p5d22.13)";
    static char fmt_1002[] = "(10x,4(12x,a,i4,9x))";
    static char fmt_1006[] = "(1x,\002ROW\002,i4,2x,1p4d28.19)";
    static char fmt_1003[] = "(10x,3(16x,a,i4,13x))";
    static char fmt_1007[] = "(1x,\002ROW\002,i4,2x,1p3d36.27)";

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    cilist ci__1;

    /* Builtin functions */
    integer s_wsfe(cilist*), e_wsfe(void), do_fio(integer*, char*, ftnlen);

    /* Local variables */
    static integer i__, j, k1, k2, lout;
    extern integer i1mach_(integer*);
    static integer ndigit;

    /* Fortran I/O blocks */
    static cilist io___6 = { 0, 0, 0, fmt_1010, 0 };
    static cilist io___8 = { 0, 0, 0, fmt_1009, 0 };
    static cilist io___10 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___11 = { 0, 0, 0, fmt_1004, 0 };
    static cilist io___12 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___13 = { 0, 0, 0, fmt_1005, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_1002, 0 };
    static cilist io___15 = { 0, 0, 0, fmt_1006, 0 };
    static cilist io___16 = { 0, 0, 0, fmt_1003, 0 };
    static cilist io___17 = { 0, 0, 0, fmt_1007, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_1009, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_1000, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_1004, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_1001, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_1005, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_1002, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_1006, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_1003, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_1007, 0 };


    /* ***BEGIN PROLOGUE  DMOUT */
    /* ***REFER TO  DBOCLS,DFC */
    /* ***ROUTINES CALLED  I1MACH */
    /* ***DESCRIPTION */

    /*     DOUBLE PRECISION MATRIX OUTPUT ROUTINE. */

    /*  INPUT.. */

    /*  M,N,LDA,A(*,*) PRINT THE DOUBLE PRECISION ARRAY A(I,J),I = 1,...,M, */
    /*                 J=1,...,N, ON OUTPUT UNIT LOUT=6. LDA IS THE DECLARED */
    /*                 FIRST DIMENSION OF A(*,*) AS SPECIFIED IN THE CALLING */
    /*                 PROGRAM. THE HEADING IN THE FORTRAN FORMAT STATEMENT */
    /*                 IFMT(*), DESCRIBED BELOW, IS PRINTED AS A FIRST STEP. */
    /*                 THE COMPONENTS A(I,J) ARE INDEXED, ON OUTPUT, IN A */
    /*                 PLEASANT FORMAT. */
    /*  IFMT(*)        A FORTRAN FORMAT STATEMENT. THIS IS PRINTED ON */
    /*                 OUTPUT UNIT LOUT=6 WITH THE VARIABLE FORMAT FORTRAN */
    /*                 STATEMENT */
    /*                       WRITE(LOUT,IFMT). */
    /*  IDIGIT         PRINT AT LEAST IABS(IDIGIT) DECIMAL DIGITS PER NUMBER. */
    /*                 THE SUBPROGRAM WILL CHOOSE THAT INTEGER 4,6,14,20 OR */
    /*                 28 WHICH WILL PRINT AT LEAST IABS(IDIGIT) NUMBER OF */
    /*                 PLACES.  IF IDIGIT.LT.0, 72 PRINTING COLUMNS ARE */
    /*                 UTILIZED TO WRITE EACH LINE OF OUTPUT OF THE ARRAY */
    /*                 A(*,*). (THIS CAN BE USED ON MOST TIME-SHARING */
    /*                 TERMINALS).  IF IDIGIT.GE.0, 133 PRINTING COLUMNS ARE */
    /*                 UTILIZED. (THIS CAN BE USED ON MOST LINE PRINTERS). */

    /*  EXAMPLE.. */

    /*  PRINT AN ARRAY CALLED (SIMPLEX TABLEAU   ) OF SIZE 10 BY 20 SHOWING */
    /*  6 DECIMAL DIGITS PER NUMBER. THE USER IS RUNNING ON A TIME-SHARING */
    /*  SYSTEM WITH A 72 COLUMN OUTPUT DEVICE. */

    /*     DOUBLE PRECISION TABLEU(20,20) */
    /*     M = 10 */
    /*     N = 20 */
    /*     LDTABL = 20 */
    /*     IDIGIT = -6 */
    /*     CALL DMOUT(M,N,LDTABL,TABLEU,21H(16H1SIMPLEX TABLEAU),IDIGIT) */



    /*     AUTHORS    JOHN A. WISNIEWSKI   SANDIA LABS ALBUQUERQUE. */
    /*                RICHARD J. HANSON    SANDIA LABS ALBUQUERQUE. */
    /*     DATE       JULY 30,1978. */
    /* ***END PROLOGUE  DMOUT */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    /* ***FIRST EXECUTABLE STATEMENT  DMOUT */
    lout = i1mach_(&c__2);
    ci__1.cierr = 0;
    ci__1.ciunit = lout;
    ci__1.cifmt = ifmt;
    s_wsfe(&ci__1);
    e_wsfe();
    if (*m <= 0 || *n <= 0 || *lda <= 0) {
        return 0;
    }
    ndigit = *idigit;
    if (*idigit == 0) {
        ndigit = 4;
    }
    if (*idigit >= 0) {
        goto L80;
    }

    ndigit = -(*idigit);
    if (ndigit > 4) {
        goto L9;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 5) {
        /* Computing MIN */
        i__2 = *n, i__3 = k1 + 4;
        k2 = min(i__2,i__3);
        io___6.ciunit = lout;
        s_wsfe(&io___6);
        i__2 = k2;
        for (i__ = k1; i__ <= i__2; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            io___8.ciunit = lout;
            s_wsfe(&io___8);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L5: */
        }
    }
    return 0;

L9:
    if (ndigit > 6) {
        goto L20;
    }

    i__2 = *n;
    for (k1 = 1; k1 <= i__2; k1 += 4) {
        /* Computing MIN */
        i__1 = *n, i__3 = k1 + 3;
        k2 = min(i__1,i__3);
        io___10.ciunit = lout;
        s_wsfe(&io___10);
        i__1 = k2;
        for (i__ = k1; i__ <= i__1; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            io___11.ciunit = lout;
            s_wsfe(&io___11);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L10: */
        }
    }
    return 0;

L20:
    if (ndigit > 14) {
        goto L40;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 2) {
        /* Computing MIN */
        i__2 = *n, i__3 = k1 + 1;
        k2 = min(i__2,i__3);
        io___12.ciunit = lout;
        s_wsfe(&io___12);
        i__2 = k2;
        for (i__ = k1; i__ <= i__2; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            io___13.ciunit = lout;
            s_wsfe(&io___13);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L30: */
        }
    }
    return 0;

L40:
    if (ndigit > 20) {
        goto L60;
    }

    i__2 = *n;
    for (k1 = 1; k1 <= i__2; k1 += 2) {
        /* Computing MIN */
        i__1 = *n, i__3 = k1 + 1;
        k2 = min(i__1,i__3);
        io___14.ciunit = lout;
        s_wsfe(&io___14);
        i__1 = k2;
        for (i__ = k1; i__ <= i__1; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            io___15.ciunit = lout;
            s_wsfe(&io___15);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L50: */
        }
    }
    return 0;

L60:
    i__1 = *n;
    for (k1 = 1; k1 <= i__1; ++k1) {
        k2 = k1;
        io___16.ciunit = lout;
        s_wsfe(&io___16);
        i__2 = k2;
        for (i__ = k1; i__ <= i__2; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            io___17.ciunit = lout;
            s_wsfe(&io___17);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L70: */
        }
    }
    return 0;

L80:
    if (ndigit > 4) {
        goto L86;
    }

    i__2 = *n;
    for (k1 = 1; k1 <= i__2; k1 += 10) {
        /* Computing MIN */
        i__1 = *n, i__3 = k1 + 9;
        k2 = min(i__1,i__3);
        io___18.ciunit = lout;
        s_wsfe(&io___18);
        i__1 = k2;
        for (i__ = k1; i__ <= i__1; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            io___19.ciunit = lout;
            s_wsfe(&io___19);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L85: */
        }
    }

L86:
    if (ndigit > 6) {
        goto L100;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 8) {
        /* Computing MIN */
        i__2 = *n, i__3 = k1 + 7;
        k2 = min(i__2,i__3);
        io___20.ciunit = lout;
        s_wsfe(&io___20);
        i__2 = k2;
        for (i__ = k1; i__ <= i__2; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            io___21.ciunit = lout;
            s_wsfe(&io___21);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L90: */
        }
    }
    return 0;

L100:
    if (ndigit > 14) {
        goto L120;
    }

    i__2 = *n;
    for (k1 = 1; k1 <= i__2; k1 += 5) {
        /* Computing MIN */
        i__1 = *n, i__3 = k1 + 4;
        k2 = min(i__1,i__3);
        io___22.ciunit = lout;
        s_wsfe(&io___22);
        i__1 = k2;
        for (i__ = k1; i__ <= i__1; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            io___23.ciunit = lout;
            s_wsfe(&io___23);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L110: */
        }
    }
    return 0;

L120:
    if (ndigit > 20) {
        goto L140;
    }

    i__1 = *n;
    for (k1 = 1; k1 <= i__1; k1 += 4) {
        /* Computing MIN */
        i__2 = *n, i__3 = k1 + 3;
        k2 = min(i__2,i__3);
        io___24.ciunit = lout;
        s_wsfe(&io___24);
        i__2 = k2;
        for (i__ = k1; i__ <= i__2; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__2 = *m;
        for (i__ = 1; i__ <= i__2; ++i__) {
            io___25.ciunit = lout;
            s_wsfe(&io___25);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L130: */
        }
    }
    return 0;

L140:
    i__2 = *n;
    for (k1 = 1; k1 <= i__2; k1 += 3) {
        /* Computing MIN */
        i__1 = *n, i__3 = k1 + 2;
        k2 = min(i__1,i__3);
        io___26.ciunit = lout;
        s_wsfe(&io___26);
        i__1 = k2;
        for (i__ = k1; i__ <= i__1; ++i__) {
            do_fio(&c__1, icol, (ftnlen)3);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
        }
        e_wsfe();
        i__1 = *m;
        for (i__ = 1; i__ <= i__1; ++i__) {
            io___27.ciunit = lout;
            s_wsfe(&io___27);
            do_fio(&c__1, (char*)&i__, (ftnlen)sizeof(integer));
            i__3 = k2;
            for (j = k1; j <= i__3; ++j) {
                do_fio(&c__1, (char*)&a[i__ + j * a_dim1], (ftnlen)sizeof(
                           doublereal));
            }
            e_wsfe();
            /* L150: */
        }
    }
    return 0;
} /* dmout_ */

