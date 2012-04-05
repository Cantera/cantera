/* dbolsm.f -- translated by f2c (version 20031025).
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

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__2 = 2;
static long int lc__2 = 2;
static integer c_n4 = -4;
static integer c__4 = 4;
static long int lc__4 = 4;

/* DECK DBOLSM */

/* 4/16/91 Harry K. Moffat Div. 1126 844-6912 */
/*         Fixed an error! */
/*         DC and DS were not declared double precision */

/* Subroutine */
int dbolsm_(doublereal* w, integer* mdw, integer* minput,
            integer* ncols, doublereal* bl, doublereal* bu, integer* ind, integer
            *iopt, doublereal* x, doublereal* rnorm, integer* mode, doublereal *
            rw, doublereal* ww, doublereal* scl, integer* ibasis, integer* ibb)
{
    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    static integer i__, j;
    static doublereal t, t1, t2, dc, sc, ds;
    static integer ip, jp, lp;
    static doublereal ss, wt, cl1, cl2, cl3, fac, big, one;
    static integer lds;
    static doublereal bou, two, beta;
    static integer jbig, jmag, ioff, jcol;
    static doublereal wbig;
    extern doublereal ddot_(integer*, doublereal*, integer*, doublereal*,
                            integer*);
    static doublereal wmag;
    static integer idum, mval, iter, nerr;
    static real rdum;
    extern /* Subroutine */ int drot_(integer*, doublereal*, integer*,
                                      doublereal*, integer*, doublereal*, doublereal*);
    static doublereal zero, xnew;
    extern doublereal dnrm2_(integer*, doublereal*, integer*);
    static real rdum2;
    static integer icase;
    static doublereal alpha;
    static integer nchar, level;
    static logical found;
    static integer nsetb;
    extern /* Subroutine */ int dcopy_(integer*, doublereal*, integer*,
                                       doublereal*, integer*);
    static integer igopr, itemp, itmax;
    extern /* Subroutine */ int drotg_(doublereal*, doublereal*, doublereal
                                       *, doublereal*), dswap_(integer*, doublereal*, integer*,
                                               doublereal*, integer*), daxpy_(integer*, doublereal*,
                                                       doublereal*, integer*, doublereal*, integer*);
    static integer lgopr;
    extern /* Subroutine */ int dmout_(integer*, integer*, integer*,
                                       doublereal*, char*, integer*, ftnlen);
    static integer jdrop;
    extern doublereal d1mach_(long int*);
    extern /* Subroutine */ int dvout_(integer*, doublereal*, char*,
                                       integer*, ftnlen), ivout_(integer*, integer*, char*, integer *
                                               , ftnlen);
    static integer mrows, jdrop1, jdrop2, jlarge;
    static doublereal colabv, colblo, wlarge;
    static integer nlevel;
    static doublereal tolind;
    static integer iprint;
    static logical constr;
    static doublereal tolsze;
    extern /* Subroutine */ int xerrwv_(char*, integer*, integer*, integer
                                        *, integer*, integer*, integer*, integer*, real*, real*,
                                        ftnlen);

    /* ***BEGIN PROLOGUE  DBOLSM */
    /* ***REFER TO  SBOCLS,SBOLS */
    /* ***ROUTINES CALLED  D1MACH,DAXPY,DCOPY,DDOT,DMOUT,DNRM2,DROT,DROTG, */
    /*                    DSWAP,DVOUT,IVOUT,XERRWV */
    /* ***DESCRIPTION */

    /*            **** Double Precision Version of SBOLSM **** */
    /*   **** All INPUT and OUTPUT real variables are DOUBLE PRECISION **** */

    /*          Solve E*X = F (least squares sense) with bounds on */
    /*            selected X values. */
    /*     The user must have dimension statements of the form: */

    /*       DIMENSION W(MDW,NCOLS+1), BL(NCOLS), BU(NCOLS), */
    /*      * X(NCOLS+NX), RW(NCOLS), WW(NCOLS), SCL(NCOLS) */
    /*       INTEGER IND(NCOLS), IOPT(1+NI), IBASIS(NCOLS), IBB(NCOLS) */

    /*     (here NX=number of extra locations required for options 1,...,7; */
    /*     NX=0 for no options; here NI=number of extra locations possibly */
    /*     required for options 1-7; NI=0 for no options; NI=14 if all the */
    /*     options are simultaneously in use.) */

    /*    INPUT */
    /*    ----- */

    /*    -------------------- */
    /*    W(MDW,*),MROWS,NCOLS */
    /*    -------------------- */
    /*     The array w(*,*) contains the matrix [E:F] on entry. The matrix */
    /*     [E:F] has MROWS rows and NCOLS+1 columns. This data is placed in */
    /*     the array W(*,*) with E occupying the first NCOLS columns and the */
    /*     right side vector F in column NCOLS+1. The row dimension, MDW, of */
    /*     the array W(*,*) must satisfy the inequality MDW .ge. MROWS. */
    /*     Other values of MDW are errors. The values of MROWS and NCOLS */
    /*     must be positive. Other values are errors. */

    /*    ------------------ */
    /*    BL(*),BU(*),IND(*) */
    /*    ------------------ */
    /*     These arrays contain the information about the bounds that the */
    /*     solution values are to satisfy. The value of IND(J) tells the */
    /*     type of bound and BL(J) and BU(J) give the explicit values for */
    /*     the respective upper and lower bounds. */

    /*    1.    For IND(J)=1, require X(J) .ge. BL(J). */
    /*    2.    For IND(J)=2, require X(J) .le. BU(J). */
    /*    3.    For IND(J)=3, require X(J) .ge. BL(J) and */
    /*                                X(J) .le. BU(J). */
    /*    4.    For IND(J)=4, no bounds on X(J) are required. */
    /*     The values of BL(*),BL(*) are modified by the subprogram. Values */
    /*     other than 1,2,3 or 4 for IND(J) are errors. In the case IND(J)=3 */
    /*     (upper and lower bounds) the condition BL(J) .gt. BU(J) is an */
    /*     error. */

    /*    ------- */
    /*    IOPT(*) */
    /*    ------- */
    /*     This is the array where the user can specify nonstandard options */
    /*     for DBOLSM( ). Most of the time this feature can be ignored by */
    /*     setting the input value IOPT(1)=99. Occasionally users may have */
    /*     needs that require use of the following subprogram options. For */
    /*     details about how to use the options see below: IOPT(*) CONTENTS. */

    /*     Option Number   Brief Statement of Purpose */
    /*     ----- ------   ----- --------- -- ------- */
    /*           1         Move the IOPT(*) processing pointer. */
    /*           2         Change rank determination tolerance. */
    /*           3         Change blow-up factor that determines the */
    /*                     size of variables being dropped from active */
    /*                     status. */
    /*           4         Reset the maximum number of iterations to use */
    /*                     in solving the problem. */
    /*           5         The data matrix is triangularized before the */
    /*                     problem is solved whenever (NCOLS/MROWS) .lt. */
    /*                     FAC. Change the value of FAC. */
    /*           6         Redefine the weighting matrix used for */
    /*                     linear independence checking. */
    /*           7         Debug output is desired. */
    /*          99         No more options to change. */

    /*    ---- */
    /*    X(*) */
    /*    ---- */
    /*     This array is used to pass data associated with options 1,2,3 and */
    /*     5. Ignore this input parameter if none of these options are used. */
    /*     Otherwise see below: IOPT(*) CONTENTS. */

    /*    ---------------- */
    /*    IBASIS(*),IBB(*) */
    /*    ---------------- */
    /*     These arrays must be initialized by the user. The values */
    /*         IBASIS(J)=J, J=1,...,NCOLS */
    /*         IBB(J)   =1, J=1,...,NCOLS */
    /*     are appropriate except when using nonstandard features. */

    /*    ------ */
    /*    SCL(*) */
    /*    ------ */
    /*     This is the array of scaling factors to use on the columns of the */
    /*     matrix E. These values must be defined by the user. To suppress */
    /*     any column scaling set SCL(J)=1.0, J=1,...,NCOLS. */

    /*    OUTPUT */
    /*    ------ */

    /*    ---------- */
    /*    X(*),RNORM */
    /*    ---------- */
    /*     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for */
    /*     the constrained least squares problem. The value RNORM is the */
    /*     minimum residual vector length. */

    /*    ---- */
    /*    MODE */
    /*    ---- */
    /*     The sign of mode determines whether the subprogram has completed */
    /*     normally, or encountered an error condition or abnormal status. */
    /*     A value of MODE .ge. 0 signifies that the subprogram has completed */
    /*     normally. The value of MODE (.ge. 0) is the number of variables */
    /*     in an active status: not at a bound nor at the value ZERO, for */
    /*     the case of free variables. A negative value of MODE will be one */
    /*     of the 20 cases -40,-39,...,-22, or -1. Values .lt. -1 correspond */
    /*     to an abnormal completion of the subprogram. To understand the */
    /*     abnormal completion codes see below: ERROR MESSAGES for DBOLSM( ) */
    /*     An approximate solution will be returned to the user only when */
    /*     max. iterations is reached, MODE=-22. */

    /*    ----------- */
    /*    RW(*),WW(*) */
    /*    ----------- */
    /*     These are working arrays each with NCOLS entries. The array RW(*) */
    /*     contains the working (scaled, nonactive) solution values. The */
    /*     array WW(*) contains the working (scaled, active) gradient vector */
    /*     values. */

    /*    ---------------- */
    /*    IBASIS(*),IBB(*) */
    /*    ---------------- */
    /*     These arrays contain information about the status of the solution */
    /*     when MODE .ge. 0. The indices IBASIS(K), K=1,...,MODE, show the */
    /*     nonactive variables; indices IBASIS(K), K=MODE+1,..., NCOLS are */
    /*     the active variables. The value (IBB(J)-1) is the number of times */
    /*     variable J was reflected from its upper bound. (normally the user */
    /*     can ignore these parameters.) */

    /*    IOPT(*) CONTENTS */
    /*    ------- -------- */
    /*     The option array allows a user to modify internal variables in */
    /*     the subprogram without recompiling the source code. A central */
    /*     goal of the initial software design was to do a good job for most */
    /*     people. Thus the use of options will be restricted to a select */
    /*     group of users. The processing of the option array proceeds as */
    /*     follows: a pointer, here called LP, is initially set to the value */
    /*     1. The value is updated as the options are processed.  At the */
    /*     pointer position the option number is extracted and used for */
    /*     locating other information that allows for options to be changed. */
    /*     The portion of the array IOPT(*) that is used for each option is */
    /*     fixed; the user and the subprogram both know how many locations */
    /*     are needed for each option. A great deal of error checking is */
    /*     done by the subprogram on the contents of the option array. */
    /*     Nevertheless it is still possible to give the subprogram optional */
    /*     input that is meaningless. For example some of the options use */
    /*     the location X(NCOLS+IOFF) for passing data. The user must manage */
    /*     the allocation of these locations when more than one piece of */
    /*     option data is being passed to the subprogram. */

    /*   1 */
    /*   - */
    /*     Move the processing pointer (either forward or backward) to the */
    /*     location IOPT(LP+1). The processing pointer is moved to location */
    /*     LP+2 of IOPT(*) in case IOPT(LP)=-1.  For example to skip over */
    /*     locations 3,...,NCOLS+2 of IOPT(*), */

    /*       IOPT(1)=1 */
    /*       IOPT(2)=NCOLS+3 */
    /*       (IOPT(I), I=3,...,NCOLS+2 are not defined here.) */
    /*       IOPT(NCOLS+3)=99 */
    /*       CALL DBOLSM( ) */

    /*     CAUTION: Misuse of this option can yield some very hard */
    /*     -to-find bugs.  Use it with care. */

    /*   2 */
    /*   - */
    /*     The algorithm that solves the bounded least squares problem */
    /*     iteratively drops columns from the active set. This has the */
    /*     effect of joining a new column vector to the QR factorization of */
    /*     the rectangular matrix consisting of the partially triangularized */
    /*     nonactive columns. After triangularizing this matrix a test is */
    /*     made on the size of the pivot element. The column vector is */
    /*     rejected as dependent if the magnitude of the pivot element is */
    /*     .le. TOL* magnitude of the column in components strictly above */
    /*     the pivot element. Nominally the value of this (rank) tolerance */
    /*     is TOL = SQRT(SRELPR), where SRELPR is relative machine */
    /*     precision. To change only the value of TOL, for example, */

    /*       X(NCOLS+1)=TOL */
    /*       IOPT(1)=2 */
    /*       IOPT(2)=1 */
    /*       IOPT(3)=99 */
    /*       CALL DBOLSM() */

    /*     Generally, if LP is the processing pointer for IOPT(*), */

    /*       X(NCOLS+IOFF)=TOL */
    /*       IOPT(LP)=2 */
    /*       IOPT(LP+1)=IOFF */
    /*        . */
    /*       CALL DBOLSM() */

    /*     The required length of IOPT(*) is increased by 2 if option 2 is */
    /*     used; The required length of X(*) is increased by 1. A value of */
    /*     IOFF .le. 0 is an error. A value of TOL .le. SRELPR gives a */
    /*     warning message; it is not considered an error. */
    /*     Here SRELPR is the relative machine precision. */

    /*   3 */
    /*   - */
    /*     A solution component is left active (not used) if, roughly */
    /*     speaking, it seems too large. Mathematically the new component is */
    /*     left active if the magnitude is .ge.((vector norm of F)/(matrix */
    /*     norm of E))/BLOWUP. Nominally the factor BLOWUP = SQRT(SRELPR) */
    /*     where SRELPR is the relative machine precision. To change only */
    /*     the value of BLOWUP, for example, */

    /*       X(NCOLS+2)=BLOWUP */
    /*       IOPT(1)=3 */
    /*       IOPT(2)=2 */
    /*       IOPT(3)=99 */
    /*       CALL DBOLSM() */

    /*     Generally, if LP is the processing pointer for IOPT(*), */

    /*       X(NCOLS+IOFF)=BLOWUP */
    /*       IOPT(LP)=3 */
    /*       IOPT(LP+1)=IOFF */
    /*        . */
    /*       CALL DBOLSM() */

    /*     The required length of IOPT(*) is increased by 2 if option 3 is */
    /*     used; the required length of X(*) is increased by 1. A value of */
    /*     IOFF .le. 0 is an error. A value of BLOWUP .le. 0.0 is an error. */

    /*   4 */
    /*   - */
    /*     Normally the algorithm for solving the bounded least squares */
    /*     problem requires between NCOLS/3 and NCOLS drop-add steps to */
    /*     converge. (this remark is based on examining a small number of */
    /*     test cases.) The amount of arithmetic for such problems is */
    /*     typically about twice that required for linear least squares if */
    /*     there are no bounds and if plane rotations are used in the */
    /*     solution method. Convergence of the algorithm, while */
    /*     mathematically certain, can be much slower than indicated. To */
    /*     avoid this potential but unlikely event ITMAX drop-add steps are */
    /*     permitted. Nominally ITMAX=5*(MAX(MROWS,NCOLS)). To change the */
    /*     value of ITMAX, for example, */

    /*       IOPT(1)=4 */
    /*       IOPT(2)=ITMAX */
    /*       IOPT(3)=99 */
    /*       CALL DBOLSM() */

    /*     Generally, if LP is the processing pointer for IOPT(*), */

    /*       IOPT(LP)=4 */
    /*       IOPT(LP+1)=ITMAX */
    /*        . */
    /*       CALL DBOLSM() */

    /*     The value of ITMAX must be .gt. 0. Other values are errors. Use */
    /*     of this option increases the required length of IOPT(*) by 2. */

    /*   5 */
    /*   - */
    /*     For purposes of increased efficiency the MROWS by NCOLS+1 data */
    /*     matrix [E:F] is triangularized as a first step whenever MROWS */
    /*     satisfies FAC*MROWS .gt. NCOLS. Nominally FAC=0.75. To change the */
    /*     value of FAC, */

    /*       X(NCOLS+3)=FAC */
    /*       IOPT(1)=5 */
    /*       IOPT(2)=3 */
    /*       IOPT(3)=99 */
    /*       CALL DBOLSM() */

    /*     Generally, if LP is the processing pointer for IOPT(*), */

    /*       X(NCOLS+IOFF)=FAC */
    /*       IOPT(LP)=5 */
    /*       IOPT(LP+1)=IOFF */
    /*        . */
    /*       CALL DBOLSM() */

    /*     The value of FAC must be nonnegative. Other values are errors. */
    /*     Resetting FAC=0.0 suppresses the initial triangularization step. */
    /*     Use of this option increases the required length of IOPT(*) by 2; */
    /*     The required length of of X(*) is increased by 1. */

    /*   6 */
    /*   - */
    /*     The norm used in testing the magnitudes of the pivot element */
    /*     compared to the mass of the column above the pivot line can be */
    /*     changed. The type of change that this option allows is to weight */
    /*     the components with an index larger than MVAL by the parameter */
    /*     WT. Normally MVAL=0 and WT=1. To change both the values MVAL and */
    /*     WT, where LP is the processing pointer for IOPT(*), */

    /*       X(NCOLS+IOFF)=WT */
    /*       IOPT(LP)=6 */
    /*       IOPT(LP+1)=IOFF */
    /*       IOPT(LP+2)=MVAL */

    /*     Use of this option increases the required length of IOPT(*) by 3. */
    /*     The length of X(*) is increased by 1. Values of MVAL must be */
    /*     nonnegative and not greater than MROWS. Other values are errors. */
    /*     The value of WT must be positive. Any other value is an error. If */
    /*     either error condition is present a message will be printed. */

    /*   7 */
    /*   - */
    /*     Debug output, showing the detailed add-drop steps for the */
    /*     constrained least squares problem, is desired. This option is */
    /*     intended to be used to locate suspected bugs. */

    /*   99 */
    /*   -- */
    /*     There are no more options to change. */

    /*     The values for options are 1,...,7,99, and are the only ones */
    /*     permitted. Other values are errors. Options -99,-1,...,-7 mean */
    /*     that the repective options 99,1,...,7 are left at their default */
    /*     values. An example is the option to modify the (rank) tolerance: */

    /*       X(NCOLS+1)=TOL */
    /*       IOPT(1)=-2 */
    /*       IOPT(2)=1 */
    /*       IOPT(3)=99 */

    /*    Error Messages for DBOLSM( ) */
    /*    ----- -------- --- --------- */
    /*    WARNING IN... */
    /*    DBOLSM(). MORE THAN (I1)=ITMAX ITERATIONS SOLVING BOUNDED LEAST */
    /*    SQUARES PROBLEM. */
    /*              IN ABOVE MESSAGE, I1=         3 */
    /*    ERROR NUMBER =        22 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM. THE OPTION NUMBER=(I1) IS NOT DEFINED. */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*    ERROR NUMBER =        23 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE OFFSET=(I1) BEYOND POSTION NCOLS=(I2) */
    /*    MUST BE POSITIVE FOR OPTION NUMBER 2. */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*              IN ABOVE MESSAGE, I2=         1 */
    /*    ERROR NUMBER =        24 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE TOLERANCE FOR RANK DETERMINATION=(R1) */
    /*    IS LESS THAN MACHINE PRECISION=(R2). */
    /*              IN ABOVE MESSAGE, R1=    0. */
    /*              IN ABOVE MESSAGE, R2=     .7105427358E-14 */
    /*    ERROR NUMBER =        25 */

    /*    WARNING IN... */
    /*    DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2) MUST */
    /*     BE POSTIVE FOR OPTION NUMBER 3. */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*              IN ABOVE MESSAGE, I2=         1 */
    /*    ERROR NUMBER =        26 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE RECIPROCAL OF THE BLOW-UP FACTOR FOR REJECTING */
    /*    VARIABLES MUST BE POSITIVE. NOW=(R1). */
    /*              IN ABOVE MESSAGE, R1=    0. */
    /*    ERROR NUMBER =        27 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE MAXIMUM NUMBER OF ITERATIONS=(I1) MUST BE POSITIVE. */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*    ERROR NUMBER =        28 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2) MUST BE */
    /*    POSTIVE FOR OPTION NUMBER 5. */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*              IN ABOVE MESSAGE, I2=         1 */
    /*    ERROR NUMBER =        29 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE FACTOR (NCOLS/MROWS) WHERE PRETRIANGULARIZING IS */
    /*    PERFORMED MUST BE NONNEGATIVE. NOW=(R1). */
    /*              IN ABOVE MESSAGE, R1=    -.2500000000E+00 */
    /*    ERROR NUMBER =        30 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE NUMBER OF ROWS=(I1) MUST BE POSITIVE. */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*    ERROR NUMBER =        31 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE NUMBER OF COLS.=(I1) MUST BE POSTIVE. */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*    ERROR NUMBER =        32 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE. THE */
    /*    NUMBER OF ROWS =(I2). */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*              IN ABOVE MESSAGE, I2=         1 */
    /*    ERROR NUMBER =        33 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). FOR J=(I1) THE CONSTRAINT INDICATOR MUST BE 1-4. */
    /*              IN ABOVE MESSAGE, I1=         1 */
    /*              IN ABOVE MESSAGE, I2=         0 */
    /*    ERROR NUMBER =        34 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). FOR J=(I1) THE LOWER BOUND=(R1) IS .GT. THE UPPER */
    /*     BOUND=(R2). */
    /*              IN ABOVE MESSAGE, I1=         1 */
    /*              IN ABOVE MESSAGE, R1=    0. */
    /*              IN ABOVE MESSAGE, R2=    -.1000000000E+01 */
    /*    ERROR NUMBER =        35 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE INPUT ORDER OF COLUMNS=(I1) IS NOT BETWEEN 1 */
    /*    AND NCOLS=(I2). */
    /*              IN ABOVE MESSAGE, I1=         0 */
    /*              IN ABOVE MESSAGE, I2=         1 */
    /*    ERROR NUMBER =        36 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE BOUND POLARITY FLAG IN COMPONENT J=(I1) MUST */
    /*    BE POSITIVE.  NOW=(I2). */
    /*              IN ABOVE MESSAGE, I1=         1 */
    /*              IN ABOVE MESSAGE, I2=         0 */
    /*    ERROR NUMBER =        37 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE ROW SEPARATOR TO APPLY WEIGHTING (I1) MUST LIE */
    /*    BETWEEN 0 AND MROWS (I2). WEIGHT (R1) MUST BE POSITIVE. */
    /*              IN ABOVE MESSAGE, I1=        -1 */
    /*              IN ABOVE MESSAGE, I2=         2 */
    /*              IN ABOVE MESSAGE, R1=    0. */
    /*    ERROR NUMBER =        38 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE OFFSET (I1) BEYOND POSITION NCOLS=(I2) MUST BE */
    /*    POSITIVE FOR OPTION NUMBER 7. */
    /*              IN ABOVE MESSAGE, I1=        -1 */
    /*              IN ABOVE MESSAGE, I2=         2 */
    /*    ERROR NUMBER =        39 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */

    /*    WARNING IN... */
    /*    DBOLSM(). THE COLUMN PIVOTING THRESHOLD FACTOR MUST BE */
    /*    POSITIVE. NOW=(R1). */
    /*              IN ABOVE MESSAGE, R1=    0. */
    /*    ERROR NUMBER =        40 */
    /*    (NORMALLY A RETURN TO THE USER TAKES PLACE FOLLOWING THIS MESSAGE.) */
    /* ***END PROLOGUE  DBOLSM */

    /*     PURPOSE */
    /*     ------- */
    /*     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE BOUNDED */
    /*     LEAST SQUARES PROBLEM.  THE PROBLEM SOLVED HERE IS: */

    /*     SOLVE E*X =  F  (LEAST SQUARES SENSE) */
    /*     WITH BOUNDS ON SELECTED X VALUES. */

    /*     REVISED 850422-0800 */
    /*     REVISED YYMMDD-HHMM */

    /*     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN */
    /*     EDITING AT THE CARD 'C++'. */
    /*     CHANGE THE SUBPROGRAM NAME TO DBOLSM AND THE STRINGS */
    /*     /SAXPY/ TO /DAXPY/, /SCOPY/ TO /DCOPY/, */
    /*     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/, */
    /*     /SROTG/ TO /DROTG/, /SROT/ TO /DROT/, /R1MACH/ TO /D1MACH/, */
    /*     /SVOUT/ TO /DVOUT/, /SMOUT/ TO /DMOUT/, */
    /*     /SSWAP/ TO /DSWAP/, /E0/ TO /D0/, */
    /*     /REAL            / TO /DOUBLE PRECISION/. */
    /* ++ */

    /*     THIS VARIABLE REMAINS TYPE REAL. */
    /* ***FIRST EXECUTABLE STATEMENT  DBOLSM */
    /* Parameter adjustments */
    w_dim1 = *mdw;
    w_offset = 1 + w_dim1;
    w -= w_offset;
    --bl;
    --bu;
    --ind;
    --iopt;
    --x;
    --rw;
    --ww;
    --scl;
    --ibasis;
    --ibb;

    /* Function Body */
    level = 1;
    /*    -> HKM 5/21/04 Put this line below in to avoid a valgrind incident */
    mrows = 0;

    /*    VERIFY THAT THE PROBLEM DIMENSIONS ARE DEFINED PROPERLY. */
    if (*minput <= 0) {
        nerr = 31;
        nchar = 51;
        xerrwv_("DBOLSM(). THE NUMBER OF ROWS=(I1) MUST BE POSITIVE.", &nchar,
                &nerr, &level, &c__1, minput, &idum, &c__0, &rdum, &rdum, (
                    ftnlen)51);
        /*     DO(RETURN TO USER PROGRAM UNIT) */
        goto L610;
    }
    if (*ncols <= 0) {
        nerr = 32;
        nchar = 51;
        xerrwv_("DBOLSM(). THE NUMBER OF COLS.=(I1) MUST BE POSTIVE.", &nchar,
                &nerr, &level, &c__1, ncols, &idum, &c__0, &rdum, &rdum, (
                    ftnlen)51);
        /*     DO(RETURN TO USER PROGRAM UNIT) */
        goto L610;
    }
    if (*mdw < *minput) {
        nerr = 33;
        nchar = 78;
        xerrwv_("DBOLSM(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE. THE N"
                "UMBER OF ROWS=(I2).", &nchar, &nerr, &level, &c__2, mdw, &
                mrows, &c__0, &rdum, &rdum, (ftnlen)78);
        /*     DO(RETURN TO USER PROGRAM UNIT) */
        goto L610;
    }

    /*     VERIFY THAT BOUND INFORMATION IS CORRECT. */
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        if (ind[j] < 1 || ind[j] > 4) {
            nerr = 34;
            nchar = 58;
            xerrwv_("DBOLSM(). FOR J=(I1) THE CONSTRAINT INDICATOR MUST BE 1"
                    "-4.", &nchar, &nerr, &level, &c__2, &j, &ind[j], &c__0, &
                    rdum, &rdum, (ftnlen)58);
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L610;
        }
        /* L10: */
    }
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        if (ind[j] == 3) {
            if (bu[j] < bl[j]) {
                nerr = 35;
                nchar = 71;
                rdum2 = (real) bl[j];
                rdum = (real) bu[j];
                xerrwv_("DBOLSM(). FOR J=(I1) THE LOWER BOUND=(R1) IS .GT. T"
                        "HE UPPER BOUND=(R2).", &nchar, &nerr, &level, &c__1, &
                        j, &idum, &c__2, &rdum2, &rdum, (ftnlen)71);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L610;
            }
        }
        /* L20: */
    }

    /*     CHECK THAT PERMUTATION AND POLARITY ARRAYS HAVE BEEN SET. */
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        if (ibasis[j] < 1 || ibasis[j] > *ncols) {
            nerr = 36;
            nchar = 74;
            xerrwv_("DBOLSM(). THE INPUT ORDER OF COLUMNS=(I1) IS NOT BETWEE"
                    "N 1 AND NCOLS=(I2).", &nchar, &nerr, &level, &c__2, &
                    ibasis[j], ncols, &c__0, &rdum, &rdum, (ftnlen)74);
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L610;
        }
        if (ibb[j] <= 0) {
            nerr = 37;
            nchar = 81;
            xerrwv_("DBOLSM(). THE BOUND POLARITY FLAG IN COMPONENT J=(I1) M"
                    "UST BE POSITIVE. NOW=(I2).", &nchar, &nerr, &level, &c__2,
                    &j, &ibb[j], &c__0, &rdum, &rdum, (ftnlen)81);
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L610;
        }
        /* L30: */
    }
    /*     DO(PROCESS OPTION ARRAY) */
    goto L580;
L40:
    /*     DO(INITIALIZE VARIABLES AND DATA VALUES) */
    goto L470;
L50:
    if (iprint > 0) {
        i__1 = *ncols + 1;
        dmout_(&mrows, &i__1, mdw, &w[w_offset], "(' PRETRI. INPUT MATRIX')",
               &c_n4, (ftnlen)25);
        dvout_(ncols, &bl[1], "(' LOWER BOUNDS')", &c_n4, (ftnlen)17);
        dvout_(ncols, &bu[1], "(' UPPER BOUNDS')", &c_n4, (ftnlen)17);
    }
L60:
    ++iter;
    if (iter <= itmax) {
        goto L80;
    }
    nerr = 22;
    nchar = 80;
    xerrwv_("DBOLSM(). MORE THAN (I1)=ITMAX ITERATIONS SOLVING BOUNDED LEAST"
            " SQUARES PROBLEM.", &nchar, &nerr, &level, &c__1, &itmax, &idum, &
            c__0, &rdum, &rdum, (ftnlen)80);
    /*     DO(RESCALE AND TRANSLATE VARIABLES) */
    igopr = 1;
    goto L130;
L70:
    /*     DO(RETURN TO USER PROGRAM UNIT) */
    goto L610;
L80:
    /*     DO(FIND A VARIABLE TO BECOME NON-ACTIVE) */
    goto L190;
L90:
    if (found) {
        goto L110;
    }
    /*     DO(RESCALE AND TRANSLATE VARIABLES) */
    igopr = 2;
    goto L130;
L100:
    *mode = nsetb;
    return 0;
L110:
    /*     DO(MAKE MOVE AND UPDATE FACTORIZATION) */
    goto L290;
L120:
    goto L60;
    /*     PROCEDURE(RESCALE AND TRANSLATE VARIABLES) */
L130:
    dcopy_(&nsetb, &x[1], &c__1, &rw[1], &c__1);
    x[1] = zero;
    dcopy_(ncols, &x[1], &c__0, &x[1], &c__1);
    i__1 = nsetb;
    for (j = 1; j <= i__1; ++j) {
        jcol = (i__2 = ibasis[j], abs(i__2));
        x[jcol] = rw[j] * (d__1 = scl[jcol], abs(d__1));
        /* L140: */
    }
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        if (ibb[j] % 2 == 0) {
            x[j] = bu[j] - x[j];
        }
        /* L150: */
    }
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        jcol = ibasis[j];
        if (jcol < 0) {
            x[-jcol] = bl[-jcol] + x[-jcol];
        }
        /* L160: */
    }
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        if (scl[j] < zero) {
            x[j] = -x[j];
        }
        /* L170: */
    }
    i__ = max(nsetb,mval);
    i__1 = mrows - i__;
    /* Computing MIN */
    i__2 = i__ + 1;
    *rnorm = dnrm2_(&i__1, &w[min(i__2,mrows) + (*ncols + 1) * w_dim1], &c__1)
             ;
    /*     END PROCEDURE */
    switch (igopr) {
    case 1:
        goto L70;
    case 2:
        goto L100;
    }
    /*     PROCEDURE(FIND A VARIABLE TO BECOME NON-ACTIVE) */
L190:

    /*     COMPUTE (NEGATIVE) OF GRADIENT VECTOR, W= */
    /*     (TRANSPOSE OF E)*(F-E*X). */
    ww[1] = zero;
    dcopy_(ncols, &ww[1], &c__0, &ww[1], &c__1);
    i__1 = *ncols;
    for (j = nsetb + 1; j <= i__1; ++j) {
        jcol = (i__2 = ibasis[j], abs(i__2));
        i__2 = mrows - nsetb;
        /* Computing MIN */
        i__3 = nsetb + 1;
        /* Computing MIN */
        i__4 = nsetb + 1;
        ww[j] = ddot_(&i__2, &w[min(i__3,mrows) + j * w_dim1], &c__1, &w[min(
                          i__4,mrows) + (*ncols + 1) * w_dim1], &c__1) * (d__1 = scl[
                                      jcol], abs(d__1));
        /* L200: */
    }
    if (iprint > 0) {
        dvout_(ncols, &ww[1], "(' GRADIENT VALUES')", &c_n4, (ftnlen)20);
        ivout_(ncols, &ibasis[1], "(' INTERNAL VARIABLE ORDER')", &c_n4, (
                   ftnlen)28);
        ivout_(ncols, &ibb[1], "(' BOUND POLARITY')", &c_n4, (ftnlen)19);
    }
L210:

    /*     IF ACTIVE SET = NUMBER OF TOTAL ROWS, QUIT. */
    if (nsetb == mrows) {
        found = FALSE_;
        /*     EXIT PROCEDURE */
        goto L90;
    }

    /*     CHOOSE AN EXTREMAL COMPONENT OF GRADIENT VECTOR */
    /*     FOR A CANDIDATE TO BECOME NON-ACTIVE. */
    wlarge = -big;
    wmag = -big;
    i__1 = *ncols;
    for (j = nsetb + 1; j <= i__1; ++j) {
        t = ww[j];
        if (t == big) {
            goto L220;
        }
        itemp = ibasis[j];
        jcol = abs(itemp);
        i__2 = mval - nsetb;
        /* Computing MIN */
        i__3 = nsetb + 1;
        t1 = dnrm2_(&i__2, &w[min(i__3,mrows) + j * w_dim1], &c__1);
        if (itemp < 0) {
            if (ibb[jcol] % 2 == 0) {
                t = -t;
            }
            if (t < zero) {
                goto L220;
            }
            if (mval > nsetb) {
                t = t1;
            }
            if (t > wlarge) {
                wlarge = t;
                jlarge = j;
            }
        } else {
            if (mval > nsetb) {
                t = t1;
            }
            if (abs(t) > wmag) {
                wmag = abs(t);
                jmag = j;
            }
        }
L220:
        ;
    }

    /*     CHOOSE MAG. OF LARGEST COMPONENT OF GRADIENT FOR CANDIDATE. */
    jbig = 0;
    wbig = zero;
    if (wlarge > zero) {
        jbig = jlarge;
        wbig = wlarge;
    }
    if (wmag >= wbig) {
        jbig = jmag;
        wbig = wmag;
    }
    if (jbig == 0) {
        found = FALSE_;
        if (iprint > 0) {
            ivout_(&c__0, &i__, "(' FOUND NO VARIABLE TO ENTER')", &c_n4, (
                       ftnlen)31);
        }
        /*     EXIT PROCEDURE */
        goto L90;
    }

    /*     SEE IF THE INCOMING COL. IS SUFFICIENTLY INDEPENDENT. */
    /*     THIS TEST IS MADE BEFORE AN ELIMINATION IS PERFORMED. */
    if (iprint > 0) {
        ivout_(&c__1, &jbig, "(' TRY TO BRING IN THIS COL.')", &c_n4, (ftnlen)
               30);
    }
    if (mval <= nsetb) {
        cl1 = dnrm2_(&mval, &w[jbig * w_dim1 + 1], &c__1);
        i__1 = nsetb - mval;
        /* Computing MIN */
        i__2 = mval + 1;
        cl2 = abs(wt) * dnrm2_(&i__1, &w[min(i__2,mrows) + jbig * w_dim1], &
                               c__1);
        i__1 = mrows - nsetb;
        /* Computing MIN */
        i__2 = nsetb + 1;
        cl3 = abs(wt) * dnrm2_(&i__1, &w[min(i__2,mrows) + jbig * w_dim1], &
                               c__1);
        drotg_(&cl1, &cl2, &dc, &ds);
        colabv = abs(cl1);
        colblo = cl3;
    } else {
        cl1 = dnrm2_(&nsetb, &w[jbig * w_dim1 + 1], &c__1);
        i__1 = mval - nsetb;
        /* Computing MIN */
        i__2 = nsetb + 1;
        cl2 = dnrm2_(&i__1, &w[min(i__2,mrows) + jbig * w_dim1], &c__1);
        i__1 = mrows - mval;
        /* Computing MIN */
        i__2 = mval + 1;
        cl3 = abs(wt) * dnrm2_(&i__1, &w[min(i__2,mrows) + jbig * w_dim1], &
                               c__1);
        colabv = cl1;
        drotg_(&cl2, &cl3, &dc, &ds);
        colblo = abs(cl2);
    }
    if (colblo <= tolind * colabv) {
        ww[jbig] = big;
        if (iprint > 0) {
            ivout_(&c__0, &i__, "(' VARIABLE IS DEPENDENT, NOT USED.')", &
                   c_n4, (ftnlen)37);
        }
        goto L210;
    }

    /*     SWAP MATRIX COLS. NSETB+1 AND JBIG, PLUS POINTER INFO., AND */
    /*     GRADIENT VALUES. */
    ++nsetb;
    if (nsetb != jbig) {
        dswap_(&mrows, &w[nsetb * w_dim1 + 1], &c__1, &w[jbig * w_dim1 + 1], &
               c__1);
        dswap_(&c__1, &ww[nsetb], &c__1, &ww[jbig], &c__1);
        itemp = ibasis[nsetb];
        ibasis[nsetb] = ibasis[jbig];
        ibasis[jbig] = itemp;
    }

    /*     ELIMINATE ENTRIES BELOW THE PIVOT LINE IN COL. NSETB. */
    if (mrows > nsetb) {
        i__1 = nsetb + 1;
        for (i__ = mrows; i__ >= i__1; --i__) {
            if (i__ == mval + 1) {
                goto L230;
            }
            drotg_(&w[i__ - 1 + nsetb * w_dim1], &w[i__ + nsetb * w_dim1], &
                   sc, &ss);
            w[i__ + nsetb * w_dim1] = zero;
            i__2 = *ncols - nsetb + 1;
            drot_(&i__2, &w[i__ - 1 + (nsetb + 1) * w_dim1], mdw, &w[i__ + (
                        nsetb + 1) * w_dim1], mdw, &sc, &ss);
L230:
            ;
        }
        if (mval >= nsetb && mval < mrows) {
            drotg_(&w[nsetb + nsetb * w_dim1], &w[mval + 1 + nsetb * w_dim1],
                   &sc, &ss);
            w[mval + 1 + nsetb * w_dim1] = zero;
            i__1 = *ncols - nsetb + 1;
            drot_(&i__1, &w[nsetb + (nsetb + 1) * w_dim1], mdw, &w[mval + 1 +
                    (nsetb + 1) * w_dim1], mdw, &sc, &ss);
        }
    }
    if (w[nsetb + nsetb * w_dim1] == zero) {
        ww[nsetb] = big;
        --nsetb;
        if (iprint > 0) {
            ivout_(&c__0, &i__, "(' PIVOT IS ZERO, NOT USED.')", &c_n4, (
                       ftnlen)29);
        }
        goto L210;
    }

    /*     CHECK THAT NEW VARIABLE IS MOVING IN THE RIGHT DIRECTION. */
    itemp = ibasis[nsetb];
    jcol = abs(itemp);
    xnew = w[nsetb + (*ncols + 1) * w_dim1] / w[nsetb + nsetb * w_dim1] / (
               d__1 = scl[jcol], abs(d__1));
    /* CONT: DO BLOCK */
    /* QUIT: DO BLOCK */
    if (itemp < 0) {
        /*     IF(WW(NSETB).GE.ZERO.AND.XNEW.LE.ZERO) EXIT(QUIT) */
        /*     IF(WW(NSETB).LE.ZERO.AND.XNEW.GE.ZERO) EXIT(QUIT) */
        if (ww[nsetb] >= zero && xnew <= zero) {
            goto L240;
        }
        if (ww[nsetb] <= zero && xnew >= zero) {
            goto L240;
        }
    }
    /*     EXIT(CONT) */
    goto L250;
    /*     END BLOCK */
L240:
    ww[nsetb] = big;
    --nsetb;
    if (iprint > 0) {
        ivout_(&c__0, &i__, "(' VARIABLE HAS BAD DIRECTION, NOT USED.')", &
               c_n4, (ftnlen)42);
    }
    goto L210;
    /*     END BLOCK */
L250:
    found = TRUE_;
    /*     EXIT PROCEDURE */
    goto L260;
L260:
    /*     END PROCEDURE */
    goto L90;
    /*     PROCEDURE(SOLVE THE TRIANGULAR SYSTEM) */
L270:
    dcopy_(&nsetb, &w[(*ncols + 1) * w_dim1 + 1], &c__1, &rw[1], &c__1);
    for (j = nsetb; j >= 1; --j) {
        rw[j] /= w[j + j * w_dim1];
        jcol = (i__1 = ibasis[j], abs(i__1));
        t = rw[j];
        if (ibb[jcol] % 2 == 0) {
            rw[j] = -rw[j];
        }
        i__1 = j - 1;
        d__1 = -t;
        daxpy_(&i__1, &d__1, &w[j * w_dim1 + 1], &c__1, &rw[1], &c__1);
        rw[j] /= (d__1 = scl[jcol], abs(d__1));
        /* L280: */
    }
    if (iprint > 0) {
        dvout_(&nsetb, &rw[1], "(' SOLN. VALUES')", &c_n4, (ftnlen)17);
        ivout_(&nsetb, &ibasis[1], "(' COLS. USED')", &c_n4, (ftnlen)15);
    }
    /*     END PROCEDURE */
    switch (lgopr) {
    case 1:
        goto L300;
    case 2:
        goto L440;
    }
    /*     PROCEDURE(MAKE MOVE AND UPDATE FACTORIZATION) */
L290:
    /*     DO(SOLVE THE TRIANGULAR SYSTEM) */
    lgopr = 1;
    goto L270;
L300:

    /*     SEE IF THE UNCONSTRAINED SOL. (OBTAINED BY SOLVING THE */
    /*     TRIANGULAR SYSTEM) SATISFIES THE PROBLEM BOUNDS. */
    alpha = two;
    beta = two;
    x[nsetb] = zero;
    i__1 = nsetb;
    for (j = 1; j <= i__1; ++j) {
        itemp = ibasis[j];
        jcol = abs(itemp);
        t1 = two;
        t2 = two;
        if (itemp < 0) {
            bou = zero;
        } else {
            bou = bl[jcol];
        }
        if (-bou != big) {
            bou /= (d__1 = scl[jcol], abs(d__1));
        }
        if (rw[j] <= bou) {
            t1 = (x[j] - bou) / (x[j] - rw[j]);
        }
        bou = bu[jcol];
        if (bou != big) {
            bou /= (d__1 = scl[jcol], abs(d__1));
        }
        if (rw[j] >= bou) {
            t2 = (bou - x[j]) / (rw[j] - x[j]);
        }

        /*     IF NOT, THEN COMPUTE A STEP LENGTH SO THAT THE */
        /*     VARIABLES REMAIN FEASIBLE. */
        if (t1 < alpha) {
            alpha = t1;
            jdrop1 = j;
        }
        if (t2 < beta) {
            beta = t2;
            jdrop2 = j;
        }
        /* L310: */
    }
    constr = alpha < two || beta < two;
    if (constr) {
        goto L320;
    }

    /*     ACCEPT THE CANDIDATE BECAUSE IT SATISFIES THE STATED BOUNDS */
    /*     ON THE VARIABLES. */
    dcopy_(&nsetb, &rw[1], &c__1, &x[1], &c__1);
    goto L120;
L320:

    /*     TAKE A STEP THAT IS AS LARGE AS POSSIBLE WITH ALL */
    /*     VARIABLES REMAINING FEASIBLE. */
    i__1 = nsetb;
    for (j = 1; j <= i__1; ++j) {
        x[j] += min(alpha,beta) * (rw[j] - x[j]);
        /* L330: */
    }
    if (alpha <= beta) {
        jdrop2 = 0;
    } else {
        jdrop1 = 0;
    }
L340:
    if (jdrop1 + jdrop2 > 0 && nsetb > 0) {
        goto L350;
    }
    goto L460;
L350:
    jdrop = jdrop1 + jdrop2;
    itemp = ibasis[jdrop];
    jcol = abs(itemp);
    if (jdrop2 > 0) {

        /*     VARIABLE IS AT AN UPPER BOUND.  SUBTRACT MULTIPLE OF THIS COL. */
        /*     FROM RIGHT HAND SIDE. */
        t = bu[jcol];
        if (itemp > 0) {
            bu[jcol] = t - bl[jcol];
            bl[jcol] = -t;
            itemp = -itemp;
            scl[jcol] = -scl[jcol];
            i__1 = jdrop;
            for (i__ = 1; i__ <= i__1; ++i__) {
                w[i__ + jdrop * w_dim1] = -w[i__ + jdrop * w_dim1];
                /* L360: */
            }
        } else {
            ++ibb[jcol];
            if (ibb[jcol] % 2 == 0) {
                t = -t;
            }
        }
        /*     VARIABLE IS AT A LOWER BOUND. */
    } else {
        if ((doublereal) itemp < zero) {
            t = zero;
        } else {
            t = -bl[jcol];
            bu[jcol] += t;
            itemp = -itemp;
        }
    }
    daxpy_(&jdrop, &t, &w[jdrop * w_dim1 + 1], &c__1, &w[(*ncols + 1) *
            w_dim1 + 1], &c__1);

    /*     MOVE CERTAIN COLS. LEFT TO ACHIEVE UPPER HESSENBERG FORM. */
    dcopy_(&jdrop, &w[jdrop * w_dim1 + 1], &c__1, &rw[1], &c__1);
    i__1 = nsetb;
    for (j = jdrop + 1; j <= i__1; ++j) {
        ibasis[j - 1] = ibasis[j];
        x[j - 1] = x[j];
        dcopy_(&j, &w[j * w_dim1 + 1], &c__1, &w[(j - 1) * w_dim1 + 1], &c__1)
        ;
        /* L370: */
    }
    ibasis[nsetb] = itemp;
    w[nsetb * w_dim1 + 1] = zero;
    i__1 = mrows - jdrop;
    dcopy_(&i__1, &w[nsetb * w_dim1 + 1], &c__0, &w[jdrop + 1 + nsetb *
            w_dim1], &c__1);
    dcopy_(&jdrop, &rw[1], &c__1, &w[nsetb * w_dim1 + 1], &c__1);

    /*     TRANSFORM THE MATRIX FROM UPPER HESSENBERG FORM TO */
    /*     UPPER TRIANGULAR FORM. */
    --nsetb;
    /* SMLL: */
    /*     *DO BLOCK */
    /* NRML: */
    /*     *DO BLOCK */
    i__1 = nsetb;
    for (i__ = jdrop; i__ <= i__1; ++i__) {

        /*     LOOK FOR SMALL PIVOTS AND AVOID MIXING WEIGHTED AND */
        /*     NONWEIGHTED ROWS. */
        if (i__ == mval) {
            t = zero;
            i__2 = nsetb;
            for (j = i__; j <= i__2; ++j) {
                jcol = (i__3 = ibasis[j], abs(i__3));
                t1 = (d__1 = w[i__ + j * w_dim1] * scl[jcol], abs(d__1));
                if (t1 > t) {
                    jbig = j;
                    t = t1;
                }
                /* L380: */
            }
            /*     EXIT(NRML) */
            goto L400;
        }
        drotg_(&w[i__ + i__ * w_dim1], &w[i__ + 1 + i__ * w_dim1], &sc, &ss);
        w[i__ + 1 + i__ * w_dim1] = zero;
        i__2 = *ncols - i__ + 1;
        drot_(&i__2, &w[i__ + (i__ + 1) * w_dim1], mdw, &w[i__ + 1 + (i__ + 1)
                * w_dim1], mdw, &sc, &ss);
        /* L390: */
    }
    /*     EXIT (SMLL) */
    goto L430;
    /*     END BLOCK */
L400:

    /*     THE TRIANGULARIZATION IS COMPLETED BY GIVING UP */
    /*     THE HESSENBERG FORM AND TRIANGULARIZING A RECTANGULAR MATRIX. */
    dswap_(&mrows, &w[i__ * w_dim1 + 1], &c__1, &w[jbig * w_dim1 + 1], &c__1);
    dswap_(&c__1, &ww[i__], &c__1, &ww[jbig], &c__1);
    dswap_(&c__1, &x[i__], &c__1, &x[jbig], &c__1);
    itemp = ibasis[i__];
    ibasis[i__] = ibasis[jbig];
    ibasis[jbig] = itemp;
    jbig = i__;
    i__1 = nsetb;
    for (j = jbig; j <= i__1; ++j) {
        i__2 = mrows;
        for (i__ = j + 1; i__ <= i__2; ++i__) {
            drotg_(&w[j + j * w_dim1], &w[i__ + j * w_dim1], &sc, &ss);
            w[i__ + j * w_dim1] = zero;
            i__3 = *ncols - j + 1;
            drot_(&i__3, &w[j + (j + 1) * w_dim1], mdw, &w[i__ + (j + 1) *
                    w_dim1], mdw, &sc, &ss);
            /* L410: */
        }
        /* L420: */
    }
    /*     END BLOCK */
L430:

    /*     SEE IF THE REMAINING COEFFICIENTS ARE FEASIBLE.  THEY SHOULD */
    /*     BE BECAUSE OF THE WAY MIN(ALPHA,BETA) WAS CHOSEN.  ANY THAT ARE */
    /*     NOT FEASIBLE WILL BE SET TO THEIR BOUNDS AND */
    /*     APPROPRIATELY TRANSLATED. */
    jdrop1 = 0;
    jdrop2 = 0;
    /*     DO(SOLVE THE TRIANGULAR SYSTEM) */
    lgopr = 2;
    goto L270;
L440:
    dcopy_(&nsetb, &rw[1], &c__1, &x[1], &c__1);
    i__1 = nsetb;
    for (j = 1; j <= i__1; ++j) {
        itemp = ibasis[j];
        jcol = abs(itemp);
        if (itemp < 0) {
            bou = zero;
        } else {
            bou = bl[jcol];
        }
        if (-bou != big) {
            bou /= (d__1 = scl[jcol], abs(d__1));
        }
        if (x[j] <= bou) {
            jdrop1 = j;
            goto L340;
        }
        bou = bu[jcol];
        if (bou != big) {
            bou /= (d__1 = scl[jcol], abs(d__1));
        }
        if (x[j] >= bou) {
            jdrop2 = j;
            goto L340;
        }
        /* L450: */
    }
    goto L340;
L460:
    /*     END PROCEDURE */
    goto L120;
    /*     PROCEDURE(INITIALIZE VARIABLES AND DATA VALUES) */
L470:

    /*     PRETRIANGULARIZE RECTANGULAR ARRAYS OF CERTAIN SIZES */
    /*     FOR INCREASED EFFICIENCY. */
    if (fac * *minput > (doublereal)(*ncols)) {
        i__1 = *ncols + 1;
        for (j = 1; j <= i__1; ++j) {
            i__2 = j + mval + 1;
            for (i__ = *minput; i__ >= i__2; --i__) {
                drotg_(&w[i__ - 1 + j * w_dim1], &w[i__ + j * w_dim1], &sc, &
                       ss);
                w[i__ + j * w_dim1] = zero;
                i__3 = *ncols - j + 1;
                drot_(&i__3, &w[i__ - 1 + (j + 1) * w_dim1], mdw, &w[i__ + (j
                        + 1) * w_dim1], mdw, &sc, &ss);
                /* L480: */
            }
            /* L490: */
        }
        mrows = *ncols + mval + 1;
    } else {
        mrows = *minput;
    }

    /*      SET THE X(*) ARRAY TO ZERO SO ALL COMPONENTS ARE DEFINED. */
    x[1] = zero;
    dcopy_(ncols, &x[1], &c__0, &x[1], &c__1);

    /*     THE ARRAYS IBASIS(*), IBB(*) ARE INITIALIZED BY THE CALLING */
    /*     PROGRAM UNIT. */
    /*     THE COL. SCALING IS DEFINED IN THE CALLING PROGRAM UNIT. */
    /*    'BIG' IS PLUS INFINITY ON THIS MACHINE. */
    big = d1mach_(&lc__2);
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        icase = ind[j];
        /*     DO CASE(ICASE,4) */
        switch (icase) {
        case 1:
            goto L500;
        case 2:
            goto L510;
        case 3:
            goto L520;
        case 4:
            goto L530;
        }
        goto L540;
        /*     CASE 1 */
L500:
        bu[j] = big;
        goto L540;
        /*     CASE 2 */
L510:
        bl[j] = -big;
        goto L540;
        /*     CASE 3 */
L520:
        goto L540;
        /*     CASE 4 */
L530:
        bl[j] = -big;
        bu[j] = big;
        /*     END CASE */
L540:
        /* L550: */
        ;
    }
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        if ((bl[j] <= zero && zero <= bu[j] &&
            (d__1 = bu[j], abs(d__1)) < (d__2 = bl[j], abs(d__2))) ||
            bu[j] < zero) {
            t = bu[j];
            bu[j] = -bl[j];
            bl[j] = -t;
            scl[j] = -scl[j];
            i__2 = mrows;
            for (i__ = 1; i__ <= i__2; ++i__) {
                w[i__ + j * w_dim1] = -w[i__ + j * w_dim1];
                /* L560: */
            }
        }

        /*     INDICES IN SET T(=TIGHT) ARE DENOTED BY NEGATIVE VALUES */
        /*     OF IBASIS(*). */
        if (bl[j] >= zero) {
            ibasis[j] = -ibasis[j];
            t = -bl[j];
            bu[j] += t;
            daxpy_(&mrows, &t, &w[j * w_dim1 + 1], &c__1, &w[(*ncols + 1) *
                    w_dim1 + 1], &c__1);
        }
        /* L570: */
    }
    nsetb = 0;
    iter = 0;
    /*     END PROCEDURE */
    goto L50;
    /*     PROCEDURE(PROCESS OPTION ARRAY) */
L580:
    zero = 0.;
    fac = .75;
    one = 1.;
    two = 2.;
    tolind = sqrt(d1mach_(&lc__4));
    tolsze = sqrt(d1mach_(&lc__4));
    itmax = max(mrows,*ncols) * 5;
    wt = one;
    mval = 0;
    iprint = 0;

    /*     CHANGES TO SOME PARAMETERS CAN OCCUR THROUGH THE OPTION */
    /*     ARRAY, IOPT(*).  PROCESS THIS ARRAY LOOKING CAREFULLY */
    /*     FOR INPUT DATA ERRORS. */
    lp = 0;
    lds = 0;
L590:
    lp += lds;

    /*     TEST FOR NO MORE OPTIONS. */
    ip = iopt[lp + 1];
    jp = abs(ip);
    if (ip == 99) {
        goto L600;
    } else if (jp == 99) {
        lds = 1;
        goto L590;
    } else if (jp == 1) {

        /*     MOVE THE IOPT(*) PROCESSING POINTER. */
        if (ip > 0) {
            lp = iopt[lp + 2] - 1;
            lds = 0;
        } else {
            lds = 2;
        }
        goto L590;
    } else if (jp == 2) {

        /*     CHANGE TOLERANCE FOR RANK DETERMINATION. */
        if (ip > 0) {
            ioff = iopt[lp + 2];
            if (ioff <= 0) {
                nerr = 24;
                nchar = 89;
                xerrwv_("DBOLSM(). THE OFFSET=(I1) BEYOND POSTION NCOLS=(I2)"
                        " MUST BE POSITIVE FOR OPTION NUMBER 2.", &nchar, &
                        nerr, &level, &c__2, &ioff, ncols, &c__0, &rdum, &
                        rdum, (ftnlen)89);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L610;
            }
            tolind = x[*ncols + ioff];
            if (tolind < d1mach_(&lc__4)) {
                nerr = 25;
                nlevel = 0;
                nchar = 88;
                rdum2 = (real) tolind;
                rdum = (real) d1mach_(&lc__4);
                xerrwv_("DBOLSM(). THE TOLERANCE FOR RANK DETERMINATION=(R1)"
                        " IS LESS THAN MACHINE PRECISION=(R2).", &nchar, &nerr,
                        &nlevel, &c__0, &idum, &idum, &c__2, &rdum2, &rdum, (
                            ftnlen)88);
            }
        }
        lds = 2;
        goto L590;
    } else if (jp == 3) {

        /*     CHANGE BLOWUP FACTOR FOR ALLOWING VARIABLES TO BECOME */
        /*     INACTIVE. */
        if (ip > 0) {
            ioff = iopt[lp + 2];
            if (ioff <= 0) {
                nerr = 26;
                nchar = 89;
                xerrwv_("DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2"
                        ") MUST BE POSTIVE FOR OPTION NUMBER 3.", &nchar, &
                        nerr, &level, &c__2, &ioff, ncols, &c__0, &rdum, &
                        rdum, (ftnlen)89);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L610;
            }
            tolsze = x[*ncols + ioff];
            if (tolsze <= zero) {
                nerr = 27;
                rdum2 = (real) tolsze;
                xerrwv_("DBOLSM(). THE RECIPROCAL OF THE BLOW-UP FACTOR FOR "
                        "REJECTING VARIABLES MUST BE POSITIVE. NOW=(R1).", &
                        nchar, &nerr, &level, &c__0, &idum, &idum, &c__1, &
                        rdum2, &rdum, (ftnlen)98);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L610;
            }
        }
        lds = 2;
        goto L590;
    } else if (jp == 4) {

        /*     CHANGE THE MAX. NO. OF ITERATIONS ALLOWED. */
        if (ip > 0) {
            itmax = iopt[lp + 2];
            if (itmax <= 0) {
                nerr = 28;
                nchar = 65;
                xerrwv_("DBOLSM(). THE MAXIMUM NUMBER OF ITERATIONS=(I1) MUS"
                        "T BE POSITIVE.", &nchar, &nerr, &level, &c__1, &itmax,
                        &idum, &c__0, &rdum, &rdum, (ftnlen)65);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L610;
            }
        }
        lds = 2;
        goto L590;
    } else if (jp == 5) {

        /*     CHANGE THE FACTOR FOR PRETRIANGULARIZING THE DATA MATRIX. */
        if (ip > 0) {
            ioff = iopt[lp + 2];
            if (ioff <= 0) {
                nerr = 29;
                nchar = 89;
                xerrwv_("DBOLSM(). THE OFFSET=(I1) BEYOND POSITION NCOLS=(I2"
                        ") MUST BE POSTIVE FOR OPTION NUMBER 5.", &nchar, &
                        nerr, &level, &c__2, &ioff, ncols, &c__0, &rdum, &
                        rdum, (ftnlen)89);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L610;
            }
            fac = x[*ncols + ioff];
            if (fac < zero) {
                nerr = 30;
                nlevel = 0;
                nchar = 104;
                rdum2 = (real) fac;
                xerrwv_("DBOLSM(). THE FACTOR (NCOLS/MROWS) WHERE PRE-TRIANG"
                        "ULARIZING IS PERFORMED MUST BE NONNEGATIVE. NOW=(R1)."
                        , &nchar, &nerr, &nlevel, &c__0, &idum, &idum, &c__1,
                        &rdum2, &rdum, (ftnlen)104);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L610;
            }
        }
        lds = 2;
        goto L590;
    } else if (jp == 6) {

        /*     CHANGE THE WEIGHTING FACTOR (FROM ONE) TO APPLY TO COMPONENTS */
        /*     NUMBERED .GT. MVAL (INITIALLY SET TO 1.)  THIS TRICK IS NEEDED */
        /*     FOR APPLICATIONS OF THIS SUBPROGRAM TO THE HEAVILY WEIGHTED */
        /*     LEAST SQUARES PROBLEM THAT COME FROM EQUALITY CONSTRAINTS. */
        if (ip > 0) {
            ioff = iopt[lp + 2];
            mval = iopt[lp + 3];
            wt = x[*ncols + ioff];
        }
        if (mval < 0 || mval > *minput || wt <= zero) {
            nerr = 38;
            nlevel = 0;
            nchar = 116;
            rdum2 = (real) wt;
            xerrwv_("DBOLSM(). THE ROW SEPARATOR TO APPLY WEIGHTING (I1) MUS"
                    "T LIE BETWEEN 0 AND MROWS (I2). WEIGHT (R1) MUST BE POSI"
                    "TIVE.", &nchar, &nerr, &nlevel, &c__2, &mval, minput, &
                    c__1, &rdum2, &rdum, (ftnlen)116);
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L610;
        }
        lds = 3;
        goto L590;

        /*     TURN ON DEBUG OUTPUT. */
    } else if (jp == 7) {
        if (ip > 0) {
            iprint = 1;
        }
        lds = 2;
        goto L590;
    } else {
        nerr = 23;
        nchar = 46;
        xerrwv_("DBOLSM. THE OPTION NUMBER=(I1) IS NOT DEFINED.", &nchar, &
                nerr, &level, &c__1, &ip, &idum, &c__0, &rdum, &rdum, (ftnlen)
                46);
        /*     DO(RETURN TO USER PROGRAM UNIT) */
        goto L610;
    }
L600:
    /*     END PROCEDURE */
    goto L40;
    /*     PROCEDURE(RETURN TO USER PROGRAM UNIT) */
L610:
    *mode = -nerr;
    return 0;
    /*     END PROCEDURE */
    /*     END PROGRAM */
} /* dbolsm_ */

