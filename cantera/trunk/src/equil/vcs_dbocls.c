/* dbocls.f -- translated by f2c (version 20031025).
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
static integer c__4 = 4;
static long int lc__4 = 4;

/* DECK DBOCLS */
/* Subroutine */
int dbocls_(doublereal* w, integer* mdw, integer* mcon,
            integer* mrows, integer* ncols, doublereal* bl, doublereal* bu,
            integer* ind, integer* iopt, doublereal* x, doublereal* rnormc,
            doublereal* rnorm, integer* mode, doublereal* rw, integer* iw)
{
    /* Initialized data */

    static integer igo = 0;

    /* System generated locals */
    integer w_dim1, w_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Local variables */
    static integer i__, j;
    static doublereal t, t1, t2;
    static integer ip, jp, lp;
    static doublereal wt;
    static integer llb;
    static doublereal one;
    static integer lds, iiw, liw, llx, irw, lrw;
    extern doublereal ddot_(integer*, doublereal*, integer*, doublereal*,
                            integer*);
    static integer idum, lbou, lmdw, lndw, mdwl, nerr;
    static real rdum;
    static integer lenx, lliw, mnew, jopt[5], lopt;
    static doublereal zero;
    static integer mopt, llrw, mout;
    extern doublereal dnrm2_(integer*, doublereal*, integer*);
    static real rdum2;
    static integer icase;
    extern /* Subroutine */ int dscal_(integer*, doublereal*, doublereal*,
                                       integer*);
    static integer modec;
    static logical accum;
    static integer nchar;
    extern /* Subroutine */ int dbols_(doublereal*, integer*, integer*,
                                       integer*, doublereal*, doublereal*, integer*, integer*,
                                       doublereal*, doublereal*, integer*, doublereal*, integer*);
    static integer level;
    extern doublereal dasum_(integer*, doublereal*, integer*);
    static doublereal anorm, cnorm;
    static integer lboum;
    extern /* Subroutine */ int dcopy_(integer*, doublereal*, integer*,
                                       doublereal*, integer*);
    static integer liopt;
    extern doublereal d1mach_(long int *);
    static integer locacc;
    static logical checkl;
    static integer iscale, locdim;
    static logical filter;
    static doublereal drelpr;
    static logical pretri;
    static integer inrows;
    extern /* Subroutine */ int xerrwv_(char*, integer*, integer*, integer
                                        *, integer*, integer*, integer*, integer*, real*, real*,
                                        ftnlen);

    /* ***BEGIN PROLOGUE  DBOCLS */
    /* ***DATE WRITTEN   821220   (YYMMDD) */
    /* ***REVISION DATE  870803   (YYMMDD) */
    /* ***CATEGORY NO.  K1A2A,G2E,G2H1,G2H2 */
    /* ***KEYWORDS  LIBRARY=SLATEC,TYPE=DOUBLE PRECISION(SBOCLS-S DBOCLS-D), */
    /*             BOUNDS,CONSTRAINTS,INEQUALITY,LEAST SQUARES,LINEAR */
    /* ***AUTHOR  HANSON, R. J., SNLA */
    /* ***PURPOSE  Solve the bounded and constrained least squares */
    /*            problem consisting of solving the equation */
    /*                      E*X = F  (in the least squares sense) */
    /*             subject to the linear constraints */
    /*                            C*X = Y. */
    /* ***DESCRIPTION */

    /*            **** Double Precision Version of SBOCLS **** */
    /*   **** All INPUT and OUTPUT real variables are DOUBLE PRECISION **** */

    /*     This subprogram solves the bounded and constrained least squares */
    /*     problem. The problem statement is: */

    /*     Solve E*X = F (least squares sense), subject to constraints */
    /*     C*X=Y. */

    /*     In this formulation both X and Y are unknowns, and both may */
    /*     have bounds on any of their components.  This formulation */
    /*     of the problem allows the user to have equality and inequality */
    /*     constraints as well as simple bounds on the solution components. */

    /*     This constrained linear least squares subprogram solves E*X=F */
    /*     subject to C*X=Y, where E is MROWS by NCOLS, C is MCON by NCOLS. */

    /*      The user must have dimension statements of the form */

    /*      DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON), BU(NCOLS+MCON), */
    /*     * X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON) */
    /*       INTEGER IND(NCOLS+MCON), IOPT(17+NI), IW(2*(NCOLS+MCON)) */

    /*     (here NX=number of extra locations required for the options; NX=0 */
    /*     if no options are in use. Also NI=number of extra locations */
    /*     for options 1-9.) */

    /*    INPUT */
    /*    ----- */

    /*    ------------------------- */
    /*    W(MDW,*),MCON,MROWS,NCOLS */
    /*    ------------------------- */
    /*     The array W contains the (possibly null) matrix [C:*] followed by */
    /*     [E:F].  This must be placed in W as follows: */
    /*          [C  :  *] */
    /*     W  = [       ] */
    /*          [E  :  F] */
    /*     The (*) after C indicates that this data can be undefined. The */
    /*     matrix [E:F] has MROWS rows and NCOLS+1 columns. The matrix C is */
    /*     placed in the first MCON rows of W(*,*) while [E:F] */
    /*     follows in rows MCON+1 through MCON+MROWS of W(*,*). The vector F */
    /*     is placed in rows MCON+1 through MCON+MROWS, column NCOLS+1. The */
    /*     values of MDW and NCOLS must be positive; the value of MCON must */
    /*     be nonnegative. An exception to this occurs when using option 1 */
    /*     for accumulation of blocks of equations. In that case MROWS is an */
    /*     OUTPUT variable only, and the matrix data for [E:F] is placed in */
    /*     W(*,*), one block of rows at a time. See IOPT(*) contents, option */
    /*     number 1, for further details. The row dimension, MDW, of the */
    /*     array W(*,*) must satisfy the inequality: */

    /*     If using option 1, */
    /*                     MDW .ge. MCON + max(max. number of */
    /*                     rows accumulated, NCOLS) + 1. */
    /*     If using option 8, */
    /*                     MDW .ge. MCON + MROWS. */
    /*     Else */
    /*                     MDW .ge. MCON + max(MROWS, NCOLS). */

    /*     Other values are errors, but this is checked only when using */
    /*     option=2.  The value of MROWS is an output parameter when */
    /*     using option number 1 for accumulating large blocks of least */
    /*     squares equations before solving the problem. */
    /*     See IOPT(*) contents for details about option 1. */

    /*    ------------------ */
    /*    BL(*),BU(*),IND(*) */
    /*    ------------------ */
    /*     These arrays contain the information about the bounds that the */
    /*     solution values are to satisfy. The value of IND(J) tells the */
    /*     type of bound and BL(J) and BU(J) give the explicit values for */
    /*     the respective upper and lower bounds on the unknowns X and Y. */
    /*     The first NVARS entries of IND(*), BL(*) and BU(*) specify */
    /*     bounds on X; the next MCON entries specify bounds on Y. */

    /*    1.    For IND(J)=1, require X(J) .ge. BL(J); */
    /*          IF J.gt.NCOLS,        Y(J-NCOLS) .ge. BL(J). */
    /*          (the value of BU(J) is not used.) */
    /*    2.    For IND(J)=2, require X(J) .le. BU(J); */
    /*          IF J.gt.NCOLS,        Y(J-NCOLS) .le. BU(J). */
    /*          (the value of BL(J) is not used.) */
    /*    3.    For IND(J)=3, require X(J) .ge. BL(J) and */
    /*                                X(J) .le. BU(J); */
    /*          IF J.gt.NCOLS,        Y(J-NCOLS) .ge. BL(J) and */
    /*                                Y(J-NCOLS) .le. BU(J). */
    /*          (to impose equality constraints have BL(J)=BU(J)= */
    /*          constraining value.) */
    /*    4.    For IND(J)=4, no bounds on X(J) or Y(J-NCOLS) are required. */
    /*          (the values of BL(J) and BU(J) are not used.) */

    /*     Values other than 1,2,3 or 4 for IND(J) are errors. In the case */
    /*     IND(J)=3 (upper and lower bounds) the condition BL(J) .gt. BU(J) */
    /*     is  an  error.   The values BL(J), BU(J), J .gt. NCOLS, will be */
    /*     changed.  Significant changes mean that the constraints are */
    /*     infeasible.  (Users must make this decision themselves.) */
    /*     The new values for BL(J), BU(J), J .gt. NCOLS, define a */
    /*     region such that the perturbed problem is feasible.  If users */
    /*     know that their problem is feasible, this step can be skipped */
    /*     by using option number 8 described below. */
    /*     See IOPT(*) description. */


    /*    ------- */
    /*    IOPT(*) */
    /*    ------- */
    /*     This is the array where the user can specify nonstandard options */
    /*     for DBOCLS( ). Most of the time this feature can be ignored by */
    /*     setting the input value IOPT(1)=99. Occasionally users may have */
    /*     needs that require use of the following subprogram options. For */
    /*     details about how to use the options see below: IOPT(*) CONTENTS. */

    /*     Option Number   Brief Statement of Purpose */
    /*     ------ ------   ----- --------- -- ------- */
    /*           1         Return to user for accumulation of blocks */
    /*                     of least squares equations.  The values */
    /*                     of IOPT(*) are changed with this option. */
    /*                     The changes are updates to pointers for */
    /*                     placing the rows of equations into position */
    /*                     for processing. */
    /*           2         Check lengths of all arrays used in the */
    /*                     subprogram. */
    /*           3         Column scaling of the data matrix, [C]. */
    /*                                                        [E] */
    /*           4         User provides column scaling for matrix [C]. */
    /*                                                             [E] */
    /*           5         Provide option array to the low-level */
    /*                     subprogram SBOLS( ). */
    /*           6         Provide option array to the low-level */
    /*                     subprogram SBOLSM( ). */
    /*           7         Move the IOPT(*) processing pointer. */
    /*           8         Do not preprocess the constraints to */
    /*                     resolve infeasibilities. */
    /*           9         Do not pretriangularize the least squares matrix. */
    /*          99         No more options to change. */

    /*    ---- */
    /*    X(*) */
    /*    ---- */
    /*     This array is used to pass data associated with options 4,5 and */
    /*     6. Ignore this parameter (on input) if no options are used. */
    /*     Otherwise see below: IOPT(*) CONTENTS. */


    /*    OUTPUT */
    /*    ------ */

    /*    ----------------- */
    /*    X(*),RNORMC,RNORM */
    /*    ----------------- */
    /*     The array X(*) contains a solution (if MODE .ge.0 or .eq.-22) for */
    /*     the constrained least squares problem. The value RNORMC is the */
    /*     minimum residual vector length for the constraints C*X - Y = 0. */
    /*     The value RNORM is the minimum residual vector length for the */
    /*     least squares equations. Normally RNORMC=0, but in the case of */
    /*     inconsistent constraints this value will be nonzero. */
    /*     The values of X are returned in the first NVARS entries of X(*). */
    /*     The values of Y are returned in the last MCON entries of X(*). */

    /*    ---- */
    /*    MODE */
    /*    ---- */
    /*     The sign of MODE determines whether the subprogram has completed */
    /*     normally, or encountered an error condition or abnormal status. A */
    /*     value of MODE .ge. 0 signifies that the subprogram has completed */
    /*     normally. The value of mode (.ge. 0) is the number of variables */
    /*     in an active status: not at a bound nor at the value zero, for */
    /*     the case of free variables. A negative value of MODE will be one */
    /*     of the cases (-57)-(-41), (-37)-(-22), (-19)-(-2). Values .lt. -1 */
    /*     correspond to an abnormal completion of the subprogram. These */
    /*     error messages are in groups for the subprograms DBOCLS(), */
    /*     SBOLSM(), and SBOLS().  An approximate solution will be returned */
    /*     to the user only when max. iterations is reached, MODE=-22. */

    /*    ----------- */
    /*    RW(*),IW(*) */
    /*    ----------- */
    /*     These are working arrays.  (normally the user can ignore the */
    /*     contents of these arrays.) */

    /*    IOPT(*) CONTENTS */
    /*    ------- -------- */
    /*     The option array allows a user to modify some internal variables */
    /*     in the subprogram without recompiling the source code. A central */
    /*     goal of the initial software design was to do a good job for most */
    /*     people. Thus the use of options will be restricted to a select */
    /*     group of users. The processing of the option array proceeds as */
    /*     follows: a pointer, here called LP, is initially set to the value */
    /*     1. At the pointer position the option number is extracted and */
    /*     used for locating other information that allows for options to be */
    /*     changed. The portion of the array IOPT(*) that is used for each */
    /*     option is fixed; the user and the subprogram both know how many */
    /*     locations are needed for each option. The value of LP is updated */
    /*     for each option based on the amount of storage in IOPT(*) that is */
    /*     required. A great deal of error checking is done by the */
    /*     subprogram on the contents of the option array. Nevertheless it */
    /*     is still possible to give the subprogram optional input that is */
    /*     meaningless. For example option 4 uses the locations */
    /*     X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) for passing scaling data. */
    /*     The user must manage the allocation of these locations. */

    /*   1 */
    /*   - */
    /*     This option allows the user to solve problems with a large number */
    /*     of rows compared to the number of variables. The idea is that the */
    /*     subprogram returns to the user (perhaps many times) and receives */
    /*     new least squares equations from the calling program unit. */
    /*     Eventually the user signals "that's all" and a solution is then */
    /*     computed. The value of MROWS is an output variable when this */
    /*     option is used. Its value is always in the range 0 .le. MROWS */
    /*     .le. NCOLS+1. It is the number of rows after the */
    /*     triangularization of the entire set of equations. If LP is the */
    /*     processing pointer for IOPT(*), the usage for the sequential */
    /*     processing of blocks of equations is */


    /*        IOPT(LP)=1 */
    /*         Move block of equations to W(*,*) starting at */
    /*         the first row of W(*,*). */
    /*        IOPT(LP+3)=# of rows in the block; user defined */

    /*     The user now calls DBOCLS( ) in a loop. The value of IOPT(LP+1) */
    /*     directs the user's action. The value of IOPT(LP+2) points to */
    /*     where the subsequent rows are to be placed in W(*,*). Both of */
    /*     these values are first defined in the subprogram. The user */
    /*     changes the value of IOPT(LP+1) (to 2) as a signal that all of */
    /*     the rows have been processed. */


    /*      .<LOOP */
    /*      . CALL DBOCLS( ) */
    /*      . IF(IOPT(LP+1) .EQ. 1) THEN */
    /*      .    IOPT(LP+3)=# OF ROWS IN THE NEW BLOCK; USER DEFINED */
    /*      .    PLACE NEW BLOCK OF IOPT(LP+3) ROWS IN */
    /*      .    W(*,*) STARTING AT ROW MCON + IOPT(LP+2). */
    /*      . */
    /*      .    IF( THIS IS THE LAST BLOCK OF EQUATIONS ) THEN */
    /*      .       IOPT(LP+1)=2 */
    /*      .<------CYCLE LOOP */
    /*      .    ELSE IF (IOPT(LP+1) .EQ. 2) THEN */
    /*      <-------EXIT LOOP SOLUTION COMPUTED IF MODE .GE. 0 */
    /*      . ELSE */
    /*      . ERROR CONDITION; SHOULD NOT HAPPEN. */
    /*      .<END LOOP */

    /*     Use of this option adds 4 to the required length of IOPT(*). */

    /*   2 */
    /*   - */
    /*     This option is useful for checking the lengths of all arrays used */
    /*     by DBOCLS( ) against their actual requirements for this problem. */
    /*     The idea is simple: the user's program unit passes the declared */
    /*     dimension information of the arrays. These values are compared */
    /*     against the problem-dependent needs within the subprogram. If any */
    /*     of the dimensions are too small an error message is printed and a */
    /*     negative value of MODE is returned, -41 to -47. The printed error */
    /*     message tells how long the dimension should be. If LP is the */
    /*     processing pointer for IOPT(*), */

    /*        IOPT(LP)=2 */
    /*        IOPT(LP+1)=Row dimension of W(*,*) */
    /*        IOPT(LP+2)=Col. dimension of W(*,*) */
    /*        IOPT(LP+3)=Dimensions of BL(*),BU(*),IND(*) */
    /*        IOPT(LP+4)=Dimension of X(*) */
    /*        IOPT(LP+5)=Dimension of RW(*) */
    /*        IOPT(LP+6)=Dimension of IW(*) */
    /*        IOPT(LP+7)=Dimension of IOPT(*) */
    /*         . */
    /*        CALL DBOCLS( ) */

    /*     Use of this option adds 8 to the required length of IOPT(*). */

    /*   3 */
    /*   - */
    /*     This option can change the type of scaling for the data matrix. */
    /*     Nominally each nonzero column of the matrix is scaled so that the */
    /*     magnitude of its largest entry is equal to the value ONE. If LP */
    /*     is the processing pointer for IOPT(*), */

    /*        IOPT(LP)=3 */
    /*        IOPT(LP+1)=1,2 or 3 */
    /*            1= Nominal scaling as noted; */
    /*            2= Each nonzero column scaled to have length ONE; */
    /*            3= Identity scaling; scaling effectively suppressed. */
    /*         . */
    /*        CALL DBOCLS( ) */

    /*     Use of this option adds 2 to the required length of IOPT(*). */

    /*   4 */
    /*   - */
    /*     This options allows the user to provide arbitrary (positive) */
    /*     column scaling for the matrix. If LP is the processing pointer */
    /*     for IOPT(*), */

    /*        IOPT(LP)=4 */
    /*        IOPT(LP+1)=IOFF */
    /*        X(NCOLS+IOFF),...,X(NCOLS+IOFF+NCOLS-1) */
    /*        = Positive scale factors for cols. of E. */
    /*         . */
    /*        CALL DBOCLS( ) */

    /*     Use of this option adds 2 to the required length of IOPT(*) */
    /*     and NCOLS to the required length of X(*). */

    /*   5 */
    /*   - */
    /*     This option allows the user to provide an option array to the */
    /*     low-level subprogram SBOLS( ). If LP is the processing pointer */
    /*     for IOPT(*), */

    /*        IOPT(LP)=5 */
    /*        IOPT(LP+1)= Position in IOPT(*) where option array */
    /*                    data for SBOLS( ) begins. */
    /*         . */
    /*        CALL DBOCLS( ) */

    /*     Use of this option adds 2 to the required length of IOPT(*). */

    /*   6 */
    /*   - */
    /*     This option allows the user to provide an option array to the */
    /*     low-level subprogram SBOLSM( ). If LP is the processing pointer */
    /*     for IOPT(*), */

    /*        IOPT(LP)=6 */
    /*        IOPT(LP+1)= Position in IOPT(*) where option array */
    /*                    data for SBOLSM( ) begins. */
    /*         . */
    /*        CALL DBOCLS( ) */

    /*     Use of this option adds 2 to the required length of IOPT(*). */

    /*   7 */
    /*   - */
    /*     Move the processing pointer (either forward or backward) to the */
    /*     location IOPT(LP+1). The processing pointer moves to locations */
    /*     LP+2 if option number 7 is used with the value -7.  For */
    /*     example to skip over locations 3,...,NCOLS+2, */

    /*       IOPT(1)=7 */
    /*       IOPT(2)=NCOLS+3 */
    /*       (IOPT(I), I=3,...,NCOLS+2 are not defined here.) */
    /*       IOPT(NCOLS+3)=99 */
    /*       CALL DBOCLS( ) */

    /*     CAUTION: Misuse of this option can yield some very hard-to-find */
    /*     bugs. Use it with care. It is intended to be used for passing */
    /*     option arrays to other subprograms. */

    /*   8 */
    /*   - */
    /*     This option allows the user to suppress the algorithmic feature */
    /*     of DBOCLS( ) that processes the constraint equations C*X = Y and */
    /*     resolves infeasibilities. The steps normally done are to solve */
    /*     C*X - Y = 0 in a least squares sense using the stated bounds on */
    /*     both X and Y. Then the "reachable" vector Y = C*X is computed */
    /*     using the solution X obtained. Finally the stated bounds for Y are */
    /*     enlarged to include C*X. To suppress the feature: */


    /*       IOPT(LP)=8 */
    /*         . */
    /*       CALL DBOCLS( ) */

    /*     Use of this option adds 1 to the required length of IOPT(*). */

    /*   9 */
    /*   - */
    /*     This option allows the user to suppress the pretriangularizing */
    /*     step of the least squares matrix that is done within DBOCLS( ). */
    /*     This is primarily a means of enhancing the subprogram efficiency */
    /*     and has little effect on accuracy. To suppress the step, set: */

    /*       IOPT(LP)=9 */
    /*         . */
    /*       CALL DBOCLS( ) */

    /*     Use of this option adds 1 to the required length of IOPT(*). */

    /*   99 */
    /*   -- */
    /*     There are no more options to change. */

    /*     Only option numbers -99, -9,-8,...,-1, 1,2,...,9, and 99 are */
    /*     permitted. Other values are errors. Options -99,-1,...,-9 mean */
    /*     that the respective options 99,1,...,9 are left at their default */
    /*     values. An example is the option to suppress the preprocessing of */
    /*     contraints: */

    /*       IOPT(1)=-8 Option is recognized but not changed */
    /*       IOPT(2)=99 */
    /*       CALL DBOCLS( ) */

    /*    Error Messages for DBOCLS() */
    /*    ----- -------- --- -------- */

    /* WARNING in... */
    /* DBOCLS(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .GE. THE NUMBER */
    /* OF EFFECTIVE ROWS=(I2). */
    /*           IN ABOVE MESSAGE, I1=         1 */
    /*           IN ABOVE MESSAGE, I2=         2 */
    /* ERROR NUMBER =        41 */

    /* WARNING IN... */
    /* DBOCLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE .GE. NCOLS+ */
    /* MCON+1=(I2). */
    /*           IN ABOVE MESSAGE, I1=         2 */
    /*           IN ABOVE MESSAGE, I2=         3 */
    /* ERROR NUMBER =        42 */

    /* WARNING IN... */
    /* DBOCLS(). THE DIMENSIONS OF THE ARRAYS BL(),BU(), AND IND()=(I1) */
    /* MUST BE .GE. NCOLS+MCON=(I2). */
    /*           IN ABOVE MESSAGE, I1=         1 */
    /*           IN ABOVE MESSAGE, I2=         2 */
    /* ERROR NUMBER =        43 */

    /* WARNING IN... */
    /* DBOCLS(). THE DIMENSION OF X()=(I1) MUST BE */
    /* .GE. THE REQD.LENGTH=(I2). */
    /*           IN ABOVE MESSAGE, I1=         1 */
    /*           IN ABOVE MESSAGE, I2=         2 */
    /* ERROR NUMBER =        44 */

    /* WARNING IN... */
    /* DBOCLS(). THE . */
    /* DBOCLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*NCOLS+2*MCON=(I2). */
    /*           IN ABOVE MESSAGE, I1=         1 */
    /*           IN ABOVE MESSAGE, I2=         4 */
    /* ERROR NUMBER =        46 */

    /* WARNING IN... */
    /* DBOCLS(). THE DIMENSION OF IOPT()=(I1) MUST BE .GE. THE REQD. */
    /* LEN.=(I2). */
    /*           IN ABOVE MESSAGE, I1=        16 */
    /*           IN ABOVE MESSAGE, I2=        18 */
    /* ERROR NUMBER =        47 */

    /* WARNING IN... */
    /* DBOCLS(). ISCALE OPTION=(I1) MUST BE 1-3. */
    /*           IN ABOVE MESSAGE, I1=         0 */
    /* ERROR NUMBER =        48 */

    /* WARNING IN... */
    /* DBOCLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PROVIDED COLUMN SCALING */
    /* MUST BE POSITIVE. */
    /*           IN ABOVE MESSAGE, I1=         0 */
    /* ERROR NUMBER =        49 */

    /* WARNING IN... */
    /* DBOCLS(). EACH PROVIDED COL. SCALE FACTOR MUST BE POSITIVE. */
    /*  COMPONENT (I1) NOW = (R1). */
    /*           IN ABOVE MESSAGE, I1=         1 */
    /*           IN ABOVE MESSAGE, R1=    0. */
    /* ERROR NUMBER =        50 */

    /* WARNING IN... */
    /* DBOCLS(). THE OPTION NUMBER=(I1) IS NOT DEFINED. */
    /*           IN ABOVE MESSAGE, I1=      1001 */
    /* ERROR NUMBER =        51 */

    /* WARNING IN... */
    /* DBOCLS(). NO. OF ROWS=(I1) MUST BE .GE. 0 .AND. .LE. MDW-MCON=(I2). */
    /*           IN ABOVE MESSAGE, I1=         2 */
    /*           IN ABOVE MESSAGE, I2=         1 */
    /* ERROR NUMBER =        52 */

    /* WARNING IN... */
    /* DBOCLS(). MDW=(I1) MUST BE POSITIVE. */
    /*           IN ABOVE MESSAGE, I1=         0 */
    /* ERROR NUMBER =        53 */

    /* WARNING IN... */
    /* DBOCLS(). MCON=(I1) MUST BE NONNEGATIVE. */
    /*           IN ABOVE MESSAGE, I1=        -1 */
    /* ERROR NUMBER =        54 */

    /* WARNING IN... */
    /* DBOCLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSITIVE. */
    /*           IN ABOVE MESSAGE, I1=         0 */
    /* ERROR NUMBER =        55 */

    /* WARNING IN... */
    /* DBOCLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4. */
    /*           IN ABOVE MESSAGE, I1=         1 */
    /*           IN ABOVE MESSAGE, I2=         0 */
    /* ERROR NUMBER =        56 */

    /* WARNING IN... */
    /* DBOCLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. BU(J)=(R2). */
    /*           IN ABOVE MESSAGE, I1=         1 */
    /*           IN ABOVE MESSAGE, R1=     .1000000000E+01 */
    /*           IN ABOVE MESSAGE, R2=    0. */
    /* ERROR NUMBER =        57 */
    /*           LINEAR CONSTRAINTS, SNLA REPT. SAND82-1517, AUG., (1982). */
    /* ***REFERENCES  HANSON, R. J. LINEAR LEAST SQUARES WITH BOUNDS AND */
    /*                 LINEAR CONSTRAINTS, SNLA REPT. SAND82-1517, AUG.,1982 */
    /* ***ROUTINES CALLED  D1MACH,DASUM,DBOLS,DCOPY,DDOT,DNRM2,DSCAL,XERRWV */
    /* ***END PROLOGUE  DBOCLS */
    /*     REVISED 850604-0900 */
    /*     REVISED YYMMDD-HHMM */

    /*    PURPOSE */
    /*    ------- */
    /*     THIS IS THE MAIN SUBPROGRAM THAT SOLVES THE LEAST SQUARES */
    /*     PROBLEM CONSISTING OF LINEAR CONSTRAINTS */

    /*              C*X = Y */

    /*     AND LEAST SQUARES EQUATIONS */

    /*              E*X = F */

    /*     IN THIS FORMULATION THE VECTORS X AND Y ARE BOTH UNKNOWNS. */
    /*     FURTHER, X AND Y MAY BOTH HAVE USER-SPECIFIED BOUNDS ON EACH */
    /*     COMPONENT.  THE USER MUST HAVE DIMENSION STATEMENTS OF THE */
    /*     FORM */

    /*     DIMENSION W(MDW,NCOLS+MCON+1), BL(NCOLS+MCON),BU(NCOLS+MCON), */
    /*               X(2*(NCOLS+MCON)+2+NX), RW(6*NCOLS+5*MCON) */

    /*     INTEGER IND(NCOLS+MCON), IOPT(16+NI), IW(2*(NCOLS+MCON)) */

    /*     TO CHANGE THIS SUBPROGRAM FROM SINGLE TO DOUBLE PRECISION BEGIN */
    /*     EDITING AT THE CARD 'C++'. */
    /*     CHANGE THIS SUBPROGRAM TO DBOCLS AND THE STRINGS */
    /*     /SDOT/ TO /DDOT/, /SNRM2/ TO /DNRM2/, /SRELPR/ TO /DRELPR/, */
    /*     /R1MACH/ TO /D1MACH/, /E0/ TO /D0/, /SCOPY/ TO /DCOPY/, */
    /*     /SSCAL/ TO /DSCAL/, /SASUM/ TO /DASUM/, /SBOLS/ TO /DBOLS/, */
    /*     /REAL            / TO /DOUBLE PRECISION/. */
    /* ++ */
    /*     THIS VARIABLE REMAINS TYPED REAL. */
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
    --iw;

    /* Function Body */
    /* ***FIRST EXECUTABLE STATEMENT  DBOCLS */
    nerr = 0;
    *mode = 0;
    level = 1;
    if (igo == 0) {
        /*     DO(CHECK VALIDITY OF INPUT DATA) */
        /*     PROCEDURE(CHECK VALIDITY OF INPUT DATA) */

        /*     SEE THAT MDW IS .GT.0. GROSS CHECK ONLY. */
        if (*mdw <= 0) {
            nerr = 53;
            nchar = 36;
            xerrwv_("DBOCLS(). MDW=(I1) MUST BE POSITIVE.", &nchar, &nerr, &
                    level, &c__1, mdw, &idum, &c__0, &rdum, &rdum, (ftnlen)36)
            ;
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L260;
        }

        /*     SEE THAT NUMBER OF CONSTRAINTS IS NONNEGATIVE. */
        if (*mcon < 0) {
            nerr = 54;
            nchar = 40;
            xerrwv_("DBOCLS(). MCON=(I1) MUST BE NONNEGATIVE.", &nchar, &nerr,
                    &level, &c__1, mcon, &idum, &c__0, &rdum, &rdum, (ftnlen)
                    40);
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L260;
        }

        /*     SEE THAT NUMBER OF UNKNOWNS IS POSITIVE. */
        if (*ncols <= 0) {
            nerr = 55;
            nchar = 59;
            xerrwv_("DBOCLS(). NCOLS=(I1) THE NO. OF VARIABLES MUST BE POSIT"
                    "IVE.", &nchar, &nerr, &level, &c__1, ncols, &idum, &c__0,
                    &rdum, &rdum, (ftnlen)59);
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L260;
        }

        /*     SEE THAT CONSTRAINT INDICATORS ARE ALL WELL-DEFINED. */
        i__1 = *ncols + *mcon;
        for (j = 1; j <= i__1; ++j) {
            if (ind[j] < 1 || ind[j] > 4) {
                nerr = 56;
                nchar = 46;
                xerrwv_("DBOCLS(). FOR J=(I1), IND(J)=(I2) MUST BE 1-4.", &
                        nchar, &nerr, &level, &c__2, &j, &ind[j], &c__0, &
                        rdum, &rdum, (ftnlen)46);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            /* L10: */
        }

        /*     SEE THAT BOUNDS ARE CONSISTENT. */
        i__1 = *ncols + *mcon;
        for (j = 1; j <= i__1; ++j) {
            if (ind[j] == 3) {
                if (bl[j] > bu[j]) {
                    nerr = 57;
                    nchar = 58;
                    rdum2 = (real) bl[j];
                    rdum = (real) bu[j];
                    xerrwv_("DBOCLS(). FOR J=(I1), BOUND BL(J)=(R1) IS .GT. "
                            "BU(J)=(R2).", &nchar, &nerr, &level, &c__1, &j, &
                            idum, &c__2, &rdum2, &rdum, (ftnlen)58);
                    /*     DO(RETURN TO USER PROGRAM UNIT) */
                    goto L260;
                }
            }
            /* L20: */
        }
        /*     END PROCEDURE */
        /*     DO(PROCESS OPTION ARRAY) */
        /*     PROCEDURE(PROCESS OPTION ARRAY) */
        zero = 0.;
        one = 1.;
        drelpr = d1mach_(&lc__4);
        checkl = FALSE_;
        filter = TRUE_;
        /*lenx = (*ncols + *mcon << 1) + 2;*/
        lenx = ((*ncols + *mcon) << 1) + 2;
        iscale = 1;
        igo = 1;
        accum = FALSE_;
        pretri = TRUE_;
        lopt = 0;
        mopt = 0;
        lp = 0;
        lds = 0;
        /*     DO FOREVER */
L30:
        lp += lds;
        ip = iopt[lp + 1];
        jp = abs(ip);

        /*     TEST FOR NO MORE OPTIIONS TO CHANGE. */
        if (ip == 99) {
            if (lopt == 0) {
                lopt = -(lp + 2);
            }
            if (mopt == 0) {
                mopt = -(abs(lopt) + 7);
            }
            if (lopt < 0) {
                lbou = abs(lopt);
            } else {
                lbou = lopt - 15;
            }

            /*     SEND COL. SCALING TO DBOLS(). */
            iopt[lbou] = 4;
            iopt[lbou + 1] = 1;

            /*     PASS AN OPTION ARRAY FOR DBOLSM(). */
            iopt[lbou + 2] = 5;

            /*     LOC. OF OPTION ARRAY FOR DBOLSM( ). */
            iopt[lbou + 3] = 8;

            /*     SKIP TO START OF USER-GIVEN OPTION ARRAY FOR DBOLS(). */
            iopt[lbou + 4] = 6;
            iopt[lbou + 6] = 99;
            if (lopt > 0) {
                iopt[lbou + 5] = lopt - lbou + 1;
            } else {
                iopt[lbou + 4] = -iopt[lbou + 4];
            }
            if (mopt < 0) {
                lboum = abs(mopt);
            } else {
                lboum = mopt - 8;
            }

            /*     CHANGE PRETRIANGULARIZATION FACTOR IN DBOLSM(). */
            iopt[lboum] = 5;
            iopt[lboum + 1] = *ncols + *mcon + 1;

            /*     PASS WEIGHT TO DBOLSM() FOR RANK TEST. */
            iopt[lboum + 2] = 6;
            iopt[lboum + 3] = *ncols + *mcon + 2;
            iopt[lboum + 4] = *mcon;

            /*     SKIP TO USER-GIVEN OPTION ARRAY FOR DBOLSM( ). */
            iopt[lboum + 5] = 1;
            iopt[lboum + 7] = 99;
            if (mopt > 0) {
                iopt[lboum + 6] = mopt - lboum + 1;
            } else {
                iopt[lboum + 5] = -iopt[lboum + 5];
            }
            /*     EXIT FOREVER */
            goto L50;
        } else if (jp == 99) {
            lds = 1;
            /*     CYCLE FOREVER */
            goto L50;
        } else if (jp == 1) {
            if (ip > 0) {

                /*     SET UP DIRECTION FLAG LOCATION, ROW STACKING POINTER */
                /*     LOCATION, AND LOCATION FOR NUMBER OF NEW ROWS. */
                locacc = lp + 2;

                /*                  IOPT(LOCACC-1)=OPTION NUMBER FOR SEQ. ACCUMULATION. */
                /*     CONTENTS..   IOPT(LOCACC  )=USER DIRECTION FLAG, 1 OR 2. */
                /*                  IOPT(LOCACC+1)=ROW STACKING POINTER. */
                /*                  IOPT(LOCACC+2)=NUMBER OF NEW ROWS TO PROCESS. */
                /*     USER ACTION WITH THIS OPTION.. */
                /*      (SET UP OPTION DATA FOR SEQ. ACCUMULATION IN IOPT(*).) */
                /*      (MOVE BLOCK OF EQUATIONS INTO W(*,*)  STARTING AT FIRST */
                /*       ROW OF W(*,*) BELOW THE ROWS FOR THE CONSTRAINT MATRIX C. */
                /*       SET IOPT(LOCACC+2)=NO. OF LEAST SQUARES EQUATIONS IN BLOCK. */
                /*              LOOP */
                /*              CALL DBOCLS() */

                /*                  IF(IOPT(LOCACC) .EQ. 1) THEN */
                /*                      STACK EQUAS. INTO W(*,*), STARTING AT */
                /*                      ROW IOPT(LOCACC+1). */
                /*                       INTO W(*,*). */
                /*                       SET IOPT(LOCACC+2)=NO. OF EQUAS. */
                /*                      IF LAST BLOCK OF EQUAS., SET IOPT(LOCACC)=2. */
                /*                  ELSE IF IOPT(LOCACC) .EQ. 2) THEN */
                /*                      (PROCESS IS OVER. EXIT LOOP.) */
                /*                  ELSE */
                /*                      (ERROR CONDITION. SHOULD NOT HAPPEN.) */
                /*                  END IF */
                /*              END LOOP */
                iopt[locacc + 1] = *mcon + 1;
                accum = TRUE_;
                iopt[locacc] = igo;
            }
            lds = 4;
            /*     CYCLE FOREVER */
            goto L30;
        } else if (jp == 2) {
            if (ip > 0) {

                /*     GET ACTUAL LENGTHS OF ARRAYS FOR CHECKING AGAINST NEEDS. */
                locdim = lp + 2;

                /*     LMDW.GE.MCON+MAX(MOUT,NCOLS), IF MCON.GT.0 .AND FILTER */
                /*     LMDW.GE.MCON+MOUT, OTHERWISE */

                /*     LNDW.GE.NCOLS+MCON+1 */
                /*     LLB .GE.NCOLS+MCON */
                /*     LLX .GE.2*(NCOLS+MCON)+2+EXTRA REQD. IN OPTIONS. */
                /*     LLRW.GE.6*NCOLS+5*MCON */
                /*     LLIW.GE.2*(NCOLS+MCON) */
                /*     LIOP.GE. AMOUNT REQD. FOR OPTION ARRAY. */
                lmdw = iopt[locdim];
                lndw = iopt[locdim + 1];
                llb = iopt[locdim + 2];
                llx = iopt[locdim + 3];
                llrw = iopt[locdim + 4];
                lliw = iopt[locdim + 5];
                liopt = iopt[locdim + 6];
                checkl = TRUE_;
            }
            lds = 8;
            /*     CYCLE FOREVER */
            goto L30;

            /*     OPTION TO MODIFY THE COLUMN SCALING. */
        } else if (jp == 3) {
            if (ip > 0) {
                iscale = iopt[lp + 2];

                /*     SEE THAT ISCALE IS 1 THRU 3. */
                if (iscale < 1 || iscale > 3) {
                    nerr = 48;
                    nchar = 41;
                    xerrwv_("DBOCLS(). ISCALE OPTION=(I1) MUST BE 1-3.", &
                            nchar, &nerr, &level, &c__1, &iscale, &idum, &
                            c__0, &rdum, &rdum, (ftnlen)41);
                    /*     DO(RETURN TO USER PROGRAM UNIT) */
                    goto L260;
                }
            }
            lds = 2;
            /*     CYCLE FOREVER */
            goto L30;

            /*     IN THIS OPTION THE USER HAS PROVIDED SCALING.  THE */
            /*     SCALE FACTORS FOR THE COLUMNS BEGIN IN X(NCOLS+IOPT(LP+2)). */
        } else if (jp == 4) {
            if (ip > 0) {
                iscale = 4;
                if (iopt[lp + 2] <= 0) {
                    nerr = 49;
                    nchar = 86;
                    xerrwv_("DBOCLS(). OFFSET PAST X(NCOLS) (I1) FOR USER-PR"
                            "OVIDED COLUMN SCALING MUST BE POSITIVE.", &nchar,
                            &nerr, &level, &c__1, &iopt[lp + 2], &idum, &c__0,
                            &rdum, &rdum, (ftnlen)86);
                    /*     DO(RETURN TO USER PROGRAM UNIT) */
                    goto L260;
                }
                dcopy_(ncols, &x[*ncols + iopt[lp + 2]], &c__1, &rw[1], &c__1)
                ;
                lenx += *ncols;
                i__1 = *ncols;
                for (j = 1; j <= i__1; ++j) {
                    if (rw[j] <= zero) {
                        nerr = 50;
                        nchar = 84;
                        rdum2 = (real) rw[j];
                        xerrwv_("DBOCLS(). EACH PROVIDED COL. SCALE FACTOR M"
                                "UST BE POSITIVE. COMP. (I1)   NOW = (R1).", &
                                nchar, &nerr, &level, &c__1, &j, &idum, &c__1,
                                &rdum2, &rdum, (ftnlen)84);
                        /*     DO(RETURN TO USER PROGRAM UNIT) */
                        goto L260;
                    }
                    /* L40: */
                }
            }
            lds = 2;
            /*     CYCLE FOREVER */
            goto L30;

            /*     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLS(). */
        } else if (jp == 5) {
            if (ip > 0) {
                lopt = iopt[lp + 2];
            }
            lds = 2;
            /*     CYCLE FOREVER */
            goto L30;

            /*     IN THIS OPTION AN OPTION ARRAY IS PROVIDED TO DBOLSM(). */
        } else if (jp == 6) {
            if (ip > 0) {
                mopt = iopt[lp + 2];
            }
            lds = 2;
            /*     CYCLE FOREVER */
            goto L30;

            /*     THIS OPTION USES THE NEXT LOC OF IOPT(*) AS A */
            /*     POINTER VALUE TO SKIP TO NEXT. */
        } else if (jp == 7) {
            if (ip > 0) {
                lp = iopt[lp + 2] - 1;
                lds = 0;
            } else {
                lds = 2;
            }
            /*     CYCLE FOREVER */
            goto L30;

            /*     THIS OPTION AVOIDS THE CONSTRAINT RESOLVING PHASE FOR */
            /*     THE LINEAR CONSTRAINTS C*X=Y. */
        } else if (jp == 8) {
            filter = !(ip > 0);
            lds = 1;
            /*     CYCLE FOREVER */
            goto L30;

            /*     THIS OPTION SUPPRESSES PRETIRANGULARIZATION OF THE LEAST */
            /*     SQUARES EQATIONS. */
        } else if (jp == 9) {
            pretri = !(ip > 0);
            lds = 1;
            /*     CYCLE FOREVER */
            goto L30;

            /*     NO VALID OPTION NUMBER WAS NOTED. THIS IS AN ERROR CONDITION. */
        } else {
            nerr = 51;
            nchar = 48;
            rdum = (real) idum;
            xerrwv_("DBOCLS(). THE OPTION NUMBER=(I1) IS NOT DEFINED.", &
                    nchar, &nerr, &level, &c__1, &jp, &idum, &c__0, &rdum, &
                    rdum, (ftnlen)48);
            /*     DO(RETURN TO USER PROGRAM UNIT) */
            goto L260;
        }
        /*     END FOREVER */
        /*     END PROCEDURE */
L50:
        if (checkl) {
            /*     DO(CHECK LENGTHS OF ARRAYS) */
            /*     PROCEDURE(CHECK LENGTHS OF ARRAYS) */

            /*     THIS FEATURE ALLOWS THE USER TO MAKE SURE THAT THE */
            /*     ARRAYS ARE LONG ENOUGH FOR THE INTENDED PROBLEM SIZE AND USE. */
            if (filter && ! accum) {
                mdwl = *mcon + max(*mrows,*ncols);
            } else {
                mdwl = *mcon + *ncols + 1;
            }
            if (lmdw < mdwl) {
                nerr = 41;
                nchar = 88;
                xerrwv_("DBOCLS(). THE ROW DIMENSION OF W(,)=(I1) MUST BE .G"
                        "E. THE NUMBER OF EFFECTIVE ROWS=(I2).", &nchar, &nerr,
                        &level, &c__2, &lmdw, &mdwl, &c__0, &rdum, &rdum, (
                            ftnlen)88);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            if (lndw < *ncols + *mcon + 1) {
                nerr = 42;
                nchar = 75;
                i__1 = *ncols + *mcon + 1;
                xerrwv_("DBOCLS(). THE COLUMN DIMENSION OF W(,)=(I1) MUST BE"
                        " .GE. NCOLS+MCON+1=(I2).", &nchar, &nerr, &level, &
                        c__2, &lndw, &i__1, &c__0, &rdum, &rdum, (ftnlen)75);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            if (llb < *ncols + *mcon) {
                nerr = 43;
                nchar = 94;
                i__1 = *ncols + *mcon;
                xerrwv_("DBOCLS(). THE DIMENSIONS OF THE ARRAYS BL(),BU(), A"
                        "ND IND()=(I1) MUST BE .GE. NCOLS+MCON=(I2).", &nchar,
                        &nerr, &level, &c__2, &llb, &i__1, &c__0, &rdum, &
                        rdum, (ftnlen)94);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            if (llx < lenx) {
                nerr = 44;
                nchar = 71;
                xerrwv_("DBOCLS(). THE DIMENSION OF X()=(I1) MUST BE .GE. TH"
                        "E REQD. LENGTH=(I2).", &nchar, &nerr, &level, &c__2, &
                        llx, &lenx, &c__0, &rdum, &rdum, (ftnlen)71);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            if (llrw < *ncols * 6 + *mcon * 5) {
                nerr = 45;
                nchar = 70;
                i__1 = *ncols * 6 + *mcon * 5;
                xerrwv_("DBOCLS(). THE DIMENSION OF RW()=(I1) MUST BE .GE. 6"
                        "*NCOLS+5*MCON=(I2).", &nchar, &nerr, &level, &c__2, &
                        llrw, &i__1, &c__0, &rdum, &rdum, (ftnlen)70);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            if (lliw < (*ncols << 1) + (*mcon << 1)) {
                nerr = 46;
                nchar = 69;
                i__1 = (*ncols << 1) + (*mcon << 1);
                xerrwv_("DBOCLS() THE DIMENSION OF IW()=(I1) MUST BE .GE. 2*"
                        "NCOLS+2*MCON=(I2).", &nchar, &nerr, &level, &c__2, &
                        lliw, &i__1, &c__0, &rdum, &rdum, (ftnlen)69);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            if (liopt < lp + 17) {
                nerr = 47;
                nchar = 72;
                i__1 = lp + 17;
                xerrwv_("DBOCLS(). THE DIMENSION OF IOPT()=(I1) MUST BE .GE."
                        " THE REQD. LEN.=(I2).", &nchar, &nerr, &level, &c__2,
                        &liopt, &i__1, &c__0, &rdum, &rdum, (ftnlen)72);
                /*     DO(RETURN TO USER PROGRAM UNIT) */
                goto L260;
            }
            /*     END PROCEDURE */
        }
    }

    /*     OPTIONALLY GO BACK TO THE USER FOR ACCUMULATION OF LEAST SQUARES */
    /*     EQUATIONS AND DIRECTIONS FOR PROCESSING THESE EQUATIONS. */
    /*     DO(ACCUMULATE LEAST SQUARES EQUATIONS) */
    /*     PROCEDURE(ACCUMULATE LEAST SQUARES EQUATIONS) */
    if (accum) {
        *mrows = iopt[locacc + 1] - 1 - *mcon;
        inrows = iopt[locacc + 2];
        mnew = *mrows + inrows;
        if (mnew < 0 || mnew + *mcon > *mdw) {
            nerr = 52;
            nchar = 66;
            i__1 = *mdw - *mcon;
            xerrwv_("DBOCLS(). NO. OF ROWS=(I1) MUST BE .GE. 0 .AND. .LE.MDW"
                    "-MCON=(I2)", &nchar, &nerr, &level, &c__2, &mnew, &i__1, &
                    c__0, &rdum, &rdum, (ftnlen)65);
            /*    (RETURN TO USER PROGRAM UNIT) */
            goto L260;
        }
    }

    /*     USE THE SOFTWARE OF DBOLS( ) FOR THE TRIANGULARIZATION OF THE */
    /*     LEAST SQUARES MATRIX.  THIS MAY INVOLVE A SYSTALTIC INTERCHANGE */
    /*     OF PROCESSING POINTERS BETWEEN THE CALLING AND CALLED (DBOLS()) */
    /*     PROGRAM UNITS. */
    jopt[0] = 1;
    jopt[1] = 2;
    jopt[3] = *mrows;
    jopt[4] = 99;
    irw = *ncols + 1;
    iiw = 1;
    if (accum || pretri) {
        dbols_(&w[*mcon + 1 + w_dim1], mdw, &mout, ncols, &bl[1], &bu[1], &
               ind[1], jopt, &x[1], rnorm, mode, &rw[irw], &iw[iiw]);
    } else {
        mout = *mrows;
    }
    if (accum) {
        accum = iopt[locacc] == 1;
        iopt[locacc + 1] = jopt[2] + *mcon;
        /* Computing MIN */
        i__1 = *ncols + 1;
        *mrows = min(i__1,mnew);
    }
    /*     END PROCEDURE */
    if (accum) {
        return 0;
    }
    /*     DO(SOLVE CONSTRAINED AND BOUNDED LEAST SQUARES PROBLEM) */
    /*     PROCEDURE(SOLVE CONSTRAINED AND BOUNDED LEAST SQUARES PROBLEM) */

    /*     MOVE RIGHT HAND SIDE OF LEAST SQUARES EQUATIONS. */
    dcopy_(&mout, &w[*mcon + 1 + (*ncols + 1) * w_dim1], &c__1, &w[*mcon + 1
            + (*ncols + *mcon + 1) * w_dim1], &c__1);
    if (*mcon > 0 && filter) {

        /*     PROJECT THE LINEAR CONSTRAINTS INTO A REACHABLE SET. */
        i__1 = *mcon;
        for (i__ = 1; i__ <= i__1; ++i__) {
            dcopy_(ncols, &w[i__ + w_dim1], mdw, &w[*mcon + 1 + (*ncols + i__)
                                                    * w_dim1], &c__1);
            /* L60: */
        }

        /*      PLACE (-)IDENTITY MATRIX AFTER CONSTRAINT DATA. */
        i__1 = *ncols + *mcon + 1;
        for (j = *ncols + 1; j <= i__1; ++j) {
            w[j * w_dim1 + 1] = zero;
            dcopy_(mcon, &w[j * w_dim1 + 1], &c__0, &w[j * w_dim1 + 1], &c__1)
            ;
            /* L70: */
        }
        w[(*ncols + 1) * w_dim1 + 1] = -one;
        i__1 = *mdw + 1;
        dcopy_(mcon, &w[(*ncols + 1) * w_dim1 + 1], &c__0, &w[(*ncols + 1) *
                w_dim1 + 1], &i__1);

        /*     OBTAIN A 'FEASIBLE POINT' FOR THE LINEAR CONSTRAINTS. */
        jopt[0] = 99;
        irw = *ncols + 1;
        iiw = 1;
        i__1 = *ncols + *mcon;
        dbols_(&w[w_offset], mdw, mcon, &i__1, &bl[1], &bu[1], &ind[1], jopt,
               &x[1], rnormc, &modec, &rw[irw], &iw[iiw]);

        /*     ENLARGE THE BOUNDS SET, IF REQUIRED, TO INCLUDE POINTS THAT */
        /*     CAN BE REACHED. */
        i__1 = *ncols + *mcon;
        for (j = *ncols + 1; j <= i__1; ++j) {
            icase = ind[j];
            if (icase < 4) {
                t = ddot_(ncols, &w[*mcon + 1 + j * w_dim1], &c__1, &x[1], &
                          c__1);
            }
            switch (icase) {
            case 1:
                goto L80;
            case 2:
                goto L90;
            case 3:
                goto L100;
            case 4:
                goto L110;
            }
            goto L120;
            /*     CASE 1 */
L80:
            /* Computing MIN */
            d__1 = t, d__2 = bl[j];
            bl[j] = min(d__1,d__2);
            goto L120;
            /*     CASE 2 */
L90:
            /* Computing MAX */
            d__1 = t, d__2 = bu[j];
            bu[j] = max(d__1,d__2);
            goto L120;
            /*     CASE 3 */
L100:
            /* Computing MIN */
            d__1 = t, d__2 = bl[j];
            bl[j] = min(d__1,d__2);
            /* Computing MAX */
            d__1 = t, d__2 = bu[j];
            bu[j] = max(d__1,d__2);
            goto L120;
            /*     CASE 4 */
L110:
L120:
            /* L130: */
            ;
        }

        /*     MOVE CONSTRAINT DATA BACK TO THE ORIGINAL AREA. */
        i__1 = *ncols + *mcon;
        for (j = *ncols + 1; j <= i__1; ++j) {
            dcopy_(ncols, &w[*mcon + 1 + j * w_dim1], &c__1, &w[j - *ncols +
                    w_dim1], mdw);
            /* L140: */
        }
    }
    if (*mcon > 0) {
        i__1 = *ncols + *mcon;
        for (j = *ncols + 1; j <= i__1; ++j) {
            w[*mcon + 1 + j * w_dim1] = zero;
            dcopy_(&mout, &w[*mcon + 1 + j * w_dim1], &c__0, &w[*mcon + 1 + j
                    * w_dim1], &c__1);
            /* L150: */
        }

        /*     PUT IN (-)IDENTITY MATRIX (POSSIBLY) ONCE AGAIN. */
        i__1 = *ncols + *mcon + 1;
        for (j = *ncols + 1; j <= i__1; ++j) {
            w[j * w_dim1 + 1] = zero;
            dcopy_(mcon, &w[j * w_dim1 + 1], &c__0, &w[j * w_dim1 + 1], &c__1)
            ;
            /* L160: */
        }
        w[(*ncols + 1) * w_dim1 + 1] = -one;
        i__1 = *mdw + 1;
        dcopy_(mcon, &w[(*ncols + 1) * w_dim1 + 1], &c__0, &w[(*ncols + 1) *
                w_dim1 + 1], &i__1);
    }

    /*     COMPUTE NOMINAL COLUMN SCALING FOR THE UNWEIGHTED MATRIX. */
    cnorm = zero;
    anorm = zero;
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        t1 = dasum_(mcon, &w[j * w_dim1 + 1], &c__1);
        t2 = dasum_(&mout, &w[*mcon + 1 + w_dim1], &c__1);
        t = t1 + t2;
        if (t == zero) {
            t = one;
        }
        cnorm = max(cnorm,t1);
        anorm = max(anorm,t2);
        x[*ncols + *mcon + j] = one / t;
        /* L170: */
    }
    switch (iscale) {
    case 1:
        goto L180;
    case 2:
        goto L190;
    case 3:
        goto L210;
    case 4:
        goto L220;
    }
    goto L230;
    /*     CASE 1 */
L180:
    goto L230;
    /*     CASE 2 */

    /*     SCALE COLS. (BEFORE WEIGHTING) TO HAVE LENGTH ONE. */
L190:
    i__1 = *ncols;
    for (j = 1; j <= i__1; ++j) {
        i__2 = *mcon + mout;
        t = dnrm2_(&i__2, &w[j * w_dim1 + 1], &c__1);
        if (t == zero) {
            t = one;
        }
        x[*ncols + *mcon + j] = one / t;
        /* L200: */
    }
    goto L230;
    /*     CASE 3 */

    /*     SUPPRESS SCALING (USE UNIT MATRIX). */
L210:
    x[*ncols + *mcon + 1] = one;
    dcopy_(ncols, &x[*ncols + *mcon + 1], &c__0, &x[*ncols + *mcon + 1], &
           c__1);
    goto L230;
    /*     CASE 4 */

    /*     THE USER HAS PROVIDED SCALING. */
L220:
    dcopy_(ncols, &rw[1], &c__1, &x[*ncols + *mcon + 1], &c__1);
L230:
    i__1 = *ncols + *mcon;
    for (j = *ncols + 1; j <= i__1; ++j) {
        x[*ncols + *mcon + j] = one;
        /* L240: */
    }

    /*     WEIGHT THE LEAST SQUARES EQUATIONS. */
    wt = drelpr;
    if (anorm > zero) {
        wt /= anorm;
    }
    if (cnorm > zero) {
        wt *= cnorm;
    }
    i__1 = mout;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dscal_(ncols, &wt, &w[i__ + *mcon + w_dim1], mdw);
        /* L250: */
    }
    dscal_(&mout, &wt, &w[*mcon + 1 + (*mcon + *ncols + 1) * w_dim1], &c__1);
    lrw = 1;
    liw = 1;

    /*     SET THE NEW TRIANGULARIZATION FACTOR. */
    /*x[(*ncols + *mcon << 1) + 1] = zero;*/
    x[((*ncols + *mcon) << 1) + 1] = zero;

    /*     SET THE WEIGHT TO USE IN COMPONENTS .GT. MCON, */
    /*     WHEN MAKING LINEAR INDEPENDENCE TEST. */
    /*x[(*ncols + *mcon << 1) + 2] = one / wt; */
    x[((*ncols + *mcon) << 1) + 2] = one / wt;
    i__1 = mout + *mcon;
    i__2 = *ncols + *mcon;
    dbols_(&w[w_offset], mdw, &i__1, &i__2, &bl[1], &bu[1], &ind[1], &iopt[
               lbou], &x[1], rnorm, mode, &rw[lrw], &iw[liw]);
    *rnorm /= wt;
    /*     END PROCEDURE */
    /*     PROCEDURE(RETURN TO USER PROGRAM UNIT) */
L260:
    if (*mode >= 0) {
        *mode = -nerr;
    }
    igo = 0;
    return 0;
    /*     END PROGRAM */
} /* dbocls_ */

