      SUBROUTINE DDASPK (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL)
C
C***BEGIN PROLOGUE  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  910624   (Added HMAX test at 525 in main driver.)
C***REVISION DATE  920929   (CJ in RES call, RES counter fix.)
C***REVISION DATE  921215   (Warnings on poor iteration performance)
C***REVISION DATE  921216   (NRMAX as optional input)
C***REVISION DATE  930315   (Name change: DDINI to DDINIT)
C***REVISION DATE  940822   (Replaced initial condition calculation)
C***REVISION DATE  941101   (Added linesearch in I.C. calculations)
C***REVISION DATE  941220   (Misc. corrections throughout)
C***REVISION DATE  950125   (Added DINVWT routine)
C***REVISION DATE  950714   (Misc. corrections throughout)
C***REVISION DATE  950802   (Default NRMAX = 5, based on tests.)
C***REVISION DATE  950808   (Optional error test added.)
C***REVISION DATE  950814   (Added I.C. constraints and INFO(14))
C***REVISION DATE  950828   (Various minor corrections.)
C***REVISION DATE  951006   (Corrected WT scaling in DFNRMK.)
C***REVISION DATE  951030   (Corrected history update at end of DDASTP.)
C***REVISION DATE  960129   (Corrected RL bug in DLINSD, DLINSK.)
C***REVISION DATE  960301   (Added NONNEG to SAVE statement.)
C***REVISION DATE  000512   (Removed copyright notices.)
C***REVISION DATE  000622   (Corrected LWM value using NCPHI.)
C***REVISION DATE  000628   (Corrected I.C. stopping tests when index = 0.)
C***REVISION DATE  000628   (Fixed alpha test in I.C. calc., Krylov case.)
C***REVISION DATE  000628   (Improved restart in I.C. calc., Krylov case.)
C***REVISION DATE  000628   (Minor corrections throughout.)
C***REVISION DATE  000711   (Fixed Newton convergence test in DNSD, DNSK.) 
C***REVISION DATE  000712   (Fixed tests on TN - TOUT below 420 and 440.)
C***CATEGORY NO.  I1A2
C***KEYWORDS  DIFFERENTIAL/ALGEBRAIC, BACKWARD DIFFERENTIATION FORMULAS,
C             IMPLICIT DIFFERENTIAL SYSTEMS, KRYLOV ITERATION
C***AUTHORS   Linda R. Petzold, Peter N. Brown, Alan C. Hindmarsh, and
C                  Clement W. Ulrich
C             Center for Computational Sciences & Engineering, L-316
C             Lawrence Livermore National Laboratory
C             P.O. Box 808,
C             Livermore, CA 94551
C***PURPOSE  This code solves a system of differential/algebraic 
C            equations of the form 
C               G(t,y,y') = 0 , 
C            using a combination of Backward Differentiation Formula 
C            (BDF) methods and a choice of two linear system solution 
C            methods: direct (dense or band) or Krylov (iterative).
C            This version is in double precision.
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C *Usage:
C
C      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR(*)
C      DOUBLE PRECISION T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*),
C         RWORK(LRW), RPAR(*)
C      EXTERNAL  RES, JAC, PSOL
C
C      CALL DDASPK (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
C     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC, PSOL)
C
C  Quantities which may be altered by the code are:
C     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL, IDID, RWORK(*), IWORK(*)
C
C
C *Arguments:
C
C  RES:EXT          This is the name of a subroutine which you
C                   provide to define the residual function G(t,y,y')
C                   of the differential/algebraic system.
C
C  NEQ:IN           This is the number of equations in the system.
C
C  T:INOUT          This is the current value of the independent 
C                   variable.
C
C  Y(*):INOUT       This array contains the solution components at T.
C
C  YPRIME(*):INOUT  This array contains the derivatives of the solution
C                   components at T.
C
C  TOUT:IN          This is a point at which a solution is desired.
C
C  INFO(N):IN       This is an integer array used to communicate details
C                   of how the solution is to be carried out, such as
C                   tolerance type, matrix structure, step size and
C                   order limits, and choice of nonlinear system method.
C                   N must be at least 20.
C
C  RTOL,ATOL:INOUT  These quantities represent absolute and relative
C                   error tolerances (on local error) which you provide
C                   to indicate how accurately you wish the solution to
C                   be computed.  You may choose them to be both scalars
C                   or else both arrays of length NEQ.
C
C  IDID:OUT         This integer scalar is an indicator reporting what
C                   the code did.  You must monitor this variable to
C                   decide what action to take next.
C
C  RWORK:WORK       A real work array of length LRW which provides the
C                   code with needed storage space.
C
C  LRW:IN           The length of RWORK.
C
C  IWORK:WORK       An integer work array of length LIW which provides
C                   the code with needed storage space.
C
C  LIW:IN           The length of IWORK.
C
C  RPAR,IPAR:IN     These are real and integer parameter arrays which
C                   you can use for communication between your calling
C                   program and the RES, JAC, and PSOL subroutines.
C
C  JAC:EXT          This is the name of a subroutine which you may
C                   provide (optionally) for calculating Jacobian 
C                   (partial derivative) data involved in solving linear
C                   systems within DDASPK.
C
C  PSOL:EXT         This is the name of a subroutine which you must
C                   provide for solving linear systems if you selected
C                   a Krylov method.  The purpose of PSOL is to solve
C                   linear systems involving a left preconditioner P.
C
C *Overview
C
C  The DDASPK solver uses the backward differentiation formulas of
C  orders one through five to solve a system of the form G(t,y,y') = 0
C  for y = Y and y' = YPRIME.  Values for Y and YPRIME at the initial 
C  time must be given as input.  These values should be consistent, 
C  that is, if T, Y, YPRIME are the given initial values, they should 
C  satisfy G(T,Y,YPRIME) = 0.  However, if consistent values are not
C  known, in many cases you can have DDASPK solve for them -- see INFO(11).
C  (This and other options are described in more detail below.)
C
C  Normally, DDASPK solves the system from T to TOUT.  It is easy to
C  continue the solution to get results at additional TOUT.  This is
C  the interval mode of operation.  Intermediate results can also be
C  obtained easily by specifying INFO(3).
C
C  On each step taken by DDASPK, a sequence of nonlinear algebraic  
C  systems arises.  These are solved by one of two types of
C  methods:
C    * a Newton iteration with a direct method for the linear
C      systems involved (INFO(12) = 0), or
C    * a Newton iteration with a preconditioned Krylov iterative 
C      method for the linear systems involved (INFO(12) = 1).
C
C  The direct method choices are dense and band matrix solvers, 
C  with either a user-supplied or an internal difference quotient 
C  Jacobian matrix, as specified by INFO(5) and INFO(6).
C  In the band case, INFO(6) = 1, you must supply half-bandwidths
C  in IWORK(1) and IWORK(2).
C
C  The Krylov method is the Generalized Minimum Residual (GMRES) 
C  method, in either complete or incomplete form, and with 
C  scaling and preconditioning.  The method is implemented
C  in an algorithm called SPIGMR.  Certain options in the Krylov 
C  method case are specified by INFO(13) and INFO(15).
C
C  If the Krylov method is chosen, you may supply a pair of routines,
C  JAC and PSOL, to apply preconditioning to the linear system.
C  If the system is A*x = b, the matrix is A = dG/dY + CJ*dG/dYPRIME
C  (of order NEQ).  This system can then be preconditioned in the form
C  (P-inverse)*A*x = (P-inverse)*b, with left preconditioner P.
C  (DDASPK does not allow right preconditioning.)
C  Then the Krylov method is applied to this altered, but equivalent,
C  linear system, hopefully with much better performance than without
C  preconditioning.  (In addition, a diagonal scaling matrix based on
C  the tolerances is also introduced into the altered system.)
C
C  The JAC routine evaluates any data needed for solving systems
C  with coefficient matrix P, and PSOL carries out that solution.
C  In any case, in order to improve convergence, you should try to
C  make P approximate the matrix A as much as possible, while keeping
C  the system P*x = b reasonably easy and inexpensive to solve for x,
C  given a vector b.
C
C
C *Description
C
C------INPUT - WHAT TO DO ON THE FIRST CALL TO DDASPK-------------------
C
C
C  The first call of the code is defined to be the start of each new
C  problem.  Read through the descriptions of all the following items,
C  provide sufficient storage space for designated arrays, set
C  appropriate variables for the initialization of the problem, and
C  give information about how you want the problem to be solved.
C
C
C  RES -- Provide a subroutine of the form
C
C             SUBROUTINE RES (T, Y, YPRIME, CJ, DELTA, IRES, RPAR, IPAR)
C
C         to define the system of differential/algebraic
C         equations which is to be solved. For the given values
C         of T, Y and YPRIME, the subroutine should return
C         the residual of the differential/algebraic system
C             DELTA = G(T,Y,YPRIME)
C         DELTA is a vector of length NEQ which is output from RES.
C
C         Subroutine RES must not alter T, Y, YPRIME, or CJ.
C         You must declare the name RES in an EXTERNAL
C         statement in your program that calls DDASPK.
C         You must dimension Y, YPRIME, and DELTA in RES.
C
C         The input argument CJ can be ignored, or used to rescale
C         constraint equations in the system (see Ref. 2, p. 145).
C         Note: In this respect, DDASPK is not downward-compatible
C         with DDASSL, which does not have the RES argument CJ.
C
C         IRES is an integer flag which is always equal to zero
C         on input.  Subroutine RES should alter IRES only if it
C         encounters an illegal value of Y or a stop condition.
C         Set IRES = -1 if an input value is illegal, and DDASPK
C         will try to solve the problem without getting IRES = -1.
C         If IRES = -2, DDASPK will return control to the calling
C         program with IDID = -11.
C
C         RPAR and IPAR are real and integer parameter arrays which
C         you can use for communication between your calling program
C         and subroutine RES. They are not altered by DDASPK. If you
C         do not need RPAR or IPAR, ignore these parameters by treat-
C         ing them as dummy arguments. If you do choose to use them,
C         dimension them in your calling program and in RES as arrays
C         of appropriate length.
C
C  NEQ -- Set it to the number of equations in the system (NEQ .GE. 1).
C
C  T -- Set it to the initial point of the integration. (T must be
C       a variable.)
C
C  Y(*) -- Set this array to the initial values of the NEQ solution
C          components at the initial point.  You must dimension Y of
C          length at least NEQ in your calling program.
C
C  YPRIME(*) -- Set this array to the initial values of the NEQ first
C               derivatives of the solution components at the initial
C               point.  You must dimension YPRIME at least NEQ in your
C               calling program. 
C
C  TOUT - Set it to the first point at which a solution is desired.
C         You cannot take TOUT = T.  Integration either forward in T
C         (TOUT .GT. T) or backward in T (TOUT .LT. T) is permitted.
C
C         The code advances the solution from T to TOUT using step
C         sizes which are automatically selected so as to achieve the
C         desired accuracy.  If you wish, the code will return with the
C         solution and its derivative at intermediate steps (the
C         intermediate-output mode) so that you can monitor them,
C         but you still must provide TOUT in accord with the basic
C         aim of the code.
C
C         The first step taken by the code is a critical one because
C         it must reflect how fast the solution changes near the
C         initial point.  The code automatically selects an initial
C         step size which is practically always suitable for the
C         problem.  By using the fact that the code will not step past
C         TOUT in the first step, you could, if necessary, restrict the
C         length of the initial step.
C
C         For some problems it may not be permissible to integrate
C         past a point TSTOP, because a discontinuity occurs there
C         or the solution or its derivative is not defined beyond
C         TSTOP.  When you have declared a TSTOP point (see INFO(4)
C         and RWORK(1)), you have told the code not to integrate past
C         TSTOP.  In this case any tout beyond TSTOP is invalid input.
C
C  INFO(*) - Use the INFO array to give the code more details about
C            how you want your problem solved.  This array should be
C            dimensioned of length 20, though DDASPK uses only the 
C            first 15 entries.  You must respond to all of the following
C            items, which are arranged as questions.  The simplest use
C            of DDASPK corresponds to setting all entries of INFO to 0.
C
C       INFO(1) - This parameter enables the code to initialize itself.
C              You must set it to indicate the start of every new 
C              problem.
C
C          **** Is this the first call for this problem ...
C                yes - set INFO(1) = 0
C                 no - not applicable here.
C                      See below for continuation calls.  ****
C
C       INFO(2) - How much accuracy you want of your solution
C              is specified by the error tolerances RTOL and ATOL.
C              The simplest use is to take them both to be scalars.
C              To obtain more flexibility, they can both be arrays.
C              The code must be told your choice.
C
C          **** Are both error tolerances RTOL, ATOL scalars ...
C                yes - set INFO(2) = 0
C                      and input scalars for both RTOL and ATOL
C                 no - set INFO(2) = 1
C                      and input arrays for both RTOL and ATOL ****
C
C       INFO(3) - The code integrates from T in the direction of TOUT
C              by steps.  If you wish, it will return the computed
C              solution and derivative at the next intermediate step
C              (the intermediate-output mode) or TOUT, whichever comes
C              first.  This is a good way to proceed if you want to
C              see the behavior of the solution.  If you must have
C              solutions at a great many specific TOUT points, this
C              code will compute them efficiently.
C
C          **** Do you want the solution only at
C               TOUT (and not at the next intermediate step) ...
C                yes - set INFO(3) = 0
C                 no - set INFO(3) = 1 ****
C
C       INFO(4) - To handle solutions at a great many specific
C              values TOUT efficiently, this code may integrate past
C              TOUT and interpolate to obtain the result at TOUT.
C              Sometimes it is not possible to integrate beyond some
C              point TSTOP because the equation changes there or it is
C              not defined past TSTOP.  Then you must tell the code
C              this stop condition.
C
C           **** Can the integration be carried out without any
C                restrictions on the independent variable T ...
C                 yes - set INFO(4) = 0
C                  no - set INFO(4) = 1
C                       and define the stopping point TSTOP by
C                       setting RWORK(1) = TSTOP ****
C
C       INFO(5) - used only when INFO(12) = 0 (direct methods).
C              To solve differential/algebraic systems you may wish
C              to use a matrix of partial derivatives of the
C              system of differential equations.  If you do not
C              provide a subroutine to evaluate it analytically (see
C              description of the item JAC in the call list), it will
C              be approximated by numerical differencing in this code.
C              Although it is less trouble for you to have the code
C              compute partial derivatives by numerical differencing,
C              the solution will be more reliable if you provide the
C              derivatives via JAC.  Usually numerical differencing is
C              more costly than evaluating derivatives in JAC, but
C              sometimes it is not - this depends on your problem.
C
C           **** Do you want the code to evaluate the partial deriv-
C                atives automatically by numerical differences ...
C                 yes - set INFO(5) = 0
C                  no - set INFO(5) = 1
C                       and provide subroutine JAC for evaluating the
C                       matrix of partial derivatives ****
C
C       INFO(6) - used only when INFO(12) = 0 (direct methods).
C              DDASPK will perform much better if the matrix of
C              partial derivatives, dG/dY + CJ*dG/dYPRIME (here CJ is
C              a scalar determined by DDASPK), is banded and the code
C              is told this.  In this case, the storage needed will be
C              greatly reduced, numerical differencing will be performed
C              much cheaper, and a number of important algorithms will
C              execute much faster.  The differential equation is said 
C              to have half-bandwidths ML (lower) and MU (upper) if 
C              equation i involves only unknowns Y(j) with
C                             i-ML .le. j .le. i+MU .
C              For all i=1,2,...,NEQ.  Thus, ML and MU are the widths
C              of the lower and upper parts of the band, respectively,
C              with the main diagonal being excluded.  If you do not
C              indicate that the equation has a banded matrix of partial
C              derivatives the code works with a full matrix of NEQ**2
C              elements (stored in the conventional way).  Computations
C              with banded matrices cost less time and storage than with
C              full matrices if  2*ML+MU .lt. NEQ.  If you tell the
C              code that the matrix of partial derivatives has a banded
C              structure and you want to provide subroutine JAC to
C              compute the partial derivatives, then you must be careful
C              to store the elements of the matrix in the special form
C              indicated in the description of JAC.
C
C          **** Do you want to solve the problem using a full (dense)
C               matrix (and not a special banded structure) ...
C                yes - set INFO(6) = 0
C                 no - set INFO(6) = 1
C                       and provide the lower (ML) and upper (MU)
C                       bandwidths by setting
C                       IWORK(1)=ML
C                       IWORK(2)=MU ****
C
C       INFO(7) - You can specify a maximum (absolute value of)
C              stepsize, so that the code will avoid passing over very
C              large regions.
C
C          ****  Do you want the code to decide on its own the maximum
C                stepsize ...
C                 yes - set INFO(7) = 0
C                  no - set INFO(7) = 1
C                       and define HMAX by setting
C                       RWORK(2) = HMAX ****
C
C       INFO(8) -  Differential/algebraic problems may occasionally
C              suffer from severe scaling difficulties on the first
C              step.  If you know a great deal about the scaling of 
C              your problem, you can help to alleviate this problem 
C              by specifying an initial stepsize H0.
C
C          ****  Do you want the code to define its own initial
C                stepsize ...
C                 yes - set INFO(8) = 0
C                  no - set INFO(8) = 1
C                       and define H0 by setting
C                       RWORK(3) = H0 ****
C
C       INFO(9) -  If storage is a severe problem, you can save some
C              storage by restricting the maximum method order MAXORD.
C              The default value is 5.  For each order decrease below 5,
C              the code requires NEQ fewer locations, but it is likely 
C              to be slower.  In any case, you must have 
C              1 .le. MAXORD .le. 5.
C          ****  Do you want the maximum order to default to 5 ...
C                 yes - set INFO(9) = 0
C                  no - set INFO(9) = 1
C                       and define MAXORD by setting
C                       IWORK(3) = MAXORD ****
C
C       INFO(10) - If you know that certain components of the
C              solutions to your equations are always nonnegative
C              (or nonpositive), it may help to set this
C              parameter.  There are three options that are
C              available:
C              1.  To have constraint checking only in the initial
C                  condition calculation.
C              2.  To enforce nonnegativity in Y during the integration.
C              3.  To enforce both options 1 and 2.
C
C              When selecting option 2 or 3, it is probably best to try the
C              code without using this option first, and only use
C              this option if that does not work very well.
C
C          ****  Do you want the code to solve the problem without
C                invoking any special inequality constraints ...
C                 yes - set INFO(10) = 0
C                  no - set INFO(10) = 1 to have option 1 enforced 
C                  no - set INFO(10) = 2 to have option 2 enforced
C                  no - set INFO(10) = 3 to have option 3 enforced ****
C
C                  If you have specified INFO(10) = 1 or 3, then you
C                  will also need to identify how each component of Y
C                  in the initial condition calculation is constrained.
C                  You must set:
C                  IWORK(40+I) = +1 if Y(I) must be .GE. 0,
C                  IWORK(40+I) = +2 if Y(I) must be .GT. 0,
C                  IWORK(40+I) = -1 if Y(I) must be .LE. 0, while
C                  IWORK(40+I) = -2 if Y(I) must be .LT. 0, while
C                  IWORK(40+I) =  0 if Y(I) is not constrained.
C
C       INFO(11) - DDASPK normally requires the initial T, Y, and
C              YPRIME to be consistent.  That is, you must have
C              G(T,Y,YPRIME) = 0 at the initial T.  If you do not know
C              the initial conditions precisely, in some cases
C              DDASPK may be able to compute it.
C
C              Denoting the differential variables in Y by Y_d
C              and the algebraic variables by Y_a, DDASPK can solve
C              one of two initialization problems:
C              1.  Given Y_d, calculate Y_a and Y'_d, or
C              2.  Given Y', calculate Y.
C              In either case, initial values for the given
C              components are input, and initial guesses for
C              the unknown components must also be provided as input.
C
C          ****  Are the initial T, Y, YPRIME consistent ...
C
C                 yes - set INFO(11) = 0
C                  no - set INFO(11) = 1 to calculate option 1 above,
C                    or set INFO(11) = 2 to calculate option 2 ****
C
C                  If you have specified INFO(11) = 1, then you
C                  will also need to identify  which are the
C                  differential and which are the algebraic
C                  components (algebraic components are components
C                  whose derivatives do not appear explicitly
C                  in the function G(T,Y,YPRIME)).  You must set:
C                  IWORK(LID+I) = +1 if Y(I) is a differential variable
C                  IWORK(LID+I) = -1 if Y(I) is an algebraic variable,
C                  where LID = 40 if INFO(10) = 0 or 2 and LID = 40+NEQ
C                  if INFO(10) = 1 or 3.
C
C       INFO(12) - Except for the addition of the RES argument CJ,
C              DDASPK by default is downward-compatible with DDASSL,
C              which uses only direct (dense or band) methods to solve 
C              the linear systems involved.  You must set INFO(12) to
C              indicate whether you want the direct methods or the
C              Krylov iterative method.
C          ****   Do you want DDASPK to use standard direct methods
C                 (dense or band) or the Krylov (iterative) method ...
C                   direct methods - set INFO(12) = 0.
C                   Krylov method  - set INFO(12) = 1,
C                       and check the settings of INFO(13) and INFO(15).
C
C       INFO(13) - used when INFO(12) = 1 (Krylov methods).  
C              DDASPK uses scalars MAXL, KMP, NRMAX, and EPLI for the
C              iterative solution of linear systems.  INFO(13) allows 
C              you to override the default values of these parameters.  
C              These parameters and their defaults are as follows:
C              MAXL = maximum number of iterations in the SPIGMR 
C                 algorithm (MAXL .le. NEQ).  The default is 
C                 MAXL = MIN(5,NEQ).
C              KMP = number of vectors on which orthogonalization is 
C                 done in the SPIGMR algorithm.  The default is 
C                 KMP = MAXL, which corresponds to complete GMRES 
C                 iteration, as opposed to the incomplete form.  
C              NRMAX = maximum number of restarts of the SPIGMR 
C                 algorithm per nonlinear iteration.  The default is
C                 NRMAX = 5.
C              EPLI = convergence test constant in SPIGMR algorithm.
C                 The default is EPLI = 0.05.
C              Note that the length of RWORK depends on both MAXL 
C              and KMP.  See the definition of LRW below.
C          ****   Are MAXL, KMP, and EPLI to be given their
C                 default values ...
C                  yes - set INFO(13) = 0
C                   no - set INFO(13) = 1,
C                        and set all of the following:
C                        IWORK(24) = MAXL (1 .le. MAXL .le. NEQ)
C                        IWORK(25) = KMP  (1 .le. KMP .le. MAXL)
C                        IWORK(26) = NRMAX  (NRMAX .ge. 0)
C                        RWORK(10) = EPLI (0 .lt. EPLI .lt. 1.0) ****
C
C        INFO(14) - used with INFO(11) > 0 (initial condition 
C               calculation is requested).  In this case, you may
C               request control to be returned to the calling program
C               immediately after the initial condition calculation,
C               before proceeding to the integration of the system
C               (e.g. to examine the computed Y and YPRIME).
C               If this is done, and if the initialization succeeded
C               (IDID = 4), you should reset INFO(11) to 0 for the
C               next call, to prevent the solver from repeating the 
C               initialization (and to avoid an infinite loop). 
C          ****   Do you want to proceed to the integration after
C                 the initial condition calculation is done ...
C                 yes - set INFO(14) = 0
C                  no - set INFO(14) = 1                        ****
C
C        INFO(15) - used when INFO(12) = 1 (Krylov methods).
C               When using preconditioning in the Krylov method,
C               you must supply a subroutine, PSOL, which solves the
C               associated linear systems using P.
C               The usage of DDASPK is simpler if PSOL can carry out
C               the solution without any prior calculation of data.
C               However, if some partial derivative data is to be
C               calculated in advance and used repeatedly in PSOL,
C               then you must supply a JAC routine to do this,
C               and set INFO(15) to indicate that JAC is to be called
C               for this purpose.  For example, P might be an
C               approximation to a part of the matrix A which can be
C               calculated and LU-factored for repeated solutions of
C               the preconditioner system.  The arrays WP and IWP
C               (described under JAC and PSOL) can be used to
C               communicate data between JAC and PSOL.
C          ****   Does PSOL operate with no prior preparation ...
C                 yes - set INFO(15) = 0 (no JAC routine)
C                  no - set INFO(15) = 1
C                       and supply a JAC routine to evaluate and
C                       preprocess any required Jacobian data.  ****
C
C         INFO(16) - option to exclude algebraic variables from
C               the error test.  
C          ****   Do you wish to control errors locally on
C                 all the variables...
C                 yes - set INFO(16) = 0
C                  no - set INFO(16) = 1
C                       If you have specified INFO(16) = 1, then you
C                       will also need to identify  which are the
C                       differential and which are the algebraic
C                       components (algebraic components are components
C                       whose derivatives do not appear explicitly
C                       in the function G(T,Y,YPRIME)).  You must set:
C                       IWORK(LID+I) = +1 if Y(I) is a differential 
C                                      variable, and
C                       IWORK(LID+I) = -1 if Y(I) is an algebraic
C                                      variable,
C                       where LID = 40 if INFO(10) = 0 or 2 and 
C                       LID = 40 + NEQ if INFO(10) = 1 or 3.
C
C       INFO(17) - used when INFO(11) > 0 (DDASPK is to do an 
C              initial condition calculation).
C              DDASPK uses several heuristic control quantities in the
C              initial condition calculation.  They have default values,
C              but can  also be set by the user using INFO(17).
C              These parameters and their defaults are as follows:
C              MXNIT  = maximum number of Newton iterations
C                 per Jacobian or preconditioner evaluation.
C                 The default is:
C                 MXNIT =  5 in the direct case (INFO(12) = 0), and
C                 MXNIT = 15 in the Krylov case (INFO(12) = 1).
C              MXNJ   = maximum number of Jacobian or preconditioner
C                 evaluations.  The default is:
C                 MXNJ = 6 in the direct case (INFO(12) = 0), and
C                 MXNJ = 2 in the Krylov case (INFO(12) = 1).
C              MXNH   = maximum number of values of the artificial
C                 stepsize parameter H to be tried if INFO(11) = 1.
C                 The default is MXNH = 5.
C                 NOTE: the maximum number of Newton iterations
C                 allowed in all is MXNIT*MXNJ*MXNH if INFO(11) = 1,
C                 and MXNIT*MXNJ if INFO(11) = 2.
C              LSOFF  = flag to turn off the linesearch algorithm
C                 (LSOFF = 0 means linesearch is on, LSOFF = 1 means
C                 it is turned off).  The default is LSOFF = 0.
C              STPTOL = minimum scaled step in linesearch algorithm.
C                 The default is STPTOL = (unit roundoff)**(2/3).
C              EPINIT = swing factor in the Newton iteration convergence
C                 test.  The test is applied to the residual vector,
C                 premultiplied by the approximate Jacobian (in the
C                 direct case) or the preconditioner (in the Krylov
C                 case).  For convergence, the weighted RMS norm of
C                 this vector (scaled by the error weights) must be
C                 less than EPINIT*EPCON, where EPCON = .33 is the
C                 analogous test constant used in the time steps.
C                 The default is EPINIT = .01.
C          ****   Are the initial condition heuristic controls to be 
C                 given their default values...
C                  yes - set INFO(17) = 0
C                   no - set INFO(17) = 1,
C                        and set all of the following:
C                        IWORK(32) = MXNIT (.GT. 0)
C                        IWORK(33) = MXNJ (.GT. 0)
C                        IWORK(34) = MXNH (.GT. 0)
C                        IWORK(35) = LSOFF ( = 0 or 1)
C                        RWORK(14) = STPTOL (.GT. 0.0)
C                        RWORK(15) = EPINIT (.GT. 0.0)  ****
C
C         INFO(18) - option to get extra printing in initial condition 
C                calculation.
C          ****   Do you wish to have extra printing...
C                 no  - set INFO(18) = 0
C                 yes - set INFO(18) = 1 for minimal printing, or
C                       set INFO(18) = 2 for full printing.
C                       If you have specified INFO(18) .ge. 1, data
C                       will be printed with the error handler routines.
C                       To print to a non-default unit number L, include
C                       the line  CALL XSETUN(L)  in your program.  ****
C
C   RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL)
C               error tolerances to tell the code how accurately you
C               want the solution to be computed.  They must be defined
C               as variables because the code may change them.
C               you have two choices --
C                     Both RTOL and ATOL are scalars (INFO(2) = 0), or
C                     both RTOL and ATOL are vectors (INFO(2) = 1).
C               In either case all components must be non-negative.
C
C               The tolerances are used by the code in a local error
C               test at each step which requires roughly that
C                        abs(local error in Y(i)) .le. EWT(i) ,
C               where EWT(i) = RTOL*abs(Y(i)) + ATOL is an error weight 
C               quantity, for each vector component.
C               (More specifically, a root-mean-square norm is used to
C               measure the size of vectors, and the error test uses the
C               magnitude of the solution at the beginning of the step.)
C
C               The true (global) error is the difference between the
C               true solution of the initial value problem and the
C               computed approximation.  Practically all present day
C               codes, including this one, control the local error at
C               each step and do not even attempt to control the global
C               error directly.
C
C               Usually, but not always, the true accuracy of
C               the computed Y is comparable to the error tolerances.
C               This code will usually, but not always, deliver a more
C               accurate solution if you reduce the tolerances and
C               integrate again.  By comparing two such solutions you 
C               can get a fairly reliable idea of the true error in the
C               solution at the larger tolerances.
C
C               Setting ATOL = 0. results in a pure relative error test
C               on that component.  Setting RTOL = 0. results in a pure
C               absolute error test on that component.  A mixed test
C               with non-zero RTOL and ATOL corresponds roughly to a
C               relative error test when the solution component is
C               much bigger than ATOL and to an absolute error test
C               when the solution component is smaller than the
C               threshold ATOL.
C
C               The code will not attempt to compute a solution at an
C               accuracy unreasonable for the machine being used.  It
C               will advise you if you ask for too much accuracy and
C               inform you as to the maximum accuracy it believes
C               possible.
C
C  RWORK(*) -- a real work array, which should be dimensioned in your
C               calling program with a length equal to the value of
C               LRW (or greater).
C
C  LRW -- Set it to the declared length of the RWORK array.  The
C               minimum length depends on the options you have selected,
C               given by a base value plus additional storage as described
C               below.
C
C               If INFO(12) = 0 (standard direct method), the base value is
C               base = 50 + max(MAXORD+4,7)*NEQ.
C               The default value is MAXORD = 5 (see INFO(9)).  With the
C               default MAXORD, base = 50 + 9*NEQ.
C               Additional storage must be added to the base value for
C               any or all of the following options:
C                 if INFO(6) = 0 (dense matrix), add NEQ**2
C                 if INFO(6) = 1 (banded matrix), then
C                    if INFO(5) = 0, add (2*ML+MU+1)*NEQ + 2*(NEQ/(ML+MU+1)+1),
C                    if INFO(5) = 1, add (2*ML+MU+1)*NEQ,
C                 if INFO(16) = 1, add NEQ.
C
C              If INFO(12) = 1 (Krylov method), the base value is
C              base = 50 + (MAXORD+5)*NEQ + (MAXL+3+MIN0(1,MAXL-KMP))*NEQ +
C                      + (MAXL+3)*MAXL + 1 + LENWP.
C              See PSOL for description of LENWP.  The default values are:
C              MAXORD = 5 (see INFO(9)), MAXL = min(5,NEQ) and KMP = MAXL 
C              (see INFO(13)).
C              With the default values for MAXORD, MAXL and KMP,
C              base = 91 + 18*NEQ + LENWP.
C              Additional storage must be added to the base value for
C              any or all of the following options:
C                if INFO(16) = 1, add NEQ.
C
C
C  IWORK(*) -- an integer work array, which should be dimensioned in
C              your calling program with a length equal to the value
C              of LIW (or greater).
C
C  LIW -- Set it to the declared length of the IWORK array.  The
C             minimum length depends on the options you have selected,
C             given by a base value plus additional storage as described
C             below.
C
C             If INFO(12) = 0 (standard direct method), the base value is
C             base = 40 + NEQ.
C             IF INFO(10) = 1 or 3, add NEQ to the base value.
C             If INFO(11) = 1 or INFO(16) =1, add NEQ to the base value.
C
C             If INFO(12) = 1 (Krylov method), the base value is
C             base = 40 + LENIWP.
C             See PSOL for description of LENIWP.
C             IF INFO(10) = 1 or 3, add NEQ to the base value.
C             If INFO(11) = 1 or INFO(16) = 1, add NEQ to the base value.
C
C
C  RPAR, IPAR -- These are arrays of double precision and integer type,
C             respectively, which are available for you to use
C             for communication between your program that calls
C             DDASPK and the RES subroutine (and the JAC and PSOL
C             subroutines).  They are not altered by DDASPK.
C             If you do not need RPAR or IPAR, ignore these
C             parameters by treating them as dummy arguments.
C             If you do choose to use them, dimension them in
C             your calling program and in RES (and in JAC and PSOL)
C             as arrays of appropriate length.
C
C  JAC -- This is the name of a routine that you may supply
C         (optionally) that relates to the Jacobian matrix of the
C         nonlinear system that the code must solve at each T step.
C         The role of JAC (and its call sequence) depends on whether
C         a direct (INFO(12) = 0) or Krylov (INFO(12) = 1) method 
C         is selected.
C
C         **** INFO(12) = 0 (direct methods):
C           If you are letting the code generate partial derivatives
C           numerically (INFO(5) = 0), then JAC can be absent
C           (or perhaps a dummy routine to satisfy the loader).
C           Otherwise you must supply a JAC routine to compute
C           the matrix A = dG/dY + CJ*dG/dYPRIME.  It must have
C           the form
C
C           SUBROUTINE JAC (T, Y, YPRIME, PD, CJ, RPAR, IPAR)
C
C           The JAC routine must dimension Y, YPRIME, and PD (and RPAR
C           and IPAR if used).  CJ is a scalar which is input to JAC.
C           For the given values of T, Y, and YPRIME, the JAC routine
C           must evaluate the nonzero elements of the matrix A, and 
C           store these values in the array PD.  The elements of PD are 
C           set to zero before each call to JAC, so that only nonzero
C           elements need to be defined.
C           The way you store the elements into the PD array depends
C           on the structure of the matrix indicated by INFO(6).
C           *** INFO(6) = 0 (full or dense matrix) ***
C               Give PD a first dimension of NEQ.  When you evaluate the
C               nonzero partial derivatives of equation i (i.e. of G(i))
C               with respect to component j (of Y and YPRIME), you must
C               store the element in PD according to
C                  PD(i,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
C           *** INFO(6) = 1 (banded matrix with half-bandwidths ML, MU
C                            as described under INFO(6)) ***
C               Give PD a first dimension of 2*ML+MU+1.  When you 
C               evaluate the nonzero partial derivatives of equation i 
C               (i.e. of G(i)) with respect to component j (of Y and 
C               YPRIME), you must store the element in PD according to 
C                  IROW = i - j + ML + MU + 1
C                  PD(IROW,j) = dG(i)/dY(j) + CJ*dG(i)/dYPRIME(j).
C
C          **** INFO(12) = 1 (Krylov method):
C            If you are not calculating Jacobian data in advance for use
C            in PSOL (INFO(15) = 0), JAC can be absent (or perhaps a
C            dummy routine to satisfy the loader).  Otherwise, you may
C            supply a JAC routine to compute and preprocess any parts of
C            of the Jacobian matrix  A = dG/dY + CJ*dG/dYPRIME that are
C            involved in the preconditioner matrix P.
C            It is to have the form
C
C            SUBROUTINE JAC (RES, IRES, NEQ, T, Y, YPRIME, REWT, SAVR,
C                            WK, H, CJ, WP, IWP, IER, RPAR, IPAR)
C
C           The JAC routine must dimension Y, YPRIME, REWT, SAVR, WK,
C           and (if used) WP, IWP, RPAR, and IPAR.
C           The Y, YPRIME, and SAVR arrays contain the current values
C           of Y, YPRIME, and the residual G, respectively.  
C           The array WK is work space of length NEQ.  
C           H is the step size.  CJ is a scalar, input to JAC, that is
C           normally proportional to 1/H.  REWT is an array of 
C           reciprocal error weights, 1/EWT(i), where EWT(i) is
C           RTOL*abs(Y(i)) + ATOL (unless you supplied routine DDAWTS
C           instead), for use in JAC if needed.  For example, if JAC
C           computes difference quotient approximations to partial
C           derivatives, the REWT array may be useful in setting the
C           increments used.  The JAC routine should do any
C           factorization operations called for, in preparation for
C           solving linear systems in PSOL.  The matrix P should
C           be an approximation to the Jacobian,
C           A = dG/dY + CJ*dG/dYPRIME.
C
C           WP and IWP are real and integer work arrays which you may
C           use for communication between your JAC routine and your
C           PSOL routine.  These may be used to store elements of the 
C           preconditioner P, or related matrix data (such as factored
C           forms).  They are not altered by DDASPK.
C           If you do not need WP or IWP, ignore these parameters by
C           treating them as dummy arguments.  If you do use them,
C           dimension them appropriately in your JAC and PSOL routines.
C           See the PSOL description for instructions on setting 
C           the lengths of WP and IWP.
C
C           On return, JAC should set the error flag IER as follows..
C             IER = 0    if JAC was successful,
C             IER .ne. 0 if JAC was unsuccessful (e.g. if Y or YPRIME
C                        was illegal, or a singular matrix is found).
C           (If IER .ne. 0, a smaller stepsize will be tried.)
C           IER = 0 on entry to JAC, so need be reset only on a failure.
C           If RES is used within JAC, then a nonzero value of IRES will
C           override any nonzero value of IER (see the RES description).
C
C         Regardless of the method type, subroutine JAC must not
C         alter T, Y(*), YPRIME(*), H, CJ, or REWT(*).
C         You must declare the name JAC in an EXTERNAL statement in
C         your program that calls DDASPK.
C
C PSOL --  This is the name of a routine you must supply if you have
C         selected a Krylov method (INFO(12) = 1) with preconditioning.
C         In the direct case (INFO(12) = 0), PSOL can be absent 
C         (a dummy routine may have to be supplied to satisfy the 
C         loader).  Otherwise, you must provide a PSOL routine to 
C         solve linear systems arising from preconditioning.
C         When supplied with INFO(12) = 1, the PSOL routine is to 
C         have the form
C
C         SUBROUTINE PSOL (NEQ, T, Y, YPRIME, SAVR, WK, CJ, WGHT,
C                          WP, IWP, B, EPLIN, IER, RPAR, IPAR)
C
C         The PSOL routine must solve linear systems of the form 
C         P*x = b where P is the left preconditioner matrix.
C
C         The right-hand side vector b is in the B array on input, and
C         PSOL must return the solution vector x in B.
C         The Y, YPRIME, and SAVR arrays contain the current values
C         of Y, YPRIME, and the residual G, respectively.  
C
C         Work space required by JAC and/or PSOL, and space for data to
C         be communicated from JAC to PSOL is made available in the form
C         of arrays WP and IWP, which are parts of the RWORK and IWORK
C         arrays, respectively.  The lengths of these real and integer
C         work spaces WP and IWP must be supplied in LENWP and LENIWP,
C         respectively, as follows..
C           IWORK(27) = LENWP = length of real work space WP
C           IWORK(28) = LENIWP = length of integer work space IWP.
C
C         WK is a work array of length NEQ for use by PSOL.
C         CJ is a scalar, input to PSOL, that is normally proportional
C         to 1/H (H = stepsize).  If the old value of CJ
C         (at the time of the last JAC call) is needed, it must have
C         been saved by JAC in WP.
C
C         WGHT is an array of weights, to be used if PSOL uses an
C         iterative method and performs a convergence test.  (In terms
C         of the argument REWT to JAC, WGHT is REWT/sqrt(NEQ).)
C         If PSOL uses an iterative method, it should use EPLIN
C         (a heuristic parameter) as the bound on the weighted norm of
C         the residual for the computed solution.  Specifically, the
C         residual vector R should satisfy
C              SQRT (SUM ( (R(i)*WGHT(i))**2 ) ) .le. EPLIN
C
C         PSOL must not alter NEQ, T, Y, YPRIME, SAVR, CJ, WGHT, EPLIN.
C
C         On return, PSOL should set the error flag IER as follows..
C           IER = 0 if PSOL was successful,
C           IER .lt. 0 if an unrecoverable error occurred, meaning
C                 control will be passed to the calling routine,
C           IER .gt. 0 if a recoverable error occurred, meaning that
C                 the step will be retried with the same step size
C                 but with a call to JAC to update necessary data,
C                 unless the Jacobian data is current, in which case
C                 the step will be retried with a smaller step size.
C           IER = 0 on entry to PSOL so need be reset only on a failure.
C
C         You must declare the name PSOL in an EXTERNAL statement in
C         your program that calls DDASPK.
C
C
C  OPTIONALLY REPLACEABLE SUBROUTINE:
C
C  DDASPK uses a weighted root-mean-square norm to measure the 
C  size of various error vectors.  The weights used in this norm
C  are set in the following subroutine:
C
C    SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, EWT, RPAR, IPAR)
C    DIMENSION RTOL(*), ATOL(*), Y(*), EWT(*), RPAR(*), IPAR(*)
C
C  A DDAWTS routine has been included with DDASPK which sets the
C  weights according to
C    EWT(I) = RTOL*ABS(Y(I)) + ATOL
C  in the case of scalar tolerances (IWT = 0) or
C    EWT(I) = RTOL(I)*ABS(Y(I)) + ATOL(I)
C  in the case of array tolerances (IWT = 1).  (IWT is INFO(2).)
C  In some special cases, it may be appropriate for you to define
C  your own error weights by writing a subroutine DDAWTS to be 
C  called instead of the version supplied.  However, this should 
C  be attempted only after careful thought and consideration. 
C  If you supply this routine, you may use the tolerances and Y 
C  as appropriate, but do not overwrite these variables.  You
C  may also use RPAR and IPAR to communicate data as appropriate.
C  ***Note: Aside from the values of the weights, the choice of 
C  norm used in DDASPK (weighted root-mean-square) is not subject
C  to replacement by the user.  In this respect, DDASPK is not
C  downward-compatible with the original DDASSL solver (in which
C  the norm routine was optionally user-replaceable).
C
C
C------OUTPUT - AFTER ANY RETURN FROM DDASPK----------------------------
C
C  The principal aim of the code is to return a computed solution at
C  T = TOUT, although it is also possible to obtain intermediate
C  results along the way.  To find out whether the code achieved its
C  goal or if the integration process was interrupted before the task
C  was completed, you must check the IDID parameter.
C
C
C   T -- The output value of T is the point to which the solution
C        was successfully advanced.
C
C   Y(*) -- contains the computed solution approximation at T.
C
C   YPRIME(*) -- contains the computed derivative approximation at T.
C
C   IDID -- reports what the code did, described as follows:
C
C                     *** TASK COMPLETED ***
C                Reported by positive values of IDID
C
C           IDID = 1 -- a step was successfully taken in the
C                   intermediate-output mode.  The code has not
C                   yet reached TOUT.
C
C           IDID = 2 -- the integration to TSTOP was successfully
C                   completed (T = TSTOP) by stepping exactly to TSTOP.
C
C           IDID = 3 -- the integration to TOUT was successfully
C                   completed (T = TOUT) by stepping past TOUT.
C                   Y(*) and YPRIME(*) are obtained by interpolation.
C
C           IDID = 4 -- the initial condition calculation, with
C                   INFO(11) > 0, was successful, and INFO(14) = 1.
C                   No integration steps were taken, and the solution
C                   is not considered to have been started.
C
C                    *** TASK INTERRUPTED ***
C                Reported by negative values of IDID
C
C           IDID = -1 -- a large amount of work has been expended
C                     (about 500 steps).
C
C           IDID = -2 -- the error tolerances are too stringent.
C
C           IDID = -3 -- the local error test cannot be satisfied
C                     because you specified a zero component in ATOL
C                     and the corresponding computed solution component
C                     is zero.  Thus, a pure relative error test is
C                     impossible for this component.
C
C           IDID = -5 -- there were repeated failures in the evaluation
C                     or processing of the preconditioner (in JAC).
C
C           IDID = -6 -- DDASPK had repeated error test failures on the
C                     last attempted step.
C
C           IDID = -7 -- the nonlinear system solver in the time integration
C                     could not converge.
C
C           IDID = -8 -- the matrix of partial derivatives appears
C                     to be singular (direct method).
C
C           IDID = -9 -- the nonlinear system solver in the time integration
C                     failed to achieve convergence, and there were repeated 
C                     error test failures in this step.
C
C           IDID =-10 -- the nonlinear system solver in the time integration 
C                     failed to achieve convergence because IRES was equal 
C                     to -1.
C
C           IDID =-11 -- IRES = -2 was encountered and control is
C                     being returned to the calling program.
C
C           IDID =-12 -- DDASPK failed to compute the initial Y, YPRIME.
C
C           IDID =-13 -- unrecoverable error encountered inside user's
C                     PSOL routine, and control is being returned to
C                     the calling program.
C
C           IDID =-14 -- the Krylov linear system solver could not 
C                     achieve convergence.
C
C           IDID =-15,..,-32 -- Not applicable for this code.
C
C                    *** TASK TERMINATED ***
C                reported by the value of IDID=-33
C
C           IDID = -33 -- the code has encountered trouble from which
C                   it cannot recover.  A message is printed
C                   explaining the trouble and control is returned
C                   to the calling program.  For example, this occurs
C                   when invalid input is detected.
C
C   RTOL, ATOL -- these quantities remain unchanged except when
C               IDID = -2.  In this case, the error tolerances have been
C               increased by the code to values which are estimated to
C               be appropriate for continuing the integration.  However,
C               the reported solution at T was obtained using the input
C               values of RTOL and ATOL.
C
C   RWORK, IWORK -- contain information which is usually of no interest
C               to the user but necessary for subsequent calls. 
C               However, you may be interested in the performance data
C               listed below.  These quantities are accessed in RWORK 
C               and IWORK but have internal mnemonic names, as follows..
C
C               RWORK(3)--contains H, the step size h to be attempted
C                        on the next step.
C
C               RWORK(4)--contains TN, the current value of the
C                        independent variable, i.e. the farthest point
C                        integration has reached.  This will differ 
C                        from T if interpolation has been performed 
C                        (IDID = 3).
C
C               RWORK(7)--contains HOLD, the stepsize used on the last
C                        successful step.  If INFO(11) = INFO(14) = 1,
C                        this contains the value of H used in the
C                        initial condition calculation.
C
C               IWORK(7)--contains K, the order of the method to be 
C                        attempted on the next step.
C
C               IWORK(8)--contains KOLD, the order of the method used
C                        on the last step.
C
C               IWORK(11)--contains NST, the number of steps (in T) 
C                        taken so far.
C
C               IWORK(12)--contains NRE, the number of calls to RES 
C                        so far.
C
C               IWORK(13)--contains NJE, the number of calls to JAC so
C                        far (Jacobian or preconditioner evaluations).
C
C               IWORK(14)--contains NETF, the total number of error test
C                        failures so far.
C
C               IWORK(15)--contains NCFN, the total number of nonlinear
C                        convergence failures so far (includes counts
C                        of singular iteration matrix or singular
C                        preconditioners).
C
C               IWORK(16)--contains NCFL, the number of convergence
C                        failures of the linear iteration so far.
C
C               IWORK(17)--contains LENIW, the length of IWORK actually
C                        required.  This is defined on normal returns 
C                        and on an illegal input return for
C                        insufficient storage.
C
C               IWORK(18)--contains LENRW, the length of RWORK actually
C                        required.  This is defined on normal returns 
C                        and on an illegal input return for
C                        insufficient storage.
C
C               IWORK(19)--contains NNI, the total number of nonlinear
C                        iterations so far (each of which calls a
C                        linear solver).
C
C               IWORK(20)--contains NLI, the total number of linear
C                        (Krylov) iterations so far.
C
C               IWORK(21)--contains NPS, the number of PSOL calls so
C                        far, for preconditioning solve operations or
C                        for solutions with the user-supplied method.
C
C               Note: The various counters in IWORK do not include 
C               counts during a call made with INFO(11) > 0 and
C               INFO(14) = 1.
C
C
C------INPUT - WHAT TO DO TO CONTINUE THE INTEGRATION  -----------------
C              (CALLS AFTER THE FIRST)
C
C     This code is organized so that subsequent calls to continue the
C     integration involve little (if any) additional effort on your
C     part.  You must monitor the IDID parameter in order to determine
C     what to do next.
C
C     Recalling that the principal task of the code is to integrate
C     from T to TOUT (the interval mode), usually all you will need
C     to do is specify a new TOUT upon reaching the current TOUT.
C
C     Do not alter any quantity not specifically permitted below.  In
C     particular do not alter NEQ, T, Y(*), YPRIME(*), RWORK(*), 
C     IWORK(*), or the differential equation in subroutine RES.  Any 
C     such alteration constitutes a new problem and must be treated 
C     as such, i.e. you must start afresh.
C
C     You cannot change from array to scalar error control or vice
C     versa (INFO(2)), but you can change the size of the entries of
C     RTOL or ATOL.  Increasing a tolerance makes the equation easier
C     to integrate.  Decreasing a tolerance will make the equation
C     harder to integrate and should generally be avoided.
C
C     You can switch from the intermediate-output mode to the
C     interval mode (INFO(3)) or vice versa at any time.
C
C     If it has been necessary to prevent the integration from going
C     past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C     code will not integrate to any TOUT beyond the currently
C     specified TSTOP.  Once TSTOP has been reached, you must change
C     the value of TSTOP or set INFO(4) = 0.  You may change INFO(4)
C     or TSTOP at any time but you must supply the value of TSTOP in
C     RWORK(1) whenever you set INFO(4) = 1.
C
C     Do not change INFO(5), INFO(6), INFO(12-17) or their associated
C     IWORK/RWORK locations unless you are going to restart the code.
C
C                    *** FOLLOWING A COMPLETED TASK ***
C
C     If..
C     IDID = 1, call the code again to continue the integration
C                  another step in the direction of TOUT.
C
C     IDID = 2 or 3, define a new TOUT and call the code again.
C                  TOUT must be different from T.  You cannot change
C                  the direction of integration without restarting.
C
C     IDID = 4, reset INFO(11) = 0 and call the code again to begin
C                  the integration.  (If you leave INFO(11) > 0 and
C                  INFO(14) = 1, you may generate an infinite loop.)
C                  In this situation, the next call to DASPK is 
C                  considered to be the first call for the problem,
C                  in that all initializations are done.
C
C                    *** FOLLOWING AN INTERRUPTED TASK ***
C
C     To show the code that you realize the task was interrupted and
C     that you want to continue, you must take appropriate action and
C     set INFO(1) = 1.
C
C     If..
C     IDID = -1, the code has taken about 500 steps.  If you want to
C                  continue, set INFO(1) = 1 and call the code again.
C                  An additional 500 steps will be allowed.
C
C
C     IDID = -2, the error tolerances RTOL, ATOL have been increased
C                  to values the code estimates appropriate for
C                  continuing.  You may want to change them yourself.
C                  If you are sure you want to continue with relaxed
C                  error tolerances, set INFO(1) = 1 and call the code
C                  again.
C
C     IDID = -3, a solution component is zero and you set the
C                  corresponding component of ATOL to zero.  If you
C                  are sure you want to continue, you must first alter
C                  the error criterion to use positive values of ATOL 
C                  for those components corresponding to zero solution
C                  components, then set INFO(1) = 1 and call the code
C                  again.
C
C     IDID = -4  --- cannot occur with this code.
C
C     IDID = -5, your JAC routine failed with the Krylov method.  Check
C                  for errors in JAC and restart the integration.
C
C     IDID = -6, repeated error test failures occurred on the last
C                  attempted step in DDASPK.  A singularity in the
C                  solution may be present.  If you are absolutely
C                  certain you want to continue, you should restart
C                  the integration.  (Provide initial values of Y and
C                  YPRIME which are consistent.)
C
C     IDID = -7, repeated convergence test failures occurred on the last
C                  attempted step in DDASPK.  An inaccurate or ill-
C                  conditioned Jacobian or preconditioner may be the
C                  problem.  If you are absolutely certain you want
C                  to continue, you should restart the integration.
C
C
C     IDID = -8, the matrix of partial derivatives is singular, with
C                  the use of direct methods.  Some of your equations
C                  may be redundant.  DDASPK cannot solve the problem
C                  as stated.  It is possible that the redundant
C                  equations could be removed, and then DDASPK could
C                  solve the problem.  It is also possible that a
C                  solution to your problem either does not exist
C                  or is not unique.
C
C     IDID = -9, DDASPK had multiple convergence test failures, preceded
C                  by multiple error test failures, on the last
C                  attempted step.  It is possible that your problem is
C                  ill-posed and cannot be solved using this code.  Or,
C                  there may be a discontinuity or a singularity in the
C                  solution.  If you are absolutely certain you want to
C                  continue, you should restart the integration.
C
C     IDID = -10, DDASPK had multiple convergence test failures
C                  because IRES was equal to -1.  If you are
C                  absolutely certain you want to continue, you
C                  should restart the integration.
C
C     IDID = -11, there was an unrecoverable error (IRES = -2) from RES
C                  inside the nonlinear system solver.  Determine the
C                  cause before trying again.
C
C     IDID = -12, DDASPK failed to compute the initial Y and YPRIME
C                  vectors.  This could happen because the initial 
C                  approximation to Y or YPRIME was not very good, or
C                  because no consistent values of these vectors exist.
C                  The problem could also be caused by an inaccurate or
C                  singular iteration matrix, or a poor preconditioner.
C
C     IDID = -13, there was an unrecoverable error encountered inside 
C                  your PSOL routine.  Determine the cause before 
C                  trying again.
C
C     IDID = -14, the Krylov linear system solver failed to achieve
C                  convergence.  This may be due to ill-conditioning
C                  in the iteration matrix, or a singularity in the
C                  preconditioner (if one is being used).
C                  Another possibility is that there is a better
C                  choice of Krylov parameters (see INFO(13)).
C                  Possibly the failure is caused by redundant equations
C                  in the system, or by inconsistent equations.
C                  In that case, reformulate the system to make it
C                  consistent and non-redundant.
C
C     IDID = -15,..,-32 --- Cannot occur with this code.
C
C                       *** FOLLOWING A TERMINATED TASK ***
C
C     If IDID = -33, you cannot continue the solution of this problem.
C                  An attempt to do so will result in your run being
C                  terminated.
C
C  ---------------------------------------------------------------------
C
C***REFERENCES
C  1.  L. R. Petzold, A Description of DASSL: A Differential/Algebraic
C      System Solver, in Scientific Computing, R. S. Stepleman et al.
C      (Eds.), North-Holland, Amsterdam, 1983, pp. 65-68.
C  2.  K. E. Brenan, S. L. Campbell, and L. R. Petzold, Numerical 
C      Solution of Initial-Value Problems in Differential-Algebraic
C      Equations, Elsevier, New York, 1989.
C  3.  P. N. Brown and A. C. Hindmarsh, Reduced Storage Matrix Methods
C      in Stiff ODE Systems, J. Applied Mathematics and Computation,
C      31 (1989), pp. 40-91.
C  4.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Using Krylov
C      Methods in the Solution of Large-Scale Differential-Algebraic
C      Systems, SIAM J. Sci. Comp., 15 (1994), pp. 1467-1488.
C  5.  P. N. Brown, A. C. Hindmarsh, and L. R. Petzold, Consistent
C      Initial Condition Calculation for Differential-Algebraic
C      Systems, SIAM J. Sci. Comp. 19 (1998), pp. 1495-1512.
C
C***ROUTINES CALLED
C
C   The following are all the subordinate routines used by DDASPK.
C
C   DDASIC computes consistent initial conditions.
C   DYYPNW updates Y and YPRIME in linesearch for initial condition
C          calculation.
C   DDSTP  carries out one step of the integration.
C   DCNSTR/DCNST0 check the current solution for constraint violations.
C   DDAWTS sets error weight quantities.
C   DINVWT tests and inverts the error weights.
C   DDATRP performs interpolation to get an output solution.
C   DDWNRM computes the weighted root-mean-square norm of a vector.
C   D1MACH provides the unit roundoff of the computer.
C   XERRWD/XSETF/XSETUN/IXSAV is a package to handle error messages. 
C   DDASID nonlinear equation driver to initialize Y and YPRIME using
C          direct linear system solver methods.  Interfaces to Newton
C          solver (direct case).
C   DNSID  solves the nonlinear system for unknown initial values by
C          modified Newton iteration and direct linear system methods.
C   DLINSD carries out linesearch algorithm for initial condition
C          calculation (direct case).
C   DFNRMD calculates weighted norm of preconditioned residual in
C          initial condition calculation (direct case).
C   DNEDD  nonlinear equation driver for direct linear system solver
C          methods.  Interfaces to Newton solver (direct case).
C   DMATD  assembles the iteration matrix (direct case).
C   DNSD   solves the associated nonlinear system by modified
C          Newton iteration and direct linear system methods.
C   DSLVD  interfaces to linear system solver (direct case).
C   DDASIK nonlinear equation driver to initialize Y and YPRIME using
C          Krylov iterative linear system methods.  Interfaces to
C          Newton solver (Krylov case).
C   DNSIK  solves the nonlinear system for unknown initial values by
C          Newton iteration and Krylov iterative linear system methods.
C   DLINSK carries out linesearch algorithm for initial condition
C          calculation (Krylov case).
C   DFNRMK calculates weighted norm of preconditioned residual in
C          initial condition calculation (Krylov case).
C   DNEDK  nonlinear equation driver for iterative linear system solver
C          methods.  Interfaces to Newton solver (Krylov case).
C   DNSK   solves the associated nonlinear system by Inexact Newton
C          iteration and (linear) Krylov iteration.
C   DSLVK  interfaces to linear system solver (Krylov case).
C   DSPIGM solves a linear system by SPIGMR algorithm.
C   DATV   computes matrix-vector product in Krylov algorithm.
C   DORTH  performs orthogonalization of Krylov basis vectors.
C   DHEQR  performs QR factorization of Hessenberg matrix.
C   DHELS  finds least-squares solution of Hessenberg linear system.
C   DGEFA, DGESL, DGBFA, DGBSL are LINPACK routines for solving 
C          linear systems (dense or band direct methods).
C   DAXPY, DCOPY, DDOT, DNRM2, DSCAL are Basic Linear Algebra (BLAS)
C          routines.
C
C The routines called directly by DDASPK are:
C   DCNST0, DDAWTS, DINVWT, D1MACH, DDWNRM, DDASIC, DDATRP, DDSTP,
C   XERRWD
C
C***END PROLOGUE DDASPK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      LOGICAL DONE, LAVL, LCFN, LCFL, LWARN
      DIMENSION Y(*),YPRIME(*)
      DIMENSION INFO(20)
      DIMENSION RWORK(LRW),IWORK(LIW)
      DIMENSION RTOL(*),ATOL(*)
      DIMENSION RPAR(*),IPAR(*)
      CHARACTER MSG*80
      EXTERNAL  RES, JAC, PSOL, DDASID, DDASIK, DNEDD, DNEDK
C
C     Set pointers into IWORK.
C
      PARAMETER (LML=1, LMU=2, LMTYPE=4, 
     *   LIWM=1, LMXORD=3, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
     *   LNS=9, LNSTL=10, LNST=11, LNRE=12, LNJE=13, LETF=14, LNCFN=15,
     *   LNCFL=16, LNIW=17, LNRW=18, LNNI=19, LNLI=20, LNPS=21,
     *   LNPD=22, LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26, LLNWP=27,
     *   LLNIWP=28, LLOCWP=29, LLCIWP=30, LKPRIN=31,
     *   LMXNIT=32, LMXNJ=33, LMXNH=34, LLSOFF=35, LICNS=41)
C
C     Set pointers into RWORK.
C
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4, LCJ=5, LCJOLD=6,
     *   LHOLD=7, LS=8, LROUND=9, LEPLI=10, LSQRN=11, LRSQRN=12,
     *   LEPCON=13, LSTOL=14, LEPIN=15,
     *   LALPHA=21, LBETA=27, LGAMMA=33, LPSI=39, LSIGMA=45, LDELTA=51)
C
      SAVE LID, LENID, NONNEG, NCPHI
C
C
C***FIRST EXECUTABLE STATEMENT  DDASPK
C
C
      IF(INFO(1).NE.0) GO TO 100
C
C-----------------------------------------------------------------------
C     This block is executed for the initial call only.
C     It contains checking of inputs and initializations.
C-----------------------------------------------------------------------
C
C     First check INFO array to make sure all elements of INFO
C     Are within the proper range.  (INFO(1) is checked later, because
C     it must be tested on every call.) ITEMP holds the location
C     within INFO which may be out of range.
C
      DO 10 I=2,9
         ITEMP = I
         IF (INFO(I) .NE. 0 .AND. INFO(I) .NE. 1) GO TO 701
 10      CONTINUE
      ITEMP = 10
      IF(INFO(10).LT.0 .OR. INFO(10).GT.3) GO TO 701
      ITEMP = 11
      IF(INFO(11).LT.0 .OR. INFO(11).GT.2) GO TO 701
      DO 15 I=12,17
         ITEMP = I
         IF (INFO(I) .NE. 0 .AND. INFO(I) .NE. 1) GO TO 701
 15      CONTINUE
      ITEMP = 18
      IF(INFO(18).LT.0 .OR. INFO(18).GT.2) GO TO 701

C
C     Check NEQ to see if it is positive.
C
      IF (NEQ .LE. 0) GO TO 702
C
C     Check and compute maximum order.
C
      MXORD=5
      IF (INFO(9) .NE. 0) THEN
         MXORD=IWORK(LMXORD)
         IF (MXORD .LT. 1 .OR. MXORD .GT. 5) GO TO 703
         ENDIF
      IWORK(LMXORD)=MXORD
C
C     Set and/or check inputs for constraint checking (INFO(10) .NE. 0).
C     Set values for ICNFLG, NONNEG, and pointer LID.
C
      ICNFLG = 0
      NONNEG = 0
      LID = LICNS
      IF (INFO(10) .EQ. 0) GO TO 20
      IF (INFO(10) .EQ. 1) THEN
         ICNFLG = 1
         NONNEG = 0
         LID = LICNS + NEQ
      ELSEIF (INFO(10) .EQ. 2) THEN
         ICNFLG = 0
         NONNEG = 1
      ELSE
         ICNFLG = 1
         NONNEG = 1
         LID = LICNS + NEQ
      ENDIF
C
 20   CONTINUE
C
C     Set and/or check inputs for Krylov solver (INFO(12) .NE. 0).
C     If indicated, set default values for MAXL, KMP, NRMAX, and EPLI.
C     Otherwise, verify inputs required for iterative solver.
C
      IF (INFO(12) .EQ. 0) GO TO 25
C
      IWORK(LMITER) = INFO(12)
      IF (INFO(13) .EQ. 0) THEN
         IWORK(LMAXL) = MIN(5,NEQ)
         IWORK(LKMP) = IWORK(LMAXL)
         IWORK(LNRMAX) = 5
         RWORK(LEPLI) = 0.05D0
      ELSE
         IF(IWORK(LMAXL) .LT. 1 .OR. IWORK(LMAXL) .GT. NEQ) GO TO 720
         IF(IWORK(LKMP) .LT. 1 .OR. IWORK(LKMP) .GT. IWORK(LMAXL))
     1      GO TO 721
         IF(IWORK(LNRMAX) .LT. 0) GO TO 722
         IF(RWORK(LEPLI).LE.0.0D0 .OR. RWORK(LEPLI).GE.1.0D0)GO TO 723
         ENDIF
C
 25   CONTINUE
C
C     Set and/or check controls for the initial condition calculation
C     (INFO(11) .GT. 0).  If indicated, set default values.
C     Otherwise, verify inputs required for iterative solver.
C
      IF (INFO(11) .EQ. 0) GO TO 30
      IF (INFO(17) .EQ. 0) THEN
        IWORK(LMXNIT) = 5
        IF (INFO(12) .GT. 0) IWORK(LMXNIT) = 15
        IWORK(LMXNJ) = 6
        IF (INFO(12) .GT. 0) IWORK(LMXNJ) = 2
        IWORK(LMXNH) = 5
        IWORK(LLSOFF) = 0
        RWORK(LEPIN) = 0.01D0
      ELSE
        IF (IWORK(LMXNIT) .LE. 0) GO TO 725
        IF (IWORK(LMXNJ) .LE. 0) GO TO 725
        IF (IWORK(LMXNH) .LE. 0) GO TO 725
        LSOFF = IWORK(LLSOFF)
        IF (LSOFF .LT. 0 .OR. LSOFF .GT. 1) GO TO 725
        IF (RWORK(LEPIN) .LE. 0.0D0) GO TO 725
        ENDIF
C
 30   CONTINUE
C
C     Below is the computation and checking of the work array lengths
C     LENIW and LENRW, using direct methods (INFO(12) = 0) or
C     the Krylov methods (INFO(12) = 1).
C
      LENIC = 0
      IF (INFO(10) .EQ. 1 .OR. INFO(10) .EQ. 3) LENIC = NEQ
      LENID = 0
      IF (INFO(11) .EQ. 1 .OR. INFO(16) .EQ. 1) LENID = NEQ
      IF (INFO(12) .EQ. 0) THEN
C
C        Compute MTYPE, etc.  Check ML and MU.
C
         NCPHI = MAX(MXORD + 1, 4)
         IF(INFO(6).EQ.0) THEN 
            LENPD = NEQ**2
            LENRW = 50 + (NCPHI+3)*NEQ + LENPD
            IF(INFO(5).EQ.0) THEN
               IWORK(LMTYPE)=2
            ELSE
               IWORK(LMTYPE)=1
            ENDIF
         ELSE
            IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
            IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
            LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
            IF(INFO(5).EQ.0) THEN
               IWORK(LMTYPE)=5
               MBAND=IWORK(LML)+IWORK(LMU)+1
               MSAVE=(NEQ/MBAND)+1
               LENRW = 50 + (NCPHI+3)*NEQ + LENPD + 2*MSAVE
            ELSE
               IWORK(LMTYPE)=4
               LENRW = 50 + (NCPHI+3)*NEQ + LENPD
            ENDIF
         ENDIF
C
C        Compute LENIW, LENWP, LENIWP.
C
         LENIW = 40 + LENIC + LENID + NEQ
         LENWP = 0
         LENIWP = 0
C
      ELSE IF (INFO(12) .EQ. 1)  THEN
         NCPHI = MXORD + 1
         MAXL = IWORK(LMAXL)
         LENWP = IWORK(LLNWP)
         LENIWP = IWORK(LLNIWP)
         LENPD = (MAXL+3+MIN0(1,MAXL-IWORK(LKMP)))*NEQ
     1         + (MAXL+3)*MAXL + 1 + LENWP
         LENRW = 50 + (MXORD+5)*NEQ + LENPD
         LENIW = 40 + LENIC + LENID + LENIWP
C
      ENDIF
      IF(INFO(16) .NE. 0) LENRW = LENRW + NEQ
C
C     Check lengths of RWORK and IWORK.
C
      IWORK(LNIW)=LENIW
      IWORK(LNRW)=LENRW
      IWORK(LNPD)=LENPD
      IWORK(LLOCWP) = LENPD-LENWP+1
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
C
C     Check ICNSTR for legality.
C
      IF (LENIC .GT. 0) THEN
        DO 40 I = 1,NEQ
          ICI = IWORK(LICNS-1+I)
          IF (ICI .LT. -2 .OR. ICI .GT. 2) GO TO 726
 40       CONTINUE
        ENDIF
C
C     Check Y for consistency with constraints.
C
      IF (LENIC .GT. 0) THEN
        CALL DCNST0(NEQ,Y,IWORK(LICNS),IRET)
        IF (IRET .NE. 0) GO TO 727
        ENDIF
C
C     Check ID for legality and set INDEX = 0 or 1.
C
      INDEX = 1
      IF (LENID .GT. 0) THEN
        INDEX = 0
        DO 50 I = 1,NEQ
          IDI = IWORK(LID-1+I)
          IF (IDI .NE. 1 .AND. IDI .NE. -1) GO TO 724
          IF (IDI .EQ. -1) INDEX = 1
 50       CONTINUE
        ENDIF
C
C     Check to see that TOUT is different from T.
C
      IF(TOUT .EQ. T)GO TO 719
C
C     Check HMAX.
C
      IF(INFO(7) .NE. 0) THEN
         HMAX = RWORK(LHMAX)
         IF (HMAX .LE. 0.0D0) GO TO 710
         ENDIF
C
C     Initialize counters and other flags.
C
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
      IWORK(LETF)=0
      IWORK(LNCFN)=0
      IWORK(LNNI)=0
      IWORK(LNLI)=0
      IWORK(LNPS)=0
      IWORK(LNCFL)=0
      IWORK(LKPRIN)=INFO(18)
      IDID=1
      GO TO 200
C
C-----------------------------------------------------------------------
C     This block is for continuation calls only.
C     Here we check INFO(1), and if the last step was interrupted,
C     we check whether appropriate action was taken.
C-----------------------------------------------------------------------
C
100   CONTINUE
      IF(INFO(1).EQ.1)GO TO 110
      ITEMP = 1
      IF(INFO(1).NE.-1)GO TO 701
C
C     If we are here, the last step was interrupted by an error
C     condition from DDSTP, and appropriate action was not taken.
C     This is a fatal error.
C
      MSG = 'DASPK--  THE LAST STEP TERMINATED WITH A NEGATIVE'
      CALL XERRWD(MSG,49,201,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  VALUE (=I1) OF IDID AND NO APPROPRIATE'
      CALL XERRWD(MSG,47,202,0,1,IDID,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  ACTION WAS TAKEN. RUN TERMINATED'
      CALL XERRWD(MSG,41,203,1,0,0,0,0,0.0D0,0.0D0)
      RETURN
110   CONTINUE
C
C-----------------------------------------------------------------------
C     This block is executed on all calls.
C
C     Counters are saved for later checks of performance.
C     Then the error tolerance parameters are checked, and the
C     work array pointers are set.
C-----------------------------------------------------------------------
C
200   CONTINUE
C
C     Save counters for use later.
C
      IWORK(LNSTL)=IWORK(LNST)
      NLI0 = IWORK(LNLI)
      NNI0 = IWORK(LNNI)
      NCFN0 = IWORK(LNCFN)
      NCFL0 = IWORK(LNCFL)
      NWARN = 0
C
C     Check RTOL and ATOL.
C
      NZFLG = 0
      RTOLI = RTOL(1)
      ATOLI = ATOL(1)
      DO 210 I=1,NEQ
         IF (INFO(2) .EQ. 1) RTOLI = RTOL(I)
         IF (INFO(2) .EQ. 1) ATOLI = ATOL(I)
         IF (RTOLI .GT. 0.0D0 .OR. ATOLI .GT. 0.0D0) NZFLG = 1
         IF (RTOLI .LT. 0.0D0) GO TO 706
         IF (ATOLI .LT. 0.0D0) GO TO 707
210      CONTINUE
      IF (NZFLG .EQ. 0) GO TO 708
C
C     Set pointers to RWORK and IWORK segments.
C     For direct methods, SAVR is not used.
C
      IWORK(LLCIWP) = LID + LENID
      LSAVR = LDELTA
      IF (INFO(12) .NE. 0) LSAVR = LDELTA + NEQ
      LE = LSAVR + NEQ
      LWT = LE + NEQ
      LVT = LWT
      IF (INFO(16) .NE. 0) LVT = LWT + NEQ
      LPHI = LVT + NEQ
      LWM = LPHI + NCPHI*NEQ
      IF (INFO(1) .EQ. 1) GO TO 400
C
C-----------------------------------------------------------------------
C     This block is executed on the initial call only.
C     Set the initial step size, the error weight vector, and PHI.
C     Compute unknown initial components of Y and YPRIME, if requested.
C-----------------------------------------------------------------------
C
300   CONTINUE
      TN=T
      IDID=1
C
C     Set error weight array WT and altered weight array VT.
C
      CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      CALL DINVWT(NEQ,RWORK(LWT),IER)
      IF (IER .NE. 0) GO TO 713
      IF (INFO(16) .NE. 0) THEN
        DO 305 I = 1, NEQ
 305      RWORK(LVT+I-1) = MAX(IWORK(LID+I-1),0)*RWORK(LWT+I-1)
        ENDIF
C
C     Compute unit roundoff and HMIN.
C
      UROUND = D1MACH(4)
      RWORK(LROUND) = UROUND
      HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
C
C     Set/check STPTOL control for initial condition calculation.
C     
      IF (INFO(11) .NE. 0) THEN
        IF( INFO(17) .EQ. 0) THEN
          RWORK(LSTOL) = UROUND**.6667D0
        ELSE
          IF (RWORK(LSTOL) .LE. 0.0D0) GO TO 725
          ENDIF
        ENDIF
C
C     Compute EPCON and square root of NEQ and its reciprocal, used
C     inside iterative solver.
C
      RWORK(LEPCON) = 0.33D0
      FLOATN = NEQ
      RWORK(LSQRN) = SQRT(FLOATN)
      RWORK(LRSQRN) = 1.D0/RWORK(LSQRN)
C
C     Check initial interval to see that it is long enough.
C
      TDIST = ABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
C
C     Check H0, if this was input.
C
      IF (INFO(8) .EQ. 0) GO TO 310
         H0 = RWORK(LH)
         IF ((TOUT - T)*H0 .LT. 0.0D0) GO TO 711
         IF (H0 .EQ. 0.0D0) GO TO 712
         GO TO 320
310    CONTINUE
C
C     Compute initial stepsize, to be used by either
C     DDSTP or DDASIC, depending on INFO(11).
C
      H0 = 0.001D0*TDIST
      YPNORM = DDWNRM(NEQ,YPRIME,RWORK(LVT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5D0/H0) H0 = 0.5D0/YPNORM
      H0 = SIGN(H0,TOUT-T)
C
C     Adjust H0 if necessary to meet HMAX bound.
C
320   IF (INFO(7) .EQ. 0) GO TO 330
         RH = ABS(H0)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H0 = H0/RH
C
C     Check against TSTOP, if applicable.
C
330   IF (INFO(4) .EQ. 0) GO TO 340
         TSTOP = RWORK(LTSTOP)
         write(*,*) 'tstop = ',tstop
         IF ((TSTOP - T)*H0 .LT. 0.0D0) GO TO 715
         IF ((T + H0 - TSTOP)*H0 .GT. 0.0D0) H0 = TSTOP - T
         IF ((TSTOP - TOUT)*H0 .LT. 0.0D0) GO TO 709
C
340   IF (INFO(11) .EQ. 0) GO TO 370
C
C     Compute unknown components of initial Y and YPRIME, depending
C     on INFO(11) and INFO(12).  INFO(12) represents the nonlinear
C     solver type (direct/Krylov).  Pass the name of the specific 
C     nonlinear solver, depending on INFO(12).  The location of the work
C     arrays SAVR, YIC, YPIC, PWK also differ in the two cases.
C     For use in stopping tests, pass TSCALE = TDIST if INDEX = 0.
C
      NWT = 1
      EPCONI = RWORK(LEPIN)*RWORK(LEPCON)
      TSCALE = 0.0D0
      IF (INDEX .EQ. 0) TSCALE = TDIST
350   IF (INFO(12) .EQ. 0) THEN
         LYIC = LPHI + 2*NEQ
         LYPIC = LYIC + NEQ
         LPWK = LYPIC
         CALL DDASIC(TN,Y,YPRIME,NEQ,INFO(11),IWORK(LID),
     *     RES,JAC,PSOL,H0,TSCALE,RWORK(LWT),NWT,IDID,RPAR,IPAR,
     *     RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *     RWORK(LYIC),RWORK(LYPIC),RWORK(LPWK),RWORK(LWM),IWORK(LIWM),
     *     RWORK(LROUND),RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *     EPCONI,RWORK(LSTOL),INFO(15),ICNFLG,IWORK(LICNS),DDASID)
      ELSE IF (INFO(12) .EQ. 1) THEN
         LYIC = LWM
         LYPIC = LYIC + NEQ
         LPWK = LYPIC + NEQ
         CALL DDASIC(TN,Y,YPRIME,NEQ,INFO(11),IWORK(LID),
     *     RES,JAC,PSOL,H0,TSCALE,RWORK(LWT),NWT,IDID,RPAR,IPAR,
     *     RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *     RWORK(LYIC),RWORK(LYPIC),RWORK(LPWK),RWORK(LWM),IWORK(LIWM),
     *     RWORK(LROUND),RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *     EPCONI,RWORK(LSTOL),INFO(15),ICNFLG,IWORK(LICNS),DDASIK)
      ENDIF
C
      IF (IDID .LT. 0) GO TO 600
C
C     DDASIC was successful.  If this was the first call to DDASIC,
C     update the WT array (with the current Y) and call it again.
C
      IF (NWT .EQ. 2) GO TO 355
      NWT = 2
      CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      CALL DINVWT(NEQ,RWORK(LWT),IER)
      IF (IER .NE. 0) GO TO 713
      GO TO 350
C
C     If INFO(14) = 1, return now with IDID = 4.
C
355   IF (INFO(14) .EQ. 1) THEN
        IDID = 4
        H = H0
        IF (INFO(11) .EQ. 1) RWORK(LHOLD) = H0
        GO TO 590
      ENDIF
C
C     Update the WT and VT arrays one more time, with the new Y.
C
      CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      CALL DINVWT(NEQ,RWORK(LWT),IER)
      IF (IER .NE. 0) GO TO 713
      IF (INFO(16) .NE. 0) THEN
        DO 357 I = 1, NEQ
 357      RWORK(LVT+I-1) = MAX(IWORK(LID+I-1),0)*RWORK(LWT+I-1)
        ENDIF
C
C     Reset the initial stepsize to be used by DDSTP.
C     Use H0, if this was input.  Otherwise, recompute H0,
C     and adjust it if necessary to meet HMAX bound.
C
      IF (INFO(8) .NE. 0) THEN
         H0 = RWORK(LH)
         GO TO 360
         ENDIF
C
      H0 = 0.001D0*TDIST
      YPNORM = DDWNRM(NEQ,YPRIME,RWORK(LVT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5D0/H0) H0 = 0.5D0/YPNORM
      H0 = SIGN(H0,TOUT-T)
C
360   IF (INFO(7) .NE. 0) THEN
         RH = ABS(H0)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H0 = H0/RH
         ENDIF
C
C     Check against TSTOP, if applicable.
C
      IF (INFO(4) .NE. 0) THEN
         TSTOP = RWORK(LTSTOP)
         write(*,*) 'tstop = ',tstop
         IF ((T + H0 - TSTOP)*H0 .GT. 0.0D0) H0 = TSTOP - T
         ENDIF
C
C     Load H and RWORK(LH) with H0.
C
370   H = H0
      RWORK(LH) = H
C
C     Load Y and H*YPRIME into PHI(*,1) and PHI(*,2).
C
      ITEMP = LPHI + NEQ
      DO 380 I = 1,NEQ
         RWORK(LPHI + I - 1) = Y(I)
380      RWORK(ITEMP + I - 1) = H*YPRIME(I)
C
      GO TO 500
C
C-----------------------------------------------------------------------
C     This block is for continuation calls only.
C     Its purpose is to check stop conditions before taking a step.
C     Adjust H if necessary to meet HMAX bound.
C-----------------------------------------------------------------------
C
400   CONTINUE
      UROUND=RWORK(LROUND)
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
         RH = ABS(H)/RWORK(LHMAX)
         IF(RH .GT. 1.0D0) H = H/RH
410   CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0D0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 490
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
420   IF((TN-T)*H .LE. 0.0D0) GO TO 490
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 425
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
425   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
430   IF(INFO(3) .EQ. 1) GO TO 440
      TSTOP=RWORK(LTSTOP)
         write(*,*) 'tstop = ',tstop
      IF((TN-TSTOP)*H.GT.0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 450
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
440   TSTOP = RWORK(LTSTOP)
         write(*,*) 'tstop = ',tstop
      IF((TN-TSTOP)*H .GT. 0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709
      IF((TN-T)*H .LE. 0.0D0) GO TO 450
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 445
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
445   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
450   CONTINUE
C
C     Check whether we are within roundoff of TSTOP.
C
      IF(ABS(TN-TSTOP).GT.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 460
      CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
460   TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490
      H=TSTOP-TN
      RWORK(LH)=H
C
490   IF (DONE) GO TO 590
C
C-----------------------------------------------------------------------
C     The next block contains the call to the one-step integrator DDSTP.
C     This is a looping point for the integration steps.
C     Check for too many steps.
C     Check for poor Newton/Krylov performance.
C     Update WT.  Check for too much accuracy requested.
C     Compute minimum stepsize.
C-----------------------------------------------------------------------
C
500   CONTINUE
C
C     Check for too many steps.
C
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500) GO TO 505
           IDID=-1
           GO TO 527
C
C Check for poor Newton/Krylov performance.
C
505   IF (INFO(12) .EQ. 0) GO TO 510
      NSTD = IWORK(LNST) - IWORK(LNSTL)
      NNID = IWORK(LNNI) - NNI0
      IF (NSTD .LT. 10 .OR. NNID .EQ. 0) GO TO 510
      AVLIN = REAL(IWORK(LNLI) - NLI0)/REAL(NNID)
      RCFN = REAL(IWORK(LNCFN) - NCFN0)/REAL(NSTD)
      RCFL = REAL(IWORK(LNCFL) - NCFL0)/REAL(NNID)
      FMAXL = IWORK(LMAXL)
      LAVL = AVLIN .GT. FMAXL
      LCFN = RCFN .GT. 0.9D0
      LCFL = RCFL .GT. 0.9D0
      LWARN = LAVL .OR. LCFN .OR. LCFL
      IF (.NOT.LWARN) GO TO 510
      NWARN = NWARN + 1
      IF (NWARN .GT. 10) GO TO 510
      IF (LAVL) THEN
        MSG = 'DASPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 501, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
        MSG = '      at T = R1. Average no. of linear iterations = R2  '
        CALL XERRWD (MSG, 56, 501, 0, 0, 0, 0, 2, TN, AVLIN)
        ENDIF
      IF (LCFN) THEN
        MSG = 'DASPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 502, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
        MSG = '      at T = R1. Nonlinear convergence failure rate = R2'
        CALL XERRWD (MSG, 56, 502, 0, 0, 0, 0, 2, TN, RCFN)
        ENDIF
      IF (LCFL) THEN
        MSG = 'DASPK-- Warning. Poor iterative algorithm performance   '
        CALL XERRWD (MSG, 56, 503, 0, 0, 0, 0, 0, 0.0D0, 0.0D0)
        MSG = '      at T = R1. Linear convergence failure rate = R2   '
        CALL XERRWD (MSG, 56, 503, 0, 0, 0, 0, 2, TN, RCFL)
        ENDIF
C
C     Update WT and VT, if this is not the first call.
C
510   CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),RWORK(LWT),
     *            RPAR,IPAR)
      CALL DINVWT(NEQ,RWORK(LWT),IER)
      IF (IER .NE. 0) THEN
        IDID = -3
        GO TO 527
        ENDIF
      IF (INFO(16) .NE. 0) THEN
        DO 515 I = 1, NEQ
 515      RWORK(LVT+I-1) = MAX(IWORK(LID+I-1),0)*RWORK(LWT+I-1)
        ENDIF
C
C     Test for too much accuracy requested.
C
      R = DDWNRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*100.0D0*UROUND
      IF (R .LE. 1.0D0) GO TO 525
C
C     Multiply RTOL and ATOL by R and return.
C
      IF(INFO(2).EQ.1)GO TO 523
           RTOL(1)=R*RTOL(1)
           ATOL(1)=R*ATOL(1)
           IDID=-2
           GO TO 527
523   DO 524 I=1,NEQ
           RTOL(I)=R*RTOL(I)
524        ATOL(I)=R*ATOL(I)
      IDID=-2
      GO TO 527
525   CONTINUE
C
C     Compute minimum stepsize.
C
      HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
C
C     Test H vs. HMAX
      IF (INFO(7) .NE. 0) THEN
         RH = ABS(H)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H = H/RH
         ENDIF
C
C     Call the one-step integrator.
C     Note that INFO(12) represents the nonlinear solver type.
C     Pass the required nonlinear solver, depending upon INFO(12).
C
      IF (INFO(12) .EQ. 0) THEN
         CALL DDSTP(TN,Y,YPRIME,NEQ,
     *      RES,JAC,PSOL,H,RWORK(LWT),RWORK(LVT),INFO(1),IDID,RPAR,IPAR,
     *      RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *      RWORK(LWM),IWORK(LIWM),
     *      RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *      RWORK(LPSI),RWORK(LSIGMA),
     *      RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),RWORK(LS),HMIN,
     *      RWORK(LROUND), RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *      RWORK(LEPCON), IWORK(LPHASE),IWORK(LJCALC),INFO(15),
     *      IWORK(LK), IWORK(LKOLD),IWORK(LNS),NONNEG,INFO(12),
     *      DNEDD)
      ELSE IF (INFO(12) .EQ. 1) THEN
         CALL DDSTP(TN,Y,YPRIME,NEQ,
     *      RES,JAC,PSOL,H,RWORK(LWT),RWORK(LVT),INFO(1),IDID,RPAR,IPAR,
     *      RWORK(LPHI),RWORK(LSAVR),RWORK(LDELTA),RWORK(LE),
     *      RWORK(LWM),IWORK(LIWM),
     *      RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *      RWORK(LPSI),RWORK(LSIGMA),
     *      RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),RWORK(LS),HMIN,
     *      RWORK(LROUND), RWORK(LEPLI),RWORK(LSQRN),RWORK(LRSQRN),
     *      RWORK(LEPCON), IWORK(LPHASE),IWORK(LJCALC),INFO(15),
     *      IWORK(LK), IWORK(LKOLD),IWORK(LNS),NONNEG,INFO(12),
     *      DNEDK)
      ENDIF
C
527   IF(IDID.LT.0)GO TO 600
C
C-----------------------------------------------------------------------
C     This block handles the case of a successful return from DDSTP
C     (IDID=1).  Test for stop conditions.
C-----------------------------------------------------------------------
C
      IF(INFO(4).NE.0)GO TO 540
           IF(INFO(3).NE.0)GO TO 530
             IF((TN-TOUT)*H.LT.0.0D0)GO TO 500
             CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
530          IF((TN-TOUT)*H.GE.0.0D0)GO TO 535
             T=TN
             IDID=1
             GO TO 580
535          CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
540   IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 542
         CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         T=TOUT
         IDID=3
         GO TO 580
542   IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 545
      TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
      H=TSTOP-TN
      GO TO 500
545   CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
550   IF((TN-TOUT)*H.GE.0.0D0)GO TO 555
      IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
552   CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
555   CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
580   CONTINUE
C
C-----------------------------------------------------------------------
C     All successful returns from DDASPK are made from this block.
C-----------------------------------------------------------------------
C
590   CONTINUE
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     This block handles all unsuccessful returns other than for
C     illegal input.
C-----------------------------------------------------------------------
C
600   CONTINUE
      ITEMP = -IDID
      GO TO (610,620,630,700,655,640,650,660,670,675,
     *  680,685,690,695), ITEMP
C
C     The maximum number of steps was taken before
C     reaching tout.
C
610   MSG = 'DASPK--  AT CURRENT T (=R1)  500 STEPS'
      CALL XERRWD(MSG,38,610,0,0,0,0,1,TN,0.0D0)
      MSG = 'DASPK--  TAKEN ON THIS CALL BEFORE REACHING TOUT'
      CALL XERRWD(MSG,48,611,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Too much accuracy for machine precision.
C
620   MSG = 'DASPK--  AT T (=R1) TOO MUCH ACCURACY REQUESTED'
      CALL XERRWD(MSG,47,620,0,0,0,0,1,TN,0.0D0)
      MSG = 'DASPK--  FOR PRECISION OF MACHINE. RTOL AND ATOL'
      CALL XERRWD(MSG,48,621,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  WERE INCREASED TO APPROPRIATE VALUES'
      CALL XERRWD(MSG,45,622,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     WT(I) .LE. 0.0D0 for some I (not at start of problem).
C
630   MSG = 'DASPK--  AT T (=R1) SOME ELEMENT OF WT'
      CALL XERRWD(MSG,38,630,0,0,0,0,1,TN,0.0D0)
      MSG = 'DASPK--  HAS BECOME .LE. 0.0'
      CALL XERRWD(MSG,28,631,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Error test failed repeatedly or with H=HMIN.
C
640   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,640,0,0,0,0,2,TN,H)
      MSG='DASPK--  ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN'
      CALL XERRWD(MSG,57,641,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Nonlinear solver failed to converge repeatedly or with H=HMIN.
C
650   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,650,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  NONLINEAR SOLVER FAILED TO CONVERGE'
      CALL XERRWD(MSG,44,651,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  REPEATEDLY OR WITH ABS(H)=HMIN'
      CALL XERRWD(MSG,40,652,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     The preconditioner had repeated failures.
C
655   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,655,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  PRECONDITIONER HAD REPEATED FAILURES.'
      CALL XERRWD(MSG,46,656,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     The iteration matrix is singular.
C
660   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,660,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  ITERATION MATRIX IS SINGULAR.'
      CALL XERRWD(MSG,38,661,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Nonlinear system failure preceded by error test failures.
C
670   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,670,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  NONLINEAR SOLVER COULD NOT CONVERGE.'
      CALL XERRWD(MSG,45,671,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  ALSO, THE ERROR TEST FAILED REPEATEDLY.'
      CALL XERRWD(MSG,49,672,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Nonlinear system failure because IRES = -1.
C
675   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,675,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  NONLINEAR SYSTEM SOLVER COULD NOT CONVERGE'
      CALL XERRWD(MSG,51,676,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  BECAUSE IRES WAS EQUAL TO MINUS ONE'
      CALL XERRWD(MSG,44,677,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Failure because IRES = -2.
C
680   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2)'
      CALL XERRWD(MSG,40,680,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  IRES WAS EQUAL TO MINUS TWO'
      CALL XERRWD(MSG,36,681,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Failed to compute initial YPRIME.
C
685   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,685,0,0,0,0,0,0.0D0,0.0D0)
      MSG = 'DASPK--  INITIAL (Y,YPRIME) COULD NOT BE COMPUTED'
      CALL XERRWD(MSG,49,686,0,0,0,0,2,TN,H0)
      GO TO 700
C
C     Failure because IER was negative from PSOL.
C
690   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2)'
      CALL XERRWD(MSG,40,690,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  IER WAS NEGATIVE FROM PSOL'
      CALL XERRWD(MSG,35,691,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C     Failure because the linear system solver could not converge.
C
695   MSG = 'DASPK--  AT T (=R1) AND STEPSIZE H (=R2) THE'
      CALL XERRWD(MSG,44,695,0,0,0,0,2,TN,H)
      MSG = 'DASPK--  LINEAR SYSTEM SOLVER COULD NOT CONVERGE.'
      CALL XERRWD(MSG,50,696,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 700
C
C
700   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     This block handles all error returns due to illegal input,
C     as detected before calling DDSTP.
C     First the error message routine is called.  If this happens
C     twice in succession, execution is terminated.
C-----------------------------------------------------------------------
C
701   MSG = 'DASPK--  ELEMENT (=I1) OF INFO VECTOR IS NOT VALID'
      CALL XERRWD(MSG,50,1,0,1,ITEMP,0,0,0.0D0,0.0D0)
      GO TO 750
702   MSG = 'DASPK--  NEQ (=I1) .LE. 0'
      CALL XERRWD(MSG,25,2,0,1,NEQ,0,0,0.0D0,0.0D0)
      GO TO 750
703   MSG = 'DASPK--  MAXORD (=I1) NOT IN RANGE'
      CALL XERRWD(MSG,34,3,0,1,MXORD,0,0,0.0D0,0.0D0)
      GO TO 750
704   MSG='DASPK--  RWORK LENGTH NEEDED, LENRW (=I1), EXCEEDS LRW (=I2)'
      CALL XERRWD(MSG,60,4,0,2,LENRW,LRW,0,0.0D0,0.0D0)
      GO TO 750
705   MSG='DASPK--  IWORK LENGTH NEEDED, LENIW (=I1), EXCEEDS LIW (=I2)'
      CALL XERRWD(MSG,60,5,0,2,LENIW,LIW,0,0.0D0,0.0D0)
      GO TO 750
706   MSG = 'DASPK--  SOME ELEMENT OF RTOL IS .LT. 0'
      CALL XERRWD(MSG,39,6,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
707   MSG = 'DASPK--  SOME ELEMENT OF ATOL IS .LT. 0'
      CALL XERRWD(MSG,39,7,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
708   MSG = 'DASPK--  ALL ELEMENTS OF RTOL AND ATOL ARE ZERO'
      CALL XERRWD(MSG,47,8,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
709   MSG='DASPK--  INFO(4) = 1 AND TSTOP (=R1) BEHIND TOUT (=R2)'
      CALL XERRWD(MSG,54,9,0,0,0,0,2,TSTOP,TOUT)
      GO TO 750
710   MSG = 'DASPK--  HMAX (=R1) .LT. 0.0'
      CALL XERRWD(MSG,28,10,0,0,0,0,1,HMAX,0.0D0)
      GO TO 750
711   MSG = 'DASPK--  TOUT (=R1) BEHIND T (=R2)'
      CALL XERRWD(MSG,34,11,0,0,0,0,2,TOUT,T)
      GO TO 750
712   MSG = 'DASPK--  INFO(8)=1 AND H0=0.0'
      CALL XERRWD(MSG,29,12,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
713   MSG = 'DASPK--  SOME ELEMENT OF WT IS .LE. 0.0'
      CALL XERRWD(MSG,39,13,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
714   MSG='DASPK-- TOUT (=R1) TOO CLOSE TO T (=R2) TO START INTEGRATION'
      CALL XERRWD(MSG,60,14,0,0,0,0,2,TOUT,T)
      GO TO 750
715   MSG = 'DASPK--  INFO(4)=1 AND TSTOP (=R1) BEHIND T (=R2)'
      CALL XERRWD(MSG,49,15,0,0,0,0,2,TSTOP,T)
      GO TO 750
717   MSG = 'DASPK--  ML (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ'
      CALL XERRWD(MSG,52,17,0,1,IWORK(LML),0,0,0.0D0,0.0D0)
      GO TO 750
718   MSG = 'DASPK--  MU (=I1) ILLEGAL. EITHER .LT. 0 OR .GT. NEQ'
      CALL XERRWD(MSG,52,18,0,1,IWORK(LMU),0,0,0.0D0,0.0D0)
      GO TO 750
719   MSG = 'DASPK--  TOUT (=R1) IS EQUAL TO T (=R2)'
      CALL XERRWD(MSG,39,19,0,0,0,0,2,TOUT,T)
      GO TO 750
720   MSG = 'DASPK--  MAXL (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. NEQ'
      CALL XERRWD(MSG,54,20,0,1,IWORK(LMAXL),0,0,0.0D0,0.0D0)
      GO TO 750
721   MSG = 'DASPK--  KMP (=I1) ILLEGAL. EITHER .LT. 1 OR .GT. MAXL'
      CALL XERRWD(MSG,54,21,0,1,IWORK(LKMP),0,0,0.0D0,0.0D0)
      GO TO 750
722   MSG = 'DASPK--  NRMAX (=I1) ILLEGAL. .LT. 0'
      CALL XERRWD(MSG,36,22,0,1,IWORK(LNRMAX),0,0,0.0D0,0.0D0)
      GO TO 750
723   MSG = 'DASPK--  EPLI (=R1) ILLEGAL. EITHER .LE. 0.D0 OR .GE. 1.D0'
      CALL XERRWD(MSG,58,23,0,0,0,0,1,RWORK(LEPLI),0.0D0)
      GO TO 750
724   MSG = 'DASPK--  ILLEGAL IWORK VALUE FOR INFO(11) .NE. 0'
      CALL XERRWD(MSG,48,24,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
725   MSG = 'DASPK--  ONE OF THE INPUTS FOR INFO(17) = 1 IS ILLEGAL'
      CALL XERRWD(MSG,54,25,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
726   MSG = 'DASPK--  ILLEGAL IWORK VALUE FOR INFO(10) .NE. 0'
      CALL XERRWD(MSG,48,26,0,0,0,0,0,0.0D0,0.0D0)
      GO TO 750
727   MSG = 'DASPK--  Y(I) AND IWORK(40+I) (I=I1) INCONSISTENT'
      CALL XERRWD(MSG,49,27,0,1,IRET,0,0,0.0D0,0.0D0)
      GO TO 750
750   IF(INFO(1).EQ.-1) GO TO 760
      INFO(1)=-1
      IDID=-33
      RETURN
760   MSG = 'DASPK--  REPEATED OCCURRENCES OF ILLEGAL INPUT'
      CALL XERRWD(MSG,46,701,0,0,0,0,0,0.0D0,0.0D0)
770   MSG = 'DASPK--  RUN TERMINATED. APPARENT INFINITE LOOP'
      CALL XERRWD(MSG,47,702,1,0,0,0,0,0.0D0,0.0D0)
      RETURN
C
C------END OF SUBROUTINE DDASPK-----------------------------------------
      END
      SUBROUTINE DDASIC (X, Y, YPRIME, NEQ, ICOPT, ID, RES, JAC, PSOL,
     *   H, TSCALE, WT, NIC, IDID, RPAR, IPAR, PHI, SAVR, DELTA, E,
     *   YIC, YPIC, PWK, WM, IWM, UROUND, EPLI, SQRTN, RSQRTN,
     *   EPCONI, STPTOL, JFLG, ICNFLG, ICNSTR, NLSIC)
C
C***BEGIN PROLOGUE  DDASIC
C***REFER TO  DDASPK
C***DATE WRITTEN   940628   (YYMMDD)
C***REVISION DATE  941206   (YYMMDD)
C***REVISION DATE  950714   (YYMMDD)
C***REVISION DATE  000628   TSCALE argument added.
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DDASIC is a driver routine to compute consistent initial values
C     for Y and YPRIME.  There are two different options:  
C     Denoting the differential variables in Y by Y_d, and
C     the algebraic variables by Y_a, the problem solved is either:
C     1.  Given Y_d, calculate Y_a and Y_d', or
C     2.  Given Y', calculate Y.
C     In either case, initial values for the given components
C     are input, and initial guesses for the unknown components
C     must also be provided as input.
C
C     The external routine NLSIC solves the resulting nonlinear system.
C
C     The parameters represent
C
C     X  --        Independent variable.
C     Y  --        Solution vector at X.
C     YPRIME --    Derivative of solution vector.
C     NEQ --       Number of equations to be integrated.
C     ICOPT     -- Flag indicating initial condition option chosen.
C                    ICOPT = 1 for option 1 above.
C                    ICOPT = 2 for option 2.
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if option 1 is chosen.
C                    ID(i) = +1 if Y_i is a differential variable,
C                    ID(i) = -1 if Y_i is an algebraic variable. 
C     RES --       External user-supplied subroutine to evaluate the
C                  residual.  See RES description in DDASPK prologue.
C     JAC --       External user-supplied routine to update Jacobian
C                  or preconditioner information in the nonlinear solver
C                  (optional).  See JAC description in DDASPK prologue.
C     PSOL --      External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  See PSOL in DDASPK prologue.
C     H --         Scaling factor in iteration matrix.  DDASIC may 
C                  reduce H to achieve convergence.
C     TSCALE --    Scale factor in T, used for stopping tests if nonzero.
C     WT --        Vector of weights for error criterion.
C     NIC --       Input number of initial condition calculation call 
C                  (= 1 or 2).
C     IDID --      Completion code.  See IDID in DDASPK prologue.
C     RPAR,IPAR -- Real and integer parameter arrays that
C                  are used for communication between the
C                  calling program and external user routines.
C                  They are not altered by DNSK
C     PHI --       Work space for DDASIC of length at least 2*NEQ.
C     SAVR --      Work vector for DDASIC of length NEQ.
C     DELTA --     Work vector for DDASIC of length NEQ.
C     E --         Work vector for DDASIC of length NEQ.
C     YIC,YPIC --  Work vectors for DDASIC, each of length NEQ.
C     PWK --       Work vector for DDASIC of length NEQ.
C     WM,IWM --    Real and integer arrays storing
C                  information required by the linear solver.
C     EPCONI --    Test constant for Newton iteration convergence.
C     ICNFLG --    Flag showing whether constraints on Y are to apply.
C     ICNSTR --    Integer array of length NEQ with constraint types.
C
C     The other parameters are for use internally by DDASIC.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DCOPY, NLSIC
C
C***END PROLOGUE  DDASIC
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),ID(*),WT(*),PHI(NEQ,*)
      DIMENSION SAVR(*),DELTA(*),E(*),YIC(*),YPIC(*),PWK(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*), ICNSTR(*)
      EXTERNAL RES, JAC, PSOL, NLSIC
C
      PARAMETER (LCFN=15)
      PARAMETER (LMXNH=34)
C
C The following parameters are data-loaded here:
C     RHCUT  = factor by which H is reduced on retry of Newton solve.
C     RATEMX = maximum convergence rate for which Newton iteration
C              is considered converging.
C
      SAVE RHCUT, RATEMX
      DATA RHCUT/0.1D0/, RATEMX/0.8D0/
C
C
C-----------------------------------------------------------------------
C     BLOCK 1.
C     Initializations.
C     JSKIP is a flag set to 1 when NIC = 2 and NH = 1, to signal that
C     the initial call to the JAC routine is to be skipped then.
C     Save Y and YPRIME in PHI.  Initialize IDID, NH, and CJ.
C-----------------------------------------------------------------------
C
      MXNH = IWM(LMXNH)
      IDID = 1
      NH = 1
      JSKIP = 0
      IF (NIC .EQ. 2) JSKIP = 1
      CALL DCOPY (NEQ, Y, 1, PHI(1,1), 1)
      CALL DCOPY (NEQ, YPRIME, 1, PHI(1,2), 1)
C
      IF (ICOPT .EQ. 2) THEN
        CJ = 0.0D0 
      ELSE
        CJ = 1.0D0/H
      ENDIF
C
C-----------------------------------------------------------------------
C     BLOCK 2
C     Call the nonlinear system solver to obtain
C     consistent initial values for Y and YPRIME.
C-----------------------------------------------------------------------
C
 200  CONTINUE
      CALL NLSIC(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JAC,PSOL,H,TSCALE,WT,
     *   JSKIP,RPAR,IPAR,SAVR,DELTA,E,YIC,YPIC,PWK,WM,IWM,CJ,UROUND,
     *   EPLI,SQRTN,RSQRTN,EPCONI,RATEMX,STPTOL,JFLG,ICNFLG,ICNSTR,
     *   IERNLS)
C
      IF (IERNLS .EQ. 0) RETURN
C
C-----------------------------------------------------------------------
C     BLOCK 3
C     The nonlinear solver was unsuccessful.  Increment NCFN.
C     Return with IDID = -12 if either
C       IERNLS = -1: error is considered unrecoverable,
C       ICOPT = 2: we are doing initialization problem type 2, or
C       NH = MXNH: the maximum number of H values has been tried.
C     Otherwise (problem 1 with IERNLS .GE. 1), reduce H and try again.
C     If IERNLS > 1, restore Y and YPRIME to their original values.
C-----------------------------------------------------------------------
C
      IWM(LCFN) = IWM(LCFN) + 1
      JSKIP = 0
C
      IF (IERNLS .EQ. -1) GO TO 350
      IF (ICOPT .EQ. 2) GO TO 350
      IF (NH .EQ. MXNH) GO TO 350
C
      NH = NH + 1
      H = H*RHCUT
      CJ = 1.0D0/H
C
      IF (IERNLS .EQ. 1) GO TO 200
C
      CALL DCOPY (NEQ, PHI(1,1), 1, Y, 1)
      CALL DCOPY (NEQ, PHI(1,2), 1, YPRIME, 1)
      GO TO 200
C
 350  IDID = -12
      RETURN
C
C------END OF SUBROUTINE DDASIC-----------------------------------------
      END
      SUBROUTINE DYYPNW (NEQ, Y, YPRIME, CJ, RL, P, ICOPT, ID, 
     *                   YNEW, YPNEW)
C
C***BEGIN PROLOGUE  DYYPNW
C***REFER TO  DLINSK
C***DATE WRITTEN   940830   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DYYPNW calculates the new (Y,YPRIME) pair needed in the
C     linesearch algorithm based on the current lambda value.  It is
C     called by DLINSK and DLINSD.  Based on the ICOPT and ID values,
C     the corresponding entry in Y or YPRIME is updated.
C
C     In addition to the parameters described in the calling programs,
C     the parameters represent
C
C     P      -- Array of length NEQ that contains the current
C               approximate Newton step.
C     RL     -- Scalar containing the current lambda value.
C     YNEW   -- Array of length NEQ containing the updated Y vector.
C     YPNEW  -- Array of length NEQ containing the updated YPRIME
C               vector.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED (NONE)
C
C***END PROLOGUE  DYYPNW
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(*), YPRIME(*), YNEW(*), YPNEW(*), ID(*), P(*)
C
      IF (ICOPT .EQ. 1) THEN
         DO 10 I=1,NEQ
            IF(ID(I) .LT. 0) THEN
               YNEW(I) = Y(I) - RL*P(I)
               YPNEW(I) = YPRIME(I)
            ELSE
               YNEW(I) = Y(I)
               YPNEW(I) = YPRIME(I) - RL*CJ*P(I)
            ENDIF
 10      CONTINUE
      ELSE
         DO 20 I = 1,NEQ
            YNEW(I) = Y(I) - RL*P(I)
            YPNEW(I) = YPRIME(I)
 20      CONTINUE
      ENDIF
      RETURN
C----------------------- END OF SUBROUTINE DYYPNW ----------------------
      END
      SUBROUTINE DDSTP(X,Y,YPRIME,NEQ,RES,JAC,PSOL,H,WT,VT,
     *  JSTART,IDID,RPAR,IPAR,PHI,SAVR,DELTA,E,WM,IWM,
     *  ALPHA,BETA,GAMMA,PSI,SIGMA,CJ,CJOLD,HOLD,S,HMIN,UROUND,
     *  EPLI,SQRTN,RSQRTN,EPCON,IPHASE,JCALC,JFLG,K,KOLD,NS,NONNEG,
     *  NTYPE,NLS)
C
C***BEGIN PROLOGUE  DDSTP
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940909   (YYMMDD) (Reset PSI(1), PHI(*,2) at 690)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DDSTP solves a system of differential/algebraic equations of 
C     the form G(X,Y,YPRIME) = 0, for one step (normally from X to X+H).
C
C     The methods used are modified divided difference, fixed leading 
C     coefficient forms of backward differentiation formulas.  
C     The code adjusts the stepsize and order to control the local error
C     per step.
C
C
C     The parameters represent
C     X  --        Independent variable.
C     Y  --        Solution vector at X.
C     YPRIME --    Derivative of solution vector
C                  after successful step.
C     NEQ --       Number of equations to be integrated.
C     RES --       External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JAC --       External user-supplied routine to update
C                  Jacobian or preconditioner information in the
C                  nonlinear solver.  See JAC description in DDASPK
C                  prologue.
C     PSOL --      External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  (This is optional).  See PSOL in DDASPK prologue.
C     H --         Appropriate step size for next step.
C                  Normally determined by the code.
C     WT --        Vector of weights for error criterion used in Newton test.
C     VT --        Masked vector of weights used in error test.
C     JSTART --    Integer variable set 0 for
C                  first step, 1 otherwise.
C     IDID --      Completion code returned from the nonlinear solver.
C                  See IDID description in DDASPK prologue.
C     RPAR,IPAR -- Real and integer parameter arrays that
C                  are used for communication between the
C                  calling program and external user routines.
C                  They are not altered by DNSK
C     PHI --       Array of divided differences used by
C                  DDSTP. The length is NEQ*(K+1), where
C                  K is the maximum order.
C     SAVR --      Work vector for DDSTP of length NEQ.
C     DELTA,E --   Work vectors for DDSTP of length NEQ.
C     WM,IWM --    Real and integer arrays storing
C                  information required by the linear solver.
C
C     The other parameters are information
C     which is needed internally by DDSTP to
C     continue from step to step.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   NLS, DDWNRM, DDATRP
C
C***END PROLOGUE  DDSTP
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),VT(*)
      DIMENSION PHI(NEQ,*),SAVR(*),DELTA(*),E(*)
      DIMENSION WM(*),IWM(*)
      DIMENSION PSI(*),ALPHA(*),BETA(*),GAMMA(*),SIGMA(*)
      DIMENSION RPAR(*),IPAR(*)
      EXTERNAL  RES, JAC, PSOL, NLS
C
      PARAMETER (LMXORD=3)
      PARAMETER (LNST=11, LETF=14, LCFN=15)
C
C
C-----------------------------------------------------------------------
C     BLOCK 1.
C     Initialize.  On the first call, set
C     the order to 1 and initialize
C     other variables.
C-----------------------------------------------------------------------
C
C     Initializations for all calls
C
      XOLD=X
      NCF=0
      NEF=0
      IF(JSTART .NE. 0) GO TO 120
C
C     If this is the first step, perform
C     other initializations
C
      K=1
      KOLD=0
      HOLD=0.0D0
      PSI(1)=H
      CJ = 1.D0/H
      IPHASE = 0
      NS=0
120   CONTINUE
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 2
C     Compute coefficients of formulas for
C     this step.
C-----------------------------------------------------------------------
200   CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN0(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
C
      BETA(1)=1.0D0
      ALPHA(1)=1.0D0
      TEMP1=H
      GAMMA(1)=0.0D0
      SIGMA(1)=1.0D0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210      CONTINUE
      PSI(KP1)=TEMP1
230   CONTINUE
C
C     Compute ALPHAS, ALPHA0
C
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0D0/I
        ALPHA0 = ALPHA0 - ALPHA(I)
240     CONTINUE
C
C     Compute leading coefficient CJ
C
      CJLAST = CJ
      CJ = -ALPHAS/H
C
C     Compute variable stepsize error coefficient CK
C
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
C
C     Change PHI to PHI STAR
C
      IF(KP1 .LT. NSP1) GO TO 280
      DO 270 J=NSP1,KP1
         DO 260 I=1,NEQ
260         PHI(I,J)=BETA(J)*PHI(I,J)
270      CONTINUE
280   CONTINUE
C
C     Update time
C
      X=X+H
C
C     Initialize IDID to 1
C
      IDID = 1
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 3
C     Call the nonlinear system solver to obtain the solution and
C     derivative.
C-----------------------------------------------------------------------
C
      CALL NLS(X,Y,YPRIME,NEQ,
     *   RES,JAC,PSOL,H,WT,JSTART,IDID,RPAR,IPAR,PHI,GAMMA,
     *   SAVR,DELTA,E,WM,IWM,CJ,CJOLD,CJLAST,S,
     *   UROUND,EPLI,SQRTN,RSQRTN,EPCON,JCALC,JFLG,KP1,
     *   NONNEG,NTYPE,IERNLS)
C
      IF(IERNLS .NE. 0)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 4
C     Estimate the errors at orders K,K-1,K-2
C     as if constant stepsize was used. Estimate
C     the local error at order K and test
C     whether the current step is successful.
C-----------------------------------------------------------------------
C
C     Estimate errors at orders K,K-1,K-2
C
      ENORM = DDWNRM(NEQ,E,VT,RPAR,IPAR)
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1)GO TO 430
      DO 405 I = 1,NEQ
405     DELTA(I) = PHI(I,KP1) + E(I)
      ERKM1=SIGMA(K)*DDWNRM(NEQ,DELTA,VT,RPAR,IPAR)
      TERKM1 = K*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5*TERK)GO TO 420
      GO TO 430
410   CONTINUE
      DO 415 I = 1,NEQ
415     DELTA(I) = PHI(I,K) + DELTA(I)
      ERKM2=SIGMA(K-1)*DDWNRM(NEQ,DELTA,VT,RPAR,IPAR)
      TERKM2 = (K-1)*ERKM2
      IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
C
C     Lower the order
C
420   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C
C     Calculate the local error for the current step
C     to see if the step was successful
C
430   CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. 1.0D0)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 5
C     The step is successful. Determine
C     the best order and stepsize for
C     the next step. Update the differences
C     for the next step.
C-----------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C
C     Estimate the error at order K+1 unless
C        already decided to lower order, or
C        already using maximum order, or
C        stepsize not constant, or
C        order raised in previous step
C
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NEQ
510      DELTA(I)=E(I)-PHI(I,KP2)
      ERKP1 = (1.0D0/(K+2))*DDWNRM(NEQ,DELTA,VT,RPAR,IPAR)
      TERKP1 = (K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550
      GO TO 530
520   IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
C
C     Raise order
C
530   K=KP1
      EST = ERKP1
      GO TO 550
C
C     Lower order
C
540   K=KM1
      EST = ERKM1
      GO TO 550
C
C     If IPHASE = 0, increase order by one and multiply stepsize by
C     factor two
C
545   K = KP1
      HNEW = H*2.0D0
      H = HNEW
      GO TO 575
C
C
C     Determine the appropriate stepsize for
C     the next step.
C
550   HNEW=H
      TEMP2=K+1
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF(R .LT. 2.0D0) GO TO 555
      HNEW = 2.0D0*H
      GO TO 560
555   IF(R .GT. 1.0D0) GO TO 560
      R = MAX(0.5D0,MIN(0.9D0,R))
      HNEW = H*R
560   H=HNEW
C
C
C     Update differences for next step
C
575   CONTINUE
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
      DO 580 I=1,NEQ
580      PHI(I,KP2)=E(I)
585   CONTINUE
      DO 590 I=1,NEQ
590      PHI(I,KP1)=PHI(I,KP1)+E(I)
      DO 595 J1=2,KP1
         J=KP1-J1+1
         DO 595 I=1,NEQ
595      PHI(I,J)=PHI(I,J)+PHI(I,J+1)
      JSTART = 1
      RETURN
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 6
C     The step is unsuccessful. Restore X,PSI,PHI
C     Determine appropriate stepsize for
C     continuing the integration, or exit with
C     an error flag if there have been many
C     failures.
C-----------------------------------------------------------------------
600   IPHASE = 1
C
C     Restore X,PHI,PSI
C
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0D0/BETA(J)
         DO 610 I=1,NEQ
610         PHI(I,J)=TEMP1*PHI(I,J)
620      CONTINUE
630   CONTINUE
      DO 640 I=2,KP1
640      PSI(I-1)=PSI(I)-H
C
C
C     Test whether failure is due to nonlinear solver
C     or error test
C
      IF(IERNLS .EQ. 0)GO TO 660
      IWM(LCFN)=IWM(LCFN)+1
C
C
C     The nonlinear solver failed to converge.
C     Determine the cause of the failure and take appropriate action.
C     If IERNLS .LT. 0, then return.  Otherwise, reduce the stepsize
C     and try again, unless too many failures have occurred.
C
      IF (IERNLS .LT. 0) GO TO 675
      NCF = NCF + 1
      R = 0.25D0
      H = H*R
      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
      IF (IDID .EQ. 1) IDID = -7
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
C
C
C     The nonlinear solver converged, and the cause
C     of the failure was the error estimate
C     exceeding the tolerance.
C
660   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) GO TO 665
C
C     On first error test failure, keep current order or lower
C     order by one.  Compute new stepsize based on differences
C     of the solution.
C
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = MAX(0.25D0,MIN(0.9D0,R))
      H = H*R
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     On second error test failure, use the current order or
C     decrease order by one.  Reduce the stepsize by a factor of
C     one quarter.
C
665   IF (NEF .GT. 2) GO TO 670
      K = KNEW
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     On third and subsequent error test failures, set the order to
C     one, and reduce the stepsize by a factor of one quarter.
C
670   K = 1
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C
C
C
C     For all crashes, restore Y to its last value,
C     interpolate to find YPRIME at last X, and return.
C
C     Before returning, verify that the user has not set
C     IDID to a nonnegative value.  If the user has set IDID
C     to a nonnegative value, then reset IDID to be -7, indicating
C     a failure in the nonlinear system solver.
C
675   CONTINUE
      CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      JSTART = 1
      IF (IDID .GE. 0) IDID = -7
      RETURN
C
C
C     Go back and try this step again.  
C     If this is the first step, reset PSI(1) and rescale PHI(*,2).
C
690   IF (KOLD .EQ. 0) THEN
        PSI(1) = H
        DO 695 I = 1,NEQ
695       PHI(I,2) = R*PHI(I,2)
        ENDIF
      GO TO 200
C
C------END OF SUBROUTINE DDSTP------------------------------------------
      END
      SUBROUTINE DCNSTR (NEQ, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
C
C***BEGIN PROLOGUE  DCNSTR
C***DATE WRITTEN   950808   (YYMMDD)
C***REVISION DATE  950814   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This subroutine checks for constraint violations in the proposed 
C new approximate solution YNEW.
C If a constraint violation occurs, then a new step length, TAU,
C is calculated, and this value is to be given to the linesearch routine
C to calculate a new approximate solution YNEW.
C
C On entry:
C
C   NEQ    -- size of the nonlinear system, and the length of arrays
C             Y, YNEW and ICNSTR.
C
C   Y      -- real array containing the current approximate y.
C
C   YNEW   -- real array containing the new approximate y.
C
C   ICNSTR -- INTEGER array of length NEQ containing flags indicating
C             which entries in YNEW are to be constrained.
C             if ICNSTR(I) =  2, then YNEW(I) must be .GT. 0,
C             if ICNSTR(I) =  1, then YNEW(I) must be .GE. 0,
C             if ICNSTR(I) = -1, then YNEW(I) must be .LE. 0, while
C             if ICNSTR(I) = -2, then YNEW(I) must be .LT. 0, while
C             if ICNSTR(I) =  0, then YNEW(I) is not constrained.
C
C   RLX    -- real scalar restricting update, if ICNSTR(I) = 2 or -2,
C             to ABS( (YNEW-Y)/Y ) < FAC2*RLX in component I.
C
C   TAU    -- the current size of the step length for the linesearch.
C
C On return
C
C   TAU    -- the adjusted size of the step length if a constraint
C             violation occurred (otherwise, it is unchanged).  it is
C             the step length to give to the linesearch routine.
C
C   IRET   -- output flag.
C             IRET=0 means that YNEW satisfied all constraints.
C             IRET=1 means that YNEW failed to satisfy all the
C                    constraints, and a new linesearch step
C                    must be computed.
C
C   IVAR   -- index of variable causing constraint to be violated.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(NEQ), YNEW(NEQ), ICNSTR(NEQ)
      SAVE FAC, FAC2, ZERO
      DATA FAC /0.6D0/, FAC2 /0.9D0/, ZERO/0.0D0/
C-----------------------------------------------------------------------
C Check constraints for proposed new step YNEW.  If a constraint has
C been violated, then calculate a new step length, TAU, to be
C used in the linesearch routine.
C-----------------------------------------------------------------------
      IRET = 0
      RDYMX = ZERO
      IVAR = 0
      DO 100 I = 1,NEQ
C
         IF (ICNSTR(I) .EQ. 2) THEN
            RDY = ABS( (YNEW(I)-Y(I))/Y(I) )
            IF (RDY .GT. RDYMX) THEN
               RDYMX = RDY
               IVAR = I
            ENDIF
            IF (YNEW(I) .LE. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ELSEIF (ICNSTR(I) .EQ. 1) THEN
            IF (YNEW(I) .LT. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ELSEIF (ICNSTR(I) .EQ. -1) THEN
            IF (YNEW(I) .GT. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ELSEIF (ICNSTR(I) .EQ. -2) THEN
            RDY = ABS( (YNEW(I)-Y(I))/Y(I) )
            IF (RDY .GT. RDYMX) THEN
               RDYMX = RDY
               IVAR = I
            ENDIF
            IF (YNEW(I) .GE. ZERO) THEN
               TAU = FAC*TAU
               IVAR = I
               IRET = 1
               RETURN
            ENDIF
C
         ENDIF
 100  CONTINUE

      IF(RDYMX .GE. RLX) THEN
         TAU = FAC2*TAU*RLX/RDYMX
         IRET = 1
      ENDIF
C
      RETURN
C----------------------- END OF SUBROUTINE DCNSTR ----------------------
      END
      SUBROUTINE DCNST0 (NEQ, Y, ICNSTR, IRET)
C
C***BEGIN PROLOGUE  DCNST0
C***DATE WRITTEN   950808   (YYMMDD)
C***REVISION DATE  950808   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This subroutine checks for constraint violations in the initial 
C approximate solution u.
C
C On entry
C
C   NEQ    -- size of the nonlinear system, and the length of arrays
C             Y and ICNSTR.
C
C   Y      -- real array containing the initial approximate root.
C
C   ICNSTR -- INTEGER array of length NEQ containing flags indicating
C             which entries in Y are to be constrained.
C             if ICNSTR(I) =  2, then Y(I) must be .GT. 0,
C             if ICNSTR(I) =  1, then Y(I) must be .GE. 0,
C             if ICNSTR(I) = -1, then Y(I) must be .LE. 0, while
C             if ICNSTR(I) = -2, then Y(I) must be .LT. 0, while
C             if ICNSTR(I) =  0, then Y(I) is not constrained.
C
C On return
C
C   IRET   -- output flag.
C             IRET=0    means that u satisfied all constraints.
C             IRET.NE.0 means that Y(IRET) failed to satisfy its
C                       constraint.
C
C-----------------------------------------------------------------------
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(NEQ), ICNSTR(NEQ)
      SAVE ZERO
      DATA ZERO/0.D0/
C-----------------------------------------------------------------------
C Check constraints for initial Y.  If a constraint has been violated,
C set IRET = I to signal an error return to calling routine.
C-----------------------------------------------------------------------
      IRET = 0
      DO 100 I = 1,NEQ
         IF (ICNSTR(I) .EQ. 2) THEN
            IF (Y(I) .LE. ZERO) THEN
               IRET = I
               RETURN
            ENDIF
         ELSEIF (ICNSTR(I) .EQ. 1) THEN
            IF (Y(I) .LT. ZERO) THEN
               IRET = I
               RETURN
            ENDIF 
         ELSEIF (ICNSTR(I) .EQ. -1) THEN
            IF (Y(I) .GT. ZERO) THEN
               IRET = I
               RETURN
            ENDIF 
         ELSEIF (ICNSTR(I) .EQ. -2) THEN
            IF (Y(I) .GE. ZERO) THEN
               IRET = I
               RETURN
            ENDIF 
        ENDIF
 100  CONTINUE
      RETURN
C----------------------- END OF SUBROUTINE DCNST0 ----------------------
      END
      SUBROUTINE DDAWTS(NEQ,IWT,RTOL,ATOL,Y,WT,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DDAWTS
C***REFER TO  DDASPK
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***END PROLOGUE  DDAWTS
C-----------------------------------------------------------------------
C     This subroutine sets the error weight vector,
C     WT, according to WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
C     I = 1 to NEQ.
C     RTOL and ATOL are scalars if IWT = 0,
C     and vectors if IWT = 1.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION RTOL(*),ATOL(*),Y(*),WT(*)
      DIMENSION RPAR(*),IPAR(*)
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 20 I=1,NEQ
         IF (IWT .EQ.0) GO TO 10
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
10         WT(I)=RTOLI*ABS(Y(I))+ATOLI
20         CONTINUE
      RETURN
C
C------END OF SUBROUTINE DDAWTS-----------------------------------------
      END
      SUBROUTINE DINVWT(NEQ,WT,IER)
C
C***BEGIN PROLOGUE  DINVWT
C***REFER TO  DDASPK
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   950125   (YYMMDD)
C***END PROLOGUE  DINVWT
C-----------------------------------------------------------------------
C     This subroutine checks the error weight vector WT, of length NEQ,
C     for components that are .le. 0, and if none are found, it
C     inverts the WT(I) in place.  This replaces division operations
C     with multiplications in all norm evaluations.
C     IER is returned as 0 if all WT(I) were found positive,
C     and the first I with WT(I) .le. 0.0 otherwise.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION WT(*)
C
      DO 10 I = 1,NEQ
        IF (WT(I) .LE. 0.0D0) GO TO 30
 10     CONTINUE
      DO 20 I = 1,NEQ
 20     WT(I) = 1.0D0/WT(I)
      IER = 0
      RETURN
C
 30   IER = I
      RETURN
C
C------END OF SUBROUTINE DINVWT-----------------------------------------
      END
      SUBROUTINE DDATRP(X,XOUT,YOUT,YPOUT,NEQ,KOLD,PHI,PSI)
C
C***BEGIN PROLOGUE  DDATRP
C***REFER TO  DDASPK
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***END PROLOGUE  DDATRP
C
C-----------------------------------------------------------------------
C     The methods in subroutine DDSTP use polynomials
C     to approximate the solution.  DDATRP approximates the
C     solution and its derivative at time XOUT by evaluating
C     one of these polynomials, and its derivative, there.
C     Information defining this polynomial is passed from
C     DDSTP, so DDATRP cannot be used alone.
C
C     The parameters are
C
C     X     The current time in the integration.
C     XOUT  The time at which the solution is desired.
C     YOUT  The interpolated approximation to Y at XOUT.
C           (This is output.)
C     YPOUT The interpolated approximation to YPRIME at XOUT.
C           (This is output.)
C     NEQ   Number of equations.
C     KOLD  Order used on last successful step.
C     PHI   Array of scaled divided differences of Y.
C     PSI   Array of past stepsize history.
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION YOUT(*),YPOUT(*)
      DIMENSION PHI(NEQ,*),PSI(*)
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
10       YPOUT(I)=0.0D0
      C=1.0D0
      D=0.0D0
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
         DO 20 I=1,NEQ
            YOUT(I)=YOUT(I)+C*PHI(I,J)
20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
30       CONTINUE
      RETURN
C
C------END OF SUBROUTINE DDATRP-----------------------------------------
      END
      DOUBLE PRECISION FUNCTION DDWNRM(NEQ,V,RWT,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DDWNRM
C***ROUTINES CALLED  (NONE)
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***END PROLOGUE  DDWNRM
C-----------------------------------------------------------------------
C     This function routine computes the weighted
C     root-mean-square norm of the vector of length
C     NEQ contained in the array V, with reciprocal weights
C     contained in the array RWT of length NEQ.
C        DDWNRM=SQRT((1/NEQ)*SUM(V(I)*RWT(I))**2)
C-----------------------------------------------------------------------
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION V(*),RWT(*)
      DIMENSION RPAR(*),IPAR(*)
      DDWNRM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
        IF(ABS(V(I)*RWT(I)) .GT. VMAX) VMAX = ABS(V(I)*RWT(I))
10    CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)*RWT(I))/VMAX)**2
      DDWNRM = VMAX*SQRT(SUM/NEQ)
30    CONTINUE
      RETURN
C
C------END OF FUNCTION DDWNRM-------------------------------------------
      END
      SUBROUTINE DDASID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JACD,PDUM,H,TSCALE,
     *  WT,JSDUM,RPAR,IPAR,DUMSVR,DELTA,R,YIC,YPIC,DUMPWK,WM,IWM,CJ,
     *  UROUND,DUME,DUMS,DUMR,EPCON,RATEMX,STPTOL,JFDUM,
     *  ICNFLG,ICNSTR,IERNLS)
C
C***BEGIN PROLOGUE  DDASID
C***REFER TO  DDASPK
C***DATE WRITTEN   940701   (YYMMDD)
C***REVISION DATE  950808   (YYMMDD)
C***REVISION DATE  951110   Removed unreachable block 390.
C***REVISION DATE  000628   TSCALE argument added.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C
C     DDASID solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
C     the initial conditions.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     ICOPT     -- Initial condition option chosen (1 or 2).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied subroutine to evaluate the
C                  residual.  See RES description in DDASPK prologue.
C     JACD      -- External user-supplied routine to evaluate the
C                  Jacobian.  See JAC description for the case
C                  INFO(12) = 0 in the DDASPK prologue.
C     PDUM      -- Dummy argument.
C     H         -- Scaling factor for this initial condition calc.
C     TSCALE    -- Scale factor in T, used for stopping tests if nonzero.
C     WT        -- Vector of weights for error criterion.
C     JSDUM     -- Dummy argument.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     DUMSVR    -- Dummy argument.
C     DELTA     -- Work vector for NLS of length NEQ.
C     R         -- Work vector for NLS of length NEQ.
C     YIC,YPIC  -- Work vectors for NLS, each of length NEQ.
C     DUMPWK    -- Dummy argument.
C     WM,IWM    -- Real and integer arrays storing matrix information
C                  such as the matrix of partial derivatives,
C                  permutation vector, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
C     UROUND    -- Unit roundoff.
C     DUME      -- Dummy argument.
C     DUMS      -- Dummy argument.
C     DUMR      -- Dummy argument.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     JFDUM     -- Dummy argument.
C     STPTOL    -- Tolerance used in calculating the minimum lambda
C                  value allowed.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length 
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0   ==> nonlinear solver converged.
C                   1,2 ==> recoverable error inside nonlinear solver.
C                           1 => retry with current Y, YPRIME
C                           2 => retry with original Y, YPRIME
C                  -1   ==> unrecoverable error in nonlinear solver.
C
C     All variables with "DUM" in their names are dummy variables
C     which are not used in this routine.
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   RES, DMATD, DNSID
C
C***END PROLOGUE  DDASID
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),ID(*),WT(*),ICNSTR(*)
      DIMENSION DELTA(*),R(*),YIC(*),YPIC(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      EXTERNAL  RES, JACD
C
      PARAMETER (LNRE=12, LNJE=13, LMXNIT=32, LMXNJ=33)
C
C
C     Perform initializations.
C
      MXNIT = IWM(LMXNIT)
      MXNJ = IWM(LMXNJ)
      IERNLS = 0
      NJ = 0
C
C     Call RES to initialize DELTA.
C
      IRES = 0
      IWM(LNRE) = IWM(LNRE) + 1
      CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 370
C
C     Looping point for updating the Jacobian.
C
300   CONTINUE
C
C     Initialize all error flags to zero.
C
      IERJ = 0
      IRES = 0
      IERNEW = 0
C
C     Reevaluate the iteration matrix, J = dG/dY + CJ*dG/dYPRIME,
C     where G(X,Y,YPRIME) = 0.
C
      NJ = NJ + 1
      IWM(LNJE)=IWM(LNJE)+1
      CALL DMATD(NEQ,X,Y,YPRIME,DELTA,CJ,H,IERJ,WT,R,
     *              WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR)
      IF (IRES .LT. 0 .OR. IERJ .NE. 0) GO TO 370
C
C     Call the nonlinear Newton solver for up to MXNIT iterations.
C
      CALL DNSID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,WT,RPAR,IPAR,DELTA,R,
     *     YIC,YPIC,WM,IWM,CJ,TSCALE,EPCON,RATEMX,MXNIT,STPTOL,
     *     ICNFLG,ICNSTR,IERNEW)
C
      IF (IERNEW .EQ. 1 .AND. NJ .LT. MXNJ) THEN
C
C        MXNIT iterations were done, the convergence rate is < 1,
C        and the number of Jacobian evaluations is less than MXNJ.
C        Call RES, reevaluate the Jacobian, and try again.
C
         IWM(LNRE)=IWM(LNRE)+1
         CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) GO TO 370
         GO TO 300
         ENDIF
C
      IF (IERNEW .NE. 0) GO TO 380

      RETURN
C
C
C     Unsuccessful exits from nonlinear solver.
C     Compute IERNLS accordingly.
C
370   IERNLS = 2
      IF (IRES .LE. -2) IERNLS = -1
      RETURN
C
380   IERNLS = MIN(IERNEW,2)
      RETURN
C
C------END OF SUBROUTINE DDASID-----------------------------------------
      END
      SUBROUTINE DNSID(X,Y,YPRIME,NEQ,ICOPT,ID,RES,WT,RPAR,IPAR,
     *   DELTA,R,YIC,YPIC,WM,IWM,CJ,TSCALE,EPCON,RATEMX,MAXIT,STPTOL,
     *   ICNFLG,ICNSTR,IERNEW)
C
C***BEGIN PROLOGUE  DNSID
C***REFER TO  DDASPK
C***DATE WRITTEN   940701   (YYMMDD)
C***REVISION DATE  950713   (YYMMDD)
C***REVISION DATE  000628   TSCALE argument added.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSID solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME
C     in the initial conditions.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     ICOPT     -- Initial condition option chosen (1 or 2).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied subroutine to evaluate the
C                  residual.  See RES description in DDASPK prologue.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     DELTA     -- Residual vector on entry, and work vector of
C                  length NEQ for DNSID.
C     WM,IWM    -- Real and integer arrays storing matrix information
C                  such as the matrix of partial derivatives,
C                  permutation vector, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
C     TSCALE    -- Scale factor in T, used for stopping tests if nonzero.
C     R         -- Array of length NEQ used as workspace by the 
C                  linesearch routine DLINSD.
C     YIC,YPIC  -- Work vectors for DLINSD, each of length NEQ.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     STPTOL    -- Tolerance used in calculating the minimum lambda
C                  value allowed.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length 
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> failed to converge, but RATE .le. RATEMX.
C                   2  ==> failed to converge, RATE .gt. RATEMX.
C                   3  ==> other recoverable error (IRES = -1, or
C                          linesearch failed).
C                  -1  ==> unrecoverable error (IRES = -2).
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DSLVD, DDWNRM, DLINSD, DCOPY
C
C***END PROLOGUE  DNSID
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),R(*)
      DIMENSION ID(*),DELTA(*), YIC(*), YPIC(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION ICNSTR(*)
      EXTERNAL  RES
C
      PARAMETER (LNNI=19, LLSOFF=35)
C
C
C     Initializations.  M is the Newton iteration counter.
C
      LSOFF = IWM(LLSOFF)
      M = 0
      RATE = 1.0D0
      RLX = 0.4D0
C
C     Compute a new step vector DELTA by back-substitution.
C
      CALL DSLVD (NEQ, DELTA, WM, IWM)
C
C     Get norm of DELTA.  Return now if norm(DELTA) .le. EPCON.
C
      DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
      FNRM = DELNRM
      IF (TSCALE .GT. 0.0D0) FNRM = FNRM*TSCALE*ABS(CJ)
      IF (FNRM .LE. EPCON) RETURN
C
C     Newton iteration loop.
C
 300  CONTINUE
      IWM(LNNI) = IWM(LNNI) + 1
C
C     Call linesearch routine for global strategy and set RATE
C
      OLDFNM = FNRM
C
      CALL DLINSD (NEQ, Y, X, YPRIME, CJ, TSCALE, DELTA, DELNRM, WT,
     *             LSOFF, STPTOL, IRET, RES, IRES, WM, IWM, FNRM, ICOPT,
     *             ID, R, YIC, YPIC, ICNFLG, ICNSTR, RLX, RPAR, IPAR)
C
      RATE = FNRM/OLDFNM
C
C     Check for error condition from linesearch.
      IF (IRET .NE. 0) GO TO 390
C
C     Test for convergence of the iteration, and return or loop.
C
      IF (FNRM .LE. EPCON) RETURN
C
C     The iteration has not yet converged.  Update M.
C     Test whether the maximum number of iterations have been tried.
C
      M = M + 1
      IF (M .GE. MAXIT) GO TO 380
C
C     Copy the residual to DELTA and its norm to DELNRM, and loop for
C     another iteration.
C
      CALL DCOPY (NEQ, R, 1, DELTA, 1)
      DELNRM = FNRM      
      GO TO 300
C
C     The maximum number of iterations was done.  Set IERNEW and return.
C
 380  IF (RATE .LE. RATEMX) THEN
         IERNEW = 1
      ELSE
         IERNEW = 2
      ENDIF
      RETURN
C
 390  IF (IRES .LE. -2) THEN
         IERNEW = -1
      ELSE
         IERNEW = 3
      ENDIF
      RETURN
C
C
C------END OF SUBROUTINE DNSID------------------------------------------
      END
      SUBROUTINE DLINSD (NEQ, Y, T, YPRIME, CJ, TSCALE, P, PNRM, WT,
     *                   LSOFF, STPTOL, IRET, RES, IRES, WM, IWM,
     *                   FNRM, ICOPT, ID, R, YNEW, YPNEW, ICNFLG,
     *                   ICNSTR, RLX, RPAR, IPAR)
C
C***BEGIN PROLOGUE  DLINSD
C***REFER TO  DNSID
C***DATE WRITTEN   941025   (YYMMDD)
C***REVISION DATE  941215   (YYMMDD)
C***REVISION DATE  960129   Moved line RL = ONE to top block.
C***REVISION DATE  000628   TSCALE argument added.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DLINSD uses a linesearch algorithm to calculate a new (Y,YPRIME)
C     pair (YNEW,YPNEW) such that 
C
C     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME) ,
C
C     where 0 < RL <= 1.  Here, f(y,y') is defined as
C
C      f(y,y') = (1/2)*norm( (J-inverse)*G(t,y,y') )**2 ,
C
C     where norm() is the weighted RMS vector norm, G is the DAE
C     system residual function, and J is the system iteration matrix
C     (Jacobian).
C
C     In addition to the parameters defined elsewhere, we have
C
C     TSCALE  --  Scale factor in T, used for stopping tests if nonzero.
C     P       -- Approximate Newton step used in backtracking.
C     PNRM    -- Weighted RMS norm of P.
C     LSOFF   -- Flag showing whether the linesearch algorithm is
C                to be invoked.  0 means do the linesearch, and
C                1 means turn off linesearch.
C     STPTOL  -- Tolerance used in calculating the minimum lambda
C                value allowed.
C     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations
C                in the proposed new approximate solution will be
C                checked for, and the maximum step length will be
C                adjusted accordingly.
C     ICNSTR  -- Integer array of length NEQ containing flags for
C                checking constraints.
C     RLX     -- Real scalar restricting update size in DCNSTR.
C     YNEW    -- Array of length NEQ used to hold the new Y in
C                performing the linesearch.
C     YPNEW   -- Array of length NEQ used to hold the new YPRIME in
C                performing the linesearch.
C     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW).
C     YPRIME  -- Array of length NEQ containing the new YPRIME 
C                (i.e.,=YPNEW).
C     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the
C                current (Y,YPRIME) on input and output.
C     R       -- Work array of length NEQ, containing the scaled 
C                residual (J-inverse)*G(t,y,y') on return.
C     IRET    -- Return flag.
C                IRET=0 means that a satisfactory (Y,YPRIME) was found.
C                IRET=1 means that the routine failed to find a new
C                       (Y,YPRIME) that was sufficiently distinct from
C                       the current (Y,YPRIME) pair.
C                IRET=2 means IRES .ne. 0 from RES.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DFNRMD, DYYPNW, DCNSTR, DCOPY, XERRWD
C
C***END PROLOGUE  DLINSD
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL  RES
      DIMENSION Y(*), YPRIME(*), WT(*), R(*), ID(*)
      DIMENSION WM(*), IWM(*)
      DIMENSION YNEW(*), YPNEW(*), P(*), ICNSTR(*)
      DIMENSION RPAR(*), IPAR(*)
      CHARACTER MSG*80
C
      PARAMETER (LNRE=12, LKPRIN=31)
C
      SAVE ALPHA, ONE, TWO
      DATA ALPHA/1.0D-4/, ONE/1.0D0/, TWO/2.0D0/
C
      KPRIN=IWM(LKPRIN)
C
      F1NRM = (FNRM*FNRM)/TWO
      RATIO = ONE
      IF (KPRIN .GE. 2) THEN
        MSG = '------ IN ROUTINE DLINSD-- PNRM = (R1)'
        CALL XERRWD(MSG, 38, 901, 0, 0, 0, 0, 1, PNRM, 0.0D0)
        ENDIF
      TAU = PNRM
      RL = ONE
C-----------------------------------------------------------------------
C Check for violations of the constraints, if any are imposed.
C If any violations are found, the step vector P is rescaled, and the 
C constraint check is repeated, until no violations are found.
C-----------------------------------------------------------------------
      IF (ICNFLG .NE. 0) THEN
 10      CONTINUE
         CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
         CALL DCNSTR (NEQ, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
         IF (IRET .EQ. 1) THEN
            RATIO1 = TAU/PNRM
            RATIO = RATIO*RATIO1
            DO 20 I = 1,NEQ
 20           P(I) = P(I)*RATIO1
            PNRM = TAU
            IF (KPRIN .GE. 2) THEN
              MSG = '------ CONSTRAINT VIOL., PNRM = (R1), INDEX = (I1)'
              CALL XERRWD(MSG, 50, 902, 0, 1, IVAR, 0, 1, PNRM, 0.0D0)
              ENDIF
            IF (PNRM .LE. STPTOL) THEN
              IRET = 1
              RETURN
              ENDIF
            GO TO 10
            ENDIF
         ENDIF
C
      SLPI = (-TWO*F1NRM)*RATIO
      RLMIN = STPTOL/PNRM
      IF (LSOFF .EQ. 0 .AND. KPRIN .GE. 2) THEN
        MSG = '------ MIN. LAMBDA = (R1)'
        CALL XERRWD(MSG, 25, 903, 0, 0, 0, 0, 1, RLMIN, 0.0D0)
        ENDIF
C-----------------------------------------------------------------------
C Begin iteration to find RL value satisfying alpha-condition.
C If RL becomes less than RLMIN, then terminate with IRET = 1.
C-----------------------------------------------------------------------
 100  CONTINUE
      CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
      CALL DFNRMD (NEQ, YNEW, T, YPNEW, R, CJ, TSCALE, WT, RES, IRES,
     *              FNRMP, WM, IWM, RPAR, IPAR)
      IWM(LNRE) = IWM(LNRE) + 1
      IF (IRES .NE. 0) THEN
        IRET = 2
        RETURN
        ENDIF
      IF (LSOFF .EQ. 1) GO TO 150
C
      F1NRMP = FNRMP*FNRMP/TWO
      IF (KPRIN .GE. 2) THEN
        MSG = '------ LAMBDA = (R1)'
        CALL XERRWD(MSG, 20, 904, 0, 0, 0, 0, 1, RL, 0.0D0)
        MSG = '------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)'
        CALL XERRWD(MSG, 43, 905, 0, 0, 0, 0, 2, F1NRM, F1NRMP)
        ENDIF
      IF (F1NRMP .GT. F1NRM + ALPHA*SLPI*RL) GO TO 200
C-----------------------------------------------------------------------
C Alpha-condition is satisfied, or linesearch is turned off.
C Copy YNEW,YPNEW to Y,YPRIME and return.
C-----------------------------------------------------------------------
 150  IRET = 0
      CALL DCOPY (NEQ, YNEW, 1, Y, 1)
      CALL DCOPY (NEQ, YPNEW, 1, YPRIME, 1)
      FNRM = FNRMP
      IF (KPRIN .GE. 1) THEN
        MSG = '------ LEAVING ROUTINE DLINSD, FNRM = (R1)'
        CALL XERRWD(MSG, 42, 906, 0, 0, 0, 0, 1, FNRM, 0.0D0)
        ENDIF
      RETURN
C-----------------------------------------------------------------------
C Alpha-condition not satisfied.  Perform backtrack to compute new RL
C value.  If no satisfactory YNEW,YPNEW can be found sufficiently 
C distinct from Y,YPRIME, then return IRET = 1.
C-----------------------------------------------------------------------
 200  CONTINUE
      IF (RL .LT. RLMIN) THEN
        IRET = 1
        RETURN
        ENDIF
C
      RL = RL/TWO
      GO TO 100
C
C----------------------- END OF SUBROUTINE DLINSD ----------------------
      END
      SUBROUTINE DFNRMD (NEQ, Y, T, YPRIME, R, CJ, TSCALE, WT,
     *                   RES, IRES, FNORM, WM, IWM, RPAR, IPAR)
C
C***BEGIN PROLOGUE  DFNRMD
C***REFER TO  DLINSD
C***DATE WRITTEN   941025   (YYMMDD)
C***REVISION DATE  000628   TSCALE argument added.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DFNRMD calculates the scaled preconditioned norm of the nonlinear
C     function used in the nonlinear iteration for obtaining consistent
C     initial conditions.  Specifically, DFNRMD calculates the weighted
C     root-mean-square norm of the vector (J-inverse)*G(T,Y,YPRIME),
C     where J is the Jacobian matrix.
C
C     In addition to the parameters described in the calling program
C     DLINSD, the parameters represent
C
C     R      -- Array of length NEQ that contains
C               (J-inverse)*G(T,Y,YPRIME) on return.
C     TSCALE -- Scale factor in T, used for stopping tests if nonzero.
C     FNORM  -- Scalar containing the weighted norm of R on return.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   RES, DSLVD, DDWNRM
C
C***END PROLOGUE  DFNRMD
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RES
      DIMENSION Y(*), YPRIME(*), WT(*), R(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
C-----------------------------------------------------------------------
C     Call RES routine.
C-----------------------------------------------------------------------
      IRES = 0
      CALL RES(T,Y,YPRIME,CJ,R,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
C-----------------------------------------------------------------------
C     Apply inverse of Jacobian to vector R.
C-----------------------------------------------------------------------
      CALL DSLVD(NEQ,R,WM,IWM)
C-----------------------------------------------------------------------
C     Calculate norm of R.
C-----------------------------------------------------------------------
      FNORM = DDWNRM(NEQ,R,WT,RPAR,IPAR)
      IF (TSCALE .GT. 0.0D0) FNORM = FNORM*TSCALE*ABS(CJ)
C
      RETURN
C----------------------- END OF SUBROUTINE DFNRMD ----------------------
      END
      SUBROUTINE DNEDD(X,Y,YPRIME,NEQ,RES,JACD,PDUM,H,WT,
     *   JSTART,IDID,RPAR,IPAR,PHI,GAMMA,DUMSVR,DELTA,E,
     *   WM,IWM,CJ,CJOLD,CJLAST,S,UROUND,DUME,DUMS,DUMR,
     *   EPCON,JCALC,JFDUM,KP1,NONNEG,NTYPE,IERNLS)
C
C***BEGIN PROLOGUE  DNEDD
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNEDD solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JACD      -- External user-supplied routine to evaluate the
C                  Jacobian.  See JAC description for the case
C                  INFO(12) = 0 in the DDASPK prologue.
C     PDUM      -- Dummy argument.
C     H         -- Appropriate step size for next step.
C     WT        -- Vector of weights for error criterion.
C     JSTART    -- Indicates first call to this routine.
C                  If JSTART = 0, then this is the first call,
C                  otherwise it is not.
C     IDID      -- Completion flag, output by DNEDD.
C                  See IDID description in DDASPK prologue.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     PHI       -- Array of divided differences used by
C                  DNEDD.  The length is NEQ*(K+1),where
C                  K is the maximum order.
C     GAMMA     -- Array used to predict Y and YPRIME.  The length
C                  is MAXORD+1 where MAXORD is the maximum order.
C     DUMSVR    -- Dummy argument.
C     DELTA     -- Work vector for NLS of length NEQ.
C     E         -- Error accumulation vector for NLS of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Parameter always proportional to 1/H.
C     CJOLD     -- Saves the value of CJ as of the last call to DMATD.
C                  Accounts for changes in CJ needed to
C                  decide whether to call DMATD.
C     CJLAST    -- Previous value of CJ.
C     S         -- A scalar determined by the approximate rate
C                  of convergence of the Newton iteration and used
C                  in the convergence test for the Newton iteration.
C
C                  If RATE is defined to be an estimate of the
C                  rate of convergence of the Newton iteration,
C                  then S = RATE/(1.D0-RATE).
C
C                  The closer RATE is to 0., the faster the Newton
C                  iteration is converging; the closer RATE is to 1.,
C                  the slower the Newton iteration is converging.
C
C                  On the first Newton iteration with an up-dated
C                  preconditioner S = 100.D0, Thus the initial
C                  RATE of convergence is approximately 1.
C
C                  S is preserved from call to call so that the rate
C                  estimate from a previous step can be applied to
C                  the current step.
C     UROUND    -- Unit roundoff.
C     DUME      -- Dummy argument.
C     DUMS      -- Dummy argument.
C     DUMR      -- Dummy argument.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     JCALC     -- Flag used to determine when to update
C                  the Jacobian matrix.  In general:
C
C                  JCALC = -1 ==> Call the DMATD routine to update
C                                 the Jacobian matrix.
C                  JCALC =  0 ==> Jacobian matrix is up-to-date.
C                  JCALC =  1 ==> Jacobian matrix is out-dated,
C                                 but DMATD will not be called unless
C                                 JCALC is set to -1.
C     JFDUM     -- Dummy argument.
C     KP1       -- The current order(K) + 1;  updated across calls.
C     NONNEG    -- Flag to determine nonnegativity constraints.
C     NTYPE     -- Identification code for the NLS routine.
C                   0  ==> modified Newton; direct solver.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0  ==> nonlinear solver converged.
C                   1  ==> recoverable error inside nonlinear solver.
C                  -1  ==> unrecoverable error inside nonlinear solver.
C
C     All variables with "DUM" in their names are dummy variables
C     which are not used in this routine.
C
C     Following is a list and description of local variables which
C     may not have an obvious usage.  They are listed in roughly the
C     order they occur in this subroutine.
C
C     The following group of variables are passed as arguments to
C     the Newton iteration solver.  They are explained in greater detail
C     in DNSD:
C        TOLNEW, MULDEL, MAXIT, IERNEW
C
C     IERTYP -- Flag which tells whether this subroutine is correct.
C               0 ==> correct subroutine.
C               1 ==> incorrect subroutine.
C 
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DDWNRM, RES, DMATD, DNSD
C
C***END PROLOGUE  DNEDD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*)
      DIMENSION DELTA(*),E(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION PHI(NEQ,*),GAMMA(*)
      EXTERNAL  RES, JACD
C
      PARAMETER (LNRE=12, LNJE=13)
C
      SAVE MULDEL, MAXIT, XRATE
      DATA MULDEL/1/, MAXIT/4/, XRATE/0.25D0/
C
C     Verify that this is the correct subroutine.
C
      IERTYP = 0
      IF (NTYPE .NE. 0) THEN
         IERTYP = 1
         GO TO 380
         ENDIF
C
C     If this is the first step, perform initializations.
C
      IF (JSTART .EQ. 0) THEN
         CJOLD = CJ
         JCALC = -1
         ENDIF
C
C     Perform all other initializations.
C
      IERNLS = 0
C
C     Decide whether new Jacobian is needed.
C
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) S = 100.D0
C
C-----------------------------------------------------------------------
C     Entry point for updating the Jacobian with current
C     stepsize.
C-----------------------------------------------------------------------
300   CONTINUE
C
C     Initialize all error flags to zero.
C
      IERJ = 0
      IRES = 0
      IERNEW = 0
C
C     Predict the solution and derivative and compute the tolerance
C     for the Newton iteration.
C
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
      PNORM = DDWNRM (NEQ,Y,WT,RPAR,IPAR)
      TOLNEW = 100.D0*UROUND*PNORM
C     
C     Call RES to initialize DELTA.
C
      IWM(LNRE)=IWM(LNRE)+1
      CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
C
C     If indicated, reevaluate the iteration matrix 
C     J = dG/dY + CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0).
C     Set JCALC to 0 as an indicator that this has been done.
C
      IF(JCALC .EQ. -1) THEN
         IWM(LNJE)=IWM(LNJE)+1
         JCALC=0
         CALL DMATD(NEQ,X,Y,YPRIME,DELTA,CJ,H,IERJ,WT,E,WM,IWM,
     *              RES,IRES,UROUND,JACD,RPAR,IPAR)
         CJOLD=CJ
         S = 100.D0
         IF (IRES .LT. 0) GO TO 380
         IF(IERJ .NE. 0)GO TO 380
      ENDIF
C
C     Call the nonlinear Newton solver.
C
      TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
      CALL DNSD(X,Y,YPRIME,NEQ,RES,PDUM,WT,RPAR,IPAR,DUMSVR,
     *          DELTA,E,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON,S,TEMP1,
     *          TOLNEW,MULDEL,MAXIT,IRES,IDUM,IERNEW)
C
      IF (IERNEW .GT. 0 .AND. JCALC .NE. 0) THEN
C
C        The Newton iteration had a recoverable failure with an old
C        iteration matrix.  Retry the step with a new iteration matrix.
C
         JCALC = -1
         GO TO 300
      ENDIF
C
      IF (IERNEW .NE. 0) GO TO 380
C
C     The Newton iteration has converged.  If nonnegativity of
C     solution is required, set the solution nonnegative, if the
C     perturbation to do it is small enough.  If the change is too
C     large, then consider the corrector iteration to have failed.
C
375   IF(NONNEG .EQ. 0) GO TO 390
      DO 377 I = 1,NEQ
377      DELTA(I) = MIN(Y(I),0.0D0)
      DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. EPCON) GO TO 380
      DO 378 I = 1,NEQ
378      E(I) = E(I) - DELTA(I)
      GO TO 390
C
C
C     Exits from nonlinear solver.
C     No convergence with current iteration
C     matrix, or singular iteration matrix.
C     Compute IERNLS and IDID accordingly.
C
380   CONTINUE
      IF (IRES .LE. -2 .OR. IERTYP .NE. 0) THEN
         IERNLS = -1
         IF (IRES .LE. -2) IDID = -11
         IF (IERTYP .NE. 0) IDID = -15
      ELSE
         IERNLS = 1
         IF (IRES .LT. 0) IDID = -10
         IF (IERJ .NE. 0) IDID = -8
      ENDIF
C
390   JCALC = 1
      RETURN
C
C------END OF SUBROUTINE DNEDD------------------------------------------
      END
      SUBROUTINE DNSD(X,Y,YPRIME,NEQ,RES,PDUM,WT,RPAR,IPAR,
     *   DUMSVR,DELTA,E,WM,IWM,CJ,DUMS,DUMR,DUME,EPCON,
     *   S,CONFAC,TOLNEW,MULDEL,MAXIT,IRES,IDUM,IERNEW)
C
C***BEGIN PROLOGUE  DNSD
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  950126   (YYMMDD)
C***REVISION DATE  000711   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSD solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     PDUM      -- Dummy argument.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     DUMSVR    -- Dummy argument.
C     DELTA     -- Work vector for DNSD of length NEQ.
C     E         -- Error accumulation vector for DNSD of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Parameter always proportional to 1/H (step size).
C     DUMS      -- Dummy argument.
C     DUMR      -- Dummy argument.
C     DUME      -- Dummy argument.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     S         -- Used for error convergence tests.
C                  In the Newton iteration: S = RATE/(1 - RATE),
C                  where RATE is the estimated rate of convergence
C                  of the Newton iteration.
C                  The calling routine passes the initial value
C                  of S to the Newton iteration.
C     CONFAC    -- A residual scale factor to improve convergence.
C     TOLNEW    -- Tolerance on the norm of Newton correction in
C                  alternative Newton convergence test.
C     MULDEL    -- A flag indicating whether or not to multiply
C                  DELTA by CONFAC.
C                  0  ==> do not scale DELTA by CONFAC.
C                  1  ==> scale DELTA by CONFAC.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     IRES      -- Error flag returned from RES.  See RES description
C                  in DDASPK prologue.  If IRES = -1, then IERNEW
C                  will be set to 1.
C                  If IRES < -1, then IERNEW will be set to -1.
C     IDUM      -- Dummy argument.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> recoverable error inside Newton iteration.
C                  -1  ==> unrecoverable error inside Newton iteration.
C
C     All arguments with "DUM" in their names are dummy arguments
C     which are not used in this routine.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DSLVD, DDWNRM, RES
C
C***END PROLOGUE  DNSD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),DELTA(*),E(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      EXTERNAL  RES
C
      PARAMETER (LNRE=12, LNNI=19)
C
C     Initialize Newton counter M and accumulation vector E. 
C
      M = 0
      DO 100 I=1,NEQ
100     E(I)=0.0D0
C
C     Corrector loop.
C
300   CONTINUE
      IWM(LNNI) = IWM(LNNI) + 1
C
C     If necessary, multiply residual by convergence factor.
C
      IF (MULDEL .EQ. 1) THEN
         DO 320 I = 1,NEQ
320        DELTA(I) = DELTA(I) * CONFAC
        ENDIF
C
C     Compute a new iterate (back-substitution).
C     Store the correction in DELTA.
C
      CALL DSLVD(NEQ,DELTA,WM,IWM)
C
C     Update Y, E, and YPRIME.
C
      DO 340 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         E(I)=E(I)-DELTA(I)
340      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     Test for convergence of the iteration.
C
      DELNRM=DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (M .EQ. 0) THEN
        OLDNRM = DELNRM
        IF (DELNRM .LE. TOLNEW) GO TO 370
      ELSE
        RATE = (DELNRM/OLDNRM)**(1.0D0/M)
        IF (RATE .GT. 0.9D0) GO TO 380
        S = RATE/(1.0D0 - RATE)
      ENDIF
      IF (S*DELNRM .LE. EPCON) GO TO 370
C
C     The corrector has not yet converged.
C     Update M and test whether the
C     maximum number of iterations have
C     been tried.
C
      M=M+1
      IF(M.GE.MAXIT) GO TO 380
C
C     Evaluate the residual,
C     and go back to do another iteration.
C
      IWM(LNRE)=IWM(LNRE)+1
      CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 300
C
C     The iteration has converged.
C
370   RETURN
C
C     The iteration has not converged.  Set IERNEW appropriately.
C
380   CONTINUE
      IF (IRES .LE. -2 ) THEN
         IERNEW = -1
      ELSE
         IERNEW = 1
      ENDIF
      RETURN
C
C
C------END OF SUBROUTINE DNSD-------------------------------------------
      END
      SUBROUTINE DMATD(NEQ,X,Y,YPRIME,DELTA,CJ,H,IER,EWT,E,
     *                 WM,IWM,RES,IRES,UROUND,JACD,RPAR,IPAR)
C
C***BEGIN PROLOGUE  DMATD
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940701   (YYMMDD) (new LIPVT)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     This routine computes the iteration matrix
C     J = dG/dY+CJ*dG/dYPRIME (where G(X,Y,YPRIME)=0).
C     Here J is computed by:
C       the user-supplied routine JACD if IWM(MTYPE) is 1 or 4, or
C       by numerical difference quotients if IWM(MTYPE) is 2 or 5.
C
C     The parameters have the following meanings.
C     X        = Independent variable.
C     Y        = Array containing predicted values.
C     YPRIME   = Array containing predicted derivatives.
C     DELTA    = Residual evaluated at (X,Y,YPRIME).
C                (Used only if IWM(MTYPE)=2 or 5).
C     CJ       = Scalar parameter defining iteration matrix.
C     H        = Current stepsize in integration.
C     IER      = Variable which is .NE. 0 if iteration matrix
C                is singular, and 0 otherwise.
C     EWT      = Vector of error weights for computing norms.
C     E        = Work space (temporary) of length NEQ.
C     WM       = Real work space for matrices.  On output
C                it contains the LU decomposition
C                of the iteration matrix.
C     IWM      = Integer work space containing
C                matrix information.
C     RES      = External user-supplied subroutine
C                to evaluate the residual.  See RES description
C                in DDASPK prologue.
C     IRES     = Flag which is equal to zero if no illegal values
C                in RES, and less than zero otherwise.  (If IRES
C                is less than zero, the matrix was not completed).
C                In this case (if IRES .LT. 0), then IER = 0.
C     UROUND   = The unit roundoff error of the machine being used.
C     JACD     = Name of the external user-supplied routine
C                to evaluate the iteration matrix.  (This routine
C                is only used if IWM(MTYPE) is 1 or 4)
C                See JAC description for the case INFO(12) = 0
C                in DDASPK prologue.
C     RPAR,IPAR= Real and integer parameter arrays that
C                are used for communication between the
C                calling program and external user routines.
C                They are not altered by DMATD.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   JACD, RES, DGEFA, DGBFA
C
C***END PROLOGUE  DMATD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),DELTA(*),EWT(*),E(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      EXTERNAL  RES, JACD
C
      PARAMETER (LML=1, LMU=2, LMTYPE=4, LNRE=12, LNPD=22, LLCIWP=30)
C
      LIPVT = IWM(LLCIWP)
      IER = 0
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C
C     Dense user-supplied matrix.
C
100   LENPD=IWM(LNPD)
      DO 110 I=1,LENPD
110      WM(I)=0.0D0
      CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR)
      GO TO 230
C
C
C     Dense finite-difference-generated matrix.
C
200   IRES=0
      NROW=0
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),
     *     ABS(1.D0/EWT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         IWM(LNRE)=IWM(LNRE)+1
         CALL RES(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         DO 220 L=1,NEQ
220        WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C
C     Do dense-matrix LU decomposition on J.
C
230      CALL DGEFA(WM,NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
C
C
C     Dummy section for IWM(MTYPE)=3.
C
300   RETURN
C
C
C     Banded user-supplied matrix.
C
400   LENPD=IWM(LNPD)
      DO 410 I=1,LENPD
410      WM(I)=0.0D0
      CALL JACD(X,Y,YPRIME,WM,CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C
C     Banded finite-difference-generated matrix.
C
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN0(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=IWM(LNPD)
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
        DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(1.D0/EWT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
        IWM(LNRE)=IWM(LNRE)+1
        CALL RES(X,Y,YPRIME,CJ,E,IRES,RPAR,IPAR)
        IF (IRES .LT. 0) RETURN
        DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),
     *      ABS(1.D0/EWT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0D0/DEL
          I1=MAX0(1,(N-IWM(LMU)))
          I2=MIN0(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530     CONTINUE
540   CONTINUE
C
C
C     Do LU decomposition of banded J.
C
550   CALL DGBFA (WM,MEBAND,NEQ,IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C
C------END OF SUBROUTINE DMATD------------------------------------------
      END
      SUBROUTINE DSLVD(NEQ,DELTA,WM,IWM)
C
C***BEGIN PROLOGUE  DSLVD
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940701   (YYMMDD) (new LIPVT)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     This routine manages the solution of the linear
C     system arising in the Newton iteration.
C     Real matrix information and real temporary storage
C     is stored in the array WM.
C     Integer matrix information is stored in the array IWM.
C     For a dense matrix, the LINPACK routine DGESL is called.
C     For a banded matrix, the LINPACK routine DGBSL is called.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DGESL, DGBSL
C
C***END PROLOGUE  DSLVD
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DELTA(*),WM(*),IWM(*)
C
      PARAMETER (LML=1, LMU=2, LMTYPE=4, LLCIWP=30)
C
      LIPVT = IWM(LLCIWP)
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
C
C     Dense matrix.
C
100   CALL DGESL(WM,NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
C
C     Dummy section for MTYPE=3.
C
300   CONTINUE
      RETURN
C
C     Banded matrix.
C
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
      CALL DGBSL(WM,MEBAND,NEQ,IWM(LML),
     *  IWM(LMU),IWM(LIPVT),DELTA,0)
      RETURN
C
C------END OF SUBROUTINE DSLVD------------------------------------------
      END
      SUBROUTINE DDASIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,JACK,PSOL,H,TSCALE,
     *   WT,JSKIP,RPAR,IPAR,SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,UROUND,
     *   EPLI,SQRTN,RSQRTN,EPCON,RATEMX,STPTOL,JFLG,
     *   ICNFLG,ICNSTR,IERNLS)
C
C***BEGIN PROLOGUE  DDASIK
C***REFER TO  DDASPK
C***DATE WRITTEN   941026   (YYMMDD)
C***REVISION DATE  950808   (YYMMDD)
C***REVISION DATE  951110   Removed unreachable block 390.
C***REVISION DATE  000628   TSCALE argument added.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C
C     DDASIK solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
C     the initial conditions.
C
C     An initial value for Y and initial guess for YPRIME are input.
C
C     The method used is a Newton scheme with Krylov iteration and a
C     linesearch algorithm.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector at x.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of equations to be integrated.
C     ICOPT     -- Initial condition option chosen (1 or 2).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JACK     --  External user-supplied routine to update
C                  the preconditioner.  (This is optional).
C                  See JAC description for the case
C                  INFO(12) = 1 in the DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning.
C                  (This is optional).  See explanation inside DDASPK.
C     H         -- Scaling factor for this initial condition calc.
C     TSCALE    -- Scale factor in T, used for stopping tests if nonzero.
C     WT        -- Vector of weights for error criterion.
C     JSKIP     -- input flag to signal if initial JAC call is to be
C                  skipped.  1 => skip the call, 0 => do not skip call.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     SAVR      -- Work vector for DDASIK of length NEQ.
C     DELTA     -- Work vector for DDASIK of length NEQ.
C     R         -- Work vector for DDASIK of length NEQ.
C     YIC,YPIC  -- Work vectors for DDASIK, each of length NEQ.
C     PWK       -- Work vector for DDASIK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information for linear system
C                  solvers, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
C     UROUND    -- Unit roundoff.  Not used here.
C     EPLI      -- convergence test constant.
C                  See DDASPK prologue for more details.
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     JFLG      -- Flag showing whether a Jacobian routine is supplied.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length 
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0   ==> nonlinear solver converged.
C                   1,2 ==> recoverable error inside nonlinear solver.
C                           1 => retry with current Y, YPRIME
C                           2 => retry with original Y, YPRIME
C                  -1   ==> unrecoverable error in nonlinear solver.
C
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   RES, JACK, DNSIK, DCOPY
C
C***END PROLOGUE  DDASIK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),ID(*),WT(*),ICNSTR(*)
      DIMENSION SAVR(*),DELTA(*),R(*),YIC(*),YPIC(*),PWK(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      EXTERNAL RES, JACK, PSOL
C
      PARAMETER (LNRE=12, LNJE=13, LLOCWP=29, LLCIWP=30)
      PARAMETER (LMXNIT=32, LMXNJ=33)
C
C
C     Perform initializations.
C
      LWP = IWM(LLOCWP)
      LIWP = IWM(LLCIWP)
      MXNIT = IWM(LMXNIT)
      MXNJ = IWM(LMXNJ)
      IERNLS = 0
      NJ = 0
      EPLIN = EPLI*EPCON
C
C     Call RES to initialize DELTA.
C
      IRES = 0
      IWM(LNRE) = IWM(LNRE) + 1
      CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 370
C
C     Looping point for updating the preconditioner.
C
 300  CONTINUE
C
C     Initialize all error flags to zero.
C
      IERPJ = 0
      IRES = 0
      IERNEW = 0
C
C     If a Jacobian routine was supplied, call it.
C
      IF (JFLG .EQ. 1 .AND. JSKIP .EQ. 0) THEN
        NJ = NJ + 1
        IWM(LNJE)=IWM(LNJE)+1
        CALL JACK (RES, IRES, NEQ, X, Y, YPRIME, WT, DELTA, R, H, CJ,
     *     WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR)
        IF (IRES .LT. 0 .OR. IERPJ .NE. 0) GO TO 370
        ENDIF
      JSKIP = 0
C
C     Call the nonlinear Newton solver for up to MXNIT iterations.
C
      CALL DNSIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,PSOL,WT,RPAR,IPAR,
     *   SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,TSCALE,SQRTN,RSQRTN,
     *   EPLIN,EPCON,RATEMX,MXNIT,STPTOL,ICNFLG,ICNSTR,IERNEW)
C
      IF (IERNEW .EQ. 1 .AND. NJ .LT. MXNJ .AND. JFLG .EQ. 1) THEN
C
C       Up to MXNIT iterations were done, the convergence rate is < 1,
C       a Jacobian routine is supplied, and the number of JACK calls
C       is less than MXNJ.  
C       Copy the residual SAVR to DELTA, call JACK, and try again.
C
        CALL DCOPY (NEQ,  SAVR, 1, DELTA, 1)
        GO TO 300
        ENDIF
C
      IF (IERNEW .NE. 0) GO TO 380
      RETURN
C
C
C     Unsuccessful exits from nonlinear solver.
C     Set IERNLS accordingly.
C
 370  IERNLS = 2
      IF (IRES .LE. -2) IERNLS = -1
      RETURN
C
 380  IERNLS = MIN(IERNEW,2)
      RETURN
C
C----------------------- END OF SUBROUTINE DDASIK-----------------------
      END
      SUBROUTINE DNSIK(X,Y,YPRIME,NEQ,ICOPT,ID,RES,PSOL,WT,RPAR,IPAR,
     *   SAVR,DELTA,R,YIC,YPIC,PWK,WM,IWM,CJ,TSCALE,SQRTN,RSQRTN,EPLIN,
     *   EPCON,RATEMX,MAXIT,STPTOL,ICNFLG,ICNSTR,IERNEW)
C
C***BEGIN PROLOGUE  DNSIK
C***REFER TO  DDASPK
C***DATE WRITTEN   940701   (YYMMDD)
C***REVISION DATE  950714   (YYMMDD)
C***REVISION DATE  000628   TSCALE argument added.
C***REVISION DATE  000628   Added criterion for IERNEW = 1 return.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSIK solves a nonlinear system of algebraic equations of the
C     form G(X,Y,YPRIME) = 0 for the unknown parts of Y and YPRIME in
C     the initial conditions.
C
C     The method used is a Newton scheme combined with a linesearch
C     algorithm, using Krylov iterative linear system methods.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     ICOPT     -- Initial condition option chosen (1 or 2).
C     ID        -- Array of dimension NEQ, which must be initialized
C                  if ICOPT = 1.  See DDASIC.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  See explanation inside DDASPK.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     SAVR      -- Work vector for DNSIK of length NEQ.
C     DELTA     -- Residual vector on entry, and work vector of
C                  length NEQ for DNSIK.
C     R         -- Work vector for DNSIK of length NEQ.
C     YIC,YPIC  -- Work vectors for DNSIK, each of length NEQ.
C     PWK       -- Work vector for DNSIK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Matrix parameter = 1/H (ICOPT = 1) or 0 (ICOPT = 2).
C     TSCALE    -- Scale factor in T, used for stopping tests if nonzero.
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPLIN     -- Tolerance for linear system solver.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     RATEMX    -- Maximum convergence rate for which Newton iteration
C                  is considered converging.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     STPTOL    -- Tolerance used in calculating the minimum lambda
C                  value allowed.
C     ICNFLG    -- Integer scalar.  If nonzero, then constraint
C                  violations in the proposed new approximate solution
C                  will be checked for, and the maximum step length
C                  will be adjusted accordingly.
C     ICNSTR    -- Integer array of length NEQ containing flags for
C                  checking constraints.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> failed to converge, but RATE .lt. 1, or the
C                          residual norm was reduced by a factor of .1.
C                   2  ==> failed to converge, RATE .gt. RATEMX.
C                   3  ==> other recoverable error.
C                  -1  ==> unrecoverable error inside Newton iteration.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DFNRMK, DSLVK, DDWNRM, DLINSK, DCOPY
C
C***END PROLOGUE  DNSIK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),ID(*),DELTA(*),R(*),SAVR(*)
      DIMENSION YIC(*),YPIC(*),PWK(*),WM(*),IWM(*), RPAR(*),IPAR(*)
      DIMENSION ICNSTR(*)
      EXTERNAL RES, PSOL
C
      PARAMETER (LNNI=19, LNPS=21, LLOCWP=29, LLCIWP=30)
      PARAMETER (LLSOFF=35, LSTOL=14)
C
C
C     Initializations.  M is the Newton iteration counter.
C
      LSOFF = IWM(LLSOFF)
      M = 0
      RATE = 1.0D0
      LWP = IWM(LLOCWP)
      LIWP = IWM(LLCIWP)
      RLX = 0.4D0
C
C     Save residual in SAVR.
C
      CALL DCOPY (NEQ, DELTA, 1, SAVR, 1)
C
C     Compute norm of (P-inverse)*(residual).
C
      CALL DFNRMK (NEQ, Y, X, YPRIME, SAVR, R, CJ, TSCALE, WT,
     *   SQRTN, RSQRTN, RES, IRES, PSOL, 1, IER, FNRM, EPLIN,
     *   WM(LWP), IWM(LIWP), PWK, RPAR, IPAR)
      IWM(LNPS) = IWM(LNPS) + 1
      IF (IER .NE. 0) THEN
        IERNEW = 3
        RETURN
      ENDIF
C
C     Return now if residual norm is .le. EPCON.
C
      IF (FNRM .LE. EPCON) RETURN
C
C     Newton iteration loop.
C
      FNRM0 = FNRM
300   CONTINUE
      IWM(LNNI) = IWM(LNNI) + 1
C
C     Compute a new step vector DELTA.
C
      CALL DSLVK (NEQ, Y, X, YPRIME, SAVR, DELTA, WT, WM, IWM,
     *   RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, RHOK,
     *   RPAR, IPAR)
      IF (IRES .NE. 0 .OR. IERSL .NE. 0) GO TO 390
C
C     Get norm of DELTA.  Return now if DELTA is zero.
C
      DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM .EQ. 0.0D0) RETURN
C
C     Call linesearch routine for global strategy and set RATE.
C
      OLDFNM = FNRM
C
      CALL DLINSK (NEQ, Y, X, YPRIME, SAVR, CJ, TSCALE, DELTA, DELNRM,
     *   WT, SQRTN, RSQRTN, LSOFF, STPTOL, IRET, RES, IRES, PSOL,
     *   WM, IWM, RHOK, FNRM, ICOPT, ID, WM(LWP), IWM(LIWP), R, EPLIN,
     *   YIC, YPIC, PWK, ICNFLG, ICNSTR, RLX, RPAR, IPAR)
C
      RATE = FNRM/OLDFNM
C
C     Check for error condition from linesearch.
      IF (IRET .NE. 0) GO TO 390
C
C     Test for convergence of the iteration, and return or loop.
C
      IF (FNRM .LE. EPCON) RETURN
C
C     The iteration has not yet converged.  Update M.
C     Test whether the maximum number of iterations have been tried.
C
      M = M + 1
      IF(M .GE. MAXIT) GO TO 380
C
C     Copy the residual SAVR to DELTA and loop for another iteration.
C
      CALL DCOPY (NEQ,  SAVR, 1, DELTA, 1)
      GO TO 300
C
C     The maximum number of iterations was done.  Set IERNEW and return.
C
380   IF (RATE .LE. RATEMX .OR. FNRM .LE. 0.1D0*FNRM0) THEN
         IERNEW = 1
      ELSE
         IERNEW = 2
      ENDIF
      RETURN
C
390   IF (IRES .LE. -2 .OR. IERSL .LT. 0) THEN
         IERNEW = -1
      ELSE
         IERNEW = 3
         IF (IRES .EQ. 0 .AND. IERSL .EQ. 1 .AND. M .GE. 2 
     1       .AND. RATE .LT. 1.0D0) IERNEW = 1
      ENDIF
      RETURN
C
C
C----------------------- END OF SUBROUTINE DNSIK------------------------
      END
      SUBROUTINE DLINSK (NEQ, Y, T, YPRIME, SAVR, CJ, TSCALE, P, PNRM,
     *   WT, SQRTN, RSQRTN, LSOFF, STPTOL, IRET, RES, IRES, PSOL,
     *   WM, IWM, RHOK, FNRM, ICOPT, ID, WP, IWP, R, EPLIN, YNEW, YPNEW,
     *   PWK, ICNFLG, ICNSTR, RLX, RPAR, IPAR)
C
C***BEGIN PROLOGUE  DLINSK
C***REFER TO  DNSIK
C***DATE WRITTEN   940830   (YYMMDD)
C***REVISION DATE  951006   (Arguments SQRTN, RSQRTN added.)
C***REVISION DATE  960129   Moved line RL = ONE to top block.
C***REVISION DATE  000628   TSCALE argument added.
C***REVISION DATE  000628   RHOK*RHOK term removed in alpha test.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DLINSK uses a linesearch algorithm to calculate a new (Y,YPRIME)
C     pair (YNEW,YPNEW) such that 
C
C     f(YNEW,YPNEW) .le. (1 - 2*ALPHA*RL)*f(Y,YPRIME)
C
C     where 0 < RL <= 1, and RHOK is the scaled preconditioned norm of
C     the final residual vector in the Krylov iteration.  
C     Here, f(y,y') is defined as
C
C      f(y,y') = (1/2)*norm( (P-inverse)*G(t,y,y') )**2 ,
C
C     where norm() is the weighted RMS vector norm, G is the DAE
C     system residual function, and P is the preconditioner used
C     in the Krylov iteration.
C
C     In addition to the parameters defined elsewhere, we have
C
C     SAVR    -- Work array of length NEQ, containing the residual
C                vector G(t,y,y') on return.
C     TSCALE  -- Scale factor in T, used for stopping tests if nonzero.
C     P       -- Approximate Newton step used in backtracking.
C     PNRM    -- Weighted RMS norm of P.
C     LSOFF   -- Flag showing whether the linesearch algorithm is
C                to be invoked.  0 means do the linesearch, 
C                1 means turn off linesearch.
C     STPTOL  -- Tolerance used in calculating the minimum lambda
C                value allowed.
C     ICNFLG  -- Integer scalar.  If nonzero, then constraint violations
C                in the proposed new approximate solution will be
C                checked for, and the maximum step length will be
C                adjusted accordingly.
C     ICNSTR  -- Integer array of length NEQ containing flags for
C                checking constraints.
C     RHOK    -- Weighted norm of preconditioned Krylov residual.
C     RLX     -- Real scalar restricting update size in DCNSTR.
C     YNEW    -- Array of length NEQ used to hold the new Y in
C                performing the linesearch.
C     YPNEW   -- Array of length NEQ used to hold the new YPRIME in
C                performing the linesearch.
C     PWK     -- Work vector of length NEQ for use in PSOL.
C     Y       -- Array of length NEQ containing the new Y (i.e.,=YNEW).
C     YPRIME  -- Array of length NEQ containing the new YPRIME 
C                (i.e.,=YPNEW).
C     FNRM    -- Real scalar containing SQRT(2*f(Y,YPRIME)) for the
C                current (Y,YPRIME) on input and output.
C     R       -- Work space length NEQ for residual vector.
C     IRET    -- Return flag.
C                IRET=0 means that a satisfactory (Y,YPRIME) was found.
C                IRET=1 means that the routine failed to find a new
C                       (Y,YPRIME) that was sufficiently distinct from
C                       the current (Y,YPRIME) pair.
C                IRET=2 means a failure in RES or PSOL.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   DFNRMK, DYYPNW, DCNSTR, DCOPY, XERRWD
C
C***END PROLOGUE  DLINSK
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      EXTERNAL  RES, PSOL
      DIMENSION Y(*), YPRIME(*), P(*), WT(*), SAVR(*), R(*), ID(*)
      DIMENSION WM(*), IWM(*), YNEW(*), YPNEW(*), PWK(*), ICNSTR(*)
      DIMENSION WP(*), IWP(*), RPAR(*), IPAR(*)
      CHARACTER MSG*80
C
      PARAMETER (LNRE=12, LNPS=21, LKPRIN=31)
C
      SAVE ALPHA, ONE, TWO
      DATA ALPHA/1.0D-4/, ONE/1.0D0/, TWO/2.0D0/
C
      KPRIN=IWM(LKPRIN)
      F1NRM = (FNRM*FNRM)/TWO
      RATIO = ONE
C
      IF (KPRIN .GE. 2) THEN
        MSG = '------ IN ROUTINE DLINSK-- PNRM = (R1)'
        CALL XERRWD(MSG, 38, 921, 0, 0, 0, 0, 1, PNRM, 0.0D0)
        ENDIF
      TAU = PNRM
      RL = ONE
C-----------------------------------------------------------------------
C Check for violations of the constraints, if any are imposed.
C If any violations are found, the step vector P is rescaled, and the 
C constraint check is repeated, until no violations are found.
C-----------------------------------------------------------------------
      IF (ICNFLG .NE. 0) THEN
 10      CONTINUE
         CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
         CALL DCNSTR (NEQ, Y, YNEW, ICNSTR, TAU, RLX, IRET, IVAR)
         IF (IRET .EQ. 1) THEN
            RATIO1 = TAU/PNRM
            RATIO = RATIO*RATIO1
            DO 20 I = 1,NEQ
 20           P(I) = P(I)*RATIO1
            PNRM = TAU
            IF (KPRIN .GE. 2) THEN
              MSG = '------ CONSTRAINT VIOL., PNRM = (R1), INDEX = (I1)'
              CALL XERRWD(MSG, 50, 922, 0, 1, IVAR, 0, 1, PNRM, 0.0D0)
              ENDIF
            IF (PNRM .LE. STPTOL) THEN
              IRET = 1
              RETURN
              ENDIF
            GO TO 10
            ENDIF
         ENDIF
C
      SLPI = -TWO*F1NRM*RATIO
      RLMIN = STPTOL/PNRM
      IF (LSOFF .EQ. 0 .AND. KPRIN .GE. 2) THEN
        MSG = '------ MIN. LAMBDA = (R1)'
        CALL XERRWD(MSG, 25, 923, 0, 0, 0, 0, 1, RLMIN, 0.0D0)
        ENDIF
C-----------------------------------------------------------------------
C Begin iteration to find RL value satisfying alpha-condition.
C Update YNEW and YPNEW, then compute norm of new scaled residual and
C perform alpha condition test.
C-----------------------------------------------------------------------
 100  CONTINUE
      CALL DYYPNW (NEQ,Y,YPRIME,CJ,RL,P,ICOPT,ID,YNEW,YPNEW)
      CALL DFNRMK (NEQ, YNEW, T, YPNEW, SAVR, R, CJ, TSCALE, WT,
     *   SQRTN, RSQRTN, RES, IRES, PSOL, 0, IER, FNRMP, EPLIN,
     *   WP, IWP, PWK, RPAR, IPAR)
      IWM(LNRE) = IWM(LNRE) + 1
      IF (IRES .GE. 0) IWM(LNPS) = IWM(LNPS) + 1
      IF (IRES .NE. 0 .OR. IER .NE. 0) THEN
        IRET = 2
        RETURN
        ENDIF
      IF (LSOFF .EQ. 1) GO TO 150
C
      F1NRMP = FNRMP*FNRMP/TWO
      IF (KPRIN .GE. 2) THEN
        MSG = '------ LAMBDA = (R1)'
        CALL XERRWD(MSG, 20, 924, 0, 0, 0, 0, 1, RL, 0.0D0)
        MSG = '------ NORM(F1) = (R1),  NORM(F1NEW) = (R2)'
        CALL XERRWD(MSG, 43, 925, 0, 0, 0, 0, 2, F1NRM, F1NRMP)
        ENDIF
      IF (F1NRMP .GT. F1NRM + ALPHA*SLPI*RL) GO TO 200
C-----------------------------------------------------------------------
C Alpha-condition is satisfied, or linesearch is turned off.
C Copy YNEW,YPNEW to Y,YPRIME and return.
C-----------------------------------------------------------------------
 150  IRET = 0
      CALL DCOPY(NEQ, YNEW, 1, Y, 1)
      CALL DCOPY(NEQ, YPNEW, 1, YPRIME, 1)
      FNRM = FNRMP
      IF (KPRIN .GE. 1) THEN
        MSG = '------ LEAVING ROUTINE DLINSK, FNRM = (R1)'
        CALL XERRWD(MSG, 42, 926, 0, 0, 0, 0, 1, FNRM, 0.0D0)
        ENDIF
      RETURN
C-----------------------------------------------------------------------
C Alpha-condition not satisfied.  Perform backtrack to compute new RL
C value.  If RL is less than RLMIN, i.e. no satisfactory YNEW,YPNEW can
C be found sufficiently distinct from Y,YPRIME, then return IRET = 1.
C-----------------------------------------------------------------------
 200  CONTINUE
      IF (RL .LT. RLMIN) THEN
        IRET = 1
        RETURN
        ENDIF
C
      RL = RL/TWO
      GO TO 100
C
C----------------------- END OF SUBROUTINE DLINSK ----------------------
      END
      SUBROUTINE DFNRMK (NEQ, Y, T, YPRIME, SAVR, R, CJ, TSCALE, WT,
     *                   SQRTN, RSQRTN, RES, IRES, PSOL, IRIN, IER,
     *                   FNORM, EPLIN, WP, IWP, PWK, RPAR, IPAR)
C
C***BEGIN PROLOGUE  DFNRMK
C***REFER TO  DLINSK
C***DATE WRITTEN   940830   (YYMMDD)
C***REVISION DATE  951006   (SQRTN, RSQRTN, and scaling of WT added.)
C***REVISION DATE  000628   TSCALE argument added.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DFNRMK calculates the scaled preconditioned norm of the nonlinear
C     function used in the nonlinear iteration for obtaining consistent
C     initial conditions.  Specifically, DFNRMK calculates the weighted
C     root-mean-square norm of the vector (P-inverse)*G(T,Y,YPRIME),
C     where P is the preconditioner matrix.
C
C     In addition to the parameters described in the calling program
C     DLINSK, the parameters represent
C
C     TSCALE -- Scale factor in T, used for stopping tests if nonzero.
C     IRIN   -- Flag showing whether the current residual vector is
C               input in SAVR.  1 means it is, 0 means it is not.
C     R      -- Array of length NEQ that contains
C               (P-inverse)*G(T,Y,YPRIME) on return.
C     FNORM  -- Scalar containing the weighted norm of R on return.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   RES, DCOPY, DSCAL, PSOL, DDWNRM
C
C***END PROLOGUE  DFNRMK
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL RES, PSOL
      DIMENSION Y(*), YPRIME(*), WT(*), SAVR(*), R(*), PWK(*)
      DIMENSION WP(*), IWP(*), RPAR(*), IPAR(*)
C-----------------------------------------------------------------------
C     Call RES routine if IRIN = 0.
C-----------------------------------------------------------------------
      IF (IRIN .EQ. 0) THEN
        IRES = 0
        CALL RES (T, Y, YPRIME, CJ, SAVR, IRES, RPAR, IPAR)
        IF (IRES .LT. 0) RETURN
        ENDIF
C-----------------------------------------------------------------------
C     Apply inverse of left preconditioner to vector R.
C     First scale WT array by 1/sqrt(N), and undo scaling afterward.
C-----------------------------------------------------------------------
      CALL DCOPY(NEQ, SAVR, 1, R, 1)
      CALL DSCAL (NEQ, RSQRTN, WT, 1)
      IER = 0
      CALL PSOL (NEQ, T, Y, YPRIME, SAVR, PWK, CJ, WT, WP, IWP,
     *           R, EPLIN, IER, RPAR, IPAR)
      CALL DSCAL (NEQ, SQRTN, WT, 1)
      IF (IER .NE. 0) RETURN
C-----------------------------------------------------------------------
C     Calculate norm of R.
C-----------------------------------------------------------------------
      FNORM = DDWNRM (NEQ, R, WT, RPAR, IPAR)
      IF (TSCALE .GT. 0.0D0) FNORM = FNORM*TSCALE*ABS(CJ)
C
      RETURN
C----------------------- END OF SUBROUTINE DFNRMK ----------------------
      END
      SUBROUTINE DNEDK(X,Y,YPRIME,NEQ,RES,JACK,PSOL,
     *   H,WT,JSTART,IDID,RPAR,IPAR,PHI,GAMMA,SAVR,DELTA,E,
     *   WM,IWM,CJ,CJOLD,CJLAST,S,UROUND,EPLI,SQRTN,RSQRTN,
     *   EPCON,JCALC,JFLG,KP1,NONNEG,NTYPE,IERNLS)
C
C***BEGIN PROLOGUE  DNEDK
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940701   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNEDK solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a matrix-free Newton scheme.
C
C     The parameters represent
C     X         -- Independent variable.
C     Y         -- Solution vector at x.
C     YPRIME    -- Derivative of solution vector
C                  after successful step.
C     NEQ       -- Number of equations to be integrated.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     JACK     --  External user-supplied routine to update
C                  the preconditioner.  (This is optional).
C                  See JAC description for the case
C                  INFO(12) = 1 in the DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  (This is optional).  See explanation inside DDASPK.
C     H         -- Appropriate step size for this step.
C     WT        -- Vector of weights for error criterion.
C     JSTART    -- Indicates first call to this routine.
C                  If JSTART = 0, then this is the first call,
C                  otherwise it is not.
C     IDID      -- Completion flag, output by DNEDK.
C                  See IDID description in DDASPK prologue.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     PHI       -- Array of divided differences used by
C                  DNEDK.  The length is NEQ*(K+1), where
C                  K is the maximum order.
C     GAMMA     -- Array used to predict Y and YPRIME.  The length
C                  is K+1, where K is the maximum order.
C     SAVR      -- Work vector for DNEDK of length NEQ.
C     DELTA     -- Work vector for DNEDK of length NEQ.
C     E         -- Error accumulation vector for DNEDK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information for linear system
C                  solvers, and various other information.
C     CJ        -- Parameter always proportional to 1/H.
C     CJOLD     -- Saves the value of CJ as of the last call to DITMD.
C                  Accounts for changes in CJ needed to
C                  decide whether to call DITMD.
C     CJLAST    -- Previous value of CJ.
C     S         -- A scalar determined by the approximate rate
C                  of convergence of the Newton iteration and used
C                  in the convergence test for the Newton iteration.
C
C                  If RATE is defined to be an estimate of the
C                  rate of convergence of the Newton iteration,
C                  then S = RATE/(1.D0-RATE).
C
C                  The closer RATE is to 0., the faster the Newton
C                  iteration is converging; the closer RATE is to 1.,
C                  the slower the Newton iteration is converging.
C
C                  On the first Newton iteration with an up-dated
C                  preconditioner S = 100.D0, Thus the initial
C                  RATE of convergence is approximately 1.
C
C                  S is preserved from call to call so that the rate
C                  estimate from a previous step can be applied to
C                  the current step.
C     UROUND    -- Unit roundoff.  Not used here.
C     EPLI      -- convergence test constant.
C                  See DDASPK prologue for more details.
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     JCALC     -- Flag used to determine when to update
C                  the Jacobian matrix.  In general:
C
C                  JCALC = -1 ==> Call the DITMD routine to update
C                                 the Jacobian matrix.
C                  JCALC =  0 ==> Jacobian matrix is up-to-date.
C                  JCALC =  1 ==> Jacobian matrix is out-dated,
C                                 but DITMD will not be called unless
C                                 JCALC is set to -1.
C     JFLG      -- Flag showing whether a Jacobian routine is supplied.
C     KP1       -- The current order + 1;  updated across calls.
C     NONNEG    -- Flag to determine nonnegativity constraints.
C     NTYPE     -- Identification code for the DNEDK routine.
C                   1 ==> modified Newton; iterative linear solver.
C                   2 ==> modified Newton; user-supplied linear solver.
C     IERNLS    -- Error flag for nonlinear solver.
C                   0 ==> nonlinear solver converged.
C                   1 ==> recoverable error inside non-linear solver.
C                  -1 ==> unrecoverable error inside non-linear solver.
C
C     The following group of variables are passed as arguments to
C     the Newton iteration solver.  They are explained in greater detail
C     in DNSK:
C        TOLNEW, MULDEL, MAXIT, IERNEW
C
C     IERTYP -- Flag which tells whether this subroutine is correct.
C               0 ==> correct subroutine.
C               1 ==> incorrect subroutine.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   RES, JACK, DDWNRM, DNSK
C
C***END PROLOGUE  DNEDK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*)
      DIMENSION PHI(NEQ,*),SAVR(*),DELTA(*),E(*)
      DIMENSION WM(*),IWM(*)
      DIMENSION GAMMA(*),RPAR(*),IPAR(*)
      EXTERNAL  RES, JACK, PSOL
C
      PARAMETER (LNRE=12, LNJE=13, LLOCWP=29, LLCIWP=30)
C
      SAVE MULDEL, MAXIT, XRATE
      DATA MULDEL/0/, MAXIT/4/, XRATE/0.25D0/
C
C     Verify that this is the correct subroutine.
C
      IERTYP = 0
      IF (NTYPE .NE. 1) THEN
         IERTYP = 1
         GO TO 380
         ENDIF
C
C     If this is the first step, perform initializations.
C
      IF (JSTART .EQ. 0) THEN
         CJOLD = CJ
         JCALC = -1
         S = 100.D0
         ENDIF
C
C     Perform all other initializations.
C
      IERNLS = 0
      LWP = IWM(LLOCWP)
      LIWP = IWM(LLCIWP)
C
C     Decide whether to update the preconditioner.
C
      IF (JFLG .NE. 0) THEN
         TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
         TEMP2 = 1.0D0/TEMP1
         IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
         IF (CJ .NE. CJLAST) S = 100.D0
      ELSE
         JCALC = 0
         ENDIF
C
C     Looping point for updating preconditioner with current stepsize.
C
300   CONTINUE
C
C     Initialize all error flags to zero.
C
      IERPJ = 0
      IRES = 0
      IERSL = 0
      IERNEW = 0
C
C     Predict the solution and derivative and compute the tolerance
C     for the Newton iteration.
C
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
      EPLIN = EPLI*EPCON
      TOLNEW = EPLIN
C
C     Call RES to initialize DELTA.
C
      IWM(LNRE)=IWM(LNRE)+1
      CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
C
C
C     If indicated, update the preconditioner.
C     Set JCALC to 0 as an indicator that this has been done.
C
      IF(JCALC .EQ. -1)THEN
         IWM(LNJE) = IWM(LNJE) + 1
         JCALC=0
         CALL JACK (RES, IRES, NEQ, X, Y, YPRIME, WT, DELTA, E, H, CJ,
     *      WM(LWP), IWM(LIWP), IERPJ, RPAR, IPAR)
         CJOLD=CJ
         S = 100.D0
         IF (IRES .LT. 0)  GO TO 380
         IF (IERPJ .NE. 0) GO TO 380
      ENDIF
C
C     Call the nonlinear Newton solver.
C
      CALL DNSK(X,Y,YPRIME,NEQ,RES,PSOL,WT,RPAR,IPAR,SAVR,
     *   DELTA,E,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON,
     *   S,TEMP1,TOLNEW,MULDEL,MAXIT,IRES,IERSL,IERNEW)
C
      IF (IERNEW .GT. 0 .AND. JCALC .NE. 0) THEN
C
C     The Newton iteration had a recoverable failure with an old
C     preconditioner.  Retry the step with a new preconditioner.
C
         JCALC = -1
         GO TO 300
      ENDIF
C
      IF (IERNEW .NE. 0) GO TO 380
C
C     The Newton iteration has converged.  If nonnegativity of
C     solution is required, set the solution nonnegative, if the
C     perturbation to do it is small enough.  If the change is too
C     large, then consider the corrector iteration to have failed.
C
      IF(NONNEG .EQ. 0) GO TO 390
      DO 360 I = 1,NEQ
 360    DELTA(I) = MIN(Y(I),0.0D0)
      DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. EPCON) GO TO 380
      DO 370 I = 1,NEQ
 370    E(I) = E(I) - DELTA(I)
      GO TO 390
C
C
C     Exits from nonlinear solver.
C     No convergence with current preconditioner.
C     Compute IERNLS and IDID accordingly.
C
380   CONTINUE
      IF (IRES .LE. -2 .OR. IERSL .LT. 0 .OR. IERTYP .NE. 0) THEN
         IERNLS = -1
         IF (IRES .LE. -2) IDID = -11
         IF (IERSL .LT. 0) IDID = -13
         IF (IERTYP .NE. 0) IDID = -15
      ELSE
         IERNLS =  1
         IF (IRES .EQ. -1) IDID = -10
         IF (IERPJ .NE. 0) IDID = -5
         IF (IERSL .GT. 0) IDID = -14
      ENDIF
C
C
390   JCALC = 1
      RETURN
C
C------END OF SUBROUTINE DNEDK------------------------------------------
      END
      SUBROUTINE DNSK(X,Y,YPRIME,NEQ,RES,PSOL,WT,RPAR,IPAR,
     *   SAVR,DELTA,E,WM,IWM,CJ,SQRTN,RSQRTN,EPLIN,EPCON,
     *   S,CONFAC,TOLNEW,MULDEL,MAXIT,IRES,IERSL,IERNEW)
C
C***BEGIN PROLOGUE  DNSK
C***REFER TO  DDASPK
C***DATE WRITTEN   891219   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  950126   (YYMMDD)
C***REVISION DATE  000711   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     DNSK solves a nonlinear system of
C     algebraic equations of the form
C     G(X,Y,YPRIME) = 0 for the unknown Y.
C
C     The method used is a modified Newton scheme.
C
C     The parameters represent
C
C     X         -- Independent variable.
C     Y         -- Solution vector.
C     YPRIME    -- Derivative of solution vector.
C     NEQ       -- Number of unknowns.
C     RES       -- External user-supplied subroutine
C                  to evaluate the residual.  See RES description
C                  in DDASPK prologue.
C     PSOL      -- External user-supplied routine to solve
C                  a linear system using preconditioning. 
C                  See explanation inside DDASPK.
C     WT        -- Vector of weights for error criterion.
C     RPAR,IPAR -- Real and integer arrays used for communication
C                  between the calling program and external user
C                  routines.  They are not altered within DASPK.
C     SAVR      -- Work vector for DNSK of length NEQ.
C     DELTA     -- Work vector for DNSK of length NEQ.
C     E         -- Error accumulation vector for DNSK of length NEQ.
C     WM,IWM    -- Real and integer arrays storing
C                  matrix information such as the matrix
C                  of partial derivatives, permutation
C                  vector, and various other information.
C     CJ        -- Parameter always proportional to 1/H (step size).
C     SQRTN     -- Square root of NEQ.
C     RSQRTN    -- reciprical of square root of NEQ.
C     EPLIN     -- Tolerance for linear system solver.
C     EPCON     -- Tolerance to test for convergence of the Newton
C                  iteration.
C     S         -- Used for error convergence tests.
C                  In the Newton iteration: S = RATE/(1.D0-RATE),
C                  where RATE is the estimated rate of convergence
C                  of the Newton iteration.
C
C                  The closer RATE is to 0., the faster the Newton
C                  iteration is converging; the closer RATE is to 1.,
C                  the slower the Newton iteration is converging.
C
C                  The calling routine sends the initial value
C                  of S to the Newton iteration.
C     CONFAC    -- A residual scale factor to improve convergence.
C     TOLNEW    -- Tolerance on the norm of Newton correction in
C                  alternative Newton convergence test.
C     MULDEL    -- A flag indicating whether or not to multiply
C                  DELTA by CONFAC.
C                  0  ==> do not scale DELTA by CONFAC.
C                  1  ==> scale DELTA by CONFAC.
C     MAXIT     -- Maximum allowed number of Newton iterations.
C     IRES      -- Error flag returned from RES.  See RES description
C                  in DDASPK prologue.  If IRES = -1, then IERNEW
C                  will be set to 1.
C                  If IRES < -1, then IERNEW will be set to -1.
C     IERSL     -- Error flag for linear system solver.
C                  See IERSL description in subroutine DSLVK.
C                  If IERSL = 1, then IERNEW will be set to 1.
C                  If IERSL < 0, then IERNEW will be set to -1.
C     IERNEW    -- Error flag for Newton iteration.
C                   0  ==> Newton iteration converged.
C                   1  ==> recoverable error inside Newton iteration.
C                  -1  ==> unrecoverable error inside Newton iteration.
C-----------------------------------------------------------------------
C
C***ROUTINES CALLED
C   RES, DSLVK, DDWNRM
C
C***END PROLOGUE  DNSK
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION Y(*),YPRIME(*),WT(*),DELTA(*),E(*),SAVR(*)
      DIMENSION WM(*),IWM(*), RPAR(*),IPAR(*)
      EXTERNAL  RES, PSOL
C
      PARAMETER (LNNI=19, LNRE=12)
C
C     Initialize Newton counter M and accumulation vector E.
C
      M = 0
      DO 100 I=1,NEQ
100     E(I) = 0.0D0
C
C     Corrector loop.
C
300   CONTINUE
      IWM(LNNI) = IWM(LNNI) + 1
C
C     If necessary, multiply residual by convergence factor.
C
      IF (MULDEL .EQ. 1) THEN
        DO 320 I = 1,NEQ
320       DELTA(I) = DELTA(I) * CONFAC
        ENDIF
C
C     Save residual in SAVR.
C
      DO 340 I = 1,NEQ
340     SAVR(I) = DELTA(I)
C
C     Compute a new iterate.  Store the correction in DELTA.
C
      CALL DSLVK (NEQ, Y, X, YPRIME, SAVR, DELTA, WT, WM, IWM,
     *   RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, RHOK,
     *   RPAR, IPAR)
      IF (IRES .NE. 0 .OR. IERSL .NE. 0) GO TO 380
C
C     Update Y, E, and YPRIME.
C
      DO 360 I=1,NEQ
         Y(I) = Y(I) - DELTA(I)
         E(I) = E(I) - DELTA(I)
360      YPRIME(I) = YPRIME(I) - CJ*DELTA(I)
C
C     Test for convergence of the iteration.
C
      DELNRM = DDWNRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (M .EQ. 0) THEN
        OLDNRM = DELNRM
        IF (DELNRM .LE. TOLNEW) GO TO 370
      ELSE
        RATE = (DELNRM/OLDNRM)**(1.0D0/M)
        IF (RATE .GT. 0.9D0) GO TO 380
        S = RATE/(1.0D0 - RATE)
      ENDIF
      IF (S*DELNRM .LE. EPCON) GO TO 370
C
C     The corrector has not yet converged.  Update M and test whether
C     the maximum number of iterations have been tried.
C
      M = M + 1
      IF (M .GE. MAXIT) GO TO 380
C
C     Evaluate the residual, and go back to do another iteration.
C
      IWM(LNRE) = IWM(LNRE) + 1
      CALL RES(X,Y,YPRIME,CJ,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 300
C
C     The iteration has converged.
C
370    RETURN
C
C     The iteration has not converged.  Set IERNEW appropriately.
C
380   CONTINUE
      IF (IRES .LE. -2 .OR. IERSL .LT. 0) THEN
         IERNEW = -1
      ELSE
         IERNEW = 1
      ENDIF
      RETURN
C
C
C------END OF SUBROUTINE DNSK-------------------------------------------
      END
      SUBROUTINE DSLVK (NEQ, Y, TN, YPRIME, SAVR, X, EWT, WM, IWM,
     *   RES, IRES, PSOL, IERSL, CJ, EPLIN, SQRTN, RSQRTN, RHOK,
     *   RPAR, IPAR)
C
C***BEGIN PROLOGUE  DSLVK
C***REFER TO  DDASPK
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940928   Removed MNEWT and added RHOK in call list.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C DSLVK uses a restart algorithm and interfaces to DSPIGM for
C the solution of the linear system arising from a Newton iteration.
C
C In addition to variables described elsewhere,
C communication with DSLVK uses the following variables..
C WM    = Real work space containing data for the algorithm
C         (Krylov basis vectors, Hessenberg matrix, etc.).
C IWM   = Integer work space containing data for the algorithm.
C X     = The right-hand side vector on input, and the solution vector
C         on output, of length NEQ.
C IRES  = Error flag from RES.
C IERSL = Output flag ..
C         IERSL =  0 means no trouble occurred (or user RES routine
C                    returned IRES < 0)
C         IERSL =  1 means the iterative method failed to converge
C                    (DSPIGM returned IFLAG > 0.)
C         IERSL = -1 means there was a nonrecoverable error in the
C                    iterative solver, and an error exit will occur.
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DSCAL, DCOPY, DSPIGM
C
C***END PROLOGUE  DSLVK
C
      INTEGER NEQ, IWM, IRES, IERSL, IPAR
      DOUBLE PRECISION Y, TN, YPRIME, SAVR, X, EWT, WM, CJ, EPLIN,
     1   SQRTN, RSQRTN, RHOK, RPAR
      DIMENSION Y(*), YPRIME(*), SAVR(*), X(*), EWT(*), 
     1  WM(*), IWM(*), RPAR(*), IPAR(*)
C
      INTEGER IFLAG, IRST, NRSTS, NRMAX, LR, LDL, LHES, LGMR, LQ, LV,
     1        LWK, LZ, MAXLP1, NPSL
      INTEGER NLI, NPS, NCFL, NRE, MAXL, KMP, MITER
      EXTERNAL  RES, PSOL
C    
      PARAMETER (LNRE=12, LNCFL=16, LNLI=20, LNPS=21) 
      PARAMETER (LLOCWP=29, LLCIWP=30)
      PARAMETER (LMITER=23, LMAXL=24, LKMP=25, LNRMAX=26)
C
C-----------------------------------------------------------------------
C IRST is set to 1, to indicate restarting is in effect.
C NRMAX is the maximum number of restarts.
C-----------------------------------------------------------------------
      DATA IRST/1/
C
      LIWP = IWM(LLCIWP)
      NLI = IWM(LNLI)
      NPS = IWM(LNPS)
      NCFL = IWM(LNCFL)
      NRE = IWM(LNRE)
      LWP = IWM(LLOCWP)
      MAXL = IWM(LMAXL) 
      KMP = IWM(LKMP)
      NRMAX = IWM(LNRMAX) 
      MITER = IWM(LMITER)
      IERSL = 0
      IRES = 0
C-----------------------------------------------------------------------
C Use a restarting strategy to solve the linear system
C P*X = -F.  Parse the work vector, and perform initializations.
C Note that zero is the initial guess for X.
C-----------------------------------------------------------------------
      MAXLP1 = MAXL + 1
      LV = 1
      LR = LV + NEQ*MAXL
      LHES = LR + NEQ + 1
      LQ = LHES + MAXL*MAXLP1
      LWK = LQ + 2*MAXL
      LDL = LWK + MIN0(1,MAXL-KMP)*NEQ
      LZ = LDL + NEQ
      CALL DSCAL (NEQ, RSQRTN, EWT, 1)
      CALL DCOPY (NEQ, X, 1, WM(LR), 1)
      DO 110 I = 1,NEQ
 110     X(I) = 0.D0
C-----------------------------------------------------------------------
C Top of loop for the restart algorithm.  Initial pass approximates
C X and sets up a transformed system to perform subsequent restarts
C to update X.  NRSTS is initialized to -1, because restarting
C does not occur until after the first pass.
C Update NRSTS; conditionally copy DL to R; call the DSPIGM
C algorithm to solve A*Z = R;  updated counters;  update X with
C the residual solution.
C Note:  if convergence is not achieved after NRMAX restarts,
C then the linear solver is considered to have failed.
C-----------------------------------------------------------------------
      NRSTS = -1
 115  CONTINUE
      NRSTS = NRSTS + 1
      IF (NRSTS .GT. 0) CALL DCOPY (NEQ, WM(LDL), 1, WM(LR),1)
      CALL DSPIGM (NEQ, TN, Y, YPRIME, SAVR, WM(LR), EWT, MAXL, MAXLP1,
     1   KMP, EPLIN, CJ, RES, IRES, NRES, PSOL, NPSL, WM(LZ), WM(LV),
     2   WM(LHES), WM(LQ), LGMR, WM(LWP), IWM(LIWP), WM(LWK),
     3   WM(LDL), RHOK, IFLAG, IRST, NRSTS, RPAR, IPAR)
      NLI = NLI + LGMR
      NPS = NPS + NPSL
      NRE = NRE + NRES
      DO 120 I = 1,NEQ
 120     X(I) = X(I) + WM(LZ+I-1) 
      IF ((IFLAG .EQ. 1) .AND. (NRSTS .LT. NRMAX) .AND. (IRES .EQ. 0))
     1   GO TO 115
C-----------------------------------------------------------------------
C The restart scheme is finished.  Test IRES and IFLAG to see if
C convergence was not achieved, and set flags accordingly.
C-----------------------------------------------------------------------
      IF (IRES .LT. 0) THEN
         NCFL = NCFL + 1
      ELSE IF (IFLAG .NE. 0) THEN
         NCFL = NCFL + 1
         IF (IFLAG .GT. 0) IERSL = 1 
         IF (IFLAG .LT. 0) IERSL = -1 
      ENDIF
C-----------------------------------------------------------------------
C Update IWM with counters, rescale EWT, and return.
C-----------------------------------------------------------------------
      IWM(LNLI)  = NLI
      IWM(LNPS)  = NPS
      IWM(LNCFL) = NCFL
      IWM(LNRE)  = NRE
      CALL DSCAL (NEQ, SQRTN, EWT, 1)
      RETURN
C
C------END OF SUBROUTINE DSLVK------------------------------------------
      END
      SUBROUTINE DSPIGM (NEQ, TN, Y, YPRIME, SAVR, R, WGHT, MAXL,
     *   MAXLP1, KMP, EPLIN, CJ, RES, IRES, NRE, PSOL, NPSL, Z, V,
     *   HES, Q, LGMR, WP, IWP, WK, DL, RHOK, IFLAG, IRST, NRSTS,
     *   RPAR, IPAR)
C
C***BEGIN PROLOGUE  DSPIGM
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C***REVISION DATE  940927   Removed MNEWT and added RHOK in call list.
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine solves the linear system A * Z = R using a scaled
C preconditioned version of the generalized minimum residual method.
C An initial guess of Z = 0 is assumed.
C
C      On entry
C
C          NEQ = Problem size, passed to PSOL.
C
C           TN = Current Value of T.
C
C            Y = Array Containing current dependent variable vector.
C
C       YPRIME = Array Containing current first derivative of Y.
C
C         SAVR = Array containing current value of G(T,Y,YPRIME).
C
C            R = The right hand side of the system A*Z = R.
C                R is also used as work space when computing
C                the final approximation and will therefore be
C                destroyed.
C                (R is the same as V(*,MAXL+1) in the call to DSPIGM.)
C
C         WGHT = The vector of length NEQ containing the nonzero
C                elements of the diagonal scaling matrix.
C
C         MAXL = The maximum allowable order of the matrix H.
C
C       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES.
C
C          KMP = The number of previous vectors the new vector, VNEW,
C                must be made orthogonal to.  (KMP .LE. MAXL.)
C
C        EPLIN = Tolerance on residuals R-A*Z in weighted rms norm.
C
C           CJ = Scalar proportional to current value of 
C                1/(step size H).
C
C           WK = Real work array used by routine DATV and PSOL.
C
C           DL = Real work array used for calculation of the residual
C                norm RHO when the method is incomplete (KMP.LT.MAXL)
C                and/or when using restarting.
C
C           WP = Real work array used by preconditioner PSOL.
C
C          IWP = Integer work array used by preconditioner PSOL.
C
C         IRST = Method flag indicating if restarting is being
C                performed.  IRST .GT. 0 means restarting is active,
C                while IRST = 0 means restarting is not being used.
C
C        NRSTS = Counter for the number of restarts on the current
C                call to DSPIGM.  If NRSTS .GT. 0, then the residual
C                R is already scaled, and so scaling of R is not
C                necessary.
C
C
C      On Return
C
C         Z    = The final computed approximation to the solution
C                of the system A*Z = R.
C
C         LGMR = The number of iterations performed and
C                the current order of the upper Hessenberg
C                matrix HES.
C
C         NRE  = The number of calls to RES (i.e. DATV)
C
C         NPSL = The number of calls to PSOL.
C
C         V    = The neq by (LGMR+1) array containing the LGMR
C                orthogonal vectors V(*,1) to V(*,LGMR).
C
C         HES  = The upper triangular factor of the QR decomposition
C                of the (LGMR+1) by LGMR upper Hessenberg matrix whose
C                entries are the scaled inner-products of A*V(*,I)
C                and V(*,K).
C
C         Q    = Real array of length 2*MAXL containing the components
C                of the givens rotations used in the QR decomposition
C                of HES.  It is loaded in DHEQR and used in DHELS.
C
C         IRES = Error flag from RES.
C
C           DL = Scaled preconditioned residual, 
C                (D-inverse)*(P-inverse)*(R-A*Z). Only loaded when
C                performing restarts of the Krylov iteration.
C
C         RHOK = Weighted norm of final preconditioned residual.
C
C        IFLAG = Integer error flag..
C                0 Means convergence in LGMR iterations, LGMR.LE.MAXL.
C                1 Means the convergence test did not pass in MAXL
C                  iterations, but the new residual norm (RHO) is
C                  .LT. the old residual norm (RNRM), and so Z is
C                  computed.
C                2 Means the convergence test did not pass in MAXL
C                  iterations, new residual norm (RHO) .GE. old residual
C                  norm (RNRM), and the initial guess, Z = 0, is
C                  returned.
C                3 Means there was a recoverable error in PSOL
C                  caused by the preconditioner being out of date.
C               -1 Means there was an unrecoverable error in PSOL.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   PSOL, DNRM2, DSCAL, DATV, DORTH, DHEQR, DCOPY, DHELS, DAXPY
C
C***END PROLOGUE  DSPIGM
C
      INTEGER NEQ,MAXL,MAXLP1,KMP,IRES,NRE,NPSL,LGMR,IWP,
     1   IFLAG,IRST,NRSTS,IPAR
      DOUBLE PRECISION TN,Y,YPRIME,SAVR,R,WGHT,EPLIN,CJ,Z,V,HES,Q,WP,WK,
     1   DL,RHOK,RPAR
      DIMENSION Y(*), YPRIME(*), SAVR(*), R(*), WGHT(*), Z(*),
     1   V(NEQ,*), HES(MAXLP1,*), Q(*), WP(*), IWP(*), WK(*), DL(*),
     2   RPAR(*), IPAR(*)
      INTEGER I, IER, INFO, IP1, I2, J, K, LL, LLP1
      DOUBLE PRECISION RNRM,C,DLNRM,PROD,RHO,S,SNORMW,DNRM2,TEM
      EXTERNAL  RES, PSOL
C
      IER = 0
      IFLAG = 0
      LGMR = 0
      NPSL = 0
      NRE = 0
C-----------------------------------------------------------------------
C The initial guess for Z is 0.  The initial residual is therefore
C the vector R.  Initialize Z to 0.
C-----------------------------------------------------------------------
      DO 10 I = 1,NEQ
 10     Z(I) = 0.0D0
C-----------------------------------------------------------------------
C Apply inverse of left preconditioner to vector R if NRSTS .EQ. 0.
C Form V(*,1), the scaled preconditioned right hand side.
C-----------------------------------------------------------------------
      IF (NRSTS .EQ. 0) THEN
         CALL PSOL (NEQ, TN, Y, YPRIME, SAVR, WK, CJ, WGHT, WP, IWP,
     1      R, EPLIN, IER, RPAR, IPAR)
         NPSL = 1
         IF (IER .NE. 0) GO TO 300
         DO 30 I = 1,NEQ
 30         V(I,1) = R(I)*WGHT(I)
      ELSE
         DO 35 I = 1,NEQ
 35         V(I,1) = R(I)
      ENDIF
C-----------------------------------------------------------------------
C Calculate norm of scaled vector V(*,1) and normalize it
C If, however, the norm of V(*,1) (i.e. the norm of the preconditioned
C residual) is .le. EPLIN, then return with Z=0.
C-----------------------------------------------------------------------
      RNRM = DNRM2 (NEQ, V, 1)
      IF (RNRM .LE. EPLIN) THEN
        RHOK = RNRM
        RETURN
        ENDIF
      TEM = 1.0D0/RNRM
      CALL DSCAL (NEQ, TEM, V(1,1), 1)
C-----------------------------------------------------------------------
C Zero out the HES array.
C-----------------------------------------------------------------------
      DO 65 J = 1,MAXL
        DO 60 I = 1,MAXLP1
 60       HES(I,J) = 0.0D0
 65     CONTINUE
C-----------------------------------------------------------------------
C Main loop to compute the vectors V(*,2) to V(*,MAXL).
C The running product PROD is needed for the convergence test.
C-----------------------------------------------------------------------
      PROD = 1.0D0
      DO 90 LL = 1,MAXL
        LGMR = LL
C-----------------------------------------------------------------------
C Call routine DATV to compute VNEW = ABAR*V(LL), where ABAR is
C the matrix A with scaling and inverse preconditioner factors applied.
C Call routine DORTH to orthogonalize the new vector VNEW = V(*,LL+1).
C call routine DHEQR to update the factors of HES.
C-----------------------------------------------------------------------
        CALL DATV (NEQ, Y, TN, YPRIME, SAVR, V(1,LL), WGHT, Z,
     1     RES, IRES, PSOL, V(1,LL+1), WK, WP, IWP, CJ, EPLIN,
     1     IER, NRE, NPSL, RPAR, IPAR)
        IF (IRES .LT. 0) RETURN
        IF (IER .NE. 0) GO TO 300
        CALL DORTH (V(1,LL+1), V, HES, NEQ, LL, MAXLP1, KMP, SNORMW)
        HES(LL+1,LL) = SNORMW
        CALL DHEQR (HES, MAXLP1, LL, Q, INFO, LL)
        IF (INFO .EQ. LL) GO TO 120
C-----------------------------------------------------------------------
C Update RHO, the estimate of the norm of the residual R - A*ZL.
C If KMP .LT. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not
C necessarily orthogonal for LL .GT. KMP.  The vector DL must then
C be computed, and its norm used in the calculation of RHO.
C-----------------------------------------------------------------------
        PROD = PROD*Q(2*LL)
        RHO = ABS(PROD*RNRM)
        IF ((LL.GT.KMP) .AND. (KMP.LT.MAXL)) THEN
          IF (LL .EQ. KMP+1) THEN
            CALL DCOPY (NEQ, V(1,1), 1, DL, 1)
            DO 75 I = 1,KMP
              IP1 = I + 1
              I2 = I*2
              S = Q(I2)
              C = Q(I2-1)
              DO 70 K = 1,NEQ
 70             DL(K) = S*DL(K) + C*V(K,IP1)
 75           CONTINUE
            ENDIF
          S = Q(2*LL)
          C = Q(2*LL-1)/SNORMW
          LLP1 = LL + 1
          DO 80 K = 1,NEQ
 80         DL(K) = S*DL(K) + C*V(K,LLP1)
          DLNRM = DNRM2 (NEQ, DL, 1)
          RHO = RHO*DLNRM
          ENDIF
C-----------------------------------------------------------------------
C Test for convergence.  If passed, compute approximation ZL.
C If failed and LL .LT. MAXL, then continue iterating.
C-----------------------------------------------------------------------
        IF (RHO .LE. EPLIN) GO TO 200
        IF (LL .EQ. MAXL) GO TO 100
C-----------------------------------------------------------------------
C Rescale so that the norm of V(1,LL+1) is one.
C-----------------------------------------------------------------------
        TEM = 1.0D0/SNORMW
        CALL DSCAL (NEQ, TEM, V(1,LL+1), 1)
 90     CONTINUE
 100  CONTINUE
      IF (RHO .LT. RNRM) GO TO 150
 120  CONTINUE
      IFLAG = 2
      DO 130 I = 1,NEQ
 130     Z(I) = 0.D0
      RETURN
 150  IFLAG = 1
C-----------------------------------------------------------------------
C The tolerance was not met, but the residual norm was reduced.
C If performing restarting (IRST .gt. 0) calculate the residual vector
C RL and store it in the DL array.  If the incomplete version is 
C being used (KMP .lt. MAXL) then DL has already been calculated.
C-----------------------------------------------------------------------
      IF (IRST .GT. 0) THEN
         IF (KMP .EQ. MAXL) THEN
C
C           Calculate DL from the V(I)'s.
C
            CALL DCOPY (NEQ, V(1,1), 1, DL, 1)
            MAXLM1 = MAXL - 1
            DO 175 I = 1,MAXLM1
               IP1 = I + 1
               I2 = I*2
               S = Q(I2)
               C = Q(I2-1)
               DO 170 K = 1,NEQ
 170              DL(K) = S*DL(K) + C*V(K,IP1)
 175        CONTINUE
            S = Q(2*MAXL)
            C = Q(2*MAXL-1)/SNORMW
            DO 180 K = 1,NEQ
 180           DL(K) = S*DL(K) + C*V(K,MAXLP1)
         ENDIF
C
C        Scale DL by RNRM*PROD to obtain the residual RL.
C
         TEM = RNRM*PROD
         CALL DSCAL(NEQ, TEM, DL, 1)
      ENDIF
C-----------------------------------------------------------------------
C Compute the approximation ZL to the solution.
C Since the vector Z was used as work space, and the initial guess
C of the Newton correction is zero, Z must be reset to zero.
C-----------------------------------------------------------------------
 200  CONTINUE
      LL = LGMR
      LLP1 = LL + 1
      DO 210 K = 1,LLP1
 210    R(K) = 0.0D0
      R(1) = RNRM
      CALL DHELS (HES, MAXLP1, LL, Q, R)
      DO 220 K = 1,NEQ
 220    Z(K) = 0.0D0
      DO 230 I = 1,LL
        CALL DAXPY (NEQ, R(I), V(1,I), 1, Z, 1)
 230    CONTINUE
      DO 240 I = 1,NEQ
 240    Z(I) = Z(I)/WGHT(I)
C Load RHO into RHOK.
      RHOK = RHO
      RETURN
C-----------------------------------------------------------------------
C This block handles error returns forced by routine PSOL.
C-----------------------------------------------------------------------
 300  CONTINUE
      IF (IER .LT. 0) IFLAG = -1
      IF (IER .GT. 0) IFLAG = 3
C
      RETURN
C
C------END OF SUBROUTINE DSPIGM-----------------------------------------
      END
      SUBROUTINE DATV (NEQ, Y, TN, YPRIME, SAVR, V, WGHT, YPTEM, RES,
     *   IRES, PSOL, Z, VTEM, WP, IWP, CJ, EPLIN, IER, NRE, NPSL,
     *   RPAR,IPAR)
C
C***BEGIN PROLOGUE  DATV
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine computes the product
C
C   Z = (D-inverse)*(P-inverse)*(dF/dY)*(D*V),
C
C where F(Y) = G(T, Y, CJ*(Y-A)), CJ is a scalar proportional to 1/H,
C and A involves the past history of Y.  The quantity CJ*(Y-A) is
C an approximation to the first derivative of Y and is stored
C in the array YPRIME.  Note that dF/dY = dG/dY + CJ*dG/dYPRIME.
C
C D is a diagonal scaling matrix, and P is the left preconditioning
C matrix.  V is assumed to have L2 norm equal to 1.
C The product is stored in Z and is computed by means of a
C difference quotient, a call to RES, and one call to PSOL.
C
C      On entry
C
C          NEQ = Problem size, passed to RES and PSOL.
C
C            Y = Array containing current dependent variable vector.
C
C       YPRIME = Array containing current first derivative of y.
C
C         SAVR = Array containing current value of G(T,Y,YPRIME).
C
C            V = Real array of length NEQ (can be the same array as Z).
C
C         WGHT = Array of length NEQ containing scale factors.
C                1/WGHT(I) are the diagonal elements of the matrix D.
C
C        YPTEM = Work array of length NEQ.
C
C         VTEM = Work array of length NEQ used to store the
C                unscaled version of V.
C
C         WP = Real work array used by preconditioner PSOL.
C
C         IWP = Integer work array used by preconditioner PSOL.
C
C           CJ = Scalar proportional to current value of 
C                1/(step size H).
C
C
C      On return
C
C            Z = Array of length NEQ containing desired scaled
C                matrix-vector product.
C
C         IRES = Error flag from RES.
C
C          IER = Error flag from PSOL.
C
C         NRE  = The number of calls to RES.
C
C         NPSL = The number of calls to PSOL.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   RES, PSOL
C
C***END PROLOGUE  DATV
C
      INTEGER NEQ, IRES, IWP, IER, NRE, NPSL, IPAR
      DOUBLE PRECISION Y, TN, YPRIME, SAVR, V, WGHT, YPTEM, Z, VTEM,
     1   WP, CJ, RPAR
      DIMENSION Y(*), YPRIME(*), SAVR(*), V(*), WGHT(*), YPTEM(*),
     1   Z(*), VTEM(*), WP(*), IWP(*), RPAR(*), IPAR(*)
      INTEGER I
      DOUBLE PRECISION EPLIN
      EXTERNAL  RES, PSOL
C
      IRES = 0
C-----------------------------------------------------------------------
C Set VTEM = D * V.
C-----------------------------------------------------------------------
      DO 10 I = 1,NEQ
 10     VTEM(I) = V(I)/WGHT(I)
      IER = 0
C-----------------------------------------------------------------------
C Store Y in Z and increment Z by VTEM.
C Store YPRIME in YPTEM and increment YPTEM by VTEM*CJ.
C-----------------------------------------------------------------------
      DO 20 I = 1,NEQ
        YPTEM(I) = YPRIME(I) + VTEM(I)*CJ
 20     Z(I) = Y(I) + VTEM(I)
C-----------------------------------------------------------------------
C Call RES with incremented Y, YPRIME arguments
C stored in Z, YPTEM.  VTEM is overwritten with new residual.
C-----------------------------------------------------------------------
      CONTINUE
      CALL RES(TN,Z,YPTEM,CJ,VTEM,IRES,RPAR,IPAR)
      NRE = NRE + 1
      IF (IRES .LT. 0) RETURN
C-----------------------------------------------------------------------
C Set Z = (dF/dY) * VBAR using difference quotient.
C (VBAR is old value of VTEM before calling RES)
C-----------------------------------------------------------------------
      DO 70 I = 1,NEQ
 70     Z(I) = VTEM(I) - SAVR(I)
C-----------------------------------------------------------------------
C Apply inverse of left preconditioner to Z.
C-----------------------------------------------------------------------
      CALL PSOL (NEQ, TN, Y, YPRIME, SAVR, YPTEM, CJ, WGHT, WP, IWP,
     1   Z, EPLIN, IER, RPAR, IPAR)
      NPSL = NPSL + 1
      IF (IER .NE. 0) RETURN
C-----------------------------------------------------------------------
C Apply D-inverse to Z and return.
C-----------------------------------------------------------------------
      DO 90 I = 1,NEQ
 90     Z(I) = Z(I)*WGHT(I)
      RETURN
C
C------END OF SUBROUTINE DATV-------------------------------------------
      END
      SUBROUTINE DORTH (VNEW, V, HES, N, LL, LDHES, KMP, SNORMW)
C
C***BEGIN PROLOGUE  DORTH
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This routine orthogonalizes the vector VNEW against the previous
C KMP vectors in the V array.  It uses a modified Gram-Schmidt
C orthogonalization procedure with conditional reorthogonalization.
C
C      On entry
C
C         VNEW = The vector of length N containing a scaled product
C                OF The Jacobian and the vector V(*,LL).
C
C         V    = The N x LL array containing the previous LL
C                orthogonal vectors V(*,1) to V(*,LL).
C
C         HES  = An LL x LL upper Hessenberg matrix containing,
C                in HES(I,K), K.LT.LL, scaled inner products of
C                A*V(*,K) and V(*,I).
C
C        LDHES = The leading dimension of the HES array.
C
C         N    = The order of the matrix A, and the length of VNEW.
C
C         LL   = The current order of the matrix HES.
C
C          KMP = The number of previous vectors the new vector VNEW
C                must be made orthogonal to (KMP .LE. MAXL).
C
C
C      On return
C
C         VNEW = The new vector orthogonal to V(*,I0),
C                where I0 = MAX(1, LL-KMP+1).
C
C         HES  = Upper Hessenberg matrix with column LL filled in with
C                scaled inner products of A*V(*,LL) and V(*,I).
C
C       SNORMW = L-2 norm of VNEW.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DDOT, DNRM2, DAXPY 
C
C***END PROLOGUE  DORTH
C
      INTEGER N, LL, LDHES, KMP
      DOUBLE PRECISION VNEW, V, HES, SNORMW
      DIMENSION VNEW(*), V(N,*), HES(LDHES,*)
      INTEGER I, I0
      DOUBLE PRECISION ARG, DDOT, DNRM2, SUMDSQ, TEM, VNRM
C
C-----------------------------------------------------------------------
C Get norm of unaltered VNEW for later use.
C-----------------------------------------------------------------------
      VNRM = DNRM2 (N, VNEW, 1)
C-----------------------------------------------------------------------
C Do Modified Gram-Schmidt on VNEW = A*V(LL).
C Scaled inner products give new column of HES.
C Projections of earlier vectors are subtracted from VNEW.
C-----------------------------------------------------------------------
      I0 = MAX0(1,LL-KMP+1)
      DO 10 I = I0,LL
        HES(I,LL) = DDOT (N, V(1,I), 1, VNEW, 1)
        TEM = -HES(I,LL)
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
 10     CONTINUE
C-----------------------------------------------------------------------
C Compute SNORMW = norm of VNEW.
C If VNEW is small compared to its input value (in norm), then
C Reorthogonalize VNEW to V(*,1) through V(*,LL).
C Correct if relative correction exceeds 1000*(unit roundoff).
C Finally, correct SNORMW using the dot products involved.
C-----------------------------------------------------------------------
      SNORMW = DNRM2 (N, VNEW, 1)
      IF (VNRM + 0.001D0*SNORMW .NE. VNRM) RETURN
      SUMDSQ = 0.0D0
      DO 30 I = I0,LL
        TEM = -DDOT (N, V(1,I), 1, VNEW, 1)
        IF (HES(I,LL) + 0.001D0*TEM .EQ. HES(I,LL)) GO TO 30
        HES(I,LL) = HES(I,LL) - TEM
        CALL DAXPY (N, TEM, V(1,I), 1, VNEW, 1)
        SUMDSQ = SUMDSQ + TEM**2
 30     CONTINUE
      IF (SUMDSQ .EQ. 0.0D0) RETURN
      ARG = MAX(0.0D0,SNORMW**2 - SUMDSQ)
      SNORMW = SQRT(ARG)
      RETURN
C
C------END OF SUBROUTINE DORTH------------------------------------------
      END
      SUBROUTINE DHEQR (A, LDA, N, Q, INFO, IJOB)
C
C***BEGIN PROLOGUE  DHEQR
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C     This routine performs a QR decomposition of an upper
C     Hessenberg matrix A.  There are two options available:
C
C          (1)  performing a fresh decomposition
C          (2)  updating the QR factors by adding a row and A
C               column to the matrix A.
C
C     DHEQR decomposes an upper Hessenberg matrix by using Givens
C     rotations.
C
C     On entry
C
C        A       DOUBLE PRECISION(LDA, N)
C                The matrix to be decomposed.
C
C        LDA     INTEGER
C                The leading dimension of the array A.
C
C        N       INTEGER
C                A is an (N+1) by N Hessenberg matrix.
C
C        IJOB    INTEGER
C                = 1     Means that a fresh decomposition of the
C                        matrix A is desired.
C                .GE. 2  Means that the current decomposition of A
C                        will be updated by the addition of a row
C                        and a column.
C     On return
C
C        A       The upper triangular matrix R.
C                The factorization can be written Q*A = R, where
C                Q is a product of Givens rotations and R is upper
C                triangular.
C
C        Q       DOUBLE PRECISION(2*N)
C                The factors C and S of each Givens rotation used
C                in decomposing A.
C
C        INFO    INTEGER
C                = 0  normal value.
C                = K  If  A(K,K) .EQ. 0.0.  This is not an error
C                     condition for this subroutine, but it does
C                     indicate that DHELS will divide by zero
C                     if called.
C
C     Modification of LINPACK.
C     Peter Brown, Lawrence Livermore Natl. Lab.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED (NONE)
C
C***END PROLOGUE  DHEQR
C
      INTEGER LDA, N, INFO, IJOB
      DOUBLE PRECISION A(LDA,*), Q(*)
      INTEGER I, IQ, J, K, KM1, KP1, NM1
      DOUBLE PRECISION C, S, T, T1, T2
C
      IF (IJOB .GT. 1) GO TO 70
C-----------------------------------------------------------------------
C A new factorization is desired.
C-----------------------------------------------------------------------
C
C     QR decomposition without pivoting.
C
      INFO = 0
      DO 60 K = 1, N
         KM1 = K - 1
         KP1 = K + 1
C
C           Compute Kth column of R.
C           First, multiply the Kth column of A by the previous
C           K-1 Givens rotations.
C
            IF (KM1 .LT. 1) GO TO 20
            DO 10 J = 1, KM1
              I = 2*(J-1) + 1
              T1 = A(J,K)
              T2 = A(J+1,K)
              C = Q(I)
              S = Q(I+1)
              A(J,K) = C*T1 - S*T2
              A(J+1,K) = S*T1 + C*T2
   10         CONTINUE
C
C           Compute Givens components C and S.
C
   20       CONTINUE
            IQ = 2*KM1 + 1
            T1 = A(K,K)
            T2 = A(KP1,K)
            IF (T2 .NE. 0.0D0) GO TO 30
              C = 1.0D0
              S = 0.0D0
              GO TO 50
   30       CONTINUE
            IF (ABS(T2) .LT. ABS(T1)) GO TO 40
              T = T1/T2
              S = -1.0D0/SQRT(1.0D0+T*T)
              C = -S*T
              GO TO 50
   40       CONTINUE
              T = T2/T1
              C = 1.0D0/SQRT(1.0D0+T*T)
              S = -C*T
   50       CONTINUE
            Q(IQ) = C
            Q(IQ+1) = S
            A(K,K) = C*T1 - S*T2
            IF (A(K,K) .EQ. 0.0D0) INFO = K
   60 CONTINUE
      RETURN
C-----------------------------------------------------------------------
C The old factorization of A will be updated.  A row and a column
C has been added to the matrix A.
C N by N-1 is now the old size of the matrix.
C-----------------------------------------------------------------------
  70  CONTINUE
      NM1 = N - 1
C-----------------------------------------------------------------------
C Multiply the new column by the N previous Givens rotations.
C-----------------------------------------------------------------------
      DO 100 K = 1,NM1
        I = 2*(K-1) + 1
        T1 = A(K,N)
        T2 = A(K+1,N)
        C = Q(I)
        S = Q(I+1)
        A(K,N) = C*T1 - S*T2
        A(K+1,N) = S*T1 + C*T2
 100    CONTINUE
C-----------------------------------------------------------------------
C Complete update of decomposition by forming last Givens rotation,
C and multiplying it times the column vector (A(N,N),A(NP1,N)).
C-----------------------------------------------------------------------
      INFO = 0
      T1 = A(N,N)
      T2 = A(N+1,N)
      IF (T2 .NE. 0.0D0) GO TO 110
        C = 1.0D0
        S = 0.0D0
        GO TO 130
 110  CONTINUE
      IF (ABS(T2) .LT. ABS(T1)) GO TO 120
        T = T1/T2
        S = -1.0D0/SQRT(1.0D0+T*T)
        C = -S*T
        GO TO 130
 120  CONTINUE
        T = T2/T1
        C = 1.0D0/SQRT(1.0D0+T*T)
        S = -C*T
 130  CONTINUE
      IQ = 2*N - 1
      Q(IQ) = C
      Q(IQ+1) = S
      A(N,N) = C*T1 - S*T2
      IF (A(N,N) .EQ. 0.0D0) INFO = N
      RETURN
C
C------END OF SUBROUTINE DHEQR------------------------------------------
      END
      SUBROUTINE DHELS (A, LDA, N, Q, B)
C
C***BEGIN PROLOGUE  DHELS
C***DATE WRITTEN   890101   (YYMMDD)
C***REVISION DATE  900926   (YYMMDD)
C
C
C-----------------------------------------------------------------------
C***DESCRIPTION
C
C This is similar to the LINPACK routine DGESL except that
C A is an upper Hessenberg matrix.
C
C     DHELS solves the least squares problem
C
C           MIN (B-A*X,B-A*X)
C
C     using the factors computed by DHEQR.
C
C     On entry
C
C        A       DOUBLE PRECISION (LDA, N)
C                The output from DHEQR which contains the upper
C                triangular factor R in the QR decomposition of A.
C
C        LDA     INTEGER
C                The leading dimension of the array  A .
C
C        N       INTEGER
C                A is originally an (N+1) by N matrix.
C
C        Q       DOUBLE PRECISION(2*N)
C                The coefficients of the N givens rotations
C                used in the QR factorization of A.
C
C        B       DOUBLE PRECISION(N+1)
C                The right hand side vector.
C
C
C     On return
C
C        B       The solution vector X.
C
C
C     Modification of LINPACK.
C     Peter Brown, Lawrence Livermore Natl. Lab.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED
C   DAXPY 
C
C***END PROLOGUE  DHELS
C
      INTEGER LDA, N
      DOUBLE PRECISION A(LDA,*), B(*), Q(*)
      INTEGER IQ, K, KB, KP1
      DOUBLE PRECISION C, S, T, T1, T2
C
C        Minimize (B-A*X,B-A*X).
C        First form Q*B.
C
         DO 20 K = 1, N
            KP1 = K + 1
            IQ = 2*(K-1) + 1
            C = Q(IQ)
            S = Q(IQ+1)
            T1 = B(K)
            T2 = B(KP1)
            B(K) = C*T1 - S*T2
            B(KP1) = S*T1 + C*T2
   20    CONTINUE
C
C        Now solve R*X = Q*B.
C
         DO 40 KB = 1, N
            K = N + 1 - KB
            B(K) = B(K)/A(K,K)
            T = -B(K)
            CALL DAXPY (K-1, T, A(1,K), 1, B(1), 1)
   40    CONTINUE
      RETURN
C
C------END OF SUBROUTINE DHELS------------------------------------------
      END
