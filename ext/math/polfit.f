*DECK POLFIT
      SUBROUTINE POLFIT (N, X, Y, W, MAXDEG, NDEG, EPS, R, IERR, A)
C***BEGIN PROLOGUE  POLFIT
C***PURPOSE  Fit discrete data in a least squares sense by polynomials
C            in one variable.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A2
C***TYPE      SINGLE PRECISION (POLFIT-S, DPOLFT-D)
C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
C***AUTHOR  Shampine, L. F., (SNLA)
C           Davenport, S. M., (SNLA)
C           Huddleston, R. E., (SNLL)
C***DESCRIPTION
C
C     Abstract
C
C     Given a collection of points X(I) and a set of values Y(I) which
C     correspond to some function or measurement at each of the X(I),
C     subroutine  POLFIT  computes the weighted least-squares polynomial
C     fits of all degrees up to some degree either specified by the user
C     or determined by the routine.  The fits thus obtained are in
C     orthogonal polynomial form.  Subroutine  PVALUE  may then be
C     called to evaluate the fitted polynomials and any of their
C     derivatives at any point.  The subroutine  PCOEF  may be used to
C     express the polynomial fits as powers of (X-C) for any specified
C     point C.
C
C     The parameters for  POLFIT  are
C
C     Input --
C         N -      the number of data points.  The arrays X, Y and W
C                  must be dimensioned at least  N  (N .GE. 1).
C         X -      array of values of the independent variable.  These
C                  values may appear in any order and need not all be
C                  distinct.
C         Y -      array of corresponding function values.
C         W -      array of positive values to be used as weights.  If
C                  W(1) is negative,  POLFIT  will set all the weights
C                  to 1.0, which means unweighted least squares error
C                  will be minimized.  To minimize relative error, the
C                  user should set the weights to:  W(I) = 1.0/Y(I)**2,
C                  I = 1,...,N .
C         MAXDEG - maximum degree to be allowed for polynomial fit.
C                  MAXDEG  may be any non-negative integer less than  N.
C                  Note -- MAXDEG  cannot be equal to  N-1  when a
C                  statistical test is to be used for degree selection,
C                  i.e., when input value of  EPS  is negative.
C         EPS -    specifies the criterion to be used in determining
C                  the degree of fit to be computed.
C                  (1)  If  EPS  is input negative,  POLFIT  chooses the
C                       degree based on a statistical F test of
C                       significance.  One of three possible
C                       significance levels will be used:  .01, .05 or
C                       .10.  If  EPS=-1.0 , the routine will
C                       automatically select one of these levels based
C                       on the number of data points and the maximum
C                       degree to be considered.  If  EPS  is input as
C                       -.01, -.05, or -.10, a significance level of
C                       .01, .05, or .10, respectively, will be used.
C                  (2)  If  EPS  is set to 0.,  POLFIT  computes the
C                       polynomials of degrees 0 through  MAXDEG .
C                  (3)  If  EPS  is input positive,  EPS  is the RMS
C                       error tolerance which must be satisfied by the
C                       fitted polynomial.  POLFIT  will increase the
C                       degree of fit until this criterion is met or
C                       until the maximum degree is reached.
C
C     Output --
C         NDEG -   degree of the highest degree fit computed.
C         EPS -    RMS error of the polynomial of degree  NDEG .
C         R -      vector of dimension at least NDEG containing values
C                  of the fit of degree  NDEG  at each of the  X(I) .
C                  Except when the statistical test is used, these
C                  values are more accurate than results from subroutine
C                  PVALUE  normally are.
C         IERR -   error flag with the following possible values.
C             1 -- indicates normal execution, i.e., either
C                  (1)  the input value of  EPS  was negative, and the
C                       computed polynomial fit of degree  NDEG
C                       satisfies the specified F test, or
C                  (2)  the input value of  EPS  was 0., and the fits of
C                       all degrees up to  MAXDEG  are complete, or
C                  (3)  the input value of  EPS  was positive, and the
C                       polynomial of degree  NDEG  satisfies the RMS
C                       error requirement.
C             2 -- invalid input parameter.  At least one of the input
C                  parameters has an illegal value and must be corrected
C                  before  POLFIT  can proceed.  Valid input results
C                  when the following restrictions are observed
C                       N .GE. 1
C                       0 .LE. MAXDEG .LE. N-1  for  EPS .GE. 0.
C                       0 .LE. MAXDEG .LE. N-2  for  EPS .LT. 0.
C                       W(1)=-1.0  or  W(I) .GT. 0., I=1,...,N .
C             3 -- cannot satisfy the RMS error requirement with a
C                  polynomial of degree no greater than  MAXDEG .  Best
C                  fit found is of degree  MAXDEG .
C             4 -- cannot satisfy the test for significance using
C                  current value of  MAXDEG .  Statistically, the
C                  best fit found is of order  NORD .  (In this case,
C                  NDEG will have one of the values:  MAXDEG-2,
C                  MAXDEG-1, or MAXDEG).  Using a higher value of
C                  MAXDEG  may result in passing the test.
C         A -      work and output array having at least 3N+3MAXDEG+3
C                  locations
C
C     Note - POLFIT  calculates all fits of degrees up to and including
C            NDEG .  Any or all of these fits can be evaluated or
C            expressed as powers of (X-C) using  PVALUE  and  PCOEF
C            after just one call to  POLFIT .
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  PVALUE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
C***END PROLOGUE  POLFIT
      DOUBLE PRECISION TEMD1,TEMD2
      DIMENSION X(*), Y(*), W(*), R(*), A(*)
      DIMENSION CO(4,3)
      SAVE CO
      DATA  CO(1,1), CO(2,1), CO(3,1), CO(4,1), CO(1,2), CO(2,2),
     1      CO(3,2), CO(4,2), CO(1,3), CO(2,3), CO(3,3),
     2  CO(4,3)/-13.086850,-2.4648165,-3.3846535,-1.2973162,
     3          -3.3381146,-1.7812271,-3.2578406,-1.6589279,
     4          -1.6282703,-1.3152745,-3.2640179,-1.9829776/
C***FIRST EXECUTABLE STATEMENT  POLFIT
      M = ABS(N)
      IF (M .EQ. 0) GO TO 30
      IF (MAXDEG .LT. 0) GO TO 30
      A(1) = MAXDEG
      MOP1 = MAXDEG + 1
      IF (M .LT. MOP1) GO TO 30
      IF (EPS .LT. 0.0  .AND.  M .EQ. MOP1) GO TO 30
      XM = M
      ETST = EPS*EPS*XM
      IF (W(1) .LT. 0.0) GO TO 2
      DO 1 I = 1,M
        IF (W(I) .LE. 0.0) GO TO 30
 1      CONTINUE
      GO TO 4
 2    DO 3 I = 1,M
 3      W(I) = 1.0
 4    IF (EPS .GE. 0.0) GO TO 8
C
C DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
C CHOOSING DEGREE OF POLYNOMIAL FIT
C
      IF (EPS .GT. (-.55)) GO TO 5
      IDEGF = M - MAXDEG - 1
      KSIG = 1
      IF (IDEGF .LT. 10) KSIG = 2
      IF (IDEGF .LT. 5) KSIG = 3
      GO TO 8
 5    KSIG = 1
      IF (EPS .LT. (-.03)) KSIG = 2
      IF (EPS .LT. (-.07)) KSIG = 3
C
C INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
C
 8    K1 = MAXDEG + 1
      K2 = K1 + MAXDEG
      K3 = K2 + MAXDEG + 2
      K4 = K3 + M
      K5 = K4 + M
      DO 9 I = 2,K4
 9      A(I) = 0.0
      W11 = 0.0
      IF (N .LT. 0) GO TO 11
C
C UNCONSTRAINED CASE
C
      DO 10 I = 1,M
        K4PI = K4 + I
        A(K4PI) = 1.0
 10     W11 = W11 + W(I)
      GO TO 13
C
C CONSTRAINED CASE
C
 11   DO 12 I = 1,M
        K4PI = K4 + I
 12     W11 = W11 + W(I)*A(K4PI)**2
C
C COMPUTE FIT OF DEGREE ZERO
C
 13   TEMD1 = 0.0D0
      DO 14 I = 1,M
        K4PI = K4 + I
        TEMD1 = TEMD1 + DBLE(W(I))*DBLE(Y(I))*DBLE(A(K4PI))
 14     CONTINUE
      TEMD1 = TEMD1/DBLE(W11)
      A(K2+1) = TEMD1
      SIGJ = 0.0
      DO 15 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = TEMD1*DBLE(A(K4PI))
        R(I) = TEMD2
        A(K5PI) = TEMD2 - DBLE(R(I))
 15     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
      J = 0
C
C SEE IF POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
C
      IF (EPS) 24,26,27
C
C INCREMENT DEGREE
C
 16   J = J + 1
      JP1 = J + 1
      K1PJ = K1 + J
      K2PJ = K2 + J
      SIGJM1 = SIGJ
C
C COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
C
      IF (J .GT. 1) A(K1PJ) = W11/W1
C
C COMPUTE NEW A COEFFICIENT
C
      TEMD1 = 0.0D0
      DO 18 I = 1,M
        K4PI = K4 + I
        TEMD2 = A(K4PI)
        TEMD1 = TEMD1 + DBLE(X(I))*DBLE(W(I))*TEMD2*TEMD2
 18     CONTINUE
      A(JP1) = TEMD1/DBLE(W11)
C
C EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
C
      W1 = W11
      W11 = 0.0
      DO 19 I = 1,M
        K3PI = K3 + I
        K4PI = K4 + I
        TEMP = A(K3PI)
        A(K3PI) = A(K4PI)
        A(K4PI) = (X(I)-A(JP1))*A(K3PI) - A(K1PJ)*TEMP
 19     W11 = W11 + W(I)*A(K4PI)**2
C
C GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
C PRECISION
C
      TEMD1 = 0.0D0
      DO 20 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = DBLE(W(I))*DBLE((Y(I)-R(I))-A(K5PI))*DBLE(A(K4PI))
 20     TEMD1 = TEMD1 + TEMD2
      TEMD1 = TEMD1/DBLE(W11)
      A(K2PJ+1) = TEMD1
C
C UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
C ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
C COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
C THE MOST SIGNIFICANT BITS ARE STORED IN  R(I) , AND THE LEAST
C SIGNIFICANT BITS ARE IN  A(K5PI) .
C
      SIGJ = 0.0
      DO 21 I = 1,M
        K4PI = K4 + I
        K5PI = K5 + I
        TEMD2 = DBLE(R(I)) + DBLE(A(K5PI)) + TEMD1*DBLE(A(K4PI))
        R(I) = TEMD2
        A(K5PI) = TEMD2 - DBLE(R(I))
 21     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
C
C SEE IF DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
C MAXDEG  HAS BEEN REACHED
C
      IF (EPS) 23,26,27
C
C COMPUTE F STATISTICS  (INPUT EPS .LT. 0.)
C
 23   IF (SIGJ .EQ. 0.0) GO TO 29
      DEGF = M - J - 1
      DEN = (CO(4,KSIG)*DEGF + 1.0)*DEGF
      FCRIT = (((CO(3,KSIG)*DEGF) + CO(2,KSIG))*DEGF + CO(1,KSIG))/DEN
      FCRIT = FCRIT*FCRIT
      F = (SIGJM1 - SIGJ)*DEGF/SIGJ
      IF (F .LT. FCRIT) GO TO 25
C
C POLYNOMIAL OF DEGREE J SATISFIES F TEST
C
 24   SIGPAS = SIGJ
      JPAS = J
      NFAIL = 0
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
C
C POLYNOMIAL OF DEGREE J FAILS F TEST.  IF THERE HAVE BEEN THREE
C SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
C
 25   NFAIL = NFAIL + 1
      IF (NFAIL .GE. 3) GO TO 29
      IF (MAXDEG .EQ. J) GO TO 32
      GO TO 16
C
C RAISE THE DEGREE IF DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
C EPS = 0.)
C
 26   IF (MAXDEG .EQ. J) GO TO 28
      GO TO 16
C
C SEE IF RMS ERROR CRITERION IS SATISFIED  (INPUT EPS .GT. 0.)
C
 27   IF (SIGJ .LE. ETST) GO TO 28
      IF (MAXDEG .EQ. J) GO TO 31
      GO TO 16
C
C RETURNS
C
 28   IERR = 1
      NDEG = J
      SIG = SIGJ
      GO TO 33
 29   IERR = 1
      NDEG = JPAS
      SIG = SIGPAS
      GO TO 33
 30   IERR = 2
      CALL XERMSG ('SLATEC', 'POLFIT', 'INVALID INPUT PARAMETER.', 2,
     +   1)
      GO TO 37
 31   IERR = 3
      NDEG = MAXDEG
      SIG = SIGJ
      GO TO 33
 32   IERR = 4
      NDEG = JPAS
      SIG = SIGPAS
C
 33   A(K3) = NDEG
C
C WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
C ALL THE DATA POINTS IF  R  DOES NOT ALREADY CONTAIN THESE VALUES
C
      IF(EPS .GE. 0.0  .OR.  NDEG .EQ. MAXDEG) GO TO 36
      NDER = 0
      DO 35 I = 1,M
        CALL PVALUE (NDEG,NDER,X(I),R(I),YP,A)
 35     CONTINUE
 36   EPS = SQRT(SIG/XM)
 37   RETURN
      END
