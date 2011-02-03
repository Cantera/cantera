*DECK DPCOEF
      SUBROUTINE DPCOEF (L, C, TC, A)
C***BEGIN PROLOGUE  DPCOEF
C***PURPOSE  Convert the DPOLFT coefficients to Taylor series form.
C***LIBRARY   SLATEC
C***CATEGORY  K1A1A2
C***TYPE      DOUBLE PRECISION (PCOEF-S, DPCOEF-D)
C***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
C***AUTHOR  Shampine, L. F., (SNLA)
C           Davenport, S. M., (SNLA)
C***DESCRIPTION
C
C     Abstract
C
C     DPOLFT  computes the least squares polynomial fit of degree  L  as
C     a sum of orthogonal polynomials.  DPCOEF  changes this fit to its
C     Taylor expansion about any point  C , i.e. writes the polynomial
C     as a sum of powers of (X-C).  Taking  C=0.  gives the polynomial
C     in powers of X, but a suitable non-zero  C  often leads to
C     polynomials which are better scaled and more accurately evaluated.
C
C     The parameters for  DPCOEF  are
C
C     INPUT -- All TYPE REAL variables are DOUBLE PRECISION
C         L -      Indicates the degree of polynomial to be changed to
C                  its Taylor expansion.  To obtain the Taylor
C                  coefficients in reverse order, input  L  as the
C                  negative of the degree desired.  The absolute value
C                  of L  must be less than or equal to NDEG, the highest
C                  degree polynomial fitted by  DPOLFT .
C         C -      The point about which the Taylor expansion is to be
C                  made.
C         A -      Work and output array containing values from last
C                  call to  DPOLFT .
C
C     OUTPUT -- All TYPE REAL variables are DOUBLE PRECISION
C         TC -     Vector containing the first LL+1 Taylor coefficients
C                  where LL=ABS(L).  If  L.GT.0 , the coefficients are
C                  in the usual Taylor series order, i.e.
C                    P(X) = TC(1) + TC(2)*(X-C) + ... + TC(N+1)*(X-C)**N
C                  If L .LT. 0, the coefficients are in reverse order,
C                  i.e.
C                    P(X) = TC(1)*(X-C)**N + ... + TC(N)*(X-C) + TC(N+1)
C
C***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
C                 Curve fitting by polynomials in one variable, Report
C                 SLA-74-0270, Sandia Laboratories, June 1974.
C***ROUTINES CALLED  DP1VLU
C***REVISION HISTORY  (YYMMDD)
C   740601  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   891006  Cosmetic changes to prologue.  (WRB)
C   891006  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C***END PROLOGUE  DPCOEF
C
      INTEGER I,L,LL,LLP1,LLP2,NEW,NR
      DOUBLE PRECISION A(*),C,FAC,SAVE,TC(*)
C***FIRST EXECUTABLE STATEMENT  DPCOEF
      LL = ABS(L)
      LLP1 = LL + 1
      CALL DP1VLU (LL,LL,C,TC(1),TC(2),A)
      IF (LL .LT. 2) GO TO 2
      FAC = 1.0D0
      DO 1 I = 3,LLP1
        FAC = FAC*(I-1)
 1      TC(I) = TC(I)/FAC
 2    IF (L .GE. 0) GO TO 4
      NR = LLP1/2
      LLP2 = LL + 2
      DO 3 I = 1,NR
        SAVE = TC(I)
        NEW = LLP2 - I
        TC(I) = TC(NEW)
 3      TC(NEW) = SAVE
 4    RETURN
      END
