      SUBROUTINE DLASQ2( M, Q, E, QQ, EE, EPS, TOL2, SMALL2, SUP, KEND,
     $                   INFO )
*
*  -- LAPACK routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, KEND, M
      DOUBLE PRECISION   EPS, SMALL2, SUP, TOL2
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   E( * ), EE( * ), Q( * ), QQ( * )
*     ..
*
*     Purpose
*     =======
*
*     DLASQ2 computes the singular values of a real N-by-N unreduced
*     bidiagonal matrix with squared diagonal elements in Q and
*     squared off-diagonal elements in E. The singular values are
*     computed to relative accuracy TOL, barring over/underflow or
*     denormalization.
*
*     Arguments
*     =========
*
*  M       (input) INTEGER
*          The number of rows and columns in the matrix. M >= 0.
*
*  Q       (output) DOUBLE PRECISION array, dimension (M)
*          On normal exit, contains the squared singular values.
*
*  E       (workspace) DOUBLE PRECISION array, dimension (M)
*
*  QQ      (input/output) DOUBLE PRECISION array, dimension (M)
*          On entry, QQ contains the squared diagonal elements of the
*          bidiagonal matrix whose SVD is desired.
*          On exit, QQ is overwritten.
*
*  EE      (input/output) DOUBLE PRECISION array, dimension (M)
*          On entry, EE(1:N-1) contains the squared off-diagonal
*          elements of the bidiagonal matrix whose SVD is desired.
*          On exit, EE is overwritten.
*
*  EPS     (input) DOUBLE PRECISION
*          Machine epsilon.
*
*  TOL2    (input) DOUBLE PRECISION
*          Desired relative accuracy of computed eigenvalues
*          as defined in DLASQ1.
*
*  SMALL2  (input) DOUBLE PRECISION
*          A threshold value as defined in DLASQ1.
*
*  SUP     (input/output) DOUBLE PRECISION
*          Upper bound for the smallest eigenvalue.
*
*  KEND    (input/output) INTEGER
*          Index where minimum d occurs.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm did not converge;  i
*                specifies how many superdiagonals did not converge.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      DOUBLE PRECISION   FOUR, HALF
      PARAMETER          ( FOUR = 4.0D+0, HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            ICONV, IPHASE, ISP, N, OFF, OFF1
      DOUBLE PRECISION   QEMAX, SIGMA, XINF, XX, YY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASQ3
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, NINT, SQRT
*     ..
*     .. Executable Statements ..
      N = M
*
*     Set the default maximum number of iterations
*
      OFF = 0
      OFF1 = OFF + 1
      SIGMA = ZERO
      XINF = ZERO
      ICONV = 0
      IPHASE = 2
*
*     Try deflation at the bottom
*
*     1x1 deflation
*
   10 CONTINUE
      IF( N.LE.2 )
     $   GO TO 20
      IF( EE( N-1 ).LE.MAX( QQ( N ), XINF, SMALL2 )*TOL2 ) THEN
         Q( N ) = QQ( N )
         N = N - 1
         IF( KEND.GT.N )
     $      KEND = N
         SUP = MIN( QQ( N ), QQ( N-1 ) )
         GO TO 10
      END IF
*
*     2x2 deflation
*
      IF( EE( N-2 ).LE.MAX( XINF, SMALL2,
     $    ( QQ( N ) / ( QQ( N )+EE( N-1 )+QQ( N-1 ) ) )*QQ( N-1 ) )*
     $    TOL2 ) THEN
         QEMAX = MAX( QQ( N ), QQ( N-1 ), EE( N-1 ) )
         IF( QEMAX.NE.ZERO ) THEN
            IF( QEMAX.EQ.QQ( N-1 ) ) THEN
               XX = HALF*( QQ( N )+QQ( N-1 )+EE( N-1 )+QEMAX*
     $              SQRT( ( ( QQ( N )-QQ( N-1 )+EE( N-1 ) ) /
     $              QEMAX )**2+FOUR*EE( N-1 ) / QEMAX ) )
            ELSE IF( QEMAX.EQ.QQ( N ) ) THEN
               XX = HALF*( QQ( N )+QQ( N-1 )+EE( N-1 )+QEMAX*
     $              SQRT( ( ( QQ( N-1 )-QQ( N )+EE( N-1 ) ) /
     $              QEMAX )**2+FOUR*EE( N-1 ) / QEMAX ) )
            ELSE
               XX = HALF*( QQ( N )+QQ( N-1 )+EE( N-1 )+QEMAX*
     $              SQRT( ( ( QQ( N )-QQ( N-1 )+EE( N-1 ) ) /
     $              QEMAX )**2+FOUR*QQ( N-1 ) / QEMAX ) )
            END IF
            YY = ( MAX( QQ( N ), QQ( N-1 ) ) / XX )*
     $           MIN( QQ( N ), QQ( N-1 ) )
         ELSE
            XX = ZERO
            YY = ZERO
         END IF
         Q( N-1 ) = XX
         Q( N ) = YY
         N = N - 2
         IF( KEND.GT.N )
     $      KEND = N
         SUP = QQ( N )
         GO TO 10
      END IF
*
   20 CONTINUE
      IF( N.EQ.0 ) THEN
*
*         The lower branch is finished
*
         IF( OFF.EQ.0 ) THEN
*
*         No upper branch; return to DLASQ1
*
            RETURN
         ELSE
*
*         Going back to upper branch
*
            XINF = ZERO
            IF( EE( OFF ).GT.ZERO ) THEN
               ISP = NINT( EE( OFF ) )
               IPHASE = 1
            ELSE
               ISP = -NINT( EE( OFF ) )
               IPHASE = 2
            END IF
            SIGMA = E( OFF )
            N = OFF - ISP + 1
            OFF1 = ISP
            OFF = OFF1 - 1
            IF( N.LE.2 )
     $         GO TO 20
            IF( IPHASE.EQ.1 ) THEN
               SUP = MIN( Q( N+OFF ), Q( N-1+OFF ), Q( N-2+OFF ) )
            ELSE
               SUP = MIN( QQ( N+OFF ), QQ( N-1+OFF ), QQ( N-2+OFF ) )
            END IF
            KEND = 0
            ICONV = -3
         END IF
      ELSE IF( N.EQ.1 ) THEN
*
*     1x1 Solver
*
         IF( IPHASE.EQ.1 ) THEN
            Q( OFF1 ) = Q( OFF1 ) + SIGMA
         ELSE
            Q( OFF1 ) = QQ( OFF1 ) + SIGMA
         END IF
         N = 0
         GO TO 20
*
*     2x2 Solver
*
      ELSE IF( N.EQ.2 ) THEN
         IF( IPHASE.EQ.2 ) THEN
            QEMAX = MAX( QQ( N+OFF ), QQ( N-1+OFF ), EE( N-1+OFF ) )
            IF( QEMAX.NE.ZERO ) THEN
               IF( QEMAX.EQ.QQ( N-1+OFF ) ) THEN
                  XX = HALF*( QQ( N+OFF )+QQ( N-1+OFF )+EE( N-1+OFF )+
     $                 QEMAX*SQRT( ( ( QQ( N+OFF )-QQ( N-1+OFF )+EE( N-
     $                 1+OFF ) ) / QEMAX )**2+FOUR*EE( OFF+N-1 ) /
     $                 QEMAX ) )
               ELSE IF( QEMAX.EQ.QQ( N+OFF ) ) THEN
                  XX = HALF*( QQ( N+OFF )+QQ( N-1+OFF )+EE( N-1+OFF )+
     $                 QEMAX*SQRT( ( ( QQ( N-1+OFF )-QQ( N+OFF )+EE( N-
     $                 1+OFF ) ) / QEMAX )**2+FOUR*EE( N-1+OFF ) /
     $                 QEMAX ) )
               ELSE
                  XX = HALF*( QQ( N+OFF )+QQ( N-1+OFF )+EE( N-1+OFF )+
     $                 QEMAX*SQRT( ( ( QQ( N+OFF )-QQ( N-1+OFF )+EE( N-
     $                 1+OFF ) ) / QEMAX )**2+FOUR*QQ( N-1+OFF ) /
     $                 QEMAX ) )
               END IF
               YY = ( MAX( QQ( N+OFF ), QQ( N-1+OFF ) ) / XX )*
     $              MIN( QQ( N+OFF ), QQ( N-1+OFF ) )
            ELSE
               XX = ZERO
               YY = ZERO
            END IF
         ELSE
            QEMAX = MAX( Q( N+OFF ), Q( N-1+OFF ), E( N-1+OFF ) )
            IF( QEMAX.NE.ZERO ) THEN
               IF( QEMAX.EQ.Q( N-1+OFF ) ) THEN
                  XX = HALF*( Q( N+OFF )+Q( N-1+OFF )+E( N-1+OFF )+
     $                 QEMAX*SQRT( ( ( Q( N+OFF )-Q( N-1+OFF )+E( N-1+
     $                 OFF ) ) / QEMAX )**2+FOUR*E( N-1+OFF ) /
     $                 QEMAX ) )
               ELSE IF( QEMAX.EQ.Q( N+OFF ) ) THEN
                  XX = HALF*( Q( N+OFF )+Q( N-1+OFF )+E( N-1+OFF )+
     $                 QEMAX*SQRT( ( ( Q( N-1+OFF )-Q( N+OFF )+E( N-1+
     $                 OFF ) ) / QEMAX )**2+FOUR*E( N-1+OFF ) /
     $                 QEMAX ) )
               ELSE
                  XX = HALF*( Q( N+OFF )+Q( N-1+OFF )+E( N-1+OFF )+
     $                 QEMAX*SQRT( ( ( Q( N+OFF )-Q( N-1+OFF )+E( N-1+
     $                 OFF ) ) / QEMAX )**2+FOUR*Q( N-1+OFF ) /
     $                 QEMAX ) )
               END IF
               YY = ( MAX( Q( N+OFF ), Q( N-1+OFF ) ) / XX )*
     $              MIN( Q( N+OFF ), Q( N-1+OFF ) )
            ELSE
               XX = ZERO
               YY = ZERO
            END IF
         END IF
         Q( N-1+OFF ) = SIGMA + XX
         Q( N+OFF ) = YY + SIGMA
         N = 0
         GO TO 20
      END IF
      CALL DLASQ3( N, Q( OFF1 ), E( OFF1 ), QQ( OFF1 ), EE( OFF1 ), SUP,
     $             SIGMA, KEND, OFF, IPHASE, ICONV, EPS, TOL2, SMALL2 )
      IF( SUP.LT.ZERO ) THEN
         INFO = N + OFF
         RETURN
      END IF
      OFF1 = OFF + 1
      GO TO 20
*
*     End of DLASQ2
*
      END
