      SUBROUTINE SIMP2(A,M,N,MP,NP,L2,NL2,IP,KP,Q1)
      implicit double precision (a-h,o-z)
      PARAMETER (EPS=1.d-6)
      DIMENSION A(MP,NP),L2(MP)
      IP=0
      IF(NL2.LT.1)RETURN
      DO 11 I=1,NL2
        IF(A(L2(I)+1,KP+1).LT.-EPS)GO TO 2
11    CONTINUE
      RETURN
2     Q1=-A(L2(I)+1,1)/A(L2(I)+1,KP+1)
      IP=L2(I)
      IF(I+1.GT.NL2)RETURN
      DO 13 I=I+1,NL2
        II=L2(I)
        IF(A(II+1,KP+1).LT.-EPS)THEN
          Q=-A(II+1,1)/A(II+1,KP+1)
          IF(Q.LT.Q1)THEN
            IP=II
            Q1=Q
          ELSE IF (Q.EQ.Q1) THEN
            DO 12 K=1,N
              QP=-A(IP+1,K+1)/A(IP+1,KP+1)
              Q0=-A(II+1,K+1)/A(II+1,KP+1)
              IF(Q0.NE.QP)GO TO 6
12          CONTINUE
6           IF(Q0.LT.QP)IP=II
          ENDIF
        ENDIF
13    CONTINUE
      RETURN
      END
