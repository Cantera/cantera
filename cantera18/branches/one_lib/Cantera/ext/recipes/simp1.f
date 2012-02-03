      SUBROUTINE SIMP1(A,MP,NP,MM,LL,NLL,IABF,KP,BMAX)
      implicit double precision (a-h,o-z)
      DIMENSION A(MP,NP),LL(NP)
      KP=LL(1)
      BMAX=A(MM+1,KP+1)
      IF(NLL.LT.2)RETURN
      DO 11 K=2,NLL
        IF(IABF.EQ.0)THEN
          TEST=A(MM+1,LL(K)+1)-BMAX
        ELSE
          TEST=ABS(A(MM+1,LL(K)+1))-ABS(BMAX)
        ENDIF
        IF(TEST.GT.0.)THEN
          BMAX=A(MM+1,LL(K)+1)
          KP=LL(K)
        ENDIF
11    CONTINUE
      RETURN
      END
