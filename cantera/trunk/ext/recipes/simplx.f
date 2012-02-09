      SUBROUTINE SIMPLX(A,M,N,MP,NP,M1,M2,M3,ICASE,IZROV,IPOSV)
      implicit double precision (a-h,o-z)
      PARAMETER(MMAX=1000,EPS=1.d-6)
      DIMENSION A(MP,NP),IZROV(N),IPOSV(M),L1(MMAX),L2(MMAX),L3(MMAX)
      IF(M.NE.M1+M2+M3)PAUSE 'Bad input constraint counts.'
      NL1=N
      DO 11 K=1,N
        L1(K)=K
        IZROV(K)=K
11    CONTINUE
      NL2=M
      DO 12 I=1,M
        IF(A(I+1,1).LT.0.d0) then
c          write(*,*) 'The A matrix input to SIMPLX is invalid.'
          PAUSE 'Bad input tableau.'
        end if
        L2(I)=I
        IPOSV(I)=N+I
12    CONTINUE
      DO 13 I=1,M2
        L3(I)=1
13    CONTINUE
      IR=0
      IF(M2+M3.EQ.0)GO TO 30
      IR=1
      DO 15 K=1,N+1
        Q1=0.
        DO 14 I=M1+1,M
          Q1=Q1+A(I+1,K)
14      CONTINUE
        A(M+2,K)=-Q1
15    CONTINUE
10    CALL SIMP1(A,MP,NP,M+1,L1,NL1,0,KP,BMAX)
      IF(BMAX.LE.EPS.AND.A(M+2,1).LT.-EPS)THEN
        ICASE=-1
        RETURN
      ELSE IF(BMAX.LE.EPS.AND.A(M+2,1).LE.EPS)THEN
        M12=M1+M2+1
        IF(M12.LE.M)THEN
          DO 16 IP=M12,M
            IF(IPOSV(IP).EQ.IP+N)THEN
              CALL SIMP1(A,MP,NP,IP,L1,NL1,1,KP,BMAX)
              IF(BMAX.GT.0.)GO TO 1
            ENDIF
16        CONTINUE
        ENDIF
        IR=0
        M12=M12-1
        IF(M1+1.GT.M12)GO TO 30
        DO 18 I=M1+1,M12
          IF(L3(I-M1).EQ.1)THEN
            DO 17 K=1,N+1
              A(I+1,K)=-A(I+1,K)
17          CONTINUE
          ENDIF
18      CONTINUE
        GO TO 30
      ENDIF
      CALL SIMP2(A,M,N,MP,NP,L2,NL2,IP,KP,Q1)
      IF(IP.EQ.0)THEN
        ICASE=-1
        RETURN
      ENDIF
1     CALL SIMP3(A,MP,NP,M+1,N,IP,KP)
      IF(IPOSV(IP).GE.N+M1+M2+1)THEN
        DO 19 K=1,NL1
          IF(L1(K).EQ.KP)GO TO 2
19      CONTINUE
2       NL1=NL1-1
        DO 21 IS=K,NL1
          L1(IS)=L1(IS+1)
21      CONTINUE
      ELSE
        IF(IPOSV(IP).LT.N+M1+1)GO TO 20
        KH=IPOSV(IP)-M1-N
        IF(L3(KH).EQ.0)GO TO 20
        L3(KH)=0
      ENDIF
      A(M+2,KP+1)=A(M+2,KP+1)+1.
      DO 22 I=1,M+2
        A(I,KP+1)=-A(I,KP+1)
22    CONTINUE
20    IS=IZROV(KP)
      IZROV(KP)=IPOSV(IP)
      IPOSV(IP)=IS
      IF(IR.NE.0)GO TO 10
30    CALL SIMP1(A,MP,NP,0,L1,NL1,0,KP,BMAX)
      IF(BMAX.LE.0.)THEN
        ICASE=0
        RETURN
      ENDIF
      CALL SIMP2(A,M,N,MP,NP,L2,NL2,IP,KP,Q1)
      IF(IP.EQ.0)THEN
        ICASE=1
        RETURN
      ENDIF
      CALL SIMP3(A,MP,NP,M,N,IP,KP)
      GO TO 20
      END
