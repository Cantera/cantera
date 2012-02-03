      SUBROUTINE SIMP3(A,MP,NP,I1,K1,IP,KP)
      implicit double precision (a-h,o-z)
      DIMENSION A(MP,NP)
      PIV=1.d0/A(IP+1,KP+1)
      IF(I1.GE.0)THEN
        DO 12 II=1,I1+1
          IF(II-1.NE.IP)THEN
            A(II,KP+1)=A(II,KP+1)*PIV
            DO 11 KK=1,K1+1
              IF(KK-1.NE.KP)THEN
                A(II,KK)=A(II,KK)-A(IP+1,KK)*A(II,KP+1)
              ENDIF
11          CONTINUE
          ENDIF
12      CONTINUE
      ENDIF
      DO 13 KK=1,K1+1
        IF(KK-1.NE.KP)A(IP+1,KK)=-A(IP+1,KK)*PIV
13    CONTINUE
      A(IP+1,KP+1)=PIV
      RETURN
      END
