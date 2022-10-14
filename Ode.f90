MODULE ODE
contains

SUBROUTINE MERSON1(X,XEND,Y,ACC,H,HMIN,IFAIL,DIFF)
IMPLICIT NONE
REAL*8 X,XEND,Y,ACC,H,HMIN,DIFF
INTEGER ISWH
REAL*8 YZ,A,B,F,W,Z,ZEND,ZMIN,BCC,S,HSW,COF,RZERO,HT
LOGICAL IFAIL
DATA RZERO/1.E-6/ 
EXTERNAL DIFF
!
    IFAIL = .FALSE. ; W=Y ; Z=X ; ZEND=XEND ; BCC=ACC ; ZMIN=HMIN ; S=H ;ISWH=0
    2 HSW = S ; COF = ZEND - Z
      IF(ABS(S).LT.ABS(COF))  GO TO 8
      S = COF
      IF(ABS(COF/HSW).LT.RZERO)  GO TO 50
      ISWH = 1
    8 YZ = W
      HT = S/3 ; F = DIFF(Z,W) ; Z = Z + HT ; A = HT*F ; W = A + YZ
      F = DIFF(Z,W) ; A = 0.5D0*A ; W = 0.5D0*HT*F+A+YZ      
	  F = DIFF(Z,W); Z = Z + 0.5D0*HT ; B = 4.5D0*HT*F; W = 0.25D0*B+0.75D0*A+YZ
	  F = DIFF(Z,W) ; Z = Z + 0.5*S ; A = 2*HT*F+A ; W = 3*A-B+YZ
	  F = DIFF(Z,W) ; B = -0.5*HT*F-B +2*A ; W = W-B ; A = ABS(5*BCC*W) ; B = ABS(B)
      IF( ABS(W).GT.RZERO .AND. B .GT. A ) goto 60
      IF(ISWH.EQ.1)  GO TO 50
      IF(B.GT.0.03125*A)  GO TO 2
      S = S + S ; GO TO 2
   50 H = HSW ;  X = Z ; Y = W
      RETURN
   60 COF = S/2
      IF(ABS(COF).GE.ZMIN)   GO TO 80
      goto 84
      S = ZMIN ; IF(HSW.LT.0.)   S = -S
      IF(ISWH.EQ.1) GO TO 50
      GO TO 2
   80 W = YZ ; Z = Z-S ;  S = COF ; ISWH = 0
      GO TO 2
   84 CONTINUE
      IFAIL = .TRUE. ; GO TO 50
      END SUBROUTINE MERSON1
END MODULE ODE
