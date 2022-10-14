MODULE POLAR
IMPLICIT NONE

PRIVATE
REAL*8, PARAMETER:: Pi=3.14159265358979324
REAL*8 r0,x0,y0,z0, Vk(3), Vp(3)
REAL*8 xp(3),yp(3),zp(3)
REAL*8 u(3),v(3),w(3)
REAL*8:: Xmax=1.D3
REAL*8  Fac,alpha,Ca,Sa, EkeV
REAL*8::  ACC=0.01, H=0.1, HMIN = 0.001
INTEGER pol

PUBLIC:: SUBPOL

CONTAINS

SUBROUTINE SUBPOL(r,Ct,Phi, kx,ky,kz,QS,VS,US,spol,EkV)  
INTEGER spol
REAL*8 r,Theta,Phi, kx,ky,kz 
REAL*8 QS,VS,US 
REAL*8 Ct,Vb(3),Babs, EkV
REAL*8 tti
Theta = acos(Ct)
EkeV = EkV ; pol = -spol
Vk(1) = kx ; Vk(2) = ky ; Vk(3) = kz
Vp(1) = r*SIN(Theta)*COS(Phi) 
Vp(2) = r*SIN(Theta)*SIN(Phi)
Vp(3) = r*Ct ! COS(Theta)
x0=Vp(1) ; y0=Vp(2) ; z0=Vp(3) ; r0=r
CALL START( )
CALL CampoB(Vb,Babs)
CALL BASI(Vb)
CALL STOKES(QS, VS, US)
END SUBROUTINE SUBPOL

! *********************  MAIN STOKES **************

 	  SUBROUTINE STOKES(Qf,Vf,Uf)
	  IMPLICIT NONE
      INTEGER,PARAMETER:: N=3, Npt=2000
      INTEGER i  ! ,Nnew 
      REAL*8 X,XEND,Y(3)
      REAL*8 Qf,Uf,Vf,Q1,ds 
      INTEGER:: JTEST=0
	  LOGICAL IFAIL
!
      CALL ANGOLO( )
      Uf = 0.D0 ; Vf = 0.D0
      Qf = pol*sqrt(1.D0-Uf**2-Vf**2)
!
	  Q1 = Qf
      Qf = Q1*COS(2*alpha) + Uf*SIN(2*alpha)
      Uf = -Q1*SIN(2*alpha) + Uf*COS(2*alpha)        

 !
 5    continue
 !    return
      if(r0.GT.500.) return
	  X = 0.D0 ; Y(1) = Qf ; Y(2) = Uf ; Y(3) = Vf  
!	  Zn(1) = X ; Q(1) = Y(1) ; U(1) = Y(2) ; V(1) = Y(3)  
      ds = 0.5D0
	  DO i = 1, Npt-1
		if(X.GT.Xmax) exit 
        Xend = X + ds  
		CALL MERSON(X,Xend,Y,N,ACC,H,HMIN,JTEST,IFAIL,SDIFS)
        X = Xend
        ds = X	 
        Qf = Y(1) ; Uf = Y(2) ; Vf = Y(3) 
	  ENDDO
!      Nnew = i-1 
	  END SUBROUTINE STOKES

      SUBROUTINE SDIFS(s,Y,DY)
	  REAL*8 s,Y(3),DY(3),Vb(3), Bn(3)
	  REAL*8 B11,Babs,bx,by, delta, Mc,Nc,Pc
      CALL PATH(s)
      CALL CampoB(Vb,Babs) ; B11 = Babs/1.0D11
	  CALL BCOMP(Vb,Bn) ! , flag)
	  bx = Bn(1) ; by = Bn(2) ;	delta = 1.D-10*B11**2 
	  CALL PARAMETRI(bx,by,delta, Mc,Nc,Pc)
      DY(1) = -2*Pc*Y(3)
      DY(3) =  2*Pc*Y(1) + (Nc-Mc)*Y(2)
      DY(2) = -(Nc - Mc)*Y(3)
      Fac = 3.45D3*EkeV*B11**2
      DY = Fac*DY
	  END SUBROUTINE SDIFS

SUBROUTINE START( )
INTEGER i
REAL*8 s, Babs,B11,Vb(3),r0p,errea,Baipoli11
!
s = 0.D0
r0p=r0
do i = 1, 5000
  CALL PATH(s) 
  CALL CampoB(Vb,Babs) ; B11 = Babs/1.0D11
  Fac = 3.45D3*EkeV*B11**2
  if(Fac.LT.50.D0) exit
  s = s+1.D0/10.D0
enddo
CALL PATH(s) ; x0 = Vp(1) ; y0 = Vp(2) ; z0 = Vp(3)
r0 = sqrt(x0**2+y0**2+z0**2)
END	SUBROUTINE START

SUBROUTINE PATH(s)
REAL*8 s  
Vp(1) = x0+Vk(1)*s ; Vp(2) = y0+Vk(2)*s ; Vp(3) = z0+Vk(3)*s ; 
END SUBROUTINE PATH

SUBROUTINE BASEX(Vb)
REAL*8 Vb(3),V(3)
CALL CROSSP(Vk, Vb, V)
CALL NORMA(V) ; xp = V
CALL CROSSP(Vk, xp, yp)
zp = vk
END	SUBROUTINE BASEX

SUBROUTINE BASEU( )
REAL*8 UNITA(3),up(3)                   
DATA UNITA/0.D0,0.D0,1.D0/
CALL CROSSP(Vk,UNITA,up)
CALL NORMA(up) ; u = up
CALL CROSSP(Vk,u,v) ; w = Vk
END	SUBROUTINE BASEU

SUBROUTINE BASI(Vb )
REAL*8 Vb(3)
CALL BASEX(Vb) ; CALL BASEU( )
END	SUBROUTINE BASI

SUBROUTINE BCOMPX(b,bn)
REAL*8 b(3), bn(3)
bn(1) = b(1)*xp(1) + b(2)*xp(2) + b(3)*xp(3)
bn(2) = b(1)*yp(1) + b(2)*yp(2) + b(3)*yp(3)
bn(3) = b(1)*zp(1) + b(2)*zp(2) + b(3)*zp(3)
END	SUBROUTINE BCOMPX

SUBROUTINE BCOMPU(b,bn)
REAL*8 b(3),bn(3)
bn(1) = b(1)*u(1) + b(2)*u(2) + b(3)*u(3)
bn(2) = b(1)*v(1) + b(2)*v(2) + b(3)*v(3)
bn(3) = b(1)*w(1) + b(2)*w(2) + b(3)*w(3)
END	SUBROUTINE BCOMPU

SUBROUTINE BCOMP(b,bn)
REAL*8 b(3), bn(3) 
   CALL BCOMPU(b,bn)
END	SUBROUTINE BCOMP

SUBROUTINE CROSSP(V1,V2,PV)
REAL*8 V1(3), V2(3), PV(3)
PV(1) =  V1(2)*V2(3)-V2(2)*V1(3)
PV(2) =  V1(3)*V2(1)-V2(3)*V1(1)
PV(3) =  V1(1)*V2(2)-V2(1)*V1(2) 
END SUBROUTINE CROSSP


SUBROUTINE ANGOLO( )
ca = u(1)*xp(1) + u(2)*xp(2) + u(3)*xp(3)
if(Vk(1)*(xp(2)*u(3) - xp(3)*u(2)).LT.0.D0) then
   alpha = -ACOS(Ca)  
else
   alpha = ACOS(Ca)  
endif
Ca = COS(alpha) ; Sa  = SIN(alpha)
END SUBROUTINE ANGOLO 


SUBROUTINE PARAMETRI(bx,by,delta, Mc,Nc,Pc)
REAL*8 bx,by,delta, Mc, Nc, Pc
REAL*8 muxx,muyy,muxy,exx,eyy,exy,BQ,Q
!  
muxx = (1.D0 - 2*delta) - 4*delta*bx**2
muyy = (1.D0 - 2*delta) - 4*delta*by**2
muxy = -4*delta*bx*by
exx = (1.D0 - 2*delta) - 7*delta*bx**2
eyy = (1.D0 - 2*delta) - 7*delta*by**2
exy = -7*delta*bx*by
Bq = bx**2*by**2 ; Q = muxx*muyy - muxy*muxy
Mc = (7*bx**2 + 4*by**2)*muxx - 12*delta*Bq
Pc = bx*by*(3*muxx - 4*delta*(7*by**2 + 4*bx**2))
Nc = (7*by**2 + 4*bx**2)*muyy - 12*delta*Bq
Mc = Mc/Q ; Nc = Nc/Q ; Pc = Pc/Q
END	SUBROUTINE PARAMETRI

SUBROUTINE NORMA(V)
REAL*8 V(3)
V = V/SQRT(SUM(V**2))
END SUBROUTINE NORMA

SUBROUTINE Trasforma(Theta,Phi, Dr,Dt,Dp, Vc)
REAL*8 Theta,Phi, Dr,Dt,Dp, Vc(3) 
  Vc(1) =  cos(phi)*(sin(theta)*Dr + cos(theta)*Dt) - sin(phi)*Dp
  Vc(2) =  sin(phi)*(sin(theta)*Dr + cos(theta)*Dt) + cos(phi)*Dp
  Vc(3)  = cos(theta)*Dr - sin(theta)*Dt
end SUBROUTINE Trasforma


! *******************  MERSON ********************

      SUBROUTINE MERSON(X,XEND,Y,N,ACC,H,HMIN,ITEST,IFAIL,DIFF)
      INTEGER K,N,ITEST,ISWH
	  REAL*8 YZ(N),A(N),B(N),F(N),W(N)
      REAL*8  Z,ZEND,ZMIN,S,HSW,COF,HT,RZERO
      DATA RZERO/1.D-6/    !!
      REAL*8 X,XEND,Y(N),ACC,H,HMIN
      LOGICAL IFAIL
	  EXTERNAL DIFF
!
      IFAIL = .TRUE.
      W = Y ; Z = X ; ZEND = XEND ;	ZMIN = HMIN ; S = H ; ISWH = 0
    2 HSW = S ; COF = ZEND - Z
      IF(ABS(S).LT.ABS(COF))  GO TO 8
      S = COF ; IF(ABS(COF/HSW).LT.RZERO)  GO TO 50
      ISWH = 1
    8 CONTINUE 
      YZ = W
      HT = S/3 ; CALL DIFF(Z,W,F) ; Z = Z + HT
      A = HT*F ; W = A + YZ
      CALL DIFF(Z,W,F)
      A = 0.5D0*A ; W = 0.5D0*HT*F + A + YZ
      CALL DIFF(Z,W,F)
      Z = Z + 0.5D0*HT
      B = 4.5D0*HT*F ; W = 0.25D0*B + 0.75D0*A + YZ
      CALL DIFF(Z,W,F)
      Z = Z + 0.5D0*S
      A = 2*HT*F + A ; W = 3*A - B + YZ
      CALL DIFF(Z,W,F) 
      B = -0.5D0*HT*F - B + 2*A ; W = W - B
      A = ABS(5*ACC*W) ; B = ABS(B)
      DO K = 1,N
         IF(ABS(W(K)).LE.RZERO)  EXIT 
         IF(B(K).GT.A(K))  GO TO 60
	  ENDDO
      IF(ISWH.EQ.1)  GO TO 50
      DO   K = 1,N
        IF(B(K).GT.0.03125D0*A(K))  GO TO 2
	  ENDDO
      S = S + S ; GO TO 2
   50 H = HSW ;  X = Z
      Y = W
      RETURN
   60 COF = S/2
      IF(ABS(COF).GE.ZMIN)   GO TO 80
      IF(ITEST.EQ.0)   GO TO 84
      S = ZMIN
      IF(HSW.LT.0.D0)   S = -S
      IF(ISWH.EQ.1)   GO TO 50
      GO TO 2
   80 W = YZ 
      Z = Z - S ;  S = COF ; ISWH = 0
      GO TO 2
   84 PRINT 88, S, ZMIN, Z, Fac
      IFAIL = .FALSE. ; GO TO 50
   88 FORMAT(//1X,'** MERSON ERROR **',13E10.2//)
      END SUBROUTINE MERSON

! ******************************************

SUBROUTINE CampoB(Vb,Babs)
USE MAGFIELD
REAL*8 Vb(3), Babs, bpol(3)
REAL*8 r,x,y,z,theta,phi, Ct, Nebeta
x = Vp(1) ; y = Vp(2) ; z = Vp(3) ; r = sqrt(x**2+ y**2+z**2)
Phi = atan2(y,x)
Ct = z/r 
CALL b_fieldm(r,Ct, bpol,Babs,Nebeta)
theta = acos(Ct)
CALL Trasforma(Theta,Phi,bpol(1),bpol(2),bpol(3), Vb)
END SUBROUTINE CampoB

SUBROUTINE CampoBstart(Vb,Babs,erre)
USE MAGFIELD
REAL*8 Vb(3), Babs, bpol(3)
REAL*8 r,x,y,z,theta,phi, Ct, Nebeta,erre
x = Vp(1) ; y = Vp(2) ; z = Vp(3) ; r = sqrt(x**2+ y**2+z**2)
erre=r
Phi = atan2(y,x)
Ct = z/r 
CALL b_fieldm(r,Ct, bpol,Babs,Nebeta)
theta = acos(Ct)
CALL Trasforma(Theta,Phi,bpol(1),bpol(2),bpol(3), Vb)
END SUBROUTINE CampoBstart

END MODULE POLAR

