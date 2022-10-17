MODULE FUNZIONI

USE CONSTANTS

IMPLICIT NONE

contains

SUBROUTINE STAR_EMIS(ct,Phi)
REAL*8 cp,sp,ct,st,Phi, kxe,kye,kze 
  CALL RAN_DIR(kxe,kye,kze)
  sp = sin(Phi)  ; cp = cos(Phi) 
  if(ct.GT.0.9999999999999) ct = 0.9999999999*ct/abs(ct)
  st = sqrt(1.D0-ct**2)
  kx = kze*st*cp + kxe*ct*cp - kye*sp
  ky = kze*st*sp + kxe*ct*sp + kye*cp
  kz = kze*ct - kxe*st
END	SUBROUTINE STAR_EMIS

SUBROUTINE RAN_DIR(kxe,kye,kze)
REAL*8 kxe,kye,kze,st,ct,Phi 
   Phi = 2*Pi*(1.D0 - 2*Randomu( ))
   ct = Randomu( )**Beam ; st = sqrt(1.- ct**2)
   kxe = st*cos(phi) ; kye = st*sin(phi) ; kze = ct
END	SUBROUTINE RAN_DIR

SUBROUTINE PATCHES( )
REAL*8 DPhi, DTheta
INTEGER i
DPhi = 2*Pi/Nphi  ; DTheta = 2.D0/NTheta
DO i=1, Nphi
  PatPhi(i) = (i-1)*Dphi+Dphi/2
ENDDO
DO i=1, NTheta
  PatTh(i) = 1.D0 - (i-1)*DTheta - DTheta/2
ENDDO
END SUBROUTINE PATCHES

REAL*8 FUNCTION BBODY( )
INTEGER i
REAL*8 u,res, q,s
u = Randomu( ) ; u = u*Randomu( )
res = -log(u*Randomu( ))
q = 1.202057D0*Randomu( )
BBODY = res
if (q.LE.1.D0) RETURN
s = 0.D0
do i=1, 30
    s = s + 1.D0/DBLE(i)**3 ; if(q.LT.s) EXIT
enddo
BBODY = res/i
END FUNCTION BBODY 

REAL*8 FUNCTION Randomu()
REAL*8  x
CALL RANDOM_NUMBER(x)
Randomu = x  
END  FUNCTION Randomu

REAl*8 FUNCTION Fe(beta)
REAL*8 beta, beta2, gamma,gammaT, gammaR 
   beta2 = beta**2 ; Fe=0.D0 
   if(beta2 .GT. 0.9999) RETURN
   if (beta2 .LT. 1.D-2) then
     gamma = 1.D0+beta2/2+0.375D0*beta2**2+0.3125D0*beta2**3
   else 
     gamma = 1.D0/sqrt(1.D0-beta2)
   endif
   gammaR = gamma*GammaB*(1.D0 - BetaB*beta)
   gammaT = (gammaR-1.D0)/Tel
   if (gammaT .GT. 80.D0) RETURN
   Fe =  exp(-gammaT)/Bes1/2.D0/GammaB
END FUNCTION Fe

REAL*8 FUNCTION BESK2N(X)
   REAL*8 X
   BESK2N = BESSK0N(X) + 2.D0*BESSK1N(X)/X
END FUNCTION BESK2N

REAL*8 FUNCTION BESSK0N(X)
      REAL*8 P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
	  REAL*8 X,Y 
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,-0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK0N = exp(x)*(-LOG(X/2.0)*BESSI0(X) + P1 + Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=(2.0/X)
        BESSK0N=(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))/SQRT(X)
      ENDIF
END FUNCTION BESSK0N
    
REAL*8 FUNCTION BESSK1N(X)
      REAL*8 P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7
	  REAL*8 X,Y 
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,-0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK1N = exp(x)*( LOG(X/2.0)*BESSI1(X)+(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))/X)
      ELSE
        Y=2.0/X
        BESSK1N= (Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))/SQRT(X)
      ENDIF
 END FUNCTION BESSK1N

 REAL*8 FUNCTION BESSI0(X)
      REAL*8 P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
	  REAL*8 X,Y, AX
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D0,0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,&
         0.2635537D-1,-0.1647633D-1,0.392377D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI0 = P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=ABS(X) ; Y=3.75/AX
        BESSI0 = (EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
 END FUNCTION BESSI0

 REAL*8 FUNCTION BESSI1(X)
      REAL*8 P1,P2,P3,P4,P5,P6,P7,Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
	  REAL*8 X,Y, AX
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,-0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,	 &
          -0.2895312D-1,0.1787654D-1,-0.420059D-2/
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X) ; Y=3.75/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
 END FUNCTION BESSI1

 SUBROUTINE locate_(xx,n,x,j)
     INTEGER j,n
     REAL*8 x,xx(n)
     INTEGER jl,jm,ju
     jl=0 ; ju=n+1
10   if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
        goto 10
     endif
     if(x.eq.xx(1))then
        j=1
     else if(x.eq.xx(n))then
        j=n-1
     else
        j=jl
     endif
 END SUBROUTINE locate_

REAL*8 FUNCTION INTERPOL(V,X,N,U)
   integer s,N
   REAL*8 V(N), X(N),U
   CALL LOCATE(X,N,U,s)
   INTERPOL = (u-x(s))*(v(s+1)-v(s))/(x(s+1) - x(s)) + v(s)
END FUNCTION INTERPOL


SUBROUTINE locate(xx,n,x,j)
    INTEGER j,n
    REAL*8 x,xx(n)
    if( xx(n) .ge. xx(1) ) then
	  call locatec(xx,n,x,j) 
	else
      call located(xx,n,x,j)
    endif
    j = max(1,j) 
END SUBROUTINE locate

SUBROUTINE located(xx,n,x,j)
    INTEGER j,n
    REAL*8 x,xx(n)
    INTEGER jl,jm,ju
    jl=0 ; ju=n+1
10  if(ju-jl.gt.1)then
      jm=(ju+jl)/2
      if( x.lt.xx(jm) )  then
        jl=jm
      else
        ju=jm
      endif
      goto 10
    endif
    if(x.eq.xx(1))then
      j=1
    else if(x.eq.xx(n))then
      j=n-1
    else
      j=jl
    endif
END SUBROUTINE located

SUBROUTINE locatec(xx,n,x,j)
    INTEGER j,n
    REAL*8 x,xx(n)
    INTEGER jl,jm,ju
    jl=0 ; ju=n+1
10  if(ju-jl.gt.1)then
      jm=(ju+jl)/2
      if( x.ge.xx(jm) )  then
        jl=jm
      else
        ju=jm
      endif
      goto 10
    endif
    if(x.eq.xx(1))then
      j=1
    else if(x.eq.xx(n))then
      j=n-1
    else
      j=jl
    endif
END SUBROUTINE locatec

REAL FUNCTION BINTERPOLA(Z,X,Y,M,N,x0,y0)
INTEGER M,N, i,j
REAL*8 X(M), Y(N), Z(M,N), p1,p2,p, a,a1,b,b1,x0,y0

  CALL LOCATE(X,M,x0,i) 
  i = MAX(1,i) ; i = MIN(i,M-1)
  a = (x0-X(i))/(X(i+1)-X(i)) ; a1 = 1.D0 - a
  if (y0 .LT. Y(1)) then 
    p1 = (Z(i,2) - Z(i,1))/(Y(2)-Y(1))
    p2 = (Z(i+1,2) - Z(i+1,1))/(Y(2)-Y(1))
    p = a*p2 + a1*p1
    BINTERPOLA = a1*Z(i,1) + a*Z(i+1,1) - p*(Y(1)-y0)
    RETURN
  endif
  if(y0 .LT. Y(N)) then 
    CALL LOCATE(Y,N,y0,j) 
    b = (y0-Y(j))/(Y(j+1)-Y(j)) ; b1 = 1.D0 - b
    BINTERPOLA = b1*( a1*Z(i,j)+a*Z(i+1,j) ) + b*( a1*Z(i,j+1)+a*Z(i+1,j+1) )
  else 
    N = N-1
    p1 = (Z(i,N) - Z(i,N-1))/(Y(N)-Y(N-1))
    p2 = (Z(i+1,N) - Z(i+1,N-1))/(Y(n)-Y(N-1))
    p = a*p2 + a1*p1
    BINTERPOLA = a1*Z(i,N) + a*Z(i+1,N) - p*(Y(N)-y0)
  endif
RETURN
END	FUNCTION BINTERPOLA

END MODULE FUNZIONI
