MODULE INTEGRATE
USE PATH
USE FUNZIONI
IMPLICIT NONE
INTEGER Jmax, bind, Maxs, Nscatt, FlagLoss, FlagStep, NTScatt
REAL*8  Rad,s1,s2,Res, Rapp,Rapm
contains

SUBROUTINE MainInteg(Jma,ip,jp)
USE FILES
USE POLAR

INTEGER i, J,Jma, ind, Ien, ip,jp, spol
REAL Sec1,Sec2
REAL*8  rand, mu, EnkeV, Leps, ell, Nebeta, s, QS,VS,US
REAL*8  omega_c, eps,Sth,Cth,Phi,b(3),bmod,beta,betap,betam
REAL*8  r, Ctin, Phin, Cthe, Phie, kbx,kby,kbz, Ctph, Phiph, smin, sman, Cthzero,Phizero, bimp
CHARACTER*45 StrP
!
   Ctin = PatTh(ip) ; Phin = PatPhi(jp) ; Jmax = jma
   Counts = 0 ; Ven = 0 ; Nord = 0 ; Next = 0 ; Nlost = 0
   QSTOKES=0.D0 ; USTOKES=0.D0 ; VSTOKES=0.D0
   if(FlagP.GT.0) Strp = ' ---- Initial ordinary photons ---- NR'
   if(FlagP.LT.0) Strp = ' ---- Initial extraordinary photons ---- NR'
   if(kinv.LT.0) StrP =  ' ---- Initial Unpolarized photons ------ NR'
   if(Printout.GT.0) write(*,'(A)') Strp
   CALL TIMES(Sec1) ; NTscatt = 0 ; Maxs = 0
   DO j = 1, Jmax
      FlagP = kinv*FlagP ; Flag = FlagP
1     FlagLoss = 0 ; Nscatt = 0

      r = 1.D0 ; Cth0 = Ctin ; Phi0 = Phin ; Sth0=sqrt(1.D0-cth0**2)
	  X0=r*Sth0*cos(Phi0) ; Y0=r*Sth0*sin(Phi0) ; Z0=r*Cth0
      eps=Tstar*BBody( )
      omega=1.519255D18*eps
      CALL STAR_EMIS(Cth0,Phi0)
      DO i = 1, 1000
         rand = Randomu( ) 	   
		 IF (Flag.EQ.1) then  
           CALL TINTEGRATE(ell, rand,SDIF1)
	     else
           CALL TINTEGRATE(ell, rand,SDIF2)
         endif 

         CALL EPSMU(ell,mu,eps,Nebeta,omega_c,r,Cth,Phi,b,bmod)
		 Sth = sqrt(1.D0-Cth**2) ; X0=r*Sth*cos(Phi) ; Y0=r*Sth*sin(Phi) ; Z0=r*Cth
         if (FlagLoss .eq.1 ) goto 10
         if (FlagLoss .eq. -1 ) goto 9
!
         IF(S1+S2.LT.1.D-30)  goto 9
         CALL BetaPM(eps,mu,betap,betam,bind)
!
         CALL DECISIONS(betap,betam,beta,ind) 
		 CALL CARTESIAN(b(1),b(2),b(3),Cth,Phi, kbx,kby,kbz) 
		 s = SIGN(1.D0,beta) ; kbx = s*kbx ; kby=s*kby ; kbz=s*kbz
!
	     CALL VANGOLI(kbx,kby,kbz, Cthe,Phie) !  
         CALL SScattering(beta,Cthe,Phie,eps,mu,ind)
         Nscatt = Nscatt+1 ;  omega=eps*omega_c 
	  ENDDO
 9    NLost = NLost+1 ; goto 1
!
10    NTscatt = NTscatt+Nscatt 
      Enkev=omega/1.519255D18 ; Leps = log10(EnkeV) 
      Ien = NINT((Nin-1)*(Leps-Emin)/(Emax-Emin))+1
      if (Ien .LE. Nin .and. Ien .GT. 0) Ven(Ien)=Ven(Ien)+1
 	  CALL VANGOLI(kx,ky,kz, Ctph,Phiph) 
	  Maxs = MAX(Maxs,Nscatt)
	  if(Flag.EQ.1) then ; 
	    Nord=Nord+1 ;  spol=1
	  else ;
	    Next =Next+1 ; spol=-1
	  endif
!
smin=-(kx*r*sqrt(1.D0-Cth**2)*cos(Phi)+ky*r*sqrt(1.D0-Cth**2)*sin(Phi)+kz*r*Cth)
bimp=sqrt((kx*smin+r*sqrt(1.D0-Cth**2)*cos(Phi))**2+(ky*smin+r*sqrt(1.D0-Cth**2)*sin(Phi))**2+(kz*smin+r*Cth)**2)
Cthzero=(kz*smin+r*Cth)/bimp
Phizero=atan2((ky*smin+r*sqrt(1.D0-Cth**2)*sin(Phi)),(kx*smin+r*sqrt(1.D0-Cth**2)*cos(Phi)))
if(bimp.LE.1.D0) then 
   bimp=1.D0
   sman=smin*(1.D0-sqrt(1.D0+(1.D0-r**2)/(smin**2)))
   Cthzero=(kz*sman+r*Cth)
   Phizero=atan2((ky*sman+r*sqrt(1.D0-Cth**2)*sin(Phi)),(kx*sman+r*sqrt(1.D0-Cth**2)*cos(Phi)))
endif
r=bimp ; Cth=Cthzero ; Phi=Phizero
      CALL SUBPOL(r,Cth,Phi, kx,ky,kz,QS,VS,US, spol, EnkeV) 
      CALL photon_collect(Ctph,Phiph,Leps,QS,US,VS)
!
      if(Printout.GE.500) then
	     if(MOD(j,Printout) .eq. 0) write(*,*) j  
	  endif 
   ENDDO
100 continue 
    if(Printout.GT.0) then
       CALL TIMES(Sec2) 
       write(*,'(2X,A23,F10.2)') 'Elapsed time (seconds)=', Sec2-Sec1
       write(*,'(2X,A23,I10)')   'rate (photons/second) =', NINT(DBLE(j)/(Sec2-Sec1))
       write(*,'(2X,A23,F10.2)') 'mean # scatts/photon  =', REAL(NTscatt)/j
       write(*,'(2X,A26,I7)')    'Max #  scatts/photon  =   ', Maxs
       CALL WRAP( )     ! Print ratio ord/extr
       write(*,'(2X,A23,2I10)')  'Lost photons          =',  NLost
       write(*,'(2X,A15,I2,1H,,I2,$)')   'Completed Patch ', ip,jp
	endif
    CALL SALVA( ) 
    CALL SAVEDATA(ip,jp)
    END SUBROUTINE MainInteg 

SUBROUTINE WRAP( )
CHARACTER*15 SDum1,Sdum2
WRITE(Sdum1,'(1X,I10)') Nord ; Sdum1 = adjustL(Sdum1)
WRITE(Sdum2,'(1X,I10)') Next ; Sdum2 = adjustL(Sdum2)
Sdum1 = TRIM(Sdum1)//'/'//TRIM(Sdum2)
WRITE(Sdum2,'(1X,F5.1)')  REAL(Nord)/(REAL(Next)+0.1) ; Sdum2 = adjustL(Sdum2)
write(*,'(A)')  '  ordinary/extra photons=     '//TRIM(Sdum1)//'   (='//TRIM(Sdum2)//')'
END SUBROUTINE WRAP

SUBROUTINE TINTEGRATE(ell,ran,SDIF)
USE ODE
INTEGER, PARAMETER:: Nloop=3000
REAL*8, PARAMETER:: ACC = 1.D-2
REAL*8  dum, SDIF, r, ran,taus,ell,ellp,ellQ, Cth,Phi
REAL*8  Y,YP, YQ,S, H,Hmin, XEND, H0, Y0,ell0,YDum
INTEGER i
LOGICAL IFAIL 
EXTERNAL SDIF

    taus = -log(ran) ; ell=0.D0 ; Y=0.D0 ; flagStep=0 
    H0=0.9D0 ; Hmin=H0/1.D2 ; H = H0/10. ; FlagLoss = 0
    i = ESCAPE(ell)
    DO i = 0, Nloop
      ellp=ell ; YP=Y
      s = SDIF(ell,YDum) ; ell = ell+H0
      IF (ESCAPE(ell) .lt. 0) RETURN
      if(s .GT. 1.D-8) goto 10  
    ENDDO
    FlagLoss = 1 ; ell=0 ; RETURN
10  Y0 = YP ; ell0 = ellp
	DO WHILE (H0 .GT. HMIN)
       ell=ellp ; Y=YP ; H = H0
       DO i = 1, Nloop
	     ellp=ell ; YP=Y ; XEND = ell+H
         CALL MERSON1(ell,XEND,Y,ACC,H,HMIN,IFAIL,SDIF)
		 IF(.NOT.IFAIL) H=MIN(2*H,H0)
		 ellQ = ell ; YQ = Y
         IF (ESCAPE(ell) .lt. 0) RETURN
         if(Y .GT. taus) goto 20
       ENDDO
	   FlagLoss = 1 ; RETURN
20     H0 = H0/10
    ENDDO
    if(abs(YP-YQ) .GT. 1.D-4) ell=ellp+(taus-YP)*(ellQ-ellp)/(YQ-YP)
    CALL POSITION(ell,r,CTh,Phi)
    if (r .LT. 0.99)  FlagLoss=-1 
	if(Flag.EQ.1) then
       dum = SDIF1(ell,Y)
	else
       dum = SDIF2(ell,Y)
    endif
END SUBROUTINE TINTEGRATE

INTEGER FUNCTION ESCAPE(ell)
    real*8 ell,mu,epsp, mup, test
    REAL*8 r, eps, Nebeta,omega_c, Ct, Phi, b(3), bmod
    if (flagStep .EQ. 0) then 
      flagStep = 1 ; ESCAPE = 1 ; FlagLoss = 0
      CALL EPSMU(ell, mu,eps,Nebeta,omega_c,r,Ct,Phi,b,bmod)
      epsp=eps ; mup=mu ; RETURN
    endif
    CALL EPSMU(ell, mu,eps,Nebeta,omega_c,r,Ct,Phi,b,bmod)
    ESCAPE = -1
    if (r .GT. 150.) then ; FlagLoss= 1 ; RETURN ; endif
    if (r .LT. 0.99) then ; FlagLoss=-1 ; RETURN ; endif
    ESCAPE = 1
    if(epsp-eps .GT. 0.D0) goto 100 
    Rad = 1.D0/eps**2 + mu**2 - 1.D0
    if(Rad .GE. 0.D0)  goto 100 
    test = abs(mup*epsp - mu*eps)
    if(test .LT. abs(epsp-eps)) then ; ESCAPE=-1 ; FlagLoss=1 ; endif
100 continue                     
    epsp = eps ; mup = mu
END FUNCTION ESCAPE

REAL*8 FUNCTION SDIF1(ell, Y)
    REAL*8  r, ell, mu, betap,betam, Nebeta, Y
    REAL*8  omega_c, eps, Ct, Phi, b(3), bmod
    SDIF1 = 0.D0 ;  Y=Y ; s1 = 0.D0 ; s2 = 0.D0
    CALL EPSMU(ell,mu,eps,Nebeta,omega_c,r,Ct,Phi,b,bmod)
    if(r .LT. 0.999) FlagLoss = -1
    CALL BetaPM(eps,mu,betap,betam,bind)
    if(bind .EQ. 0) RETURN
    Rapp = abs(mu-betap)/(1.D0-mu*betap)
    s1 = Rapp*Fe(betap)*nebeta/vmean
    Rapm = abs(mu-betam)/(1.D0 - mu*betam)
    s2 = Rapm*Fe(betam)*nebeta/vmean
    SDIF1 = omega_c*CSigma*(s1+s2)/omega**2 
 END FUNCTION SDIF1

 REAL*8 FUNCTION SDIF2(ell, Y)
    REAL*8  r, ell, mu, betap,betam,Nebeta,Y 
    REAL*8  omega_c, eps, Ct, Phi, b(3), bmod
    SDIF2 = 0.D0 ; Y = Y ; s1 = 0.D0 ; s2 = 0.D0
    CALL EPSMU(ell,mu,eps,Nebeta,omega_c,r,Ct,Phi,b,bmod)
    if(r .LT. 0.999) FlagLoss = -1
    CALL BetaPM(eps,mu,betap,betam,bind)
    if(bind .EQ. 0) RETURN
    if(abs(mu-betap) .LT. 1.D-2) RETURN
    Rapp = (1.D0-mu*betap)/abs(mu-betap)
    s1 = Rapp*Fe(betap)*nebeta/vmean
    if(abs(mu-betam) .LT. 1.D-2) RETURN
    Rapm = (1.D0-mu*betam)/abs(mu-betam)
    s2 = nebeta*Rapm*Fe(betam)/vmean
    SDIF2 = omega_c*CSigma*(s1+s2)/omega**2
 END FUNCTION SDIF2
 
 SUBROUTINE BetaPM(eps,mu,betap,betam, bind)
 REAL*8 eps, mu, betap, betam, e2,mu2
 integer bind
 e2 = eps**2 ; mu2 = mu**2 ; Rad = 1.D0/e2 + mu2 - 1.D0
 if (Rad .LE. 0.D0) then 
    bind = 0
    betap = e2*mu/(1.D0+e2*mu2)
    betam = betap
 else 
    bind = 2
    betap = e2*(mu + sqrt(Rad)/eps)/(1.D0 + e2*mu2)
    betam = e2*(mu - sqrt(Rad)/eps)/(1.D0 + e2*mu2)
	if(betap .GT.0.9999) bind = 0
 endif
 END SUBROUTINE BetaPM

 SUBROUTINE EPSMU(ell,mu,eps,Nebeta,omega_c, r,Ct,Phi,b,bmod)
 USE MAGFIELD
 REAL*8 omega_c, bmod, b(3), eps, ell
 REAL*8 r, Ct, Phi, Nebeta, mu, vx,vy,vz
 CALL POSITION(ell,r,Ct,Phi)
 CALL b_fieldm(r,Ct,b,bmod,Nebeta)
 omega_c = bmod*omega_p
 eps=omega/omega_c
 CALL CARTESIAN(b(1),b(2),b(3),Ct,Phi, vx,vy,vz)
 mu = vx*kx+vy*ky+vz*kz
 END SUBROUTINE EPSMU 

SUBROUTINE DECISIONS(betap,betam, beta,ind)
INTEGER ind
REAL*8  R12,Rbm,beta,betap,betam 
!
R12 = Randomu( )
if (Flag.EQ.1) then
  ind = 1
  if (R12 .GT. 0.25D0) then
      Flag = -1 ; ind = 2
  endif
else                        
  ind = 2
  if (R12 .GT. 0.75D0) then
     Flag = 1 ; ind = 1
  endif
endif
beta = betap
Rbm = Randomu( )  
if(Rbm .lt. s1/(s1+s2)) then ; beta=betap ; else ; beta=betam ; endif
END SUBROUTINE DECISIONS

SUBROUTINE SScattering(beta,Cthe,Phie,eps,mu,ind)
INTEGER ind
REAL*8 beta, Cthe, Phie, eps,mu, Phi,Cphi,Sphi, cth2, ct2, om
REAL*8 u, gam2, sthe, chi, cosom, ct,sdf, sinom	, st, cdf,df
    cth2 = cthe**2 ; if(cth2 .GT. 0.99999999) cth2 = 0.99999999
    gam2 = 1.D0/(1.D0-beta**2) ; sthe = sqrt(1.d0-cth2)
    chi = 2*Pi*Randomu( ); cosom = 1.D0-2*Randomu( )   	
	if (ind .EQ. 2)  goto 10
!
    u = cosom ; cosom = ABS(cosom)**0.33333333333333D0
    if (u .LT. 0.d0) cosom = -cosom 
 10 continue
    eps = gam2*eps*(1.D0-beta*mu)*(1.D0+beta*cosom)
    cosom = (cosom+beta)/(1.D0+beta*cosom)
	if(abs(cosom) .GT. 0.99999999) then
	 cosom = 0.9999999*cosom/abs(cosom)
	endif
    om = acos(cosom) ; sinom = sin(om)
    ct = cthe*cosom + sthe*sinom*cos(chi) 
    ct2 = ct**2 ; if(ct2 .GT. 0.99999999) ct2 = 0.99999999
    st = sqrt(1.D0-ct2)
    sdf = sin(chi)*sinom/st ; cdf=(cosom-cthe*ct)/(sthe*st) 
    if(abs(cdf).gt.0.9999999) then
      df = 0.D0
    else
      df = acos(cdf)
      if (sdf .LT. 0.D0) df = 2*pi-df
    endif
    phi = phie-df
    DO WHILE (phi .LT. 0.0D0) ; phi=phi+2*Pi ; ENDDO
    cphi=cos(phi) ; sphi=sin(phi) ; kx=st*cphi ; ky=st*sphi ; kz=ct
END SUBROUTINE SScattering

SUBROUTINE TIMES(time)
character(4) S4
!character*10 ST$
character*10 ST
integer h,m,s,ms 
real  time
!CALL DATE_AND_TIME (time = ST$)
CALL DATE_AND_TIME (time = ST)
!
S4 = ST(1:2)  ; read(S4,'(i2)') h
S4 = ST(3:4)  ; read(S4,'(i2)') m
S4 = ST(5:6)  ; read(S4,'(i2)') s
S4 = ST(8:10) ; read(S4,'(i3)') ms
time = 3600*DBLE(h)+ 60*DBLE(m) + DBLE(s) + DBLE(ms)/1000
END SUBROUTINE TIMES
END MODULE INTEGRATE
