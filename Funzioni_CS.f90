MODULE FUNZIONI

USE CONSTANTS

IMPLICIT NONE

contains

SUBROUTINE STAR_EMIS(prob,temp,ct,Phi,energy,Flagr)
REAL*8 cp,sp,ct,st,Phi, kxe,kye,kze, energy, prob,temp
INTEGER Flagr
  CALL RAN_DIR(prob,temp, kxe,kye,kze,energy,Flagr)
  sp = sin(Phi)  ; cp = cos(Phi) 
  if(ct.GT.0.9999999999999) ct = 0.9999999999*ct/abs(ct)
  st = sqrt(1.D0-ct**2)
  kx = kze*st*cp + kxe*ct*cp - kye*sp
  ky = kze*st*sp + kxe*ct*sp + kye*cp
  kz = kze*ct - kxe*st
END	SUBROUTINE STAR_EMIS

SUBROUTINE RAN_DIR(prob,temp, kxe,kye,kze,energy,Flagr)
REAL*8 kxe,kye,kze,st,ct,Phi,energy,prob,temp
INTEGER Flagr
!
   CALL BARE_REJC(prob,temp,energy,ct,Phi,Flagr)
   st = sqrt(1.D0- ct**2)
   kxe = st*cos(phi) ; kye = st*sin(phi) ; kze = ct
END	SUBROUTINE RAN_DIR

SUBROUTINE BARE_REJC(prob,temp,energy,ct,phi,Flagr)
REAL*8  emis, ekT, y, theta, ct, phi, emistot, emis1, emiss1, temp, prob, ener, energy
REAL*8  seedpol, seed0, seedth, seedph
REAL*8  B13,dens,epe,eci,ece,ec,ecfx,n0
INTEGER  Flagr

COMMON/bare/B13,dens,epe,eci,ece,ec,ecfx,n0

emis = 0.D0
y = 1.D0
CALL RANDOM_NUMBER(seedpol)
IF (seedpol .LT. prob) THEN
   Flagr = -1
ELSE
   Flagr = 1
ENDIF
DO WHILE (emis .LT. y)
   CALL RANDOM_NUMBER(seed0)
   CALL RANDOM_NUMBER(seedth)
   CALL RANDOM_NUMBER(seedph)
   ekT = BBODY()
   theta = (pi/2.D0)*seedth
   phi = (2.D0*pi)*seedph
   y = seed0
   ener = ekT*temp
   emistot = baremiss(ener,theta,phi)
   emis1 = baremiss1(ener,theta,phi)
   IF (Flagr .LT. 0) THEN
      emis = emis1
   ELSE
      emis = 2.D0*emistot - emis1
   ENDIF
ENDDO

energy = ener
ct = cos(theta)

END SUBROUTINE BARE_REJC

SUBROUTINE PROBTOT(prob,temp)
INTEGER, PARAMETER ::  nenerf=100, nmuf=100, nphif=100
REAL*8  mukmin, mukmax, xmin, xmax, temp
REAL*8  muk(nmuf), phk(nphif), xv(nenerf), xvl(nenerf)
REAL*8  preint(nenerf,nmuf,nphif), preint1(nenerf,nmuf,nphif)
REAL*8  preintphi(nenerf,nmuf,nphif-1), preintphi1(nenerf,nmuf,nphif-1)
REAL*8  workphi(nphif-1), workphi1(nphif-1), workphis, workphi1s
REAL*8  INTPHIF(nenerf,nmuf), INTPHIF1(nenerf,nmuf) 
REAL*8  preintmu(nenerf,nmuf-1), preintmu1(nenerf,nmuf-1)
REAL*8  workmu(nmuf-1), workmu1(nmuf-1), workmus, workmu1s
REAL*8  INTMUF(nenerf), INTMUF1(nenerf)
REAL*8  preinten(nenerf-1), preinten1(nenerf-1)
REAL*8  preintang(nmuf,nphif), preintang1(nmuf,nphif), intang(nenerf), intang1(nenerf)
REAL*8  preintango1(nmuf,nphif)
REAL*8  ENNE, ENNEUNO, prob, thetaB
REAL*8  B13,dens,epe,eci,ece,ec,ecfx,n0
INTEGER  i, j, k

COMMON/bare/B13,dens,epe,eci,ece,ec,ecfx,n0
COMMON thetaB

mukmin=0.001D0 ; mukmax=0.999D0
DO j = 1, nmuf
   muk(j) = mukmin + ((mukmax-mukmin)*(j-1))/float(nmuf-1)
ENDDO
DO k = 1, nphif
   phk(k) = 2*pi*(k-1)/float(nphif-1)
ENDDO
xmin=0.01D0 ; xmax=10.D0
DO i = 1, nenerf
   xv(i) = xmin + ((xmax-xmin)*(i-1))/float(nenerf-1)
ENDDO
xvl = log(xv)

INTPHIF = 0.D0 ; INTPHIF1 = 0.D0
INTMUF = 0.D0 ; INTMUF1 = 0.D0

DO i = 1, nenerf
   DO j = 1, nmuf
      DO k = 1, nphif
         preint(i,j,k) = (baremiss(xv(i)*temp,acos(muk(j)),phk(k)))*bb(xv(i))
         preint1(i,j,k) = (baremiss1(xv(i)*temp,acos(muk(j)),phk(k)))*bb(xv(i))
      ENDDO
   ENDDO
ENDDO
!
DO i = 1, nenerf
   DO j = 1, nmuf
      DO k = 2, nphif
         preintphi(i,j,k-1) = (0.5D0*(preint(i,j,k)+preint(i,j,k-1)))*(phk(k)-phk(k-1))
         preintphi1(i,j,k-1) = (0.5D0*(preint1(i,j,k)+preint1(i,j,k-1)))*(phk(k)-phk(k-1))
      ENDDO
   ENDDO
ENDDO
DO i = 1, nenerf
   DO j = 1, nmuf
      workphi = preintphi(i,j,1:nphif-1)
      workphi1 = preintphi1(i,j,1:nphif-1)
      workphis = SUM(workphi)
      workphi1s = SUM(workphi1)
      INTPHIF(i,j) = workphis
      INTPHIF1(i,j) = workphi1s
   ENDDO
ENDDO
!
DO i = 1, nenerf
   DO j = 2, nmuf
      preintmu(i,j-1) = (0.5D0*(INTPHIF(i,j)+INTPHIF(i,j-1)))*(muk(j)-muk(j-1))
      preintmu1(i,j-1) = (0.5D0*(INTPHIF1(i,j)+INTPHIF1(i,j-1)))*(muk(j)-muk(j-1))
   ENDDO
ENDDO
DO i = 1, nenerf
   workmu = preintmu(i,1:nmuf-1)
   workmu1 = preintmu1(i,1:nmuf-1)
   workmus = SUM(workmu)
   workmu1s = SUM(workmu1)
   INTMUF(i) = workmus
   INTMUF1(i) = workmu1s
ENDDO
!
DO i = 2, nenerf
   preinten(i-1) = (0.5D0*(INTMUF(i)+INTMUF(i-1)))*(xv(i)*temp-xv(i-1)*temp)
   preinten1(i-1) = (0.5D0*(INTMUF1(i)+INTMUF1(i-1)))*(xv(i)*temp-xv(i-1)*temp)
ENDDO
ENNE = SUM(preinten)
ENNEUNO = SUM(preinten1)
!
prob = ENNEUNO/(2.D0*ENNE)

ENDSUBROUTINE PROBTOT

SUBROUTINE FLUSSOTOT(flusso,temp,thetaB,B13)
INTEGER, PARAMETER ::  nenerf=100, nmuf=100, nphif=100
REAL*8  mukmin, mukmax, xmin, xmax, temp
REAL*8  muk(nmuf), phk(nphif), xv(nenerf), xvl(nenerf)
REAL*8  preint1(nenerf,nphif,nmuf), preint2(nenerf,nphif,nmuf)
REAL*8  preintmu1(nenerf,nphif,nmuf-1), preintmu2(nenerf,nphif,nmuf-1)
REAL*8  workmu1(nmuf-1), workmu2(nmuf-1), workmu1s, workmu2s
REAL*8  INTMUF1(nenerf,nphif), INTMUF2(nenerf,nphif) 
REAL*8  preintphi1(nenerf,nphif-1), preintphi2(nenerf,nphif-1)
REAL*8  workphi1(nphif-1), workphi2(nphif-1), workphi1s, workphi2s
REAL*8  INTPHIF1(nenerf), INTPHIF2(nenerf)
REAL*8  preinten1(nenerf-1), preinten2(nenerf-1)
REAL*8  ENNEUNO, ENNEDUE, flusso
REAL*8  thetaB, B13
INTEGER  i, j, k

ENNEUNO = 0.D0
ENNEDUE = 0.D0
mukmin=0.001D0 ; mukmax=0.999D0
DO j = 1, nmuf
   muk(j) = mukmin + ((mukmax-mukmin)*(j-1))/float(nmuf-1)
ENDDO
DO k = 1, nphif
   phk(k) = 2*pi*(k-1)/float(nphif-1)
ENDDO
xmin=0.01D0 ; xmax=20.D0
DO i = 1, nenerf
   xv(i) = xmin + ((xmax-xmin)*(i-1))/float(nenerf-1)
ENDDO
xvl = log(xv)

INTPHIF1 = 0.D0 ; INTPHIF2 = 0.D0
INTMUF2 = 0.D0 ; INTMUF2 = 0.D0

DO i = 1, nenerf
   DO j = 1, nphif
      DO k = 1, nmuf
         preint1(i,j,k) = ((baremiss1n(xv(i)*temp,acos(muk(k)),phk(j),thetaB,B13))*bb(xv(i)))*muk(k)
         preint2(i,j,k) = ((baremiss2n(xv(i)*temp,acos(muk(k)),phk(j),thetaB,B13))*bb(xv(i)))*muk(k)
      ENDDO
   ENDDO
ENDDO
!
DO i = 1, nenerf
   DO j = 1, nphif
      DO k = 2, nmuf
         preintmu1(i,j,k-1) = (0.5D0*(preint1(i,j,k)+preint1(i,j,k-1)))*(muk(k)-muk(k-1))
         preintmu2(i,j,k-1) = (0.5D0*(preint2(i,j,k)+preint2(i,j,k-1)))*(muk(k)-muk(k-1))
      ENDDO
   ENDDO
ENDDO
DO i = 1, nenerf
   DO j = 1, nphif
      workmu1 = preintmu1(i,j,1:nmuf-1)
      workmu2 = preintmu2(i,j,1:nmuf-1)
      workmu1s = SUM(workmu1)
      workmu2s = SUM(workmu2)
      INTMUF1(i,j) = workmu1s
      INTMUF2(i,j) = workmu2s
   ENDDO
ENDDO
!
DO i = 1, nenerf
   DO j = 2, nphif
      preintphi1(i,j-1) = (0.5D0*(INTMUF1(i,j)+INTMUF1(i,j-1)))*(phk(j)-phk(j-1))
      preintphi2(i,j-1) = (0.5D0*(INTMUF2(i,j)+INTMUF2(i,j-1)))*(phk(j)-phk(j-1))
   ENDDO
ENDDO
DO i = 1, nenerf
   workphi1 = preintphi1(i,1:nphif-1)
   workphi2 = preintphi2(i,1:nphif-1)
   workphi1s = SUM(workphi1)
   workphi2s = SUM(workphi2)
   INTPHIF1(i) = workphi1s
   INTPHIF2(i) = workphi2s
ENDDO
!
DO i = 2, nenerf
   preinten1(i-1) = (0.5D0*(INTPHIF1(i)+INTPHIF1(i-1)))*(xv(i)*temp-xv(i-1)*temp)
   preinten2(i-1) = (0.5D0*(INTPHIF2(i)+INTPHIF2(i-1)))*(xv(i)*temp-xv(i-1)*temp)
ENDDO
ENNEUNO = SUM(preinten1)
ENNEDUE = SUM(preinten2)

flusso = ENNEUNO+ENNEDUE

ENDSUBROUTINE FLUSSOTOT

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

real*8 function baremiss(ene,thetak,phik)
real*8 thetaB,thetB,dens
real*8 al,ali,alr,ctb,stb,stk,ctk,epet,ectfx,eci,epe,ene,thetak,phik,ece,ec,ecfx,B13,n0
real*8 cali,calr,sgncali,sgncalr,alphai,alphar,alpha
real*8 a1,aa1
real*8 n0p,n0m,r0p,r0m,j0,j00,a,ja,ja1
real*8 ject,jeci,pfx,jbfx,ntfx,jcfx,el,wlfx,xfx,ellfx,jecti,ject1,p1,rellfx,fLfx,jb1fx,jeci1
real*8 ypie1i,xpie1i,ypie2i,xpie2i,ypre1r,xpre1r,ypre2r,xpre2r
real*8 emisfx,emiss1fx,emis1fx,emis2fx,emisXfx,emisOfx
COMMON/bare/B13,dens,epe,eci,ece,ec,ecfx,n0
common thetaB
!
thetB=thetaB
if (cos(thetaB).LT.0.D0) then
   thetB=4.D0*atan(1.D0)-thetaB
   phik=4.D0*atan(1.D0)-phik
endif
stb=sin(thetB)
ctb=cos(thetB)
stk=sin(thetak)
ctk=cos(thetak)
if (ctk.gt.0.9999D0) then
   ctk = 0.9999D0
   stk = sqrt(1.-ctk**2)
endif
cali = stb*stk*cos(phik)-ctb*ctk
calr = stb*stk*cos(phik)+ctb*ctk
sgncali = cali/abs(cali)
sgncalr = calr/abs(calr)
if (abs(cali).gt.0.9999D0) then 
   cali = 0.9999D0*sgncali
endif
if (abs(calr).gt.0.9999D0) then 
   calr = 0.9999D0*sgncalr
endif
ali=acos(cali)
alr=acos(calr)
al=min(ali,alr)
alphai=ali
alphar=alr
alpha=al
!
epet=epe*sqrt(3.d0-2.d0*ctk)
ectfx=epet**2/ece
!
ypie1i = (cos(thetB)*sin(thetak)+sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
xpie1i = (sin(thetB)*sin(phik))/sin(alphai)
ypie2i = (sin(thetB)*sin(phik))/sin(alphai)
xpie2i = (-cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
!
ypre1r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
xpre1r = (sin(thetB)*sin(phik))/sin(alphar)
ypre2r = (-sin(thetB)*sin(phik))/sin(alphar)
xpre2r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
!
ject=0.5d0+0.05d0*(1.d0+ctb*stk)/(1.d0+B13)-0.15d0*(1.d0-ctb)*sin(al)
pfx=0.1d0*(1.d0+stb)/(1.d0+B13)
jbfx=ject/(1.d0-pfx+pfx*(ectfx/ene)**0.6)
ject1= 0.5d0 + 0.05d0/(1.d0 + B13) + stb/4.d0
jb1fx=ject1/(0.1d0+0.9d0*(ectfx/ene)**0.4)
!
if (ene.le.ectfx) then
        emisfx=jbfx
        emiss1fx=jb1fx
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then 
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss=(emisXfx+emisOfx)/2.D0
        return
endif

if (ene.gt.ectfx) then
        el=epe*(1.d0+1.2d0*(1.d0-ctk)**1.5)*(1.d0-stb**2/3.d0)
        ntfx=sqrt(1.d0-epet**2/(ece*(ene)))
        jcfx=4.d0*ntfx/(1.d0+ntfx)**2
        wlfx=0.8d0*(ectfx/epe)**0.2*sqrt(sin(al/2.d0))*(1.d0+stb**2)
        xfx=(ene-el)/(1.d0-ctk)/2.d0/epe/wlfx
        ellfx=stk**2*wlfx*(0.17d0*epe/ecfx/(1.d0+xfx**4)+0.21d0*exp(-(ene/epe)**2))
        rellfx=stb**0.25*(2.d0-(sin(al))**4)*ellfx/(1.d0+ellfx)
        fLfx = 1./(1.+exp(5*((EL-ene)/(EL-ectfx))))
        emisfx=jbfx*(1.d0-jcfx)+jcfx/(1.d0+ellfx)
        emiss1fx=jb1fx*(1.d0-jcfx)+jcfx*(1.d0-rellfx)
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss=(emisXfx+emisOfx)/2.D0
        return
endif
end function baremiss

real*8 function baremiss1(ene,thetak,phik)
real*8 thetaB,thetB,dens
real*8 al,ali,alr,ctb,stb,stk,ctk,epet,ectfx,eci,epe,ene,thetak,phik,ece,ec,ecfx,B13,n0
real*8 cali,calr,sgncali,sgncalr,alphai,alphar,alpha
real*8 a1,aa1
real*8 n0p,n0m,r0p,r0m,j0,j00,a,ja,ja1
real*8 ject,jeci,pfx,jbfx,ntfx,jcfx,el,wlfx,xfx,ellfx,jecti,ject1,p1,rellfx,fLfx,jb1fx,jeci1
real*8 ypie1i,xpie1i,ypie2i,xpie2i,ypre1r,xpre1r,ypre2r,xpre2r
real*8 emisfx,emiss1fx,emis1fx,emis2fx,emisXfx,emisOfx
COMMON/bare/B13,dens,epe,eci,ece,ec,ecfx,n0
common thetaB
!
thetB=thetaB
if (cos(thetaB).LT.0.D0) then
   thetB=4.D0*atan(1.D0)-thetaB
   phik=4.D0*atan(1.D0)-phik
endif
stb=sin(thetB)
ctb=cos(thetB)
stk=sin(thetak)
ctk=cos(thetak)
if (ctk.gt.0.9999D0) then
   ctk = 0.9999D0
   stk = sqrt(1.-ctk**2)
endif
cali = stb*stk*cos(phik)-ctb*ctk
calr = stb*stk*cos(phik)+ctb*ctk
sgncali = cali/abs(cali)
sgncalr = calr/abs(calr)
if (abs(cali).gt.0.9999D0) then 
   cali = 0.9999D0*sgncali
endif
if (abs(calr).gt.0.9999D0) then 
   calr = 0.9999D0*sgncalr
endif
ali=acos(cali)
alr=acos(calr)
al=min(ali,alr)
alphai=ali
alphar=alr
alpha=al
!
epet=epe*sqrt(3.d0-2.d0*ctk)
ectfx=epet**2/ece
!
ypie1i = (cos(thetB)*sin(thetak)+sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
xpie1i = (sin(thetB)*sin(phik))/sin(alphai)
ypie2i = (sin(thetB)*sin(phik))/sin(alphai)
xpie2i = (-cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
!
ypre1r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
xpre1r = (sin(thetB)*sin(phik))/sin(alphar)
ypre2r = (-sin(thetB)*sin(phik))/sin(alphar)
xpre2r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
!
ject=0.5d0+0.05d0*(1.d0+ctb*stk)/(1.d0+B13)-0.15d0*(1.d0-ctb)*sin(al)
pfx=0.1d0*(1.d0+stb)/(1.d0+B13)
jbfx=ject/(1.d0-pfx+pfx*(ectfx/ene)**0.6)
ject1= 0.5d0 + 0.05d0/(1.d0 + B13) + stb/4.d0
jb1fx=ject1/(0.1d0+0.9d0*(ectfx/ene)**0.4)
!
if (ene.le.ectfx) then
        emisfx=jbfx
        emiss1fx=jb1fx
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then 
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss1=emisXfx
        return
endif

if (ene.gt.ectfx) then
        el=epe*(1.d0+1.2d0*(1.d0-ctk)**1.5)*(1.d0-stb**2/3.d0)
        ntfx=sqrt(1.d0-epet**2/(ece*(ene)))
        jcfx=4.d0*ntfx/(1.d0+ntfx)**2
        wlfx=0.8d0*(ectfx/epe)**0.2*sqrt(sin(al/2.d0))*(1.d0+stb**2)
        xfx=(ene-el)/(1.d0-ctk)/2.d0/epe/wlfx
        ellfx=stk**2*wlfx*(0.17d0*epe/ecfx/(1.d0+xfx**4)+0.21d0*exp(-(ene/epe)**2))
        rellfx=stb**0.25*(2.d0-(sin(al))**4)*ellfx/(1.d0+ellfx)
        fLfx = 1./(1.+exp(5*((EL-ene)/(EL-ectfx))))
        emisfx=jbfx*(1.d0-jcfx)+jcfx/(1.d0+ellfx)
        emiss1fx=jb1fx*(1.d0-jcfx)+jcfx*(1.d0-rellfx)
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss1=emisXfx
        return
endif
end function baremiss1

real*8 function baremiss2(ene,thetak,phik)
real*8 thetaB,thetB,dens
real*8 al,ali,alr,ctb,stb,stk,ctk,epet,ectfx,eci,epe,ene,thetak,phik,ece,ec,ecfx,B13,n0
real*8 cali,calr,sgncali,sgncalr,alphai,alphar,alpha
real*8 a1,aa1
real*8 n0p,n0m,r0p,r0m,j0,j00,a,ja,ja1
real*8 ject,jeci,pfx,jbfx,ntfx,jcfx,el,wlfx,xfx,ellfx,jecti,ject1,p1,rellfx,fLfx,jb1fx,jeci1
real*8 ypie1i,xpie1i,ypie2i,xpie2i,ypre1r,xpre1r,ypre2r,xpre2r
real*8 emisfx,emiss1fx,emis1fx,emis2fx,emisXfx,emisOfx
COMMON/bare/B13,dens,epe,eci,ece,ec,ecfx,n0
common thetaB
!
thetB=thetaB
if (cos(thetaB).LT.0.D0) then
   thetB=4.D0*atan(1.D0)-thetaB
   phik=4.D0*atan(1.D0)-phik
endif
stb=sin(thetB)
ctb=cos(thetB)
stk=sin(thetak)
ctk=cos(thetak)
if (ctk.gt.0.9999D0) then
   ctk = 0.9999D0
   stk = sqrt(1.-ctk**2)
endif
cali = stb*stk*cos(phik)-ctb*ctk
calr = stb*stk*cos(phik)+ctb*ctk
sgncali = cali/abs(cali)
sgncalr = calr/abs(calr)
if (abs(cali).gt.0.9999D0) then 
   cali = 0.9999D0*sgncali
endif
if (abs(calr).gt.0.9999D0) then 
   calr = 0.9999D0*sgncalr
endif
ali=acos(cali)
alr=acos(calr)
al=min(ali,alr)
alphai=ali
alphar=alr
alpha=al
!
epet=epe*sqrt(3.d0-2.d0*ctk)
ectfx=epet**2/ece
!
ypie1i = (cos(thetB)*sin(thetak)+sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
xpie1i = (sin(thetB)*sin(phik))/sin(alphai)
ypie2i = (sin(thetB)*sin(phik))/sin(alphai)
xpie2i = (-cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
!
ypre1r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
xpre1r = (sin(thetB)*sin(phik))/sin(alphar)
ypre2r = (-sin(thetB)*sin(phik))/sin(alphar)
xpre2r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
!
ject=0.5d0+0.05d0*(1.d0+ctb*stk)/(1.d0+B13)-0.15d0*(1.d0-ctb)*sin(al)
pfx=0.1d0*(1.d0+stb)/(1.d0+B13)
jbfx=ject/(1.d0-pfx+pfx*(ectfx/ene)**0.6)
ject1= 0.5d0 + 0.05d0/(1.d0 + B13) + stb/4.d0
jb1fx=ject1/(0.1d0+0.9d0*(ectfx/ene)**0.4)
!
if (ene.le.ectfx) then
        emisfx=jbfx
        emiss1fx=jb1fx
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then 
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss2=emisOfx
        return
endif

if (ene.gt.ectfx) then
        el=epe*(1.d0+1.2d0*(1.d0-ctk)**1.5)*(1.d0-stb**2/3.d0)
        ntfx=sqrt(1.d0-epet**2/(ece*(ene)))
        jcfx=4.d0*ntfx/(1.d0+ntfx)**2
        wlfx=0.8d0*(ectfx/epe)**0.2*sqrt(sin(al/2.d0))*(1.d0+stb**2)
        xfx=(ene-el)/(1.d0-ctk)/2.d0/epe/wlfx
        ellfx=stk**2*wlfx*(0.17d0*epe/ecfx/(1.d0+xfx**4)+0.21d0*exp(-(ene/epe)**2))
        rellfx=stb**0.25*(2.d0-(sin(al))**4)*ellfx/(1.d0+ellfx)
        fLfx = 1./(1.+exp(5*((EL-ene)/(EL-ectfx))))
        emisfx=jbfx*(1.d0-jcfx)+jcfx/(1.d0+ellfx)
        emiss1fx=jb1fx*(1.d0-jcfx)+jcfx*(1.d0-rellfx)
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss2=emisOfx
        return
endif
end function baremiss2

real*8 function baremiss1n(ene,thetak,phik,thetaB,B13)
real*8 thetaB,thetB,dens
real*8 al,ali,alr,ctb,stb,stk,ctk,epet,ectfx,eci,epe,ene,thetak,phik,ece,ec,ecfx,B13,n0
real*8 cali,calr,sgncali,sgncalr,alphai,alphar,alpha
real*8 a1,aa1
real*8 n0p,n0m,r0p,r0m,j0,j00,a,ja,ja1
real*8 ject,jeci,pfx,jbfx,ntfx,jcfx,el,wlfx,xfx,ellfx,jecti,ject1,p1,rellfx,fLfx,jb1fx,jeci1
real*8 ypie1i,xpie1i,ypie2i,xpie2i,ypre1r,xpre1r,ypre2r,xpre2r
real*8 emisfx,emiss1fx,emis1fx,emis2fx,emisXfx,emisOfx
REAL*8, parameter ::  Z=26.D0, capA=56.D0, eta=1.D0

dens=8.9D3*eta*capA*Z**(-0.6)*B13**(1.2)
epe=0.0288D0*sqrt(dens*Z/capA)
eci=0.0635D0*(Z/capA)*B13
ece=115.77D0*B13
ec=eci+(epe**2)/ece
ecfx=(epe**2)/ece
n0=sqrt(1.D0+(epe**2)/(2.D0*ece*eci))
!
thetB=thetaB
if (cos(thetaB).LT.0.D0) then
   thetB=4.D0*atan(1.D0)-thetaB
   phik=4.D0*atan(1.D0)-phik
endif
stb=sin(thetB)
ctb=cos(thetB)
stk=sin(thetak)
ctk=cos(thetak)
if (ctk.gt.0.9999D0) then
   ctk = 0.9999D0
   stk = sqrt(1.-ctk**2)
endif
cali = stb*stk*cos(phik)-ctb*ctk
calr = stb*stk*cos(phik)+ctb*ctk
sgncali = cali/abs(cali)
sgncalr = calr/abs(calr)
if (abs(cali).gt.0.9999D0) then 
   cali = 0.9999D0*sgncali
endif
if (abs(calr).gt.0.9999D0) then 
   calr = 0.9999D0*sgncalr
endif
ali=acos(cali)
alr=acos(calr)
al=min(ali,alr)
alphai=ali
alphar=alr
alpha=al
!
epet=epe*sqrt(3.d0-2.d0*ctk)
ectfx=epet**2/ece
!
ypie1i = (cos(thetB)*sin(thetak)+sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
xpie1i = (sin(thetB)*sin(phik))/sin(alphai)
ypie2i = (sin(thetB)*sin(phik))/sin(alphai)
xpie2i = (-cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
!
ypre1r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
xpre1r = (sin(thetB)*sin(phik))/sin(alphar)
ypre2r = (-sin(thetB)*sin(phik))/sin(alphar)
xpre2r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
!
ject=0.5d0+0.05d0*(1.d0+ctb*stk)/(1.d0+B13)-0.15d0*(1.d0-ctb)*sin(al)
pfx=0.1d0*(1.d0+stb)/(1.d0+B13)
jbfx=ject/(1.d0-pfx+pfx*(ectfx/ene)**0.6)
ject1= 0.5d0 + 0.05d0/(1.d0 + B13) + stb/4.d0
jb1fx=ject1/(0.1d0+0.9d0*(ectfx/ene)**0.4)
!
if (ene.le.ectfx) then
        emisfx=jbfx
        emiss1fx=jb1fx
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then 
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss1n=emisXfx
        return
endif

if (ene.gt.ectfx) then
        el=epe*(1.d0+1.2d0*(1.d0-ctk)**1.5)*(1.d0-stb**2/3.d0)
        ntfx=sqrt(1.d0-epet**2/(ece*(ene)))
        jcfx=4.d0*ntfx/(1.d0+ntfx)**2
        wlfx=0.8d0*(ectfx/epe)**0.2*sqrt(sin(al/2.d0))*(1.d0+stb**2)
        xfx=(ene-el)/(1.d0-ctk)/2.d0/epe/wlfx
        ellfx=stk**2*wlfx*(0.17d0*epe/ecfx/(1.d0+xfx**4)+0.21d0*exp(-(ene/epe)**2))
        rellfx=stb**0.25*(2.d0-(sin(al))**4)*ellfx/(1.d0+ellfx)
        fLfx = 1./(1.+exp(5*((EL-ene)/(EL-ectfx))))
        emisfx=jbfx*(1.d0-jcfx)+jcfx/(1.d0+ellfx)
        emiss1fx=jb1fx*(1.d0-jcfx)+jcfx*(1.d0-rellfx)
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss1n=emisXfx
        return
endif
end function baremiss1n

real*8 function baremiss2n(ene,thetak,phik,thetaB,B13)
real*8 thetaB,thetB,dens
real*8 al,ali,alr,ctb,stb,stk,ctk,epet,ectfx,eci,epe,ene,thetak,phik,ece,ec,ecfx,B13,n0
real*8 cali,calr,sgncali,sgncalr,alphai,alphar,alpha
real*8 a1,aa1
real*8 n0p,n0m,r0p,r0m,j0,j00,a,ja,ja1
real*8 ject,jeci,pfx,jbfx,ntfx,jcfx,el,wlfx,xfx,ellfx,jecti,ject1,p1,rellfx,fLfx,jb1fx,jeci1
real*8 ypie1i,xpie1i,ypie2i,xpie2i,ypre1r,xpre1r,ypre2r,xpre2r
real*8 emisfx,emiss1fx,emis1fx,emis2fx,emisXfx,emisOfx
REAL*8, parameter ::  Z=26.D0, capA=56.D0, eta=1.D0
   
dens=8.9D3*eta*capA*Z**(-0.6)*B13**(1.2)
epe=0.0288D0*sqrt(dens*Z/capA)
eci=0.0635D0*(Z/capA)*B13
ece=115.77D0*B13
ec=eci+(epe**2)/ece
ecfx=(epe**2)/ece
n0=sqrt(1.D0+(epe**2)/(2.D0*ece*eci))
!
thetB=thetaB
if (cos(thetaB).LT.0.D0) then
   thetB=4.D0*atan(1.D0)-thetaB
   phik=4.D0*atan(1.D0)-phik
endif
stb=sin(thetB)
ctb=cos(thetB)
stk=sin(thetak)
ctk=cos(thetak)
if (ctk.gt.0.9999D0) then
   ctk = 0.9999D0
   stk = sqrt(1.-ctk**2)
endif
cali = stb*stk*cos(phik)-ctb*ctk
calr = stb*stk*cos(phik)+ctb*ctk
sgncali = cali/abs(cali)
sgncalr = calr/abs(calr)
if (abs(cali).gt.0.9999D0) then 
   cali = 0.9999D0*sgncali
endif
if (abs(calr).gt.0.9999D0) then 
   calr = 0.9999D0*sgncalr
endif
ali=acos(cali)
alr=acos(calr)
al=min(ali,alr)
alphai=ali
alphar=alr
alpha=al
!
epet=epe*sqrt(3.d0-2.d0*ctk)
ectfx=epet**2/ece
!
ypie1i = (cos(thetB)*sin(thetak)+sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
xpie1i = (sin(thetB)*sin(phik))/sin(alphai)
ypie2i = (sin(thetB)*sin(phik))/sin(alphai)
xpie2i = (-cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphai)
!
ypre1r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
xpre1r = (sin(thetB)*sin(phik))/sin(alphar)
ypre2r = (-sin(thetB)*sin(phik))/sin(alphar)
xpre2r = (cos(thetB)*sin(thetak)-sin(thetB)*cos(thetak)*cos(phik))/sin(alphar)
!
ject=0.5d0+0.05d0*(1.d0+ctb*stk)/(1.d0+B13)-0.15d0*(1.d0-ctb)*sin(al)
pfx=0.1d0*(1.d0+stb)/(1.d0+B13)
jbfx=ject/(1.d0-pfx+pfx*(ectfx/ene)**0.6)
ject1= 0.5d0 + 0.05d0/(1.d0 + B13) + stb/4.d0
jb1fx=ject1/(0.1d0+0.9d0*(ectfx/ene)**0.4)
!
if (ene.le.ectfx) then
        emisfx=jbfx
        emiss1fx=jb1fx
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then 
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss2n=emisOfx
        return
endif

if (ene.gt.ectfx) then
        el=epe*(1.d0+1.2d0*(1.d0-ctk)**1.5)*(1.d0-stb**2/3.d0)
        ntfx=sqrt(1.d0-epet**2/(ece*(ene)))
        jcfx=4.d0*ntfx/(1.d0+ntfx)**2
        wlfx=0.8d0*(ectfx/epe)**0.2*sqrt(sin(al/2.d0))*(1.d0+stb**2)
        xfx=(ene-el)/(1.d0-ctk)/2.d0/epe/wlfx
        ellfx=stk**2*wlfx*(0.17d0*epe/ecfx/(1.d0+xfx**4)+0.21d0*exp(-(ene/epe)**2))
        rellfx=stb**0.25*(2.d0-(sin(al))**4)*ellfx/(1.d0+ellfx)
        fLfx = 1./(1.+exp(5*((EL-ene)/(EL-ectfx))))
        emisfx=jbfx*(1.d0-jcfx)+jcfx/(1.d0+ellfx)
        emiss1fx=jb1fx*(1.d0-jcfx)+jcfx*(1.d0-rellfx)
        if (emiss1fx.lt.2.D0*emisfx-1.D0) then
           emis1fx = 2.D0*emisfx-1.D0
        else
           if (emiss1fx.gt.2.D0*emisfx) then
              emis1fx = 2.D0*emisfx
           else
              if (emiss1fx.gt.1.D0) then
                 emis1fx = 1.D0
              else
                 emis1fx = emiss1fx
              endif
           endif
        endif
        emis2fx = 2.D0*emisfx-emis1fx
        
        emisXfx = xpre1r**2*emis2fx+xpre2r**2*emis1fx
        emisOfx = ypre1r**2*emis2fx+ypre2r**2*emis1fx
        
        if (emisXfx.lt.0.D0) then 
           emisXfx=0.D0
        endif
        if (emisOfx.lt.0.D0) then 
           emisOfx=0.D0
        endif
        if (ISNAN(emisXfx)) then 
           emisXfx=0.D0
        endif
        if (ISNAN(emisOfx)) then 
           emisOfx=0.D0
        endif
        
        baremiss2n=emisOfx
        return
endif
end function baremiss2n

real*8 function bb(x)
real*8 x
!
if (x.lt.0.01d0) then 
   bb=x 
else 
   bb=x**2/(exp(x)-1.d0)
endif
return
end function bb

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
