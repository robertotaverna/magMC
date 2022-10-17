MODULE PATH
IMPLICIT NONE

contains

SUBROUTINE POSITION(ell,r,cth,Phi)
USE CONSTANTS
REAL*8 ell,x,y,z,sth,cp,Phi,r, ct2, cth
   x = x0+kx*ell ; y=y0+ky*ell ; z=z0+kz*ell ; r=sqrt(x**2+y**2+z**2)
   cth = z/r  ; ct2 = cth**2 
   if(ct2.GT.0.999999) then
      Phi = 0.D0
   else
      sth = sqrt(1.D0-ct2) ;  cp = x/(r*sth) ;  Phi = acos(cp)
      if (y .LT. 0.D0)  Phi = 2*Pi - Phi
   endif
END	SUBROUTINE POSITION

SUBROUTINE CARTESIAN(vr,vtheta,vphi,Ct,Phi, vx,vy,vz)
!
REAL*8 vr,vtheta,vphi,Ct,St, Phi, vx,vy,vz, Cphi, Sphi
   St = sqrt(abs(1.D0-Ct**2)) ; Cphi = COS(Phi) ; Sphi = SIN(Phi)
   vx = vr*St*Cphi + vtheta*Ct*Cphi - vphi*Sphi
   vy = vr*St*Sphi + vtheta*Ct*Sphi + vphi*Cphi
   vz = vr*Ct - vtheta*St
END SUBROUTINE CARTESIAN

SUBROUTINE POLAR(vx,vy,vz, Theta, Phi, vr,vtheta, vphi)
REAL*8 vr,vtheta,vphi, Ct,St, Phi, vx,vy,vz, Cphi, Sphi, Theta
   Ct=COS(theta) ; St=SIN(theta) ; Cphi=COS(Phi) ; Sphi=SIN(Phi)
   vr = St*(vx*Cphi + vy*Sphi) + vz*Ct
   vtheta = Ct*(vx*Cphi + vy*Sphi) - vz*St
   vphi = vy*Cphi - vx*Sphi
END	SUBROUTINE POLAR

SUBROUTINE VANGOLI(vx,vy,vz, Ct,Phi)
REAL*8,parameter:: Pi=3.14159265358793
REAL*8 vx,vy,vz, Ct,St, Phi, sp, cp, ct2
  Ct = vz ; ct2 = ct**2 ; if(ct2 .GT. 0.99999999) ct2=0.99999999
  st = sqrt(1.D0-ct**2)
  cp = vx/st ; sp=vy/st ; Phi=acos(cp)
  if (sp .LT. 0.D0) Phi = 2*Pi - Phi
END SUBROUTINE VANGOLI

END MODULE PATH
