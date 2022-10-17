MODULE CONSTANTS
 INTEGER, parameter:: NPhi=6, NTheta=7
 INTEGER, parameter:: Nmuc=10,Nphic=10,Nener=50
 INTEGER, PARAMETER:: Nin=80
 INTEGER  Flag, FlagP, kinv
 INTEGER:: Counts(Nmuc,Nphic,Nener)=0
 REAL*8:: QSTOKES(Nmuc,Nphic,Nener)=0.D0,USTOKES(Nmuc,Nphic,Nener)=0.D0,VSTOKES(Nmuc,Nphic,Nener) = 0.D0
 REAL*8,  parameter:: c_=2.997925d10,h_=6.626196d-27,k_B=1.380622d-16, Pi=3.14159265358793
 REAL*8,  parameter:: m_e=9.109558d-28,e_=4.80325d-10,r_e=2.817939d-13 
 REAL*8,  parameter:: Rstar=1.D06, omega_p=e_/(m_e *c_)
 REAL*8,  parameter:: CSigma = 2*Pi**2*r_e*c_*Rstar 
 REAL*8,  parameter:: temperaturek = 0.450D0
 REAL*8   PatPhi(NPhi), PatTh(NTheta)
 REAL*8   kx,ky,kz, omega, Bes1
 REAL*8   Beam, BetaB, GammaB, Tel, TStar, Vmean
 real*8   cth0, sth0, Phi0, X0, Y0, Z0
 REAL*8   Bp
 REAL*8   EMIN, EMAX
 REAL*8   mugrd(Nmuc), Phigrd(Nphic), egrd(Nener)

 REAL*8   Xen(Nin), Xeno(Nin)
 INTEGER  VEN(Nin), Veno(Nin)

 CHARACTER*20::  Directory='output/' 
 INTEGER   Nlost,Nord,Next   
 INTEGER:: itmin=1,itmax=NTheta,ipmin=1,ipmax=NPhi
 INTEGER   Printout

 COMMON/FIELD/Bp

END MODULE CONSTANTS

