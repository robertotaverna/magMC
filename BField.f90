MODULE MAGFIELD
IMPLICIT NONE

PRIVATE
INTEGER,PARAMETER:: npt = 100
REAL*8  p, coeff,CDensity
REAL*8  Bmu(npt), Fmu(npt), Dfmu(npt)

PUBLIC:: b_fieldm,init_field 

contains

SUBROUTINE b_fieldm(r,Ct,b,bmod,Nebeta)	
USE FUNZIONI
     real*8 r,Ct,St,mubar, dfint, fint
	 REAL*8 B(3), bmod, Nebeta, RapB
	 REAL*8 smu
     mubar= Ct ; smu = SIGN(1.D0,mubar) ; mubar = abs(mubar) 
	 if(mubar.GT.0.999999D0) then
	 	b(1) = 1.D0 ; b(2) = 0.D0 ; b(3) = 0.D0
		Bmod = Bp*r**(-2.d0-p) ; Nebeta = 0.D0 ; RETURN
	 endif
	 St = sqrt(1.D0 - Ct**2)
     fint = interpol(Fmu,Bmu,npt,mubar) ; dfint = smu*interpol(Dfmu,Bmu,npt,mubar)
     B(1) = -dfint
     B(2) = p*fint/St
     RapB = sqrt(coeff/p/(p+1.d0))*fint**(1.D0/p)
     B(3) = RapB*B(2)
     B = (B*Bp/2)*r**(-2.d0-p)
	 Bmod = SQRT(SUM(B**2))
     B = B/Bmod
     Nebeta = CDensity*Bmod*RapB/r
     end SUBROUTINE b_fieldm

     SUBROUTINE init_field(p0)
USE CONSTANTS
USE FUNZIONI
	 integer k, n, err 
     integer npa, npta
	 REAL*8, ALLOCATABLE:: f(:,:), df(:,:), pv(:), cc(:)
	 REAL*8 p0 
     open(unit=1,file='Ftwist.out',status='old')
     read(1,*) npa,npta
     ALLOCATE(cc(1:npa), f(1:npa,1:npta), df(1:npa,1:npta), pv(1:npa), STAT=err)
	 read(1,*) (Bmu(k), k=1,npta)
	 DO k=1,npa
        read(1,*) pv(k)
	write(*,*) pv(k)
        DO n=1,npta
		  read(1,*) f(k,n), df(k,n), cc(k)
        ENDDO
     endDO
	 close(1)
     call locate(pv,npa,p0,n)
	 p=pv(n) ; coeff = cc(n)
	 do k=1, npta ; Fmu(k) = f(n,k) ; Dfmu(k) = df(n,k) ; enddo
     CDensity = (1.D0+p)/(4*Pi*Rstar*e_)
	 DEALLOCATE(df,STAT=err) ; DEALLOCATE(pv,STAT=err) ;
	 DEALLOCATE(cc,STAT=err) ; DEALLOCATE(f,STAT=err)    
     end SUBROUTINE init_field
END MODULE MAGFIELD

