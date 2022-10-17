!   ************* VERSIONE NON RELATIVISTICA **************
!       ===== ORDINE DI COMPILAZIONE DEI MODULI: =====
!  MODULE  CONSTANTS
!  MODULE  FUNZIONI
!  MODULE  MAGFIELD
!  MODULE  FILES
!  MODULE  PATH
!  MODULE  ODE
!  MODULE  POLAR
!  MODULE  INTEGRATE
!  Program MAGNETARS   <==== Last
! 
! PER CREARE UN UNICO FILE ASSEMBLARE I MODULI NELLO STESSO ORDINE ====
!

 INCLUDE  'Constants.f90'
 INCLUDE  'Funzioni.f90'
 INCLUDE  'BField.f90'
 INCLUDE  'Files.f90'
 INCLUDE  'Path.f90'
 INCLUDE  'Ode.f90'
 INCLUDE  'Polar.f90'
 INCLUDE  'Integrate.f90'

Program MAGNETARS
USE MAGFIELD
USE INTEGRATE
USE FILES
USE CONSTANTS
!
USE MAGFIELD
!
IMPLICIT NONE
integer ip,jp, Nphot
integer:: Npa = 0
REAL*8  p0
!
REAL*8  BS(3),BmodS,NebetaS,thetaBi, Ctin, B13i
REAL*8  temperature(Ntheta),Nphottruek(NTheta),flussoi(Ntheta),flussoscal,flussoref
INTEGER i,Nphottrue
!
COMMON/TWIST/p0
  CALL LEGGI(ip,jp,p0,Nphot)

  GammaB=1.D0/sqrt(1.D0-BetaB**2) ; Bes1=BESSK1N(1.D0/Tel)
  VMean = BetaB
  call init_field(p0)            
  call PATCHES( )
  !
  DO i=itmin,itmax
     Ctin=PatTh(i)
     CALL b_fieldm(1.D0,Ctin,BS,BmodS,NebetaS)
     thetaBi=acos(BS(1)) ; B13i=BmodS/1.D13
     temperature(i)=temperaturek
     CALL FLUSSOTOT(flussoscal,temperaturek,thetaBi,B13i)
     flussoi(i)=flussoscal
  ENDDO
  flussoref=flussoi(MINLOC(flussoi,1))
  DO i=itmin,itmax
     Nphottruek(i)=(flussoi(i)/flussoref)*Nphot
  ENDDO
  CALL RANDOM_SEED( )
  DO ip=itmin,itmax
     DO jp=ipmin,ipmax
        !
        Nphottrue=int(Nphottruek(ip))
        !
        write(*,'(A)') ' '
        CALL MainInteg(Nphottrue,ip,jp)  ; Npa = Npa+1
		  write(*,'(1X,A15,I3,A7,I3)') '- This is Patch ',Npa,' out of',(itmax-itmin+1)*(ipmax-ipmin+1)
		  write(*,*) 'Total Nphot=',Nphottrue,'temperature(keV)=',temperature(ip)
     ENDDO
  ENDDO
  PAUSE ' PREMI UN TASTO...........'
END Program MAGNETARS
