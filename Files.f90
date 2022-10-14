MODULE FILES
USE CONSTANTS

IMPLICIT NONE
character*50 Filename
REAL*8    Dum1,Dum2,Dum3
INTEGER   Nordo, Nexto
contains

SUBROUTINE LEGGI(ip,jp,p0,Nphot)
INTEGER i,ip,jp, Nphot
REAL*8 p0
    open(unit=1,file='Mag.inp',status='old',err=100)
    read(1,*,err=100)  ip, jp, FlagP
    read(1,*,err=100)  Nphot
    read(1,*,err=100)  p0, Bp
    read(1,*,err=100)  Tstar, Tel
    read(1,*,err=100)  BetaB
    read(1,*,err=100)  beam
    read(1,*,err=100)  Emin, Emax
    read(1,*,err=100)  Printout
    close(1)
    ip = MIN(ip,NTheta) ; jp = MIN(jp,NPhi)
    DO i=1,Nin
     Xen(i) = Emin + (i-1)*(Emax-Emin)/(Nin-1)
    ENDDO
    Tel = Tel/511.D0         ! T electrons in mc^2 units
    if(ip .GT. 0) then ; itmin = ip ; itmax = ip ; endif
    if(jp .GT. 0) then ; ipmin = jp ; ipmax = jp ; endif
    if(FlagP.EQ.0) then ; FlagP=1 ; kinv=-1 ; else ; kinv=1 ; endif
    RETURN
100	write(*,*) ' Reading error in Mag.inp ...... '
END SUBROUTINE LEGGI

SUBROUTINE SALVA( ) 
    INTEGER k, Nold
    LOGICAL IFAIL
    Filename = TRIM(Directory)//'AllData.out'
    CALL READLIST(Filename,IFAIL)
    if(IFAIL) then
      open(unit=1,file=Filename,status='new')
      goto 50
    endif
    open(unit=1,file=Filename,status='old')
    read(1,*) Nold
    if(Nold.NE.Nin) goto 100
    read(1,'(5F10.4)')  (Xeno(k), k = 1, Nold)
    read(1,'(5I10)')    (Veno(k), k = 1, Nold)
    read(1,'(3F10.5)')  Dum1, Dum2, Dum3
    read(1,'(5X,2I10)')  Nordo, Nexto
    Ven = Ven + Veno
    Nord = Nord+Nordo ; Next = Next+Nexto
    REWIND(1)
50  write(1,*) Nin
    write(1,'(5F10.4)')  (Xen(k), k = 1, Nin)
    write(1,'(5I10)')  (Ven(k), k = 1, Nin)
    write(1,'(3F10.5)')  BetaB, Tel*511, Tstar
    write(1,'(5X,2I10)')  Nord, Next
100 close(1)
END SUBROUTINE SALVA

subroutine READLIST(filename,IFAIL)
character*(*) filename
   LOGICAL IFAIL
   IFAIL = .FALSE.
   open(unit=1,file=filename,status='old',ERR=10)
   RETURN
10 IFAIL = .TRUE. 
end subroutine READLIST

SUBROUTINE photon_collect(mu,phi,ener,QS,US,VS)
USE FUNZIONI
integer locmu,locphi,locen, i  ! ,reset
REAL*8  mu, phi, ener, QS,US,VS
LOGICAL FIRST
DATA FIRST/.FALSE./
if(.NOT.FIRST) then
  FIRST = .TRUE.
  do i=1, Nmuc
    mugrd(i) = -1.D0 + 2.D0*(i-1)/(Nmuc-1)
  enddo
  do i=1, Nphic
    Phigrd(i) =  2.D0*Pi*(i-1)/(NPhic-1)
  enddo
  do i=1, Nener
    egrd(i) = Emin + (i-1)*(Emax-Emin)/(Nener-1)
  enddo
endif
CALL locate(mugrd,Nmuc,mu,locmu)
CALL locate(phigrd,Nphic,phi,locphi)
CALL locate(egrd,Nener,ener,locen)
if (locen .lt. 1 .or. locen .gt. Nener) then 
!   lostphot=lostphot+1 
else 
 counts(locmu,locphi,locen) = counts(locmu,locphi,locen) + 1
 QSTOKES(locmu,locphi,locen)= QSTOKES(locmu,locphi,locen) + QS
 USTOKES(locmu,locphi,locen)= USTOKES(locmu,locphi,locen) + US
 VSTOKES(locmu,locphi,locen)= VSTOKES(locmu,locphi,locen) + VS
endif
end	SUBROUTINE photon_collect

SUBROUTINE SAVEDATA(ip,jp)
   INTEGER, ALLOCATABLE:: Co(:,:,:)
   REAL*8, ALLOCATABLE:: Qso(:,:,:),Uso(:,:,:),Vso(:,:,:) 
   character*10 SA, SB, SC
   character*12:: SFORM= '(1X,10F13.3)'
   integer i,j,k, ip,jp, Ntot, er, Nt
   REAL*8  Dum1,Dum2,p0
   COMMON/TWIST/p0
   Nt = 0
   SA = ' ' ; SB = ' ' ; SC = ' ' 
   write(SA,'(1X,I2.2)') ip ; SA = adjustL(SA) 
   write(SB,'(1X,I2.2)') jp ; SB = adjustL(SB) 
   write(SC,'(1X,I2.2)') nint(10.*p0) ; SC = adjustL(SC) 
   Filename = TRIM(Directory)//'DATA'//TRIM(SA)//'_'//TRIM(SB)//'_'//TRIM(SC)//'.out'
   CLOSE(1)
   open(unit=1,file=Filename,status='old',ERR=10)
   ALLOCATE(Co(1:Nmuc,1:Nphic,1:Nener), STAT=er)
   ALLOCATE(Qso(1:Nmuc,1:Nphic,1:Nener), STAT=er)
   ALLOCATE(Uso(1:Nmuc,1:Nphic,1:Nener), STAT=er)
   ALLOCATE(Vso(1:Nmuc,1:Nphic,1:Nener), STAT=er)

   Co = 0 ; Qso = 0.D0 ; Uso = 0.D0 ; Vso = 0.D0
   read(1,*)  i, j, k, Ntot
   read(1,*)  Dum1,Dum2
   do i=1,Nmuc-1  !  <===== ATTENZIONE
     do j=1,Nphic-1
       read(1,'(10I6)')  (Co(i,j,k), k = 1, Nener)
     enddo				
   enddo
   read(1,'(3F10.5)')  Dum1, Dum2, Dum3
   read(1,'(5X,2I10)')  Nordo, Nexto
!
   do i=1,Nmuc-1
     do j=1,Nphic-1
       read(1,SFORM)  (Qso(i,j,k), k = 1, Nener)
     enddo				
   enddo
   do i=1,Nmuc-1
     do j=1,Nphic-1
       read(1,SFORM)  (Uso(i,j,k), k = 1, Nener)
     enddo				
   enddo
   do i=1,Nmuc-1
     do j=1,Nphic-1
       read(1,SFORM)  (Vso(i,j,k), k = 1, Nener)
     enddo				
   enddo

   Nord = Nord + Nordo ; Next = Next + Nexto
   Counts = Counts + Co
   QSTOKES = QSTOKES + Qso    
   USTOKES = USTOKES + Uso    
   VSTOKES = VSTOKES + Vso    
   rewind(1)
   DEALLOCATE(Co,STAT=er)
   DEALLOCATE(Qso,STAT=er)
   DEALLOCATE(Uso,STAT=er)
   DEALLOCATE(Vso,STAT=er)
   goto 15
10 continue
   open(unit=1,file=Filename,status='new')
15 continue
   Ntot = SUM(Counts)
   write(1,*) Nmuc,Nphic,Nener,Ntot
   write(1,*) Emin,Emax
   do i=1,Nmuc-1
     do j=1,Nphic-1
       write(1,'(10I6)')  (Counts(i,j,k), k = 1, Nener)
     enddo				
   enddo
   write(1,'(3F10.5)')  BetaB, Tel*511, Tstar
   write(1,'(5X,2I10)')  Nord, Next
!
   do i=1,Nmuc-1
     do j=1,Nphic-1
       write(1,SFORM)  (QSTOKES(i,j,k), k = 1, Nener)
     enddo				
   enddo
   do i=1,Nmuc-1
     do j=1,Nphic-1
       write(1,SFORM)  (USTOKES(i,j,k), k = 1, Nener)
     enddo				
   enddo
   do i=1,Nmuc-1
     do j=1,Nphic-1
       write(1,SFORM)  (VSTOKES(i,j,k), k = 1, Nener)
     enddo				
   enddo
   close(1)
   END SUBROUTINE SAVEDATA
END MODULE FILES
