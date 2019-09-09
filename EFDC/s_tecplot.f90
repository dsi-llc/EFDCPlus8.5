SUBROUTINE S_TECPLOT
  
  USE GLOBAL
  
  IMPLICIT NONE
  
  INTEGER :: COUNTER,M
  REAL :: UPVEL,DNVEL,ATVEL,ZMIDLAYER
  INTEGER :: I,J,L,LE,K,ITEMPMSK,ILL,IDOWN
  INTEGER :: ICOUNT(TCOUNT,3)
  REAL :: UTMPS,VTMPS,UTMP,VTMP
  REAL :: MAG,MCHANGE
  REAL,DIMENSION(5) :: MB
  REAL,DIMENSION(TCOUNT,2) :: ESTORAGE
  REAL,DIMENSION(LC) :: AVGSED
  REAL,DIMENSION(KC) :: CTEMP1
  REAL,DIMENSION(IC-4) :: MAGREF1,MAGREF2
  REAL :: RHOH2O, HTCAP, HTCONT, HTEMP              !VB HEAT CONTENT VARIABLES
  REAL :: WQV2, WQV3, WQV19, WQV22, WQV16,WQV17, WQV6, TIMESTEP      !VB TEMPORARY VARIABLES FOR TIMESERIES O/P
  REAL :: WTEMP, WVOL, VOLTEMP, WTMP, WQV10, WQV14, WQV15
  REAL :: VELPDF !velocity PDF variable for ocean model
  REAL,SAVE :: ELAST,TLAST
  INTEGER,SAVE :: nstep
  LOGICAL,SAVE :: FIRSTTIME=.FALSE.
  INTEGER :: LL1,LL2
  REAL,DIMENSION(IC,JC) :: PUPSTREAM,PDNSTREAM,PUPKIN,PDNKIN,PUPPOT,PDNPOT,DELTA
  REAL,DIMENSION(IC,JC,KC) :: MDOTU,MDOTD
  REAL,DIMENSION(IC) :: PUP,PDN,KUP,KDN,PPUP,PPDN
  CHARACTER(LEN=66) :: timeline
  
  DOUBLE PRECISION,SAVE,ALLOCATABLE,DIMENSION(:)     :: CBLTOT        !(LCM)
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)     :: THCK          !(LCM)
   DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: UVEL          !(LCM,KC)
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: VVEL          !(LCM,KC)
   ALLOCATE(CBLTOT(LC))
    ALLOCATE(THCK(LC))
    ALLOCATE(UVEL(LC,KC))
    ALLOCATE(VVEL(LC,KC))
  CBLTOT = 0.0  ! tecplot only - move
  THCK=0.0        !(LCM)
    UVEL=0.0        !(LCM,KC)
  VVEL=0.0        !(LCM,KC)

  
  timeline='TEXT X=9, Y=84, T="&(ZONENAME[0001])", F=TIMES, CS=FRAME, H=3, ZN='
  IF(  .NOT. FIRSTTIME )THEN
! This opens the Tecplot output file
    IF( MAXVAL(MVEGL(2:LA))>90 )THEN !MHK devices exist
      OPEN(UNIT=222,FILE=OUTDIR//'powerout.dat')
      ITEMPMSK=1
      DO ILL=2,LA
        IF( MVEGL(ILL)>90 )THEN
          ICOUNT(ITEMPMSK,1)=ITEMPMSK
          ICOUNT(ITEMPMSK,2)=IL(ILL)
          ICOUNT(ITEMPMSK,3)=JL(ILL)
          ITEMPMSK=ITEMPMSK+1
        ENDIF
      ENDDO
      WRITE(222,'("TURBINE",100(I6,6X))')(ICOUNT(ILL,1),ILL=1,TCOUNT)
      WRITE(222,'(7X,100(3X,I3,3X,I3))') (ICOUNT(ILL,2),ICOUNT(ILL,3),ILL=1,TCOUNT)
    ENDIF
    OPEN (UNIT=111,FILE=OUTDIR//'tecplot2d.dat')
    IF( ISTRAN(8) == 1 )THEN !WQ data
      WRITE(111,'(A36)')'TITLE = "EFDC 2D Tecplot Algae Data"'
      WRITE(111,'(A77)')'VARIABLES="X","Y","U","V","SPEED","CHG(g)","P4D","NHX","NOX","DO(g)","CO2(g)"' !,"HEAT(kJ)","VOL(M3)"'
    ELSEIF( LSEDZLJ )THEN !SEDZLJ data
      WRITE(111,*)'TITLE = "EFDC 2D Tecplot Sediment Data"'
      WRITE(111,*)'VARIABLES = "I","J","X","Y","TAU","D50","CBL","SED","U","V","THICK1","SPEED"'
      OPEN (UNIT = 112,FILE = 'massbal.dat')
      WRITE(112,*)'TITLE = "EFDC Mass Balance Data"'
      WRITE(112,*)'VARIABLES = "Time","MB1","MB2","MB3","MB4","MB5","ERATE","D50"'
    ELSE !Flow data
      WRITE(111,'(A36)')'TITLE = "EFDC 2D Tecplot Flow Data"'
      WRITE(111,'(A77)')'VARIABLES="X","Y","U","V","SPEED","SHEAR","DEPTH"' !,"HEAT(kJ)","VOL(M3)"'
    ENDIF
    IF( OUTPUTFLAG == 4 )THEN   
      OPEN(444,FILE=OUTDIR//'CALIBRATION.DAT')
      WRITE(444,'("TIME        VEL-FED DEFECT  UPVEL   ATVEL   DNVEL")')
    ELSEIF( OUTPUTFLAG == 3 )THEN
      OPEN(864,FILE=OUTDIR//'TIDALREF.DAT')
      WRITE(864,'("3 ROWS OF AVERAGE VELOCITY, TOP VELOCITY, AND DEPTH")')
    ELSEIF( OUTPUTFLAG == 2 )THEN
      OPEN(468,FILE=OUTDIR//'VELPDF.DAT')
      WRITE(468,'("TIME AVERAGE VELOCITY ACROSS NARROWS I=60?")')
      OPEN(579,FILE=OUTDIR//'ZPROF.DAT')
      WRITE(579,'("TIME VELOCITIES FOR EACH LAYER FOR PROFILES")')
    ELSEIF( OUTPUTFLAG == 1 )THEN
      OPEN(765,FILE=OUTDIR//'POWERDIF.DAT')
      WRITE(765,'("TIME   POWER DIFFERENCES ACROSS DIFFERENT I-CELL COLUMNS")')
      OPEN(567,FILE=OUTDIR//'POTENTIAL.DAT')
      WRITE(567,'("TIME   POTENTIAL DIFFERENCES ACROSS DIFFERENT I-CELL COLUMNS")')
      OPEN(678,FILE=OUTDIR//'KINETIC.DAT')
      WRITE(678,'("TIME   KINETIC DIFFERENCES ACROSS DIFFERENT I-CELL COLUMNS")')
    ENDIF
    FIRSTTIME=.TRUE.
  ENDIF
  TIMESTEP=tbegin+float(N)*dt/86400.0
!  write(765,*)TIMESTEP,SUM(WQV(3,1:KC,3)*DZC(L,1:KC))
  ITEMPMSK=1
  DO ILL=2,LA
    IF( MVEGL(ILL)>90 )THEN
      ESTORAGE(ITEMPMSK,1)=SUM(ESUP(:,ILL))
      ESTORAGE(ITEMPMSK,2)=SUM(EMHK(:,ILL))
      ITEMPMSK=ITEMPMSK+1
    ENDIF
  ENDDO  
  IF( MAXVAL(MVEGL(2:LA))>90)WRITE(222,'(F7.3,100(F7.4,1X))')TIMESTEP,(ESTORAGE(ILL,1),ESTORAGE(ILL,2),ILL=1,TCOUNT)
  nstep=nstep+1
  IF( nstep>9999)PRINT*,'Tecplot timestamp is greater than 9999, increase field width or reduce writing frequency'                
  WRITE(timeline(31:34),'(I4.4)')nstep
  WRITE(111,'(A66,I4.4)')timeline,nstep
  WRITE(111,*)'ZONE T="',TIMESTEP,'" I= ' ,IC-4,' J= ' ,JC-4,' F=POINT'
  IF( ISTRAN(8) == 1 )THEN

    !VB HEAT CONTENT VARIABLE
    RHOH2O=1000     !KG/M3
    HTCAP=4.187     !kJ/KG-K
    HTCONT=0.0
    ! 2 Dimensional Output !VB REWRITTEN TO OUTPUT TIME SERIES PLOTS
    DO J=3,JC-2
      DO I=3,IC-2
        IF( LIJ(I,J) > 0 )THEN
!VB       TWATER=SUM(TEM(LIJ(I,J),1:KC)*DZC(L,1:KC))
!        TALT=MAXVAL(TWQ(LIJ(I,J))
          L=LIJ(I,J)
        UTMPS=SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
        VTMPS=SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
        MAG=SQRT((UTMPS)**2+(VTMPS)**2)
        WTEMP=SUM((TEM(L,1:KC)+273.15)*DZC(L,1:KC)*DXYP(L)*HP(L))
        WTMP=WTMP+WTEMP
        HTEMP=RHOH2O*HTCAP*WTEMP
        HTCONT=HTCONT+HTEMP            !HEAT CONTENT IN KILOJOULES                    
        WQV2=SUM(WQV(L,1:KC,2)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        WQV3=SUM(WQV(L,1:KC,3)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        WQV6=SUM(WQV(L,1:KC,6)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV10=SUM(WQV(L,1:KC,10)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV14=SUM(WQV(L,1:KC,14)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV15=SUM(WQV(L,1:KC,15)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV16=SUM(WQV(L,1:KC,16)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV17=SUM(WQV(L,1:KC,17)*DZC(L,1:KC))  !*DXYP(LIJ(I,J))*HP(L))
        WQV19=SUM(WQV(L,1:KC,19)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        WQV22=SUM(WQV(L,1:KC,22)*DZC(L,1:KC)) !*DXYP(LIJ(I,J))*HP(L))
        WRITE(111,'(5(1X,F11.2), 7(1X,F11.4))')DLON(L),DLAT(L),UTMPS,VTMPS,MAG,WQV3,WQV10,WQV14,WQV15,WQV19,WQV22
      ENDIF
      ENDDO
    ENDDO
    DO L=2, LA
      VOLTEMP = SUM(DZC(L,1:KC)*DXYP(L)*HP(L))
      WVOL = WVOL + VOLTEMP
    ENDDO
  ELSEIF( LSEDZLJ )THEN      
! OUTPUT FOR TECPLOT, SEDIMENT DATA
    MB(1) = 0.0 !Mass in bedload (kg)
    MB(2) = 0.0 !Mass in suspension (kg)
    !MB(3)     !Mass exiting in suspension (stored)
    !MB(4)     !Mass exiting as bedload (stored)
    MB(5) = 0.0 !Mass eroded from bed (kg)
    FORALL(L=2:LA) !added to convert THCK from g/cm^2 to cm
      TSEDT(L)=SUM(TSED(1:KB,L)/BULKDENS(1:KB,L))
    ENDFORALL
    FORALL(L = 2:LA) THCK(L) = TSEDT(L) - TSET0T(L)
    DO J = 3,JC-2
      DO I = 3,IC-2
        L=LIJ(I,J)
        CBLTOT(L) = 10.0*SUM(CBL(L,1:NSCM))*DXYP(L) !g/cm^3*cm*m^2*(0.001*100*100)
        CTEMP1(1:KC) = 0.001*SUM(SED(L,1:KC,1:NSCM))
        MB(1) = MB(1) + CBLTOT(L)
        DO K = 1,KC
          MB(2) = MB(2) + CTEMP1(K)*DZC(L,K)*DXYP(L)*HP(L) !g/m^3*m^2*m*(0.001)
          UVEL(L,K) = U(L,K)*CUE(L) + V(L,K)*CVE(L)
          VVEL(L,K) = U(L,K)*CUN(L) + V(L,K)*CVN(L)
          IF( I == IC-3 )THEN
            !MB(3) = MB(3) - DT*VVEL(LIJ(I,J),K)*DYP(LIJ(I,J))*DZC(L,K)*HP(LIJ(I,J))*SUM(SED(LIJ(I,J),K,1:NSCM))*0.001 !Dt*m/s*m*m*g/m^3*(0.001)
            MB(3) = MB(3) - 0.25*DT*VVEL(L,K)*DYP(L)*DZC(L,K)*(HP(L) + HP(LIJ(I+1,J)))*(SUM(SED(L,K,1:NSCM)) + SUM(SED(LIJ(I+1,J),K,1:NSCM)))*0.001 !Dt*m/s*m*m*g/m^3*(0.001)
          ENDIF
        ENDDO
        IF( I == IC-3 )THEN
          !MB(4) = MB(4) - 0.05*DT*DYP(L)*SUM((UBL(L,1:NSCM)*CUN(L) + VBL(L,1:NSCM)*CVN(L))*(CBL(L,1:NSCM) + CBL(LIJ(I+1,J),1:NSCM),1:NSCM)) !Dt*m*cm/s*g/cm^3*cm*(0.001*100)
        ENDIF
        MB(5) = MB(5) - 10.0*THCK(L)*DXYP(L) !g/cm^2*m^2*(0.001*100*100)
        AVGSED(L) = SUM(CTEMP1(1:KC)*DZC(L,1:KC))*DXYP(L)*HP(L)
        !WRITE(111,'(I4,1X,I4,1X,10(E13.4,1X))')I, J, DLON(LIJ(I,J)), DLAT(LIJ(I,J)), TAU(LIJ(I,J)), D50AVG(LIJ(I,J)), CBLTOT(LIJ(I,J)), AVGSED(LIJ(I,J)), SUM((U(LIJ(I,J),1:KC)*CUE(LIJ(I,J)) + V(LIJ(I,J),1:KC)*CVE(LIJ(I,J)))*DZC(L,1:KC)), SUM((U(LIJ(I,J),1:KC)*CUN(LIJ(I,J)) + V(LIJ(I,J),1:KC)*CVN(LIJ(I,J)))*DZC(L,1:KC)), SUM(TSED(1:KB,LIJ(I,J)))/1.4625, THCK(LIJ(I,J))/1.4625
        UTMPS=SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
        VTMPS=SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
        MAG=SQRT((UTMPS)**2+(VTMPS)**2)
        WRITE(111,'(I4,1X,I4,1X,10(E13.4,1X))')I, J, DLON(L), DLAT(L), TAU(L), D50AVG(L), CBLTOT(L), AVGSED(L), SUM((U(L,1:KC)*CUE(L) + V(L,1:KC)*CVE(L))*DZC(L,1:KC)), SUM((U(L,1:KC)*CUN(L) + V(L,1:KC)*CVN(L))*DZC(L,1:KC)), THCK(L), MAG
      ENDDO
    ENDDO
    !PRINT*,TBEGIN + FLOAT(N-1)*DT,MB(1),MB(2),MB(3),MB(4),MB(5),SUM(TAU(2:LA))/FLOAT(LA-1),MAXVAL(TAU(2:LA))
    !WRITE(112,'(8(E13.4,1X))')TBEGIN + FLOAT(N-1)*DT,MB(1),MB(2),MB(3),MB(4),MB(5),SUM(TAU(2:LA))/FLOAT(LA-1),MAXVAL(TAU(2:LA))
    WRITE(112,'(8(E13.4,1X))')TIMESTEP,MB(1),MB(2),MB(3),MB(4),MB(5),(MB(5)-(MB(1) + MB(2))-ELAST)/((TBEGIN + FLOAT(N-1)*DT)-TLAST),SUM(D50AVG(:))/(FLOAT(LA-1))
    ELAST = MB(5)-(MB(1) + MB(2)) !Erosion is saved (as mass eroded, minus the mass in the water column)
    TLAST = TBEGIN + FLOAT(N-1)*DT
    !print *,SNGL(CBL(LIJ(13,4),1:7))
    !print *,SNGL(SED(LIJ(13,4),KC,1:7))
    !print *, maxval(tau(2:LA)), minval(tau(2:LA))
  ELSE
    DO J=3,JC-2
      DO I=3,IC-2
        IF( LIJ(I,J) > 0 )THEN
          L=LIJ(I,J)
          UTMPS=SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
          VTMPS=SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
          MAG=SQRT((UTMPS)**2+(VTMPS)**2)
          WRITE(111,'(5(1X,F11.2), 7(1X,F11.4))')DLON(L),DLAT(L),UTMPS,VTMPS,MAG,TAUB(L),HP(L)
        ENDIF
      ENDDO
    ENDDO
  ENDIF	
  IF( OUTPUTFLAG == 1 )THEN !MHK look at power up- and down-stream of the MHK device
    PUPSTREAM=0.0;PDNSTREAM=0.0;PUPPOT=0.0;PDNPOT=0.0;PUPKIN=0.0;PDNKIN=0.0;PUP=0.0;PDN=0.0;KUP=0.0;KDN=0.0;PPUP=0.0;PPDN=0.0;MDOTU=0.0;MDOTD=0.0
    DO I=1,1
      DO J=3,JC-2
        LL1=LIJ(5-I,J) !looking where whatever integer is in this expression and the next
        LL2=LIJ(I+6,J)
        DO K=1,KC  
          MDOTU(I,J,K)=1024.0*U(LL1,K)*DYU(LL1)*DZC(L,K)*HU(LL1)
          MDOTD(I,J,K)=1024.0*U(LL2,K)*DYU(LL2)*DZC(L,K)*HU(LL2)
          PUPKIN(I,J)=PUPKIN(I,J)+0.5*MDOTU(I,J,K)*U(LL1,K)**2
          PDNKIN(I,J)=PDNKIN(I,J)+0.5*MDOTD(I,J,K)*U(LL2,K)**2
          PUPPOT(I,J)=PUPPOT(I,J)+MDOTU(I,J,K)*(FLOAT(K)-0.5)/DZI*G*HU(LL1)*DZC(L,K)
          PDNPOT(I,J)=PDNPOT(I,J)+MDOTD(I,J,K)*(FLOAT(K)-0.5)/DZI*G*HU(LL2)*DZC(L,K)
          continue
        ENDDO
        DELTA(I,J)=SUM(MDOTU(I,J,1:KC))-SUM(MDOTD(I,J,1:KC))
 !      EUPSTREAM(I,J)=0.5*DXP(LL1)*HP(LL1)*1024.*SUM((U(LL1,1:KC)+U(LL1+1,1:KC))*DZC(L,1:KC)*(0.5*(0.5*(U(LL1,1:KC)+U(LL1+1,1:KC)))**2+9.8106*HP(LL1)))
  !     EDNSTREAM(I,J)=0.5*DXP(LL2)*HP(LL2)*1024.*SUM((U(LL2,1:KC)+U(LL2+1,1:KC))*DZC(L,1:KC)*(0.5*(0.5*(U(LL2,1:KC)+U(LL2+1,1:KC)))**2+9.8106*HP(LL2)))
        PUPSTREAM(I,J)=PUPPOT(I,J)+PUPKIN(I,J)
        PDNSTREAM(I,J)=PDNPOT(I,J)+PDNKIN(I,J)
      ENDDO
      MCHANGE=SUM(MDOTU(I,3:JC-2,1:KC))-SUM(MDOTD(I,3:JC-2,1:KC))
      PUP(I)=SUM(PUPPOT(I,3:JC-2))
      PDN(I)=SUM(PDNPOT(I,3:JC-2))
      KUP(I)=SUM(PUPKIN(I,3:JC-2))
      KDN(I)=SUM(PDNKIN(I,3:JC-2))
      PPUP(I)=SUM(PUPSTREAM(I,3:JC-2))
      PPDN(I)=SUM(PDNSTREAM(I,3:JC-2))
      continue
    ENDDO
    WRITE(765,'(E10.3,8(1X,F10.1))')TIMESTEP,(SUM(PUPSTREAM(I,3:JC-2))-SUM(PDNSTREAM(I,3:JC-2)),I=2,8)
    WRITE(567,'(E10.3,16(1X,F10.0))')TIMESTEP,(SUM(PUPPOT(I,3:JC-2)),I=2,8),(SUM(PDNPOT(I,3:JC-2)),I=2,8)
    WRITE(678,'(E10.3,16(1X,F10.1))')TIMESTEP,(SUM(PUPKIN(I,3:JC-2)),I=2,8),(SUM(PDNKIN(I,3:JC-2)),I=2,8)
  ENDIF
 !OUTPUTFLAG=2 !Tidal reference model average and z-profile velocities
  IF( OUTPUTFLAG == 2 )THEN
    I=60
    VELPDF=0.0
    DO J=3,JC-2
      L=LIJ(I,J)
      UTMPS=SUM((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))
      VTMPS=SUM((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))
      MAG=SQRT((UTMPS)**2+(VTMPS)**2)
      VELPDF=VELPDF+MAG
    ENDDO
    VELPDF=VELPDF/FLOAT(JC-4)
    WRITE(468,*)TIMESTEP,VELPDF
    WRITE(579,'(E10.3,10(1X,F6.3))')TIMESTEP,((U(L,K)*CUE(L)+V(L,K)*CVE(L)),K=1,KC)
  ENDIF
 IF( OUTPUTFLAG == 3 )THEN !average and surface velocities for tidal reference model
   J=20
   DO I=3,IC-2
     L=LIJ(I,J)
     LE=LIJ(I+1,J)
     UTMPS=SUM(0.5*((U(L,1:KC)*CUE(L)+V(L,1:KC)*CVE(L))*DZC(L,1:KC))+0.5*((U(LE,1:KC)*CUE(LE)+V(LE,1:KC)*CVE(LE))*DZC(L,1:KC)))
     VTMPS=SUM(0.5*((U(L,1:KC)*CUN(L)+V(L,1:KC)*CVN(L))*DZC(L,1:KC))+0.5*((U(LE,1:KC)*CUN(LE)+V(LE,1:KC)*CVN(LE))*DZC(L,1:KC)))
     MAGREF1(I-2)=VTMPS
     UTMP=     0.5*(U(L,KC) *CUE(L) +V(L,KC) *CVE(L))
     UTMP=UTMP+0.5*(U(LE,KC)*CUE(LE)+V(LE,KC)*CVE(LE))
     UTMPS=UTMP
     VTMP=     0.5*(U(L,KC) *CUN(L)+ V(L,KC)* CVN(L))
     VTMP=VTMP+0.5*(U(LE,KC)*CUN(LE)+V(LE,KC)*CVN(LE))
     VTMPS=VTMP
     MAGREF2(I-2)=VTMPS
   ENDDO
   WRITE(864,'(84(F6.3,1X))')(MAGREF1(I-2),I=3,IC-2)
   WRITE(864,'(84(F6.3,1X))')(MAGREF2(I-2),I=3,IC-2)
   WRITE(864,'(84(F6.3,1X))')(HP(I-2),I=3,IC-2)
 ENDIF
 IF( OUTPUTFLAG == 4 )THEN !this interrogates the W-2-E straight flow channel for wake calibration
   UPVEL=0.0;DNVEL=0.0;ATVEL=0.0;COUNTER=0
   DO L=2,LA
     IF( MVEGL(L)>90 )THEN
       M=MVEGL(L)-90
       DO K=1,KC
         ZMIDLAYER=HP(L)*(SUM(DZC(L,1:K))-0.5*DZC(L,K))+BELV(L) !midlayer height
         IF( ZMIDLAYER >= ZMINMHK(M,L) .AND. ZMIDLAYER <= ZMAXMHK(M,L) )THEN !MHK device exists in this layer
           COUNTER=COUNTER+1 
           IDOWN=NINT(20.*WIDTHMHK(M)/DXP(L))
           ATVEL=ATVEL+U(LEC(L),K)    !1 cell downstream
           UPVEL=UPVEL+U(L-IDOWN/4,K) !20 cells upstream
           DNVEL=DNVEL+U(L+IDOWN,K)   !80 cells downstream for a 4-m turbine
         ENDIF
       ENDDO
     ENDIF
   ENDDO
   ATVEL=ATVEL/FLOAT(COUNTER)
   UPVEL=UPVEL/FLOAT(COUNTER)
   DNVEL=DNVEL/FLOAT(COUNTER)
   WRITE(444,'(e10.2,1x,5(f7.4,1x))')TIMESTEP,1.0-DNVEL/UPVEL,1.0-ATVEL/UPVEL,UPVEL,ATVEL,DNVEL !calculate velocity ratio, we are looking for 90% recovery at 20D downstream
 ENDIF
 RETURN
END SUBROUTINE S_TECPLOT
