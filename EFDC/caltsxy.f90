SUBROUTINE CALTSXY  

  ! ***SUBROUTINE CALTSXY UPDATES TIME VARIABLE SURFACE WIND STRESS  
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------!
  !    2015-01       Paul M. Craig     Added fully coupled Ice Sub-model
  !                  Dang H Chung        
  !    2014-12       Paul M. Craig     Added the new Evaporation Approach
  !    2014-10       D H CHUNG         ADD ICE THICKNESS CALCULATION USING CE-QUAL-W2
  !    2011-03       PAUL M. CRAIG     MERGED LATEST CODES AND RESTRUCTURED, ADDED OMP

  USE GLOBAL
  USE HEAT_MODULE, ONLY:EQUILIBRIUM_TEMPERATURE
  
  IMPLICIT NONE
  
  INTEGER :: ND,LF,LL,LP,L,LS,LN,I,J,NA,NW,M1,M2,NI,ICECOVL,IWMP,ITMP,IMAP,IICE
  REAL    :: TIME,TDIFF,WTM1,WTM2,DEGM1,DEGM2,WINDS1,WINDS2,WINDE1,WINDE2,WINDN1,WINDN2,TAUICE,CONVRT2
  REAL    :: WNDFAC,C2,TSEAST,TSNORT,WINDXX,WINDYY,TMPVAL,SVPWET,CLEVAPTMP,CCNHTTTMP,TMPVL1,U10,CD10
  REAL    :: VAPORP, WSPD, ET, CSHE
  CHARACTER(1) :: VE1,VC1
  CHARACTER(6) :: WS4(20),AS4(20)
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: CLOUDTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: EVAPTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PATMTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: RAINTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: RHATT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SOLSWRTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SVPATT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TATMTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TDEWTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TWETTT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: VPATT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WINDE  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WINDN  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WINDSXX
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WINDSXY
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WINDSYX
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WINDSYY
  
  IF(  .NOT. ALLOCATED(CLOUDTT) )THEN
    ALLOCATE(CLOUDTT(NASERM))  
    ALLOCATE(EVAPTT(NASERM))  
    ALLOCATE(PATMTT(NASERM))  
    ALLOCATE(RAINTT(NASERM))  
    ALLOCATE(RHATT(NASERM))  
    ALLOCATE(SOLSWRTT(NASERM))  
    ALLOCATE(SVPATT(NASERM))  
    ALLOCATE(TATMTT(NASERM))  
    ALLOCATE(TDEWTT(NASERM))  
    ALLOCATE(TWETTT(NASERM))  
    ALLOCATE(VPATT(NASERM))  
    ALLOCATE(WINDE(NWSERM))  
    ALLOCATE(WINDN(NWSERM))  
    ALLOCATE(WINDSXX(LCM))  
    ALLOCATE(WINDSXY(LCM))  
    ALLOCATE(WINDSYX(LCM))  
    ALLOCATE(WINDSYY(LCM))  

    CLOUDTT=0.0   
    EVAPTT=0.0   
    PATMTT=0.0   
    RAINTT=0.0   
    RHATT=0.0 
    SOLSWRTT=0.0   
    SVPATT=0.0  
    TATMTT=0.0 
    TDEWTT=0.0
    TWETTT=0.0   
    VPATT=0.0  
    WINDE=0.0   
    WINDN=0.0 
    WINDSXX=0.0
    WINDSXY=0.0
    WINDSYX=0.0
    WINDSYY=0.0

    IF( ISVHEAT > 0 )THEN
      DO L=2,LA  
        CLEVAP(L) = SVREVC(L)
        CCNHTT(L) = SVRCHC(L) 
      ENDDO 
    ELSE
      DO L=2,LA  
        CLEVAP(L) = 0.001*ABS(REVC)
        CCNHTT(L) = 0.001*ABS(RCHC)  
      ENDDO 
    ENDIF
    
    IF ( NASER == 0 ) TATMT(2:LA) = TEM(2:LA,KC)  ! *** ASSIGNED INITIAL TEMP.
    
    ! *** WRITE TO LOG  FILE
    IF( ISTRAN(2) > 0 .AND. ISVHEAT > 0 )THEN
      OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
      WRITE(8,'(//A)') 'SURFACE HEAT PARAMETERS'
      FORALL(I=1:NWSER) WS4(I) = 'WSER'//CHAR(48+I)
      FORALL(I=1:NASER) AS4(I) = 'ASER'//CHAR(48+I)
      WRITE(8,'(A5,3A10,"  | ",20A8)') 'L','REVC','RCHC','BACKKE', WS4(1:NWSER), AS4(1:NASER)
      DO L = 2,LA
        IF( LSVHTWINDE(L) )THEN
          VE1='+'
        ELSE
          VE1=' '
        ENDIF
        IF( LSVHTWINDC(L) )THEN
          VC1='+'
        ELSE
          VC1=' '
        ENDIF
        IF( NWSER > 1 .AND. NASER > 1 )THEN
          WRITE(8,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), WNDWHT(1:NWSER,L,1), ATMWHT(1:NASER,L,1)
        ELSEIF( NWSER > 1 )THEN
          WRITE(8,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), WNDWHT(1:NWSER,L,1), 1.
        ELSEIF( NASER > 1 )THEN
          WRITE(8,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), 1., ATMWHT(1:NASER,L,1)
        ELSE
          WRITE(8,'(L5,2(F9.3,A1),F10.3,"  | ",20F8.4)') L, 1000.*SVREVC(L),VE1, 1000.*SVRCHC(L),VC1, SVKEBACK(L), 1., 1.
        ENDIF
      ENDDO
    ENDIF
    
  ENDIF

  !**********************************************************************!
  ! INITIALIZE WIND SHELTERED SURFACE GAS TRANSFER
  IF( N == -1 .AND. NWSER > 0 )THEN
    IF( DEBUG )THEN
      OPEN(1,FILE=OUTDIR//'WINDSHELT.OUT')  
      CLOSE(1,STATUS='DELETE')  
      OPEN(1,FILE=OUTDIR//'WINDSHELT.OUT')  
    ENDIF
    DO L=2,LA  
      I=IL(L)  
      J=JL(L)  
      IF( WINDSTKA(L) > 0.0 )THEN  
        ! ** IF WINDSTKA > 0 BOTH X AND Y COMPONENTS ARE APPLIED  
        WINDSXX(L)=CVN(L)  
        WINDSXY(L)=-CVE(L)  
        WINDSYX(L)=-CUN(L)  
        WINDSYY(L)=CUE(L)  
      ELSE  
        ! ** IF WINDSTKA < 0 SLECTIVELY APPLY X AND Y COMPONENTS  
        ! ** FIRST CASE IS FULLY OPEN WATER  
        WINDSXX(L)=CVN(L)  
        WINDSXY(L)=-CVE(L)  
        WINDSYX(L)=-CUN(L)  
        WINDSYY(L)=CUE(L)  
        LS=LSC(L)  
        LN=LNC(L)  
        ! ** SECOND CASE IS 1D CHANNEL IN COMP X DIRECTION  
        IF( SVB(L) < 0.5 .AND. IJCT(I,J-1)/=5 )THEN  
          IF( SVB(LN) < 0.5 .AND. IJCT(I,J+1)/=5 )THEN  
            WINDSXX(L)=CVN(L)  
            WINDSXY(L)=-CVE(L)  
            WINDSYX(L)=-1000.  
            WINDSYY(L)=0.  
          ENDIF  
        ENDIF  
        ! ** THIRD CASE IS 1D CHANNEL IN COMP Y DIRECTION  
        IF( SUB(L) < 0.5 .AND. IJCT(I-1,J)/=5 )THEN  
          IF( SUB(LEC(L)) < 0.5 .AND. IJCT(I+1,J)/=5 )THEN  
            WINDSXX(L)=0.  
            WINDSXY(L)=-1000.  
            WINDSYX(L)=-CUN(L)  
            WINDSYY(L)=CUE(L)  
          ENDIF  
        ENDIF  
      ENDIF
      IF( DEBUG)WRITE(1,1111)IL(L),JL(L),WINDSTKA(L),WINDSXX(L),WINDSXY(L),WINDSYX(L),WINDSYY(L)  
    ENDDO  
    IF( DEBUG)CLOSE(1)  
    1111 FORMAT(2I5,10F10.6)  
  ENDIF  

  LDAYLIGHT = .FALSE.
  LCHECKICE = .FALSE.

  !**********************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)

  !**********************************************************************
  ! *** INTERPOLATE THE WINDS FOR EACH SERIES TO THE CURRENT TIME
  IF( NWSER > 0 )THEN  
    ! *** UPDATE THE FORCING WIND DATA TO THE CURRENT TIME
    !$OMP SINGLE
    DO NW=1,NWSER  
      TIME=TIMESEC/TCWSER(NW)  
      M2=MWTLAST(NW)
      DO WHILE (TIME > TSWND(NW).TIM(M2))
        M2=M2+1
        IF( M2 > TSWND(NW).NREC )THEN
          M2 = TSWND(NW).NREC   !** THIS ALLOWS USING EXTRAPOLATION !!!
          EXIT
        ENDIF
      ENDDO
      MWTLAST(NW)=M2
      M1 = M2-1
            
      TDIFF=TSWND(NW).TIM(M2)-TSWND(NW).TIM(M1)  
      WTM1=(TSWND(NW).TIM(M2)-TIME)/TDIFF  
      WTM2=(TIME-TSWND(NW).TIM(M1))/TDIFF  
      DEGM1=90.-TSWND(NW).VAL(M1,2)
      DEGM2=90.-TSWND(NW).VAL(M2,2)  
      WINDS1=WTM1*TSWND(NW).VAL(M1,1)+WTM2*TSWND(NW).VAL(M2,1)  
      WINDS2=WTM1*TSWND(NW).VAL(M1,1)+WTM2*TSWND(NW).VAL(M2,1)  
      WINDE1=TSWND(NW).VAL(M1,1)*COS(DEGM1/57.29578)  
      WINDN1=TSWND(NW).VAL(M1,1)*SIN(DEGM1/57.29578)  
      WINDE2=TSWND(NW).VAL(M2,1)*COS(DEGM2/57.29578)  
      WINDN2=TSWND(NW).VAL(M2,1)*SIN(DEGM2/57.29578)  
      WINDE(NW)=WTM1*WINDE1+WTM2*WINDE2  
      WINDN(NW)=WTM1*WINDN1+WTM2*WINDN2  
      
      ! *** EE7.3 CONVERT TO 2 METERS FOR ALL CALCULATIONS (0.003 IS OPEN GRASSLAND Z0)
      CONVRT2  = LOG(2.0/0.003)/LOG(WINDH(NW)/0.003)
      WINDE(NW)=CONVRT2*WINDE(NW)
      WINDN(NW)=CONVRT2*WINDN(NW)

    ENDDO  
    
    ! --------------------------------------------------------------------------
    ! *** CALCULATE THE WIND STRESS
    ! *** HANDLE TIME VARIABLE WIND MAPS
    IF( NWSER > 1 )THEN      
      IWMP=1
      IF( NWNDMAP > 1 )THEN
        DO ITMP=1,NWNDMAP
          IF( TIMEDAY >= TWNDMAPBEG(ITMP) .AND. TIMEDAY < TWNDMAPEND(ITMP) )THEN
            IWMP=ITMP
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDIF

    ! *** CONVERSION FACTOR FROM 2 METERS TO 10 METERS (0.003 IS OPEN GRASSLAND Z0)
    CONVRT2  = LOG(10.0/0.003)/LOG(2./0.003)
    !$OMP END SINGLE
    
    !$OMP DO PRIVATE(ND,LF,LL,L,WNDFAC,C2,TSEAST,TSNORT,WSPD,WINDXX,WINDYY,U10,CD10)
    DO ND=1,NDM  
      LF=2+(ND-1)*LDM  
      LL=MIN(LF+LDM-1,LA)

      ! *** GET 10 METER WIND VECTORS FOR EACH CELL
      IF( NWSER > 1 )THEN
        ! *** MULTIPLE STATIONS
        DO L=LF,LL 
          WNDVELE(L) = SUM(WNDWHT(1:NWSER,L,IWMP)*WINDE(1:NWSER))*CONVRT2
          WNDVELN(L) = SUM(WNDWHT(1:NWSER,L,IWMP)*WINDN(1:NWSER))*CONVRT2
        ENDDO              
      ELSEIF( NWSER == 1 )THEN
        ! *** SINGLE STATION
        DO L=LF,LL
          WNDVELE(L)=WINDE(1)*CONVRT2
          WNDVELN(L)=WINDN(1)*CONVRT2
        ENDDO 
      ENDIF

      DO L=LF,LL
        ! ** CASE 0 MAGNITUDE SHELTERING AND NO DIRECTIONAL SHELTERING  
        IF( WINDSTKA(L) > 0.0 )THEN  
          WNDFAC = WINDSTKA(L)
          WNDVELE(L) = WNDFAC*WNDVELE(L)  
          WNDVELN(L) = WNDFAC*WNDVELN(L)  
          WSPD = SQRT( WNDVELE(L)*WNDVELE(L) + WNDVELN(L)*WNDVELN(L) )  
          U10 = MAX(WSPD,0.25)
          IF( IWDRAG < 2 )THEN
            ! *** ORIGINAL EFDC WIND DRAG
            IF( U10 < 5. )THEN
              C2 = 1./U10
              CD10 = 3.83111E-005*(C2**3) - 0.000308715*(C2**2) + 0.00116012*C2 + 0.000899602
            ELSEIF( U10 >= 5.0 .AND. U10 <= 7. )THEN
              CD10 = -5.37642E-006*(U10**3) + 0.000112556*(U10**2) - 0.000721203*U10 + 0.00259657
            ELSEIF( U10 >= 7. )THEN
              CD10 = -3.99677E-007*(U10**2) + 7.32937E-005*U10 + 0.000726716
            ENDIF
          ELSEIF( IWDRAG == 2 )THEN
            ! *** USE HERSBACH 2011, EUROPEAN CENTRE FOR MEDIUM-RANGE WEATHER FORECASTS (ECMWF) 
            CD10 = (0.00103 + 0.00004 * U10**1.48) / U10**0.21            
          ELSEIF( IWDRAG == 3 )THEN
            ! *** USE COARE3.6 (EDSON, ET.AL. 2013) SIMPLIFICATION AT NUETRAL BOUYANCY
            IF( U10 > 20. )THEN
              CD10 = 8.394E-05*U10 + 8.285E-04
            ELSE
              CD10 = 1.376E-08*U10**4 - 9.841E-07*U10**3 + 2.324E-05*U10**2 - 1.012E-04*U10 + 8.971E-04
            ENDIF
          ENDIF
          IF( IWDRAG > 0 )THEN
            ! *** COMPUTE WIND SPEED RELATIVE TO WATER
            TSEAST = WINDSXX(L)*WNDVELE(L) + WINDSXY(L)*WNDVELN(L) - U(L,KC)
            TSNORT = WINDSYX(L)*WNDVELE(L) + WINDSYY(L)*WNDVELN(L) - V(L,KC)
            ! *** 1.225E-3 = RHOa/RHOw (DIMENSIONLESS)
            TSX(L) = 1.225E-3*CD10*U10*TSEAST             ! *** TSX IS THE WIND SHEAR IN THE U DIRECTION (M2/S2)
            TSY(L) = 1.225E-3*CD10*U10*TSNORT             ! *** TSY IS THE WIND SHEAR IN THE V DIRECTION (M2/S2)
          ELSE
            TSEAST = 1.225E-3*CD10*U10*WNDVELE(L)
            TSNORT = 1.225E-3*CD10*U10*WNDVELN(L)
            TSX(L) = WINDSXX(L)*TSEAST + WINDSXY(L)*TSNORT
            TSY(L) = WINDSYX(L)*TSEAST + WINDSYY(L)*TSNORT
          ENDIF
          WINDCD10(L) = CD10
          
          ! *** WINDST IS WIND MAGNITUDE AT 2 METERS (HEAT EXCHANGE)
          WINDST(L) = WSPD/CONVRT2

        ELSEIF( WINDSTKA(L) < 0.0 )THEN  
          ! ** CASE 1 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, OPEN WATER  
          IF( WINDSYX(L) > -99.0 .AND. WINDSXY(L) > -99.0 )THEN  
            WNDFAC = ABS(WINDSTKA(L))  
            WNDVELE(L) = WNDFAC*WNDVELE(L)  
            WNDVELN(L) = WNDFAC*WNDVELN(L)  
            WINDST(L) = SQRT( WNDVELE(L)*WNDVELE(L)+WNDVELN(L)*WNDVELN(L) )
            C2 = 1.2E-6*(0.8+0.065*WINDST(L))
            TSEAST = C2*WINDST(L)*WNDVELE(L)  
            TSNORT = C2*WINDST(L)*WNDVELN(L)  
            TSX(L) = WINDSXX(L)*TSEAST + WINDSXY(L)*TSNORT  
            TSY(L) = WINDSYX(L)*TSEAST + WINDSYY(L)*TSNORT  
          ENDIF  
        
          ! ** CASE 2 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, X CHANNEL  
          IF( WINDSYX(L) < -99.0 )THEN  
            WINDXX=WINDSXX(L)*WNDVELE(L)+WINDSXY(L)*WNDVELN(L)  
            WNDFAC=ABS(WINDSTKA(L))  
            WINDXX=WNDFAC*WNDVELE(L)  
            WINDST(L) = ABS(WINDXX)  
            TSX(L) = 1.2E-6*(0.8+0.065*WINDST(L))*WINDST(L)*WINDXX  
            TSY(L) = 0.  
          ENDIF  

          ! ** CASE 3 MAGNITUDE SHELTERING AND DIRECTIONAL SHELTERING, Y CHANNEL  
          IF( WINDSXY(L) < -99.0 )THEN  
            WINDYY=WINDSYX(L)*WNDVELE(L)+WINDSYY(L)*WNDVELN(L)  
            WNDFAC=ABS(WINDSTKA(L))  
            WINDYY=WNDFAC*WINDYY  
            WINDST(L)=ABS(WINDYY)  
            TSX(L)=0  
            TSY(L)=1.2E-6*(0.8+0.065*WINDST(L))*WINDST(L)*WINDYY  
          ENDIF  
        ENDIF  
      ENDDO  
      
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO
    
  ENDIF    ! *** END OF NWSER>0

  ! *************************************************************************
  ! *** UPDATE ASER DATA AND COMPUTE SURFACE EXCHANGE PARAMETERS
  IF( NASER > 0 )THEN  
    ! *** UPDATE ATMOSPHERIC DATA TO THE CURRENT TIME
    !$OMP SINGLE
    DO NA=1,NASER  
      TIME=TIMESEC/TCASER(NA)  
      M2=MATLAST(NA)   
      DO WHILE (TIME > TSATM(NA).TIM(M2))
        M2=M2+1
        IF( M2 > TSATM(NA).NREC )THEN
          M2 = TSATM(NA).NREC
          EXIT
        ENDIF
      ENDDO
      MATLAST(NA) = M2 
      M1 = M2-1
                  
      TDIFF=TSATM(NA).TIM(M2)-TSATM(NA).TIM(M1)  
      WTM1=(TSATM(NA).TIM(M2)-TIME)/TDIFF  
      WTM2=(TIME-TSATM(NA).TIM(M1))/TDIFF  
      PATMTT(NA)=WTM1*TSATM(NA).VAL(M1,1)+WTM2*TSATM(NA).VAL(M2,1)  
      TATMTT(NA)=WTM1*TSATM(NA).VAL(M1,2)+WTM2*TSATM(NA).VAL(M2,2)  
      TWETTT(NA)=WTM1*TSATM(NA).VAL(M1,3)+WTM2*TSATM(NA).VAL(M2,3)  
      RAINTT(NA)=WTM1*TSATM(NA).VAL(M1,4)+WTM2*TSATM(NA).VAL(M2,4)  
      EVAPTT(NA)=WTM1*TSATM(NA).VAL(M1,5)+WTM2*TSATM(NA).VAL(M2,5)  

      SOLSWRTT(NA)=WTM1*TSATM(NA).VAL(M1,6)+WTM2*TSATM(NA).VAL(M2,6)  
      SOLSWRTT(NA)=SOLSWRTT(NA)*(1.-REFL)     ! *** ADJUST INCIDENT LIGHT FOR REFLECTANCE
      IF( SOLSWRTT(NA) > 1.0E-4 ) LDAYLIGHT = .TRUE.
    
      CLOUDTT(NA)=WTM1*TSATM(NA).VAL(M1,7)+WTM2*TSATM(NA).VAL(M2,7)  
      SVPATT(NA)=10.**((0.7859+0.03477*TATMTT(NA))/(1.+0.00412*TATMTT(NA)))  
      IF( IRELH(NA) == 0 )THEN  
        ! *** COMPUTE RHA FROM WET BULB
        TMPVAL = 0.00066*(1.0+0.00115*TWETTT(NA))
        SVPWET = 10.**((0.7859+0.03477*TWETTT(NA))/(1.+0.00412*TWETTT(NA)))  
        TMPVL1 = SVPWET-TMPVAL*PATMTT(NA)*(TATMTT(NA)-TWETTT(NA))
        RHATT(NA)=MAX(TMPVL1/ SVPATT(NA),.01)
      ELSE
        ! *** DIRECT ENTRY OF RELATIVE HUMIDITY
        RHATT(NA)=TWETTT(NA)  
      ENDIF 

      ! *** PREVENT LOG OF ZERO
      RHATT(NA) = MAX(RHATT(NA),.0001) 
      
      ! *** Compute Dew Point Temp in C.  From Jensen et al. (1990) ASCE Manual No. 70
      ! *** Ambient vapor pressure in kPa 
      VAPORP = RHATT(NA)* 0.611*EXP(17.27*TATMTT(NA)/(TATMTT(NA)+237.3)) 
      TDEWTT(NA) = (116.9+237.3*LOG(VAPORP))/(16.78-LOG(VAPORP))
      
      ! *** ACTUAL VAPOR PRESSURE (MB)
      VPATT(NA)=RHATT(NA)*SVPATT(NA)  
    ENDDO
        
    ! *** HANDLE TIME VARIABLE ATMOSPHERIC MAPS
    IF( NASER > 1 )THEN
      IMAP=1
      IF( NATMMAP > 1 )THEN
        DO ITMP=1,NATMMAP
          IF( TIMEDAY >= TATMMAPBEG(ITMP) .AND. TIMEDAY < TATMMAPEND(ITMP) )THEN
            IMAP=ITMP
            EXIT
          ENDIF
        ENDDO
      ENDIF
      
    ENDIF
    !$OMP END SINGLE
    
    ! *** UPDATE CELL ATMOSPHERIC PARAMETERS FOR ALL CELLS NOT JUST WET
    !$OMP DO PRIVATE(ND,LF,LL,L)
    DO ND=1,NDM  
      LF=2+(ND-1)*LDM  
      LL=MIN(LF+LDM-1,LA)
        
      ! *************************************************************************
      ! *** ADJUST PARAMETERS BASED ON ATM WEIGHTING, IF NEEDED
      IF( NASER > 1 )THEN  
        DO L=LF,LL
          PATMT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*PATMTT(1:NASER))
          TATMT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*TATMTT(1:NASER)) 
          RAINT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*RAINTT(1:NASER))  
          EVAPT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*EVAPTT(1:NASER))  
          SOLSWRT(L) = SUM(ATMWHT(1:NASER,L,IMAP)*SOLSWRTT(1:NASER)) 
          CLOUDT(L)  = SUM(ATMWHT(1:NASER,L,IMAP)*CLOUDTT(1:NASER))  
          SVPAT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*SVPATT(1:NASER)) 
          RHAT(L)    = SUM(ATMWHT(1:NASER,L,IMAP)*RHATT(1:NASER))  
          TDEWT(L)   = SUM(ATMWHT(1:NASER,L,IMAP)*TDEWTT(1:NASER))  
          VPAT(L)    = SUM(ATMWHT(1:NASER,L,IMAP)*VPATT(1:NASER))  
        ENDDO  
        
        IF( IATMP > 0 )THEN
          DO L=LF,LL
            ! *** CONVERT ATM PRESSURE (MB) TO METERS OF WATER * G FOR UNITS OF M2/S2 
            ATMP(L) = PATMT(L)*0.0101974*G 
          ENDDO  
        ENDIF
        
      ELSEIF( NASER == 1 )THEN
        DO L=LF,LL
          PATMT(L)   = PATMTT(1)  
          TATMT(L)   = TATMTT(1)  
          RAINT(L)   = RAINTT(1)  
          EVAPT(L)   = EVAPTT(1)  
          SOLSWRT(L) = SOLSWRTT(1)  
          CLOUDT(L)  = CLOUDTT(1)  
          SVPAT(L)   = SVPATT(1)  
          RHAT(L)    = RHATT(1)  
          TDEWT(L)   = TDEWTT(1)  
          VPAT(L)    = VPATT(1)    
        ENDDO
      ENDIF
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,CLEVAPTMP,CCNHTTTMP)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
        
      ! *** EVAPORATION WIND FUNCTION (FW)
      IF( IEVAP == 2 .OR. ISTOPT(2) == 1 .OR. ISTOPT(2) == 5 )THEN
        IF( ISVHEAT > 0 )THEN
          DO LP=LF,LL
            L=LWET(LP)  
            IF( LSVHTWINDE(L) )THEN
              CLEVAPTMP = 1.E-3*(0.8+0.065*WINDST(L))  
              CLEVAP(L) = MAX(SVREVC(L),CLEVAPTMP)  
            ENDIF
          ENDDO
        ELSE
          IF( REVC < 0. )THEN  
            CLEVAPTMP = 0.001*ABS(REVC)
            DO LP=LF,LL
              L=LWET(LP)  
              CLEVAP(L) = 1.E-3*(0.8+0.065*WINDST(L))  
              CLEVAP(L) = MAX(CLEVAP(L),CLEVAPTMP)  
            ENDDO  
          ENDIF
        ENDIF
        
      ENDIF

      ! *** CONDUCTIVE HEAT TRANSFER WIND FUNCTION
      IF(  ISTOPT(2) == 1 .OR. ISTOPT(2) == 5 )THEN
        IF( ISVHEAT > 0 )THEN
          DO LP=LF,LL
            L=LWET(LP)  
            IF( LSVHTWINDC(L) )THEN
              CCNHTTTMP = 1.E-3*(0.8+0.065*WINDST(L))  
              CCNHTT(L) = MAX(SVRCHC(L),CCNHTTTMP)  
            ENDIF
          ENDDO
        ELSE
          IF( RCHC < 0. )THEN  
            CCNHTTTMP=0.001*ABS(RCHC)
            DO LP=LF,LL
              L=LWET(LP)  
              CCNHTT(L)=1.E-3*(0.8+0.065*WINDST(L))  
              CCNHTT(L)=MAX(CCNHTT(L),CCNHTTTMP)  
            ENDDO  
          ENDIF
        ENDIF
      ENDIF
      
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO
    !$OMP SINGLE
  
    ! *** CHECK FOR POSSIBLE ICE CONDITIONS
    IF( ISICE > 2 .AND. NWSER > 0 .AND. NASER > 0 )THEN
      DO NA=1,NASER
        ! *** COMPUTE EQUILIBRIUM TEMPERATURE
        WSPD = SQRT(WINDE(1)*WINDE(1) + WINDN(1)*WINDN(1))
        CALL EQUILIBRIUM_TEMPERATURE(WSPD, TDEWTT(NA), TATMTT(NA), SOLSWRTT(NA), ET, CSHE) 
        IF( ET < 1.0 ) LCHECKICE = .TRUE.
      ENDDO
    ENDIF
    !$OMP END SINGLE

  ENDIF    ! *** END OF NASER>0 SECTION

  !**********************************************************************C
  IF( ISICE > 0 )THEN
    IF( (ISICE == 1 .OR. ISICE == 2) .AND. NISER>0 )THEN
      !$OMP SINGLE
      DO NI=1,NISER    
        TIME=TIMESEC/TCISER(NI)         ! ** THE SAME TIME UNIT IN ISER.INP
        M2=MITLAST(NI)
        DO WHILE (TIME > TSICE(NI).TIM(M2))
          M2=M2+1
          IF( M2 > TSICE(NI).NREC )THEN 
            M2 = TSICE(NI).NREC
            EXIT
          ENDIF
        ENDDO
        MITLAST(NI)=M2           
        M1 = M2-1
        
        IF( ISICE == 1 )THEN
          TDIFF=TSICE(NI).TIM(M2)-TSICE(NI).TIM(M1)
          WTM1=(TSICE(NI).TIM(M2)-TIME)/TDIFF
          WTM2=(TIME-TSICE(NI).TIM(M1))/TDIFF
          RICECOVT(NI)=WTM1*TSICE(NI).VAL(M1,2)+WTM2*TSICE(NI).VAL(M2,2)
          RICETHKT(NI)=WTM1*TSICE(NI).VAL(M1,1)+WTM2*TSICE(NI).VAL(M2,1)
          
        ELSEIF( ISICE == 2 )THEN
          RICECOVT(NI)=TSICE(NI).VAL(M1,2)
          RICETHKT(NI)=TSICE(NI).VAL(M1,1)
        ENDIF
        
      ENDDO
      !$OMP END SINGLE
                   
      IF( NISER > 1 )THEN
        ! ** ONLY FOR ISICE=1 & NISER>1
        ! *** HANDLE TIME VARIABLE ICE MAPS
        !$OMP SINGLE
        IICE=1
        IF( NICEMAP > 1 )THEN
          DO ITMP=1,NICEMAP
            IF( TIMEDAY >= TICEMAPBEG(ITMP) .AND. TIMEDAY < TICEMAPEND(ITMP) )THEN
              IICE=ITMP
              EXIT
            ENDIF
          ENDDO
        ENDIF      
        !$OMP END SINGLE

        !$OMP DO PRIVATE(ND,LF,LL,LP,L,CLEVAPTMP,CCNHTTTMP)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
        
          DO LP=LF,LL
            L=LWET(LP)  
            ICECOVER(L) = SUM(RICEWHT(IICE,L,1:NISER)*RICECOVT(1:NISER))
            ICETHICK(L) = SUM(RICEWHT(IICE,L,1:NISER)*RICETHKT(1:NISER))
          ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO
        
      ELSEIF( NISER == 1 )THEN
        !$OMP DO PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
        
          DO LP=LF,LL
            L=LWET(LP)  
            ICECOVER(L) = RICECOVT(1) 
            ICETHICK(L) = RICETHKT(1)
          ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      ENDIF

      !$OMP DO PRIVATE(ND,LF,LL,LP,L,ICECOVL)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        
        DO LP=LF,LL
          L=LWET(LP)  
          ICECOVL=NINT(ICECOVER(L))
          ICECOVER(L)=MIN(FLOAT(ICECOVL),1.0)
          IF( ICECOVL == 0 )THEN
            ICETHICK(L)=0.0
            ICECELL(L)=.FALSE.
          ELSE
            ICECELL(L)=.TRUE.
          ENDIF
        ENDDO
      ENDDO   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    
    ENDIF

    !**********************************************************************
    !$OMP DO PRIVATE(ND,LF,LL,L,TAUICE)
    DO ND=1,NDM  
      LF=2+(ND-1)*LDM  
      LL=MIN(LF+LDM-1,LA)

      DO L=LF,LL
        IF( ICECELL(L) )THEN
          ! *** SHEAR ON TOP WATER LAYER DUE TO DRAG ON BOTTOM OF ICE (CDICE=0.001)
          TAUICE = -CDICE*SQRT( U(L,KC)*U(L,KC) + V(L,KC)*V(L,KC) )
          TAUICE = TAUICE*(1. + 9.17*ICETHICK(L)/(HP(L) + ICETHICK(L)) )
          IF( TAUICE < -0.075 )THEN
            TAUICE = -0.075
            IF( HP(L) < 0.5 ) TAUICE = HP(L)*TAUICE
          ENDIF
          TSX(L) = TAUICE*U(L,KC)
          TSY(L) = TAUICE*V(L,KC)

        ENDIF
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO
    
  ENDIF     ! *** END OF ISICE>0
  
  !**********************************************************************C
  !$OMP END PARALLEL
  
  RETURN  
  
END  

