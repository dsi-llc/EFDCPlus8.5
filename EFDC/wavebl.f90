SUBROUTINE WAVEBL 
  ! ** WAVEBL subroutine is used to specify information for wave boundary layer only        
  ! ** NTSWV  = Number of time steps for gradual introduction of wave forcing
  ! ** WVLCAL = 1: To calculate wave length / 0: wave length is from wave.inp
  ! ** ISDZBR = Write diagnostics for effective wave current boundary layer roughness
  ! ** IFWAVE = 0 To use wave.inp
  ! **          1 To directly input SWAN formatted files 
  ! ** SWANGRP= 1 To use SWAN output for whole grid / 0 SWAN output for locations (x,y)
      
  ! CHANGE RECORD 
  ! DATE MODIFIED     BY               DESCRIPTION        
  !-- ------------------------------------------------------------------
  ! 2017-05-23        Dang Chung       IMPROVED WAVE CALCULATION
  ! 2012-08-04        Dang Chung       IMPROVED WAVE INPUT & CALCULATION 
  ! 2011-08-08        PAUL M. CRAIG &  RESTRUCTURED TO F90 AND CORRECTED CODE
  !                   Dang Chung      
  

  USE GLOBAL
  USE GETSWANMOD     ! ** DIRECTLY READ SWAN OUTPUT
  USE WAVELENGTH
  
  IMPLICIT NONE
  INTEGER :: NWV,IWV,JWV,NW,L,K
  INTEGER :: IOS,ND,LF,LL
  
  REAL(RKD) :: WDEP1,RRLS,WPRD

  REAL   :: WWVH,WANGLE,WWPRDP,WVLEN,DISPTMP
  REAL   :: AEXTMP,UWORBIT,REYWAVE
  REAL   :: RA,FCW,CDTMP,VISMUDD
  REAL, EXTERNAL :: CSEDVIS
  
  JSWRPH=0
  IF( JSWAVE == 0 )THEN
    ! *** INITIALIZE AND INPUT WAVE INFORMATION       
    JSWRPH=1
    DO L=1,LC
      HMPW(L)=0.
      HMUW(L)=0.
      HMVW(L)=0.
      WVKHC(L)=0.
      WVKHU(L)=0.
      WVKHV(L)=0.
      WVTMP1(L)=0.
      WVTMP2(L)=0.
      WVTMP3(L)=0.
      WVTMP4(L)=0.
      WVENEP(L)=0.
      UWVSQ(L)=0.
      QQWC(L)=1.E-12
      QQWCR(L)=1.E-12
      QQWV1(L)=1.E-12   ! *** BOUNDARY LAYER STRESS DUE TO WAVES
      QQWV2(L)=1.E-12   ! *** WATER COLUMN STRESS DUE TO WAVES
      QQWV3(L)=1.E-12
    ENDDO
    DO K=1,KC
      DO L=1,LC
        WVHUU(L,K)=0.
        WVHVV(L,K)=0.
        WVHUV(L,K)=0.
        WVPP(L,K)=0.
        WVPU(L,K)=0.
        WVPV(L,K)=0.
        FXWAVE(L,K)=0.
        FYWAVE(L,K)=0.
      ENDDO
    ENDDO

    IF( ISSTEAD == 0 )THEN
      ! ** UNSTEADY WAVE
      WRITE(*,'(A)')'WAVE: READING WAVETIME.INP'     
      CALL GETWAVEDAY
      IF( WAVEDAY(1) > TIMEDAY .OR. WAVEDAY(NWVTIM) < TIMEDAY) STOP 'TIMEDAY IS OUT OF THE WAVE-TIME RANGE!'
      DO NW=1,NWVTIM-1  
        IF( WAVEDAY(NW+1) > TIMEDAY .AND. WAVEDAY(NW) <= TIMEDAY )THEN
          IWVCOUNT = NW
          EXIT
        ENDIF
      ENDDO
    ELSE
      ! ** STEADY WAVE
      IWVCOUNT = 1
      NWVTIM   = 1
      ALLOCATE(WAVEDAY(2))
      WAVEDAY(1:2) = TIMEDAY
    ENDIF
    
    ! *** READ FIRST WAVE PERIOD FIELD
    IF( IFWAVE == 0 )THEN
      WRITE(*,'(A)')'READING WAVE.INP'
      OPEN(WUNI,FILE='wave.inp')
      CALL SKIPCOM(WUNI,'*')
      DO NW=1,IWVCOUNT-1
        DO NWV=2,LA
          READ(WUNI,*,IOSTAT=IOS) IWV,JWV,WWVH,WANGLE,WWPRDP,WVLEN,DISPTMP
          IF( IOS > 0 )THEN
            WRITE(*,'(A45,I5,A12,I5)') '***  READ ERROR ON FILE WAVE.INP AT CELL L = ',NWV,'BLOCK NO = ',NW
            STOP
          ENDIF
        ENDDO
      ENDDO
      
      ! ** NEXT WAVE RECORD
      WRITE(*,'(A36,F12.4)')'WAVE: READING WAVE.INP AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
      CALL GETWAVEINP
      IF( IWVCOUNT == NWVTIM )THEN
        CLOSE(WUNI)
        WRITE(*,'(A)') '** WAVE: END OF WAVE RECORD'
      ENDIF      
      
    ELSEIF (IFWAVE == 1 )THEN
      WRITE(*,'(A38,F12.4)')'WAVE: READING SWAN OUPUT AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
      IF( SWANGRP == 1 )THEN
        CALL GETSWAN_GRP
      ELSE
        CALL GETSWAN_LOC
      ENDIF
    ENDIF    
    JSWAVE=1
    
  ELSEIF( JSWAVE == 1 .AND. IWVCOUNT < NWVTIM .AND. TIMEDAY >= WAVEDAY(IWVCOUNT+1) )THEN
    ! *** UPDATE WAVE PARAMETERS
    IWVCOUNT = IWVCOUNT + 1   
    IF( IFWAVE == 0 )THEN
      WRITE(*,'(A36,F12.4)')'WAVE: READING WAVE.INP AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
      CALL GETWAVEINP
      IF( IWVCOUNT == NWVTIM )THEN
        CLOSE(WUNI)
        WRITE(*,'(A)') '** WAVE: END OF WAVE RECORD'
      ENDIF      
      
    ELSEIF( IFWAVE == 1 )THEN
      WRITE(*,'(A38,F12.4)')'WAVE: READING SWAN OUPUT AT WAVE-TIME:',REAL(WAVEDAY(IWVCOUNT))
      IF( SWANGRP == 1 )THEN
        CALL GETSWAN_GRP
      ELSE
        CALL GETSWAN_LOC
      ENDIF
    ENDIF
  ENDIF 
    
  ! ****************************************************************************
  ! *** FINISHED READING DATA, NOW SETUP THE WAVE TABLE AND COMPUTE WAVE PARAMETERS
  ! ***  GENERATE WAVE TABLE
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(ND,L,LF,LL,WDEP1,WPRD,RRLS)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
        
    DO L=LF,LL
      WV(L).HEIGHT = MIN(0.75*HP(L),WV(L).HEISIG)                     ! *** INCLUDING BREAKING WAVE
      IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRYWAV )THEN
        WDEP1 = HP(L)
        WPRD  = 2.*PI/WV(L).FREQ
        IF( WVLCAL == 1 )THEN
          CALL BISEC(DISRELATION,WLMIN,WLMAX,EPS,WPRD,WDEP1,0._8,0._8,RRLS)
          WV(L).LENGTH = RRLS
        ENDIF
        WV(L).K   = MAX( 2.*PI/WV(L).LENGTH, 0.01 )                   ! *** ANGULAR WAVE NUMBER (RAD/M)
        WV(L).KHP = MIN(WV(L).K*HP(L),SHLIM)
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO
  
  ! **  HDRYWAV IS USED TO LIMIT WAVE HEIGHT INSTEAD OF HP(L) = 0.55
  ! **  THE WAVE TURBULENT INTENSITY, QQWV
  ! **  AND SQUARED HORIZONTAL WAVE OBRITAL VELOCITY MAGNITUDE

  ! *** WRITE DIAGNOSTICS
  !$OMP SINGLE
  IF( ISDZBR > 0 )THEN
    OPEN(1,FILE=OUTDIR//'WAVEBL_DIA.OUT')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'WAVEBL_DIA.OUT')
  ENDIF
  !$OMP END SINGLE

  !$OMP DO PRIVATE(ND,L,LF,LL,AEXTMP,UWORBIT,VISMUDD,REYWAVE,RA,FCW,CDTMP)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
    DO L= LF,LL  !2,LA
      ! *** SET ZBRE AS N IKURADSE ROUGHNESS
      IF( ISTRAN(7) > 0 )THEN
        ! *** BASE ON NIKURADSE ROUGHNESS (APPROXIMATE 2.5*D50)
        ZBRE(L)=MAX(SEDDIA50(L,KBT(L)),1E-6)*2.5
      ELSE
        ZBRE(L)=KSW    !Z0 = KSW/30
      ENDIF

      IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRYWAV )THEN
        AEXTMP   = 0.5*WV(L).HEIGHT/SINH(WV(L).KHP)
        UWORBIT  = AEXTMP*WV(L).FREQ
        UWVSQ(L) = UWORBIT*UWORBIT
        IF( UWORBIT < 1.E-6 )THEN
          UWVSQ(L)=0.    
          QQWV1(L)=0.
          CYCLE
        ENDIF
        VISMUDD=1.36E-6      !DHC: 2010-05-06
        IF( ISMUD >= 1 ) VISMUDD=CSEDVIS(SED(L,KSZ(L),1))
        REYWAVE=UWORBIT*AEXTMP/VISMUDD
        RA= AEXTMP/ZBRE(L)       !Relative roughness = A/Ksw

        ! *** COMPUTE FRICTION FACTOR DUE TO WAVE
        IF( REYWAVE <= 5D5 )THEN
          !LAMINAR
          FCW  = 2*REYWAVE**(-0.5)

        ELSEIF( REYWAVE>5D5 .AND. RA>1.57 )THEN
          !** TURBULENT SMOOTH WAVE BOUNDARY LAYER
          FCW = 0.09*REYWAVE**(-0.2)

        ELSEIF( REYWAVE>5D5 .AND. RA <= 1.57 )THEN
          !** TURBULENT ROUGH WAVE BOUNDARY LAYER
          FCW = EXP(5.2*RA**(-0.19)-6)    ! *** Baird's paper
          FCW = MIN(FCW,0.3)
        ENDIF
        CDTMP=0.5*FCW
        QQWV1(L)=MIN(CDTMP*UWORBIT*UWORBIT,QQMAX)

      ELSE
        QQWV1(L)=0.
        UWVSQ(L)=0.
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL
  
  IF( ISDZBR > 0 )THEN
    WRITE(1,'(A)') '**  L    I    J    WVWHA(L)   WVFRQL(L)    WVLEN(L)    WVKHP(L)    UWVSQ(L)    QQWV1(L)     ZBRE(L)'
    DO L=2,LA
        WRITE(1,'(3I5,7f12.5)') L,IL(L),JL(L),WV(L).HEIGHT ,WV(L).FREQ ,WV(L).LENGTH ,WV(L).K ,UWVSQ(L) ,QQWV1(L) ,ZBRE(L)
    ENDDO
    CLOSE(1)
  ENDIF
  
END SUBROUTINE
