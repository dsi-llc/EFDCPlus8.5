MODULE WINDWAVE

! CHANGE RECORD  
! DATE MODIFIED     BY               DESCRIPTION
!----------------------------------------------------------------------!
! 2017-05-23        Dang Chung       IMPROVED WAVE CALCULATION 
! 2011-08-19        Paul M. Craig    Added OMP
! 2010-04-26        Dang Chung &     Built the WindWave Sub-Model
!                   Paul M. Craig    

USE GLOBAL
USE INFOMOD,ONLY:SKIPCOM

IMPLICIT NONE

INTEGER(IK4) ,PARAMETER :: UFET=214       ! *** FETCH.OUT
INTEGER(IK4) ,PARAMETER :: UWIN=215       ! *** LIJXY.OUT
INTEGER(IK4) ,PARAMETER :: UTAU=216       ! *** TAUW.OUT

INTEGER(IK4),PRIVATE,PARAMETER :: NZONE=16      ! *** NUMBER OF ZONES

REAL(RKD)             :: ROTAT            ! *** COUNTER-CLOCKWISE ROTATION OF DOMAIN [0,360]

REAL(RKD) ,PRIVATE,PARAMETER :: FETANG(NZONE) = (/0.0,22.5,45.0,67.5,90.0,112.5,135.0,157.5,180.0,202.5,225.0,247.5,270.0,292.5,315.0,337.5/)

PRIVATE              :: FETZONE

CONTAINS

  SUBROUTINE WINDWAVETUR 
  ! THIS OPTION OF WIND WAVE NEVER CALL THE FOLLOWING SUBROUTINES:
  ! WAVEBL & WAVESXY
    USE GLOBAL  
    INTEGER   :: L,ND,LF,LL,LP
    REAL      :: RA,CDTMP,AEXTMP,UWORBIT,VISMUDD,REYWAVE,WVFF
    REAL      :: TMPVAL,TAUTMP,CORZBR,CDRGTMP
    REAL,EXTERNAL :: CSEDVIS 

    CALL WINDWAVECAL         !OUTPUT: WV(L).HEIGHT,UDEL,WDIR,WV.FREQ,RLS
    
    !**  INITIALIZE WAVE-CURRENT BOUNDARY LAYER MODEL CALCULATING  
    !**  THE WAVE TURBULENT INTENSITY, QQWV  
    !**  AND SQUARED HORIZONTAL WAVE OBRITAL VELOCITY MAGNITUDE  
    !    WAVE REYNOLD NUMBER MUST BE USED TO DISTINGUISH HYRODYNAMIC REGIME
    
    !$OMP PARALLEL DEFAULT(SHARED) 
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)  &
    !$OMP    PRIVATE(UWORBIT,AEXTMP,VISMUDD,REYWAVE,RA,CDTMP,WVFF,TMPVAL,TAUTMP,CORZBR,CDRGTMP) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
    
      DO LP=LF,LL  
        L=LWET(LP)  
        IF( LWVMASK(L) )THEN
          ! *** SET ZBRE AS NIKURADSE ROUGHNESS
          IF( ISTRAN(7) > 0 )THEN
            ! *** BASE ON NIKURADSE ROUGHNESS (APPROXIMATE 2.5*D50)
            ZBRE(L)=MAX(SEDDIA50(L,KBT(L)),1D-6)*2.5
          ELSE
            ZBRE(L)=KSW 
          ENDIF
          IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRYWAV )THEN
            UWORBIT  = WV(L).UDEL        
            AEXTMP   = UWORBIT/WV(L).FREQ 
            UWVSQ(L) = UWORBIT*UWORBIT  
            
            IF( UWORBIT < 1.E-6 )THEN
              UWVSQ(L)=0.    
              QQWV1(L)=0.
              QQWV2(L)=0.
              CYCLE
            ENDIF
            
            IF( ISWAVE == 3 )THEN
              VISMUDD=1.36D-6                          !DHC: 2010-05-06
              IF( ISMUD >= 1 ) VISMUDD=CSEDVIS(SED(L,KSZ(L),1))  
              REYWAVE=UWORBIT*AEXTMP/VISMUDD           
              RA= AEXTMP/ZBRE(L)

              ! *** COMPUTE FRICTION FACTOR DUE TO WAVE: WVFF
              IF( REYWAVE <= 5D5 )THEN   
                !** LAMINAR
                WVFF  = 2*REYWAVE**(-0.5)
              ELSEIF( REYWAVE>5D5 .AND. RA>1.57 )THEN
                !** TURBULENT SMOOTH WAVE BOUNDARY LAYER  
                WVFF = 0.09*REYWAVE**(-0.2)
              ELSEIF( REYWAVE>5D5 .AND. RA <= 1.57 )THEN  
                !** TURBULENT ROUGH WAVE BOUNDARY LAYER  
                WVFF = EXP(5.2*RA**(-0.19)-6)    ! *** Baird's paper
                WVFF = MIN(WVFF,0.3)
              ENDIF 
              CDTMP=0.5*WVFF
              QQWV1(L)=MIN(CDTMP*UWORBIT*UWORBIT,QQMAX)
            
            ELSEIF( ISWAVE == 4 )THEN
              ! *** CHECK FOR NON-COHESIVE TRANSPORT
              IF( ISTRAN(7) > 0 )THEN
                TMPVAL=UWORBIT*SQRT( AEXTMP/(30.*ZBRE(L)) )  
                TAUTMP=TMPVAL/TAUR(NSED+1)
                CORZBR=1.+1.2*TAUTMP/(1.+0.2*TAUTMP)
                ZBRE(L)=ZBRE(L)*CORZBR
              ENDIF
              CDRGTMP=(30.*ZBRE(L)/AEXTMP)**0.2
              CDRGTMP=5.57*CDRGTMP-6.13
              CDRGTMP=EXP(CDRGTMP)
              CDRGTMP=MIN(CDRGTMP,0.22)
              TMPVAL=0.5*CDRGTMP*UWVSQ(L)
              QQWV2(L)=MIN(CTURB2*TMPVAL,QQMAX)
            ENDIF
            
          ELSE
            QQWV1(L)=QQLMIN
            QQWV2(L)=QQLMIN
            UWVSQ(L)=0.
          ENDIF
    
          WV(L).TWX=RHO*QQWV1(L)*WV(L).TWX   !True East
          WV(L).TWY=RHO*QQWV1(L)*WV(L).TWY   !True North
        ENDIF
      ENDDO  
    ENDDO    ! *** END OF DOMAIN
    !$OMP END DO
    !$OMP END PARALLEL
    
    !WAVE TURBULENT INTENSITY: QQWV2 SIMILAR TO ISWAVE=2
    !TO DECREASE TURBULENCE DECREASE CTURB (DEFAULT:16.6)
    
    IF( ISWAVE == 4) CALL W_WAVESXY

    IF( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) .AND. DEBUG )THEN
      OPEN(UTAU,FILE=OUTDIR//'TAUW.OUT',FORM='BINARY')
      WRITE(UTAU) (WV(L).HEIGHT,REAL(WV(L).PERIOD,4),REAL(WV(L).LENGTH,4),WV(L).TWX,WV(L).TWY,L=2,LA)
      CLOSE(UTAU)
    ENDIF  

  END SUBROUTINE

  SUBROUTINE WINDWAVECAL
    !CALCULATING WAVE PARAMETERS FOR EVERY CELL
    !BASED ON COMPUTED WIND PARAMETERS FROM WSER.INP AND SHELTERING
    !INPUT:
    !WNDVELE(L),WNDVELN(L),HP(L)
    !OUTPUT:
    !  WV(L).HEIGHT, WV(L).FREQ, WV(L).DIR, WV(L).UDEL
    !  WV(L).HEIGHT  - WAVE HEIGHT (M)
    !  WV(L).DIR - WAVE ANGLE (RADIANS)
    !  WV(L).FREQ - WAVE FREQENCY (SEC)
    !  WV(L).TWX, WV(L).TWY
    INTEGER :: L,ZONE,ND,LF,LL,LP
    REAL(RKD)  :: AVEDEP,WVEL2,FC1,FC2,FC3
    REAL(RKD)  :: WDIR           !WIND DIRECTION IN DEG [0,360]
    REAL(RKD)  :: WINX,WINY      !IN CURVI-LINEAR SYS
    REAL(RKD)  :: WVEL           !INTERPOLATED WIND VELOCITY

    !CALCULATING WAVE HEIGHT,PERIOD,ORBITAL VELOCITY AND LENGTH
    AVEDEP=SUM(HP(2:LA))/FLOAT(LA-1)
    
    !$OMP PARALLEL DEFAULT(SHARED) 

    !$OMP DO PRIVATE(ND,LF,LL,LP,L)  &
    !$OMP    PRIVATE(WINX,WINY,WVEL2,WVEL,WDIR,ZONE,FC1,FC2,FC3) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
    
      DO LP=LF,LL  
        L=LWET(LP)  
        IF( LWVMASK(L) )THEN
          ! *** 10 METER WINDS WNDVELE, WNDVELN ALREADY ADJUSTED FOR WIND SHELTERING
          WINX = WNDVELE(L)    ! *** X IS TRUE EAST
          WINY = WNDVELN(L)    ! *** Y IS TRUE NORTH
          
          WVEL2 = WINX*WINX + WINY*WINY
          WVEL  = SQRT(WVEL2)

          IF( HP(L) > HDRY .AND. WVEL > 1D-3 )THEN
            WV(L).TWX = WINX/WVEL
            WV(L).TWY = WINY/WVEL
            !AVEDEP=HP(L)
            IF( WINX >= 0 )THEN
              WDIR  = ACOS(WV(L).TWY)*180./PI     !DEG. (NORTH,WIND TO)
            ELSE
              WDIR  = 360-ACOS(WV(L).TWY)*180./PI
            ENDIF
            ZONE = FETZONE(WDIR)

            ! *** WAVE HEIGHT
            FC3 = TANH(0.530*(9.81*AVEDEP/WVEL2)**0.75)
            FC1 = WVEL2/9.81*0.283*FC3  
            FC2 = TANH(0.0125*(9.81*FWDIR(L,ZONE)/WVEL2)**0.42/FC3)   
            WV(L).HEIGHT = MIN(0.75*HP(L),0.7*FC1*FC2)                 ! *** WAVE HEIGHT INCLUDING BREAKING WAVE
            
            ! *** VEGETATION EFFECT
            IF( ISVEG > 0 )THEN
              IF( MVEGL(L) /= MVEGOW ) THEN
                WV(L).HEIGHT = 0.
              ELSE
                IF( (MVEGL(LWC(L)) /= MVEGOW) .OR. &
                    (MVEGL(LEC(L)) /= MVEGOW) .OR. &
                    (MVEGL(LSC(L)) /= MVEGOW) .OR. &
                    (MVEGL(LNC(L)) /= MVEGOW)      )  WV(L).HEIGHT=0.5*WV(L).HEIGHT
              ENDIF
            ENDIF
            
            ! *** WAVE FREQUENCY
            FC3 = TANH(0.833*(9.81*AVEDEP/WVEL2)**0.375)
            FC1 = (WVEL/9.81)*7.54*FC3
            FC2 = TANH(0.077*(9.81*FWDIR(L,ZONE)/WVEL2)**0.25/FC3)   
            WV(L).PERIOD = MAX(1D-3,1.1*FC1*FC2)       
            WV(L).FREQ   = 2.0*PI/WV(L).PERIOD   

            FC1 = WV(L).FREQ**2*HP(L)*GI
            FC2 = FC1 + 1.0/( 1.0 + 0.6522*FC1 + 0.4622*FC1**2 + 0.0864*FC1**4 + 0.0675*FC1**5 )
            WV(L).LENGTH = WV(L).PERIOD*SQRT(G*HP(L)/FC2)   
            
            WV(L).K = MAX( 2.*PI/WV(L).LENGTH, 0.01 )
            WV(L).KHP = MIN(WV(L).K*HP(L),SHLIM)
            WV(L).UDEL = MAX(1D-6,PI*WV(L).HEIGHT/(WV(L).PERIOD*SINH(WV(L).KHP)))        ! *** ORBITAL VELOCITY
            
            ! *** WAVE DIRECTION (RADIANS) COUNTER-CLOCKWISE (CELL-EAST AXIS,WAVE)
            ! WV(L).DIR=(90-WDIR-ROTAT)*PI/180._8  
            WINX =  CVN(L)*WNDVELE(L) - CVE(L)*WNDVELN(L)                                ! *** CURVI-LINEAR SYS
            WINY = -CUN(L)*WNDVELE(L) + CUE(L)*WNDVELN(L)
            WV(L).DIR= ATAN2(WINY,WINX)                    ![-pi,pi] no rotation
            
          ELSE
            WV(L).HEIGHT = 0.
            WV(L).HEISIG = 0.
            WV(L).FREQ = 1.
            WV(L).DIR = 0.
            WV(L).LENGTH = 0.
            WV(L).UDEL = 0.
            WV(L).PERIOD = 0.
            WV(L).K  = 0.01
          ENDIF
        ENDIF  ! *** END OF VALID WAVE CELL
        
      ENDDO
    ENDDO    ! *** END OF DOMAIN
    !$OMP END DO
    
    !$OMP END PARALLEL
    
  END SUBROUTINE

  FUNCTION FETZONE8(WDIR) RESULT(ZONE)
    !DETERMINING FETCH ZONE AND FETCH MAIN ANGLE
    !BASED ON THE GIVEN WIND DIRECTION WDIR
    !WDIR     : INTERPOLATED WIND DIRECTION FROM WSER.INP
    !UNIT     : [0,360]
    !FORMATION: ANGLE BY (NORTH,WIND TO)IN CLOCKWISE DIRECTION
    !ZONE 1: NORTH       >337.5 OR  <= 22.5
    !ZONE 2: NORTH-EAST  >22.5  OR  <= 67.5
    !ZONE 3: EAST        >
    !ZONE 4: SOUTH-EAST
    !ZONE 5: SOUTH
    !ZONE 6: SOUTH-WEST
    !ZONE 7: WEST
    !ZONE 8: NORTH-WEST
    REAL(RKD) ,INTENT(IN ) :: WDIR   ![0,360]
    INTEGER :: ZONE

    IF     (WDIR>337.5 .OR. WDIR <= 22.5 )THEN
     ZONE = 1
    ELSEIF (WDIR>22.5 .AND. WDIR <= 67.5 )THEN
     ZONE = 2
    ELSEIF (WDIR>67.5 .AND. WDIR <= 112.5 )THEN
     ZONE = 3
    ELSEIF (WDIR>112.5 .AND. WDIR <= 157.5 )THEN
     ZONE = 4
    ELSEIF (WDIR>157.5 .AND. WDIR <= 202.5 )THEN
     ZONE = 5
    ELSEIF (WDIR>202.5 .AND. WDIR <= 247.5 )THEN
     ZONE = 6
    ELSEIF (WDIR>247.5 .AND. WDIR <= 292.5 )THEN
     ZONE = 7
    ELSEIF (WDIR>292.5 .AND. WDIR <= 337.5 )THEN
     ZONE = 8
    ENDIF
  END FUNCTION
 
  FUNCTION FETZONE(WDIR) RESULT(ZONE)
    !DETERMINING FETCH ZONE AND FETCH MAIN ANGLE
    !BASED ON THE GIVEN WIND DIRECTION WDIR
    !WDIR     : INTERPOLATED WIND DIRECTION FROM WSER.INP
    !UNIT     : [0,360]
    !FORMATION: ANGLE BY (NORTH,WIND TO)IN CLOCKWISE DIRECTION
    REAL(RKD) ,INTENT(IN ) :: WDIR   ![0,360]
    INTEGER :: ZONE

    IF     (WDIR>348.75 .OR. WDIR <= 11.25 )THEN
     ZONE = 1
    ELSEIF (WDIR>11.25 .AND. WDIR <= 33.75 )THEN
     ZONE = 2
    ELSEIF (WDIR>33.75 .AND. WDIR <= 56.25 )THEN
     ZONE = 3
    ELSEIF (WDIR>56.25 .AND. WDIR <= 78.75 )THEN
     ZONE = 4
    ELSEIF (WDIR>78.75 .AND. WDIR <= 101.25 )THEN
     ZONE = 5
    ELSEIF (WDIR>101.25 .AND. WDIR <= 123.75 )THEN
     ZONE = 6
    ELSEIF (WDIR>123.75 .AND. WDIR <= 146.25 )THEN
     ZONE = 7
    ELSEIF (WDIR>146.25 .AND. WDIR <= 168.75 )THEN
     ZONE = 8
    ELSEIF (WDIR>168.75 .AND. WDIR <= 191.25 )THEN
     ZONE = 9
    ELSEIF (WDIR>191.25 .AND. WDIR <= 213.75 )THEN
     ZONE = 10
    ELSEIF (WDIR>213.75 .AND. WDIR <= 236.25 )THEN
     ZONE = 11
    ELSEIF (WDIR>236.25 .AND. WDIR <= 258.75 )THEN
     ZONE = 12
    ELSEIF (WDIR>258.75 .AND. WDIR <= 281.25 )THEN
     ZONE = 13
    ELSEIF (WDIR>281.25 .AND. WDIR <= 303.75 )THEN
     ZONE = 14
    ELSEIF (WDIR>303.75 .AND. WDIR <= 326.25 )THEN
     ZONE = 15
    ELSEIF (WDIR>326.25 .AND. WDIR <= 348.75 )THEN
     ZONE = 16
    ENDIF
  END FUNCTION
  
  SUBROUTINE FETCH
    !DETERMINING THE FETCHES OF ALL CELLS:
    !OUTPUT: FWDIR(2:LA,1:NZONE) IN M
    USE DRIFTER,ONLY:INSIDECELL
    REAL(RKD) :: AL(NZONE),RL,XM,YM,RL0,DOTX,DOTY
    INTEGER :: I,J,L,NZ,IM,JM,LM,STATUS,LS,LN,ND,LF,LL,LP
    LOGICAL(4) :: ULOG,VLOG

    FWDIR = 0    
    AL  = (180+90-FETANG-ROTAT)*PI/180._8   !COUNTER-CLOCKWISE (TRUE EAST,WIND FR) 
    RL0 = 0.25*MIN(MINVAL(DXP(2:LA)),MINVAL(DYP(2:LA)))

    !$OMP PARALLEL DEFAULT(SHARED) 

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,I,J,NZ,IM,JM,LM,STATUS,LS,LN)  &
    !$OMP    PRIVATE(RL,XM,YM,DOTX,DOTY,ULOG,VLOG)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
    
      DO LP=LF,LL  
        L=LWET(LP)  
        IF( LWVMASK(L) )THEN
          DO NZ=1,NZONE
            RL=0
            IM=IL(L)
            JM=JL(L)
            
            LOOP1:DO WHILE(1)
              STATUS=0
              RL=RL+RL0                        !SEARCH FETCH FOR MIN (DX,DY)/4
              XM = XCOR(L,5)+RL*COS(AL(NZ))    !UPWIND DISTANCE=FETCH
              YM = YCOR(L,5)+RL*SIN(AL(NZ))
              LOOP2:DO J=JM-1,JM+1
               DO I=IM-1,IM+1
                 LM=LIJ(I,J)
                 IF( LM<2) CYCLE
                 IF( INSIDECELL(LM,XM,YM) )THEN
                   DOTX = CUE(L)*COS(AL(NZ))+CUN(L)*SIN(AL(NZ))
                   DOTY = CVE(L)*COS(AL(NZ))+CVN(L)*SIN(AL(NZ))
                   IF( ABS(DOTX)<1D-3 )THEN
                     DOTX = 0
                     DOTY = 1
                   ENDIF
                   IF( ABS(DOTX)>0.999 )THEN
                     DOTX = 1
                     DOTY = 0
                   ENDIF
                   IF( ABS(DOTY)<1D-3 )THEN
                     DOTX = 1
                     DOTY = 0
                   ENDIF
                   IF( ABS(DOTY)>0.999 )THEN
                     DOTX = 0
                     DOTY = 1
                   ENDIF
                   LS = LSC(LM)
                   LN = LNC(LM)
                   IF( NZ == 1 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LM) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM) == 1 .OR. UMASK(LM+1) == 1)
                   ELSEIF( NZ>1 .AND. NZ <= 4 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LM) == 1 .OR. VMASK(LM-1) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM) == 1 .OR. UMASK(LS) == 1)                                  
                   ELSEIF( NZ == 5 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LM) == 1 .OR. VMASK(LN  ) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM) == 1)
                   ELSEIF( NZ>5 .AND. NZ <= 8 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LN) == 1 .OR. VMASK(LN-1) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM) == 1 .OR. UMASK(LN) == 1)
                   ELSEIF( NZ == 9 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LN) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM) == 1 .OR. UMASK(LM+1) == 1)
                   ELSEIF( NZ>9 .AND. NZ <= 12 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LN) == 1 .OR. VMASK(LN+1) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM+1) == 1 .OR. UMASK(LN+1) == 1)           
                   ELSEIF( NZ == 13 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LM) == 1 .OR. VMASK(LN  ) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM+1) == 1)
                   ELSEIF( NZ>13 .AND. NZ <= 16 )THEN
                     ULOG = ABS(DOTX)/=1 .AND. (VMASK(LM) == 1 .OR. VMASK(LM+1) == 1)
                     VLOG = ABS(DOTY)/=1 .AND. (UMASK(LM+1) == 1 .OR. UMASK(LS+1) == 1)                 
                   ENDIF
                   
                   IF( ULOG .OR. VLOG )THEN
                     EXIT LOOP1
                   ELSE
                     STATUS=1
                   ENDIF
                   EXIT LOOP2
                 ENDIF
               ENDDO
             ENDDO LOOP2
             IF( STATUS == 0) EXIT LOOP1
             IM=IL(LM)
             JM=JL(LM)
            ENDDO LOOP1
            FWDIR(L,NZ)=RL
          ENDDO
        ENDIF
      ENDDO  !
    ENDDO    ! *** END OF DOMAIN
    !$OMP END DO
    
    !$OMP END PARALLEL

    ! *** SAVE THE FETCH OUT TO A FILE
    OPEN(UFET,FILE=OUTDIR//'FETCH.OUT')
    DO L=2,LA
      WRITE(UFET,'(2I6,16F15.4)') IL(L),JL(L),(FWDIR(L,NZ),NZ=1,NZONE)
    ENDDO
    CLOSE(UFET)
    
  END SUBROUTINE 
 
  SUBROUTINE WINDWAVEINIT   
    ! *** INITIALIZES WAVE VARIABLES AND GENERATES FETCH.OUT
    USE GLOBAL  
    INTEGER   :: L,K
    
    RSWRSR=FLOAT(ISWRSR)
    RSWRSI=FLOAT(ISWRSI)
    ROTAT   = 0
    WACCWE  = 0
    FWVTP   = 0

    ! *** Wave Computational Cell Defaults
    NWVCELLS = LA-1
    DO L=2,LA
      LWVCELL(L-1) = L
      LWVMASK(L) = .TRUE.
    ENDDO

    KSW = MAX(1D-6,KSW)
    WRITE(*,'(A)')'COMPUTING FETCH'
    CALL FETCH
    OPEN(UWIN,FILE=OUTDIR//'LIJXY.OUT',ACTION='WRITE')
    DO L=2,LA
      WRITE(UWIN,'(3I10,2F15.5)') L,IL(L),JL(L),DLON(L),DLAT(L)
    ENDDO
    CLOSE(UWIN)
    
    JSWRPH=1  
    
    CALL ZEROWAVE
   
    ITWCBL1=1  
    ITWCBL2=0  
    ITWCBL3=0  
    
    !****************************************************************
    !OPEN(1,FILE='WAVE.INP',STATUS='UNKNOWN')
    !READ(1,*,IOSTAT=ISO)NWVCELLS,WVPRD,ISWCBL,ISWRSR,ISWRSI,NWUPDT, &
    !    NTSWV,WVDISV,WVDISH,WVLSH,WVLSX,ISWVSD,ISDZBR

    !---------------------------------------------------------------
    ! **  INITIALIZE VERTICAL DISTRIBUTION OF WAVE DISSIPATION AS SOURCE                                                    
    ! **  TO VERTICAL TKE CLOSURE                                                                                           
    IF( KC == 2 )THEN
      WVDTKEM(1)=WVDISV
      WVDTKEP(1)=WVDISV
    ENDIF
    IF( KC == 3 )THEN
      WVDTKEM(1)=WVDISV
      WVDTKEP(1)=0.5*WVDISV
      WVDTKEM(2)=0.5*WVDISV
      WVDTKEP(2)=WVDISV
    ENDIF
    IF( KC >= 4 )THEN
      WVDTKEM(1)=WVDISV
      WVDTKEP(1)=0.5*WVDISV
      WVDTKEM(KS)=0.5*WVDISV
      WVDTKEP(KS)=WVDISV
      DO K=2,KS-1
        WVDTKEM(K)=0.5*WVDISV
        WVDTKEP(K)=0.5*WVDISV
      ENDDO
    ENDIF

    DO L=1,LC
      ZBRE(L)=KSW  ! *** INPUT NIKURADSE ROUGHNESS
    ENDDO
    
  END SUBROUTINE
  
  SUBROUTINE ZEROWAVE
  INTEGER :: L,K
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
    QQWC(L)=1.E-12   ! *** NOT USED
    QQWCR(L)=1.E-12  
    QQWV1(L)=1.E-12  ! *** BED TURBULENT INTENSITY DUE TO WAVES ONLY 
    QQWV2(L)=0.      ! *** WATER COLUMN TURBULENT INTENSITY DUE TO WAVES  
    QQWV3(L)=1.E-12  ! *** WATER COLUMN TURBULENT INTENSITY DUE TO WAVES MODIFIED FOR NON-COHESIVE MOVING BED
    WV(L).HEIGHT=0.  
    WV(L).K=0.  
    WV(L).DIR=0.  
  ENDDO  
  DO K=1,KC  
    DO L=1,LC  
     WVHUU(L,K)=0.  
     WVHVV(L,K)=0.  
     WVHUV(L,K)=0.  
     WVPP(L,K)=0.  
     WVPU(L,K)=0.  
     WVPV(L,K)=0.  
     WV(L).DISSIPA(K)=0.  
     FXWAVE(L,K)=0.  
     FYWAVE(L,K)=0.  
    ENDDO  
  ENDDO  
END SUBROUTINE
 
! **********************************************************************************************
! *** Subroutine to compute wind wave generated currents
! ***
SUBROUTINE W_WAVESXY
  REAL :: DISPTMP
  REAL :: TMPVAL,WVWHAUT,WVWHAVT
  REAL(RKD) :: SNHTOPU,SNHBOTU,SNHTOPV,SNHBOTV,TMPPU,TMPPV,TMP
  REAL(RKD) :: WG,WN,WK,RKHM1,RKHM2,TMPP1,TMPP2
  REAL(RKD) :: COSH3,RATIO,ZTOP,ZBOT,COSHTOP,COSHBOT,SINHTOP,SINHBOT,SINH2
  
  INTEGER :: K,L,LS,LW,LE,LSW,LN,LNW,LSE,ND,LF,LL,LP

  !$OMP PARALLEL DEFAULT(SHARED) 
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,DISPTMP,WK,WG,WN)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    DO LP=LF,LL  
      L=LWET(LP)  
      IF( LWVMASK(L) )THEN
        IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRYWAV )THEN
          WVENEP(L)= G*WV(L).HEIGHT**2./16.                       ! *** ENERGY/RHO (M3/S2) FOR RANDOM WAVE
          DISPTMP =  0.25*G*WV(L).HEIGHT**2./WV(L).PERIOD         ! *** ENERGY DISSIPATION DUE TO BREAKING WAVE   
          WV(L).DISSIPA(KC)= 0.01*DISPTMP
          WK = 2.0*WV(L).KHP                                      
          WG = WK/SINH(WK)                                  
          WN = (WG+1)/2                                     
          WVHUU(L,KC) = WVENEP(L)*(WG/2+WN*COS(WV(L).DIR)**2)     ! *** SXXTMP (M3/S2)
          WVHVV(L,KC) = WVENEP(L)*(WG/2+WN*SIN(WV(L).DIR)**2)     ! *** SYYTMP (M3/S2)
          WVHUV(L,KC) = WVENEP(L)*WN/2*SIN(2*WV(L).DIR)           ! *** SXYTMP (M3/S2)
          HMPW(L) = HP(L)+WV(L).HEIGHT
        ELSE
          WVENEP(L)=0
          WV(L).DISSIPA(KC)=0
          WVHUU(L,KC) =0
          WVHVV(L,KC) =0
          WVHUV(L,KC) =0
          HMPW(L) =HP(L)
        ENDIF
      ENDIF
    ENDDO
  ENDDO    ! *** END OF DOMAIN
  !$OMP END DO

  ! *** COMPUTE THE DERIVED WATER SURFACES                                                                                
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,LS,LW,LSW)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    DO LP=LF,LL  
      L=LWET(LP)  
      IF( LWVMASK(L) )THEN
        LS=LSC(L)
        LW=LWC(L)
        LSW=LSWC(L)
        HMUW(L) = 0.5*(DXYP(L)*HMPW(L)+DXYP(LW)*HMPW(LW))/(DXU(L)*DYU(L))
        HMVW(L) = 0.5*(DXYP(L)*HMPW(L)+DXYP(LS)*HMPW(LS))/(DXV(L)*DYV(L))
      ENDIF
    ENDDO
  ENDDO    ! *** END OF DOMAIN
  !$OMP END DO
  
  ! *** DISTRIBUTE VALUES ACROSS KC                                                                                       
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K) &
  !$OMP    PRIVATE(RKHM1,RKHM2,SINH2,COSH3,RATIO,ZTOP,ZBOT)  &
  !$OMP    PRIVATE(SINHTOP,SINHBOT,COSHTOP,COSHBOT,TMPVAL,TMP,TMPP1,TMPP2)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    DO LP=LF,LL  
      L=LWET(LP)  
      IF( LWVMASK(L) )THEN
        IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRY )THEN
          RKHM1 = WV(L).KHP
          RKHM2 = 2.*RKHM1
          SINH2 = SINH(RKHM2)
          COSH3 = COSH(RKHM1)
          RATIO = 0.5+(RKHM1/SINH2)
          DO K=KSZ(L),KC
            ZTOP = Z(L,K)
            ZBOT = Z(L,K-1)
            ! *** START Moved to after ZTOP is defined (2009_09_03)
            SINHTOP = SINH(RKHM2*ZTOP)
            SINHBOT = SINH(RKHM2*ZBOT)
            COSHTOP = COSH(RKHM1*ZTOP)
            COSHBOT = COSH(RKHM1*ZBOT)
            TMPVAL = (RKHM2*(ZTOP-ZBOT)+SINHTOP-SINHBOT)/(RKHM2+SINH2)
            
            ! *** APPLY FACTOR
            WVHUU(L,K) = TMPVAL*WVHUU(L,KC)
            WVHVV(L,K) = TMPVAL*WVHVV(L,KC)
            WVHUV(L,K) = TMPVAL*WVHUV(L,KC)
            WV(L).DISSIPA(K) = TMPVAL*WV(L).DISSIPA(KC)

            TMPP1 = -0.5*(ZTOP-ZBOT)+(ZTOP*COSHTOP-ZBOT*COSHBOT)/COSH3
            TMP = SINH2-2.
            IF( ABS(TMP) < 1D-3 ) TMP = SIGN(1D-3,TMP)
            TMPP2 = (RATIO-1.)*(SINHTOP-SINHBOT-2.*(ZTOP-ZBOT))/TMP
            
            ! *** LIMIT RANGE WHEN WV.K~0.72
            IF( ABS(TMPP1) > 0.5 )THEN
              TMPP1 = SIGN(0.5,TMPP1)
            ENDIF
            IF( ABS(TMPP2) > 0.5 )THEN
              TMPP2 = SIGN(0.5,TMPP2)
            ENDIF
            WVPP(L,K) = WVENEP(L)*(TMPP1+TMPP2)
          ENDDO
        ELSE
          WVHUU(L,1:KC) = 0
          WVHVV(L,1:KC) = 0
          WVHUV(L,1:KC) = 0
          WV(L).DISSIPA(1:KC) = 0
          WVPP(L,1:KC)=0
        ENDIF
      ENDIF
    ENDDO
  ENDDO    ! *** END OF DOMAIN
  !$OMP END DO

  ! **  The Updated of QQWV2 was moved to CALTBXY
  
  !IF(ISRESTI/=0 )THEN
  !  !$OMP SINGLE
  !  WRITE(*,'(A)')'WAVE: WVQWCP.INP'
  !  OPEN(1,FILE='wvqwcp.inp',STATUS='UNKNOWN')
  !  DO L=2,LA
  !    READ(1,*)IDUM,JDUM,QQWV1(L),QQWV2(L),QQWV2(L),QQWC(L),QQWCR(L)
  !  ENDDO
  !  !$OMP END SINGLE
  !ENDIF

  ! **  COMPUTE CELL FACE QUANTITIES WVPU,WVPV                                                                            
  !$OMP DO PRIVATE(ND,LF,LL,LP,L)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    DO LP=LF,LL  
      L=LWET(LP)  
      IF( LWVMASK(L) )THEN
        IF( WV(L).HEIGHT >= WHMI .AND. HP(L)>HDRY )THEN
          WVKHU(L) = MIN(HMUW(L)*WV(L).K,SHLIM)
          WVKHV(L) = MIN(HMVW(L)*WV(L).K,SHLIM)
        ELSE
          WVKHU(L) = 1.D-12
          WVKHV(L) = 1.D-12
        ENDIF
      ENDIF
    ENDDO
  ENDDO    ! *** END OF DOMAIN
  !$OMP END DO

  !$OMP DO PRIVATE(ND,LF,LL,LP,L,LS)  &
  !$OMP    PRIVATE(TMPVAL,WVWHAUT,WVWHAVT)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    DO LP=LF,LL  
      L=LWET(LP)  
      IF( LWVMASK(L) )THEN
        TMPVAL    = 0.5*WV(L).FREQ*WV(L).FREQ
        WVTMP1(L) = MAX(SINH(WVKHU(L)),1D-6)
        WVWHAUT   = (WV(L).HEIGHT + SUB(L)*WV(LWC(L)).HEIGHT)/(1.+SUB(L))
        WVTMP2(L) = TMPVAL*WVWHAUT*WVWHAUT/(WVTMP1(L)*WVTMP1(L))
        WVWHAVT   = (WV(L).HEIGHT + SVB(L)*WV(LSC(L)).HEIGHT)/(1.+SVB(L))
        WVTMP3(L) = MAX(SINH(WVKHV(L)),1D-6)
        WVTMP4(L) = TMPVAL*WVWHAVT*WVWHAVT/(WVTMP3(L)*WVTMP3(L))
      ENDIF
    ENDDO
  ENDDO    ! *** END OF DOMAIN
  !$OMP END DO

  !$OMP DO PRIVATE(ND,K,LP,L)  &
  !$OMP    PRIVATE(ZTOP,ZBOT,SNHTOPU,SNHBOTU,SNHTOPV,SNHBOTV,TMPPU,TMPPV)
  DO ND=1,NDM  
    DO K=1,KC
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IF( LWVMASK(L) )THEN
          ZTOP=Z(L,K)
          ZBOT=Z(L,K-1)
          SNHTOPU=SINH(WVKHU(L)*ZTOP)
          SNHBOTU=SINH(WVKHU(L)*ZBOT)
          SNHTOPV=SINH(WVKHV(L)*ZTOP)
          SNHBOTV=SINH(WVKHV(L)*ZBOT)
          TMPPU=(1.-ZTOP)*SNHTOPU*(ZTOP*WVTMP1(L)-SNHTOPU)-(1.-ZBOT)*SNHBOTU*(ZBOT*WVTMP1(L)-SNHBOTU)
          TMPPV=(1.-ZTOP)*SNHTOPV*(ZTOP*WVTMP3(L)-SNHTOPV)-(1.-ZBOT)*SNHBOTV*(ZBOT*WVTMP3(L)-SNHBOTV)
          WVPU(L,K)=WVTMP2(L)*TMPPU
          WVPV(L,K)=WVTMP4(L)*TMPPV
        ENDIF
      ENDDO
    ENDDO
  ENDDO    ! *** END OF DOMAIN
  !$OMP END DO

  ! **  CALCULATE THE NET X AND  Y WAVE REYNOLDS STRESS FORCINGS                                                          
  !$OMP DO PRIVATE(ND,K,LP,L,LW,LE,LS,LN,LNW,LSE)  
  DO ND=1,NDM  
    DO K=1,KC    
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IF( LWVMASK(L) )THEN
          LW=LWC(L)
          LE=LEC(L)
          LS=LSC(L)
          LN=LNC(L)
          LNW=LNWC(L)
          LSE=LSEC(L)
          FXWAVE(L,K)=DZIC(L,K)*SUB(L)*SPB(L) * &
                      ( RSWRSI*(DYU(L)*(WVPP(L,K) - WVPP(LW,K))         + DYU(L)*WVPU(L,K)*(HMPW(L) - HMPW(LW)))   &
                       +RSWRSR*(DYP(L)*WVHUU(L,K) - DYP(LW)*WVHUU(LW,K) + 0.5*(DXV(LN) + DXV(LNW))*WVHUV(LN,K)     &
                       -0.5*(DXV(L)+DXV(LW))*WVHUV(L,K)) )
          FYWAVE(L,K)=DZIC(L,K)*SVB(L)*SPB(L) * &
                      ( RSWRSI*(DXV(L)*(WVPP(L,K) - WVPP(LS,K))         + DXV(L)*WVPV(L,K)*(HMPW(L) - HMPW(LS )))  &
                       +RSWRSR*(DXP(L)*WVHVV(L,K) - DXP(LS)*WVHVV(LS,K) + 0.5*(DYU(LE) + DYU(LSE))*WVHUV(LE,K)     &
                       -0.5*(DYU(L)+DYU(LS))*WVHUV(L,K)) )         
        ENDIF
      ENDDO
    ENDDO
  ENDDO    ! *** END OF DOMAIN
  !$OMP END DO
  !$OMP END PARALLEL
  
END SUBROUTINE

SUBROUTINE READWAVECELLS
  ! *** READS THE FILE WAVECELLS.INP TO DEFINE A SUBSET OF CELLS TO USE FOR WAVE COMPUTATIONS
  !
  INTEGER :: L, IWV, JWV, ISO, NWV

  ! *** DEFAULT IS OFF
  NWVCELLS = 0
  DO NWV = 1,LA
    LWVMASK(NWV) = .FALSE.
    LWVCELL(NWV) = 0
  ENDDO
  
  ! *** READ THE WAVE COMPUTATIONAL CELL LIST
  WRITE(*,'(A)')'READING WAVECELLS.INP'
  OPEN(1,FILE='wavecells.inp',STATUS='OLD')

  ! *** SKIP HEADERS, IF ANY
  CALL SKIPCOM(1,'*')
  
  ! *** LOOP OVER THE CELLS
  DO NWV=2,LA
    READ(1,*,IOSTAT=ISO,END=1000) IWV, JWV
    IF( ISO > 0 ) GOTO 9000
    
    L=LIJ(IWV,JWV)
    IF( L > 1 )THEN
      NWVCELLS = NWVCELLS+1
      LWVCELL(NWVCELLS) = L
      LWVMASK(L) = .TRUE.
    ENDIF
  ENDDO

1000 CLOSE(1)  
  RETURN
  
  9000 STOP ' ERROR READING WAVECELLS.INP.'
  
END SUBROUTINE 
 

END MODULE
