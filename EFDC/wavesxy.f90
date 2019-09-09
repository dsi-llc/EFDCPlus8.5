SUBROUTINE WAVESXY 
  
  ! *** WAVESXY subroutine links and processes externally generated wave data for
  
  ! *** * First Line: Initial Data Settings (1 line)
  ! *** NWVCELLS  - Number Of Cells With Wave Data                                                                      
  ! *** WVPRD   - Wave Period (sec)                                                                                   
  ! *** ISWCBL  - 1 Activates Wave-Current BL Model                                                                   
  ! *** ISWRSR  - 1 Activates Inclusion Of Rotational Component Of Rad Stress                                        
  ! *** ISWRSI  - 1 Activates Inclusion Of Irrotational Component Of Rad Stress                                      
  ! *** NTSWV   - Number Of Time Steps For Gradual Introduction Of Wave Forcing                                      
  ! *** WVDISV  - Fraction Of Wave Dissipation As Source In Vertical TKE Closure                                     
  ! *** WVDISH  - Fraction Of Wave Dissipation As Source In Horiz Smagorinky's Subgrid Closure  (NOT USED)                       
  ! *** WVLSH   - Weight For Depth As The Horiz SSG Eddy Viscosity Length Scale,       ISHMDF>0                                      
  ! *** WVLSX   - Weight For Sqrt(Dxdy) As The Horiz SSG Eddy Viscosity Length Scale,  ISHMDF>0
  ! *** ISWVSD  - 1 Include Nondiverg Wave Stokes Drift In Mass Transport                                             
  ! *** ISDZBR  - 1 Write Diagnostics For Effect Wave Current Bndry Layer Roughness                                  
  ! ***
  ! *** * Second Line: (repeated NWVCELLS times)                                                                        
  ! *** IWV     - I Of The Wave Cell                                                                                  
  ! *** JWV     - J Of The Wave Cell                                                                                  
  ! *** ENETMP  - Wave Energy, 0.5*G*Abs(Amplitude)**2  (M3/S2)                                                       
  ! *** SXXTMP  - Radiation Stresses- XX Component (kg/s^2)                                                           
  ! *** SYYTMP  - Radiation Stresses- YY Component (kg/s^2)                                                           
  ! *** SXYTMP  - Radiation Stresses- XY Component (kg/s^2)                                                           
  ! *** WVDISP  - Energy Dissipation  (M3/S3)                        [WV.DISSIPA]                                                                       
  ! *** WANGLE  - Wave Angle Measured from East positive CCW (deg)   
                                                                                   
  ! CHANGE RECORD                                                                                                     
  ! DATE MODIFIED     BY               DESCRIPTION                                                                    
  !-- -------------------------------------------------------------------!   
  ! 2017-05-23        Dang Chung       IMPROVED WAVE CALCULATION
  ! 2011-04-26        Dang Chung       Corrected Wave Formulation & Orbital Velocity                                 
  ! 2011-03-02        PAUL M. CRAIG    RESTRUCTURED AND CORRECTED CODE     
  ! 2012-08-04        Dang Chung       IMPROVED WAVE INPUT & CALCULATION 

  USE GLOBAL
  USE GETSWANMOD
  USE WAVELENGTH
  
  IMPLICIT NONE
  
  REAL :: UWORBIT,AEXTMP
  REAL :: WWVH, WANGLE, WVLEN,DISPTMP
  REAL :: TAUTMP,CORZBR,CDRGTMP
  REAL :: WVWHAUT,WVWHAVT,DZITMP
  REAL(RKD) :: SNHBOTU,SNHTOPV,SNHBOTV,TMPPU,TMPPV,SNHTOPU
  REAL(RKD) :: WG,WN,WK,WE,WDEP1,RRLS,WPRD
  REAL(RKD) :: SXXTMP,SXYTMP,SYYTMP,RATIO,RKHM1,RKHM2
  REAL(RKD) :: ZTOP,ZBOT,SINHTOP,SINHBOT,SINH2,COSH3
  REAL(RKD) :: COSHTOP,COSHBOT,TMPVAL,TMPP1,TMPP2,TMP
  REAL, EXTERNAL :: CSEDVIS

 
  INTEGER :: NW,L,K,IOS,NWV,IWV,JWV,LS,LW,LSW,LE
  INTEGER :: LN,LNW,LSE,ND,LF,LL
  
  ! **  INPUT WAVE INFORMATION                                                                                            
  IF( JSWAVE == 0 )THEN
    JSWRPH=1
    RSWRSR=FLOAT(ISWRSR)
    RSWRSI=FLOAT(ISWRSI)
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
      QQWV1(L)=1.E-12  ! *** BED TURBULENT INTENSITY DUE TO WAVES ONLY
      QQWV2(L)=1.E-12  ! *** WATER COLUMN TURBULENT INTENSITY DUE TO WAVES
      QQWV3(L)=1.E-12  ! *** WATER COLUMN TURBULENT INTENSITY DUE TO WAVES MODIFIED FOR NON-COHESIVE MOVING BED
    ENDDO
    DO K=1,KC
      DO L=1,LC
        WVHUU(L,K) =0.
        WVHVV(L,K) =0.
        WVHUV(L,K) =0.
        WVPP(L,K)  =0.
        WVPU(L,K)  =0.
        WVPV(L,K)  =0.
        FXWAVE(L,K)=0.
        FYWAVE(L,K)=0.
      ENDDO
    ENDDO
    
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
        WVDTKEM(K)=0.5*WVDISV  ! *** BOTTOM
        WVDTKEP(K)=0.5*WVDISV  ! *** TOP
      ENDDO
    ENDIF
    
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
       
    ! *** INPUT THE WAVE FILE: FIRST READ 
    IF( IFWAVE == 0 )THEN
      WRITE(*,'(A)')'WAVE: READING WAVE.INP'
      OPEN(WUNI,FILE='wave.inp',STATUS='UNKNOWN')
      CALL SKIPCOM(WUNI,'*')
      DO NW=1,IWVCOUNT-1
        DO NWV=2,LA
          READ(WUNI,*,IOSTAT=IOS)IWV,JWV,WWVH,WANGLE,WVPRD,WVLEN,DISPTMP
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
    CALL HDEP0_PLSWAVE
    
  ELSEIF( JSWAVE == 1 .AND. IWVCOUNT<NWVTIM .AND. TIMEDAY >= WAVEDAY(IWVCOUNT+1) )THEN
    ! ** WAVE UPDATE
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
    CALL HDEP0_PLSWAVE
  ENDIF
  ! ** END OF WAVE UPDATE *************************************************
  
  ! **  DISTRIBUTE WVHUU, WVHVV, AND WV.DISSIPA OVER DEPTH            
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(ND,L,LF,LL,WDEP1,WPRD,RRLS,WK,WE,WG,WN,SXXTMP,SXYTMP,SYYTMP)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
    
    DO L=LF,LL
      WV(L).HEIGHT = MIN(0.75*HP(L),WV(L).HEISIG)            ! *** INCLUDING BREAKING WAVE
      IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRYWAV )THEN
        WDEP1 = HP(L)
        WPRD  = 2.*PI/WV(L).FREQ
        IF( WVLCAL == 1 )THEN
          CALL BISEC(DISRELATION,WLMIN,WLMAX,EPS,WPRD,WDEP1,0._8,0._8,RRLS)
          WV(L).LENGTH = RRLS
        ENDIF
        WV(L).KHP = MIN(WV(L).K*HP(L),SHLIM)
        WK = 2.0*WV(L).KHP                             ! *** 4*PI/WV(L).LENGTH*HP(L) = 2KH
        WE = 9.81*1000.*WV(L).HEIGHT**2 /16.           ! *** TOTAL WAVE ENERGY : KG/S2  
        WVENEP(L) = WE/1000.                           ! *** WAVE ENERGY/RHO: M3/S2
        WG = WK/SINH(WK)                               ! *** 2KH/SINH(2KH)
        WN = (WG+1)/2                               
        SXXTMP = WE*(WG/2+WN*COS(WV(L).DIR)**2)        ! *** RADIATION SHEAR STRESS [Kg/S2] 
        SXYTMP = WE*WN/2*SIN(2*WV(L).DIR)              ! *** RADIATION SHEAR STRESS [Kg/S2] 
        SYYTMP = WE*(WG/2+WN*SIN(WV(L).DIR)**2)        ! *** RADIATION SHEAR STRESS [Kg/S2] 
        WVHUU(L,KC) = SXXTMP/1000.
        WVHVV(L,KC) = SYYTMP/1000.
        WVHUV(L,KC) = SXYTMP/1000.
      ELSE
        WVHUU(L,KC)=0.
        WVHVV(L,KC)=0.
        WVHUV(L,KC)=0.    
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO
  
  ! *** DISTRIBUTE VALUES ACROSS KC    
  !$OMP DO PRIVATE(ND,L,LF,LL,K,RKHM1,RKHM2,SINH2,COSH3,RATIO,ZTOP,ZBOT) &
  !$OMP    PRIVATE(SINHTOP,SINHBOT,COSHTOP,COSHBOT,TMPVAL,TMPP1,TMP,TMPP2)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
    DO L=LF,LL 
      IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRYWAV )THEN
        RKHM1 = WV(L).KHP
        RKHM2 = 2.*RKHM1
        SINH2 = SINH(RKHM2)
        COSH3 = COSH(RKHM1)
        RATIO = 0.5+(RKHM1/SINH2)
        DO K=1,KC
          ZTOP = Z(L,K)
          ZBOT = Z(L,K-1)
          SINHTOP = SINH(RKHM2*ZTOP)
          SINHBOT = SINH(RKHM2*ZBOT)
          COSHTOP = COSH(RKHM1*ZTOP)
          COSHBOT = COSH(RKHM1*ZBOT)
          TMPVAL  = (RKHM2*(ZTOP-ZBOT)+SINHTOP-SINHBOT)/(RKHM2+SINH2)

          ! *** APPLY FACTOR
          WVHUU(L,K) = TMPVAL*WVHUU(L,KC)
          WVHVV(L,K) = TMPVAL*WVHVV(L,KC)
          WVHUV(L,K) = TMPVAL*WVHUV(L,KC)
          WV(L).DISSIPA(K) = TMPVAL*WV(L).DISSIPA(KC)  ! FROM GETSWANMOD
          
          TMPP1 = -0.5*(ZTOP-ZBOT)+(ZTOP*COSHTOP-ZBOT*COSHBOT)/COSH3
          TMP   = SINH2-2.
          IF( ABS(TMP)<1D-3) TMP = SIGN(1D-3,TMP)
          TMPP2 = (RATIO-1.)*(SINHTOP-SINHBOT-2.*(ZTOP-ZBOT))/TMP

          ! *** LIMIT RANGE WHEN WV.K~0.72
          IF( ABS(TMPP1)>0.5 )THEN
            TMPP1=SIGN(0.5,TMPP1)
          ENDIF
          IF( ABS(TMPP2)>0.5 )THEN
            TMPP2=SIGN(0.5,TMPP2)
          ENDIF
          WVPP(L,K) = WVENEP(L)*(TMPP1+TMPP2)
        ENDDO
      ELSE
        WVHUU(L,1:KC) = 1D-6
        WVHVV(L,1:KC) = 1D-6
        WVHUV(L,1:KC) = 1D-6
        WV(L).DISSIPA(1:KC) = 1D-6
        WVPP(L,1:KC)  = 1D-6
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO
  
  ! **  INITIALIZE WAVE-CURRENT BOUNDARY LAYER MODEL CALCULATING                                                          
  ! **  THE WAVE TURBULENT INTENSITY, QQWV                                                                                
  ! **  AND SQUARED HORIZONTAL WAVE ORBITAL VELOCITY MAGNITUDE            
  !$OMP DO PRIVATE(ND,L,LF,LL,AEXTMP,UWORBIT) &
  !$OMP    PRIVATE(TMPVAL,TAUTMP,CORZBR,CDRGTMP)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
    
    DO L=LF,LL
      ! *** SET ZBRE AS GRAIN/SKIN ROUGHNESS (M)
      IF( ISTRAN(7) > 0 )THEN
        ! *** BASE ON NIKURADSE ROUGHNESS (APPROXIMATE 2.5*D50)
        !ZBRE(L)=SEDDIA50(L,KBT(L)) !*2.5/30.
        ZBRE(L) = MAX(SEDDIA50(L,KBT(L)),1E-6)*2.5
      ELSE
        ZBRE(L) = KSW
      ENDIF

      IF( WV(L).HEIGHT >= WHMI .AND. HP(L) > HDRYWAV )THEN

        ! *** WAVE HEIGHT / HYPERBOLIC SINE OF WAVE NUMBER
        AEXTMP = 0.5*WV(L).HEIGHT/SINH(WV(L).KHP)

        ! *** SQUARED HORIZONTAL WAVE ORBITAL VELOCITY MAGNITUDE
        UWORBIT  = AEXTMP*WV(L).FREQ  !WVFRQ
        UWVSQ(L) = UWORBIT*UWORBIT

        IF( UWORBIT >= 1.E-6 )THEN
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
        ELSE
          QQWV2(L)=QQLMIN
        ENDIF
      ELSE
        UWVSQ(L)=0
        WV(L).FREQ=1
        QQWV2(L)=QQLMIN
      ENDIF
    ENDDO
  ENDDO
  !$OMP END DO
  
  ! **  COMPUTE CELL FACE QUANTITIES WVPU,WVPV         
  !$OMP DO PRIVATE(ND,L,LF,LL)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
    DO L=LF,LL
      IF( WV(L).HEIGHT >= WHMI ) WVKHU(L)= MIN(HMUW(L)*WV(L).K,SHLIM)     !VALKH(HFFDG)
      IF( WV(L).HEIGHT >= WHMI ) WVKHV(L)= MIN(HMVW(L)*WV(L).K,SHLIM)     !VALKH(HFFDG)
    ENDDO
  ENDDO
  !$OMP END DO
  
  !$OMP DO PRIVATE(ND,L,LF,LL,LS,TMPVAL,WVWHAUT,WVWHAVT)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
    DO L=LF,LL
      LS=LSC(L)
      TMPVAL    = 0.5*WV(L).FREQ*WV(L).FREQ     !WVFRQ*WVFRQ
      WVTMP1(L) = MAX(SINH(WVKHU(L)),1D-6)
      WVWHAUT   = (WV(L).HEIGHT + SUB(L)*WV(LWC(L)).HEIGHT)/(1.+SUB(L))
      WVTMP2(L) = TMPVAL*WVWHAUT*WVWHAUT /(WVTMP1(L)*WVTMP1(L))
      WVWHAVT   = (WV(L).HEIGHT + SVB(L)*WV(LSC(L)).HEIGHT)/(1.+SVB(L))
      WVTMP3(L) = MAX(SINH(WVKHV(L)),1D-6)
      WVTMP4(L) = TMPVAL*WVWHAVT*WVWHAVT /(WVTMP3(L)*WVTMP3(L))
    ENDDO
  ENDDO
  !$OMP END DO
  
  !$OMP DO PRIVATE(ND,LF,LL,K,L,ZTOP,ZBOT,SNHTOPU,SNHBOTU,SNHTOPV,SNHBOTV,TMPPU,TMPPV)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)
    DO K=1,KC
      DO L=LF,LL
        ZTOP = Z(L,K)
        ZBOT = Z(L,K-1)
        SNHTOPU = SINH(WVKHU(L)*ZTOP)
        SNHBOTU = SINH(WVKHU(L)*ZBOT)
        SNHTOPV = SINH(WVKHV(L)*ZTOP)
        SNHBOTV = SINH(WVKHV(L)*ZBOT)
        TMPPU = (1.-ZTOP)*SNHTOPU*(ZTOP*WVTMP1(L)-SNHTOPU) - (1.-ZBOT)*SNHBOTU*(ZBOT*WVTMP1(L)-SNHBOTU)
        TMPPV = (1.-ZTOP)*SNHTOPV*(ZTOP*WVTMP3(L)-SNHTOPV) - (1.-ZBOT)*SNHBOTV*(ZBOT*WVTMP3(L)-SNHBOTV)
        WVPU(L,K) = WVTMP2(L)*TMPPU
        WVPV(L,K) = WVTMP4(L)*TMPPV
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! **  CALCULATE THE INTERNAL MODE NET X AND Y WAVE REYNOLDS STRESS FORCINGS                                             
  !$OMP DO PRIVATE(ND,LF,LL,K,L,LS,LN,LW,LE,LNW,LSE,DZITMP)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)

    DO K=1,KC
      DO L=LF,LL
        IF( LKSZ(L,K) )CYCLE
        DZITMP = 1.*DZIC(L,K)
        LS = LSC(L)
        LN = LNC(L)
        LW = LWC(L)
        LE = LEC(L)
        LNW = LNWC(L)
        LSE = LSEC(L)
        FXWAVE(L,K) = DZITMP*SUB(L)*SPB(L)*( RSWRSI*(DYU(L)*(WVPP(L,K) - WVPP(LW,K)) + DYU(L)*WVPU(L,K)*(HMPW(L)-HMPW(LW)))    &
                                            +RSWRSR*(DYP(L)*WVHUU(L,K) - DYP(LW)*WVHUU(LW,K)                                   &
                                                    + 0.5*(DXV(LN)+DXV(LNW))*WVHUV(LN,K) - 0.5*(DXV(L)+DXV(LW))*WVHUV(L,K)) )
      
        FYWAVE(L,K) = DZITMP*SVB(L)*SPB(L)*( RSWRSI*(DXV(L)*(WVPP(L,K) - WVPP(LS,K)) + DXV(L)*WVPV(L,K)*(HMPW(L)-HMPW(LS )))   &
                                            +RSWRSR*(DXP(L)*WVHVV(L,K) - DXP(LS)*WVHVV(LS,K)                                   &
                                                    + 0.5*(DYU(LE)+DYU(LSE))*WVHUV(LE,K) - 0.5*(DYU(L)+DYU(LS))*WVHUV(L,K)) )
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  !$OMP END PARALLEL

CONTAINS

! *** 
SUBROUTINE HDEP0_PLSWAVE
  ! *** COMPUTE THE DERIVED WATER DEPTHS AT U/V CELL FACES
  DO L=2,LA
    LS=LSC(L)
    LW=LWC(L)
    LSW=LSWC(L)
    HMPW(L) = HP(L) + WV(L).HEIGHT
    
    ! HEIGHT U AND V FACE
    HMUW(L) = 0.5*( DXYP(L)*HMPW(L) + DXYP(LW)*HMPW(LW) )/(DXU(L)*DYU(L))
    HMVW(L) = 0.5*( DXYP(L)*HMPW(L) + DXYP(LS)*HMPW(LS) )/(DXV(L)*DYV(L))
  ENDDO
  END SUBROUTINE
  
END SUBROUTINE

SUBROUTINE GETWAVEINP
  
  USE GLOBAL
  
  ! *** READ WAVE PARAMS FROM WAVE.INP
  INTEGER :: L,IOS,NWV,IWV,JWV
  REAL   :: WVCX,WVCY,WVDX,WVDY,DISPTMP
  REAL   :: WWVH,WANGLE,WVLEN
    
  NWVCELLS = 0
  DO NWV=2,LA
    READ(WUNI,*,IOSTAT=IOS)IWV,JWV,WWVH,WANGLE,WVPRD,WVLEN,DISPTMP
    IF( IOS > 0 )THEN
      WRITE(*,'(A45,I5)') '***  READ ERROR ON FILE WAVE.INP AT CELL L = ',NWV
      STOP
    ENDIF
    L=LIJ(IWV,JWV)
    IF( L > 1 .AND. L <= LA )THEN
      NWVCELLS = NWVCELLS+1
      LWVCELL(NWVCELLS) = L   ! *** PMC - 2015-12 - This was pointing to L-1 but LWVCELL is used as an L array reference.  This has to be L not L-1.
      LWVMASK(L) = .TRUE.

      ! *** WAVE HEIGHT (M) THEN ZERO FOR VEGETATION EFFECT
      WV(L).HEISIG = WWVH 
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
      
      IF( WV(L).HEISIG >= WHMI .AND. WVPRD > 0. )THEN      
        WVDX= COS(WANGLE*PI/180)
        WVDY= SIN(WANGLE*PI/180)       
        WVCX =  CVN(L)*WVDX - CVE(L)*WVDY       
        WVCY = -CUN(L)*WVDX + CUE(L)*WVDY
        WV(L).DIR    = ATAN2(WVCY,WVCX)        
        WV(L).LENGTH = WVLEN     
        WV(L).FREQ   = 2.*PI/WVPRD
        WV(L).DISSIPA(KC) = DISPTMP                        ! *** ENERGY DISSIPATION DUE TO BREAKING WAVE (M3/S3)
        WV(L).K = MAX( 2.*PI/WV(L).LENGTH, 0.01 )          ! *** ANGULAR WAVE NUMBER (RAD/M)
      ELSE
        WV(L).DIR    = 0.
        WV(L).LENGTH = 0.
        WV(L).FREQ   = 1.
        WV(L).DISSIPA(KC) = 0.
        WV(L).HEISIG = 0.
        WV(L).K      = 1.
      ENDIF          
    ENDIF
  ENDDO
  
  WV(2:LA).HEIGHT = WV(2:LA).HEISIG  ! *** ASSIGN INCIDENT SIGNIFICANT WAVE HEIGHT TO APPLIED WAVE HEIGHT
  
END SUBROUTINE 
  
