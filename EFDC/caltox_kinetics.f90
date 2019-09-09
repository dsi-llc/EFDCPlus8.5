SUBROUTINE CALTOX_KINETICS

  ! **  SUBROUTINE CALTOX_KINETICS CALCULATES TOXICS DECAY FOR THE WATER COLUMN AND
  ! **  THE SEDIMENT BED AND IS CALLED FROM SSEDTOX
  !
  !--------------------------------------------------------------------------------
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !--------------------------------------------------------------------------------
  ! 2012-12-05        PAUL M. CRAIG    RESTRUCTURED AND ADDED OMP
  ! 2017-09-01        PAUL M. CRAIG    REWRITTEN TO UPDATE FRAMEWORK TO ADD
  !                                    MULTIPLE KINETIC PROCESSES.
  !                                    ADDED BIODEGRADATAION AND VOLATILIZATION 
  !********************************************************************************

  ! *** ITOXKIN(1,NT) = BULK DECAY,         0-DO NOT USE, 1-USE
  ! *** ITOXKIN(2,NT) = BIODEGRADATION,     0-DO NOT USE, 1-USE
  ! *** ITOXKIN(3,NT) = VOLATILIZATION,     0-DO NOT USE, 1-SIMPLE, 2-COMPUTED
  ! *** ITOXKIN(4,NT) = PHOTOLYSIS,         0-DO NOT USE, 1-SIMPLE  (NOT IMPLEMENTED)
  ! *** ITOXKIN(5,NT) = HYDROLYSIS,         0-DO NOT USE, 1-SIMPLE  (NOT IMPLEMENTED)
  ! *** ITOXKIN(6,NT) = DAUGHTER PRODUCTS,  0-DO NOT USE, 1-SIMPLE  (NOT IMPLEMENTED)
  
  USE GLOBAL
  USE CALCSERMOD,ONLY:LIN_INTER_COEF
  
  IMPLICIT NONE
  
  INTEGER :: K,L,NT,ND,LF,LL,LP,IDECAYB,IDECAYW,KBOT,KINC,NS,M1,M2
  REAL :: COEFF,DENA,DENW,DA,DW,CDRAG,HE,SCA,SCW,TKA,TKW,TOXMW,VISCA,VISCW,USTR
  REAL :: U_AVG, V_AVG, VEL_AVG, K_L, K_G, KV, HPU, TOXDIS, VOLTERM, DEPTH, WTM1, WTM2, TIME
  REAL,SAVE :: TOXTIME, RA, RCONST
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: CDECAYB
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: CDECAYW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)     :: TLAYER
  !REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:) :: KINTERM                        ! *** TO BE USED FOR DAUGHTER PRODUCTS
  
  REAL,EXTERNAL :: FUNDEN

  IF(  .NOT. ALLOCATED(CDECAYB) )THEN
    ALLOCATE(CDECAYB(LCM,KBM))  
    ALLOCATE(CDECAYW(LCM,KCM))  
    ALLOCATE(TLAYER(KCM))
    !ALLOCATE(KINTERM(LCM,KCM,NTXM))
  
    CDECAYB=0.0
    CDECAYW=0.0
    TLAYER=0.0
    
    TOXTIME = 0.0
    RCONST = 8.206E-5                                                ! *** Universal gas constant (atm-m3/mole K)
    CALL VISC_TABLE(1,20.,0.,VISCW)                                  ! *** Initialize viscosity table (kg/m/s)
    TOXSTEPW = TOXSTEPW - DELT/10.
    
    IF( ITOXTEMP == 1 )THEN
      ! *** SET CONSTANT TEMPERATURE
      DO L=2,LA
        DO K=1,KC
          TEM(L,K) = TOXTEMP
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  
  TOXTIME = TOXTIME + DTSED
  IF( TOXTIME < TOXSTEPW ) RETURN
  
  IF( ITOXTEMP > 1 )THEN
    ! *** ASSIGN TIME VARYING TEMPERATURES
    NS = ITOXTEMP - 1
    TIME = TIMEDAY
    CALL LIN_INTER_COEF(TIME,TSTEM(NS).TIM,NS,2,M1,M2,WTM1,WTM2)
    DO K=1,KC
      TLAYER(K) = WTM1*TSTEM(NS).VAL(M1,K) + WTM2*TSTEM(NS).VAL(M2,K)
    ENDDO      
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LP,L,K) 
    DO ND=1,NDM  
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          TEM(L,K) = TLAYER(K)
        ENDDO
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF
    
  ! *** SET LAYER ORIENTATION
  IF( LSEDZLJ )THEN
    KINC  = 1
    KBOT  = KB
  ELSE
    KINC  = -1
    KBOT  = 1
  ENDIF

  !KINTERM = 0.0    
  !RCONST = 8.3144
  IDECAYW = 0
  IDECAYB = 0

  ! *** LOOP OVER EACH TOXIC AND APPLY LOSS TERMS
  DO NT=1,NTOX
    ! *** INITIALIZE TOXIC CONSTANTS
    IF( TOX_MW(NT) > 0. ) TOXMW = 1./TOX_MW(NT)**0.66667

    !$OMP PARALLEL DEFAULT(SHARED)
  
    IF( (ITOXKIN(1,NT)+ITOXKIN(2,NT) ) > 0 )THEN
      ! *** ZERO LOSS COEFFICIENTS
      !$OMP DO PRIVATE(ND,LP,L,K) 
      DO ND=1,NDM  
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            CDECAYW(L,K) = 0.0
          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO
    
      ! *** ZERO LOSS COEFFICIENTS
      !$OMP DO PRIVATE(ND,LL,LF,L,K) 
      DO ND=1,NDM  
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          DO K=KBOT,KBT(L),-KINC
            CDECAYB(L,K) = 0.0
          ENDDO
        ENDDO
      ENDDO   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    ENDIF
    
    ! *** BULK DECAY
    IF( ITOXKIN(1,NT) > 0 )THEN
      ! *** COMPUTE DECAY COEFFICIENTS FOR EVERY TIMESTEP IF USING DYNAMIC TIME STEPPING
  
      ! *** BULK DECAY COEFFICIENT: WATER
      IF( TOX_BLK_KW(NT) > 0. )THEN
        IDECAYW = 1
        !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,COEFF) 
        DO ND=1,NDM  
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              CDECAYW(L,K) = TOX_BLK_KW(NT)
            ENDDO
          ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      ENDIF
      
      ! *** BULK DECAY COEFFICIENT: SEDIMENT BED
      IF( TOX_BLK_KB(NT) > 0. )THEN
        IDECAYB = 1
        !$OMP DO PRIVATE(ND,LF,LL,L,K,COEFF,DEPTH) 
        DO ND=1,NDM  
          LF=2+(ND-1)*LDM
          LL=MIN(LF+LDM-1,LA)
          DO L=LF,LL
            DEPTH = 0.0
            DO K=KBT(L),KBOT,KINC
              DEPTH = DEPTH + HBED(L,K)
              IF( DEPTH > TOX_BLK_MXD(NT) )CYCLE
              CDECAYB(L,K) = TOX_BLK_KB(NT)
            ENDDO
          ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      ENDIF
    ENDIF
  
    ! *** BIODEGRADATION
    IF( ITOXKIN(2,NT) > 0 .AND. ( ISTRAN(2) > 0 .OR. ITOXTEMP > 0 ) )THEN
      ! *** DECAY COEFFICIENT: WATER
      IF( TOX_BIO_KW(NT) > 0. )THEN
        IDECAYW = 1
        !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,COEFF) 
        DO ND=1,NDM  
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              COEFF = TOX_BIO_KW(NT)*TOX_BIO_Q10W(NT)**( 0.1*(TEM(L,K) - TOX_BIO_TW(NT)) )
              CDECAYW(L,K) = CDECAYW(L,K) + COEFF
            ENDDO
          ENDDO
        ENDDO   ! *** END OF DOMAIN LOOP
        !$OMP END DO
      ENDIF
      
      ! *** DECAY COEFFICIENT: SEDIMENT BED
      IF(  TOX_BIO_KB(NT) > 0. )THEN
        IDECAYB = 1
        IF( (HTBED1 + HTBED2) > 0. .AND. ITOXTEMP == 0 )THEN
          ! *** BED TEMPERATURES ARE SIMULATED, SO USE
          !$OMP DO PRIVATE(ND,LF,LL,L,K,COEFF,DEPTH) 
          DO ND=1,NDM  
            LF=2+(ND-1)*LDM
            LL=MIN(LF+LDM-1,LA)
            DO L=LF,LL
              DEPTH = 0.0
              DO K=KBT(L),KBOT,KINC
                DEPTH = DEPTH + HBED(L,K)
                IF( DEPTH > TOX_BIO_MXD(NT) )CYCLE
                COEFF = TOX_BIO_KB(NT)*TOX_BIO_Q10B(NT)**( 0.1*(TEMB(L) - TOX_BIO_TB(NT)) )
                CDECAYB(L,K) = CDECAYB(L,K) + COEFF
              ENDDO
            ENDDO
          ENDDO   ! *** END OF DOMAIN LOOP
          !$OMP END DO
        ELSE
          ! *** BED TEMPERATURES ARE NOT SIMULATED, USE BOTTOM LAYER WATER TEMPERATURES
          !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,COEFF,DEPTH) 
          DO ND=1,NDM  
            LF=2+(ND-1)*LDM
            LL=MIN(LF+LDM-1,LA)
            DO L=LF,LL
              DEPTH = 0.0
              DO K=KBT(L),KBOT,KINC
                DEPTH = DEPTH + HBED(L,K)
                IF( DEPTH > TOX_BIO_MXD(NT) )CYCLE
                COEFF = TOX_BIO_KB(NT)*TOX_BIO_Q10B(NT)**( 0.1*(TEM(L,KSZ(L)) - TOX_BIO_TB(NT)) )
                CDECAYB(L,K) = CDECAYB(L,K) + COEFF
              ENDDO
            ENDDO
          ENDDO   ! *** END OF DOMAIN LOOP
          !$OMP END DO
        ENDIF
      ENDIF  ! *** END OF BED DEGRADATION BLOCK  
    ENDIF    ! *** END OF BIODEGRADATION BLOCK

    IF( IDECAYW == 1 )THEN
      ! *** APPLY DEGRADATION TO THE WATER COLUMN
      !$OMP DO PRIVATE(ND,LP,L,K) 
      DO ND=1,NDM  
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            TOX(L,K,NT) = TOX(L,K,NT)*(1.0 - TOXTIME*CDECAYW(L,K))
          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF

    IF( IDECAYB == 1 )THEN
      ! *** APPLY DEGRADATION TO THE SEDIMENT BED
      !$OMP DO PRIVATE(ND,LL,LF,L,K) 
      DO ND=1,NDM  
        LF=2+(ND-1)*LDM
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL
          DO K=KBOT,KBT(L),-KINC
            TOXB(L,K,NT) = TOXB(L,K,NT)*(1.0 - TOXTIME*CDECAYB(L,K))
          ENDDO
        ENDDO
      ENDDO   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    ENDIF
    
    IF( ITOXKIN(3,NT) > 0 .AND. ( ISTRAN(2) > 0 .OR. ITOXTEMP > 0 ) )THEN
      ! ***    VOLATILIZATION
      ! ***      TOX_MW(NT)=MOLECULAR WEIGHT
  
      !! *** UPDATE DRAG COEFF, USTAR
      !IF( NWSER > 1 .OR. ISHELTERVARY > 0 )THEN      
      !  !$OMP DO PRIVATE(ND,LF,LL,L,UWIND,CDRAG) 
      !  DO ND=1,NDM  
      !    LF=2+(ND-1)*LDM
      !    LL=MIN(LF+LDM-1,LA)
      !    DO L=LF,LL
      !      UWIND = WINDST(L)
      !      CDRAG = (6.1+0.63*UWIND)*1.0E-04
      !      USTR(L) = UWIND*CDRAG**0.5                                   ! m/s
      !    ENDDO
      !  ENDDO   ! *** END OF DOMAIN LOOP
      !  !$OMP END DO
      !ELSEIF( NWSER == 1 )THEN
      !  CDRAG = (6.1+0.63*WINDST(2))*1.0E-04
      !  UWIND = WINDST(2)*CDRAG**0.5
      !  !$OMP DO PRIVATE(ND,LF,LL,L) 
      !  DO ND=1,NDM  
      !    LF=2+(ND-1)*LDM
      !    LL=MIN(LF+LDM-1,LA)
      !    DO L=LF,LL
      !      USTR(L) = UWIND                                             ! m/s
      !    ENDDO
      !  ENDDO   ! *** END OF DOMAIN LOOP
      !  !$OMP END DO
      !ENDIF
      
      ! *** UPDATE MOLECULAR DIFFUSIVITY,SCHMIDT NUMBER, HENRY'S CONST (SURFACE ONLY)
      !$OMP DO PRIVATE(ND,K,LP,L,U_AVG,V_AVG,VEL_AVG,HPU,TKW,VISCW,DENW,DW,SCW,TKA,RA,VISCA,DENA,DA,SCA,CDRAG,USTR,K_L,K_G,HE,KV,TOXDIS,VOLTERM) 
      DO ND=1,NDM  
        DO LP=1,LLWET(KC,ND)
          L = LKWET(LP,KC,ND) 
            
          ! *** CALCULATE AVERAGE VELOCITY
          U_AVG = STCUV(L)*UHE(L)/HU(L)
          V_AVG = STCUV(L)*VHE(L)/HV(L)
          VEL_AVG = SQRT( U_AVG*U_AVG + V_AVG*V_AVG )
          
          TKW = TEM(L,KC) + 273.15                                       ! *** Temperature (Kelvin)
          IF( VEL_AVG > TOX_VEL_MAX .OR. HP(L) < TOX_DEP_MAX )THEN
            ! *** RIVER/STREAM CONDITIONS: WATER TURBULENCE CONTROLED GAS TRANSFER
            HPU = 13.584*VEL_AVG**2.9135
            
            IF( HP(L) < 0.61 )THEN
              ! *** OWENS FORMULA FOR K_L (K_A)
              K_L = 5.349*VEL_AVG**0.667 * HPI(L)**1.85                  ! *** Compute liquid transfer coefficient (m/day)
              K_L = K_L/86400.                                           ! *** Convert from m/day to m/s
            ELSEIF( VEL_AVG < 0.518 .OR. HP(L) > HPU )THEN
              ! *** O'CONNOR-DOBBINS
              DW = 22.E-9*TOXMW                                          ! *** Diffusivity of toxic in water (m^2/s)
              K_L = DW*VEL_AVG**0.5 * HPI(L)**1.5                        ! *** Compute liquid transfer coefficient (m/s)
            ELSE
              ! *** CHURCHILL
              K_L = 5.049*VEL_AVG**0.969 * HPI(L)**1.673                 ! *** Compute liquid transfer coefficient (m/day)
              K_L = K_L/86400.                                           ! *** Convert from m/day to m/s
            ENDIF
            K_G = 0.0011574                                              ! *** Set gas transfer coefficient (m/s)  [100 m/day]
          ELSE
            ! *** LAKE/QUIESCENT CONDITIONS: AIR TURBULENCE CONTROLED GAS TRANSFER

            ! *** SCHMIDT NUMBER FOR WATER (dimensionless)
            CALL VISC_TABLE(0,TEM(L,KC),SAL(L,KC),VISCW)                 ! *** Viscosity of water (kg/m/s)
            DENW = FUNDEN(SAL(L,KC),0.0,TEM(L,KC))                       ! *** kg/m^3
              
            DW = 22.E-9*TOXMW                                            ! *** Diffusivity of toxic in water (m^2/s)
            SCW = VISCW/DENW/DW                                          ! *** Schmidt Number 
            
            ! *** SCHMIDT NUMBER FOR AIR (dimensionless)
            TKA = TATMT(L) + 273.15                                      ! *** Temperature (Kelvin)
            RA = 287.05                                                  ! *** Specific Gas Constant (m^2/s^2/K) (J/kg/K)  (J=kg*m^2/s^2)
            DENA = PATMT(L)*100./RA/TKA                                  ! *** Density of air (kg/m^3)

            VISCA = 1.458E-6*TKA**1.5/(TKA+110.4)                        ! *** Viscosity of air (kg/m/s)
            DA = 1.9E-4*TOXMW                                            ! *** Diffusivity of toxic in air (m^2/s)
            SCA = VISCA/DENA/DA                                          ! *** Schmidt Number 
            
            ! *** COMPUTE GAS AND LIQUID TRANSFER COEFFICIENTS (M/S)
            CDRAG = ( 6.1 + 0.63*WINDST(L) )*1.0E-04
            CDRAG = AMAX1(CDRAG,0.0011)
            USTR = WINDST(L)*CDRAG**0.5                                  ! *** Shear velocity (m/s)
            IF( ITOXKIN(3,NT) == 2 )THEN
              ! *** CALCULATE LIQUID FILM COEFF AS PER MACKAY AND YEUN (1983) (WIND SPEED DRIVEN)
              IF( USTR < 0.3 )THEN
                K_L = 0.0144*USTR**2.2/SQRT(SCW) + 1E-6
              ELSE
                K_L = 0.0341*USTR/SQRT(SCW) + 1E-6
              END IF
              K_G = 0.0462*USTR/SCA**0.6667                              ! *** .0462 = k^0.33/LAMDA
            ELSE              
              K_L = 0.0462*USTR*(DENA/DENW)**.5/SCW**0.6667              ! *** .0462 = k^0.33/LAMDA
              K_G = 0.0462*USTR/SCA**0.6667                              ! *** .0462 = k^0.33/LAMDA
            ENDIF
            K_G = AMAX1(K_G,1.0E-06)                                     ! *** Gas mass transfer (m/s)
          ENDIF
          K_L = AMAX1(K_L,1.0E-9)                                        ! *** Liquid mass transfer (m/s)
            
          HE = TOX_HE(NT)/RCONST/TKW                                     ! *** Henry's term (dimensionless)
          KV = 1.0 / ( 1.0/K_L + 1.0/(K_G*HE) )                          ! *** Mass transfer rate at 20 degC (m/s)
          KV = KV*TOX_KV_TCOEFF(NT)**(TEM(L,KC)-20.)                     ! *** Mass transfer rate at water temperature (m/s)
  
          ! *** DISSOLVED FRACTION
          TOXDIS = 1./(1.+TOXPFTW(L,KC,NT))
          TOXDIS = TOX(L,KC,NT)*TOXDIS
            
          ! *** COMPUTE VOLATIZATION FLUX TERM (ug-m/l-s or mg/m^2-s)
          !VOLTERM = KV*HPKI(L,KC)*( TOXDIS - TOXATM(L,NT)/HE )
          VOLTERM = TOX_ADJ(NT)*KV*HPKI(L,KC)*( TOXDIS - TOX_ATM(NT)/HE )
          TOX(L,KC,NT) = TOX(L,KC,NT) - VOLTERM*TOXTIME
          
        ENDDO
      ENDDO   ! *** END OF DOMAIN LOOP
      !$OMP END DO
    ENDIF     ! *** END OF VOLATILIZATION
  
    IF( ITOXKIN(4,NT) > 0 )THEN
      ! ***    PHOTOLOSIS
      ! ***      RKTOXP(NT)=BASE RATE
      ! ***      SKTOXP(NT)=SOLAR RADIATION AT BASE RATE
    ENDIF
  
    !$OMP END PARALLEL
  
  ENDDO   ! *** END OF NTOX LOOP
    
  TOXTIME = 0.0
  
  RETURN

  END

  !****************************************************************
  !  SUBROUTINE FOR CALCULATING WATER VISCOSITY FROM TEMPERATURE
  SUBROUTINE VISC_TABLE(INIT_FLAG,TEMP,SAL,VISC)

    ! *** THIS SUBROUTINE CALCULATES AN INTERPOLATED VALUE OF WATER
    ! *** VISCOSITY

    ! *** TEMP=TEMPERATURE AT WHICH TO CALCULATE WATER VISCOSITY
    ! *** VISC=CALCULATED WATER VISCOSITY (cp)

    INTEGER, INTENT(IN) :: INIT_FLAG
    REAL, INTENT(IN)  :: TEMP,SAL
    REAL, INTENT(OUT) :: VISC
    INTEGER :: NX
    REAL,SAVE, DIMENSION(11) :: TEMP_VALS,VISC_VALS

    ! *** ONLY DEFINE REFERENCE VALUES FIRST TIME THE SUBROUTINE IS CALLED
    IF(INIT_FLAG == 1 )THEN
      ! *** REFERENCE: CRC HANDBOOK OF PHYSICS AND CHEMISTRY, 1995 (STUDENT EDITION)
      TEMP_VALS = [   0.,  10.,  20.,  30.,  40., 50.,   60.,  70.,  80.,  90., 100.]    ! *** deg.C.
      VISC_VALS = [1793.,1307.,1002.,797.7,653.2,547.0,466.5,404.0,354.4,314.5,281.8]    ! *** Pa*s*1E6  (kg/m/s)
    
      ! *** CONVERT TO centipoise (10^-2 g/cm*s)
      !VISC_VALS = VISC_VALS/1000.
    END IF

    ! *** LINEARLY INTERPOLATE BETWEEN TABLE VALUES
    IF( TEMP < 0. )THEN
      ! *** CALHEAT HANDLES THIS CASE, SO JUST APPY MINIMUM
      VISC = VISC_VALS(1)
    ELSEIF( TEMP >= 100. )THEN
      WRITE(220,*) 'WATER TEMPERATURE ABOVE 100 deg.C. ENCOUNTERED!'
      WRITE(220,*) 'SUBROUTINE VISC_TABLE STOPPED SIMULATION'
      STOP
    ELSE
      DO NX=1,9
        IF( (TEMP >= TEMP_VALS(NX)) .AND. (TEMP < TEMP_VALS(NX+1)) )THEN
          VISC = VISC_VALS(NX) + (TEMP-TEMP_VALS(NX))*( VISC_VALS(NX+1)-VISC_VALS(NX) )/( TEMP_VALS(NX+1)-TEMP_VALS(NX) )
        ENDIF
      ENDDO
    ENDIF
!
END SUBROUTINE VISC_TABLE
