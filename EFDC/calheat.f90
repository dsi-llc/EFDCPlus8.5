MODULE HEAT_MODULE

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2015-01       Paul M. Craig     Added fully coupled Ice Sub-model
  !                  Dang H Chung        
  ! 2014-05-08       Paul M. Craig     Converted to Module and fixed OMP issue with equilibrium temperature

  USE GLOBAL
  USE INFOMOD,ONLY:READSTR
  
  IMPLICIT NONE

  ! ******* PARAMETER DECLARATIONS
  REAL, PARAMETER, PRIVATE :: MPS_TO_MPH          = 2.23714
  REAL, PARAMETER, PRIVATE :: W_M2_TO_BTU_FT2_DAY = 7.60796
  REAL, PARAMETER, PRIVATE :: FLUX_BR_TO_FLUX_SI  = 0.23659
  REAL, PARAMETER, PRIVATE :: BTU_FT2_DAY_TO_W_M2 = 0.1314
  REAL, PARAMETER, PRIVATE :: BCONV               = 1.520411

  REAL,PARAMETER,PRIVATE   :: LHF     = 333507.0    ! *** LATENT HEAT OF FUSION FOR ICE (J/KG)
  REAL,PARAMETER,PRIVATE   :: CP      = 4179.0      ! *** SPECIFIC HEAT (J/KG/degC)
  
  REAL,SAVE,PRIVATE        :: REFICE
  REAL,SAVE,PRIVATE        :: RHOWCP
  REAL,SAVE,PRIVATE        :: RHOWCPI
  REAL,SAVE,PRIVATE        :: RHOICP
  REAL,SAVE,PRIVATE        :: RHOICPI
  REAL,SAVE,PRIVATE        :: CREFLI
  REAL,SAVE,PRIVATE        :: RHOILHFI

CONTAINS

  SUBROUTINE CALHEAT(ISTL_)  

    !   Subroutine CALHEAT takes the information from the atmospheric boundary
    !   file and the wind forcing file and calculates the net heat flux across
    !   the water surface boundary.  The heat flux is then used to update the
    !   water temperature either in the surface cells, or distributed across
    !   the cells in the vertical and into the bottom.  The subroutine has
    !   four options.  These are:
    !
    !    ISOPT(2)=1: Full surface and internal heat transfer calculation
    !                using meteorologic data from input stream.
    !                IASWRAD=0:  ADSORB SW SOLR RAD TO ALL LAYERS AND BED
    !                IASWRAD=1:  ADSORB SW SOLR RAD TO TO SURFACE LAYER
    !    ISOPT(2)=2: Transient equilibrium surface heat transfer calculation
    !                using external equilibrium temperature and heat transfer
    !                coefficient data from the meteorologic input data.
    !    ISOPT(2)=3: Equilibrium surface heat transfer calculation using constant
    !                equilibrium temperature and heat transfer coefficients
    !                HEQT and TEMO read in through input stream.
    !    ISOPT(2)=4: Equilibrium surface heat transfer calculation using algorithm
    !                from CE-QUAL-W2.
    !
    !   The heat flux terms are derived from a paper by Rosati
    !   and Miyakoda (1988) entitled "A General Circulation Model for Upper Ocean
    !   Simulation".  The heat flux is prescribed by term for the following
    !   influxes and outfluxes:
    !
    !     - Short Wave Incoming Radiation (+)
    !     - Net Long Wave Radiation (+/-)
    !     - Sensible Heat Flux (convection -)
    !     - Latent Heat Flux (evaporation +/-)
    !
    !   Two formulations of the Latent Heat Flux are provided.  The first is from
    !   the Rosati and Miyakoda paper, the second is an alternate formulation by
    !   Brady, Graves, and Geyer (1969).  The second formulation was taken from
    !   "Hydrodynamics and Transport for Water Quality Modeling" (Martin and
    !   McCutcheon, 1999).  The Rosati and Miyakoda formulation will have zero
    !   evaporative cooling or heating if wind speed goes to zero.  The Brady,
    !   Graves, and Geyer formulation provides for a minimum evaporative cooling
    !   under zero wind speed.
    !
    !
    ! VARIABLE LIST:
    !
    !   CLOUDT  = Cloud Cover (0 to 10)
    !   HCON    = Sensible Heat Flux (W/m2)
    !   HLAT    = Latent Heat Flux (W/m2)
    !   HLWBR   = Net Longwave Radiation (Atmospheric Long Wave plus Back Radiation, W/m2)
    !   SOLSWRT = Total Short Wave Incoming Radiation (W/m2)
    !   SVPW    = Saturation Vapor Pressure in mb Based on the Water Surface Temperature (mb)
    !   TATMT   = Temperature of Air Above Water Surface (deg C)
    !   TEM     = Water Temperature in Cell (deg C)
    !   VPAT    = Vapor Pressure of Air at Near Surface Air Temperature (mb)
    !   WINDST  = Wind Speed at 2 Meters Over Cell Surface (m/s)
    !
    ! MODIFICATION HISTORY:
    !
    !   Date       Author             Comments*6
    !   ---------- ------------------ -----------------------------------------------------------------
    !   2015-01    Paul M. Craig      Added Fully Heat Coupled Ice Sub-model
    !              Dang Chung        
    !   2014-12    Paul M. Craig      Added the new Evaporation Approach
    !   2012-10    Paul M. Craig      Added DOC to the light extinction and 
    !                                 standardized for all other routines
    !   2011-03    Paul M. Craig      Converted to F90, added OMP
    !   11/01/2005 Paul M. Craig      Added Option 4, to use the Equilibrium Temperature 
    !                                 algorithym from CE-QUAL-W2.  Also added the sub-option
    !                                 under this option to couple or decouple the bottom temperature
    !                                 to the water column temperatures.
    !                                 Added the ability to input spatially variable bed temps and
    !                                 thermally active bed thicknesses.
    !                                 Also cleaned up the code and added structure.
    !   11/07/2000 Steven Peene       Cleaned code, provided more detailed
    !                                 descriptions, added alternate formulation
    !                                 for latent heat flux, separated out
    !                                 individual heat flux terms
    !   06/01/1992 John M. Hamrick    Orignial author

    USE INFOMOD,ONLY:SKIPCOM

    INTEGER, INTENT(IN) :: ISTL_ 
    INTEGER :: L, K, I, J, L1, ND, LF, LL, LP, IOS
    
    INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ICOUNT
    
    REAL       :: T1, T2, SVPW1, DTHEQT, WCKESS, TFLUX
    REAL       :: TFAST, TFAST1, TSLOW, TSLOW1, RSN, C1, C2, UBED, VBED
    REAL       :: USPD, TMPKC, SRO, SRON
    REAL       :: THICK, TMIN
    REAL       :: ET, CSHE, TEMP, RICETHKL0, ET0, CSHE0, WSHLTR0, PSHADE0
    REAL       :: HBLW,HBCD,HBCV,HBEV,HBNT
    REAL,SAVE  :: TIMEDRY
    
    REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TBEDTHK
    REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: FLUXTB
    REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: PSHADE_OLD
    REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WSHLTR_OLD
    REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TMIND

    CHARACTER(200) :: STR
    
    IF(  .NOT. ALLOCATED(TBEDTHK) )THEN
      ALLOCATE(TBEDTHK(LCM))
      ALLOCATE(FLUXTB(LCM))
      ALLOCATE(PSHADE_OLD(NDM))
      ALLOCATE(WSHLTR_OLD(NDM))
      ALLOCATE(ICOUNT(NDM))
      ALLOCATE(TMIND(NDM))
        
      FLUXTB     = 0.0
      PSHADE_OLD = 0.0
      WSHLTR_OLD = 1.0
      
      ! *** Ininitialze Heat Exchange Terms
      WRITE(*,'(A)')'CALHEAT: INITIALIZING'

      RPEMMINT = 0.
      RADNET  = 0.
      TBEDTHK = ABS(DABEDT)
      TEMB    = ABS(TBEDIT)
      TIMEDRY = TBEGIN + 0.25
      
      IF( DABEDT < 0. .AND. .NOT. ( ISRESTI > 0 .AND. ISCI(2) == 1 ) ) THEN
        ! *** READ IN THE SPATIALLY VARYING INIT T AND BED THICKNESS (DABEDT) 
        WRITE(*,'(A)')' CALHEAT: READ IN THE SPATIALLY VARYING INIT T AND BED THICKNESS: TEMB.INP'
        OPEN(1001,FILE='temb.inp',ERR=1000,STATUS='OLD')  
        STR = READSTR(1001)
        DO L1=2,LA
          READ(1001,*,END=1000) I,J,T1,T2
          L=LIJ(I,J)
          TEMB(L)    = T1
          TBEDTHK(L) = T2
        ENDDO  
        1000 CLOSE(1001)  
      ENDIF

      PSHADE = 1.0          ! *** DEFAULT SHADE FACTOR

      IF( ISTOPT(2) == 1 .OR. ISTOPT(2) >= 4 .OR. ISICE > 2 )THEN 
        IF( USESHADE )THEN
          ! *** READ IN SPATIALLY VARYING SHADING FACTORS
          WRITE(*,'(A)')' CALHEAT: READ IN SPATIALLY VARYING SHADE: PSHADE.INP'
          OPEN(1001,FILE='pshade.inp',ERR=1010,STATUS='OLD')
 
          STR = READSTR(1001)

          DO L1 = 2, LA
            READ(1001,*,END=1010) I,J,T1
            L=LIJ(I,J)
            PSHADE(L)=T1
          ENDDO  
          1010 CLOSE(1001)  

        ELSE
          WRITE(*,'(A)')' CALHEAT: SETTING CONSTANT SHADE TO: 1.0 (NO SHADE)'
        ENDIF
        
        ! ** READING ICEINIT.INP FILE FOR THE INITIAL ICE THICKNESS
        IF( ISICE > 2 )THEN
          IF( ISRESTI == 0 )THEN
            WRITE(*,'(A)')' CALHEAT: READ IN SPATIALLY VARYING ICE CONDITIONS: ICE.INP'
            OPEN(1001,FILE='ice.inp',STATUS='OLD')      
            CALL SKIPCOM(1001,'*')
            DO LL=2,LA
              READ(1001,*,IOSTAT=IOS) I,J,RICETHKL0
              IF( IOS > 0 ) STOP 'ICE.INP: READING ERROR'
              L = LIJ(I,J)      
              ICETHICK(L) = RICETHKL0        ! INITIAL ICE THICKNESS IN M
            ENDDO
            CLOSE(1001)
          ENDIF
          
          DO L=2,LA
            IF( ICETHICK(L) >= MINICETHICK )THEN
              ICECELL(L) = .TRUE.
              ICECOVER(L) = 1.0            ! RATIO OF ICE COVER: 0/1
              LFRAZIL = .TRUE.
            ELSE
              ICECELL(L) = .FALSE.
              ICEVOL(L) = 0.0
            ENDIF       
          ENDDO
        ENDIF
        
        ICOUNT  = 1  ! *** FORCE INITIAL ZEROING
      ENDIF

      ! *** INITIALIZE EXTINCTION COEFFICIENTS
      DO K=1,KC
        DO L=2,LA
          WCKESS = GET_EXT_COEFF(L,K)
        ENDDO
      ENDDO
      RHOWCP  = 1000.0*CP             ! *** CP-SPECIFIC HEAT (J/KG/degC) * RHO (KG/M3)
      RHOWCPI = 1./RHOWCP

      IF( DEBUG )THEN
        WRITE(*,'(A,F8.1)')' CALHEAT: Bed Temp(L=2):', TEMB(2)
        OPEN(77,FILE=OUTDIR//'calheat.dia',status='unknown')
        CLOSE(77,status='DELETE')
        OPEN(77,FILE=OUTDIR//'calheat.dia',status='NEW')
        WRITE(77,998)'TIMEDAY','SRON','ET','TD_C','TA_C','TDEW_F','TAIR_F','FW'
        998 FORMAT(A11,8A9)
      ENDIF
    
    ENDIF    ! *** END OF CALHEAT INITIALIZATION BLOCK

    IF( ISTL_ == 2 )THEN  
      IF( ISDYNSTP == 0 )THEN  
        DELT=DT  
      ELSE  
        DELT=DTDYN  
      END IF  
    ELSE
      DELT=DT2 
    ENDIF  
    TMIND = 9999.

    ! *** OVERWRITE THE INPUT SOLAR RAD WITH A COMPUTED ONE
    IF( COMPUTESOLRAD )THEN
      CALL SHORT_WAVE_RADIATION(CLOUDT(2), SRO, SRON)
      IF( SRON > 1.0E-4 )THEN
        LDAYLIGHT = .TRUE.
      ELSE
        LDAYLIGHT = .FALSE.
      ENDIF
    ENDIF

    ! *** PREPARE SOLAR RADIATION AND SURFACE PROPERTIES BEFORE 
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, LP, L, K)                       
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      IF( COMPUTESOLRAD )THEN
        ! *** OVERWRITE THE INPUT SOLAR RAD WITH A COMPUTED ONE
        DO LP=LF,LL  
          L=LWET(LP)  
          SOLSWRT(L)=SRON
        ENDDO
      ENDIF

      IF( LDAYLIGHT .OR. IWQSUN /= 2 )THEN
        IF( ISTRAN(8) > 0 .AND. (ISTOPT(2) >= 4 .OR. WQKECHL < 0.0) )THEN
          ! *** If using WQ then use the WQ Coefficients for light extinction
    
          ! *** Water Quality is Active so account for Chlorophyll
          ! *** Compute WQCHL (Chlorophyll) Using Algal Biomass & factors
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              WQCHL(L,K) = WQV(L,K,1)*WQCHLC + WQV(L,K,2)*WQCHLD + WQV(L,K,3)*WQCHLG  
            ENDDO
          ENDDO
        ENDIF

        IF( USESHADE )THEN
          DO LP=LF,LL  
            L=LWET(LP)  
            ! *** APPLY PSHADE FACTORS
            SOLSWRT(L) = SOLSWRT(L)*PSHADE(L)
          ENDDO
        ENDIF

        ! *** SET SOLAR RADIATION AT SURFACE AND FOR LAYERS
        IF( ISTOPT(2) == 1 .OR. ISTOPT(2) >= 4 .OR. (ISTRAN(8) > 0 .AND. WQKECHL < 0.0 ) )THEN 
          ICOUNT(ND) = ICOUNT(ND) + 1
          DO LP=LF,LL  
            L=LWET(LP)  
            CALL SET_LIGHT(L)
          ENDDO
        ELSEIF( ISTRAN(8) > 0 )THEN
          ICOUNT(ND) = ICOUNT(ND) + 1
          DO LP=LF,LL  
            L=LWET(LP)  
            RADTOP(L,KC) = SOLSWRT(L)
          ENDDO
        ENDIF

      ELSEIF( ICOUNT(ND) > 0 )THEN
        ! *** NO LIGHT CASE (ONLY ZERO IF LAST TIMESTEP WAS DAYLIGHT)
        ICOUNT(ND) = 0
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            RADBOT(L,K)=0.
            RADTOP(L,K)=0.
            RADNET(L,K)=0.
          ENDDO
        ENDDO
        DO LP=LF,LL  
          L=LWET(LP)  
          RADTOP(L,0)=0.
        ENDDO
      ELSE
        ! *** ALWAYS ZERO TOP NET RADIATION 
        K = KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          RADBOT(L,K)=0.
          RADTOP(L,K)=0.
          RADNET(L,K)=0.
        ENDDO
      ENDIF   ! *** END OF LDAYLIGHT BLOCK

    ENDDO   ! *** END OF DOMAIN
    !$OMP END PARALLEL DO

    ! *** CONDUCT HEAT TRANSFER AT SURFACE, BOTTOM AND SOLAR RADIATION EXTINCTION, BASED ON SELECTED METHOD
    IF( ISTOPT(2) == 1 )THEN  
      ! *** FULL HEAT BALANCE
      
      !$OMP PARALLEL DEFAULT(SHARED) 
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,C2,SVPW1) 
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)

        ! *** SURFACE HEAT FLUX WITHOUT SOLAR RADIATION (m * degC)
        DO LP=LF,LL 
          L=LWET(LP)  
          C2=-DELT*DZIC(L,KC)
          IF( .NOT. ICECELL(L) )THEN
            ! *** CELL WITHOUT ICE
            SVPW1=SVPW(L)
            RADNET(L,KC) = C2*( 1.312E-14*((TEM(L,KC)+273.)**4)*(0.39-0.05*SQRT(VPAT(L)))*(1.-.8*CLOUDT(L))  &
                                +5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))                        &
                                +CCNHTT(L)*0.288E-3*WINDST(L)*(TEM(L,KC)-TATMT(L))                           &
                                +CLEVAP(L)*0.445*WINDST(L)*(SVPW1-VPAT(L))/PATMT(L) )
          ELSE
            RADNET(L,KC) = 0.
          ENDIF
        ENDDO
      ENDDO   ! *** END OF DOMAIN
      !$OMP END DO

      ! *** NET SHORTWAVE SOLAR RADIATION
      IF( IASWRAD == 0. )THEN  

        ! *** ADSORB SW SOLR RAD TO ALL LAYERS AND BED  
        IF( LDAYLIGHT )THEN

          !$OMP DO PRIVATE(ND,LF,LL,LP,L,TFAST,TFAST1,TSLOW,TSLOW1,C1,RSN)
          DO ND=1,NDM  
            LF=(ND-1)*LDMWET+1  
            LL=MIN(LF+LDMWET-1,LAWET)

            ! *** SURFACE LAYER
            DO LP=LF,LL  
              L=LWET(LP)  
              TFAST  = SWRATNF*(Z(L,KC)  -1.)
              TFAST1 = SWRATNF*(Z(L,KC-1)-1.)
              TSLOW  = SWRATNS*(Z(L,KC)  -1.)
              TSLOW1 = SWRATNS*(Z(L,KC-1)-1.)
              C1 = DELT*DZIC(L,KC)*RHOWCPI
              RSN = RADTOP(L,KC)*( FSWRATF   *EXP(TFAST*HP(L))     &
                                +(1.-FSWRATF)*EXP(TSLOW*HP(L))     &
                                - FSWRATF    *EXP(TFAST1*HP(L))    &
                                -(1.-FSWRATF)*EXP(TSLOW1*HP(L)) )
              RADNET(L,KC) = RADNET(L,KC) + RSN*C1

              ! *** SAVE THE SOLAR RADIATION
              RADBOT(L,KC) = RADTOP(L,KC) - RSN
            ENDDO
          ENDDO   ! *** END OF DOMAIN
          !$OMP END DO

          ! *** ALL REMAINING LAYERS
          IF( KC > 1 )THEN
            !$OMP DO PRIVATE(ND,K,LP,L,TFAST,TFAST1,TSLOW,TSLOW1,C1,RSN)
            DO ND=1,NDM  
              DO K=KS,1,-1
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  

                  TFAST  = SWRATNF*(Z(L,K)-1.)
                  TFAST1 = SWRATNF*(Z(L,K-1)-1.)
                  TSLOW  = SWRATNS*(Z(L,K)-1.)
                  TSLOW1 = SWRATNS*(Z(L,K-1)-1.)
                  C1 = DELT*DZIC(L,K)*0.2393E-6
                  RSN = RADTOP(L,KC)*( FSWRATF     *EXP(TFAST*HP(L))     &  
                                      +(1.-FSWRATF)*EXP(TSLOW*HP(L))     &
                                      - FSWRATF    *EXP(TFAST1*HP(L))    &
                                      -(1.-FSWRATF)*EXP(TSLOW1*HP(L)) )
                  RADNET(L,K) = RSN*C1

                  ! *** SAVE THE SOLAR RADIATION
                  RADTOP(L,K) = RADBOT(L,K+1)   
                  RADBOT(L,K) = RADTOP(L,K)-RSN
                ENDDO  
              ENDDO
            ENDDO   ! *** END OF DOMAIN
            !$OMP END DO
          ENDIF
          
          ! *** NOW FINALIZE THE TEMPERATURE 
          !$OMP DO PRIVATE(ND,K,LP,L)
          DO ND=1,NDM  
            DO K=1,KC  
              DO LP=1,LLWET(K,ND)
                L=LKWET(LP,K,ND)  
                TEM(L,K) = TEM(L,K) + HPI(L)*RADNET(L,K)  
              ENDDO  
            ENDDO
          ENDDO   ! *** END OF DOMAIN
          !$OMP END DO
       
        ELSE
          ! *** NOW FINALIZE THE TEMPERATURE FOR THE TOP LAYER
          K = KC
          !$OMP DO PRIVATE(ND,LP,L)
          DO ND=1,NDM  
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              TEM(L,K) = TEM(L,K) + HPI(L)*RADNET(L,K)  
            ENDDO  
          ENDDO
          !$OMP END DO
          
        ENDIF   ! *** END OF DAYLIGHT BLOCK

      ELSE   ! IF( IASWRAD == 1 )THEN  

        ! *** ADSORB SOLAR RAD TO TO SURFACE LAYER  
        !$OMP DO PRIVATE(ND,LF,LL,LP,L,C1,RSN,UBED,VBED,USPD)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)

          IF( LDAYLIGHT )THEN
            DO LP=LF,LL  
              L=LWET(LP)  
              C1 = DELT*DZIC(L,KC)*0.2385E-6
              RSN = RADTOP(L,KC)
              RADNET(L,KC) = RADNET(L,KC)+RSN*C1
              RADBOT(L,KC) = 0.
            ENDDO
          ENDIF
        
          ! *** NOW FINALIZE THE TEMPERATURE 
          DO LP=LF,LL  
            L=LWET(LP)  
            TEM(L,KC) = TEM(L,KC)+HPI(L)*RADNET(L,KC)  
          ENDDO  
        ENDDO   ! *** END OF DOMAIN
        !$OMP END DO
      ENDIF  
      !$OMP END PARALLEL

    ELSEIF( ISTOPT(2) == 2 )THEN  
      ! *** IMPLEMENT EXTERNALLY SPECIFIED EQUILIBRIUM TEMPERATURE FORMULATION  
      ! *** THIS OPTION USES CLOUDT IN A SPECIALIZED MANOR TO CONTAIN THE SURFACE EXCHANGE COEFFICIENT.  IT IS NOT FRACTION OF CLOUD COVER.
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL  
          L=LWET(LP)  
          TMPKC=DELT*DZIC(L,KC)
          TEM(L,KC)=TEM(L,KC)-TMPKC*CLOUDT(L)*HPI(L)*(TEM(L,KC)-TATMT(L))          
          TMIND(ND) = MIN(TEM(L,KC),TMIND(ND))
        ENDDO  
      ENDDO  

    ELSEIF( ISTOPT(2) == 3 )THEN  
      ! *** IMPLEMENT CONSTANT COEFFICIENT EQUILIBRIUM TEMPERATURE FORMULATION  
      DTHEQT=DELT*HEQT*DZI  
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL  
          L=LWET(LP)  
          TEM(L,KC)=TEM(L,KC)-DTHEQT*HPI(L)*(TEM(L,KC)-TEMO)  
          TMIND(ND) = MIN(TEM(L,KC),TMIND(ND))
        ENDDO  
      ENDDO  

    ELSEIF( ISTOPT(2) == 4 )THEN  
      ! *** IMPLEMENT W2 EQUILIBRIUM TEMPERATURE FORMULATION 
    
      IF(  .NOT. COMPUTESOLRAD )THEN
        ! *** MUST MAKE AT LEAST ONE CALL TO THIS TO INITIALIZE VARIABLES
        CALL SHORT_WAVE_RADIATION(CLOUDT(2),SRO,SRON)
      ENDIF
      
      ! *** SWRATNF - Background/Clear Water Extinction Coefficient
      ! *** SWRATNS - Light Extinction for TSS (1/m per g/m3)
      ! *** FSWRATF - Minimum Fraction of Solar Rad Absobed in the Surface Layer
      ! *** HTBED2  - Bottom Heat Exchange Coefficient (W/m2/s)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K)  &
      !$OMP                             PRIVATE(ET,CSHE,THICK,TFLUX,TEMP,C1,ET0,CSHE0,WSHLTR0,PSHADE0) 
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)

        ! *** INITIALIZE DOMAIN ET
        L=LWET(LF)
        CALL EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE) 
        PSHADE_OLD(ND) = PSHADE(L)
        PSHADE0 = PSHADE(L)
        WSHLTR_OLD(ND) = WINDSTKA(L)
        WSHLTR0 = WINDSTKA(L)
        ET0 = ET
        CSHE0 = CSHE
        
        ! *** SURFACE HEAT EXCHANGE PROCESSES
        DO LP=LF,LL  
          L=LWET(LP)  
          IF( .NOT. ICECELL(L) )THEN
            IF( NASER > 1 )THEN
              CALL EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE) 

            ELSEIF( PSHADE_OLD(ND) /= PSHADE(L) .OR. WSHLTR_OLD(ND) /= WINDSTKA(L) )THEN
              IF( WINDSTKA(L) == WSHLTR0 .AND. (PSHADE(L) == PSHADE0 .OR. SOLSWRT(L) == 0.0) )THEN
                ! *** USE PREVIOUS CONDITIONS
                CALL SWAP(ET,ET0)
                CALL SWAP(CSHE,CSHE0)
                CALL SWAP(PSHADE0,PSHADE_OLD(ND))
                CALL SWAP(WSHLTR0,WSHLTR_OLD(ND))
              ELSE
                ET0 = ET
                CSHE0 = CSHE
                PSHADE0 = PSHADE_OLD(ND)
                WSHLTR0 = WSHLTR_OLD(ND)
                CALL EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE) 
              ENDIF
              PSHADE_OLD(ND) = PSHADE(L)
              WSHLTR_OLD(ND) = WINDSTKA(L) 
            ENDIF
            
            ! *** SURFACE HEAT FLUX
            THICK = HPK(L,KC)
            TFLUX = CSHE*(ET-TEM(L,KC))/THICK*DELT
            TEM(L,KC) = TEM(L,KC)+TFLUX 
            TMIND(ND) = MIN(TEM(L,KC),TMIND(ND))
          ENDIF  ! *** END OF ICE COVER CHECK FOR CELL
        ENDDO

        ! *** Distribute Solar Radiation Across Water Column
        IF( LDAYLIGHT .AND. KC > 1 )THEN  
          ! *** NOW FINALIZE THE TEMPERATURE 
          DO K=1,KS
            ! *** RHO = 1000.0  Density (kg / m^3)  
            ! *** CP  = 4179.0  Specific Heat (J / kg / degC)
            ! *** RHOWCPI = 1/RHO/CP = 0.2393E-6 (m^3 degC)/(W s)
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              C1=DELT*DZIC(L,K)*RHOWCPI
              TEMP=TEM(L,K)+HPI(L)*C1*RADNET(L,K)
              TEM(L,K)=TEMP
            ENDDO  
          ENDDO  
      
        ENDIF    ! *** END OF ACTIVE SOLAR RADIATION

      ENDDO   ! *** END OF DOMAIN
      !$OMP END PARALLEL DO

    ELSEIF( ISTOPT(2) == 5 )THEN  
      ! *** FULL HEAT BALANCE WITH SPATIALLY AND TEMPORALLY VARIABLE EXTINCTION COEFFICIENTS
      C2 = -DELT
      C1 = DELT*RHOWCPI
      !$OMP PARALLEL DEFAULT(SHARED) 
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,SVPW1,HBLW,HBCD,HBCV,HBEV,HBNT,TEMP) 
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)

        ! *** SURFACE HEAT FLUX WITHOUT SOLAR RADIATION (m * degC)
        DO LP=LF,LL 
          L=LWET(LP)  
          
          ! *** UPDATE NET NETRAD FOR SOLAR RADIATION ADSORBTION. CONVERT TO M*degC FROM W/M2
          RADNET(L,KC) = RADNET(L,KC)*C1                                                          ! *** m*degC
          
          IF( .NOT. ICECELL(L) )THEN
            ! *** CELL WITHOUT ICE
            SVPW1=SVPW(L)
            HBLW = 1.312E-14*((TEM(L,KC)+273.)**4)*(0.39-0.05*SQRT(VPAT(L)))*(1.-.8*CLOUDT(L))    ! *** m*degC/s
            HBCD = 5.248E-14*((TEM(L,KC)+273.)**3)*(TEM(L,KC)-TATMT(L))                           ! *** m*degC/s
            HBCV = CCNHTT(L)*0.288E-3*WINDST(L)*(TEM(L,KC)-TATMT(L))                              ! *** m*degC/s
            HBEV = CLEVAP(L)*0.445*WINDST(L)*(SVPW1-VPAT(L))/PATMT(L)                             ! *** m*degC/s
            HBNT = C2*( HBLW + HBCD + HBCV + HBEV )                                               ! *** m*degC
            
            PMCTESTX(1,L) = RHOWCP*HBLW   !  W/M2   SAVE FOR USE WITH EE_ARRAYS
            PMCTESTX(2,L) = RHOWCP*HBCD   !  W/M2
            PMCTESTX(3,L) = RHOWCP*HBCV   !  W/M2
            PMCTESTX(4,L) = RHOWCP*HBEV   !  W/M2
            PMCTESTX(5,L) = RHOWCP*( HBLW + HBCD + HBCV + HBEV )    !  W/M2

            RADNET(L,KC) = RADNET(L,KC) + HBNT
            
            ! *** FINALIZE TEMPERATURE
            TEMP = TEM(L,KC) + HPKI(L,KC)*RADNET(L,KC)
            TEM(L,KC) = TEMP
            TMIND(ND) = MIN(TEM(L,KC),TMIND(ND))
          ELSE
            ! *** FINALIZE TEMPERATURE
            TEMP = TEM(L,KC) + HPKI(L,KC)*RADNET(L,KC)
            TEM(L,KC) = TEMP
          ENDIF
        ENDDO
      ENDDO   ! *** END OF DOMAIN
      !$OMP END DO
      
      IF( LDAYLIGHT .AND. KC > 1 )THEN
        ! *** INCORPORATE WATER COLUMN SOLAR RADIATION ON TEMPERATURE, RADNET (W/M2) WAS COMPUTED IN SET_LIGHT
        C1 = DELT*RHOWCPI
        !$OMP DO PRIVATE(ND,K,LP,L,TEMP) 
        DO ND=1,NDM  
          DO K=1,KS
            ! *** RHOWCPI = 1/RHO/CP = 0.2393E-6
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TEMP = TEM(L,K) + HPKI(L,K)*C1*RADNET(L,K)
              TEM(L,K) = TEMP
            ENDDO  
          ENDDO  
        ENDDO   ! *** END OF DOMAIN
        !$OMP END DO
      ENDIF          
      !$OMP END PARALLEL
    ENDIF    ! *** END OF HEAT CALCULATION OPTIONS

    ! *** BED/WATER COLUMN INTERFACE HEAT FLUX
    !$OMP PARALLEL DEFAULT(SHARED) 
    IF( NASER > 0 .AND. ISTOPT(2) /= 0 .AND. (HTBED1 + HTBED2) > 0. )THEN
      !$OMP DO PRIVATE(ND,LF,LL,LP,L,UBED,VBED,USPD,TFLUX,THICK,TEMP) 
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL  
          L=LWET(LP)  
          UBED = 0.5*( U(L,KSZU(L)) + U(LEC(L),KSZU(LEC(L))) )  
          VBED = 0.5*( V(L,KSZV(L)) + V(LNC(L),KSZV(LNC(L))) )  
          USPD = SQRT( UBED*UBED+VBED*VBED )  
          TFLUX = ( HTBED1*USPD + HTBED2 )*( TEMB(L) - TEM(L,KSZ(L)) )*DELT  

          ! *** ACCUMULATE BOTTOM HEAT FLUXTB (TO AVOID TRUNCATION ERROR, APPLY MINIMUM FLUXTB)
          FLUXTB(L) = FLUXTB(L) + TFLUX
          THICK = HPK(L,KSZ(L))  
          IF( ABS(FLUXTB(L))/THICK > 0.001 )THEN
            TEM(L,KSZ(L)) = TEM(L,KSZ(L)) + FLUXTB(L)/THICK 

            ! *** UPDATE BED TEMPERATURE
            IF( TBEDIT > 0. )THEN 
              TEMB(L) = TEMB(L) - FLUXTB(L)/TBEDTHK(L)
            ENDIF
            FLUXTB(L) = 0.
          ENDIF
        ENDDO
      ENDDO   ! *** END OF DOMAIN
      !$OMP END DO
    ENDIF  

    ! *** GET THE MINIMUM SURFACE TEMPERATURE
    IF( ISICE < 3 )THEN
      !$OMP DO PRIVATE(ND,K,LP,L)
      DO ND=1,NDM  
        K = KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          TMIND(ND) = MIN(TEM(L,K),TMIND(ND))
        ENDDO
      ENDDO   ! *** END OF DOMAIN
      !$OMP END DO
    ENDIF
    !$OMP END PARALLEL
    
    ! *** APPLY DRY CELL CORRECTIONS
    IF( ISDRY > 0 .AND. LADRY > 0 )THEN
      DO K=1,KC  
        DO LP=1,LADRY  
          L=LDRY(LP)

          ! *** ZERO HEAT ARRAYS
          RADBOT(L,K)=0.
          RADTOP(L,K)=0.
          RADNET(L,K)=0.
        ENDDO  
      ENDDO  
    ENDIF

    ! *** LIMIT WATER TEMPERATURES WHEN ICE FREEZE/MELT IS NOT COUPLED TO THE HEAT SUB-MODEL
    TMIN = MINVAL(TMIND)
    RPEMMINT = TMIN 
    IF( TMIN < 1.0 .AND. ISICE < 3 )THEN
      IF( ISICE == 0 )THEN
        ! *** LIMIT WATER TEMPERATURES TO > FREEZING TEMPERATURE IF ICE IS NOT SIMULATED
        K=KC
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)

          DO LP=LF,LL
            L=LWET(LP)  
            ! *** LIMIT THE WATER TEMPERATURE
            IF( ISTRAN(1) > 0 )THEN
              IF( TEM(L,KC) < -1.3 )THEN
                TEM(L,KC) = -1.3*(SAL(L,KC)/35.)
              ENDIF
            ELSE
              IF( TEM(L,KC) < 0.001 )THEN
                ! *** LIMIT THE WATER TEMPERATURE
                TEM(L,KC) = 0.001
              ENDIF
            ENDIF
          ENDDO
        ENDDO   ! *** END OF DOMAIN
        !$OMP END PARALLEL DO
      
      ELSEIF( ISICE == 1 .OR. ISICE == 2 )THEN
        ! *** FORCE A LIMITATION OF TEMPERATURE TO ICE TEMP FOR USER SPECIFIED ICE COVER
        K=KC
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)

          DO LP=LF,LL
            L=LWET(LP)  
            ! *** LIMIT THE WATER TEMPERATURE
            IF( TEM(L,K) < TEMPICE )THEN
              TEM(L,K) = TEMPICE
            ELSE
              TEM(L,K) = ICECOVER(L)*TEMPICE + (1.-ICECOVER(L))*TEM(L,K)
            ENDIF
          ENDDO
        ENDDO   ! *** END OF DOMAIN
        !$OMP END PARALLEL DO
      ENDIF
    ENDIF   ! *** END OF TMIN < 1.0 BLOCK
    
    ! *** CHECK NON-ICE "DRY" CELLS
    IF( LAWET < LA-1 .AND. ISICE > 2 .AND. TIMEDAY > TIMEDRY )THEN
      TIMEDRY = TIMEDRY + 0.25
      DO L=2,LA
        IF( .NOT. LMASKDRY(L) )THEN
          IF( ICETHICK(L) == 0.0 )THEN
            IF( TEM(L,KC) < -5 .OR. TEM(L,KC) > 30. )THEN
              TEM(L,:) = MAX(TATMT(L),0.0)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    
    RETURN  

  END SUBROUTINE CALHEAT

  SUBROUTINE ICECOMP(ISTL_, DELTD2, ICESTEP)
  
    ! *************************************************************************************************
    ! *** HEAT COUPLED ICE COMPUTATIONS WITH MASS BALANCE IN FREEZE/MELT
    ! *** 
    ! *** ISICE = 0 - ICE NOT INCLUDED
    ! *** ISICE = 1 - USER SPECIFIED TEMPORALLY/SPATIALLY VARYING ICE COVER 
    ! *** ISICE = 2 - USER SPECIFIED CTIME VARYING COMPLETE ICE COVER USING BINARY ON/OFF 
    ! *** ISICE = 3 - FULLY HEAT COUPLED 
    ! *** ISICE = 4 - FULLY HEAT COUPLED WITH FRAZIL ICE TRANSPORT
  
    ! CHANGE RECORD  
    ! DATE MODIFIED     BY               DESCRIPTION
    !------------------------------------------------------------------------------------------------!
    ! 2015-05           PAUL M. CRAIG    ADDED FULLY HEAT COUPLED ICE FORMATION/MELT
    !                   DANG H CHUNG     
  
    INTEGER, INTENT(IN) :: ISTL_ 
    REAL, INTENT(IN)    :: DELTD2,ICESTEP
    
    INTEGER :: L, K, ND, LF, LL, IFILE
    INTEGER,SAVE :: ITERMAX, ITERMAX2

    REAL       :: THICKMIN, MINICE, TMPVAL, ET, CSHE
    REAL       :: WTEMP,RHOA
    REAL       :: DEL,RANLW,RN,TF,TFS, RB,RC,RE,RT,TICEBOT
    REAL       :: RATIODENS,DICETHI,DICETHW,HICE
    
    REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: ICETEM_OLD

    REAL,EXTERNAL :: FUNDEN

    CHARACTER*9 :: FMODE
    
    IF( .NOT. ALLOCATED(ICETEM_OLD) )THEN
      ALLOCATE(ICETEM_OLD(LCM))
      ICETEM_OLD = -9999.

      ! *** RHOI = Density of Ice (kg / m^3)
      ! *** CP  = 4179.0  Specific Heat (J / kg / degC)
      ! *** RHOWICPI = 1/RHOI/CP
      RHOICP   = RHOI*CP
      RHOICPI  = 1./RHOICP
      RHOILHFI = 1./RHOI/LHF
      CREFLI   = 1./(1.-REFL)*(1.0-ALBEDOI)

      ITERMAX = 500
      ITERMAX2= 0
    ENDIF
    
    IFILE = -1
    IF( ( NCTBC /= NTSTBC .OR. IS2TL > 0 ) .AND. (LCHECKICE .OR. LFRAZIL) )THEN
      ! *** DENSITY OF WATER AT 0 DEG C = 999.8426 KG/M3

      IF( ISDYNSTP == 0 )THEN
        ! *** FIXED DELTA T
        DELT=DT
      ELSE
        ! *** DYNAMIC DELTA T
        DELT=DTDYN   
      ENDIF
      
      RATIODENS=RHOI/999.8426
      LFRAZIL = .FALSE.
        
      ! *** NEED TO LOOP OVER ALL THE CELLS BECAUSE SOME ICE COVERED CELLS MAY BE "DRY"
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND, LF, LL, L, K, ITER)                     &                       
      !$OMP                             PRIVATE(DICETHW, DICETHI, MINICE, TF, TFS, TMPVAL)  &
      !$OMP                             PRIVATE(WTEMP, RHOA, THICKMIN, ET, CSHE, TICEBOT)            &
      !$OMP                             PRIVATE(RANLW, RT, RB, RC, RE, RN, DEL, HICE)
      DO ND=1,NDM  
        LF=2+(ND-1)*LDM  
        LL=MIN(LF+LDM-1,LA)
        DO L=LF,LL  
          DICETHW   = 0.0   ! *** CHANGE IN ICE THICKNESS DUE TO WATER TEMPERATURES FOR OPEN WATER (NO ICE)/ICE WATER INTERFACE (WITH ICE COVER)  (M)
          DICETHI   = 0.0   ! *** CHANGE IN ICE THICKNESS DUE TO ICE TEMPERATURES, I.E. BOTTOM GROWTH (M)
          
          IF( ISTL_ == 3 .OR. IS2TL > 0 ) ICETHICK1(L) = ICETHICK(L)

          ! *** HANDLE SHALLOW CELLS
          MINICE = MINICETHICK
          IF( HP(L) <= HDRY*1.2 .AND. ISDRY > 0 )THEN
            MINICE = MINICETHICK*DZ
          ENDIF
          
          ! *** FREEZING TEMPERATURE OF WATER
          IF( ISTRAN(1) > 0 )THEN                            
            IF( SAL(L,KC) < 35. )THEN
              TF = -0.0545*SAL(L,KC)          
            ELSE
              TF = -0.31462-0.04177*SAL(L,KC)-0.000166*SAL(L,KC)*SAL(L,KC)    
            ENDIF                                                                                
          ELSE
            TF=0.0
          ENDIF
          TFS = TF - 0.01  ! *** SUPERCOOLED STATE

          ! *** DETERMINE ICE FORMATION AND HEAT FLUX DUE TO ICE FORMATION.  ALLOW ONLY FOR "WET" CELLS
          IF( ISICE == 4 .AND. .NOT. ICECELL(L) )THEN
            ! *** COUPLED HEAT/ICE MDOEL WITH FRAZIL ICE TRANSPORT
            DO K = KSZ(L),KC
              IF( TEM(L,K) <= TFS )THEN
                IF( HP(L) >= HDRYICE )THEN
                  ! *** FRAZIL ICE GROWTH
                  TMPVAL = -(TF-TEM(L,K))*RHOW(L,K)*CP*HP(L)*DZC(L,K)
                  TMPVAL = -TMPVAL*RHOILHFI
                  FRAZILICE(L,K) = FRAZILICE(L,K) + TMPVAL
                
                  ! *** REMOVE THE VOLUME OF WATER FROM THE WATER COLUMN 
                  ICEVOL(L) = ICEVOL(L) - TMPVAL*RATIODENS*DXYP(L)     

                  ! *** ICE FORMATION LATENT HEAT OF FUSION ADDS HEAT
                  IF( IS2TL == 0 )THEN
                    IF( ISTL_ == 3 )THEN
                      TEM1(L,K) = TF
                    ELSE
                      TEM(L,K) = TF                            
                    ENDIF
                  ELSE
                    TEM1(L,K) = TF
                    TEM(L,K)  = TF
                  ENDIF

                  ! ***  AMBIENT DENSITY
                  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
                    RHOW(L,K)=FUNDEN(SAL(L,KC),SEDT(L,KC),TF)
                  ELSE
                    RHOW(L,K)=FUNDEN(SAL(L,KC),0.,TF)
                  ENDIF
                  LFRAZIL = .TRUE.
                ENDIF

              ELSEIF( TEM(L,K) > 0.0 )THEN
                ! *** CHECK FOR MELTING FRAZIL ICE
                IF( FRAZILICE(L,K) > 0. )THEN
                  TMPVAL     = -DELT*HWI*(TEM(L,K)-TF)*RHOILHFI 
                  FRAZILICE(L,K) = FRAZILICE(L,K) + TMPVAL
                  FRAZILICE(L,K) = MAX(FRAZILICE(L,K),0.)
                  IF( FRAZILICE(L,K) > 0. )LFRAZIL = .TRUE.
                    
                  ! *** ICE MELT ADDS WATER AT T=0.0 DEGC
                  ICEVOL(L) = ICEVOL(L) - TMPVAL*RATIODENS*DXYP(L)
                ENDIF
                
                ! *** ACCOUNT FOR ICE "COVER" MELT
                IF( ICETHICK(L) > 0. .AND. K == KC )THEN
                  TMPVAL     = -DELT*HWI*(TEM(L,K)-TF)*RHOILHFI 
                  ICETHICK(L) = ICETHICK(L) + TMPVAL
                  ICETHICK(L) = MAX(ICETHICK(L),0.)
                  IF( ICETHICK(L) > 0. )LFRAZIL = .TRUE.
                    
                  ! *** ICE MELT ADDS WATER AT T=0.0 DEGC
                  ICEVOL(L) = ICEVOL(L) - TMPVAL*RATIODENS*DXYP(L)
                ENDIF

              ENDIF
            
            ENDDO
                
          ELSEIF( ICETHICK(L) < 1E-6 )THEN   ! *** ISICE = 3
              
            ! *** COUPLED HEAT MODEL WITHOUT FRAZIL ICE TRANSPORT
            IF( TEM(L,KC) < TFS )THEN

              ! *** GET ANY "FREEZING" WATER THAT MAY HAVE SETTLED/TRANSPORTED IN PRIOR TIMESTEPS WHERE ICE THICKNESS < MINICETHICK
              WTEMP  = 0.
              DO K=KSZ(L),KC
                IF( TEM(L,K) <= TF )THEN
                  WTEMP  = WTEMP  + (TF-TEM(L,K))*DZC(L,K)
                ENDIF
              ENDDO
              WTEMP = WTEMP*DZIC(L,KC)

              ! ***  AMBIENT DENSITY
              IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
                RHOA=FUNDEN(SAL(L,KC),SEDT(L,KC),TF)
              ELSE
                RHOA=FUNDEN(SAL(L,KC),0.,TF)
              ENDIF
                
              ! *** INITIAL ICE FORMATION FOR OPEN WATER
              THICKMIN = HP(L)*DZC(L,KC)
              IF( THICKMIN < HDRY )THEN
                THICKMIN = 0.
                DO K = KSZ(L),KC
                  IF( TEM(L,K) <= TF )THEN
                    THICKMIN  = THICKMIN  + HP(L)*DZC(L,K)
                  ENDIF
                ENDDO
              ENDIF
              TMPVAL  = WTEMP*RHOA*CP*THICKMIN
              DICETHW = TMPVAL*RHOILHFI
                
              ! *** CHECK IF SUFFICIENT WATER EXISTS
              IF( DICETHW > 0.1*HP(L) .AND. HP(L) < HDRYICE )THEN
                DICETHW = 0.0
              ENDIF
              
              ! *** UPDATE ALL THE TEMPERATURES
              DO K = KSZ(L),KC
                IF( TEM(L,K) <= TF )THEN
                  TEM1(L,K) = TF
                  TEM(L,K)  = TF                            ! *** ICE FORMATION LATENT HEAT OF FUSION ADDS HEAT

                  ! ***  AMBIENT DENSITY
                  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
                    RHOW(L,K)=FUNDEN(SAL(L,K),SEDT(L,K),TF)
                  ELSE
                    RHOW(L,K)=FUNDEN(SAL(L,K),0.,TF)
                  ENDIF
                ENDIF
              ENDDO
              IF( IS2TL == 0 )THEN
                TMPVAL = ICEVOL(L)/DXYP(L)
                H1P(L) = H1P(L) - DICETHW*RHOI/999.8426 + TMPVAL
                DO K=KSZ(L),KC  
                  H1PK(L,K) = H1P(L)*DZC(L,K)
                ENDDO
                IF( ISTL_ == 3 )HP(L) = HP(L) - DICETHW*RHOI/999.8426 + TMPVAL
                ! ICEVOL(L) should be zeroed here?
              ELSE
                ICEVOL(L) = ICEVOL(L) - DICETHW*RATIODENS*DXYP(L)    ! *** REMOVE THE VOLUME OF WATER FROM THE WATER COLUMN 
              ENDIF
                
            ENDIF       ! *** END BLOCK TEM(L,KC) < TF
          ENDIF         ! *** ISICE OPTION
        
          ! *** ICE BALANCE, IF ICE COVER ALREADY EXISTS FOR THE CELL. ALLOW ICE MELT EVEN FOR "DRY" CELLS
          IF( ICECELL(L) )THEN
            ! ******************************************************************************
            ! *** DETERMINE ICE TEMPERATURE AT TOP OF ICE COVER (Ts) AFTER REACHING 
            ! *** QUASI STEADY STATE CONDITIONS
            IF( DELTD2 >= ICESTEP )THEN

              IF( .NOT. LMASKDRY(L) ) CALL SET_LIGHT(L)   ! *** UPDATE DRY CELL SOLAR RADIATION CONDITIONS

              IF( ICETEM_OLD(L) == -9999. )THEN
                ICETEMP(L) = TATMT(L)  
              ELSE
                ICETEMP(L) = ICETEM_OLD(L) + SIGN(0.005,TATMT(L))
              ENDIF
    
              ! ** INCIDENT LONG WAVE RADIATION
              IF( TATMT(L) >= 5.0 )THEN
                RANLW = 5.31E-13*(273.15+TATMT(L))**6*(1.0+0.17*CLOUDT(L)**2)*0.97
              ELSE
                RANLW = 5.62E-8*(273.15+TATMT(L))**4*(1.-0.261*EXP(-7.77E-4*TATMT(L)**2))*(1.0+0.17*CLOUDT(L)**2)*0.97
              ENDIF         
        
              ! *** RT = TOTAL RADIATION AT TOP OF ICE, SHORTWAVE PLUS LONGWAVE
              RT = SOLSWRT(L)*CREFLI + RANLW

              ITER = 0
              DO WHILE ( ITER <= ITERMAX )     
                ITER = ITER+1
                CALL SURFACE_TERMS (ICETEMP(L),L,RB,RC,RE)          ! *** OUT: RB,RC,RE 
                RN = RT-RB-RE-RC                                    ! *** W/M2      
                DEL = RN + ICEK*(TF-ICETEMP(L))/ICETHICK(L)         ! *** W/M2
                TMPVAL = DEL*ICETHICK(L)/10.
                IF( ABS(TMPVAL) > 1. )THEN
                  TMPVAL = SIGN(1.0,TMPVAL)
                ENDIF
                ICETEMP(L) = ICETEMP(L) + TMPVAL
                IF( ABS(TMPVAL) < 0.01 ) EXIT
              ENDDO

              IF( ITER > ITERMAX2 )THEN
                ! PROVIDE SOME INFORMATION ON THE ITERATIONS
              
                IF( IFILE == -1 )THEN
                  IFILE = 8
                  OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
                ENDIF
                ITERMAX2 = ITER
              
                INQUIRE(FILE=OUTDIR//'EFDCLOG.OUT', ACTION = FMODE)
                IF( FMODE == 'WRITE    ' .OR. FMODE == 'READWRITE' )THEN
                  WRITE(8, "(A,2I5,2F8.3,F10.4)" )'ICE TEMPERATURE: NEW MAX ITERATION (MAX L ICETHK ICET TIME): ',ITER,L,ICETHICK(L),ICETEMP(L),TIMEDAY
              
                  IF( ITER > ITERMAX )THEN
                    WRITE(8,"(A,2I5,2F8.3,F10.4)")'ERROR: ICE TEMPERATURE ITERATIONS EXCEEDED MAX  ALLOWABLE: (MAX L ICETHK ICET TIME): ',ITER,L,ICETHICK(L),ICETEMP(L),TIMEDAY
                  ENDIF
                ENDIF
                IF( ITER > ITERMAX )THEN
                  PRINT '(A,2I5)','ERROR: ICE TEMPERATURE ITERATIONS EXCEEDED MAX  ALLOWABLE: (MAX L ICETHK ICET TIME): ',ITER,L
                  IF( ICETEM_OLD(L) == -9999. )THEN
                    ICETEMP(L) = MIN(0.5*TATMT(L),0.0)
                  ELSE
                    ICETEMP(L) = ICETEM_OLD(L)  ! *** USE LAST GOOD ICE TEMPERATURE
                  ENDIF
                ENDIF
              ENDIF
              ICETEM_OLD(L) = ICETEMP(L) 
            
              ! ******************************************************************************
              ! ** SOLAR RADIATION ATTENUATION THROUGH ICE IS HANDLED IN SET_LIGHT FUNCTION       
            
              ! *** GET INTERACTION AT THE ICE INTERFACE DUE TO ICE TEMPERATURE
              IF( ICETEMP(L) > 0.0 )THEN
                ! *** ICE MELT DUE TO ICE WARMING FROM AIR/ICE INTERFACE.  ICE IS PURE WATER SO, TF=0.0
                ICETEMP(L) = MIN(ICETEMP(L),0.001/ICETHICK(L))  
                DICETHI = -CP*ICETEMP(L)*ICETHICK(L)*MELTFACTOR/LHF*ICESTEP/DELT
                
                ! *** ICE MELT ADDS WATER AT T=0.0 DEGC
                ICEVOL(L) = ICEVOL(L) - DICETHI*RATIODENS*DXYP(L)
              ENDIF
            ENDIF

            IF( ICETEMP(L) <= 0.0 )THEN
              ! *** ICE GROWTH AT ICE/WATER INTERFACE (ICE BOTTOM)
              IF( LMASKDRY(L) .AND. ICETHICK(L) <= ICETHMX )THEN
                TICEBOT = MIN(ICETEMP(L) + TEM(L,KC),0.0)       ! *** ADJUST ICE TEMPERATURE FOR WATER CONTACT
                HICE    = ICEK*(TF-TICEBOT)/ICETHICK(L)         ! *** W/M2
                DICETHI = DELT*HICE*RHOILHFI

                ! *** ICE FORMATION LATENT HEAT OF FUSION CONSUMES THE HEAT
                ICEVOL(L) = ICEVOL(L) - DICETHI*RATIODENS*DXYP(L)
                
                IF( TEM(L,KC) > 1.0 )THEN
                  TEM(L,KC)  = TEM(L,KC)*MAX(1.-DICETHI*RATIODENS*HPKI(L,KC),0.0)
                  TEM1(L,KC) = TEM(L,KC)
                ENDIF
              ENDIF
            ELSE
              ICETEMP(L) =  0.0
            ENDIF
          ENDIF
          
          IF( ICETHICK(L) >= 1E-6 )THEN
            ! *** ICE MELT/FREEZE FROM WATER-ICE INTERFACE
            IF( TEM(L,KC) > 0.0 )THEN
              ! *** MELT
              HICE   = -HWI*TEM(L,KC)*RHOILHFI 
              DICETHW = DELT*HICE

              ! *** ICE MELT ADDS WATER AT T=0.0 DEGC
              ICEVOL(L) = ICEVOL(L) - DICETHW*RATIODENS*DXYP(L)
            ELSEIF( TEM(L,KC) < TFS )THEN
              ! *** FREEZE
              WTEMP = TF - TEM(L,KC)
              IF( HP(L) < HDRY )THEN
                THICKMIN = HP(L)
              ELSE
                THICKMIN = HP(L)*DZC(L,KC)
              ENDIF
              TMPVAL  = WTEMP*CP*THICKMIN*999.842
              DICETHW = TMPVAL*RHOILHFI
              
              IF( IS2TL == 0 )THEN
                TMPVAL = ICEVOL(L)/DXYP(L)
                H1P(L) = H1P(L) - DICETHW*RHOI/999.8426 + TMPVAL
                DO K=KSZ(L),KC  
                  H1PK(L,K) = H1P(L)*DZC(L,K)
                ENDDO
                IF( ISTL_ == 3 ) HP(L) = HP(L) - DICETHW*RHOI/999.8426 + TMPVAL
                TEM1(L,KC) = TF
                ! ICEVOL(L) should be zeroed here?
              ELSE
                ICEVOL(L) = ICEVOL(L) - DICETHW*RATIODENS*DXYP(L)    ! *** REMOVE THE VOLUME OF WATER FROM THE WATER COLUMN 
                TEM1(L,KC) = TF
                TEM(L,KC)  = TF
              ENDIF
            ENDIF
          
          ENDIF   ! *** END BLOCK ICECELL(L)
        
          ! ** TOTAL ICE THICKNESS AT THE TOP OF THE WATER COLUMN (M)
          ICETHICK(L) = ICETHICK(L) + (DICETHI + DICETHW)
          
          IF( ICETHICK(L) >= MINICE )THEN
            ICECELL(L) = .TRUE.
            ICECOVER(L) = 1.
            LFRAZIL = .TRUE.
          ELSE
            ICETEM_OLD(L) = -9999.
            IF( ISICE == 3 )THEN
              IF( ICETHICK(L) > 0.0 .AND. ICETHICK(L) < 1E-6 )THEN
                CALL EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE)
                IF( ET > 1.0 )THEN
                  ! *** COMPLETE MELT
                  IF( IS2TL == 0 )THEN
                    TMPVAL = ICEVOL(L)/DXYP(L)
                    H1P(L) = H1P(L) + ICETHICK(L)*RHOI/999.8426 + TMPVAL
                    DO K=KSZ(L),KC  
                      H1PK(L,K) = H1P(L)*DZC(L,K)
                    ENDDO
                    IF( ISTL_ == 3 ) HP(L) = HP(L) + ICETHICK(L)*RHOI/999.8426 + TMPVAL
                    TMPVAL = 0.
                    ICETHICK(L)=0. 
                  ELSE
                    TMPVAL = ICETHICK(L)*RATIODENS*DXYP(L)
                  ENDIF
                  ICEVOL(L) = ICEVOL(L) + TMPVAL
                  ICETHICK(L) = 0.
                  ICECOVER(L) = 0.
                ENDIF
              ENDIF
              
              ! *** ACCUMULATE ANY REMAINING ICE VOLUME AND ADD IT BACK TO THE WATER COLUMN
              IF( ICEVOL(L) /= 0. ) LFRAZIL = .TRUE.   ! ACCOUNT FOR FINAL "MELT"
            ELSE
              ICECOVER(L) = ICETHICK(L)/MINICETHICK
              IF( FRAZILICE(L,KC) > 0. .OR. ICETHICK(L) > 0. )LFRAZIL = .TRUE.   ! ALLOW FRAZIL ICE TRANSPORT EVEN IF NO ICE COVER OR NEW FORMATION
            ENDIF

            ! *** CHECK FOR SPECIAL CASES
            IF( ICECELL(L) )THEN
              ! *** COMPLETE MELT
              IF( ISDRY > 0 .AND. HP(L) < 3.*HDRY )THEN
                ! *** COMPUTE EQUILIBRIUM TEMPERATURE
                CALL EQUILIBRIUM_TEMPERATURE(WINDST(L), TDEWT(L), TATMT(L), SOLSWRT(L), ET, CSHE)
                ET = MAX(ET,0.1)
                IF( ET+2. < TEM(L,KC) )THEN
                  DO K=KSZ(L),KC
                    TEM(L,K) = ET
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
            
            ICECELL(L) = .FALSE.
          ENDIF
        ENDDO
      ENDDO   ! *** END OF DOMAIN
      !$OMP END PARALLEL DO

    ENDIF   ! *** END OF ISICE > 2 BLOCK

    IF( IFILE == 8 ) CLOSE(8)

    RETURN  

  END SUBROUTINE ICECOMP

  SUBROUTINE SHORT_WAVE_RADIATION(CLD,SRO,SRON) 
    ! ************************************************************************
    ! **                 S H O R T  W A V E  R A D I A T I O N              **
    ! **                      FROM CE-QUAL-W2 (VER 3.1)                     **
    ! ************************************************************************

    ! ******* Type declaration
    REAL, INTENT(IN)    ::  CLD 
    REAL, INTENT(OUT)   :: SRO, SRON
    INTEGER   :: IDAY 
    REAL      :: JDAY
    REAL      :: STANDARD, THOUR, PMC1, EQTNEW, H, DECL, SINAL, A0, CLD10
  
    ! ******* Input Conversions
    CLD10 = CLD*10.                ! *** CONVERT FROM 0-1 (EFDC) TO 0-10 (W2)

    ! ******* Shortwave Radiation
    STANDARD = 15.0*INT(DS_LONG/15.0)

    ! *** Day of the Year
    THOUR    = (TIMEDAY-INT(TIMEDAY))*24.0
    IDAY     =  TIMEDAY-INT(TIMEDAY/365.)*365.
    IDAY     =  IDAY+INT(INT(TIMEDAY/365.)/4.)
    JDAY     = REAL(IDAY)
    PMC1     = (2.*PI*(JDAY-1.))/365.
    EQTNEW   =  0.170*SIN(4.*PI*(JDAY-80.)/373.)-0.129*SIN(2.*PI*(JDAY-8.)/355.)
    H     = 0.2618*(THOUR-(DS_LONG-STANDARD)*0.066667+EQTNEW-12.0)
    DECL =  0.006918-0.399912*COS(PMC1)  +0.070257*SIN(PMC1)-0.006758*COS(2*PMC1)+0.000907*SIN(2*PMC1)-0.002697*COS(3*PMC1)+0.001480*SIN(3*PMC1)
    SINAL = SIN(DS_LAT*.017453)*SIN(DECL)+COS(DS_LAT*.017453)*COS(DECL)*COS(H)
    A0    = 57.2958*ASIN(SINAL)
    
    IF( A0 > 0.0 )THEN
      SRO  = 2.044*A0+0.1296*A0**2-1.941E-3*A0**3+7.591E-6*A0**4
      SRO  = (1.0-0.0065*CLD10**2)*SRO*24.0
      SRO  = SRO*BTU_FT2_DAY_TO_W_M2
      SRON = SRO*(1.-REFL)                    ! *** Adjust for surface reflection
    ELSE
      SRO  = 0.0
      SRON = 0.0
    END IF
       
    RETURN
    
  END SUBROUTINE SHORT_WAVE_RADIATION

  SUBROUTINE EQUILIBRIUM_TEMPERATURE(WSPD, TD, TAIR, SRON, ET, CSHE) 
    !************************************************************************
    !**             E Q U I L I B R I U M  T E M P E R A T U R E           **
    !**                      FROM CE-QUAL-W2 (VER 3.1)                     **
    !************************************************************************

    ! ******* Type declaration
    REAL, INTENT(IN)    :: WSPD, TD, TAIR, SRON 
    REAL, INTENT(OUT)   :: ET, CSHE
    REAL      :: TDEW_F, TAIR_F, WIND_2M 
    REAL      :: WIND_MPH
    REAL      :: SRO_BR, TSTAR, BETA, FW, RA, ETP
  
    INTEGER :: J  

    ! ******* Input Conversions

    ! ******* British units
    TDEW_F   = TD*1.8+32.0
    TAIR_F   = TAIR*1.8+32.0
    WIND_MPH = WSPD*MPS_TO_MPH
    WIND_2M  = WIND_MPH         ! *** EE7.3 AND LATER APPLY CORRECTION TO ALL WINDS
    WIND_2M = MAX(WIND_2M, 1.)  ! PMC - APPLY A MINUMUM TO ALLOW SURFACE HEAT EXCHANGE EVEN WITH CALM CONDITIONS

    ! *** SRON Should already be adjusted for Shading & Reflection
    SRO_BR   = SRON*W_M2_TO_BTU_FT2_DAY 

    ! ******* Equilibrium temperature and heat exchange coefficient
    
    ET    = TDEW_F
    TSTAR = (ET+TDEW_F)*0.5
    BETA  = 0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
    FW    = W_M2_TO_BTU_FT2_DAY*AFW+BCONV*BFW*WIND_2M**CFW
    CSHE  = 15.7+(0.26+BETA)*FW
    RA    = 3.1872E-08*(TAIR_F+459.67)**4
    ETP   = (SRO_BR+RA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE*(0.26+BETA))
    J     = 0
    DO WHILE (ABS(ETP-ET) > 0.05 .AND. J < 10)
      ET    = ETP
      TSTAR = (ET+TDEW_F)*0.5
      BETA  = 0.255-(8.5E-3*TSTAR)+(2.04E-4*TSTAR*TSTAR)
      CSHE  = 15.7+(0.26+BETA)*FW
      ETP   = (SRO_BR+RA-1801.0)/CSHE+(CSHE-15.7)*(0.26*TAIR_F+BETA*TDEW_F)/(CSHE*(0.26+BETA))
      J     = J+1
    END DO

    ! ******* SI units

    ! *** RHOWCPI = 1/RHO/CP = 0.2393E-6
    ET   = (ET-32.0)*5.0/9.0                 ! *** DEG C
    CSHE = CSHE*FLUX_BR_TO_FLUX_SI*RHOWCPI   ! *** M/S

  END SUBROUTINE EQUILIBRIUM_TEMPERATURE

  REAL FUNCTION GET_EXT_COEFF(L,K)
  
    ! *** COMPUTE THE LIGHT EXTINCTION COEFFICIENT
    INTEGER, INTENT(IN) :: L, K
    REAL    :: WCKESS, CHLKE, TSS

    IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
      TSS = SEDT(L,K)+SNDT(L,K) 
    ELSE
      TSS = 0.
    ENDIF
    
    IF( ISTRAN(8) > 0 )THEN
      ! *** If using WQ then use the WQ Coefficients
      WCKESS = WQKEB(IWQZMAP(L,K))                       ! *** BACKGROUND EXTINCTION (BY WQ ZONE)
      WCKESS = WCKESS + WQKETSS*TSS                      ! *** INORGANIC SOLIDS COMPONENT
      WCKESS = WCKESS + WQKEPOC*(WQV(L,K,4)+WQV(L,K,5))  ! *** PARTICULATE ORGANIC MATTER
      WCKESS = WCKESS + WQKEDOM*WQV(L,K,6)               ! *** DISSOLVED ORGANIC MATTER

      IF( WQKECHL < 0.0 )THEN                            ! *** CHLOROPHYLL
        ! *** Compute Extinction Factor as a fn(Chla)
        CHLKE = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)
      ELSE
        CHLKE = WQKECHL*WQCHL(L,K)**WQKECHLE              
      ENDIF  
      WCKESS = WCKESS + CHLKE  

    ELSE
      WCKESS = SVKEBACK(L) + WQKETSS*TSS
    ENDIF
    RADKE(L,K) = WCKESS
    GET_EXT_COEFF = WCKESS

  END FUNCTION GET_EXT_COEFF

  SUBROUTINE SET_LIGHT(L)
    ! *** SETS THE LIGHT EXTINCTION FACTORS AND SOLAR RADIATION FOR EACH LAYER
  
    ! *** RADBOT - SOLAR RADIATION AT THE BOTTOM OF THE LAYER (W/M2)
    ! *** RADNET - NET SOLAR RADIATION ABSORBED IN THE LAYER (W/M2)
    ! *** RADTOP - SOLAR RADIATION AT THE TOP OF THE LAYER (W/M2)
    ! *** RADKE  - EXTINCTION COEFFICIENT AT THE MIDPOINT OF THE LAYER

    INTEGER, INTENT(IN) :: L
    INTEGER :: K
    REAL    :: WCKESS, BOT, FRACLAYER, SOLBOT

    ! *** INITIALIZE CELL

    ! *** SOLSWRT HAS BEEN ADJUSTED BEFORE HERE TO INCLUDE REFLECTANCE AND ICE (IF ANY)
    IF( ICECELL(L) )THEN
      ! *** ADJUST INCIDENT SOLAR RADIATION AT THE TOP OF THE WATER DUE TO ICE
      SOLBOT = EXP(-GAMMAI*ICETHICK(L))
      BETAI = 1.0 - SOLBOT
      REFICE  = (1. - ALBEDOI)/(1.-REFL)*SOLBOT
      RADTOP(L,KC) = SOLSWRT(L)*REFICE
    ELSE
      RADTOP(L,KC) = SOLSWRT(L)
    ENDIF
    
    ! *** LOOP OVER THE LAYERS
    IF( ISTOPT(2) == 4 )THEN
      ! *** EQUILIBRIUM TEMP APPROACH
      DO K = KC,KSZ(L),-1
    
        ! *** Extinction Coefficient
        WCKESS = GET_EXT_COEFF(L,K)

        ! *** FRACTION OF LIGHT AT THE BOTTOM OF THE CURRENT LAYER
        BOT=MAX(-WCKESS*HPK(L,K),-40.0)
        FRACLAYER=EXP(BOT)
        IF( K == KC )THEN
          ! *** ENSURE AT LEAST THE FSWRATF FRACTION OF SRO IS ATTENUATED IN THE TOP LAYER
          IF( FRACLAYER > (1.-FSWRATF) .AND. HP(L) > 2.*HDRY ) FRACLAYER = 1.-FSWRATF
        ENDIF

        RADBOT(L,K) = RADTOP(L,K)*FRACLAYER

        ! *** Compute Net Energy (W/M2)
        RADNET(L,K) = RADTOP(L,K) - RADBOT(L,K)   

        ! *** UPDATE LAYER VARIABLES
        RADTOP(L,K-1) = RADBOT(L,K)
      ENDDO  
      
    ELSEIF( ISTOPT(2) == 5 )THEN
      ! *** FULL HEAT BALANCE APPROACH
      DO K = KC,KSZ(L),-1
    
        ! *** Extinction Coefficient
        WCKESS = GET_EXT_COEFF(L,K)

        ! *** FRACTION OF LIGHT AT THE BOTTOM OF THE CURRENT LAYER
        BOT=MAX(-WCKESS*HPK(L,K),-40.0)
        FRACLAYER=EXP(BOT)
        RADBOT(L,K) = RADTOP(L,K)*FRACLAYER

        ! *** Compute Net Energy (W/M2)
        RADNET(L,K) = RADTOP(L,K) - RADBOT(L,K)   

        ! *** UPDATE LAYER VARIABLES
        RADTOP(L,K-1) = RADBOT(L,K)
      ENDDO 
      
    ELSEIF( ISTRAN(8) > 0 )THEN
      ! *** SET EXTINCTION FOR WATER QUALITY FOR CASES WITHOUT SURFACE HEAT EXCHANGE (I.E. TEST CASES)
      DO K = KC,KSZ(L),-1
        WCKESS = GET_EXT_COEFF(L,K)
      ENDDO
    ENDIF
    
  END SUBROUTINE

  SUBROUTINE SURFACE_TERMS(TSUR,L,RB,RC,RE)
    ! ** CALCULATION OF HEAT EXCHANGE TERMS ON WATER SURFACE
    ! ** OUTPUT:
    ! ** FW,RE,RB,RC
  
    INTEGER,INTENT(IN)    :: L
    REAL,   INTENT(IN)    :: TSUR
    REAL,   INTENT(INOUT) :: RB,RC,RE
    REAL::EA,ES,TAIRV,DTV,DTVL,BOWEN_CONSTANT,FW
    DATA BOWEN_CONSTANT /0.47/
  
    ! *** ATMOSPHERIC CONDITIONS: VAPOR PRESSURE (mmHg)
    EA = EXP(2.3026*(7.5*TDEWT(L)/(TDEWT(L)+237.3)+0.6609))

    ! *** SATURATION VAPOR PRESSURE AT WATER TEMPERATURE
    IF( TSUR<0.0 )THEN
      ES = EXP(2.3026*(9.5*TSUR/(TSUR+265.5)+0.6609))
    ELSE
      ES = EXP(2.3026*(7.5*TSUR/(TSUR+237.3)+0.6609))
    ENDIF

    ! ** EVAPORATIVE WIND SPEED FUNCTION, W/M2/mmHG
    IF( ISRHEVAP == 1 )THEN
      TAIRV = (TATMT(L)+273.0)/(1.0-0.378*EA/760.0)    ! *** Tav
      DTV   = (TSUR+273.0)/(1.0-0.378*ES/760.0)-TAIRV  ! *** Tsv-Tav
      DTVL  =  0.0084*WINDST(L)**3
      IF( DTV < DTVL) DTV = DTVL
      FW = (3.59*DTV**0.3333+4.26*WINDST(L))           ! *** Rayan-Harleman,1974
    ELSE
      FW = AFWI+BFWI*WINDST(L)**CFWI
    END IF

    ! ** EVAPORATIVE HEAT LOSS, (W/M2)
    RE = FW*(ES-EA)

    ! ** HEAT CONDUCTION, (W/M2)
    RC = FW*BOWEN_CONSTANT*(TSUR-TATMT(L))

    ! ** BACK RADIATION FROM WATER SURFACE, (W/M2)
    RB = 5.51E-8*(TSUR+273.15)**4
  
  END SUBROUTINE
  
  SUBROUTINE SWAP (X1,X2)
    ! *** FUNCTION TO SWAP THE VALUES OF TWO VARIABLES
    REAL, INTENT(INOUT) :: X1,X2
    REAL                :: X
    X = X1
    X1 = X2
    X2 = X
    RETURN
  END SUBROUTINE
  
END MODULE
