SUBROUTINE WQSKE1  

  ! *** ORGINALLY CODED BY K.-Y. PARK 
  ! *** WQSKE1 - CE-QUAL-ICM KINETICS WITH UPDATES
  
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2013-03           Paul M. Craig    Added OMP
  ! 2012-08           Paul M. Craig    Added DOM component to light extinction WQKEDOM
  ! 2012-06           Paul M. Craig    Removed thickness from concentration based WQKETSS & WQKEPOC
  ! 2011-09           Paul M. Craig    Fixed Reaeration
  ! 2011-07           Paul M. Craig    Rewritten to F90
  ! 2008              SCOTT JAMES      ADDED CARBON DIOXIDE
  ! 2006-01-12        PAUL M. CRAIG    MAJOR REWRITE

  USE GLOBAL  

  IMPLICIT NONE

  INTEGER NQ,NW,NS,IZ,IMWQZ,NSTPTMP,IOBC
  INTEGER ND,LF,LL,LP,L,K,IFLAG
  
  ! *** DEBUG DECLARATIONS
  !INTEGER, STATIC :: LDBG
  !REAL   , STATIC :: TDBG

  REAL WQAVGIO,CNS1,RMULTMP,TIME,RLIGHT1,RLIGHT2
  REAL WQGNC,WQGND,WQGNG,WQGNM,WQGPM,WQF1NM,WQGPC,WQGPD,WQGPG
  REAL WQF1NC,WQF1ND,WQF1NG,WQKESS,XMRM,YMRM,WQTT1
  REAL WQFDI0,WQHTT,WQTTT
  REAL WQF2IC,WQF2ID,WQF2IG,SADWQ,WQGSD,WQTTB,WQISM,WQFDM,WQF2IM
  REAL UMRM,VMRM,WQVEL,WQLVF,WQF4SC,WQKDOC,WQKHP,WQTTS
  REAL WQKHN,WQTTM,TVAL1,TVAL2,TVAL3,TVAL4,TVAL5
  REAL RLNSAT1,RLNSAT2,XNUMER,XDENOM,WQLDF,WQTTC,WQTTD,WQTTG
  REAL WINDREA,WQWREA,WQVREA,WQA1C,WQVA1C,WQR1C,WQA2D
  REAL WQR2D,WQA3G,WQR3G,WQB4,WQA4,WQR4,WQC5,WQA5,WQR5
  REAL WQD6,WQA6C,WQA6D,WQA6G,WQA6,WQA6M,WQR6
  REAL WQE7,WQA7C,WQA7D,WQA7G,WQA7,WQR7
  REAL WQF8,WQA8C,WQA8D,WQA8G,WQA8,WQR8
  REAL WQF9,WQA9C,WQA9D,WQA9G,WQA9,WQR9
  REAL WQA10C,WQA10D,WQA10G,WQR10,WQKKL
  REAL WQI11,WQA11C,WQA11D,WQA11G,WQA11,WQR11
  REAL WQJ12,WQA12C,WQA12D,WQA12G,WQA12,WQR12
  REAL WQF13,WQA13C,WQA13D,WQA13G,WQA13,WQR13
  REAL WQR14,WQF14,WQA14C,WQA14D,WQA14G,WQA14
  REAL WQR15,WQA15C,WQA15D,WQA15G,WQA15,WQB15
  REAL WQM16,WQA16D,WQR16,WQR17,WQR18
  REAL TEMFAC,DTWQxH,DTWQxH2,WQA19C,WQA19D,WQA19G
  REAL WQA19,WQA19A,WQSUM,WQRea,WQPOC,WQDOC,WQNH3,WQCOD
  REAL WQT20,WQR21,TIMTMP,WQTAMD
  REAL PPCDO,TMP22, WQA22, WQA22C, WQA22D, WQA22G, WQCDDOC
  REAL WQCDREA, WQCDSUM
  REAL EXPA0,EXPA1                          ! VARIABLES FOR LIGHT EXTINCTION
  REAL WQGCO2M,WQGCO2C,WQGCO2G,WQGCO2D      !  CO2 Limitation Consts added by AA

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DZCHP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQISC
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQISD
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQISG
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQIBOT
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: WQI0TOP
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:)   :: WQO
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:,:) :: WQOLD

  ! ***  1) CHC - cyanobacteria
  ! ***  2) CHD - diatom algae
  ! ***  3) CHG - green algae
  ! ***  4) ROC - refractory particulate organic carbon
  ! ***  5) LOC - labile particulate organic carbon
  ! ***  6) DOC - dissolved organic carbon
  ! ***  7) ROP - refractory particulate organic phosphorus
  ! ***  8) LOP - labile particulate organic phosphorus
  ! ***  9) DOP - dissolved organic phosphorus
  ! *** 10) P4D - total phosphate
  ! *** 11) RON - refractory particulate organic nitrogen 23) macroalgae
  ! *** 12) LON - labile particulate organic nitrogen
  ! *** 13) DON - dissolved organic nitrogen
  ! *** 14) NHX - ammonia nitrogen
  ! *** 15) NOX - nitrate nitrogen
  ! *** 16) SUU - particulate biogenic silica
  ! *** 17) SAA - dissolved available silica
  ! *** 18) COD - chemical oxygen demand
  ! *** 19) DOX - dissolved oxygen
  ! *** 20) TAM - total active metal
  ! *** 21) FCB - fecal coliform bacteria
  ! *** 22) CO2 - dissolved carbon dioxide
  ! *** 23) macroalgae

  ! *** DTWQ - Water quality time step, which is in units of days
  ! *** DTWQO2 = DTWQ*0.5 

  ! *** WQCHL   = Chlorophyll a (ug/l)
  ! *** WQCHLC  = carbon-to-chlorophyll ratio for cyanobacteria (mg C / ug Chl)
  ! *** WQCHLD  = carbon-to-chlorophyll ratio for algae diatoms (mg C / ug Chl)
  ! *** WQCHLG  = carbon-to-chlorophyll ratio for algae greens (mg C / ug Chl)
  ! *** WQKECHL = Light Extinction Coeff for CHLa (1/m per mg/l)
  ! *** WQKETSS = Light Extinction Coeff for TSS (1/m per mg/l)
  ! *** WQKEPOC = Light Extinction Coeff for POC (1/m per mg/l)
  ! *** WQKEDOM = Light Extinction Coeff for DOC (1/m per mg/l)
  ! *** WQKETOT(L,K) = Total Light Extinction
  ! *** WQDOPG  = Optimal Depth for Growth - Green Algae

  ! *** RNH4WQ(L)  = Ammonia (for Current Layer)
  ! *** RNO3WQ(L)  = Nitrate (for Current Layer)
  ! *** PO4DWQ(L)  = Phosphate (for Current Layer)
  ! *** RNH4NO3(L) = Total Inorganic Nitrogen (for Current Layer)

  ! *** WQKHNG = Nitrogen half-saturation for Algae-Greens (mg/L)
  ! *** WQKHPG = Phosphorus half-saturation for Algae-Greens (mg/L)
  ! *** WQKHCO2G = CO2 half-saturation for Algae-Greens (mg/L)

  ! *** XLIMIG = Rate Limiting Factor - Light
  ! *** XLIMTG = Rate Limiting Factor - Temperature   (Lookup Table: WQTDGG)
  ! *** XLIMNG = Rate Limiting Factor - Nitrogen      (Local-WQGNG)
  ! *** XLIMPG = Rate Limiting Factor - Phosphorus    (Local-WQGPG)
  ! *** XLIMCO2G = Rate Limiting Factor - CO2    (Local-WQGPG)
  ! *** WQF1NG = Rate Limiting Factor, Minimum of N & P

  ! *** WQPMG  = Maximum Growth Rate for Algae-Greens (1/d)
  ! *** WQPG   = Current Growth Rate for Algae-Greens (1/d)
  ! *** WQBMG  = Current Basal Metabolism Rate for Algae-Greens (1/d)
  ! *** WQPRG  = Current Predation Metabolism Rate for Algae-Greens (1/d)

  ! *** WQBMRG   = Basal Metabolism Rate for Algae-Greens (1/d)
  ! *** WQPRRG   = Predation Rate for Algae-Greens (1/d)
  ! *** WQTDRG   = Lookup Table for Temperature Rate Effect - Algae-Greens

  ! *** WQPC   = Final Net Growth Rate - Cyanobacteria
  ! *** WQPD   = Final Net Growth Rate - Diatoms Algae  
  ! *** WQPG   = Final Net Growth Rate - Green Algae
  ! *** WQPM   = Final Net Growth Rate - Macroalgae

  ! *** WQPNC = Preference for ammonium uptake - Cyanobacteria
  ! *** WQPND = Preference for ammonium uptake - Diatoms Algae
  ! *** WQPNG = Preference for ammonium uptake - Green Algae

  ! *** WQOBTOT  = Total Algal Biomass (mg/l)
  ! *** WQKRC    = Minimum Dissolution Rate of Refractory POC (1/day)
  ! *** WQKLC    = Minimum Dissolution Rate of Labile POC (1/day)
  ! *** WQKLCALG = Constant Refractory POC Dissolution Rate
  ! *** WQTDHDR  = Lookup Table for Temperature Rate Effect for Hydrolysis
  ! *** WQKRPC   = Current Dissolution Rate for POC

  ! *** WQI0   = SOLAR RADIATION for Current Time
  ! *** WQI1   = SOLAR RADIATION ON PREVIOUS DAY  
  ! *** WQI2   = SOLAR RADIATION TWO DAYS AGO  
  ! *** WQI3   = SOLAR RADIATION THREE DAYS AGO  

  ! *** WQKHR  = DOC Heterotrophic Respiration Rate

  ! *** WQWSSET = Water quality settling speed L;(L:1) is for top water layer; (L,2) is for lower water layers
  ! *** WQTTM   = Temporary concentration variable
  IF(  .NOT. ALLOCATED(DZCHP) )THEN
    ALLOCATE(DZCHP(LCM))  
    ALLOCATE(WQISC(LCM))  
    ALLOCATE(WQISD(LCM))  
    ALLOCATE(WQISG(LCM))  
    ALLOCATE(WQIBOT(LCM))                ! *** Solar Radiation at the bottom of the current layer
    ALLOCATE(WQI0TOP(LCM))               ! *** Solar Radiation at the surface, adjusted to shade
    ALLOCATE(WQO(LCM,NWQVM))
    ALLOCATE(WQOLD(NBCSOP,KCM,0:NWQVM))  ! *** ALLOCATE MEMORY FOR VARIABLE TO STORE CONCENTRATIONS AT OPEN BOUNDARIES
    DZCHP=0.0
    WQISC=0.0
    WQISD=0.0
    WQISG=0.0
    WQIBOT=0.0
    WQI0TOP=0.0
    WQO=0.0
    WQOLD=0.0
    
    !TDBG=TIMEDAY+ 0.5
    !LDBG=1137
    !OPEN(11,FILE=OUTDIR//'LIGHT.LOG',STATUS='UNKNOWN')
    !WRITE(11,'(A)')'   TIMEDAY    K   SOLSWRT    WQKESS      WQI0    WQITOP    WQIBOT  ICETHICK    WQF2IG'
  ENDIF
  
  ! ***  COMPUTE THE CURRENT OPTIMAL LIGHT INTENSITY
  IF( IWQSUN  ==  2 )THEN  
    WQAVGIO = WQCIA*WQI1 + WQCIB*WQI2 + WQCIC*WQI3  
  ELSE
    WQAVGIO = WQCIA*WQI0 + WQCIB*WQI1 + WQCIC*WQI2
  ENDIF 

  ! *** HARDWIRE TO SET OPTIMAL GROWTH TO A FIXED VALUE USED BY CHAPRA (Algal Growth Example (Chapra 33.2 with Fixed IS250).
  ! *** SET TO FALSE FOR NORMAL CASE
  IF( .FALSE. )THEN
      WQAVGIO = 250.        ! *** Is
      !WQFD = 0.5
      !DO L=LF,LL
      !  WQI0TOP(L) = 400.   ! *** Ia
      !ENDDO
  ENDIF

  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO NQ=0,NWQVM
    DO K=1,KC
      DO IOBC=1,NBCSOP  
        L=LOBCS(IOBC)
        WQOLD(IOBC,K,NQ) = WQV(L,K,NQ)
      ENDDO
    ENDDO  
  ENDDO  
  
  CNS1=2.718  
  NS=1  
  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NQ,NW,IZ,IMWQZ)  &
  !$OMP             PRIVATE(WQGNC,WQGND,WQGNG,WQGNM,WQGPM,WQF1NM,WQGPC,WQGPD,WQGPG)  &
  !$OMP             PRIVATE(WQF1NC,WQF1ND,WQF1NG,WQKESS,XMRM,YMRM,WQTT1)  &
  !$OMP             PRIVATE(WQFDI0,WQHTT,WQTTT)  &
  !$OMP             PRIVATE(WQF2IC,WQF2ID,WQF2IG,SADWQ,WQGSD,WQTTB,WQISM,WQFDM,WQF2IM)  &
  !$OMP             PRIVATE(UMRM,VMRM,WQVEL,WQLVF,WQF4SC,WQKDOC,WQKHP,WQTTS)  &
  !$OMP             PRIVATE(WQKHN,WQTTM,TVAL1,TVAL2,TVAL3,TVAL4,TVAL5)  &
  !$OMP             PRIVATE(RLNSAT1,RLNSAT2,XNUMER,XDENOM,WQLDF,WQTTC,WQTTD,WQTTG)  &
  !$OMP             PRIVATE(WINDREA,WQWREA,WQVREA,WQA1C,WQVA1C,WQR1C,WQA2D)  &
  !$OMP             PRIVATE(WQR2D,WQA3G,WQR3G,WQB4,WQA4,WQR4,WQC5,WQA5,WQR5)  &
  !$OMP             PRIVATE(WQD6,WQA6C,WQA6D,WQA6G,WQA6,WQA6M,WQR6)  &
  !$OMP             PRIVATE(WQE7,WQA7C,WQA7D,WQA7G,WQA7,WQR7)  &
  !$OMP             PRIVATE(WQF8,WQA8C,WQA8D,WQA8G,WQA8,WQR8)  &
  !$OMP             PRIVATE(WQF9,WQA9C,WQA9D,WQA9G,WQA9,WQR9)  &
  !$OMP             PRIVATE(WQA10C,WQA10D,WQA10G,WQR10,WQKKL)  &
  !$OMP             PRIVATE(WQI11,WQA11C,WQA11D,WQA11G,WQA11,WQR11)  &
  !$OMP             PRIVATE(WQJ12,WQA12C,WQA12D,WQA12G,WQA12,WQR12)  &
  !$OMP             PRIVATE(WQF13,WQA13C,WQA13D,WQA13G,WQA13,WQR13)  &
  !$OMP             PRIVATE(WQR14,WQF14,WQA14C,WQA14D,WQA14G,WQA14)  &
  !$OMP             PRIVATE(WQR15,WQA15C,WQA15D,WQA15G,WQA15,WQB15)  &
  !$OMP             PRIVATE(WQM16,WQA16D,WQR16,WQR17,WQR18)  &
  !$OMP             PRIVATE(TEMFAC,DTWQxH,DTWQxH2,WQA19C,WQA19D,WQA19G)  &
  !$OMP             PRIVATE(WQA19,WQA19A,WQSUM,WQRea,WQPOC,WQDOC,WQNH3,WQCOD)  &
  !$OMP             PRIVATE(WQT20,WQR21,TIMTMP,WQTAMD)  &
  !$OMP             PRIVATE(PPCDO,TMP22, WQA22, WQA22C, WQA22D, WQA22G, WQCDDOC)  &
  !$OMP             PRIVATE(WQCDREA, WQCDSUM, EXPA0,EXPA1)                 &
  !$OMP             PRIVATE(WQGCO2M,WQGCO2C,WQGCO2G,WQGCO2D)
  DO ND=1,NDM 
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
      
    ! COMPUTE WQCHL,WQTAMP,WQPO4D,WQSAD AT A NEW TIME STEP: WQCHLX=1/WQCHLX  

    ! ***************************************************************************
    ! *** INITIALIZE SOLAR RADIATION AND OPTIMAL LIGHT
    IF( LDAYLIGHT .AND. IWQSUN == 2 )THEN
      ! *** INITIAL SOLAR RADIATION AT TOP OF SURFACE LAYER (SHADING AND ICECOVER ALREADY ACCOUNTED FOR)
      DO LP=LF,LL
        L=LWET(LP)
        WQI0TOP(L) = PARADJ*2.065*RADTOP(L,KC)   ! *** SOLAR RADIATION IN LANGLEYS/DAY 
      ENDDO
      !PRINT '(F10.4,2(4X,3F9.2))',TIMEDAY,WQI0TOP(1221),RADTOP(1221,KC),SOLSWRT(1221),WQI0TOP(1045),RADTOP(1045,KC),SOLSWRT(1045)
    ELSE
      ! *** INITIAL SOLAR RADIATION AT TOP OF SURFACE LAYER (SHADING AND ICECOVER ALREADY ACCOUNTED FOR)
      DO LP=LF,LL
        L=LWET(LP)
        WQI0TOP(L) = WQI0                        ! *** SOLAR RADIATION IN LANGLEYS/DAY 
      ENDDO
    ENDIF

    ! ***************************************************************************
    ! *** DZWQ=1/H (for a layer), VOLWQ=1/VOL (m^-3)
    DO K=KC,1,-1  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        TWQ(L)    = TEM(L,K)              !layer temperature for WQ calcs
        SWQ(L)    = MAX(SAL(L,K), 0.0)    !layer salinity for WQ calcs
        DZCHP(L)  = HPK(L,K)              !layer thickness of a cell in meters
        DZWQ(L)   = HPKI(L,K)             !inverse layer thickness
        VOLWQ(L)  = DZWQ(L)*DXYIP(L)      !inverse volume of each cell in a layer
        IMWQZT(L) = IWQZMAP(L,K)          !WQ Zone Map for current layer
      ENDDO  
      
      ! *** ZERO WQWPSL IF FLOWS ARE NEGATIVE.  THESE ARE HANDLED IN CALFQC (PMC)
      IF( IWQPSL /= 2 )THEN
        DO NQ=1,NQSIJ  
          IF( (QSERCELL(K,NQ)+QSS(K,NQ)) <= 0.0 )THEN
            ! *** ZERO THE FLUX
            L=LQS(NQ)  
            DO NW=1,NWQV
              WQWPSL(L,K,NW)=0.0
            ENDDO
          ENDIF
        ENDDO
      ENDIF
            
      ! *** DETERMINE THE RATE OF ALGAE LEAVING THE CELL THROUGH SETTLING OR FLOATING
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        WQBCSET(L,1) = 0.0
        
        ! *** CYANOBACTERIA
        IF( WQWSC(IMWQZT(L)) < 0 )THEN
          IF( K == KC )THEN      
            WQBCSET(L,1) = 0.0                         ! *** ALGAE AT THE WATER SURFACE CAN'T LEAVE CELL
          ELSE
            WQBCSET(L,1) = -WQWSC(IMWQZT(L))*DZWQ(L)   ! *** CYANOBACTERIA  !NEEDS TO BE A POSITIVE QTY  
          ENDIF
        ELSE
          WQBCSET(L,1) = WQWSC(IMWQZT(L))*DZWQ(L)   
        ENDIF

        ! *** DIATOMS
        IF( WQWSD(IMWQZT(L)) < 0 )THEN                  ! *** PERMITS DIATOMS TO FLOAT AND/OR SETTLE
          IF( K == KC )THEN
            WQBDSET(L,1) = 0.0
          ELSE
            WQBDSET(L,1) = -WQWSD(IMWQZT(L))*DZWQ(L)    
          ENDIF
        ELSE
          WQBDSET(L,1) = WQWSD(IMWQZT(L))*DZWQ(L) 
        ENDIF
  
        ! *** GREEN ALGAE  
        IF( WQWSG(IMWQZT(L)) < 0 )THEN                  ! ***  PERMITS GREEN ALGAE TO FLOAT AND/OR SETTLE
          IF( K == KC )THEN
            WQBGSET(L,1) = 0.0
          ELSE
            WQBGSET(L,1) = -WQWSG(IMWQZT(L))*DZWQ(L)  
          ENDIF
        ELSE
          WQBGSET(L,1) = WQWSG(IMWQZT(L))*DZWQ(L)
        ENDIF

        ! *** ZONE SPECIFIC SETTING VELOCITIES for POM, (m/day)   
        WQRPSET(L,1) = WQWSRP(IMWQZT(L))*DZWQ(L)  ! *** Refractory POM 
        WQLPSET(L,1) = WQWSLP(IMWQZT(L))*DZWQ(L)  ! *** Labile POM 
      ENDDO

      ! *** SET SETTLING FOR TAM SORPTION: CURRENT LAYER  
      IF( IWQSRP == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQWSSET(L,1) = WQWSS(IMWQZT(L))*DZWQ(L)  
        ENDDO  
      ENDIF  

      ! *** SET CURRENT LAYER WQ ZONE FOR LAYERS ABOVE AND BELOW
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IF( K /= KC )THEN 
          IMWQZT1(L)=IWQZMAP(L,K+1)   ! *** IMWQZT1 - LAYER ABOVE CURRENT LAYER
        ENDIF
        IF( K /= KSZ(L) )THEN
          IMWQZT2(L)=IWQZMAP(L,K-1)   ! *** IMWQZT2 - LAYER BELOW CURRENT LAYER
        ENDIF
      ENDDO 

      ! *** COMPUTE THE MATERIAL COMING INTO THE CURRENT LAYER FROM LAYER ABOVE
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        ! *** ALL LAYERS EXCEPT THE TOP LAYER
        IF( K /= KC )THEN 
          IF( WQWSC(IMWQZT1(L)) < 0 )THEN
            WQBCSET(L,2) = 0.0  
          ELSE
            WQBCSET(L,2) = WQWSC(IMWQZT1(L))*DZWQ(L)  
          ENDIF
          IF( WQWSD(IMWQZT1(L)) < 0 )THEN
            WQBDSET(L,2) = 0.0  
          ELSE
            WQBDSET(L,2) = WQWSD(IMWQZT1(L))*DZWQ(L)  
          ENDIF
          IF( WQWSG(IMWQZT1(L)) < 0 )THEN
            WQBGSET(L,2) = 0.0  
          ELSE
            WQBGSET(L,2) = WQWSG(IMWQZT1(L))*DZWQ(L)  
          ENDIF
        ENDIF

        ! *** ALL LAYERS EXCEPT THE BOTTOM LAYER
        IF( K /= KSZ(L) )THEN
          IF( WQWSC(IMWQZT2(L)) < 0 )THEN
            WQBCSET(L,2) = WQBCSET(L,2)-WQWSC(IMWQZT2(L))*DZWQ(L)  
          ENDIF
          IF( WQWSD(IMWQZT2(L)) < 0 )THEN
            WQBDSET(L,2) = WQBDSET(L,2)-WQWSD(IMWQZT2(L))*DZWQ(L) 
          ENDIF
          IF( WQWSG(IMWQZT2(L)) < 0 )THEN
            WQBGSET(L,2) = WQBGSET(L,2)-WQWSG(IMWQZT2(L))*DZWQ(L)  
          ENDIF
        ENDIF
      ENDDO

      IF( K /= KC )THEN
        ! *** Flux of particulate into the layer from the layer above
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQRPSET(L,2) = WQWSRP(IMWQZT1(L))*DZWQ(L)  
          WQLPSET(L,2) = WQWSLP(IMWQZT1(L))*DZWQ(L)  
        ENDDO  
      ENDIF

      ! *** Set settling for tam sorption: One layer up
      IF( IWQSRP == 1 .AND. K /= KC )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQWSSET(L,2) = WQWSS(IMWQZT1(L))*DZWQ(L)  
        ENDDO  
      ENDIF

      ! *** FIND AN INDEX FOR LOOK-UP TABLE FOR TEMPERATURE DEPENDENCY  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        
        IWQT(L)=NINT((TWQ(L)-WQTDMIN)/WQTDINC)+1  
        IF( IWQT(L) < 1 .OR. IWQT(L) > NWQTD )THEN  
          OPEN(1,FILE=OUTDIR//'ERROR.LOG',POSITION='APPEND',STATUS='UNKNOWN')  
          WRITE(1,*)' *** ERROR IN WQSKE1:TEMPERATURE LOOKUP TABLE'
          WRITE(1,911) TIMEDAY, L, IL(L), JL(L), K, TWQ(L),TEM1(L,K), HP(L), H1P(L)
          WRITE(6,600)IL(L),JL(L),K,TWQ(L),HP(L)
            
          WRITE (6,'(A)')'SURROUNDING DEPTHS'
          WRITE (6,'(2X,A14,I5,4F14.4)')'HP  ',L,HP(LWC(L)), HP(LEC(L)), HP(LSC(L)), HP(LNC(L))
          WRITE (6,'(2X,A14,I5,4F14.4)')'H1P ',L,H1P(LWC(L)),H1P(LEC(L)),H1P(LSC(L)),H1P(LNC(L))
          WRITE (6,'(A)')'FLUX TERMS'
          WRITE (6,'(2X,A14,I5,4E14.6)')'UHDYE/VHDXE'  ,L,UHDYE(L), UHDYE(LEC(L)), VHDXE(L), VHDXE(LNC(L))
          WRITE (6,'(2X,A14,I5,4E14.6)')'UHDY1E/VHDX1E',L,UHDY1E(L),UHDY1E(LEC(L)),VHDX1E(L),VHDX1E(LNC(L))

          IWQT(L)=MAX(IWQT(L),1)  
          IWQT(L)=MIN(IWQT(L),NWQTD)  
          CLOSE(1,STATUS='KEEP')
        ENDIF  
      ENDDO  

      600 FORMAT(' I,J,K,TEM,HP = ',3I5,4E12.4)  
      911 FORMAT('ERROR: TIME, L, I, J, K, TWQ, TEM, HP, H1P = ',F10.5, 4I4, 4E12.4,/)  

      !C NOTE: MRM 04/29/99  ADDED ARRAYS TO KEEP TRACK OF  
      !C       NITROGEN, PHOSPHORUS, LIGHT, AND TEMPERATURE LIMITS  
      !C       FOR ALGAE GROWTH FOR CYANOBACTERIA, DIATOMS, GREENS,  
      !C       AND MACROALGAE.  THESE ARE THE ARRAYS:  
      !C        XLIMNX(L,K) = NITROGEN    LIMITATION FOR ALGAE GROUP X  
      !C        XLIMPX(L,K) = PHOSPHORUS  LIMITATION FOR ALGAE GROUP X  
      !C        XLIMIX(L,K) = LIGHT       LIMITATION FOR ALGAE GROUP X  
      !C        XLIMTX(L,K) = TEMPERATURE LIMITATION FOR ALGAE GROUP X  
   
      ! *** BEGIN HORIZONTAL LOOP FOR ALGAE PARMETERS  

      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        RNH4WQ(L) = MAX (WQV(L,K,14), 0.0)  ! *** Ammonia
        RNO3WQ(L) = MAX (WQV(L,K,15), 0.0)  ! *** Nitrate
        PO4DWQ(L) = MAX (WQPO4D(L,K), 0.0)  ! *** Phosphate
        RNH4NO3(L) = RNH4WQ(L) + RNO3WQ(L)  ! *** Total Inorganic Nitrogen
        IF( LMASKDRY(L) )THEN
          WQGNC = RNH4NO3(L) / (WQKHNC+RNH4NO3(L)+ 1.E-18)  
          WQGND = RNH4NO3(L) / (WQKHND+RNH4NO3(L)+ 1.E-18)  
          WQGNG = RNH4NO3(L) / (WQKHNG+RNH4NO3(L)+ 1.E-18)  
          WQGPC = PO4DWQ(L) / (WQKHPC+PO4DWQ(L)+ 1.E-18)  
          WQGPD = PO4DWQ(L) / (WQKHPD+PO4DWQ(L)+ 1.E-18)  
          WQGPG = PO4DWQ(L) / (WQKHPG+PO4DWQ(L)+ 1.E-18)  
          IF( ISTRWQ(22) > 0 )THEN
            CO2WQ(L)  = MAX (WQV(L,K,22), 0.0)  ! *** CO2       !AA Added
            WQGCO2C = CO2WQ(L) / (WQKHCO2C+CO2WQ(L)+ 1.E-18)    !AA
            WQGCO2D = CO2WQ(L) / (WQKHCO2D+CO2WQ(L)+ 1.E-18)    !AA
            WQGCO2G = CO2WQ(L) / (WQKHCO2G+CO2WQ(L)+ 1.E-18)    !AA
          ENDIF

          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN            
            WQGNM = RNH4NO3(L) / (WQKHNM+RNH4NO3(L) + 1.E-18)  
            WQGPM = PO4DWQ(L)  / (WQKHPM+PO4DWQ(L)  + 1.E-18) 
            IF( ISTRWQ(22) > 0 )THEN                
              WQGCO2M = CO2WQ(L) / (WQKHCO2M+CO2WQ(L) + 1.E-18)    !AA
              WQF1NM = MIN(WQGNM, WQGPM, WQGCO2M)                  !AA  Minimum of the N/P/CO2 Limit: Mats
            ELSE                          
              WQF1NM = MIN(WQGNM, WQGPM)                           !AA  Minimum of the N/P     Limit: Mats
            ENDIF
          ENDIF

          ! *** Limit: Cyanobacteria
          IF( ISTRWQ(22) > 0 )THEN 
            WQF1NC = MIN(WQGNC, WQGPC, WQGCO2C)        !AA  Minimum of the N/P/CO2     
          ELSE
            WQF1NC = MIN(WQGNC, WQGPC)                 ! *** Minimum of the N/P     
          ENDIF

          ! *** Limit: Diatoms
          IF( IWQSI == 1 )THEN  
            SADWQ = MAX (WQSAD(L,K), 0.0)  
            WQGSD = SADWQ / (WQKHS+SADWQ+ 1.E-18)
            WQF1ND = MIN(WQGND, WQGPD, WQGSD)          ! *** Minimum of the N/P/S
          ELSE
            WQF1ND = MIN(WQGND, WQGPD)                 ! *** Minimum of the N/P
          ENDIF 
          IF( ISTRWQ(22) > 0 )WQF1ND = MIN(WQF1ND, WQGCO2D)     ! ***  Minimum of the CO2

          ! *** Limit: Greens
          IF( ISTRWQ(22)  > 0 )THEN  
            WQF1NG = MIN(WQGNG, WQGPG, WQGCO2G)         ! ***  Minimum of the N/P/CO2         
          ELSE
            WQF1NG = MIN(WQGNG, WQGPG)                  ! ***  Minimum of the N/P 
          ENDIF
                                      
          ! *** IN C&C, F2IC=F2IC/FCYAN, FACTOR TO ALLOW CYANOBACTERIA MAT FORMATION  
          ! *** LIGHT EXTINCTION (THIS WILL ALWAYS BE TRUE EXCEPT FOR IWQSUN=2)
          IF( WQI0 > 0.1 )THEN
            ! *** GET THE EXTINCTION COEFFICIENT
            WQKESS = RADKE(L,K)
            
            ! *** OPTIMAL LIGHT INTENSITY AT OPTIMAL DEPTH
            IF( K == KC )THEN
              WQISC(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPC),WQISMIN)  
              WQISD(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPD),WQISMIN)  
              WQISG(L) = MAX( WQAVGIO*EXP(-WQKESS*WQDOPG),WQISMIN) 
            ENDIF

            ! *** CURRENT LIGHT GROWTH LIMITING FACTOR. 
            ! *** WQITOP is solar radiation at the TOP of layer K
            ! *** WQIBOT is solar radiation at the BOTTOM of layer K

            IF( K == KC )THEN    
              WQITOP(L,K) = WQI0TOP(L)
              WQIBOT(L)   = WQI0TOP(L)*EXP(-WQKESS*DZCHP(L)) 
            ELSE            
              WQITOP(L,K) = WQIBOT(L)
              WQIBOT(L)   = WQITOP(L,K)*EXP(-WQKESS*DZCHP(L)) 
            ENDIF !SEE DiTORO ET AL (1971, EQNS. (11)&(12)) 

            WQTT1 = EXP(1.0)*WQFD
            
            EXPA0=EXP(-WQITOP(L,K)/WQISC(L))
            EXPA1=EXP(-WQIBOT(L)  /WQISC(L))
            WQF2IC=WQTT1/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0)
          
            EXPA0=EXP(-WQITOP(L,K)/WQISD(L))
            EXPA1=EXP(-WQIBOT(L)  /WQISD(L))
            WQF2ID=WQTT1/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0)
          
            EXPA0=EXP(-WQITOP(L,K)/WQISG(L))
            EXPA1=EXP(-WQIBOT(L)  /WQISG(L))
            WQF2IG=WQTT1/(DZCHP(L)*WQKESS)*(EXPA1-EXPA0)

          ELSE
            ! *** No Light Case
            WQIBOT(L)=0.
            WQITOP(L,K)=0.
            WQF2IC=0.0
            WQF2ID=0.0
            WQF2IG=0.0
            WQKESS = 0 
          ENDIF

          ! *** UPDATE SOLAR RADIATION AT BOTTOM OF THIS LAYER  !AA
          ! *** MACROALGAE SUBMODEL - COMPUTE IF HAVE SOLAR RADIATION
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) .AND. WQITOP(L,K) > 1.0E-18 )THEN  

            IZ=IWQZMAP(L,K)  
          
            ! *** Light Limitation
            WQFDI0 = - WQIBOT(L) / (WQFD + 1.E-18)
            WQISM = MAX( WQAVGIO*EXP(-WQKESS*WQDOPM(IZ)), WQISMIN )  ! *** Optimal Light
            WQFDM = WQFDI0 / (WQISM + 1.E-18)                        ! *** Ratio of Actual to Optimal Light
            WQHTT = WQHT(K) * HP(L)  
            WQTTB = EXP( -WQKESS * (WQHTT+1.0/DZWQ(L)) )  
            WQTTT = EXP( -WQKESS * WQHTT )  
            WQTT1 = (CNS1 * WQFD * DZWQ(L)) / WQKESS  
            WQF2IM = WQTT1 * (EXP(WQFDM*WQTTB) - EXP(WQFDM*WQTTT))   ! *** Light Based Macroalgae Growth Limiting Factor
          
            ! *** Velocity Limitation
            UMRM = 0.5*( U(L,K) + U(LEC(L),K) )  
            VMRM = 0.5*( V(L,K) + V(LNC(L),K) )  
            WQVEL = SQRT(UMRM*UMRM + VMRM*VMRM)  
            WQLVF = 1.0  
      
            !C OPTION 1 FOR VELOCITY LIMITATION ASSUMES MACROALGAE GROWTH  
            !C IS LIMITED AT LOW VELOCITIES DUE TO REDUCED AVAILABILITY OF  
            !C NUTRIENTS REACHING THE ALGAE BIOMASS.  USES A MICHAELIS-MENTON  
            !C TYPE OF EQUATION.  
            IF( IWQVLIM  ==  1 )THEN  
              IF( WQVEL  > WQKMVMIN(L) )THEN  
                WQLVF = WQVEL / (WQKMV(L) + WQVEL)  
              ELSE  
                WQLVF = WQKMVMIN(L) / (WQKMV(L) + WQKMVMIN(L))  
              ENDIF  
            ENDIF  
      
            !C OPTION 2 FOR VELOCITY LIMITATION APPLIES A FIVE-PARAMETER LOGISTIC  
            !C FUNCTION THAT CAN BE ADJUSTED TO LIMIT MACROALGAE GROWTH FOR  
            !C EITHER LOW OR HIGH (SCOUR) VELOCITIES.  IN STREAMS WITH LOW NUTRIENTS,  
            !C THE LOW VELOCITY WILL LIKELY BE LIMITING SINCE AMPLE NUTRIENTS MAY  
            !C NOT REACH THE ALGAE BIOMASS DUE TO REDUCED FLOW.  IN STREAMS WITH  
            !C ABUNDANT NUTRIENTS, LOW VELOCITIES WILL NOT LIMIT MACROALGAE GROWTH,  
            !C INSTEAD, HIGH VELOCITIES WILL LIKELY SCOUR THE MACROALGAE AND DETACH  
            !C IT FROM THE SUBSTRATE.  
            IF( IWQVLIM  == 2 )THEN  
              XNUMER = WQKMVA(L) - WQKMVD(L)  
              XDENOM = 1.0 + (WQVEL/WQKMVC(L))**WQKMVB(L)  
              WQLVF = WQKMVD(L) + ( XNUMER / (XDENOM**WQKMVE(L)) )  
            ENDIF  
      
            ! *** USE THE MORE SEVERELY LIMITING OF VELOCITY OR NUTRIENT FACTORS:  
            WQF1NM = MIN(WQLVF, WQF1NM)  
      
            ! *** FIRST CONVERT FROM MACROALGAE FROM A CONCENTRATION: WQV (MG C/M3) TO A DENSITY: XMRM (MG C/M2).  
            XMRM = WQV(L,K,IDNOTRVA)*DZCHP(L)  
            WQLDF = WQKBP(L) / (WQKBP(L) + XMRM)  
            WQPM(L)= WQPMM(IMWQZT(L))*WQF1NM*WQF2IM*WQTDGM(IWQT(L))*WQLDF   ! *** Final Macroalgae Growth Limiting Factor
          ELSE
            WQPM(L) = 0.
          ENDIF  
      
          ! *** Compute the Growth Rate based on Maximums & Limiting Factors
          IF( IWQSTOX == 1 )THEN  
            WQF4SC = WQSTOX / (WQSTOX + SWQ(L)*SWQ(L)+1.E-12)  
            WQPC(L)=WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L))*WQF4SC  
          ELSE  
            WQPC(L) = WQPMC(IMWQZT(L))*WQF1NC*WQF2IC*WQTDGC(IWQT(L))  
          ENDIF  
          WQPD(L) = WQPMD(IMWQZT(L))*WQF1ND*WQF2ID*WQTDGD(IWQT(L))  
          WQPG(L) = WQPMG(IMWQZT(L))*WQF1NG*WQF2IG*WQTDGG(IWQT(L))  
          
          ! ALGAL BASAL METABOLISM & PREDATION  
          WQBMC(L) = WQBMRC(IMWQZT(L)) * WQTDRC(IWQT(L))  
          WQPRC(L) = WQPRRC(IMWQZT(L)) * WQTDRC(IWQT(L))  

          ! THE VARIABLE WQTDGP ADJUSTS PREDATION AND BASAL METABOLISM BASED ON A  
          ! LOWER/UPPER OPTIMUM TEMPERATURE FUNCTION.  THIS WILL ALLOW DIATOMS TO  
          ! BLOOM IN WINTER IF WQTDGP IS CLOSE TO ZERO.  
          WQBMD(L) = WQBMRD(IMWQZT(L))*WQTDRD(IWQT(L))*WQTDGP(IWQT(L))  
          WQPRD(L) = WQPRRD(IMWQZT(L))*WQTDRD(IWQT(L))*WQTDGP(IWQT(L)) 
          
          WQBMG(L) = WQBMRG(IMWQZT(L)) * WQTDRG(IWQT(L))  
          WQPRG(L) = WQPRRG(IMWQZT(L)) * WQTDRG(IWQT(L))  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQBMM(L) = WQBMRM(IMWQZT(L)) * WQTDRM(IWQT(L))  
            WQPRM(L) = WQPRRM(IMWQZT(L)) * WQTDRM(IWQT(L))  
          ENDIF  

        ELSE
          ! *** DRY CELL BYPASS
          WQPC(L) = 0.
          WQPD(L) = 0.
          WQPG(L) = 0.
          WQBMC(L) = 0.
          WQPRC(L) = 0.
          WQBMD(L) = 0.
          WQPRD(L) = 0.
          WQBMG(L) = 0.
          WQPRG(L) = 0.
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQPM(L) = 0.  
            WQBMM(L) = 0.
            WQPRM(L) = 0.
            WQPRM(L) = 0.
          ENDIF  
        ENDIF    ! *** DRY CELL CHECK

        ! *** DSI NOTE  PMC
        !IF(TIMEDAY > TDBG .AND. L == LDBG )THEN
        !  WRITE(11,'(F10.2,I5,5F10.1,2F10.4)')TIMEDAY,K, SOLSWRT(L), WQKESS, WQI0, WQITOP(L,K), WQIBOT(L), ICETHICK(L), WQF2IG
        !  IF( K == 1 )TDBG = TDBG+1.0
        !ENDIF

      ENDDO      ! *** END HORIZONTAL LOOP FOR ALGAE PARMETERS

      ! *** 
      XMRM = 0.0  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IZ=IWQZMAP(L,K)  
        WQOBTOT(L) = WQV(L,K,1) + WQV(L,K,2) + WQV(L,K,3)  
        WQKRPC(L) = (WQKRC + WQKRCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))  
        WQKLPC(L) = (WQKLC + WQKLCALG*WQOBTOT(L)) * WQTDHDR(IWQT(L))  
        IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
          XMRM = WQKDCALM(IZ) * WQV(L,K,IDNOTRVA)  
        ENDIF  

        ! M. MORTON 08/28/99: ADDED SPATIALLY VARIABLE DOC HYDROLYSIS RATE WQKDC  
        !    TO ACHIEVE BETTER CONTROL IN SYSTEMS WITH A COMBINATION OF FRESHWAT  
        !    STREAMS AND TIDAL RIVERS WITH DIFFERENT CHARACTERISTICS.  

        WQKDOC = (WQKDC(IZ) + WQKDCALG*WQOBTOT(L) + XMRM)*WQTDMNL(IWQT(L))  
        O2WQ(L) = MAX(WQV(L,K,19), 0.0)  
        WQTT1 = WQKDOC / (WQKHORDO + O2WQ(L)+ 1.E-18)  
        WQKHR(L)   = WQTT1*O2WQ(L)  
        WQDENIT(L) = WQTT1*WQAANOX*RNO3WQ(L)/(WQKHDNN+RNO3WQ(L)+ 1.E-18)  
      ENDDO  

      ! ***********************************
      ! 7-10 PHOSPHORUS  
      ! *** HYDROLYSIS  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        WQAPC(L) = 1.0/(WQCP1PRM+WQCP2PRM*EXP(-WQCP3PRM*PO4DWQ(L)))  
        WQKHP = (WQKHPC+WQKHPD+WQKHPG) / 3.0  
        WQTT1 = WQKHP / (WQKHP+PO4DWQ(L)+ 1.E-18) * WQOBTOT(L)  
        WQKRPP(L) = (WQKRP + WQKRPALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** RPOP--> PO4
        WQKLPP(L) = (WQKLP + WQKLPALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** LPOP--> DOP
        WQKDOP(L) = (WQKDP + WQKDPALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** DOP --> PO4
      ENDDO
      
      ! *** PHOSPHATE SETTLING   
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IF( IWQSRP == 1 )THEN  
          WQTTM = WQKPO4P*WQTAMP(L,K)  
          WQH10(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)  
          IF( K /= KC )THEN  
            WQTTM = WQKPO4P*WQTAMP(L,K+1)  
            WQT10(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)  
          ENDIF  
        ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN  
          WQTTS = WQKPO4P*SEDT(L,K)  
          WQH10(L) = - WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
          IF( K /= KC )THEN  
            WQTTS = WQKPO4P*SEDT(L,K)  
            WQT10(L) = WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
          ENDIF  
        ELSE  
          WQH10(L) = 0.0  
          WQT10(L) = 0.0  
        ENDIF 
        WQH10(L) = WQH10(L)*DTWQO2 
      ENDDO  

      ! ***********************************
      ! 11-15 NITROGEN  
      ! *** HYDROLYSIS  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        WQKHN = (WQKHNC+WQKHND+WQKHNG) / 3.0  
        WQTT1 = WQKHN / (WQKHN+RNH4NO3(L)+ 1.E-18) * WQOBTOT(L)  
        WQKRPN(L) = (WQKRN + WQKRNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** RPON-->NH3
        WQKLPN(L) = (WQKLN + WQKLNALG*WQTT1) * WQTDHDR(IWQT(L))  ! *** LPON-->DON
        WQKDON(L) = (WQKDN + WQKDNALG*WQTT1) * WQTDMNL(IWQT(L))  ! *** DON -->NH3
      ENDDO  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IF( RNH4NO3(L) == 0.0 )THEN  
          WQPNC(L)=0.0  
          WQPND(L)=0.0  
          WQPNG(L)=0.0  
          WQPNM(L)=0.0  
        ELSE  
          WQTTC = RNH4WQ(L)/(WQKHNC+RNO3WQ(L)+ 1.E-18)  
          WQTTD = RNH4WQ(L)/(WQKHND+RNO3WQ(L)+ 1.E-18)  
          WQTTG = RNH4WQ(L)/(WQKHNG+RNO3WQ(L)+ 1.E-18)  
          WQTTM = RNH4WQ(L)/(WQKHNM+RNO3WQ(L)+ 1.E-18)  
          WQPNC(L) = (RNO3WQ(L)/(WQKHNC+RNH4WQ(L)+ 1.E-18) + WQKHNC/(RNH4NO3(L)+ 1.E-18)) * WQTTC  
          WQPND(L) = (RNO3WQ(L)/(WQKHND+RNH4WQ(L)+ 1.E-18) + WQKHND/(RNH4NO3(L)+ 1.E-18)) * WQTTD  
          WQPNG(L) = (RNO3WQ(L)/(WQKHNG+RNH4WQ(L)+ 1.E-18) + WQKHNG/(RNH4NO3(L)+ 1.E-18)) * WQTTG  
          WQPNM(L) = (RNO3WQ(L)/(WQKHNM+RNH4WQ(L)+ 1.E-18) + WQKHNM/(RNH4NO3(L)+ 1.E-18)) * WQTTM  
        ENDIF  
        WQNIT(L) = WQTDNIT(IWQT(L)) * O2WQ(L) / (WQKHNDO + O2WQ(L) + 1.E-18) * RNH4WQ(L) / (WQKHNN + RNH4WQ(L) + 1.E-18)
      ENDDO  
      IF( IWQSI == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( IWQSRP == 1 )THEN  
            WQTTM = WQKSAP*WQTAMP(L,K)  
            WQN17(L) = - WQWSSET(L,1) * WQTTM / (1.0+WQTTM)  
            IF( K /= KC )THEN  
              WQTTM = WQKSAP*WQTAMP(L,K+1)  
              WQT17(L) = WQWSSET(L,2) * WQTTM / (1.0+WQTTM)  
            ENDIF  
          ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN  
            WQTTS = WQKSAP*SEDT(L,K)  
            WQN17(L) = - WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
            IF( K /= KC )THEN  
              WQTTS = WQKSAP*SEDT(L,K+1)  
              WQT17(L) = WSEDO(NS) * WQTTS * DZWQ(L) / (1.0+WQTTS)  
            ENDIF  
          ELSE  
            WQN17(L) = 0.0  
            WQT17(L) = 0.0  
          ENDIF  
        ENDDO  
        WQN17(L) = WQN17(L)*DTWQO2 
      ENDIF  

      ! ***********************************
      ! *** DISSOLVED OXYGEN
      PPCDO=-3.45  !PARTIAL PRES OF CO2 IN 10^ppcdo ATM; TEMPORARILY DECLARED HERE. SHLD BE READ IN FROM INPUT FILE
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IZ=IWQZMAP(L,K)  
        WQO18(L)= -DTWQO2*WQKCOD(IWQT(L),IZ)*O2WQ(L) / (WQKHCOD(IZ) + O2WQ(L) + 1.E-18)  

        ! *** DO Saturation, Modified by SCJ, see Garcia and Gordon, Limnology and Oceanography 37(6), 1992, Eqn. 8 and Table 1
        TVAL1=log((298.15-TWQ(L))/(273.15+TWQ(L)))
        TVAL2=TVAL1*TVAL1
        TVAL3=TVAL1*TVAL2
        TVAL4=TVAL1*TVAL3
        TVAL5=TVAL1*TVAL4
        RLNSAT1=5.80818+3.20684*TVAL1+4.11890*TVAL2+4.93845*TVAL3+1.01567*TVAL4+1.41575*TVAL5
        RLNSAT2=SWQ(L)*(-7.01211E-3-7.25958E-3*TVAL1-7.93334E-3*TVAL2-5.54491E-3*TVAL3)-1.32412E-7*SWQ(L)*SWQ(L)
        WQDOS(L) = EXP(RLNSAT1+RLNSAT2)*32E-3 !32E-3 approximately converts micromol/L to mg/L or g/m^3
        XDOSAT(L,K) = XDOSAT(L,K) + WQDOS(L)*DTWQ*DZCHP(L)
      
        !************* CO2 parameters
        !VB COMPUTING THE pK FOR SAT CONC OF CO2; K - HENRY'S CONST
        CDOSATIDX(L) = -2385.73/(TWQ(L) + 273.15) -  0.0152642 * (TWQ(L) + 273.15) + 14.0184
        !          K * MOL WT OF CO2 * PARTAL PRES OF CO2 IN ATM
        WQCDOS(L) = 10.**(-CDOSATIDX(L)+PPCDO) * (44.* 1000.) !VB EVALUATING CONC OF CO2 IN G/M^3 
        !************* CO2 parameters

        ! *** Compute Reaeration
        IF( K == KC )THEN 
       
          ! DO NOT ALLOW WIND SPEEDS ABOVE 11 M/SEC IN THE FOLLOWING EQUATION
          WINDREA = MIN(WINDST(L),11.)  
          WQWREA = 0.728*SQRT(WINDREA) + (0.0372*WINDREA-0.317)*WINDREA  
    
          IF( IWQKA(IZ)  ==  0 )THEN
            ! *** Constant  
            WQVREA = WQKRO(IZ)  
            WQWREA = 0.0  
          ELSEIF( IWQKA(IZ)  ==  1 )THEN  
            ! *** Constant plus Wind
            WQVREA = WQKRO(IZ)  
          ELSEIF( IWQKA(IZ)  ==  2 )THEN
            ! *** OCONNOR-DOBBINS REAERATION FORMULA
            UMRM = 0.5*(U(L,K)+U(LEC(L),K))  
            VMRM = 0.5*(V(L,K)+V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            ! *** WQKRO = 3.933 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.5 / HP(L)**1.5  
          ELSEIF( IWQKA(IZ)  ==  3 )THEN
            ! *** OWENS & GIBBS (1964) REAERATION FORMULA
            UMRM = MAX(U(L,K), U(LEC(L),K))  
            VMRM = MAX(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            ! *** WQKRO = 5.32 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.67 / HP(L)**1.85  
          ELSEIF( IWQKA(IZ)  ==  4 )THEN  
            ! *** MODIFIED OWENS AND GIBBS REAERATION EQUATION:  
            ! *** NOTE: NORMALIZED TO A DEPTH OF 1.0 FT, I.E., THIS EQUATION GIVES THE  
            ! ***       SAME REAERATION AS OWENS & GIBBS AT 1.0 FT DEPTH; AT HIGHER  
            ! ***       DEPTHS IT GIVES LARGER REAERATION THAN OWENS & GIBBS.  
            UMRM = MAX(U(L,K), U(LEC(L),K))  
            VMRM = MAX(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            YMRM = HP(L)*3.0*(1.0 - HP(L)/(HP(L)+0.1524))  
            ! *** WQKRO = 5.32 TYPICALLY  
            WQVREA = WQKRO(IZ) * XMRM**0.67 / YMRM**1.85  
          ELSEIF( IWQKA(IZ)  ==  5 )THEN  
            UMRM = MAX(U(L,K), U(LEC(L),K))  
            VMRM = MAX(V(L,K), V(LNC(L),K))  
            XMRM = SQRT(UMRM*UMRM + VMRM*VMRM)  
            WQVREA = 3.7*XMRM  
          ENDIF  

          ! *** NOW COMBINE REAERATION DUE TO WATER VELOCITY AND WIND STRESS
          WQVREA = WQVREA * REAC(IZ)  ! *** Diffusive Flux
          WQWREA = WQWREA * REAC(IZ)  ! *** Wind Component
          WQP19(L) = - (WQVREA + WQWREA) * DZWQ(L)* WQTDKR(IWQT(L),IZ)  
          WQKRDOS(L) = -WQP19(L)*WQDOS(L)
          WQP22(L) = WQP19(L)*((32./44.)**0.25)     !VB Kr FOR CO2 ANALOGOUS TO WQP19 ; 44 = MOL WT OF CO2
          WQKRCDOS(L) = -WQP22(L) * WQCDOS(L)       !VB EVALUATING Kr*SAT CONC OF CO2
        ELSE
          WQKRDOS(L)  = 0.0
          WQKRCDOS(L) = 0.0  
          WQP19(L)    = 0.0  
          WQP22(L)    = 0.0              !VB Kr FOR CO2 IS ZERO FOR CELLS NOT AT THE SURFACE
        ENDIF  
      ENDDO  
    
      ! *** ICE
      IF( ISICE > 0 .AND. K == KC )THEN
        ! *** REDUCE REAERATION AND SUFACE EXCHANGE DUE TO ICE
        DO LP=LF,LL
          L=LWET(LP)
          WQP19(L)   = (1.-ICECOVER(L))*WQP19(L)
          WQKRDOS(L) = (1.-ICECOVER(L))*WQKRDOS(L)
        ENDDO
      ENDIF
        
      ! ***********************************
      IF( IWQSRP == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQR20(L) = WQWPSL(L,K,20)*VOLWQ(L) + (WQV(L,K,20) - WQTAMP(L,K)) * WQWSSET(L,1)  
        ENDDO

        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( K == KSZ(L) ) WQR20(L) = WQR20(L) + WQTDTAM(IWQT(L))*DZWQ(L)/(WQKHBMF+O2WQ(L)+ 1.E-18)
        ENDDO

        IF( K /= KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQR20(L) = WQR20(L) + (WQV(L,K+1,20) - WQTAMP(L,K+1)) * WQWSSET(L,2)  
          ENDDO 
        ELSE     ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQR20(L)=WQR20(L)+(WQWDSL(L,KC,20)+WQATML(L,KC,20))*VOLWQ(L)
          ENDDO
        ENDIF 
      ENDIF  
    
      ! ***********************************
      ! *** MACROPHYTE
      IF( IDNOTRVA > 0 )THEN
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            WQA1C = (WQPM(L)-WQBMM(L)-WQPRM(L)-WQWSM*DZWQ(L)) * DTWQO2  
            WQVA1C = 1.0 / (1.0 - WQA1C)
            WQV(L,K,IDNOTRVA) = (WQV(L,K,IDNOTRVA)+WQA1C*WQV(L,K,IDNOTRVA))*WQVA1C*SMAC(L)  
            WQV(L,K,IDNOTRVA) = MAX(WQV(L,K,IDNOTRVA),WQMCMIN)*SMAC(L)  
            WQO(L,IDNOTRVA) = WQVO(L,K,IDNOTRVA)+WQV(L,K,IDNOTRVA)  
          ENDIF
        ENDDO
      ENDIF  
  
      !******************************************************************************
      ! ***
      ! *** NOW COMPUTE KINETICS FOR EACH CONSTITUENT

      ! ****  PARAM 01  CHC - cyanobacteria
      IF( ISTRWQ(1) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
          WQA1C=(WQPC(L)-WQBMC(L)-WQPRC(L)-WQBCSET(L,1))*DTWQO2 !production per unit time multiplied by half time step
          WQKK(L) = 1.0 / (1.0 - WQA1C) 

          ! ***   PT_SRC_LOADS    VOLUME  
          WQR1C = WQWPSL(L,K,1) * VOLWQ(L)  !point source load rate multiplied by inverse cell volume  g/m^3/t
          WQRR(L) = WQV(L,K,1) + DTWQ*WQR1C + WQA1C*WQV(L,K,1)   !transported biomass conc. (CALWQC) + point source load rate X time step + growth rate X previous biomass conc.
        ENDDO  

        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQBCSET(L,2)*WQO(L,1) !biomass conc. + DtX(1/t)* biomass conc.
          ENDDO
        ELSE  ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***      ATM DRY DEP      ATM WET DEP    VOLUME  
            WQR1C = (WQWDSL(L,KC,1) + WQATML(L,KC,1))*VOLWQ(L)  !atmospheric loading mass per time / cell volume
            WQRR(L) = WQRR(L) + DTWQ*WQR1C  !biomass conc. + Dt*loading rate per unit volume
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,1)=SCB(L)*(WQRR(L)*WQKK(L))+(1.0-SCB(L))*WQV(L,K,1)  !boundary condition implementation
          WQV(L,K,1) = WQRR(L)*WQKK(L)
          WQO(L,1)   = WQVO(L,K,1) + WQV(L,K,1) !depth totaled biomass conc = old biomass conc in cell + biomass conc from this iteration
        ENDDO  
      ENDIF  
  
      ! ****  PARAM 02  CHD - diatom algae
      IF( ISTRWQ(2) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
          WQA2D=(WQPD(L)-WQBMD(L)-WQPRD(L)-WQBDSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQA2D)  

          ! ***   PT_SRC_LOADS    VOLUME  
          WQR2D = WQWPSL(L,K,2) * VOLWQ(L)  
          WQRR(L) = WQV(L,K,2) + DTWQ*WQR2D + WQA2D*WQV(L,K,2)  
        ENDDO  

        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQBDSET(L,2)*WQO(L,2)
          ENDDO 
        ELSE  ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR2D = (WQWDSL(L,KC,2)+WQATML(L,KC,2))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR2D  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,2)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,2)  
          WQV(L,K,2) = WQRR(L)*WQKK(L)
          WQO(L,2)   = WQVO(L,K,2)+WQV(L,K,2)
        ENDDO  
      ENDIF  
  
      ! ****  PARAM 03  CHG - green algae
      IF( ISTRWQ(3) == 1 )THEN
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! *** GROWTH BASAL_METAB PREDATION SETTLING  TIME STEP  
          WQA3G=(WQPG(L)-WQBMG(L)-WQPRG(L)-WQBGSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQA3G)  

          ! ***   PT_SRC_LOADS    VOLUME  
          WQR3G = WQWPSL(L,K,3) * VOLWQ(L) 
          ! ***                   External      Internal 
          WQRR(L) = WQV(L,K,3) + DTWQ*WQR3G + WQA3G*WQV(L,K,3)  
        ENDDO
        
        IF( K /= KC )THEN
          ! *** Add the Algae settled in from the cell above  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQBGSET(L,2)*WQO(L,3)
          ENDDO  
        ELSE  ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP     VOLUME  
            WQR3G = (WQWDSL(L,KC,3)+WQATML(L,KC,3))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR3G  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,3)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,3)  
          WQV(L,K,3) = WQRR(L)*WQKK(L)
          WQO(L,3)   = WQVO(L,K,3)+WQV(L,K,3)
        ENDDO  
      ENDIF  

      ! ****  PARAM 04  ROC - refractory particulate organic carbon
      IF( ISTRWQ(4) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS    SETTLING  
          WQB4 = -( WQKRPC(L) + WQRPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQB4)  

          ! ***  ALGAE PREDATION SOURCE OF RPOC  
          WQA4 = WQFCRP * (WQPRC(L)*WQO(L,1) + WQPRD(L)*WQO(L,2) + WQPRG(L)*WQO(L,3))  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA4 = WQA4 + WQFCRPM*WQPRM(L)*WQVO(L,K,IDNOTRVA)  
          ENDIF  

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR4 = WQWPSL(L,K,4) * VOLWQ(L)  
          WQRR(L) = WQV(L,K,4) + DTWQ*WQR4 + DTWQO2*WQA4 + WQB4*WQV(L,K,4)  
        ENDDO

        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQO(L,4)  
          ENDDO  
        ELSE  ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR4 = (WQWDSL(L,KC,4)+WQATML(L,KC,4))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR4  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,4)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,4)  
          WQV(L,K,4) = WQRR(L)*WQKK(L)
          WQO(L,4)   = WQVO(L,K,4)+WQV(L,K,4)
        ENDDO  
      ENDIF  

      ! ****  PARAM 05  LOC - labile particulate organic carbon
      IF( ISTRWQ(5) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS    SETTLING  
          WQC5 = - (WQKLPC(L)  + WQLPSET(L,1))*DTWQO2 
          WQKK(L) = 1.0 / (1.0 - WQC5)
                                                !Predation    
          WQA5 = WQFCLP * (WQPRC(L)*WQO(L,1) + WQPRD(L)*WQO(L,2) + WQPRG(L)*WQO(L,3))    

          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA5 = WQA5 + WQFCLPM*WQPRM(L)*WQVO(L,K,IDNOTRVA)  
          ENDIF  

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR5 = WQWPSL(L,K,5) * VOLWQ(L)

          WQRR(L) = WQV(L,K,5) + DTWQ*WQR5 + DTWQO2*WQA5 + WQC5*WQV(L,K,5)
        ENDDO  
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQO(L,5)  
          ENDDO  
        ELSE  ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR5 = (WQWDSL(L,K,5)+WQATML(L,KC,5))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR5  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,5)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,5)  
          WQV(L,K,5) = WQRR(L)*WQKK(L)
          WQO(L,5)   = WQVO(L,K,5)+WQV(L,K,5)
        ENDDO  
      ENDIF  

      ! ****  PARAM 06  DOC - dissolved organic carbon
      IF( ISTRWQ(6) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***    RESPIRATION  DENITRIFICATION
          WQD6 = - ( WQKHR(L) +   WQDENIT(L)) *DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQD6)  
          WQA6C=WQFCDC + CFCDCWQ*( WQKHRC/(WQKHRC+O2WQ(L)+ 1.E-18) )  
          WQA6D=WQFCDD + CFCDDWQ*( WQKHRD/(WQKHRD+O2WQ(L)+ 1.E-18) )  
          WQA6G=WQFCDG + CFCDGWQ*( WQKHRG/(WQKHRG+O2WQ(L)+ 1.E-18) )  
          WQA6 = ( WQA6C*WQBMC(L) + WQFCDP*WQPRC(L) )*WQO(L,1) + ( WQA6D*WQBMD(L) + WQFCDP*WQPRD(L) )*WQO(L,2) + ( WQA6G*WQBMG(L) + WQFCDP*WQPRG(L) )*WQO(L,3)  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            IZ=IWQZMAP(L,K)  
            WQA6M = (WQFCDM+(1.-WQFCDM)*WQKHRM(IZ) / (WQKHRM(IZ) + O2WQ(L) + 1.E-18))*WQBMM(L)  
            WQA6 = WQA6 + (WQA6M+ WQFCDPM*WQPRM(L))*WQVO(L,K,IDNOTRVA)  
          ENDIF  

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR6 = WQWPSL(L,K,6) * VOLWQ(L)

          WQRR(L) = WQV(L,K,6) + DTWQ*WQR6 + WQD6*WQV(L,K,6) + DTWQO2*(WQA6 +WQKLPC(L)*WQO(L,5))
        ENDDO

        IF( K == KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR6 = (WQWDSL(L,K,6)+WQATML(L,KC,6))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR6  
          ENDDO
        ENDIF
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,6)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,6)  
          WQV(L,K,6) = WQRR(L)*WQKK(L)
          WQO(L,6)   = WQVO(L,K,6)+WQV(L,K,6)
        ENDDO  
      ENDIF  

      ! ****  PARAM 07  ROP - refractory particulate organic phosphorus
      IF( ISTRWQ(7) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQE7 = - (WQKRPP(L)+WQRPSET(L,1))*DTWQO2 
          WQKK(L) = 1.0 / (1.0 - WQE7)  
          WQA7C = (WQFPRC*WQBMC(L) + WQFPRP*WQPRC(L)) * WQO(L,1)  
          WQA7D = (WQFPRD*WQBMD(L) + WQFPRP*WQPRD(L)) * WQO(L,2)  
          WQA7G = (WQFPRG*WQBMG(L) + WQFPRP*WQPRG(L)) * WQO(L,3)  
          WQA7 = (WQA7C+WQA7D+WQA7G) * WQAPC(L)  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA7 = WQA7 + (WQFPRM*WQBMM(L) + WQFPRPM*WQPRM(L))* WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM  
          ENDIF  

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR7 = WQWPSL(L,K,7) * VOLWQ(L)  
          WQRR(L) = WQV(L,K,7) + DTWQ*WQR7 + DTWQO2*WQA7 + WQE7*WQV(L,K,7)   
        ENDDO  
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQO(L,7)
          ENDDO  
        ELSE
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR7 = (WQWDSL(L,K,7)+WQATML(L,KC,7))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR7  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,7)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,7)  
          WQV(L,K,7) = WQRR(L)*WQKK(L)
          WQO(L,7)   = WQVO(L,K,7)+WQV(L,K,7)
        ENDDO  
      ENDIF  

      ! ****  PARAM 08  LOP - labile particulate organic phosphorus
      IF( ISTRWQ(8) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***    HYDROLYSIS  SETTLING
          WQF8 = - (WQKLPP(L)+WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQF8)  
          WQA8C = (WQFPLC*WQBMC(L) + WQFPLP*WQPRC(L)) * WQO(L,1)  
          WQA8D = (WQFPLD*WQBMD(L) + WQFPLP*WQPRD(L)) * WQO(L,2)  
          WQA8G = (WQFPLG*WQBMG(L) + WQFPLP*WQPRG(L)) * WQO(L,3)  
          WQA8 = (WQA8C+WQA8D+WQA8G) * WQAPC(L)  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA8 = WQA8 + (WQFPLM*WQBMM(L) + WQFPLPM*WQPRM(L))* WQVO(L,K,IDNOTRVA)* WQAPC(L)*WQAPCM  
          ENDIF  

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR8 = WQWPSL(L,K,8) * VOLWQ(L)  
          WQRR(L) = WQV(L,K,8) + DTWQ*WQR8 + DTWQO2*WQA8 + WQF8*WQV(L,K,8)  
        ENDDO  
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQO(L,8)
          ENDDO  
        ELSE
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP   ATM WET DEP    VOLUME  
            WQR8 = (WQWDSL(L,K,8)+WQATML(L,KC,8))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR8  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,8)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,8)  
          WQV(L,K,8) = WQRR(L)*WQKK(L)
          WQO(L,8)   = WQVO(L,K,8)+WQV(L,K,8)
        ENDDO  
      ENDIF  

      ! ****  PARAM 09  DOP - dissolved organic phosphorus
      IF( ISTRWQ(9) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQF9 = - DTWQO2*WQKDOP(L)  
          WQKK(L) = 1.0 / (1.0 - WQF9)
          WQA9C = (WQFPDC*WQBMC(L) + WQFPDP*WQPRC(L)) * WQO(L,1)  
          WQA9D = (WQFPDD*WQBMD(L) + WQFPDP*WQPRD(L)) * WQO(L,2)  
          WQA9G = (WQFPDG*WQBMG(L) + WQFPDP*WQPRG(L)) * WQO(L,3)  
          WQA9 = (WQA9C+WQA9D+WQA9G) * WQAPC(L)  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA9 = WQA9 + (WQFPDM*WQBMM(L) + WQFPDPM*WQPRM(L)) * WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM  
          ENDIF  

          ! ***  PT_SRC_LOADS    VOLUME  
          WQR9 = WQWPSL(L,K,9) * VOLWQ(L)  
          WQRR(L) = WQV(L,K,9) + DTWQ*WQR9 + WQF9*WQV(L,K,9) + DTWQO2*(WQA9 + WQKLPP(L)*WQO(L,8) )
        ENDDO

        IF( K == KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***    ATM DRY DEP    ATM WET DEP    VOLUME  
            WQR9 = (WQWDSL(L,KC,9)+WQATML(L,KC,9))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR9  
          ENDDO
        ENDIF
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,9)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,9)  
          WQV(L,K,9) = WQRR(L)*WQKK(L)
        ENDDO  
      ENDIF  

      ! ****  PARAM 10  P4D - total phosphate
      IF( ISTRWQ(10) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQA10C=(WQFPIC*WQBMC(L)+WQFPIP*WQPRC(L)-WQPC(L))*WQO(L,1)  
          WQA10D=(WQFPID*WQBMD(L)+WQFPIP*WQPRD(L)-WQPD(L))*WQO(L,2)  
          WQA10G=(WQFPIG*WQBMG(L)+WQFPIP*WQPRG(L)-WQPG(L))*WQO(L,3)  
          WQKK(L) = (WQA10C+WQA10D+WQA10G) * WQAPC(L)  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQKK(L) = WQKK(L) + (WQFPIM*WQBMM(L)+WQFPIP*WQPRM(L)-WQPM(L))*WQVO(L,K,IDNOTRVA) * WQAPC(L)*WQAPCM  
          ENDIF  

          ! ***    PT_SRC_LOADS    VOLUME  
          WQRR(L) = WQWPSL(L,K,10) * VOLWQ(L)  
        ENDDO  

        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
              WQRR(L) = WQRR(L) + WQBFPO4D(L)*DZWQ(L) ! *** Add in Benthic Flux
          ENDIF
        ENDDO  

        IF( K == KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***      ATM DRY DEP    ATM WET DEP     VOLUME  
            WQR10 = (WQWDSL(L,KC,10)+WQATML(L,KC,10))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR10  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQRR(L) = WQV(L,K,10) + DTWQ*WQRR(L) + WQH10(L)*WQV(L,K,10) + DTWQO2*(WQKK(L) + WQKRPP(L)*WQO(L,7) + WQKDOP(L)*WQO(L,9))
        ENDDO  
      
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQT10(L)*WQO(L,10)
          ENDDO  
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQKKL = 1.0 / (1.0 - WQH10(L)) 
          !WQV(L,K,10)=SCB(L)*(WQRR(L)*WQKKL)+(1.-SCB(L))*WQV(L,K,10)  
          WQV(L,K,10) = WQRR(L)*WQKKL
          WQO(L,10)   = WQVO(L,K,10)+WQV(L,K,10)
        ENDDO  
      ENDIF  

      ! ****  PARAM 11  RON - refractory particulate organic nitrogen
      IF( ISTRWQ(11) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS     SETTLING
          WQI11 = - (WQKRPN(L) + WQRPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQI11)  
          WQA11C=(WQFNRC*WQBMC(L)+WQFNRP*WQPRC(L))*WQANCC*WQO(L,1)  
          WQA11D=(WQFNRD*WQBMD(L)+WQFNRP*WQPRD(L))*WQANCD*WQO(L,2)  
          WQA11G=(WQFNRG*WQBMG(L)+WQFNRP*WQPRG(L))*WQANCG*WQO(L,3)  
          WQA11 = WQA11C+WQA11D+WQA11G  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA11 = WQA11 + (WQFNRM*WQBMM(L)+WQFNRPM*WQPRM(L))*WQANCM*WQVO(L,K,IDNOTRVA)  
          ENDIF  

          ! ***    PT_SRC_LOADS    VOLUME  
          WQR11 = WQWPSL(L,K,11) * VOLWQ(L)

          WQRR(L) = WQV(L,K,11) + DTWQ*WQR11 + DTWQO2*WQA11 + WQI11*WQV(L,K,11) 
        ENDDO 

        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQRPSET(L,2)*WQO(L,11)
          ENDDO  
        ELSE   ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
            WQR11 = (WQWDSL(L,KC,11)+WQATML(L,KC,11))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR11  
          ENDDO
        ENDIF
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,11)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,11)  
          WQV(L,K,11) = WQRR(L)*WQKK(L)
          WQO(L,11)   = WQVO(L,K,11)+WQV(L,K,11)
        ENDDO  
      ENDIF  

      ! ****  PARAM 12  LON - labile particulate organic nitrogen
      IF( ISTRWQ(12) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***     HYDROLYSIS     SETTLING
          WQJ12 = - (WQKLPN(L)+WQLPSET(L,1))*DTWQO2
          WQKK(L) = 1.0 / (1.0 - WQJ12)  
          WQA12C=(WQFNLC*WQBMC(L)+WQFNLP*WQPRC(L))*WQANCC*WQO(L,1)  
          WQA12D=(WQFNLD*WQBMD(L)+WQFNLP*WQPRD(L))*WQANCD*WQO(L,2)  
          WQA12G=(WQFNLG*WQBMG(L)+WQFNLP*WQPRG(L))*WQANCG*WQO(L,3)  
          WQA12 = WQA12C+WQA12D+WQA12G  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA12 = WQA12 + (WQFNLM*WQBMM(L)+WQFNLPM*WQPRM(L)) *WQANCM*WQVO(L,K,IDNOTRVA)  
          ENDIF  

          ! ***    PT_SRC_LOADS    VOLUME  
          WQR12 = WQWPSL(L,K,12) * VOLWQ(L)  

          WQRR(L) = WQV(L,K,12) + DTWQ*WQR12 + DTWQO2*WQA12 + WQJ12*WQV(L,K,12)
        ENDDO  
        IF( K /= KC )THEN  
          ! *** Add in settling from above
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQRR(L) + DTWQO2*WQLPSET(L,2)*WQO(L,12)
          ENDDO  
        ELSE   ! K == KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
            WQR12 = (WQWDSL(L,KC,12)+WQATML(L,KC,12))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR12  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,12)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,12)  
          WQV(L,K,12) = WQRR(L)*WQKK(L)
          WQO(L,12)   = WQVO(L,K,12)+WQV(L,K,12)
        ENDDO  
      ENDIF  

      ! ****  PARAM 13  DON - dissolved organic nitrogen
      IF( ISTRWQ(13) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQF13 = - DTWQO2*WQKDON(L)  
          WQKK(L) = 1.0 / (1.0 - WQF13)
          WQA13C=(WQFNDC*WQBMC(L)+WQFNDP*WQPRC(L))*WQANCC*WQO(L,1)  
          WQA13D=(WQFNDD*WQBMD(L)+WQFNDP*WQPRD(L))*WQANCD*WQO(L,2)  
          WQA13G=(WQFNDG*WQBMG(L)+WQFNDP*WQPRG(L))*WQANCG*WQO(L,3)  
          WQA13 = WQA13C+WQA13D+WQA13G  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA13 = WQA13 + (WQFNDM*WQBMM(L)+WQFNDPM*WQPRM(L)) *WQANCM*WQVO(L,K,IDNOTRVA)  
          ENDIF  

          ! ***    PT_SRC_LOADS    VOLUME  
          WQR13 = WQWPSL(L,K,13) * VOLWQ(L)

          WQRR(L) = WQV(L,K,13) + DTWQ*WQR13 + WQF13*WQV(L,K,13) + DTWQO2*( WQA13 + WQKLPN(L)*WQO(L,12))
        ENDDO

        IF( K == KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
            WQR13 = (WQWDSL(L,KC,13)+WQATML(L,KC,13))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + DTWQ*WQR13  
          ENDDO
        ENDIF  
      
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          !WQV(L,K,13)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,13)  
          WQV(L,K,13) = WQRR(L)*WQKK(L)
          WQO(L,13)   = WQVO(L,K,13)+WQV(L,K,13)
        ENDDO  
      ENDIF  

      ! ****  PARAM 14  NHX - ammonia nitrogen
      IF( ISTRWQ(14) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***      PT_SRC_LOADS    VOLUME  
          WQRR(L) = WQWPSL(L,K,14) * VOLWQ(L)  
        ENDDO  

        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            WQRR(L) = WQRR(L) + WQBFNH4(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF  
        ENDDO  

        IF( K == KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
            WQR14 = (WQWDSL(L,KC,14)+WQATML(L,KC,14))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR14  
          ENDDO
        ENDIF  

        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQF14 = - DTWQO2*WQNIT(L)  
          WQKK(L) = 1.0 / (1.0 - WQF14) 
          WQA14C=WQFNIC*WQBMC(L)+WQFNIP*WQPRC(L)-WQPNC(L)*WQPC(L)  
          WQA14D=WQFNID*WQBMD(L)+WQFNIP*WQPRD(L)-WQPND(L)*WQPD(L)  
          WQA14G=WQFNIG*WQBMG(L)+WQFNIP*WQPRG(L)-WQPNG(L)*WQPG(L)  
          WQA14 = WQA14C*WQANCC*WQO(L,1) + WQA14D*WQANCD*WQO(L,2) + WQA14G*WQANCG*WQO(L,3)  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA14 = WQA14 + (WQFNIM*WQBMM(L)+WQFNIPM*WQPRM(L)-WQPNM(L)*WQPM(L))*WQANCM*WQVO(L,K,IDNOTRVA)  
          ENDIF  
          WQRR(L) = WQV(L,K,14) + DTWQ*WQRR(L) + WQF14*WQV(L,K,14) + DTWQO2*( WQA14 + WQKRPN(L)*WQO(L,11) + WQKDON(L)*WQO(L,13) )                 
          !WQV(L,K,14)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,14)  
          WQV(L,K,14) = WQRR(L)*WQKK(L)
          WQO(L,14)   = WQVO(L,K,14)+WQV(L,K,14)
        ENDDO  
      ENDIF  

      ! ****  PARAM 15  NOX - nitrate nitrogen
      IF( ISTRWQ(15) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          ! ***      PT_SRC_LOADS    VOLUME  
          WQRR(L) = WQWPSL(L,K,15) * VOLWQ(L)  
        ENDDO  
        
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            WQRR(L) = WQRR(L) + WQBFNO3(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF  
        ENDDO  
        
        IF( K == KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
            WQR15 = (WQWDSL(L,KC,15)+WQATML(L,KC,15))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR15  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQA15C = (WQPNC(L)-1.0)*WQPC(L) * WQANCC * WQO(L,1)  
          WQA15D = (WQPND(L)-1.0)*WQPD(L) * WQANCD * WQO(L,2)  
          WQA15G = (WQPNG(L)-1.0)*WQPG(L) * WQANCG * WQO(L,3)  
          WQA15 = WQA15C + WQA15D + WQA15G  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            WQA15 = WQA15 + (WQPNM(L)-1.0)*WQPM(L)*WQANCM*WQVO(L,K,IDNOTRVA)  
          ENDIF  
          WQB15 = WQV(L,K,15) + DTWQ*WQRR(L) + DTWQO2*( WQA15 - WQANDC*WQDENIT(L)*WQO(L,6) + WQNIT(L)*WQO(L,14))

          !WQV(L,K,15)=SCB(L)*WQB15 + (1.-SCB(L))*WQV(L,K,15)  
          WQV(L,K,15) = WQB15
          WQO(L,15)   = WQVO(L,K,15)+WQV(L,K,15)
        ENDDO  
      ENDIF  

      ! ****  PARAM 16  SUU - particulate biogenic silica
      IF( ISTRWQ(16) == 1 )THEN  
        IF( IWQSI == 1 )THEN  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQM16 = - (WQKSUA(IWQT(L)) + WQBDSET(L,1)) * DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQM16)  
            WQA16D = (WQFSPD*WQBMD(L) + WQFSPP*WQPRD(L)) * WQASCD * WQO(L,2)  
            ! ***    PT_SRC_LOADS    VOLUME  
            WQR16 = WQWPSL(L,K,16) * VOLWQ(L)

            WQRR(L) = WQV(L,K,16) + DTWQ*WQR16 + DTWQO2*WQA16D + WQM16*WQV(L,K,16)  
          ENDDO  
          IF( K /= KC )THEN  
            ! *** Add in settling from above
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQBDSET(L,2)*WQO(L,16)     ! *** PMC
            ENDDO  
          ELSE
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR16 = (WQWDSL(L,KC,16)+WQATML(L,KC,16))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR16  
            ENDDO
          ENDIF  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            !WQV(L,K,16)=SCB(L)*( WQRR(L)*WQKK(L) ) + (1.-SCB(L))*WQV(L,K,16)  
            WQV(L,K,16) = WQRR(L)*WQKK(L)
            WQO(L,16)   = WQVO(L,K,16)+WQV(L,K,16)
          ENDDO  
        ENDIF  
      ENDIF  

      ! ****  PARAM 17  SAA - dissolved available silica
      IF( ISTRWQ(17) == 1 )THEN  
        IF( IWQSI == 1 )THEN  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQKK(L) = (WQFSID*WQBMD(L) + WQFSIP*WQPRD(L) - WQPD(L)) * WQASCD * WQO(L,2)  
            ! ***      PT_SRC_LOADS    VOLUME  
            WQRR(L) = WQWPSL(L,K,17) * VOLWQ(L)  
          ENDDO  
          
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            IF( K == KSZ(L) )THEN
              WQRR(L) = WQRR(L) + WQBFSAD(L)*DZWQ(L)   ! *** Add in Benthic Flux
            ENDIF  
          ENDDO  

          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQRR(L) = WQV(L,K,17) + DTWQ*WQRR(L) +WQN17(L)*WQV(L,K,17) + DTWQO2*( WQKK(L) + WQKSUA(IWQT(L))*WQO(L,16))  
          ENDDO  

          IF( K /= KC )THEN  
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQT17(L)*WQO(L,17)  
            ENDDO  
          ELSE
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
              WQR17 = (WQWDSL(L,KC,17)+WQATML(L,KC,17))*VOLWQ(L)  
              WQRR(L) = WQRR(L) + DTWQ*WQR17  
            ENDDO
          ENDIF  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - WQN17(L))
            !WQV(L,K,17)=SCB(L)*( WQRR(L)*WQKK(L) ) +(1.-SCB(L))*WQV(L,K,17)  
            WQV(L,K,17) = WQRR(L)*WQKK(L)
            WQO(L,17)   = WQVO(L,K,17)+WQV(L,K,17)
          ENDDO  
        ENDIF  
      ENDIF  

      ! ****  PARAM 18  COD - chemical oxygen demand
      IF( ISTRWQ(18) == 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQKK(L) = 1.0 / (1.0 - WQO18(L))  
            ! ***    PT_SRC_LOADS    VOLUME  
          WQRR(L) = WQWPSL(L,K,18) * VOLWQ(L)  
        ENDDO  
        
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            WQRR(L) = WQRR(L) + WQBFCOD(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF  
        ENDDO  

        IF( K == KC )THEN
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            ! ***     ATM DRY DEP     ATM WET DEP    VOLUME  
            WQR18 = (WQWDSL(L,KC,18)+WQATML(L,KC,18))*VOLWQ(L)  
            WQRR(L) = WQRR(L) + WQR18  
          ENDDO
        ENDIF  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQRR(L)=WQV(L,K,18)+DTWQ*WQRR(L)+WQO18(L)*WQV(L,K,18)  
          !WQV(L,K,18)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,18)  
          WQV(L,K,18) = WQRR(L)*WQKK(L)
          WQO(L,18)   = WQVO(L,K,18)+WQV(L,K,18)  
        ENDDO  
      ENDIF

      ! ****  PARAM 19  DOX - dissolved oxygen

      ! ***  1) CHC - cyanobacteria
      ! ***  2) CHD - diatom algae
      ! ***  3) CHG - green algae
      ! ***  4) ROC - refractory particulate organic carbon
      ! ***  5) LOC - labile particulate organic carbon
      ! ***  6) DOC - dissolved organic carbon
      ! ***  7) ROP - refractory particulate organic phosphorus
      ! ***  8) LOP - labile particulate organic phosphorus
      ! ***  9) DOP - dissolved organic phosphorus
      ! *** 10) P4D - total phosphate
      ! *** 11) RON - refractory particulate organic nitrogen 
      ! *** 12) LON - labile particulate organic nitrogen
      ! *** 13) DON - dissolved organic nitrogen
      ! *** 14) NHX - ammonia nitrogen
      ! *** 15) NOX - nitrate nitrogen
      ! *** 16) SUU - particulate biogenic silica
      ! *** 17) SAA - dissolved available silica
      ! *** 18) COD - chemical oxygen demand
      ! *** 19) DOX - dissolved oxygen
      ! *** 20) TAM - total active metal
      ! *** 21) FCB - fecal coliform bacteria
      ! *** 22) CO2 - dissolved carbon dioxide
      ! *** 23) MAC - macroalgae

      ! *** EE7.2 REMOVED THE DO COMPONENT ANALYSIS

      IF( ISTRWQ(19) == 1 )THEN 
        ! *** Point Source Loading 
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQRR(L) = WQWPSL(L,K,19) * VOLWQ(L)  
        ENDDO
      
        ! *** Handle Surface Processes
        IF( K == KC )THEN  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP19(L)) 
            ! ***                  ATM DRY DEP     ATM WET DEP    VOLUME  
            WQRR(L) = WQRR(L) + (WQWDSL(L,KC,19)+WQATML(L,KC,19))*VOLWQ(L)
            ! *** Reaeration - Atm to water flux  (Offset later in O2 Mass Balance by WQRea term)
            WQRR(L) = WQRR(L) + WQKRDOS(L)
          ENDDO
        ELSE  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQKK(L) = 1.0
          ENDDO
        ENDIF 

        ! *** Bottom Processes 
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( K == KSZ(L) )THEN
            TEMFAC=STEMFAC**(TEM(L,KSZ(L))-20.)
            WQRR(L) = WQRR(L) + TEMFAC*WQBFO2(L)*DZWQ(L)   ! *** Add in Benthic Flux
          ENDIF
        ENDDO
  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          DTWQxH = DTWQ*DZCHP(L)  
          DTWQxH2= DTWQO2*DZCHP(L)

          IF( WQI0  <=  0.001 )THEN  
            WQTTC = 0.0
            WQTTD = 0.0
            WQTTG = 0.0
          ELSE
            WQTTC = (1.3 - 0.3*WQPNC(L)) * WQPC(L) 
            WQTTD = (1.3 - 0.3*WQPND(L)) * WQPD(L)  
            WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)

            ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
            !TMP19 = WQAOCR*DTWQxH2*(WQTTC*WQO(L,1)+WQTTD*WQO(L,2)+WQTTG*WQO(L,3))
          ENDIF

          ! *** RESPIRATION OF TOTAL CHLOROPHYLL - CYANOBACTERIA 
          XMRM = CFCDCWQ*O2WQ(L)*WQBMC(L)/(WQKHRC+O2WQ(L)+ 1.E-18)  
          WQA19C = WQTTC - XMRM

          ! *** RESPIRATION OF TOTAL CHLOROPHYLL - DIATOMS 
          XMRM = CFCDDWQ*O2WQ(L)*WQBMD(L)/(WQKHRD+O2WQ(L)+ 1.E-18)  
          WQA19D = WQTTD - XMRM

          ! *** RESPIRATION OF TOTAL CHLOROPHYLL -  GREENS
          XMRM = CFCDGWQ*O2WQ(L)*WQBMG(L)/(WQKHRG+O2WQ(L)+ 1.E-18)  
          WQA19G = WQTTG - XMRM

          ! *** TOTAL NET RESPIRATION/PHOTOSYNTHESIS
          WQA19=(WQA19C*WQO(L,1) + WQA19D*WQO(L,2) + WQA19G*WQO(L,3)) * WQAOCR
        
          ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
          ! ***   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
          IF( IDNOTRVA > 0 .AND. K == KSZ(L) )THEN  
            ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
            WQO(L,IDNOTRVA) = WQVO(L,K,IDNOTRVA) + WQV(L,K,IDNOTRVA)
            IZ = IWQZMAP(L,K)  
            WQTTM = (1.3 - 0.3*WQPNM(L)) * WQPM(L)  
            XMRM=(1.0-WQFCDM)*O2WQ(L)*WQBMM(L)/(WQKHRM(IZ)+O2WQ(L)+1.E-18)  
            WQA19A = WQTTM*WQO(L,IDNOTRVA)*WQAOCRPM - XMRM*WQO(L,IDNOTRVA)*WQAOCRRM  
            WQA19 = WQA19 + WQA19A
          ENDIF  

          ! *** O2 Mass Balance
          ! WQA19                           ! *** Total Net Respiration/Photosynthesis
          WQSUM=DTWQ*WQRR(L)                ! *** Sum of Loadings/Demands
          WQRea=WQP19(L)*WQV(L,K,19)        ! *** Reaeration - Offsetting Flux (WQP19 is negative)
          WQPOC=WQAOCR*WQKRPC(L)*WQO(L,4)   ! *** POC
          WQDOC=WQAOCR*WQKHR(L) *WQO(L,6)   ! *** DOC
          WQNH3=WQAONT*WQNIT(L) *WQO(L,14)  ! *** Ammonia
          WQCOD=WQO18(L)*WQO(L,18)          ! *** COD
          WQRR(L) = WQV(L,K,19) + WQSUM  + WQCOD  + DTWQO2*(WQA19 - WQPOC - WQDOC - WQNH3 + WQRea)

          !WQV(L,K,19)=SCB(L)*(WQRR(L)*WQKK(L))+(1.-SCB(L))*WQV(L,K,19)
          WQV(L,K,19) = WQRR(L)*WQKK(L)
          WQV(L,K,19) = MAX(WQV(L,K,19), 0.0)  
        ENDDO  
      
      ENDIF  ! *** END OF ISTRWQ(19)

      ! ****  PARAM 20  TAM - total active metal
      IF( ISTRWQ(20) == 1 )THEN  
        IF( IWQSRP == 1 )THEN  
          DO LP=1,LLWET(K,ND)   
            L=LKWET(LP,K,ND)  
            WQT20 = - DTWQ*WQWSSET(L,1)    ! *** DTWQO2
            WQKK(L) = 1.0 / (1.0 - WQT20)  
            WQRR(L)=WQV(L,K,20)+DTWQ*WQR20(L)+WQT20*WQV(L,K,20)  
          ENDDO  
          IF( K /= KC )THEN  
            ! *** Add in settling from above
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              WQRR(L) = WQRR(L) + DTWQO2*WQWSSET(L,2)*WQO(L,20)
            ENDDO  
          ENDIF  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQV(L,K,20) = WQRR(L)*WQKK(L)
            WQO(L,20)   = WQVO(L,K,20)+WQV(L,K,20)
          ENDDO  
        ENDIF  
      ENDIF  

      ! ****  PARAM 21  FCB - fecal coliform bacteria
      IF( ISTRWQ(21) == 1 )THEN  
        IF( IWQFCB == 1 )THEN  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQKK(L) = WQTD2FCB(IWQT(L))  
    
            ! ***      ATM DRY DEP       LOADS          VOLUME  
            WQR21= (WQWDSL(L,K,21)+WQWPSL(L,K,21))*VOLWQ(L)      !VB CHANGED NWQV TO 21
            ! WQR21= (WQWDSL(L,K,NWQV)+WQWPSL(L,K,NWQV))*VOLWQ(L)  
    
            IF( K == KC .AND. LMASKDRY(L) )THEN
              ! ***                   ATM WET DEP      VOLUME  
              WQR21 = WQR21 + (WQATML(L,KC,21)*VOLWQ(L))  
            ENDIF
          
            ! WQRR(L) = WQV(L,K,NWQV)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21    !VB CHANGED NWQV TO 21
            WQRR(L) = WQV(L,K,21)*WQTD1FCB(IWQT(L)) + DTWQ*WQR21
            !WQV(L,K,21)=SCB(L)*( WQRR(L)*WQKK(L) ) + (1.-SCB(L))*WQV(L,K,21)  
            WQV(L,K,21) = WQRR(L)*WQKK(L)
          ENDDO  
        ENDIF  
      ENDIF

      ! ****  PARAM 22 DISSOLVED CARBON DIOXIDE

      !C THE FOLLOWING ARRAYS WERE ADDED TO KEEP TRACK OF THE VARIOUS COMPONENT  
      !C OF DISSOLVED CARBON DIOXIDE.  
      !C THE ARRAY DESCRIPTIONS ARE:  
      !C  XCDOKAR(L,K) = CDO. COMPONENT FOR REAERATION  
      !C  XCDODOC(L,K) = CDO. COMPONENT FOR DISS. ORG. CARBON DECAY  
      !C  XCDOPPB(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
      !C  XCDORRB(L,K) = CDO. COMPONENT FOR RESPIRATION OF TOTAL CHLOROPHYLL  
      !C  XCDOPPM(L,K) = CDO. COMPONENT FOR PHOTOSYNTHESIS OF MACROALGAE  
      !C  XCDORRM(L,K) = CDO. COMPONENT FOR RESPIRATION OF MACROALGAE  
      !C  XCDOALL(L,K) = SUM OF THE ABOVE 6 CDO. COMPONENTS  

      IF( ISTRWQ(22) == 1 )THEN  
        ! *** CO2
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQRR(L) = WQWPSL(L,K,22) * VOLWQ(L)  
          TMP22=WQRR(L)*DTWQ*DZCHP(L) 
        ENDDO  
      
        ! *** Handle Surface Processes
        IF( K == KC )THEN  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQKK(L) = 1.0 / (1.0 - DTWQO2*WQP22(L)) 
            ! ***             ATM DRY DEP    ATM WET DEP    VOLUME  
            WQRR(L)=WQRR(L)+(WQWDSL(L,KC,22)+WQATML(L,KC,22)*VOLWQ(L))  
            ! *** Reaeration
            WQRR(L) = WQRR(L) + WQKRCDOS(L) 
          ENDDO
        ELSE
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            WQKK(L) = 1.0
          ENDDO
        ENDIF   
      
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          DTWQxH = DTWQ*DZCHP(L)  
          DTWQxH2= DTWQO2*DZCHP(L)
          IF( WQI0  <=  0.001 )THEN  
            WQTTC = 0.0
            WQTTD = 0.0
            WQTTG = 0.0
          ELSE
            ! *** PHOTOSYNTHESIS OF TOTAL CHLOROPHYLL  
            WQTTC = (1.3 - 0.3*WQPNC(L)) * WQPC(L) 
            WQTTD = (1.3 - 0.3*WQPND(L)) * WQPD(L)  
            WQTTG = (1.3 - 0.3*WQPNG(L)) * WQPG(L)
          ENDIF

          ! *** RESPIRATION OF TOTAL CHLOROPHYLL - CYANOBACTERIA 
          XMRM = CFCDCWQ*O2WQ(L)*WQBMC(L)/(WQKHRC+O2WQ(L)+ 1.E-18)  
          WQA22C = WQTTC - XMRM

          ! *** RESPIRATION OF TOTAL CHLOROPHYLL - DIATOMS 
          XMRM = CFCDDWQ*O2WQ(L)*WQBMD(L)/(WQKHRD+O2WQ(L)+ 1.E-18)  
          WQA22D = WQTTD - XMRM

          ! *** RESPIRATION OF TOTAL CHLOROPHYLL -  GREENS
          XMRM = CFCDGWQ*O2WQ(L)*WQBMG(L)/(WQKHRG+O2WQ(L)+ 1.E-18)  
          WQA22G = WQTTG - XMRM
        
          ! *** TOTAL NET RESPIRATION/PHOTOSYNTHESIS
          WQA22=3.67*(WQA22C*WQO(L,1)+WQA22D*WQO(L,2)+WQA22G*WQO(L,3))  !VB 3.67 CONVERTS g CARBON TO g CO2

          ! *** MODIFIED BY MRM 05/23/99 TO ALLOW DIFFERENT AOCR CONSTANTS TO BE APPLIED  
          ! ***   TO PHOTOSYNTHESIS AND RESPIRATION TERMS FOR MACROALGAE  
          ! *** TRAPEZOIDAL AVERAGE CONCENTRATIONS
          ! *** CO2 Mass Balance
          ! WQA22                                          ! *** Total Net Respiration/Photosynthesis
          WQCDSUM=DTWQ*WQRR(L)                             ! *** Sum of Loadings/Demands
          WQCDRea=WQP22(L)*WQV(L,K,22)                     ! *** Reaeration
          WQCDDOC=(WQKHR(L)+WQDENIT(L))*WQO(L,6)*3.67      ! *** DOC FROM HYDROLYSIS AND DENITRIFICATION 3.67 CONVERTS G CARBON TO G CO2    

          WQRR(L) = WQV(L,K,22) + WQCDSUM + DTWQO2*(-WQA22 + WQCDDOC + WQCDRea)
          WQV(L,K,22) = WQRR(L)*WQKK(L)
       
          ! *** WQV(L,K,22) After this WQV(L,K,22) can not be < 0.
          WQV(L,K,22) = MAX(WQV(L,K,22), 0.0)  
        ENDDO  
      ENDIF  
    ENDDO  ! *** END OF THE KC LOOP

    ! *** PHOSPHATE AND SILICA SORPTION OPTIONS  
    IF( IWQSRP == 1 )THEN  
      ! *** Sorption Option: TAM
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          O2WQ(L) = MAX(WQV(L,K,19), 0.0)  
          WQTAMD = MIN( WQTAMDMX*EXP(-WQKDOTAM*O2WQ(L)), WQV(L,K,20) )  
          WQTAMP(L,K) = WQV(L,K,20) - WQTAMD  
          WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*WQTAMP(L,K))  
          WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*WQTAMP(L,K))  
        ENDDO  
      ENDDO  
    ELSE IF( IWQSRP == 2 .AND. ISTRAN(6) > 0 )THEN
      ! *** Sorption Option: Sediments
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQPO4D(L,K) = WQV(L,K,10) / (1.0 + WQKPO4P*SEDT(L,K))  
          WQSAD(L,K)  = WQV(L,K,17) / (1.0 + WQKSAP*SEDT(L,K))  
        ENDDO  
      ENDDO  
    ELSE  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          WQPO4D(L,K) = WQV(L,K,10)  
          WQSAD(L,K)  = WQV(L,K,17)  
        ENDDO  
      ENDDO  
    ENDIF  

    ! ***************************************************************************
    ! COUPLING TO SEDIMENT MODEL  
    ! EVALUATE DEP. FLUX USING NEW VALUES CAUSE IMPLICIT SCHEME IS USED IN SPM
    IF( IWQBEN == 1 )THEN  
      DO LP=LF,LL
        L=LWET(LP)
        IMWQZ = IWQZMAP(L,KSZ(L))  
        WQDFBC(L) = SCB(L)*WQWSC(IMWQZ)*WQV(L,KSZ(L),1)  
        WQDFBD(L) = SCB(L)*WQWSD(IMWQZ)*WQV(L,KSZ(L),2)  
        WQDFBG(L) = SCB(L)*WQWSG(IMWQZ)*WQV(L,KSZ(L),3) + WQWSM*DZWQ(L)*WQV(L,KSZ(L),IDNOTRVA)  
        WQDFRC(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,KSZ(L),4)  
        WQDFLC(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,KSZ(L),5)  
        WQDFRP(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,KSZ(L),7)  
        WQDFLP(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,KSZ(L),8)  
        WQDFRN(L) = SCB(L)*WQWSRP(IMWQZ)*WQV(L,KSZ(L),11)  
        WQDFLN(L) = SCB(L)*WQWSLP(IMWQZ)*WQV(L,KSZ(L),12)  
        IF( IWQSI == 1 ) WQDFSI(L) = SCB(L)*WQWSD(IMWQZ)*WQV(L,KSZ(L),16)  
      ENDDO  
      IF( IWQSRP == 1 )THEN  
        DO LP=LF,LL
          L=LWET(LP)
          IMWQZ = IWQZMAP(L,KSZ(L))  
          WQDFLP(L) = SCB(L)*( WQDFLP(L) + WQWSS(IMWQZ)*( WQV(L,KSZ(L),10)-WQPO4D(L,1) ) )  
          IF( IWQSI == 1 ) WQDFSI(L) = SCB(L)*( WQDFSI(L) + WQWSS(IMWQZ)*( WQV(L,KSZ(L),17)-WQSAD(L,1) ) )  
        ENDDO  
      ELSE IF( IWQSRP == 2 )THEN  
        DO LP=LF,LL
          L=LWET(LP)
          WQDFLP(L) = SCB(L)*( WQDFLP(L)+WSEDO(NS)*( WQV(L,KSZ(L),10)-WQPO4D(L,1) ) )  
          IF( IWQSI == 1 ) WQDFSI(L) = SCB(L)*( WQDFSI(L) + WSEDO(NS)*( WQV(L,KSZ(L),17)-WQSAD(L,1) ) )  
        ENDDO  
      ENDIF  
    ENDIF  
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO

  ! *** ZERO THIN LAYERS WHEN DELTA T WQ IS LARGE
  IF( KC > 1 .AND. ISDRY > 0 )THEN
    IFLAG = 0
    !$OMP PARALLEL DO PRIVATE(ND,LF,LL,LP,L,NW,K)
    DO ND=1,NDM 
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      DO LP=LF,LL
        L=LWET(LP)
        IF( HP(L) < 2.*HDRY )THEN
          DO NW = 1,NWQV
            IF( ISTRWQ(NW) > 0 )THEN
              DO K=KSZ(L),KC
                IF( WQV(L,K,NW) < 0.0 )THEN
                  WQV(L,K,NW) = 0.0
                  IFLAG = IFLAG + 1
                ENDIF
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    
    IF( IFLAG > 0 )THEN
      OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  
      WRITE(8,'(A,F12.4,I10)')' WARNING!  NEGATIVE WQ CONCENTRATIONS: TIMEDAY, # OF VALUES < 0.0:',TIMEDAY,IFLAG
      CLOSE(8)
    ENDIF
  ENDIF
  
  ! *** RESTORE PREVIOUS WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO NQ=0,NWQVM
    DO K=1,KC
      DO IOBC=1,NBCSOP  
        L=LOBCS(IOBC)
        WQV(L,K,NQ) = WQOLD(IOBC,K,NQ)
      ENDDO
    ENDDO  
  ENDDO  
  
  ! ***************************************************************************
  ! DIURNAL DO ANALYSIS  
  IF( NDDOAVG >= 1 .AND. DEBUG )THEN  
    OPEN(1,FILE=OUTDIR//'DIURNDO.OUT',POSITION='APPEND')  
    NDDOCNT=NDDOCNT+1  
    NSTPTMP=NDDOAVG*NTSPTC/2  
    RMULTMP=1./FLOAT(NSTPTMP)  
    DO K=1,KC  
      DO L=2,LA  
        DDOMAX(L,K)=MAX(DDOMAX(L,K),WQV(L,K,19))  
        DDOMIN(L,K)=MIN(DDOMIN(L,K),WQV(L,K,19))  
      ENDDO  
    ENDDO  
    IF( NDDOCNT == NSTPTMP )THEN  
      NDDOCNT=0  
      IF( ISDYNSTP == 0 )THEN  
        TIME=DT*FLOAT(N)+TCON*TBEGIN  
        TIME=TIME/TCON  
      ELSE  
        TIME=TIMESEC/TCON  
      ENDIF  
      WRITE(1,1111)N,TIME  
      DO L=2,LA  
        WRITE(1,1112)IL(L),JL(L),(DDOMIN(L,K),K=1,KC),(DDOMAX(L,K),K=1,KC)  
      ENDDO  
      DO K=1,KC  
        DO L=2,LA  
          DDOMAX(L,K)=-1.E6  
          DDOMIN(L,K)=1.E6  
        ENDDO  
      ENDDO  
    ENDIF  
    CLOSE(1)  
  ENDIF  

  ! ***************************************************************************
  ! *** LIGHT EXTINCTION ANALYSIS  
  IF( NDLTAVG >= 1 )THEN  
    OPEN(1,FILE=OUTDIR//'LIGHT.OUT',POSITION='APPEND')  
    NDLTCNT=NDLTCNT+1  
    NSTPTMP=NDLTAVG*NTSPTC/2  
    RMULTMP=1./FLOAT(NSTPTMP)  
    DO K=1,KC  
      DO L=2,LA 
        IF ( ISTRAN(6) > 0 ) THEN
          RLIGHT1=WQKEB(IMWQZT(L))+WQKETSS*SEDT(L,K) 
        ELSE
          RLIGHT1=WQKEB(IMWQZT(L))
        ENDIF
        XMRM = WQKECHL*WQCHL(L,K)  
        IF( WQKECHL  < 0.0 )THEN  
          XMRM = 0.054*WQCHL(L,K)**0.6667 + 0.0088*WQCHL(L,K)  
        ENDIF  
        RLIGHT2 = XMRM  
        RLIGHTT(L,K)=RLIGHTT(L,K)+RLIGHT1  
        RLIGHTC(L,K)=RLIGHTC(L,K)+RLIGHT1+RLIGHT2  
      ENDDO  
    ENDDO  
    IF( NDLTCNT == NSTPTMP )THEN  
      NDLTCNT=0  
      IF( ISDYNSTP == 0 )THEN  
        TIME=DT*FLOAT(N)+TCON*TBEGIN  
        TIME=TIME/TCON  
      ELSE  
        TIME=TIMESEC/TCON  
      ENDIF  
      DO K=1,KC  
        DO L=LF,LL  
          RLIGHTT(L,K)=RMULTMP*RLIGHTT(L,K)  
          RLIGHTC(L,K)=RMULTMP*RLIGHTC(L,K)  
        ENDDO  
      ENDDO  
      WRITE(1,1111)N,TIME  
      DO L=2,LA  
        WRITE(1,1113)IL(L),JL(L),(RLIGHTT(L,K),K=1,KC),(RLIGHTC(L,K),K=1,KC)  
      ENDDO  
      DO K=1,KC  
        DO L=2,LA  
          RLIGHTT(L,K)=0.  
          RLIGHTC(L,K)=0.  
        ENDDO  
      ENDDO  
    ENDIF  
    CLOSE(1)  
  ENDIF  

  ! ***************************************************************************
  ! INCREMENT COUNTER FOR LIMITATION AND XDOXXX DO COMPONENT ARRAYS:  
  IF( ISDYNSTP == 0 )THEN  
    TIMTMP=DT*FLOAT(N)+TCON*TBEGIN  
    TIMTMP=TIMTMP/TCTMSR  
  ELSE  
    TIMTMP=TIMESEC/TCTMSR  
  ENDIF  
  TIMESUM3 = TIMESUM3 + TIMTMP  
  
  1111 FORMAT(I12,F10.4)  
  1112 FORMAT(2I5,12F7.2)  
  1113 FORMAT(2I5,12E12.4)  
  1414 FORMAT(I12,11E12.4)  

  RETURN  

END  

