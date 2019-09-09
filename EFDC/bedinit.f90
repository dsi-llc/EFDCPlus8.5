SUBROUTINE BEDINIT  
  !  
  ! ***  SUBROUTINE BEDINIT INITIALIZES SEDIMENT AND TOXIC VARIABLES  
  ! ***  IN SEDIMENT BED FOR HOT AND COLD START CONDITIONS  
  !     CHANGE RECORD  
  !     ADDED ADDITIONAL DIAGNOSTIC OUTPUT  
  !     MOVED TOXIC INITIALIZATIONS FROM SSEDTOX  
  !  
  
  USE GLOBAL  
  USE INFOMOD,ONLY:SKIPCOM,READSTR

  IMPLICIT NONE
  
  INTEGER :: K,L,NS,NX,NT,KTOPP1,IVAL,KTOPTP,IHOTSTRT
  INTEGER :: NSKIP,ISSTYPE,LD,ID,JD,LCORE,IDUM,NL,KT1,KT2
  INTEGER :: ISO,NSITES,NCOHSEDS,NCOHSEDL
  CHARACTER*80 STR*200
  
  INTEGER,ALLOCATABLE,DIMENSION(:) :: LSSCOHSED 

  REAL :: CSEDTAUS,CSEDRESS,CSEDTAUB,CSEDRESB,TMPCVT
  REAL :: SURF,FVOLSSD,FVOLSED,FVOLSND,TMP1,SXD
  REAL :: RMULADJ,RUMLADJC,FRACT1,FRACT2,FRACT2C,HBEDP,TTHICK

  REAL,ALLOCATABLE,DIMENSION(:) :: FRACACT  
  REAL,ALLOCATABLE,DIMENSION(:) :: FRACPAR 
  REAL,ALLOCATABLE,DIMENSION(:) :: RADJCOHSEDS 
  REAL,ALLOCATABLE,DIMENSION(:) :: TAUDSS 
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TAURSS 
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TAUNSS 
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TEXPSS 
  REAL,ALLOCATABLE,DIMENSION(:,:) :: DEPBBSS 
  REAL,ALLOCATABLE,DIMENSION(:,:) :: WRSPOSS 
  REAL,ALLOCATABLE,DIMENSION(:,:) :: SEDBALL

  ALLOCATE(LSSCOHSED(LCM))
  ALLOCATE(FRACACT(LCM))  
  ALLOCATE(FRACPAR(LCM))  
  ALLOCATE(RADJCOHSEDS(LCM))
  ALLOCATE(SEDBALL(LCM,KBM))

  LSSCOHSED=0
  FRACACT=0.0
  FRACPAR=0.0
  RADJCOHSEDS=0.0
  SEDBALL=0.0

  ! *** SEDIMENT PROPERTIES ARE IN MASS PER UNIT AREA
  ! ***
  ! *** SEDB - COHESIVE SEDIMENTS (G/M2)
  ! *** SNDB - NONCOHESIVE SEDIMENTS (G/M2)
  ! *** TOXB - TOXICS (MG/M2)
  
  ! *** ZERO LOCAL ARRAYS
  FRACACT=0.0
  FRACPAR=0.0

  ! *** SITE SPECIFIC RESUSPENSION INFORMATION BASED ON HAMRICK'S
  ! *** ANALYSIS OF SEDFLUME CORES JANUARY 2008
  IF( NSED > 0 .AND. IWRSP(1) >= 99 .AND.  .NOT. LSEDZLJ )THEN
    WRITE(*,'(A)')'READING SSCOHSEDPROP.INP'
    OPEN(99,FILE='sscohsedprop.inp')
    STR=READSTR(99)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
    READ(99,*,IOSTAT=ISO)NCOHSEDS,NCOHSEDL,RMULADJ
    ALLOCATE(TAUDSS(NCOHSEDS))
    ALLOCATE(WRSPOSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(TAURSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(TAUNSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(TEXPSS(NCOHSEDL,NCOHSEDS))
    ALLOCATE(DEPBBSS(NCOHSEDL,NCOHSEDS))
    DO NSITES=1,NCOHSEDS
      READ(99,*,IOSTAT=ISO)IDUM,RUMLADJC,TAUDSS(NSITES)
      READ(99,*,IOSTAT=ISO)IDUM,(WRSPOSS(NL,NSITES),NL=1,NCOHSEDL)
      DO NL=1,NCOHSEDL
        WRSPOSS(NL,NSITES) = RUMLADJC*RMULADJ*WRSPOSS(NL,NSITES)
      ENDDO
      READ(99,*,IOSTAT=ISO)IDUM,(TAURSS(NL,NSITES),NL=1,NCOHSEDL)
      READ(99,*,IOSTAT=ISO)IDUM,(TAUNSS(NL,NSITES),NL=1,NCOHSEDL)
      READ(99,*,IOSTAT=ISO)IDUM,(TEXPSS(NL,NSITES),NL=1,NCOHSEDL)
      READ(99,*,IOSTAT=ISO)IDUM,(DEPBBSS(NL,NSITES),NL=1,NCOHSEDL)
    ENDDO
    CLOSE(99)
    IF( ISO > 0 )THEN
      STOP 'ERROR READING SSCOHSEDPROP.INP FILE'
    ENDIF
  ENDIF

  ! SET BOTTOM LAYER NUMBER: ORIGINAL:KBB=1, SEDZLJ:KBB=KB
  IF( LSEDZLJ )THEN
    KBB = KB
  ELSE
    KBB = 1
  ENDIF
  
  ! ***  DETERMINE START UP MODE  
  IHOTSTRT=0  
  IF( ISRESTI /= 0 )THEN  
    IF( ISCI(6) /= 0 .OR. ISCI(7) /= 0 )THEN  
      IHOTSTRT=1  
    ENDIF  
  ENDIF  

  ! ***********************************************************************************
  ! ***  HOT START INITIALIZATION.  SEDZLJ HOTSTART IS HANDLED IN S_SEDIC
  IF( IHOTSTRT /= 0 .AND. .NOT. LSEDZLJ )THEN

    ! ***  SET POROSITY  
    DO K=1,KB  
      DO L=2,LA  
        PORBED(L,K) = VDRBED(L,K)/(1.0+VDRBED(L,K))  
        PORBED1(L,K) = VDRBED1(L,K)/(1.0+VDRBED1(L,K))  
      ENDDO  
    ENDDO  

    ! ***  SET BULK DENSITY  
    DO K=1,KB  
      DO L=2,LA  
        SEDBT(L,K)=0.0
        SNDBT(L,K)=0.0
      ENDDO  
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            SEDBT(L,K) = SEDBT(L,K)+SEDB(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NS=1,NSND  
        DO K=1,KB  
          DO L=2,LA  
            SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF
    
    DO K=1,KB  
      DO L=2,LA  
        IF( HBED(L,K) > 0.0 )THEN
          ! *** COMPUTE TOTAL/WET DENSITY
          BDENBED(L,K) = 1000.0*PORBED(L,K)+0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)  
        ELSE  
          BDENBED(L,K) = 0.0
        ENDIF  
      ENDDO  
    ENDDO
    DO K=1,KB  
      DO L=2,LA  
        SEDBT(L,K)=0.0
        SNDBT(L,K)=0.0
      ENDDO  
    ENDDO  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            SEDBT(L,K) = SEDBT(L,K)+SEDB1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NS=1,NSND  
        DO K=1,KB  
          DO L=2,LA  
            SNDBT(L,K) = SNDBT(L,K)+SNDB1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    DO K=1,KB  
      DO L=2,LA  
        IF( HBED1(L,K) > 0.0 )THEN  
          BDENBED1(L,K) = 1000.0*PORBED1(L,K)+0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED1(L,K)  
          ELSE  
          BDENBED1(L,K) = 0.0
        ENDIF
      ENDDO
    ENDDO

    ! ***  SET TOP BED LAYER  
    DO L=2,LA  
      KBT(L)=1  
    ENDDO  
    DO K=1,KB  
      DO L=2,LA  
        IF( HBED(L,K) > 0.)KBT(L)=K  
      ENDDO  
    ENDDO  

    ! ***  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES  
    IF( ISTRAN(6) >= 1 )THEN
      IF( IWRSP(1) == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURS(L,K)=TAUR(1)  
            WRSPS(L,K)=WRSPO(1)  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( IWRSPB(1) == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURB(L,K)=1.E6  
            WRSPB(L,K)=0.0  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( IWRSP(1) >= 1 .AND. IWRSP(1) < 99 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)  
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1))  
          ENDDO  
        ENDDO  
      ENDIF 

      ! *** READ IN THE SPATIALLY ASSIGNED COHESIVE CRITICAL BED SHEAR STRESS AND SURFACE EROSION RATE
      IF( IWRSP(1) >= 99 .AND. ISHOUSATONIC == 0 )THEN
        WRITE(*,'(A)')'READING SSCOHSEDPMAP.INP'
        OPEN(1,FILE='sscohsedpmap.inp')
        OPEN(2,FILE=OUTDIR//'SSCOHSEDPMAP.OUT')
        STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
        READ(1,*)ISSTYPE
        IF( ISSTYPE == 0 )THEN
          DO L=2,LA
            READ(1,*)LD,ID,JD,LSSCOHSED(L)
            RADJCOHSEDS(L)=1.
          ENDDO
        ELSE
          DO L=2,LA
            READ(1,*)LD,ID,JD,LSSCOHSED(L),RADJCOHSEDS(L)
          ENDDO
        ENDIF
        DO L=2,LA
          LCORE=LSSCOHSED(L)
          TAUDS(L)=TAUDSS(LCORE)
        ENDDO
        IF( NCOHSEDL == 1 )THEN
          DO K=1,KB
            DO L=2,LA
              LCORE=LSSCOHSED(L)
              TAURS(L,K) = TAURSS(1,LCORE)
              TAUNS(L,K) = TAUNSS(1,LCORE)
              WRSPS(L,K) = RADJCOHSEDS(L)*WRSPOSS(1,LCORE)
              TEXPS(L,K) = TEXPSS(1,LCORE)
            ENDDO
          ENDDO
        ELSE
          DO K=1,KB
            DO L=2,LA
              LCORE = LSSCOHSED(L)
              TAURS(L,K) = TAURSS(K,LCORE)
              !TAUNS(L,K) = TAUNSS(K,LCORE)
              WRSPS(L,K) = RADJCOHSEDS(L)*WRSPOSS(K,LCORE)
              !TEXPS(L,K) = TEXPSS(K,LCORE)
            ENDDO
          ENDDO
          ENDIF
        DO L=2,LA
          K=KBT(L)
          WRITE(2,*)L,IL(L),JL(L),TAURS(L,K), WRSPS(L,K) 
        ENDDO
        CLOSE(1)
        CLOSE(2)
      ENDIF
      
      IF( IWRSPB(1) >= 1 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))  
            WRSPB(L,K) = CSEDRESB
          ENDDO  
        ENDDO  
      ENDIF  
    ENDIF  

    ! ***  SET SEDIMENT VOLUME FRACTIONS  
    DO K=1,KB  
      DO L=2,LA  
        BEDLINIT(L,K)=0.0
        BEDDINIT(L,K)=0.0
      ENDDO  
    ENDDO
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)  
            VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)  
            VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  

    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)  
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)  
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA
            ! *** BEGIN DS-INTL
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K) 
            ELSE  
              VFRBED(L,K,NS)=0.0
              BEDLINIT(L,K) =0.0  
            ENDIF  
            IF( BEDDINIT(L,K) > 0.0 )THEN  
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)  
            ELSE  
              VFRBED1(L,K,NS)=0.0  
              BEDDINIT(L,K) =0.0  
            ENDIF  
            ! *** END DS-INTL
          ENDDO
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K) 
            ELSE  
              VFRBED(L,K,NS)=0.0
              BEDLINIT(L,K) =0.0  
            ENDIF  
            IF( BEDDINIT(L,K) > 0.0 )THEN  
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)  
            ELSE  
              VFRBED1(L,K,NS)=0.0  
              BEDDINIT(L,K) =0.0  
            ENDIF  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  

    ! ***  INITIALIZE BED BOTTOM ELEVATION  
    DO L=2,LA  
      HBEDA(L)=0.0
    ENDDO  
    DO L=2,LA  
      DO K=1,KBT(L)  
        HBEDA(L) = HBEDA(L)+HBED(L,K)  
      END DO  
    ENDDO  
    DO L=2,LA  
      ZELBEDA(L) = BELV(L)-HBEDA(L)  
    ENDDO  

    ! ***  INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA  
    DO K=1,KB  
      DO L=2,LA  
        SEDBT(L,K)=0.0
        SNDBT(L,K)=0.0
      ENDDO  
    ENDDO  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            SEDBT(L,K) = SEDBT(L,K)+SEDB(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NS=1,NSND  
        DO K=1,KB  
          DO L=2,LA  
            SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    
    IF( LSEDZLJ )THEN
      DO L=2,LA
        FORALL(K=1:KB) SEDDIA50(L,K) = SUM(PERSED(1:NSCM,K,L)*D50(1:NSCM))
      ENDDO
    ENDIF

    GOTO 1000    ! *** JUMP TO FINAL DIAGNOSTICS PRINTOUT
 
  ENDIF  
  ! ***  END HOT START INITIALIZATION  

  ! ***********************************************************************************
  ! ***  COLD START INITIALIZATION: IBMECH=0  
  IF( (IBMECH == 0 .OR. ISEDINT <= 1) .AND. .NOT. LSEDZLJ )THEN

    ! ***  SET POROSITY AND VOID RATIO  
    DO K=1,KB  
      DO L=2,LA  
        PORBED(L,K) = BEDPORC
        PORBED1(L,K) = BEDPORC
        VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))
        VDRBED1(L,K) = PORBED1(L,K)/(1.0-PORBED1(L,K))
        HBED(L,K) = 0.0
        HBED1(L,K) = 0.0
        KBT(L) = 1
      ENDDO  
    ENDDO  

    ! ***  UNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS  
    ! ***  CALCULATE LAYER THICKNESS AND BULK DENSITY  
    IF( ISEDINT <= 1 )THEN  
      DO K=1,KB  
        DO L=2,LA  
          SEDBT(L,K)=0.0
          SNDBT(L,K)=0.0
        ENDDO  
      ENDDO  
      IF( ISTRAN(6) >= 1 )THEN  
        DO NS=1,NSED  
          DO K=1,KB  
            DO L=2,LA
              HBED(L,K) = HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)
              SEDBT(L,K) = SEDBT(L,K)+SEDB(L,K,NS)
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(7) >= 1 )THEN  
        DO NX=1,NSND  
          NS=NSED+NX  
          DO K=1,KB  
            DO L=2,LA  
              HBED(L,K) = HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)  
              SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NX)  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF 
      DO K=1,KB  
        DO L=2,LA  
          HBED(L,K) = (1.0+VDRBED(L,K))*HBED(L,K)  
          IF( HBED(L,K) > 0.0) KBT(L)=K  
        ENDDO
      ENDDO
      DO K=1,KB  
        DO L=2,LA  
          IF( HBED(L,K) > 0.0 )THEN  
            ! *** COMPUTE TOTAL/WET DENSITY
            BDENBED(L,K) = 1000.*PORBED(L,K) + 0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)  
           ELSE  
            BDENBED(L,K) = 0.0
          ENDIF
        ENDDO
      ENDDO
      DO K=1,KB  
        DO L=2,LA  
          HBED1(L,K) = HBED(L,K)  
          BDENBED1(L,K) = BDENBED(L,K)  
        ENDDO  
      ENDDO  

      ! ***  HANDLE ACTIVE ARMORING LAYER IF NOT PRESENT IN INITIAL OR RESTART  
      IF( ISNDAL == 2 .AND. IALSTUP > 0 .AND. KB > 1 )THEN  
        DO L=2,LA  
          KBT(L)=KBT(L)-1
          IF( KBT(L) < 1)KBT(L)=1
          HBED(L,KBT(L)+1)  = 0.
          HBED1(L,KBT(L)+1) = 0.
        ENDDO
      ENDIF
      
    ENDIF  

    ! ***  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS  
    ! ***  AND INITIAL CONDITIONS IN SEDB.INP ARE IN MASS PER UNIT AREA  
    ! ***  CALCULATE LAYER THICKNESS AND BULK DENSITY  
    IF( ISEDINT >= 2 )THEN  
      IF( ISEDBINT == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            SEDBT(L,K) = 0.0
            SNDBT(L,K) = 0.0
          ENDDO  
        ENDDO  
        IF( ISTRAN(6) >= 1 )THEN  
          DO NS=1,NSED  
            DO K=1,KB  
              DO L=2,LA  
                HBED(L,K) = HBED(L,K)+SDEN(NS)*SEDB(L,K,NS)  
                SEDBT(L,K) = SEDBT(L,K)+SEDB(L,K,NS)  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
        IF( ISTRAN(7) >= 1 )THEN  
          DO NX=1,NSND  
            NS=NSED+NX  
            DO K=1,KB  
              DO L=2,LA  
                HBED(L,K) = HBED(L,K)+SDEN(NS)*SNDB(L,K,NX)  
                SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NX)  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
        DO K=1,KB  
          DO L=2,LA  
            HBED(L,K) = (1.+VDRBED(L,K))*HBED(L,K)  
            IF( HBED(L,K) > 0.0) KBT(L)=K
          ENDDO  
        ENDDO  
        DO K=1,KB  
          DO L=2,LA  
            IF( HBED(L,K) > 0.0 )THEN  
              ! *** COMPUTE TOTAL/WET DENSITY
              BDENBED(L,K) = 1000.*PORBED(L,K)  +0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)  
            ELSE  
              BDENBED(L,K) = 0.0
            ENDIF  
          ENDDO  
        ENDDO  
        DO K=1,KB  
          DO L=2,LA  
            HBED1(L,K) = HBED(L,K)  
            BDENBED1(L,K) = BDENBED(L,K)  
          ENDDO  
        ENDDO  
      ENDIF  
    ENDIF  

    ! ***  NONUNIFORM SEDIMENT MASS PER UNIT AREA ALL CELLS, ALL BED LAYERS  
    ! ***  AND INITIAL CONDITIONS  IN SEDB.INP ARE IN MASS FRACTION  
    ! ***  CALCULATE LAYER THICKNESS AND BULK DENSITY  
    ! ***  THIS OPTION REQUIRES INITIAL LAYER THICKNESSES  
    IF( ISEDINT >= 2 )THEN  
      IF( ISEDBINT == 1 )THEN  
        IF( IBEDLAYU == 1 )THEN  
          DO K=1,KB  
            DO L=2,LA  
              BEDLINIT(L,K) = 0.1*BEDLINIT(L,K)  
            ENDDO  
          ENDDO  
        ENDIF  
        DO K=1,KB  
          DO L=2,LA  
            HBED(L,K) = BEDLINIT(L,K)  
            HBED1(L,K) = BEDLINIT(L,K)  
            IF( HBED(L,K) > 0.0)KBT(L)=K  
          ENDDO  
        ENDDO  
        DO K=1,KB  
          DO L=2,LA  
            BDENBED(L,K) = 0.0
          ENDDO  
        ENDDO  
        
        ! *** FIRST COMPUTE DRY DENSITY (G/M2)
        IF( ISTRAN(6) >= 1 )THEN  
          DO NS=1,NSED  
            DO K=1,KB  
              DO L=2,LA  
                BDENBED(L,K) = BDENBED(L,K) + 1000.0*SSG(NS)*SEDB(L,K,NS)  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
        IF( ISTRAN(7) >= 1 )THEN  
          DO NX=1,NSND  
            NS=NSED+NX  
            DO K=1,KB  
              DO L=2,LA  
                BDENBED(L,K) = BDENBED(L,K) + 1000.0*SSG(NS)*SNDB(L,K,NX)  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
        
        ! *** CONVERT DRY DENSITY TO TOTAL/WET DENSITY
        DO K=1,KB  
          DO L=2,LA  
            BDENBED(L,K) = 1000.0*PORBED(L,K) + (1.0-PORBED(L,K))*BDENBED(L,K)  
            BDENBED1(L,K) = BDENBED(L,K)  
          ENDDO  
        ENDDO  
        
        ! *** TOTAL SEDIMENT FOR ALL CLASSES (G/M2)
        DO K=1,KB  
          DO L=2,LA  
            SEDBT(L,K) = 1000.0*HBED(L,K)*(BDENBED(L,K) - 1000.0*PORBED(L,K))  
          ENDDO  
        ENDDO  
        
        ! *** SPLIT TOTAL SEDIMENT INTO CLASSES
        IF( ISTRAN(6) >= 1 )THEN  
          DO NS=1,NSED  
            DO K=1,KB  
              DO L=2,LA  
                SEDB(L,K,NS) = SEDB(L,K,NS)*SEDBT(L,K)  
                SEDB1(L,K,NS) = SEDB(L,K,NS)  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
        IF( ISTRAN(7) >= 1 )THEN  
          DO NX=1,NSND  
            DO K=1,KB  
              DO L=2,LA  
                SNDB(L,K,NX) = SNDB(L,K,NX)*SEDBT(L,K)  
                SNDB1(L,K,NX) = SNDB(L,K,NX)  
              ENDDO  
            ENDDO  
          ENDDO  
        ENDIF  
      ENDIF  
    ENDIF  

    ! ***  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES  
    IF( ISTRAN(6) >= 1 .AND. .NOT. LSEDZLJ )THEN  
      IF( IWRSP(1) == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURS(L,K) = TAUR(1)  
            WRSPS(L,K) = WRSPO(1)  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( IWRSPB(1) == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURB(L,K) = 1.E6  
            WRSPB(L,K) = 0.0  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( IWRSP(1) >= 1 .AND. IWRSP(1) < 99 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)  
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1))  
          ENDDO  
        ENDDO  
      ENDIF  
      
      ! *** READ IN THE SPATIALLY ASSIGNED COHESIVE CRITICAL BED SHEAR STRESS AND SURFACE EROSION RATE
      IF( IWRSP(1) >= 99 .AND. ISHOUSATONIC == 0 )THEN
        WRITE(*,'(A)')'READING SSCOHSEDPMAP.INP'
        OPEN(1,FILE='sscohsedpmap.inp')
        OPEN(2,FILE=OUTDIR//'SSCOHSEDPMAP.OUT')
        STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
        READ(1,*)ISSTYPE
        IF( ISSTYPE == 0 )THEN
          DO L=2,LA
            READ(1,*)LD,ID,JD,LSSCOHSED(L)
            RADJCOHSEDS(L) = 1.
          ENDDO
        ELSE
          DO L=2,LA
            READ(1,*)LD,ID,JD,LSSCOHSED(L),RADJCOHSEDS(L)
          ENDDO
        ENDIF
        DO L=2,LA
          LCORE = LSSCOHSED(L)
          TAUDS(L) = TAUDSS(LCORE)
        ENDDO
        IF( NCOHSEDL == 1 )THEN
          DO K=1,KB
            DO L=2,LA
              LCORE = LSSCOHSED(L)
              TAURS(L,K) = TAURSS(1,LCORE)
              TAUNS(L,K) = TAUNSS(1,LCORE)
              WRSPS(L,K) = RADJCOHSEDS(L)*WRSPOSS(1,LCORE)
              TEXPS(L,K) = TEXPSS(1,LCORE)
            ENDDO
          ENDDO
        ELSE
          DO K=1,KB
            DO L=2,LA
              LCORE = LSSCOHSED(L)
              TAURS(L,K) = TAURSS(K,LCORE)
              TAUNS(L,K) = TAUNSS(K,LCORE)
              WRSPS(L,K) = RADJCOHSEDS(L)*WRSPOSS(K,LCORE)
              TEXPS(L,K) = TEXPSS(K,LCORE)
            ENDDO
          ENDDO
        ENDIF
        DO L=2,LA
          K=KBT(L)
          WRITE(2,*)L,IL(L),JL(L),TAURS(L,K),WRSPS(L,K)
        ENDDO
        CLOSE(1)
        CLOSE(2)
      ENDIF
      IF( IWRSPB(1) >= 1 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))  
            WRSPB(L,K) = CSEDRESB
          ENDDO  
        ENDDO  
      ENDIF  
    ENDIF  

    ! ***  SET SEDIMENT VOLUME FRACTIONS  
    DO K=1,KB  
      DO L=2,LA  
        BEDLINIT(L,K) = 0.0
        BEDDINIT(L,K) = 0.0
      ENDDO  
    ENDDO  
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)  
            VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)  
            VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(6) >= 1 )THEN
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)  
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)  
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K) 
            ELSE  
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0  
            ENDIF  
            IF( BEDDINIT(L,K) > 0.0 )THEN  
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)  
            ELSE  
              VFRBED1(L,K,NS) = 0.0  
              BEDDINIT(L,K)  = 0.0  
            ENDIF  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K) 
            ELSE  
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0  
            ENDIF  
            IF( BEDDINIT(L,K) > 0.0 )THEN  
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)  
            ELSE  
              VFRBED1(L,K,NS) = 0.0  
              BEDDINIT(L,K)  = 0.0  
            ENDIF  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
  ENDIF  
  ! ***  END COLD START INITIALIZATION: IBMECH=0  

  ! ***********************************************************************************
  ! ***  COLD START INITIALIZATION: IBMECH >= 1  
  IF( IBMECH >= 1 .AND. ISEDINT > 1 .AND. .NOT. LSEDZLJ )THEN

    ! ***  CONVERT AND INITIALIZE BED LAYER THICKNESS AND DEFINE  
    ! ***  INITIAL TOP LAYER  
    IF( IBEDLAYU == 1 ) TMPCVT=0.001    ! *** THICKNESS IN MM
    IF( IBEDLAYU == 2 ) TMPCVT=0.01     ! *** THICKNESS IN CM 
    IF( IBEDLAYU == 3 ) TMPCVT=1.0      ! *** THICKNESS IN M
    IF( IBEDLAYU >= 1 )THEN  
      DO K=1,KB  
        DO L=2,LA  
          BEDLINIT(L,K) = TMPCVT*BEDLINIT(L,K)  
        ENDDO  
      ENDDO  
    ENDIF  
    DO L=2,LA  
      KBT(L)=1  
    ENDDO  
    DO K=1,KB  
      DO L=2,LA  
        HBED(L,K) = BEDLINIT(L,K)  
        HBED1(L,K) = BEDLINIT(L,K)  
        IF( HBED(L,K) > 0.0) KBT(L)=K  
      ENDDO  
    ENDDO  

    ! ***  CONVERT AND INITIALIZE BED BULK DENSITY  
    ! ***   IBEDBDNU=0 BEDBINIT IS NOT BULK DENSITY  
    ! ***   IBEDBDNU=1 BEDBINIT BULK DENSITY IN KG/M**3  
    ! ***   IBEDBDNU=3 BEDBINIT BULK DENSITY IN GM/CM**3  
    IF( IBEDBDNU >= 1 )THEN  
      IF( IBEDBDNU == 2 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            BEDBINIT(L,K) = 1000.*BEDBINIT(L,K)  
          ENDDO  
        ENDDO  
      ENDIF  
      DO K=1,KB  
        DO L=2,LA  
          BDENBED(L,K) = BEDBINIT(L,K)  
          BDENBED1(L,K) = BEDBINIT(L,K)  
        ENDDO  
      ENDDO  
    ENDIF  

    ! ***  CONVERT AND DRY DENSITY OF BED  
    ! ***  IBEDDDNU=0, 1 ACTUAL DRY DENSITY, = 2 POROSITY, = 3 VOID RATIO  
    IF( IBEDDDNU == 1 )THEN  
      DO K=1,KB  
        DO L=2,LA  
          BEDDINIT(L,K) = 1000.*BEDDINIT(L,K)  
        ENDDO  
      ENDDO  
    ENDIF  

    ! ***  CALCULATE POROSITY AND VOID RATIO  
    IF( IBEDDDNU <= 1 )THEN  
      DO K=1,KB  
        DO L=2,LA  
          PORBED(L,K) = 0.001*(BEDBINIT(L,K)-BEDDINIT(L,K))
          VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( IBEDDDNU == 2 )THEN  
      DO K=1,KB  
        DO L=2,LA  
          PORBED(L,K) = BEDDINIT(L,K)  
          VDRBED(L,K) = PORBED(L,K)/(1.-PORBED(L,K))  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( IBEDDDNU == 3 )THEN  
      DO K=1,KB  
        DO L=2,LA  
          VDRBED(L,K) = BEDDINIT(L,K)  
          PORBED(L,K) = VDRBED(L,K)/(1.+VDRBED(L,K))  
        ENDDO  
      ENDDO  
    ENDIF  
    DO K=1,KB  
      DO L=2,LA  
        VDRBED1(L,K) = VDRBED(L,K)  
        PORBED1(L,K) = PORBED(L,K)  
      ENDDO  
    ENDDO  

    ! ***  INITIALIZE BED SEDIMENT FOR MASS FRACTION INPUT BY CACLULATING  
    ! ***  AND STORING TOTAL MASS OF SED/AREA IN BEDDINIT(L,K)  
    DO K=1,KB  
      DO L=2,LA  
        ! ***              M         KG/M3             KG/M3
        BEDDINIT(L,K) = HBED(L,K)*(BDENBED(L,K) - 1000.*PORBED(L,K))     ! *** BEDDINIT is Dry Density in KG/M2
      ENDDO  
    ENDDO  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        IF( ISEDBU(NS) == 1 )THEN  
          DO K=1,KB  
            DO L=2,LA  
              ! *** G/M2     G/KG    FRAC (NO DIM)     KG/M2
              SEDB(L,K,NS) = 1000.*SEDBINIT(L,K,NS)*BEDDINIT(L,K)        
              SEDB1(L,K,NS) = SEDB(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDIF  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NS=1,NSND  
        IF( ISNDBU(NS) == 1 )THEN  
          DO K=1,KB  
            DO L=2,LA  
              SNDB(L,K,NS) = 1000.*SNDBINIT(L,K,NS)*BEDDINIT(L,K)  
              SNDB1(L,K,NS) = SNDB(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDIF  
      ENDDO  
    ENDIF  

    ! ***  SET COHESIVE BED CRITICAL STRESSES AND RESUSPENSION RATES  
    IF( ISTRAN(6) >= 1 )THEN  
      IF( IWRSP(1) == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURS(L,K) = TAUR(1)  
            WRSPS(L,K) = WRSPO(1)  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( IWRSPB(1) == 0 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURB(L,K) = 1.E6  
            WRSPB(L,K) = 0.0  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( IWRSP(1) >= 1 .AND. IWRSP(1) < 99 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURS(L,K) = CSEDTAUS(BDENBED(L,K),TAUR(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1),L)  
            WRSPS(L,K) = CSEDRESS(BDENBED(L,K),WRSPO(1),VDRRSPO(1),  VDRBED(L,K),VDRBED(L,K),IWRSP(1))  
          ENDDO  
        ENDDO  
      ENDIF  
      
      ! *** READ IN THE SPATIALLY ASSIGNED COHESIVE CRITICAL BED SHEAR STRESS AND SURFACE EROSION RATE
      IF( IWRSP(1) >= 99 .AND. ISHOUSATONIC == 0 )THEN
        WRITE(*,'(A)')'READING SSCOHSEDPMAP.INP'
        OPEN(1,FILE='sscohsedpmap.inp')
        OPEN(2,FILE=OUTDIR//'SSCOHSEDPMAP.OUT')
        STR=READSTR(1)  ! *** SKIP OVER TITLE AND AND HEADER LINES 
        READ(1,*)ISSTYPE
        IF( ISSTYPE == 0 )THEN
          DO L=2,LA
            READ(1,*)LD,ID,JD,LSSCOHSED(L)
            RADJCOHSEDS(L)=1.
          ENDDO
        ELSE
          DO L=2,LA
            READ(1,*)LD,ID,JD,LSSCOHSED(L),RADJCOHSEDS(L)
          ENDDO
        ENDIF
        DO L=2,LA
          LCORE = LSSCOHSED(L)
          TAUDS(L) = TAUDSS(LCORE)
        ENDDO
        IF( NCOHSEDL == 1 )THEN
          DO K=1,KB
            DO L=2,LA
              LCORE = LSSCOHSED(L)
              TAURS(L,K) = TAURSS(1,LCORE)
              TAUNS(L,K) = TAUNSS(1,LCORE)
              WRSPS(L,K) = RADJCOHSEDS(L)*WRSPOSS(1,LCORE)
              TEXPS(L,K) = TEXPSS(1,LCORE)
            ENDDO
          ENDDO
        ELSE
          DO K=1,KB
            DO L=2,LA
              LCORE = LSSCOHSED(L)
              TAURS(L,K) = TAURSS(K,LCORE)
              TAUNS(L,K) = TAUNSS(K,LCORE)
              WRSPS(L,K) = RADJCOHSEDS(L)*WRSPOSS(K,LCORE)
              TEXPS(L,K) = TEXPSS(K,LCORE)
            ENDDO
          ENDDO
        ENDIF
        DO L=2,LA
          K=KBT(L)
          WRITE(2,*)L,IL(L),JL(L),TAURS(L,K), TAUNS(L,K),WRSPS(L,K), TEXPS(L,K)
        ENDDO
        CLOSE(1)
        CLOSE(2)
      ENDIF

      IF( IWRSPB(1) >= 1 )THEN  
        DO K=1,KB  
          DO L=2,LA  
            TAURB(L,K) = CSEDTAUB(BDENBED(L,K),IWRSPB(1))  
            WRSPB(L,K) = CSEDRESB
          ENDDO  
        ENDDO  
      ENDIF  
    ENDIF    ! *** END OF SETTING RESUSPENSION RATES FOR COHESIVES

    ! ***  SET SEDIMENT VOLUME FRACTIONS  
    DO K=1,KB  
      DO L=2,LA  
        BEDLINIT(L,K) = 0.  
        BEDDINIT(L,K) = 0.  
      ENDDO  
    ENDDO  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            VFRBED(L,K,NS) = SDEN(NS)*SEDB(L,K,NS)  
            VFRBED1(L,K,NS) = SDEN(NS)*SEDB1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            VFRBED(L,K,NS) = SDEN(NS)*SNDB(L,K,NX)  
            VFRBED1(L,K,NS) = SDEN(NS)*SNDB1(L,K,NX)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)  
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            BEDLINIT(L,K) = BEDLINIT(L,K)+VFRBED(L,K,NS)  
            BEDDINIT(L,K) = BEDDINIT(L,K)+VFRBED1(L,K,NS)  
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            ! *** BEGIN DS-INTL
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K) 
            ELSE  
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0  
            ENDIF  
            IF( BEDDINIT(L,K) > 0.0 )THEN  
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)  
            ELSE  
              VFRBED1(L,K,NS) = 0.0  
              BEDDINIT(L,K)  = 0.0  
            ENDIF  
            ! *** END DS-INTL
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NX=1,NSND  
        NS=NSED+NX  
        DO K=1,KB  
          DO L=2,LA  
            ! *** BEGIN DS-INTL
            IF( BEDLINIT(L,K) > 0.0 )THEN
              VFRBED(L,K,NS) = VFRBED(L,K,NS)/BEDLINIT(L,K) 
            ELSE  
              VFRBED(L,K,NS) = 0.0
              BEDLINIT(L,K)  = 0.0  
            ENDIF  
            IF( BEDDINIT(L,K) > 0.0 )THEN  
              VFRBED1(L,K,NS) = VFRBED1(L,K,NS)/BEDDINIT(L,K)  
            ELSE  
              VFRBED1(L,K,NS) = 0.0  
              BEDDINIT(L,K)  = 0.0  
            ENDIF  
            ! *** END DS-INTL
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
  ENDIF  
  ! ***  END COLD START INITIALIZATION: IBMECH >= 1  

  ! ***********************************************************************************
  ! ***  INITIALIZE BED BOTTOM ELEVATION  
  IF( .NOT. LSEDZLJ )THEN
    DO L=2,LA  
      HBEDA(L)=0.0
    ENDDO  
    DO L=2,LA  
      DO K=1,KBT(L)  
        HBEDA(L) = HBEDA(L)+HBED(L,K)  
      END DO  
    ENDDO  
    DO L=2,LA  
      ZELBEDA(L) = BELV(L) - HBEDA(L)    ! *** HARD BOTTOM ELEVATION
    ENDDO  

    ! ***  INITIALIZE TOTAL SEDIMENT MASS PER UNIT AREA  
    DO K=1,KB  
      DO L=2,LA  
        SEDBT(L,K) = 0.0
        SNDBT(L,K) = 0.0
      ENDDO  
    ENDDO  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO K=1,KB  
          DO L=2,LA  
            IF( SEDB(L,K,NS) > 0.0 )THEN
              SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)  
            ELSE
              SEDB(L,K,NS) = 0.0
            ENDIF
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NS=1,NSND  
        DO K=1,KB  
          DO L=2,LA
            IF( SNDB(L,K,NS) > 0.0 )THEN
              SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NS)  
            ELSE
              SNDB(L,K,NS) = 0.0
            ENDIF
          ENDDO  
        ENDDO  
      ENDDO  
    ENDIF  

    ! *** GET TOTAL THICKNESSES AND BEDROCK ELEVATION
    DO L=2,LA  
      HBEDA(L)=0.0
    ENDDO  
    DO L=2,LA  
      DO K=1,KB
        IF( SEDBT(L,K) <= 1E-6 )THEN
          SEDB(L,K,1:NSED) = 0.0 
          SEDBT(L,K) = 0.0
        ENDIF
        IF( SNDBT(L,K) <= 1E-6 )THEN
          SNDB(L,K,1:NSND) = 0.0 
          SNDBT(L,K) = 0.0
        ENDIF
        IF( (SEDBT(L,K) + SNDBT(L,K)) <= 1E-6 )THEN
          HBED(L,K) = 0.0
          SEDB(L,K,1:NSED) = 0.0 
          SNDB(L,K,1:NSND) = 0.0 
        ENDIF
        HBEDA(L) = HBEDA(L) + HBED(L,K)  
        
        ! *** QC
        IF( HBED(L,K) > 0. .AND. PORBED(L,K) > 0. )THEN
          TTHICK = 0.
          DO NS=1,NSED  
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SEDB(L,K,NS)*DSEDGMM
          ENDDO  
          DO NX=1,NSND  
            NS=NSED+NX  
            DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
            TTHICK = TTHICK + SNDB(L,K,NX)*DSEDGMM
          ENDDO 
          TTHICK = TTHICK/(1.0 - PORBED(L,K))
          IF( ABS(HBED(L,K) - TTHICK) > 1E-4 )THEN
            WRITE(7,'(" *** WARNING - BED LAYER THICKNESS IS NOT CONSISTENT.  ADJUSTING THICKNESS (L, K, HOLD, HNEW): ",2I5,2F12.5)') L,K,HBED(L,K),TTHICK
            HBED(L,K) = TTHICK
          ENDIF
        ELSE
          HBED(L,K) = 0.0
          PORBED(L,K) = BEDPORC
          VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))
          SEDB(L,K,:) = 0.
          SNDB(L,K,:) = 0.
        ENDIF
      END DO  
      IF( HBEDA(L) <= 1E-6 ) KBT(L) = 1
    ENDDO  
    DO L=2,LA  
      ZELBEDA(L) = BELV(L) - HBEDA(L)    ! *** HARD BOTTOM ELEVATION
    ENDDO  

  ENDIF
  
  ! ***********************************************************************************
  ! ***  IF N=1 AND ISTRAN(5)=1 CHECK INITIAL TOXIC CONCENTRATIONS IN  
  ! ***  BED AND REINITILIZE IF NECESSARY  
  IF( ISTRAN(5) >= 1 )THEN  
    ! *** GET EACH SOLIDS CLASS FOR EACH TOXIC: NSP2
    DO NT=1,NTOX  
      NSP2(NT) = NSED + NSND                      ! *** Kd  Approach ISTOC(NT)=1
      IF( ISTOC(NT) == 1 ) NSP2(NT) = NSP2(NT)+2  ! *** DOC AND POC (POC NON-SEDIMENT RELATED)  (3 Phase)
      IF( ISTOC(NT) == 2 ) NSP2(NT) = NSP2(NT)+1  ! *** DOC AND POC FRACTIONALLY DISTRIBUTED    (3 Phase)
      ! *** ISTOC(NT)=0 and ISTOC(NT)=3               POC fOC*SED/SND BASED ONLY              (2 Phase) 
    END DO  

    IF( ISRESTI == 0 .OR. ISCI(5) == 0 )THEN  

      ! ***  CALCULATE TOTAL SEDIMENT IN THE BED  
      DO K=1,KB  
        DO L=1,LC  
          SEDBALL(L,K) = 0.  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=1,LC  
          SEDBALL(L,K) = SEDBT(L,K)+SNDBT(L,K)  
        ENDDO  
      ENDDO  

      !**********************************************************************C
      !
      ! **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
      !
      ! **  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED  (DIMENSIONLESS)
      ! **  TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED  (DIMENSIONLESS)
      CALL CALTOXB_FRACTIONS(2,LA)

      ! *** COMPUTE FINAL FRACTION (DIMENSIONLESS) ADSORBED ON ALL SOLIDS AND DOC (I.E. NOT DISSOLVED)
      DO NT=1,NTOX  
        DO K=1,KB  
          DO L=2,LA
            IF( SEDBALL(L,K) > 0.0 )THEN
              ! ***                 M           (   POREWATER (M)        SOLIDS (M)   )   
              TOXPFTB(L,K,NT) = TOXPFTB(L,K,NT)/(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))  
            ELSE  
              TOXPFTB(L,K,NT) = 1.0
            ENDIF  
          ENDDO  
        ENDDO  
      ENDDO  

      ! ***  CONVERT MASS TOX/MASS SED INITIAL CONDITION TO TOTAL TOXIC  
      ! ***  CONCENTRATION IN BED 0.001 CONVERTS TOXINTB UNITS OF MG/KG  
      ! ***  TO TOXB UNITS OF OF MG/M**2  
      DO NT=1,NTOX  
        IF( ITXBDUT(NT) == 0 )THEN
          ! *** INPUT UNITS OF UG/L ( = MG/M3)
          DO K=1,KB  
            DO L=2,LA  
              ! *** MG/M2   =   M          MG/M3
              TOXB(L,K,NT)  = HBED(L,K)*TOXB(L,K,NT)  
              TOXB1(L,K,NT) = TOXB(L,K,NT)  
            ENDDO  
          ENDDO  
        ENDIF  
        IF( ITXBDUT(NT) == 1 )THEN  
          ! *** PMC - BEFORE 2017-06 EFDC ASSUMED ONLY PARTICLE MASS IN FILE.  TYPICAL LABRATORY ANALYSIS ARE TOTAL MASS TOXIC/MASS OF SEDIMENT.
          ! ***       ITXBDUT(NT)=1 NOW REPRESENTS THE INPUT TOTAL TOXB AS MG/KG
          ! *** INPUT UNITS OF MG/KG
          DO K=1,KB  
            DO L=2,LA
              !TOXB(L,K,NT)  = 0.001*TOXB(L,K,NT)*(SEDBT(L,K) + SNDBT(L,K))/TOXPFTB(L,K,NT)   DEPRECATED
              ! *** MG/M2   =  KG/G    MG/KG     (          G/M2         )         
              TOXB(L,K,NT)  = 0.001*TOXB(L,K,NT)*(SEDBT(L,K) + SNDBT(L,K))
              TOXB1(L,K,NT) = TOXB(L,K,NT)  
            ENDDO  
          ENDDO  
        ENDIF  
      ENDDO  

      ! *** DIAGNOSTICS OF INITIALIZATION  
      OPEN(2,FILE=OUTDIR//'TOXBED.DIA')  
      CLOSE(2,STATUS='DELETE')  
      OPEN(2,FILE=OUTDIR//'TOXBED.DIA') 
      DO NT=1,NTOX
        WRITE(2,*) 'TOXIC BED FOR CLASS: ',NT
        DO L=2,LA
          DO K=1,KB
            IF( HBED(L,K) > 0.0 )THEN
              TMP1 = TOXB(L,K,NT)/HBED(L,K)   ! *** MG/M^3
            ELSE
              IF( TOXB(L,K,NT) > 0.0 )THEN
                WRITE(8,*) 'WARNING - TOXIC IC HAS TOXB > 0 FOR ZERO THICKNESS LAYER: ',L,K,NT,TOXB(L,K,NT)
                TOXB(L,K,NT) = 0.0
              ENDIF
              TMP1 = 0.0
            ENDIF
            WRITE(2,2222) IL(L),JL(L),K,TOXPFTB(L,K,NT),TOXB(L,K,NT),TMP1,TOX(L,KSZ(L),NT)
          ENDDO
        ENDDO  
      ENDDO
      CLOSE(2)  
    ENDIF
  ENDIF  
  2222 FORMAT(3I5,7E13.4)  

  ! ***  INITIALIZE FRACTION OF PARTICULATE ORGANIC CARBON IN BED BY FUNCTION
  IF( ISTPOCB == 4 )THEN
    IVAL=0  
    DO NT=1,NTOX  
      IF( ISTOC(NT) >= 2)IVAL=1  
    ENDDO  
    IF( IVAL == 1 )CALL SETFPOCB(0)  
  ENDIF
  
  ! ***  CALCULATE COHESIVE AND NONCOHESIVE VOID RATIOS  
  IF( .NOT. LSEDZLJ )THEN
    DO K=1,KB  
      DO L=2,LA  
        IF( K <= KBT(L) )THEN  
          FVOLSSD = 1.0/(1.0+VDRBED(L,K))  
          FVOLSED = 0.0  
          DO NS=1,NSED  
            FVOLSED = FVOLSED+VFRBED(L,K,NS)  
          ENDDO  
          FVOLSND = 0.0  
          DO NX=1,NSND  
            NS = NSED+NX  
            FVOLSND = FVOLSND+VFRBED(L,K,NS)  
          ENDDO  
          FVOLSED = FVOLSSD*FVOLSED  
          FVOLSND = FVOLSSD*FVOLSND  
          VDRBEDSND(L,K) = SNDVDRD
          IF( FVOLSED > 1.0E-18 )THEN  
            VDRBEDSED(L,K) = ((FVOLSED+FVOLSND)*VDRBED(L,K)-  FVOLSND*SNDVDRD)/FVOLSED  
           ELSE
            VDRBEDSED(L,K) = 0.0  
          ENDIF
         ELSE  
          VDRBEDSND(L,K) = 0.0  
          VDRBEDSED(L,K) = 0.0  
        ENDIF  
      ENDDO  
    ENDDO 

    ! ***  ADD ACTIVE ARMORING LAYER IF NOT PRESENT IN INITIAL OR RESTART  
    IF( ISNDAL == 2 .AND. IALSTUP > 0 .AND. KB > 1 )THEN  
    
      DO L=2,LA  
        KTOPTP = KBT(L)     ! *** NEW PARENT LAYER
        KTOPP1 = KBT(L)+1   ! *** NEW ACTIVE LAYER

        TTHICK=0.
        DO K=1,KB
          TTHICK = TTHICK+HBED(L,K)
        ENDDO

        ! *** MAKE SURE THERE IS SEDIMENT IN THE CELL
        IF( TTHICK > 0. )THEN
          IF( KTOPTP == KB )THEN
            ! *** CASE KBT = KB

            IF( HBED(L,KTOPTP) < HBEDAL )THEN
              ! *** ACTIVE LAYER IS NOT THICK ENOUGH. REPARTITIION KBT AND KBT-1
          
              KT1 = KBT(L)     ! *** ACTIVE LAYER
              KT2 = KBT(L)-1   ! *** PARENT LAYER
              IF( (HBED(L,KT1)+HBED(L,KT2))<HBEDAL )THEN
                WRITE(*,'(A,I5)')'BEDINIT: ISNDAL=2 AND IALSTUP>0: BED LAYER IS NOT THICK ENOUGH.  L = ', L
                STOP 
              ENDIF
              HBEDP = HBEDAL
          
              FRACT1  = HBEDP/HBED(L,KT1)  
              FRACT2  = (HBED(L,KT2)-(HBEDP-HBED(L,KT1)))/HBED(L,KT2)  
              FRACT2C = 1.-FRACT2
          
              ! *** INCREASE TOP LAYER TO HBEDP
              HBED(L,KT1 ) = HBEDP
              SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)  
              SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)  
              DO NS=1,NSED  
                SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)  
                SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)  
              ENDDO
              DO NS=1,NSND  
                SNDB(L,KT1,NS)  = FRACT1*SNDB(L,KT1,NS)  
                SNDB1(L,KT1,NS) = FRACT1*SNDB1(L,KT1,NS)  
              ENDDO
              DO NT=1,NTOX  
                TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)  
                TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)  
              ENDDO
            
              ! *** REDUCE LAYER BELOW                                    
              HBED(L,KT2)  = FRACT2*HBED(L,KT2)  
              SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)  
              SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)  
              DO NS=1,NSED  
                SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)  
                SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)  
              ENDDO
              DO NS=1,NSND  
                SNDB(L,KT2,NS)  = FRACT2*SNDB(L,KT2,NS)  
                SNDB1(L,KT2,NS) = FRACT2*SNDB1(L,KT2,NS)  
              ENDDO
              DO NT=1,NTOX  
                TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)  
                TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)  
              ENDDO
            ELSE
              ! *** ACTIVE LAYER IS TOO THICK. REPARTITIION KBT AND KBT-1
              KT1 = KBT(L)     ! *** ACTIVE LAYER
              KT2 = KBT(L)-1   ! *** PARENT LAYER
          
              HBEDP = HBEDAL
          
              FRACT1  = HBEDP/HBED(L,KT1)  
              FRACT2  = (HBED(L,KT2)+(HBED(L,KT1)-HBEDP))/HBED(L,KT2)  
              FRACT2C = 1.-FRACT2
          
              ! *** DECREASE TOP LAYER TO HBEDP
              HBED(L,KT1)  = HBEDP
              SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)  
              SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)  
              DO NS=1,NSED  
                SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)  
                SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)  
              ENDDO
              DO NS=1,NSND  
                SNDB(L,KT1,NS)  = FRACT1*SNDB(L,KT1,NS)  
                SNDB1(L,KT1,NS) = FRACT1*SNDB1(L,KT1,NS)  
              ENDDO
              DO NT=1,NTOX  
                TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)  
                TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)  
              ENDDO
          
              ! *** INCREASE THE LAYER BELOW                                    
              HBED(L,KT2)  = FRACT2*HBED(L,KT2)  
              SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)  
              SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)  
              DO NS=1,NSED  
                SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)  
                SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)  
              ENDDO
              DO NS=1,NSND  
                SNDB(L,KT2,NS)  = FRACT2*SNDB(L,KT2,NS)  
                SNDB1(L,KT2,NS) = FRACT2*SNDB1(L,KT2,NS)  
              ENDDO
              DO NT=1,NTOX  
                TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)  
                TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)  
              ENDDO
        
            ENDIF
      
          ELSE
      
            ! *** STANDARD CASE OF KBT<KB
            IF( HBED(L,KTOPTP) < HBEDAL )THEN
              ! *** PARENT LAYER IS NOT THICK ENOUGH. REPARTITIION KBT AND KBT-1

              KT1 = KBT(L)     ! *** PARENT LAYER
              KT2 = KBT(L)-1   ! *** LAYER BELOW PARENT
              IF( KT2 > 0 )THEN
                HBEDP = MIN(HBEDAL*2.,HBED(L,KT2))
            
                FRACT1  = HBEDP/HBED(L,KT1)  
                FRACT2  = (HBED(L,KT2)-(HBEDP-HBED(L,KT1)))/HBED(L,KT2)  
                FRACT2C = 1.-FRACT2
            
                ! *** INCREASE TOP LAYER TO HBEDP
                HBED(L,KT1) = HBEDP
                IF( ISTRAN(6) > 0 )THEN
                  SEDBT(L,KT1) = FRACT1*SEDBT(L,KT1)  
                  DO NS=1,NSED  
                    SEDB(L,KT1,NS)  = FRACT1*SEDB(L,KT1,NS)  
                    SEDB1(L,KT1,NS) = FRACT1*SEDB1(L,KT1,NS)  
                  ENDDO
                ENDIF
                IF( ISTRAN(7) > 0 )THEN
                  SNDBT(L,KT1) = FRACT1*SNDBT(L,KT1)  
                  DO NS=1,NSND  
                    SNDB(L,KT1,NS)  = FRACT1*SNDB(L,KT1,NS)  
                    SNDB1(L,KT1,NS) = FRACT1*SNDB1(L,KT1,NS)  
                  ENDDO
                ENDIF
                IF( ISTRAN(5) > 0 )THEN
                  DO NT=1,NTOX  
                    TOXB(L,KT1,NT)  = FRACT1*TOXB(L,KT1,NT)  
                    TOXB1(L,KT1,NT) = FRACT1*TOXB1(L,KT1,NT)  
                  ENDDO
                ENDIF
                
                ! *** REDUCE LAYER BELOW                                    
                HBED(L,KT2) = FRACT2*HBED(L,KT2)  
                IF( ISTRAN(6) > 0 )THEN
                  SEDBT(L,KT2) = FRACT2*SEDBT(L,KT2)  
                  DO NS=1,NSED  
                    SEDB(L,KT2,NS)  = FRACT2*SEDB(L,KT2,NS)  
                    SEDB1(L,KT2,NS) = FRACT2*SEDB1(L,KT2,NS)  
                  ENDDO
                ENDIF
                IF( ISTRAN(7) > 0 )THEN
                  SNDBT(L,KT2) = FRACT2*SNDBT(L,KT2)  
                  DO NS=1,NSND  
                    SNDB(L,KT2,NS)  = FRACT2*SNDB(L,KT2,NS)  
                    SNDB1(L,KT2,NS) = FRACT2*SNDB1(L,KT2,NS)  
                  ENDDO
                ENDIF
                IF( ISTRAN(5) > 0 )THEN
                  DO NT=1,NTOX  
                    TOXB(L,KT2,NT)  = FRACT2*TOXB(L,KT2,NT)  
                    TOXB1(L,KT2,NT) = FRACT2*TOXB1(L,KT2,NT)  
                  ENDDO
                ENDIF
              ENDIF
            ENDIF
        
            ! *** SPLIT THE PARENT LAYER INTO THE PARENT (KTOPTP) AND THE ACTIVE (KTOPP1)
            FRACACT(L) = HBEDAL/HBED(L,KTOPTP)  
            FRACPAR(L) = (HBED(L,KTOPTP)-HBEDAL)/HBED(L,KTOPTP)  
            HBED(L,KTOPP1) = FRACACT(L)*HBED(L,KTOPTP)  ! ACTIVE
            HBED(L,KTOPTP) = FRACPAR(L)*HBED(L,KTOPTP)  ! PARENT
            PORBED(L,KTOPP1) = PORBED(L,KTOPTP)  
            PORBED1(L,KTOPP1) = PORBED1(L,KTOPTP)  
            VDRBED(L,KTOPP1) = VDRBED(L,KTOPTP)  
            VDRBED1(L,KTOPP1) = VDRBED1(L,KTOPTP)  
            BDENBED(L,KTOPP1) = BDENBED(L,KTOPTP)  
            BDENBED1(L,KTOPP1) = BDENBED1(L,KTOPTP)  
            SEDBT(L,KTOPP1) = FRACACT(L)*SEDBT(L,KTOPTP)  
            SEDBT(L,KTOPTP) = FRACPAR(L)*SEDBT(L,KTOPTP)  
            SNDBT(L,KTOPP1) = FRACACT(L)*SNDBT(L,KTOPTP)  
            SNDBT(L,KTOPTP) = FRACPAR(L)*SNDBT(L,KTOPTP)  
            STDOCB(L,KTOPP1) = STDOCB(L,KTOPTP)  
            STPOCB(L,KTOPP1) = STPOCB(L,KTOPTP)  

            IF( ISTRAN(6) > 0 )THEN
              DO NS=1,NSED  
                SEDB(L,KTOPP1,NS) = FRACACT(L)*SEDB(L,KTOPTP,NS)  
                SEDB1(L,KTOPP1,NS) = FRACACT(L)*SEDB1(L,KTOPTP,NS)  
                SEDB(L,KTOPTP,NS) = FRACPAR(L)*SEDB(L,KTOPTP,NS)  
                SEDB1(L,KTOPTP,NS) = FRACPAR(L)*SEDB1(L,KTOPTP,NS)  
                STFPOCB(L,KTOPP1,NS) = STFPOCB(L,KTOPTP,NS)  
              ENDDO  
            ENDIF
            IF( ISTRAN(7) > 0 )THEN
              DO NS=1,NSND  
                NX=NSED+NS  
                SNDB(L,KTOPP1,NS) = FRACACT(L)*SNDB(L,KTOPTP,NS)  
                SNDB1(L,KTOPP1,NS) = FRACACT(L)*SNDB1(L,KTOPTP,NS)  
                SNDB(L,KTOPTP,NS) = FRACPAR(L)*SNDB(L,KTOPTP,NS)  
                SNDB1(L,KTOPTP,NS) = FRACPAR(L)*SNDB1(L,KTOPTP,NS)  
                STFPOCB(L,KTOPP1,NX) = STFPOCB(L,KTOPTP,NX)  
              ENDDO  
            ENDIF
            IF( ISTRAN(5) > 0 )THEN
              DO NT=1,NTOX  
                TOXB(L,KTOPP1,NT) = FRACACT(L)*TOXB(L,KTOPTP,NT)  
                TOXB1(L,KTOPP1,NT) = FRACACT(L)*TOXB1(L,KTOPTP,NT)  
                TOXB(L,KTOPTP,NT) = FRACPAR(L)*TOXB(L,KTOPTP,NT)  
                TOXB1(L,KTOPTP,NT) = FRACPAR(L)*TOXB1(L,KTOPTP,NT)  
                TOXPFTB(L,KTOPP1,NT) = TOXPFTB(L,KTOPTP,NT)  
              ENDDO  
              DO NT=1,NTOX  
                DO NS=1,NSED+NSND+2  
                  TOXPFB(L,KTOPP1,NS,NT) = TOXPFB(L,KTOPTP,NS,NT)  
                ENDDO  
              ENDDO
            ENDIF
            
            ! *** UPDATE TOP LAYER
            KBT(L) = KBT(L)+1          
      
          ENDIF  ! *** END OF KBT<KB CHECK
        ENDIF    ! *** END OF TOTAL THICKNESS CHECK: TTHICK>0
      ENDDO  ! *** END OF LA LOOP
    
    ENDIF    ! *** END OF REPARTITIONING ACTIVE ARMORING LAYER AT STARTUP IF REQUESTED

    ! ***********************************************************************!
    ! ***  ADJUST POROSITY AND VOID RATIO FOR IBMECH == 99
    IF( IBMECH == 99 .AND. .NOT. LSEDZLJ )THEN
      DO K=1,KB
        DO L=2,LA
          FRACCOH(L,K) = 0.0
          FRACNON(L,K) = 0.0
        ENDDO
      ENDDO

      DO NS=1,NSED
        DO K=1,KB
          DO L=2,LA
            IF( K <= KBT(L) )THEN
              FRACCOH(L,K) = FRACCOH(L,K)+VFRBED(L,K,NS)
              FRACNON(L,K) = FRACNON(L,K)+VFRBED1(L,K,NS)
            ENDIF
          ENDDO
        ENDDO
      ENDDO

      DO K=1,KB
        DO L=2,LA
          IF( K <= KBT(L) )THEN
            PORBED(L,K) = BMECH1*(FRACCOH(L,K)**BMECH2)+BMECH3
            PORBED1(L,K) = BMECH1*(FRACNON(L,K)**BMECH2)+BMECH3
            VDRBED(L,K) = PORBED(L,K)/(1.-PORBED(L,K))
            VDRBED1(L,K) = PORBED1(L,K)/(1.-PORBED1(L,K))
          ENDIF
        ENDDO
      ENDDO

    ENDIF

  ENDIF   ! *** END OF LSEDZLJ BYPASS

  IF( .NOT. LSEDZLJ )THEN
    !**********************************************************************C                                                
    ! *** SET MEAN D50 AND D90                                                                                              
    IF( ISTRAN(7) >= 1 )THEN
      DO K=1,KB
        DO L=2,LA
          SEDDIA50(L,K)=0.
          SEDDIA90(L,K)=0.
          SNDBT(L,K)=0.
        ENDDO
      ENDDO
      DO NX=1,NSND
        NS=NSED+NX
        DO K=1,KB
          DO L=2,LA
            SEDDIA50(L,K)=SEDDIA50(L,K)+SNDB(L,K,NX)*LOG(SEDDIA(NS))
            SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)
          ENDDO
        ENDDO
      ENDDO
      DO K=1,KB
        DO L=2,LA
          IF( SNDBT(L,K) > 0. )THEN
            SEDDIA50(L,K)=SEDDIA50(L,K)/SNDBT(L,K)
          ENDIF
        ENDDO
      ENDDO
                            
      DO K=1,KB
        DO L=2,LA
          IF( SNDBT(L,K) > 0.1 )THEN  ! *** Handle cases when very little sand but not zero (SNDBT is in g/m2)
            SEDDIA50(L,K)=EXP(SEDDIA50(L,K))
          ELSE
            SEDDIA50(L,K)=SEDDIA(NSED+1)
          ENDIF
        ENDDO
      ENDDO
                            
    ENDIF  ! *** END OF NON-COHESIVE GRAIN SIZE CALCS
  ENDIF    ! *** END OF LSEDZLJ BYPASS

  ! ***********************************************************************************************************
  1000 CONTINUE     ! JUMP FROM THE HOTSTART (SEDZLJ INITIALIZATION CONTINUES FROM THE HOTSTART SECTION)

  DO K=1,KB
    DO L=2,LA
      IF( (K <= KBT(L) .AND. .NOT. LSEDZLJ) .OR. (K >= KBT(L) .AND. LSEDZLJ) )THEN
        SDENAVG(L,K) = (BDENBED(L,K)-1000.0*PORBED(L,K))/(1.0-PORBED(L,K))
      ELSE
        SDENAVG(L,K) = (BDENBED(L,KBT(L))-1000.0*PORBED(L,KBT(L)))/(1.0-PORBED(L,KBT(L)))
      ENDIF
      IF( SDENAVG(L,K) <= 0. )SDENAVG(L,K) = 0.0
    ENDDO
  ENDDO
  
  ! ***********************************************************************************************************
  ! ***  WRITE DIAGNOSTIC FILES FOR BED INITIALIZATION  
  OPEN(1,FILE=OUTDIR//'BEDINIT.SED')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.SED')  
  WRITE(1,111)  
  DO L=2,LA  
    WRITE(1,101)IL(L),JL(L),(SEDB(L,K,1),K=1,KB)  
    IF( NSED > 1 )THEN  
      DO NX=2,NSED  
        WRITE(1,102)(SEDB(L,K,NX),K=1,KB)  
      END DO  
    ENDIF  
  ENDDO  
  CLOSE(1)
  
  OPEN(1,FILE=OUTDIR//'BEDINIT.SND')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.SND')  
  WRITE(1,112)  
  DO L=2,LA  
    WRITE(1,101)IL(L),JL(L),(SNDB(L,K,1),K=1,KB)  
    IF( NSND > 1 )THEN  
      DO NX=2,NSND  
        WRITE(1,102)(SNDB(L,K,NX),K=1,KB)  
      END DO  
    ENDIF  
  ENDDO  
  
  CLOSE(1)  
  OPEN(1,FILE=OUTDIR//'BEDINIT.VDR')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.VDR')  
  WRITE(1,113)  
  DO L=2,LA  
    WRITE(1,101)IL(L),JL(L),(VDRBED(L,K),K=1,KB)  
  ENDDO  
  CLOSE(1)  

  OPEN(1,FILE=OUTDIR//'BEDINIT.POR')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.POR')  
  WRITE(1,114)  
  DO L=2,LA  
    WRITE(1,101)IL(L),JL(L),(PORBED(L,K),K=1,KB)  
  ENDDO  
  CLOSE(1) 
     
  OPEN(1,FILE=OUTDIR//'BEDINIT.ZHB')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.ZHB')  
  WRITE(1,115)  
  DO L=2,LA  
    WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),(HBED(L,K),K=1,KB)  
  ENDDO  
  CLOSE(1)  
    
  OPEN(1,FILE=OUTDIR//'BEDINIT.BDN')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.BDN')  
  WRITE(1,116)  
  DO L=2,LA  
    WRITE(1,101)IL(L),JL(L),(BDENBED(L,K),K=1,KB)  
  ENDDO  
  CLOSE(1)  
    
  OPEN(1,FILE=OUTDIR//'BEDINIT.ELV')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.ELV')  
  WRITE(1,117)  
  DO L=2,LA  
    SURF=HP(L)+BELV(L)  
    WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),BELV(L),HP(L),SURF  
  ENDDO  
  CLOSE(1)  
    
  OPEN(1,FILE=OUTDIR//'BEDINIT.TOX')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.TOX')  
  DO NT=1,NTOX  
    WRITE(1,118)NT  
    DO L=2,LA  
      WRITE(1,101)IL(L),JL(L),(TOXB(L,K,NT),K=1,KB)  
    ENDDO  
  ENDDO  
  CLOSE(1)  
    
  OPEN(1,FILE=OUTDIR//'BEDINIT.VRS')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINIT.VRS')  
  DO L=2,LA  
    WRITE(1,191)IL(L),JL(L),(VDRBED(L,K),K=1,KB)  
    WRITE(1,192,IOSTAT=ISO) (VDRBEDSED(L,K),K=1,KB)  
  ENDDO  
  CLOSE(1)  
    
  OPEN(1,FILE=OUTDIR//'BEDINITC.VVF')  
  CLOSE(1,STATUS='DELETE')  
  OPEN(1,FILE=OUTDIR//'BEDINITC.VVF')  
  OPEN(2,FILE=OUTDIR//'BEDINITF.VVF')  
  CLOSE(2,STATUS='DELETE')  
  OPEN(2,FILE=OUTDIR//'BEDINITF.VVF')  
  DO L=2,LA  
    K=KBT(L)  
    IF( HMP(L) > 0.05 )THEN  
      WRITE(1,191,IOSTAT=ISO)IL(L),JL(L),VFRBED(L,K,1),PORBED(L,K),VDRBED(L,K),  VDRBEDSED(L,K)  
    ELSE  
      WRITE(2,191)IL(L),JL(L),VFRBED(L,K,1),PORBED(L,K),VDRBED(L,K),  VDRBEDSED(L,K)  
    ENDIF  
  ENDDO 
  ! *** END DIAGNOSTIC FILE DUMP

  ! *** SAVE STARTUP BED PARAMETERS FOR LATER  
  DO K=1,KB  
    DO L=2,LA  
      IF( K <= KBT(L) )THEN  
        VDRBED2(L,K)=VDRBED(L,K)  
      ELSE  
        VDRBED2(L,K)=VDRBED(L,KBT(L))  
      ENDIF  
    ENDDO  
  ENDDO  

  DO L=2,LA  
    DO K=1,KB  
      BEDTHKSV(L,K)=HBED(L,K)  
      BEDPORSV(L,K)=PORBED(L,K)  
      BEDVDRSV(L,K)=VDRBED(L,K)  
      BEDBKDSV(L,K)=BDENBED(L,K)  
    ENDDO  
  ENDDO  

  CLOSE(1)  
  CLOSE(2)  
  
  191 FORMAT(2I5,18F10.3)  
  192 FORMAT(10X,18F10.3)  
  101 FORMAT(2I5,18E13.5)  
  102 FORMAT(10X,18E13.5)  
  111 FORMAT('   IL   JL    SEDBT(K=1,KB)')  
  112 FORMAT('   IL   JL    SNDBT(K=1,KB)')  
  113 FORMAT('   IL   JL    VRDBED(K=1,KB)')  
  114 FORMAT('   IL   JL    PORBED(K=1,KB)')  
  115 FORMAT('   IL   JL    ZBEDB        HBEDT        HBED(K=1,KB)')  
  116 FORMAT('   IL   JL    BDENBED(K=1,KB)')  
  117 FORMAT('   IL   JL    ZBEDB        HBEDT        BELV',  '        HWCOL        SELV')  
  118 FORMAT('   IL   JL    TOXB(K=1,KB,NT)  NT = ',I5)  
  RETURN  

END  

