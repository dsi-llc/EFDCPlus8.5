SUBROUTINE CALTOX

  ! ***  SUBROUTINE CALTOX CALCULATES CONTAMINANT TRANSPORT.
  ! ***  IT IS CALLED FROM SSEDTOX
  ! ***  USED BY BOTH THE STANDARD EFDC SEDTRAN MODULE AND SEDZLJ
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2014-08-26        PAUL M. CRAIG    Updated to allow for any mix of ISTOC
  ! 2012-12-05        PAUL M. CRAIG    Updated OMP
  ! 2012-10-02        Dang H Chung     Added OMP
  !**********************************************************************!

  USE GLOBAL
  IMPLICIT NONE

  INTEGER :: ND,LF,LL,LP,L,IVAL,NT,NS,NX,K,NFD,LUTMP,LDTMP,KTOPM1,NP
  INTEGER :: LBANK,LCHAN,KBANK,LE,LW,LS,LN,NSB,KCHAN,KTOPTP
  INTEGER, SAVE :: NDUMP
  
  REAL :: TMPEXP,TMPVAL,TMPTOXB,TMPTOXC,TMPTOXE,TMPTOXW,TMPTOXN,TMPTOXS,TOXFLUX,CBLTOXTMP
  REAL :: AA11,BB11,BB22,SNDFEFF,CLEFT,CRIGHT,WVEL
  REAL :: FTPOS,FTNEG,HBEDMIN0,SORBMIN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TOXFPA
  
  IF(  .NOT. ALLOCATED(TOXFPA) )THEN
    ALLOCATE(TOXFPA(LCM))
    TOXFPA=0.0
    NDUMP = 0
  ENDIF

  IF( LSEDZLJ )THEN
      HBEDMIN0 = 1E-9
  ELSE
      HBEDMIN0 = 1E-4
  ENDIF

  !**********************************************************************CC
  ! **  UPDATE FRACTION OF PARTICULATE ORGANIC CARBON IN BED
  IF( ISTPOCB == 4 )THEN
    IVAL=0
    DO NT=1,NTOX
      IF( ISTOC(NT) >= 2)IVAL=1
    ENDDO
  
    IF( IVAL == 1 )THEN
      ! *** HOUSATONIC
      CALL SETFPOCB(1)
    ENDIF
  ENDIF
  
  !**********************************************************************C
  ! **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN WATER COLUMN
  ! **
  ! **  TOXPFW(L,K,NS,NT) = PARTICULATE FRACTION IN WATER COLUMN
  ! **  TOXPARW(NS,NT) = PARTITION COEFFICIENT IN WATER COLUMN
  ! **  TOXPFTW(L,K,NT) = TOTAL PARTICULATE FRACTION IN WATER COLUMN
  !                       USED AS TEMPORARY VARIBLE IN THIS AND
  !                       FOLLOWING CODE BLOCK
  NFD=NSED+NSND+1

  DO NT=1,NTOX
    TOXBLB(NT)=0.0
  ENDDO

  NP = PRECISION(SORBMIN)
  SORBMIN = 0.99  ! 1. - 1./10**NP
  
  ! *** OMP LOOP OVER THE TOXICS COMPUTATIONS (SETTLING, DEPOSTION AND EROSION ONLY)
  !$OMP PARALLEL DEFAULT(SHARED)
 
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,NT,NS,NX,LE,LW,LN,LS,TMPEXP,TMPVAL,SNDFEFF)
  DO ND=1,NDM
    LF=2+(ND-1)*LDM
    LL=MIN(LF+LDM-1,LA)

    !**********************************************************************C
    ! **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN WATER COLUMN
    ! **
    ! **  TOXPFW(L,K,NS,NT) = PARTICULATE FRACTION IN WATER COLUMN
    ! **  TOXPARW(NS,NT)    = PARTITION COEFFICIENT IN WATER COLUMN
    ! **  TOXPFTW(L,K,NT)   = TOTAL PARTICULATE FRACTION IN WATER COLUMN

    DO NT=1,NTOX
      DO NS=1,NSP2(NT)
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TOXPFW(L,K,NS,NT)=0.
          ENDDO
        ENDDO
      ENDDO
    ENDDO

    DO NT=1,NTOX
      ! *** PARTITION TO COHESIVE
      IF( ISTRAN(6) >= 1 )THEN
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO NS=1,NSED
            IF( ITXPARW(NS,NT) == 0 )THEN   ! *** NON-SOLIDS BASED PARTITIONING
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  ! ***            =      G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SED(L,K,NS)*STFPOCW(L,K,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF
            IF( ITXPARW(NS,NT) == 1 )THEN   ! *** SOLIDS BASED PARTITIONING
              TMPEXP=CONPARW(NS,NT)
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  TMPVAL=1.
                  IF( SED(L,K,NS) > 0.) TMPVAL=SED(L,K,NS)**TMPEXP
                  ! ***             =           G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SED(L,K,NS)*STFPOCW(L,K,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN
          ! *** Kd APPROACH
          DO NS=1,NSED
            IF( ITXPARW(NS,NT) == 0 )THEN
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  ! ***             =     G/M3      L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SED(L,K,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF
            IF( ITXPARW(NS,NT) == 1 )THEN
              TMPEXP=CONPARW(NS,NT)
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  TMPVAL=1.
                  IF( SED(L,K,NS) > 0.) TMPVAL=SED(L,K,NS)**TMPEXP
                  ! ***             =          G/M3       L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SED(L,K,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      ! *** PARTITION TO NONCOHESIVE
      IF( ISTRAN(7) >= 1 )THEN
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO NX=1,NSND
            NS=NX+NSED
            IF( ITXPARW(NS,NT) == 0 )THEN
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  ! ***             =     G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SND(L,K,NX)*STFPOCW(L,K,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF

            IF( ITXPARW(NS,NT) == 1 )THEN
              TMPEXP=CONPARW(NS,NT)
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  TMPVAL=1.
                  IF( SND(L,K,NX) > 0.) TMPVAL=SND(L,K,NX)**TMPEXP
                  ! ***             =            G/M3                     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SND(L,K,NX) *STFPOCW(L,K,NS)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN
          ! *** Kd APPROACH
          DO NX=1,NSND
            NS=NX+NSED
            IF( ITXPARW(NS,NT) == 0 )THEN
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  ! ***             =     G/M3     L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = SND(L,K,NX)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF

            IF( ITXPARW(NS,NT) == 1 )THEN
              TMPEXP=CONPARW(NS,NT)
              DO K=1,KC
                DO LP=1,LLWET(K,ND)
                  L=LKWET(LP,K,ND)  
                  TMPVAL=1.
                  IF( SND(L,K,NX) > 0.) TMPVAL=SND(L,K,NX)**TMPEXP
                  ! ***             =           G/M3      L/MG (M3/G)   
                  TOXPFW(L,K,NS,NT) = TMPVAL*SND(L,K,NX)*TOXPARW(NS,NT)
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      ! *** PARTITION (COMPLEX TO DISSOLVED ORGANIC CARBON)
      IF( ISTOC(NT) == 1 .OR. ISTOC(NT) == 2 )THEN
        NS=1+NSED+NSND
        IF( ITXPARWC(1,NT) == 0 )THEN
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              ! ***                DOC(G/M3)   L/MG (M3/G)
              TOXPFW(L,K,NS,NT) = STDOCW(L,K)*TOXPARWC(1,NT)
            ENDDO
          ENDDO
        ENDIF

        IF( ITXPARWC(1,NT) == 1 )THEN
          TMPEXP=CONPARWC(1,NT)
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              TMPVAL=1.
              IF( STDOCW(L,K) > 0.) TMPVAL=STDOCW(L,K)**TMPEXP
              ! ***                       DOC(G/M3)   L/MG (M3/G)
              TOXPFW(L,K,NS,NT) = TMPVAL*STDOCW(L,K)*TOXPARWC(1,NT)
            ENDDO
          ENDDO
        ENDIF
      ENDIF

      ! *** PARTITION TO PARTICULATE ORGANIC CARBON
      IF( ISTOC(NT) == 1 )THEN
        NS=2+NSED+NSND
        IF( ITXPARWC(2,NT) == 0 )THEN
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              ! ***                POC(G/M3)   L/MG (M3/G)
              TOXPFW(L,K,NS,NT) = STPOCW(L,K)*TOXPARWC(2,NT)
            ENDDO
          ENDDO
        ENDIF

        IF( ITXPARW(NS,NT) == 1 )THEN
          TMPEXP=CONPARW(NS,NT)
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              TMPVAL=1.
              IF( STPOCW(L,K) > 0.) TMPVAL=STPOCW(L,K)**TMPEXP
              ! ***                       POC(G/M3)   L/MG (M3/G)
              TOXPFW(L,K,NS,NT) = TMPVAL*STPOCW(L,K)*TOXPARWC(2,NT)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDDO   ! *** NTOX

    ! ** TOXPFTW IS TEMPORARILY USED TO STORE TOTAL SORBED (DIMENSIONLESS)
    DO NT=1,NTOX
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          TOXPFTW(L,K,NT)=0.
        ENDDO
      ENDDO
      DO NS=1,NSP2(NT)
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TOXPFTW(L,K,NT) = TOXPFTW(L,K,NT) + TOXPFW(L,K,NS,NT)
          ENDDO
        ENDDO
      ENDDO
    ENDDO   ! *** NTOX
    
    ! ** COMPUTE THE SORBED FRACTION (TOXPFW) (DIMENSIONLESS)
    DO NT=1,NTOX
      DO NS=1,NSP2(NT)
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TOXPFW(L,K,NS,NT) = TOXPFW(L,K,NS,NT)/(1.+TOXPFTW(L,K,NT))
          ENDDO
        ENDDO
      ENDDO
    ENDDO   ! *** NTOX

    IF( ISTMSR >= 1 .OR. ISFDCH > 0 )THEN
      ! *** ONLY NEEDED FOR TIME SERIES USING TMSR OR FOODCHAIN
      DO NT=1,NTOX
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TOXFDFW(L,K,NT) = 1./(1.+TOXPFTW(L,K,NT))
            TOXCDFW(L,K,NT) = TOXPFW(L,K,NFD,NT)
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    !**********************************************************************C
    !
    ! **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
    !
    ! **  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED 
    ! **  TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED (INTERIM UNITS: METERS)
    CALL CALTOXB_FRACTIONS(LF,LL)

    ! ******************************************************************************
    ! ******************************************************************************
    !
    ! *** CALCULATE PARTICULATE TOXIC CONTAMINANT SETTLING AND BED EXCHANGE FLUXES

    ! *** TOXF(L,1:KS,NT) = TOXIC CONTAMINANT SETTLING AND BED EXCHANGE FLUX 
    ! ***                    (+) UPWARD FLUX,  (-) DOWNWARD FLUX   (M/S)
    ! *** TOXF(L,0,NT)    = TOXIC CONTAMINANT DEPOSITIONAL FLUX    (M/S)
    ! *** TOXFB(L,KBT,NT) = TOXIC CONTAMINANT EROSIONAL FLUX       (FINAL 1/S, INTERIM M/S)
    
    DO NT=1,NTOX
      DO K=0,KC
        DO L=LF,LL
          TOXF(L,K,NT)=0.
        ENDDO
      ENDDO
    ENDDO

    ! ******************************************************************************
    ! *** TOXICS VERTICAL FLUX (TOXF) IN THE WATER COLUMN (K=1,KS)
    ! *** SEDF (G/M2/S) IS ALWAYS NEGATIVE IN THE WATER COLUMN
    IF( KC >= 2 )THEN
      DO NT=1,NTOX
        ! *** PARTICLE COHESIVE FLUX
        IF( ISTRAN(6) >= 1 )THEN
          IF( ISTOC(NT) > 1 )THEN
            ! *** fPOC BASED
            DO NS=1,NSED
              IF( ITXPARW(NS,NT) == 0 )THEN
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND)  
                    ! *** M/S          M/S           G/M2/S                         M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SEDF(L,K,NS)*STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ELSEIF( ITXPARW(NS,NT) == 1 )THEN
                TMPEXP=CONPARW(NS,NT)
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND)  
                    TMPVAL = 1.
                    IF( SED(L,K+1,NS) > 0.) TMPVAL = SED(L,K+1,NS)**TMPEXP
                    ! *** M/S          M/S                  G/M2/S                         M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SEDF(L,K,NS)*STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSEIF( ISTOC(NT) == 0 )THEN
            ! *** Kd APPROACH
            DO NS=1,NSED
              IF( ITXPARW(NS,NT) == 0 )THEN
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND) 
                    ! *** M/S          M/S            G/M2/S       M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SEDF(L,K,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ELSEIF( ITXPARW(NS,NT) == 1 )THEN
                TMPEXP=CONPARW(NS,NT)
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND)  
                    TMPVAL = 1.
                    IF( SED(L,K+1,NS) > 0.) TMPVAL = SED(L,K+1,NS)**TMPEXP
                    ! *** M/S          M/S                   G/M2/S       M3/G
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SEDF(L,K,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF            
        ENDIF

        ! *** PARTICLE NONCOHESIVE FLUX
        IF( ISTRAN(7) >= 1 )THEN
          IF( ISTOC(NT) > 1 )THEN
            ! *** fPOC BASED
            DO NX=1,NSND
              NS=NX+NSED
              IF( ITXPARW(NS,NT) == 0 )THEN
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND)  
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SNDF(L,K,NX)*STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
              IF( ITXPARW(NS,NT) == 1 )THEN
                TMPEXP=CONPARW(NS,NT)
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND)  
                    TMPVAL = 1.
                    IF( SND(L,K+1,NX) > 0.) TMPVAL = SND(L,K+1,NX)**TMPEXP
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SNDF(L,K,NX)*STFPOCW(L,K+1,NS)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ELSEIF( ISTOC(NT) == 0 )THEN
            ! *** Kd APPROACH
            DO NX=1,NSND
              NS=NX+NSED
              IF( ITXPARW(NS,NT) == 0 )THEN
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND)  
                    TOXF(L,K,NT) = TOXF(L,K,NT) + SNDF(L,K,NX)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
              IF( ITXPARW(NS,NT) == 1 )THEN
                TMPEXP=CONPARW(NS,NT)
                DO K=1,KS
                  DO LP=1,LLWET(K,ND)
                    L=LKWET(LP,K,ND)  
                    TMPVAL = 1.
                    IF( SND(L,K+1,NX) > 0.) TMPVAL = SND(L,K+1,NX)**TMPEXP
                    TOXF(L,K,NT) = TOXF(L,K,NT) + TMPVAL*SNDF(L,K,NX)*TOXPARW(NS,NT)
                  ENDDO
                ENDDO
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO   ! *** NTOX

      DO NT=1,NTOX
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TOXF(L,K,NT) = TOXF(L,K,NT)/(1.+TOXPFTW(L,K+1,NT))
          ENDDO
        ENDDO
      ENDDO
    ENDIF    ! *** KC >= 2

    ! ******************************************************************************
    ! *** DEPOSITIONAL TOXICS FLUX (TOXF) AT THE SEDIMENT BED/WATER COLUMN INTERFACE 
    ! *** FROM THE WATER COLUMMN ONLY
    DO NT=1,NTOX
      ! *** PARTICLE COHESIVE DEPOSITIONAL FLUX,  (SEDF) < 0 
      IF( ISTRAN(6) >= 1 )THEN
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO NS=1,NSED
            IF( ITXPARW(NS,NT) == 0 )THEN
              DO L=LF,LL
                ! *** M/S         M/S                G/M2/S                                M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SEDF(L,0,NS),0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
            IF( ITXPARW(NS,NT) == 1 )THEN
              TMPEXP=CONPARW(NS,NT)
              DO L=LF,LL
                TMPVAL=1.
                IF( SED(L,KSZ(L),NS) > 0.) TMPVAL = SED(L,KSZ(L),NS)**TMPEXP
                ! *** M/S         M/S                     G/M2/S                               M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SEDF(L,0,NS),0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN
          ! *** Kd APPROACH
          DO NS=1,NSED
            IF( ITXPARW(NS,NT) == 0 )THEN
              DO L=LF,LL
                ! *** M/S         M/S                G/M2/S           M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SEDF(L,0,NS),0.)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
            IF( ITXPARW(NS,NT) == 1 )THEN
              TMPEXP=CONPARW(NS,NT)
              DO L=LF,LL
                TMPVAL=1.
                IF( SED(L,KSZ(L),NS) > 0.) TMPVAL = SED(L,KSZ(L),NS)**TMPEXP
                ! *** M/S         M/S                     G/M2/S           M3/G
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SEDF(L,0,NS),0.)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      ! *** PARTICLE NONCOHESIVE DEPOSITIONAL FLUX,  (SNDF-SNDFBL) < 0 
      IF( ISTRAN(7) >= 1 )THEN
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO NX=1,NSND
            NS=NX+NSED
            IF( ITXPARW(NS,NT) == 0 )THEN
              DO L=LF,LL
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SNDFEFF,0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
            IF( ITXPARW(NS,NT) == 1 )THEN
              TMPEXP=CONPARW(NS,NT)
              DO L=LF,LL
                TMPVAL = 1.
                IF( SND(L,KSZ(L),NX) > 0.) TMPVAL = SND(L,KSZ(L),NX)**TMPEXP
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SNDFEFF,0.)*STFPOCW(L,KSZ(L),NS)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN
          ! *** Kd APPROACH
          DO NX=1,NSND
            NS=NX+NSED
            IF( ITXPARW(NS,NT) == 0 )THEN
              DO L=LF,LL
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + MIN(SNDFEFF,0.)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
            IF( ITXPARW(NS,NT) == 1 )THEN
              TMPEXP=CONPARW(NS,NT)
              DO L=LF,LL
                TMPVAL = 1.
                IF( SND(L,KSZ(L),NX) > 0.) TMPVAL = SND(L,KSZ(L),NX)**TMPEXP
                SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
                TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*MIN(SNDFEFF,0.)*TOXPARW(NS,NT)
              ENDDO
            ENDIF
          ENDDO

        ENDIF
      ENDIF
    ENDDO   ! *** NTOX

    ! *** FINALIZE THE DEPOSITIONAL TOXICS FLUX (M/S)
    DO NT=1,NTOX
      DO L=LF,LL
        TOXF(L,0,NT) = TOXF(L,0,NT)/(1.+TOXPFTW(L,KSZ(L),NT))
      ENDDO
    ENDDO

    ! ******************************************************************************
    ! *** EROSIONAL TOXICS FLUX (TOXFB) AT THE SEDIMENT BED/WATER COLUMN INTERFACE
    ! *** WHEN SOLIDS FLUX IS > 0
    DO NT=1,NTOX
      DO L=LF,LL
        TOXFB(L,KBT(L),NT)=0.   ! *** INTERIM UNITS: M/S, FINAL UNITS: 1/S
      ENDDO
    ENDDO

    DO NT=1,NTOX
      ! *** PARTICLE COHESIVE EROSIONAL FLUX,  (SEDF) > 0
      IF( ISTRAN(6) >= 1 )THEN
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO NS=1,NSED
            DO L=LF,LL
              ! ***     M/S               M/S                 G/M2/S                               M3/G
              TOXFB(L,KBT(L),NT) = TOXFB(L,KBT(L),NT) + MAX(SEDF(L,0,NS),0.)*STFPOCB(L,KBT(L),NS)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN
          ! *** Kd APPROACH
          DO NS=1,NSED
            DO L=LF,LL
              ! ***     M/S               M/S                  G/M2/S          M3/G
              TOXFB(L,KBT(L),NT) = TOXFB(L,KBT(L),NT) + MAX(SEDF(L,0,NS),0.)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      
      ! *** PARTICLE NONCOHESIVE EROSIONAL FLUX,  (SNDF-SNDFBL) > 0
      IF( ISTRAN(7) >= 1 )THEN
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO NX=1,NSND
            NS=NX+NSED
            DO L=LF,LL
              SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
              TOXFB(L,KBT(L),NT) = TOXFB(L,KBT(L),NT) + MAX(SNDFEFF,0.)*STFPOCB(L,KBT(L),NS)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN
          ! *** Kd APPROACH
          DO NX=1,NSND
            NS=NX+NSED
            DO L=LF,LL
              SNDFEFF = SNDF(L,0,NX) - SNDFBL(L,NX)
              TOXFB(L,KBT(L),NT) = TOXFB(L,KBT(L),NT) + MAX(SNDFEFF,0.)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ENDIF      
      ENDIF
    ENDDO

    ! *** FINALIZE EROSIONAL FLUX TOXFB (1/S)
    DO NT=1,NTOX
      DO L=LF,LL
        IF( HBED(L,KBT(L)) > 1E-7 )THEN
          ! ***                   M                           M
          TMPVAL = PORBED(L,KBT(L))*HBED(L,KBT(L)) + TOXPFTB(L,KBT(L),NT)    ! *** UNITS: M
          ! *** 1/S                 M/S              M
          TOXFB(L,KBT(L),NT) = TOXFB(L,KBT(L),NT)/TMPVAL
        ELSE
          TOXFB(L,KBT(L),NT) = 0.
        ENDIF
      ENDDO
    ENDDO

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! ******************************************************************************
  ! *** BANK EROSION
  IF( ISBKERO >= 1 )THEN
    !$OMP SINGLE
    DO NT=1,NTOX
      DO L=2,LA
        TOXFBEBKB(L,NT)=0.
        TOXFBECHB(L,NT)=0.
        TOXFBECHW(L,NT)=0.
      ENDDO
    ENDDO

    DO NT=1,NTOX
      ! PARTICLE cohesive flux
      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            KBANK=KBT(LBANK)
            TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)+SEDFBEBKB(NP,NS)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(NS,NT)
            TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)+SEDFBECHB(NP,NS)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(NS,NT)
            TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)+SEDFBECHW(NP,NS)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(NS,NT)
          ENDDO
        ENDDO
      ENDIF

      ! PARTICLE noncohesive flux
      IF( ISTRAN(7) >= 1 )THEN
        DO NX=1,NSND
          NS=NX+NSED
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            KBANK=KBT(LBANK)
            TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)+SNDFBEBKB(NP,NX)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(NS,NT)
            TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)+SNDFBECHB(NP,NX)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(NS,NT)
            TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)+SNDFBECHW(NP,NX)*STFPOCB(LBANK,KBANK,NS)*TOXPARB(NS,NT)
          ENDDO
        ENDDO
      ENDIF
    ENDDO

    DO NT=1,NTOX
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        KBANK=KBT(LBANK)
        IF( HBED(LBANK,KBANK) > 0.0 )THEN
          TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)/(PORBED(LBANK,KBANK)*HBED(LBANK,KBANK)+TOXPFTB(LBANK,KBANK,NT))
          TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)/(PORBED(LBANK,KBANK)*HBED(LBANK,KBANK)+TOXPFTB(LBANK,KBANK,NT))
          TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)/(PORBED(LBANK,KBANK)*HBED(LBANK,KBANK)+TOXPFTB(LBANK,KBANK,NT))
        ELSE
          TOXFBEBKB(NP,NT)=0.
          TOXFBECHB(NP,NT)=0.
          TOXFBECHW(NP,NT)=0.
        ENDIF
      ENDDO
    ENDDO
    !$OMP END SINGLE
  ENDIF 
  ! *** END BANK EROSION

  !$OMP DO PRIVATE(ND,LF,LL,NT,K,LP,L)
  DO ND=1,NDM
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    !**********************************************************************CC
    ! **  CALCULATE TOTAL PARTICULATE FRACTIONS IN WATER COLUMN AND BED
    ! **  NOTING THAT TO THIS POINT TOXPFTW AND TOXPFTB HAVE BEEN USED
    ! **  TO TEMPORARILY STORE THE SORBED PORTION
    DO NT=1,NTOX
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          ! ***                                                          SORBED TO DOC
          TOXPFTW(L,K,NT) = ( TOXPFTW(L,K,NT)/(1.+TOXPFTW(L,K,NT)) ) - TOXPFW(L,K,NFD,NT)
          TOXPFTW(L,K,NT) = MIN(TOXPFTW(L,K,NT),SORBMIN)
        ENDDO
      ENDDO
    ENDDO
    
    ! *** COMPUTE MASS FRACTION (TOXPFTB) OF SORBED TOXICS TO PARTICULATES (DIMENSIONLESS)
    DO NT=1,NTOX
      DO K=1,KB
        DO LP=LF,LL
          L=LWET(LP)
          IF( HBED(L,K) >= HBEDMIN0 .AND. PORBED(L,K) > 0. )THEN
            ! ***                                                                                 SORBED TO DOC
            TOXPFTB(L,K,NT) = ( TOXPFTB(L,K,NT)/( PORBED(L,K)*HBED(L,K) + TOXPFTB(L,K,NT) ) ) - TOXPFB(L,K,NFD,NT)
            TOXPFTB(L,K,NT) = MIN(TOXPFTB(L,K,NT),SORBMIN)
          ELSE
            TOXPFTB(L,K,NT)=0.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO       ! *** END OF DOMAIN LOOP
  !$OMP END DO

  !**********************************************************************C
  ! **  DETERMINE TOXIC FLUX FROM BED LOAD SORBED MATERIAL
  IF( ISBDLDBC > 0 .AND. .NOT. LSEDZLJ )THEN  
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,NT,NS,NX,LE,LW,LN,LS) &
    !$OMP    PRIVATE(TMPTOXC,TMPTOXE,TMPTOXW,TMPTOXN,TMPTOXS)
    DO ND=1,NDM
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO NT=1,NTOX
        DO LP=LF,LL
          L=LWET(LP)
          TOXFBL(L,NT)=0.0     ! *** UNITS:  MG/M2/S
        ENDDO
      ENDDO

      DO NT=1,NTOX
        DO NX=1,NSND
          NS=NX+NSED
          DO LP=LF,LL
            L=LWET(LP)
            LE=LEC(L)
            LW=LWC(L)
            LN=LNC(L)
            LS=LSC(L)
            TMPTOXC=0.0
            TMPTOXE=0.0
            TMPTOXW=0.0
            TMPTOXN=0.0
            TMPTOXS=0.0
            K=KBT(L)
            IF( SNDB1(L,K,NX) > 0.0 )THEN
              ! ***         DIMENSIONLESS    MG/M2        M          M2/G      (ORIGINAL CODE)
              !TMPTOXC = TOXPFB(L,K,NS,NT)*TOXB(L,K,NT)*HBED(L,K)/SNDB(L,K,NX)
              ! MG/G         DIMENSIONLESS     MG/M2        M2/G 
              TMPTOXC = TOXPFB(L,K,NS,NT)*TOXB(L,K,NT)/SNDB1(L,K,NX) 
            ENDIF
            K=KBT(LE)
            IF( SUB(LE) > 0.5 .AND. SNDB1(LE,K,NX) > 0.0 )THEN
              TMPTOXE = TOXPFB(LE,K,NS,NT)*TOXB(LE,K,NT)/SNDB1(LE,K,NX)
            ENDIF
            K=KBT(LW)
            IF( SUB(L) > 0.5 .AND. SNDB1(LW,K,NX) > 0.0 )THEN
              TMPTOXW = TOXPFB(LW,K,NS,NT)*TOXB(LW,K,NT)/SNDB1(LW,K,NX)
            ENDIF
            K=KBT(LN)
            IF( SVB(LN) > 0.5 .AND. SNDB1(LN,K,NX) > 0.0 )THEN
              TMPTOXN = TOXPFB(LN,K,NS,NT)*TOXB(LN,K,NT)/SNDB1(LN,K,NX)
            ENDIF
            K=KBT(LS)
            IF( SVB(L) > 0.5 .AND. SNDB1(LS,K,NX) > 0.0 )THEN
              TMPTOXS = TOXPFB(LS,K,NS,NT)*TOXB(LS,K,NT)/SNDB1(LS,K,NX)
            ENDIF

            ! ***             MG/M2/S       1/M2                        G/S           MG/G
            TOXFBL(L,NT) = TOXFBL(L,NT) + DXYIP(L)*( SUB(LE)*MAX(QSBDLDX(LE,NX),0.)*TMPTOXC + SUB(LE)*MIN(QSBDLDX(LE,NX),0.)*TMPTOXE  &
                                                   - SUB(L) *MAX(QSBDLDX(L,NX),0.) *TMPTOXW - SUB(L) *MIN(QSBDLDX(L,NX),0.) *TMPTOXC  &
                                                   + SVB(LN)*MAX(QSBDLDY(LN,NX),0.)*TMPTOXC + SVB(LN)*MIN(QSBDLDY(LN,NX),0.)*TMPTOXN  &
                                                   - SVB(L) *MAX(QSBDLDY(L,NX),0.) *TMPTOXS - SVB(L) *MIN(QSBDLDY(L,NX),0.) *TMPTOXC )
          ENDDO
        ENDDO
      ENDDO   ! *** NTOX
    ENDDO       ! *** END OF DOMAIN LOOP
    !$OMP END DO

    !$OMP SINGLE
    ! *** ADJUST FOR BEDLOAD TRANSPORT AT BOUNDARY CELLS
    IF( NSBDLDBC > 0 )THEN
      ! *** UPSTREAM AND DOWNSTREAM 
      DO NSB=1,NSBDLDBC
        LUTMP=LSBLBCU(NSB)
        LDTMP=LSBLBCD(NSB)
        DO NT=1,NTOX
          DO NX=1,NSND
            NS=NX+NSED
            TMPTOXC=0.0
            K=KBT(LUTMP)
            IF( SNDB1(LUTMP,K,NX) > 0.0 )THEN
              TMPTOXC=TOXPFB(LUTMP,K,NS,NT)*TOXB(LUTMP,K,NT)/SNDB1(LUTMP,K,NX)
            ENDIF  

            ! *** MG/M2/S            MG/M2/S          1/M2          G/S            MG/G          
            TOXFBL(LUTMP,NT) = TOXFBL(LUTMP,NT) - DXYIP(LUTMP)*QSBDLDOT(LUTMP,NX)*TMPTOXC
            IF( LDTMP /= 0 ) TOXFBL(LDTMP,NT) = TOXFBL(LDTMP,NT) + DXYIP(LDTMP)*QSBDLDOT(LUTMP,NX)*TMPTOXC
          ENDDO
        ENDDO
      ENDDO
    ENDIF   ! *** END OF NSBDLDBC > 0 BLOCK
    !$OMP END SINGLE
  ENDIF
  
  ! *** HANDLE BEDLOAD FOR SEDZLJ
  IF( NCALC_BL > 0 .AND. LSEDZLJ )THEN
    !$OMP DO PRIVATE(ND,LF,LL,NT,LP,L,NS,LE,LW,LN,LS) &
    !$OMP    PRIVATE(TMPTOXB,TMPTOXC,TMPTOXE,TMPTOXW,TMPTOXN,TMPTOXS,TMPVAL,TOXFLUX,CBLTOXTMP)
    DO ND=1,NDM
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO NT=1,NTOX
        DO LP=LF,LL
          L=LWET(LP)
          TOXFBL(L,NT)=0.0     ! *** UNITS:  MG/M2/S  (+ is Erosional, - is Depositional)
        ENDDO
      ENDDO

      ! *** Compute bedload toxics concentration
      DO NT=1,NTOX
        DO NS=NNONCO,NSCM
          DO LP=LF,LL
            L=LWET(LP)
            LE=LEC(L)
            LW=LWC(L)
            LN=LNC(L)
            LS=LSC(L)
            IF( SEDB1(L,KBT(L),NS) > 1.E-6 )THEN
              ! MG/G          DIMENSIONLESS          MG/M2            M2/G
              TMPTOXB = TOXPFB(L,KBT(L),NS,NT)*TOXB1(L,KBT(L),NT)/SEDB1(L,KBT(L),NS) 
            ELSE
              TMPTOXB=0.0
            ENDIF
            ! MG/G                 DIMENSIONLESS           MG/G
            TMPTOXC =         TOXPFB(L,KBT(L),NS,NT)  *CBLTXCON(L,NS,NT)
            TMPTOXE = SUB(LE)*TOXPFB(LE,KBT(LE),NS,NT)*CBLTXCON(LE,NS,NT)
            TMPTOXW = SUB(L)* TOXPFB(LW,KBT(LW),NS,NT)*CBLTXCON(LW,NS,NT)
            TMPTOXN = SVB(LN)*TOXPFB(LN,KBT(LN),NS,NT)*CBLTXCON(LN,NS,NT)
            TMPTOXS = SVB(L)* TOXPFB(LS,KBT(LS),NS,NT)*CBLTXCON(LS,NS,NT)

            ! *** NET TOXIC FLUX IN/OUT OF THE BED (MG/M2)
            TOXFLUX = EBL(L,NS)*TMPTOXB - DBL(L,NS)*TMPTOXC

            ! *** CONCENTRATION OF TOXIC IN BEDLOAD LAYER (MG/M2)
            TMPVAL = DTSED*DXYIP(L)
            CBLTOXTMP = CBLTOX(L,NT) + TMPVAL*(  MAX(QSBDLDX(L,NS),0.)*TMPTOXW - MIN(QSBDLDX(LE,NS),0.)*TMPTOXE   &   ! *** EAST/WEST: IN 
                                               + MIN(QSBDLDX(L,NS),0.)*TMPTOXC - MAX(QSBDLDX(LE,NS),0.)*TMPTOXC   &   ! *** EAST/WEST: OUT 
                                               + MAX(QSBDLDY(L,NS),0.)*TMPTOXS - MIN(QSBDLDY(LN,NS),0.)*TMPTOXN   &   ! *** NORTH/SOUTH: IN
                                               + MIN(QSBDLDY(L,NS),0.)*TMPTOXC - MAX(QSBDLDY(LN,NS),0.)*TMPTOXC ) &   ! *** NORTH/SOUTH: OUT
                                               + TOXFLUX
            
            IF( CBLTOXTMP >= 0. )THEN
              CBLTOX(L,NT) = CBLTOXTMP
            ELSE
              ! *** LIMIT FLUX TO AVAILABLE MASS
              CBLTOX(L,NT) = 0.0
              TOXFBL(L,NT) = TOXFBL(L,NT) - CBLTOXTMP*DELTI
              EXIT
            ENDIF
            
            TOXFBL(L,NT) = TOXFBL(L,NT) + TOXFLUX*DELTI                                                                   ! *** TOXFBL FLUX RATE   (MG/M2/S)
          ENDDO
        ENDDO
      ENDDO   ! *** NTOX

      ! *** CHECK FOR DEPLETED BEDLOAD MASS
      DO LP=LF,LL
        L=LWET(LP)
        IF( SUM(CBL(L,NNONCO:NSCM)) <= 1.E-8 )THEN
          ! *** MOVE ANY REMAINING MASS INTO BED
          TOXB(L,KBT(L),1:NTOX) = TOXB(L,KBT(L),1:NTOX) + CBLTOX(L,1:NTOX)
          CBLTOX(L,1:NTOX) = 0.0
        ENDIF
      ENDDO
      
    ENDDO       ! *** END OF DOMAIN LOOP
    !$OMP END DO

    ! *** ADJUST FOR BEDLOAD TRANSPORT AT BOUNDARY CELLS
    IF( NSBDLDBC > 0 )THEN
      !$OMP SINGLE
      ! *** UPSTREAM AND DOWNSTREAM 
      DO NSB=1,NSBDLDBC
        LUTMP=LSBLBCU(NSB)
        LDTMP=LSBLBCD(NSB)
        DO NT=1,NTOX
          TOXFBL(LUTMP,NT) = 0.0
          DO NS=NNONCO,NSCM
            TMPTOXC = TOXPFB(LUTMP,KBT(LUTMP),NS,NT)*CBLTXCON(LUTMP,NS,NT)

            ! *** MG/M2/S            MG/M2/S          1/M2          G/S            MG/G
            !TOXFBL(LUTMP,NT) = TOXFBL(LUTMP,NT) - DXYIP(LUTMP)*QSBDLDOT(LUTMP,NS)*TMPTOXC
            
            ! *** MG/M2              MG/M2/S          1/M2          G/S            MG/G    S
            !CBLTOX(LUTMP,NT) = CBLTOX(LUTMP,NT) - DXYIP(LUTMP)*QSBDLDOT(LUTMP,NS)*TMPTOXC*DELT
            IF( LDTMP /= 0 )THEN
              ! *** MG/M2/S            MG/M2/S          1/M2          G/S            MG/G
              TOXFBL(LDTMP,NT) = TOXFBL(LDTMP,NT) - DXYIP(LDTMP)*QSBDLDOT(LUTMP,NS)*TMPTOXC
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      !$OMP END SINGLE
    ENDIF   ! *** END OF NSBDLDBC > 0 BLOCK
  ENDIF     ! *** END OF ISBDLDBC > 0 BLOCK

  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,NT,NS,NX,TMPVAL)
  DO ND=1,NDM
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    ! ***
    IF( (NCALC_BL > 0 .AND. LSEDZLJ) .OR. (ISBDLDBC > 0 .AND. .NOT. LSEDZLJ) ) THEN
      DO NT=1,NTOX
        TOXFBLT(NT) = 0.
        DO LP=LF,LL
          L=LWET(LP)
          TOXFBLT(NT) = TOXFBLT(NT) + DXYP(L)*TOXFBL(L,NT)
        ENDDO
      ENDDO
    ENDIF
    
    !**********************************************************************C
    ! **  ADJUST TOXIC FLUXES ACROSS WATER COLUMN - BED INTERFACE TO
    ! **  INCLUDE WATER ENTRAINMENT AND EXPULSION ASSOCIATED WITH
    ! **  DEPOSITION AND RESUSPENSION
    
    ! *** DEPSITIONAL FLUX ADJUSTMENT DUE TO ENTRAPMENT (M/S)
    DO NT=1,NTOX
      DO LP=LF,LL
        L=LWET(LP)
        TMPVAL = ( MIN(QSBDTOP(L),0.0) + MIN(QWBDTOP(L),0.0) )    ! *** TOTAL DEPOSITIONAL VOLUMETRIC RATE (M/S)
        TOXF(L,0,NT) = TOXF(L,0,NT) + TMPVAL*(1.-TOXPFTW(L,KSZ(L),NT))
      ENDDO
    ENDDO
  
    ! *** EROSONAL FLUX ADJUSTMENT DUE TO POREWATER EXPULSION (1/S)
    DO NT=1,NTOX
      DO LP=LF,LL
        L=LWET(LP)
        K=KBT(L)
        IF( HBED(L,K) > 1.E-6 )THEN
          TMPVAL = ( MAX(QSBDTOP(L),0.0) + MAX(QWBDTOP(L),0.0) )/HBED(L,K)      ! *** TOTAL EROSIONAL VOLUMETRIC RATE (M/S)/BED THICKNESS (M)
        ELSE
          TMPVAL = 0.
        ENDIF       
        TOXFB(L,K,NT) = TOXFB(L,K,NT) + TMPVAL*(1.-TOXPFTB(L,K,NT))
      ENDDO
    ENDDO
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! *** START BANK EROSION
  IF( ISBKERO >= 1 )THEN
    !$OMP SINGLE
    DO NT=1,NTOX
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        KBANK=KBT(LBANK)
        TMPVAL=(1.-TOXPFTB(LBANK,KBANK,NT))/HBED(LBANK,KBANK)
        TOXFBEBKB(NP,NT)=TOXFBEBKB(NP,NT)+TMPVAL*(QSBDTOPBEBKB(NP)+QWBDTOPBEBKB(NP))
        TOXFBECHB(NP,NT)=TOXFBECHB(NP,NT)+TMPVAL*(QSBDTOPBECHB(NP)+QWBDTOPBECHB(NP))
        TOXFBECHW(NP,NT)=TOXFBECHW(NP,NT)+TMPVAL*(QSBDTOPBECHW(NP)+QWBDTOPBECHW(NP))
      ENDDO
    ENDDO
    !$OMP END SINGLE
  ENDIF
  !---END BANK EROSION-------------------------

  ! ***********************************************************************************
  ! ***  UPDATE WATER COLUMN BOTTOM LAYER AND TOP SEDIMENT LAYER FOR TOXIC CONTAMINANT
  ! ***  PMC - NEW APPROACH FOR BETTER MASS BALANCE 2017-07
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,NT,NS,NX,LE,LW,LN,LS)  &
  !$OMP    PRIVATE(AA11,BB11,BB22,TMPVAL,CLEFT,CRIGHT,WVEL)
  DO ND=1,NDM
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    ! ----------------------------------------------------------------
    ! *** KC=1 (SINGLE LAYER IN VERTICAL)
    IF( KC == 1 )THEN
      DO NT=1,NTOX
        ! *** WATER COLUMN/SEDIMENT BED INTERFACE
        DO LP=LF,LL
          L=LWET(LP)

          ! *** MASS BALANCE AROUND BOTTOM WATER COLUMN LAYER
          BB11 = TOXFB(L,KBT(L),NT)*TOXB(L,KBT(L),NT)     ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          BB22 = TOXF(L,0,NT)*TOX(L,KSZ(L),NT)            ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S
          
          ! *** LIMIT DEPOSITING MASS TO AVAILABLE MASS
          TMPVAL = TOX(L,KSZ(L),NT) + DTSED*( BB11 + BB22 )*HPKI(L,KSZ(L))
          IF( TMPVAL < 0. )THEN
            BB11 = -TOX(L,KSZ(L),NT)*HPK(L,KSZ(L))*DELTI
            BB22 = 0.0
            TOX(L,KSZ(L),NT) = 0.0
          ELSE
            TOX(L,KSZ(L),NT) = TMPVAL
          ENDIF
          
          ! *** MASS BALANCE AROUND TOP SEDIMENT LAYER
          !TOXB1(L,KBT(L),NT) = S3TL*TOXB(L,KBT(L),NT) + S2TL*TOXB1(L,KBT(L),NT)   ! NOT USED
          TOXB(L,KBT(L),NT) = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          TOXB(L,KBT(L),NT) = MAX(TOXB(L,KBT(L),NT),0.0)
        ENDDO
      ENDDO   ! *** NTOX
    ENDIF     ! *** KC=1

    ! ----------------------------------------------------------------
    ! *** KC=2 (TWO LAYERS IN VERTICAL)
    IF( KC == 2 )THEN
      DO NT=1,NTOX
        ! *** WATER COLUMN
        K = KC
        DO L=LF,LL
          WVEL=DELTI*HPK(L,K)
          CLEFT=WVEL-TOXF(L,K-1,NT)
          CRIGHT=WVEL*TOX(L,K,NT)
          TOX(L,K,NT)=CRIGHT/CLEFT
        ENDDO
       
        ! *** WATER COLUMN/SEDIMENT BED INTERFACE
        DO LP=LF,LL
          L=LWET(LP)

          ! *** MASS BALANCE AROUND BOTTOM WATER COLUMN LAYER
          AA11 = TOXF(L,KSZ(L),NT)*TOX(L,KSZ(L)+1,NT)     ! *** SETTLING FROM TOP OF BOTTOM LAYER,        MG/M2/S
          BB11 = TOXFB(L,KBT(L),NT)*TOXB(L,KBT(L),NT)     ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          BB22 = TOXF(L,0,NT)*TOX(L,KSZ(L),NT)            ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S
          
          ! *** LIMIT DEPOSITING MASS TO AVAILABLE MASS
          TMPVAL = TOX(L,KSZ(L),NT) + DTSED*( BB11 + BB22 - AA11 )*HPKI(L,KSZ(L))
          IF( TMPVAL < 0. )THEN
            BB11 = -TOX(L,KSZ(L),NT)*HPK(L,KSZ(L))*DELTI
            BB22 = 0.0
            TOX(L,KSZ(L),NT) = 0.0
          ELSE
            TOX(L,KSZ(L),NT) = TMPVAL
          ENDIF
          
          ! *** MASS BALANCE AROUND TOP SEDIMENT LAYER
          !TOXB1(L,KBT(L),NT) = S3TL*TOXB(L,KBT(L),NT) + S2TL*TOXB1(L,KBT(L),NT)   ! NOT USED
          TOXB(L,KBT(L),NT) = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          TOXB(L,KBT(L),NT) = MAX(TOXB(L,KBT(L),NT),0.0)
        ENDDO
      ENDDO   ! *** NTOX
    ENDIF     ! *** KC = 2

    ! ----------------------------------------------------------------
    ! *** KC >= 3 (THREE OR MORE LAYERS IN VERTICAL)
    IF( KC >= 3 )THEN
      DO NT=1,NTOX
        ! *** WATER COLUMN
        K=KC
        DO LP=LF,LL
          L=LWET(LP)
          WVEL = DELTI*HPK(L,K)
          CLEFT = WVEL-TOXF(L,K-1,NT)
          CRIGHT = WVEL*TOX(L,K,NT)
          TOX(L,K,NT) = CRIGHT/CLEFT
        ENDDO
       
        DO K=KS,2,-1 
          DO LP=1,LLWET(K-1,ND)
            L = LKWET(LP,K,ND)  
            WVEL = DELTI*HPK(L,K)
            CLEFT = WVEL - TOXF(L,K-1,NT)
            CRIGHT = WVEL*TOX(L,K,NT) - TOXF(L,K,NT)*TOX(L,K+1,NT)
            TOX(L,K,NT) = CRIGHT/CLEFT
          ENDDO
        ENDDO

        ! *** WATER COLUMN/SEDIMENT BED INTERFACE
        DO LP=LF,LL
          L=LWET(LP)

          ! *** MASS BALANCE AROUND BOTTOM WATER COLUMN LAYER
          AA11 = TOXF(L,KSZ(L),NT)*TOX(L,KSZ(L)+1,NT)     ! *** SETTLING FROM TOP OF BOTTOM LAYER,        MG/M2/S
          BB11 = TOXFB(L,KBT(L),NT)*TOXB(L,KBT(L),NT)     ! *** EROSIONAL FLUX TO THE WATER COLUMN,       MG/M2/S
          BB22 = TOXF(L,0,NT)*TOX(L,KSZ(L),NT)            ! *** DEPOSITIONAL FLUX FROM THE WATER COLUMN,  MG/M2/S
          
          ! *** LIMIT DEPOSITING MASS TO AVAILABLE MASS
          TMPVAL = TOX(L,KSZ(L),NT) + DTSED*( BB11 + BB22 - AA11 )*HPKI(L,KSZ(L))
          IF( TMPVAL < 0. )THEN
            BB11 = -TOX(L,KSZ(L),NT)*HPK(L,KSZ(L))*DELTI
            BB22 = 0.0
            TOX(L,KSZ(L),NT) = 0.0
          ELSE
            TOX(L,KSZ(L),NT) = TMPVAL
          ENDIF
          
          ! *** MASS BALANCE AROUND TOP SEDIMENT LAYER
          !TOXB1(L,KBT(L),NT) = S3TL*TOXB(L,KBT(L),NT) + S2TL*TOXB1(L,KBT(L),NT)   ! NOT USED
          TOXB(L,KBT(L),NT) = TOXB(L,KBT(L),NT) - DTSED*( BB22 + BB11 + TOXFBL(L,NT) )
          TOXB(L,KBT(L),NT) = MAX(TOXB(L,KBT(L),NT),0.0)
        ENDDO
      ENDDO   ! *** NTOX
    ENDIF     ! *** KC >= 3
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  ! *** START BANK EROSION
  IF( ISBKERO >= 1 )THEN
    !$OMP SINGLE
    DO NP=1,NBEPAIR
      LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
      LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
      KBANK=KBT(LBANK)
      KCHAN=KBT(LCHAN)
      TOXFBEBKB(NP,NT) = TOXFBEBKB(NP,NT)*TOXB(LBANK,KBANK,NT)
      TOXFBECHB(NP,NT) = TOXFBECHB(NP,NT)*TOXB(LBANK,KBANK,NT)
      TOXFBECHW(NP,NT) = TOXFBECHW(NP,NT)*TOXB(LBANK,KBANK,NT)
      TOXB(LBANK,KBANK,NT) = TOXB(LBANK,KBANK,NT)-DTSED*TOXFBEBKB(NP,NT)
      TOXB(LCHAN,KCHAN,NT) = TOXB(LCHAN,KCHAN,NT)-DTSED*TOXFBECHB(NP,NT)
      TOX(LCHAN,KSZ(LCHAN),NT) = TOX(LCHAN,KSZ(LCHAN),NT) + DTSED*HPKI(L,KSZ(LCHAN))*TOXFBECHW(NP,NT)
    ENDDO
    !$OMP END SINGLE
  ENDIF

  IF( .NOT. LSEDZLJ )THEN
    !$OMP DO PRIVATE(ND,LF,LL,L,K,NT,NS,NX,LE,LW,LN,LS,KTOPTP,KTOPM1)  &
    !$OMP    PRIVATE(FTPOS,FTNEG)
    DO ND=1,NDM
      LF=2+(ND-1)*LDM
      LL=MIN(LF+LDM-1,LA)

      !**********************************************************************C
      !
      ! **  ADD PARENT TO ACTIVE LAYER TOXIC TRANSPORT
      !
      IF( ISNDAL >= 2 )THEN
        DO NT=1,NTOX
          DO L=LF,LL
            TOXFPA(L)=0.0
          ENDDO
          IF( ISTRAN(6) > 0 )THEN
            DO NS=1,NSED
              DO L=LF,LL
                KTOPTP=KBT(L)
                KTOPM1=KBT(L)-1
                IF( KTOPM1 < 1 ) CYCLE
                FTPOS=0.0
                FTNEG=0.0
                IF( SEDFPA(L,NS) > 0.0 .AND. SEDB(L,KTOPM1,NS) > 0.0 ) FTPOS = SEDFPA(L,NS)*TOXPFB(L,KTOPM1,NS,NT)/SEDB(L,KTOPM1,NS)
                IF( SEDFPA(L,NS) < 0.0 .AND. SEDB(L,KTOPTP,NS) > 0.0 ) FTNEG = SEDFPA(L,NS)*TOXPFB(L,KTOPTP,NS,NT)/SEDB(L,KTOPTP,NS)
                TOXFPA(L) = TOXFPA(L) + FTPOS*TOXB(L,KTOPM1,NT) + FTNEG*TOXB(L,KTOPTP,NT)
              ENDDO
            ENDDO
          ENDIF
          
          IF( ISTRAN(7) > 0 )THEN
            DO NX=1,NSND
              NS=NX+NSED
              DO L=LF,LL
                KTOPTP=KBT(L)
                KTOPM1=KBT(L)-1
                IF( KTOPM1 < 1 ) CYCLE
                FTPOS=0.0
                FTNEG=0.0
                IF( SNDFPA(L,NX) > 0.0 .AND. SNDB(L,KTOPM1,NX) > 0.0 ) FTPOS = SNDFPA(L,NX)*TOXPFB(L,KTOPM1,NS,NT)/SNDB(L,KTOPM1,NX)
                IF( SNDFPA(L,NX) < 0.0 .AND. SNDB(L,KTOPTP,NX) > 0.0 ) FTNEG = SNDFPA(L,NX)*TOXPFB(L,KTOPTP,NS,NT)/SNDB(L,KTOPTP,NX)
                TOXFPA(L) = TOXFPA(L) + FTPOS*TOXB(L,KTOPM1,NT) + FTNEG*TOXB(L,KTOPTP,NT)
              ENDDO
            ENDDO
          ENDIF
          
          DO L=LF,LL
            KTOPTP=KBT(L)
            KTOPM1=KBT(L)-1
            IF( KTOPM1 < 1 ) CYCLE
            FTPOS=0.0
            FTNEG=0.0
            IF( QWATPA(L) > 0.0 ) FTPOS=QSSDPA(L)*(1.-TOXPFTB(L,KTOPM1,NT))
            IF( QWATPA(L) < 0.0 ) FTNEG=QSSDPA(L)*(1.-TOXPFTB(L,KTOPTP,NT))
            TOXFPA(L) = TOXFPA(L) + FTPOS*TOXB(L,KTOPM1,NT) + FTNEG*TOXB(L,KTOPTP,NT)
          ENDDO
          DO L=LF,LL
            KTOPTP=KBT(L)
            KTOPM1=KBT(L)-1
            IF( KTOPM1 < 1 ) CYCLE
            TOXB(L,KTOPTP,NT) = TOXB(L,KTOPTP,NT)+DTSED*TOXFPA(L)
            TOXB(L,KTOPM1,NT) = TOXB(L,KTOPM1,NT)-DTSED*TOXFPA(L)
          ENDDO
        ENDDO
      ENDIF   ! *** ISNDAL >= 2
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO
  ENDIF

  !$OMP END PARALLEL
  
  RETURN

END

