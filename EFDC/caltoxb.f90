SUBROUTINE CALTOXB

  ! *** SUBROUTINE CALTOXB CALCULATES CONTAMINANT TRANSPORT WITHIN THE SEDIMENT BED
  ! *** AND THE FLUX OF CONTAMINANTS AT THE SEDIMENT/WATER INTERFACE
  ! *** 
  ! *** USED FOR BOTH STANDARD EFDC SEDIMENT TRANSPORT AND SEDZLJ 
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2017-05-15        PAUL M. CRAIG    Updated to work for Original and SEDZLJ
  ! 2014-08-26        PAUL M. CRAIG    Updated to allow for any mix of ISTOC
  ! 2012-12-05        PAUL M. CRAIG    Updated OMP
  ! 2012-10-02        Dang H Chung     Added OMP
  !**********************************************************************!

  USE GLOBAL

  IMPLICIT NONE

  INTEGER :: ND,LF,LL,L,NT,K,KBTM1
  INTEGER :: NFD,KM,KK,KBTP1,KBOT,KINC,KBOT2
  REAL   :: CELLMASS,DEPINBED,DIFBWFAC,BETTMP
  REAL   :: ERRT,ERRB,ERRW,HBEDMIN0,SORBMIN,TOXTIMEI
  REAL,SAVE :: TOXTIME
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: DIFTOXBW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: PARTDIF
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DERRB
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TOXBBALN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TOXBBALO

  IF(  .NOT. ALLOCATED(DIFTOXBW) )THEN
    ALLOCATE(DIFTOXBW(LCM))
    ALLOCATE(PARTDIF(LCM,KBM))
    ALLOCATE(DERRB(KBM,NDM))  
    ALLOCATE(TOXBBALN(LCM))  
    ALLOCATE(TOXBBALO(LCM))  
    DIFTOXBW=0.0
    PARTDIF=0.0
    DERRB=0.0
    TOXBBALN=0.0
    TOXBBALO=0.0
    TOXTIME = 0.0
    TOXSTEPB = TOXSTEPB - DTSED/10.
  ENDIF
  NFD=NSED+NSND+1

  TOXTIME = TOXTIME + DTSED
  IF( TOXTIME < TOXSTEPB ) RETURN
  TOXTIMEI = 1./TOXTIME
  
  ! *** SET LAYER ORIENTATION
  IF( LSEDZLJ )THEN
    KINC  = 1
    KBOT  = KB
    KBOT2 = KB-1
    HBEDMIN0 = 1E-4  ! *** ACCOUNT FOR ACTIVE LAYERS
  ELSE
    KINC  = -1
    KBOT  = 1
    KBOT2 = 2
    HBEDMIN0 = 1E-4
  ENDIF
  
  ND = PRECISION(SORBMIN)
  SORBMIN = 0.99  ! 1. - 1./10**ND
  
  !**********************************************************************C
  !
  ! **  UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
  !
  !$OMP PARALLEL DEFAULT(SHARED) 
  !$OMP DO PRIVATE(ND,LF,LL,L,K,NT)
  DO ND=1,NDM
    LF=2+(ND-1)*LDM
    LL=MIN(LF+LDM-1,LA)

    !**********************************************************************C
    !
    ! **  CALCULATE TOXIC CONTAMINANT PARTICULATE FRACTIONS IN SEDIMENT BED
    !
    ! **  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED 
    ! **  TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED
    CALL CALTOXB_FRACTIONS(LF,LL)
    
    ! *** COMPUTE MASS FRACTION (TOXPFTB) OF SORBED TOXICS TO PARTICULATES (DIMENSIONLESS)
    DO NT=1,NTOX
      DO K=1,KB
        DO L=LF,LL
          IF( HBED(L,K) >= HBEDMIN0 .AND. PORBED(L,K) > 0. )THEN
            ! ***                                                                             SORBED TO DOC
            TOXPFTB(L,K,NT) = ( TOXPFTB(L,K,NT)/(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT)) ) - TOXPFB(L,K,NFD,NT)
            TOXPFTB(L,K,NT) = MIN(TOXPFTB(L,K,NT),SORBMIN)
          ELSE
            TOXPFTB(L,K,NT) = 0.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  !$OMP DO PRIVATE(ND,LF,LL,L,NT,K,KM,KK,DEPINBED,KBTP1,KBTM1,DIFBWFAC,BETTMP,ERRT,CELLMASS,ERRB,ERRW)
  DO ND=1,NDM
    LF=2+(ND-1)*LDM
    LL=MIN(LF+LDM-1,LA)

    !**********************************************************************C
    !
    ! **  ADVECT AND DIFFUSE TOXICS IN BED AND INTO BOTTOM WATER
    ! **  COLUMN LAYER

    ! **  ADD PARTICLE MIXING AND SCALE PARTICLE DIFFUSION FOR SOLUTION
    DO NT=1,NTOX
      DO L=LF,LL
        ! *** BYPASS FOR HARDBOTTOM
        !IF( BEDMAP(L) == 2 )CYCLE
          
        DEPINBED=0.
        PARTDIF(L,KBT(L)) = 0.0
        IF( KBT(L) /= KBOT2 )THEN
          DO K=KBT(L),KBOT2,KINC
            KM = K+KINC
            DEPINBED = DEPINBED+HBED(L,K)
            PARTDIF(L,KM) = 0.0
            IF( DEPINBED < DPDIFTOX(NT) )THEN
              PARTDIF(L,KM) = 2.*PDIFTOX(NT)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
            ENDIF
          ENDDO
        ELSE
          K  = KBOT2
          KM = K+KINC
          DEPINBED = DEPINBED+HBED(L,K)
          PARTDIF(L,KM) = 0.0
          IF( DEPINBED < DPDIFTOX(NT) )THEN
            PARTDIF(L,KM) = 2.*PDIFTOX(NT)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT))
          ENDIF
        ENDIF
      ENDDO

      !DO L=LF,LL
      !  DIFTOXBW(L) = 0.0
      !ENDDO
    
      DO L=LF,LL
        DIFTOXBW(L) = DIFTOXS(NT)
      ENDDO
    
      DO L=LF,LL
        IF( HBED(L,KBT(L)) < HBEDMIN0 ) CYCLE
        
        DIFBWFAC=2./HBED(L,KBT(L))
        IF( ISDIFBW(NT) == 1 ) DIFBWFAC=1.0
        TOXBBALO(L) = 0.
        KBTP1 = KBT(L)-KINC   ! *** EMPTY SEDIMENT LAYER TO ACCUMULATE MASS FLUX TO WATER COLUMN
        KBTM1 = KBT(L)+KINC
        ALOW(L,KBOT) = 0. 
        CUPP(L,KBTP1) = 0.

        DO K=KBOT,KBTM1,-KINC
          CUPP(L,K) = MIN(QWTRBED(L,K),0.) - (DIFTOX(NT)+PARTDIF(L,K))*(PORBED(L,K)+PORBED(L,K-KINC))/(HBED(L,K)+HBED(L,K-KINC))
        ENDDO
        CUPP(L,KBT(L)) = MIN(QWTRBED(L,KBT(L)),0.) - DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))
        
        DO K=KBOT2,KBT(L),-KINC
          ALOW(L,K) = -MAX(QWTRBED(L,K+KINC),0.) - (DIFTOX(NT)+PARTDIF(L,K+KINC))*(PORBED(L,K+KINC)+PORBED(L,K))/(HBED(L,K+KINC)+HBED(L,K))
        ENDDO
        ALOW(L,KBTP1) = -MAX(QWTRBED(L,KBT(L)),0.) - DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))
        
        DO K=KBOT,KBT(L),-KINC
          BMNN(L,K) = TOXTIMEI*HBED(L,K)*PORBED(L,K)/(1.-TOXPFTB(L,K,NT))
        ENDDO
        BMNN(L,KBTP1) = TOXTIMEI*HPK(L,KSZ(L))/(1.-TOXPFTW(L,KSZ(L),NT))  ! *** BOTTOM LAYER OF WATER COLUMN
        
        ! *** BOTTOM LAYER
        BMNN(L,KBOT)  = BMNN(L,KBOT) - MIN(QWTRBED(L,0),0.)
        BMNN(L,KBOT)  = BMNN(L,KBOT) + MAX(QWTRBED(L,KBOT),0.) + (DIFTOX(NT)+PARTDIF(L,KBOT))*(PORBED(L,KBOT2)+PORBED(L,KBOT))/(HBED(L,KBOT2)+HBED(L,KBOT))
        
        ! *** MIDDLE LAYERS
        DO K=KBOT2,KBTM1,-KINC
          BMNN(L,K) = BMNN(L,K)+MAX(QWTRBED(L,K),0.)     +(DIFTOX(NT)+PARTDIF(L,K))     *(PORBED(L,K-KINC)+PORBED(L,K))/(HBED(L,K-KINC)+HBED(L,K)) &
                               -MIN(QWTRBED(L,K+KINC),0.)+(DIFTOX(NT)+PARTDIF(L,K+KINC))*(PORBED(L,K+KINC)+PORBED(L,K))/(HBED(L,K+KINC)+HBED(L,K))
        ENDDO
        
        ! *** TOP LAYERS
        K=KBT(L)
        IF( K == KBOT )THEN
          BMNN(L,K) = BMNN(L,K)+MAX(QWTRBED(L,K),0.) +DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))
        ELSE
          BMNN(L,K) = BMNN(L,K)+MAX(QWTRBED(L,K),0.)      + DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L)) &
                               -MIN(QWTRBED(L,K+KINC),0.) + (DIFTOX(NT)+PARTDIF(L,K+KINC))*(PORBED(L,K+KINC)+PORBED(L,K))/(HBED(L,K+KINC)+HBED(L,K))  
        ENDIF
      
        BMNN(L,KBTP1) = BMNN(L,KBTP1) - MIN(QWTRBED(L,KBT(L)),0.) + DIFBWFAC*DIFTOXBW(L)*PORBED(L,KBT(L))

        ! *** COMPUTE RRHS KG/S
        DO K=KBOT,KBT(L),-KINC
          RRHS(L,K) =  TOXTIMEI*TOXB(L,K,NT)
          TOXBBALO(L) = TOXBBALO(L)+TOXB(L,K,NT)               ! *** TOTAL MASS OF TOXICS IN BED BEFORE DIFFUSION STEPS (MG/M2)
        ENDDO
        TOXWBALO(L) = HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)           ! *** TOTAL MASS OF TOXICS IN BOTTOM LAYER BEFORE DIFFUSION STEPS (MG/M2)
        
        ! *** ADD FLUX OF GROUNDWATER INTO THE BOTTOM SEDIMENT LAYER
        RRHS(L,KBOT) = RRHS(L,KBOT) + MAX(QWTRBED(L,0),0.)*CONGW(L,NT+4)
        RRHS(L,KBTP1) = TOXTIMEI*TOXWBALO(L)                   ! *** RRHS - MASS FLUX RATE (MG/M2/S) 
      ENDDO

      ! **  TRI-DIAGONAL SOLVER
      DO L=LF,LL
        IF( HBED(L,KBT(L)) < HBEDMIN0 .OR. ABS(BMNN(L,KBOT)) < 1.0E-10 ) CYCLE
        KBTP1 = KBT(L)-KINC
        BETTMP = BMNN(L,KBOT)
        TOXTMP(L,KBOT) = RRHS(L,KBOT)/BETTMP
        DO KK=KBOT2,KBTP1,-KINC
          GAMTMP(L,KK) = CUPP(L,KK+KINC)/BETTMP
          BETTMP = BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)
          TOXTMP(L,KK) = (RRHS(L,KK) - ALOW(L,KK)*TOXTMP(L,KK+KINC))/BETTMP
        ENDDO
        DO KK=KBT(L),KBOT,KINC
          TOXTMP(L,KK) = TOXTMP(L,KK) - GAMTMP(L,KK-KINC)*TOXTMP(L,KK-KINC)
        ENDDO
      ENDDO

      ! **  CONVERT SCALED SOLUTION VARIABLES AND CALCULATE FINAL MASS
      DO L=LF,LL
        IF( HBED(L,KBT(L)) < HBEDMIN0 .OR. ABS(BMNN(L,KBOT)) < 1.0E-10 ) CYCLE

        TOXBBALN(L) = 0.0
        DO K=KBOT,KBT(L),-KINC
          TOXB(L,K,NT) = HBED(L,K)*PORBED(L,K)*TOXTMP(L,K)/(1.0-TOXPFTB(L,K,NT))
          TOXBBALN(L) = TOXBBALN(L) + TOXB(L,K,NT)                               ! *** TOTAL MASS OF TOXICS IN BED AFTER DIFFUSION STEPS (MG/M2)
        ENDDO
        
        KBTP1 = KBT(L)-KINC
        TOX(L,KSZ(L),NT) = TOXTMP(L,KBTP1)/(1.-TOXPFTW(L,KSZ(L),NT))
        TOXWBALN(L) = HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)                             ! *** TOTAL MASS OF TOXICS IN BOTTOM LAYER AFTER DIFFUSION STEPS (MG/M2)
      ENDDO
  
      ! **  CORRECT MASS ERROR AND DETERMINE NET FLUX FROM BED TO WATER COLUMN
      DO L=LF,LL
        IF( HBED(L,KBT(L)) < HBEDMIN0 .OR. ABS(BMNN(L,KBOT)) < 1.0E-10) CYCLE
        CELLMASS = TOXBBALN(L)+TOXWBALN(L)
        ERRT     = CELLMASS - (TOXBBALO(L)+TOXWBALO(L))
        
        ! *** HANDLE ZERO SEDIMENT AND/OR WATER CONCENTRATIONS
        IF( TOXBBALN(L) > 1.E-12 .AND. ERRT /= 0. )THEN
          ERRB = ERRT*TOXBBALN(L)/CELLMASS
          ERRW = ERRT-ERRB
          TOX(L,KSZ(L),NT) = TOX(L,KSZ(L),NT) - ERRW/HPK(L,KSZ(L))
          DO K=KBOT,KBT(L),-KINC
            DERRB(K,ND) = TOXB(L,K,NT)/TOXBBALN(L)
          ENDDO
          DO K=KBOT,KBT(L),-KINC
            TOXB(L,K,NT) = TOXB(L,K,NT) - DERRB(K,ND)*ERRB
          ENDDO
          TADFLUX(L,NT) = (HPK(L,KSZ(L))*TOX(L,KSZ(L),NT)-TOXWBALO(L))/TOXTIME
        ELSE
          TADFLUX(L,NT) = 0.
        ENDIF
      ENDDO
    ENDDO
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP END PARALLEL

  TOXTIME = 0.0
  
  RETURN

END

SUBROUTINE CALTOXB_FRACTIONS(LF,LL)

  !**********************************************************************C
  !
  ! ***  UPDATE TOTAL PARTICULATE FRACTION OF EACH TOXIC IN THE BED
  !

  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2014-08-26        PAUL M. CRAIG    Updated to allow for any mix of ISTOC

  ! ***  TOXPFB(L,NS,NT) = PARTICULATE FRACTION IN SEDIMENT BED (DIMENSIONLESS, INTERMEDIATE UNITS METERS) 
  ! ***  TOXPARB(NS,NT)  = PARTITION COEFFICIENT IN SEDIMENT BED (L/MG = M3/G)
  ! ***  TOXPFTB(L,NT)   = TOTAL PARTICULATE FRACTION IN SEDIMENT BED (DIMENSIONLESS)
  ! ***  STFPOCB(L,K,NS) = FRACTION OF OC ON SEDIMENT CLASS NS (DIMENSIONLESS)
  ! ***  SEDB(L,K)       = COHESIVE SEDIMENTS (G/M2)
  ! ***  SNDB(L,K)       = NON-COHESIVE SEDIMENTS (G/M2)
  ! ***  STDOCB(L,K)     = DOC CONCENTRATION IN SEDIMENTS (G/M3)
  ! ***  STPOCB(L,K)     = POC CONCENTRATION IN SEDIMENTS, NON SEDIMENT COMPONENT (E.G. ALGAE) (G/M3)

  USE GLOBAL
  
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: LL,LF
  INTEGER             :: L,K,NT,NS,NX,NFD

  REAL                :: HBEDMIN0,TMPVAL

  IF( LSEDZLJ )THEN
      HBEDMIN0 = 1E-9
  ELSE
      HBEDMIN0 = 1E-4
  ENDIF
  NFD = NSED + NSND + 1

  ! *** ZERO THE TOXIC SORBED FRACTIONS
  DO NT=1,NTOX
    DO NS=1,NSP2(NT)
      DO K=1,KB
        DO L=LF,LL
          TOXPFB(L,K,NS,NT) = 0.
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  DO NT=1,NTOX

    ! *** PARTITION TO COHESIVES
    IF( ISTRAN(6) >= 1 )THEN
      IF( ISTOC(NT) > 1 )THEN
        ! *** fPOC BASED
        DO NS=1,NSED
          DO K=1,KB
            DO L=LF,LL
              ! ***     M       =    G/M2                       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SEDB(L,K,NS)*STFPOCB(L,K,NS)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ENDDO
      ELSEIF( ISTOC(NT) == 0 )THEN
        ! *** Kd APPROACH
        DO NS=1,NSED
          DO K=1,KB
            DO L=LF,LL
              ! ***     M       =    G/M2       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SEDB(L,K,NS)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    ! *** PARTITION TO NONCOHESIVES
    IF( ISTRAN(7) >= 1 )THEN
      DO NX=1,NSND
        NS=NX+NSED
        IF( ISTOC(NT) > 1 )THEN
          ! *** fPOC BASED
          DO K=1,KB
            DO L=LF,LL
              ! ***     M       =    G/M2                       L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SNDB(L,K,NX)*STFPOCB(L,K,NS)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ELSEIF( ISTOC(NT) == 0 )THEN
          ! *** Kd APPROACH
          DO K=1,KB
            DO L=LF,LL
              ! ***     M       =    G/M2      L/MG (M3/G)   
              TOXPFB(L,K,NS,NT) = SNDB(L,K,NX)*TOXPARB(NS,NT)
            ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    ! *** PARTITION (COMPLEX) TO DOC
    IF( ISTOC(NT) == 1 .OR. ISTOC(NT) == 2 )THEN
      NS = 1 + NSED + NSND
      DO K=1,KB
        DO L=LF,LL
          ! ***   M                      M          DOC(G/M3)   L/MG (M3/G)
          TOXPFB(L,K,NS,NT) = PORBED(L,K)*HBED(L,K)*STDOCB(L,K)*TOXPARBC(1,NT)
        ENDDO
      ENDDO
    ENDIF
      
    ! *** POC SORPTION (NON-SEDIMENT RELATED) COMPONENT
    IF( ISTOC(NT) == 1 )THEN
      NS = 2 + NSED + NSND
      DO K=1,KB
        DO L=LF,LL
          ! ***    M        =   M        POC(G/M3)   L/MG (M3/G)
          TOXPFB(L,K,NS,NT) = HBED(L,K)*STPOCB(L,K)*TOXPARBC(2,NT)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  
  ! ** TOXPFTB (M) IS TEMPORARILY USED TO STORE TOTAL SORBED PER TOXIC
  DO NT=1,NTOX
    DO K=1,KB
      DO L=LF,LL
        TOXPFTB(L,K,NT) = 0.
      ENDDO
    ENDDO
    DO NS=1,NSP2(NT)
      DO K=1,KB
        DO L=LF,LL
          ! *** M
          TOXPFTB(L,K,NT) = TOXPFTB(L,K,NT)+TOXPFB(L,K,NS,NT)
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  
  ! *** COMPUTE MASS FRACTION (TOXPFB) OF SORBED TOXICS FOR EACH SEDIMENT CLASS (DIMENSIONLESS)
  DO NT=1,NTOX
    DO NS=1,NSP2(NT)
      DO K=1,KB
        DO L=LF,LL
          TMPVAL = PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT)
          IF( TMPVAL > 0. )THEN
            TOXPFB(L,K,NS,NT) = TOXPFB(L,K,NS,NT)/TMPVAL
          ELSE
            TOXPFB(L,K,NS,NT) = 1.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  IF( ISTMSR >= 1 .OR. ISFDCH > 0 )THEN
    ! *** ONLY NEEDED FOR TIME SERIES USING TMSR OR FOODCHAIN
    DO NT=1,NTOX
      DO K=1,KB
        DO L=LF,LL
          IF( HBED(L,K) > 1E-12 )THEN
            TOXFDFB(L,K,NT) = PORBED(L,K)*HBED(L,K) /(PORBED(L,K)*HBED(L,K)+TOXPFTB(L,K,NT))
            TOXCDFB(L,K,NT) = TOXPFB(L,K,NFD,NT)
          ELSE
            TOXFDFB(L,K,NT) = 0.
            TOXCDFB(L,K,NT) = 0.
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

END 
