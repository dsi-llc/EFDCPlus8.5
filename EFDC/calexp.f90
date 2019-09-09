SUBROUTINE CALEXP (ISTL_)

  !**********************************************************************!
  ! **  SUBROUTINE CALEXP2T CALCULATES EXPLICIT MOMENTUM EQUATION TERMS  
  ! **  USING A TWO TIME LEVEL SCHEME  
  
  ! *** VARIABLES   DESCRIPTION                                 UNITS
  ! *** FUHU,FVHU   GROSS MOMENTUM, U COMPONENTS                M4/S2
  ! *** FVHV,FUHV   GROSS MOMENTUM, V COMPONENTS                M4/S2
  ! *** FX, FY      INTERNAL MODE FORCING BY LAYER              M4/S2
  ! *** FBBX, FBBY  INTERNAL MODE BOUYANCY FORCING BY LAYER     M4/S2
  ! *** FCAX,FCAY   CORIOLIS FORCING BY LAYER                   M4/S2
  ! *** DU, DV      INTERNAL SHEARS BY LAYER                    M2/S2
  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! ** 2018-01     PAUL M. CRAIG      CORRECTED WITHDRAWAL/RETURN MOMENTUM OPTION
  ! **                                  ADDED MOMENTUM OPTION FOR STANDARD FLOW BC
  ! ** 2016-02     PAUL M. CRAIG      UPDATED SIGMA-Z (SGZ) FOR EE8.0 
  ! ** 2015-12     PAUL M. CRAIG      ADOPTED AQEA ISHDMF>0 FOR 3TL
  ! ** 2015-06     PAUL M. CRAIG      IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  ! ** 2015-02     Paul M. Craig      UPDATED OMP AND ADDED LDRY/LWET
  ! ** 2014-01     Paul M. Craig      Fixed the TOT FXVEGE/FYVEGE when partial vegetation penetration using VEGK
  ! ** 2012-09     Dang H Chung       Added OMP
  ! ** 2011-03     John Hamrick/Scott James
  ! **                                Fixed the Vegetative Resistance from AVG to TOT using FKC
  ! ** 2010-10     Scott James        Added MHK
  !
  !----------------------------------------------------------------------C
  !
  !**********************************************************************C
  !
  USE GLOBAL

  IMPLICIT NONE
  
  INTEGER :: LU, NS, ID,  JD, KD, LD, K,  L,  LL, LW, LNW, LSE, LE, LP, IT                                                                        
  INTEGER :: ND, LF, NWR, IU, JU, KU, LN, LS, ISTL_   
  
  REAL :: QMF, QUMF, WUU, VTMPATU, UTMPATV, UMAGTMP, VMAGTMP, CACSUMT                                                                      
  REAL :: WVFACT, WVV, CFEFF, TMPVAL, DZPU, DZPV, FMDUY_TMP, FMDVX_TMP                                                       
  REAL :: VHC, VHB, DELTD2, UHC, UHB, QWRABS, VDIR                                                            

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: CACSUM  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DZPC
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: FUHJ
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: FVHJ

  IF(  .NOT. ALLOCATED(DZPC) )THEN
    ALLOCATE(CACSUM(NTHREADS))  
    ALLOCATE(FUHJ(LCM,KCM))
    ALLOCATE(FVHJ(LCM,KCM))
    ALLOCATE(DZPC(LCM,KCM))
    CACSUM=0.
    FUHJ=0.
    FVHJ=0.
    DZPC=0.
  ENDIF
  
  !**********************************************************************C
  DELT=DT2
  DELTD2=DT
  IF( ISTL_ == 2 )THEN
    DELT=DT
    DELTD2=0.5*DT
  ENDIF

  DELTI=1./DELT

  IF( N == 1 .AND. DEBUG )THEN  
    OPEN(1,FILE=OUTDIR//'MFLUX.DIA')  
    CLOSE(1,STATUS='DELETE')  
  ENDIF  

  ! *** WAVE RAMPUP FACTOR
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN
    IF( N < NTSWV )THEN  
      TMPVAL = FLOAT(N)/FLOAT(NTSWV)  
      WVFACT = 0.5-0.5*COS(PI*TMPVAL)  
    ELSE  
      WVFACT = 1.0  
    ENDIF  
  ENDIF
  CACSUMT = 0.

  !**********************************************************************!  
  ! *** ZERO NEWLY DRY CELL SHEARS/MOMENTUM
  IF( LADRY > 0 )THEN
    DO LP=1,LADRY
      L=LDRY(LP)  
      FCAXE(L)=0.  
      FCAYE(L)=0.  
      FXE(L)=0.  
      FYE(L)=0.  
    ENDDO

    DO K=1,KC
      DO LP=1,LADRY
        L=LDRY(LP)  
        FUHU(L,K)=0.  
        FVHU(L,K)=0.  
        FUHV(L,K)=0.  
        FVHV(L,K)=0.  
        FWU(L,K)=0.
        FWV(L,K)=0.
        CAC(L,K)=0.0
        FCAX(L,K)=0.
        FCAY(L,K)=0.
        FX(L,K)=0.
        FY(L,K)=0.
        FBBX(L,K)=0.
        FBBY(L,K)=0.
        DU(L,K)=0.0  
        DV(L,K)=0.0  
        
        ! *** TWO LAYER ROTATIONAL EFFECTS OR WITHDRAWAL/RETURN
        FUHJ(L,K)=0.  
        FVHJ(L,K)=0.
        
        TVAR2E(L,K)=0.
        TVAR2N(L,K)=0.
      ENDDO
    ENDDO
    
    IF( ISVEG > 0 )THEN
      DO LP=1,LADRY
        L=LDRY(LP)  
        FXVEGE(L)=0.
        FYVEGE(L)=0.
      ENDDO
      DO K=1,KC
        DO LP=1,LADRY
          L=LDRY(LP)  
          FXVEG(L,K)=0.0
          FYVEG(L,K)=0.0
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  !**********************************************************************!  
  !**********************************************************************!  
  ! *** INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS
  !----------------------------------------------------------------------!
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    DO LP=LF,LL
      L=LWET(LP)  
      FCAXE(L)=0.  
      FCAYE(L)=0.  
      FXE(L)=0.  
      FYE(L)=0.  
    ENDDO  
  ENDDO
  !$OMP END PARALLEL DO
  
  !**********************************************************************C
  ! *** SELECT ADVECTIVE FLUX FORM
  ! ***
  IF( ISTL_ == 2 )THEN
  
    ! *** THREE TIME LEVEL CORRECTOR STEP
    ! *** CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION
    ! *** AVERAGED BETWEEN (N) AND (N+1) AND ADVECTED FIELD AT N

    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          TVAR2E(L,K) = 0.5*(UHDY(L,K) + UHDY1(L,K))
          TVAR2N(L,K) = 0.5*(VHDX(L,K) + VHDX1(L,K))
        ENDDO
      ENDDO
    ENDDO
    !$OMP END DO

    !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW,UHC,UHB,VHC,VHB)
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LE=LEC(L)
          LN=LNC(L)
          LS=LSC(L)
          LW=LWC(L)

          ! *** U COMPONENTS  
          UHB = 0.5*(TVAR2E(L,K)+TVAR2E(LE,K))
          VHC = 0.5*(TVAR2N(L,K)+TVAR2N(LW,K))

          ! ***       |-- EAST FLOWING --|    |-- WEST FLOWING --|
          FUHU(L,K) = (MAX(UHB,0.)*U1(L,K)  + MIN(UHB,0.)*U1(LE,K))
          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHU(L,K) = (MAX(VHC,0.)*U1(LS,K) + MIN(VHC,0.)*U1(L,K))

          ! *** V COMPONENTS
          VHB = 0.5*(TVAR2N(L,K)+TVAR2N(LN,K))
          UHC = 0.5*(TVAR2E(L,K)+TVAR2E(LS,K))
          
          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHV(L,K) = (MAX(VHB,0.)*V1(L,K)  + MIN(VHB,0.)*V1(LN, K))
          ! ***       |-- EAST  FLOWING --|   |-- WEST  FLOWING --|
          FUHV(L,K) = (MAX(UHC,0.)*V1(LW,K) + MIN(UHC,0.)*V1(L,K))
        ENDDO
      ENDDO
    ENDDO
    !$OMP END DO

    IF( KC > 1 )THEN
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,WUU,WVV)
      DO ND=1,NDM
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LW=LWC(L)
            LS=LSC(L)
            WUU = 0.25*DXYU(L)*(W2(L,K)+W2(LW,K))
            WVV = 0.25*DXYV(L)*(W2(L,K)+W2(LS,K))
            FWU(L,K) = MAX(WUU,0.)*U1(L,K) + MIN(WUU,0.)*U1(L,K+1)
            FWV(L,K) = MAX(WVV,0.)*V1(L,K) + MIN(WVV,0.)*V1(L,K+1)
          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF
    !$OMP END PARALLEL
  
  ELSEIF( ISTL_ /= 2 .AND. ISCDMA == 0 )THEN

    ! *** THREE TIME LEVEL (LEAP-FROG) STEP
    ! *** WITH TRANSPORT AT (N) AND TRANSPORTED FIELD AT (N-1)
    ! *** UPWIND DIFFERENCE MOMENTUM ADVECTION 

    !$OMP PARALLEL DEFAULT(SHARED)
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW,UHC,UHB,VHC,VHB)
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LE=LEC(L)
          LN=LNC(L)
          LS=LSC(L)
          LW=LWC(L)

          ! *** U COMPONENTS  
          UHB = 0.5*(UHDY(L,K)+UHDY(LE,K))
          VHC = 0.5*(VHDX(L,K)+VHDX(LW,K))
          
          ! ***       |-- EAST FLOWING --|    |-- WEST FLOWING --|
          FUHU(L,K) = (MAX(UHB,0.)*U1(L,K)  + MIN(UHB,0.)*U1(LE,K))  
          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHU(L,K) = (MAX(VHC,0.)*U1(LS,K) + MIN(VHC,0.)*U1(L,K))

          ! *** V COMPONENTS
          VHB = 0.5*(VHDX(L,K)+VHDX(LN,K))
          UHC = 0.5*(UHDY(L,K)+UHDY(LS,K))

          ! ***       |-- NORTH FLOWING --|   |-- SOUTH FLOWING --|
          FVHV(L,K) = (MAX(VHB,0.)*V1(L,K)  + MIN(VHB,0.)*V1(LN, K))
          ! ***       |-- EAST  FLOWING --|   |-- WEST  FLOWING --|
          FUHV(L,K) = (MAX(UHC,0.)*V1(LW,K) + MIN(UHC,0.)*V1(L,K))
        ENDDO
      ENDDO
    ENDDO
    !$OMP END DO

    !----------------------------------------------------------------------!  
    ! *** COMPUTE VERTICAL ACCELERATIONS
    IF( KC > 1 )THEN
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,WUU,WVV)
      DO ND=1,NDM

        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LW=LWC(L)
            LS=LSC(L)
            WUU = 0.5*DXYU(L)*(W(L,K)+W(LW,K))
            WVV = 0.5*DXYV(L)*(W(L,K)+W(LS,K))
            FWU(L,K) = MAX(WUU,0.)*U1(L,K) + MIN(WUU,0.)*U1(L,K+1)
            FWV(L,K) = MAX(WVV,0.)*V1(L,K) + MIN(WVV,0.)*V1(L,K+1)
          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF
    !$OMP END PARALLEL
  
  ELSEIF( ISTL_ /= 2 .AND. ISCDMA == 1 )THEN
    
    ! *** THREE TIME LEVEL (LEAP-FROG) STEP
    ! *** AT (N) AND TRANSPORTED FIELD AT (N)
    ! *** CENTRAL DIFFERENCE MOMENTUM ADVECTION
  
    !$OMP PARALLEL DEFAULT(SHARED) 
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW)
    DO ND=1,NDM

      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LE=LEC(L)
          LN=LNC(L)
          LS=LSC(L)
          LW=LWC(L)
          
          ! *** U COMPONENTS
          FUHU(L,K) = 0.25*(UHDY(L,K) + UHDY(LE,K))*(U(L,K)+U(LE,K))
          FVHU(L,K) = 0.25*(VHDX(L,K) + VHDX(LW,K))*(U(L,K)+U(LS,K))

          ! *** V COMPONENTS
          FVHV(L,K) = 0.25*(VHDX(L,K) + VHDX(LN,K))*(V(L,K)+V(LN,K))
          FUHV(L,K) = 0.25*(UHDY(L,K) + UHDY(LS,K))*(V(L,K)+V(LW,K))
        ENDDO
      ENDDO

      IF( KC > 1 )THEN
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LW=LWC(L)
            LS=LSC(L)
            FWU(L,K) = 0.25*DXYU(L)*(W(L,K)+W(LW,K))*(U(L,K+1)+U(L,K))
            FWV(L,K) = 0.25*DXYV(L)*(W(L,K)+W(LS,K))*(V(L,K+1)+V(L,K))
          ENDDO
        ENDDO
      ENDIF
    ENDDO   ! ***  END OF DOMAIN
    !$OMP END DO
    !$OMP END PARALLEL   
  
  ELSEIF( ISTL_ /= 2 .AND. ISCDMA == 2 )THEN
    
    ! *** THREE TIME LEVEL (LEAP-FROG) STEP
    ! *** FIRST HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE
    ! *** WITH TRANSPORT AT (N-1/2) AND TRANSPORTED FIELD AT (N-1)
    ! *** SECOND HALF STEP CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE
    ! *** WITH TRANSPORT AT (N+1/2) AND TRANSPORTED FIELD AT (N)
    
    !$OMP PARALLEL DEFAULT(SHARED) 
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW,UHC,UHB,VHC,VHB,WUU,WVV)
    DO ND=1,NDM

      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          U2(L,K) = U1(L,K)+U(L,K)
          V2(L,K) = V1(L,K)+V(L,K)
        ENDDO
      ENDDO

      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LE=LEC(L)
          LN=LNC(L)
          LS=LSC(L)
          LW=LWC(L)

          ! *** U COMPONENTS
          UHB = 0.25*(UHDY(L,K)+UHDY(LE,K))
          VHC = 0.25*(VHDX(L,K)+VHDX(LW,K))
          
          FUHU(L,K) = MAX(UHB,0.)*U2(L,K)  + MIN(UHB,0.)*U2(LE,K)
          FVHU(L,K) = MAX(VHC,0.)*U2(LS,K) + MIN(VHC,0.)*U2(L,K)

          ! *** V COMPONENTS
          VHB = 0.25*(VHDX(L,K)+VHDX(LN,K))
          UHC = 0.25*(UHDY(L,K)+UHDY(LS,K))

          FVHV(L,K) = MAX(VHB,0.)*V2(L,K)  + MIN(VHB,0.)*V2(LN,K)
          FUHV(L,K) = MAX(UHC,0.)*V2(LW,K) + MIN(UHC,0.)*V2(L,K)
        ENDDO
      ENDDO
  
      IF( KC > 1 )THEN
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LW=LWC(L)
            LS=LSC(L)
            WUU = 0.25*DXYU(L)*(W(L,K)+W(LW,K))
            WVV = 0.25*DXYV(L)*(W(L,K)+W(LS,K))
            FWU(L,K) = MAX(WUU,0.)*U2(L,K) + MIN(WUU,0.)*U2(L,K+1)
            FWV(L,K) = MAX(WVV,0.)*V2(L,K) + MIN(WVV,0.)*V2(L,K+1)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL
 
  ENDIF ! *** NUMERICAL SCHEMES
  
  ! *** ADD WITHDRAWAL/RETURN FLOW MOMENTUM FLUXES
  DO NWR=1,NQWR  
    IF( ABS(NQWRMFU(NWR)) > 0 )THEN  
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NS = NQWRSERQ(NWR)  
      IF( QWRSERT(NS) >= 0. )THEN
        ! *** Original Withdrawal/Return
        IU=IQWRU(NWR)  
        JU=JQWRU(NWR)  
        KU=KQWRU(NWR)
        VDIR = 1.
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        IU=IQWRD(NWR)  
        JU=JQWRD(NWR)  
        KU=KQWRD(NWR) 
        VDIR = -1.
      ENDIF
      LU=LIJ(IU,JU)  

      QWRABS = ABS(QWRSERT(NS))
      QMF = QWR(NWR) + QWRABS 
      QUMF = VDIR*QMF*( QMF/( HPK(LU,KU)*BQWRMFU(NWR) ) )   ! *** M4/S2
      IF( NQWRMFU(NWR) ==  1 ) FUHJ(LU     ,KU) = -QUMF  
      IF( NQWRMFU(NWR) ==  2 ) FVHJ(LU     ,KU) = -QUMF  
      IF( NQWRMFU(NWR) ==  3 ) FUHJ(LEC(LU),KU) = -QUMF  
      IF( NQWRMFU(NWR) ==  4 ) FVHJ(LNC(LU),KU) = -QUMF  
      IF( NQWRMFU(NWR) == -1 ) FUHJ(LU     ,KU) = QUMF  
      IF( NQWRMFU(NWR) == -2 ) FVHJ(LU     ,KU) = QUMF  
      IF( NQWRMFU(NWR) == -3 ) FUHJ(LEC(LU),KU) = QUMF  
      IF( NQWRMFU(NWR) == -4 ) FVHJ(LNC(LU),KU) = QUMF  
    ENDIF  
    IF( ABS(NQWRMFD(NWR)) > 0 )THEN  
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NS = NQWRSERQ(NWR)  
      IF( QWRSERT(NS) >= 0. )THEN
        ! *** Original Withdrawal/Return
        ID=IQWRD(NWR)  
        JD=JQWRD(NWR)  
        KD=KQWRD(NWR)  
        VDIR = 1.
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        ID=IQWRU(NWR)  
        JD=JQWRU(NWR)  
        KD=KQWRU(NWR) 
        VDIR = -1.
      ENDIF
      LD=LIJ(ID,JD)

      QWRABS = ABS(QWRSERT(NS))
      QMF = QWR(NWR) + QWRABS 
      QUMF = VDIR*QMF*( QMF/( HPK(LD,KD)*BQWRMFD(NWR) ) )   ! *** M4/S2
      IF( NQWRMFD(NWR) ==  1 ) FUHJ(LD     ,KD) = -QUMF  
      IF( NQWRMFD(NWR) ==  2 ) FVHJ(LD     ,KD) = -QUMF  
      IF( NQWRMFD(NWR) ==  3 ) FUHJ(LEC(LD),KD) = -QUMF  
      IF( NQWRMFD(NWR) ==  4 ) FVHJ(LNC(LD),KD) = -QUMF  
      IF( NQWRMFD(NWR) == -1 ) FUHJ(LD     ,KD) = QUMF  
      IF( NQWRMFD(NWR) == -2 ) FVHJ(LD     ,KD) = QUMF  
      IF( NQWRMFD(NWR) == -3 ) FUHJ(LEC(LD),KD) = QUMF  
      IF( NQWRMFD(NWR) == -4 ) FVHJ(LNC(LD),KD) = QUMF  
    ENDIF  
  ENDDO  
  
  ! *** ADD QSER MOMENTUM FLUXES
  DO LL=1,NQSIJ
    IF( ABS(NQSMF(LL)) > 0 )THEN  
      L=LQS(LL)
      DO K=KSZ(L),KC
        ! *** Handle reversing flows in/out of domain
        IF( QSERCELL(K,LL) >= 0. )THEN
          VDIR = 1.
        ELSE
          VDIR = -1.
        ENDIF
        !LD = LIJ(IQS(LL),JQS(LL))
        
        QMF = ABS(QSERCELL(K,LL))
        QUMF = VDIR*QMF*( QMF/( HPK(L,K)*QWIDTH(LL) ) )   ! *** M4/S2
        IF( NQSMF(LL) ==  1 ) FUHJ(L     ,K) = -QUMF  
        IF( NQSMF(LL) ==  2 ) FVHJ(L     ,K) = -QUMF  
        IF( NQSMF(LL) ==  3 ) FUHJ(LEC(L),K) = -QUMF  
        IF( NQSMF(LL) ==  4 ) FVHJ(LNC(L),K) = -QUMF  
        IF( NQSMF(LL) == -1 ) FUHJ(L     ,K) = QUMF  
        IF( NQSMF(LL) == -2 ) FVHJ(L     ,K) = QUMF  
        IF( NQSMF(LL) == -3 ) FUHJ(LEC(L),K) = QUMF  
        IF( NQSMF(LL) == -4 ) FVHJ(LNC(L),K) = QUMF  
      ENDDO
    ENDIF  
  ENDDO

  !$OMP PARALLEL DEFAULT(SHARED)     
  ! **********************************************************************!  
  ! ** BLOCK MOMENTUM FLUX ON LAND SIDE OF TRIANGULAR CELLS  
  IF( ITRICELL > 0 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          FUHU(L,K) = STCUV(L)*FUHU(L,K)
          FVHV(L,K) = STCUV(L)*FVHV(L,K)
        ENDDO
      ENDDO
    ENDDO   ! ***  END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  !**********************************************************************C
  !
  ! *** CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS
  !$OMP SINGLE
  CACSUM=0. 
  CFMAX=CF  
  !$OMP END SINGLE

  IF( ISCURVATURE )THEN
    IF( ISDCCA == 0 )THEN
      ! *** STANDARD CALCULATIONS, NO DIAGNOSTICS
      !$OMP DO PRIVATE(ND,K,LP,L,LE,LN)
      DO ND=1,NDM  
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LE=LEC(L)
            LN=LNC(L)
            CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)
            CACSUM(ND) = CACSUM(ND)+CAC(L,K)
          ENDDO
        ENDDO
      ENDDO   ! ***  END OF DOMAIN
      !$OMP END DO
      
    ELSE
      ! *** STANDARD CALCULATIONS, WITH DIAGNOSTICS
      !$OMP SINGLE
      IT = 1
      DO K=1,KC
        DO L=2,LA
          LE = LEC(L)
          LN = LNC(L)
          CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)
          CFEFF = ABS(CAC(L,K))*DXYIP(L)*HPI(L)
          CFMAX = MAX(CFMAX,CFEFF)
          CACSUM(IT) = CACSUM(IT)+CAC(L,K)
        ENDDO
      ENDDO
      !$OMP END SINGLE
    ENDIF

    ! *** ENSURE FCAY & FCAX ARE RESET
    !$OMP SINGLE
    CACSUMT = ABS(SUM(CACSUM(:)))
    !$OMP END SINGLE

    IF( CACSUMT < 1.E-7 )THEN
      !$OMP DO PRIVATE(ND,K,LP,L)
      DO ND=1,NDM  
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            FCAX(L,K) = 0.
            FCAY(L,K) = 0.
          ENDDO
        ENDDO
      ENDDO  
      !$OMP END DO
    ENDIF
  ENDIF

  !**********************************************************************C
  !
  ! *** CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS
  !
  !----------------------------------------------------------------------C

  ! **  STANDARD CALCULATION
  IF( CACSUMT > 1.E-7 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LE,LN,LS,LW,LNW,LSE)
    DO ND=1,NDM

      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)  
          LS = LSC(L)  
          LW = LWC(L)
          LNW = LNWC(L)  
          LSE = LSEC(L)  
          FCAX(L,K) = ROLD*FCAX(L,K) + 0.25*RNEW*SCAX(L)*(CAC(L,K)*(V(LN,K)+V(L,K)) + CAC(LW,K)*(V(LNW,K)+V(LW,K)))
          FCAY(L,K) = ROLD*FCAY(L,K) + 0.25*RNEW*SCAY(L)*(CAC(L,K)*(U(LE,K)+U(L,K)) + CAC(LS,K)*(U(LSE,K)+U(LS,K)))
        ENDDO
      ENDDO
    ENDDO
    !$OMP END DO

    !----------------------------------------------------------------------C
    !
    ! *** MODIFICATION FOR TYPE 2 OPEN BOUNDARIES
    !$OMP SINGLE
    DO LL=1,NPBW
      IF( ISPBW(LL) == 2 )THEN
        L = LPBW(LL)+1
        LN = LNC(L)
        DO K=KSZ(L),KC
          FCAX(L,K) = 0.5*SCAX(L)*CAC(L,K)*(V(LN,K)+V(L,K))
        ENDDO
      ENDIF
    ENDDO
  
    DO LL=1,NPBE
      IF( ISPBE(LL) == 2 )THEN
        L = LPBE(LL)
        LNW = LNWC(L)
        DO K=KSZ(L),KC
          FCAX(L,K) = 0.5*SCAX(L)*CAC(LWC(L),K)*(V(LNW,K)+V(LWC(L),K))
        ENDDO
      ENDIF
    ENDDO
  
    DO LL=1,NPBS
      IF( ISPBS(LL) == 2 )THEN
        L = LNC(LPBS(LL))
        DO K=KSZ(L),KC
          FCAY(L,K) = 0.5*SCAY(L)*CAC(L,K)*(U(LEC(L),K)+U(L,K))
        ENDDO
      ENDIF
    ENDDO
  
    DO LL=1,NPBN
      IF( ISPBN(LL) == 2 )THEN
        L = LPBN(LL)
        LS = LSC(L)
        LSE = LSEC(L)
        DO K=KSZ(L),KC
          FCAY(L,K) = 0.5*SCAY(L)*CAC(LS,K)*(U(LSE,K)+U(LS,K))
        ENDDO
      ENDIF
    ENDDO
    !$OMP END SINGLE

  ENDIF    ! *** END OF CACSUMT > 1.E-7

  !----------------------------------------------------------------------C
  !  INITIALIZE EXTERNAL FORCINGS WITH GROSS MOMENTUM
  !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW)
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        LE = LEC(L)
        LN = LNC(L)
        LS = LSC(L)
        LW = LWC(L)
        FX(L,K) = FSGZU(L,K)*( FUHU(L,K)-FUHU(LW,K) + FVHU(LN,K)-FVHU(L,K) ) + FUHJ(L,K)   ! ***  M4/S2
        FY(L,K) = FSGZV(L,K)*( FVHV(L,K)-FVHV(LS,K) + FUHV(LE,K)-FUHV(L,K) ) + FVHJ(L,K)   ! ***  M4/S2
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! *** TREAT BC'S NEAR EDGES
  !$OMP SINGLE
  DO LL=1,NBCS
    ! *** BC CELL
    L = LBCS(LL)
    DO K=KSZ(L),KC
      FX(L,K) = SAAX(L)*FX(L,K)
      FY(L,K) = SAAY(L)*FY(L,K)
    ENDDO

    ! *** EAST/WEST ADJACENT CELL
    L = MAX(1,LBERC(LL))
    DO K=KSZ(L),KC
      FX(L,K) = SAAX(L)*FX(L,K)
    ENDDO

    ! *** NORTH/SOUTH ADJACENT CELL
    L = MAX(1,LBNRC(LL))
    DO K=KSZ(L),KC
      FY(L,K) = SAAY(L)*FY(L,K)
    ENDDO
  ENDDO
  !$OMP END SINGLE

  !**********************************************************************!
  !
  ! *** ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS
  !
  !----------------------------------------------------------------------!
  IF( ISVEG > 0 )THEN
    ! *** ADD IN MHK DEVICES, IF NEEDED
    !$OMP SINGLE
    IF( LMHK )CALL MHKPWRDIS
    !$OMP END SINGLE

    IF( ISDRY > 0 .AND. LADRY > 0 )THEN
      ! *** UPDATE ACTIVE CELL BY LAYER LIST
      !$OMP DO PRIVATE(ND,K,LN,LP,L)
      DO ND=1,NDM  
        DO K=1,KC  
          LN=0
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            IF( LVEG(L) )THEN
              LN = LN+1
              LKVEG(LN,K,ND) = L
            ENDIF
          ENDDO
          LLVEG(K,ND)=LN
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF
    
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LW,LE,LS,LN,LNW,LSE) &
    !$OMP    PRIVATE(UTMPATV,VTMPATU,UMAGTMP,VMAGTMP)
    DO ND=1,NDM  

      DO LP=1,LLVEG(KC,ND)
        L=LKVEG(LP,KC,ND) 
        FXVEGE(L) = 0.  
        FYVEGE(L) = 0.  
      ENDDO  

      DO K=1,KC  
        DO LP=1,LLVEG(K,ND)
          L=LKVEG(LP,K,ND)
          LW = LWC(L)    !west cell
          LE = LEC(L)    !east cell
          LS = LSC(L)    !south cell
          LN = LNC(L)    !north cell
          LNW = LNWC(L)  !northwest cell
          LSE = LSEC(L)  !southeast cell
          UTMPATV = 0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))       !u-velocity at v face
          VTMPATU = 0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))       !v-velocity at u face
          UMAGTMP = SQRT( U(L,K)*U(L,K) +   VTMPATU*VTMPATU )    !u-face velocity vector magnitude
          VMAGTMP = SQRT( UTMPATV*UTMPATV + V(L,K)*V(L,K) )      !v-face velocity vector magnitude

          !FXVEG/FYVEG come from CALTBXY unitless, but they are really just a form of the drag coefficient with terms accounting for the area density
          !FXVEG/FYVEG only change inasmuch as the water depth changes and are zero in layers not penetrated by vegetation
          !FXVEG/FYVEG are C_d(N/L^2)
          !FXVEG/FYVEG are now multiplied by the cell area and cell-averaged velocity
          !FXVEG/FYVEG are C_d(N/L^2)A|q|
          FXVEG(L,K) = UMAGTMP*SUB3D(L,K)*FXVEG(L,K)  ![m/s] q_xC_d
          FYVEG(L,K) = VMAGTMP*SVB3D(L,K)*FYVEG(L,K)  ![m/s] q_yC_d

          !FXVEG/FXVEGE are multiplied by the local velocity to yield units of [m^4/s^2]
          !FXVEG/FXVEGE are added to the body forces as C_d(N/L^2)A|q|q
          FXVEGE(L) = FXVEGE(L) + FXVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FUHDXE
          FYVEGE(L) = FYVEGE(L) + FYVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FVHDYE
        ENDDO
      ENDDO

      ! *** ADD VEGETATIVE DRAG TO INTERNAL MODE SHEAR AND EXTERNAL MODE MOMENTUM
      DO K=1,KC  
        DO LP=1,LLVEG(K,ND)
          L=LKVEG(LP,K,ND)
          IF( (K-KSZ(L)+1) > INT(VEGK(L)+1.) )CYCLE
          FX(L,K) = FX(L,K) + (FXVEG(L,K)-FXVEGE(L))*U(L,K)*DXYU(L) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) FXVEGE goes into FUHDXE for momentum conservation
          FY(L,K) = FY(L,K) + (FYVEG(L,K)-FYVEGE(L))*V(L,K)*DXYV(L) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) FYVEGE goes into FVHDYE for momentum conservation
        ENDDO
      ENDDO
    
      ! *** CONVERT THE AVG FXVEGE/FYVEGE TO TOTAL FXVEGE/FYVEGE
      DO LP=1,LLVEG(KC,ND)
        L=LKVEG(LP,KC,ND) 
        FXVEGE(L) = FXVEGE(L)*HUI(L)*VEGK(L) !Calculate vegetative dissipation for FUHDYE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
        FYVEGE(L) = FYVEGE(L)*HVI(L)*VEGK(L) !Calculate vegetative dissipation for FVHDXE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
      ENDDO
  
      IF( LMHK )THEN
        DO LP=1,LLVEG(KC,ND)
          L=LKVEG(LP,KC,ND) 
          FXVEGE(L) = FXVEGE(L) + FXMHKE(L)   ! Add MHK to vegetative dissipation in FUHDYE for momentum conservation in CALPUV multiply by HUI took place in MHKPWRDIS
          FYVEGE(L) = FYVEGE(L) + FYMHKE(L)   ! Add MHK to vegetative dissipation in FVHDXE for momentum conservation in CALPUV multiply by HVI took place in MHKPWRDIS
          FXVEGE(L) = FXVEGE(L) + FXSUPE(L)   ! Add MHK support to vegetative dissipation in FUHDYE for momentum conservation in CALPUV
          FYVEGE(L) = FYVEGE(L) + FYSUPE(L)   ! Add MHK support to vegetative dissipation in FVHDXE for momentum conservation in CALPUV
        ENDDO
      ENDIF
      
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF

  !**********************************************************************!
  !
  ! *** ADD HORIZONTAL MOMENTUM DIFFUSION TO ADVECTIVE ACCELERATIONS
  !
  !----------------------------------------------------------------------!
  IF( ISHDMF >= 1 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LS,LN,LW,FMDUY_TMP,FMDVX_TMP)
    DO ND=1,NDM  
      DO K=1,KC
        DO LP=1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          ! *** THE FOLLOWING IS FROM AQEA 2015 CODE
          LS=LSC(L)
          LN=LNC(L)  
          LW=LWC(L)
          LE=LEC(L)
          ! ** CALCULATE DIFFUSIVE FLUX ON THE NORTHERN EXTERNAL BOUNDARY  
          IF( SVB(LN) < 0.5 ) THEN
            ! ** DXV1 IS 0 ON NORTHERN CLOSED BOUNDARY
            ! ** DYU1 IS CALCULATED USING THE U1 CLOSEST TO THE WALL AND 0
            FMDUY_TMP = DXU(L)*H1C(L)*AHC(L,K)*(-1.*U1(L,K)/DYU(L))
            FX(L,K) = FX(L,K) - SUB3D(L,K)*SDX(L)*( FMDUX(L,K) - FMDUX(LW,K) + FMDUY_TMP   - FMDUY(L,K) )
          ELSE
            FX(L,K) = FX(L,K) - SUB3D(L,K)*SDX(L)*( FMDUX(L,K) - FMDUX(LW,K) + FMDUY(LN,K) - FMDUY(L,K) )
          END IF
            
          ! ** CALCULATE DIFFUSIVE FLUX ON THE EASTERN EXTERNAL BOUNDARY  
          IF( SUB(LE) < 0.5 )THEN
            ! ** DXV1 IS 0 IS CALCUALATED USING V1 CLOSEST TO THE WALL AND 0
            ! ** DYU1 IS 0 ON EASTERN CLOSED BOUNDARY
            FMDVX_TMP = (DYU(L))*H1C(L)*AHC(L,K)*((-1.*V1(L,K)/DXU(L)))
            FY(L,K) = FY(L,K) - SVB3D(L,K)*SDY(L)*( FMDVY(L,K) - FMDVY(LS,K) + FMDVX_TMP   - FMDVX(L,K) )
          ELSE 
            FY(L,K) = FY(L,K) - SVB3D(L,K)*SDY(L)*( FMDVY(L,K) - FMDVY(LS,K) + FMDVX(LE,K) - FMDVX(L,K) )
          END IF
        ENDDO  
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  !**********************************************************************C
  !
  ! *** ADD BODY FORCE TO ADVECTIVE ACCELERATIONS
  ! *** DISTRIBUTE UNIFORMLY OVER ALL LAYERS IF ISBODYF=1
  ! *** DISTRIBUTE OVER SURFACE LAYER IF ISBODYF=2
  !
  !----------------------------------------------------------------------C
  IF( ISBODYF == 1 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          FX(L,K) = FX(L,K) - DYU(L)*HU(L)*FBODYFX(L,K)
          FY(L,K) = FY(L,K) - DXV(L)*HV(L)*FBODYFY(L,K)
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  IF( ISBODYF == 2 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      DO LP=LF,LL
        L=LWET(LP)
        FX(L,KC) = FX(L,KC) - DZIC(L,KC)*DYU(L)*HU(L)*FBODYFX(L,KC)
        FY(L,KC) = FY(L,KC) - DZIC(L,KC)*DXV(L)*HV(L)*FBODYFY(L,KC)
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  !**********************************************************************C
  !
  ! *** ADD EXPLICIT NONHYDROSTATIC PRESSURE
  !
  IF( KC > 1 .AND. ISPNHYDS >= 1 )THEN
    !$OMP DO PRIVATE(ND,L,LP,K,TMPVAL) 
    DO ND=1,NDM  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        TMPVAL    = 2./( DZC(L,KSZ(L)) + DZC(L,KSZ(L)+1) )
        DZPC(L,1) = TMPVAL*(PNHYDS(L,2)-PNHYDS(L,1))
      ENDDO
    ENDDO  
    !$OMP END DO
  
    !$OMP DO PRIVATE(ND,L,LP,TMPVAL) 
    DO ND=1,NDM  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        TMPVAL     = 2./( DZC(L,KC) + DZC(L,KC-1) )
        DZPC(L,KC) = TMPVAL*(PNHYDS(L,KC)-PNHYDS(L,KC-1))
      ENDDO
    ENDDO  
    !$OMP END DO

    IF( KC >= 3 )THEN
      !$OMP DO PRIVATE(ND,K,LP,L,TMPVAL)
      DO ND=1,NDM  
        DO K=2,KS
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            TMPVAL    = 2./( DZC(L,K+1)+2.*DZC(L,K)+DZC(L,K-1) )
            DZPC(L,K) = TMPVAL*(PNHYDS(L,K+1)-PNHYDS(L,K-1))
          ENDDO
        ENDDO
      ENDDO  
      !$OMP END DO
    ENDIF
  
   !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,DZPU,DZPV)
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LW=LWC(L)
          LS=LSC(L)
          DZPU=0.5*(DZPC(L,K)+DZPC(LW,K))
          DZPV=0.5*(DZPC(L,K)+DZPC(LS,K))
          FX(L,K) = FX(L,K) + SUB3D(L,K)*DYU(L)*( HU(L)*(PNHYDS(L,K)-PNHYDS(LW,K) ) - ( BELV(L)-BELV(LW) + ZZ(L,K) *(HP(L)-HP(LW)) )*DZPU )
          FY(L,K) = FY(L,K) + SVB3D(L,K)*DXV(L)*( HV(L)*(PNHYDS(L,K)-PNHYDS(LS,K) ) - ( BELV(L)-BELV(LS) + ZZ(LS,K)*(HP(L)-HP(LS)) )*DZPV )
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  IF( IATMP > 0 )THEN
    ! *** ADD AIR PRESSURE GRADIENT
    !$OMP DO PRIVATE(ND,K,LP,L,LS,LW)
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LW=LWC(L)
          LS=LSC(L)
          !                                M       M         M2/S2        
          FX(L,K) = FX(L,K) + SUB3D(L,K)*DYU(L)*HU(L)*( ATMP(L)-ATMP(LW) )
          FY(L,K) = FY(L,K) + SVB3D(L,K)*DXV(L)*HV(L)*( ATMP(L)-ATMP(LS) )
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF

  !----------------------------------------------------------------------!
  !
  ! *** ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.
  !
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN
  
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L) 
    DO ND=1,NDM  
      LF=(ND-1)*NWVCELLS+1  
      LL=MIN(LF+NWVCELLS-1,NWVCELLS)
      
      DO K=1,KC  
        DO LP=LF,LL
          L = LWVCELL(LP)
          IF( LKSZ(L,K) )CYCLE  
          FX(L,K) = FX(L,K) + SUB3D(L,K)*WVFACT*SAAX(L)*FXWAVE(L,K)
          FY(L,K) = FY(L,K) + SVB3D(L,K)*WVFACT*SAAY(L)*FYWAVE(L,K)
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  
  ENDIF
  
  !**********************************************************************!
  !
  ! *** CALCULATE TOTAL EXTERNAL ACCELERATIONS
  !
  !----------------------------------------------------------------------!
  !$OMP DO PRIVATE(ND,K,LP,L)
  DO ND=1,NDM  
      
    DO K=1,KC  
      DO LP=1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        FCAXE(L) = FCAXE(L) + FCAX(L,K)*SGZU(L,K)
        FCAYE(L) = FCAYE(L) + FCAY(L,K)*SGZV(L,K)
        FXE(L) = FXE(L) + FX(L,K)*SGZU(L,K)
        FYE(L) = FYE(L) + FY(L,K)*SGZV(L,K)
      ENDDO
    ENDDO
  ENDDO   ! *** END OF DOMAIN
  !$OMP END DO

  !**********************************************************************!  
  ! **  COMPLETE CALCULATION OF INTERNAL MODE ADVECTIVE ACCELERATIONS  
  IF( KC > 1 )THEN
    ! *** LIMIT THE VERTICAL MOMENTUM AT THE KSZ DIFFERENCES
  
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          FX(L,K) = FX(L,K) + SAAX(L)*(FWU(L,K)-FWU(L,K-1))*DZIC(L,K)
          FY(L,K) = FY(L,K) + SAAY(L)*(FWV(L,K)-FWV(L,K-1))*DZIC(L,K)
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  !**********************************************************************C
  !
  ! *** CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR
  ! *** THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP
  ! *** SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV
  !
  !----------------------------------------------------------------------C
  IF( BSC > 1.E-6 .AND. KC > 1 )THEN

    IF( IGRIDV == 1 )THEN
      ! *** SIGMA-ZED BOUYANCY SHEARS
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      DO ND=1,NDM  
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LS=LSC(L) 
            LW=LWC(L)
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SUB3D(L,K)*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*SGZU(L,K+1) + (B(L,K)-B(LW,K))*SGZU(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SVB3D(L,K)*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*SGZV(L,K+1) + (B(L,K)-B(LS,K))*SGZV(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO
      
    ELSEIF( IGRIDV == 2 )THEN
      ! *** SIGMA-ZED BOUYANCY SHEARS
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      DO ND=1,NDM
        ! *** ALL LAYERS
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND) 
            LS = LSC(L) 
            LW = LWC(L)
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SUB3D(L,K)*SBX(L)*GP*HU(L)*( BW(L,K+1)*HPW(L)*SGZW(L,K+1) - BE(LW,K+1)*HPE(LW)*SGZE(LW,K+1) + BW(L,K)*HPW(L)*SGZW(L,K) - BE(LW,K)*HPE(LW)*SGZE(LW,K)   &
                                                    - (BW(L,K+1)-BW(L,K)+BE(LW,K+1)-BE(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SVB3D(L,K)*SBY(L)*GP*HV(L)*( BS(L,K+1)*HPS(L)*SGZS(L,K+1) - BN(LS,K+1)*HPN(LS)*SGZN(LS,K+1) + BS(L,K)*HPS(L)*SGZS(L,K) - BN(LS,K)*HPN(LS)*SGZN(LS,K)   &
                                                    - (BS(L,K+1)-BS(L,K)+BN(LS,K+1)-BN(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          ENDDO  
        ENDDO
      ENDDO
      !$OMP END DO
      

    ELSEIF( IINTPG == 0 )THEN  
      ! *** IINTPG=0  
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LS,LW) 
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO K=1,KS  
          DO LP=LF,LL
            L = LWET(LP)  
            LS = LSC(L)
            LW = LWC(L)
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*DZCK(K+1) + (B(L,K)-B(LW,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*(BELV(L)-BELV(LW)+Z(L,K)*(HP(L)-HP(LW))) )
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZCK(K+1) + (B(L,K)-B(LS,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*(BELV(L)-BELV(LS)+Z(L,K)*(HP(L)-HP(LS))) )
          ENDDO
        ENDDO
      ENDDO
      !$OMP END DO  
  
    ELSEIF( IINTPG == 1 )THEN
      ! *** JACOBIAN
      !$OMP SINGLE
      K=1
      DO L=2,LA
        LW=LWC(L)
        LS=LSC(L)
        FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( 0.5*HU(L)*( (B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K  )-B(LW,K  ))*DZC(L,K  )+(B(L,K  )-B(LW,K  ))*DZC(L,K  ) ) &
                   -0.5*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))*(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.5*(B(L,K  )-B(L,K  )+B(LW,K  )-B(LW,K  ))*(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
  
        FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( 0.5*HV(L)*( (B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K  )-B(LS ,K  ))*DZC(L,K  )+(B(L,K  )-B(LS ,K  ))*DZC(L,K  ) ) &
                   -0.5*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.5*(B(L,K  )-B(L,K  )+B(LS ,K  )-B(LS ,K  ))*(BELV(L)-BELV(LS )+Z(L,K-1)*(HP(L)-HP(LS ))) )
      ENDDO
  
      IF( KC > 2 )THEN
        K=KS
        DO L=2,LA
          LW=LWC(L)
          LS=LSC(L)
          FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( 0.5*HU(L)*( (B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K  )-B(LW,K  ))*DZC(L,K  )+(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) ) &
                     -0.5*(B(L,K+1)-B(L,K+1)+B(LW,K+1)-B(LW,K+1))*(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.5*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1))*(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
          FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( 0.5*HV(L)*( (B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K  )-B(LS ,K  ))*DZC(L,K  )+(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) ) &
                     -0.5*(B(L,K+1)-B(L,K+1)+B(LS ,K+1)-B(LS ,K+1))*(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.5*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*(BELV(L)-BELV(LS )+Z(L,K-1)*(HP(L)-HP(LS ))) )
        ENDDO
      ENDIF
  
      IF( KC > 3 )THEN
        DO K=1,KS
          DO L=2,LA
            LW=LWC(L)
            LS=LSC(L)
            FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( 0.5*HU(L)*( (B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K  )-B(LW,K  ))*DZC(L,K  )+(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) ) &
                       -0.5*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))*(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.5*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1))*(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
            FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( 0.5*HV(L)*( (B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+(B(L,K  )-B(LS ,K  ))*DZC(L,K  )+(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) ) &
                       -0.5*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))*(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.5*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))*(BELV(L)-BELV(LS )+Z(L,K-1)*(HP(L)-HP(LS ))) )
          ENDDO
        ENDDO
      ENDIF
      !$OMP END SINGLE
        
    ELSEIF( IINTPG == 2 )THEN
      ! *** FINITE VOLUME
      !$OMP SINGLE
      DO K=1,KS
        DO L=2,LA
          LW=LWC(L)
          LS=LSC(L)
          FBBX(L,K) = ROLD*FBBX(L,K) + RNEW*SBX(L)*GP*HU(L)*( ( HP(L)*B(L,K+1)-HP(LW)*B(LW,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LW)*B(LW,K  ) )*DZC(L,K  ) )-RNEW*SBX(L)*GP*(BELV(L)-BELV(LW)) &
                     *( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LW)*B(LW,K+1)-HP(LW)*B(LW,K) )   - RNEW*SBX(L)*GP*(HP(L)-HP(LW))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K)+HP(LW)*ZZ(L,K+1)*B(LW,K+1)-HP(LW)*ZZ(L,K)*B(LW,K) )
          FBBY(L,K) = ROLD*FBBY(L,K) + RNEW*SBY(L)*GP*HV(L)*( ( HP(L)*B(L,K+1)-HP(LS )*B(LS ,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LS )*B(LS ,K  ) )*DZC(L,K  ) )-RNEW*SBY(L)*GP*(BELV(L)-BELV(LS )) &
                     *( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LS)*B(LS ,K+1)-HP(LS)*B(LS ,K) ) - RNEW*SBY(L)*GP*(HP(L)-HP(LS ))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K)+HP(LS)*ZZ(L,K+1)*B(LS ,K+1)-HP(LS)*ZZ(L,K)*B(LS ,K) )
        ENDDO
      ENDDO
      !$OMP END SINGLE
    ENDIF

  ENDIF  ! *** END OF BOUYANCY

  !**********************************************************************!
  !
  ! **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
  !
  !----------------------------------------------------------------------!
  IF( KC > 1 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K)
    DO ND=1,NDM  
      ! *** COMPUTE THE INTERNAL SHEARS FOR THE LOWER LAYERS
      DO K=1,KS
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          DU(L,K) = CDZFU(L,K)*( H1U(L)*(U1(L,K+1)-U1(L,K))*DELTI + DXYIU(L)*(FCAX(L,K+1)-FCAX(L,K)   + FBBX(L,K) + SNLT*(FX(L,K)-FX(L,K+1))) )
          DV(L,K) = CDZFV(L,K)*( H1V(L)*(V1(L,K+1)-V1(L,K))*DELTI + DXYIV(L)*(FCAY(L,K)  -FCAY(L,K+1) + FBBY(L,K) + SNLT*(FY(L,K)-FY(L,K+1))) )
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF

  ! *** ADD WIND SHEAR TO THE KC/KS INTERFACE
  IF( ISTL_ == 2 .AND. NWSER > 0 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    DO ND=1,NDM
      DO LP=LF,LL
        L=LWET(LP)  
        DU(L,KS) = DU(L,KS) - CDZUU(L,KS)*TSX(L)
        DV(L,KS) = DV(L,KS) - CDZUV(L,KS)*TSY(L)
      ENDDO
    ENDDO    ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  !$OMP END PARALLEL
  
  !**********************************************************************!

  RETURN
END
