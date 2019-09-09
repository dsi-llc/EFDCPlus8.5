SUBROUTINE CALTRAN (ISTL_,IS2TL_,MVAR,MO,CON,CON1,IT)  

  ! **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE  
  ! **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO  
  ! **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES  
  ! **  THE NUMBER OF TIME LEVELS IN THE STEP  

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2015-01       PAUL M. CRAIG     Added fully coupled Ice Sub-model with Frazil Ice Transport
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH

  USE GLOBAL  
  !  USE OMP_LIB
  
  IMPLICIT NONE  
    
  INTEGER, INTENT(IN) :: ISTL_,IS2TL_,MVAR,MO,IT  
  REAL    :: CON(LCM,KCM),CON1(LCM,KCM)  
  
  REAL    :: BSMALL, DDELT, DDELTA, DDELTD2  
  REAL    :: CTMP, CBT, AUHU, AVHV, UTERM, VTERM, WTERM  
  REAL    :: CBSTMP, CBWTMP, CBETMP, CBNTMP, UHU, VHV, AWW, WW  
  REAL    :: CWMAX, CEMAX, CSMAX, CNMAX, CMAXT  
  REAL    :: CWMIN, CEMIN, CSMIN, CNMIN, CMINT  
  
  INTEGER :: M
  INTEGER :: ISUD, K, NSID, IOBC, NMNLOD
  INTEGER :: LP,L, LN, LS, LE, LW, LSE, LNW, LL, ITRANFLOC  

  ! *** SET UP FLOC TRANSPORT
  ITRANFLOC=0
  
  BSMALL=1.0E-6  
  ISUD=1  
  IF( ISDYNSTP == 0 )THEN 
    ! *** FIXED DELTA T
    DDELT=DT2  
    DDELTA=DT2  
    IF( ISCDCA(MVAR) == 2 ) DDELTA=DT   ! *** Central Differencing (3TL) [EXPERIMENTAL]
    DDELTD2=DT
    IF( ISCDCA(MVAR) == 1 ) ISUD=0 
    IF( ISTL_/=3 )THEN  
      DDELT=DT  
      DDELTA=DT  
      DDELTD2=0.5*DT  
      IF( IS2TIM == 0 ) ISUD=0          ! *** 3TL CORRECTOR TIME STEP (ISTL=2)
    ENDIF  
  ELSE  
    ! *** DYNAMIC DELTA T
    DDELT=DTDYN  
    DDELTA=DTDYN  
    DDELTD2=0.5*DTDYN  
  END IF  
    
  M=MO  
  IF( IS2TL_ == 1 )THEN  
    ! *** ADVANCE CONCENTRATIONS BEFORE THE 2TL TRANSPORT CALCULATIONS
    ! *** SKIP UPDATING VARIABLES IF ALREADY COMPLETED BEFORE THIS STEP
    IF( MVAR /= 8 )THEN  
      CON1=CON 
    ENDIF  
  ENDIF  

  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC)  
    DO K=1,KC  
      WQBCCON(IOBC,K,IT)  = CON(L,K)  
      WQBCCON1(IOBC,K,IT) = CON1(L,K)  
    ENDDO  
  ENDDO  
  
  ! **  CALCULATED EXTERNAL SOURCES AND SINKS  
  CALL CALFQC (ISTL_,IS2TL_,MVAR,M,CON,CON1,IT)  

  
  ! **  SELECT TRANSPORT OPTION, ISPLIT=1 FOR HORIZONTAL-VERTICAL  
  ! **  OPERATOR SPLITTING  
  ! **  BEGIN COMBINED ADVECTION SCHEME  
  ! **  ADVECTIVE FLUX CALCULATION  

  IF( ISTL_ == 2 ) GOTO 300  
  IF( ISCDCA(MVAR) == 0 ) GOTO 300   ! *** Upwind Differencing  (3TL)
  IF( ISCDCA(MVAR) == 1 ) GOTO 400   ! *** Central Differencing (3TL)  
  IF( ISCDCA(MVAR) == 2 ) GOTO 350   ! *** Upwind Differencing  (3TL) [EXPERIMENTAL]
  ! *** UHDY2 AND VHDX2 ARE LAYER FLOWS (SIGMA-Z VERSION)

  ! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION  
  ! **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED  
  ! **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY  

  300 CONTINUE  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      FUHUD(L,K,IT) = UHDY2(L,K)*CON1(LUPU(L,K),K)  
      FVHUD(L,K,IT) = VHDX2(L,K)*CON1(LUPV(L,K),K)  
    ENDDO  
  ENDDO
  IF( KC > 1 )THEN  
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        FWUU(L,K,IT) = W2(L,K)*CON1(L,KUPW(L,K))  
      ENDDO  
    ENDDO  
  ENDIF  
  GOTO 500  

  ! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION  
  ! **  AVERAGED BETWEEN  (N-1) AND (N+1) AND ADVECTED FIELD AVERAGED  
  ! **  BETWEEN AT (N-1) AND (N) IF ISTL 3 ONLY  

  350 CONTINUE  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      CONTD(L,K,IT) = 0.5*(CON(L,K)+CON1(L,K)) + DDELT*0.5*FQC(L,K,IT)*DXYIP(L)/H2P(L)  
    ENDDO  
  ENDDO  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      FUHUD(L,K,IT) = UHDY2(L,K)*CONTD(LUPU(L,K),K,IT)  
      FVHUD(L,K,IT) = VHDX2(L,K)*CONTD(LUPV(L,K),K,IT)  
    ENDDO  
  ENDDO  
  IF( KC > 1 )THEN  
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        FWUU(L,K,IT) = W2(L,K)*CONTD(L,KUPW(L,K),IT)  
      ENDDO  
    ENDDO  
  ENDIF  
  GOTO 500  

  ! **  CALCULATE ADVECTIVE FLUXES BY CENTRAL DIFFERENCE WITH TRANSPORT  
  ! **  AVERAGED BETWEEN (N+1) AND (N-1) AND TRANSPORTED FIELD AT (N)  

  400 CONTINUE  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      LS=LSC(L)  
      FUHUD(L,K,IT) = 0.5*UHDY2(L,K)*(CON(L,K)+CON(LWC(L),K))  
      FVHUD(L,K,IT) = 0.5*VHDX2(L,K)*(CON(L,K)+CON(LS,K))  
    ENDDO  
  ENDDO  
  DO K=1,KC  
    DO LL=1,NCBS  
      L=LCBS(LL)  
      LN=LNC(L)  
      IF( VHDX2(LN,K) < 0.) FVHUD(LN,K,IT) = VHDX2(LN,K)*CON1(LN,K)  
    ENDDO  
    DO LL=1,NCBW  
      L=LCBW(LL)  
      IF( UHDY2(LEC(L),K) < 0.) FUHUD(LEC(L),K,IT) = UHDY2(LEC(L),K)*CON1(LEC(L),K)  
    ENDDO  
    DO LL=1,NCBE  
      L=LCBE(LL)  
      IF( UHDY2(L,K) > 0.) FUHUD(L,K,IT) = UHDY2(L,K)*CON1(LWC(L),K)  
    ENDDO  
    DO LL=1,NCBN  
      L=LCBN(LL)  
      LS =LSC(L)  
      IF( VHDX2(L,K) > 0.) FVHUD(L,K,IT) = VHDX2(L,K)*CON1(LS,K)  
    ENDDO  
  ENDDO  
  DO K=1,KS  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0) 
      FWUU(L,K,IT) = 0.5*W2(L,K)*(CON(L,K+1)+CON(L,K))  
    ENDDO  
  ENDDO  

  ! **  STANDARD ADVECTION CALCULATION  
  500 CONTINUE  

  ! *** CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX (PMC MOVED)  
  IF( ISHDMF == 2 ) CALL CALDIFF (CON1,IT)
    
  ! *** BEGIN IF ON TRANSPORT OPTION CHOICE  
  ! *** IF ISACAC EQ 0 INCLUDE FQC MASS SOURCES IN UPDATE  
  IF( ISCDCA(MVAR) == 0 )THEN
    ! *** Upwind Differencing (3TL & 2TL)  
    IF( ISTL_ == 2 )THEN  
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          CD(L,K,IT) = CON1(L,K)*H1PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                       FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT) ) * DXYIP(L)   &
                                                    + (FWUU(L,K-1,IT)-FWUU(L,K,IT)) )  
        ENDDO  
      ENDDO  

      IF( ISFCT(MVAR) >= 1 .AND. ISADAC(MVAR) > 0 )THEN 
        DO K=1,KC  
          DO LP=1,LLWET(K,0)
            L=LKWET(LP,K,0)  
            CON2(L,K,IT) = MAX(CON1(L,K),0.0)
          ENDDO  
        ENDDO  
      ENDIF  
    
    ELSE  ! *** IF ISTL = 3
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          CD(L,K,IT) = CON1(L,K)*H2PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                       FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT) ) * DXYIP(L)   &
                                                    + (FWUU(L,K-1,IT)-FWUU(L,K,IT))  )  
        ENDDO  
      ENDDO
      
      IF( ISFCT(MVAR) >= 1 .AND. ISADAC(MVAR) > 0 )THEN
        DO K=1,KC  
          DO LP=1,LLWET(K,0)
            L=LKWET(LP,K,0)  
            CON2(L,K,IT) = MAX(CON(L,K),0.0)
          ENDDO  
        ENDDO  
      ENDIF  
    ENDIF     ! *** ENDIF ON TIME LEVEL CHOICE FOR ISCDCA=0  
  
    IF( IS2TL_ == 0 .AND. ISUD == 1 )THEN  
      ! *** ADVANCE CON1 TO CON (3TL)
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          CON1(L,K) = CON(L,K)  
        ENDDO  
      ENDDO  
    ENDIF  
      
    ! *** UPDATE NEW CONCENTRATIONS  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)
        CON(L,K) = CD(L,K,IT)*HPKI(L,K)
      ENDDO
    ENDDO  
    
  ELSE  ! *** ELSE ON TRANSPORT OPTION CHOICE: ISCDCA(MVAR)/=0
    
    ! *** Central Differencing (3TL)
    ! *** Experimental Upwind Differencing (3TL)  
    IF( ISTL_ == 2 )THEN  
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          CD(L,K,IT) = CON1(L,K)*H1PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                    &
                                                     FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT) ) * DXYIP(L)   &
                                                  + (FWUU(L,K-1,IT)-FWUU(L,K,IT)) )  
        ENDDO  
      ENDDO
      
    ELSE   ! *** ISTL == 3  
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0) 
          CD(L,K,IT) = CON1(L,K)*H2PK(L,K) + DDELT*( ( FQC(L,K,IT) +                                                                      &
                                                       FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT) ) * DXYIP(L)   &
                                                    + (FWUU(L,K-1,IT)-FWUU(L,K,IT)) )  
        ENDDO  
      ENDDO  
    ENDIF   ! *** ENDIF ON TIME LEVEL CHOICE FOR ISCDCA/=0  
  
    ! *** SAVE CONCENTRATION FOR FLUX CORRECTOR
    IF( ISFCT(MVAR) >= 1 )THEN  
      DO K=1,KC
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          CON2(L,K,IT) = MAX(CON1(L,K),0.0)
        ENDDO
      ENDDO
    ENDIF  

    IF( ISUD == 1 )THEN  
      ! *** ADVANCE CON1 TO CON (3TL SOLUTION)
      DO K=1,KC  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          CON1(L,K) = CON(L,K)  
        ENDDO  
      ENDDO  
    ENDIF  
      
    ! *** COMPUTE THE CURRENT CONCENTRATIONS
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        CON(L,K) = CD(L,K,IT)*HPKI(L,K)
      ENDDO  
    ENDDO  
    
  ENDIF ! *** END OF TRANSPORT OPTION CHOICE  

  ! *** RESTORE ORIGINAL CONCENTRATIONS PRIOR TO APPLYING OPEN BC'S - 2TL & 3TL
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC)  
    DO K=1,KC  
      CON1(L,K) = WQBCCON1(IOBC,K,IT)
    ENDDO  
  ENDDO  

  ! *** ALL OTHER WATER CONSTITUENTS  
  IF( MVAR == 8 )THEN  ! .AND. IWQPSL == 2 )THEN
    M=4+NTOX+NSED+NSND+MO  
  ENDIF
  
  ! ******************************************************************************************
  ! *** APPLY OPEN BOUNDARY CONDITIONS, BASED ON DIRECTION OF FLOW  

  ! *** SOUTH OPEN BC, WITHOUT FLOCS
  IF( NCBS > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBS  
        NSID=NCSERS(LL,M)  
        L=LCBS(LL) 
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        LN=LNC(L)  
        IF( VHDX2(LN,K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
                                    CTMP = CON1(L,K) + DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2 ) CTMP = CON1(L,K) + DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPKI(L,K)  
            IF( ISCDCA(MVAR) == 2 ) CTMP = 0.5*(CON1(L,K)+CON(L,K)) + 0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)  &  
                                           + DDELT*(0.5*VHDX2(LN,K)*(CON1(L,K)+CON(L,K))-FVHUD(LN,K,IT))*DXYIP(L)*HPKI(L,K)
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          IF( M == 1 )THEN
            ! *** LIMIT CONCENTRATIONS TO MAXIMUM BC CONCENTRATIONS AT BOTTOM LAYER (SALINITY ONLY)
            CBSTMP = CBS(LL,1,M) + CSERT(1,NSID,M)  
            IF( CON(L,K) > CBSTMP )THEN
              CON(L,K) = CBSTMP  
            ENDIF
          ENDIF  
          CLOS(LL,K,M)=CON(L,K)  
          NLOS(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT=WTCI(K,1)*CBS(LL,1,M)+WTCI(K,2)*CBS(LL,2,M)+CSERT(K,NSID,M)  
          NMNLOD=NITER-NLOS(LL,K,M)  
          IF( NMNLOD >= NTSCRS(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBSTMP = CLOS(LL,K,M) + (CBT-CLOS(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRS(LL))  
            CON(L,K) = MAX(CBSTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** WEST OPEN BC, WITHOUT FLOCS  
  IF( NCBW > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBW  
        NSID=NCSERW(LL,M)  
        L=LCBW(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        IF( UHDY2(LEC(L),K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR)/=2) CTMP=CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPKI(L,K)
            IF( ISCDCA(MVAR) == 2 ) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)         &  
                               +DDELT*(0.5*UHDY2(LEC(L),K)*(CON1(L,K)+CON(L,K))-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPKI(L,K) 
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBWTMP=CBW(LL,1,M)+CSERT(1,NSID,M)  
          IF( M == 1 .AND. CON(L,K) > CBWTMP) CON(L,K) = CBWTMP  
          CLOW(LL,K,M)=CON(L,K)  
          NLOW(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT=WTCI(K,1)*CBW(LL,1,M)+WTCI(K,2)*CBW(LL,2,M)+CSERT(K,NSID,M)  
          NMNLOD=NITER-NLOW(LL,K,M)  
          IF( NMNLOD >= NTSCRW(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBWTMP=CLOW(LL,K,M)+(CBT-CLOW(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRW(LL))  
            CON(L,K) = MAX(CBWTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** EAST OPEN BC, WITHOUT FLOCS  
  IF( NCBE > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBE  
        NSID=NCSERE(LL,M)  
        L=LCBE(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) ) CYCLE 
        IF( UHDY2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR) /= 2) CTMP=CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)  
            IF( ISCDCA(MVAR) == 2 ) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2P(L)*HPKI(L,K)+DDELT*(FUHUD(L,K,IT) &  
                                      -0.5*UHDY2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPKI(L,K)
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBETMP = CBE(LL,1,M)+CSERT(1,NSID,M)  
          IF( M == 1 .AND. CON(L,K) > CBETMP) CON(L,K) = CBETMP  
          CLOE(LL,K,M)=CON(L,K)  
          NLOE(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT = WTCI(K,1)*CBE(LL,1,M) + WTCI(K,2)*CBE(LL,2,M) + CSERT(K,NSID,M)  
          NMNLOD = NITER-NLOE(LL,K,M)  
          IF( NMNLOD >= NTSCRE(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBETMP = CLOE(LL,K,M) + (CBT-CLOE(LL,K,M)) * FLOAT(NMNLOD)/FLOAT(NTSCRE(LL))  
            CON(L,K) = MAX(CBETMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF
    
  ! *** NORTH OPEN BC, WITHOUT FLOCS  
  IF( NCBN > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBN  
        NSID=NCSERN(LL,M)  
        L=LCBN(LL)  
        IF( LKSZ(L,K) .OR. .NOT. LMASKDRY(L) )CYCLE 
        LS=LSC(L)  
        IF( VHDX2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K)
          ELSE  
            IF( ISCDCA(MVAR)/=2) CTMP=CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPKI(L,K) 
            IF( ISCDCA(MVAR) == 2 ) CTMP=0.5*(CON1(L,K)+CON(L,K))+0.5*(CON1(L,K)-CON(L,K))*H2PK(L,K)*HPKI(L,K) + DDELT*(FVHUD(L,K,IT) &  
                                      -0.5*VHDX2(L,K)*(CON1(L,K)+CON(L,K)))*DXYIP(L)*HPKI(L,K) 
          ENDIF  
          CON(L,K) = MAX(CTMP  ,0.)
          CBNTMP=CBN(LL,1,M)+CSERT(1,NSID,M)  
          IF( M == 1 .AND. CON(L,K) > CBNTMP) CON(L,K) = CBNTMP  
          CLON(LL,K,M)=CON(L,K)  
          NLON(LL,K,M)=NITER  
        ELSE  
          ! *** FLOWING INTO DOMAIN  
          CBT=WTCI(K,1)*CBN(LL,1,M)+WTCI(K,2)*CBN(LL,2,M)+CSERT(K,NSID,M)  
          NMNLOD=NITER-NLON(LL,K,M)  
          IF( NMNLOD >= NTSCRN(LL) )THEN  
            CON(L,K) = CBT  
          ELSE  
            CBNTMP=CLON(LL,K,M)+(CBT-CLON(LL,K,M))*FLOAT(NMNLOD)/FLOAT(NTSCRN(LL))  
            CON(L,K) = MAX(CBNTMP,0.)
          ENDIF  
        ENDIF  
        IF( ISUD == 1 ) CON1(L,K) = CON(L,K)
      ENDDO  
    ENDDO  
  ENDIF    

  ! ****************************************************************************************
  ! **  ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATIONS WITH FLUX CORRECTOR  
  IF( ISADAC(MVAR) == 0 ) GOTO 2000  ! *** SKIP IF NOT USING ANTI-DIFFUSION
  IF( ISCDCA(MVAR) == 1 ) GOTO 2000  ! *** SKIP IF CENTRAL DIFFERENCE (3TL)
  
  ! ----------------------------------------------------------------------------------------
  ! **  STANDARD ANTI-DIFFUSIVE ADVECTIVE FLUX CALCULATION  
  
  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  DO IOBC=1,NBCSOP  
    L=LOBCS(IOBC) 
    DO K=1,KC  
      WQBCCON(IOBC,K,IT)  = CON(L,K)  
    ENDDO  
  ENDDO  

  ! *** GET ONLY POSITIVE CONCENTRATIONS FOR ANTI-DIFUSION CALCULATIONS
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      POS(L,K,IT) = MAX(CON(L,K),0.) 
    ENDDO  
  ENDDO  
   
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      UUUU(L,K,IT) = U2(L,K)*( POS(L,K,IT) - POS(LWC(L),K,IT) )*DXIU(L)  
      VVVV(L,K,IT) = V2(L,K)*( POS(L,K,IT) - POS(LSC(L),K,IT) )*DYIV(L)  
    ENDDO  
  ENDDO

  IF( KC > 1 )THEN  
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        WWWW(L,K,IT) = W2(L,K)*( POS(L,K+1,IT) - POS(L,K,IT) )*HPI(L)*DZIG(L,K)
      ENDDO  
    ENDDO  
    
    ! *** ASSIGN FOR ZERO BOTTOM GRADIENT FOR VARIABLE SGZ LAYERS
    IF( IGRIDV > 0 .AND. .FALSE. )THEN
      DO LP=1,LAWET
        L = LWET(LP)
        IF( KSZ(L) > 1 )WWWW(L,KSZ(L)-1,IT) = -WWWW(L,KSZ(L),IT)
      ENDDO
    ENDIF
  ENDIF

  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      LE=LEC(L)
      LW=LWC(L)
      LN=LNC(L)  
      LS=LSC(L)  
      LNW=LNWC(L)  
      LSE=LSEC(L)
      
      ! *** U COMPONENTS
      AUHU=ABS(UHDY2(L,K))  
      UTERM=AUHU*( POS(L,K,IT) - POS(LW,K,IT) )  
      IF( UHDY2(L,K) >= 0.0 )THEN  
        ! *** QUANTITIES FROM WEST CELL
        UTERM = UTERM - 0.5*DDELTA*UHDY2(L,K)*( VVVV(LNW,K,IT)+VVVV(LW,K,IT) + WWWW(LW,K,IT)+WWWW(LW,K-1,IT) + UUUU(L,K,IT)+UUUU(LW,K,IT) )  
      ELSE  
        ! *** QUANTITIES FROM CURRENT CELL
        UTERM = UTERM - 0.5*DDELTA*UHDY2(L,K)*( VVVV(LN,K,IT) +VVVV(L,K,IT)  + WWWW(L,K,IT) +WWWW(L,K-1,IT)  + UUUU(L,K,IT)+UUUU(LE,K,IT) )  
      ENDIF  
      UHU = UTERM/( POS(L,K,IT) + POS(LW,K,IT) + BSMALL )  
      FUHUD(L,K,IT) = MAX(UHU,0.)*POS(LW,K,IT) + MIN(UHU,0.)*POS(L,K,IT)

      ! *** V COMPONENTS
      AVHV = ABS(VHDX2(L,K))  
      VTERM = AVHV*( POS(L,K,IT) - POS(LS,K, IT) )  
      IF( VHDX2(L,K) >= 0.0 )THEN
        ! *** QUANTITIES FROM SOUTH CELL
        VTERM = VTERM - 0.5*DDELTA*VHDX2(L,K)*( UUUU(LS,K,IT)+UUUU(LSE,K,IT) + WWWW(LS,K,IT)+WWWW(LS,K-1,IT) + VVVV(LS,K,IT)+VVVV(L,K,IT) )  
      ELSE  
        ! *** QUANTITIES FROM CURRENT CELL
        VTERM = VTERM - 0.5*DDELTA*VHDX2(L,K)*( UUUU(L,K,IT) +UUUU(LE,K,IT)  + WWWW(L,K,IT) +WWWW(L ,K-1,IT) + VVVV(LN,K,IT)+VVVV(L,K,IT) )  
      ENDIF  
      VHV = VTERM/( POS(L,K,IT) + POS(LS,K ,IT) + BSMALL )  
      FVHUD(L,K,IT) = MAX(VHV,0.)*POS(LS,K ,IT) + MIN(VHV,0.)*POS(L,K,IT)  
    ENDDO  
  ENDDO  

  IF( KC > 1 )THEN
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)
        LE=LEC(L)
        LN=LNC(L)  
        AWW=ABS(W2(L,K))  
        WTERM=AWW*(POS(L,K+1,IT)-POS(L,K,IT))  
        IF( W2(L,K) >= 0.0 )THEN  
          WTERM=WTERM-0.5*DDELTA*W2(L,K)*(UUUU(L,K,IT)  +UUUU(LE,K,IT)   + VVVV(L,K,IT)  +VVVV(LN,K,IT)   + WWWW(L,K,IT)+WWWW(L,K-1,IT))
        ELSE  
          WTERM=WTERM-0.5*DDELTA*W2(L,K)*(UUUU(L,K+1,IT)+UUUU(LE,K+1,IT) + VVVV(L,K+1,IT)+VVVV(LN,K+1,IT) + WWWW(L,K,IT)+WWWW(L,K+1,IT))
        ENDIF  
        WW = WTERM/( POS(L,K+1,IT)+POS(L,K,IT) + BSMALL )  
        FWUU(L,K,IT) = MAX(WW,0.)*POS(L,K,IT)+MIN(WW,0.)*POS(L,K+1,IT)  
      ENDDO  
    ENDDO  
  ENDIF

  ! ----------------------------------------------------------------------------------------
  ! ** ZERO ANTIDIFFUSIVE FLUXES FOR SELECTED CELLS
  
  ! ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR SOURCE CELLS  
  IF( ISADAC(MVAR) == 1 )THEN  
    ! *** ANTI-DIFFUSION TURNED OFF FOR SOURCE CELLS
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        IF( ABS(QSUM(L,K)) > 1.E-8 )THEN
          LN=LNC(L)  
          FUHUD(L  ,K,IT) = 0.  
          FUHUD(LEC(L),K,IT) = 0.  
          FVHUD(L  ,K,IT) = 0.  
          FVHUD(LN ,K,IT) = 0.  
          FWUU(L,K  ,IT)  = 0.  
          FWVV(L,K-1,IT)  = 0.  
        ENDIF  
      ENDDO  
    ENDDO  
  ENDIF  
  
  ! ** SET ANTIDIFFUSIVE FLUXES TO ZERO FOR OPEN BOUNDARY CELLS  
  DO K=1,KC  
    DO LL=1,NCBS  
      L=LCBS(LL)  
      LN=LNC(L)  
      FVHUD(LN,K,IT) = 0.0  
    ENDDO  
    DO LL=1,NCBW  
      L=LCBW(LL)  
      FUHUD(LEC(L),K,IT) = 0.0  
    ENDDO  
    DO LL=1,NCBE  
      L=LCBE(LL)  
      FUHUD(L,K,IT) = 0.0  
    ENDDO  
    DO LL=1,NCBN  
      L=LCBN(LL)  
      FVHUD(L,K,IT) = 0.0  
    ENDDO  
  ENDDO  
  
  ! ----------------------------------------------------------------------------------------
  ! **  CALCULATE AND APPLY FLUX CORRECTED TRANSPORT LIMITERS  
  IF( ISFCT(MVAR) /= 0 )THEN  
  
    ! **  DETERMINE MAX AND MIN CONCENTRATIONS  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        CONTMX(L,K,IT) = MAX(POS(L,K,IT),ABS(CON2(L,K,IT)))  
        CONTMN(L,K,IT) = MIN(POS(L,K,IT),ABS(CON2(L,K,IT)))  
      ENDDO  
    ENDDO
    
    DO LP=1,LAWET
      L=LWET(LP) 
      CMAX(L,KSZ(L),IT) = MAX(CONTMX(L,KSZ(L),IT), CONTMX(L,KSZ(L)+1,IT))  
      CMAX(L,KC    ,IT) = MAX(CONTMX(L,KS    ,IT), CONTMX(L,KC      ,IT))  
      CMIN(L,KSZ(L),IT) = MIN(CONTMN(L,KSZ(L),IT), CONTMN(L,KSZ(L)+1,IT))  
      CMIN(L,KC    ,IT) = MIN(CONTMN(L,KS    ,IT), CONTMN(L,KC      ,IT))  
    ENDDO  
    
    DO K=2,KS  
      DO LP=1,LLWET(K-1,0)
        L=LKWET(LP,K-1,0)  
        CMAXT        = MAX(CONTMX(L,K-1,IT),CONTMX(L,K+1,IT))
        CMAX(L,K,IT) = MAX(CONTMX(L,K,IT),CMAXT)  
        CMINT        = MIN(CONTMN(L,K-1,IT),CONTMN(L,K+1,IT))  
        CMIN(L,K,IT) = MIN(CONTMN(L,K,IT),CMINT)  
      ENDDO  
    ENDDO  
    
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        LE=LEC(L)
        LS=LSC(L)  
        LN=LNC(L)  
        LW=LWC(L)
        
        CWMAX = SUB3D(L,K) *CONTMX(LW,K,IT)  
        CEMAX = SUB3D(LE,K)*CONTMX(LE,K,IT)  
        CSMAX = SVB3D(L,K) *CONTMX(LS,K,IT)  
        CNMAX = SVB3D(LN,K)*CONTMX(LN,K,IT)  
        CMAXT = MAX(CNMAX,CEMAX)  
        CMAXT = MAX(CMAXT,CSMAX)  
        CMAXT = MAX(CMAXT,CWMAX)  
        CMAX(L,K,IT) = MAX(CMAX(L,K,IT),CMAXT)  
        
        CWMIN = SUB3D(L,K) *CONTMN(LW,K,IT) + 1.E+6*(1.-SUB3D(L,K))
        CEMIN = SUB3D(LE,K)*CONTMN(LE,K,IT) + 1.E+6*(1.-SUB3D(LE,K))  
        CSMIN = SVB3D(L,K) *CONTMN(LS,K,IT) + 1.E+6*(1.-SVB3D(L,K))
        CNMIN = SVB3D(LN,K)*CONTMN(LN,K,IT) + 1.E+6*(1.-SVB3D(LN,K))  
        CMINT=MIN(CNMIN,CEMIN)  
        CMINT=MIN(CMINT,CSMIN)  
        CMINT=MIN(CMINT,CWMIN)  
        CMIN(L,K,IT) = MIN(CMIN(L,K,IT),CMINT)  
      ENDDO  
    ENDDO  
  
    ! **  SEPARATE POSITIVE AND NEGATIVE FLUXES  
    ! **  PUT NEGATIVE FLUXES INTO FUHVD, FVHVD, AND FWVV  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        FUHVD(L,K,IT) = MIN(FUHUD(L,K,IT),0.)
        FUHUD(L,K,IT) = MAX(FUHUD(L,K,IT),0.)
        FVHVD(L,K,IT) = MIN(FVHUD(L,K,IT),0.)
        FVHUD(L,K,IT) = MAX(FVHUD(L,K,IT),0.)  
      ENDDO  
    ENDDO

    IF( KC > 1 )THEN  
      DO K=1,KS  
        DO LP=1,LLWET(K,0)
          L=LKWET(LP,K,0)  
          FWVV(L,K,IT) = MIN(FWUU(L,K,IT),0.)  
          FWUU(L,K,IT) = MAX(FWUU(L,K,IT),0.)  
        ENDDO  
      ENDDO  
    ENDIF

    ! **  CALCULATE INFLUX AND OUTFLUX IN CONCENTRATION UNITS  
    ! **  LOAD INTO DUU AND DVV, THEN ZERO VALUES AT BOUNDARIES  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        LE=LEC(L)
        LN=LNC(L)
        ! ***                              MAX            MIN              MAX            MIN                 MAX            MIN
        DUU(L,K,IT) = DDELT*( DXYIP(L)*( FUHUD(L,K,IT) -FUHVD(LE,K,IT) + FVHUD(L,K,IT) -FVHVD(LN,K,IT) ) + ( FWUU(L,K-1,IT)-FWVV(L,K,IT))   )*HPKI(L,K)
        DVV(L,K,IT) = DDELT*( DXYIP(L)*( FUHUD(LE,K,IT)-FUHVD(L,K,IT)  + FVHUD(LN,K,IT)-FVHVD(L,K,IT)  ) + ( FWUU(L,K,IT)  -FWVV(L,K-1,IT)) )*HPKI(L,K)
      ENDDO  
    ENDDO  

    ! *** ZERO SELECTED FLUX CORRECTORS
    DO K=1,KC  
      DO IOBC=1,NBCSOP  
        L=LOBCS(IOBC)  
        DUU(L,K,IT) = 0.  
        DVV(L,K,IT) = 0.  
      ENDDO  
    END DO  
  
    ! **  CALCULATE BETA COEFFICIENTS WITH BETAUP AND BETADOWN IN DUU AND DVV  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        IF( DUU(L,K,IT) > 0. ) DUU(L,K,IT) = (CMAX(L,K,IT)-POS(L,K,IT))/(DUU(L,K,IT)+BSMALL)  
        DUU(L,K,IT) = MIN(DUU(L,K,IT),1.)  
        IF( DVV(L,K,IT) > 0. ) DVV(L,K,IT) = (POS(L,K,IT)-CMIN(L,K,IT))/(DVV(L,K,IT)+BSMALL)  
        DVV(L,K,IT) = MIN(DVV(L,K,IT),1.)  
      ENDDO  
    ENDDO  
  
    ! **  LIMIT FLUXES  
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        LS=LSC(L)
        LW=LWC(L)
        FUHUD(L,K,IT) = MIN(DVV(LW,K,IT),DUU(L,K,IT))*FUHUD(L,K,IT) + MIN(DUU(LW,K,IT),DVV(L,K,IT))*FUHVD(L,K,IT)  
        FVHUD(L,K,IT) = MIN(DVV(LS,K,IT),DUU(L,K,IT))*FVHUD(L,K,IT) + MIN(DUU(LS,K,IT),DVV(L,K,IT))*FVHVD(L,K,IT)  
      ENDDO  
    ENDDO  

    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        FWUU(L,K,IT) = MIN(DVV(L,K,IT),DUU(L,K+1,IT))*FWUU(L,K,IT) + MIN(DUU(L,K,IT),DVV(L,K+1,IT))*FWVV(L,K,IT)  
      ENDDO  
    ENDDO  
    
    ! *** APPLY OPEN BOUNDARYS 
    DO LL=1,NBCSOP
      L=LOBCS(LL)
      DO K=1,KS  
        FWUU(L,K,IT) = 0.0
      ENDDO  
    ENDDO 
    DO LL=1,NBCSOP2
      L=LOBCS2(LL)
      DO K=1,KS  
        FWUU(L,K,IT) = 0.0
      ENDDO  
    ENDDO 

  ENDIF  ! *** END OF FLUX CORRECTOR SECTION
  
  ! *** APPLY THE ANTI-DIFFUSIVE ADVECTION CALCULATION
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      CD(L,K,IT) = CON(L,K)*HPK(L,K) + DDELT*( ( FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT) ) * DXYIP(L)   &
                                              + ( FWUU(L,K-1,IT)-FWUU(L,K,IT) )  )  
      CON(L,K) = CD(L,K,IT)*HPKI(L,K)
    ENDDO  
  ENDDO  
  
  ! *** RESET OPEN BC CONCENTRATIONS  
  DO K = 1,KC
    DO IOBC=1,NBCSOP
      L=LOBCS(IOBC)
      CON(L,K) = WQBCCON(IOBC,K,IT)
    ENDDO
  ENDDO
  
  ! ----------------------------------------------------------------------------------------
  ! *** CALTRAN EXIT 
 2000 CONTINUE  
 
  ! *** ZERO HEAT FLUXES
  IF( MVAR == 2 )THEN        
    ! *** ZERO EVAP/RAINFALL
    DO L=1,LC  
      FQC(L,KC,IT) = 0.  
    ENDDO  
    IF( ISADAC(MVAR) >= 2 )THEN
      DO L=1,LC  
        FQCPAD(L,KC,IT) = 0.  
      ENDDO  
    ENDIF
    IF( ISADAC(MVAR) > 0 )THEN
      DO L=1,LC  
        QSUMPAD(L,KC,IT) = 0.  
      ENDDO  
    ENDIF
  ENDIF

  !  $ print *, ' CALTRAN-Leaving: Thread, MVAR, MO',IT,MVAR,MO  

  RETURN  
END  
