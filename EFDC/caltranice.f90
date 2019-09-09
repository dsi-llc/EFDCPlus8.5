SUBROUTINE CALTRANICE(ISTL_,CON,CON1,IT)  

  ! **  SUBROUTINE CALTRAN CALCULATES THE ADVECTIVE  
  ! **  TRANSPORT OF FRAZIL ICE OR ANY OTHER ICE IMPACTED TRANSPORT
  ! **  A NEW VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL INDICATES  
  ! **  THE NUMBER OF TIME LEVELS IN THE STEP  

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------------------------------!
  ! ** 2015-01     PAUL M. CRAIG      Added fully coupled Ice Sub-model with Frazil Ice Transport
  ! ** 2014-09     PAUL M. CRAIG      ADDED THE LWET BYPASS APPROACH

  USE GLOBAL  
  !  USE OMP_LIB
  
  IMPLICIT NONE  
    
  INTEGER, INTENT(IN) :: ISTL_,IT  
  REAL    :: CON(LCM,KCM),CON1(LCM,KCM)  
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: FLUX

  REAL    :: DDELT  
  REAL    :: RDZIC, CTMP, WVEL, CLEFT, CRIGHT 
  
  INTEGER :: LP,L, LN, LS, LL, K
  
  IF(  .NOT. ALLOCATED(FLUX) )THEN
    ALLOCATE(FLUX(LCM,KCM))
    FLUX = 0.0
  ENDIF

  !ISUD=1  
  IF( ISDYNSTP == 0 )THEN
    ! *** FIXED DELTA T
    IF( ISTL_ == 3 )THEN  
      DDELT=DT2
    ELSE
      DDELT=DT
    ENDIF
  ELSE  
    DDELT=DTDYN   ! *** DYNAMIC DELTA T
  END IF  
    
 ! *** ADVANCE THE CONSTITUENT
  CON1=CON 
    
  ! *** SKIP BOUNDARY LOADINGS FOR FRAZIL TRANSPORT (MAY BE ADDED FOR FUTURE VERSIONS)
  !IF( .FALSE. )THEN
  !  ! *** SAVE OLD WQ CONCENTRATIONS FOR OPEN BOUNDARY CELLS
  !  DO IOBC=1,NBCSOP  
  !    L=LOBCS(IOBC)  
  !    DO K=1,KC  
  !      WQBCCON1(IOBC,K,IT) = CON1(L,K)  
  !      WQBCCON(IOBC,K,IT)  = CON(L,K)  
  !    ENDDO  
  !  ENDDO  
  !
  !  ! **  CALCULATED EXTERNAL SOURCES AND SINKS  
  !  CALL CALFQC (ISTL_,IS2TL_,MVAR,10,CON,CON1,IT)  
  !ENDIF

  ! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH ADVECTION  
  ! **  AVERAGED BETWEEN (N) AND (N+1) OR (N-1) AND (N+1) AND ADVECTED  
  ! **  AT (N) OR (N-1) IF ISTL EQUALS 2 OR 3 RESPECTIVELY  

  ! *** COMPUTE FLUXES
  DO K=1,KC
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      FUHUD(L,K,IT)=UHDY2(L,K)*CON1(LUPU(L,K),K)
      FVHUD(L,K,IT)=VHDX2(L,K)*CON1(LUPV(L,K),K)  
    ENDDO  
  ENDDO
  
  IF( KC > 1 )THEN  
    DO K=1,KS  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        FWUU(L,K,IT)=W2(L,K)*CON1(L,KUPW(L,K))  
      ENDDO  
    ENDDO  
  ENDIF  
  
  ! *** CALCULATE AND ADD HORIZONTAL DIFFUSION FLUX (PMC MOVED)  
  !IF(ISHDMF == 2 ) CALL CALDIFF (CON1,IT)  
    
  ! *** Upwind Differencing (3TL & 2TL)  
  IF( ISTL_ == 2 )THEN  
    DO K=1,KC  
      RDZIC=DZIC(L,K)  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        CD(L,K,IT)=CON1(L,K)*H1P(L)+DDELT*( ( FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT))*DXYIP(L) &
                                          + ( FWUU(L,K-1,IT)-FWUU(L,K,IT))*RDZIC )  
      ENDDO  
    ENDDO    
  ELSE  ! *** IF ISTL = 3
    DO K=1,KC  
      RDZIC=DZIC(L,K)  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        CD(L,K,IT)=CON1(L,K)*H2P(L)+DDELT*( ( FUHUD(L,K,IT)-FUHUD(LEC(L),K,IT) + FVHUD(L,K,IT)-FVHUD(LNC(L),K,IT))*DXYIP(L) &
                                          + ( FWUU(L,K-1,IT)-FWUU(L,K,IT))*RDZIC )  
      ENDDO  
    ENDDO
  ENDIF     ! *** ENDIF ON TIME LEVEL CHOICE FOR ISCDCA=0  
  
  IF( ISTL_ == 3 .AND. .FALSE. )THEN  
    ! *** ADVANCE CON1 TO CON
    DO K=1,KC  
      DO LP=1,LLWET(K,0)
        L=LKWET(LP,K,0)  
        CON1(L,K)=CON(L,K)  
      ENDDO
      
  !    ! *** RESET OPEN BC CONCENTRATIONS  
  !    IF( ISUD == 1 )THEN
  !      DO IOBC=1,NBCSOP  
  !        L=LOBCS(IOBC)  
  !        CON1(L,K)=WQBCCON1(IOBC,K,IT)  
  !      ENDDO  
  !    ENDIF
    ENDDO  
  ENDIF  
      
  ! *** UPDATE NEW CONCENTRATIONS  
  DO K=1,KC  
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      CON(L,K)=CD(L,K,IT)*HPI(L) 
    ENDDO

    ! *** RESET OPEN BC CONCENTRATIONS  
    !DO IOBC=1,NBCSOP  
    !  L=LOBCS(IOBC)  
    !  CON(L,K)=WQBCCON(IOBC,K,IT)  
    !ENDDO  
  ENDDO  

  !IF( ISUD == 1 )THEN
  !  ! *** RESTORE ORIGINAL CONCENTRATIONS PRIOR TO APPLYING OPEN BC'S  
  !  DO K=1,KC  
  !    DO IOBC=1,NBCSOP  
  !      L=LOBCS(IOBC)  
  !      CON1(L,K)=WQBCCON(IOBC,K,IT)  
  !    ENDDO  
  !  ENDDO  
  !ENDIF  

  ! ******************************************************************************************
  ! *** COMPUTE THE FRAZIL ICE RISE
  
  ! *** BOTTOM LAYER
  DO LP=1,LAWET
    K = KSZ(L)
    L=LWET(LP) 
    WVEL   = DDELT*HPI(L)*DZIC(L,K)
    CLEFT  = 1.+RISEVEL*WVEL
    CRIGHT = CON(L,K)
    CON(L,K) = CRIGHT/CLEFT
    FLUX(L,K) = RISEVEL*CON(L,K) 
  ENDDO

  ! *** ALL OTHER LAYERS
  DO K=2,KC
    DO LP=1,LLWET(K,0)
      L=LKWET(LP,K,0)  
      WVEL   = DDELT*HPI(L)*DZIC(L,K)
      CLEFT  = 1.+RISEVEL*WVEL
      CRIGHT = CON(L,K)+FLUX(L,K-1)*WVEL
      CON(L,K)=CRIGHT/CLEFT
      FLUX(L,K) = RISEVEL*CON(L,K)
    ENDDO
  ENDDO

  ! *** ACCUMULATE ICE COVER
  DO LP=1,LAWET
    L=LWET(LP)
    ICETHICK(L) = ICETHICK(L) + FLUX(L,KC)*DDELT/(HP(L)*DZC(L,KC))
  ENDDO
  
  ! ******************************************************************************************
  ! *** APPLY OPEN BOUNDARY CONDITIONS, BASED ON DIRECTION OF FLOW  

  ! *** SOUTH OPEN BC
  IF( NCBS > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBS  
        L=LCBS(LL)  
        LN=LNC(L)  
        IF( VHDX2(LN,K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPI(L)  
          ELSE  
            CTMP=CON1(L,K)+DDELT*(VHDX2(LN,K)*CON1(L,K)-FVHUD(LN,K,IT))*DXYIP(L)*HPI(L)  
            CON1(L,K)=CON(L,K)  
          ENDIF  
          CON(L,K)=MAX(CTMP  ,0.)
        ENDIF  
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** WEST OPEN BC 
  IF( NCBW > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBW  
        L=LCBW(LL)  
        IF( UHDY2(LEC(L),K) <= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPI(L)  
          ELSE  
            CTMP=CON1(L,K)+DDELT*(UHDY2(LEC(L),K)*CON1(L,K)-FUHUD(LEC(L),K,IT))*DXYIP(L)*HPI(L)  
            CON1(L,K)=CON(L,K)  
          ENDIF  
          CON(L,K)=MAX(CTMP  ,0.)
        ENDIF  
      ENDDO  
    ENDDO  
  ENDIF    

  ! *** EAST OPEN BC
  IF( NCBE > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBE  
        L=LCBE(LL)  
        IF( UHDY2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
          ELSE  
            CTMP=CON1(L,K)+DDELT*(FUHUD(L,K,IT)-UHDY2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
            CON1(L,K)=CON(L,K)  
          ENDIF  
          CON(L,K)=MAX(CTMP  ,0.)
        ENDIF  
      ENDDO  
    ENDDO  
  ENDIF
    
  ! *** NORTH OPEN BC 
  IF( NCBN > 0 )THEN
    DO K=1,KC  
      DO LL=1,NCBN  
        L=LCBN(LL)  
        LS=LSC(L)  
        IF( VHDX2(L,K) >= 0. )THEN  
          ! *** FLOWING OUT OF DOMAIN  
          IF( ISTL_ == 2 )THEN  
            CTMP=CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
          ELSE  
            CTMP=CON1(L,K)+DDELT*(FVHUD(L,K,IT)-VHDX2(L,K)*CON1(L,K))*DXYIP(L)*HPI(L)  
            CON1(L,K)=CON(L,K)  
          ENDIF  
          CON(L,K)=MAX(CTMP  ,0.)
        ENDIF  
      ENDDO  
    ENDDO  
  ENDIF

  ! ----------------------------------------------------------------------------------------
  ! *** CALTRANICE EXIT 

  !  $ print *, ' CALTRAN-Leaving: Thread, MVAR, MO',IT,MVAR,MO  

  RETURN  
END  
