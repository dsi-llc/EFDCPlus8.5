SUBROUTINE CALCONC (ISTL_,IS2TL_)  

  ! *** SUBROUTINE CALCULATES THE CONCENTRATION OF DISSOLVED AND  
  ! *** SUSPENDED CONSTITUTENTS, INCLUDING SALINITY, TEMPERATURE, DYE AND  
  ! *** AND SUSPENDED SEDIMENT AT TIME LEVEL (N+1). THE VALUE OF ISTL  
  ! *** INDICATES THE NUMBER OF TIME LEVELS IN THE STEP  
  !  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2015-01       PAUL M. CRAIG     Added fully coupled Ice Sub-model with Frazil Ice Transport
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2014-08       D H CHUNG         SET EXPLICIT PRECISIONS OF INTEGER & REAL
  !    2011-03       PAUL M. CRAIG     MERGED LATEST CODES AND RESTRUCTURED, ADDED OMP
  !    2002-05       John Hamrick      Modified calls to calbal and budget subroutines
  !                                     added calls to bal2t2, bal2t3
  !------------------------------------------------------------------------------------------------!

  USE GLOBAL  
  USE OMP_LIB
  USE HEAT_MODULE, ONLY:CALHEAT
  
  IMPLICIT NONE
  
  INTEGER :: ISTL_,IS2TL_
  INTEGER :: K,L,LP,NT,NS,ND,IT,M
  INTEGER :: NTMP,LF,LL,LE,LN
  INTEGER, SAVE :: NICE
  
  REAL      :: RCDZKMK, RCDZKK
  REAL      :: DELTD2,CDYETMP,DAGE
  REAL,SAVE :: SEDTIME
  
  REAL(RKD), EXTERNAL :: DSTIME 
  REAL(RKD)           :: CCUBTMP, CCMBTMP     ! VERTICAL DIFFUSION TEMPORARY VARIABLES
  REAL(RKD)           :: TTDS                 ! MODEL TIMING TEMPORARY VARIABLE

  ! *** VERTICAL DIFFUSION VARIABLES
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:)   :: CCLBTMP 
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:)   :: EEB     
  REAL(RKD),SAVE,ALLOCATABLE,DIMENSION(:,:) :: VCU         

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: TOXASM
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SEDASM
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SNDASM

  IF(  .NOT. ALLOCATED(EEB) )THEN
    ALLOCATE(CCLBTMP(LCM))
    ALLOCATE(EEB(LCM))
    ALLOCATE(VCU(LCM,KCM))
    ALLOCATE(TOXASM(NTXM))
    ALLOCATE(SEDASM(NSCM))
    ALLOCATE(SNDASM(NSNM))
    EEB=0.0 
    CCLBTMP=0.0 
    VCU=0.
    TOXASM=0.0 
    SEDASM=0.0 
    SNDASM=0.0 
    NICE=0
    SEDTIME=0.0
  ENDIF

  DELT=DT2  
  IF( ISTL_ == 2 )THEN  
    IF( ISDYNSTP == 0 )THEN  
      DELT=DT  
    ELSE  
      DELT=DTDYN  
    END IF  
  ENDIF  
  DELTD2=DELT  
  
  ! *** MASS BALANCE
  IF( IS2TIM >= 1 )THEN  
    IF( ISBAL >= 1 )THEN  
      CALL BAL2T3A  
    ENDIF  
  ENDIF  
  IT=1
  
  ! ***************************************************************************************
  ! *** 3D ADVECTI0N TRANSPORT CALCULATION, STANDARD TRANSPORT
  TTDS=DSTIME(0)  

  !$OMP PARALLEL DEFAULT(SHARED)

  ! *** PRESPECIFY THE UPWIND CELLS FOR 3D ADVECTION  
  !$OMP DO PRIVATE(ND,K,LP,L,LE,LN)
  DO ND=1,NDM  
    DO K=1,KC  
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        IF( UHDY2(L,K) >= 0.0 )THEN  
          LUPU(L,K)=LWC(L)  
        ELSE  
          LUPU(L,K)=L  
        END IF  
        IF( VHDX2(L,K) >= 0.0 )THEN  
          LUPV(L,K)=LSC(L)  
        ELSE  
          LUPV(L,K)=L  
        END IF  
      ENDDO  
    ENDDO  
    
    IF( KC > 1 )THEN  
      DO K=1,KS  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( W2(L,K) >= 0. )THEN
            KUPW(L,K)=K  
          ELSE  
            KUPW(L,K)=K+1
          END IF  
        ENDDO  
      ENDDO  
    ENDIF  
  ENDDO
  !$OMP END DO

  ! *** EE7.2 - MUST ZERO ALL THREADS FOR EVERY INSTANCE DO TO POTENTIAL OF DIFFERENT THREADS MAY NOT BE ZEROED FOR CURRENT TIME
  ! *** ZERO DRY CELL FLUXES 
  IF( LADRY > 0 )THEN
    !$OMP SINGLE
    DO ND=1,NDM  
      DO K=1,KC
        DO LP=1,LADRY
          L=LDRY(LP)  
          FUHUD(L,K,ND)=0.  
          FVHUD(L,K,ND)=0.  
          FUHVD(L,K,ND)=0.  
          FVHVD(L,K,ND)=0.  
          UUUU(L,K,ND) =0.
          VVVV(L,K,ND) =0.
          DUU(L,K,ND)  =0.
          DVV(L,K,ND)  =0. 
          POS(L,K,ND)  =0.
          WWWW(L,K,ND) =0.
          FWUU(L,K,ND) =0. 

          LN=LNC(L)
          LE=LEC(L)
          FUHUD(LE,K,ND)=0.  
          FVHUD(LN,K,ND)=0.  
          FUHVD(LE,K,ND)=0.  
          FVHVD(LN,K,ND)=0.  
          UUUU(LE,K,ND) =0.
          VVVV(LN,K,ND) =0.
          DUU(LE,K,ND)  =0.
          DVV(LN,K,ND)  =0. 
        ENDDO
      ENDDO
    ENDDO
    !$OMP END SINGLE
  ENDIF
  
  IF( ISTRAN(1) == 1 .AND. ISCDCA(1) < 4 )THEN
    !$OMP SINGLE PRIVATE(IT)
    !$  IT = OMP_GET_THREAD_NUM() + 1
    CALL CALTRAN(ISTL_,IS2TL_,1,1,SAL,SAL1,IT)
    !$OMP END SINGLE NOWAIT
  ENDIF

  IF( ISTRAN(2) == 1 .AND. ISCDCA(2) < 4 )THEN
    !$OMP SINGLE PRIVATE(IT)
    !$  IT = OMP_GET_THREAD_NUM() + 1
    CALL CALTRAN(ISTL_,IS2TL_,2,2,TEM,TEM1,IT)  
    !$OMP END SINGLE NOWAIT 

    IF( ISICE == 4 .AND. (IS2TL > 0 .OR. (IS2TL == 0 .AND. NCTBC /= NTSTBC)) )THEN
      IF( LFRAZIL )THEN
        ! *** FRAZIL ICE TRANSPORT
        !$OMP SINGLE PRIVATE(IT)
        !$  IT = OMP_GET_THREAD_NUM() + 1
        CALL CALTRANICE(ISTL_,FRAZILICE,FRAZILICE1,IT)
        NICE = NICE+1
        !$OMP END SINGLE NOWAIT

      ELSEIF( NICE > 0 )THEN
        ! *** ADVANCE THE FRAZIL ICE VARIABLE
        !$OMP SINGLE
        FRAZILICE1 = FRAZILICE
        NICE = 0
        !$OMP END SINGLE NOWAIT
        
      ENDIF
    ENDIF

  ENDIF

  IF( ISTRAN(3) == 1 .AND. ISCDCA(3) < 4 )THEN
    !$OMP SINGLE PRIVATE(IT)
    !$  IT = OMP_GET_THREAD_NUM() + 1
    CALL CALTRAN(ISTL_,IS2TL_,3,3,DYE,DYE1,IT)  
    !$OMP END SINGLE NOWAIT
  ENDIF

  IF( ISTRAN(5) == 1 .AND. ISCDCA(5) < 4 )THEN
    !$OMP DO SCHEDULE(STATIC,1) PRIVATE(NT,M,IT)
    DO NT=1,NTOX  
      M=MSVTOX(NT)  
      !$  IT = OMP_GET_THREAD_NUM() + 1
      CALL CALTRAN(ISTL_,IS2TL_,5,M,TOX(1,1,NT),TOX1(1,1,NT),IT)       
    ENDDO  
    !$OMP END DO NOWAIT
  ENDIF  
  
  IF( ISTRAN(6) == 1 .AND. ISCDCA(6) < 4 )THEN  
    !$OMP DO SCHEDULE(STATIC,1) PRIVATE(NS,M,IT)
    DO NS=1,NSED  
      M=MSVSED(NS)  
      !$  IT = OMP_GET_THREAD_NUM() + 1
      CALL CALTRAN(ISTL_,IS2TL_,6,M,SED(1,1,NS),SED1(1,1,NS),IT)  
    ENDDO  
    !$OMP END DO NOWAIT
  ENDIF  
  
  IF( ISTRAN(7) == 1 .AND. ISCDCA(7) < 4 )THEN  
    !$OMP DO SCHEDULE(STATIC,1) PRIVATE(NS,M,IT)
    DO NS=1,NSND  
      M=MSVSND(NS)  
      !$  IT = OMP_GET_THREAD_NUM() + 1
      CALL CALTRAN(ISTL_,IS2TL_,7,M,SND(1,1,NS),SND1(1,1,NS),IT)  
    ENDDO  
    !$OMP END DO NOWAIT
  ENDIF  
  !$OMP END PARALLEL
  TSADV=TSADV+(DSTIME(0)-TTDS)
  
  ! *** COSMIC TRANSPORT REMOVED 2014-07 PMC
  
  ! *** 1D ADVECTI0N TRANSPORT CALCULATION  
  ! *** REMOVED 2004-09-19  PMC  

  ! ******************************************************************************************
  ! *** VERTICAL DIFFUSION IMPLICIT HALF STEP CALCULATION  
  IF( KC == 1 ) GOTO 1500
  TTDS=DSTIME(0)  

  IF( LADRY > 0 )THEN
    DO LP=1,LADRY
      L=LDRY(LP)  
      EEB(L)=1.
      CCLBTMP(L)=0. 
    ENDDO

    DO K=1,KC
      DO LP=1,LADRY
        L=LDRY(LP)  
        VCU(L,K)=0.
      ENDDO
    ENDDO
  ENDIF

  !$OMP PARALLEL DEFAULT(SHARED)
  
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,CCUBTMP,CCMBTMP,RCDZKK,RCDZKMK,NT,NS) SCHEDULE(DYNAMIC,1)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    ! *** BOTTOM LAYER
    DO LP=1,LLWET(KS,ND)
      L=LKWET(LP,KS,ND) 
      RCDZKK=-DELTD2*CDZKK(L,KSZ(L))  
      CCUBTMP=RCDZKK*HPI(L)*AB(L,KSZ(L))  
      CCMBTMP=1._8-CCUBTMP  
      EEB(L)=1._8/CCMBTMP  
      VCU(L,KSZ(L))=CCUBTMP*EEB(L)  
    ENDDO  
    IF( ISTRAN(1) >= 1 )THEN  
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        SAL(L,KSZ(L))=SAL(L,KSZ(L))*EEB(L)  
      ENDDO  
    ENDIF  
    IF( ISTRAN(2) >= 1 )THEN  
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        TEM(L,KSZ(L))=TEM(L,KSZ(L))*EEB(L)  
      ENDDO  
    ENDIF  
    IF( ISTRAN(3) >= 1 )THEN  
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        DYE(L,KSZ(L))=DYE(L,KSZ(L))*EEB(L)  
      ENDDO  
    ENDIF  
    IF( ISTRAN(5) >= 1 )THEN  
      DO NT=1,NTOX  
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND) 
          TOX(L,KSZ(L),NT)=TOX(L,KSZ(L),NT)*EEB(L)  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND) 
          SED(L,KSZ(L),NS)=SED(L,KSZ(L),NS)*EEB(L)  
        ENDDO  
      ENDDO  
    ENDIF
    IF( ISTRAN(7) >= 1 )THEN  
      DO NS=1,NSND  
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND) 
          SND(L,KSZ(L),NS)=SND(L,KSZ(L),NS)*EEB(L)  
        ENDDO  
      ENDDO  
    ENDIF  

    ! *** MIDDLE LAYERS
    DO K=2,KS  
      DO LP=1,LLWET(K-1,ND)
        L=LKWET(LP,K-1,ND) 
        RCDZKMK=-DELTD2*CDZKMK(L,K)  
        RCDZKK =-DELTD2*CDZKK(L,K)  
        CCLBTMP(L)=RCDZKMK*HPI(L)*AB(L,K-1)  
        CCUBTMP=RCDZKK*HPI(L)*AB(L,K)  
        CCMBTMP=1._8-CCLBTMP(L)-CCUBTMP  
        EEB(L)=1._8/(CCMBTMP-CCLBTMP(L)*VCU(L,K-1))  
        VCU(L,K)=CCUBTMP*EEB(L)  
      ENDDO  
      IF( ISTRAN(1) >= 1 )THEN  
        DO LP=1,LLWET(K-1,ND)
          L=LKWET(LP,K-1,ND) 
          SAL(L,K)=(SAL(L,K)-CCLBTMP(L)*SAL(L,K-1))*EEB(L)  
        ENDDO  
      ENDIF  
      IF( ISTRAN(2) >= 1 )THEN  
        DO LP=1,LLWET(K-1,ND)
          L=LKWET(LP,K-1,ND) 
          TEM(L,K)=(TEM(L,K)-CCLBTMP(L)*TEM(L,K-1))*EEB(L)  
        ENDDO  
      ENDIF  
      IF( ISTRAN(3) >= 1 )THEN  
        DO LP=1,LLWET(K-1,ND)
          L=LKWET(LP,K-1,ND) 
          DYE(L,K)=(DYE(L,K)-CCLBTMP(L)*DYE(L,K-1))*EEB(L)  
        ENDDO  
      ENDIF  
      IF( ISTRAN(5) >= 1 )THEN  
        DO NT=1,NTOX  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(6) >= 1 )THEN  
        DO NS=1,NSED  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            SED(L,K,NS)=(SED(L,K,NS)-CCLBTMP(L)*SED(L,K-1,NS))*EEB(L)  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(7) >= 1 )THEN  
        DO NS=1,NSND  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            SND(L,K,NS)=(SND(L,K,NS)-CCLBTMP(L)*SND(L,K-1,NS))*EEB(L)  
          ENDDO  
        ENDDO  
      ENDIF  
    ENDDO  

    ! *** TOP LAYER
    K=KC  
    DO LP=1,LLWET(KS,ND)
      L=LKWET(LP,KS,ND) 
      RCDZKMK=-DELTD2*CDZKMK(L,K)  
      CCLBTMP(L)=RCDZKMK*HPI(L)*AB(L,K-1)  
      CCMBTMP=1._8-CCLBTMP(L) 
      EEB(L)=1._8/(CCMBTMP-CCLBTMP(L)*VCU(L,K-1))  
    ENDDO  
    IF( ISTRAN(1) >= 1 )THEN  
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        SAL(L,K)=(SAL(L,K)-CCLBTMP(L)*SAL(L,K-1))*EEB(L)  
      ENDDO  
    ENDIF  
    IF( ISTRAN(2) >= 1 )THEN  
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        TEM(L,K)=(TEM(L,K)-CCLBTMP(L)*TEM(L,K-1))*EEB(L)  
      ENDDO  
    ENDIF  
    IF( ISTRAN(3) >= 1 )THEN  
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND) 
        DYE(L,K)=(DYE(L,K)-CCLBTMP(L)*DYE(L,K-1))*EEB(L)  
      ENDDO  
    ENDIF  
    IF( ISTRAN(5) >= 1 )THEN  
      DO NT=1,NTOX  
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND) 
          TOX(L,K,NT)=(TOX(L,K,NT)-CCLBTMP(L)*TOX(L,K-1,NT))*EEB(L)  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(6) >= 1 )THEN  
      DO NS=1,NSED  
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND) 
          SED(L,K,NS)=(SED(L,K,NS)-CCLBTMP(L)*SED(L,K-1,NS))*EEB(L)  
        ENDDO  
      ENDDO  
    ENDIF  
    IF( ISTRAN(7) >= 1 )THEN  
      DO NS=1,NSND  
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND) 
          SND(L,K,NS)=(SND(L,K,NS)-CCLBTMP(L)*SND(L,K-1,NS))*EEB(L)  
        ENDDO  
      ENDDO  
    ENDIF  

    ! *** FINAL PASS
    DO K=KS,1,-1  
      IF( ISTRAN(1) >= 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          SAL(L,K)=SAL(L,K)-VCU(L,K)*SAL(L,K+1)  
        ENDDO  
      ENDIF  
      IF( ISTRAN(2) >= 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          TEM(L,K)=TEM(L,K)-VCU(L,K)*TEM(L,K+1)  
        ENDDO  
      ENDIF  
      IF( ISTRAN(3) >= 1 )THEN  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          DYE(L,K)=DYE(L,K)-VCU(L,K)*DYE(L,K+1)  
        ENDDO  
      ENDIF  
      IF( ISTRAN(5) >= 1 )THEN  
        DO NT=1,NTOX  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TOX(L,K,NT)=TOX(L,K,NT)-VCU(L,K)*TOX(L,K+1,NT)  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(6) >= 1 )THEN  
        DO NS=1,NSED  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SED(L,K,NS)=SED(L,K,NS)-VCU(L,K)*SED(L,K+1,NS)  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(7) >= 1 )THEN  
        DO NS=1,NSND  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SND(L,K,NS)=SND(L,K,NS)-VCU(L,K)*SND(L,K+1,NS)  
          ENDDO  
        ENDDO  
      ENDIF  
    ENDDO  
  ENDDO     ! *** END OF DOMAIN
  !$OMP END DO
  !$OMP END PARALLEL

  TVDIF=TVDIF+(DSTIME(0)-TTDS)  
  ! *** END OF VERTICAL DIFFUSION STEP
  
1500 CONTINUE  

  ! ***  COMPUTE TOTALS FOR SEDIMENT TRANSPORT (AFTER ADVECTION/DIFFUSION)
  IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
    !$OMP PARALLEL DEFAULT(SHARED)
    ! *** ZERO SED SEDIMENT ACCUMULATION ARRAYS
    IF( ISTRAN(6) > 0 )THEN
      !$OMP DO PRIVATE(ND,L,K,LF,LL,LP,NS) SCHEDULE(DYNAMIC)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
      
        ! *** SEDIMENT BED  PMC - 2016-12-04 - NOT USED.  COMMENTED OUT
        !DO K=1,KB  
        !  DO LP=LF,LL
        !    L=LWET(LP)  
        !    SEDBT(L,K)=0.  
        !  ENDDO  
        !ENDDO  

        ! *** WATER COLUMN
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SEDT(L,K) = 0.  
          ENDDO  
        ENDDO  

        ! *** SEDIMENT BED  PMC - 2016-12-04 - NOT USED.  COMMENTED OUT
        !DO K=1,KB  
        !  DO NS=1,NSED  
        !    DO LP=LF,LL
        !      L=LWET(LP)  
        !      SEDBT(L,K)=SEDBT(L,K)+SEDB(L,K,NS)  
        !    ENDDO  
        !  ENDDO  
        !ENDDO  
        DO NS=1,NSED  
          DO K=1,KC  
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              SEDT(L,K) = SEDT(L,K) + SED(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO
    ENDIF
  
    ! *** ZERO SND SEDIMENT ACCUMULATION ARRAYS
    IF( ISTRAN(7) > 0 )THEN
      !$OMP DO PRIVATE(ND,L,K,LF,LL,LP,NS) SCHEDULE(DYNAMIC)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)

        ! *** SEDIMENT BED  PMC - 2016-12-04 - NOT USED.  COMMENTED OUT
        !DO K=1,KB  
        !  DO LP=LF,LL
        !    L=LWET(LP)  
        !    SNDBT(L,K)=0.  
        !  ENDDO  
        !ENDDO  

        ! *** WATER COLUMN
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SNDT(L,K) = 0.  
          ENDDO  
        ENDDO  

        ! *** SEDIMENT BED  PMC - 2016-12-04 - NOT USED.  COMMENTED OUT
        !DO NS=1,NSND  
        !  DO K=1,KB  
        !    DO LP=LF,LL
        !      L=LWET(LP)  
        !      SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NS)  
        !    ENDDO  
        !  ENDDO  
        !ENDDO  

        DO NS=1,NSND  
          DO K=1,KC  
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              SNDT(L,K) = SNDT(L,K) + SND(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO
    ENDIF
    !$OMP END PARALLEL
  ENDIF

  ! ***************************************************************************************
  ! *** SURFACE AND INTERNAL HEAT SOURCE-SINK CALCULATION  
  IF( ISTRAN(2) >= 1 )THEN
    TTDS=DSTIME(0)
    CALL CALHEAT(ISTL_)
    THEAT=THEAT+(DSTIME(0)-TTDS)  
  ENDIF

  ! *** FULL IMPLICIT DYE AND TOXIC CONTAMINANT DECAY/GROWTH CALCULATION  
  IF( ISTRAN(3) >= 1 )THEN
    IF( RKDYE == 1000.0 )THEN   
      ! *** Age of Water
      DAGE=DELT/86400.
      !$OMP PARALLEL DO PRIVATE(ND,K,LP,L)
      DO ND=1,NDM  

        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            DYE(L,K)=DYE(L,K)+DAGE
          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
  
    ELSEIF( RKDYE /= 0.0 )THEN
      IF( RKDYE < 0.0 )THEN
        CDYETMP=EXP(-RKDYE*DELT)    
      ELSE
        CDYETMP=1./(1.+DELT*RKDYE)  
      ENDIF

      !$OMP PARALLEL DO PRIVATE(ND,K,LP,L)
      DO ND=1,NDM  

        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            DYE(L,K)=CDYETMP*DYE(L,K)  
          ENDDO
        ENDDO  
      ENDDO 
      !$OMP END PARALLEL DO
    ENDIF 
  ENDIF  

  ! **************************************************************************************************
  ! *** BOTTOM AND INTERNAL SEDIMENT AND TOXIC CONTAMINATION SOURCE-SINK CALCULATION  

  ! *** SEDIMENT AND TOXICS SETTLING,DEPOSITION,RESUSPENSION,ETC  
  IF( ( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 ) .AND. TIMEDAY >= SEDSTART )THEN 
    IF( IS2TIM >= 1 )THEN  
      ! *** FOR TWO TIME LEVEL SIMULATION  
      SEDTIME = SEDTIME + DELT
      IF( SEDTIME >= SEDSTEP )THEN 
        DTSED = SEDTIME
        CALL SSEDTOX  
        SEDTIME = 0.0  
      ENDIF  
    ELSE 
      ! *** FOR THREE TIME LEVEL SIMULATION  
      IF( NCTBC == 1 )THEN
        DTSED = FLOAT(NTSTBC)*DT
        CALL SSEDTOX  
        SEDTIME = 0.0  
      ENDIF
    ENDIF  
  ENDIF  
  
  ! *** OPTIONAL MASS BALANCE CALCULATION  
  IF( IS2TIM == 0 )THEN  
    IF( ISTL_/=2 .AND. ISBAL >= 1 )THEN  
      CALL CALBAL2  
      CALL CALBAL3  
      NTMP=MOD(N,2)  
      IF( NTMP == 0 )THEN  
        CALL CBALEV2  
        CALL CBALEV3  
      ELSE  
        CALL CBALOD2  
        CALL CBALOD3  
      ENDIF  
    ENDIF  
  ENDIF  

  ! *** CALLS TO TWO-TIME LEVEL BALANCES  
  !  
  IF( IS2TIM >= 1 )THEN  
    IF( ISBAL >= 1 )THEN  
      CALL BAL2T2
      CALL BAL2T3B(1)  
    ENDIF  
  ENDIF  

  ! *** SEDIMENT BUDGET CALCULATION    (DLK 10/15)  
  !  
  IF( IS2TIM == 0 )THEN  
    IF( ISTL_/=2 .AND. ISSBAL >= 1 )THEN  
      CALL BUDGET2  
      CALL BUDGET3  
    ENDIF  
  ENDIF  
  
  ! *** 2014-09 - REMOVED DATA ASSIMILATION CODE (PMC)

  RETURN
  
END  
      
