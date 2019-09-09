SUBROUTINE HDMT

  ! *** SUBROUTINE HDMT EXECUTES THE FULL HYDRODYNAMIC AND MASS TRANSPORT
  ! *** TIME INTERGATION
  !
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2014-08-12    D H CHUNG         SET EXPLICIT PRECISIONS OF INTEGER & REAL
  !    2012-09       Chung Dang        OMP
  !    05/01/2002    John Hamrick      Modified calls to calbal and budget subroutines
  !    09-22-2004    Paul M. Craig     Merged DS and TT versions with the 06-04-2004 TT code
  !
  !
  !**********************************************************************C
  !
  USE GLOBAL
  USE DRIFTER
  USE WINDWAVE ,ONLY:WINDWAVEINIT,WINDWAVETUR,READWAVECELLS
  USE HIFREQOUT
  USE RESTART_MODULE
  USE EFDCOUT
  USE CALCSERMOD,ONLY: CALCSER
#ifdef NCOUT
  USE MOD_NETCDF
#endif

  IMPLICIT NONE
  INTEGER(IK4) :: NS
  INTEGER :: L,ND,NTMP1,NTMP2,NTMP,K,IMAX,JMAX,KMAX
  INTEGER :: IMIN,JMIN,KMIN,NMD,ICALLTP,LS,LP
  INTEGER :: ILOGC
  INTEGER :: LN,LNW,LSE,LF,LL,LE,LW
  
  REAL :: SALMIN,HPPTMP,WTM,WTMP,TMP
  REAL :: DELVOL,SALMAX,TAUB2,DELTD2,DZDDELT
  REAL :: TAUBC,TAUBC2,UTMP,VTMP,CURANG,USGZ,VSGZ
  REAL :: CTIM
  REAL(RKD), EXTERNAL :: DSTIME
  REAL(RKD)           :: TTDS, T1TMP, TIMELAST         ! MODEL TIMING TEMPORARY VARIABLES
  REAL(RKD)           :: TIMEHARM

#ifdef _WIN  
  LOGICAL, EXTERNAL :: KEY_PRESSED
  LOGICAL, EXTERNAL :: ISEXIT
#endif
  REAL(8)           :: DEL !RESTTIME
  LOGICAL           :: RES

  WRITE(*,'(A)')'STARTING HDMT 3TL'
  T1TMP=DSTIME(0)
  FOURDPI=4./PI
  RES = .TRUE.

  ! *** INITIALIZE TIME FOR INITIAL CONDITION COMPUTATIONS
  TIMESEC = DBLE(TCON)*DBLE(TBEGIN)
  TIMEDAY = DBLE(TCON)*DBLE(TBEGIN)/86400._8
  TIMEHARM = DBLE(TCON)*DBLE(TBEGIN) + 300.   ! *** 05 MINUTE INCREMENTS

  !**********************************************************************C
  !
  ! *** INITIALIZE COURNT NUMBER DIAGNOSTICS
  !
  !**********************************************************************!  
  ! *** REINITIALIZE VARIABLES  
  IF( ISRESTI == 0 )THEN 
    DO L=2,LA  
      H2P(L)=HP(L)  
      H1P(L)=HP(L)  
      H1U(L)=HU(L)  
      H1UI(L)=HUI(L)  
      H1V(L)=HV(L)  
      H1VI(L)=HVI(L)  
      UHDY1E(L)=UHDYE(L)  
      VHDX1E(L)=VHDXE(L)  
    ENDDO  

    DO K=1,KC  
      DO L=2,LA  
        HPK(L,K) = HP(L)*DZC(L,K)
        H2PK(L,K) = HPK(L,K)
        H1PK(L,K) = HPK(L,K)
        U1(L,K)=U(L,K)  
        V1(L,K)=V(L,K)  
        UHDY1(L,K)=UHDY(L,K)  
        VHDX1(L,K)=VHDX(L,K)  
      ENDDO  
    ENDDO  
  ENDIF

  !**********************************************************************!  
  ! *** INITIALIZE COURANT NUMBER DIAGNOSTICS  
  IF( ISINWV == 1 )THEN
    DO K=1,KC
      DO L=2,LA
        CFLUUU(L,K)=0.
        CFLVVV(L,K)=0.
        CFLWWW(L,K)=0.
        CFLCAC(L,K)=0.
      ENDDO
    ENDDO
  ENDIF
  ILOGC=0

  !**********************************************************************C
  ! *** CALCULATE U AT V AND V AT U USING ENERGY CONSERVING WEIGHTING
  ! *** CALCULATE VELOCITY GRADIENTS
  !----------------------------------------------------------------------C
  DO L=2,LA
    LN=LNC(L)
    LS=LSC(L)
    LE=LEC(L)
    LW=LWC(L)
    LNW=LNWC(L)
    LSE=LSEC(L)

    UV(L)  = 0.25*( HP(LS) *(U(LSE,KSZV(L)) +U(LS,KSZV(L))  ) + HP(L) *(U(LE,KSZV(L)) +U(L,KSZV(L))  ) )*HVI(L)
    U1V(L) = 0.25*( H1P(LS)*(U1(LSE,KSZV(L))+U1(LS,KSZV(L)) ) + H1P(L)*(U1(LE,KSZV(L))+U1(L,KSZV(L)) ) )*H1VI(L)  
    VU(L)  = 0.25*( HP(LW) *(V(LNW,KSZU(L)) +V(LW,KSZU(L))  ) + HP(L) *(V(LN,KSZU(L)) +V(L,KSZU(L))  ) )*HUI(L)  
    V1U(L) = 0.25*( H1P(LW)*(V1(LNW,KSZU(L))+V1(LW,KSZU(L)) ) + H1P(L)*(V1(LN,KSZU(L))+V1(L,KSZU(L)) ) )*H1UI(L)  

  ENDDO

  !**********************************************************************C
  ! *** CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
  IF( ISWAVE == 1 ) CALL WAVEBL
  IF( ISWAVE == 2 ) CALL WAVESXY
  IF( ISWAVE >= 3 .AND. NWSER > 0 )THEN
    CALL WINDWAVEINIT
    CALL WINDWAVETUR   !DHC FIRST CALL

    ! *** READ IN WAVE COMPUTATIONAL CELL LIST
    IF( IUSEWVCELLS /= 0 )THEN
      CALL READWAVECELLS
    ENDIF
  ENDIF

#ifdef NCOUT
  ! ** NETCDF INIT
  NC_DATESTAMP = ''
  IF (NCDFOUT > 0) THEN
    CALL READCORN
    IF( TIMEDAY >= TBEGNCDF ) CALL NETCDF_WRITE(NC_ID)
  ENDIF

  IF( HFREOUT == 1 )THEN
    DO NS=1,NSUBSET
      CALL HFREHYOUT(1,NS)
      CALL HFREWCOUT(1,NS)
      CALL HFREWQOUT(1,NS)
      CALL HFRERPEMOUT(1,NS)
    ENDDO
  ENDIF
#endif  

  !**********************************************************************C
  ! *** FIRST CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS
  ISTL=3
  IS2TL=0
  CALL CALTBXY(ISTL,IS2TL)

  !**********************************************************************C
  ! *** CALCULATE HORIZONTAL VISCOSITY AND DIFFUSIVE MOMENTUM FLUXES
  IF( ISHDMF >= 1 ) CALL CALHDMF3(ISTL)

  !**********************************************************************C
  ! *** CALCULATE BOTTOM AND SURFACE STRESS AT TIME LEVEL (N-1) AND N
  !----------------------------------------------------------------------C
  N=-1
  CALL CALTSXY

  DO L=2,LA
    TBX1(L)=(AVCON1*H1UI(L) + STBX(L)*SQRT(V1U(L)*V1U(L) + U1(L,KSZ(L))*U1(L,KSZ(L)))) * U1(L,KSZ(L))
    TBY1(L)=(AVCON1*H1VI(L) + STBY(L)*SQRT(U1V(L)*U1V(L) + V1(L,KSZ(L))*V1(L,KSZ(L)))) * V1(L,KSZ(L))
    TSX1(L)=TSX(L)
    TSY1(L)=TSY(L)
  ENDDO

  !**********************************************************************!
  ! *** SECOND CALL TO INITIALIZE BOTTOM STRESS COEFFICIENTS
  CALL CALTBXY(ISTL,IS2TL)

  !**********************************************************************C
  ! *** SET BOTTOM AND SURFACE STRESSES
  DO L=2,LA
    USGZ=U(L,KSZU(L))
    VSGZ=V(L,KSZV(L))
    TBX(L)=(AVCON1*HUI(L)+STBX(L)*SQRT(VU(L)*VU(L)+USGZ*USGZ))*USGZ
    TBY(L)=(AVCON1*HVI(L)+STBY(L)*SQRT(UV(L)*UV(L)+VSGZ*VSGZ))*VSGZ
  ENDDO

  N=0
  CALL CALTSXY
 
  ! *** SET DEPTH DEVIATION FROM UNIFORM FLOW ON FLOW FACES
  DO L=2,LA
    HDFUFX(L)=1.
    HDFUFY(L)=1.
    HDFUF(L)=1.
  ENDDO
 
  ! *** SET CORNER CORRECTIONS
  DO L=2,LA
    WCOREST(L)=1.
    WCORWST(L)=1.
    WCORNTH(L)=1.
    WCORSTH(L)=1.
  ENDDO
 
  !**********************************************************************!
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED
  IF( ISWAVE == 0 )THEN
    DO L=2,LA
      TVAR3W(L)=TSX1(LEC(L))
      TVAR3S(L)=TSY1(LNC(L))
      TVAR3E(L)=TBX1(LEC(L))
      TVAR3N(L)=TBY1(LNC(L))
    ENDDO

    DO L=2,LA
      QQ1(L,0)  = 0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX1(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY1(L))**2)
      QQ1(L,KC) = 0.5*CTURB2*SQRT((RSSBCE(L)*TVAR3W(L)+RSSBCW(L)*TSX1(L))**2 + (RSSBCN(L)*TVAR3S(L)+RSSBCS(L)*TSY1(L))**2)
    ENDDO

    DO L=2,LA
      TVAR3W(L) = TSX(LEC(L))
      TVAR3S(L) = TSY(LNC(L))
      TVAR3E(L) = TBX(LEC(L))
      TVAR3N(L) = TBY(LNC(L))
    ENDDO

    DO L=2,LA
      QQ(L,0)  = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L))**2  + (RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L))**2 )  
      QQ(L,KC) = 0.5*CTURB2*SQRT( (RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L))**2  + (RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L))**2 )  
      QQSQR(L,0) = SQRT(QQ(L,0))  
    ENDDO
  ENDIF    ! *** END OF ISWAVE=0

  !**********************************************************************!  
  ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED  
  !----------------------------------------------------------------------!  
  IF( ISWAVE >= 1 )THEN  
    DO L=2,LA
      TVAR3S(L)=TSY1(LNC(L))
      TVAR3W(L)=TSX1(LEC(L))
      TVAR3E(L)=TBX1(LEC(L))
      TVAR3N(L)=TBY1(LNC(L))
    ENDDO

    DO L=2,LA
      TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2
      TAUBC=0.5*SQRT(TAUBC2)
      CTAUC(L)=TAUBC
      UTMP = 0.5*STCUV(L)*( U1(LEC(L),KSZ(L)) + U1(L,KSZ(L)) )+1.E-12
      VTMP = 0.5*STCUV(L)*( V1(LNC(L),KSZ(L)) + V1(L,KSZ(L)) )
      CURANG=ATAN2(VTMP,UTMP)
      TAUB2 = TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
      TAUB2 = MAX(TAUB2,0.)
      QQ1(L,0)  = CTURB2*SQRT(TAUB2)
      QQ1(L,KC) = 0.5*CTURB2*SQRT( (TVAR3W(L)+TSX1(L))**2 + (TVAR3S(L)+TSY1(L))**2 )
    ENDDO

    DO L=2,LA
      TVAR3S(L)=TSY(LNC(L))
      TVAR3W(L)=TSX(LEC(L))
      TVAR3E(L)=TBX(LEC(L)   )
      TVAR3N(L)=TBY(LNC(L))
    ENDDO

    DO L=2,LA
      TAUBC2 = 0.25*( (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2 + (RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2 )
      TAUBC = SQRT(TAUBC2)
      CTAUC(L) = TAUBC
      UTMP = 0.5*STCUV(L)*(U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)))+1.E-12
      VTMP = 0.5*STCUV(L)*(V(LN,KSZV(LNC(L)))     + V(L,KSZV(L)))
      CURANG = ATAN2(VTMP,UTMP)
      TAUB2 =  TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
      TAUB2 =  MAX(TAUB2,0.)
      QQ(L,0)  = CTURB2*SQRT(TAUB2)
      QQ(L,KC) = 0.5*CTURB2*SQRT( (TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2 )
      QQSQR(L,0) = SQRT(QQ(L,0))            
    ENDDO
  ENDIF

  ! ***  SET GRAIN STRESS
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    DO L=2,LA
      TAUBSED(L)=QQ(L,0)/CTURB2
      TAUBSND(L)=QQ(L,0)/CTURB2
    ENDDO
  ENDIF
  
  !**********************************************************************C
  !
  ! ***  SET SWITCHES FOR THREE TIME LEVEL STEP
  !
  ISTL=3
  IS2TL=0
  DELT=DT2
  DELTD2=DT
  DZDDELT=DZ/DELT
  ROLD=0.
  RNEW=1.

  !**********************************************************************C
  ! *** BEGIN TIME LOOP FOR FULL HYDRODYNAMIC AND MASS TRANSPORT
  ! *** CALCULATION

  ! *** SET CYCLE COUNTER AND CALL TIMER
  !**********************************************************************
  NTIMER=0
  N=0
  TIMELAST = TIMEDAY + TIMERST
  
  ! *** INITIALIZE & RECORD TIME
  TIMEDAY=TCON*TBEGIN/86400.
  !RESTTIME = TIMEDAY+HREST/24.
  CALL TIMELOG(N,TIMEDAY,OUTDIR)

  NTIMER=1

  !***************************************************************************
  !***************************************************************************
  ! *** BEGINNING OF THE MAIN TIME ITERATION LOOP FOR TWO TIME LEVEL SOLUTION
  !***************************************************************************
  !***************************************************************************
  DO 1000 N=1,NTS

    TIMESEC = DBLE(DT)*DBLE(N)+DBLE(TCON)*DBLE(TBEGIN)
    TIMEDAY = TIMESEC/86400._8
    NITER = N
    
    ILOGC=ILOGC+1

    IF( N <= NLTS )THEN
      SNLT=0.
    ELSEIF( N > NLTS .AND. N <= NTTS )THEN
      NTMP1=N-NLTS
      NTMP2=NTTS-NLTS+1
      SNLT=FLOAT(NTMP1)/FLOAT(NTMP2)
    ELSE
      SNLT=1.
    ENDIF

    ! *** SET THE GRAVITY ACCELERATION
    IF( N <= NTSVB )THEN
      GP=GPO*(FLOAT(N)/FLOAT(NTSVB))
    ELSE
      GP=GPO
    ENDIF

    !----------------------------------------------------------------------C
    !
    ! *** INITIALIZE VOLUME, MASS, MOMENTUM, AND ENERGY BALANCE
    !
    IF( NCTBC /= NTSTBC .AND. ISBAL >= 1 )THEN
      CALL CALBAL1
      NTMP=MOD(N,2)
      IF( NTMP == 0 )THEN
        CALL CBALEV1
      ELSE
        CALL CBALOD1
      ENDIF
    ENDIF

    !  ** INITIALIZE SEDIMENT BUDGET CALCULATION   (DLK 10/15)
    IF( NCTBC /= NTSTBC .AND. ISSBAL >= 1 )THEN
      CALL BUDGET1
    ENDIF

    !----------------------------------------------------------------------C
    ! *** REENTER HERE FOR TWO TIME LEVEL CORRECTION
    500 CONTINUE
    
    !**********************************************************************C
    ! *** CALCULATE VERTICAL VISCOSITY AND DIFFUSIVITY AT TIME LEVEL (N)
    TTDS=DSTIME(0)
    IF( KC > 1 )THEN
      CALL CALAVB
    ENDIF
    TAVB=TAVB+(DSTIME(0)-TTDS)

    !**********************************************************************C
    ! *** CALCULATE WAVE BOUNDARY LAYER AND WAVE REYNOLDS STRESS FORCINGS
    IF( ISTL == 3 )THEN
      TTDS=DSTIME(0)
      IF( ISWAVE == 1 ) CALL WAVEBL
      IF( ISWAVE == 2 ) CALL WAVESXY
      IF( ISWAVE >= 3 .AND. NWSER > 0 ) CALL WINDWAVETUR   !DHC NEXT CALL
      TTBXY=TTBXY+(DSTIME(0)-TTDS)  
    ENDIF

    !**********************************************************************C
    ! *** CALCULATE EXPLICIT MOMENTUM EQUATION TERMS
    ! *** CALEXP   -  PRODUCTION VERSION WITH HORIZONTAL MOMENTUM SOURCE
    ! ***             AND 3D IMPLICIT VEGETATION DRAG
    TTDS=DSTIME(0)
    CALL CALEXP (ISTL)
    TCEXP=TCEXP+(DSTIME(0)-TTDS)

    !**********************************************************************C
    ! *** UPDATE TIME VARIABLE VOLUME SOURCES AND SINKS, CONCENTRATIONS,
    ! *** VEGETATION CHARACTERISTICS AND SURFACE ELEVATIONS
    CALL CALCSER (ISTL)
    CALL CALVEGSER
    CALL CALQVS (ISTL)
    PSERT(0)=0.
    IF( NPSER >= 1 ) CALL CALPSER

    !**********************************************************************C
    ! *** SOLVE EXTERNAL MODE EQUATIONS FOR P, UHDYE, AND VHDXE
    TTDS=DSTIME(0)
    CALL CALPUV9C(ISTL)
    TPUV=TPUV+(DSTIME(0)-TTDS)

    !**********************************************************************C
    ! *** ADVANCE TIME SNAPSHOT INTERNAL VARIABLES FOR THREE TIME LEVEL STEP
    IF( ISTL == 3 )THEN
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,LP,L)
      DO ND=1,NDM  
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            UHDYF2(L,K)=UHDYF1(L,K)   ! *** FLOW BY LAYER FOR THE ENTIRE DEPTH
            UHDYF1(L,K)=UHDYF(L,K)
            VHDXF2(L,K)=VHDXF1(L,K)
            VHDXF1(L,K)=VHDXF(L,K)
            
            UHDY2(L,K)=UHDY1(L,K)     ! *** FLOW BY LAYER 
            UHDY1(L,K)=UHDY(L,K)
            VHDX2(L,K)=VHDX1(L,K)
            VHDX1(L,K)=VHDX(L,K)
            
            U2(L,K)=U1(L,K)           ! *** N-2 TIME, NOT INTERVAL AVERAGE
            V2(L,K)=V1(L,K)           ! *** N-2 TIME, NOT INTERVAL AVERAGE
            U1(L,K)=U(L,K)
            V1(L,K)=V(L,K)
            W2(L,K)=W1(L,K)           ! *** N-2 TIME, NOT INTERVAL AVERAGE
            W1(L,K)=W(L,K)
          ENDDO
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ENDIF

    !**********************************************************************C
    ! *** ADVANCE TIME VARIABLE SURFACE WIND STRESS AND LOAD INTO INTERNAL
    ! *** MODE FORCING
    !----------------------------------------------------------------------C
    IF( ISTL == 3 )THEN
      TTDS=DSTIME(0)

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL
          L=LWET(LP)  
          TSX1(L)=TSX(L)
          TSY1(L)=TSY(L)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      CALL CALTSXY

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL
          L=LWET(LP)  
          DU(L,KS) = DU(L,KS) - CDZUU(L,KS)*TSX(L)
          DV(L,KS) = DV(L,KS) - CDZUV(L,KS)*TSY(L)
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

      TTBXY=TTBXY+(DSTIME(0)-TTDS)
    ENDIF

    !**********************************************************************C
    ! *** SOLVE INTERNAL SHEAR MODE EQUATIONS FOR U, UHDY, V, VHDX, AND W
    ! *** USING THE INTERNAL SHEARS STORED IN DU & DV 
    TTDS=DSTIME(0)
    IF( KC > 1 )THEN
      CALL CALUVW (ISTL)
    ELSE
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL
          L=LWET(LP)  
          UHDYF(L,1) = UHDYE(L)
          UHDY(L,1)  = UHDYE(L)
          U(L,1)=UHDYE(L)*HUI(L)*DYIU(L)
          VHDXF(L,1) = VHDXE(L)
          VHDX(L,1)  = VHDXE(L)
          V(L,1)=VHDXE(L)*HVI(L)*DXIV(L)
          W(L,1)=0.
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO
      CALL CALUVW (ISTL)
    ENDIF
    TUVW=TUVW+(DSTIME(0)-TTDS)

    !**********************************************************************C
    ! *** CALCULATE SALINITY, TEMPERATURE, DYE AND SEDIMENT CONCENTRATIONS
    ! *** AT TIME LEVEL (N+1)
    IF( ISTRANACTIVE > 0 )  CALL CALCONC (ISTL,IS2TL)

    !----------------------------------------------------------------------C
    ! *** CHECK RANGE OF SALINITY AND DYE CONCENTRATION
    IF( ISMMC == 1 .AND. DEBUG )THEN

      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
        DO L=2,LA
          IF( SAL(L,K) > SALMAX )THEN
            SALMAX=SAL(L,K)
            IMAX=IL(L)
            JMAX=JL(L)
            KMAX=K
          ENDIF
          IF( SAL(L,K) < SALMIN )THEN
            SALMIN=SAL(L,K)
            IMIN=IL(L)
            JMIN=JL(L)
            KMIN=K
          ENDIF
        ENDDO
      ENDDO

      WRITE(6,6001)N
      WRITE(6,6002)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6003)SALMIN,IMIN,JMIN,KMIN

      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
        DO L=2,LA
          IF( DYE(L,K) > SALMAX )THEN
            SALMAX=DYE(L,K)
            IMAX=IL(L)
            JMAX=JL(L)
            KMAX=K
          ENDIF
          IF( DYE(L,K) < SALMIN )THEN
            SALMIN=DYE(L,K)
            IMIN=IL(L)
            JMIN=JL(L)
            KMIN=K
          ENDIF
        ENDDO
      ENDDO

      WRITE(6,6004)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6005)SALMIN,IMIN,JMIN,KMIN

      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
        DO L=2,LA
          IF( SFL(L,K) > SALMAX )THEN
            SALMAX=SFL(L,K)
            IMAX=IL(L)
            JMAX=JL(L)
            KMAX=K
          ENDIF
          IF( SFL(L,K) < SALMIN )THEN
            SALMIN=SFL(L,K)
            IMIN=IL(L)
            JMIN=JL(L)
            KMIN=K
          ENDIF
        ENDDO
      ENDDO

      WRITE(6,6006)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6007)SALMIN,IMIN,JMIN,KMIN

    ENDIF

    IF( ISMMC == 2 .AND. DEBUG )THEN
      SALMAX=-100000.
      SALMIN=100000.
      DO K=1,KC
        DO L=2,LA
          IF( TEM(L,K) > SALMAX )THEN
            SALMAX=TEM(L,K)
            IMAX=IL(L)
            JMAX=JL(L)
            KMAX=K
          ENDIF
          IF( TEM(L,K) < SALMIN )THEN
            SALMIN=TEM(L,K)
            IMIN=IL(L)
            JMIN=JL(L)
            KMIN=K
          ENDIF
        ENDDO
      ENDDO

      WRITE(6,6001)N
      WRITE(6,6008)SALMAX,IMAX,JMAX,KMAX
      WRITE(6,6009)SALMIN,IMIN,JMIN,KMIN

    ENDIF

   6001 FORMAT('  N=',I10)
   6002 FORMAT('  SALMAX=',F14.4,5X,'I,J,K=',(3I10))
   6003 FORMAT('  SALMIN=',F14.4,5X,'I,J,K=',(3I10))
   6004 FORMAT('  DYEMAX=',F14.4,5X,'I,J,K=',(3I10))
   6005 FORMAT('  DYEMIN=',F14.4,5X,'I,J,K=',(3I10))
   6006 FORMAT('  SFLMAX=',F14.4,5X,'I,J,K=',(3I10))
   6007 FORMAT('  SFLMIN=',F14.4,5X,'I,J,K=',(3I10))
   6008 FORMAT('  TEMMAX=',F14.4,5X,'I,J,K=',(3I10))
   6009 FORMAT('  TEMMIN=',F14.4,5X,'I,J,K=',(3I10))

    !**********************************************************************C
    !
    ! *** CALCULATE SHELL FISH LARVAE AND/OR WATER QUALITY CONSTITUENT
    ! *** CONCENTRATIONS AT TIME LEVEL (N+1) AFTER SETTING DOUBLE TIME
    ! *** STEP TRANSPORT FIELD
    !
    !----------------------------------------------------------------------C
    IF( ISWQFLUX == 1 )THEN
      NTMP=MOD(N,2)
      IF( NTMP == 0 .AND. ISTL == 3 )THEN

        ! *** CALCULATE CONSERVATION OF VOLUME FOR THE WATER QUALITY ADVECTION
        !$OMP PARALLEL DEFAULT(SHARED)
        !$OMP DO PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)  
            !UHDY2E(L)=0.
            !VHDX2E(L)=0.
            HWQ(L) = 0.25*(H2P(L) + 2.*H1P(L)+HP(L))
          ENDDO
        ENDDO
        !$OMP END DO

        !$OMP DO PRIVATE(ND,K,LP,L)
        DO ND=1,NDM
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              UHDYWQ(L,K) = 0.5*(UHDY1(L,K)+UHDY2(L,K))
              VHDXWQ(L,K) = 0.5*(VHDX1(L,K)+VHDX2(L,K))
              !UHDY2E(L) = UHDY2E(L)+UHDYWQ(L,K)
              !VHDX2E(L) = VHDX2E(L)+VHDXWQ(L,K)
              UWQ(L,K) = 0.5*(U1(L,K)+U2(L,K))
              VWQ(L,K) = 0.5*(V1(L,K)+V2(L,K))
            ENDDO
          ENDDO
        ENDDO
        !$OMP END DO

        !$OMP DO PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)  
            TVAR3E(L) = UHDY2E(LEC(L))
            TVAR3N(L) = VHDX2E(LNC(L))
          ENDDO
        ENDDO
        !$OMP END DO

        !$OMP DO PRIVATE(ND,K,LP,L)
        DO ND=1,NDM
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              TVAR2E(L,K) = UHDYWQ(LEC(L),K)
              TVAR2N(L,K) = VHDXWQ(LNC(L),K)
            ENDDO
          ENDDO
        ENDDO
        !$OMP END DO

        !$OMP DO PRIVATE(ND,K,LP,L)
        DO ND=1,NDM
          DO K=1,KS
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              WWQ(L,K) = SWB(L)*( WWQ(L,K-1)-DZC(L,K)*(TVAR2E(L,K)-UHDYWQ(L,K)-TVAR3E(L)+UHDY2E(L)             &
                                                  +TVAR2N(L,K)-VHDXWQ(L,K)-TVAR3N(L)+VHDX2E(L))*DXYIP(L) ) &
                       + SWB(L)*( QSUM(L,K)-DZC(L,K)*QSUME(L) )*DXYIP(L)
            ENDDO
          ENDDO
        ENDDO
        !$OMP END DO

        IF( ISICM >= 1 )THEN
          !$OMP DO PRIVATE(ND,LF,LL,LP,L,HPPTMP)
          DO ND=1,NDM
            LF=(ND-1)*LDMWET+1  
            LL=MIN(LF+LDMWET-1,LAWET)
            DO LP=LF,LL
              L=LWET(LP)  
              HPPTMP=H2WQ(L)+DT2*DXYIP(L)*( QSUME(L) - (TVAR3E(L)-UHDY2E(L)+TVAR3N(L)-VHDX2E(L)) )
              HWQ(L)=SPB(L)*HPPTMP+(1.-SPB(L))*HWQ(L)
            ENDDO
          ENDDO
          !$OMP END DO
        ENDIF
        !$OMP END PARALLEL

        ! *** ADD CHANNEL INTERACTIONS  
        IF( MDCHH >= 1 )THEN  
          DO NMD=1,MDCHH  
            IF( MDCHTYP(NMD) == 1 )THEN  
              HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD))  +DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD))  
              HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD))  -DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))  
            ENDIF  
            IF( MDCHTYP(NMD) == 2 )THEN  
              HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))  
              HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD)) - DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))  
            ENDIF  
            IF( MDCHTYP(NMD) == 3 )THEN  
              HWQ(LMDCHH(NMD))=HWQ(LMDCHH(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANU(NMD)) + DT2*DXYIP(LMDCHH(NMD))*(QCHANV(NMD))  
              HWQ(LMDCHU(NMD))=HWQ(LMDCHU(NMD)) - DT2*DXYIP(LMDCHU(NMD))*(QCHANU(NMD))  
              HWQ(LMDCHV(NMD))=HWQ(LMDCHV(NMD)) - DT2*DXYIP(LMDCHV(NMD))*(QCHANV(NMD))  
            ENDIF  
          ENDDO  
        ENDIF      ! *** END ADD CHANNEL INTERACTIONS  

        ! *** CALL WATER QAULITY KINETICS AND TRANSPORT 
        IF( ISTRAN(8) >= 1 ) CALL WQ3D(ISTL,IS2TL)
        IF( ISTRAN(4) >= 1 ) CALL CALSFT(ISTL,IS2TL)

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)  
            H2WQ(L)=HWQ(L)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO

      ENDIF
    ENDIF         ! *** END OF WQ SECTION

    !**********************************************************************C
    ! *** UPDATE BUOYANCY AND CALCULATE NEW BUOYANCY USING AN EQUATION OF STATE
    IF( BSC > 1.E-6 )THEN
      IF( ISTL == 3 )THEN
        CALL CALBUOY(.TRUE.)
      ELSE
        CALL CALBUOY(.FALSE.)
      ENDIF
    ENDIF

    IF( NCTBC /= NTSTBC .AND. ISBAL >= 1 )THEN
      CALL CALBAL4
      NTMP=MOD(N,2)
      IF( NTMP == 0 )THEN
        CALL CBALEV4
      ELSE
        CALL CBALOD4
      ENDIF
    ENDIF

    !**********************************************************************C
    !
    ! *** CALCULATE U AT V AND V AT U AT TIME LEVEL (N+1)

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,LN,LS,LNW,LSE,LE,LW)            &
    !$OMP             SHARED(NDM,LDMWET,LAWET,LWET,LEC,LNC,LSC,LWC,LNWC,LSEC,LSWC,LKSZ)   &
    !$OMP             SHARED(SUB,SVB,UV,VU,HP,U,V,U1V,V1U,HVI,HUI,KSZ,KSZU,KSZV,HPW,HPE,HPN,HPS) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      DO LP=LF,LL
        L=LWET(LP)  
        LN=LNC(L)  
        LS=LSC(L) 
        LE=LEC(L) 
        LW=LWC(L)
        LNW=LNWC(L)  
        LSE=LSEC(L)  
        
        U1V(L) = UV(L)  
        V1U(L) = VU(L)  

        IF( LKSZ(LS,KSZ(L)) )THEN
          UV(L)=0.5*( U(LE,KSZ(L)) + U(L,KSZ(L)) ) 
        ELSE
          UV(L)= 0.25*( HP(LS)*( U(LSE,KSZU(LSE)) + U(LS,KSZU(LS)) ) + HP(L)*( U(LE,KSZU(LE)) + U(L,KSZU(L)) ) )*HVI(L)
        ENDIF
        IF( LKSZ(LW,KSZ(L)) )THEN
          VU(L)= 0.5*( V(LN,KSZ(L)) + V(L,KSZ(L)) )
        ELSE
          VU(L)= 0.25*( HP(LW)*( V(LNW,KSZV(LNW)) + V(LW,KSZV(LW)) ) + HP(L)*( V(LN,KSZV(LN)) + V(L,KSZV(L)) ) )*HUI(L)
        ENDIF
      ENDDO  
    ENDDO
    !$OMP END PARALLEL DO

    !**********************************************************************C
    ! *** CALCULATE HORIZONTAL VISCOSITY AND MOMENTUM DIFFUSION FLUXES
    ! *** AT TIME LEVEL (N)
    IF( ISHDMF >= 1 )THEN
      TTDS=DSTIME(0)
      CALL CALHDMF3(ISTL)
      THMDF=THMDF+(DSTIME(0)-TTDS)  
    ENDIF

    !**********************************************************************C
    !
    ! *** UPDATE BOTTOM STRESSES AND SURFACE AND BOTTOM TURBULENT
    ! *** INTENSITIES
    IF( ISTL == 3 )THEN
      IF( ISCDMA == 2 )THEN
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)  
            TBX1(L)   = TBX(L)
            TBY1(L)   = TBY(L)
            QQ2(L,0)  = QQ(L,0)+QQ1(L,0)
            QQ2(L,KC) = QQ(L,KC)+QQ1(L,KC)
            QQ1(L,0)  = QQ(L,0)
            QQ1(L,KC) = QQ(L,KC)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO

      ELSE  ! *** IF( ISCDMA < 2 )THEN

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)  
            TBX1(L)   = TBX(L)
            TBY1(L)   = TBY(L)
            QQ2(L,0)  = QQ1(L,0)+QQ1(L,0)
            QQ2(L,KC) = QQ1(L,KC)+QQ1(L,KC)
            QQ1(L,0)  = QQ(L,0)
            QQ1(L,KC) = QQ(L,KC)
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
  
      ENDIF
    ENDIF

    !**********************************************************************C
    !
    ! *** CALCULATE BOTTOM STRESS AT LEVEL (N+1)

    TTDS=DSTIME(0)
    CALL CALTBXY(ISTL,IS2TL)

    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,USGZ,VSGZ)       &
    !$OMP             SHARED(NDM,LDMWET,LAWET,LWET,KSZ,KSZU,KSZV)          &
    !$OMP             SHARED(AVCON1, TBX,HUI,STBX,VU,U, TBY,HVI,STBY,UV,V)
    DO ND=1,NDM
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        USGZ=U(L,KSZU(L))
        VSGZ=V(L,KSZV(L))
        TBX(L) = ( AVCON1*HUI(L) + STBX(L)*SQRT( VU(L)*VU(L) + USGZ*USGZ ) )*USGZ  
        TBY(L) = ( AVCON1*HVI(L) + STBY(L)*SQRT( UV(L)*UV(L) + VSGZ*VSGZ ) )*VSGZ  
      ENDDO
    ENDDO  
    !$OMP END PARALLEL DO

    TTBXY=TTBXY+(DSTIME(0)-TTDS)

    !**********************************************************************C
    !
    ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
    !
    !----------------------------------------------------------------------C
    TTDS=DSTIME(0)
    IF( ISWAVE == 0 )THEN

      !*** STANDARD CALCULATIONS - NO CORNER CORRECTS
      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,TMP)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)

        DO LP=LF,LL
          L=LWET(LP)  
          TVAR3W(L)=TSX(LEC(L))
          TVAR3S(L)=TSY(LNC(L))
          TVAR3E(L)=TBX(LEC(L))
          TVAR3N(L)=TBY(LNC(L))
        ENDDO

        DO LP=LF,LL
          L=LWET(LP)
          TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2  
          QQ(L,0)  = 0.5*CTURB2*SQRT(TMP)
          !QQ(L,0)  = 0.5*CTURB2*SQRT( (TVAR3E(L)+TBX(L))**2 + (TVAR3N(L)+TBY(L))**2 )

          TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2  
          QQ(L,KC) = 0.5*CTURB2*SQRT(TMP) 
          !QQ(L,KC) = 0.5*CTURB2*SQRT( (TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2 )
           
          QQSQR(L,0) = SQRT(QQ(L,0))
        ENDDO
      ENDDO
      !$OMP END PARALLEL DO

    ENDIF

    !**********************************************************************C
    !
    ! *** SET BOTTOM AND SURFACE TURBULENT INTENSITY SQUARED AT (N+1)
    !
    !----------------------------------------------------------------------C
    IF( ISWAVE >= 1 )THEN

      !$OMP PARALLEL DEFAULT(SHARED)
      !$OMP DO PRIVATE(ND,LF,LL,LP,L)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL
          L=LWET(LP)  
          TVAR3W(L)=TSX(LEC(L))
          TVAR3S(L)=TSY(LNC(L))
          TVAR3E(L)=TBX(LEC(L))
          TVAR3N(L)=TBY(LNC(L))
        ENDDO
      ENDDO
      !$OMP END DO

      !$OMP DO PRIVATE(ND,LF,LL,LP,L,TAUBC,TAUBC2,UTMP,VTMP,CURANG,TAUB2,TMP)
      DO ND=1,NDM
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        DO LP=LF,LL
          L=LWET(LP)  
          IF( LWVMASK(L) )THEN
            TAUBC2 = (RSSBCE(L)*TVAR3E(L)+RSSBCW(L)*TBX(L))**2  +(RSSBCN(L)*TVAR3N(L)+RSSBCS(L)*TBY(L))**2  
            TAUBC=0.5*SQRT(TAUBC2)   ! *** CURRENT ONLY
            CTAUC(L)=TAUBC 
            UTMP=0.5*STCUV(L)*( U(LEC(L),KSZU(LEC(L))) + U(L,KSZU(L)) )+1.E-12  
            VTMP=0.5*STCUV(L)*( V(LN ,KSZV(LN) ) + V(L,KSZV(L)) )  
            CURANG=ATAN2(VTMP,UTMP)  
            TAUB2=TAUBC*TAUBC +      (QQWV3(L)*QQWV3(L)) + 2.     *TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)  
            TAUB2=MAX(TAUB2,0.)              ! *** CURRENT & WAVE      
            QQ(L,0 )   = CTURB2*SQRT(TAUB2)  ! *** CELL CENTERED TURBULENT INTENSITY DUE TO CURRENTS & WAVES
          
            QQ(L,KC)   = 0.5*CTURB2*SQRT((TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2)  
            QQSQR(L,0) = SQRT(QQ(L,0))
          ELSE
            TMP = ( RSSBCE(L)*TVAR3E(L) + RSSBCW(L)*TBX(L) )**2 + ( RSSBCN(L)*TVAR3N(L) + RSSBCS(L)*TBY(L) )**2  
            QQ(L,0 ) = 0.5*CTURB2*SQRT(TMP)
            !QQ(L,0 ) = 0.5*CTURB2*SQRT( (TVAR3E(L)+TBX(L))**2 + (TVAR3N(L)+TBY(L))**2 )

            TMP = ( RSSBCE(L)*TVAR3W(L) + RSSBCW(L)*TSX(L) )**2 + ( RSSBCN(L)*TVAR3S(L) + RSSBCS(L)*TSY(L) )**2  
            QQ(L,KC) = 0.5*CTURB2*SQRT(TMP) 
            !QQ(L,KC) = 0.5*CTURB2*SQRT( (TVAR3W(L)+TSX(L))**2 + (TVAR3S(L)+TSY(L))**2 )
           
            QQSQR(L,0)=SQRT(QQ(L,0))
          ENDIF
        ENDDO
      ENDDO
      !$OMP END DO
      !$OMP END PARALLEL
    ENDIF

    !**********************************************************************C
    !
    ! *** CALCULATE TURBULENT INTENSITY SQUARED
    !
    IF( KC > 1 )THEN
      CALL CALQQ1 (ISTL)
    ENDIF
    TQQQ=TQQQ+(DSTIME(0)-TTDS)

    !**********************************************************************C
    ! *** CALCULATE MEAN MASS TRANSPORT FIELD
    IF( ISSSMMT /= 2 )THEN
      IF( ISICM >= 1 )THEN
        NTMP=MOD(N,2)
        IF( ISTL == 3 .AND. NTMP == 0 ) CALL CALMMT
      ENDIF
    ENDIF

    !**********************************************************************C
    !
    ! *** HYDRODYNAMIC CALCULATIONS FOR THIS TIME STEP ARE COMPLETED
    ! *** IF NCTBC EQ NTSTBC APPLY TRAPEZOIDAL CORRECTION
    !
    !----------------------------------------------------------------------C
    IF( NCTBC == NTSTBC )THEN
      NCTBC=0
      ISTL=2
      DELT=DT          ! *** NOT USED IN HDMT BUT INDICATIVE
      DELTD2=0.5*DT
      DZDDELT=DZ/DELT
      ROLD=0.5
      RNEW=0.5
      GOTO 500
    ELSE
      NCTBC=NCTBC+1
      ISTL=3
      DELT=DT2         ! *** NOT USED IN HDMT BUT INDICATIVE
      DELTD2=DT
      DZDDELT=DZ/DELT
      ROLD=0.
      RNEW=1.
    ENDIF

    !**********************************************************************C
    !
    ! *** WRITE TO TIME SERIES FILES
    !
    CTIM=DT*FLOAT(N)+TCON*TBEGIN
    CTIM=CTIM/TCON

    ICALLTP=0
    IF( ISTMSR >= 1 )THEN
      IF( NCTMSR >= NWTMSR )THEN
        CALL TMSR
        ICALLTP=1
        NCTMSR=1
      ELSE
        NCTMSR=NCTMSR+1
      ENDIF
    ENDIF

    !**********************************************************************C
    ! *** OUTPUT ZERO DIMENSION VOLUME BALANCE
    IF( ISDRY >= 1 .AND. ICALLTP == 1 .AND. DEBUG )THEN
      OPEN(1,FILE=OUTDIR//'ZVOLBAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')
      DO LS=1,LORMAX
        IF( VOLZERD >= VOLSEL(LS) .AND. VOLZERD < VOLSEL(LS+1) )THEN
          WTM=VOLSEL(LS+1)-VOLZERD
          WTMP=VOLZERD-VOLSEL(LS)
          DELVOL=VOLSEL(LS+1)-VOLSEL(LS)
          WTM=WTM/DELVOL
          WTMP=WTMP/DELVOL
          SELZERD=WTM*BELSURF(LS)+WTMP*BELSURF(LS+1)
          ASFZERD=WTM*ASURFEL(LS)+WTMP*ASURFEL(LS+1)
        ENDIF
      ENDDO
      CTIM=(DT*FLOAT(N)+TCON*TBEGIN)/TCTMSR
      WRITE(1,5304)CTIM,SELZERD,ASFZERD,VOLZERD,VETZERD
      CLOSE(1)
     5304 FORMAT(2X,F10.4,2X,F10.5,3(2X,E12.4))
    ENDIF
    ICALLTP=0


    !**********************************************************************C
    ! *** CALCULATE MEAN MASS TRANSPORT FIELD
    IF( ISSSMMT /= 2 .AND. NTSMMT > 0 )THEN
      IF( ISICM == 0 ) CALL CALMMT
    ENDIF

    !**********************************************************************C
    !
    ! *** ADVANCE NEUTRALLY BUOYANT PARTICLE DRIFTER TRAJECTORIES
    !
    IF( ISPD >= 2 )THEN
      IF( TIMEDAY >= LA_BEGTI0 .AND. TIMEDAY <= LA_ENDTI0 )THEN
        TTDS=DSTIME(0)
        CALL DRIFTERC
        TLRPD=TLRPD+(DSTIME(0)-TTDS)
      ENDIF
    ENDIF

    !**********************************************************************C
    !
    ! *** CALCULATE VOLUME MASS, MOMENTUM AND ENERGY BALANCES
    !
    IF( ISBAL >= 1 )THEN
      CALL CALBAL5
      NTMP=MOD(N,2)
      IF( NTMP == 0 )THEN
        CALL CBALEV5
      ELSE
        CALL CBALOD5
      ENDIF
    ENDIF

    !   SEDIMENT BUDGET CALCULATION     (DLK 10/15)
    IF( ISSBAL >= 1 )THEN
      CALL BUDGET5
    ENDIF

    !**********************************************************************C
    ! *** PERFORM AN M2 TIDE HARMONIC ANALYSIS EVERY 2 M2 PERIODS
    IF( ISHTA == 1 ) CALL CALHTA

    !**********************************************************************C
    ! *** CALCULATE DISPERSION COEFFICIENTS
    IF( N >= NDISP .AND. NCTBC == 1 )THEN
      IF( ISDISP == 2 ) CALL CALDISP2
      IF( ISDISP == 3 ) CALL CALDISP3
    ENDIF

    !**********************************************************************C
    ! *** PERFORM LEAST SQUARES HARMONIC ANALYSIS AT SELECTED LOCATIONS
    IF( ISLSHA == 1 .AND. TIMESEC >= TIMEHARM )THEN  
      CALL LSQHARM  
      TIMEHARM = TIMEHARM + 300.
    ENDIF  


  !**********************************************************************!  
    ! *** WRITE TO TIME VARYING GRAPHICS FILES
    IF( HFREOUT == 1 )THEN
      DO NS=1,NSUBSET
        DEL = ABS(84600*(TIMEDAY-HFREDAY(NS)))
        IF( DEL <= DELT .AND. TIMEDAY <= HFREDAYEN(NS) )THEN  
          CALL HFREHYOUT(0,NS)
          CALL HFREWCOUT(0,NS)
          CALL HFREWQOUT(0,NS)
          CALL HFRERPEMOUT(0,NS)
          HFREDAY(NS) = HFREDAY(NS) + HFREMIN(NS)/1440
        ENDIF
      ENDDO
    ENDIF

#ifdef NCOUT
    ! ** NETCDF OUTPUT 
    IF (NCDFOUT > 0) THEN
      IF( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) .AND. (TIMEDAY >= TBEGNCDF .AND. TIMEDAY <= TENDNCDF )) THEN
        CALL NETCDF_WRITE(NC_ID)    
      ELSEIF (TIMEDAY > TENDNCDF ) THEN
        CALL NC_CLOSE_FILE(NC_ID)
      ENDIF
    ENDIF
#endif  		
    !**********************************************************************C
    ! *** WRITE EFDC EXPLORER FORMAT OUTPUT
    IF( ISPPH == 1 )THEN
      IF( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )THEN
        CALL EE_LINKAGE(0)
      ENDIF
    ENDIF
    IF( TIMEDAY >= SNAPSHOTS(NSNAPSHOTS) )THEN
      NSNAPSHOTS=NSNAPSHOTS+1
    ENDIF

    !**********************************************************************C
    ! *** WRITE TO TIME VARYING 3D HDF GRAPHICS FILES
    IF( N == NC3DO .AND. IS3DO == 1 )THEN
      CALL OUT3D
      NC3DO=NC3DO+(NTSPTC/NP3DO)
    ENDIF

    !**********************************************************************C
    ! *** WRITE RESTART FILE EVERY ISRESTO M2 REFERENCE PERIODS 
    IF( ISRESTO >= 1 )THEN  
      IF( TIMEDAY >= TIMELAST )THEN  
        CALL RESTOUT(0)  
        IF( ISTRAN(8) >= 1 )THEN  
          IF( IWQRST == 1 ) CALL WQWCRST  
          IF( IWQBEN == 1 ) CALL WQSDRST  
          IF( ISRPEM > 0  ) CALL WQRPEMRST  
        ENDIF
        TIMELAST = TIMEDAY + TIMERST
      ENDIF  
    ENDIF  

    !**********************************************************************C
    ! *** RECORD TIME

    ! *** DTIME AND FLUSH ARE SUPPORTED ON SUN SYSTEMS, BUT MAY NOT BE
    ! *** SUPPORTED ON OTHER SYSTEMS.
    !
    IF( NTIMER == NTSPTC )THEN
      ! *** EE BEGIN BLOCK
      CALL TIMELOG(N,TIMEDAY,OUTDIR)
      ! *** EE END BLOCK
      NTIMER=1
    ELSE
      NTIMER=NTIMER+1
    ENDIF

    !**********************************************************************C
    IF( ISHOW > 0 ) CALL SHOWVAL

    !**********************************************************************C
#ifdef _WIN  
    IF( KEY_PRESSED() )THEN
      IF( ISEXIT())GOTO 1001
    ENDIF
#endif

  1000 CONTINUE  

  !**********************************************************************C
  !**********************************************************************C
  !
  ! *** TIME LOOP COMPLETED
  !
  1001 THDMT=THDMT+(DSTIME(0)-T1TMP)

  !**********************************************************************C
  !
  2000 CONTINUE

  !**********************************************************************C

  ! *** WRITE RESTART FILE
  IF( ABS(ISRESTO) > 0 .AND. TIMEDAY > TIMELAST-TIMERST*0.5 )THEN
    WRITE(6,'(A,F12.4)') 'FINAL RESTART FILE: ',TIMEDAY
    CALL RESTOUT(0)  
    IF( ISTRAN(8) >= 1 )THEN  
      IF( IWQRST == 1 ) CALL WQWCRST  
      IF( IWQBEN == 1 ) CALL WQSDRST  
      IF( ISRPEM > 0  ) CALL WQRPEMRST  
    ENDIF
  ENDIF
  IF( ISRESTO == -2 )THEN
    CALL RESTMOD
  ENDIF

  !**********************************************************************C
  ! *** COMPLETE LEAST SQUARES HARMONIC ANALYSIS
  LSLSHA=1
  IF( ISLSHA == 1 ) CALL LSQHARM

  !**********************************************************************C
  ! *** OUTPUT COSMETIC VOLUME LOSSES FORM DRY CELLS  
  IF( NDRYSTP < 0 )THEN  
    OPEN(1,FILE=OUTDIR//'DRYLOSS.OUT')  
    CLOSE(1,STATUS='DELETE')  
    OPEN(1,FILE=OUTDIR//'DRYLOSS.OUT')  

    WRITE(1,'(3A5,A14,A10)')'L','I','J','VOL WASTED','EQUIV HP'
    DO L=2,LA  
      WRITE(1,1993)L,IL(L),JL(L),VDWASTE(L),VDWASTE(L)/DXYP(L)
    ENDDO  
    1993 FORMAT(3I5,E14.6,F10.3)  

    CLOSE(1)  
  ENDIF  

  !**********************************************************************C
  ! *** OUTPUT COURANT NUMBER DIAGNOSTICS
  IF( ISINWV == 1 .AND. DEBUG )THEN
    OPEN(1,FILE=OUTDIR//'CFLMAX.OUT')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'CFLMAX.OUT')

    DO L=2,LA
      WRITE(1,1991)IL(L),JL(L),(CFLUUU(L,K),K=1,KC)
      WRITE(1,1992)(CFLVVV(L,K),K=1,KC)
      WRITE(1,1992)(CFLWWW(L,K),K=1,KC)
      WRITE(1,1992)(CFLCAC(L,K),K=1,KC)
    ENDDO

    CLOSE(1)
     1991 FORMAT(2I5,12F7.2)
     1992 FORMAT(10X,12F7.2)
  ENDIF

  !**********************************************************************C
  !
  ! *** OUTPUT FINAL FOOD CHAIN AVERAGING PERIOD  
  IF( ISTRAN(5) >= 1 .AND. ISFDCH >= 1 )CALL FOODCHAIN(1)  
#ifdef NCOUT
  IF (NCDFOUT > 0) CALL NC_CLOSE_FILE(NC_ID)  
#endif  
END
