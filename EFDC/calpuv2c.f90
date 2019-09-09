SUBROUTINE CALPUV2C 

  !**********************************************************************C
  !
  ! ** SUBROUTINE CALPUV2C CALCULATES THE EXTERNAL SOLUTION FOR P, UHDYE,
  ! ** AND VHDXE, FOR FREE SURFACE FLOWS FOR THE 2TL SOLUTION.
  ! ** WITH PROVISIONS FOR WETTING AND DRYING OF CELLS
  !
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP
  !    2002-02-28    John Hamrick      Modified Drying And Wetting Scheme.
  !                                      The Old Formulation Remains  
  !       See (ISDRY > 0 .AND. ISDRY < 98). The New Formulation Is Activated  
  !       By ( ISDRY == 99 ). Also Added Option To Waste Water From Essentially  
  !       Dry Cells Having Water Depths Greater Than HDRY.  Ie The High And  
  !       Wet Cells Blocked By Dry Cells. This Is Actived By A Negative Value  
  !       Of NDRYSTP Parameter Is The Efdc.Inp File  
  !       Added Save Of Old Values Of Horizontal Flow Face Switches SUB1 & SVB1  
  !       And Transport Bypass Mask, LMASKDRY For Dry Cells.
  !       Added QDWASTE(L) To Save Source Equivalent Of Volume Loss Rate  
  !       For Reducing Depth Of High/Dry Cells.  Also Added Concentration  
  !       Adjustment  

  USE GLOBAL 
  USE EFDCOUT
  
  IMPLICIT NONE
  
  INTEGER :: NMD,ITERHP,NCORDRY,ICORDRY,ND,LF,LL,L,LP,LS,LN,LW,LE,LHOST,LCHNU,LCHNV
  INTEGER :: IUW,IUE,IVS,IVN,IFACE,LMAX,LMIN,NTMP,IMAX,IMIN,JMAX,JMIN,I,K,NNEGFLG
  INTEGER,SAVE :: NOPTIMAL,LDMOPT,INOTICE,NCORDRYMAX,NCORDRYAVG,ITERMAX,ITERAVG,NITERAVG,NCOUNT,NMLWC
  
  REAL :: DELTD2,RLAMN,RLAMO,TMPX,TMPY,C1,TMPVAL,HDRY2,BELVAVG,SVPW1
  REAL :: SUBW,SUBE,SVBS,SVBN,DHPDT,DHPDT2,RDRY,CCMNM,CCMNMI,HDRY90
  REAL :: RVAL,HOLDTMP,RNPORI,DIVEXMX,DIVEXMN,DIVEX,ETGWTMP,ETGWAVL
  REAL(RKD),SAVE :: DAYOLD,HPMLWC
  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: IACTIVE  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: ICORDRYD
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: NNATDRY
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:,:) :: LNATDRY

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: QCHANUT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: QCHANVT  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: QSUMTMP 
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SUB1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SVB1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: CCMNMD

  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: FSGZUDXYPI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: FSGZVDXYPI
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: HPOLD

  IF( .NOT. ALLOCATED(IACTIVE) )THEN
    ! *** SET THE OPTIMAL NUMBER OF THREADS.  USE 100 CELLS PER THREAD AS GENERAL RULE
    NOPTIMAL = MIN(NTHREADS,8,MAX(LA/100,1))
    LDMOPT   = INT(FLOAT(LA-1)/FLOAT(NOPTIMAL))+1
    WRITE(*,'(A,I5)')'FIRST CALL TO 2TL PRESSURE SOLUTION.  CALPUV THREADS:',NOPTIMAL

    ALLOCATE(IACTIVE(NCHANM))  
    ALLOCATE(QCHANUT(NCHANM))  
    ALLOCATE(QCHANVT(NCHANM))  
    ALLOCATE(QSUMTMP(LCM))  
    ALLOCATE(SUB1(0:LCM))  
    ALLOCATE(SVB1(0:LCM))
    ALLOCATE(CCMNMD(NOPTIMAL))
    ALLOCATE(ICORDRYD(NOPTIMAL))
    ALLOCATE(FSGZUDXYPI(LCM))
    ALLOCATE(FSGZVDXYPI(LCM))
    ALLOCATE(HPOLD(LCM))
    ALLOCATE(LNATDRY(NOPTIMAL,(INT(LCM/NOPTIMAL)+1)))
    ALLOCATE(NNATDRY(NOPTIMAL))
    
    IACTIVE=0
    QCHANUT=0.
    QCHANVT=0.
    QSUMTMP=0.
    SUB1=SUB
    SVB1=SVB
    ISCDRY=0
    CCMNMD=0
    ICORDRYD=0
    RCX=1.0
    RCY=1.0
    RCX(1)=0.  
    RCY(1)=0.  
    RCX(LC)=0.  
    RCY(LC)=0.  
    LNATDRY=0
    NNATDRY=0
    INOTICE=-999
    NCORDRYMAX=0
    NCORDRYAVG=0
    ITERMAX=0
    ITERAVG=0
    NITERAVG=0
    NCOUNT=0
    DAYOLD=INT(TIMEDAY)
    NMLWC = 0
    HPMLWC = 0.
    
    ! INITIALIZE DIAGONAL
    CC=1.0

    FSGZUDXYPI=1.
    FSGZVDXYPI=1.
    DO L=2,LA
      LS=LSC(L)
      LW=LWC(L)
      FSGZUDXYPI(L) =  0.5/DXYU(L)
      FSGZVDXYPI(L) =  0.5/DXYV(L)
    ENDDO
    DO L=2,LA
      HPOLD(L) = HP(L)
    ENDDO
    
    IF( ISDRY > 0 .AND. NTCNB > 98 )THEN
      PRINT '(A)','MCCLELLAND LAKE WETLAND COMPLEX HARDWIRE!'
    ENDIF

  ENDIF

  IF( ISDYNSTP == 0 )THEN  
    DELT=DT  
    DELTD2=0.5*DT  
    DELTI=1./DELT  
  ELSE  
    DELT=DTDYN  
    DELTD2=0.5*DTDYN  
    DELTI=1./DELT  
  ENDIF  
  ISTL=2  

  ! *** INITIALIZE SUBGRID SCALE CHANNEL INTERACTIONS  
  IF( MDCHH >= 1 )THEN  
    DO NMD=1,MDCHH  
      QCHANUT(NMD)=QCHANU(NMD)  
      QCHANVT(NMD)=QCHANV(NMD)  
    ENDDO  
  ENDIF  

  ! *** SET SWITCHES FOR DRYING AND WETTING  
  ITERHP=0  
  NCORDRY=0  
  NNEGFLG = 0

  IF( ISDRY > 0 )THEN
    DO LL=1,NBCSOP
      L=LOBCS(LL)
      LOPENBCDRY(L) = .FALSE.
    ENDDO
  ENDIF

  ! ***************************************************************************
  ! *** CALCULATE EXTERNAL BUOYANCY INTEGRALS AT TIME LEVEL (N) 
  IF( BSC > 1.E-6 )CALL CALEBI  

  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LE,LS,LN,LW,K,TMPX,TMPY,C1)
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)

    ! ***
    
    ! *** SET SWITCHES FOR DRYING AND WETTING  
    IF( ISDRY > 0 )THEN
      DO L=LF,LL  
        ISCDRY(L)=0
        SUB1(L)=SUB(L)
        SVB1(L)=SVB(L)
        OLDMASK(L)=LMASKDRY(L)
      ENDDO
  
      ! *** ZERO VOLUME WASTING ARRAY COUNTERS
      IF( NDRYSTP < 0 )THEN
        NNATDRY(ND)=0
        IF( ISBAL > 0 )THEN
          DO L=LF,LL  
            QDWASTE(L)=0.
          ENDDO
        ENDIF
      ENDIF

    ENDIF

    IF( BSC > 1.E-6 )THEN
      ! *** CALCULATE EXPLICIT EXTERNAL DENSITY GRADIENTS  
      IF( IGRIDV == 0 )THEN
        DO L=LF,LL  
          LW=LWC(L) 
          FPGXE(L) = -SBX(L)*HU(L)*GP*((BI2W(L)+BI2W(LW))*(HP(L)-HP(LW))+2.0*HU(L)*(BI1W(L)-BI1W(LW))+(BEW(L)+BEW(LW))*(BELV(L)-BELV(LW)))  

          LS=LSC(L)  
          FPGYE(L) = -SBY(L)*HV(L)*GP*((BI2S(L)+BI2S(LS))*(HP(L)-HP(LS))+2.0*HV(L)*(BI1S(L)-BI1S(LS))+(BES(L)+BES(LS))*(BELV(L)-BELV(LS)))  
        ENDDO  
      ELSE
        DO L=LF,LL  
          LW=LWC(L) 
          FPGXE(L) = -SBX(L)*HU(L)*GP*( (BI2W(L)+BI2E(LW))*(HPW(L)-HPE(LW)) + 2.0*HU(L)*(BI1W(L)-BI1E(LW)) + (BEW(L)+BEE(LW))*(BELVW(L)-BELVE(LW)) )  
          
          LS=LSC(L)  
          FPGYE(L) = -SBY(L)*HV(L)*GP*( (BI2S(L)+BI2N(LS))*(HPS(L)-HPN(LS)) + 2.0*HV(L)*(BI1S(L)-BI1N(LS)) + (BES(L)+BEN(LS))*(BELVS(L)-BELVN(LS)) )  
        ENDDO  
      ENDIF
    ENDIF
    
    ! *** CALCULATE EXPLICIT EXTERNAL UHDYE AND VHDXE EQUATION TERMS  
    DO L=LF,LL
      LS=LSC(L)
      LW=LWC(L) 
      FUHDYE(L) = UHDYE(L) - DELTD2*SUB(L)*HRUO(L)*HU(L)*(P(L)-P(LW)) + SUB(L)*DELT*DXIU(L)*(DXYU(L)*(TSX(L)-RITB1*TBX(L)) + FCAXE(L) + FPGXE(L) - SNLT*FXE(L))  

      FVHDXE(L) = VHDXE(L) - DELTD2*SVB(L)*HRVO(L)*HV(L)*(P(L)-P(LS)) + SVB(L)*DELT*DYIV(L)*(DXYV(L)*(TSY(L)-RITB1*TBY(L)) - FCAYE(L) + FPGYE(L) - SNLT*FYE(L)) 
    ENDDO  

    ! *** SET IMPLICIT BOTTOM AND VEGETATION DRAG AS APPROPRIATE  
    IF( ISITB >= 1 )THEN  
      ! *** IMPLICIT BOTTOM DRAG WITH VEGETATION  
      DO L=LF,LL
        TMPX=1.0
        TMPY=1.0
        IF( UHE(L) /= 0.0) TMPX = U(L,KSZU(L))*HU(L)/UHE(L)
        IF( VHE(L) /= 0.0) TMPY = V(L,KSZV(L))*HV(L)/VHE(L)
        RCX(L) = 1./( 1.+TMPX*RITB*DELT*HUI(L)*STBX(L)*SQRT(VU(L)*VU(L)+U(L,KSZU(L))*U(L,KSZU(L)))+DELT*FXVEGE(L) )
        RCY(L) = 1./( 1.+TMPY*RITB*DELT*HVI(L)*STBY(L)*SQRT(UV(L)*UV(L)+V(L,KSZV(L))*V(L,KSZV(L)))+DELT*FYVEGE(L) )
        FUHDYE(L) = FUHDYE(L)*RCX(L)
        FVHDXE(L) = FVHDXE(L)*RCY(L)
      ENDDO
    ELSEIF( ISVEG > 0 )THEN
      ! *** IMPLICIT VEGETATION DRAG ONLY.  REDUCE BY THE TOTAL OF ENERGY
      DO L=LF,LL
        RCX(L) = 1./( 1.+DELT*FXVEGE(L) )
        RCY(L) = 1./( 1.+DELT*FYVEGE(L) )
        FUHDYE(L) = FUHDYE(L)*RCX(L)
        FVHDXE(L) = FVHDXE(L)*RCY(L)
      ENDDO
    ENDIF

    ! *** RESET BOUNDARY CONDITIONS SWITCHES  
    IF( ISDRY > 0 )THEN  
      DO L=LF,LL  
        SUB(L) = SUBO(L)  
        SVB(L) = SVBO(L)  
        SBX(L) = SBXO(L)  
        SBY(L) = SBYO(L)  
      ENDDO  
    ENDIF

    ! *** ADVANCE EXTERNAL VARIABLES  
    DO L=LF,LL  
      UHDY1E(L)=UHDYE(L)  
      VHDX1E(L)=VHDXE(L)  
      P1(L)=P(L)  
      H1U(L)=HU(L)  
      H1V(L)=HV(L)  
      H1UI(L)=HUI(L)  
      H1VI(L)=HVI(L)  
      H2P(L)=H1P(L)
      H1P(L)=HP(L)
    ENDDO

    DO K=1,KC  
      DO L=LF,LL
        H1PK(L,K) = HPK(L,K)   
        UHDY1EK(L,K) = UHDYEK(L,K)
        VHDX1EK(L,K) = VHDXEK(L,K) 
      ENDDO
    ENDDO

    IF( ISGWIE >= 1 )THEN
      DO L=LF,LL  
        AGWELV2(L)=AGWELV1(L)  
        AGWELV1(L)=AGWELV(L)  
      ENDDO
    ENDIF  

    ! *** SET OLD TIME LEVEL TERMS IN CONTINUITY EQUATION FOR NON BOUNDARY POINTS  
    ! *** DXYP = DXP*DYP (M^2),  G-9.82 (M/S^2), P (M2/S2), FP1 (M4/S3)
    C1=0.5*G
    DO L=LF,LL
      LE=LEC(L)
      LN=LNC(L)
      FP1(L) = DELTI*DXYP(L)*P(L) - C1*( UHDYE(LE)-UHDYE(L)+VHDXE(LN)-VHDXE(L) )
    ENDDO  

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO
  
  ! *************************************************************************
  ! ***  SET NEW TIME LEVEL TERMS IN CONTINUITY EQUATION INCLUDING
  ! ***  HOST-GUEST CHANNEL INTERACTION FOR NON BOUNDARY POINTS 

  !----------------------------------------------------------------------C
  ! ***  REENTER AT 1000 FOR WETTING-DRYING CORRECTION AND CHANNEL 
  ! ***  INTERACTION
  1000 CONTINUE  
  !----------------------------------------------------------------------C
  
  ! *** SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM  
  CCMNMD=1.E+18 
  
  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DEFAULT(SHARED)
  !$OMP DO PRIVATE(ND,LF,LL,L,LE,LS,LN,C1) 
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)

    C1=0.5*G
    DO L=LF,LL
      LE=LEC(L)
      LN=LNC(L)
      ! ***  USE THE SUB & SVB SWITCHES FOR MOMENTUM FOR THE CURRENT WET/DRY ITERATION
      FP(L) = FP1(L) - C1*( SUB(LE)*FUHDYE(LE) - SUB(L)*FUHDYE(L) + SVB(LN)*FVHDXE(LN) - SVB(L)*FVHDXE(L) - 2.0*QSUME(L) )
    ENDDO  

    IF( ISGWIE >= 1 )THEN  
      DO L=LF,LL  
        FP(L) = FP(L) - G*SPB(L)*(EVAPSW(L)-QGW(L))   ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
      ENDDO  
    ENDIF  

    C1=-0.5*DELTD2*G  
    DO L=LF,LL  
      LE=LEC(L)
      LN=LNC(L)  
      CS(L) = C1*SVB(L) *HRVO(L) *RCY(L) *HV(L)  
      CW(L) = C1*SUB(L) *HRUO(L) *RCX(L) *HU(L)  
      CE(L) = C1*SUB(LE)*HRUO(LE)*RCX(LE)*HU(LE)  
      CN(L) = C1*SVB(LN)*HRVO(LN)*RCY(LN)*HV(LN)  
    ENDDO  

    ! *** SET THE CENTER
    DO L=LF,LL  
      CC(L) = DELTI*DXYP(L) - CS(L)-CW(L)-CE(L)-CN(L)  
    ENDDO  

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP SINGLE
  
  ! *** APPLY THE OPEN BOUNDARY CONDITIONS
  IF( NBCSOP > 0 ) CALL SETOPENBC(DELTD2,DELTI,HU,HV)    

  ! *** INSERT IMPLICT SUB-GRID SCALE CHANNEL INTERACTIONS  
  IF( MDCHH >= 1 )THEN
    RLAMN=QCHERR  
    RLAMO=1.-RLAMN  
    CALL SUBCHAN(QCHANUT,QCHANVT,IACTIVE,DELT)  
  ENDIF
  !$OMP END SINGLE

  ! *** SCALE COEFFICIENTS IN EXTERNAL MODEL LINEAR EQUATION SYSTEM
  !$OMP DO PRIVATE(ND,LF,LL,L) 
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)
    DO L=LF,LL
      CCMNMD(ND)=MIN(CCMNMD(ND),CC(L))
      FPTMP(L)=FP(L)  
    ENDDO  
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  !$OMP SINGLE
  CCMNM =1.E+18 
  DO ND=1,NOPTIMAL 
    CCMNM=MIN(CCMNM,CCMNMD(ND))
  ENDDO
  CCMNMI=1./CCMNM  

  ! *** APPLY THE OPEN BOUNDARY CONDITIONS FOR ADJACENT CELLS
  IF( NBCSOP > 0 ) CALL SETOPENBC2   
  !$OMP END SINGLE

  ! *** SCALE BY MINIMUM DIAGONAL  
  ! *** BEGIN DOMAIN LOOP
  !$OMP DO PRIVATE(ND,LF,LL,L) 
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)
    DO L=LF,LL  
      CCS(L)=CS(L)*CCMNMI  
      CCW(L)=CW(L)*CCMNMI  
      CCE(L)=CE(L)*CCMNMI  
      CCN(L)=CN(L)*CCMNMI  
      CCC(L)=CC(L)*CCMNMI  
      FPTMP(L)=FPTMP(L)*CCMNMI  
      CCCI(L)=1./CCC(L)  
    ENDDO  
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP END PARALLEL

  IF( MDCHH >= 1 )THEN  
    DO NMD=1,MDCHH  
      CCCCHH(NMD)=CCCCHH(NMD)*CCMNMI  
    ENDDO  
  ENDIF  
  
  ! *************************************************************************
  ! *** CALL THE PRECONDITIONED CONJUGATE GRADIENT SOLVER
  IF( MDCHH == 0 ) CALL CONGRAD(NOPTIMAL,LDMOPT)
  IF( MDCHH >= 1 ) CALL CONGRADC  
  ! *************************************************************************
  ITERMAX = MAX(ITERMAX,IDRYTBP)
  NITERAVG = NITERAVG + 1
  ITERAVG = ITERAVG + IDRYTBP
  
  !$OMP PARALLEL DEFAULT(SHARED) 
  !$OMP DO PRIVATE(ND,LF,LL,L,LS,LW) 
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)
    
    ! *** CELL FACE DISCHARGE (M3/S)
    DO L=LF,LL
      LS=LSC(L)  
      LW=LWC(L) 
      UHDYE(L) = SUB(L)*( FUHDYE(L) - DELTD2*HRUO(L)*RCX(L)*HU(L)*(P(L)-P(LW)) )  
      VHDXE(L) = SVB(L)*( FVHDXE(L) - DELTD2*HRVO(L)*RCY(L)*HV(L)*(P(L)-P(LS)) )  
    ENDDO  
    
    ! *** UNIT WIDTH DISCHARGE AT CELL FACE (M2/S)  
    DO L=LF,LL  
      UHE(L) = UHDYE(L)*DYIU(L)  
      VHE(L) = VHDXE(L)*DXIV(L)  
    ENDDO  

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  ! *** CALCULATE NEW SUB-GRID SCALE CHANNEL EXCHANGE FLOWS  
  !$OMP SINGLE
  IF( MDCHH >= 1 )THEN  
    DO NMD=1,MDCHH  
      IF( IACTIVE(NMD) > 0 )THEN  
        LHOST=LMDCHH(NMD)  
        LCHNU=LMDCHU(NMD)  
        LCHNV=LMDCHV(NMD)  
        IF( MDCHTYP(NMD) == 1 )THEN  
          QCHANU(NMD)=CCCCHU(NMD)*QCHANUT(NMD)-RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNU))-RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNU))  
          QCHANV(NMD)=0.  
        ENDIF 
        IF( MDCHTYP(NMD) == 2 )THEN  
          QCHANU(NMD)=0.  
          QCHANV(NMD)=CCCCHU(NMD)*QCHANVT(NMD)-RLAMN*CCCCHU(NMD)*CCCCHV(NMD)*(P(LHOST)-P(LCHNV))-RLAMO*CCCCHU(NMD)*CCCCHV(NMD)*(P1(LHOST)-P1(LCHNV))  
        ENDIF  
      ELSE  
        QCHANV(NMD)=0.  
        QCHANVN(NMD)=0.  
        QCHANU(NMD)=0.  
        QCHANUN(NMD)=0.  
      ENDIF  
    ENDDO  
  ENDIF  
  !$OMP END SINGLE
  
  ! **  CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL TRANSPORTS AT (N+1)
  !$OMP DO PRIVATE(ND,LF,LL,L,LE,LS,LN) 
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)

    ! *** CALCULATE REVISED CELL DEPTHS BASED ON NEW HORIZONTAL TRANSPORTS AT (N+1)  
    DO L=LF,LL  
      LE=LEC(L)
      LN=LNC(L)  
      HP(L) = H1P(L) + DELTD2*DXYIP(L)*( 2.*QSUME(L) - ( UHDYE(LE)+UHDY1E(LE)-UHDYE(L)-UHDY1E(L)  &
                                                       + VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L) ) )  
    ENDDO  

    IF( ISGWIE >= 1 )THEN  
      DO L=LF,LL  
        HP(L) = HP(L) - DELT*DXYIP(L)*(EVAPSW(L)-QGW(L))   ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
      ENDDO  
    ENDIF  
    
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  
  ! *** APPLY OPEN BOUNDARYS 
  !$OMP SINGLE
  DO LL=1,NBCSOP
    L=LOBCS(LL)
    HP(L)=GI*P(L)-BELV(L)  

    ! *** CHECK FOR CONDITION OF PSERT<BELV
    IF( ISDRY > 0 )THEN
      
      ! ***                 ** HP CHECKING IS FOR RADIATION BC'S **
      IF( (LOPENBCDRY(L) .OR. (H1P(L) < 0.9*HDRY .AND. HP(L) < H1P(L)) .OR. HP(L) < 0. ) .AND. ISCDRY(L) == 0 )THEN
        IF( HP(L) <= 0. .AND. LOPENBCDRY(L) )THEN
          PRINT '(A,I5,3F10.4,F12.4)',' WARNING!  NEG DEPTH AT OPEN BC: L,HP,H1P,H2P,TIMEDAY',L,HP(L),H1P(L),H2P(L),TIMEDAY
          OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  
          WRITE(8,'(A,I5,3F10.4,F12.4)')' WARNING!  NEG DEPTH AT OPEN BC: L,HP,H1P,H2P,TIMEDAY',L,HP(L),H1P(L),H2P(L),TIMEDAY
          CLOSE(8)
          HP(L) = 0.2*HDRY
        ELSE
          HP(L) = MIN(MAX(H1P(L),0.1*HDRY),0.9*HDRY)
        ENDIF
        
        ISCDRY(L)=1  
        ICORDRY=-999  

        LE=LEC(L)
        LN=LNC(L)
        
        SUB(L)=0.
        SVB(L)=0.
        SUB(LE)=0.
        SVB(LN)=0.
        SBX(L)=0.  
        SBY(L)=0.  
        SBX(LE)=0.  
        SBY(LN)=0.  
        UHDYE(L)=0.
        UHDYE(LE)=0.
        VHDXE(L)=0.
        VHDXE(LN)=0.
        
        LOPENBCDRY(L) = .TRUE.
        CC(L)  = DELTI*DXYP(L)
        P(L)   = (HP(L)+BELV(L))*G
        FP1(L) = DELTI*DXYP(L)*P(L)
      ENDIF
    ENDIF

  ENDDO 
  
  ! *** ADD CHANNEL INTERACTION EXCHANGES  
  IF( MDCHH >= 1 )THEN  
    DO NMD=1,MDCHH  
      IF( IACTIVE(NMD) > 0 )THEN  
        LHOST=LMDCHH(NMD)  
        LCHNU=LMDCHU(NMD)  
        LCHNV=LMDCHV(NMD)  
        IF( MDCHTYP(NMD) == 1 )THEN  
          TMPVAL=DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUT(NMD))  
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)  
          HP(LCHNU)=HP(LCHNU)-TMPVAL*DXYIP(LCHNU)  
        ENDIF  
        IF( MDCHTYP(NMD) == 2 )THEN  
          TMPVAL=DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVT(NMD))  
          HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)  
          HP(LCHNV)=HP(LCHNV)-TMPVAL*DXYIP(LCHNV)  
        ENDIF  
      ENDIF  
    ENDDO  
  ENDIF  
  !$OMP END SINGLE

  ! *** PERFORM INTERMEDIATE UPDATES OF P
  IF( ISDRY > 0 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,L) 
    DO ND=1,NOPTIMAL  
      LF=2+(ND-1)*LDMOPT  
      LL=MIN(LF+LDMOPT-1,LA)
      DO L=LF,LL  
        P(L) = G*(HP(L)+BELV(L))  
      ENDDO  
    ENDDO
    !$OMP END DO
  ENDIF
  !$OMP END PARALLEL
  
  ! *** CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY  
  IF( ISDRY > 0 .AND. ISDRY < 98 )THEN  
    ICORDRY=0
    ICORDRYD=0  
    DO L=2,LA  
      LE=LEC(L)  
      LN=LNC(L)  
      IF( HP(L) <= HDRY )THEN  
        IF( ISCDRY(L) == 0 )THEN  
          ISCDRY(L)=1  
          ICORDRY=1  
        ENDIF  
        SUB(L)=0.  
        SVB(L)=0.  
        SUB(LE)=0.  
        SVB(LN)=0.  
        SBX(L)=0.  
        SBY(L)=0.  
        SBX(LE)=0.  
        SBY(LN)=0.  
      ENDIF  
    ENDDO  
    IF( ICORDRY == 1 )THEN  
      NCORDRY=NCORDRY+1  
      IF (NCORDRY < 500) THEN
        GOTO 1000
      ELSE
        STOP '*** NCORDRY > 500 WHEN 0 < ISDRY < 98'
      ENDIF
    ENDIF 
  ENDIF  

  ! *** CHECK FOR DRYING AND RESOLVE EQUATIONS IF NECESSARY  
  IF( ISDRY == 99 )THEN  
    HDRY2=2.*HDRY  
    HDRY90 = 0.9*HDRY
    ICORDRY=0
    ICORDRYD=0  
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LS,LN,LW,LE)  &
    !$OMP             PRIVATE(SUBW,SUBE,SVBS,SVBN,DHPDT,DHPDT2,RDRY,TMPVAL)
    DO ND=1,NOPTIMAL  
      LF=2+(ND-1)*LDMOPT  
      LL=MIN(LF+LDMOPT-1,LA)

      DO L=LF,LL
        LW=LWC(L) 
        LE=LEC(L) 
        LS=LSC(L)  
        LN=LNC(L)  
          
        IF( HP(L) <= HDRY )THEN
          SUBW=SUB(L)  
          SUBE=SUB(LE)  
          SVBS=SVB(L)  
          SVBN=SVB(LN)  
          DHPDT = ( HP(L) - H1P(L) )
          
          ! *** ALLOW RE-WETTING
          IF( DHPDT > 0.0 )THEN 
            SUB(L)=0.0  
            SUB(LE)=0.0  
            SVB(L)=0.0  
            SVB(LN)=0.0  
            SBX(L)=0.0  
            SBX(LE)=0.0  
            SBY(L)=0.0  
            SBY(LN)=0.0  
            
            ! *** RAISING WATER, IS IT FAST ENOUGH TO STAY WET
            DHPDT2 = (HDRY - H1P(L))*0.01
            IF( DHPDT > DHPDT2 )THEN
              IF( SUBO(L) > 0.5 )THEN  
                IF( UHDYE(L) > 0.0 .AND. HP(LW) > HDRY2 )THEN 
                  SUB(L)=SUBO(L)
                  SBX(L)=SBXO(L)
                ENDIF  
              ENDIF  
              IF( SUBO(LE) > 0.5 )THEN  
                IF( UHDYE(LE) < 0.0 .AND. HP(LE) > HDRY2 )THEN
                  SUB(LE)=SUBO(LE)
                  SBX(LE)=SBXO(LE)  
                ENDIF  
              ENDIF  
              IF( SVBO(L) > 0.5 )THEN  
                IF( VHDXE(L) > 0.0 .AND. HP(LS) > HDRY2 )THEN
                  SVB(L)=SVBO(L)
                  SBY(L)=SBYO(L)
                ENDIF  
              ENDIF  
              IF( SVBO(LN) > 0.5 )THEN  
                IF( VHDXE(LN) < 0.0 .AND. HP(LN) > HDRY2 )THEN
                  SVB(LN)=SVBO(LN)
                  SBY(LN)=SBYO(LN)
                ENDIF  
              ENDIF  
              RDRY=SUB(L)+SUB(LE)+SVB(L)+SVB(LN)  
              IF( RDRY < 0.5 )THEN  
                ISCDRY(L)=1  
              ELSE  
                ISCDRY(L)=0  
              ENDIF  
              TMPVAL=ABS(SUB(L)-SUBW)  
              IF( TMPVAL > 0.5 )THEN
                ICORDRYD(ND)=1  
              ELSE
                TMPVAL=ABS(SUB(LE)-SUBE)  
                IF( TMPVAL > 0.5 )THEN
                  ICORDRYD(ND)=1  
                ELSE
                  TMPVAL=ABS(SVB(L)-SVBS)  
                  IF( TMPVAL > 0.5 )THEN
                    ICORDRYD(ND)=1  
                  ELSE
                    TMPVAL=ABS(SVB(LN)-SVBN)  
                    IF( TMPVAL > 0.5)THEN
                      ICORDRYD(ND)=1  
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ELSE
              ! *** CASE: HP < HDRY BUT RISING, JUST NOT FAST ENOUGH
              IF( ISCDRY(L) == 0 )THEN
                ISCDRY(L)=1  
                ICORDRYD(ND)=1  
              ENDIF  
            ENDIF     ! *** END OF REWETTING SECTION, DHPDT > DHPDT2
            
          ELSEIF( HP(L) < HDRY90 .OR. H1P(L) < HDRY )THEN  
            ! *** HP < HDRY.  SET SWITCHES TO DRY
            SUB(L)=0.0
            SUB(LE)=0.0  
            SVB(L)=0.0  
            SVB(LN)=0.0  
            SBX(L)=0.0  
            SBX(LE)=0.0  
            SBY(L)=0.0  
            SBY(LN)=0.0  
            IF( ISCDRY(L) == 0 )THEN  
              ISCDRY(L)=1  
              ICORDRYD(ND)=1  
            ENDIF  
          ENDIF  
        ENDIF  
        
      ENDDO  
    
    ENDDO    ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
    
    ICORDRY=SUM(ICORDRYD(:))
    IF( ICORDRY > 0 )THEN  
      NCORDRY=NCORDRY+1  
      IF (NCORDRY < 500 ) THEN
        GOTO 1000
      ELSE
        ! *** WET/DRY MAXIMUM ITERATIONS EXCEEDED
        OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  

        ! *** SAVE A SNAPSHOT FOR EE
        IF( ISPPH == 1 )THEN
          WRITE (6,*)'THE LATEST MODEL RESULTS HAVE BEEN SAVED TO THE EE LINKAGE'
          WRITE(8,*) 'CALPUV: WRITING EE LINKAGE',TIMEDAY
          CALL EE_LINKAGE(-1)  
        ENDIF
        
        ! *** WRITE TO LOG THE CELLS THAT HAVE CHANGED
        DO L=2,LA
          ! ***    WAS WET BUT NOW DRY                           WAS DRY BUT NOW WET
          IF( (OLDMASK(L) .AND. ISCDRY(L) == 1) .OR. (.NOT. OLDMASK(L) .AND. ISCDRY(L) == 0) ) THEN
            WRITE(8,'(A,F10.3,I10,I5,2F10.4)') 'TIMEDAY,N,L,H1P(L),HP(L) = ',TIMEDAY,N,L,H1P(L),HP(L)
          ENDIF  
        ENDDO  
        CLOSE(8)
        
        STOP '*** NCORDRY > 500 WHEN ISDRY > 0'
      ENDIF
    ENDIF  

    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***

  ENDIF    ! *** END OF WET/DRY SECTION  ISDRY>0

  NCORDRYMAX = MAX(NCORDRY,NCORDRYMAX)
  NCORDRYAVG = NCORDRYAVG+NCORDRY
  NCOUNT = NCOUNT+1
  
  !**********************************************************************C
  ! *** FINISHED WITH WETTING/DRYING ITERATIONS

  IF( INT(TIMEDAY) /= DAYOLD .OR. TIMEDAY >= (TIMEEND-DELT/86400.) )THEN
    OPEN(888,FILE=OUTDIR//'CALPUV.LOG',POSITION='APPEND')  
    WRITE(888, '(I15,I10,F12.3,2(10X,2I10))') N,NITER,TIMEDAY,NINT(FLOAT(NCORDRYAVG)/FLOAT(NCOUNT)),NCORDRYMAX,NINT(FLOAT(ITERAVG)/FLOAT(NITERAVG)),ITERMAX
    CLOSE(888)
    NCORDRYAVG = 0
    NCORDRYAVG = 0
    ITERMAX = 0
    NITERAVG = 0
    ITERAVG = 0
    NCOUNT = 0
    DAYOLD = INT(TIMEDAY)
    
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
    ! ***
  ENDIF
  
  ! *** REPORT CELLS THAT HAVE JUST REWETTED
  IF( ISDRY > 0 .AND. DEBUG )THEN
    OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  
    DO L=2,LA
      IF( .NOT. LMASKDRY(L) )THEN 
        IF( HP(L) >= HDRY )THEN
          ! *** PREVIOUSLY CELL WAS DRY, NOW IT IS WET
          WRITE(8,'(A,2I5,F12.5,I5,L5,3F10.5)') 'REWETTED CELL: ',N,2,TIMEDAY,L,LMASKDRY(L),H1P(L),HP(L)
        ENDIF
      ENDIF
    ENDDO
    CLOSE(8)
  ENDIF
      
  ! *** BEGIN DOMAIN LOOP
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,L,LS,LN,LE,LW,NTMP,IUW,IUE,IVS,IVN,IFACE,K)  &
  !$OMP          PRIVATE(RDRY,BELVAVG,RVAL,HOLDTMP,TMPVAL,SVPW1,RNPORI,ETGWTMP,ETGWAVL)
  DO ND=1,NOPTIMAL  
    LF=2+(ND-1)*LDMOPT  
    LL=MIN(LF+LDMOPT-1,LA)

    ! *** COUNT THE NUMBER TO TIME STEPS A CELL IS ISOLATED, AND IF IT HAS BEEN
    ! *** ISOLATED FOR MORE THAN ABS(NDRYSTP), AND ITS BOTTOM ELEVATION IS HIGHER
    ! *** THAN THE SURROUNDING DRY CELLS, THEN REDUCE ITS DEPTH BELOW THE 
    ! *** DRYING DEPTH IF NECESSARY.  SAVE VOLUME REDUCTION RATE AS QDWASTE
    ! *** DEFINED AS POSITIVE OUT. 
    IF( ISDRY > 0 )THEN  
      IF( NDRYSTP < 0 )THEN
        NTMP=ABS(NDRYSTP)
        DO L=LF,LL
          IF( HP(L) >= HDRY )THEN
            ! *** WET CELL, DETERMINE IF ISOLATED
            LW=LWC(L)
            LE=LEC(L)
            LS=LSC(L)
            LN=LNC(L)
            RDRY = SUB(L) + SUB(LE) + SVB(L) + SVB(LN)
            IF( RDRY > 0.5 )THEN
              NATDRY(L)=0                         ! *** CELL IS NOT ISOLATED
            ELSE
              NATDRY(L) = NATDRY(L)+1             ! *** CELL IS ISOLATED
            
              IF( NATDRY(L) > NTMP )THEN
                ! *** EXCEEDED THE NUMBER OF ISOLATED STEPS SO DETERMINE IF NEED TO DRY OUT THE CELL
                BELVAVG=0.0
                RVAL=0.0
                IF( SUBO(LE) > 0.5 .AND. BELV(LE) < BELV(L) )THEN
                  RVAL=RVAL+1.
                ENDIF
                IF( SUBO(L)  > 0.5 .AND. BELV(LW)<BELV(L) )THEN
                  RVAL=RVAL+1.
                ENDIF
                IF( SVBO(LN) > 0.5 .AND. BELV(LN)<BELV(L) )THEN
                  RVAL=RVAL+1.
                ENDIF
                IF( SVBO(L)  > 0.5 .AND. BELV(LS)<BELV(L) )THEN
                  RVAL=RVAL+1.
                ENDIF
                IF( RVAL > 0. .OR. HP(L) < HDRY*2. )THEN
                  HOLDTMP=HP(L)
                  H1P(L) = HP(L)
                  HP(L)  = 0.90*HDRY
                  DO K=KSZ(L),KC 
                    HPK(L,K)  = HP(L)*DZC(L,K)
                    HPKI(L,K) = 1./HPK(L,K)
                    H1PK(L,K) = H1P(L)*DZC(L,K)
                  ENDDO

                  NATDRY(L) = 0
                  QDWASTE(L) = DELTI*DXYP(L)*(HOLDTMP-HP(L))
                  VDWASTE(L) = VDWASTE(L) + DXYP(L)*(HOLDTMP-HP(L))
                  NNATDRY(ND) = NNATDRY(ND)+1
                  LNATDRY(ND,NNATDRY(ND)) = L

                ENDIF
              ENDIF
            ENDIF
          ENDIF  ! *** END OF WET CELL TEST
        ENDDO
      END IF
    END IF            

    !**********************************************************************C
    ! *** PERFORM FINAL UPDATES OF HU, AND HV  
    IF( ISDRY == 0 )THEN
      DO L=LF,LL        
        P(L) = G*(HP(L)+BELV(L))  
      ENDDO  
    ENDIF
    IF( IGRIDV == 0 )THEN
      DO L=LF,LL  
        LS=LSC(L)
        LW=LWC(L)
        HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )*FSGZUDXYPI(L)
        HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )*FSGZVDXYPI(L)
      ENDDO  

      DO L=LF,LL  
        HPI(L)=1./HP(L)  
        HUI(L)=1./HU(L)  
        HVI(L)=1./HV(L)  
      ENDDO  
    ENDIF

    ! *** SET TRANSPORT MASK FOR DRY CELLS  
    IF( ISDRY > 0 )THEN  
      DO L=LF,LL  
        LMASKDRY(L)=.TRUE.  
      END DO  
      DO L=LF,LL  
        IF( HP(L) < HDRY )THEN
          LE=LEC(L)
          LN=LNC(L) 
          IUW=0  
          IUE=0  
          IVS=0  
          IVN=0  
          ! *** THIS REQUIRES THE CELL HAVE HP<HDRY FOR 2 ITERATIONS
          IF( SUB1(L)   < 0.5 .AND. SUB(L)   < 0.5 )IUE=1  
          IF( SUB1(LE)  < 0.5 .AND. SUB(LE)  < 0.5 )IUW=1  
          IF( SVB1(L)   < 0.5 .AND. SVB(L)   < 0.5 )IVS=1  
          IF( SVB1(LN)  < 0.5 .AND. SVB(LN)  < 0.5 )IVN=1  
          IFACE=IUW+IUE+IVS+IVN  
          IF( IFACE == 4 )THEN  
            LMASKDRY(L) = .FALSE.  
          END IF  
        ENDIF
      END DO  
    END IF  

    ! *** PERFORM UPDATE ON GROUNDWATER ELEVATION  
    IF( ISGWIE >= 1 )THEN  
      DO L=LF,LL  
        QSUM(L,KC)     = QSUM(L,KC)     - EVAPSW(L)  
        QSUM(L,KSZ(L)) = QSUM(L,KSZ(L)) + QGW(L)   ! *** 2018-10-24 PMC CHANGED SIGN CONVENTION FOR SEEPAGE +(IN), -(OUT) 
      ENDDO  

      ! *** INFILTRATION STEP  
      RNPORI=1./RNPOR  
      IF( ISTL == 3 )THEN  
        DO L=LF,LL  
          AGWELV(L) = AGWELV2(L) - RNPORI*DELT*DXYIP(L)*QGW(L)  
        ENDDO  
      ELSE  
        DO L=LF,LL  
          AGWELV(L) = AGWELV1(L) - RNPORI*DELT*DXYIP(L)*QGW(L)  
        ENDDO  
      ENDIF  
      DO L=LF,LL  
        AGWELV(L) = MIN(AGWELV(L),BELV(L))  
      ENDDO  

      ! *** ET STEP  
      DO L=LF,LL  
        IF( IEVAP > 1 )THEN  
          SVPW1 = (10.**((0.7859+0.03477*TEM(L,KC))/(1.+0.00412*TEM(L,KC))))  
          EVAPT(L) = CLEVAP(L)*0.7464E-3*WINDST(L)*(SVPW1-VPAT(L))/PATMT(L)  
        ENDIF  
        ETGWTMP = EVAPT(L) - EVAPSW(L)*DXYIP(L)      ! *** EXCESS EVAPORATION
        ETGWTMP = MAX(ETGWTMP,0.0)  
        ETGWAVL = RNPOR*DELTI*(AGWELV(L)-BELAGW(L))  
        ETGWAVL = MAX(ETGWAVL,0.0)  
        ETGWTMP = MIN(ETGWTMP,ETGWAVL)  
        EVAPGW(L) = ETGWTMP*DXYP(L)                  ! *** TRANSPIRATION
      ENDDO  
      DO L=LF,LL  
        AGWELV(L)=AGWELV(L)-RNPORI*DELT*DXYIP(L)*EVAPGW(L)  
      ENDDO  
      DO L=LF,LL  
        AGWELV(L)=MAX(AGWELV(L),BELAGW(L))  
      ENDDO  
    ENDIF  

    IF( ISDRY > 0 )THEN
      DO K=1,KC
        DO L=LF,LL
          IF( LKSZ(L,K) )CYCLE
          SUB3D(L,K) = SUB(L)*SUB3DO(L,K)
          SVB3D(L,K) = SVB(L)*SVB3DO(L,K)
        ENDDO
      ENDDO
    ENDIF

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END PARALLEL DO
  
  ! *** IF ANY CELLS WASTED WATER THEN REPORT
  NTMP = SUM(NNATDRY)
  IF( NTMP > 0 )THEN
    IF( INOTICE /= INT(TIMEDAY) )THEN    ! *** ONLY DISPLAY WARNING ONCE PER DAY
      PRINT '(A,F14.5,I6)',' QDWASTED CELLS. SEE EFDCLOG.OUT FOR DETAILS. [TIMEDAY, # OF CELLS]: ',TIMEDAY,NTMP
      INOTICE = INT(TIMEDAY)
    ENDIF

    OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  
    DO ND=1,NOPTIMAL
      DO I=1,NNATDRY(ND)
        L=LNATDRY(ND,I)
        
        TMPVAL=QDWASTE(L)/DXYP(L)
        WRITE(8,8888)TIMEDAY,ND,IL(L),JL(L),TIMEDAY,H1P(L),HP(L),QDWASTE(L),TMPVAL
        
        ! *** UPDATE H1P FOR THE DRY/WASTED CELL
        H1P(L) = 0.90*HDRY
        DO K=KSZ(L),KC 
          H1PK(L,K) = H1P(L)*DZC(L,K)
        ENDDO
        
        ! *** ZERO QDWASTE IF NOT CONDUCTING MASS BALANCE
        IF( ISBAL == 0 ) QDWASTE(L) = 0.0
      ENDDO
    ENDDO
    CLOSE(8)
    8888 FORMAT(' QDWASTE ',F12.5,3I6,F12.4,2F10.4,E14.6,F10.4)
  ENDIF

  ! *** CHECK FOR NEGATIVE DEPTHS  
  CALL NEGDEP(NOPTIMAL,LDMOPT,QCHANUT,QCHANVT,2,SUB1,SVB1,NNEGFLG)  

  ! *** CALCULATE THE EXTERNAL DIVERGENCE  
  IF( ISDIVEX > 0 )THEN  
    DIVEXMX=0.  
    DIVEXMN=1000000.  
    DO L=2,LA  
      IF( SPB(L) /= 0 )THEN  
        LN=LNC(L)  
        DIVEX=SPB(L)*(DXYP(L)*(HP(L)-H1P(L))*DELTI+0.5*(UHDYE(LEC(L))+UHDY1E(LEC(L))-UHDYE(L)-UHDY1E(L)+VHDXE(LN)+VHDX1E(LN)-VHDXE(L)-VHDX1E(L))-QSUME(L)-QGW(L)+EVAPSW(L))  
        IF( ISDIVEX == 2 ) DIVEX = DIVEX*DXYIP(L)*HPI(L)  !  *** RELATIVE DIVERGENCE
        IF( DIVEX > DIVEXMX )THEN  
          DIVEXMX=DIVEX  
          LMAX=L  
        ENDIF  
        IF( DIVEX < DIVEXMN )THEN  
          DIVEXMN=DIVEX  
          LMIN=L  
        ENDIF  
      ENDIF  
    ENDDO  
    IMAX=IL(LMAX)  
    JMAX=JL(LMAX)  
    IMIN=IL(LMIN)  
    JMIN=JL(LMIN)  
    WRITE(6,6628)DIVEXMX,IMAX,JMAX  
    WRITE(6,6629)DIVEXMN,IMIN,JMIN  
    6628 FORMAT('  DIVEXMX=',E13.5,5X,2I10)  
    6629 FORMAT('  DIVEXMN=',E13.5,5X,2I10)  
  ENDIF  

  ! *** DETERMINE THE WET AND DRY CELL LIST
  IF( ISDRY == 0 )THEN
    LAWET = LA-1
    LADRY = 0
  ELSE
    LAWET = 0
    LADRY = 0
    DO L=2,LA
      IF( LMASKDRY(L) )THEN
        LAWET = LAWET+1
        LWET(LAWET)=L
      ELSEIF( OLDMASK(L) .OR. N < NTSTBC+1 )THEN
        ! *** ONLY FLAG NEWLY DRY CELLS
        LADRY = LADRY+1
        LDRY(LADRY)=L
      ENDIF
    ENDDO
  ENDIF
  LDMWET = INT(FLOAT(LAWET)/FLOAT(NDM))+1
  LDMDRY = INT(FLOAT(LADRY)/FLOAT(NDM))+1

  ! *** GET CELL LIST FOR ACTIVE LAYERS = 1
  IF( IGRIDV > 0 .AND. KMINV == 1 )THEN
    LASGZ1 = 0
    DO L=2,LA
      IF( KSZ(L) == KC )THEN
        IF( LMASKDRY(L) )THEN
          LASGZ1 = LASGZ1+1
          LSGZ1(LASGZ1)=L
        ENDIF
      ENDIF
    ENDDO
    LDMSGZ1 = INT(LASGZ1/NDM)+1
  ENDIF
  
  ! *** COMPUTATIONAL CELL LIST FOR EACH SUB-DOMAIN
  IF( ISDRY > 0 )THEN
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LDMWET,LAWET,LWET,LKSZ,LLWET,LKWET,KS,KC) PRIVATE(ND,LF,LL,K,LN,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
    
      DO K=1,KC
        LN=0
        DO LP=LF,LL
          L=LWET(LP)  
          IF( LKSZ(L,K) )CYCLE
          LN = LN+1
          LKWET(LN,K,ND)=L
        ENDDO
        LLWET(K,ND)=LN    ! *** NUMBER OF WET CELLS FOR THE CURRENT LAYER
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO

    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LLWET,LLWETZ,LKWET,LKWETZ,KS,KC) PRIVATE(ND,K,LP,L)
    DO ND=1,NDM
      DO K=1,KS
        LLWETZ(K,ND) = LLWET(K,ND)
        DO LP=1,LLWET(K,ND)
          LKWETZ(LP,K,ND) = LKWET(LP,K,ND)  
        ENDDO 
      ENDDO
      
      LLWETZ(KC,ND) = LLWET(KS,ND)
      DO LP=1,LLWET(KS,ND)
        LKWETZ(LP,KC,ND) = LKWET(LP,KS,ND)  
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    
  ENDIF

  ! *** THIRD PASS CELL CONSTANTS
  IF( IGRIDV > 0 )THEN
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NOPTIMAL,LDMOPT,LA,LWC,LSC,KSZ,DZC,HP,HPK,HU,HV,HPI,HUI,HVI)        &
    !$OMP                           SHARED(DXYP,FSGZUDXYPI,FSGZVDXYPI)  &
    !$OMP                           PRIVATE(ND,LF,LL,L,LW,LS) 
    DO ND=1,NOPTIMAL  
      LF=2+(ND-1)*LDMOPT  
      LL=MIN(LF+LDMOPT-1,LA)

      ! *** UPDATE HU & HV FOR ALL CELLS
      DO L=LF,LL
        LW=LWC(L)
        LS=LSC(L)

        IF( KSZ(LW) > KSZ(L) )THEN
          HU(L) = MAX( 0.5*HPK(L,KSZ(LW)), HP(LW)*(1.+DZC(L,KSZ(LW))*0.1) )
        ELSEIF( KSZ(LW) < KSZ(L) )THEN
          HU(L) = MAX( 0.5*HPK(LW,KSZ(L)), HP(L)*(1.+DZC(LW,KSZ(L))*0.1) )
        ELSE
          HU(L) = ( DXYP(L)*HP(L) + DXYP(LW)*HP(LW) )*FSGZUDXYPI(L)
        ENDIF

        IF( KSZ(LS) > KSZ(L) )THEN
          HV(L) = MAX( 0.5*HPK(L,KSZ(LS)), HP(LS)*(1.+DZC(L,KSZ(LS))*0.1) )
        ELSEIF( KSZ(LS) < KSZ(L) )THEN
          HV(L) = MAX( 0.5*HPK(LS,KSZ(L)), HP(L)*(1.+DZC(LS,KSZ(L))*0.1) )
        ELSE
          HV(L) = ( DXYP(L)*HP(L) + DXYP(LS)*HP(LS) )*FSGZVDXYPI(L)
        ENDIF
      ENDDO
  
      DO L=LF,LL
        HPI(L)=1./HP(L)  
        HUI(L)=1./HU(L)  
        HVI(L)=1./HV(L)  
      ENDDO  

    ENDDO
    !$OMP END PARALLEL DO
  ENDIF

  ! *** COMPUTATIONAL CELL LIST FOR ENTIRE DOMAIN
  IF( ISDRY > 0 )THEN
    DO K=1,KC
      LN=0
      DO L=2,LA
        IF( LKSZ(L,K) )CYCLE
        IF( LMASKDRY(L) )THEN
          LN = LN+1
          LKWET(LN,K,0)=L   ! *** Wet Cell for Layer K
        ENDIF
      ENDDO
      LLWET(K,0) =LN        ! *** Total Wet Cells for Layer K
    ENDDO
  ENDIF
  
  RETURN  

END  

