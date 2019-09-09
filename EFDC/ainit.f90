SUBROUTINE AINIT

  ! CHANGE RECORD
  !
  !  ALL ZEROING OF ARRAYS MOVED TO VARZEROINT AND VARZEROREAL

  USE GLOBAL
  IMPLICIT NONE  
                                                                                                        
  INTEGER :: L,I,J,LS,NT,LCHNV,IVAL,NS,K,NMD,LHOST,LCHNU,NV,NX,LW                                                       
  INTEGER :: NTMPC,NTMPN                                                                                                   

  ! **  INITIALIZE ARRAYS                                                                                                 
  ZBR(1)=ZBRADJ
  ZBRE(1)=ZBRADJ
  HMP(1)=HMIN
  HMU(1)=HMIN
  HMV(1)=HMIN
  HWQ(1)=HMIN
  H2WQ(1)=HMIN
  DXP(1)=DX
  DYP(1)=DY
  DXU(1)=DX
  DYU(1)=DY
  DXV(1)=DX
  DYV(1)=DY
  DXYP(1)=DX*DY
  MVEGL(1)=1
  BELV(1)=BELV(2)
  
  ZBR(LC)=ZBRADJ
  ZBRE(LC)=ZBRADJ
  HMP(LC)=HMIN
  HMU(LC)=HMIN
  HMV(LC)=HMIN
  HWQ(LC)=HMIN
  H2WQ(LC)=HMIN
  DXP(LC)=DX
  DYP(LC)=DY
  DXU(LC)=DX
  DYU(LC)=DY
  DXV(LC)=DX
  DYV(LC)=DY
  DXYP(LC)=DX*DY
  MVEGL(LC)=1
  BELV(LC)=BELV(LA)
  BELV0=BELV
  
  IF( ISGWIE == 0 ) DAGWZ=0.
  DO L=2,LA
    I=IL(L)
    J=JL(L)
    KBT(L)=1
    BELAGW(L)=BELV(L)-DAGWZ
    ZBRE(L)=ZBR(L)
    DLON(L)=CDLON1+(CDLON2*FLOAT(I)+CDLON3)/60.
    DLAT(L)=CDLAT1+(CDLAT2*FLOAT(J)+CDLAT3)/60.
    CUE(L)=1.
    CVE(L)=0.
    CUN(L)=0.
    CVN(L)=1.
  ENDDO
  KBT(1)=1
  KBT(LC)=1
  
  DO L=2,LA
    LW=LWC(L)
    LS=LSC(L)
    DXU(L)=0.5*(DXP(L)+DXP(LW))
    DYU(L)=0.5*(DYP(L)+DYP(LW))
    DXV(L)=0.5*(DXP(L)+DXP(LS))
    DYV(L)=0.5*(DYP(L)+DYP(LS))
  ENDDO
  
  DO L=2,LA
    LW=LWC(L)
    LS=LSC(L)
    HMU(L)=0.5*(DXP(L)*DYP(L)*HMP(L)+DXP(LW)*DYP(LW)*HMP(LW))/(DXU(L)*DYU(L))
    HMV(L)=0.5*(DXP(L)*DYP(L)*HMP(L)+DXP(LS)*DYP(LS)*HMP(LS))/(DXV(L)*DYV(L))
  ENDDO

  HMU(1)=HMU(2)    ! *** PMC
  HMV(1)=HMV(2)    ! *** PMC
  HMU(LC)=HMU(LA)  ! *** PMC
  HMV(LC)=HMV(LA)  ! *** PMC
  DO L=1,LC
    CC(L)=1.
    CCC(L)=1.
    P(L)=G*(HMP(L)+BELV(L))
    P1(L)=G*(HMP(L)+BELV(L))
    HP(L)=HMP(L)
    HU(L)=HMU(L)
    HV(L)=HMV(L)
    HPI(L)=1./HP(L)
    HUI(L)=1./HU(L)
    HVI(L)=1./HV(L)
    HWQ(L)=HMP(L)
    H1P(L)=HMP(L)
    H1U(L)=HMU(L)
    H1V(L)=HMV(L)
    H1UI(L)=1./H1U(L)
    H1VI(L)=1./H1V(L)
    H2WQ(L)=HMP(L)
    SCB(L)=1.
    SPB(L)=1.
    SUB(L)=1.
    SVB(L)=1.
    SWB(L)=1.
    STCUV(L)=1.
    STCAP(L)=1.
    STBX(L)=1.
    STBY(L)=1.
    SAAX(L)=1.
    SAAY(L)=1.
    SCAX(L)=1.
    SCAY(L)=1.
    SBX(L)=1.
    SBY(L)=1.
    SDX(L)=1.
    SDY(L)=1.
    LMASKDRY(L)=.TRUE.
  ENDDO
  
  LOPENBCDRY = .FALSE.
  IF( ISTRAN(8) > 0 )THEN
    RADKE = WQKEB(1)
  ELSE
    RADKE = SWRATNF
  ENDIF
  
  ! *** OPEN WATER DEFAULT SETTINGS
  IF( ISVEG > 0 )THEN
    NV=0
    PVEGX(NV)=1.
    PVEGY(NV)=1.
    PVEGZ(NV)=1.
  ENDIF
  
  IF( ISTRAN(5) > 0 )THEN
    DO NT=1,NTOX
      DO K=1,KB
        DO L=1,LC
          TOXB(L,K,NT)=TOXBINIT(L,K,NT)
          TOXB1(L,K,NT)=TOXBINIT(L,K,NT)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  IF( ISTRAN(6) > 0 .AND.  .NOT. LSEDZLJ )THEN
    DO NS=1,NSED
      DO K=1,KB
        DO L=1,LC
          SEDB(L,K,NS)=SEDBINIT(L,K,NS)
          SEDB1(L,K,NS)=SEDBINIT(L,K,NS)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  IF( ISTRAN(7) > 0 )THEN
    SNDDMX=0.
    DO NS=1,NSND
      NX=NS+NSED
      DO K=1,KB
        DO L=1,LC
          SNDB(L,K,NS)=SNDBINIT(L,K,NS)
          SNDB1(L,K,NS)=SNDBINIT(L,K,NS)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  
  !      IF( IS1DCHAN == 1 )THEN
  !        DO L=1,LC
  !          FADYP(L)=1.
  !          FADYP1(L)=1.
  !          FADYP2(L)=1.
  !          WPDYP(L)=1.
  !          WPDYP1(L)=1.
  !          FADXP(L)=1.
  !          FADXP1(L)=1.
  !          FADXP2(L)=1.
  !          WPDXP(L)=1.
  !          WPDXP1(L)=1.
  !          FADYU(L)=1.
  !          FADYU1(L)=1.
  !          WPDYU(L)=1.
  !          WPDYU1(L)=1.
  !          FADXV(L)=1.
  !          FADXV1(L)=1.
  !          WPDXV(L)=1.
  !          WPDXV1(L)=1.
  !          DADH(L)=1.
  !          DADH1(L)=1.
  !          SRFXP(L)=0.
  !          SRFYP(L)=0.  C
  !          SRFXP1(L)=0.
  !          SRFYP1(L)=0.
  !          SRFXV(L)=0.
  !          SRFYU(L)=0.
  !          SRFXV1(L)=0.
  !          SRFYU1(L)=0.
  !        ENDDO
  !      ENDIF

  DO K=1,KS
    DO L=1,LC
      AV(L,K)=AVO
      AVVI(L,K)=1./AVO
      AVUI(L,K)=1./AVO
      AB(L,K)=ABO
      QQL(L,K)=QQLMIN
      QQL1(L,K)=QQLMIN
      QQL2(L,K)=QQLMIN
      DML(L,K)=DMLMIN
      ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
    ENDDO
  ENDDO
  DO K=1,KC
    DO L=1,LC
      AH(L,K)=AHO
      !AHU(L,K)=AHO
      !AHV(L,K)=AHO
      AHC(L,K)=AHO
      AQ(L,K)=AVO
      ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      CTURBB1(L,K)=CTURB
      CTURBB2(L,K)=CTURB2B
      ! *** TEMPERATURE INITIATION
      TEM(L,K)=TEMO
      TEM1(L,K)=TEMO
    ENDDO
  ENDDO
  IF( ISWQFLUX == 1 )THEN
    DO K=1,KC
      DO L=1,LC
        AHULPF(L,K)=AHO
        AHVLPF(L,K)=AHO
      ENDDO
    ENDDO
  ENDIF
  
  NTMPC=MAX(NSED,1)
  DO NS=1,NTMPC
    DO K=1,KC
      DO L=1,LC
        SED(L,K,NS)=SEDO(NS)
        SED1(L,K,NS)=SEDO(NS)
  ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      ENDDO
    ENDDO
  ENDDO

  NTMPN=MAX(NSND,1)
  DO NX=1,NTMPN
    NS=NX+NTMPC
    DO K=1,KC
      DO L=1,LC
        SND(L,K,NX)=SEDO(NS)
        SND1(L,K,NX)=SEDO(NS)
  ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      ENDDO
    ENDDO
  ENDDO
  DO NT=1,NTOX
    DO K=1,KC
      DO L=1,LC
        TOX(L,K,NT)=TOXINTW(NT)
        TOX1(L,K,NT)=TOXINTW(NT)
        ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      ENDDO
    ENDDO
  ENDDO

  DO K=0,KC
    DO L=1,LC
      ! *** ALL ZEROING OF ARRAYS MOVED TO ZERO
      QQ(L,K)=QQMIN
      QQ1(L,K)=QQMIN
      QQ2(L,K)=QQMIN
      QQSQR(L,K)=SQRT(QQMIN)
    ENDDO
  ENDDO

  IF( MDCHH >= 1 )THEN
    DO NMD=1,MDCHH
      LHOST=LMDCHH(NMD)
      LCHNU=LMDCHU(NMD)
      LCHNV=LMDCHV(NMD)
  !
  !         SET HOST DRYING DEPTH
  !
      IF( PMDCH(NMD) < 0.0) PMDCH(NMD)=HWET
  !
  !         X-DIRECTION CHANNEL
  !
      IF( MDCHTYP(NMD) == 1 )THEN
        IF( CHANLEN(NMD) < 0.0 )THEN
          CHANLEN(NMD)=0.25*DYP(LHOST)
        ELSE
          CHANLEN(NMD)=CHANLEN(NMD)-0.5*DYP(LCHNU)
        ENDIF
      ENDIF
  !
  !         Y-DIRECTION CHANNEL
  !
      IF( MDCHTYP(NMD) == 2 )THEN
        IF( CHANLEN(NMD) < 0.0 )THEN
          CHANLEN(NMD)=0.25*DXP(LHOST)
        ELSE
          CHANLEN(NMD)=CHANLEN(NMD)-0.5*DXP(LCHNV)
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  ! ** INITIALIZE ORGANIC CARBON VARIABLES OF SEDIMENT-TOXICS
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) > 0 )IVAL=1
  ENDDO

  IF( IVAL == 0 .AND. ISTRAN(5) > 0 )THEN
    ! *** Kd APPROACH ONLY USED IN MODEL
    DO NS=1,NSED+NSND
      DO K=1,KB
        DO L=1,LC
          STFPOCB(L,K,NS)=1.0
        ENDDO
      ENDDO
    ENDDO
    DO NS=1,NSED+NSND
      DO K=1,KC
        DO L=1,LC
          STFPOCW(L,K,NS)=1.0
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF( IVAL == 1 .AND. ISTRAN(5) > 0 )THEN
    ! *** IF ANY fPOC, POC AND DOC OPTIONS USED INITIALIZE HERE
    
    ! *** SEDIMENT BED
    IF( ISTDOCB == 0 )THEN  !   ISTDOCB == 1 STDOCB IS INITIALIZED FROM DOCB.INP
      DO K=1,KB
        DO L=1,LC
          STDOCB(L,K)=STDOCBC
        ENDDO
      ENDDO
    ENDIF
    IF( ISTPOCB == 0 )THEN  !   ISTPOCB == 1 STPOCB IS INITIALIZED FROM POCB.INP
      DO K=1,KB
        DO L=1,LC
          STPOCB(L,K)=STPOCBC
        ENDDO
      ENDDO
    ENDIF
    IF( ISTPOCB /= 3 )THEN
      ! *** SPATIALLY CONSTANT FPOC.  ISTPOCB == 3 STFPOCB IS INITIALIZED FROM FPOCB.INP
      DO NS=1,NSED+NSND
        DO K=1,KB
          DO L=1,LC
            STFPOCB(L,K,NS)=FPOCBST(NS,1)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    
    ! *** WATER COLUMN
    IF( ISTDOCW == 0 )THEN
      DO K=1,KC
        DO L=1,LC
          STDOCW(L,K)=STDOCWC
        ENDDO
      ENDDO
    ENDIF
    IF( ISTPOCW == 0 )THEN
      IF( STPOCWC <= 0.0 ) STPOCWC = 1E-12
      DO K=1,KC
        DO L=1,LC
          STPOCW(L,K)=STPOCWC
        ENDDO
      ENDDO
    ENDIF
    IF( ISTPOCW /= 3 )THEN
      ! *** SPATIALLY CONSTANT FPOC.  ISTPOCW == 3 STFPOCB IS INITIALIZED FROM FPOCW.INP
      DO NS=1,NSED+NSND
        DO K=1,KC
          DO L=1,LC
            STFPOCW(L,K,NS)=FPOCWST(NS,1)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  ! *** NEW VARIABLES FOR QCTL NQCTYP=3 & 4 
  LOWCHORDU=-9999.
  LOWCHORDV=-9999.
  NLOWCHORD=0

  RETURN
END

