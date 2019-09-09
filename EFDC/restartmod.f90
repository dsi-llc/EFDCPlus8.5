MODULE RESTART_MODULE
  ! ** THIS MODULE INCLUDES THE SUBROUTINES FOR 
  ! ** WRITTING AND READING RESTART FILES:
  ! ** RESTOUT, WSMRST, WQWCRST, RESTMOD
  ! ** RESTIN1, RESTIN2, RESTIN10, RSMRST, WQWCRST_IN
  ! ** ALL THE BINARY OUTPUT FILES WILL BE APPENDED
  ! ** UPDATED: 2013-06-23 BY DH CHUNG
  USE GLOBAL
  USE JULIANMOD
  USE IFPORT
  USE INFOMOD,ONLY: NUMCOL
  
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE RESTOUT(IRSTYP)
       
    ! ** GENERATE A RESTART FILE FOR USE BY BY EFDC_DSI
  
    INTEGER,INTENT(IN) :: IRSTYP
    INTEGER :: K,L,LL,M,NT,NX,IDUM,JDUM,NS,NMD,NHR,VER,NCTL,IMASK,NP
    
    REAL :: RVAL,TMP1,TMP2,TMP3,TMP4,SURF
    CHARACTER(2) :: SNHR

    IF( ISGREGOR == 1 )CALL TOGREGOR(SDATE,TIMEDAY)

    IF( IRSTYP == 0 )THEN
      PRINT *,'Restart Snapshot @ Timeday: ',SNGL(TIMEDAY)
      IF( ISGREGOR == 1 )THEN
        NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
        WRITE(SNHR,'(I2.2)') NHR
        RESFILE = 'RESTART'//'_'//TRIM(SDATE)//'_'//SNHR//'.OUT'  
      ELSE
        RESFILE = 'RESTART.OUT'  
      ENDIF
      OPEN(99,FILE=OUTDIR//TRIM(RESFILE),STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE=OUTDIR//TRIM(RESFILE),STATUS='UNKNOWN')
    ENDIF
    IF( IRSTYP == 1 )THEN
      OPEN(99,FILE=OUTDIR//'CRASHST.OUT',STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE=OUTDIR//'CRASHST.OUT',STATUS='UNKNOWN')
    ENDIF
    IF( ISRESTO == -11 )THEN
      DO L=1,LC
        HP(L)=-BELV(L)
        H1P(L)=-BELV(L)
        HWQ(L)=-BELV(L)
        H2WQ(L)=-BELV(L)
        UHDYE(L)=0.
        UHDY1E(L)=0.
        VHDXE(L)=0.
        VHDX1E(L)=0.
      ENDDO
      DO K=0,KC
        DO L=1,LC
          QQ(L,K)=QQMIN
          QQSQR(L,K)=SQRT(QQ(L,K))  ! *** DS-INTL
          QQ1(L,K)=QQMIN
          QQL(L,K)=QQLMIN
          QQL1(L,K)=QQLMIN
          DML(L,K)=DMLMIN
        ENDDO
      ENDDO
      DO K=1,KC
        DO L=1,LC
          U(L,K)=0.
          U1(L,K)=0.
          V(L,K)=0.
          V1(L,K)=0.
        ENDDO
      ENDDO
      DO K=1,KC
        DO L=1,LC
          SAL(L,K)=MAX(SAL(L,K),0.)
          SAL1(L,K)=MAX(SAL1(L,K),0.)
        ENDDO
      ENDDO
    ENDIF

    ! *** BEGIN RESTART.OUT 
    VER = 720
    WRITE(99,'(I20,4X,F12.4,2I10)')N,TIMEDAY,VER    
    WRITE(99,'(9I6)') ISCO(0:8)
    DO L=2,LA
      WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L),BELV(L)
      WRITE(99,906)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
      WRITE(99,906)(U(L,K),K=1,KC)
      WRITE(99,906)(U1(L,K),K=1,KC)
      WRITE(99,906)(V(L,K),K=1,KC)
      WRITE(99,906)(V1(L,K),K=1,KC)
      WRITE(99,906)(QQ(L,K),K=0,KC)
      WRITE(99,906)(QQ1(L,K),K=0,KC)
      WRITE(99,906)(QQL(L,K),K=0,KC)
      WRITE(99,906)(QQL1(L,K),K=0,KC)
      WRITE(99,906)(DML(L,K),K=0,KC)
      IF( ISTRAN(2) >= 1 .AND. ISICE > 2 )THEN
        TMP1 = ICETHICK(L) + ICEVOL(L)*999.8426/RHOI*DXYIP(L)
        IF( ISICE == 4 ) TMP1 = TMP1 + SUM(FRAZILICE(L,KSZ(L):KC))
        WRITE(99,906)TMP1,ICETEMP(L)                ! ** EFDC_7.3
      ENDIF        
      IF( ISTRAN(1) >= 1 .AND. ISCO(1) == 1 )THEN
        WRITE(99,906)(SAL(L,K),K=1,KC)
        WRITE(99,906)(SAL1(L,K),K=1,KC)
      ENDIF
      IF( ISTRAN(2) >= 1 .AND. ISCO(2) == 1 )THEN
        WRITE(99,906) TEMB(L),(TEM(L,K),K=1,KC)     ! ** EFDC_7.2
        WRITE(99,906)(TEM1(L,K),K=1,KC)
      ENDIF
      IF( ISTRAN(3) >= 1 .AND. ISCO(3) == 1 )THEN
        WRITE(99,906)(DYE(L,K),K=1,KC)
        WRITE(99,906)(DYE1(L,K),K=1,KC)
      ENDIF
      IF( ISTRAN(4) >= 1 .AND. ISCO(4) == 1 )THEN
        WRITE(99,906)SFLSBOT(L),(SFL(L,K),K=1,KC)
        WRITE(99,906)SFLSBOT(L),(SFL2(L,K),K=1,KC)
      ENDIF
      IF( ISTRAN(5) >= 1 .AND. ISCO(5) == 1 )THEN
        DO NT=1,NTXM
          WRITE(99,906)(TOXB(L,K,NT),K=1,KB)
          WRITE(99,906)(TOX(L,K,NT),K=1,KC)
          WRITE(99,906)(TOXB1(L,K,NT),K=1,KB)
          WRITE(99,906)(TOX1(L,K,NT),K=1,KC)
        ENDDO
      ENDIF
      IF( ISTRAN(6) >= 1 .AND. ISCO(6) == 1 )THEN
        DO NS=1,NSCM
          WRITE(99,906)(SEDB(L,K,NS),K=1,KB)
          WRITE(99,906)(SED(L,K,NS),K=1,KC)
          WRITE(99,906)(SEDB1(L,K,NS),K=1,KB)
          WRITE(99,906)(SED1(L,K,NS),K=1,KC)
        ENDDO
      ENDIF
      IF( ISTRAN(7) >= 1 .AND. ISCO(7) == 1 )THEN
        DO NS=1,NSNM
          WRITE(99,906)(SNDB(L,K,NS),K=1,KB)
          WRITE(99,906)(SND(L,K,NS),K=1,KC)
          WRITE(99,906)(SNDB1(L,K,NS),K=1,KB)
          WRITE(99,906)(SND1(L,K,NS),K=1,KC)
        ENDDO
      ENDIF
      IF( (ISTRAN(6) >= 1 .AND. ISCO(6) == 1 ) .OR. (ISTRAN(7) >= 1 .AND. ISCO(7) == 1 ) )THEN
        WRITE(99,906)(HBED(L,K),K=1,KB)
        WRITE(99,906)(HBED1(L,K),K=1,KB)
        WRITE(99,906)(VDRBED(L,K),K=1,KB)
        WRITE(99,906)(VDRBED1(L,K),K=1,KB)
      ENDIF
    ENDDO

    ! *** BOUNDARY CONDITIONS
    DO M=1,4
      IF( ISTRAN(M) >= 1 .AND. ISCO(M) == 1 )THEN
        DO LL=1,NCBS
          DO K=1,KC
            NLOS(LL,K,M)=MAX( NLOS(LL,K,M)-N , -NTSCRS(LL) )
          ENDDO
          WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOS(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBW
          DO K=1,KC
            NLOW(LL,K,M)=MAX( NLOW(LL,K,M)-N , -NTSCRW(LL) )
          ENDDO
          WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOW(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBE
          DO K=1,KC
            NLOE(LL,K,M)=MAX( NLOE(LL,K,M)-N , -NTSCRE(LL) )
          ENDDO
          WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOE(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBN
          DO K=1,KC
            NLON(LL,K,M)=MAX( NLON(LL,K,M)-N , -NTSCRN(LL) )
          ENDDO
          WRITE(99,908)(NLON(LL,K,M),K=1,KC)
          WRITE(99,906)(CLON(LL,K,M),K=1,KC)
        ENDDO
      ENDIF
    ENDDO
    IF( ISTRAN(5) >= 1 .AND. ISCO(5) == 1 )THEN
      DO NT=1,NTXM
        M=MSVTOX(NT)
        DO LL=1,NCBS
          DO K=1,KC
            NLOS(LL,K,M)=MAX( NLOS(LL,K,M)-N , -NTSCRS(LL) )
          ENDDO
          WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOS(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBW
          DO K=1,KC
            NLOW(LL,K,M)=MAX( NLOW(LL,K,M)-N , -NTSCRW(LL) )
          ENDDO
          WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOW(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBE
          DO K=1,KC
            NLOE(LL,K,M)=MAX( NLOE(LL,K,M)-N , -NTSCRE(LL) )
          ENDDO
          WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOE(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBN
          DO K=1,KC
            NLON(LL,K,M)=MAX( NLON(LL,K,M)-N , -NTSCRN(LL) )
          ENDDO
          WRITE(99,908)(NLON(LL,K,M),K=1,KC)
          WRITE(99,906)(CLON(LL,K,M),K=1,KC)
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(6) >= 1 .AND. ISCO(6) == 1 )THEN
      DO NT=1,NSCM
        M=MSVSED(NT)
        DO LL=1,NCBS
          DO K=1,KC
            NLOS(LL,K,M)=MAX( NLOS(LL,K,M)-N , -NTSCRS(LL) )
          ENDDO
          WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOS(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBW
          DO K=1,KC
            NLOW(LL,K,M)=MAX( NLOW(LL,K,M)-N , -NTSCRW(LL) )
          ENDDO
          WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOW(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBE
          DO K=1,KC
            NLOE(LL,K,M)=MAX( NLOE(LL,K,M)-N , -NTSCRE(LL) )
          ENDDO
          WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOE(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBN
          DO K=1,KC
            NLON(LL,K,M)=MAX( NLON(LL,K,M)-N , -NTSCRN(LL) )
          ENDDO
          WRITE(99,908)(NLON(LL,K,M),K=1,KC)
          WRITE(99,906)(CLON(LL,K,M),K=1,KC)
        ENDDO
      ENDDO
    ENDIF
    IF( ISTRAN(7) >= 1 .AND. ISCO(7) == 1 )THEN
      DO NT=1,NSNM
        M=MSVSND(NT)
        DO LL=1,NCBS
          DO K=1,KC
            NLOS(LL,K,M)=MAX( NLOS(LL,K,M)-N , -NTSCRS(LL) )
          ENDDO
          WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOS(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBW
          DO K=1,KC
            NLOW(LL,K,M)=MAX( NLOW(LL,K,M)-N , -NTSCRW(LL) )
          ENDDO
          WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOW(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBE
          DO K=1,KC
            NLOE(LL,K,M)=MAX( NLOE(LL,K,M)-N , -NTSCRE(LL) )
          ENDDO
          WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOE(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBN
          DO K=1,KC
            NLON(LL,K,M)=MAX( NLON(LL,K,M)-N , -NTSCRN(LL) )
          ENDDO
          WRITE(99,908)(NLON(LL,K,M),K=1,KC)
          WRITE(99,906)(CLON(LL,K,M),K=1,KC)
        ENDDO
      ENDDO
    ENDIF
    DO L=2,LA
      WRITE(99,906)QSUME(L),(QSUM(L,K),K=1,KC)
    ENDDO
    IF( MDCHH >= 1 )THEN
      DO NMD=1,MDCHH
        WRITE(99,'(6I5,2X,E17.8,2X,E17.8)') IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),IMDCHV(NMD),JMDCHV(NMD),QCHANU(NMD),QCHANV(NMD)
      ENDDO
    ENDIF
    IF( ISGWIE >= 1 )THEN
      DO L=2,LA
        WRITE(99,906) AGWELV(L),AGWELV1(L)
      ENDDO
    ENDIF
    ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE WITH A LOWER CHORD ACTIVATION TOGGLE
    DO NCTL=1,NQCTL
      IF( NQCTYP(NCTL) == 3 .OR. NQCTYP(NCTL) == 4 )THEN
        ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE (4) DEPENDANT FLOWS
        WRITE(99,908)NCTL,NLOWCHORD(NCTL)
        WRITE(99,906)SAVESUB(1,NCTL),SAVESVB(1,NCTL),SAVESUB(2,NCTL),SAVESVB(2,NCTL)
        WRITE(99,906)LOWCHORDU(NCTL),LOWCHORDV(NCTL)
      ENDIF
    ENDDO

    CLOSE(99)
    
    ! *** SEDZLJ
    IF( LSEDZLJ ) THEN
      IF( ISGREGOR == 1 )THEN
        RSTFILE = 'SEDBED_HOT'//'_'//TRIM(SDATE)//'_'//SNHR//'.OUT'  
      ELSE
        RSTFILE = 'SEDBED_HOT.OUT'
      ENDIF
      OPEN(99,FILE=OUTDIR//TRIM(RSTFILE),STATUS='UNKNOWN')
      CLOSE(99, STATUS='DELETE')
      OPEN(99,FILE=OUTDIR//TRIM(RSTFILE),FORM='FORMATTED',STATUS='UNKNOWN')

      WRITE(99,34569)((LAYERACTIVE(K,L),K=1,KB),L=2,LA)
      WRITE(99,34569)(KBT(L),L=2,LA)
      WRITE(99,34567)(D50AVG(L),L=2,LA)
      WRITE(99,34568)((BULKDENS(LL,L),LL=1,KB),L=2,LA)
      WRITE(99,34568)((TSED(LL,L),LL=1,KB),L=2,LA)
      WRITE(99,34568)(((PERSED(K,LL,L),K=1,NSCM),LL=1,KB),L=2,LA)
      CLOSE(99)
34567 FORMAT(E17.9)
34568 FORMAT(6E17.9)
34569 FORMAT(8I8)
    END IF
    
    ! *** LPT RESTART
    IF( ISPD == 2 )THEN
      IF( ISGREGOR == 1 )THEN
        NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
        WRITE(SNHR,'(I2.2)') NHR
        RESFILE = 'EE_DRIFTER'//'_'//TRIM(SDATE)//'_'//SNHR//'.RST'  
      ELSE
        RESFILE = 'EE_DRIFTER.RST'  
      ENDIF
      OPEN(307,FILE=OUTDIR//TRIM(RESFILE),STATUS='UNKNOWN')  
      CLOSE(307, STATUS='DELETE')
      OPEN(307,FILE=OUTDIR//TRIM(RESFILE),STATUS='UNKNOWN')
      WRITE(307,'(50I2)') JSPD(1:NPD)
      WRITE(307,'(30I3)') KLA(1:NPD)
      WRITE(307,'(2(I8,2F15.5,F10.5))') (LLA(NP),XLA(NP),YLA(NP),ZLA(NP),NP=1,NPD)
      CLOSE(307,STATUS='KEEP')
    ENDIF
     
    ! *** SPECIAL FILES
    IF( ISWAVE >= 1 )THEN
      OPEN(1,FILE=OUTDIR//'WVQWCP.OUT',STATUS='UNKNOWN')
      CLOSE(1, STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'WVQWCP.OUT',STATUS='UNKNOWN')
      DO L=2,LA
        WRITE(1,'(2I5,2X,6E13.4)')IL(L),JL(L),QQWV1(L),QQWV2(L),QQWV3(L),QQWC(L),QQWCR(L),QQ(L,0)
      ENDDO
      CLOSE(1)
    ENDIF
    
    IF( ISDRY == 99 )THEN
      IF( ISGREGOR == 1 )THEN
        RSTFILE = 'RSTWD'//'_'//TRIM(SDATE)//'_'//SNHR//'.OUT'  
      ELSE
        RSTFILE = 'RSTWD.OUT'
      ENDIF
      OPEN(1,FILE=OUTDIR//TRIM(RSTFILE),STATUS='UNKNOWN')
      CLOSE(1, STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//TRIM(RSTFILE),STATUS='UNKNOWN')
      DO L=2,LA
        IMASK=0
        IF( .NOT. LMASKDRY(L) )IMASK=1
        WRITE(1,913)L,IL(L),JL(L),ISCDRY(L),NATDRY(L),IMASK,SUB(L),SVB(L),SUBO(L),SVBO(L)
      ENDDO
      CLOSE(1)
    ENDIF
  
    ! **  OUTPUT SALINITY AND TEMPATURE DATA ASSIMILATION
    IF( NLCDA > 0 )THEN
      ! *** DEPRECATED
    ENDIF
  
    IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0) .AND. ISDTXBUG == 1 .AND. N == NTS )THEN
      OPEN(1,FILE=OUTDIR//'BEDRST.SED')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.SED')
      WRITE(1,111)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SEDB(L,K,1),K=1,KB)
        IF( NSCM > 1 )THEN
          DO NX=2,NSCM
            WRITE(1,102)(SEDB(L,K,NX),K=1,KB)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'BEDRST.SND')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.SND')
      WRITE(1,112)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SNDB(L,K,1),K=1,KB)
        IF( NSNM > 1 )THEN
          DO NX=2,NSNM
            WRITE(1,102)(SNDB(L,K,NX),K=1,KB)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'BEDRST.VDR')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.VDR')
      WRITE(1,113)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(VDRBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'BEDRST.POR')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.POR')
      WRITE(1,114)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(PORBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'BEDRST.ZHB')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.ZHB')
      WRITE(1,115)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),(HBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'BEDRST.BDN')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.BDN')
      WRITE(1,116)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(BDENBED(L,K),K=1,KB)
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'BEDRST.ELV')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.ELV')
      WRITE(1,117)
      RVAL=0.
      TMP1=0.
      TMP2=0.
      TMP3=0.
      TMP4=0.
      DO L=2,LA
        RVAL=RVAL+1.
        TMP1=TMP1+ZELBEDA(L)
        TMP2=TMP2+HBEDA(L)
        TMP3=TMP3+BELV(L)
        TMP4=TMP4+HP(L)
        SURF=HP(L)+BELV(L)
        WRITE(1,101)IL(L),JL(L),ZELBEDA(L),HBEDA(L),BELV(L),HP(L),SURF
      ENDDO
      TMP1=TMP1/RVAL
      TMP2=TMP2/RVAL
      TMP3=TMP3/RVAL
      TMP4=TMP4/RVAL
      IDUM=0
      JDUM=0
      WRITE(1,101)IDUM,JDUM,TMP1,TMP2,TMP3,TMP4
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'WATRST.SED')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'WATRST.SED')
      WRITE(1,118)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SED(L,K,1),K=1,KC)
        IF( NSCM > 1 )THEN
          DO NX=2,NSCM
            WRITE(1,102)(SED(L,K,NX),K=1,KC)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'WATRST.SND')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'WATRST.SND')
      WRITE(1,119)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),(SND(L,K,1),K=1,KC)
        IF( NSNM > 1 )THEN
          DO NX=2,NSNM
            WRITE(1,102)(SND(L,K,NX),K=1,KC)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
      OPEN(1,FILE=OUTDIR//'BEDRST.BDL')
      CLOSE(1,STATUS='DELETE')
      OPEN(1,FILE=OUTDIR//'BEDRST.BDL')
      WRITE(1,120)
      DO L=2,LA
        WRITE(1,101)IL(L),JL(L),QSBDLDX(L,1),QSBDLDY(L,1)
        IF( NSNM > 1 )THEN
          DO NX=2,NSNM
            WRITE(1,102)QSBDLDX(L,NX),QSBDLDY(L,NX)
          END DO
        ENDIF
      ENDDO
      CLOSE(1)
      IF( ISTRAN(5) > 0 )THEN
        OPEN(1,FILE=OUTDIR//'BEDRST.TOX')
        CLOSE(1,STATUS='DELETE')
        OPEN(1,FILE=OUTDIR//'BEDRST.TOX')
        DO NT=1,NTXM
          WRITE(1,121)NT
          DO L=2,LA
            WRITE(1,101)IL(L),JL(L),(TOXB(L,K,NT),K=1,KB)
          ENDDO
        ENDDO
        CLOSE(1)
      ENDIF
    ENDIF
    101 FORMAT(2I5,50E13.5)
    102 FORMAT(10X,50E13.5)
    111 FORMAT('   IL   JL    SEDBT(K=1,KB)')
    112 FORMAT('   IL   JL    SNDBT(K=1,KB)')
    113 FORMAT('   IL   JL    VRDBED(K=1,KB)')
    114 FORMAT('   IL   JL    PORBED(K=1,KB)')
    115 FORMAT('   IL   JL    ZBEDB        HBEDT        HBED(K=1,KB)')
    116 FORMAT('   IL   JL    BDENBED(K=1,KB)')
    117 FORMAT('   IL   JL    ZBEDB        HBEDT        BELV', '        HWCOL        SELV')
    118 FORMAT('   IL   JL    SEDT(K=1,KC)')
    119 FORMAT('   IL   JL    SNDT(K=1,KC)')
    120 FORMAT('   IL   JL    QSBDLDX      QSBDLDY')
    121 FORMAT('   IL   JL    TOXB(K=1,KB,NT)  NT = ',I5)
    906 FORMAT(50E17.8)
    908 FORMAT(50I10)
    912 FORMAT(3I5,12F7.3,5(15X,12F7.3))   ! *** DS-INTL Single Line
    913 FORMAT(6I10,4F7.3)                 ! *** DS-INTL Single Line
    
  END SUBROUTINE

  SUBROUTINE WQSDRST
    ! ** WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT ISMORST.
    INTEGER :: L,NW,NHR
    CHARACTER(2) :: SNHR
    
    IF( ISGREGOR ==  1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = OUTDIR//'WQSDRST_'//TRIM(SDATE)//'_'//SNHR//'.OUT'  
    ELSE
      RESFILE = OUTDIR//'WQSDRST.OUT'  
    ENDIF
    OPEN(1,FILE=RESFILE,STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=RESFILE,STATUS='UNKNOWN')

    ! *** BEGIN WQSDRST
    WRITE(1,101) N,TIMEDAY
    WRITE(1,888)
    DO L=2,LA
      WRITE(1,90) L,(SMPON(L,NW),NW=1,NSMG), &
          (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L), &
          SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L), &
          SMBST(L),SMT(L)
    ENDDO
    CLOSE(1)

     90 FORMAT(I5, 1P, 18E12.4)
    101 FORMAT('CC  SM RESTART FILE TIME STEP, TIME = ',I10,F13.5)
    888 FORMAT('    L', &
        '       GPON1       GPON2       GPON3       GPOP1 &
        GPOP2','       GPOP3       GPOC1       GPOC2       GPOC3 &
        G1NH4','       G2NH4       G2NO3       G2PO4       G2H2S &
        GPSI','        G2SI        GBST          GT')
        
  END SUBROUTINE

  SUBROUTINE WQWCRST
    ! ** WRITE SPATIAL DISTRIBUTIONS AT THE END OF SIMULATION TO UNIT IWQORST.
    INTEGER :: NWQV0,L,K,NW,NHR,VER
    CHARACTER(2) :: SNHR
    
    IF( ISGREGOR ==  1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = OUTDIR//'WQWCRST_'//TRIM(SDATE)//'_'//SNHR//'.OUT'  
    ELSE
      RESFILE = OUTDIR//'WQWCRST.OUT'  
    ENDIF
    OPEN(1,FILE=RESFILE,STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=RESFILE,STATUS='UNKNOWN')
    
    ! *** BEGIN WQWCRST
    VER = 8400
    WRITE(1,101) VER,N,TIMEDAY
    WRITE(1,102) WQI1,WQI2,WQI3
    WRITE(1,103)
    NWQV0=NWQV
    IF( IDNOTRVA > 0 ) NWQV0=NWQV0+1
    DO L=2,LA
      DO K=1,KC
        WRITE(1,104) L,K,(WQV(L,K,NW),NW=1,NWQV0)
      ENDDO
    ENDDO
    CLOSE(1)

    101 FORMAT('CC  WQ RESTART FILE TIME STEP, TIME = ',2I10,F12.4)
    102 FORMAT('CC  PREVIOUS DAYS SOLAR RADIATION   = ',3F12.5)
    103 FORMAT('C   L    K  BC          BD          BG          ', &
        'RPOC        LPOC        DOC         ', &
        'RPOP        LPOP        DOP         PTO4        ', &
        'RPON        LPON        DON         AMN         ', &
        'NIT         SU          SA          COD         ', &
        'DO          TAM         FCB        MALG')
    104 FORMAT(2I5, 1P, 22E12.4)
        
   END SUBROUTINE

  SUBROUTINE WQRPEMRST
  
    INTEGER :: L,NHR
    CHARACTER(2) :: SNHR

    ! **********************************************************************C                                                
    ! ** RESTART OUTPUT                                                                                                     
    IF( ISGREGOR ==  1 )THEN
      NHR = NINT((TIMEDAY-INT(TIMEDAY))*24.0)
      WRITE(SNHR,'(I2.2)') NHR
      RESFILE = OUTDIR//'WQRPEMRST_'//TRIM(SDATE)//'_'//SNHR//'.OUT'  
    ELSE
      RESFILE = OUTDIR//'WQRPEMRST.OUT'  
    ENDIF
    OPEN(1,FILE=RESFILE,STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=RESFILE,STATUS='UNKNOWN')

    WRITE(1,110)TIMEDAY                                                                                                   
    DO L=2,LA
      IF( LMASKRPEM(L) )THEN                                                                                               
        WRITE(1,111) L,WQRPS(L),WQRPR(L),WQRPE(L),WQRPD(L)                                                                
      ENDIF
    ENDDO
    CLOSE(1)                                                                                                             
    110 FORMAT(F12.4)                                                                                                  
    111 FORMAT(I6,6E13.5)                                                                                                 
  
  END SUBROUTINE
  
  SUBROUTINE RESTMOD
    ! ** SUBROUTINE RESTOUT WRITES A RESTART FILE
    ! ** ANOTHER FORMAT OF RESTART.INP (ISRESTO == -2)
    INTEGER :: L,K,NNIJ,NIJMOD,ITMP,JTMP,LSMOD,M,LL
    INTEGER :: LIJMOD(100)

    RESFILE = 'RESTART.OUT'
    OPEN(99,FILE=TRIM(RESFILE),STATUS='UNKNOWN')
    CLOSE(99,STATUS='DELETE')
    OPEN(99,FILE=TRIM(RESFILE),STATUS='UNKNOWN')

    ! *** BEGIN RESTMOD
    WRITE(99,'(I20,4X,F15.5)')N,TIMEDAY
    WRITE(*,'(A)')'READING RESTART FILE: RESTMOD.INP'
    OPEN(1,FILE='restmod.inp',STATUS='UNKNOWN')
    READ(1,*)NIJMOD
    DO NNIJ=1,NIJMOD
      READ(1,*)ITMP,JTMP
      LIJMOD(NNIJ)=LIJ(ITMP,JTMP)
    ENDDO
    CLOSE(1)
    DO L=2,LA
      LSMOD=1
      DO NNIJ=1,NIJMOD
        IF( L == LIJMOD(NNIJ) ) LSMOD=0
      ENDDO
      IF( LSMOD == 1 )THEN
        WRITE(99,906)HP(L),H1P(L),HWQ(L),H2WQ(L),BELV(L)
        WRITE(99,906)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
        WRITE(99,906)(U(L,K),K=1,KC)
        WRITE(99,906)(U1(L,K),K=1,KC)
        WRITE(99,906)(V(L,K),K=1,KC)
        WRITE(99,906)(V1(L,K),K=1,KC)
        WRITE(99,906)(QQ(L,K),K=0,KC)
        WRITE(99,906)(QQ1(L,K),K=0,KC)
        WRITE(99,906)(QQL(L,K),K=0,KC)
        WRITE(99,906)(QQL1(L,K),K=0,KC)
        WRITE(99,906)(DML(L,K),K=0,KC)
        IF( ISCO(1) == 1 )THEN
          WRITE(99,906)(SAL(L,K),K=1,KC)
          WRITE(99,906)(SAL1(L,K),K=1,KC)
        ENDIF
        IF( ISCO(2) == 1 )THEN
          WRITE(99,906)(TEM(L,K),K=1,KC)
          WRITE(99,906)(TEM1(L,K),K=1,KC)
        ENDIF
        IF( ISCO(3) == 1 )THEN
          WRITE(99,906)(DYE(L,K),K=1,KC)
          WRITE(99,906)(DYE1(L,K),K=1,KC)
        ENDIF
        IF( ISCO(4) == 1 )THEN
          WRITE(99,906)SEDB(L,1,1),(SED(L,K,1),K=1,KC)
          WRITE(99,906)SEDB1(L,1,1),(SED1(L,K,1),K=1,KC)
        ENDIF
        IF( ISCO(5) == 1 )THEN
          WRITE(99,906)(SFL(L,K),K=1,KC)
          WRITE(99,906)(SFL2(L,K),K=1,KC)
        ENDIF
      ENDIF
    ENDDO
    DO M=1,5
      IF( ISCO(M) == 1 )THEN
        DO LL=1,NCBS
          DO K=1,KC
            NLOS(LL,K,M)=MAX( NLOS(LL,K,M)-N , -NTSCRS(LL) )
          ENDDO
          WRITE(99,908)(NLOS(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOS(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBW
          DO K=1,KC
            NLOW(LL,K,M)=MAX( NLOW(LL,K,M)-N , -NTSCRW(LL) )
          ENDDO
          WRITE(99,908)(NLOW(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOW(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBE
          DO K=1,KC
            NLOE(LL,K,M)=MAX( NLOE(LL,K,M)-N , -NTSCRE(LL) )
          ENDDO
          WRITE(99,908)(NLOE(LL,K,M),K=1,KC)
          WRITE(99,906)(CLOE(LL,K,M),K=1,KC)
        ENDDO
        DO LL=1,NCBN
          DO K=1,KC
            NLON(LL,K,M)=MAX( NLON(LL,K,M)-N , -NTSCRN(LL) )
          ENDDO
          WRITE(99,908)(NLON(LL,K,M),K=1,KC)
          WRITE(99,906)(CLON(LL,K,M),K=1,KC)
        ENDDO
      ENDIF
    ENDDO
    CLOSE(99)
    906 FORMAT(50E17.8)
    908 FORMAT(50I10)
  END SUBROUTINE

  SUBROUTINE RESTIN1(OPT)
    ! ** ADDED CODE TO PROPERLY INITIAL RESTART INPUT FOR DRYING AND WETTING
    ! ** SUBROUTINE RESTIN1 READS A RESTART FILE EXPORTED BY EFDC_DSI
    ! ** USING TRUE RESATRT FILE IN EFDC_DSI (ISRESTI == 1/-1/11)
    ! ** ISRESTI=1/-1,VER=720: RESTART FILES (.INP) BY EFDC_7.2 
    ! **              VER=710: RESTART FILES (.INP) BY EFDC_7.0-7.1
    INTEGER,INTENT(IN),OPTIONAL :: OPT    
    INTEGER :: ISBELVC,L,K,NT,NS,M,LL,NMD,ITMP1,JTMP1,VER,NCOL,NCTL,ISCOCHK(0:10),IMASK
    INTEGER :: ITMP2,JTMP2,ITMP3,JTMP3,LS,LN,LDUM,IDUM,JDUM,IDFLAG,ICORDRY,IU,JU,LU,ID,JD,LD
    
    REAL :: BELTMP,TMPVAL,RDZC,HDRY2,DHPDT,SUBW,SUBE,SVBS,SVBN,RDRY
    REAL :: SUBL,SVBL
    
    REAL,ALLOCATABLE,DIMENSION(:) :: TDUMMY
    CHARACTER(100) :: STR

    ALLOCATE(TDUMMY(KCM))
    TDUMMY=0.
    RSTFIRST_WS  = 0
    RSTFIRST_VEL = 0
    RSTFIRST_WC  = 0
    RSTFIRST_WQ  = 0
    RSTFIRST_SEDZLJ = 0
    RSTFIRST_RPEM = 0
    RSTFIRST_BED  = 0
    RSTFIRST_SD   = 0 
    RSTFIRST_BC   = 0
    
    WRITE(*,'(A)')'READING RESTART FILE: RESTART.INP'
    OPEN(UINP,FILE='restart.inp',STATUS='UNKNOWN')
    ISBELVC=0
    
    READ(UINP,'(A)',ERR=1000) STR
    NCOL = NUMCOL(STR)
    IF( NCOL == 3 )THEN
      READ(STR,*,ERR=1000) NRESTART,TBEGINC,VER
    ELSE
      READ(STR,*,ERR=1000) NRESTART,TBEGINC
      VER = 710
    ENDIF
            
    IF( PRESENT(OPT) )THEN
      CLOSE(UINP)
      RETURN
    ENDIF
    
    IF( VER <= 710 .AND. ISTRAN(2) > 0 .AND. ISCI(2) >= 1 )THEN
      WRITE(*,'(A)')'READING RESTART FILE: TEMP.RST'
      OPEN(111,FILE='temp.rst',STATUS='UNKNOWN')
      DO L=2,LA
        READ(111,*)LDUM,IDUM,JDUM,(TDUMMY(K),K=1,KC),TEMB(L)
      ENDDO
      CLOSE(111)
      
    ELSEIF( VER >= 720 )THEN
      READ(UINP,*) ISCOCHK(0:8)
      ! ** RESET ISCI(I) BASED ON RESTART.INP
      !ISCI = ISCOCHK
    ENDIF
    
    DO L=2,LA
      IF( ISRESTIOPT == 0 )THEN
        READ(UINP,*,ERR=1001)HP(L),H1P(L),HWQ(L),H2WQ(L),BELV(L)
      ELSE
        READ(UINP,*,ERR=1002)HP(L),H1P(L),HWQ(L),H2WQ(L),BELTMP
        IF( BELTMP /= BELV(L) )THEN
          ISBELVC=1
          WRITE(6,600)IL(L),JL(L),BELTMP,BELV(L)
          WRITE(8,600)IL(L),JL(L),BELTMP,BELV(L)
          HP(L)   = HP(L)+BELTMP-BELV(L)
          H1P(L)  = H1P(L)+BELTMP-BELV(L)
          HWQ(L)  = HWQ(L)+BELTMP-BELV(L)
          H2WQ(L) = H2WQ(L)+BELTMP-BELV(L)
        ENDIF
      ENDIF
      IF( HP(L) < 0.0 .OR. H1P(L) < 0.0 )THEN
        WRITE(6,9696)L,IL(L),JL(L),HP(L),H1P(L)
        STOP
      ELSEIF( HP(L) < 0.0001 .OR. H1P(L) < 0.0001 )THEN
        WRITE(6,9698)L,IL(L),JL(L),HP(L),H1P(L)
        HP(L) = HDRY*0.9
        H1P(L)= HDRY*0.9
      ENDIF
      READ(UINP,*,ERR=1003)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
      READ(UINP,*,ERR=1004)(U(L,K),K=1,KC)
      READ(UINP,*,ERR=1005)(U1(L,K),K=1,KC)
      READ(UINP,*,ERR=1006)(V(L,K),K=1,KC)
      READ(UINP,*,ERR=1007)(V1(L,K),K=1,KC)
      READ(UINP,*,ERR=1008)(QQ(L,K),K=0,KC)
      READ(UINP,*,ERR=1009)(QQ1(L,K),K=0,KC)
      READ(UINP,*,ERR=1010)(QQL(L,K),K=0,KC)
      READ(UINP,*,ERR=1011)(QQL1(L,K),K=0,KC)
      READ(UINP,*,ERR=1012)(DML(L,K),K=0,KC)

      IF( ISCOCHK(2) == 1 .AND. VER >= 720 .AND. ISICE > 2 )THEN
        IF( ISTRAN(2) > 0 .AND. ISCI(2) > 0 )THEN
          READ(UINP,*,ERR=1035)ICETHICK(L),ICETEMP(L)       ! ** EFDC_7.3
          IF( ICETHICK(L) > 0.0 ) ICECELL(L) = .TRUE.
        ELSE
          READ(UINP,*,ERR=1035)TMPVAL,TMPVAL
        ENDIF
      ENDIF        

      ! *** REMOVED THE ISTRAN CHECK SINCE IT IS NOT REQUIRED JUST TO LOAD THE RESTART FILE.  (EE7.2)
      IF( ISCOCHK(1) == 1 )THEN
        IF( ISTRAN(1) > 0 .AND. ISCI(1) > 0 )THEN
          READ(UINP,*,ERR=1013)(SAL(L,K),K=1,KC)
          READ(UINP,*,ERR=1014)(SAL1(L,K),K=1,KC)
        ELSE
          READ(UINP,*,ERR=1013)(TMPVAL,K=1,KC)
          READ(UINP,*,ERR=1014)(TMPVAL,K=1,KC)
        ENDIF
      ENDIF 
      IF( ISCOCHK(2) == 1 )THEN
        IF( ISTRAN(2) > 0 .AND. ISCI(2) > 0 )THEN
          IF( VER >= 720 )THEN
            READ(UINP,*,ERR=1015) TEMB(L),(TEM(L,K),K=1,KC)  ! ** EFDC_7.2
          ELSEIF (VER <= 710 )THEN      
            READ(UINP,*,ERR=1015) (TEM(L,K),K=1,KC)          ! ** EFDC_7.0-7.1
          ENDIF
          READ(UINP,*,ERR=1016)(TEM1(L,K),K=1,KC)
        ELSE
          IF( VER >= 720 )THEN
            READ(UINP,*,ERR=1015) TMPVAL,(TMPVAL,K=1,KC)
          ELSEIF (VER <= 710 )THEN      
            READ(UINP,*,ERR=1015) (TMPVAL,K=1,KC) 
          ENDIF
          READ(UINP,*,ERR=1016)(TMPVAL,K=1,KC)
        ENDIF
      ENDIF
      IF( ISCOCHK(3) == 1 )THEN
        IF( ISTRAN(3) > 0 .AND. ISCI(3) > 0 )THEN
          READ(UINP,*,ERR=1017)(DYE(L,K),K=1,KC)
          READ(UINP,*,ERR=1018)(DYE1(L,K),K=1,KC)
        ELSE
          READ(UINP,*,ERR=1017)(TMPVAL,K=1,KC)
          READ(UINP,*,ERR=1018)(TMPVAL,K=1,KC)
        ENDIF
      ENDIF
      IF( ISCOCHK(4) == 1 )THEN
        IF( ISTRAN(4) > 0 .AND. ISCI(4) > 0 )THEN
          READ(UINP,*,ERR=1021)SFLSBOT(L),(SFL(L,K),K=1,KC)
          READ(UINP,*,ERR=1022)SFLSBOT(L),(SFL2(L,K),K=1,KC)
        ELSE
          READ(UINP,*,ERR=1021)TMPVAL,(TMPVAL,K=1,KC)
          READ(UINP,*,ERR=1022)TMPVAL,(TMPVAL,K=1,KC)
        ENDIF
      ENDIF
      IF( ISCOCHK(5) == 1 )THEN
        IF( ISTRAN(5) > 0 .AND. ISCI(5) > 0 )THEN
          DO NT=1,NTXM
            READ(UINP,*,ERR=1019)(TOXB(L,K,NT),K=1,KB)
            READ(UINP,*,ERR=1019)(TOX(L,K,NT),K=1,KC)
            READ(UINP,*,ERR=1020)(TOXB1(L,K,NT),K=1,KB)
            READ(UINP,*,ERR=1020)(TOX1(L,K,NT),K=1,KC)
          ENDDO
        ELSE
          DO NT=1,NTXM
            READ(UINP,*,ERR=1019)(TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1019)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1020)(TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1020)(TMPVAL,K=1,KC)
          ENDDO
        ENDIF
      ENDIF
      IF( ISCOCHK(6) == 1 )THEN
        IF( ISTRAN(6) > 0 .AND. ISCI(6) > 0 )THEN
          DO NS=1,NSCM
            READ(UINP,*,ERR=1019)(SEDB(L,K,NS),K=1,KB)
            READ(UINP,*,ERR=1019)(SED(L,K,NS),K=1,KC)
            READ(UINP,*,ERR=1020)(SEDB1(L,K,NS),K=1,KB)
            READ(UINP,*,ERR=1020)(SED1(L,K,NS),K=1,KC)
          ENDDO
        ELSE
          DO NS=1,NSCM
            READ(UINP,*,ERR=1019)(TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1019)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1020)(TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1020)(TMPVAL,K=1,KC)
          ENDDO
        ENDIF
      ENDIF
      IF( ISCOCHK(7) == 1 )THEN
        IF( ISTRAN(7) > 0 .AND. ISCI(7) > 0 )THEN
          DO NS=1,NSNM
            READ(UINP,*,ERR=1019)(SNDB(L,K,NS),K=1,KB)
            READ(UINP,*,ERR=1019)(SND(L,K,NS),K=1,KC)
            READ(UINP,*,ERR=1020)(SNDB1(L,K,NS),K=1,KB)
            READ(UINP,*,ERR=1020)(SND1(L,K,NS),K=1,KC)
          ENDDO
        ELSE
          DO NS=1,NSNM
            READ(UINP,*,ERR=1019)(TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1019)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1020)(TMPVAL,K=1,KB)
            READ(UINP,*,ERR=1020)(TMPVAL,K=1,KC)
          ENDDO
        ENDIF
      ENDIF
      IF( (ISCOCHK(6) == 1 ) .OR. (ISCOCHK(7) == 1 ) )THEN
        IF( ( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 ) .AND. ( ISCI(6) > 0 .OR. ISCI(7) > 0 ) )THEN
          READ(UINP,*,ERR=1019)(HBED(L,K),K=1,KB)
          READ(UINP,*,ERR=1019)(HBED1(L,K),K=1,KB)
          READ(UINP,*,ERR=1019)(VDRBED(L,K),K=1,KB)
          READ(UINP,*,ERR=1019)(VDRBED1(L,K),K=1,KB)
        ELSE
          READ(UINP,*,ERR=1019)(TMPVAL,K=1,KB)
          READ(UINP,*,ERR=1019)(TMPVAL,K=1,KB)
          READ(UINP,*,ERR=1019)(TMPVAL,K=1,KB)
          READ(UINP,*,ERR=1019)(TMPVAL,K=1,KB)
        ENDIF
      ENDIF
    ENDDO

    ! *** BOUNDARY CONDITIONS
    DO M=1,4
      IF( ISCOCHK(M) == 1 )THEN
        IF( ISTRAN(M) == 1 .AND. ISCI(M) > 0 )THEN
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
          ENDDO
        ELSE
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1024)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1026)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1028)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1030)(TMPVAL,K=1,KC)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    
    IF( ISCOCHK(5) == 1 )THEN
      IF( ISTRAN(5) == 1 .AND. ISCI(5) > 0 )THEN
        DO NT=1,NTXM
          M=MSVTOX(NT)
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
          ENDDO
        ENDDO
      ELSE
        DO NT=1,NTXM
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1024)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1026)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1028)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1030)(TMPVAL,K=1,KC)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    IF( ISCOCHK(6) == 1 )THEN
      IF( ISTRAN(6) == 1 .AND. ISCI(6) > 0 )THEN
        DO NT=1,NSCM
          M=MSVSED(NT)
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
          ENDDO
        ENDDO
      ELSE
        DO NT=1,NSCM
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1024)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1026)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1028)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1030)(TMPVAL,K=1,KC)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    IF( ISCOCHK(7) == 1 )THEN
      IF( ISTRAN(7) == 1 .AND. ISCI(7) > 0 )THEN
        DO NT=1,NSNM
          M=MSVSND(NT)
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(NLOS(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1024)(CLOS(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(NLOW(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1026)(CLOW(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(NLOE(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1028)(CLOE(LL,K,M),K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(NLON(LL,K,M),K=1,KC)
            READ(UINP,*,ERR=1030)(CLON(LL,K,M),K=1,KC)
          ENDDO
        ENDDO
      ELSE
        DO NT=1,NSNM
          DO LL=1,NCBS
            READ(UINP,*,ERR=1023)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1024)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBW
            READ(UINP,*,ERR=1025)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1026)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBE
            READ(UINP,*,ERR=1027)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1028)(TMPVAL,K=1,KC)
          ENDDO
          DO LL=1,NCBN
            READ(UINP,*,ERR=1029)(TMPVAL,K=1,KC)
            READ(UINP,*,ERR=1030)(TMPVAL,K=1,KC)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    
    DO L=2,LA
      READ(UINP,*,ERR=1031)QSUME(L),(QSUM(L,K),K=1,KC)
    ENDDO
    IF( MDCHH >= 1 )THEN
      DO NMD=1,MDCHH
        READ(UINP,*,ERR=1032)ITMP1,JTMP1,ITMP2,JTMP2,ITMP3,JTMP3,QCHANU(NMD),QCHANV(NMD)
      ENDDO
    ELSE
      DO NMD=1,MDCHH
        QCHANU(NMD)=0.
        QCHANV(NMD)=0.
        QCHANUN(NMD)=0.
        QCHANVN(NMD)=0.
      ENDDO
    ENDIF
    IF( ISGWIE >= 1 )THEN
      DO L=2,LA
        READ(UINP,*,ERR=1033) AGWELV(L),AGWELV1(L)
      ENDDO
    ENDIF
    ! *** UPSTREAM OR UPSTREAM-DOWNSTREAM DEPENDENT RATING CURVE WITH A LOWER CHORD ACTIVATION TOGGLE
    DO NCTL=1,NQCTL
      IF( NQCTYP(NCTL) == 3 .OR. NQCTYP(NCTL) == 4 )THEN
        ! *** UPSTREAM DEPTH (3) OR ELEVATION/STAGE (4) DEPENDANT FLOWS
        READ(1,*,ERR=1034)LL,NLOWCHORD(NCTL)
        READ(1,*,ERR=1034)SAVESUB(1,NCTL),SAVESVB(1,NCTL),SAVESUB(2,NCTL),SAVESVB(2,NCTL)
        READ(1,*,ERR=1034)LOWCHORDU(NCTL),LOWCHORDV(NCTL)

        IU=IQCTLU(NCTL)
        JU=JQCTLU(NCTL)
        LU=LIJ(IU,JU)

        ID=IQCTLD(NCTL)
        JD=JQCTLD(NCTL)
        IF( ID == 0 .AND. JD == 0 )THEN
          ! *** INVALID SPECIFICATION
          CYCLE
        ENDIF
        LD=LIJ(ID,JD)

        IF( NLOWCHORD(NCTL) < 1 )THEN
          IF( ID > IU )THEN
            IF( SAVESUB(2,NCTL) > 0.5 )THEN
              SUB(LD)   = 1.0
              SUBO(LD)  = 1.0
            ENDIF
          ELSE
            IF( SAVESUB(1,NCTL) > 0.5 )THEN
              SUB(LU)   = 1.0
              SUBO(LU)  = 1.0
            ENDIF
          ENDIF
          IF( JD > JU )THEN
            IF( SAVESVB(2,NCTL) > 0.5 )THEN
              SVB(LD)   = 1.0
              SVBO(LD)  = 1.0
            ENDIF
          ELSE
            IF( SAVESVB(1,NCTL) > 0.5 )THEN
              SVB(LU)=1.0
              SVBO(LU)=1.0
            ENDIF
          ENDIF
        ELSE
          ! *** SET THE CELL FACE FLAGS - BLOCK FULL HYDRODYNAMICS
          IF( ID > IU )THEN
            SUB(LD)=0.0
            SUBO(LD)=0.0
            UHDYE(LD)=0.0
          ELSE
            SUB(LU)=0.0
            SUBO(LU)=0.0
            UHDYE(LU)=0.0
          ENDIF
          IF( JD > JU )THEN
            SVB(LD)=0.0
            SVBO(LD)=0.0
            VHDXE(LD)=0.0
          ELSE
            SVB(LU)=0.0
            SVBO(LU)=0.0
            VHDXE(LU)=0.0
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    ! *** CLOSE MAIN RESTART.OUT FILE
    CLOSE(UINP)

    6666 FORMAT(3I10,F12.6)
    6667 FORMAT(7I5,2X,E12.4,2X,E12.4)
    DO K=1,KC
      SAL(1,K)=0.
      TEM(1,K)=0.
      DYE(1,K)=0.
      IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
        SED(1,K,1)=0.
        SND(1,K,1)=0.
        SED1(1,K,1)=0.
        SND1(1,K,1)=0.
        SED(LC,K,1)=0.
        SND(LC,K,1)=0.
        SED1(LC,K,1)=0.
        SND1(LC,K,1)=0.
        IF( ISTRAN(5) > 0 )THEN
          TOX(1,K,1)=0.
          TOX1(1,K,1)=0.
          TOX(LC,K,1)=0.
          TOX1(LC,K,1)=0.
        ENDIF
      ENDIF
      SFL(1,K)=0.
      VHDX(1,K)=0.
      UHDY(1,K)=0.
      VHDXF(1,K)=0.
      UHDYF(1,K)=0.
      SAL1(1,K)=0.
      TEM1(1,K)=0.
      DYE1(1,K)=0.
      SFL2(1,K)=0.
      VHDX1(1,K)=0.
      UHDY1(1,K)=0.
      VHDXF1(1,K)=0.
      UHDYF1(1,K)=0.
      SAL(LC,K)=0.
      TEM(LC,K)=0.
      DYE(LC,K)=0.
      SFL(LC,K)=0.
      VHDX(LC,K)=0.
      UHDY(LC,K)=0.
      VHDXF(LC,K)=0.
      UHDYF(LC,K)=0.
      SAL1(LC,K)=0.
      TEM1(LC,K)=0.
      DYE1(LC,K)=0.
      SFL2(LC,K)=0.
      VHDX1(LC,K)=0.
      UHDY1(LC,K)=0.
      VHDXF1(LC,K)=0.
      UHDYF1(LC,K)=0.
    ENDDO
    IF( ISWQFLUX == 1 )THEN
      DO K=1,KC
        UHDYWQ(1,K)=0.
        UHDYWQ(LC,K)=0.
        VHDXWQ(1,K)=0.
        VHDXWQ(LC,K)=0.
      ENDDO
    ENDIF
    
    DO L=2,LA
      UHDYE(L)=SUB(L)*UHDYE(L)
      UHDY1E(L)=SUB(L)*UHDY1E(L)
      VHDXE(L)=SVB(L)*VHDXE(L)
      VHDX1E(L)=SVB(L)*VHDX1E(L)
    ENDDO
    DO K=1,KC
      DO L=2,LA
        U(L,K)=SUB(L)*U(L,K)
        U1(L,K)=SUB(L)*U1(L,K)
        V(L,K)=SVB(L)*V(L,K)
        V1(L,K)=SVB(L)*V1(L,K)
      ENDDO
    ENDDO
    DO L=2,LA
      LS=LSC(L)
      H1U(L)=0.5*(DXP(L)*DYP(L)*H1P(L)+DXP(LWC(L))*DYP(LWC(L))*H1P(LWC(L)))/(DXU(L)*DYU(L))
      H1V(L)=0.5*(DXP(L)*DYP(L)*H1P(L)+DXP(LS )*DYP(LWC(L))*H1P(LS ))/(DXV(L)*DYV(L))
      P1(L)=G*(H1P(L)+BELV(L))
      HU(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(LWC(L))*DYP(LWC(L))*HP(LWC(L)))/(DXU(L)*DYU(L))
      HV(L)=0.5*(DXP(L)*DYP(L)*HP(L)+DXP(LS )*DYP(LWC(L))*HP(LS ))/(DXV(L)*DYV(L))
      P(L)=G*(HP(L)+BELV(L))
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      H1UI(L)=1./H1U(L)
      H1VI(L)=1./H1V(L)
    ENDDO
    H1U(1)=H1U(2)
    H1V(1)=H1V(2)
    P1(1)=P1(2)
    HU(1)=HU(2)
    HV(1)=HV(2)
    P(1)=P(2)
    HPI(1)=1./HP(2)
    HUI(1)=1./HU(2)
    HVI(1)=1./HV(2)
    H1UI(1)=1./H1U(2)
    H1VI(1)=1./H1V(2)
    H1U(LC)=H1U(LA)
    H1V(LC)=H1V(LA)
    P1(LC)=P1(LA)
    HU(LC)=HU(LA)
    HV(LC)=HV(LA)
    P(LC)=P(LA)
    HPI(LC)=1./HP(LA)
    HUI(LC)=1./HU(LA)
    HVI(LC)=1./HV(LA)
    H1UI(LC)=1./H1U(LA)
    H1VI(LC)=1./H1V(LA)
    
    DO K=1,KC
      DO L=2,LA
        UHDYF1(L,K) = DYU(L)*H1U(L)*U1(L,K)
        VHDXF1(L,K) = DXV(L)*H1V(L)*V1(L,K)
        UHDY1(L,K)  = UHDYF1(L,K)*DZC(L,K)
        VHDX1(L,K)  = VHDXF1(L,K)*DZC(L,K)
        
        UHDYF(L,K)  = DYU(L)*HU(L)*U(L,K)
        VHDXF(L,K)  = DXV(L)*HV(L)*V(L,K)
        UHDY(L,K)   = UHDYF(L,K)*DZC(L,K)
        VHDX(L,K)   = VHDXF(L,K)*DZC(L,K)
        
        SAL(L,K)=MAX(SAL(L,K),0.)
        SAL1(L,K)=MAX(SAL1(L,K),0.)
      ENDDO
    ENDDO

    ! **  CORRECT FOR CHANGED BOTTOM ELEV
    IF( ISRESTIOPT == 0 .AND. ISBELVC == 1 )THEN
      DO L=2,LA
        UHE(L)=0.
        VHE(L)=0.
      ENDDO
      DO K=1,KC
        DO L=2,LA
          UHE(L)=UHE(L)+UHDYF1(L,K)
          VHE(L)=VHE(L)+VHDXF1(L,K)
        ENDDO
      ENDDO
      DO L=2,LA
        IF( UHE(L) /= 0. )THEN
          TMPVAL=UHDY1E(L)/UHE(L)
          DO K=1,KC
            U1(L,K)=TMPVAL*U1(L,K)
          ENDDO
        ENDIF
        IF( VHE(L) /= 0. )THEN
          TMPVAL=VHDX1E(L)/VHE(L)
          DO K=1,KC
            V1(L,K)=TMPVAL*V1(L,K)
          ENDDO
        ENDIF
      ENDDO
      DO L=2,LA
        UHE(L)=0.
        VHE(L)=0.
      ENDDO
      DO K=1,KC
        DO L=2,LA
          UHE(L)=UHE(L)+UHDYF(L,K)
          VHE(L)=VHE(L)+VHDXF(L,K)
        ENDDO
      ENDDO
      
      DO L=2,LA
        IF( UHE(L) /= 0. )THEN
          TMPVAL=UHDYE(L)/UHE(L)
          DO K=1,KC
            U(L,K)=TMPVAL*U(L,K)
          ENDDO
        ENDIF
        IF( VHE(L) /= 0. )THEN
          TMPVAL=VHDXE(L)/VHE(L)
          DO K=1,KC
            V(L,K)=TMPVAL*V(L,K)
          ENDDO
        ENDIF
      ENDDO
      DO K=1,KC
        DO L=2,LA
          UHDYF1(L,K) = DYU(L)*H1U(L)*U1(L,K)
          VHDXF1(L,K) = DXV(L)*H1V(L)*V1(L,K)
          UHDY1(L,K)  = UHDYF1(L,K)*DZC(L,K)
          VHDX1(L,K)  = VHDXF1(L,K)*DZC(L,K)

          UHDYF(L,K) = DYU(L)*HU(L)*U(L,K)
          VHDXF(L,K) = DXV(L)*HV(L)*V(L,K)
          UHDY(L,K)  = UHDYF(L,K)*DZC(L,K)
          VHDX(L,K)  = VHDXF(L,K)*DZC(L,K)
        ENDDO
      ENDDO
    ENDIF
    N=0
    
    IF( ISDRY == 0 )THEN
      CALL CALTSXY
      CALL CALQVS (2)
    ENDIF
    
    DO K=1,KS
      DO L=2,LA
        IF( LKSZ(L,K) )CYCLE
        RDZC=DZC(L,K)
        LN=LNC(L)
        W(L,K)  = SWB(L)*( W(L,K-1) - RDZC*(UHDYF(LEC(L),K)-UHDYF(L,K)-UHDYE(LEC(L))+UHDYE(L) &
                                          + VHDXF(LN,K) -VHDXF(L,K)-VHDXE(LN) +VHDXE(L))*DXYIP(L)) & 
                + SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
        
        W1(L,K) = SWB(L)*( W1(L,K-1) - RDZC*(UHDYF1(LEC(L),K)-UHDYF1(L,K)-UHDY1E(LEC(L))+UHDY1E(L) &
                                           + VHDXF1(LN,K) -VHDXF1(L,K)-VHDX1E(LN) +VHDX1E(L))*DXYIP(L)) &
                + SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
      ENDDO
    ENDDO
    !
    ! **  SET DRYING AND WETTING FLAGS
    IF( ISDRY > 0 .AND. ISDRY < 97 )THEN
      DO L=2,LA
        ISCDRY(L)=0
        LS=LSC(L)
        LN=LNC(L)
        IF( HP(L) <= HDRY )THEN
          ISCDRY(L)=1
          SUB(L)=0.
          SVB(L)=0.
          SUB(LEC(L))=0.
          SVB(LN)=0.
          SBX(L)=0.
          SBY(L)=0.
          SBX(LEC(L))=0.
          SBY(LN)=0.
        ENDIF
      ENDDO
    ENDIF
    IF( ISDRY == 98 )THEN
      HDRY2=2.*HDRY
      DO L=2,LA
        LS=LSC(L)
        LN=LNC(L)
        IF( HP(L) <= HDRY )THEN
          DHPDT=(HP(L)-H1P(L))/DT
          IF( DHPDT > 0.0 )THEN
            SUBW=SUB(L)
            SUBE=SUB(LEC(L))
            SVBS=SVB(L)
            SVBN=SVB(LN)
            SUB(L)=0.0
            SUB(LEC(L))=0.0
            SVB(L)=0.0
            SVB(LN)=0.0
            SBX(L)=0.0
            SBX(LEC(L))=0.0
            SBY(L)=0.0
            SBY(LN)=0.0
            IF( SUBO(L) > 0.5 )THEN
              IF( UHDYE(L) > 0.0 .AND. HP(LWC(L)) > HDRY2 )THEN
                SUB(L)=1.
                SBX(L)=1.
              ENDIF
            ENDIF
            IF( SUBO(LEC(L)) > 0.5 )THEN
              IF( UHDYE(LEC(L)) < 0.0 .AND. HP(LEC(L)) > HDRY2 )THEN
                SUB(LEC(L))=1.
                SBX(LEC(L))=1.
              ENDIF
            ENDIF
            IF( SVBO(L) > 0.5 )THEN
              IF( VHDXE(L) > 0.0 .AND. HP(LS) > HDRY2 )THEN
                SVB(L)=1.
                SBY(L)=1.
              ENDIF
            ENDIF
            IF( SVBO(LN) > 0.5 )THEN
              IF( VHDXE(LN) < 0.0 .AND. HP(LN) > HDRY2 )THEN
                SVB(LN)=1.
                SBY(LN)=1.
              ENDIF
            ENDIF
            RDRY=SUB(L)+SUB(LEC(L))+SVB(L)+SVB(LN)
            IF( RDRY < 0.5 )THEN
              ISCDRY(L)=1
            ELSE
              ISCDRY(L)=0
            ENDIF
            IDFLAG=0
            TMPVAL=ABS(SUB(L)-SUBW)
            IF( TMPVAL > 0.5 )ICORDRY=1
            TMPVAL=ABS(SUB(LEC(L))-SUBE)
            IF( TMPVAL > 0.5 )ICORDRY=1
            TMPVAL=ABS(SVB(L)-SVBS)
            IF( TMPVAL > 0.5 )ICORDRY=1
            TMPVAL=ABS(SVB(LN)-SVBN)
            IF( TMPVAL > 0.5 )ICORDRY=1
          ELSE
            SUB(L)=0.0
            SUB(LEC(L))=0.0
            SVB(L)=0.0
            SVB(LN)=0.0
            SBX(L)=0.0
            SBX(LEC(L))=0.0
            SBY(L)=0.0
            SBY(LN)=0.0
            IF( ISCDRY(L) == 0 )THEN
              ISCDRY(L)=1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    
    IF( ISDRY == 99 )THEN
      WRITE(*,'(A)')'READING RESTART FILE: RSTWD.INP'
      OPEN(1,FILE='rstwd.inp',STATUS='UNKNOWN')
      OPEN(2,FILE=OUTDIR//'rstwd.rck',STATUS='UNKNOWN')
      CLOSE(2,STATUS='DELETE')
      OPEN(2,FILE=OUTDIR//'rstwd.rck',STATUS='UNKNOWN')
      DO L=2,LA
        READ(1,*)LDUM,IDUM,JDUM,ISCDRY(L),NATDRY(L),IMASK,SUBL,SVBL
        WRITE(2,'(6I10,4F7.3)')LDUM,IDUM,JDUM,ISCDRY(L),NATDRY(L),IMASK,SUBL,SVBL,SUBO(L),SVBO(L)
        IF( IMASK == 0 )THEN
          LMASKDRY(L) = .TRUE.
        ELSE
          LMASKDRY(L) = .FALSE.
        ENDIF          
        
        ! *** ALLOW FOR BREAKPOINT CHECKS
        IF( SUBL /= SUB(L) )THEN
          SUB(L) = SUBL
        ENDIF
        IF( SVBL /= SVB(L) )THEN
          SVB(L) = SVBL
        ENDIF
      ENDDO
      
      CLOSE(1)
      CLOSE(2)
    ENDIF
    
    ! *** DS-INTL BEGIN BLOCK
    !DO L=2,LA
    !  IF( IMASKDRY(L) == 0 ) LMASKDRY(L)=.TRUE.
    !  IF( IMASKDRY(L) > 0 ) LMASKDRY(L)=.FALSE.
    !END DO

    RETURN

    101 FORMAT(I5)
    102 FORMAT(3I5,12F8.2)

    ! **  WRITE READ ERRORS ON RESTART
    1000 WRITE(6,2000)
    STOP

    1001 WRITE(6,2001)L
    STOP

    1002 WRITE(6,2002)L
    STOP

    1003 WRITE(6,2003)L
    STOP

    1004 WRITE(6,2004)L
    STOP

    1005 WRITE(6,2005)L
    STOP

    1006 WRITE(6,2006)L
    STOP

    1007 WRITE(6,2007)L
    STOP

    1008 WRITE(6,2008)L
    STOP

    1009 WRITE(6,2009)L
    STOP

    1010 WRITE(6,2010)L
    STOP

    1011 WRITE(6,2011)L
    STOP

    1012 WRITE(6,2012)L
    STOP

    1013 WRITE(6,2013)L
    STOP

    1014 WRITE(6,2014)L
    STOP

    1015 WRITE(6,2015)L
    STOP

    1016 WRITE(6,2016)L
    STOP

    1017 WRITE(6,2017)L
    STOP

    1018 WRITE(6,2018)L
    STOP

    1019 WRITE(6,2019)L
    STOP

    1020 WRITE(6,2020)L
    STOP

    1021 WRITE(6,2021)L
    STOP

    1022 WRITE(6,2022)L
    STOP

    1023 WRITE(6,2023)L
    STOP

    1024 WRITE(6,2024)L
    STOP

    1025 WRITE(6,2025)L
    STOP

    1026 WRITE(6,2026)L
    STOP

    1027 WRITE(6,2027)L
    STOP

    1028 WRITE(6,2028)L
    STOP

    1029 WRITE(6,2029)L
    STOP

    1030 WRITE(6,2030)L
    STOP

    1031 WRITE(6,2031)L
    STOP

    1032 WRITE(6,2032)NMD
    STOP

    1033 WRITE(6,2033)L
    STOP

    1034 WRITE(6,2034)NCTL
    STOP

    1035 WRITE(6,2035)L
    STOP

    600 FORMAT(2X,'I,J,BELVOLD,BELVNEW',2I5,2F12.2)
    906 FORMAT(5E15.7)
    907 FORMAT(12E12.4)
    908 FORMAT(12I10)
    2000 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1000')
    2001 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1001 L =',I6)
    2002 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1002 L =',I6)
    2003 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1003 L =',I6)
    2004 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1004 L =',I6)
    2005 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1005 L =',I6)
    2006 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1006 L =',I6)
    2007 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1007 L =',I6)
    2008 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1008 L =',I6)
    2009 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1009 L =',I6)
    2010 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1010 L =',I6)
    2011 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1011 L =',I6)
    2012 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1012 L =',I6)
    2013 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1013 L =',I6)
    2014 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1014 L =',I6)
    2015 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1015 L =',I6)
    2016 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1016 L =',I6)
    2017 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1017 L =',I6)
    2018 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1018 L =',I6)
    2019 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1019 L =',I6)
    2020 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1020 L =',I6)
    2021 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1021 L =',I6)
    2022 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1022 L =',I6)
    2023 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1023 L =',I6)
    2024 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1024 L =',I6)
    2025 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1025 L =',I6)
    2026 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1026 L =',I6)
    2027 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1027 L =',I6)
    2028 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1028 L =',I6)
    2029 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1029 L =',I6)
    2030 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1030 L =',I6)
    2031 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1031 L =',I6)
    2032 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1032 NMD =',I6)
    2033 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1033 L =',I6)
    2034 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1034 NCTL =',I6)
    2035 FORMAT('  READ ERROR ON FILE RESTART.INP ERR 1035 L =',I6)
    9696 FORMAT('  NEGATIVE DEPTH RESTART, L,I,J,HP,H1P = ',3I7,2F10.4)
    9698 FORMAT('  ZERO DEPTH RESTART (WARN), L,I,J,HP,H1P = ',3I7,2F10.4)

  END SUBROUTINE

  SUBROUTINE RESTIN2
    ! ** READ AN OLD RESTART.INP FILE (ISRESTI == 2 )
    ! ** SUBROUTINE RESTINP READS A RESTART FILE FOR (KC/2) LAYERS AND
    ! ** AND INITIALIZES FOR KC LAYERS
    INTEGER(IK4) :: K,M,LL,L,LS,LN
    INTEGER(IK8) :: NREST
    REAL :: RDZC

    WRITE(*,'(A)')'READING RESTIN2 FILE: RESTART.INP'
    OPEN(1,FILE='restart.inp',STATUS='UNKNOWN')
    
    READ(1,*,ERR=1000)NREST
    DO L=2,LA
      READ(1,*,ERR=1000)HP(L),H1P(L),HWQ(L),H2WQ(L)
      READ(1,*,ERR=1000)UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
      READ(1,*,ERR=1000)(U(L,K),K=2,KC,2)
      READ(1,*,ERR=1000)(U1(L,K),K=2,KC,2)
      READ(1,*,ERR=1000)(V(L,K),K=2,KC,2)
      READ(1,*,ERR=1000)(V1(L,K),K=2,KC,2)
      READ(1,*,ERR=1000)(QQ(L,K),K=0,KC,2)
      READ(1,*,ERR=1000)(QQ1(L,K),K=0,KC,2)
      READ(1,*,ERR=1000)(QQL(L,K),K=0,KC,2)
      READ(1,*,ERR=1000)(QQL1(L,K),K=0,KC,2)
      READ(1,*,ERR=1000)(DML(L,K),K=0,KC,2)
      DO K=1,KS,2
        U(L,K)=U(L,K+1)
        U1(L,K)=U1(L,K+1)
        V(L,K)=V(L,K+1)
        V1(L,K)=V1(L,K+1)
        QQ(L,K)=0.5*(QQ(L,K-1)+QQ(L,K+1))
        QQSQR(L,K)=SQRT(QQ(L,K))  ! *** DS-INTL
        QQ1(L,K)=0.5*(QQ1(L,K-1)+QQ1(L,K+1))
        QQL(L,K)=0.5*(QQL(L,K-1)+QQL(L,K+1))
        QQL1(L,K)=0.5*(QQL1(L,K-1)+QQL1(L,K+1))
        DML(L,K)=0.5*(DML(L,K-1)+DML(L,K+1))
      ENDDO
      IF( ISTRAN(1) >= 1 .AND. ISCI(1) == 1 )THEN
        READ(1,*,ERR=1000)(SAL(L,K),K=2,KC,2)
        READ(1,*,ERR=1000)(SAL1(L,K),K=2,KC,2)
        DO K=1,KS,2
          SAL(L,K)=SAL(L,K+1)
          SAL1(L,K)=SAL1(L,K+1)
        ENDDO
      ENDIF
      IF( ISTRAN(2) >= 1 .AND. ISCI(2) == 1 )THEN
        READ(1,*,ERR=1000)(TEM(L,K),K=2,KC,2)
        READ(1,*,ERR=1000)(TEM1(L,K),K=2,KC,2)
        DO K=1,KS,2
          TEM(L,K)=TEM(L,K+1)
          TEM1(L,K)=TEM1(L,K+1)
        ENDDO
      ENDIF
      IF( ISTRAN(3) >= 1 .AND. ISCI(3) == 1 )THEN
        READ(1,*,ERR=1000)(DYE(L,K),K=2,KC,2)
        READ(1,*,ERR=1000)(DYE1(L,K),K=2,KC,2)
        DO K=1,KS,2
          DYE(L,K)=DYE(L,K+1)
          DYE1(L,K)=DYE1(L,K+1)
        ENDDO
      ENDIF
      IF( ISTRAN(4) >= 1 .AND. ISCI(4) == 1 )THEN
        READ(1,*,ERR=1000)SEDB(L,1,1),(SED(L,K,1),K=2,KC,2)
        READ(1,*,ERR=1000)SEDB1(L,1,1),(SED1(L,K,1),K=2,KC,2)
        DO K=1,KS,2
          SED(L,K,1)=SED(L,K+1,1)
          SED1(L,K,1)=SED1(L,K+1,1)
        ENDDO
      ENDIF
      IF( ISTRAN(5) >= 1 .AND. ISCI(5) == 1 )THEN
        READ(1,*,ERR=1000)(SFL(L,K),K=2,KC,2)
        READ(1,*,ERR=1000)(SFL2(L,K),K=2,KC,2)
        DO K=1,KS,2
          SFL(L,K)=SFL(L,K+1)
          SFL2(L,K)=SFL2(L,K+1)
        ENDDO
      ENDIF
    ENDDO
    DO M=1,5
      IF( ISTRAN(M) >= 1 .AND. ISCI(M) == 1 )THEN
        DO LL=1,NCBS
          READ(1,*,ERR=1000)(NLOS(LL,K,M),K=2,KC,2)
          READ(1,*,ERR=1000)(CLOS(LL,K,M),K=2,KC,2)
        ENDDO
        DO LL=1,NCBW
          READ(1,*,ERR=1000)(NLOW(LL,K,M),K=2,KC,2)
          READ(1,*,ERR=1000)(CLOW(LL,K,M),K=2,KC,2)
        ENDDO
        DO LL=1,NCBE
          READ(1,*,ERR=1000)(NLOE(LL,K,M),K=2,KC,2)
          READ(1,*,ERR=1000)(CLOE(LL,K,M),K=2,KC,2)
        ENDDO
        DO LL=1,NCBN
          READ(1,*,ERR=1000)(NLON(LL,K,M),K=2,KC,2)
          READ(1,*,ERR=1000)(CLON(LL,K,M),K=2,KC,2)
        ENDDO
      ENDIF
    ENDDO
    DO M=1,5
      DO K=1,KS,2
        DO LL=1,NCBS
          NLOS(LL,K,M)=NLOS(LL,K+1,M)
          CLOS(LL,K,M)=CLOS(LL,K+1,M)
        ENDDO
        DO LL=1,NCBW
          NLOW(LL,K,M)=NLOW(LL,K+1,M)
          CLOW(LL,K,M)=CLOW(LL,K+1,M)
        ENDDO
        DO LL=1,NCBE
          NLOE(LL,K,M)=NLOE(LL,K+1,M)
          CLOE(LL,K,M)=CLOE(LL,K+1,M)
        ENDDO
        DO LL=1,NCBN
          NLON(LL,K,M)=NLON(LL,K+1,M)
          CLON(LL,K,M)=CLON(LL,K+1,M)
        ENDDO
      ENDDO
    ENDDO
    CLOSE(1)
    DO K=1,KC
      SAL(1,K)=0.
      TEM(1,K)=0.
      DYE(1,K)=0.
      SED(1,K,1)=0.
      SFL(1,K)=0.
      VHDX(1,K)=0.
      UHDY(1,K)=0.
      SAL1(1,K)=0.
      TEM1(1,K)=0.
      DYE1(1,K)=0.
      SED1(1,K,1)=0.
      SFL2(1,K)=0.
      VHDX1(1,K)=0.
      UHDY1(1,K)=0.
      VHDXWQ(1,K)=0.
      UHDYWQ(1,K)=0.
      SAL(LC,K)=0.
      TEM(LC,K)=0.
      DYE(LC,K)=0.
      SED(LC,K,1)=0.
      SFL(LC,K)=0.
      VHDX(LC,K)=0.
      UHDY(LC,K)=0.
      SAL1(LC,K)=0.
      TEM1(LC,K)=0.
      DYE1(LC,K)=0.
      SED1(LC,K,1)=0.
      SFL2(LC,K)=0.
      VHDX1(LC,K)=0.
      UHDY1(LC,K)=0.
      VHDXWQ(LC,K)=0.
      UHDYWQ(LC,K)=0.
    ENDDO
    DO L=2,LA
      LS=LSC(L)
      H1U(L)=0.5*(H1P(L)+H1P(LWC(L)))
      H1V(L)=0.5*(H1P(L)+H1P(LS))
      P1(L)=G*(H1P(L)+BELV(L))
      HU(L)=0.5*(HP(L)+HP(LWC(L)))
      HV(L)=0.5*(HP(L)+HP(LS))
      P(L)=G*(HP(L)+BELV(L))
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      H1UI(L)=1./H1U(L)
      H1VI(L)=1./H1V(L)
    ENDDO
    DO K=1,KC
      DO L=2,LA
        UHDYF1(L,K) = DYU(L)*H1U(L)*U1(L,K)
        VHDXF1(L,K) = DXV(L)*H1V(L)*V1(L,K)
        UHDY1(L,K)  = UHDYF1(L,K)*DZC(L,K)
        VHDX1(L,K)  = VHDXF1(L,K)*DZC(L,K)
        
        UHDYF(L,K)  = DYU(L)*HU(L)*U(L,K)
        VHDXF(L,K)  = DXV(L)*HV(L)*V(L,K)
        UHDY(L,K)   = UHDYF(L,K)*DZC(L,K)
        VHDX(L,K)   = VHDXF(L,K)*DZC(L,K)

        SAL(L,K)=MAX(SAL(L,K),0.)
        SAL1(L,K)=MAX(SAL1(L,K),0.)
      ENDDO
    ENDDO
    N=0
    CALL CALQVS (2)
    DO K=1,KS
      DO L=2,LA
        IF( LKSZ(L,K) )CYCLE
        RDZC=DZC(L,K)
        LN=LNC(L)
        W(L,K)=SWB(L)*(W(L,K-1) &
            -RDZC*(UHDYF(LEC(L),K)-UHDYF(L,K)-UHDYE(LEC(L))+UHDYE(L) &
            +VHDXF(LN,K)-VHDXF(L,K)-VHDXE(LN)+VHDXE(L))*DXYIP(L)) &
            +SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
        W1(L,K)=SWB(L)*(W1(L,K-1) &
            -RDZC*(UHDYF1(LEC(L),K)-UHDYF1(L,K)-UHDY1E(LEC(L))+UHDY1E(L) &
            +VHDXF1(LN,K)-VHDXF1(L,K)-VHDX1E(LN)+VHDX1E(L))*DXYIP(L)) &
            +SWB(L)*( QSUM(L,K)-RDZC*QSUME(L) )*DXYIP(L)
      ENDDO
    ENDDO
    !
    ! **  WRITE READ ERRORS ON RESTART
    !
    GOTO 1002
     1000 WRITE(6,1001)
     1001 FORMAT('  READ ERROR ON FILE RESTART.INP ')
    STOP
     1002 CONTINUE
      906 FORMAT(4E15.7)
      907 FORMAT(12E12.4)
      908 FORMAT(12I10)   

  END SUBROUTINE

  SUBROUTINE RESTIN10
    ! **  SUBROUTINE RESTINP READS A RESTART FILE GENERATED BY A
    ! **  PRE SEPTEMBER 8, 1992 VERSION OF EFDC.FOR
    INTEGER :: K,LL,L,LS
    INTEGER(IK8) :: NREST
    
    WRITE(*,'(A,I6)')'READING RESTIN10 FILE: RESTART.INP'
    OPEN(1,FILE='restart.inp',STATUS='UNKNOWN')
    READ(1,*,ERR=1000)NREST
    DO L=2,LA
      READ(1,*,ERR=1000)P(L),P1(L),UHDYE(L),UHDY1E(L),VHDXE(L),VHDX1E(L)
      READ(1,*,ERR=1000)(U(L,K),K=1,KC)
      READ(1,*,ERR=1000)(U1(L,K),K=1,KC)
      READ(1,*,ERR=1000)(V(L,K),K=1,KC)
      READ(1,*,ERR=1000)(V1(L,K),K=1,KC)
      READ(1,*,ERR=1000)(W(L,K),K=1,KS)
      READ(1,*,ERR=1000)(W1(L,K),K=1,KS)
      READ(1,*,ERR=1000)(QQ(L,K),K=0,KC)
      READ(1,*,ERR=1000)(QQ1(L,K),K=0,KC)
      READ(1,*,ERR=1000)(QQL(L,K),K=0,KC)
      READ(1,*,ERR=1000)(QQL1(L,K),K=0,KC)
      READ(1,*,ERR=1000)(DML(L,K),K=0,KC)
      IF( ISTRAN(1) >= 1 .AND. ISCI(1) == 1 )THEN
        READ(1,*,ERR=1000)(SAL(L,K),K=1,KC)
        READ(1,*,ERR=1000)(SAL1(L,K),K=1,KC)
      ENDIF
      IF( ISTRAN(2) >= 1 .AND. ISCI(2) == 1 )THEN
        READ(1,*,ERR=1000)(TEM(L,K),K=1,KC)
        READ(1,*,ERR=1000)(TEM1(L,K),K=1,KC)
      ENDIF
      IF( ISTRAN(3) >= 1 .AND. ISCI(3) == 1 )THEN
        READ(1,*,ERR=1000)(DYE(L,K),K=1,KC)
        READ(1,*,ERR=1000)(DYE1(L,K),K=1,KC)
      ENDIF
      IF( ISTRAN(4) >= 1 .AND. ISCI(4) == 1 )THEN
        READ(1,*,ERR=1000) SEDB(L,1,1),(SED(L,K,1),K=1,KC)
        READ(1,*,ERR=1000) SEDB1(L,1,1),(SED1(L,K,1),K=1,KC)
      ENDIF
      IF( ISTRAN(5) >= 1 .AND. ISCI(5) == 1 )THEN
        READ(1,*,ERR=1000)(SFL(L,K),K=1,KC)
        READ(1,*,ERR=1000)(SFL2(L,K),K=1,KC)
      ENDIF
    ENDDO
    IF( ISTRAN(1) >= 1 .AND. ISCI(1) == 1 )THEN
      DO LL=1,NCBS
        L=LCBS(LL)
        READ(1,*,ERR=1000)(NLOS(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(CLOS(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBW
        L=LCBW(LL)
        READ(1,*,ERR=1000)(NLOW(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(CLOW(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBE
        L=LCBE(LL)
        READ(1,*,ERR=1000)(NLOE(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(CLOE(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBN
        L=LCBN(LL)
        READ(1,*,ERR=1000)(NLON(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(CLON(LL,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SAL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SAL1(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
      ENDDO
    ENDIF
    IF( ISTRAN(2) >= 1 .AND. ISCI(2) == 1 )THEN
      DO LL=1,NCBS
        L=LCBS(LL)
        READ(1,*,ERR=1000)(NLOS(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(CLOS(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBW
        L=LCBW(LL)
        READ(1,*,ERR=1000)(NLOW(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(CLOW(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBE
        L=LCBE(LL)
        READ(1,*,ERR=1000)(NLOE(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(CLOE(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBN
        L=LCBN(LL)
        READ(1,*,ERR=1000)(NLON(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(CLON(LL,K,2),K=1,KC)
        READ(1,*,ERR=1000)(TEM(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(TEM1(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
      ENDDO
    ENDIF
    IF( ISTRAN(3) >= 1 .AND. ISCI(3) == 1 )THEN
      DO LL=1,NCBS
        L=LCBS(LL)
        READ(1,*,ERR=1000)(NLOS(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(CLOS(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBW
        L=LCBW(LL)
        READ(1,*,ERR=1000)(NLOW(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(CLOW(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBE
        L=LCBE(LL)
        READ(1,*,ERR=1000)(NLOE(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(CLOE(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBN
        L=LCBN(LL)
        READ(1,*,ERR=1000)(NLON(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(CLON(LL,K,3),K=1,KC)
        READ(1,*,ERR=1000)(DYE(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(DYE1(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
      ENDDO
    ENDIF
    IF( ISTRAN(4) >= 1 .AND. ISCI(4) == 1 )THEN
      DO LL=1,NCBS
        L=LCBS(LL)
        READ(1,*,ERR=1000)(NLOS(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(CLOS(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
      ENDDO
      DO LL=1,NCBW
        L=LCBW(LL)
        READ(1,*,ERR=1000)(NLOW(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(CLOW(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
      ENDDO
      DO LL=1,NCBE
        L=LCBE(LL)
        READ(1,*,ERR=1000)(NLOE(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(CLOE(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
        READ(1,*,ERR=1000)(UHDY(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDY1(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBN
        L=LCBN(LL)
        READ(1,*,ERR=1000)(NLON(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(CLON(LL,K,4),K=1,KC)
        READ(1,*,ERR=1000)(SED(LC,K,1),K=1,KC)
        READ(1,*,ERR=1000)(SED1(LC,K,1),K=1,KC)
        READ(1,*,ERR=1000)(VHDX(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDX1(LC,K),K=1,KC)
      ENDDO
    ENDIF
    IF( ISTRAN(5) >= 1 .AND. ISCI(5) == 1 )THEN
      DO LL=1,NCBS
        L=LCBS(LL)
        READ(1,*,ERR=1000)(NLOS(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(CLOS(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBW
        L=LCBW(LL)
        READ(1,*,ERR=1000)(NLOW(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(CLOW(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBE
        L=LCBE(LL)
        READ(1,*,ERR=1000)(NLOE(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(CLOE(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDYWQ(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(UHDYWQ(LC,K),K=1,KC)
      ENDDO
      DO LL=1,NCBN
        L=LCBN(LL)
        READ(1,*,ERR=1000)(NLON(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(CLON(LL,K,5),K=1,KC)
        READ(1,*,ERR=1000)(SFL(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(SFL2(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDXWQ(LC,K),K=1,KC)
        READ(1,*,ERR=1000)(VHDXWQ(LC,K),K=1,KC)
      ENDDO
    ENDIF
    CLOSE(1)
    DO K=1,KC
      SAL(1,K)=0.
      TEM(1,K)=0.
      DYE(1,K)=0.
      SED(1,K,1)=0.
      SFL(1,K)=0.
      VHDX(1,K)=0.
      UHDY(1,K)=0.
      SAL1(1,K)=0.
      TEM1(1,K)=0.
      DYE1(1,K)=0.
      SED1(1,K,1)=0.
      SFL2(1,K)=0.
      VHDX1(1,K)=0.
      UHDY1(1,K)=0.
      VHDXWQ(1,K)=0.
      UHDYWQ(1,K)=0.
      SAL(LC,K)=0.
      TEM(LC,K)=0.
      DYE(LC,K)=0.
      SED(LC,K,1)=0.
      SFL(LC,K)=0.
      VHDX(LC,K)=0.
      UHDY(LC,K)=0.
      SAL1(LC,K)=0.
      TEM1(LC,K)=0.
      DYE1(LC,K)=0.
      SED1(LC,K,1)=0.
      SFL2(LC,K)=0.
      VHDX1(LC,K)=0.
      UHDY1(LC,K)=0.
      VHDXWQ(LC,K)=0.
      UHDYWQ(LC,K)=0.
    ENDDO
    DO L=2,LA
      LS=LSC(L)
      H1U(L)=0.5*GI*(P1(L)+P1(LWC(L)))-0.5*(BELV(L)+BELV(LWC(L)))
      H1V(L)=0.5*GI*(P1(L)+P1(LS))-0.5*(BELV(L)+BELV(LS))
      H1P(L)=GI*P1(L)-BELV(L)
      HU(L)=0.5*GI*(P(L)+P(LWC(L)))-0.5*(BELV(L)+BELV(LWC(L)))
      HV(L)=0.5*GI*(P(L)+P(LS))-0.5*(BELV(L)+BELV(LS))
      HP(L)=GI*P(L)-BELV(L)
      HPI(L)=1./HP(L)
      HUI(L)=1./HU(L)
      HVI(L)=1./HV(L)
      H1UI(L)=1./H1U(L)
      H1VI(L)=1./H1V(L)
      H2WQ(L)=HP(L)
    ENDDO
    DO K=1,KC
      DO L=2,LA
        UHDYF1(L,K) = DYU(L)*H1U(L)*U1(L,K)
        VHDXF1(L,K) = DXV(L)*H1V(L)*V1(L,K)
        UHDY1(L,K)  = UHDYF1(L,K)*DZC(L,K)
        VHDX1(L,K)  = VHDXF1(L,K)*DZC(L,K)
        
        UHDYF(L,K)  = DYU(L)*HU(L)*U(L,K)
        VHDXF(L,K)  = DXV(L)*HV(L)*V(L,K)
        UHDY(L,K)   = UHDYF(L,K)*DZC(L,K)
        VHDX(L,K)   = VHDXF(L,K)*DZC(L,K)

        SAL(L,K)=MAX(SAL(L,K),0.)
        SAL1(L,K)=MAX(SAL1(L,K),0.)
      ENDDO
    ENDDO
  
    GOTO 1002

    ! **  WRITE READ ERRORS ON RESTART
      1000 WRITE(6,1001)
      1001 FORMAT('  READ ERROR ON FILE RESTART.INP ')
    STOP
      1002 CONTINUE
      907 FORMAT(12E12.4)
      908 FORMAT(12I10)

  END SUBROUTINE

  SUBROUTINE RSMRST
    ! ** READ ICS FROM RESTART FILE FROM INSMRST.
    ! ** THIS FILE IS ONLY USED WHEN iICI=2 (WQ3DSD.INP)
    ! ** CHECK FIRST TO SEE IF BINARY RESTART FILE EXISTS.  IF NOT, USE
    ! ** THE ASCII FILE INSTEAD.
    !
    LOGICAL FEXIST
    INTEGER :: NREST,M,L,NW
    REAL    :: RSTTIME

    INQUIRE(FILE='wqsdrst.bin', EXIST=FEXIST)
    IF( .NOT. FEXIST )THEN
      ! *** FORMATTED
      WRITE(*,'(A)')' WQ: READING WQSDRST.INP'
      OPEN(1,FILE='wqsdrst.inp',STATUS='UNKNOWN')
      READ(1,999)
      READ(1,999)
      DO M=2,LA
        READ(1,*) L,(SMPON(L,NW),NW=1,NSMG), &
            (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L), &
            SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L), &
            SMBST(L),SMT(L)
      ENDDO
      CLOSE(1)
    ELSE
      ! *** BINARY
      WRITE(*,'(A)')' WQ: READING WQSDRST.BIN'
      OPEN(1,FILE='wqsdrst.bin',STATUS='UNKNOWN',FORM='BINARY')
      READ(1) NREST, RSTTIME
      WRITE(0,911) NREST, RSTTIME
      911 FORMAT(' READING BINARY WQSDRST.BIN FILE ...    NN, TIME = ', I10, F11.5)
      DO M=2,LA
        READ(1) L
        READ(1) (SMPON(L,NW),NW=1,NSMG), &
            (SMPOP(L,NW),NW=1,NSMG),(SMPOC(L,NW),NW=1,NSMG),SM1NH4(L), &
            SM2NH4(L),SM2NO3(L),SM2PO4(L),SM2H2S(L),SMPSI(L),SM2SI(L), &
            SMBST(L),SMT(L)
      ENDDO
      CLOSE(1)
    ENDIF
     90 FORMAT(I5, 18E12.4)
    999 FORMAT(1X)

  END SUBROUTINE

  SUBROUTINE WQWCRST_IN
    ! ** 2013-06-09 DH CHUNG: CORRECT BINARY FORM
    ! ** INITIAL CONDITION FILE
    ! ** READ ICS FROM RESTART FILE FROM INWQRST(WQ3DWC.INP: IWQICI == 2 ).
    INTEGER :: L,K,NREST,LK,NWQV0,M,NW,LL,KK,VER,NCOL
    INTEGER(IK8) :: NNN
    REAL    :: RSTTIME
    CHARACTER(100) :: STR
    LOGICAL FEXIST

    LK=(LA-1)*KC

    ! CHECK FIRST TO SEE IF BINARY RESTART FILE EXISTS.  IF NOT, USE THE ASCII FILE INSTEAD.
    INQUIRE(FILE='wqwcrst.bin', EXIST=FEXIST)
    IF( .NOT. FEXIST )THEN
      ! *** FORMATTED
      WRITE(*,'(A)')' WQ: RESTART: WQWCRST.INP'
      OPEN(1,FILE='wqwcrst.inp',STATUS='UNKNOWN')
      
      ! *** SKIP THE FIRST LINE
      READ(1,999)

      READ(1,'(A)') STR
      NCOL = INDEX(STR,'RADIATION')
      IF( NCOL > 0 )THEN
        NCOL = INDEX(STR,'=') + 1
        !READ(STR,*) VER,NNN,RSTTIME
        READ(STR(NCOL:100),* ) WQI1,WQI2,WQI3

        ! *** SKIP THE HEADER LINE
        READ(1,999)
      ENDIF

      NWQV0=NWQV
      IF( IDNOTRVA > 0 ) NWQV0=NWQV0+1
      DO M=1,LK
        READ(1,* ) L,K,(WQV(L,K,NW),NW=1,NWQV0)
      ENDDO
      
      CLOSE(1)
    ELSE
      ! *** BINARY
      WRITE(*,'(A)')' WQ: RESTART: WQWCRST.BIN'
      OPEN(1,FILE='wqwcrst.bin',STATUS='UNKNOWN',FORM='BINARY')
      READ(1) NREST, RSTTIME
      WRITE(0,911) NREST, RSTTIME
911   FORMAT(' READING BINARY WQWCRST.BIN FILE ...    NN, TIME = ', I7, F11.5)
      
      NWQV0=NWQV
      IF( IDNOTRVA > 0 ) NWQV0=NWQV0+1
      DO LL=2,LA
        DO KK=1,KC
          READ(1) L, K  
          READ(1) (WQV(L, K, NW),NW=1,NWQV0)
        ENDDO
      ENDDO
      CLOSE(1)
    ENDIF
    
    ! INITIALIZE MACROALGAE SO BIOMASS ONLY EXISTS IN BOTTOM LAYER:
    IF( IDNOTRVA > 0 )THEN
      IF( KC > 1 )THEN
        DO K=1,KC
          DO L=2,LA
            IF( K == KSZ(L) ) CYCLE
            WQV(L,K,22)=0.
          ENDDO
        ENDDO
      ENDIF
    ENDIF
     90 FORMAT(2I5, 21E12.4)
    999 FORMAT(1X)
    
  END SUBROUTINE
 
  SUBROUTINE GENRSTINP(RESTARTF)
    ! ** EDIT EFDC.INP FOR RESTART
    ! ** COPY RESTART FILES INTO PROJECT FOLDER
    ! ** CHANGE THE FOLLOWING PARAMETERS:
    ! ** C2: ISRESTI =1 EFDC RUNS RESTART
    ! **     ICONTINUE=0 EFDC WRITES OUT FILES FROM BEGINING AS USUAL 
    ! **     ICONTINUE=1 EFDC READS OUT FILES & WRITES THE NEXT OUTPUTS AT THE TRUE
    ! **                LOCATION RIGHT AFTER THE TBEGINC
    ! ** C7  NTC CHANGES TO THE TOTAL NUMBER OF DAYS IN EFDC.INP
    ! ** C8: TBEGIN STILL KEEPS THE PREVIOUS VALUE   IN EFDC.INP
    ! ** RESET ISCI(1:8),ISCO(1:8) & SOME PARAMETERS
    CHARACTER(*) :: RESTARTF
    LOGICAL(4) :: RESLOG
    
    INTEGER :: SLEN,L
    
    CHARACTER*3  FLABEL*12
    CHARACTER(40) :: STR*200,RSTFILE1,RSTFILE2,RSTFILE3,RSTFILE4,RSTFILE5
  
    ! ** COPY FILES
    IF( ISGREGOR ==  1 )THEN
      RESTARTF = ADJUSTL(RESTARTF)
      SLEN = LEN_TRIM(RESTARTF)
      DO L=1,SLEN
        IF( RESTARTF(L:L) == '_') EXIT
      ENDDO
      FLABEL = RESTARTF(L:SLEN) 
    ELSE
      FLABEL = ''
    ENDIF
  
    RSTFILE1 = OUTDIR//'RESTART'//TRIM(FLABEL)//'.OUT'
    RSTFILE2 = OUTDIR//'RSTWD'//TRIM(FLABEL)//'.OUT'
    STR ='copy '//trim(RSTFILE1)//' restart.inp'

    WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: RESTART.INP'
    RESLOG = SYSTEMQQ(TRIM(STR))  

    IF( ISDRY > 0 )THEN
      WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: RSTWD.INP'
      STR ='copy '//trim(RSTFILE2)//' rstwd.inp'
      RESLOG = SYSTEMQQ(TRIM(STR))
    ENDIF
    
    IF( LSEDZLJ )THEN
      WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: SEDBED_HOT.SDF'
      RSTFILE5 = OUTDIR//'SEDBED_HOT'//TRIM(FLABEL)//'.OUT'
      STR ='copy '//trim(RSTFILE5)//' sedbed_hot.sdf'
      RESLOG = SYSTEMQQ(TRIM(STR))
    ENDIF
    
    ! ** WQ: THIS MUST COME AFTER READING EFDC.INP **     
    IF( ISTRAN(8) >= 1 )THEN    
      ! ** wq3dwc.inp
      IF( IWQRST == 1 )THEN
        WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: WQWCRST.INP'
        RSTFILE4 = OUTDIR//'WQWCRST'//TRIM(FLABEL)//'.OUT'   
        STR ='copy '//trim(RSTFILE4)//' wqwcrst.inp'
        RESLOG = SYSTEMQQ(trim(str))       
      ENDIF
      
      IF( IWQBEN == 1 .AND. ISMICI == 2 )THEN
        ! ** wq3dsd.inp
        WRITE(*,'(A)')'COPYING CONTINUATION FILE TO: WQSDRST.INP'
        RSTFILE3 = OUTDIR//'WQSDRST'//TRIM(FLABEL)//'.OUT'
        STR ='copy '//trim(RSTFILE3)//' wqsdrst.inp'
        RESLOG = SYSTEMQQ(TRIM(STR))       
      ENDIF
      
    ENDIF  
  
  CALL RESTIN1(1)  !READ NRESTART, TBEGINC, VER ONLY

  END SUBROUTINE
  
END MODULE
