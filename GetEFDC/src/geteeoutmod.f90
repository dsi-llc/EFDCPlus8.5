    MODULE GETEEOUTMOD
    ! ** AUTHOR: DH CHUNG

    USE GLOBALVARS
    USE EFDCPROMOD
    USE INFOMOD
    USE IFPORT
    USE TECMOD

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE GETEE_WS
    INTEGER(4)::VER,LINES,L,NT,N,I,J,L1,L2,HSIZE,BSIZE,IOS,ITMP
    REAL(4)::DTIME,TMP,DELT,DTIME1,DTIME2,MAXJTIM,MINJTIM
    REAL(8)::DTIME8
    CHARACTER(200)::SS*5

    WRITE(SS,'(I5)') 2*NLOC
    SFMT='(F12.4,'//SS//'F10.3)'
    OPEN(UOUTI(1),FILE=TSWSF,ACTION='WRITE')
    OPEN(UWSO,FILE=WSFILEO,ACTION='WRITE')
    DO N=1,NLOC
        WRITE(UOUTI(1),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
    ENDDO
    WRITE(UOUTI(1),'(A31,I3,A3)') '** TIME,(HP(L(N)),BELV(L),N=1, ',NLOC,')'
    IF (TECCHKL)  CALL TECPLOT(0.,1,1)
    IF (IGRIDV>0) THEN
        CALL TECPLOT(0.,11,1)
        CALL TECPLOT(1.,11,2)
        CLOSE(UTEC5)
    ENDIF

    OPEN(UWSI,FILE=WSFILEI,FORM='BINARY')
    READ(UWSI,IOSTAT=IOS) VER
    IF (IOS /= 0) GOTO 999

    WSVER = VER

    IF (VER>=8400) THEN  !2018-07-24, NTL:
        READ(UWSI) HSIZE,BSIZE,IC,JC,LINES
    ELSEIF (VER>7000) THEN  !2017-09-12, NTL: .AND.VER<=7300
        READ(UWSI) HSIZE,IC,JC,LINES
    ELSE  !2017-09-12, NTL: IF(VER<=7000) THEN
        READ(UWSI) IC,JC,LINES
    ENDIF
    LA = LINES+1
    NTM=0
    DO WHILE(1)
        IF (VER>=8300) THEN  !2017-09-12, NTL:
            READ (UWSI,END=100)N,DTIME8,DELT,TMP
            DTIME = REAL(DTIME8)
        ELSEIF(VER>=7300) THEN
            READ (UWSI,END=100)N,DTIME8,DELT
            DTIME = REAL(DTIME8)
        ELSE  !2017-09-12, NTL: IF (VER<7300) THEN
            READ (UWSI,END=100)N,DTIME,DELT
        ENDIF

        NTM=NTM+1
        DO L=2,LA
            READ(UWSI) TMP
        ENDDO

    ENDDO
100 REWIND(UWSI)

    JULTIME = MIN(JULTIME,DTIME)
    NTMWC = NTM+1
    DTIME2 = 0
    READ(UWSI) VER
    IF (VER>=8400) THEN  !2018-07-24, NTL:
        READ(UWSI) HSIZE,BSIZE,IC,JC,LINES
    ELSEIF (VER>7000) THEN
        READ(UWSI) HSIZE,LINES
    ELSE
        READ(UWSI) IC,JC,LINES
    ENDIF
    IF (.NOT. ALLOCATED(HPT)) ALLOCATE(HPT(LCM,NTM))
    IF (.NOT. ALLOCATED(JTIME)) ALLOCATE(JTIME(NTM))

    DO NT=1,NTM
        IF (VER>=8300) THEN   !2017-09-12, NTL:
            READ (UWSI,END=100)N,DTIME8,DELT,ITMP
            !PRINT* ,'I10,F15.5,F15.5,I10',N,DTIME8,DELT,ITMP
            DTIME = REAL(DTIME8)
        ELSEIF(VER>=7300) THEN
            READ (UWSI,END=100)N,DTIME8,DELT
            DTIME = REAL(DTIME8)
        ELSE  !2017-09-12, NTL: IF (VER<7300) THEN
            READ (UWSI,END=100)N,DTIME,DELT
        ENDIF

        DTIME1=DTIME2
        DTIME2=DTIME
        JTIME(NT) = DTIME
        DO L=2,LA
            READ(UWSI)HPT(L,NT)
        ENDDO

        !** TIME SERIES
        IF(VER>=7300) THEN
            WRITE(UOUTI(1),SFMT) DTIME8,(HPT(LIJ(ICEL(N),JCEL(N)),NT),BELV(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)
        ELSE
            WRITE(UOUTI(1),SFMT) DTIME,(HPT(LIJ(ICEL(N),JCEL(N)),NT),BELV(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)
        ENDIF
        !** TECPLOT & LAYER
        IF(ALLSNAP==0.AND.DTIME1<JULTIME.AND.JULTIME<=DTIME2) THEN
            IF (TECCHKL) CALL TECPLOT(DTIME,1,0)
            WRITE(UWSO,'(A8,I5)') '** LA = ',LA
            WRITE(UWSO,'(A)') '** WATER DEPTH [M]: (HP(L),L=2,LA)'
            IF (IGRIDV>0) THEN
                DO L=2,LA
                    IF (KSZ(L)==1) THEN
                        WRITE(UWSO,'(A27,I8)') '** CELL WITH MAX LAYER L = ',L
                    ENDIF
                ENDDO
            ENDIF
            WRITE(UWSO,'(A)') '** TIME      DELT[S]'
            WRITE(UWSO,'(F15.5,F10.5)') DTIME,DELT
            WRITE(UWSO,'(10F10.3)') (HPT(L,NT),L=2,LA)
        ENDIF

        IF (ALLSNAP==1) THEN
            IF(TECCHKL) CALL TECPLOT(DTIME,1,0)
            !OUTPUT FOR THE WHOLE DOMAIN AND SNAPSHOTS
            WRITE(UWSO,'(A)') '** WATER DEPTH [M]: (HP(L),L=2,LA)'
            IF (IGRIDV>0) THEN
                DO L=2,LA
                    IF (KSZ(L)==1) THEN
                        WRITE(UWSO,'(A27,I8)') '** CELL WITH MAX LAYER L = ',L
                    ENDIF
                ENDDO
            ENDIF
            WRITE(UWSO,'(A)') '** N:SNAPSHOT    TIME      DELT[S]'
            WRITE(UWSO,'(I10,F15.5,F10.5)') NT,DTIME,DELT
            WRITE(UWSO,'(10F10.3)') (HPT(L,NT),L=2,LA)
        ENDIF
    ENDDO
101 CLOSE(UWSI)
    CLOSE(UOUTI(1))
    CLOSE(UWSO)
    IF (TECCHKL) CLOSE(UTEC1)
    MAXJTIM = MAXVAL(JTIME)
    MINJTIM = MINVAL(JTIME)
    RETURN
999 STOP ' **** OPENING EE_WS.OUT ERROR'
    END SUBROUTINE

    SUBROUTINE GETEE_VE
    INTEGER(4)::VER,LINES,L,NT,K,NTM1,N,LA1,LN,I,J,LE,IOS,KLL,K1
    INTEGER(4)::HSIZE,BSIZE,IGRIDV,CELL3D
    REAL(4)   ::DTIME,TMP,DELT,UTMPS,VTMPS,DTIME1,DTIME2
    REAL(4)   ::ZK1,ZK2,VKT,VKB,ZABV
    REAL(8)   ::DTIME8
    CHARACTER(200)::SS*5

    OPEN(UVEI,FILE=VEFILEI,FORM='BINARY')
    WRITE(SS,'(I5)') 3*NLOC
    SFMT='(F12.4,'//SS//'F12.5)'

    IF (LAYK>=0) OPEN(UVEO,FILE=VEFILEO,ACTION='WRITE')      !LAYER
    OPEN(UOUTI(1),FILE=TSVEF,ACTION='WRITE')    !TSERIES
    WRITE(UOUTI(1),'(A12,I5)') '** LAYER K =',LAYK
    DO N=1,NLOC
        WRITE(UOUTI(1),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
    ENDDO
    ALLOCATE(UKZ(NLOC),VKZ(NLOC),WKZ(NLOC))
    IF (LAYK/=0) THEN
        !VEL AT LAYER K
        IF (TECCHKL.AND.LAYK>0) CALL TECPLOT(0.,2,1)
        WRITE(UOUTI(1),'(A)') '** TIME,(U(L(N),K),V(L(N),K),W(L(N),K),N=1, '//SSTR3//') [M/S]' !TSER
    ELSEIF(LAYK==0) THEN
        !DEPTH-AVERAGED VEL
        IF (TECCHKL) CALL TECPLOT(0.,21,1)
        WRITE(UOUTI(1),'(A)') '** TIME,(UA(L(N)),VA(L(N)),WA(L(N)),N=1, '//SSTR3//') [M/S]'    !TSER
    ENDIF
    IF (VPROFCHK.AND.LAYK>=0) THEN
        OPEN(UCVPRF(1),FILE=VELPROFILE,ACTION='WRITE')
        WRITE(UCVPRF(1),'(A)') '** VERTICAL PROFILE OF VELOCITY MAGNITUDE (M/S)'
        IF (VPROFILE==1) THEN
            WRITE(UCVPRF(1),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
        ELSEIF(VPROFILE==2) THEN
            WRITE(UCVPRF(1),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
        ENDIF
        DO N=1,NLOC
            WRITE(UCVPRF(1),'(A6,I5,A4,I5)') '** I =',ICEL(N),' J =',JCEL(N)
        ENDDO
    ENDIF

    READ(UVEI,IOSTAT=IOS) VER
    IF (IOS /= 0) GOTO 999

    IF (VER>=8400) THEN  !2018-07-24, NTL:
        READ(UVEI) HSIZE,BSIZE,IGRIDV,CELL3D,IC,JC,KC,LINES
    ELSEIF (VER>7000) THEN
        READ(UVEI) HSIZE,IC,JC,KC,LINES
    ELSE
        READ(UVEI) IC,JC,KC,LINES
    ENDIF

    LA = LINES+1
    READ(UVEI) RSSBCE,RSSBCW,RSSBCS,RSSBCN
    NTM1=0
    DO WHILE(1)
        IF (VER<7300) THEN
            READ (UVEI,END=100)N,DTIME,DELT
        ELSEIF(VER>=7300) THEN
            READ (UVEI,END=100)N,DTIME8,DELT
            DTIME = DTIME8
        ENDIF
        NTM1=NTM1+1
        DO L=2,LA
            K1 = KSZ(L)
            READ(UVEI) (TMP,TMP,TMP,K=K1,KC)
        ENDDO
    ENDDO
100 REWIND(UVEI)

    READ(UVEI) VER
    IF (VER>=8400) THEN  !2018-07-24, NTL:
        READ(UVEI) HSIZE,BSIZE,IGRIDV,CELL3D,IC,JC,KC,LINES
    ELSEIF (VER>7000) THEN
        READ(UVEI) HSIZE,IC,JC,KC,LINES
    ELSE
        READ(UVEI) IC,JC,KC,LINES
    ENDIF
    READ(UVEI) RSSBCE,RSSBCW,RSSBCS,RSSBCN

    IF (NTM>0.AND.NTM1.NE.NTM) THEN
        !SNAPSHOT NUMBERS OF WS AND VE ARE DIFFERENT!
        STOP ' **** SNAPSHOT NUMBERS OF WS AND VE ARE DIFFERENT!'
    ENDIF
    U  = 0
    V  = 0
    W  = 0
    UK = 0
    VK = 0
    DTIME2=0
    DO NT=1,NTM
        !READ (UVEI,END=101) N,DTIME,DELT
        IF (VER<7300) THEN
            READ (UVEI,END=100)N,DTIME,DELT
        ELSEIF(VER>=7300) THEN
            READ (UVEI,END=100)N,DTIME8,DELT
            DTIME = DTIME8
        ENDIF
        IF (VER>=8400) THEN  !2018-07-24, NTL:
            READ(UVEI) ((U(L,K),K=KSZ(L),KC),L=2,LA)
            READ(UVEI) ((V(L,K),K=KSZ(L),KC),L=2,LA)
            READ(UVEI) ((W(L,K),K=KSZ(L),KC),L=2,LA)
        ELSE
            DO L=2,LA
                K1 = KSZ(L)
                READ(UVEI) (U(L,K),V(L,K),W(L,K),K=K1,KC)
            ENDDO
        ENDIF
        DTIME1=DTIME2
        DTIME2=DTIME
        DO L=2,LA
            DO K=1,KC
                IF (ROT==1) THEN
                    !CONVERTING CELL CENTER VELOCITY COMPONENTS TO TRUE EAST AND NORTH
                    LN=LNC(L)
                    LE=LEC(L)
                    UTMPS=0.5*(RSSBCE(L)*U(LE,K)+RSSBCW(L)*U(L,K))
                    VTMPS=0.5*(RSSBCN(L)*V(LN ,K)+RSSBCS(L)*V(L,K))
                    UK(L,K)=CUE(L)*UTMPS+CVE(L)*VTMPS
                    VK(L,K)=CUN(L)*UTMPS+CVN(L)*VTMPS
                ELSE
                    UK(L,K) = U(L,K)
                    VK(L,K) = V(L,K)
                ENDIF
            ENDDO
            UA(L) =SUM(UK(L,1:KC)*HPT(L,NT)*DZCK(1:KC))/HPT(L,NT)
            VA(L) =SUM(VK(L,1:KC)*HPT(L,NT)*DZCK(1:KC))/HPT(L,NT)
            WA(L) =SUM(W(L,1:KC)*HPT(L,NT)*DZCK(1:KC))/HPT(L,NT)
        ENDDO

        IF (LAYK>0) THEN
            IF(ALLSNAP==0.AND.DTIME1<JULTIME.AND.JULTIME<=DTIME2) THEN
                WRITE(UVEO,'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                WRITE(UVEO,'(A)') '** VELOCITY [M/S]:(U(L,K),V(L,K),W(L,K),L=2,LA)'
                WRITE(UVEO,'(A)') '** TIME      DELT [S]'
                WRITE (UVEO,'(F15.5,F10.5)') DTIME,DELT
                WRITE (UVEO,'(10F10.5)') (UK(L,LAYK),VK(L,LAYK),W(L,LAYK),L=2,LA)
                IF (TECCHKL)  CALL TECPLOT(DTIME,2,0)
                IF (VPROFCHK) CALL VPROFOUT('VEL',DTIME,NT)
            ENDIF
            IF(ALLSNAP==1) THEN
                WRITE(UVEO,'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                WRITE(UVEO,'(A)') '** VELOCITY [M/S]:(U(L,K),V(L,K),W(L,K),L=2,LA)'
                WRITE(UVEO,'(A)') '** TIME      DELT [S]'
                WRITE (UVEO,'(F15.5,F10.5)') DTIME,DELT
                WRITE (UVEO,'(10F10.5)') (UK(L,LAYK),VK(L,LAYK),W(L,LAYK),L=2,LA)
                IF (TECCHKL) CALL TECPLOT(DTIME,2,0)
                IF (VPROFCHK) CALL VPROFOUT('VEL',DTIME,NT)
            ENDIF
            !** TSERIES
            WRITE(UOUTI(1),SFMT) DTIME,(UK(LIJ(ICEL(N),JCEL(N)),LAYK),VK(LIJ(ICEL(N),JCEL(N)),LAYK),W(LIJ(ICEL(N),JCEL(N)),LAYK),N=1,NLOC)

        ELSEIF(LAYK==0) THEN
            !DEPTH-AVERAGE
            IF(ALLSNAP==0.AND.DTIME1<JULTIME.AND.JULTIME<=DTIME2) THEN
                WRITE(UVEO,'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                WRITE(UVEO,'(A)') '** VELOCITY [M/S]: DO K=1,KC (UA(L),VA(L),WA(L),L=2,LA) ENDDO'
                WRITE(UVEO,'(A)') '**   TIME      DELT[S]'
                WRITE (UVEO,'(F15.5,F10.5)') DTIME,DELT
                WRITE (UVEO,'(10F10.5)') (UA(L),VA(L),WA(L),L=2,LA)
                IF (TECCHKL) CALL TECPLOT(DTIME,21,0)
                IF(VPROFCHK) CALL VPROFOUT('VEL',DTIME,NT)
            ENDIF

            IF(ALLSNAP==1) THEN
                WRITE(UVEO,'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                WRITE(UVEO,'(A)') '** VELOCITY [M/S]: DO K=1,KC (UA(L),VA(L),WA(L),L=2,LA) ENDDO'
                WRITE(UVEO,'(A)') '**   TIME      DELT[S]'
                WRITE (UVEO,'(F15.5,F10.5)') DTIME,DELT
                WRITE (UVEO,'(10F10.5)') (UA(L),VA(L),WA(L),L=2,LA)
                IF (TECCHKL) CALL TECPLOT(DTIME,21,0)
                IF (VPROFCHK) CALL VPROFOUT('VEL',DTIME,NT)
            ENDIF

            !** TSERIES
            WRITE(UOUTI(1),SFMT) DTIME,(UA(LIJ(ICEL(N),JCEL(N))),VA(LIJ(ICEL(N),JCEL(N))),WA(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)

        ELSEIF (LAYK==-2) THEN
            DO N=1,NLOC
                IF (ZOPT==1) THEN
                    ZABV = HPT(LIJ(ICEL(N),JCEL(N)),NT) - ZINT(N)  !CONVERT DEPTH UNDER WS TO HEIGHT ABOVE BED
                ELSEIF(ZOPT==2) THEN
                    ZABV = ZINT(N)                                !THE HEIGHT ABOVE BED
                ENDIF
                ZSIG(N) = ZABV/HPT(LIJ(ICEL(N),JCEL(N)),NT)
                DO K=1,KC
                    IF (ZSIG(N) <= SUM(DZCK(1:K))) THEN
                        EXIT
                    ENDIF
                ENDDO
                K = MIN(K,KC)
                ZK2 = SUM(DZCK(1:K))
                ZK1 = ZK2 - DZCK(K)

                VKT = UK(LIJ(ICEL(N),JCEL(N)),K)
                VKB = UK(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1))
                UKZ(N) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB

                VKT = VK(LIJ(ICEL(N),JCEL(N)),K)
                VKB = VK(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1))
                VKZ(N) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB

                VKT = W(LIJ(ICEL(N),JCEL(N)),K)
                VKB = W(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1))
                WKZ(N) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB

            ENDDO
            !** TSERIES
            WRITE(UOUTI(1),SFMT) DTIME,(UKZ(N),VKZ(N),WKZ(N),N=1,NLOC)

        ENDIF
    ENDDO
101 CLOSE(UVEI)
    CLOSE(UOUTI(1))
    IF (LAYK>=0)  CLOSE(UVEO)
    IF (TECCHKL)  CLOSE(UTEC2)
    IF (VPROFCHK) CLOSE(UCVPRF(1))
    RETURN
999 STOP ' **** OPENING EE_VEL.OUT ERROR'
    END SUBROUTINE

    SUBROUTINE GETEE_WC
    INTEGER(4)::VER,LINES,L,NT,K,NTM1,N,LA1,NSXD,NS,NBEDSTEPS,NTS,NLYR
    INTEGER(4)::NW,MW,NACTIVE,N1,NX,JSEXPLORER,I,J,IOS
    INTEGER(4)::HSIZE,BSIZE,IGRIDV,CELL3D,IEVAP,ISGWIE,ISICE
    INTEGER(4)::IWQ(40)
    REAL(4)   ::TMP,WQ,EETIME,SHEAR,TMPVAL,DTIME1,DTIME2
    REAL(4)   ::ZK1,ZK2,VKT,VKB,ZABV
    REAL(8)   :: DTIME8
    CHARACTER*8 ARRAYNAME
    CHARACTER(200)::SS*5
    LOGICAL(4)::ASW

    IF (TECCHKL.AND.LAYK>=0)  CALL TECPLOT(0.,3,1)            !ALL CONSITUENTS
    IF(ISTRAN(1)==1) THEN
        SALLAYF  = '#output\RESULT\SAL_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
        SALTSF   = '#output\RESULT\SAL_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
        OPEN(UOUTI(11),FILE=SALTSF,ACTION='WRITE')   !TSER
        DO N=1,NLOC
            WRITE(UOUTI(11),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
        ENDDO
        WRITE(UOUTI(11),'(A)') '** EETIME,(SAL(L(N),'//SSTR2//',NT),N=1,'//SSTR3//') [PPT]'
        IF (LAYK>=0) OPEN(UOUTI(1),FILE=SALLAYF,ACTION='WRITE')    !LAYER
        IF (VPROFCHK.AND.LAYK>=0) THEN
            OPEN(UCVPRF(1),FILE=SALPROFILE,ACTION='WRITE')
            WRITE(UCVPRF(1),'(A)') '** VERTICAL PROFILE OF SALINITY (PPT)'
            IF (VPROFILE==1) THEN
                WRITE(UCVPRF(1),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
            ELSEIF(VPROFILE==2) THEN
                WRITE(UCVPRF(1),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
            ENDIF
            DO N=1,NLOC
                WRITE(UCVPRF(1),'(A6,I5,A4,I5)') '** I =',ICEL(N),' J =',JCEL(N)
            ENDDO
        ENDIF
    ENDIF

    IF(ISTRAN(2)==1) THEN
        TEMLAYF  = '#output\RESULT\TEM_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
        TEMTSF   = '#output\RESULT\TEM_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
        OPEN(UOUTI(12),FILE=TEMTSF,ACTION='WRITE')                     !TSER
        DO N=1,NLOC
            WRITE(UOUTI(12),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
        ENDDO
        WRITE(UOUTI(12),'(A)') '** EETIME,(TEM(L(N),'//SSTR2//'),N=1,'//SSTR3//') [C]'
        IF (LAYK>=0) OPEN(UOUTI(2),FILE=TEMLAYF,ACTION='WRITE')                    !LAYER
        IF (VPROFCHK.AND.LAYK>=0) THEN
            OPEN(UCVPRF(2),FILE=TEMPROFILE,ACTION='WRITE')
            WRITE(UCVPRF(2),'(A)') '** VERTICAL PROFILE OF TEMPERATURE (C)'
            IF (VPROFILE==1) THEN
                WRITE(UCVPRF(2),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
            ELSEIF(VPROFILE==2) THEN
                WRITE(UCVPRF(2),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
            ENDIF
            DO N=1,NLOC
                WRITE(UCVPRF(2),'(A6,I5,A4,I5)') '** I =',ICEL(N),' J =',JCEL(N)
            ENDDO
        ENDIF
    ENDIF

    IF(ISTRAN(3)==1) THEN
        DYELAYF  = '#output\RESULT\DYE_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
        DYETSF   = '#output\RESULT\DYE_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
        OPEN(UOUTI(13),FILE=DYETSF,ACTION='WRITE')
        DO N=1,NLOC
            WRITE(UOUTI(13),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
        ENDDO
        WRITE(UOUTI(13),'(A)') '** EETIME,(DYE(L(N),'//SSTR2//'),N=1,'//SSTR3//') [MG/L]'
        IF (LAYK>=0) OPEN(UOUTI(3),FILE=DYELAYF,ACTION='WRITE')
        IF (VPROFCHK.AND.LAYK>=0) THEN
            OPEN(UCVPRF(3),FILE=DYEPROFILE,ACTION='WRITE')
            WRITE(UCVPRF(3),'(A)') '** VERTICAL PROFILE OF DYE (MG/L)'
            IF (VPROFILE==1) THEN
                WRITE(UCVPRF(3),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
            ELSEIF(VPROFILE==2) THEN
                WRITE(UCVPRF(3),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
            ENDIF
            DO N=1,NLOC
                WRITE(UCVPRF(3),'(A6,I5,A4,I5)') '** I =',ICEL(N),' J =',JCEL(N)
            ENDDO
        ENDIF
    ENDIF

    IF(ISTRAN(4)==1) THEN
        SFLLAYF  = '#output\RESULT\SFL_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
        SFLTSF   = '#output\RESULT\SFL_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
        OPEN(UOUTI(14),FILE=SFLTSF,ACTION='WRITE')
        WRITE(UOUTI(14),'(A)') '** EETIME,(SFL(L(N),'//SSTR2//'),N=1,'//SSTR3//')'
        IF (LAYK>=0) OPEN(UOUTI(4),FILE=SFLLAYF,ACTION='WRITE')
    ENDIF

    IF(ISTRAN(5)==1) THEN
        TOXLAYF  = '#output\RESULT\TOX_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
        TOXTSF   = '#output\RESULT\TOX_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
        OPEN(UOUTI(15),FILE=TOXTSF,ACTION='WRITE')
        DO N=1,NLOC
            WRITE(UOUTI(15),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
        ENDDO
        WRITE(UOUTI(15),'(A)') '** EETIME,((TOX(L,'//SSTR2//',NT),L=2,LA),NT=1,NTOX) [uG/L]'
        IF (LAYK>=0) OPEN(UOUTI(5),FILE=TOXLAYF,ACTION='WRITE')

        IF (VPROFCHK.AND.LAYK>=0) THEN
            OPEN(UCVPRF(5),FILE=TOXPROFILE,ACTION='WRITE')
            WRITE(UCVPRF(5),'(A)') '** VERTICAL PROFILE OF TOXIC (MG/L)'
            IF (VPROFILE==1) THEN
                WRITE(UCVPRF(5),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
            ELSEIF(VPROFILE==2) THEN
                WRITE(UCVPRF(5),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
            ENDIF
            DO N=1,NLOC
                WRITE(UCVPRF(5),'(A6,I5,A4,I5)') '** I =',ICEL(N),' J =',JCEL(N)
            ENDDO
        ENDIF
    ENDIF

    IF(ISTRAN(6)==1) THEN
        SEDLAYF  = '#output\RESULT\SED_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
        SEDTSF  = '#output\RESULT\SED_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
        OPEN(UOUTI(16),FILE=SEDTSF,ACTION='WRITE')
        DO N=1,NLOC
            WRITE(UOUTI(16),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
        ENDDO
        WRITE(UOUTI(16),'(A)') '** EETIME,((SED(L(N),'//SSTR2//',NS),N=1,'//SSTR3//'),NS=1,NSED) [MG/L]'

        IF (LAYK>=0) OPEN(UOUTI(6), FILE=SEDLAYF,ACTION='WRITE')

        IF (VPROFCHK.AND.LAYK>=0) THEN
            OPEN(UCVPRF(6),FILE=SEDPROFILE,ACTION='WRITE')
            WRITE(UCVPRF(6),'(A)') '** VERTICAL PROFILE OF SEDIMENT (MG/L)'
            IF (VPROFILE==1) THEN
                WRITE(UCVPRF(6),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
            ELSEIF(VPROFILE==2) THEN
                WRITE(UCVPRF(6),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
            ENDIF
            DO N=1,NLOC
                WRITE(UCVPRF(6),'(A6,I5,A4,I5)') '** I =',ICEL(N),' J =',JCEL(N)
            ENDDO
        ENDIF
    ENDIF

    IF(ISTRAN(7)==1) THEN
        SNDLAYF  = '#output\RESULT\SND_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
        SNDTSF  = '#output\RESULT\SND_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
        OPEN(UOUTI(17),FILE=SNDTSF,ACTION='WRITE')
        DO N=1,NLOC
            WRITE(UOUTI(17),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
        ENDDO
        WRITE(UOUTI(17),'(A)') '** EETIME,(SND(L(N),'//SSTR2//',NX),N=1,'//SSTR3//'),NX=1,NSND) [MG/L]'
        IF (LAYK>=0) OPEN(UOUTI(7),FILE=SNDLAYF,ACTION='WRITE')

        IF (VPROFCHK.AND.LAYK>=0) THEN
            OPEN(UCVPRF(7),FILE=SNDPROFILE,ACTION='WRITE')
            WRITE(UCVPRF(7),'(A)') '** VERTICAL PROFILE OF NON-COHESIVE SED (MG/L)'
            IF (VPROFILE==1) THEN
                WRITE(UCVPRF(7),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
            ELSEIF(VPROFILE==2) THEN
                WRITE(UCVPRF(7),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
            ENDIF
            DO N=1,NLOC
                WRITE(UCVPRF(7),'(A6,I5,A4,I5)') '** I =',ICEL(N),' J =',JCEL(N)
            ENDDO
        ENDIF
    ENDIF

    IF (KB==1) CALL BEDOUT(1)

    IF (ISICE >=3) THEN
        ALLOCATE(RICECOVL(LCM),RICETHKL(LCM))
        OPEN(UICE,FILE=TSICF,ACTION='WRITE')
        DO N=1,NLOC
            WRITE(UICE,'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
        ENDDO
        WRITE(UICE,'(A)') '** TIME,(RICETHKL(L(N)),N=1, '//TRIM(SS)//') [M]'
    ENDIF

    !READ DATA FROM EE_WC.OUT
    OPEN(UWCI,FILE=WCFILEI,FORM='BINARY')
    READ(UWCI,IOSTAT=IOS)VER
    IF (IOS /= 0) GOTO 999

    IF (VER>=8400) THEN  !2018-07-24, NTL:
        READ(UWCI) HSIZE,BSIZE,IGRIDV,CELL3D,IC,JC,KC,KB,LA1
        READ(UWCI) (ISTRAN(I),I=1,7)
        READ(UWCI) NSED,NSND,NTOX
        READ(UWCI) ISWAVE,ISBEDSTR,LSEDZLJ,ISBDLDBC,TBEDIT
        READ(UWCI) IEVAP,ISGWIE,ISICE
    ELSEIF (VER>=109) THEN
        READ(UWCI) HSIZE
        READ(UWCI) KC,KB,LA1
        READ(UWCI) (ISTRAN(I),I=1,7)
        READ(UWCI) NSED,NSND,NTOX
        READ(UWCI) ISWAVE,ISBEDSTR,LSEDZLJ,ISBDLDBC,TBEDIT
    ELSEIF (VER<109) THEN
        READ(UWCI)ISTRAN(1),ISTRAN(2),ISTRAN(3),ISTRAN(4)
        READ(UWCI)ISTRAN(5),ISTRAN(6),ISTRAN(7)
        READ(UWCI)NSED,NSND,KB,KC,NTOX
    ENDIF
    NSXD=NSED+NSND

    WRITE(SS,'(I5)') NLOC
    SFMT1='(F12.4,'//SS//'F12.5)'
    WRITE(SS,'(I5)') NLOC*NTOX
    SFMT2='(F12.4,'//SS//'F10.3)'
    WRITE(SS,'(I5)') NLOC*NSED
    SFMT3='(F12.4,'//SS//'F10.3)'
    WRITE(SS,'(I5)') NLOC*NSND
    SFMT4='(F12.4,'//SS//'F10.3)'

    ALLOCATE(SAL(LA,KC))
    ALLOCATE(SJDAY(NTMWC))
    ALLOCATE(TEM(LA,KC))
    ALLOCATE(DYE(LA,KC))
    ALLOCATE(SFL(LA,KC))
    ALLOCATE(TOX(LCM,KCM,NTXM))
    ALLOCATE(SED(LCM,KCM,NSCM))
    ALLOCATE(SND(LCM,KCM,NSNM))
    ALLOCATE(TEMB(LA))
    ALLOCATE(QQWV1(LA))
    ALLOCATE(SEDDIA(NSXD))
    ALLOCATE(KBT(LA))
    ALLOCATE(TAUBSED(LA))
    ALLOCATE(TAUBSND(LA))
    ALLOCATE(TAUB(LA))
    ALLOCATE(BEDSHR(LA))
    ALLOCATE(WVWHA(LA))
    ALLOCATE(WVFRQL(LA))
    ALLOCATE(WACCWE(LA))
    ALLOCATE(WVDISP(LA,KC))
    ALLOCATE(WVHUU(LA,KC))
    ALLOCATE(WVHVV(LA,KC))
    ALLOCATE(WVHUV(LA,KC))
    ALLOCATE(TOXB(LA,KBM,NTOX))
    ALLOCATE(HBED(LA,KBM))
    ALLOCATE(BDENBED(LA,KBM))
    ALLOCATE(PORBED(LA,KBM))
    ALLOCATE(SEDB(LA,KBM,NSED))
    ALLOCATE(SNDB(LA,KBM,NSND))
    ALLOCATE(CQBEDLOADX(LA,NSND))
    ALLOCATE(CQBEDLOADY(LA,NSND))
    ALLOCATE(VFRBED(LA,KBM,NSTM))
    ALLOCATE(SALA(LA))
    ALLOCATE(TEMA(LA))
    ALLOCATE(DYEA(LA))
    ALLOCATE(SFLA(LA))
    ALLOCATE(TOXA(LCM,NTXM))
    ALLOCATE(SEDA(LCM,NSCM))
    ALLOCATE(SNDA(LCM,NSNM))
    ALLOCATE(SALZ(NLOC),TEMZ(NLOC),DYEZ(NLOC))
    ALLOCATE(TOXZ(NLOC,NTOX),SEDZ(NLOC,NSED),SNDZ(NLOC,NSND))

    SAL = BLK
    TEM = BLK
    DYE = BLK
    TOX = BLK
    SED = BLK
    SND = BLK

    SALA = BLK
    TEMA = BLK
    DYEA = BLK
    TOXA = BLK
    SEDA = BLK
    SNDA = BLK
    BEDSHR = BLK

    DO NS=1,NSXD
        READ(UWCI) SEDDIA(NS)
    ENDDO

    NTS=0
    DTIME2=0
    DO WHILE(1)
        NTS = NTS+1
        IF(VER<=200) THEN
            READ(UWCI,END=100) EETIME,NACTIVE
        ELSEIF(VER > 200 .AND. VER <7300) THEN
            READ(UWCI,END=100) EETIME
        ELSEIF(VER >=7300) THEN
            READ(UWCI,END=100) DTIME8
            EETIME = DTIME8
        ENDIF
        SJDAY(NTS) = EETIME
        DTIME1 = DTIME2
        DTIME2 = EETIME
        ! *** READ THE TOP LAYER INDEX
        IF(ISTRAN(6) > 0 .OR. ISTRAN(7) > 0)THEN
            DO L=2,LA
                READ(UWCI) KBT(L)
            ENDDO
        ENDIF

        ! *** READ THE WATER COLUMN AND TOP LAYER OF SEDIMENT DATA, IF NEEDED
        IF (VER>=8400) THEN  !2018-07-24, NTL:
            IF(ISTRAN(6) > 0 .OR. ISTRAN(7) > 0)THEN
                IF(LSEDZLJ)THEN
                    READ(UWCI)  (BEDSHR(L),L=2,LA)   !=TAU(L) * 0.1 / WATERDENS / 1000.
                ELSEIF(ISBEDSTR >= 1)THEN
                    READ(UWCI) (BEDSHR(L),L=2,LA)
                    IF(ISBEDSTR == 1)THEN
                        READ(UWCI) (BEDSHR(L),L=2,LA)
                    ENDIF
                ELSE
                    ! *** TOTAL BED SHEAR STRESS
                    READ(UWCI) (BEDSHR(L),L=2,LA)
                ENDIF
            ELSE
                ! *** TOTAL BED SHEAR STRESS
                READ(UWCI) (BEDSHR(L),L=2,LA)
            ENDIF

            IF(ISWAVE >= 1)THEN
                ! *** Shear due to Current Only
                READ(UWCI) (QQWV1(L),L=2,LA)  ! *** Bed Shear due to Waves Only
                READ(UWCI) (SHEAR,L=2,LA)     ! *** Bed Shear due to Current Only
                IF(ISWAVE >= 3)THEN
                    READ(UWCI) (WVWHA(L),L=2,LA)
                    READ(UWCI) (WVFRQL(L),L=2,LA)
                    READ(UWCI) (WACCWE(L),L=2,LA)
                    IF(ISWAVE == 4)THEN
                        ! ***      DISSIPATION  SXX         SYY         SXY   (M3/S2)
                        READ(UWCI) (WVDISP(L,KC),L=2,LA)
                        READ(UWCI) (WVHUU(L,KC),L=2,LA)
                        READ(UWCI) (WVHVV(L,KC),L=2,LA)
                        READ(UWCI) (WVHUV(L,KC),L=2,LA)
                    ENDIF
                ENDIF
            ENDIF

            IF(ISTRAN(1) > 0) THEN
                READ(UWCI) ((SAL(L,K),K=KSZ(L),KC),L=2,LA)
                DO L=2,LA
                    NLYR = KC - KSZ(L) + 1
                    SALA(L) = SUM(SAL(L,KSZ(L):KC))/NLYR
                ENDDO
            ENDIF

            IF(ISTRAN(2) > 0)THEN
                READ(UWCI) ((TEM(L,K),K=KSZ(L),KC),L=2,LA)
                IF(TBEDIT > 0) READ(UWCI) (TEMB(L),L=2,LA)
                DO L=2,LA
                    NLYR = KC - KSZ(L) + 1
                    TEMA(L) = SUM(TEM(L,KSZ(L):KC))/NLYR
                ENDDO
                !** VER 7300
                IF (IEVAP > 1 )THEN
                    READ(UWCI) (TMPVAL,L=2,LA)  !REAL(EVAPT(L),4)
                    READ(UWCI) (TMPVAL,L=2,LA)  !REAL(RAINT(L),4)
                ENDIF

            ENDIF

            IF(ISTRAN(3) > 0) THEN
                READ(UWCI) ((DYE(L,K),K=KSZ(L),KC),L=2,LA)
                DO L=2,LA
                    NLYR = KC - KSZ(L) + 1
                    DYEA(L) = SUM(DYE(L,KSZ(L):KC))/NLYR
                ENDDO
            ENDIF

            IF(ISTRAN(4) > 0) THEN
                READ(UWCI) ((SFL(L,K),K=KSZ(L),KC),L=2,LA)
                DO L=2,LA
                    NLYR = KC - KSZ(L) + 1
                    SFLA(L) = SUM(SFL(L,KSZ(L):KC))/NLYR
                ENDDO
            ENDIF

            IF(ISTRAN(5) > 0)THEN
                READ(UWCI) ((TOXB(L,KBT(L),NT),L=2,LA),NT=1,NTOX)
                READ(UWCI) (((TOX(L,K,NT),K=KSZ(L),KC),L=2,LA),NT=1,NTOX)
                DO L=2,LA
                    NLYR = KC - KSZ(L) + 1
                    DO NT=1,NTOX
                        TOXA(L,NT) = SUM(TOX(L,KSZ(L):KC,NT))/NLYR
                    ENDDO
                ENDDO
            ENDIF

            IF(ISTRAN(6) > 0 .OR. ISTRAN(7) > 0)THEN
                READ(UWCI) (BELV(L),L=2,LA)
                READ(UWCI) (HBED(L,KBT(L)),L=2,LA)
                READ(UWCI) (BDENBED(L,KBT(L)),L=2,LA)
                READ(UWCI) (PORBED(L,KBT(L)),L=2,LA)
                IF(ISTRAN(6) > 0) READ(UWCI) ((SEDB(L,KBT(L),NS),L=2,LA),NS=1,NSED)
                IF(ISTRAN(7) > 0) READ(UWCI) ((SNDB(L,KBT(L),NS),L=2,LA),NS=1,NSND)
                READ(UWCI) ((VFRBED(L,KBT(L),NS),L=2,LA),NS=1,NSED+NSND)
                IF(ISTRAN(6) > 0)THEN
                    READ(UWCI) (((SED(L,K,NS),K=KSZ(L),KC),L=2,LA),NS=1,NSED)
                    DO L=2,LA
                        NLYR = KC - KSZ(L) + 1
                        DO NS=1,NSED
                            SEDA(L,NS) = SUM(SED(L,KSZ(L):KC,NS))/NLYR
                        ENDDO
                    ENDDO
                ENDIF
                IF(ISTRAN(7) > 0)THEN
                    READ(UWCI) (((SND(L,K,NX),K=KSZ(L),KC),L=2,LA),NX=1,NSND)
                    DO L=2,LA
                        NLYR = KC - KSZ(L) + 1
                        DO NX=1,NSND
                            SNDA(L,NX) = SUM(SND(L,KSZ(L):KC,NX))/NLYR
                        ENDDO
                    ENDDO
                    IF(ISBDLDBC > 0)THEN
                        READ(UWCI) ((CQBEDLOADX(L,NX),L=2,LA),NX=1,NSND)
                        READ(UWCI) ((CQBEDLOADY(L,NX),L=2,LA),NX=1,NSND)
                    ENDIF
                ENDIF
            ENDIF

            IF (ISGWIE > 0 ) THEN
                READ(UWCI) (EVAPSW(L),L=2,LA)
                READ(UWCI) (EVAPGW(L),L=2,LA)
                READ(UWCI) (RIFTR(L),L=2,LA)
                READ(UWCI) (AGWELV(L),L=2,LA)
            ENDIF

            IF(ISTRAN(2) > 0 .AND. ISICE >=3)THEN
                READ(UWCI) (RICETHKL(L),L=2,LA)
                READ(UWCI) (RICECOVL(L),L=2,LA)
            ENDIF
        ELSE
            DO L=2,LA
                N1=KBT(L)
                NLYR = KC - KSZ(L) + 1
                IF(ISTRAN(6) > 0 .OR. ISTRAN(7) > 0)THEN
                    IF(LSEDZLJ)THEN
                        READ(UWCI)  BEDSHR(L)   !=TAU(L) * 0.1 / WATERDENS / 1000.
                    ELSEIF(ISBEDSTR >= 1)THEN
                        READ(UWCI) BEDSHR(L)
                        IF(ISBEDSTR == 1)THEN
                            READ(UWCI) BEDSHR(L)
                        ENDIF
                    ELSE
                        ! *** TOTAL BED SHEAR STRESS
                        READ(UWCI) BEDSHR(L)
                    ENDIF
                ELSE
                    ! *** TOTAL BED SHEAR STRESS
                    READ(UWCI) BEDSHR(L)
                ENDIF

                IF(ISWAVE >= 1)THEN
                    ! *** Shear due to Current Only
                    READ(UWCI)QQWV1(L)  ! *** Bed Shear due to Waves Only
                    READ(UWCI)SHEAR     ! *** Bed Shear due to Current Only
                    IF(ISWAVE >= 3)THEN
                        READ(UWCI) WVWHA(L),WVFRQL(L),WACCWE(L)
                        IF(ISWAVE == 4)THEN
                            ! ***      DISSIPATION  SXX         SYY         SXY   (M3/S2)
                            READ(UWCI) WVDISP(L,KC),WVHUU(L,KC),WVHVV(L,KC),WVHUV(L,KC)
                        ENDIF
                    ENDIF
                ENDIF

                IF(ISTRAN(1) == 1) THEN
                    READ(UWCI)(SAL(L,K),K=KSZ(L),KC)
                    SALA(L) = SUM(SAL(L,KSZ(L):KC))/NLYR
                ENDIF

                IF(ISTRAN(2) == 1)THEN
                    READ(UWCI)(TEM(L,K),K=KSZ(L),KC)
                    IF(TBEDIT > 0) READ(UWCI) TEMB(L)
                    TEMA(L) = SUM(TEM(L,KSZ(L):KC))/NLYR
                    !** VER 7300
                    IF( VER >= 7300 .AND. IEVAP > 1 )THEN
                        READ(UWCI) TMPVAL,TMPVAL  !REAL(EVAPT(L),4),REAL(RAINT(L),4)
                    ENDIF

                ENDIF

                IF(ISTRAN(3) == 1) THEN
                    READ(UWCI)(DYE(L,K),K=KSZ(L),KC)
                    DYEA(L) = SUM(DYE(L,KSZ(L):KC))/NLYR
                ENDIF

                IF(ISTRAN(4) == 1) THEN
                    READ(UWCI)(SFL(L,K),K=KSZ(L),KC)
                    SFLA(L) = SUM(SFL(L,KSZ(L):KC))/NLYR
                ENDIF

                IF(ISTRAN(5) == 1)THEN
                    READ(UWCI)(TOXB(L,N1,NT),NT=1,NTOX)
                    READ(UWCI)((TOX(L,K,NT),K=KSZ(L),KC),NT=1,NTOX)
                    DO NT=1,NTOX
                        TOXA(L,NT) = SUM(TOX(L,KSZ(L):KC,NT))/NLYR
                    ENDDO
                ENDIF

                IF(ISTRAN(6) == 1 .OR. ISTRAN(7) >= 1)THEN
                    READ(UWCI)BELV(L),HBED(L,N1),BDENBED(L,N1),PORBED(L,N1)

                    IF(ISTRAN(6) == 1)THEN
                        READ(UWCI)(SEDB(L,N1,NS),VFRBED(L,N1,NS),NS=1,NSED)
                        READ(UWCI)((SED(L,K,NS),K=KSZ(L),KC),NS=1,NSED)
                        DO NS=1,NSED
                            SEDA(L,NS) = SUM(SED(L,KSZ(L):KC,NS))/NLYR
                        ENDDO
                    ENDIF

                    IF(ISTRAN(7) == 1)THEN
                        READ(UWCI)(SNDB(L,N1,NX),VFRBED(L,N1,NX+NSED),NX=1,NSND)
                        READ(UWCI)((SND(L,K,NX),K=KSZ(L),KC),NX=1,NSND)
                        DO NX=1,NSND
                            SNDA(L,NX) = SUM(SND(L,KSZ(L):KC,NX))/NLYR
                        ENDDO

                        IF(ISBDLDBC > 0)THEN
                            READ(UWCI)(CQBEDLOADX(L,NX),CQBEDLOADY(L,NX),NX=1,NSND)
                        ENDIF
                    ENDIF
                ENDIF

                IF(VER >= 7300 ) THEN
                    IF (ISGWIE > 0 )  READ(UWCI) EVAPSW(L),EVAPGW(L),RIFTR(L),AGWELV(L)
                    IF (ISICE >=3) READ(UWCI) RICETHKL(L),RICECOVL(L)
                ENDIF

            ENDDO
        ENDIF

        IF (LAYK>0) THEN  !----------------------------------------------------------------------------
            IF (ALLSNAP==0.AND.DTIME1<JULTIME.AND.JULTIME<=DTIME2) THEN
                !LAYER & PROFILE
                IF(TECCHKL) CALL TECPLOT(EETIME,3,0)
                IF(ISTRAN(1)==1) THEN
                    WRITE(UOUTI(1),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(1),'(A28,F15.5,A7)') '** (SAL(L,'//SSTR2//'),L=2,LA);TIME=',EETIME,'  [PPT]'
                    WRITE(UOUTI(1),'(10F10.5)') (SAL(L,LAYK),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('SAL',EETIME,NTS)
                ENDIF
                IF(ISTRAN(2)==1) THEN
                    WRITE(UOUTI(2),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(2),'(A28,F15.5,A5)') '** (TEM(L,'//SSTR2//'),L=2,LA);TIME=',EETIME,'  [C]'
                    WRITE(UOUTI(2),'(10F10.5)') (TEM(L,LAYK),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('TEM',EETIME,NTS)
                ENDIF
                IF(ISTRAN(3)==1) THEN
                    WRITE(UOUTI(3),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(3),'(A28,F15.5,A8)') '** (DYE(L,'//SSTR2//'),L=2,LA);TIME=',EETIME,'  [MG/L]'
                    WRITE(UOUTI(3),'(10F10.5)') (DYE(L,LAYK),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('DYE',EETIME,NTS)
                ENDIF
                IF(ISTRAN(4)==1) THEN
                    WRITE(UOUTI(4),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(4),'(A28,F15.5)') '** (SFL(L,'//SSTR2//'),L=2,LA);TIME=',EETIME
                    WRITE(UOUTI(4),'(10F10.5)') (SFL(L,LAYK),L=2,LA)
                ENDIF
                IF(ISTRAN(5)==1) THEN
                    WRITE(UOUTI(5),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER K = ',LAYK,'NTOX =',NTOX
                    WRITE(UOUTI(5),'(A45,F15.5,A8)') '** DO NT=1,NTOX (TOX(L,'//SSTR2//',NT),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
                    DO NT=1,NTOX
                        WRITE(UOUTI(5),'(A7,I5)') '** NT =',NT
                        WRITE(UOUTI(5),'(10F10.5)') (TOX(L,LAYK,NT),L=2,LA)
                    ENDDO
                    IF (VPROFCHK) CALL VPROFOUT('TOX',EETIME,NTS)
                ENDIF
                IF(ISTRAN(6)==1) THEN
                    WRITE(UOUTI(6),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER K = ',LAYK,'NSED =',NSED
                    WRITE(UOUTI(6),'(A50,F15.5,A8)') '** DO NS=1,NSED (SED(L,'//SSTR2//',NS),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
                    DO NS=1,NSED
                        WRITE(UOUTI(6),'(A7,I5)') '** NS =',NS
                        WRITE(UOUTI(6),'(10F10.5)') (SED(L,LAYK,NS),L=2,LA)
                    ENDDO
                    IF (VPROFCHK) CALL VPROFOUT('SED',EETIME,NTS)
                ENDIF
                IF(ISTRAN(7)==1) THEN
                    WRITE(UOUTI(7),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER K = ',LAYK,'NSND =',NSND
                    WRITE(UOUTI(7),'(A50,F15.5,A8)') '** DO NX=1,NSND (SND(L,'//SSTR2//',NX),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
                    DO NX=1,NSND
                        WRITE(UOUTI(7),'(A7,I5)') '** NX =',NX
                        WRITE(UOUTI(7),'(10F10.5)') (SND(L,LAYK,NX),L=2,LA)
                    ENDDO
                    IF (VPROFCHK) CALL VPROFOUT('SND',EETIME,NTS)
                ENDIF

                IF (KB==1) CALL BEDOUT(2,EETIME)

            ENDIF

            IF(ALLSNAP==1) THEN
                IF (TECCHKL) CALL TECPLOT(EETIME,3,0)
                IF(ISTRAN(1)==1) THEN
                    WRITE(UOUTI(1),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(1),'(A28,F15.5,A7)') '** (SAL(L,'//SSTR2//'),L=2,LA);TIME=',EETIME,'  [PPT]'
                    WRITE(UOUTI(1),'(10F10.5)') (SAL(L,LAYK),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('SAL',EETIME,NTS)
                ENDIF
                IF(ISTRAN(2)==1) THEN
                    WRITE(UOUTI(2),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(2),'(A28,F15.5,A5)') '** (TEM(L,'//SSTR2//'),L=2,LA);TIME=',EETIME,'  [C]'
                    WRITE(UOUTI(2),'(10F10.5)') (TEM(L,LAYK),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('TEM',EETIME,NTS)
                ENDIF
                IF(ISTRAN(3)==1) THEN
                    WRITE(UOUTI(3),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(3),'(A28,F15.5,A8)') '** (DYE(L,'//SSTR2//'),L=2,LA);TIME=',EETIME,'  [MG/L]'
                    WRITE(UOUTI(3),'(10F10.5)') (DYE(L,LAYK),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('DYE',EETIME,NTS)
                ENDIF
                IF(ISTRAN(4)==1) THEN
                    WRITE(UOUTI(4),'(A8,I5,2X,A9,I5)') '** LA = ',LA,'LAYER K = ',LAYK
                    WRITE(UOUTI(4),'(A28,F15.5)') '** (SFL(L,'//SSTR2//'),L=2,LA);TIME=',EETIME
                    WRITE(UOUTI(4),'(10F10.5)') (SFL(L,LAYK),L=2,LA)
                ENDIF
                IF(ISTRAN(5)==1) THEN
                    WRITE(UOUTI(5),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER K = ',LAYK,'NTOX =',NTOX
                    WRITE(UOUTI(5),'(A45,F15.5,A8)') '** DO NT=1,NTOX (TOX(L,'//SSTR2//',NT),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
                    DO NT=1,NTOX
                        WRITE(UOUTI(5),'(A7,I5)') '** NT =',NT
                        WRITE(UOUTI(5),'(10F10.3)') (TOX(L,LAYK,NT),L=2,LA)
                    ENDDO
                    IF (VPROFCHK) CALL VPROFOUT('TOX',EETIME,NTS)
                ENDIF
                IF(ISTRAN(6)==1) THEN
                    WRITE(UOUTI(6),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER K = ',LAYK,'NSED =',NSED
                    WRITE(UOUTI(6),'(A50,F15.5,A8)') '** DO NS=1,NSED (SED(L,'//SSTR2//',NS),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
                    DO NS=1,NSED
                        WRITE(UOUTI(6),'(A7,I5)') '** NS =',NS
                        WRITE(UOUTI(6),'(10F10.3)') (SED(L,LAYK,NS),L=2,LA)
                    ENDDO
                    IF (VPROFCHK) CALL VPROFOUT('SED',EETIME,NTS)
                ENDIF
                IF(ISTRAN(7)==1) THEN
                    WRITE(UOUTI(7),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER K = ',LAYK,'NSND =',NSND
                    WRITE(UOUTI(7),'(A50,F15.5,A8)') '** DO NX=1,NSND (SND(L,'//SSTR2//',NX),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
                    DO NX=1,NSND
                        WRITE(UOUTI(7),'(A7,I5)') '** NX =',NX
                        WRITE(UOUTI(7),'(10F10.3)') (SND(L,LAYK,NX),L=2,LA)
                    ENDDO
                    IF (VPROFCHK) CALL VPROFOUT('SND',EETIME,NTS)
                ENDIF
                IF (KB==1) CALL BEDOUT(2,EETIME)
            ENDIF

            !TIME SERIES
            IF(ISTRAN(1)==1) WRITE(UOUTI(11),SFMT1) EETIME,(SAL(LIJ(ICEL(N),JCEL(N)),LAYK),N=1,NLOC)
            IF(ISTRAN(2)==1) WRITE(UOUTI(12),SFMT1) EETIME,(TEM(LIJ(ICEL(N),JCEL(N)),LAYK),N=1,NLOC)
            IF(ISTRAN(3)==1) WRITE(UOUTI(13),SFMT1) EETIME,(DYE(LIJ(ICEL(N),JCEL(N)),LAYK),N=1,NLOC)
            IF(ISTRAN(4)==1) WRITE(UOUTI(14),SFMT1) EETIME,(SFL(LIJ(ICEL(N),JCEL(N)),LAYK),N=1,NLOC)

            IF(ISTRAN(5)==1) THEN
                WRITE(UOUTI(15),SFMT2) EETIME,((TOX(LIJ(ICEL(N),JCEL(N)),LAYK,NT),N=1,NLOC),NT=1,NTOX)
            ENDIF
            IF(ISTRAN(6)==1) THEN
                WRITE(UOUTI(16),SFMT3) EETIME,((SED(LIJ(ICEL(N),JCEL(N)),LAYK,NS),N=1,NLOC),NS=1,NSED)
            ENDIF
            IF(ISTRAN(7)==1) THEN
                WRITE(UOUTI(17),SFMT4) EETIME,((SND(LIJ(ICEL(N),JCEL(N)),LAYK,NX),N=1,NLOC),NX=1,NSND)
            ENDIF

        ELSEIF(LAYK==0) THEN  !-----------------------------------------------------------------------------
            IF (ALLSNAP==0.AND.DTIME1<JULTIME.AND.JULTIME<=DTIME2) THEN
                !LAYER & PROFILE FOR DEPTH-AVERAGED DAT AT ONE TIME MOMENT
                IF(TECCHKL) CALL TECPLOT(EETIME,31,0)
                IF(ISTRAN(1)==1) THEN
                    WRITE(UOUTI(1),'(A41,F15.5)') '** ((SALA(L),L=2,LA) [PPT],TIME=',EETIME
                    WRITE(UOUTI(1),'(10F10.5)') (SALA(L),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('SAL',EETIME,NTS)
                ENDIF
                IF(ISTRAN(2)==1) THEN
                    WRITE(UOUTI(2),'(A39,F15.5)') '** ((TEMA(L),L=2,LA) [C],TIME=',EETIME
                    WRITE(UOUTI(2),'(10F10.5)') (TEMA(L),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('TEM',EETIME,NTS)
                ENDIF
                IF(ISTRAN(3)==1) THEN
                    WRITE(UOUTI(3),'(A42,F15.5)') '** ((DYEA(L),L=2,LA) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(3),'(10F10.5)') (DYEA(L),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('DYE',EETIME,NTS)
                ENDIF
                IF(ISTRAN(4)==1) THEN
                    WRITE(UOUTI(4),'(A35,F15.5)') '** ((SFLA(L),L=2,LA),TIME=',EETIME
                    WRITE(UOUTI(4),'(10F10.5)') (SFLA(L),L=2,LA)
                ENDIF
                IF(ISTRAN(5)==1) THEN
                    WRITE(UOUTI(5),'(A57,F15.5)') '** ((TOXA(L,NT),L=2,LA),NT=1,NTOX) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(5),'(10F10.3)') ((TOXA(L,NT),L=2,LA),NT=1,NTOX)
                    IF (VPROFCHK) CALL VPROFOUT('TOX',EETIME,NTS)
                ENDIF
                IF(ISTRAN(6)==1) THEN
                    WRITE(UOUTI(6),'(A57,F15.5)') '** ((SEDA(L,NS),L=2,LA),NS=1,NSED) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(6),'(10F10.3)') ((SEDA(L,NS),L=2,LA),NS=1,NSED)
                    IF (VPROFCHK) CALL VPROFOUT('SED',EETIME,NTS)
                ENDIF
                IF(ISTRAN(7)==1) THEN
                    WRITE(UOUTI(7),'(A57,F15.5)') '** ((SNDA(L,NX),L=2,LA),NX=1,NSND) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(7),'(10F10.3)') ((SNDA(L,NX),L=2,LA),NX=1,NSND)
                    IF (VPROFCHK) CALL VPROFOUT('SND',EETIME,NTS)
                ENDIF
            ENDIF

            !TIME SERIES
            IF(ALLSNAP==1) THEN
                IF(TECCHKL) CALL TECPLOT(EETIME,31,0)
                IF(ISTRAN(1)==1) THEN
                    WRITE(UOUTI(1),'(A41,F15.5)') '** ((SALA(L),L=2,LA) [PPT],TIME=',EETIME
                    WRITE(UOUTI(1),'(10F10.5)') (SALA(L),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('SAL',EETIME,NTS)
                ENDIF
                IF(ISTRAN(2)==1) THEN
                    WRITE(UOUTI(2),'(A39,F15.5)') '** ((TEMA(L),L=2,LA) [C],TIME=',EETIME
                    WRITE(UOUTI(2),'(10F10.5)') (TEMA(L),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('TEM',EETIME,NTS)
                ENDIF
                IF(ISTRAN(3)==1) THEN
                    WRITE(UOUTI(3),'(A42,F15.5)') '** ((DYEA(L),L=2,LA) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(3),'(10F10.5)') (DYEA(L),L=2,LA)
                    IF (VPROFCHK) CALL VPROFOUT('DYE',EETIME,NTS)
                ENDIF
                IF(ISTRAN(4)==1) THEN
                    WRITE(UOUTI(4),'(A35,F15.5)') '** ((SFLA(L),L=2,LA),TIME=',EETIME
                    WRITE(UOUTI(4),'(10F10.5)') (SFLA(L),L=2,LA)
                ENDIF
                IF(ISTRAN(5)==1) THEN
                    WRITE(UOUTI(5),'(A57,F15.5)') '** ((TOXA(L,NT),L=2,LA),NT=1,NTOX) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(5),'(10F10.3)') ((TOXA(L,NT),L=2,LA),NT=1,NTOX)
                    IF (VPROFCHK) CALL VPROFOUT('TOX',EETIME,NTS)
                ENDIF
                IF(ISTRAN(6)==1) THEN
                    WRITE(UOUTI(6),'(A57,F15.5)') '** ((SEDA(L,NS),L=2,LA),NS=1,NSED) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(6),'(10F10.3)') ((SEDA(L,NS),L=2,LA),NS=1,NSED)
                    IF (VPROFCHK) CALL VPROFOUT('SED',EETIME,NTS)
                ENDIF
                IF(ISTRAN(7)==1) THEN
                    WRITE(UOUTI(7),'(A57,F15.5)') '** ((SNDA(L,NX),L=2,LA),NX=1,NSND) [MG/L],TIME=',EETIME
                    WRITE(UOUTI(7),'(10F10.3)') ((SNDA(L,NX),L=2,LA),NX=1,NSND)
                    IF (VPROFCHK) CALL VPROFOUT('SND',EETIME,NTS)
                ENDIF
            ENDIF

            !TIME SERIES
            IF(ISTRAN(1)==1) THEN
                !LAYER OF TIME FOR AVER
                WRITE(UOUTI(1),'(A41,F15.5)') '** ((SALA(L,K),L=2,LA),K=1,KC) [PPT],TIME=',EETIME
                WRITE(UOUTI(1),'(10F10.5)') (SALA(L),L=2,LA)
                WRITE(UOUTI(11),SFMT1) EETIME,(SALA(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)
            ENDIF
            IF(ISTRAN(2)==1) THEN
                WRITE(UOUTI(2),'(A39,F15.5)') '** ((TEMA(L),L=2,LA) [C],TIME=',EETIME
                WRITE(UOUTI(2),'(10F10.5)') (TEMA(L),L=2,LA)
                WRITE(UOUTI(12),SFMT1) EETIME,(TEMA(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)
            ENDIF
            IF(ISTRAN(3)==1) THEN
                WRITE(UOUTI(3),'(A42,F15.5)') '** ((DYEA(L),L=2,LA) [MG/L],TIME=',EETIME
                WRITE(UOUTI(3),'(10F10.5)') (DYEA(L),L=2,LA)
                WRITE(UOUTI(13),SFMT1) EETIME,(DYEA(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)
            ENDIF
            IF(ISTRAN(4)==1) THEN
                WRITE(UOUTI(4),'(A35,F15.5)') '** ((SFLA(L),L=2,LA),TIME=',EETIME
                WRITE(UOUTI(4),'(10F10.5)') (SFLA(L),L=2,LA)
                WRITE(UOUTI(14),SFMT1) EETIME,(SFLA(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)
            ENDIF
            IF(ISTRAN(5)==1) THEN
                WRITE(UOUTI(5),'(A57,F15.5)') '** (((TOXA(L,NT),L=2,LA),NT=1,NTOX) [MG/L],TIME=',EETIME
                WRITE(UOUTI(5),'(10F10.3)') ((TOXA(L,NT),L=2,LA),NT=1,NTOX)
                WRITE(UOUTI(15),SFMT2) EETIME,((TOXA(LIJ(ICEL(N),JCEL(N)),NT),N=1,NLOC),NT=1,NTOX)
            ENDIF
            IF(ISTRAN(6)==1) THEN
                WRITE(UOUTI(6),'(A57,F15.5)') '** (((SEDA(L,NS),L=2,LA),NS=1,NSED) [MG/L],TIME=',EETIME
                WRITE(UOUTI(6),'(10F10.3)') ((SEDA(L,NS),L=2,LA),NS=1,NSED)
                WRITE(UOUTI(16),SFMT3) EETIME,((SEDA(LIJ(ICEL(N),JCEL(N)),NS),N=1,NLOC),NS=1,NSED)
            ENDIF
            IF(ISTRAN(7)==1) THEN
                WRITE(UOUTI(7),'(A57,F15.5)') '** (((SNDA(L,NX),L=2,LA),NX=1,NSND) [MG/L],TIME=',EETIME
                WRITE(UOUTI(7),'(10F10.3)') ((SNDA(L,NX),L=2,LA),NX=1,NSND)
                WRITE(UOUTI(17),SFMT4) EETIME,((SNDA(LIJ(ICEL(N),JCEL(N)),NX),N=1,NLOC),NX=1,NSND)
            ENDIF

        ELSEIF (LAYK==-2) THEN

            DO N=1,NLOC
                IF (ZOPT==1) THEN
                    ZABV = HPT(LIJ(ICEL(N),JCEL(N)),NTS) - ZINT(N)  !CONVERT DEPTH UNDER WS TO HEIGHT ABOVE BED
                ELSEIF(ZOPT==2) THEN
                    ZABV = ZINT(N)                             !THE HEIGHT ABOVE BED
                ENDIF
                ZSIG(N) = ZABV/HPT(LIJ(ICEL(N),JCEL(N)),NTS)
                DO K=1,KC
                    IF (ZSIG(N) <= SUM(DZCK(1:K))) THEN
                        EXIT
                    ENDIF
                ENDDO
                K = MIN(K,KC)
                ZK2 = SUM(DZCK(1:K))
                ZK1 = ZK2 - DZCK(K)
                IF (ISTRAN(1)==1) THEN
                    VKT = SAL(LIJ(ICEL(N),JCEL(N)),K)
                    VKB = SAL(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1))
                    SALZ(N) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                ENDIF
                IF (ISTRAN(2)==1) THEN
                    VKT = TEM(LIJ(ICEL(N),JCEL(N)),K)
                    VKB = TEM(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1))
                    TEMZ(N) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                ENDIF
                IF (ISTRAN(2)==1) THEN
                    VKT = TEM(LIJ(ICEL(N),JCEL(N)),K)
                    VKB = TEM(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1))
                    DYEZ(N) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                ENDIF
                IF (ISTRAN(3)==1) THEN
                    VKT = DYE(LIJ(ICEL(N),JCEL(N)),K)
                    VKB = DYE(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1))
                    DYEZ(N) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                ENDIF
                IF (ISTRAN(5)==1) THEN
                    DO NT=1,NTOX
                        VKT = TOX(LIJ(ICEL(N),JCEL(N)),K,NT)
                        VKB = TOX(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1),NT)
                        TOXZ(N,NT) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                    ENDDO
                ENDIF
                IF (ISTRAN(6)==1) THEN
                    DO NS=1,NSED
                        VKT = SED(LIJ(ICEL(N),JCEL(N)),K,NS)
                        VKB = SED(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1),NS)
                        SEDZ(N,NS) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                    ENDDO
                ENDIF
                IF (ISTRAN(7)==1) THEN
                    DO NS=1,NSND
                        VKT = SND(LIJ(ICEL(N),JCEL(N)),K,NS)
                        VKB = SND(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1),NS)
                        SNDZ(N,NS) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                    ENDDO
                ENDIF

            ENDDO

            !TIME SERIES
            IF(ISTRAN(1)==1) WRITE(UOUTI(11),SFMT1) EETIME,(SALZ(N),N=1,NLOC)
            IF(ISTRAN(2)==1) WRITE(UOUTI(12),SFMT1) EETIME,(TEMZ(N),N=1,NLOC)
            IF(ISTRAN(3)==1) WRITE(UOUTI(13),SFMT1) EETIME,(DYEZ(N),N=1,NLOC)
            IF(ISTRAN(5)==1) WRITE(UOUTI(15),SFMT2) EETIME,((TOXZ(N,NT),N=1,NLOC),NT=1,NTOX)
            IF(ISTRAN(6)==1) WRITE(UOUTI(16),SFMT3) EETIME,((SEDZ(N,NS),N=1,NLOC),NS=1,NSED)
            IF(ISTRAN(7)==1) WRITE(UOUTI(17),SFMT4) EETIME,((SNDZ(N,NX),N=1,NLOC),NX=1,NSND)

        ENDIF

        ! ** ICE TS WRITTING
        IF (ISICE >=3) WRITE(UICE,SFMT1) EETIME,(RICETHKL(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)

    ENDDO !FINISHED EE_WC.OUT

100 CLOSE(UWCI)
    DO N=1,7
        IF(ISTRAN(N)==1) THEN
            CLOSE(UOUTI(N))
            CLOSE(UOUTI(N+10))
            CLOSE(UOUTI(N*10+1))
            IF (VPROFCHK) CLOSE(UCVPRF(N  ))
        ENDIF
    ENDDO
    IF(TECCHKL) CLOSE(UTEC3)
    INQUIRE(UOUTI(1),OPENED=ASW)
    IF (ASW) CLOSE(UOUTI(1))
    !CALL FREE_WC
    RETURN
999 STOP ' **** OPENING EE_WC.OUT ERROR'
    END SUBROUTINE

    SUBROUTINE GETEE_WQ
    INTEGER(4)::NW,L,K,N,N1,I,J,IOS,NT,NLYR
    INTEGER(4) IWQ(40), NACTIVE,LA1,HSIZE,BSIZE,IGRIDV,CELL3D
    REAL(4)::EETIME,WQ,DTIME1,DTIME2
    REAL(4)   ::ZK1,ZK2,VKT,VKB,ZABV
    REAL(8)   :: DTIME8
    CHARACTER(200)::SS*5,VPROWQF

    ! *** WATER QUALITY MODEL (HEM3D) RESULTS
    OPEN(UWQI,FILE=WQFILEI,FORM='BINARY')
    WQFILEO  = '#output\RESULT\WQ_TSK_'//TRIM(SSTR2)//'_DOM.DAT'
    IF (LAYK>=0) OPEN(UWQO,FILE=WQFILEO,ACTION='WRITE')

    !IF (LAYK>0) THEN
    !  1) CHC - cyanobacteria
    !  2) CHD - diatom algae
    !  3) CHG - green algae
    !  4) ROC - refractory particulate organic carbon
    !  5) LOC - labile particulate organic carbon
    !  6) DOC - dissolved organic carbon
    !  7) ROP - refractory particulate organic phosphorus
    !  8) LOP - labile particulate organic phosphorus
    !  9) DOP - dissolved organic phosphorus
    ! 10) P4D - total phosphate
    ! 11) RON - refractory particulate organic nitrogen 22) macroalgae
    ! 12) LON - labile particulate organic nitrogen
    ! 13) DON - dissolved organic nitrogen
    ! 14) NHX - ammonia nitrogen
    ! 15) NOX - nitrate nitrogen
    ! 16) SUU - particulate biogenic silica
    ! 17) SAA - dissolved available silica
    ! 18) COD - chemical oxygen demand
    ! 19) DOX - dissolved oxygen
    ! 20) TAM - total active metal
    ! 21) FCB - fecal coliform bacteria

    CHCTSF = '#output\RESULT\CHC_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    CHDTSF = '#output\RESULT\CHD_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    CHGTSF = '#output\RESULT\CHG_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    ROCTSF = '#output\RESULT\ROC_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    LOCTSF = '#output\RESULT\LOC_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    DOCTSF = '#output\RESULT\DOC_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    ROPTSF = '#output\RESULT\ROP_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    LOPTSF = '#output\RESULT\LOP_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    DOPTSF = '#output\RESULT\DOP_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    P4DTSF = '#output\RESULT\P4D_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    RONTSF = '#output\RESULT\RON_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    LONTSF = '#output\RESULT\LON_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    DONTSF = '#output\RESULT\DON_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    NHXTSF = '#output\RESULT\NHX_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    NOXTSF = '#output\RESULT\NOX_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    SUUTSF = '#output\RESULT\SUU_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    SAATSF = '#output\RESULT\SAA_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    CODTSF = '#output\RESULT\COD_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    DOXTSF = '#output\RESULT\DOX_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    TAMTSF = '#output\RESULT\TAM_TSK_'//TRIM(SSTR2)//'_CEL.DAT'
    FCBTSF = '#output\RESULT\FCB_TSK_'//TRIM(SSTR2)//'_CEL.DAT'

    TSWQF(1:21) = (/CHCTSF,CHDTSF,CHGTSF,ROCTSF,LOCTSF,DOCTSF,ROPTSF,LOPTSF,&
        DOPTSF,P4DTSF,RONTSF,LONTSF,DONTSF,NHXTSF,NOXTSF,SUUTSF,SAATSF,CODTSF,&
        DOXTSF,TAMTSF,FCBTSF /)

    !IF (LAYK>0) THEN
    DO N=1,21
        OPEN(UOUTI(N),FILE=TSWQF(N),ACTION='WRITE')
        WRITE(UOUTI(N),'(A12,I5)') '** LAYER K =',LAYK
        DO N1=1,NLOC
            WRITE(UOUTI(N),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N1),'J =',JCEL(N1)
        ENDDO
        IF (N<21) THEN
            WRITE(UOUTI(N),'(A22,I2,A16)') '** EETIME,(WQV(L(N),'//SSTR2//',',N,'),N=1,'//SSTR3//') [MG/L]'
        ELSEIF(N==21) THEN
            WRITE(UOUTI(N),'(A22,I2,A21)') '** EETIME,(WQV(L(N),'//SSTR2//',',N,'),N=1,'//SSTR3//') [MPN/100ML]'
        ENDIF
    ENDDO
    IF (TECCHKL) CALL TECPLOT(0.,4,1)
    !ENDIF

    WRITE(SS,'(I5)') NLOC
    SFMT='(F12.4,'//SS//'F10.3)'
    DO N=1,21
        IF (VPROFCHK) THEN
            VPROWQF = TSWQF(N)(1:19)//'_PROF.DAT'
            OPEN(UCVPRF(N),FILE=VPROWQF,ACTION='WRITE')
            WRITE(UCVPRF(N),'(A)') '** VERTICAL PROFILE OF '//TSWQF(N)(16:18)
            IF (VPROFILE==1) THEN
                WRITE(UCVPRF(N),'(A)') '** (Z_ABOVE_BOT, VALUE),N=1:NLOC)'
            ELSEIF(VPROFILE==2) THEN
                WRITE(UCVPRF(N),'(A)') '** (ELEVATION, VALUE),N=1:NLOC)'
            ENDIF
            DO N1=1,NLOC
                WRITE(UCVPRF(N),'(A6,I5,A4,I5)') '** I =',ICEL(N1),' J =',JCEL(N1)
            ENDDO
        ENDIF
    ENDDO
    ! ***************

    READ(UWQI,IOSTAT=IOS)VER_WQ
    IF (IOS /= 0) GOTO 999

    IF (VER_WQ>=8400) THEN  !2018-07-24, NTL:
        READ(UWQI) HSIZE,BSIZE,IGRIDV,CELL3D
        READ(UWQI) IC,JC,KC,KB,LA1
        READ(UWQI) NWQV
        READ(UWQI) (IWQ(NW),NW=1,NWQV)
    ELSEIF(VER_WQ>=102) THEN
        READ(UWQI) HSIZE
        READ(UWQI) KC,KB,LA1
        READ(UWQI) NWQV
        READ(UWQI) (IWQ(NW),NW=1,NWQV)
    ELSEIF (VER_WQ<102) THEN
        READ(UWQI)NWQV
        IF (VER_WQ==100) READ(UWQI)(ISTRWQ(NW),NW=1,NWQV)   !VERSION 06_OMP
        READ(UWQI)(IWQ(NW),NW=1,NWQV)
    ENDIF

    ALLOCATE(WQV(LA,KC,NTSWQVM))
    ALLOCATE(WQVA(LA,NTSWQVM))
    ALLOCATE(WQVZ(NLOC,NTSWQVM))

    WQV  = BLK
    DTIME2=0
    NT=0
    DO WHILE(1)

        IF (VER_WQ < 7300) THEN
            READ(UWQI,END=200) EETIME
        ELSE
            READ(UWQI,END=200) DTIME8
            EETIME = DTIME8
        ENDIF

        DTIME1=DTIME2
        DTIME2=EETIME
        NT=NT+1
        IF (VER_WQ>=8400) THEN  !2018-07-24, NTL:
            DO NW=1,NWQV
                IF(IWQ(NW) > 0)THEN
                    READ(UWQI) ((WQV(L,K,NW),K=KSZ(L),KC),L=2,LA)
                ELSE
                    WQV(:,:,NW) = 0.0
                ENDIF
            ENDDO

            DO L=2,LA
                NLYR = KC-KSZ(L)+1
                DO NW=1,NWQV
                    IF (IWQ(NW) == 0)THEN
                        DO K=KSZ(L),KC
                            WQV(L,K,NW)=0.0
                        ENDDO
                    ENDIF
                    WQVA(L,NW) = SUM(WQV(L,KSZ(L):KC,NW))/NLYR
                ENDDO
            ENDDO
        ELSE
            DO L=2,LA
                NLYR = KC-KSZ(L)+1
                DO K=KSZ(L),KC
                    DO NW=1,NWQV
                        IF(IWQ(NW) > 0)THEN
                            READ(UWQI) WQ
                            WQV(L,K,NW)=WQ
                        ELSE
                            WQV(L,K,NW)=0.0
                        ENDIF
                    ENDDO
                ENDDO
                DO NW=1,NWQV
                    WQVA(L,NW) = SUM(WQV(L,KSZ(L):KC,NW))/NLYR
                ENDDO
            ENDDO
        ENDIF

        IF (LAYK>0) THEN
            IF(ALLSNAP==0.AND.DTIME1<JULTIME.AND.JULTIME<=DTIME2) THEN
                IF(TECCHKL)   CALL TECPLOT(EETIME,4,0)
                WRITE(UWQO,'(A50,F15.5)') '** DO NW=1,NWQV (WQV(L,'//SSTR2//',NW),L=2,LA) ENDDO; TIME=',EETIME
                DO NW=1,NWQV
                    IF (VPROFCHK) CALL VPROFOUT_WQ(NW,EETIME,NT)
                    WRITE(UWQO,'(10F10.3)') (WQV(L,LAYK,NW),L=2,LA)
                ENDDO
            ENDIF

            IF(ALLSNAP==1) THEN
                IF(TECCHKL)  CALL TECPLOT(EETIME,4,0)
                WRITE(UWQO,'(A50,F15.5)') '** DO NW=1,NWQV (WQV(L,'//SSTR2//',NW),L=2,LA) ENDDO; TIME=',EETIME
                DO NW=1,NWQV
                    IF(VPROFCHK) CALL VPROFOUT_WQ(NW,EETIME,NT)
                    WRITE(UWQO,'(10F10.3)') (WQV(L,LAYK,NW),L=2,LA)
                ENDDO
            ENDIF

            !TIME SERIES
            DO NW=1,NWQV
                WRITE(UOUTI(NW),SFMT) EETIME,(WQV(LIJ(ICEL(N),JCEL(N)),LAYK,NW),N=1,NLOC)
            ENDDO

        ELSEIF(LAYK==0) THEN
            IF(ALLSNAP==0.AND.DTIME1<JULTIME.AND.JULTIME<=DTIME2) THEN
                IF(TECCHKL)   CALL TECPLOT(EETIME,41,0)
                WRITE(UWQO,'(A50,F15.5)') '** DO NW=1,NWQV (WQV(L,'//SSTR2//',NW),L=2,LA) ENDDO; TIME=',EETIME
                DO NW=1,NWQV
                    IF (VPROFCHK) CALL VPROFOUT_WQ(NW,EETIME,NT)
                    WRITE(UWQO,'(10F10.3)') (WQVA(L,NW),L=2,LA)
                ENDDO
            ENDIF

            IF(ALLSNAP==1) THEN
                IF(TECCHKL)  CALL TECPLOT(EETIME,41,0)
                WRITE(UWQO,'(A50,F15.5)') '** DO NW=1,NWQV (WQV(L,'//SSTR2//',NW),L=2,LA) ENDDO; TIME=',EETIME
                DO NW=1,NWQV
                    IF(VPROFCHK) CALL VPROFOUT_WQ(NW,EETIME,NT)
                    WRITE(UWQO,'(10F10.3)') (WQVA(L,NW),L=2,LA)
                ENDDO

            ENDIF

            !TIME SERIES
            DO NW=1,NWQV
                WRITE(UOUTI(NW),SFMT) EETIME,(WQVA(LIJ(ICEL(N),JCEL(N)),NW),N=1,NLOC)
            ENDDO

        ELSEIF (LAYK==-2) THEN
            DO N=1,NLOC
                IF (ZOPT==1) THEN
                    ZABV = HPT(LIJ(ICEL(N),JCEL(N)),NT) - ZINT(N)  !CONVERT DEPTH UNDER WS TO HEIGHT ABOVE BED
                ELSEIF(ZOPT==2) THEN
                    ZABV = ZINT(N)                             !THE HEIGHT ABOVE BED
                ENDIF
                ZSIG(N) = ZABV/HPT(LIJ(ICEL(N),JCEL(N)),NT)
                DO K=1,KC
                    IF (ZSIG(N) <= SUM(DZCK(1:K))) THEN
                        EXIT
                    ENDIF
                ENDDO
                K = MIN(K,KC)
                ZK2 = SUM(DZCK(1:K))
                ZK1 = ZK2 - DZCK(K)
                DO NW=1,NWQV
                    VKT = WQV(LIJ(ICEL(N),JCEL(N)),K,NW)
                    VKB = WQV(LIJ(ICEL(N),JCEL(N)),MAX(K-1,1),NW)
                    WQVZ(N,NW) = (VKT -VKB )*(ZSIG(N)-ZK1)/(ZK2-ZK1)+VKB
                ENDDO
            ENDDO
            !TIME SERIES
            DO NW=1,NWQV
                WRITE(UOUTI(NW),SFMT) EETIME,(WQVZ(N,NW),N=1,NLOC)
            ENDDO

        ENDIF
    ENDDO
200 CLOSE(UWQI)
    PRINT '(A,I8)',' *** NUMBER OF WQ VARIABLES:',NWQV
    IF (LAYK>=0) CLOSE(UWQO)
    DO N=1,21
        CLOSE(UOUTI(N))
        IF (VPROFCHK) CLOSE(UCVPRF(N))
    ENDDO
    IF (LAYK>0) CLOSE(UTEC4)
    DEALLOCATE(WQVA)
    RETURN
999 STOP ' **** OPENING EE_WS.OUT ERROR'
    END SUBROUTINE

    SUBROUTINE VPROFOUT(ITEM,DTIME,NTI)
    CHARACTER(*),INTENT(IN)::ITEM
    REAL(4),     INTENT(IN)::DTIME
    INTEGER(4),  INTENT(IN)::NTI
    INTEGER(4)::L,N,K,NT
    REAL(4)::VAL(NLOC),ZZH(NLOC,KC)
    CHARACTER(20)::SFMT,SS*4

    WRITE(SS,'(I4)') 2*NLOC
    SFMT = '('//SS//'F12.5)'
    ZZH = 0
    IF( IGRIDV > 0 )THEN
        DO N=1,NLOC
            L=LIJ(ICEL(N),JCEL(N))
            ZZH(N,1) = HPT(L,NTI)*DZCS(L,1)
            DO K=2,KC
                ZZH(N,K) = ZZH(N,K-1)+HPT(L,NTI)*DZCS(L,K)
            ENDDO
            IF (VPROFILE==2) ZZH(N,KSZ(L):KC) = ZZH(N,KSZ(L):KC)+BELV(L)
        ENDDO
    ELSE
        DO N=1,NLOC
            L=LIJ(ICEL(N),JCEL(N))
            ZZH(N,1) = HPT(L,NTI)*DZCK(1)
            DO K=2,KC
                ZZH(N,K) = ZZH(N,K-1)+HPT(L,NTI)*DZCK(K)
            ENDDO
            IF (VPROFILE==2) ZZH(N,1:KC) = ZZH(N,1:KC)+BELV(L)
        ENDDO
    ENDIF

    IF(ITEM(1:3)=='VEL') THEN
        WRITE(UCVPRF(1),'(A10,F15.5)') '** TIME = ',DTIME
        DO K=KC,1,-1
            DO N=1,NLOC
                L=LIJ(ICEL(N),JCEL(N))
                VAL(N) = SQRT(UK(L,K)**2+VK(L,K)**2+W(L,K)**2)
            ENDDO
            WRITE(UCVPRF(1),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
        ENDDO
    ELSEIF(ITEM(1:3)=='SAL') THEN
        WRITE(UCVPRF(1),'(A10,F15.5)') '** TIME = ',DTIME
        DO K=KC,1,-1
            DO N=1,NLOC
                L=LIJ(ICEL(N),JCEL(N))
                VAL(N) = SAL(L,K)
            ENDDO
            WRITE(UCVPRF(1),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
        ENDDO
    ELSEIF(ITEM(1:3)=='TEM') THEN
        WRITE(UCVPRF(2),'(A10,F15.5)') '** TIME = ',DTIME
        DO K=KC,1,-1
            DO N=1,NLOC
                L=LIJ(ICEL(N),JCEL(N))
                VAL(N) = TEM(L,K)
            ENDDO
            WRITE(UCVPRF(2),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
        ENDDO
    ELSEIF(ITEM(1:3)=='DYE') THEN
        WRITE(UCVPRF(3),'(A10,F15.5)') '** TIME = ',DTIME
        DO K=KC,1,-1
            DO N=1,NLOC
                L=LIJ(ICEL(N),JCEL(N))
                VAL(N) = DYE(L,K)
            ENDDO
            WRITE(UCVPRF(3),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
        ENDDO
    ELSEIF(ITEM(1:3)=='TOX') THEN
        WRITE(UCVPRF(5),'(A10,F15.5)') '** TIME = ',DTIME
        DO NT=1,NTOX
            WRITE(UCVPRF(5),'(A7,I5)') '** NT =',NT
            DO K=KC,1,-1
                DO N=1,NLOC
                    L=LIJ(ICEL(N),JCEL(N))
                    VAL(N) = TOX(L,K,NT)
                ENDDO
                WRITE(UCVPRF(5),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
            ENDDO
        ENDDO
    ELSEIF(ITEM(1:3)=='SED') THEN
        WRITE(UCVPRF(6),'(A10,F15.5)') '** TIME = ',DTIME
        DO NT=1,NSED
            WRITE(UCVPRF(6),'(A9,I5)') '** NSED =',NT
            DO K=KC,1,-1
                DO N=1,NLOC
                    L=LIJ(ICEL(N),JCEL(N))
                    VAL(N) = SED(L,K,NT)
                ENDDO
                WRITE(UCVPRF(6),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
            ENDDO
        ENDDO
    ELSEIF(ITEM(1:3)=='SND') THEN
        WRITE(UCVPRF(7),'(A10,F15.5)') '** TIME = ',DTIME
        DO NT=1,NSND
            WRITE(UCVPRF(7),'(A9,I5)') '** NSND =',NT
            DO K=KC,1,-1
                DO N=1,NLOC
                    L=LIJ(ICEL(N),JCEL(N))
                    VAL(N) = SND(L,K,NT)
                ENDDO
                WRITE(UCVPRF(7),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
            ENDDO
        ENDDO

    ENDIF

    END SUBROUTINE

    SUBROUTINE VPROFOUT_WQ(NW,DTIME,NTI)
    INTEGER(4),  INTENT(IN)::NW
    REAL(4),     INTENT(IN)::DTIME
    INTEGER(4),  INTENT(IN)::NTI
    INTEGER(4)::L,N,K,NT
    REAL(4)::VAL(NLOC),ZZH(NLOC,KC)
    CHARACTER(20)::SFMT,SS*4

    WRITE(SS,'(I4)') 2*NLOC
    SFMT = '('//SS//'F12.5)'

    ZZH = 0
    IF( IGRIDV > 0 )THEN
        DO N=1,NLOC
            L=LIJ(ICEL(N),JCEL(N))
            ZZH(N,1) = HPT(L,NTI)*DZCS(L,1)
            DO K=2,KC
                ZZH(N,K) = ZZH(N,K-1)+HPT(L,NTI)*DZCS(L,K)
            ENDDO
            IF (VPROFILE==2) ZZH(N,KSZ(L):KC) = ZZH(N,KSZ(L):KC)+BELV(L)
        ENDDO
    ELSE
        DO N=1,NLOC
            L=LIJ(ICEL(N),JCEL(N))
            ZZH(N,1) = HPT(L,NTI)*DZCK(1)
            DO K=2,KC
                ZZH(N,K) = ZZH(N,K-1)+HPT(L,NTI)*DZCK(K)
            ENDDO
            IF (VPROFILE==2) ZZH(N,1:KC) = ZZH(N,1:KC)+BELV(L)
        ENDDO
    ENDIF

    WRITE(UCVPRF(NW),'(A10,F15.5)') '** TIME = ',DTIME
    DO K=KC,1,-1
        DO N=1,NLOC
            L=LIJ(ICEL(N),JCEL(N))
            VAL(N) = WQV(L,K,NW)
        ENDDO
        WRITE(UCVPRF(NW),SFMT) (ZZH(N,K),VAL(N),N=1,NLOC)
    ENDDO

    END SUBROUTINE

    SUBROUTINE FREE_WC
    DEALLOCATE(SAL)
    DEALLOCATE(SJDAY)
    DEALLOCATE(TEM)
    DEALLOCATE(DYE)
    DEALLOCATE(SFL)
    DEALLOCATE(TOX)
    DEALLOCATE(SED)
    DEALLOCATE(SND)
    DEALLOCATE(TEMB)
    DEALLOCATE(QQWV1)
    DEALLOCATE(TAUBSED)
    DEALLOCATE(TAUBSND)
    DEALLOCATE(TAUB)
    DEALLOCATE(WVWHA)
    DEALLOCATE(WVFRQL)
    DEALLOCATE(WACCWE)
    DEALLOCATE(WVDISP)
    DEALLOCATE(WVHUU)
    DEALLOCATE(WVHVV)
    DEALLOCATE(WVHUV)
    DEALLOCATE(CQBEDLOADX)
    DEALLOCATE(CQBEDLOADY)
    DEALLOCATE(VFRBED)
    DEALLOCATE(SALA)
    DEALLOCATE(TEMA)
    DEALLOCATE(DYEA)
    DEALLOCATE(SFLA)
    DEALLOCATE(TOXA)
    DEALLOCATE(SEDA)
    DEALLOCATE(SNDA)
    IF (ISICE >=3) THEN
        CLOSE(UICE)
        DEALLOCATE(RICECOVL,RICETHKL)
    ENDIF
    IF (KB == 1) THEN
        DEALLOCATE(SEDDIA)
        DEALLOCATE(KBT)
        DEALLOCATE(TOXB)
        DEALLOCATE(HBED)
        DEALLOCATE(BDENBED)
        DEALLOCATE(PORBED)
        DEALLOCATE(SEDB)
        DEALLOCATE(SNDB)
    ENDIF

    END SUBROUTINE

    SUBROUTINE GETLPT
    ! ** GET LPT FROM EFDC
    ! ** INPUT:
    ! ** EE_DRIFTER.OUT
    ! ** OUTPUT:
    ! ** EE_DRIFTER.DAT
    INTEGER(4)::VER,KC,LLA,NP,NACT,N,LENDIR,NUMC,NM,K,NS
    INTEGER(4):: IX,IY,IZ,IV,ITMP,IOS,NTMP
    REAL(4)::DTIME4,ZLA,DVOL,ZETA,TMP
    REAL(8)::XLA,YLA,DTIME8
    INTEGER(4),ALLOCATABLE::ULPT(:)
    CHARACTER(200),ALLOCATABLE:: OUTFILE(:)
    CHARACTER*80 TITLE,METHOD,drive*2,dir,name,ext,SS*6,STR*200
    REAL(8),ALLOCATABLE:: TTMP(:)

    CALL GETSTRNUM(NDRIFTER,NM)  !GET NSET

    ALLOCATE(OUTFILE(NM),ULPT(NM))
    ALLOCATE(TTMP(LPTSNP))

    OPEN(ULGR,FILE=LPTFILE,FORM='BINARY',ACTION='READ')

    READ(ULGR,IOSTAT=IOS) VER
    IF (IOS /= 0) GOTO 999
    READ(ULGR) TITLE
    READ(ULGR) METHOD

    ! **
    DO N=1,NM
        WRITE(SS,'(I6.6)') NSET(N)
        OUTFILE(N) = '#output\RESULT\DRIFTER_'//SS//'.DAT'
        ULPT(N) = 100+N
        OPEN(ULPT(N),FILE=OUTFILE(N),ACTION='WRITE')
        IF (VER <720) THEN
            WRITE(ULPT(N),'(A)') '**  JTIME        LLA         X             Y            Z        WS-Z     '
        ELSEIF (VER >= 720) THEN
            WRITE(ULPT(N),'(A)') '**  JTIME        LLA         X             Y            Z        WS-Z     VOL'
        ENDIF
    ENDDO
    ! **

    IF (VER==7300) THEN
        READ(ULGR) NPD,KC,XYZSCL
    ELSE
        READ(ULGR) NPD,KC
    ENDIF

    NTMP = 0
    DO WHILE (1)
        IF (VER==7300) THEN
            READ(ULGR,END=100) DTIME8
            DTIME4 = DTIME8
        ELSE
            READ(ULGR,END=100) DTIME4
            DTIME8 = DTIME4
        ENDIF

        NTMP = NTMP + 1
        TTMP(NTMP) = DTIME8

        READ(ULGR,END=100) NACT

        IF (VER==720) THEN
            DO N=1,NACT
                READ(ULGR,END=100) NP,LLA,XLA,YLA,ZLA,DVOL
                DO NS=1,NM
                    IF (NP==NSET(NS)) THEN
                        CALL DRIFTERWDEP(DTIME4,LLA,XLA,YLA,ZETA)
                        WRITE(ULPT(NS),'(F12.4,I8,2F15.3,2F10.3,F12.3)') DTIME4,LLA,XLA,YLA,ZLA,ZETA-ZLA,DVOL
                        EXIT
                    ENDIF
                ENDDO
            ENDDO
        ELSEIF(VER<720) THEN
            DO N=1,NACT
                READ(ULGR,END=100) NP,LLA,XLA,YLA,ZLA
                DO NS=1,NM
                    IF (NP==NSET(NS)) THEN
                        CALL DRIFTERWDEP(DTIME4,LLA,XLA,YLA,ZETA)
                        WRITE(ULPT(NS),'(F12.4,I8,2F15.3,2F10.3)') DTIME4,LLA,XLA,YLA,ZLA,ZETA-ZLA
                        EXIT
                    ENDIF
                ENDDO
            ENDDO
        ELSEIF(VER==7300) THEN
            IF (ANY(ISOILSPI==1)) THEN
                DO N=1,NACT
                    READ(ULGR,END=100) NP,LLA,IX,IY,IZ,IV
                    DO NS=1,NM
                        IF (NP==NSET(NS)) THEN
                            XLA = REAL(IX)/XYZSCL
                            YLA = REAL(IY)/XYZSCL
                            ZLA = REAL(IZ)/XYZSCL
                            DVOL= REAL(IV)*1D-6
                            CALL DRIFTERWDEP(DTIME4,LLA,XLA,YLA,ZETA)
                            WRITE(ULPT(NS),'(F12.4,I8,2F15.3,2F10.3,F12.3)') DTIME4,LLA,XLA,YLA,ZLA,ZETA-ZLA,DVOL
                            EXIT
                        ENDIF
                    ENDDO
                ENDDO
            ELSE
                DO N=1,NACT
                    READ(ULGR,END=100) NP,LLA,IX,IY,IZ
                    DO NS=1,NM
                        IF (NP==NSET(NS)) THEN
                            XLA = REAL(IX)/XYZSCL
                            YLA = REAL(IY)/XYZSCL
                            ZLA = REAL(IZ)/XYZSCL
                            CALL DRIFTERWDEP(DTIME4,LLA,XLA,YLA,ZETA)
                            WRITE(ULPT(NS),'(F12.4,I8,2F15.3,2F10.3)') DTIME4,LLA,XLA,YLA,ZLA,ZETA-ZLA
                            EXIT
                        ENDIF
                    ENDDO
                ENDDO
            ENDIF
        ENDIF
    ENDDO
100 LPTSNP = NTMP
    ALLOCATE(LPTTIME(LPTSNP))
    LPTTIME(1:LPTSNP) = TTMP(1:NTMP)
    DEALLOCATE(TTMP)
    CLOSE(ULGR)
    CLOSE(UOUT)
    RETURN
999 STOP ' **** OPENING EE_WS.OUT ERROR'
    END SUBROUTINE

    SUBROUTINE DRIFTERWDEP(TIMEDAY,LNI,XLA,YLA,ZETA)
    !INTERPOLATION OF THE TOTAL WATER DEPTH
    !FOR THE DRIFTER NI AT EACH TIME INSTANT AND EACH LOCATION
    INTEGER(4),INTENT(IN)::LNI
    REAL(8),   INTENT(IN)::XLA,YLA
    REAL(4),   INTENT(IN)::TIMEDAY
    REAL(4),  INTENT(OUT)::ZETA
    INTEGER(4)::NT
    REAL(RK)  ::HPL

    DO NT=1,NTM-1
        IF(TIMEDAY>=JTIME(NT).AND.TIMEDAY<=JTIME(NT+1)) THEN
            HPL = (HPT(LNI,NT+1)-HPT(LNI,NT))*(TIMEDAY-JTIME(NT))/(JTIME(NT+1)-JTIME(NT))+HPT(LNI,NT)
            EXIT
        ENDIF
    ENDDO

    ZETA = HPL+BELV(LNI)

    END SUBROUTINE

    SUBROUTINE GETEE_BC
    INTEGER(4)::VER,LINES,L,NT,N,K,HSIZE,NS,LL,ITMP,IOS,BSIZE

    INTEGER(IK4) :: NSXD,NBCCELLS,NCELLLIST,ISTAT
    INTEGER(IK4),ALLOCATABLE :: BCCELLS(:)

    REAL(8)::DTIME8
    REAL(4)::TIMESTEP,TMP
    CHARACTER(200)::SS*5

    WRITE(SS,'(I5)') 2*NLOC
    SFMT='(F12.4,'//SS//'F10.3)'
    OPEN(UOUTI(1),FILE=TSBCF,ACTION='WRITE')
    DO N=1,NLOC
        WRITE(UOUTI(1),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
    ENDDO
    WRITE(UOUTI(1),'(A38,I3,A3)') '** TIME,(QSUM(LN,LAYK),QSUME(LN),N=1, ',NLOC,')'

    OPEN(UINP,FILE=BCFILEI,FORM='BINARY')

    !READ(UINP,IOSTAT=IOS) VER,HSIZE,NBCS,NBCCELLS,NCELLLIST
    READ(UINP,IOSTAT=IOS) VER,HSIZE
    IF(VER>=8400) THEN
        READ(UINP,IOSTAT=IOS) BSIZE,ITMP,ITMP,ITMP,ITMP
    ENDIF
    READ(UINP,IOSTAT=IOS) NBCS,NBCCELLS,NCELLLIST
    IF (IOS /= 0) GOTO 999

    ALLOCATE(BCCELLS(NBCCELLS))
    READ(UINP) (BCCELLS(L),L=1,NCELLLIST)
    IF(VER>=8400) THEN
        READ(UINP,IOSTAT=IOS) NPBS,NPBW,NPBE,NPBN
        READ(UINP,IOSTAT=IOS) NQCTL,NQWR,NQCTLSER,NQCRULES,NGWSER,ISGWIT
    ENDIF

    ISTAT = FSEEK(UINP,HSIZE,0)


    NTM=0
    DO WHILE(1)
        IF( ISTAT /= 0 ) EXIT
        READ (UINP,END=100) DTIME8
        NTM=NTM+1

        ! *** OUTPUT SELECTIVE QSUM
        IF( KC > 1 )THEN
            DO L=2,LA
                READ(UINP) TMP
            ENDDO
            DO NS=1,NBCS
                READ(UINP) (TMP,K=1,KC)
            ENDDO
            IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
                DO L=2,LA
                    READ(UINP) TMP
                ENDDO
            ENDIF
        ELSE
            ! *** SINGLE LAYER
            DO L=2,LA
                READ(UINP) TMP
            ENDDO
        ENDIF

        ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
        DO LL=1,NPBS
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        DO LL=1,NPBW
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        DO LL=1,NPBE
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        DO LL=1,NPBN
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        IF (NQCTL > 0) THEN
            DO L=1,NQCTL
                DO NS=1,2
                    READ(UINP) (TMP,K=1,KC)
                ENDDO
            ENDDO
        ENDIF
        IF(NQCTLSER > 0 .OR. NQCRULES > 0) THEN
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
        ENDIF
        IF (NQWR > 0) THEN
            READ(UINP) (TMP,L=1,NQWR)
            READ(UINP) (TMP,L=1,NQWR)
        ENDIF
    ENDDO

100 REWIND(UINP)
    ALLOCATE(QSUM(LCM,KC),QSUME(LCM))
    QSUM = 0
    QSUME = 0
    
    !READ(UINP,IOSTAT=IOS) VER,HSIZE
    !!READ(UINP) VER,HSIZE,ITMP,NBCCELLS,NCELLLIST,(BCCELLS(L),L=1,NCELLLIST)
    !IF(VER>=8400) THEN
    !    READ(UINP,IOSTAT=IOS) BSIZE,ITMP,ITMP,ITMP,ITMP
    !ENDIF
    !READ(UINP,IOSTAT=IOS) NBCS,NBCCELLS,NCELLLIST
    !IF (IOS /= 0) GOTO 999
    !
    !ALLOCATE(BCCELLS(NBCCELLS))
    !READ(UINP) (BCCELLS(L),L=1,NCELLLIST)
    !IF(VER>=8400) THEN
    !    READ(UINP,IOSTAT=IOS) NPBS,NPBW,NPBE,NPBN
    !    READ(UINP,IOSTAT=IOS) NQCTL,NQWR,NQCTLSER,NQCRULES,NGWSER,ISGWIT
    !ENDIF

    ISTAT = FSEEK(UINP,HSIZE,0)
    
    CALL TECPLOT_TMP(0.0,1)
    DO NT=1,NTM
        IF( ISTAT /= 0 ) EXIT
        READ (UINP,END=200) DTIME8
        TIMESTEP  = DTIME8

        ! *** OUTPUT SELECTIVE QSUM
        IF( KC > 1 )THEN
            DO L=2,LA
                READ(UINP) QSUM(L,KC)
            ENDDO
            DO NS=1,NBCS
                L=BCCELLS(NS)
                IF(L < 2 .OR. L > LA) THEN
                    L=1
                ENDIF
                READ(UINP) (QSUM(L,K),K=1,KC)
            ENDDO
            IF( NGWSER > 0 .OR. ISGWIT /= 0 )THEN
                DO L=2,LA
                    READ(UINP) QSUM(L,KSZ(L))
                ENDDO
            ENDIF
        ELSE
            ! *** SINGLE LAYER
            DO L=2,LA
                READ(UINP) QSUME(L)
            ENDDO
        ENDIF

        ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES (OUTPUT IN LOBCS ORDER)
        DO LL=1,NPBS
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        DO LL=1,NPBW
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        DO LL=1,NPBE
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        DO LL=1,NPBN
            READ(UINP) (TMP,K=1,KC)
        ENDDO
        IF (NQCTL > 0) THEN
            DO L=1,NQCTL
                DO NS=1,2
                    READ(UINP) (TMP,K=1,KC)
                ENDDO
            ENDDO
        ENDIF
        IF(NQCTLSER > 0 .OR. NQCRULES > 0) THEN
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
            READ(UINP) (TMP,L=1,NQCTL)
        ENDIF

        IF (NQWR > 0) THEN
            READ(UINP) (TMP,L=1,NQWR)
            READ(UINP) (TMP,L=1,NQWR)
        ENDIF
        !** TIME SERIES
        IF(LAYK > 0) WRITE(UOUTI(1),SFMT) DTIME8,(QSUM(LIJ(ICEL(N),JCEL(N)),LAYK),QSUME(LIJ(ICEL(N),JCEL(N))),N=1,NLOC)
        ARR(2:LA) = QSUM(2:LA,LAYK)
        CALL TECPLOT_TMP(TIMESTEP,2)

    ENDDO
200 CLOSE(UINP)
    CLOSE(UOUTI(1))
    CLOSE(UTEC1)
    RETURN
999 PRINT *,' **** OPENING EE_BC.OUT ERROR OR IT DOES NOT EXIST'
    END SUBROUTINE

    SUBROUTINE READ_TMP
    INTEGER(4)::NT,L,IOS
    REAL(4)   ::TIMESTEP
    INPFILE  = trim(OUTFOLDER)//'TMP.DAT'
    OPEN(10,FILE=INPFILE,ACTION='READ',STATUS='OLD')

    NT = 1
    CALL TECPLOT_TMP(TIMESTEP,NT,1)
    DO WHILE(1)
        READ(10,*,END=100,IOSTAT=IOS) TIMESTEP
        IF (IOS /= 0) GOTO 999
        READ(10,*) (ARR(L),L=2,LA)
        NT=NT+1
        CALL TECPLOT_TMP(TIMESTEP,NT)
    ENDDO
100 CLOSE(10)
    CLOSE(UOUT)
    RETURN
999 PRINT *,' **** OPENING TMP.DAT ERROR OR IT DOES NOT EXIST'
    END SUBROUTINE

    SUBROUTINE BEDOUT(ISTA,JULTME4)
    INTEGER(4),INTENT(IN),OPTIONAL::ISTA
    REAL(4),INTENT(IN),OPTIONAL   ::JULTME4
    INTEGER(IK4) :: VER,HSIZE,I,NS,L,K,NX,NT,ITMP,N1,N,NSXD,OPT,IOS
    REAL(4)    :: TMP,EETIME,DTIME1,DTIME2
    REAL(8)    :: PTIME
    CHARACTER(200)::SS*5

    WRITE(SS,'(I5)') NLOC
    SFMT1='(F12.4,'//SS//'F12.5)'
    WRITE(SS,'(I5)') NLOC*NTOX
    SFMT2='(F12.4,'//SS//'F10.3)'
    WRITE(SS,'(I5)') NLOC*NSED
    SFMT3='(F12.4,'//SS//'F10.3)'
    WRITE(SS,'(I5)') NLOC*NSND
    SFMT4='(F12.4,'//SS//'F10.3)'
    N1 = KB
    OPT = 0

    IF (PRESENT(ISTA)) OPT = ISTA

    IF (OPT < 2) THEN
        ! *** OPEN FILES & WRITING THE HEADER LINES
        IF (OPT==0) OPEN(10,FILE=BEFILEI,STATUS='OLD',ACTION='READ',FORM='BINARY')
        OPEN(UOUTI(51),FILE='#output\RESULT\TOXB_TS_CEL.DAT',ACTION='WRITE')
        OPEN(UOUTI(52),FILE='#output\RESULT\SEDB_TS_CEL.DAT',ACTION='WRITE')
        OPEN(UOUTI(53),FILE='#output\RESULT\SNDB_TS_CEL.DAT',ACTION='WRITE')
        OPEN(UOUTI(54),FILE='#output\RESULT\HBED_TS_CEL.DAT',ACTION='WRITE')
        OPEN(UOUTI(55),FILE='#output\RESULT\BDENBED_TS_CEL.DAT',ACTION='WRITE')
        OPEN(UOUTI(56),FILE='#output\RESULT\PORBED_TS_CEL.DAT',ACTION='WRITE')
        OPEN(UOUTI(57),FILE='#output\RESULT\TOXB_TS_DOM.DAT',ACTION='WRITE')
        OPEN(UOUTI(58),FILE='#output\RESULT\SEDB_TS_DOM.DAT',ACTION='WRITE')
        OPEN(UOUTI(59),FILE='#output\RESULT\SNDB_TS_DOM.DAT',ACTION='WRITE')
        OPEN(UOUTI(60),FILE='#output\RESULT\HBED_TS_DOM.DAT',ACTION='WRITE')
        OPEN(UOUTI(61),FILE='#output\RESULT\BDENBED_TS_DOM.DAT',ACTION='WRITE')
        OPEN(UOUTI(62),FILE='#output\RESULT\PORBED_TS_DOM.DAT',ACTION='WRITE')

        DO I=51,56
            DO N=1,NLOC
                WRITE(UOUTI(I),'(A6,I5,2X,A3,I5)') '** I =',ICEL(N),'J =',JCEL(N)
            ENDDO
        ENDDO

        WRITE(UOUTI(51),'(A30,I3,A6)') '** TIME  (TOXB(LN,KB), N=1, ',NLOC,'),ug/l'
        WRITE(UOUTI(52),'(A30,I3,A6)') '** TIME  (SEDB(LN,KB), N=1, ',NLOC,'),mg/l'
        WRITE(UOUTI(53),'(A30,I3,A6)') '** TIME  (SNDB(LN,KB), N=1, ',NLOC,'),mg/l'
        WRITE(UOUTI(54),'(A30,I3,A3)') '** TIME  (HBED(LN,KB), N=1, ',NLOC,'),m'
        WRITE(UOUTI(55),'(A30,I3,A3)') '** TIME  (BDENBED(LN,KB), N=1, ',NLOC,')'
        WRITE(UOUTI(56),'(A30,I3,A3)') '** TIME  (PORBED(LN,KB), N=1, ',NLOC,')'
    ENDIF

    IF (OPT == 1) THEN
        ! *** WRITING THE HEADER BLOCK
        RETURN

    ELSEIF (OPT == 0) THEN
        ! *** READING BINARY DATA & EXPORTING
        SEDDIA = 0.
        HBED = 0.
        BDENBED= 0.
        PORBED = 0.
        SEDB = 0.
        SNDB = 0.
        TOXB = 0.

        READ(10,IOSTAT=IOS) VER,HSIZE
        IF (IOS /= 0) GOTO 999

        READ(10) ITMP,ITMP,ITMP
        READ(10) (ITMP,I=1,7)
        READ(10) NSED,NSND,NTOX

        WRITE(SS,'(I5)') NLOC
        SFMT1='(F12.4,'//SS//'F12.5)'
        WRITE(SS,'(I5)') NLOC*NTOX
        SFMT2='(F12.4,'//SS//'F10.3)'
        WRITE(SS,'(I5)') NLOC*NSED
        SFMT3='(F12.4,'//SS//'F10.3)'
        WRITE(SS,'(I5)') NLOC*NSND
        SFMT4='(F12.4,'//SS//'F10.3)'

        NSXD=NSED+NSND
        DTIME2 = 0

        DO NS=1,NSXD
            READ(10) SEDDIA(NS)
        ENDDO

        DO WHILE (1)
            READ(10, END=100) PTIME
            DTIME1 = DTIME2
            DTIME2 = PTIME
            EETIME = PTIME
            DO L=2,LA
                READ(10) KBT(L)
            ENDDO
            DO L=2,LA
                DO K=1,KB
                    READ(10) HBED(L,K),BDENBED(L,K),PORBED(L,K)
                    IF( ISTRAN(6) >= 1 )THEN
                        DO NS=1,NSED
                            READ(10) SEDB(L,K,NS)
                        ENDDO
                    ENDIF
                    IF( ISTRAN(7) >= 1 )THEN
                        DO NX=1,NSND
                            NS=NSED+NX
                            READ(10) SNDB(L,K,NX)
                        ENDDO
                    ENDIF
                    IF( ISTRAN(5) >= 1 )THEN
                        DO NT=1,NTOX
                            READ(10) TOXB(L,K,NT)
                        ENDDO
                    ENDIF
                ENDDO
            ENDDO

            ! ** EXPORT
            IF ((ALLSNAP==0 .AND. DTIME1<JULTIME .AND. JULTIME<=DTIME2) .OR. ALLSNAP==1) CALL DATEXPORT

        ENDDO
100     CLOSE(10,STATUS='KEEP')
        DO I=51,62
            CLOSE(UOUTI(I))
        ENDDO
        DEALLOCATE(SEDDIA)
        DEALLOCATE(KBT)
        DEALLOCATE(TOXB)
        DEALLOCATE(HBED)
        DEALLOCATE(BDENBED)
        DEALLOCATE(PORBED)
        DEALLOCATE(SEDB)
        DEALLOCATE(SNDB)

    ELSEIF (OPT ==2) THEN
        ! ** EXPORT DATA FROM EE_WC.OUT
        EETIME = JULTME4
        CALL DATEXPORT
    ENDIF
    RETURN
999 PRINT*, ' *** OPENING EE_BC.OUT ERROR OR IT DOES NOT EXIST'

    CONTAINS
    SUBROUTINE DATEXPORT
    IF (ISTRAN(5)==1) THEN
        ! ** TIME SERIES FOR CELLS
        WRITE(UOUTI(51),SFMT2) EETIME,((TOXB(LIJ(ICEL(N),JCEL(N)),N1,NT),N=1,NLOC),NT=1,NTOX)
        ! ** TIME SERIES FOR THE WHOLE DOMAIN
        WRITE(UOUTI(57),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER B = ',N1,'NTOX =',NTOX
        WRITE(UOUTI(57),'(A45,F15.5,A8)') '** DO NT=1,NTOX (TOXB(L,'//SSTR2//',NT),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
        DO NT=1,NTOX
            WRITE(UOUTI(57),'(A7,I5)') '** NT =',NT
            WRITE(UOUTI(57),'(10F10.5)') (TOXB(L,N1,NT),L=2,LA)
        ENDDO
    ENDIF

    IF(ISTRAN(6)==1) THEN
        WRITE(UOUTI(52),SFMT3) EETIME,((SEDB(LIJ(ICEL(N),JCEL(N)),N1,NS),N=1,NLOC),NS=1,NSED)
        WRITE(UOUTI(58),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER B = ',N1,'NSED =',NSED
        WRITE(UOUTI(58),'(A50,F15.5,A8)') '** DO NS=1,NSED (SEDB(L,'//SSTR2//',NS),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
        DO NS=1,NSED
            WRITE(UOUTI(58),'(A7,I5)') '** NS =',NS
            WRITE(UOUTI(58),'(10F10.5)') (SEDB(L,N1,NS),L=2,LA)
        ENDDO
    ENDIF

    IF(ISTRAN(7)==1) THEN
        WRITE(UOUTI(53),SFMT4) EETIME,((SNDB(LIJ(ICEL(N),JCEL(N)),N1,NX),N=1,NLOC),NX=1,NSND)
        WRITE(UOUTI(59),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER B = ',N1,'NSND =',NSND
        WRITE(UOUTI(59),'(A50,F15.5,A8)') '** DO NX=1,NSND (SNDB(L,'//SSTR2//',NX),L=2,LA) ENDDO;TIME=',EETIME,'  [MG/L]'
        DO NX=1,NSND
            WRITE(UOUTI(59),'(A7,I5)') '** NX =',NX
            WRITE(UOUTI(59),'(10F10.5)') (SNDB(L,N1,NX),L=2,LA)
        ENDDO
    ENDIF

    IF (ISTRAN(6)==1 .OR. ISTRAN(7)==1) THEN
        WRITE(UOUTI(54),SFMT1) EETIME,(HBED(LIJ(ICEL(N),JCEL(N)),N1),N=1,NLOC)
        WRITE(UOUTI(55),SFMT1) EETIME,(BDENBED(LIJ(ICEL(N),JCEL(N)),N1),N=1,NLOC)
        WRITE(UOUTI(56),SFMT1) EETIME,(PORBED(LIJ(ICEL(N),JCEL(N)),N1),N=1,NLOC)

        WRITE(UOUTI(60),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER B = ',N1,'NSND =',NSND
        WRITE(UOUTI(60),'(A40,F15.5,A8)') '** (HBED(L,KB),L=2,LA); TIME=',EETIME,'  [MG/L]'
        WRITE(UOUTI(60),'(10F10.5)') (HBED(L,N1),L=2,LA)

        WRITE(UOUTI(61),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER B = ',N1,'NSND =',NSND
        WRITE(UOUTI(61),'(A40,F15.5,A8)') '** (BDENBED(L,KB),L=2,LA); TIME=',EETIME,'  [MG/L]'
        WRITE(UOUTI(61),'(10F10.5)') (BDENBED(L,N1),L=2,LA)

        WRITE(UOUTI(62),'(A8,I5,2X,A9,I5,2X,A6,I5)') '** LA = ',LA,'LAYER B = ',N1,'NSND =',NSND
        WRITE(UOUTI(62),'(A40,F15.5,A8)') '** (PORBED(L,KB),L=2,LA); TIME=',EETIME,'  [MG/L]'
        WRITE(UOUTI(62),'(10F10.5)') (PORBED(L,N1),L=2,LA)
    ENDIF

    END SUBROUTINE

    END SUBROUTINE

    SUBROUTINE GETEE_TU

    INTEGER(4):: ITMP,JTMP,KTMP,K,L,N,VER,HSIZE,NACTIVE,NN,LO,IOS
    REAL(8)   :: PTIME,DZ
    REAL(4),ALLOCATABLE:: KTUR(:,:),VTUR(:,:),TURNU(:,:),Z(:,:)

    OPEN(95,FILE=TUFILEI,FORM='BINARY')
    OPEN(1,FILE='#output\RESULT\VPROF_TUR.DAT', ACTION='WRITE')

    READ(95,IOSTAT = IOS) VER,HSIZE
    IF (IOS /=0) GOTO 999
    READ(95) ITMP,JTMP,KTMP,NACTIVE
    ALLOCATE(KTUR(LA,KC),VTUR(LA,KC),TURNU(LA,KC),Z(LA,KC))

    DO WHILE(1)
        READ (95,END=998) NN,PTIME
        DO L=2,LA
            READ(95) (Z(L,K),KTUR(L,K),VTUR(L,K),TURNU(L,K),K=1,KC)
        ENDDO
        WRITE(1,'(A1,I8,F12.2)') '*',NN,PTIME
        DO K=1,KC
            WRITE(1,'(F12.4,30F12.6)') (Z(LIJ(ICEL(N),JCEL(N)),K),KTUR(LIJ(ICEL(N),JCEL(N)),K),VTUR(LIJ(ICEL(N),JCEL(N)),K),TURNU(LIJ(ICEL(N),JCEL(N)),K),N=1,NLOC)
        ENDDO
    ENDDO
998 WRITE(*,'(A)')' *** FINISHED READING TUR'
    CLOSE(95)
    CLOSE(1)
    RETURN
999 PRINT*, ' **** OPENING EE_TUR.OUT ERROR OR IT DOES NOT EXIST'
    END SUBROUTINE

    END MODULE
