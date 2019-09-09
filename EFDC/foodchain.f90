SUBROUTINE FOODCHAIN(IFINISH)

  ! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
  !
  ! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001

  !----------------------------------------------------------------------C
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY                 DATE APPROVED    BY

  !----------------------------------------------------------------------C
  !
  ! **  SUBROUTINES OUTPUT SPACE AND TIME AVERAGE TOXICS CONCENTRATIONS
  ! **  FOR FOOD CHAIN MODEL

  !**********************************************************************C

  USE GLOBAL
  IMPLICIT NONE

  INTEGER :: IFINISH,JSFDCH,M,NT,L,K,NS,NX,KSTOP,KTMP
  REAL    :: TIMFDCH,HBEDTMP,HBSTOP,TMPVAL,PORHINV,FDCHVAL
  
  INTEGER,ALLOCATABLE,DIMENSION(:) :: KBFC

  REAL,ALLOCATABLE,DIMENSION(:,:) :: VALPOCW
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPVOLW

  REAL,ALLOCATABLE,DIMENSION(:,:) :: WTBED
  REAL,ALLOCATABLE,DIMENSION(:,:) :: VALPOCB
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPVOLB

  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPTXWF
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPTXWC
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPTXWP
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPTXBF
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPTXBC
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPTXBP
  REAL,ALLOCATABLE,DIMENSION(:,:) :: TMPTXBPD

  REAL,ALLOCATABLE,DIMENSION(:,:) :: VALBCONC

  REAL,ALLOCATABLE,DIMENSION(:) :: TMPDOCW
  REAL,ALLOCATABLE,DIMENSION(:) :: TMPPOCW
  REAL,ALLOCATABLE,DIMENSION(:) :: TMPDOCB
  REAL,ALLOCATABLE,DIMENSION(:) :: TMPPOCB
  REAL,ALLOCATABLE,DIMENSION(:) :: VOLFCW
  REAL,ALLOCATABLE,DIMENSION(:) :: VOLFCB

  LOGICAL,ALLOCATABLE,DIMENSION(:) :: LMASKFC

  REAL,ALLOCATABLE,DIMENSION(:) :: FDCHDOCW
  REAL,ALLOCATABLE,DIMENSION(:) :: FDCHPOCW
  REAL,ALLOCATABLE,DIMENSION(:) :: FDCHDOCB
  REAL,ALLOCATABLE,DIMENSION(:) :: FDCHPOCB

  REAL,ALLOCATABLE,DIMENSION(:,:) :: FDCHTXWF
  REAL,ALLOCATABLE,DIMENSION(:,:) :: FDCHTXWC
  REAL,ALLOCATABLE,DIMENSION(:,:) :: FDCHTXWP
  REAL,ALLOCATABLE,DIMENSION(:,:) :: FDCHTXBF
  REAL,ALLOCATABLE,DIMENSION(:,:) :: FDCHTXBC
  REAL,ALLOCATABLE,DIMENSION(:,:) :: FDCHTXBP
  REAL,ALLOCATABLE,DIMENSION(:,:) :: FDCHTXBD

  IF(  .NOT. ALLOCATED(KBFC) )THEN
    ALLOCATE(KBFC(LCM))
    ALLOCATE(VALPOCW(LCM,KCM))
    ALLOCATE(TMPVOLW(LCM,KCM))

    ALLOCATE(WTBED(LCM,KBM))
    ALLOCATE(VALPOCB(LCM,KBM))
    ALLOCATE(VALBCONC(LCM,KBM))

    ALLOCATE(TMPVOLB(LCM,KBM))
    ALLOCATE(TMPTXWF(NFDCHZ,NTXM))
    ALLOCATE(TMPTXWC(NFDCHZ,NTXM))
    ALLOCATE(TMPTXWP(NFDCHZ,NTXM))
    ALLOCATE(TMPTXBF(NFDCHZ,NTXM))
    ALLOCATE(TMPTXBC(NFDCHZ,NTXM))
    ALLOCATE(TMPTXBP(NFDCHZ,NTXM))
    ALLOCATE(TMPTXBPD(NFDCHZ,NTXM))

    ALLOCATE(TMPDOCW(NFDCHZ))
    ALLOCATE(TMPPOCW(NFDCHZ))
    ALLOCATE(TMPDOCB(NFDCHZ))
    ALLOCATE(TMPPOCB(NFDCHZ))
    ALLOCATE(VOLFCW(NFDCHZ))
    ALLOCATE(VOLFCB(NFDCHZ))

    ALLOCATE(FDCHDOCW(NFDCHZ))
    ALLOCATE(FDCHPOCW(NFDCHZ))
    ALLOCATE(FDCHDOCB(NFDCHZ))
    ALLOCATE(FDCHPOCB(NFDCHZ))

    ALLOCATE(FDCHTXWF(NFDCHZ,NTXM))
    ALLOCATE(FDCHTXWC(NFDCHZ,NTXM))
    ALLOCATE(FDCHTXWP(NFDCHZ,NTXM))
    ALLOCATE(FDCHTXBF(NFDCHZ,NTXM))
    ALLOCATE(FDCHTXBC(NFDCHZ,NTXM))
    ALLOCATE(FDCHTXBP(NFDCHZ,NTXM))
    ALLOCATE(FDCHTXBD(NFDCHZ,NTXM))

    ! *** ALLOCATE LOCAL ARRAYS
    KBFC=0
    VALPOCW=0.
    TMPVOLW=0.

    WTBED=0.
    VALPOCB=0.
    VALBCONC=0.

    TMPVOLB=0.
    TMPTXWF=0.
    TMPTXWC=0.
    TMPTXWP=0.
    TMPTXBF=0.
    TMPTXBC=0.
    TMPTXBP=0.
    TMPTXBPD=0.

    TMPDOCW=0.
    TMPPOCW=0.
    TMPDOCB=0.
    TMPPOCB=0.
    VOLFCW=0.
    VOLFCB=0.

    FDCHDOCW=0.
    FDCHPOCW=0.
    FDCHDOCB=0.
    FDCHPOCB=0.

    FDCHTXWF=0.
    FDCHTXWC=0.
    FDCHTXWP=0.
    FDCHTXBF=0.
    FDCHTXBC=0.
    FDCHTXBP=0.
    FDCHTXBD=0.
  ENDIF

  !**********************************************************************C

  IF( IFINISH == 1 ) GO TO 2000
  IF( JSFDCH == 0 )GO TO 1000

  !      WRITE(8,*)' FIRST ENTRY TO FOODCHAIN.FOR '

  IF( DEBUG )THEN
    OPEN(1,FILE=OUTDIR//'FOODCHAIN.OUT')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'FOODCHAIN.OUT')
    WRITE(1,121)
    WRITE(1,122)
    WRITE(1,123)
    CLOSE(1)
  ENDIF

  !     JSFDCH=0

  DO M=1,NFDCHZ
    FDCHDOCW(M)=0.
    FDCHPOCW(M)=0.
    FDCHDOCB(M)=0.
    FDCHPOCB(M)=0.
  ENDDO

  DO NT=1,NTOX
    DO M=1,NFDCHZ
      FDCHTXWF(M,NT)=0.
      FDCHTXWC(M,NT)=0.
      FDCHTXWP(M,NT)=0.
      FDCHTXBF(M,NT)=0.
      FDCHTXBC(M,NT)=0.
      FDCHTXBP(M,NT)=0.
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      FDCHTXBD(M,NT)=0.
    !####################################################################################
    ENDDO
  ENDDO

  TIMFDCH=0.0

  !**********************************************************************C

    1000 CONTINUE

  TIMFDCH = TIMFDCH + DTSED

  ! **  INITIALIZE VOLUMES AND VOLUME AVERAGES

  DO M=1,NFDCHZ
    VOLFCW(M)=0.
    VOLFCB(M)=0.
  ENDDO

  DO M=1,NFDCHZ
    TMPDOCW(M)=0.
    TMPPOCW(M)=0.
    TMPDOCB(M)=0.
    TMPPOCB(M)=0.
  ENDDO

  DO NT=1,NTOX
    DO M=1,NFDCHZ
      TMPTXWF(M,NT)=0.
      TMPTXWC(M,NT)=0.
      TMPTXWP(M,NT)=0.
      TMPTXBF(M,NT)=0.
      TMPTXBC(M,NT)=0.
      TMPTXBP(M,NT)=0.
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      TMPTXBPD(M,NT)=0.
    !####################################################################################
    ENDDO
  ENDDO

  ! **  INITIALIZE MASK
  DO L=2,LA
    LMASKFC(L)=.FALSE.
  ENDDO
  !
  DO L=2,LA
    IF( LMASKDRY(L) )THEN
      IF( MFDCHZ(L) > 0 )LMASKFC(L)=.TRUE.
    ENDIF
  ENDDO

  !----------------------------------------------------------------------C
  !
  ! **  VOLUME WEIGHTED AVERAGE OVER WATER COLUMN ZONES
  !
  !     STDOCW(L,K) HAS UNITS: MG/L OR GM/M**3
  !     STPOCW(L,K) AND VALPOCW(L,K) HAVE UNITS: MG/L OR GM/M**3
  IF( ISTPOCW <= 1 )THEN
    DO K=1,KC
      DO L=2,LA
        IF( LMASKFC(L)) VALPOCW(L,K)=STPOCW(L,K)
      ENDDO
    ENDDO
  ENDIF

  IF( ISTPOCW >= 2 )THEN
    DO K=1,KC
      DO L=2,LA
        IF( LMASKFC(L)) VALPOCW(L,K)=0.
      ENDDO
    ENDDO
    DO NS=1,NSED
      DO K=1,KC
        DO L=2,LA
          IF( LMASKFC(L)) VALPOCW(L,K)= &
            VALPOCW(L,K)+SED(L,K,NS)*STFPOCW(L,K,NS)
        ENDDO
      ENDDO
    ENDDO
    DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KC
        DO L=2,LA
          IF( LMASKFC(L)) VALPOCW(L,K)= &
            VALPOCW(L,K)+SND(L,K,NX)*STFPOCW(L,K,NS)
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  !     changed to areal weighting from voloume weighting
  DO K=1,KC
    DO L=2,LA
    !        IF( LMASKFC(L)) TMPVOLW(L,K)=DXYP(L)*HP(L)*DZC(L,K)
      IF( LMASKFC(L)) TMPVOLW(L,K)=DXYP(L)*DZC(L,K)
    ENDDO
  ENDDO

  DO K=1,KC
    DO L=2,LA
      IF( LMASKFC(L) )THEN
        M=MFDCHZ(L)
        VOLFCW(M)=VOLFCW(M)+TMPVOLW(L,K)
      ENDIF
    ENDDO
  ENDDO

  !     TMPTXWF(M,NT)/TMPVOLW HAS UNITS: UG/L OR MG/M**3
  !     TMPTXWC(M,NT)/TMPVOLW HAS UNITS: UG/L OR MG/M**3
  !     TMPTXWP(M,NT)/TMPVOLW HAS UNITS: UG/MG
  !     TMPDOCW(M,NT)/TMPVOLW AND STDOCW(L,K) HAVE UNITS: MG/L OR GM/M**3
  !     TMPPOCW(M,NT)/TMPVOLW AND VALPOCW(L,K) HAVE UNITS: MG/L OR GM/M**3

  DO K=1,KC
  DO L=2,LA
    IF( LMASKFC(L) )THEN
    M=MFDCHZ(L)
    TMPDOCW(M)=TMPDOCW(M)+TMPVOLW(L,K)*STDOCW(L,K)
    TMPPOCW(M)=TMPPOCW(M)+TMPVOLW(L,K)*VALPOCW(L,K)
    ENDIF
  ENDDO
  ENDDO

  DO NT=1,NTOX
  DO K=1,KC
  DO L=2,LA
    IF( LMASKFC(L) )THEN
      M=MFDCHZ(L)
      TMPTXWF(M,NT)=TMPTXWF(M,NT) + TMPVOLW(L,K)*TOXFDFW(L,K,NT)*TOX(L,K,NT)
      TMPTXWC(M,NT)=TMPTXWC(M,NT) + TMPVOLW(L,K)*TOXCDFW(L,K,NT)*TOX(L,K,NT)
      IF( VALPOCW(L,K) > 0.) TMPTXWP(M,NT) = TMPTXWP(M,NT) + TMPVOLW(L,K)*TOXPFTW(L,K,NT)*TOX(L,K,NT) / VALPOCW(L,K)
    ENDIF
  ENDDO
  ENDDO
  ENDDO

  DO M=1,NFDCHZ
    IF( VOLFCW(M) > 0.0 )THEN
      TMPDOCW(M)=TMPDOCW(M)/VOLFCW(M)
      TMPPOCW(M)=TMPPOCW(M)/VOLFCW(M)
    ENDIF
  ENDDO


  DO NT=1,NTOX
    DO M=1,NFDCHZ
      IF( VOLFCW(M) > 0.0 )THEN
        TMPTXWF(M,NT)=TMPTXWF(M,NT)/VOLFCW(M)
        TMPTXWC(M,NT)=TMPTXWC(M,NT)/VOLFCW(M)
        TMPTXWP(M,NT)=TMPTXWP(M,NT)/VOLFCW(M)
      ENDIF
    ENDDO
  ENDDO

  !     CONVERT PARTICULATE FROM UG/MG TO UG/GM

  DO NT=1,NTOX
    DO M=1,NFDCHZ
      TMPTXWP(M,NT)=1000.*TMPTXWP(M,NT)
    ENDDO
  ENDDO

  !----------------------------------------------------------------------C
  !
  ! **  VOLUME WEIGHTED AVERAGE OVER BED ZONES
  !
  !     STDOCB(L,K) HAS UNITS: MG/L OR GM/M**3 (MASS PER VOLUME OF PORE WATER)
  !     STPOCB(L,K) AND VALPOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER TOTAL VOLUME)

  IF( ISTPOCB <= 1 )THEN
    DO K=1,KB
      DO L=2,LA
        IF( LMASKFC(L)) VALPOCB(L,K)=STPOCB(L,K)
      ENDDO
    ENDDO
  ENDIF

  !      IF( ISTPOCB >= 2 )THEN   ! PMC  FPOCB IS NOT INIITALIZED UNTIL ISTPOCB>3
  IF( ISTPOCB >= 4 )THEN
    DO K=1,KB
      DO L=2,LA
        VALPOCB(L,K)=0.
      ENDDO
    ENDDO
    DO NS=1,NSED
      DO K=1,KB
        DO L=2,LA
          IF( LMASKFC(L) )THEN
            IF( K <= KBT(L)) VALPOCB(L,K)=VALPOCB(L,K)+SEDB(L,K,NS)*FPOCB(L,K)/HBED(L,K)
            !####################################################################################
            ! RM 05/14/04
            ! Change to average using data-based foc rather than partitioning foc
            !     &                    +SEDB(L,K,NS)*STFPOCB(L,K,NS)/HBED(L,K)
            !####################################################################################
          ENDIF
        ENDDO
      ENDDO
    ENDDO
    DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KB
        DO L=2,LA
          IF( LMASKFC(L) )THEN
            IF( K <= KBT(L)) VALPOCB(L,K)=VALPOCB(L,K)+SNDB(L,K,NX)*FPOCB(L,K)/HBED(L,K)
            !####################################################################################
            ! RM 05/14/04
            ! Change to average using data-based foc rather than partitioning foc
            !     &                    +SNDB(L,K,NX)*STFPOCB(L,K,NS)/HBED(L,K)
            !####################################################################################
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  !####################################################################################
  ! RM 05/14/04
  ! Change to average dry weight PCB
  DO K=1,KB
    DO L=2,LA
        VALBCONC(L,K)=0.
    ENDDO
  ENDDO
  DO NS=1,NSED
    DO K=1,KB
      DO L=2,LA
        IF( LMASKFC(L) )THEN
          IF( K <= KBT(L)) VALBCONC(L,K)=VALBCONC(L,K) + SEDB(L,K,NS)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  DO NX=1,NSND
    NS=NX+NSED
    DO K=1,KB
      DO L=2,LA
        IF( LMASKFC(L) )THEN
          IF( K <= KBT(L)) VALBCONC(L,K)=VALBCONC(L,K) + SNDB(L,K,NX)
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  !####################################################################################

  DO K=1,KB
    DO L=2,LA
      WTBED(L,K)=0.0
    ENDDO
  ENDDO

  DO L=2,LA
    IF( LMASKFC(L) )THEN
      KBFC(L)=KBT(L)
      HBEDTMP=0.0
      KSTOP=0
      DO K=KBT(L),1,-1
        HBEDTMP=HBEDTMP+HBED(L,K)
        IF( HBEDTMP > HBFDCH .AND. KSTOP == 0 )THEN
          KBFC(L)=K
          KSTOP=1
          HBSTOP=HBED(L,K)-HBEDTMP+HBFDCH
          WTBED(L,K)=HBSTOP/HBFDCH
        ENDIF
      ENDDO
      KTMP=KBFC(L)+1
      DO K=KTMP,KBT(L)
        !####################################################################################
        ! RM 05/14/04
        ! Weightages greater than 1 could occur with this method of depth-weighting.
        ! When the thickness of the top layer is greater than HBFDCH (0.1524
        ! meters), the weightage assigned to this layer could become greater
        ! than 1. Need to confirm this with JH.
        !####################################################################################
        WTBED(L,K)=HBED(L,K)/HBFDCH
      ENDDO
    ELSE
      KBFC(L)=0
    ENDIF
  ENDDO
  !
  IF( JSFDCH == 1 .AND. DEBUG )THEN
    OPEN(1,FILE=OUTDIR//'FOODCHAIN.DIA')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'FOODCHAIN.DIA')
    DO L=2,LA
      IF( LMASKFC(L) )THEN
        WRITE(1,111)IL(L),JL(L),KBFC(L),KBT(L),(WTBED(L,K),K=KBFC(L),KBT(L))
        WRITE(1,112)(HBED(L,K),K=KBFC(L),KBT(L))
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
  
  DO K=1,KB
  DO L=2,LA
    IF( LMASKFC(L) )THEN
      IF( K >= KBFC(L) .AND. K <= KBT(L) )THEN
        TMPVOLB(L,K)=DXYP(L)*WTBED(L,K)*HBED(L,K)
      ENDIF
    ENDIF
  ENDDO
  ENDDO

  DO K=1,KB
  DO L=2,LA
    IF( LMASKFC(L) )THEN
      M=MFDCHZ(L)
      IF( K >= KBFC(L) .AND. K <= KBT(L) )THEN
        VOLFCB(M)=VOLFCB(M)+TMPVOLB(L,K)
      ENDIF
    ENDIF
  ENDDO
  ENDDO

  !     TMPTXBF(M,NT)/TMPVOLB HAS UNITS: UG/L OR MG/M**3  (MASS PER VOLUME PORE WATER)
  !     TMPTXBC(M,NT)/TMPVOLB HAS UNITS: UG/L OR MG/M**3  (MASS PER VOLUME PORE WATER)
  !     TMPTXWP(M,NT)/TMPVOLB HAS UNITS: UG/MG
  !     TMPDOCW(M,NT)/TMPVOLB AND STDOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER VOLUME PORE WATER)
  !     TMPPOCW(M,NT)/TMPVOLB AND VALPOCB(L,K) HAVE UNITS: MG/L OR GM/M**3 (MASS PER TOTAL VOLUME)

  DO K=KB,1,-1
  DO L=2,LA
    IF( LMASKFC(L) )THEN
      M=MFDCHZ(L)
      IF( K >= KBFC(L) .AND. K <= KBT(L) )THEN
        TMPDOCB(M)=TMPDOCB(M)+TMPVOLB(L,K)*STDOCB(L,K)
        TMPPOCB(M)=TMPPOCB(M)+TMPVOLB(L,K)*VALPOCB(L,K)
      ENDIF
    ENDIF
  ENDDO
  ENDDO

  DO NT=1,NTOX
  DO K=KB,1,-1
  DO L=2,LA
    IF( LMASKFC(L) )THEN
      M=MFDCHZ(L)
      IF( K >= KBFC(L) .AND. K <= KBT(L) )THEN
        TMPVAL=HBED(L,K)*VALPOCB(L,K)
        PORHINV=1.0/(HBED(L,K)*PORBED(L,K))
        TMPTXBF(M,NT)=TMPTXBF(M,NT)+TMPVOLB(L,K)*PORHINV*TOXFDFB(L,K,NT)*TOXB(L,K,NT)
        TMPTXBC(M,NT)=TMPTXBC(M,NT)+TMPVOLB(L,K)*PORHINV*TOXCDFB(L,K,NT)*TOXB(L,K,NT)
        
        IF( TMPVAL > 0.) TMPTXBP(M,NT)=TMPTXBP(M,NT)+TMPVOLB(L,K)*TOXPFTB(L,K,NT)*TOXB(L,K,NT)/TMPVAL
        !####################################################################################
        ! RM 05/14/04
        ! Change to average dry weight PCBs
        TMPTXBPD(M,NT)=TMPTXBPD(M,NT) + TMPVOLB(L,K)*TOXPFTB(L,K,NT)*TOXB(L,K,NT)/VALBCONC(L,K)
        !####################################################################################
      ENDIF
    ENDIF
  ENDDO
  ENDDO
  ENDDO
  
  DO M=1,NFDCHZ
    IF( VOLFCB(M) > 0.0 )THEN
      TMPDOCB(M)=TMPDOCB(M)/VOLFCB(M)
      TMPPOCB(M)=TMPPOCB(M)/VOLFCB(M)
    ENDIF
  ENDDO
  !
  DO NT=1,NTOX
  DO M=1,NFDCHZ
    IF( VOLFCB(M) > 0.0 )THEN
      TMPTXBF(M,NT)=TMPTXBF(M,NT)/VOLFCB(M)
      TMPTXBC(M,NT)=TMPTXBC(M,NT)/VOLFCB(M)
      TMPTXBP(M,NT)=TMPTXBP(M,NT)/VOLFCB(M)
      !####################################################################################
      ! RM 05/14/04
      ! Change to average dry weight PCBs
          TMPTXBPD(M,NT)=TMPTXBPD(M,NT)/VOLFCB(M)
      !####################################################################################
    ENDIF
  ENDDO
  ENDDO

  !     CONVERT PARTICULATE FROM UG/MG TO UG/GM
  DO NT=1,NTOX
    DO M=1,NFDCHZ
      TMPTXBP(M,NT)=1000.*TMPTXBP(M,NT)
    ENDDO
  ENDDO

  !----------------------------------------------------------------------C
  !
  ! **  ACCUMULATE THE TIME AVERAGE
  DO M=1,NFDCHZ
    FDCHDOCW(M)=FDCHDOCW(M)+DTSED*TMPDOCW(M)
    FDCHPOCW(M)=FDCHPOCW(M)+DTSED*TMPPOCW(M)
    FDCHDOCB(M)=FDCHDOCB(M)+DTSED*TMPDOCB(M)
    FDCHPOCB(M)=FDCHPOCB(M)+DTSED*TMPPOCB(M)
  ENDDO
  !
  DO NT=1,NTOX
    DO M=1,NFDCHZ
      FDCHTXWF(M,NT)=FDCHTXWF(M,NT)+DTSED*TMPTXWF(M,NT)
      FDCHTXWC(M,NT)=FDCHTXWC(M,NT)+DTSED*TMPTXWC(M,NT)
      FDCHTXWP(M,NT)=FDCHTXWP(M,NT)+DTSED*TMPTXWP(M,NT)
      FDCHTXBF(M,NT)=FDCHTXBF(M,NT)+DTSED*TMPTXBF(M,NT)
      FDCHTXBC(M,NT)=FDCHTXBC(M,NT)+DTSED*TMPTXBC(M,NT)
      FDCHTXBP(M,NT)=FDCHTXBP(M,NT)+DTSED*TMPTXBP(M,NT)
      !####################################################################################
      ! RM 05/14/04
      ! Change to average dry weight PCBs
      FDCHTXBD(M,NT)=FDCHTXBD(M,NT)+DTSED*TMPTXBPD(M,NT)
      !####################################################################################
      ENDDO
  ENDDO

  JSFDCH=0

  IF( TIMFDCH < TFCAVG) RETURN

  !**********************************************************************C
  !
  ! **  COMPLETE AVERAGING AND OUTPUT RESULTS
  !
    2000 CONTINUE

  FDCHVAL=1./TIMFDCH
  DO M=1,NFDCHZ
    FDCHDOCW(M)=FDCHVAL*FDCHDOCW(M)
    FDCHPOCW(M)=FDCHVAL*FDCHPOCW(M)
    FDCHDOCB(M)=FDCHVAL*FDCHDOCB(M)
    FDCHPOCB(M)=FDCHVAL*FDCHPOCB(M)
  ENDDO
  !
  DO NT=1,NTOX
    DO M=1,NFDCHZ
      FDCHTXWF(M,NT)=FDCHVAL*FDCHTXWF(M,NT)
      FDCHTXWC(M,NT)=FDCHVAL*FDCHTXWC(M,NT)
      FDCHTXWP(M,NT)=FDCHVAL*FDCHTXWP(M,NT)
      FDCHTXBF(M,NT)=FDCHVAL*FDCHTXBF(M,NT)
      FDCHTXBC(M,NT)=FDCHVAL*FDCHTXBC(M,NT)
      FDCHTXBP(M,NT)=FDCHVAL*FDCHTXBP(M,NT)
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      FDCHTXBD(M,NT)=FDCHVAL*FDCHTXBD(M,NT)
    !####################################################################################
    ENDDO
  ENDDO

  IF( DEBUG )THEN
    OPEN(1,FILE=OUTDIR//'FOODCHAIN.OUT',POSITION='APPEND')
    WRITE(1,101)TIMEDAY,NTOX,NFDCHZ,TIMFDCH
    DO NT=1,NTOX
      DO M=1,NFDCHZ
        WRITE(1,102)NT,M, &
                      FDCHTXWF(M,NT),FDCHTXWC(M,NT),FDCHTXWP(M,NT), &
                      FDCHDOCW(M),FDCHPOCW(M),FDCHTXBF(M,NT), &
                      FDCHTXBC(M,NT),FDCHTXBP(M,NT),FDCHDOCB(M), &
                      FDCHPOCB(M),FDCHTXBD(M,NT)
      ENDDO
    ENDDO
    CLOSE(1)
  ENDIF

  !**********************************************************************C
  !
  ! **  INITIALIZE FOR NEXT AVERAGING PERIOD
  DO M=1,NFDCHZ
    FDCHDOCW(M)=0.
    FDCHPOCW(M)=0.
    FDCHDOCB(M)=0.
    FDCHPOCB(M)=0.
  ENDDO
  !
  DO NT=1,NTOX
    DO M=1,NFDCHZ
      FDCHTXWF(M,NT)=0.
      FDCHTXWC(M,NT)=0.
      FDCHTXWP(M,NT)=0.
      FDCHTXBF(M,NT)=0.
      FDCHTXBC(M,NT)=0.
      FDCHTXBP(M,NT)=0.
    !####################################################################################
    ! RM 05/14/04
    ! Change to average dry weight PCBs
      FDCHTXBD(M,NT)=0.
    !####################################################################################
    ENDDO
  ENDDO

  TIMFDCH=0.0

  !**********************************************************************C
    111 FORMAT(4I5,10F10.4)
    112 FORMAT(20X,10F10.4)
    101 FORMAT(F12.4,2I7,F12.3)
    102 FORMAT(1X,2I6,10E13.5)
    103 FORMAT('              TXWF         TXWC         TXWP', &
          '         DOCW         POCW         TXBF         TXBC', &
          '         TXBP (roc)   DOCB         POCB        TXBPD (r)')
    121 FORMAT('DATA: OUTPUT TIME (DAYS), NTOX, NZONES, ', &
          'AERAGING PERIOD (SECS)')
    122 FORMAT('DATA: NT    NZ   TXWF         TXWC         TXWP', &
          '         DOCW         POCW         TXBF         TXBC', &
          '         TXBP (roc)   DOCB         POCB        TXBPD (r)')
    123 FORMAT('DATA:            UG/L         UG/L         UG/GM', &
          '        MG/L         MG/L         UG/L         UG/L', &
          '         UG/GM OC     MG/L         MG/L        UG/GM Dry')
  !
  !**********************************************************************C
  !
  RETURN

END
