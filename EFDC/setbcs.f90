SUBROUTINE SETBCS

  ! CHANGE RECORD
  !  MODIFIED BOUNDARY CONDITION FLAGS FOR TYPE 2 OPEN BOUNDARIES
  !  ADDED REAL FLAGS RSSBCE(L),RSSBCW(L),RSSBCN(L),RSSBCS(L)
  !  TO MODIFIED CALCULATION OF CELL CENTER BED STRESS (STORED AS QQ(L,0))
  !  AND THE OUTPUTED CELL CENTER VELOCITY FOR CELLS HAVE SOURCE/SINKS
  ! **  SUBROUTINE SETBCS SETS BOUNDARY CONDITION SWITCHES

  USE GLOBAL

  IMPLICIT NONE

  INTEGER :: L,I,J,NPN,LE,LW,LS,LN,LL,NCTL,IU,JU,IQ,L1,ID,JD,LTMP,NWR,NJP,NMD,NDRYTMP
  REAL    :: DDYDDDX,DDXDDDY,RQDW 
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SUBEW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SVBNS

  IF(  .NOT. ALLOCATED(SUBEW) )THEN
    ALLOCATE(SUBEW(LCM))
    ALLOCATE(SVBNS(LCM))
    SUBEW=0.0
    SVBNS=0.0
  ENDIF

  ! **  SET LAND-WATER BOUNDARY SWITCHES
  ITRICELL=0

  DO L=2,LA
    I=IL(L)
    J=JL(L)
    IF( LCT(L) == 1 )THEN
      STCUV(L)=0.
      ITRICELL=1
      STCAP(L)=0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 2 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 3 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 4 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 5 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 9 ) SUB(L)=0.
      IF( IJCT(I,J-1) == 1 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 2 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 3 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 4 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 5 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 9 ) SVB(L)=0.
    ENDIF
    IF( LCT(L) == 2 )THEN
      STCUV(L)=0.
      ITRICELL=1
      STCAP(L)=0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 2 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 3 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 4 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 5 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 9 ) SUB(L)=0.
      IF( IJCT(I,J-1) == 1 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 2 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 3 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 4 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 5 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 9 ) SVB(L)=0.
    ENDIF
    IF( LCT(L) == 3 )THEN
      STCUV(L)=0.
      ITRICELL=1
      STCAP(L)=0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 2 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 3 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 4 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 9 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 5 ) SUB(L)=0.
      IF( IJCT(I,J-1) == 1 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 2 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 3 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 4 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 5 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 9 ) SVB(L)=0.
    ENDIF
    IF( LCT(L) == 4 )THEN
      STCUV(L)=0.
      ITRICELL=1
      STCAP(L)=0.5
      IF( IJCT(I-1,J) == 1 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 2 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 3 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 4 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 5 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 9 ) SUB(L)=0.
      IF( IJCT(I,J-1) == 1 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 2 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 3 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 4 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 5 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 9 ) SVB(L)=0.
    ENDIF
    IF( LCT(L) == 5 )THEN
      IF( IJCT(I-1,J) == 1 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 2 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 3 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 4 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 5 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 9 ) SUB(L)=0.
      IF( IJCT(I,J-1) == 1 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 2 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 3 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 4 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 5 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 9 ) SVB(L)=0.
    ENDIF
    IF( LCT(L) == 6 )THEN
      IF( IJCT(I-1,J) == 1 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 2 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 3 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 4 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 5 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 6 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 7 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 9 ) SUB(L)=0.
      IF( IJCT(I,J-1) == 1 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 2 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 3 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 4 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 5 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 6 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 7 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 9 ) SVB(L)=0.
    ENDIF
    IF( LCT(L) == 7 )THEN
      IF( IJCT(I-1,J) == 1 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 2 ) SUB(L)=0.
      IF( IJCT(I-1,J) == 3 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 4 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 5 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 6 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 7 ) SUB(L)=1.
      IF( IJCT(I-1,J) == 9 ) SUB(L)=0.
      IF( IJCT(I,J-1) == 1 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 2 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 3 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 4 ) SVB(L)=0.
      IF( IJCT(I,J-1) == 5 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 6 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 7 ) SVB(L)=1.
      IF( IJCT(I,J-1) == 9 ) SVB(L)=0.
    ENDIF
  ENDDO
  SUB(1)=0.
  SVB(1)=0.
  SUB(LC)=0.
  SVB(LC)=0.

  ! **  MODIFY LAND-WATER BNDRY CONDS FOR PERIOD GRID IN E-W DIRECTION
  IF( ISCONNECT >= 2 )THEN
    DO NPN=1,NPEWBP
      L=LIJ(IWPEW(NPN),JWPEW(NPN))
      SUB(L)=1.
      SUBO(L)=1.
    ENDDO
  ENDIF

  ! **  MODIFY LAND-WATER BNDRY CONDS FOR PERIOD GRID IN N-S DIRECTION
    IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    DO NPN=1,NPNSBP
      LS=LIJ(ISPNS(NPN),JSPNS(NPN))
      SVB(LS)=1.
      SVBO(LS)=1.
    ENDDO
  ENDIF

  ! **  SET WATER-WATER (P OR SURFACE ELEVATION) BOUNDARY SWITCHES
  DO LL=1,NPBW
    I=IPBW(LL)
    J=JPBW(LL)
    L=LIJ(I,J)
    LPBW(LL)=L
    SPB(L)=0.           ! *** Used for On/Off of various BC forcings at !open boundaries
    SUB(L)=0.
    SVB(L)=0.
    SVB(LNC(L))=0.
    SWB(L)=0.               ! *** Used for On/Off of Vertical Velocities
    SAAX(L)=0.              ! *** Used for On/Off of Horizontal Momentum Stresses (X Dir)
    IF( ISPBW(LL) <= 1 )THEN
      SWB(LEC(L))=0.
      SVB(LEC(L))=0.
      SVB(LNC(LEC(L)))=0.
      SCAX(LEC(L))=0.       ! *** Used for On/Off of Coriolis & Curvature Stresses
    END IF
    SAAX(LEC(L))=0.         ! *** Used for On/Off of Horizontal Momentum Stresses (X Dir)
    SDX(LEC(L))=0.          ! *** Used for On/Off of Horizontal Momentum Diffusion Stresses
  ENDDO
  DO LL=1,NPBE
    I=IPBE(LL)
    J=JPBE(LL)
    L=LIJ(I,J)
    LPBE(LL)=L
    SPB(L)=0.
    SVB(L)=0.
    SVB(LNC(L))=0.
    SWB(L)=0.
    IF( ISPBE(LL) <= 1 )THEN
      SWB(LWC(L))=0.
      SVB(LWC(L))=0.
      SVB(LNC(LWC(L)))=0.
    END IF
    SCAX(L)=0.
    SAAX(L)=0.
    SDX(L)=0.
  ENDDO
  DO LL=1,NPBS
    I=IPBS(LL)
    J=JPBS(LL)
    L=LIJ(I,J)
    LPBS(LL)=L
    LN=LNC(L)
    SPB(L)=0.
    SVB(L)=0.
    SUB(L)=0.
    SUB(LEC(L))=0.
    SWB(L)=0.
    IF( ISPBS(LL) <= 1 )THEN
      SWB(LN)=0.
      SUB(LN)=0.
      SUB(LEC(LN))=0.
      SCAY(LN)=0.
    END IF
    SAAY(LN)=0.
    SDY(LN)=0.
  ENDDO
  DO LL=1,NPBN
    I=IPBN(LL)
    J=JPBN(LL)
    L=LIJ(I,J)
    LPBN(LL)=L
    LS=LSC(L)
    SPB(L)=0.
    SUB(L)=0.
    SUB(LEC(L))=0.
    SWB(L)=0.
    SAAY(L)=0.
    SDY(L)=0.
    SCAX(L)=0.
    SCAX(LEC(L))=0.
    SAAX(L)=0.
    SAAX(LEC(L))=0.
    IF( ISPBN(LL) <= 1 )THEN
      SUB(LS)=0.
      SUB(LEC(LS))=0.
      SWB(LS)=0.
      SCAY(L)=0.
      SCAX(LS)=0.
      SCAX(LEC(LS))=0.
      SAAX(LS)=0.
      SAAX(LEC(LS))=0.
    END IF
  ENDDO
 
  ! *********************************************************************
  ! *** SET THE CELL FACES SWITCHES FOR HEAD CONTROL STRUCTURES
  ! *** UPSTREAM CONTROL
  DO NCTL=1,NQCTL
    IU=IQCTLU(NCTL)
    JU=JQCTLU(NCTL)
    L=LIJ(IU,JU)
    SAVESUB(1,NCTL) = SUB(L)
    SAVESVB(1,NCTL) = SVB(L)

    ! *** SET U FACE
    LW=LWC(L)
    DO IQ=1,NQCTL
      IF( IQ /= NCTL .AND. NQCTYP(IQ) /= 3 .AND. NQCTYP(IQ) /= 4 )THEN
        I=IQCTLU(IQ)
        J=JQCTLU(IQ)
        L1=LIJ(I,J)
        IF( L1 == LW )THEN
          SUB(L)=0.0
          EXIT
        ENDIF
      ENDIF
    ENDDO

    ! *** SET V FACE
    LS=LSC(L)
    DO IQ=1,NQCTL
      IF( IQ /= NCTL .AND. NQCTYP(IQ) /= 3 .AND. NQCTYP(IQ) /= 4 )THEN
        I=IQCTLU(IQ)
        J=JQCTLU(IQ)
        L1=LIJ(I,J)
        IF( L1 == LS )THEN
          SVB(L)=0.0
          EXIT
        ENDIF
      ENDIF
    ENDDO
  ENDDO

  ! *** DOWNSTREAM CONTROL
  DO NCTL=1,NQCTL
    ID=IQCTLD(NCTL)
    JD=JQCTLD(NCTL)
    IF( ID /= 0 .AND. JD /= 0 )THEN
      L=LIJ(ID,JD)
      SAVESUB(2,NCTL) = SUB(L)
      SAVESVB(2,NCTL) = SVB(L)

      ! *** SET U FACE
      LW=LWC(L)
      DO IQ=1,NQCTL
        IF( IQ /= NCTL .AND. NQCTYP(IQ) /= 3 .AND. NQCTYP(IQ) /= 4 )THEN
          I=IQCTLD(IQ)
          J=JQCTLD(IQ)
          IF( I > 0 .AND. J > 0 )THEN
            L1=LIJ(I,J)
            IF( L1 == LW )THEN
              SUB(L)=0.0
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDDO

      ! *** SET V FACE
      LS=LSC(L)
      DO IQ=1,NQCTL
        IF( IQ /= NCTL .AND. NQCTYP(IQ) /= 3 .AND. NQCTYP(IQ) /= 4 )THEN
          I=IQCTLD(IQ)
          J=JQCTLD(IQ)
          IF( I > 0 .AND. J > 0 )THEN
            L1=LIJ(I,J)
            IF( L1 == LS )THEN
              SVB(L)=0.0
              EXIT
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO

  ! ** RESET DXU,DYU,DXV,DYV BASED ON BOUNDARY CONDITION SWITCHES
  DO L=2,LA
    IF( SUB(L) > 0.5 )THEN
      DXU(L)=0.5*(DXP(L)+DXP(LWC(L)))
      DYU(L)=0.5*(DYP(L)+DYP(LWC(L)))
    ENDIF
    IF( SUB(L) < 0.5 .AND. SUB(LEC(L)) > 0.5 )THEN
      DXU(L)=DXP(L)
      DDYDDDX=2.*(DYP(LEC(L))-DYP(L))/(DXP(L)+DXP(LEC(L)))
      DYU(L)=DYP(L)-0.5*DXP(L)*DDYDDDX
    ENDIF
    IF( SUB(L) < 0.5 .AND. SUB(LEC(L)) < 0.5 )THEN
      DXU(L)=DXP(L)
      DYU(L)=DYP(L)
    ENDIF
  ENDDO
  DO L=2,LA
    LN=LNC(L)
    LS=LSC(L)
    IF( SVB(L) > 0.5 )THEN
      DXV(L)=0.5*(DXP(L)+DXP(LS))
      DYV(L)=0.5*(DYP(L)+DYP(LS))
    ENDIF
    IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
      DDXDDDY=2.*(DXP(LN)-DXP(L))/(DYP(L)+DYP(LN))
      DXV(L)=DXP(L)-0.5*DYP(L)*DDXDDDY
      DYV(L)=DYP(L)
    ENDIF
    IF( SVB(L) < 0.5 .AND. SVB(LN) < 0.5 )THEN
      DXV(L)=DXP(L)
      DYV(L)=DYP(L)
    ENDIF
  ENDDO

  ! **  SET THIN BARRIERS BY CALLING CELLMASK
  ! **  CALL MOVED FROM AAEFDC ON 23 JAN 2004
  IF( ISMASK == 1 ) CALL CELLMASK

  ! **  SET VOLUMETRIC & CONCENTRATION SOURCE LOCATIONS AND BED STRESS
  ! **  AND CELL CENTER BED STRESS AND VELOCITY MODIFERS
  DO LL=1,NQSIJ
    I=IQS(LL)
    J=JQS(LL)
    LTMP=LIJ(I,J)
    LQS(LL)=LTMP
    IF( NQSMUL(LL) == 0 )RQSMUL(LL)=1.
    IF( NQSMUL(LL) == 1 )RQSMUL(LL)=DYP(LTMP)
    IF( NQSMUL(LL) == 2 )RQSMUL(LL)=DXP(LTMP)
    IF( NQSMUL(LL) == 3 )RQSMUL(LL)=DXP(LTMP)+DYP(LTMP)
    IF( NQSMUL(LL) == 4 )RQSMUL(LL)=DXP(LTMP)*DYP(LTMP)
  ENDDO
  
  DO NCTL=1,NQCTL
    RQDW=1.
    IU=IQCTLU(NCTL)
    JU=JQCTLU(NCTL)
    LTMP=LIJ(IU,JU)
    IF( NQCMUL(NCTL) == 0 )RQCMUL(NCTL)=1.
    IF( NQCMUL(NCTL) == 1 )RQCMUL(NCTL)=DYP(LTMP)
    IF( NQCMUL(NCTL) == 2 )RQCMUL(NCTL)=DXP(LTMP)
    IF( NQCMUL(NCTL) == 3 )RQCMUL(NCTL)=DXP(LTMP)+DYP(LTMP)
    IF( NQCMUL(NCTL) == 4 )RQCMUL(NCTL)=DXP(LTMP)*DYP(LTMP)
  ENDDO

  ! *********************************************************************
  ! *** SET THE VELOCITY AND FLUX BOUNDARY CONDITIONS MULTIPLIERS

  ! *** DEFAULT CONDITION
  DO L=2,LA
    RSSBCE(L)=1.0
    RSSBCW(L)=1.0
    RSSBCN(L)=1.0
    RSSBCS(L)=1.0
    SUBEW(L) = SUB(L) + SUB(LEC(L))
    SVBNS(L) = SVB(L) + SVB(LNC(L))
  ENDDO

  ! *** STANDARD BORDER CELLS
  DO L=2,LA
    LE=LEC(L)
    LN=LNC(L)
    IF( SUBEW(L) < 1.5 )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
        RSSBCW(L)=0.0
        RSSBCE(L)=1.0
      ELSEIF( SUB(L) > 0.5 .AND. SUB(LE) < 0.5 )THEN
        RSSBCW(L)=1.0
        RSSBCE(L)=0.0
      ELSE
        RSSBCW(L)=0.0
        RSSBCE(L)=0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
        RSSBCS(L)=0.0
        RSSBCN(L)=1.0
      ELSEIF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
        RSSBCS(L)=1.0
        RSSBCN(L)=0.0
      ELSE
        RSSBCN(L)=0.0
        RSSBCS(L)=0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** FLOW BOUNDARY CONDITIONS
  DO LL=1,NQSIJ
    L=LQS(LL)
    LE=LEC(L)
    LN=LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) >  0.5 )THEN
        RSSBCW(L)=0.0
        RSSBCE(L)=2.0
      ENDIF
      IF( SUB(L) > 0.5 .AND. SUB(LE) <  0.5 )THEN
        RSSBCW(L)=2.0
        RSSBCE(L)=0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) >  0.5 )THEN
        RSSBCS(L)=0.0
        RSSBCN(L)=2.0
      ENDIF
      IF( SVB(L) >  0.5 .AND. SVB(LN) <  0.5 )THEN
        RSSBCS(L)=2.0
        RSSBCN(L)=0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** HYDRAULIC STRUCTURE: UPSTREAM
  DO NCTL=1,NQCTL
    IU=IQCTLU(NCTL)
    JU=JQCTLU(NCTL)
    L=LIJ(IU,JU)
    LE=LEC(L)
    LN=LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
        RSSBCW(L)=0.0
        RSSBCE(L)=2.0
      ENDIF
      IF( SUB(L) >  0.5 .AND. SUB(LE) < 0.5 )THEN
        RSSBCW(L)=2.0
        RSSBCE(L)=0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
        RSSBCS(L)=0.0
        RSSBCN(L)=2.0
      ENDIF
      IF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
        RSSBCS(L)=2.0
        RSSBCN(L)=0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** HYDRAULIC STRUCTURE: DOWNSTREAM
  DO NCTL=1,NQCTL
    ID=IQCTLD(NCTL)
    JD=JQCTLD(NCTL)
    IF( ID /= 0 .AND. JD /= 0 )THEN
      L=LIJ(ID,JD)
      LE=LEC(L)
      LN=LNC(L)
      ! *** ACCOUNT FOR DEAD END CHANNELS
      IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
        IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
          RSSBCW(L)=0.0
          RSSBCE(L)=2.0
        ENDIF
        IF( SUB(L) > 0.5 .AND. SUB(LE) < 0.5 )THEN
          RSSBCW(L)=2.0
          RSSBCE(L)=0.0
        ENDIF
      ENDIF
      IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
        IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
          RSSBCS(L)=0.0
          RSSBCN(L)=2.0
        ENDIF
        IF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
          RSSBCS(L)=2.0
          RSSBCN(L)=0.0
        ENDIF
      ENDIF
    END IF
  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: UPSTREAM
  DO NWR=1,NQWR
    IU=IQWRU(NWR)
    JU=JQWRU(NWR)
    L=LIJ(IU,JU)
    LE=LEC(L)
    LN=LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) < 0.5 .AND. SUB(LE) > 0.5 )THEN
        RSSBCW(L)=0.0
        RSSBCE(L)=2.0
      ENDIF
      IF( SUB(L) > 0.5 .AND. SUB(LE) < 0.5 )THEN
        RSSBCW(L)=2.0
        RSSBCE(L)=0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) < 1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) < 0.5 .AND. SVB(LN) > 0.5 )THEN
        RSSBCS(L)=0.0
        RSSBCN(L)=2.0
      ENDIF
      IF( SVB(L) > 0.5 .AND. SVB(LN) < 0.5 )THEN
        RSSBCS(L)=2.0
        RSSBCN(L)=0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: DOWNSTREAM
  DO NWR=1,NQWR
    ID=IQWRD(NWR)
    JD=JQWRD(NWR)
    L=LIJ(ID,JD)
    LE=LEC(L)
    LN=LNC(L)
    ! *** ACCOUNT FOR DEAD END CHANNELS
    IF( SUBEW(L) < 1.5 .AND. DXP(L) > DYP(L) .AND. SVB(L) == 0. .AND. SVB(LNC(L)) == 0. )THEN
      IF( SUB(L) <  0.5 .AND. SUB(LE) >  0.5 )THEN
        RSSBCW(L)=0.0
        RSSBCE(L)=2.0
      ENDIF
      IF( SUB(L) >  0.5 .AND. SUB(LE) <  0.5 )THEN
        RSSBCW(L)=2.0
        RSSBCE(L)=0.0
      ENDIF
    ENDIF
    IF( SVBNS(L) <  1.5 .AND. DYP(L) > DXP(L) .AND. SUB(L) == 0. .AND. SUB(LEC(L)) == 0. )THEN
      IF( SVB(L) <  0.5 .AND. SVB(LN) >  0.5 )THEN
        RSSBCS(L)=0.0
        RSSBCN(L)=2.0
      ENDIF
      IF( SVB(L) >  0.5 .AND. SVB(LN) <  0.5 )THEN
        RSSBCS(L)=2.0
        RSSBCN(L)=0.0
      ENDIF
    ENDIF
  ENDDO

  ! *** OPEN BOUNDARIES
  DO LL=1,NPBS
    I=IPBS(LL)
    J=JPBS(LL)
    L=LIJ(I,J)
    RSSBCS(L)=0.0
    RSSBCN(L)=2.0
  ENDDO
  DO LL=1,NPBW
    I=IPBW(LL)
    J=JPBW(LL)
    L=LIJ(I,J)
    RSSBCW(L)=0.0
    RSSBCE(L)=2.0
  ENDDO
  DO LL=1,NPBE
    I=IPBE(LL)
    J=JPBE(LL)
    L=LIJ(I,J)
    RSSBCW(L)=2.0
    RSSBCE(L)=0.0
  ENDDO
  DO LL=1,NPBN
    I=IPBN(LL)
    J=JPBN(LL)
    L=LIJ(I,J)
    RSSBCS(L)=2.0
    RSSBCN(L)=0.0
  ENDDO

  ! *********************************************************************
  ! *** SET BOUNDARY MOMENTUM SWITCHES FOR FLOW & HEAD CONTROL

  ! *** GLOBAL BOUNDARY CELL LIST
  NBCS=0

  ! *** FLOW BC'S
  DO LL=1,NQSIJ
    I=IQS(LL)
    J=JQS(LL)
    L=LIJ(I,J)
    NBCS=NBCS+1
    LBCS(NBCS)=L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS)=L
    IF( SUB(L) < 0.5 )THEN
      SAAX(L)=0.
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L)=0.
    ENDIF
    IF( L < LA-1 )THEN
      ! **                                                   L+2
      IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5)) THEN
        LBERC(NBCS)=LEC(L)
        SAAX(LBERC(NBCS))=BC_EDGEFACTOR
        SAAY(LBERC(NBCS))=BC_EDGEFACTOR
      ENDIF
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND. (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        SAAX(L)=BC_EDGEFACTOR
        SAAY(L)=BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS)=L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND. SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS)=LNC(L)
      SAAX(LBNRC(NBCS))=BC_EDGEFACTOR
      SAAY(LBNRC(NBCS))=BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND. (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      SAAX(L)=BC_EDGEFACTOR
      SAAY(L)=BC_EDGEFACTOR
    ENDIF

    IF( ABS(NQSMF(LL)) > 0 )THEN
      IF( QWIDTH(LL) <= 0.0 )THEN
        IF( ABS(NQSMF(LL)) == 1 .OR. ABS(NQSMF(LL)) == 3 ) QWIDTH(LL) = DYP(L)
        IF( ABS(NQSMF(LL)) == 2 .OR. ABS(NQSMF(LL)) == 4 ) QWIDTH(LL) = DXP(L)
      ENDIF
    ENDIF
  ENDDO

  ! *** HEAD CONTROL: UPSTREAM
  DO NCTL=1,NQCTL
    RQDW=1.
    IU=IQCTLU(NCTL)
    JU=JQCTLU(NCTL)
    L=LIJ(IU,JU)
    NBCS=NBCS+1
    LBCS(NBCS)=L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS)=L
    IF( SUB(L) < 0.5 )THEN
      SAAX(L)=0.
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L)=0.
    ENDIF
    IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
      LBERC(NBCS)=LEC(L)
      SAAX(LBERC(NBCS))=BC_EDGEFACTOR
      SAAY(LBERC(NBCS))=BC_EDGEFACTOR
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND. (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        SAAX(L)=BC_EDGEFACTOR
        SAAY(L)=BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS)=L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND. SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS)=LNC(L)
      SAAX(LBNRC(NBCS))=BC_EDGEFACTOR
      SAAY(LBNRC(NBCS))=BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND. (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      SAAX(L)=BC_EDGEFACTOR
      SAAY(L)=BC_EDGEFACTOR
    ENDIF

  ENDDO

  ! *** HEAD CONTROL: DOWNSTREAM
  DO NCTL=1,NQCTL
    ID=IQCTLD(NCTL)
    JD=JQCTLD(NCTL)
    IF( ID /= 0 .AND. JD /= 0 )THEN
      L=LIJ(ID,JD)
      NBCS=NBCS+1
      LBCS(NBCS)=L

      ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
      ! *** EAST/WEST MOMENTUM
      LBERC(NBCS)=L
      IF( SUB(L) < 0.5 )THEN
        SAAX(L)=0.
      ENDIF
      IF( SVB(L) < 0.5 )THEN
        SAAY(L)=0.
      ENDIF
      IF( L < LA-2 )THEN
        IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5)) THEN
          LBERC(NBCS)=LEC(L)
          SAAX(LBERC(NBCS))=BC_EDGEFACTOR
          SAAY(LBERC(NBCS))=BC_EDGEFACTOR
        ENDIF
      ENDIF
      IF( L > 2 )THEN
        IF( (SUB(L  ) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND.  &
       (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
          SAAX(L)=BC_EDGEFACTOR
          SAAY(L)=BC_EDGEFACTOR
        ENDIF
      ENDIF

      ! *** NORTH/SOUTH MOMENTUM
      LBNRC(NBCS)=L
      IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND.  &
                          SVB(LNC(LNC(L))) > 0.5) )THEN
        LBNRC(NBCS)=LNC(L)
        SAAX(LBNRC(NBCS))=BC_EDGEFACTOR
        SAAY(LBNRC(NBCS))=BC_EDGEFACTOR
      ENDIF
      IF( (SVB(L     ) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND.  &
       (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
        SAAX(L)=BC_EDGEFACTOR
        SAAY(L)=BC_EDGEFACTOR
      ENDIF
    ENDIF
  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: UPSTREAM
  DO NWR=1,NQWR
    IU=IQWRU(NWR)
    JU=JQWRU(NWR)
    L=LIJ(IU,JU)
    NBCS=NBCS+1
    LBCS(NBCS)=L

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS)=L
    IF( SUB(L) < 0.5 )THEN
      SAAX(L)=0.
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L)=0.
    ENDIF
    IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
      LBERC(NBCS)=LEC(L)
      SAAX(LBERC(NBCS))=BC_EDGEFACTOR
      SAAY(LBERC(NBCS))=BC_EDGEFACTOR
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L  ) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND.  &
       (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        SAAX(L)=BC_EDGEFACTOR
        SAAY(L)=BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS)=L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND.  &
                          SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS)=LNC(L)
      SAAX(LBNRC(NBCS))=BC_EDGEFACTOR
      SAAY(LBNRC(NBCS))=BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L     ) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND. (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      !LBNRC(NBCS)=LSC(L)
      !SAAX(LBNRC(NBCS))=BC_EDGEFACTOR
      !SAAY(LBNRC(NBCS))=BC_EDGEFACTOR
      SAAX(L)=BC_EDGEFACTOR
      SAAY(L)=BC_EDGEFACTOR
    ENDIF

  ENDDO

  ! *** WITHDRAWAL & RETURN BOUNDARY CONDITIONS: DOWNSTREAM
  DO NWR=1,NQWR
    ID=IQWRD(NWR)
    JD=JQWRD(NWR)
    L=LIJ(ID,JD)
    NBCS=NBCS+1
    LBCS(NBCS)=L
    !LBERC(NBCS)=1
    !LBNRC(NBCS)=1

    ! *** SET SAAX & SAAY FOR BOUNDARY MOMENTUM FLUXES
    ! *** EAST/WEST MOMENTUM
    LBERC(NBCS)=L
    IF( SUB(L) < 0.5 )THEN
      SAAX(L)=0.
    ENDIF
    IF( SVB(L) < 0.5 )THEN
      SAAY(L)=0.
    ENDIF
    IF( SUB(L) < 0.5 .AND. (SUB(LEC(L)) > 0.5 .AND. SUB(LEC(LEC(L))) > 0.5) )THEN
      LBERC(NBCS)=LEC(L)
      SAAX(LBERC(NBCS))=BC_EDGEFACTOR
      SAAY(LBERC(NBCS))=BC_EDGEFACTOR
    ENDIF
    IF( L > 2 )THEN
      IF( (SUB(L  ) > 0.5 .AND. SUB(LEC(L)) < 0.5) .AND.  &
       (SUB(LWC(L)) > 0.5 .AND. SUB(LWC(LWC(L))) > 0.5) )THEN
        !LBERC(NBCS)=LWC(L)
        !SAAX(LBERC(NBCS))=BC_EDGEFACTOR
        !SAAY(LBERC(NBCS))=BC_EDGEFACTOR
        SAAX(L)=BC_EDGEFACTOR
        SAAY(L)=BC_EDGEFACTOR
      ENDIF
    ENDIF

    ! *** NORTH/SOUTH MOMENTUM
    LBNRC(NBCS)=L
    IF( SVB(L) < 0.5 .AND. (SVB(LNC(L)) > 0.5 .AND. SVB(LNC(LNC(L))) > 0.5) )THEN
      LBNRC(NBCS)=LNC(L)
      SAAX(LBNRC(NBCS))=BC_EDGEFACTOR
      SAAY(LBNRC(NBCS))=BC_EDGEFACTOR
    ENDIF
    IF( (SVB(L     ) > 0.5 .AND. SVB(LNC(L)) < 0.5) .AND.  &
       (SVB(LSC(L)) > 0.5 .AND. SVB(LSC(LSC(L))) > 0.5) )THEN
      !LBNRC(NBCS)=LSC(L)
      !SAAX(LBNRC(NBCS))=BC_EDGEFACTOR
      !SAAY(LBNRC(NBCS))=BC_EDGEFACTOR
      SAAX(L)=BC_EDGEFACTOR
      SAAY(L)=BC_EDGEFACTOR
    ENDIF
  ENDDO

  ! *** OPEN BOUNDARIES
  NBCSOP=0
  NBCSOP2=0
  DO LL=1,NPBS
    I=IPBS(LL)
    J=JPBS(LL)
    L=LIJ(I,J)
    NBCSOP=NBCSOP+1
    LOBCS(NBCSOP)=L
    NBCS=NBCS+1
    LBCS(NBCS)=L
    LBERC(NBCS)=L
    LBNRC(NBCS)=LNC(L)
    IF( ISPBS(LL) <= 1 )THEN
      NBCSOP2=NBCSOP2+1
      LOBCS2(NBCSOP2)=LNC(L)
    ENDIF
  ENDDO
  DO LL=1,NPBW
    I=IPBW(LL)
    J=JPBW(LL)
    L=LIJ(I,J)
    NBCSOP=NBCSOP+1
    LOBCS(NBCSOP)=L
    NBCS=NBCS+1
    LBCS(NBCS)=L
    LBERC(NBCS)=LEC(L)
    LBNRC(NBCS)=L
    IF( ISPBW(LL) <= 1 )THEN
      NBCSOP2=NBCSOP2+1
      LOBCS2(NBCSOP2)=LEC(L)
    ENDIF
  ENDDO
  DO LL=1,NPBE
    I=IPBE(LL)
    J=JPBE(LL)
    L=LIJ(I,J)
    NBCSOP=NBCSOP+1
    LOBCS(NBCSOP)=L
    NBCS=NBCS+1
    LBCS(NBCS)=L
    LBERC(NBCS)=LWC(L)
    LBNRC(NBCS)=L
    IF( ISPBE(LL) <= 1 )THEN
      NBCSOP2=NBCSOP2+1
      LOBCS2(NBCSOP2)=LWC(L)
    ENDIF
  ENDDO
  DO LL=1,NPBN
    I=IPBN(LL)
    J=JPBN(LL)
    L=LIJ(I,J)
    NBCSOP=NBCSOP+1
    LOBCS(NBCSOP)=L
    NBCS=NBCS+1
    LBCS(NBCS)=L
    LBERC(NBCS)=L
    LBNRC(NBCS)=LSC(L)
    IF( ISPBN(LL) <= 1 )THEN
      NBCSOP2=NBCSOP2+1
      LOBCS2(NBCSOP2)=LSC(L)
    ENDIF
  ENDDO

  ! *********************************************************************
  ! ***  SET OPEN BOUNDARY FLAGS FOR CONSTITUENTS
  DO LL=1,NCBS
    I=ICBS(LL)
    J=JCBS(LL)
    LCBS(LL)=LIJ(I,J)
    L=LIJ(I,J)
    SCB(L)=0.
  ENDDO
  DO LL=1,NCBW
    I=ICBW(LL)
    J=JCBW(LL)
    LCBW(LL)=LIJ(I,J)
    L=LIJ(I,J)
    SCB(L)=0.
  ENDDO
  DO LL=1,NCBE
    I=ICBE(LL)
    J=JCBE(LL)
    LCBE(LL)=LIJ(I,J)
    L=LIJ(I,J)
    SCB(L)=0.
  ENDDO
  DO LL=1,NCBN
    I=ICBN(LL)
    J=JCBN(LL)
    LCBN(LL)=LIJ(I,J)
    L=LIJ(I,J)
    SCB(L)=0.
  ENDDO

  ! *********************************************************************
  ! ***  SET JET-PLUME VOLUMES SOURCES
  DO NJP=1,NQJPIJ
    L=LIJ(IQJP(NJP),JQJP(NJP))
    NBCS=NBCS+1
    LBCS(NBCS)=L
    LBERC(NBCS)=1
    LBNRC(NBCS)=1

    IF( ICALJP(NJP) == 2 )THEN
      ! *** WITHDRAWAL CELL
      L=LIJ(IUPCJP(NJP),JUPCJP(NJP))
      NBCS=NBCS+1
      LBCS(NBCS)=L
      LBERC(NBCS)=1
      LBNRC(NBCS)=1
    ENDIF
  ENDDO

  ! **  SET CHANNEL HOST AND GUEST LOCATION MAPPINGS
  IF( MDCHH >= 1 )THEN
    DO NMD=1,MDCHH
      L=LIJ(IMDCHH(NMD),JMDCHH(NMD))
      ! *** HOST
      LMDCHH(NMD)=L
      NBCS=NBCS+1
      LBCS(NBCS)=L
      
      ! *** DOWNSTREAM
      IF( IMDCHU(NMD) == 1 .AND. JMDCHU(NMD) == 1 )THEN
        LMDCHU(NMD)=1
      ELSE
        L=LIJ(IMDCHU(NMD),JMDCHU(NMD))
        LMDCHU(NMD)=L
      ENDIF
      IF( IMDCHV(NMD) == 1 .AND. JMDCHV(NMD) == 1 )THEN
        LMDCHV(NMD)=1
      ELSE
        L=LIJ(IMDCHV(NMD),JMDCHV(NMD))
        LMDCHV(NMD)=L
      ENDIF
      NBCS=NBCS+1
      LBCS(NBCS)=L
    ENDDO
  ENDIF

  ! **  SET CELL FACE WET DEPTHS
  HUWET(1)=HWET
  HUWET(LC)=HWET
  HVWET(1)=HWET
  HVWET(LC)=HWET
  HUDRY(1)=HDRY
  HUDRY(LC)=HDRY
  HVDRY(1)=HDRY
  HVDRY(LC)=HDRY
  DO L=2,LA
    LS=LSC(L)
    HUDRY(L)=HDRY+0.5*ABS(BELV(L)-BELV(LWC(L)))
    HVDRY(L)=HDRY+0.5*ABS(BELV(L)-BELV(LS))
    HUWET(L)=HWET+0.5*ABS(BELV(L)-BELV(LWC(L)))
    HVWET(L)=HWET+0.5*ABS(BELV(L)-BELV(LS))
  ENDDO
  
  IF( ISDRY > 0 )THEN
    NDRYTMP=MOD(ISDRY,2)
    IF( NDRYTMP /= 0 )THEN
      DO L=2,LA
        HUWET(L)=HWET
        HVWET(L)=HWET
        HUDRY(L)=HDRY
        HVDRY(L)=HDRY
      ENDDO
    ENDIF
  ENDIF

  ! *** SET PERMANENT FACE SWITCHES
  DO L=1,LC
    SUBO(L)=SUB(L)
    SVBO(L)=SVB(L)
  ENDDO
  
  ! **  DIAGNOSTIC OUTPUT
  IF( DEBUG )THEN
    OPEN(1,FILE=OUTDIR//'SETBC.DIA',STATUS='UNKNOWN')
    CLOSE(1,STATUS='DELETE')
    OPEN(1,FILE=OUTDIR//'SETBC.DIA')
    DO L=2,LA
      WRITE(1,1001)IL(L),JL(L),SUB(L),SUB(LEC(L)),SVB(L),SVB(LNC(L)),SPB(L)
    ENDDO
    CLOSE(1)
  ENDIF
  
  1001 FORMAT(2I5,8E13.4)
  
  RETURN
END

