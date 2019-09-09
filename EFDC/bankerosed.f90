SUBROUTINE BANKEROSED

  ! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION 1.0a
  !
  ! **  LAST MODIFIED BY JOHN HAMRICK ON 1 NOVEMBER 2001
  !
  !----------------------------------------------------------------------C
  !
  ! CHANGE RECORD
  ! DATE MODIFIED     BY                 DATE APPROVED    BY
  !----------------------------------------------------------------------C
  !
  !**********************************************************************C
  !
  ! **  SUBROUTINE BANKEROSED CALCULATES SEDIMENT TRANSPORT DUE TO BANK
  ! **  EROSION.  TRANSPORT IS FROM BANK BED CELL TO CHANNEL BED AND
  ! **  WATER COLUMN CELLS
  !
  !**********************************************************************C

  USE GLOBAL

  IMPLICIT NONE

  INTEGER :: L,NS,NT,LBANK,LCHAN,K,NX,NP
  REAL    :: TIME,BKEROBKB,BKEROCHW,BKEROCHB,WVEL,TMPVAL

  !**********************************************************************C
  !
  IF( ISDYNSTP == 0 )THEN
    TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON
  ELSE
    TIME=TIMESEC/TCON
  ENDIF

  !**********************************************************************C
  !
  !  INITIALIZE BANK EROSION VARIABLES
  !
  DO L=1,LC
    QSBDTOPBEBKB(L)=0.
    QSBDTOPBECHB(L)=0.
    QSBDTOPBECHW(L)=0.
    QWBDTOPBEBKB(L)=0.
    QWBDTOPBECHB(L)=0.
    QWBDTOPBECHW(L)=0.
  ENDDO
  !
  DO NS=1,NSED
  DO L=1,LC
    SEDFBEBKB(L,NS)=0.
    SEDFBECHB(L,NS)=0.
    SEDFBECHW(L,NS)=0.
  ENDDO
  ENDDO
  !
  DO NS=1,NSND
  DO L=1,LC
    SNDFBEBKB(L,NS)=0.
    SNDFBECHB(L,NS)=0.
    SNDFBECHW(L,NS)=0.
  ENDDO
  ENDDO
  !
  DO NT=1,NTOX
  DO L=1,LC
    TOXFBEBKB(L,NT)=0.
    TOXFBECHB(L,NT)=0.
    TOXFBECHW(L,NT)=0.
  ENDDO
  ENDDO
  !
  !**********************************************************************C
  !
  !  LOAD SEDIMENT FLUXES
  !
  IF( ISTRAN(6) > 0 )THEN
    DO NS=1,NSED
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        BKEROBKB=VFRBED(LBANK,K,NS)*FBESER(NP) &
                 *BESERT(NBESERN(NP))
        BKEROCHW=FWCBESERT(NBESERN(NP))
        BKEROCHB=1.-FWCBESERT(NBESERN(NP))
  !**********************************************************************C
  ! HQI Change to restrict bank erosion to cells with available sediment
  ! RM, 01/08/07
  !            SEDFBEBKB(NP,NS)=BKEROBKB*DXYIP(LBANK)
  !            SEDFBECHB(NP,NS)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
  !            SEDFBECHW(NP,NS)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        IF( (DELT*BKEROBKB*DXYIP(LBANK)) > SEDB(LBANK,KBT(LBANK),NS) )THEN

          SEDFBEBKB(NP,NS)=0.
          SEDFBECHB(NP,NS)=0.
          SEDFBECHW(NP,NS)=0.
        ELSE
          SEDFBEBKB(NP,NS)=BKEROBKB*DXYIP(LBANK)
          SEDFBECHB(NP,NS)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
          SEDFBECHW(NP,NS)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        ENDIF
  ! End HQI Change
  !**********************************************************************C
      ENDDO
    ENDDO
  ENDIF
  !
  IF( ISTRAN(7) > 0 )THEN
    DO NX=1,NSND
    NS=NSED+NX
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        BKEROBKB=VFRBED(LBANK,K,NS)*FBESER(NP) &
                 *BESERT(NBESERN(NP))
        BKEROCHW=FWCBESERT(NBESERN(NP))
        BKEROCHB=1.-FWCBESERT(NBESERN(NP))
  !**********************************************************************C
  ! HQI Change to restrict bank erosion to cells with available sediment
  ! RM, 01/08/07
  !     SNDFBEBKB(NP,NX)=BKEROBKB*DXYIP(LBANK)
  !     SNDFBECHB(NP,NX)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
  !     SNDFBECHW(NP,NX)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        IF( (DELT*BKEROBKB*DXYIP(LBANK)) > SNDB(LBANK,KBT(LBANK),NX) )THEN

          SNDFBEBKB(NP,NX)=0.
          SNDFBECHB(NP,NX)=0.
          SNDFBECHW(NP,NX)=0.
        ELSE
          SNDFBEBKB(NP,NX)=BKEROBKB*DXYIP(LBANK)
          SNDFBECHB(NP,NX)=-BKEROCHB*BKEROBKB*DXYIP(LCHAN)
          SNDFBECHW(NP,NX)=BKEROCHW*BKEROBKB*DXYIP(LCHAN)
        ENDIF
  ! End HQI Change
  !**********************************************************************C
      ENDDO
    ENDDO
  ENDIF
  !
  !**********************************************************************C
  !
  !  UPDATE BED AND WATER COLUMN SEDIMENT CONCENTRATION
  !
  IF( ISTRAN(6) > 0 )THEN
    DO NS=1,NSED
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        WVEL=DELT*HPI(LCHAN)*DZIC(LCHAN,KSZ(LCHAN))
        SEDB(LBANK,KBT(LBANK),NS)=SEDB(LBANK,KBT(LBANK),NS) &
                                 -DELT*SEDFBEBKB(NP,NS)
        SEDB(LCHAN,KBT(LCHAN),NS)=SEDB(LCHAN,KBT(LCHAN),NS) &
                                 -DELT*SEDFBECHB(NP,NS)
        SED(LCHAN,KSZ(LCHAN),NS)=SED(LCHAN,1,NS)+WVEL*SEDFBECHW(NP,NS)
      ENDDO
    ENDDO
  ENDIF
  !
  IF( ISTRAN(7) > 0 )THEN
    DO NS=1,NSND
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        WVEL=DELT*HPI(LCHAN)*DZIC(LCHAN,KSZ(LCHAN))
        SNDB(LBANK,KBT(LBANK),NS)=SNDB(LBANK,KBT(LBANK),NS) &
                                 -DELT*SNDFBEBKB(NP,NS)
        SNDB(LCHAN,KBT(LCHAN),NS)=SNDB(LCHAN,KBT(LCHAN),NS) &
                                 -DELT*SNDFBECHB(NP,NS)
        SND(LCHAN,KSZ(LCHAN),NS)=SND(LCHAN,1,NS)+WVEL*SNDFBECHW(NP,NS)
      ENDDO
    ENDDO
  ENDIF
  !
  !**********************************************************************C
  !
  !  CALCULATE SEDIMENT AND WATER VOLUME FLUXES FROM BANK
  !
  !  COHESIVE
  !
  IF( ISTRAN(6) > 0 )THEN
  IF( IBMECH == 1 .AND. SEDVRDT < 0.0 )THEN
  !
    DO NS=1,NSED
      DO NP=1,NBEPAIR
        LBANK = LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN = LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        QSBDTOPBEBKB(NP) = QSBDTOPBEBKB(NP) + 0.001*SEDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
        QWBDTOPBEBKB(NP) = QWBDTOPBEBKB(NP) + 0.001*VDRBED(LBANK,K)*SEDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
    ENDDO
    ENDDO
  !
  ELSE
  !
    DO NS=1,NSED
      DSEDGMM=1./(1.E6*SSG(NS))
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        QSBDTOPBEBKB(NP)=QSBDTOPBEBKB(NP) &
          +DSEDGMM*SEDFBEBKB(NP,NS)
        QWBDTOPBEBKB(NP)=QWBDTOPBEBKB(NP) &
          +DSEDGMM*VDRBED(LBANK,K)*SEDFBEBKB(NP,NS)
      ENDDO
    ENDDO
  !
  ENDIF
  ENDIF
  !
  !  NONCOHESIVE
  !
  IF( ISTRAN(7) > 0 )THEN
  IF( IBMECH == 1 .AND. SEDVRDT < 0.0 )THEN
  !
    DO NS=1,NSND
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        QSBDTOPBEBKB(NP)=QSBDTOPBEBKB(NP) &
        +0.001*SNDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
        QWBDTOPBEBKB(NP)=QWBDTOPBEBKB(NP) &
        +0.001*VDRBED(LBANK,K)*SNDFBEBKB(NP,NS)/SDENAVG(LBANK,K)
      ENDDO
    ENDDO
  !
  ELSE
  !
    DO NS=1,NSND
      DSEDGMM=1./(1.E6*SSG(NS+NSED))
      DO NP=1,NBEPAIR
        LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
        LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
        K=KBT(LBANK)
        QSBDTOPBEBKB(NP)=QSBDTOPBEBKB(NP) &
          +DSEDGMM*SNDFBEBKB(NP,NS)
        QWBDTOPBEBKB(NP)=QWBDTOPBEBKB(NP) &
          +DSEDGMM*VDRBED(LBANK,K)*SNDFBEBKB(NP,NS)
      ENDDO
    ENDDO
  !
  ENDIF
  ENDIF
  !
  !**********************************************************************C
  !
  !  CALCULATE SEDIMENT AND WATER VOLUME FLUXES FOR CHANNEL BED AND
  !  WATER COLUMN
  !
  DO NP=1,NBEPAIR
    LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
    LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
    BKEROCHW=FWCBESERT(NBESERN(NP))
    BKEROCHB=1.-FWCBESERT(NBESERN(NP))
    TMPVAL=DXYP(LBANK)*DXYIP(LCHAN)
    QSBDTOPBECHB(NP)=-TMPVAL*BKEROCHB*QSBDTOPBEBKB(NP)
    QWBDTOPBECHB(NP)=-TMPVAL*BKEROCHB*QWBDTOPBEBKB(NP)
    QSBDTOPBECHW(NP)=TMPVAL*BKEROCHW*QSBDTOPBEBKB(NP)
    QWBDTOPBECHW(NP)=TMPVAL*BKEROCHW*QWBDTOPBEBKB(NP)
  ENDDO
  !
  !**********************************************************************C
  !
  RETURN
END
