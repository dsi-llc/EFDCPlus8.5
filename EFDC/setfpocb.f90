SUBROUTINE SETFPOCB(ITYPE)

  ! **  SUBROUTINE SETFPOC SETS FPOC BASED ON SEDIMENT BED COMPOSITION
  ! **  FUNCTIONAL RELATIONSHIPS FOR HOUSATONIC RIVER PROVIDED BY
  ! **  PROVIDED BY HYDROQUAL DOCUMENT DATED MARCH 17, 2003
  ! **  NOTE THAT PER CENTS IN PROVIDED RELATIONS ARE CONVERTED
  ! **  TO FRACTIONS CONSISTENT WITH EFDC VARIABLES
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !**********************************************************************!

  USE GLOBAL

  IMPLICIT NONE

  INTEGER :: L,K,ITYPE
  REAL    :: TMPVAL,TMPVAL1,TMPVAL2

  !###########################################################################
  ! HQI Change to include spatially varying, but time constant bulk foc and
  ! pseudo-foc
  REAL PREDFOCB
  ! RM, 02/29/04
  !###########################################################################

  !
  !
  !**********************************************************************C
  !   NOTE LIMIT VALUES TO LESS THAN 100 %
  IF( ITYPE == 0 )THEN

    !     INITIALIZE ALL FULL BED LAYERS
    !
    !     COHESIVE SEDIMENT (CLASS 1)
    DO K=1,KB
      DO L=2,LA
        IF( K <= KBT(L) )THEN
          IF( VFRBED(L,K,1) > 0.0 )THEN
          !#######################################################################
          !     09/18/03, RM, HQI, relationship for cohesives corrected and
          !     regression based foc limited to 100%
          !     02/29/04, RM, HQI, log-space to arithmetic space conversion
          !             STFPOCB(L,K,1)=0.8657/((100.*VFRBED(L,K,1))**1.3326)
            STFPOCB(L,K,1)=0.0384*((100.*VFRBED(L,K,1))**0.0366)
            STFPOCB(L,K,1)=LOG(STFPOCB(L,K,1))
            STFPOCB(L,K,1)=STFPOCB(L,K,1) + ((0.5955**2)/2)
            STFPOCB(L,K,1)=EXP(STFPOCB(L,K,1))
            IF( STFPOCB(L,K,1) > 1.) STFPOCB(L,K,1)=1.0
          !#######################################################################
          ELSE
            STFPOCB(L,K,1)=0.0
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !     COARSER NONCOHESIVE SEDIMENT (CLASS 2)
    !#######################################################################

    DO K=1,KB
      DO L=2,LA
        IF( K <= KBT(L) )THEN
          IF( VFRBED(L,K,2) > 0.0 )THEN
            !#######################################################################
            !     09/18/03, RM, HQI, Change to account for NC4 and
            !     regression based foc limited to 100%
            !     02/29/04, RM, HQI, log-space to arithmetic space conversion
            TMPVAL2=( (100.*VFRBED(L,K,2))**0.9189)
            TMPVAL=VFRBED(L,K,1)/ &
               (VFRBED(L,K,2)+VFRBED(L,K,3)+VFRBED(L,K,4))
            TMPVAL1=TMPVAL**0.9204
            STFPOCB(L,K,2)=0.7428*TMPVAL1/TMPVAL2
            STFPOCB(L,K,2)=LOG(STFPOCB(L,K,2))
            STFPOCB(L,K,2)=STFPOCB(L,K,2) + ((1.1371**2)/2)
            STFPOCB(L,K,2)=EXP(STFPOCB(L,K,2))
            IF( STFPOCB(L,K,2) > 1.) STFPOCB(L,K,2)=1.0
            !#######################################################################
          ELSE
            STFPOCB(L,K,2)=0.0
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !     FINER NONCOHESIVE SEDIMENT (CLASS 2)
    !#######################################################################
    !     09/18/03, RM, HQI, this is Non-Cohesive 2
    !#######################################################################
    DO K=1,KB
      DO L=2,LA
        IF( K <= KBT(L) )THEN
          IF( (VFRBED(L,K,3)+VFRBED(L,K,4)) > 0.0 )THEN
            !#######################################################################
            !     09/18/03, RM, HQI, regression corrected, account for NC4 and
            !     regression based foc limited to 100%
            !     02/29/04, RM, HQI, log-space to arithmetic space conversion
            !             STFPOCB(L,K,3)=0.033*((100.*VFRBED(L,K,3))**0.152)
            STFPOCB(L,K,3)=0.4084/ &
               ((100.*(VFRBED(L,K,3)+VFRBED(L,K,4)))**1.116)
            STFPOCB(L,K,3)=LOG(STFPOCB(L,K,3))
            STFPOCB(L,K,3)=STFPOCB(L,K,3) + ((1.4436**2)/2)
            STFPOCB(L,K,3)=EXP(STFPOCB(L,K,3))
            IF( STFPOCB(L,K,3) > 1.) STFPOCB(L,K,3)=1.0
            STFPOCB(L,K,4)=STFPOCB(L,K,3)
            !#######################################################################
          ELSE
            STFPOCB(L,K,3)=0.0
            STFPOCB(L,K,4)=0.0
          ENDIF
        ENDIF
      ENDDO
    ENDDO

    !###########################################################################
    ! HQI Change to include spatially varying, but time constant bulk foc and
    ! pseudo-foc
    DO K=1,KB
      DO L=2,LA
        IF( K <= KBT(L) )THEN
          PREDFOCB = (STFPOCB(L,K,1)*VFRBED(L,K,1)) + &
                  (STFPOCB(L,K,2)*VFRBED(L,K,2)) + &
                  (STFPOCB(L,K,3)*VFRBED(L,K,3)) + &
                  (STFPOCB(L,K,4)*VFRBED(L,K,4))
          IF( PREDFOCB > 0.0 )THEN
            STFPOCB(L,K,1) = STFPOCB(L,K,1)*PFPOCB(L,K)/PREDFOCB
            STFPOCB(L,K,2) = STFPOCB(L,K,2)*PFPOCB(L,K)/PREDFOCB
            STFPOCB(L,K,3) = STFPOCB(L,K,3)*PFPOCB(L,K)/PREDFOCB
            STFPOCB(L,K,4) = STFPOCB(L,K,4)*PFPOCB(L,K)/PREDFOCB
          ENDIF
          IF( STFPOCB(L,K,1) > 1.) STFPOCB(L,K,1)=1.0
          IF( STFPOCB(L,K,2) > 1.) STFPOCB(L,K,2)=1.0
          IF( STFPOCB(L,K,3) > 1.) STFPOCB(L,K,3)=1.0
          IF( STFPOCB(L,K,4) > 1.) STFPOCB(L,K,4)=1.0
        ENDIF
      ENDDO
    ENDDO
    ! RM, 02/29/04
    !###########################################################################
  ENDIF

  !C**********************************************************************C
  IF( ITYPE == 1 )THEN
    !     UPDATE TOP LAYER OF BED
    !
    !     COHESIVE SEDIMENT (CLASS 1)
    !
    DO L=2,LA
      K=KBT(L)
      IF( VFRBED(L,K,1) > 0.0 )THEN
        !#######################################################################
        !     09/18/03, RM, HQI, relationship for cohesives corrected and
        !     regression based foc limited to 100%
        !     02/29/04, RM, HQI, log-space to arithmetic space conversion
        !           STFPOCB(L,K,1)=0.8657/((100.*VFRBED(L,K,1))**1.3326)
        STFPOCB(L,K,1)=0.0384*((100.*VFRBED(L,K,1))**0.0366)
        STFPOCB(L,K,1)=LOG(STFPOCB(L,K,1))
        STFPOCB(L,K,1)=STFPOCB(L,K,1) + ((0.5955**2)/2)
        STFPOCB(L,K,1)=EXP(STFPOCB(L,K,1))
        IF( STFPOCB(L,K,1) > 1.) STFPOCB(L,K,1)=1.0
        !#######################################################################
      ELSE
        STFPOCB(L,K,1)=0.0
      ENDIF
    ENDDO

    !     COARSER NONCOHESIVE SEDIMENT (CLASS 2)
    !#######################################################################
    !     09/18/03, RM, HQI, this is Non-Cohesive 1
    !#######################################################################

    DO L=2,LA
      K=KBT(L)
      IF( VFRBED(L,K,2) > 0.0 )THEN
        !#######################################################################
        !     09/18/03, RM, HQI, Change to account for NC4 and
        !     regression based foc limited to 100%
        !     02/29/04, RM, HQI, log-space to arithmetic space conversion
        !           TMPVAL=VFRBED(L,K,3)/(VFRBED(L,K,1)+VFRBED(L,K,2))
        TMPVAL2=( (100.*VFRBED(L,K,2))**0.9189)
        TMPVAL=VFRBED(L,K,1)/ &
               (VFRBED(L,K,2)+VFRBED(L,K,3)+VFRBED(L,K,4))
        TMPVAL1=TMPVAL**0.9204
        STFPOCB(L,K,2)=0.7428*TMPVAL1/TMPVAL2
        STFPOCB(L,K,2)=LOG(STFPOCB(L,K,2))
        STFPOCB(L,K,2)=STFPOCB(L,K,2) + ((1.1371**2)/2)
        STFPOCB(L,K,2)=EXP(STFPOCB(L,K,2))
        IF( STFPOCB(L,K,2) > 1.) STFPOCB(L,K,2)=1.0
        !#######################################################################
      ELSE
        STFPOCB(L,K,2)=0.0
      ENDIF
    ENDDO

    !     FINER NONCOHESIVE SEDIMENT (CLASS 2)
    !#######################################################################
    !     09/18/03, RM, HQI, this is Non-Cohesive 2
    !#######################################################################

    DO L=2,LA
      K=KBT(L)
      IF( VFRBED(L,K,3) > 0.0 )THEN
        !#######################################################################
        !     09/18/03, RM, HQI, regression corrected, account for NC4 and
        !     regression based foc limited to 100%
        !     02/29/04, RM, HQI, log-space to arithmetic space conversion
        !           STFPOCB(L,K,3)=0.033*((100.*VFRBED(L,K,3))**0.152)
        STFPOCB(L,K,3)=0.4084/ &
               ((100.*(VFRBED(L,K,3)+VFRBED(L,K,4)))**1.116)
        STFPOCB(L,K,3)=LOG(STFPOCB(L,K,3))
        STFPOCB(L,K,3)=STFPOCB(L,K,3) + ((1.4436**2)/2)
        STFPOCB(L,K,3)=EXP(STFPOCB(L,K,3))
        IF( STFPOCB(L,K,3) > 1.) STFPOCB(L,K,3)=1.0
        STFPOCB(L,K,4)=STFPOCB(L,K,3)
        !#######################################################################
      ELSE
        STFPOCB(L,K,3)=0.0
        STFPOCB(L,K,4)=0.0
      ENDIF
    ENDDO

    !###########################################################################
    ! HQI Change to include spatially varying, but time constant bulk foc and
    ! pseudo-foc
    DO L=2,LA
      K=KBT(L)
      PREDFOCB = (STFPOCB(L,K,1)*VFRBED(L,K,1)) + &
            (STFPOCB(L,K,2)*VFRBED(L,K,2)) + &
            (STFPOCB(L,K,3)*VFRBED(L,K,3)) + &
            (STFPOCB(L,K,4)*VFRBED(L,K,4))
      IF( PREDFOCB > 0.0 )THEN
        STFPOCB(L,K,1) = STFPOCB(L,K,1)*PFPOCB(L,K)/PREDFOCB
        STFPOCB(L,K,2) = STFPOCB(L,K,2)*PFPOCB(L,K)/PREDFOCB
        STFPOCB(L,K,3) = STFPOCB(L,K,3)*PFPOCB(L,K)/PREDFOCB
        STFPOCB(L,K,4) = STFPOCB(L,K,4)*PFPOCB(L,K)/PREDFOCB
      ENDIF
      IF( STFPOCB(L,K,1) > 1.) STFPOCB(L,K,1)=1.0
      IF( STFPOCB(L,K,2) > 1.) STFPOCB(L,K,2)=1.0
      IF( STFPOCB(L,K,3) > 1.) STFPOCB(L,K,3)=1.0
      IF( STFPOCB(L,K,4) > 1.) STFPOCB(L,K,4)=1.0
    ENDDO
    ! RM, 02/29/04
    !###########################################################################
    !
  ENDIF

  !**********************************************************************C
  !
  RETURN

END
