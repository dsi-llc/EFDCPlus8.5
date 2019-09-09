SUBROUTINE CALBED  

  !**********************************************************************!
  ! ***  SUBROUTINE CALSND CALCULATES NONCOHESIVE SEDIMENT SETTLING,  
  ! ***  DEPOSITION AND RESUSPENSION AND IS CALLED FOR SSEDTOX  
  ! ***  NOT USED FOR SEDZLJ
  !**********************************************************************!
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-03           Paul M. Craig      Rewritten to F90 and added OMP
  !----------------------------------------------------------------------!

  USE GLOBAL  

  IMPLICIT NONE

  INTEGER :: K,L,IFLAG,NS,KK,NX,KBTM1,ND,LF,LL
  REAL :: TMPVAL,WDENKGM3,WDENGMM3,TMPVALK
  REAL :: TMPVALKP,BETTMP,VOIDCON1,HBEDTMP,TMPVALO
  REAL :: TMPVALN,TMPEXP,TMPTOP,TMPBOT,FSTRSE,FDSTRSE
  REAL :: FHYDCN,DSTRESET,POROTMP

  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: SNDHYDCN  

  IF(  .NOT. ALLOCATED(SNDHYDCN) )THEN
    ALLOCATE(SNDHYDCN(LCM,KBM))  
    SNDHYDCN = 0.0 
  ENDIF

  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN  
    HBEDMIN = 1.E-4  
    IF( ISTRAN(7) >= 1 )THEN  
      HBEDMIN = MAX(HBEDMIN,SNDDMX)  
    END IF  

    ! *********************************************************************
    ! *** CONSTANT POROSITY BED  
    IF( IBMECH == 0 )THEN

      ! *** UPDATE TOP LAYER THICKNESS TO MAINTAIN CONSTANT POROSITY  
      VOIDCON1 = BEDPORC/(1.-BEDPORC)  

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,L,LF,LL)              &
      !$OMP             PRIVATE(HBEDTMP,TMPVALO,TMPVALN)
      DO ND=1,NDM  
        LF = 2 + (ND-1)*LDM  
        LL = MIN(LF+LDM-1,LA)

        DO L=LF,LL  
          K = KBT(L)  
          HBEDTMP = (1. + VOIDCON1)*HBED(L,K)/(1. + VDRBED(L,K))  
          TMPVALO = VDRBED(L,K)*HBED(L,K)/(1. + VDRBED(L,K))  
          TMPVALN = VOIDCON1*HBEDTMP/(1. + VOIDCON1)  
          QWBDTOP(L) = DELTI*(TMPVALO-TMPVALN)  
          HBED(L,K) = HBEDTMP  
          QWTRBED(L,K) = QWBDTOP(L) + QGW(L)/DXYP(L) 
          VDRBED(L,K) = VOIDCON1     ! *** PMC 2010-12-06
        ENDDO  
        DO K=0,KBT(L)-1  
          DO L=LF,LL  
            QWTRBED(L,K) = QGW(L)/DXYP(L)  
          ENDDO  
        END DO  

      ENDDO   ! *** END OF DOMAIN LOOP
      !$OMP END PARALLEL DO
      
    ENDIF  

    ! *********************************************************************
    ! *** SIMPLE CONSOLIDATING BED  
    IF( IBMECH == 1 )THEN

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,L,LF,LL)              &
      !$OMP             PRIVATE(TMPEXP,TMPTOP,TMPBOT,TMPVAL)
      DO ND=1,NDM  
        LF = 2+(ND-1)*LDM  
        LL = MIN(LF+LDM-1,LA)

        ! ***  IF SEDVRDT>0.0001 CONSOLIDATE TO SEDVRM (THE MINIMUM VOID RATIO)  
        IF( SEDVRDT > 0.0001 )THEN  
          TMPEXP = EXP(-DTSED/SEDVRDT)  
          DO K=1,KB  
            DO L=LF,LL  
              IF( K <= KBT(L) )THEN  
                VDRBED1(L,K) = VDRBED(L,K)  
                HBED1(L,K) = HBED(L,K)  
                VDRBED(L,K) = SEDVDRM + (VDRBED1(L,K)-SEDVDRM)*TMPEXP  
                TMPTOP = 1. + VDRBED(L,K)  
                TMPBOT = 1. + VDRBED1(L,K)  
                HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT  
              ENDIF  
            ENDDO  
          ENDDO  
        ENDIF  

        ! ***  IF SEDVRDT>0.0001 CONSOLIDATE TO SEDVRM INSTANTANEOUSLY  
        IF( SEDVRDT >= 0.0 .AND. SEDVRDT <= 0.0001 )THEN  
          TMPEXP = 0.0  
          DO K=1,KB  
            DO L=LF,LL  
              IF( K <= KBT(L) )THEN  
                VDRBED1(L,K) = VDRBED(L,K)  
                HBED1(L,K) = HBED(L,K)  
                VDRBED(L,K) = SEDVDRM  
                TMPTOP = 1. + VDRBED(L,K)  
                TMPBOT = 1. + VDRBED1(L,K)  
                HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT  
              ENDIF  
            ENDDO  
          ENDDO  
        ENDIF  

        ! ***  IF SEDVRDT<0.0 MAINTAIN INITIAL VOID RATIO (SAVED IN VDRBED2)  
        IF( SEDVRDT < 0.0 )THEN  
          TMPEXP = 1.0  
          DO L=LF,LL  
            K = KBT(L)  
            VDRBED1(L,K) = VDRBED(L,K)  
            HBED1(L,K) = HBED(L,K)  
            VDRBED(L,K) = VDRBED2(L,K)  
            TMPTOP = 1. + VDRBED(L,K)  
            TMPBOT = 1. + VDRBED1(L,K)  
            HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT  
          ENDDO  
          DO L=LF,LL  
            K = KBT(L)-1  
            IF( K > 0 )THEN  
              VDRBED1(L,K) = VDRBED(L,K)  
              HBED1(L,K) = HBED(L,K)  
              VDRBED(L,K) = VDRBED2(L,K)  
              TMPTOP = 1. + VDRBED(L,K)  
              TMPBOT = 1. + VDRBED1(L,K)  
              HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT  
            ENDIF  
          ENDDO  
        ENDIF  

        ! *** UPDATE POROSITY  
        DO K=1,KB  
          DO L=LF,LL  
            IF( K <= KBT(L) )THEN  
              PORBED(L,K)  = VDRBED(L,K) /(1. + VDRBED(L,K))  
              PORBED1(L,K) = VDRBED1(L,K)/(1. + VDRBED1(L,K))  
            ENDIF  
          ENDDO  
        ENDDO  

        ! *** UPDATE PORE WATER FLOWS  
        DO L=LF,LL  
          QWTRBED(L,0) = QGW(L)/DXYP(L)  
        ENDDO  
        DO K=1,KB  
          DO L=LF,LL  
            IF( K <= KBT(L) )THEN  
              TMPVAL = HBED(L,K)/(1. + VDRBED(L,K))  
              QWTRBED(L,K) = QWTRBED(L,K-1)  -DELTI*TMPVAL*(VDRBED(L,K)-VDRBED1(L,K))  
            ENDIF  
          ENDDO  
        ENDDO  

      ENDDO   ! *** END OF DOMAIN LOOP
      !$OMP END PARALLEL DO

    ENDIF  

    ! *********************************************************************
    ! *** CONSOLIDATE TO VOID RATIO/POROSITY DEFINED AS FUNCTION OF
    ! *** COHESIVE SEDIMENT FRACTION, FRACCOH(L,K)
    IF( IBMECH == 99 )THEN

      !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,L,LF,LL)              &
      !$OMP             PRIVATE(TMPEXP,TMPTOP,TMPBOT,TMPVAL,POROTMP)
      DO ND=1,NDM  
        LF = 2+(ND-1)*LDM  
        LL = MIN(LF+LDM-1,LA)
        
        DO K=1,KB
          DO L=LF,LL
            IF( K <= KBT(L) )THEN
              VDRBED1(L,K) = VDRBED(L,K)
              HBED1(L,K) = HBED(L,K)
              POROTMP = BMECH1*(FRACCOH(L,K)**BMECH2) + BMECH3  ! *** PMC - THIS APPROACH DISCONNECTS THE SEDIMENT MASS FROM THE VOID RATIO
              VDRBED(L,K) = POROTMP/(1.-POROTMP)
              TMPTOP = 1. + VDRBED(L,K)
              TMPBOT = 1. + VDRBED1(L,K)
              HBED(L,K) = TMPTOP*HBED1(L,K)/TMPBOT
            ENDIF
          ENDDO
        ENDDO

        ! *** UPDATE POROSITY
        DO K=1,KB
          DO L=LF,LL
            IF( K <= KBT(L) )THEN
              PORBED(L,K) = VDRBED(L,K)/(1. + VDRBED(L,K))
              PORBED1(L,K) = VDRBED1(L,K)/(1. + VDRBED1(L,K))
            ENDIF
          ENDDO
        ENDDO

        ! *** UPDATE PORE WATER FLOWS
        DO L=LF,LL
          QWTRBED(L,0) = QGW(L)/DXYP(L)
        ENDDO

        DO K=1,KB
          DO L=LF,LL
            IF( K <= KBT(L) )THEN
              TMPVAL = HBED(L,K)/(1. + VDRBED(L,K))
              QWTRBED(L,K) = QWTRBED(L,K-1)-DELTI*TMPVAL*(VDRBED(L,K)-VDRBED1(L,K))
            ENDIF
          ENDDO
        ENDDO

      ENDDO   ! *** END OF DOMAIN LOOP
      !$OMP END PARALLEL DO

    ENDIF     ! *** END OF IBMECH = 99


    ! *********************************************************************
    ! *** SET FLAG FOR FINITE STRAIN CONSOLIDATION  
    IFLAG = 0  
    IF( ISTRAN(6) >= 1 )IFLAG = IFLAG + 1  
    IF( ISTRAN(7) >= 1 )IFLAG = IFLAG + 1  

    ! *** FINITE STRAIN CONSOLIDATING HOMOGENEOUS BED  
    IF( IBMECH >= 2 .AND. IBMECH < 98 .AND. IFLAG == 1 )THEN
      WDENKGM3 = 1.E3  
      WDENGMM3 = 1.E6  

      ! ***  SET PHYSICAL VERTICAL COORDINATES OF THE BED  

      ! ***  ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
      ! ***  ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
      DO L=2,LA  
        ZBEDG(L,0) = ZELBEDA(L)  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDG(L,K) = ZBEDG(L,K-1) + HBED(L,K)  
          ELSE  
            ZBEDG(L,K) = HBED(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        ZBEDGT(L) = ZBEDG(L,KBT(L))  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDC(L,K) = 0.5*(ZBEDG(L,K) + ZBEDG(L,K-1))  
          ELSE  
            ZBEDC(L,K) = 0.5*ZBEDG(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  CALCULATE TRANSFORMED THICKNESS OF BED LAYERS  
      ! ***  DZBTR = TRANSFORMED THICKNESS OF BED LAYER
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            DZBTR(L,K) = HBED(L,K)/(1. + VDRBED(L,K))  
            DZBTR1(L,K) = HBED1(L,K)/(1. + VDRBED1(L,K))  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  CALCULATE WATER SPECIFIC WEIGHT NORMALIZED  
      ! ***    EFFECTIVE STRESS USING FSTRSE  
      ! ***  CALCULATE DERIVATIVE OF EFFECTIVE STRESS  
      ! ***    WITH RESPECT TO VOID RATIO, DSTRSE USING  
      ! ***    FUNCTION FDSTRSE  
      ! ***  CALCULATE HYDRAULIC CONDUCTIVITY DIVIDED BY (1+VOID),  
      ! ***    HYDCN =USING FUNCTION FHYDCN 
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            STRSE(L,K)=  FSTRSE(VDRBED(L,K),BMECH1,BMECH2,BMECH3)  
            DSTRSE(L,K)=  FDSTRSE(VDRBED(L,K),BMECH1,BMECH2,BMECH3)  
            HYDCN(L,K)=  FHYDCN(VDRBED(L,K),BMECH4,BMECH5,BMECH6,IBMECHK)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          KBTM1 = KBT(L)-1  
          IF( K <= KBTM1 )THEN  
            TMPVAL = ( DZBTR(L,K)/HYDCN(L,K) )  +( DZBTR(L,K+1)/HYDCN(L,K+1) )  
            COEFK(L,K) = (DZBTR(L,K) + DZBTR(L,K+1))/TMPVAL  
            DSTRESET = (DZBTR(L,K)*DSTRSE(L,K+1)  +DZBTR(L,K+1)*DSTRSE(L,K))  /(DZBTR(L,K) + DZBTR(L,K+1))  
            COEFSK(L,K) = DSTRESET*COEFK(L,K)  
          ENDIF  
          IF( K == KBT(L) )THEN  
            COEFK(L,K) = HYDCN(L,K)  
            COEFSK(L,K) = DSTRSE(L,K)*HYDCN(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  CALCULATE PRESSURE COMPONENTS  
      DO K=1,KB  
        DO L=2,LA  
          SGSM1(L,K) = 0.  
        ENDDO  
      ENDDO  
      IF( ISTRAN(6) > 0 )THEN  
        DO NS=1,NSED  
          DO K=1,KB  
            DO L=2,LA  
              IF( K <= KBT(L) )THEN  
                SGSM1(L,K) = SGSM1(L,K) + SSG(NS)*VFRBED(L,K,NS)  
              ENDIF  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(7) > 0 )THEN  
        DO NX=1,NSND  
          NS = NSED + NX  
          DO K=1,KB  
            DO L=2,LA  
              IF( K <= KBT(L) )THEN  
                SGSM1(L,K) = SGSM1(L,K) + SSG(NS)*VFRBED(L,K,NS)  
              ENDIF  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            SGSM1(L,K) = SGSM1(L,K)-1.  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  NEW IMPLICIT CONSOLIDATION SOLUTION BEGINS HERE  
      DO K=1,KB  
        DO L=2,LA  
          IF( K<KBT(L) )THEN  
            ACOEF(L,K) = 2.*COEFSK(L,K)/(DZBTR(L,K) + DZBTR(L,K+1))  
            TMPVALK = DZBTR(L,K)*SGSM1(L,K)  
            TMPVALKP = DZBTR(L,K+1)*SGSM1(L,K+1)  
            QCOEF(L,K) = (TMPVALK + TMPVALKP)*COEFK(L,K)/  (DZBTR(L,K) + DZBTR(L,K+1))  
          ELSE  
            ACOEF(L,K) = 0.0  
            QCOEF(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        K = KBT(L)  
        ACOEF(L,K) = 2.*COEFSK(L,K)/DZBTR(L,K)  
        QCOEF(L,K) = SGSM1(L,K)*COEFK(L,K)  
        QCOEF(L,0) = QGW(L)*DXYIP(L)  
        QWTRBED(L,0) = QGW(L)*DXYIP(L)  
      ENDDO  
      DO L=2,LA  
        ALOW(L,1) = 0.  
        CUPP(L,KBT(L)) = 0.  
        DO K=1,KBT(L)-1  
          CUPP(L,K) = -DTSED*ACOEF(L,K)/DZBTR(L,K)  
        ENDDO  
        DO K=2,KBT(L)  
          ALOW(L,K) = -DTSED*ACOEF(L,K-1)/DZBTR(L,K)  
        ENDDO  
        DO K=1,KBT(L)  
          BMNN(L,K) = 1.0-ALOW(L,K)-CUPP(L,K)  
        ENDDO  
        K = KBT(L)  
        BMNN(L,K) = BMNN(L,K) + DTSED*ACOEF(L,K)/DZBTR(L,K)  
        DO K=1,KBT(L)  
          RRHS(L,K) = VDRBED(L,K)  +DTSED*(QCOEF(L,K-1)-QCOEF(L,K))/DZBTR(L,K)  
          VDRBED1(L,K) = VDRBED(L,K)  
          HBED1(L,K) = HBED(L,K)  
        ENDDO  
      ENDDO  
      
      IF( IBMECH == 2 )THEN  
        DO L=2,LA  
          K = KBT(L)  
          RRHS(L,K) = RRHS(L,K) + DTSED*ACOEF(L,K)*VDRDEPO(1)/DZBTR(L,K)  
        ENDDO  
      ENDIF  
      IF( IBMECH == 3 )THEN  
        DO L=2,LA  
          K = KBT(L)  
          RRHS(L,K) = RRHS(L,K) + DTSED*ACOEF(L,K)*(VDRBED(L,K)  +(STRSE(L,K)/DSTRSE(L,K)))/DZBTR(L,K)  
        ENDDO  
      ENDIF  
      DO L=2,LA  
        BETTMP = BMNN(L,1)  
        TOXTMP(L,1) = RRHS(L,1)/BETTMP  
        DO KK=2,KBT(L)  
          GAMTMP(L,KK) = CUPP(L,KK-1)/BETTMP  
          BETTMP = BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)  
          TOXTMP(L,KK) = (RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/  BETTMP  
        ENDDO  
        DO KK=KBT(L)-1,1,-1  
          TOXTMP(L,KK) = TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K<KBT(L) )THEN  
            QWTRBED(L,K) = -ACOEF(L,K)*(TOXTMP(L,K+1)-TOXTMP(L,K))  +QCOEF(L,K)  
          ELSE  
            QWTRBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  
      IF( IBMECH == 2 )THEN  
        DO L=2,LA  
          K = KBT(L)  
          QWTRBED(L,K) = -ACOEF(L,K)*(SEDVDRD-TOXTMP(L,K)) + QCOEF(L,K)  
        ENDDO  
      ELSE  
        DO L=2,LA  
          K = KBT(L)  
          QWTRBED(L,K) = -ACOEF(L,K)*(VDRBED(L,K)  +(STRSE(L,K)/DSTRSE(L,K))-TOXTMP(L,K)) + QCOEF(L,K)  
        ENDDO  
      END IF  

      ! ***  CALCULATE VOID RATIOS  
      ! ***  VDRBED =  VOID RATIO OF BED LAYER 
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            VDRBED(L,K) = VDRBED1(L,K)  -DTSED*(QWTRBED(L,K)-QWTRBED(L,K-1))/DZBTR1(L,K)  
          ELSE  
            VDRBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  UPDATE LAYER THICKNESS  
      ! ***  HBED = BED LAYER THICKNESS
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            HBED(L,K) = HBED1(L,K)*(1. + VDRBED(L,K))/(1. + VDRBED1(L,K))  
          ELSE  
            HBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        QWTRBED(L,0) = QGW(L)/DXYP(L)  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            QWTRBED(L,K) = QWTRBED(L,K-1)  -DELTI*(HBED(L,K)-HBED1(L,K))  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
      ! ***  ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDG(L,K) = ZBEDG(L,K-1) + HBED(L,K)  
          ELSE  
            ZBEDG(L,K) = HBED(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        ZBEDGT(L) = ZBEDG(L,KBT(L))  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDC(L,K) = 0.5*(ZBEDG(L,K) + ZBEDG(L,K-1))  
          ELSE  
            ZBEDC(L,K) = 0.5*ZBEDG(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            DZBTR(L,K) = HBED(L,K)/(1. + VDRBED(L,K))  
          ENDIF  
        ENDDO  
      ENDDO  

      ! *** UPDATE POROSITY  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            PORBED(L,K) = VDRBED(L,K)/(1. + VDRBED(L,K))  
            PORBED1(L,K) = VDRBED1(L,K)/(1. + VDRBED1(L,K))  
          ELSE  
            PORBED(L,K) = 0.0  
            PORBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  

      ! *** ADD OR REMOVE LAYERS  
      !CALL CALBLAY  

    ENDIF  

    ! *********************************************************************
    ! *** FINITE STRAIN CONSOLIDATING NON-HOMOGENEOUS BED  
    IF( IBMECH >= 2 .AND. IBMECH<98 .AND. IFLAG == 2 )THEN
      WDENKGM3 = 1.E3  
      WDENGMM3 = 1.E6  

      ! ***  SET PHYSICAL VERTICAL COORDINATES OF THE BED  
        
      ! ***  ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
      ! ***  ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
      DO L=2,LA  
        ZBEDG(L,0) = ZELBEDA(L)  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDG(L,K) = ZBEDG(L,K-1) + HBED(L,K)  
          ELSE  
            ZBEDG(L,K) = HBED(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        ZBEDGT(L) = ZBEDG(L,KBT(L))  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDC(L,K) = 0.5*(ZBEDG(L,K) + ZBEDG(L,K-1))  
          ELSE  
            ZBEDC(L,K) = 0.5*ZBEDG(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  CALCULATE TRANSFORMED THICKNESS OF BED LAYERS  
      ! ***  DZBTR = TRANSFORMED THICKNESS OF BED LAYER
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            DZBTR(L,K) = HBED(L,K)/(1. + VDRBED(L,K))  
            DZBTR1(L,K) = HBED1(L,K)/(1. + VDRBED1(L,K))  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  CALCULATE WATER SPECIFIC WEIGHT NORMALIZED  
      ! ***    EFFECTIVE STRESS USING FSTRSE  
      ! ***  CALCULATE DERIVATIVE OF EFFECTIVE STRESS  
      ! ***    WITH RESPECT TO VOID RATIO, DSTRSE USING  
      ! ***    FUNCTION FDSTRSE  
      ! ***  CALCULATE HYDRAULIC CONDUCTIVITY DIVIED BY (1+VOID),  
      ! ***    HYDCN =USING FUNCTION FHYDCN 
      ! ***  NONCOHESIVE HYDRAULIC CONDUCTIVITY BASED ON R. R. RUMER, CHAP 3,  
      ! ***  EQ 13 AND 14, IN 'FLOW THROUGH POROUS MEDIA' ED. R. J. M. DE WIEST  
      ! ***  ACADEMIC PRESS, 1969  
      ! ***  KH = 2854*(n**2)*(d**2)/(1-n) = 2854*(e**2)*(d**2)/(1+e)
      ! ***  FOR KH IN M/S AND D IN METERS  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            STRSE(L,K)=  FSTRSE(VDRBEDSED(L,K),BMECH1,BMECH2,BMECH3)  
            DSTRSE(L,K)=  FDSTRSE(VDRBEDSED(L,K),BMECH1,BMECH2,BMECH3)  
            HYDCN(L,K)=  FHYDCN(VDRBEDSED(L,K),BMECH4,BMECH5,BMECH6,IBMECHK)  
            TMPVAL = VDRBEDSND(L,K)/(1.+VDRBEDSND(L,K))  
            SNDHYDCN(L,K) = 2854.0*TMPVAL*TMPVAL*(SEDDIA50(L,K)**2)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            COEFK(L,K) = (FRACCOH(L,K)/HYDCN(L,K))  +(FRACNON(L,K)/SNDHYDCN(L,K))  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            HYDCN(L,K) = 1./COEFK(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          COEFK(L,K) = 0.  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          KBTM1 = KBT(L)-1  
          IF( K <= KBTM1 )THEN  
            TMPVAL = ( DZBTR(L,K)/HYDCN(L,K) )  +( DZBTR(L,K+1)/HYDCN(L,K+1) )  
            COEFK(L,K) = (DZBTR(L,K) + DZBTR(L,K+1))/TMPVAL  
            DSTRESET = (DZBTR(L,K)*DSTRSE(L,K+1)  +DZBTR(L,K+1)*DSTRSE(L,K))  /(DZBTR(L,K) + DZBTR(L,K+1))  
            COEFSK(L,K) = DSTRESET*COEFK(L,K)  
          ENDIF  
          IF( K == KBT(L) )THEN  
            COEFK(L,K) = HYDCN(L,K)  
            COEFSK(L,K) = DSTRSE(L,K)*HYDCN(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  CALCULATE PRESSURE COMPONENTS  
      DO K=1,KB  
        DO L=2,LA  
          SGSM1(L,K) = 0.  
        ENDDO  
      ENDDO  
      IF( ISTRAN(6) > 0 )THEN  
        DO NS=1,NSED  
          DO K=1,KB  
            DO L=2,LA  
              IF( K <= KBT(L) )THEN  
                SGSM1(L,K) = SGSM1(L,K) + SSG(NS)*VFRBED(L,K,NS)  
              ENDIF  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(7) > 0 )THEN  
        DO NX=1,NSND  
          NS = NSED + NX  
          DO K=1,KB  
            DO L=2,LA  
              IF( K <= KBT(L) )THEN  
                SGSM1(L,K) = SGSM1(L,K) + SSG(NS)*VFRBED(L,K,NS)  
              ENDIF  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            SGSM1(L,K) = SGSM1(L,K)-1.  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  NEW IMPLICIT CONSOLIDATION SOLUTION BEGINS HERE  
      DO K=1,KB  
        DO L=2,LA  
          IF( K<KBT(L) )THEN  
            ACOEF(L,K) = 2.*COEFSK(L,K)/(DZBTR(L,K) + DZBTR(L,K+1))  
            TMPVALK = DZBTR(L,K)*SGSM1(L,K)  
            TMPVALKP = DZBTR(L,K+1)*SGSM1(L,K+1)  
            QCOEF(L,K) = (TMPVALK + TMPVALKP)*COEFK(L,K)/  (DZBTR(L,K) + DZBTR(L,K+1))  
          ELSE  
            ACOEF(L,K) = 0.0  
            QCOEF(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        K = KBT(L)  
        ACOEF(L,K) = 2.*COEFSK(L,K)/DZBTR(L,K)  
        QCOEF(L,K) = SGSM1(L,K)*COEFK(L,K)  
        QCOEF(L,0) = QGW(L)*DXYIP(L)  
        QWTRBED(L,0) = QGW(L)*DXYIP(L)  
      ENDDO  
      DO L=2,LA  
        ALOW(L,1) = 0.  
        CUPP(L,KBT(L)) = 0.  
        DO K=1,KBT(L)-1  
          CUPP(L,K) = -DTSED*ACOEF(L,K)/DZBTR(L,K)  
        ENDDO  
        DO K=2,KBT(L)  
          ALOW(L,K) = -DTSED*ACOEF(L,K-1)/DZBTR(L,K)  
        ENDDO  
        DO K=1,KBT(L)  
          BMNN(L,K) = FRACCOH(L,K)-ALOW(L,K)-CUPP(L,K)  
        ENDDO  
        K = KBT(L)  
        BMNN(L,K) = BMNN(L,K) + DTSED*ACOEF(L,K)/DZBTR(L,K)  
        DO K=1,KBT(L)  
          RRHS(L,K) = FRACCOH(L,K)*VDRBEDSED(L,K)  +DTSED*(QCOEF(L,K-1)-QCOEF(L,K))/DZBTR(L,K)  
          VDRBED1(L,K) = VDRBED(L,K)  
          HBED1(L,K) = HBED(L,K)  
        ENDDO  
      ENDDO  
      IF( IBMECH == 2 )THEN  
        DO L=2,LA  
          K = KBT(L)  
          RRHS(L,K) = RRHS(L,K) + DTSED*ACOEF(L,K)*VDRDEPO(1)/DZBTR(L,K)  
        ENDDO  
      ENDIF  
      IF( IBMECH == 3 )THEN  
        DO L=2,LA  
          K = KBT(L)  
          RRHS(L,K) = RRHS(L,K)  +DTSED*ACOEF(L,K)*(VDRBEDSED(L,K)  +(STRSE(L,K)/DSTRSE(L,K)))/DZBTR(L,K)  
        ENDDO  
      ENDIF  
      DO L=2,LA  
        BETTMP = BMNN(L,1)  
        TOXTMP(L,1) = RRHS(L,1)/BETTMP  
        DO KK=2,KBT(L)  
          GAMTMP(L,KK) = CUPP(L,KK-1)/BETTMP  
          BETTMP = BMNN(L,KK)-ALOW(L,KK)*GAMTMP(L,KK)  
          TOXTMP(L,KK) = (RRHS(L,KK)-ALOW(L,KK)*TOXTMP(L,KK-1))/  BETTMP  
        ENDDO  
        DO KK=KBT(L)-1,1,-1  
          TOXTMP(L,KK) = TOXTMP(L,KK)-GAMTMP(L,KK+1)*TOXTMP(L,KK+1)  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K<KBT(L) )THEN  
            QWTRBED(L,K) = -ACOEF(L,K)*(TOXTMP(L,K+1)-TOXTMP(L,K))  +QCOEF(L,K)  
          ELSE  
            QWTRBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  
      IF( IBMECH == 2 )THEN  
        DO L=2,LA  
          K = KBT(L)  
          QWTRBED(L,K) = -ACOEF(L,K)*(VDRDEPO(1) - TOXTMP(L,K)) + QCOEF(L,K)  
        ENDDO  
      ELSE  
        DO L=2,LA  
          K = KBT(L)  
          QWTRBED(L,K) = -ACOEF(L,K)*(VDRBEDSED(L,K) + (STRSE(L,K)/DSTRSE(L,K))-TOXTMP(L,K)) + QCOEF(L,K)  
        ENDDO  
      END IF  

      ! ***  CALCULATE VOID RATIOS  
      ! ***  VDRBED =  VOID RATIO OF BED LAYER 
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            IF( FRACCOH(L,K)>0.0 )THEN
            VDRBEDSED(L,K) = VDRBEDSED(L,K) - DTSED*(QWTRBED(L,K)-QWTRBED(L,K-1))/(FRACCOH(L,K)*DZBTR1(L,K))  
              VDRBEDSED(L,K) = MAX(VDRBEDSED(L,K),VDRBEDSND(L,K))
            ELSE
              VDRBEDSED(L,K) = 0.0
            ENDIF
          ELSE  
            VDRBEDSED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  UPDATE VOID RATIO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            VDRBED(L,K) = FRACCOH(L,K)*VDRBEDSED(L,K) + FRACNON(L,K)*VDRBEDSND(L,K)  
          ELSE  
            VDRBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  UPDATE LAYER THICKNESS  
      ! ***  HBED = BED LAYER THICKNESS
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            HBED(L,K) = HBED1(L,K)*(1. + VDRBED(L,K))/(1. + VDRBED1(L,K))  
          ELSE  
            HBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        QWTRBED(L,0) = QGW(L)/DXYP(L)  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            QWTRBED(L,K) = QWTRBED(L,K-1)  -DELTI*(HBED(L,K)-HBED1(L,K))  
          ENDIF  
        ENDDO  
      ENDDO  

      ! ***  ZBEDC = VERTICAL COORDINATE OF THE CENTER OF BED LAYER
      ! ***  ZBEDG = VERTICAL COORDINATE AT TOP OF BED LAYER
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDG(L,K) = ZBEDG(L,K-1) + HBED(L,K)  
          ELSE  
            ZBEDG(L,K) = HBED(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO L=2,LA  
        ZBEDGT(L) = ZBEDG(L,KBT(L))  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            ZBEDC(L,K) = 0.5*(ZBEDG(L,K) + ZBEDG(L,K-1))  
          ELSE  
            ZBEDC(L,K) = 0.5*ZBEDG(L,K)  
          ENDIF  
        ENDDO  
      ENDDO  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            DZBTR(L,K) = HBED(L,K)/(1. + VDRBED(L,K))  
          ENDIF  
        ENDDO  
      ENDDO  

      ! *** UPDATE POROSITY  
      DO K=1,KB  
        DO L=2,LA  
          IF( K <= KBT(L) )THEN  
            PORBED(L,K) = VDRBED(L,K)/(1. + VDRBED(L,K))  
            PORBED1(L,K) = VDRBED1(L,K)/(1. + VDRBED1(L,K))  
          ELSE  
            PORBED(L,K) = 0.0  
            PORBED(L,K) = 0.0  
          ENDIF  
        ENDDO  
      ENDDO  

      ! *** ADD OR REMOVE LAYERS  
      !CALL CALBLAY  
    ENDIF   ! *** END OF FINITE STRAIN OPION
    
    ! *********************************************************************
    ! *** UPDATE BULK DENSITY AND TOTAL SEDIMENT MASS
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,K,L,NS) 
    DO ND=1,NDM  
      LF = 2+(ND-1)*LDM  
      LL = MIN(LF+LDM-1,LA)

      DO K=1,KB
        DO L=LF,LL
          SEDBT(L,K) = 0.0
          SNDBT(L,K) = 0.0
        ENDDO
      ENDDO

      IF( ISTRAN(6) >= 1 )THEN
        DO NS=1,NSED
          DO K=1,KB
            DO L=LF,LL
              SEDBT(L,K) = SEDBT(L,K) + SEDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      IF( ISTRAN(7) >= 1 )THEN
        DO NS=1,NSND
          DO K=1,KB
            DO L=LF,LL
              SNDBT(L,K) = SNDBT(L,K) + SNDB(L,K,NS)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      
      DO K=1,KB
        DO L=LF,LL
          IF( HBED(L,K) > 0. )THEN
            ! *** UPDATE TOTAL/WET DENSITY
            BDENBED(L,K) = 1000.*PORBED(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)
          ELSE
            BDENBED(L,K) = 0.
          ENDIF
        ENDDO
      ENDDO
    ENDDO  ! END OF DOMAIN
    !$OMP END PARALLEL DO

  ENDIF

  8011 FORMAT('CBED',I12,6E13.5)    

  !**********************************************************************!

  RETURN  
  END  

