SUBROUTINE CALQQ2T()

  ! **  SUBROUTINE CALQQ CALCULATES THE TURBULENT INTENSITY SQUARED AT                   
  ! **  TIME LEVEL (N+1).  THE VALUE OF ISTL INDICATES THE NUMBER OF                     
  ! **  TIME LEVELS INVOLVED                            

  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP
  !    2010-10       Scott James       Added MHK

  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER :: L,K,LS,LN,LE,LW,LF,ND,LL,LP
  
  REAL    :: BETAVEG_P=1.0, BETAVEG_D=5.1,CE4VEG=0.9 !from Katul et al. 2003
  REAL    :: BETASUP_P=1.0, BETASUP_D=5.1,CE4SUP=0.9 !from Katul et al. 2003

  REAL :: TMPQQI,TMPQQE,BSMALL,WB,UHUW,VHVW,DELB,CTE3TMP,PQQB,PQQU,PQQ,PQQL,PQQW,PQQV
  REAL :: CLQTMP,CUQTMP,CLQLTMP,CUQLTMP,CMQTMP,CMQLTMP,EQ,EQL,QQHDH,DMLTMP,DMLMAX,FFTMP
  REAL :: TMPVAL,WVFACT,SQLDSQ

  REAL,SAVE,ALLOCATABLE :: PQQVEGI(:,:),PQQVEGE(:,:)
  REAL,SAVE,ALLOCATABLE :: PQQMHKI(:,:),PQQMHKE(:,:)
  REAL,SAVE,ALLOCATABLE :: PQQSUPI(:,:),PQQSUPE(:,:)

  IF(  .NOT. ALLOCATED(PQQVEGI) )THEN
    ALLOCATE(PQQVEGI(LCM,KCM))
    ALLOCATE(PQQVEGE(LCM,KCM))
    ALLOCATE(PQQMHKI(LCM,KCM))
    ALLOCATE(PQQMHKE(LCM,KCM))
    ALLOCATE(PQQSUPI(LCM,KCM))
    ALLOCATE(PQQSUPE(LCM,KCM))
    PQQVEGI=0.0
    PQQVEGE=0.0
    PQQMHKI=0.0
    PQQMHKE=0.0
    PQQSUPI=0.0
    PQQSUPE=0.0
  ENDIF

  IF( ISDYNSTP == 0 )THEN
    DELT=DT
    DELTI=1./DELT
  ELSE
    DELT=DTDYN
    DELTI=1./DELT
  END IF
  BSMALL=1.E-12
  
  ! *** SET WAVE RAMPUP FACTOR
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN
    IF( N < NTSWV )THEN
      TMPVAL = FLOAT(N)/FLOAT(NTSWV)
      WVFACT = 0.5-0.5*COS(PI*TMPVAL)
    ELSE
      WVFACT=1.0
    ENDIF
  ENDIF
  
  ! *** ZERO FOR INITIALLY DRY CELLS
  IF( LADRY > 0 )THEN
    DO K=1,KC
      DO LP=1,LADRY
        L=LDRY(LP)  
        LE=LEC(L)
        LN=LNC(L)
        FWQQ(L,K)=0.
        FWQQL(L,K)=0.
        FUHU(L,K)=0.
        FUHV(L,K)=0.
        FUHU(LE,K)=0.
        FUHV(LE,K)=0.
        FVHU(L,K)=0.
        FVHV(L,K)=0.
        FVHU(LN,K)=0.
        FVHV(LN,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        PQQVEGI(L,K)=0.
        PQQVEGE(L,K)=0.
        PQQMHKI(L,K)=0.
        PQQMHKE(L,K)=0.
        PQQSUPI(L,K)=0.
        PQQSUPE(L,K)=0.
        CU1(L,K)=0.
        CU2(L,K)=0.
        TVAR1W(L,K)=0.
      ENDDO  
    ENDDO 
  ENDIF
  
  ! *** SET RATIO OF LENTH SCALE*TURB_INTENSITY TO TURB_INTENSITY DIFFUSION
  SQLDSQ=1.0
  IF( ISTOPT(0) == 3 )SQLDSQ=0.377/0.628

  !$OMP PARALLEL DEFAULT(SHARED)

  ! *** ZERO ACCUMULATION ARRAYS FOR ACTIVE CELLS
  !$OMP DO PRIVATE(ND,K,LP,L)
  DO ND=1,NDM  

    DO K=1,KS
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        QQ2(L,K)  = QQ(L,K)  + QQ(L,K)
        QQL2(L,K) = QQL(L,K) + QQL(L,K)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO
  
  ! **  CALCULATE ADVECTIVE FLUXES BY UPWIND DIFFERENCE WITH TRANSPORT                   
  ! **  AVERAGED BETWEEN (N) AND (N+1) AND TRANSPORTED FIELD AT (N) OR                   
  ! **  TRANSPORT BETWEEN (N-1) AND (N+1) AND TRANSPORTED FIELD AT (N-1)                 
  ! **  FOR ISTL EQUAL TO 2 AND 3 RESPECTIVELY          

  ! *** VERTICAL FLUXES
  !$OMP DO PRIVATE(ND,K,LP,L,WB)
  DO ND=1,NDM  

    DO K=1,KC
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        WB = 0.5*DXYP(L)*(W2(L,K-1)+W2(L,K))
        FWQQ(L,K)  = MAX(WB,0.)*QQ(L,K-1)        + MIN(WB,0.)*QQ(L,K)
        FWQQL(L,K) = MAX(WB,0.)*QQL(L,K-1)*HP(L) + MIN(WB,0.)*QQL(L,K)*HP(L)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! *** HORIZONTAL FLUXES
  !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,UHUW,VHVW)
  DO ND=1,NDM  
    DO K=1,KS
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        LS=LSC(L)
        LW=LWC(L)
        UHUW = 0.5*(UHDYF2(L,K)+UHDYF2(L,K+1))
        FUHU(L,K) = MAX(UHUW,0.)*QQ(LW,K)         + MIN(UHUW,0.)*QQ(L,K)
        FUHV(L,K) = MAX(UHUW,0.)*QQL(LW,K)*HP(LW) + MIN(UHUW,0.)*QQL(L,K)*HP(L)

        VHVW = 0.5*(VHDXF2(L,K)+VHDXF2(L,K+1))
        FVHU(L,K) = MAX(VHVW,0.)*QQ(LS,K)         + MIN(VHVW,0.)*QQ(L,K)
        FVHV(L,K) = MAX(VHVW,0.)*QQL(LS,K)*HP(LS) + MIN(VHVW,0.)*QQL(L,K)*HP(L)
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! *** APPLY LAYER SPECIFIC SUB/SVB TO FLUX TERMS
  IF( IGRIDV > 0 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KS
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          FUHU(L,K) = SUB3D(L,K)*FUHU(L,K)
          FUHV(L,K) = SUB3D(L,K)*FUHV(L,K)
          FVHU(L,K) = SVB3D(L,K)*FVHU(L,K)
          FVHV(L,K) = SVB3D(L,K)*FVHV(L,K)
        ENDDO
      ENDDO
    ENDDO  ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  ! **  CALCULATE PRODUCTION, LOAD BOUNDARY CONDITIONS AND SOLVE                         
  ! **  TRANSPORT EQUATIONS                             
  ! **  FUHQQ=FUHU, FVHQQ=FVHU, FUHQQL=FUHV, FVHQQL=FVHV
  !$OMP SINGLE
  DO K=1,KC
    DO LL=1,NPBS
      L=LPBS(LL)
      LN=LNC(L)
      IF( FVHU(LN,K) > 0. )THEN
        FVHU(LN,K)=0.0
        FVHV(LN,K)=0.0
      ENDIF
    ENDDO
    DO LL=1,NPBW
      L=LPBW(LL)
      IF( FUHU(LEC(L),K) > 0. )THEN
        FUHU(LEC(L),K)=0.0
        FUHV(LEC(L),K)=0.0
      ENDIF
    ENDDO
    DO LL=1,NPBE
      L=LPBE(LL)
      IF( FUHU(L,K) < 0. )THEN
        FUHU(L,K)=0.0
        FUHV(L,K)=0.0
      ENDIF
    ENDDO
    DO LL=1,NPBN
      L=LPBN(LL)
      IF( FVHU(L,K) < 0. )THEN
        FVHU(L,K)=0.0
        FVHV(L,K)=0.0
      ENDIF
    ENDDO
  ENDDO
  !$OMP END SINGLE

  ! *** ADD VEGETATION IMPACTS ON TURBULENCE
  IF( ISVEG > 0 )THEN !SCJ vegetative/MHK impact on K-epsilon
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,TMPQQI,TMPQQE)
    DO ND=1,NDM  
      DO K=1,KS
        DO LP=1,LLVEG(K,ND)
          L=LKVEG(LP,K,ND)  
          LE=LEC(L)
          LN=LNC(L)
          TMPQQI = 0.25*BETAVEG_P
          TMPQQE = 0.25*BETAVEG_D
          PQQVEGI(L,K) = TMPQQI*( FXVEG(L,K)+FXVEG(L,K+1)+FXVEG(LE,K)+FXVEG(LE,K+1) + FYVEG(L,K)+FYVEG(L,K+1)+FYVEG(LN,K)+FYVEG(LN,K+1) )
          PQQVEGE(L,K) = TMPQQE*( FXVEG(L ,K)*U(L ,K)*U(L ,K) + FXVEG(L ,K+1)*U(L ,K+1)*U(L ,K+1) &
                                 +FXVEG(LE,K)*U(LE,K)*U(LE,K) + FXVEG(LE,K+1)*U(LE,K+1)*U(LE,K+1) &
                                 +FYVEG(L ,K)*V(L ,K)*V(L ,K) + FYVEG(L ,K+1)*V(L ,K+1)*V(L ,K+1) &
                                 +FYVEG(LN,K)*V(LN,K)*V(LN,K) + FYVEG(LN,K+1)*V(LN,K+1)*V(LN,K+1))
          IF( MVEGL(L)>90 )THEN
            TMPQQI = 0.5*BETAMHK_P
            TMPQQE = 0.5*BETAMHK_D
            PQQMHKI(L,K) = TMPQQI*(FXMHK(L,K)+FXMHK(L,K+1)+FYMHK(L,K)+FYMHK(L,K+1))
            PQQMHKE(L,K) = TMPQQE*(FXMHK(L,K)*U(L,K)*U(L,K) + FXMHK(L,K+1)*U(L,K+1)*U(L,K+1) &
                                  +FYMHK(L,K)*V(L,K)*V(L,K) + FYMHK(L,K+1)*V(L,K+1)*V(L,K+1))
            TMPQQI = 0.5*BETASUP_P
            TMPQQE = 0.5*BETASUP_D
            PQQSUPI(L,K) = TMPQQI*(FXSUP(L,K)+FXSUP(L,K+1)+FYSUP(L,K  )+FYSUP(L,K+1))
            PQQSUPE(L,K) = TMPQQE*(FXSUP(L,K)*U(L,K)*U(L,K) + FXSUP(L,K+1)*U(L,K+1)*U(L,K+1) &
                                  +FYSUP(L,K)*V(L,K)*V(L,K) + FYSUP(L,K+1)*V(L,K+1)*V(L,K+1))
          ENDIF
        ENDDO
      ENDDO
    ENDDO     ! *** END OF DOMAIN
    !$OMP END DO
    
  ENDIF

  ! *** CALCS WITHOUT INTERNAL RADIATION SHEAR STRESS DUE TO WAVE ACTION
  IF( ISWAVE <= 1 .OR. ISWAVE == 3 )THEN

    !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,DELB,CTE3TMP,PQQB,PQQU,PQQ,PQQL,PQQV)
    DO ND=1,NDM  

      DO K=1,KS
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LN=LNC(L)
          LE=LEC(L)
          UUU(L,K) = QQ(L,K)*HP(L)        + DELT*( FUHU(L,K)-FUHU(LE,K) + FVHU(L,K)-FVHU(LN,K) + (FWQQ(L,K) -FWQQ(L,K+1) )*DZIG(L,K))*DXYIP(L)
          VVV(L,K) = QQL(L,K)*HP(L)*HP(L) + DELT*( FUHV(L,K)-FUHV(LE,K) + FVHV(L,K)-FVHV(LN,K) + (FWQQL(L,K)-FWQQL(L,K+1))*DZIG(L,K))*DXYIP(L)
        ENDDO
      ENDDO

      DO K=1,KS
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)
          
          DELB=B(L,K)-B(L,K+1)
          CTE3TMP=CTE3
          IF( DELB < 0.0 )CTE3TMP=CTE1
          
          LN=LNC(L)
          LE=LEC(L)

          PQQB = AB(L,K)*GP*HP(L)*DZIG(L,K)*(B(L,K+1)-B(L,K))
          
          PQQU = AV(L,K)*DZIGSD4U(L,K)*(U(LE,K+1)-U(LE,K) + U(L,K+1)-U(L,K))**2  
          PQQV = AV(L,K)*DZIGSD4V(L,K)*(V(LN,K+1)-V(LN,K) + V(L,K+1)-V(L,K))**2  
          PQQ  = DELT*( PQQB+PQQU+PQQV+PQQVEGE(L,K)+PQQMHKE(L,K)+PQQSUPE(L,K) )
          UUU(L,K) = UUU(L,K) + 2.*PQQ
          PQQL = DELT*HP(L)*(CTE3TMP*PQQB + CTE1*(PQQU+PQQV) + CE4VEG*PQQVEGE(L,K) + CE4MHK*PQQMHKE(L,K) + CE4SUP*PQQSUPE(L,K))
          VVV(L,K) = VVV(L,K) + DML(L,K)*PQQL
        ENDDO
      ENDDO
    ENDDO     ! *** END OF DOMAIN
    !$OMP END DO

  ENDIF

  ! *** Wave Effects on Boundary and Water Column
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN

    !$OMP DO PRIVATE(ND,K,LP,L,LN,LE,LS)  &
    !$OMP    PRIVATE(DELB,CTE3TMP,PQQB,PQQU,PQQV,PQQW,PQQ,PQQL,FFTMP)
    DO ND=1,NDM  

      ! *** SUM VERTICAL WAVE DISSIPATION DUE TO TKE CLOSURE
      DO K=1,KS
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( LWVMASK(L) )THEN
            ! ***                       BOTTOM                     TOP
            TVAR1W(L,K) = WVDTKEM(K)*WV(L).DISSIPA(K) + WVDTKEP(K)*WV(L).DISSIPA(K+1)
          ELSE
            TVAR1W(L,K)=0.0
          ENDIF
        ENDDO
      ENDDO

      DO K=1,KS
        DO LP=1,LLWET(K,ND) 
          L=LKWET(LP,K,ND)  
          LN=LNC(L)
          LE=LEC(L)
          
          CTE3TMP=CTE3
          DELB=B(L,K)-B(L,K+1)  
          IF( DELB < 0.0 )CTE3TMP=CTE1
          
          PQQB = AB(L,K)*GP*HP(L)*DZIG(L,K)*(B(L,K+1)-B(L,K))
          PQQU = AV(L,K)*DZIGSD4U(L,K)*(U(LE,K+1)-U(LE,K)+U(L,K+1)-U(L,K))**2 
          PQQV = AV(L,K)*DZIGSD4V(L,K)*(V(LN,K+1)-V(LN,K)+V(L,K+1)-V(L,K))**2
          PQQW = WVFACT*TVAR1W(L,K)

          PQQ = DELT*( PQQU + PQQV + PQQB + PQQW + PQQVEGE(L,K) + PQQMHKE(L,K) + PQQSUPE(L,K))
          FFTMP = MAX( FUHU(L,K)-FUHU(LE,K) + FVHU(L,K)-FVHU(LN,K) + (FWQQ(L,K)-FWQQ(L,K+1))*DZIG(L,K),0.)
          UUU(L,K) = QQ(L,K)*HP(L) + DELT*FFTMP*DXYIP(L) + 2.*PQQ

          PQQL  = DELT*HP(L)*(CTE3TMP*PQQB + CTE1*(PQQU+PQQV+PQQW) + CE4VEG*PQQVEGE(L,K) + CE4MHK*PQQMHKE(L,K) + CE4SUP*PQQSUPE(L,K))
          FFTMP = MAX( FUHV(L,K)-FUHV(LE,K) + FVHV(L,K)-FVHV(LN,K) + (FWQQL(L,K)-FWQQL(L,K+1))*DZIG(L,K),0.)
          VVV(L,K) = QQL(L,K)*HP(L)*HP(L)+DELT*FFTMP*DXYIP(L) + DML(L,K)*PQQL
        ENDDO
      ENDDO

    ENDDO     ! *** END OF DOMAIN
    !$OMP END DO

  ENDIF

  ! *****************************************************************************
  IF( KC <= 2 )THEN
    ! *** 1 AND 2 LAYER CASE
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,LN,CLQTMP,CUQTMP,CLQLTMP,CUQLTMP,CMQTMP,CMQLTMP,EQ,EQL)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL  
        L=LWET(LP)  
        CLQTMP = -DELT*CDZKK(L,1) *AQ(L,1)*HPI(L)
        CUQTMP = -DELT*CDZKKP(L,1)*AQ(L,2)*HPI(L)
        CLQLTMP = SQLDSQ*CLQTMP
        CUQLTMP = SQLDSQ*CUQTMP
        CMQTMP  = 1. - CLQTMP-CUQTMP+2.*DELT*QQSQR(L,1)/(CTURBB1(L,1)*DML(L,1)*HP(L)) &
                  + DELT*HPI(L)*(PQQVEGI(L,1)+PQQMHKI(L,1)+PQQSUPI(L,1))
        CMQLTMP = 1. - CLQLTMP-CUQLTMP+DELT*(QQSQR(L,1)/(CTURBB1(L,1)*DML(L,1)*HP(L)))*(1.+CTE4*DML(L,1)*DML(L,1)*FPROX(L,1)) &
                  + DELT*HPI(L)*(CE4VEG*PQQVEGI(L,1)+CE4MHK*PQQMHKI(L,1)+CE4SUP*PQQSUPI(L,1))
        EQ = 1./CMQTMP
        EQL = 1./CMQLTMP
        CU1(L,1) = CUQTMP*EQ
        CU2(L,1) = CUQLTMP*EQL
        UUU(L,1) = (UUU(L,1)-CLQTMP*HP(L)*QQ(L,0)-CUQTMP*HP(L)*QQ(L,KC))*EQ
        VVV(L,1) = VVV(L,1)*EQL
      ENDDO
    ENDDO     ! *** END OF DOMAIN
    !$OMP END DO

  ELSE   ! IF( KC > 2 )THEN

    ! *** MULTI-LAYER CASE
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K,LN,CLQTMP,CUQTMP,CLQLTMP,CUQLTMP,CMQTMP,CMQLTMP,EQ,EQL)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      ! *** BOTTOM ACTIVE LAYER
      DO LP=1,LLWET(KS,ND)
        L=LKWET(LP,KS,ND)  
        K=KSZ(L)
        CLQTMP = -DELT*CDZKK(L,K) *AQ(L,K)  *HPI(L)
        CUQTMP = -DELT*CDZKKP(L,K)*AQ(L,K+1)*HPI(L)
        CLQLTMP = SQLDSQ*CLQTMP
        CUQLTMP = SQLDSQ*CUQTMP
        CMQTMP  = 1. - CLQTMP -CUQTMP+2.*DELT*QQSQR(L,K) /(CTURBB1(L,K)*DML(L,K)*HP(L)) &
                      + DELT*HPI(L)*(PQQVEGI(L,K)+PQQMHKI(L,K)+PQQSUPI(L,K))
        CMQLTMP = 1. - CLQLTMP-CUQLTMP  +DELT*(QQSQR(L,K)/(CTURBB1(L,K)*DML(L,K)*HP(L)))*(1. + CTE4*DML(L,K)*DML(L,K)*FPROX(L,K)) &
                      + DELT*HPI(L)*(CE4VEG*PQQVEGI(L,K) + CE4MHK*PQQMHKI(L,K) + CE4SUP*PQQSUPI(L,K)) 
        EQ  = 1./CMQTMP
        EQL = 1./CMQLTMP
        CU1(L,K) = CUQTMP*EQ
        CU2(L,K) = CUQLTMP*EQL
        IF( K == KS )THEN
          ! *** CASE OF WHEN BOTTOM ACTIVE LAYER = KC-1
          UUU(L,K) = ( UUU(L,K) - CLQTMP*HP(L)*QQ(L,0) - CUQTMP*HP(L)*QQ(L,KC) )*EQ
          VVV(L,K) = VVV(L,K)*EQL
        ELSE
          UUU(L,K)  = ( UUU(L,K) - CLQTMP*HP(L)*QQ(L,0) )*EQ
          VVV(L,K)  = VVV(L,K)*EQL
          CUQTMP    = -DELT*CDZKKP(L,KS)*AQ(L,KC)
          UUU(L,KS) = UUU(L,KS) - CUQTMP*QQ(L,KC)
        ENDIF
      ENDDO

      ! *** MIDDLE LAYERS
      DO K=2,KS
        DO LP=1,LLWET(K-1,ND)
          L=LKWET(LP,K-1,ND) 
          CLQTMP = -DELT*CDZKK(L,K) *AQ(L,K)  *HPI(L)
          CUQTMP = -DELT*CDZKKP(L,K)*AQ(L,K+1)*HPI(L)
          CLQLTMP = SQLDSQ*CLQTMP
          CUQLTMP = SQLDSQ*CUQTMP
          CMQTMP  = 1. - CLQTMP-CUQTMP+2.*DELT*QQSQR(L,K)/(CTURBB1(L,K)*DML(L,K)*HP(L)) &
                       + DELT*HPI(L)*(PQQVEGI(L,K)+PQQMHKI(L,K)+PQQSUPI(L,K))
          CMQLTMP = 1. - CLQLTMP-CUQLTMP+DELT*(QQSQR(L,K)/(CTURBB1(L,K)*DML(L,K)*HP(L)))*(1. + CTE4*DML(L,K)*DML(L,K)*FPROX(L,K)) &
                       + DELT*HPI(L)*(CE4VEG*PQQVEGI(L,K) + CE4MHK*PQQMHKI(L,K) + CE4SUP*PQQSUPI(L,K))
          EQ  = 1./(CMQTMP-CLQTMP*CU1(L,K-1))
          EQL = 1./(CMQLTMP-CLQLTMP*CU2(L,K-1))
          CU1(L,K) = CUQTMP*EQ
          CU2(L,K) = CUQLTMP*EQL
          UUU(L,K) = (UUU(L,K)-CLQTMP*UUU(L,K-1))*EQ
          VVV(L,K) = (VVV(L,K)-CLQLTMP*VVV(L,K-1))*EQL
        ENDDO
      ENDDO
      
      DO K=KS-1,1,-1
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          UUU(L,K)=UUU(L,K)-CU1(L,K)*UUU(L,K+1)
          VVV(L,K)=VVV(L,K)-CU2(L,K)*VVV(L,K+1)
        ENDDO
      ENDDO
    ENDDO     ! *** END OF DOMAIN
    !$OMP END DO

  ENDIF    ! *** END OF KC>2 LOOP

  ! *****************************************************************************
  ! **  ORIGINAL FORM MODIFIED FOR DIMENSIONAL LENGTH SCALE TRANSPORT                 
  !$OMP DO PRIVATE(ND,K,LP,L,QQHDH,DMLTMP,DELB,DMLMAX)
  DO ND=1,NDM  
    DO K=1,KS
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        QQHDH    = UUU(L,K)*HPI(L)
        QQ(L,K)  = MAX(QQHDH,QQMIN)
      ENDDO
    ENDDO

    ! *** DIMENSIONAL LENGTH SCALE TRANSPORT                    
    DO K=1,KS
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        QQHDH = VVV(L,K)*HPI(L)*HPI(L)
        QQHDH = MIN(QQHDH,HP(L))              ! LIMIT DML
        QQL(L,K) = MAX(QQHDH,QQLMIN)
        DMLTMP = QQL(L,K)/QQ(L,K)
        DMLTMP = MAX(DMLTMP,DMLMIN)
        DELB = B(L,K)-B(L,K+1)
        IF( DELB > 1.0E-12 .AND. ISLLIM == 2 )THEN
          DMLMAX = SQRT(RIQMAX)*SQRT(QQ(L,K)/(G*HP(L)*DZIG(L,K)*DELB))
          DML(L,K) = MIN(DMLMAX,DMLTMP)
          QQL(L,K) = QQ(L,K)*DML(L,K)
        ELSE
          DML(L,K) = DMLTMP
        ENDIF
      ENDDO
    ENDDO

  ENDDO     ! *** END OF DOMAIN
  !$OMP END DO

  ! ****************************************************************************
  ! *** CHECK FOR DEPTHS LESS THAN ZBR
  IF( ISDRY > 0 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K) 
    DO ND=1,NDM  
      DO LP=1,LLWET(KS,ND)  
        L = LKWET(LP,KS,ND)
        IF( HPK(L,KSZ(L)) < ZBR(L) )THEN
          ! *** SPECIAL CASE: LAYER 1 OR MORE THICKNESSES < Z0
          DO K=KSZ(L),KS
            IF( HP(L)*Z(L,K-1)>ZBR(L) )EXIT
            QQ(L,K)=QQMIN
            QQL(L,K)=QQLMIN
            DML(L,K)=DMLMIN
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !$OMP END DO
  ENDIF
        
  !$OMP SINGLE
  DO K=1,KS
    DO LL=1,NPBS
      L=LPBS(LL)
      LN=LNC(L)
      QQ(L,K)=QQ(LN,K)
      QQL(L,K)=QQL(LN,K)
      DML(L,K)=DML(LN,K)
    ENDDO
  ENDDO
  DO K=1,KS
    DO LL=1,NPBW
      L=LPBW(LL)
      QQ(L,K)=QQ(LEC(L),K)
      QQL(L,K)=QQL(LEC(L),K)
      DML(L,K)=DML(LEC(L),K)
    ENDDO
  ENDDO
  DO K=1,KS
    DO LL=1,NPBE
      L=LPBE(LL)
      QQ(L,K)=QQ(LWC(L),K)
      QQL(L,K)=QQL(LWC(L),K)
      DML(L,K)=DML(LWC(L),K)
    ENDDO
  ENDDO

  DO K=1,KS
    DO LL=1,NPBN
      L=LPBN(LL)
      LS=LSC(L)
      QQ(L,K)=QQ(LS,K)
      QQL(L,K)=QQL(LS,K)
      DML(L,K)=DML(LS,K)
    ENDDO
  ENDDO
  !$OMP END SINGLE

  ! *** SAVE THE SQRT OF THE TURBULENCE (M/S)
  !$OMP DO PRIVATE(ND,K,LP,L)
  DO ND=1,NDM  
    DO K=1,KS
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        QQSQR(L,K)=SQRT(QQ(L,K))  
      ENDDO
    ENDDO
  ENDDO 
  !$OMP END DO
  !$OMP END PARALLEL
  
  RETURN

END

