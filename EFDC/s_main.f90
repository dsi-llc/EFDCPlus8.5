!**********************************************************************!
SUBROUTINE SEDZLJ_MAIN

  !**********************************************************************!
  !
  ! **  SUBROUTINE CALSED CALCULATES COHESIVE SEDIMENT SETTLING,
  ! **  DEPOSITION AND RESUSPENSION According to SEDZLJ MODEL
  !
  ! ORIGINAL DATE :  May 24, 2006
  !  Craig Jones and Scott James
  ! REVISED: SIGMA-ZED AND OMP - 2016-11-07
  !  Paul M. Craig

  !**********************************************************************!

  USE GLOBAL
  IMPLICIT NONE

  REAL(RKD) :: CRNUM,SEDAVG,GRADSED,CLEFT,CRIGHT,WVEL,SMASSD,SMASSU,SMASS
  REAL(RKD),SAVE :: LASTTIME=-9999.
  REAL(RKD),DIMENSION(NSCM) :: BLFLUXU,BLFLUXD
  REAL    :: T1TMP
  INTEGER :: L,K,NS,ND,LL,LF,LP,LUTMP,LDTMP,NSB
  
  ! *** ENFORCE STARTUP FOR ENTIRE SEDIMENT PROCESS
  IF( LASTTIME == -9999. ) LASTTIME = TIMEDAY
  T1TMP=SECNDS(0.0)  
  
  SSGI(1:NSCM)=1.0/(1.0E6*SSG(1:NSCM)) 
  
  ! *** INITIALIZE CURRENT TIME STEP VARIABLES FOR SEDZLJ
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    DO LP=LF,LL
      L=LWET(LP)
      HPCM(L) = 100.0*HP(L)
    ENDDO
  ENDDO 
  !$OMP END PARALLEL DO

  ! *** Calculates the shear stresses from the velocity output of the hydraulic model.
  CALL SEDZLJ_SHEAR
  
  IF( ISSLOPE ) CALL SEDZLJ_SLOPE

  ! *** Calling bed load calculation.  BEDLOADJ subroutine provides CBL, the sediment concentration load.
  IF( NCALC_BL > 0 )THEN
    CALL BEDLOADJ
  ENDIF
  
  ! *** Calls the sedflume transport. NOTE: ISEDTIME is the time difference between when sediment
  ! *** transport begins versus the beginning of hydraulic transport.
  IF( NSEDFLUME == 1 )THEN

    !*************************************************************************************!
    !   SEDZLJ Sediment Transport
    !
    !$OMP PARALLEL DEFAULT(SHARED)

    !$OMP DO PRIVATE(ND,NS,K,LP,L,WVEL,CLEFT,CRIGHT)
    DO ND=1,NDM  
      ! *** ZERO LAYER FLUXES, INCLUDING BOTTOM LAYER
      DO NS=1,NSCM
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SEDF(L,K,NS) = 0.0
            WSETA(L,K,NS) = WSEDO(NS)
          ENDDO
          
          ! *** ZERO BED/WATER INTERFACE FLUX
          IF( K == KS )THEN
            DO LP=1,LLWET(KC,ND)
              L=LKWET(LP,KC,ND)  
              SEDF(L,0,NS) = 0.0
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    
      !-----------------------------------------------------------------------------------!
      !
      ! *** HORIZONTAL LOOPS
      ! *** These loops account for sediment flux between the water layers due to settling.
      IF( KC /= 1 )THEN
        DO NS=1,NSCM
          DO LP=1,LLWET(KC,ND)
            L=LKWET(LP,KC,ND)  
            WVEL = DTSEDJ*HPKI(L,KC)
            CLEFT = 1.0 + WSETA(L,KC-1,NS)*WVEL
            CRIGHT = MAX(SED(L,KC,NS),0.0)
            SED(L,KC,NS) = CRIGHT/CLEFT
            SEDF(L,KC-1,NS) = -WSETA(L,KC-1,NS)*SED(L,KC,NS)
          ENDDO
           
          IF( KC /= 2 )THEN
            DO K=KS,2,-1
              DO LP=1,LLWET(K-1,ND)
                L=LKWET(LP,K-1,ND) 
                WVEL = DTSEDJ*HPKI(L,K)
                CLEFT = 1.0+WSETA(L,K-1,NS)*WVEL
                CRIGHT = MAX(SED(L,K,NS),0.0) - SEDF(L,K,NS)*WVEL
                SED(L,K,NS) = CRIGHT/CLEFT
                SEDF(L,K-1,NS) = -WSETA(L,K-1,NS)*SED(L,K,NS)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    !$OMP END DO
    
    ! **********************************************************************************
    ! *** Erosion/Deposition Loop
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      DO LP=LF,LL
        L=LWET(LP)
        
        ! *** Calculate total flux from bed and determine if 
        ! *** there is enough sediment to erode from SEDZLJ
        CALL SEDZLJ(L)
        
      ENDDO
    ENDDO
    !$OMP END DO

    ! **********************************************************************************
    ! *** Compute cell centered bedload fluxes for outflow or recirculation boundary
    !$OMP SINGLE
    IF( NSBDLDBC > 0 )THEN
      BLFLUXU = 0.0
      BLFLUXD = 0.0
      DO NSB=1,NSBDLDBC
        LUTMP = LSBLBCU(NSB)
        LDTMP = LSBLBCD(NSB)

        SMASSU = TSED(1,LUTMP)
        DO NS=NNONCO,NSCM
          BLFLUXU(NS) = ( DBL(LUTMP,NS) - EBL(LUTMP,NS) )/10000.             ! *** Flux in/out of bedload CBL                   (g/cm^2)
          TSED(1,LUTMP) = TSED(1,LUTMP) - BLFLUXU(NS)                        ! *** Back out erosion/depostion due to bedload    (g/cm^2)
          QBLFLUX(LUTMP,NS)  = -DBL(LUTMP,NS)/10000.                         ! *** Flux out of bed for bedload                  (g/cm^2)
          QSBDLDOT(LUTMP,NS) =  DBL(LUTMP,NS)*DXYP(LUTMP)/DTSEDJ             ! *** Sediment flux from upstream cell             (g/s)
          EBL(LUTMP,NS) = 0.0                            
        ENDDO
        
        ! *** UPDATE PERSED WITH MACHINE PRECISION A CONSIDERATION: UPSTREAM
        SMASS = 0.0
        DO NS=1,NSCM
          SMASS = SMASS + MAX(PERSED(NS,1,LUTMP)*SMASSU - BLFLUXU(NS),0.0)
        ENDDO
        IF( SMASS > 0.0 )THEN 
          DO NS=1,NSCM
            PERSED(NS,1,LUTMP) = MAX(PERSED(NS,1,LUTMP)*SMASSU - BLFLUXU(NS),0.0)/SMASS
          ENDDO
          TSED(1,LUTMP) = SMASS
        ELSE
          TSED(1,LUTMP) = 0.0
          PERSED(1:NSCM,1,LUTMP) = 0.0
        ENDIF

        ! *** Update EFDC variables
        SEDDIA50(LUTMP,1) = SUM(PERSED(1:NSCM,1,LUTMP)*D50(1:NSCM))
        HBED(LUTMP,1)  = 0.01_8*TSED(1,LUTMP)/BULKDENS(1,LUTMP)  
        SEDBT(LUTMP,1) = TSED(1,LUTMP)*10000.
        SEDB(LUTMP,1,1:NSCM)  = SEDBT(LUTMP,1)*PERSED(1:NSCM,1,LUTMP)

        IF( LDTMP > 0 )THEN
          ! *** The flow is returning.  Add to downstream cell
          SMASSD = TSED(1,LDTMP)
          DO NS=NNONCO,NSCM
            BLFLUXD(NS)    = DBL(LUTMP,NS)/10000.*DXYP(LUTMP)*DXYIP(LDTMP)     ! *** Depostion due to bedload from upstream cell  (g/cm^2)
            TSED(1,LDTMP) = TSED(1,LDTMP) + BLFLUXD(NS)                        ! *** Add mass to downstream cell                  (g/cm^2)
          ENDDO
          
          ! *** UPDATE PERSED WITH MACHINE PRECISION A CONSIDERATION: DOWNSTREAM
          SMASS = 0.0
          DO NS=1,NSCM
            SMASS = SMASS + MAX(PERSED(NS,1,LDTMP)*SMASSD + BLFLUXD(NS),0.0)
          ENDDO
          IF( SMASS > 0.0 )THEN 
            DO NS=1,NSCM
              PERSED(NS,1,LDTMP) = MAX(PERSED(NS,1,LDTMP)*SMASSD + BLFLUXD(NS),0.0)/SMASS
            ENDDO
            TSED(1,LDTMP) = SMASS
          ELSE
            TSED(1,LUTMP) = 0.0
            PERSED(1:NSCM,1,LUTMP) = 0.0
          ENDIF
          
          ! *** Update EFDC variables
          SEDDIA50(LDTMP,1) = SUM(PERSED(1:NSCM,1,LDTMP)*D50(1:NSCM))
          HBED(LDTMP,1)  = 0.01_8*TSED(1,LDTMP)/BULKDENS(1,LDTMP)  
          SEDBT(LDTMP,1) = TSED(1,LDTMP)*10000.
          SEDB(LDTMP,1,1:NSCM)  = SEDBT(LDTMP,1)*PERSED(1:NSCM,1,LDTMP)
        ENDIF
          
      ENDDO
    ENDIF
    !$OMP END SINGLE
    
    ! **********************************************************************************
    ! ***  ANTI-DIFFUSION AND FINAL BED/WATER FLUX FOR SEDIMENTS WHEN KC > 1
    IF( KC > 1 .AND. ISTOPT(6) > 0 )THEN
      !$OMP DO PRIVATE(ND,LF,LL,NS,K,LP,L,CRNUM,GRADSED,SEDAVG)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
      
        DO NS=1,NSCM
          ! ***  ANTI-DIFFUSION OF SEDIMENT FOR KC > 1
          DO K=1,KS
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              CRNUM = 1.+DTSEDJ*WSETA(L,K,NS)*HPKI(L,K+1)
              GRADSED = ( SED(L,K+1,NS) - SED(L,K,NS) )/(DZC(L,K+1)+DZC(L,K))
              SEDAVG = 0.5*( SED(L,K+1,NS) + SED(L,K,NS) + 1.E-16 )
              WSETA(L,K,NS) = -CRNUM*DZC(L,K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG
            ENDDO
          ENDDO

          ! *** TVAR1S=LOWER DIAGONAL
          DO LP=LF,LL
            L  =LWET(LP)  
            TVAR1S(L,KSZ(L)) = 0.
          ENDDO
          DO K=2,KC
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TVAR1S(L,K) = MIN(WSETA(L,K-1,NS),0.)
            ENDDO
          ENDDO

          ! *** TVAR1N=UPPER DIAGONAL
          DO LP=LF,LL
            L = LWET(LP)  
            TVAR1N(L,KC) = 0.
          ENDDO
          DO K=1,KS
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TVAR1N(L,K) = -MAX(WSETA(L,K,NS),0.)
            ENDDO
          ENDDO
        
          ! *** TVAR1W=MAIN DIAGONAL
          DO LP=LF,LL
            L=LWET(LP)  
            TVAR1W(L,KSZ(L)) = DELTI*HPK(L,KSZ(L)) - MIN(WSETA(L,1,NS),0.)
            TVAR1W(L,KC)     = DELTI*HPK(L,KC)     + MAX(WSETA(L,KC-1,NS),0.)
          ENDDO
          DO K=2,KS
            DO LP=1,LLWET(K-1,ND)
              L = LKWET(LP,K-1,ND) 
              TVAR1W(L,K) = DELTI*HPK(L,K) + MAX(WSETA(L,K-1,NS),0.) - MIN(WSETA(L,K,NS),0.)
            ENDDO
          ENDDO

          ! *** TVAR1E=RIGHT HAND SIDE
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              TVAR1E(L,K) = DELTI*HPK(L,K)*SED(L,K,NS)
            ENDDO
          ENDDO

          ! *** TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
          DO LP=LF,LL
            L=LWET(LP)   
            TVAR3S(L) =TVAR1W(L,KSZ(L))
          ENDDO
          DO LP=LF,LL
            L=LWET(LP)  
            TVAR2N(L,KSZ(L)) = TVAR1E(L,KSZ(L))/TVAR3S(L)
          ENDDO
          DO K=2,KC
            DO LP=1,LLWET(K-1,ND)
              L = LKWET(LP,K-1,ND) 
              TVAR2S(L,K) = TVAR1N(L,K-1)/TVAR3S(L)
              TVAR3S(L)   = TVAR1W(L,K)  - TVAR1S(L,K)*TVAR2S(L,K)
              TVAR2N(L,K) = (TVAR1E(L,K) - TVAR1S(L,K)*TVAR2N(L,K-1))/TVAR3S(L)
            ENDDO
          ENDDO
          DO K=KS,1,-1
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              TVAR2N(L,K) = TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
            ENDDO
          ENDDO
          DO K=1,KC
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              SED(L,K,NS) = TVAR2N(L,K)
            ENDDO
          ENDDO
          ! *** END OF ANTI-DIFFUSION CALCULATIONS
    
          ! *** RECOMPUTE FLUXES FROM FINAL SED: KC-1 LAYER
          DO LP=LF,LL
            L=LWET(LP)  
            SEDF(L,KS,NS) = DELTI*HPK(L,KC)*( SED(L,KC,NS) - SEDS(L,KC,NS) )
          ENDDO 

          ! *** RECOMPUTE FLUXES FROM FINAL SED: MIDDLE LAYERS
          DO K=KS-1,1,-1
            DO LP=1,LLWET(K+1,ND)
              L=LKWET(LP,K+1,ND)  
              SEDF(L,K,NS) = DELTI*HPK(L,K+1)*( SED(L,K+1,NS)-SEDS(L,K+1,NS) ) + SEDF(L,K+1,NS)
            ENDDO  
          ENDDO  
        ENDDO   ! *** End do of NS=1,NSCM loop
        
      ENDDO  ! *** End of Domain Loop
      !$OMP END DO
    ENDIF    ! *** END OF KC > 1 ANTI-DIFFUSION BLOCK    
    
    !$OMP END PARALLEL
     
  ELSEIF( NSEDFLUME == 2 )THEN
    ! **********************************************************************************
    ! *** SEDZLJ Sediment and Contaminant Transport model
    ! *** PMC - 2017-08 - Deprecated as toxic transport is handled by CALTOX and CALTOXB
    STOP 'NSEDFLUME=2 Deprecated (EE8.3) as toxic transport is handled by CALTOX and CALTOXB'
  ENDIF
  ! *** End of condition  NSEDFLUME == 1 .AND. N >= ISEDTIME
  ! ************************************************************************************

  ! *** PMC deleted the MORPHJ routine as it was redundant and not as complete as the morphology update in SSEDTOX  (2016-12)

  ! *** Mass balance routine will be called every X minutes divided by 1440
  IF( ISBAL > 0 .AND. ( NDYCOUNT == 0 .OR. ( TIMEDAY >= LASTTIME+5./1440. ) ) ) THEN
    CALL MASS_BALANCE
    LASTTIME = TIMEDAY
  END IF
  NDYCOUNT=NDYCOUNT+1

  TSSEDZLJ = TSSEDZLJ+SECNDS(T1TMP)   

  RETURN

END SUBROUTINE SEDZLJ_MAIN
