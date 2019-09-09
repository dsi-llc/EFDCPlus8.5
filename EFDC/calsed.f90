SUBROUTINE CALSED  

  !**********************************************************************!
  ! **  SUBROUTINE CALSED CALCULATES COHESIVE SEDIMENT SETTLING,  
  ! **  DEPOSITION AND RESUSPENSION AND IS CALLED FOR SSEDTOX  
  !  
  ! *** STANDARD EFDC COHESIVE SEDIMENT TRANSPORT
  ! *** NOT USED FOR SEDZLJ
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-07       PAUL M. CRAIG     CHANGED COHESIVE CONCENTRATIONS FOR SETTLING
  !                                      TO TOTAL CONCENTRATION (SEDT) RATHER THAN
  !                                      BY CLASS (SED)
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  ! 2011-03-02       PAUL M. CRAIG     RESTRUCTURED AND CORRECTED CODE
  !                                      REMOVED KC DEPENDENT DUPLCIATE CODE
  !                                      EXPANDED CAPABILITIES AND ADDED OMP
  !**********************************************************************!

  USE GLOBAL
  
  IMPLICIT NONE
  
  INTEGER :: K,LP,L,NS,ND,LF,LL,LN,K1,IOBC,IFLAG
  REAL    :: TIME,STRESS,SHEAR,TAUBC,UTMP,VTMP,CURANG,TAUB2,UUSTARTMP,TAUDSS
  REAL    :: WVEL,CLEFT,CRIGHT,WESE,TAUE,TMPSTR,TMPSEDHID,TAUBHYDRO,PROBDEP,WSETMP
  REAL    :: SEDBTMP,CRNUM,GRADSED,SEDAVG,WESEMX,TAURTMP
  REAL, EXTERNAL :: CSEDSET,FPROBDEP 

  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: IFIRST

  IF(  .NOT. ALLOCATED(IFIRST) )THEN
    ALLOCATE(IFIRST(2))
    IFIRST = 0

    IF( ISEDVW == 0 )THEN
      ! *** CONSTANT (ONLY ASSIGN AT START OF RUN)
      DO NS=1,NSED
        DO K=0,KS
          ! *** USING 2,LA INTENTIONAL TO ENSURE ASSIGNMENT EVEN IF DRY CELLS EXIST
          DO L=2,LA
            WSETA(L,K,NS) = WSEDO(NS)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  
  !**********************************************************************!  
  IF( ISDYNSTP == 0 )THEN  
    TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON  
  ELSE  
    TIME=TIMESEC/TCON  
  ENDIF  

  !**********************************************************************!
  ! **
  ! **  COHESIVE SEDIMENT FLUX

  IF( LADRY > 0 )THEN
    DO K=1,KC
      DO LP=1,LADRY
        L=LDRY(LP)  
        TVAR2E(L,K)=0.
        TVAR2N(L,K)=0.
      ENDDO
    ENDDO
  ENDIF
  
  DO NS=1,NSED
    DSEDGMM = 1./(1.E6*SSG(NS))    ! *** SPECIFIC VOLUME (M**3/G)

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,K1,LN,IOBC)  &
    !$OMP             PRIVATE(STRESS,SHEAR,TAUBC,UTMP,VTMP,CURANG,TAUB2,UUSTARTMP,TAUDSS) &
    !$OMP             PRIVATE(WVEL,CLEFT,CRIGHT,WESE,TAUE,TMPSTR,TMPSEDHID,TAUBHYDRO,PROBDEP,WSETMP) &
    !$OMP             PRIVATE(SEDBTMP,CRNUM,GRADSED,SEDAVG,WESEMX,TAURTMP)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      
      !----------------------------------------------------------------------!
      ! ***  SET SETTLING VELOCITIES  ( ISEDVW == 0 SET ABOVE )
      IF( ISEDVW == 1 )THEN
        ! *** SIMPLE CONCENTRATION DEPENDENT (HUANG AND METHA)
        DO K=0,KS
          DO LP=1,LLWET(K+1,ND)
            L = LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),0.0,ISEDVW)
          ENDDO
        ENDDO
      ENDIF

      IF( ISEDVW == 2 )THEN
        ! *** CONCENTRATION AND SHEAR/TURBULENCE DEPENDENT COHESIVE SETTLING  (SHRESTA AND ORLOB)
        IF( ISWAVE == 0 )THEN
          DO K=0,KS
            DO LP=1,LLWET(K+1,ND)
              L=LKWET(LP,K+1,ND)  
              STRESS=QQ(L,0)/CTURB2
              SHEAR=2.*HPKI(L,K+1)*SQRT(STRESS)/VKC
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),SHEAR,ISEDVW)
            ENDDO
          ENDDO
        ELSE
          K=0
          DO LP=LF,LL
            L=LWET(LP)  
            K1=KSZ(L)
            IF( LWVMASK(L) )THEN
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(LEC(L),   K1) + U(L,K1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),K1) + V(L,K1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              SHEAR=2.*HPKI(L,K1)*SQRT(STRESS)/VKC
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),SHEAR,ISEDVW)
            ELSE
              STRESS=QQ(L,0)/CTURB2
              SHEAR=2.*HPKI(L,K1)*SQRT(STRESS)/VKC
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),SHEAR,ISEDVW)
            ENDIF
          ENDDO
          IF( KC > 1 )THEN
            DO K=1,KS
              DO LP=1,LLWET(K+1,ND)
                L=LKWET(LP,K+1,ND)  
                IF( LWVMASK(L) )THEN
                  LN=LNC(L)
                  SHEAR=HPI(L)*SQRT( DZIGSD4U(L,K) )*SQRT( (U(LEC(L),K+1)-U(LEC(L),K)+U(L,K+1)-U(L,K))**2  &
                                                         + (V(LN ,K+1)-V(LN ,K)+V(L,K+1)-V(L,K))**2 )
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),SHEAR,ISEDVW)
                ELSE
                  STRESS=QQ(L,0)/CTURB2
                  SHEAR=2.*HPKI(L,K+1)*SQRT(STRESS)/VKC
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),SHEAR,ISEDVW)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF  
      ENDIF

      IF( ISEDVW >= 3 .AND. ISEDVW <= 4 )THEN
        ! *** CONCENTRATION AND SHEAR/TURBULENCE DEPENDENT COHESIVE SETTLING  (ZIEGLER AND NESBIT)
        IF( ISWAVE == 0 )THEN
          DO K=0,KS
            DO LP=1,LLWET(K+1,ND)
              L=LKWET(LP,K+1,ND)  
              STRESS=0.5*QQ(L,0)/CTURB2
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
            ENDDO
          ENDDO
        ELSE
          K=0
          DO LP=LF,LL
            L=LWET(LP)  
            K1=KSZ(L)
            IF( LWVMASK(L) )THEN
              TAUBC=QQ(L,0)/CTURB2
              UTMP=0.5*STCUV(L)*(U(LEC(L),   K1) + U(L,K1))+1.E-12
              VTMP=0.5*STCUV(L)*(V(LNC(L),K1) + V(L,K1))
              CURANG=ATAN2(VTMP,UTMP)
              TAUB2=TAUBC*TAUBC + (QQWV3(L)*QQWV3(L)) + 2.*TAUBC*QQWV3(L)*COS(CURANG-WV(L).DIR)
              TAUB2=MAX(TAUB2,0.)
              STRESS=SQRT(TAUB2)
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),STRESS,ISEDVW)
            ELSE
              STRESS=0.5*QQ(L,0)/CTURB2
              WSETA(L,K,NS) = CSEDSET(SEDT(L,K1),STRESS,ISEDVW)
            ENDIF
          ENDDO

          IF( KC > 1 )THEN
            DO K=1,KS
              DO LP=1,LLWET(K+1,ND)
                L=LKWET(LP,K+1,ND)  
                IF( LWVMASK(L) )THEN
                  LN=LNC(L)
                  STRESS = AV(L,K)*SQRT( DZIGSD4U(L,K)*(U(LEC(L),K+1)-U(LEC(L),K) + U(L,K+1)-U(L,K))**2  &
                                      +  DZIGSD4V(L,K)*(V(LN ,K+1)-V(LN ,K) + V(L,K+1)-V(L,K))**2 )
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
                ELSE
                  STRESS = 0.5*QQ(L,0)/CTURB2
                  WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
                ENDIF
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDIF

      IF( ISEDVW == 5 )THEN
        ! *** CONCENTRATION AND SHEAR/TURBULENCE DEPENDENT COHESIVE SETTLING  (HOUSATONIC)
        DO K=0,KS
          DO LP=1,LLWET(K+1,ND)
            L=LKWET(LP,K+1,ND)  
            UUSTARTMP=QCELLCTR(L)*SQRT(QQ(L,0)/CTURB2)
            STRESS = SQRT(HPI(L)*HPI(L)*UUSTARTMP)
            WSETA(L,K,NS) = CSEDSET(SEDT(L,K+1),STRESS,ISEDVW)
          ENDDO
        ENDDO
      ENDIF

      ! *** HANDLE LAYER 0 FOR SIGMA-ZED GRIDS
      IF( IGRIDV > 0 )THEN
        ! *** ASSIGN LAYER 0 WHEN KSZ(L) > 1
        DO LP = 1,LLWET(KC,ND)
          L = LKWET(LP,KC,ND)  
          IF( KSZ(L) > 1 )THEN
            WSETA(L,0,NS) = WSETA(L,KSZ(L)-1,NS) 
          ENDIF
        ENDDO
      ENDIF
      
      ! *** ZERO SETTLING VELOCITIES AT OPEN BOUNDARIES
      DO K=0,KS
        DO IOBC=1,NBCSOP  
          L=LOBCS(IOBC)  
          WSETA(L,K,NS)=0.0
        ENDDO  
      ENDDO  

      !----------------------------------------------------------------------!
      !
      ! **  HORIZONTAL LOOPS

      ! *** SET FLUX FOR THE BOTTOM OF THE TOP LAYER
      IF( KC > 1 )THEN
        K=KC
        DO LP=LF,LL
          L=LWET(LP)  
          SEDF(L,K,NS)=0.
          WVEL=DTSED*HPKI(L,K)
          CLEFT=1.+WSETA(L,K-1,NS)*WVEL
          CRIGHT=MAX(SED(L,K,NS),0.)
          SED(L,K,NS)=CRIGHT/CLEFT
          SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
        ENDDO

        ! *** SET FLUX FOR MIDDLE LAYERS
        DO K=KS,2,-1
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND)
            WVEL=DTSED*HPKI(L,K)
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL
            CRIGHT=MAX(SED(L,K,NS),0.)-SEDF(L,K,NS)*WVEL
            SED(L,K,NS)=CRIGHT/CLEFT
            SEDF(L,K-1,NS)=-WSETA(L,K-1,NS)*SED(L,K,NS)
          ENDDO
        ENDDO
      ENDIF
      
      ! *** SET CONSTANTS
      TMPSEDHID=1.0
      IF( IWRSP(1) == 4 )THEN
        TMPSTR=0.0
      ELSE
        TMPSTR=1.0
      ENDIF

      ! *** UPDATE SEDIMENT BED MASS & BOTTOM LAYER CONCENTRATION
      DO LP=LF,LL
        L=LWET(LP)
        
        ! *** SET MAXIMUM EROSION RATE
        WESEMX = DELTI*SEDB(L,KBT(L),NS)

        IF( TAUBSED(L) > TAURB(L,KBT(L)) )THEN
          ! *** MASS/BULK EROSION
          WESE = WRSPB(L,KBT(L))*VFRBED(L,KBT(L),NS)
          WESE = MIN(WESE,WESEMX)
        ELSE
          TAUE=0.
          IF( TAUBSED(L) > TAURS(L,KBT(L)) )THEN
            ! *** SURFACE EROSION
            WESE = WRSPS(L,KBT(L))*VFRBED(L,KBT(L),NS)

            ! *** SET NORMALIZING SHEAR BASED ON INPUT OPTION
            IF( IWRSP(1) >= 99 )THEN
              ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU FROM SSCOHSEDPMAP
              TAURTMP=TAUNS(L,KBT(L))
            ELSEIF ((IWRSP(1) >= 2) .AND. (IWRSP(1) < 4) )THEN
              ! *** 1 HWANG AND METHA - LAKE OKEECHOBEE
              ! *** 2 HAMRICK'S MODIFICATION OF SANFORD AND MAA
              ! *** 3 SAME AS 2 EXCEPT VOID RATIO OF COHESIVE SEDIMENT FRACTION IS USED
              TAURTMP=TAUR(1)
            ELSE
              ! *** USE DIRECTLY INPUT PROPERTIES
              TAURTMP=TAURS(L,KBT(L))
            ENDIF

            TAUE=(TAUBSED(L)-TMPSTR*TAURS(L,KBT(L)))/TAURTMP
            TAUE=MAX(TAUE,0.0)

            ! *** SET NON-COHESIVE HIDING FACTOR
            IF( ISTRAN(7) >= 1 .AND. COSEDHID(1)/=0.0 )THEN
              TMPSEDHID=(FRACCOH(L,KBT(L)))**COSEDHID(1)
            ENDIF
              
            ! *** FINALIZE SURFACE EROSION RATE
            IF( IWRSP(1) < 99 )THEN
              WESE=TMPSEDHID*WESE*( TAUE**TEXP(NS) )
            ELSE
              ! *** SITE SPECIFIC/SPACIALLY VARIABLE TAU EXPONENT FROM SSCOHSEDPMAP
              WESE=TMPSEDHID*WESE*( TAUE**TEXPS(L,KBT(L)) )
            ENDIF
            WESE=MIN(WESE,WESEMX)
          ELSE
            ! *** NO EROSION 
            WESE=0.0
          ENDIF
        ENDIF

        ! *** SET PROBABILITY OF DEPOSITION 
        IF( IWRSP(1) < 99 )THEN
          TAUDSS=TAUD(NS)
        ELSE
          TAUDSS=TAUDS(L)
        ENDIF
        TAUBHYDRO=QQ(L,0)/CTURB2
        PROBDEP=0.
        IF( ISPROBDEP(NS) == 0 )THEN
          ! *** KRONE PROBABILITY OF DEPOSITION USING COHESIVE GRAIN STRESS
          IF( TAUBSED(L) < TAUDSS) PROBDEP=(TAUDSS-TAUBSED(L))/TAUDSS
        ELSEIF( ISPROBDEP(NS) == 1 )THEN
          !*** KRONE PROBABILITY OF DEPOSITION USING TOTAL BED STRESS
          IF( TAUBHYDRO < TAUDSS) PROBDEP=(TAUDSS-TAUBHYDRO)/TAUDSS
        ELSEIF( ISPROBDEP(NS) == 2 )THEN
          ! *** PARTHENIADES PROBABILITY OF DEPOSITION USING COHESIVE GRAIN STRESS
          IF( TAUBSED(L) <= TAUDSS )THEN
            PROBDEP=1.0
          ELSE
            PROBDEP = FPROBDEP(TAUDSS,TAUBSED(L))
          ENDIF
        ELSEIF( ISPROBDEP(NS) == 3 )THEN
          ! *** PARTHENIADES PROBABILITY OF DEPOSITION USING TOTAL BED STRESS
          IF( TAUBHYDRO <= TAUDSS )THEN
            PROBDEP=1.0
          ELSE
            PROBDEP=FPROBDEP(TAUDSS,TAUBHYDRO)
          ENDIF
        ENDIF
        IF( SED(L,KSZ(L),NS) > SEDMDGM) PROBDEP=1.
          
        WSETMP = PROBDEP*WSETA(L,0,NS)
        WVEL   = DTSED*HPKI(L,KSZ(L))
        CLEFT  = 1. + WSETMP*WVEL
        CRIGHT = MAX(SED(L,KSZ(L),NS),0.) + (WESE-SEDF(L,KSZ(L),NS))*WVEL
        SED(L,KSZ(L),NS) = CRIGHT/CLEFT
        SEDF(L,0,NS)     = -WSETMP*SED(L,KSZ(L),NS) + WESE
        SEDBTMP          = SEDB(L,KBT(L),NS) - DTSED*SEDF(L,0,NS)

        ! *** Limit Erosion to Available Mass
        IF( SEDBTMP < 0.0 )THEN
          SEDF(L,0,NS) = DELTI*SEDB(L,KBT(L),NS)
          SEDBTMP = 0.0
          SED(L,KSZ(L),NS) = SEDS(L,KSZ(L),NS) + (SEDF(L,0,NS)-SEDF(L,KSZ(L),NS))*WVEL
        ENDIF
        SEDB1(L,KBT(L),NS) = SEDB(L,KBT(L),NS)
        SEDB(L,KBT(L),NS) = SEDBTMP
        QSBDTOP(L) = QSBDTOP(L) + DSEDGMM*SEDF(L,0,NS)
        IF( IBMECH == 1 .AND. SEDVRDT < 0.00001 )THEN
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*VDRBED(L,KBT(L))*SEDF(L,0,NS)  ! *** IF EITHER CONSTANT OR INSTANTLY CONSOLIDATING, USE BED VR  
        ELSE
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*( VDRBED(L,KBT(L))*MAX(SEDF(L,0,NS),0.) + VDRDEPO(NS)*MIN(SEDF(L,0,NS),0.) )
        ENDIF
      ENDDO      ! *** END OF UPDATE SEDIMENT BED MASS & BOTTOM LAYER LOOP

      !----------------------------------------------------------------------!
      ! **  ANTI-DIFFUSION OF COHESIVE SEDIMENT  KC > 1
      IF( ISTOPT(6) == 1 .AND. KC > 1 )THEN

        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            CRNUM=1.+DTSED*WSETA(L,K,NS)*HPKI(L,K+1)
            GRADSED=(SED(L,K+1,NS)-SED(L,K,NS))/(DZC(L,K+1)+DZC(L,K))
            SEDAVG=0.5*(SED(L,K+1,NS)+SED(L,K,NS)+1.E-16)
            WSETA(L,K,NS)=-CRNUM*DZC(L,K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG
          ENDDO
        ENDDO

        ! *** TVAR1S=LOWER DIAGONAL
        DO LP=LF,LL
          L=LWET(LP)  
          TVAR1S(L,KSZ(L))=0.
        ENDDO
        DO K=2,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TVAR1S(L,K)=MIN(WSETA(L,K-1,NS),0.)
          ENDDO
        ENDDO

        ! *** TVAR1N=UPPER DIAGONAL
        DO LP=LF,LL
          L=LWET(LP)  
          TVAR1N(L,KC)=0
        ENDDO
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TVAR1N(L,K)=-MAX(WSETA(L,K,NS),0.)
          ENDDO
        ENDDO
        
        ! *** TVAR1W=MAIN DIAGONAL
        DO LP=LF,LL
          L=LWET(LP)  
          TVAR1W(L,KSZ(L)) = DELTI*HPK(L,KSZ(L)) - MIN(WSETA(L,KSZ(L),NS),0.)
          TVAR1W(L,KC)     = DELTI*HPK(L,KC)     + MAX(WSETA(L,KC-1,NS),0.)
        ENDDO
        DO K=2,KS
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            TVAR1W(L,K)=DELTI*HPK(L,K)+MAX(WSETA(L,K-1,NS),0.)-MIN(WSETA(L,K,NS),0.)
          ENDDO
        ENDDO

        ! *** TVAR1E=RIGHT HAND SIDE
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TVAR1E(L,K)=DELTI*HPK(L,K)*SED(L,K,NS)
          ENDDO
        ENDDO

        ! *** TVAR3S=BET,TVAR2N=U,TVAR2S=GAM ARE WORKING ARRAYS
        DO LP=LF,LL
          L=LWET(LP)  
          TVAR3S(L)=TVAR1W(L,KSZ(L))
        ENDDO
        DO LP=LF,LL
          L=LWET(LP)  
          TVAR2N(L,KSZ(L))=TVAR1E(L,KSZ(L))/TVAR3S(L)
        ENDDO
        DO K=2,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TVAR2S(L,K)=TVAR1N(L,K-1)/TVAR3S(L)
            TVAR3S(L)=TVAR1W(L,K)-TVAR1S(L,K)*TVAR2S(L,K)
            TVAR2N(L,K)=(TVAR1E(L,K)-TVAR1S(L,K)*TVAR2N(L,K-1))/TVAR3S(L)
          ENDDO
        ENDDO
        DO K=KS,1,-1
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)
          ENDDO
        ENDDO
        DO K=1,KC
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SED(L,K,NS)=TVAR2N(L,K)
          ENDDO
        ENDDO
      ENDIF   ! *** END ANTI-DIFFUSION

      !----------------------------------------------------------------------!
      !
      ! **  FINAL FLUX
      IF( KC > 1 .AND. ( ISTRAN(5) > 0 .OR. ISBAL > 0 )  )THEN
        ! *** KC-1 LAYER
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND)  
          SEDF(L,KS,NS)=DELTI*DZC(L,KC)*HP(L)*(SED(L,KC,NS)-SEDS(L,KC,NS))
        ENDDO 

        ! *** MIDDLE LAYERS
        DO K=KS-1,1,-1
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SEDF(L,K,NS)=DELTI*DZC(L,K+1)*HP(L)*(SED(L,K+1,NS)-SEDS(L,K+1,NS))+SEDF(L,K+1,NS)
          ENDDO  
        ENDDO  
      ENDIF
      
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
    
  ENDDO     ! *** END OF SEDIMENT CLASS LOOP

  !**********************************************************************!
  IF( ISDTXBUG > 0 )THEN
    IFLAG = 0  
    DO NS=1,NSED  
      DO L=2,LA
        K = KSZ(L)
        IF( SED(L,K,NS) < 0. )THEN
          IF( IFLAG == 0 )THEN  
            OPEN(1,FILE=OUTDIR//'NEGSEDSND.OUT',POSITION='APPEND')  
            IFLAG = 1  
          ENDIF  
          WRITE(1,107) TIME,NS,IL(L),JL(L),K,SED(L,K,NS)  
          SED(L,K,NS) = 0.0    ! *** Continue with warning
        ENDIF  
      ENDDO  
    ENDDO  

    DO NS=1,NSED
      DO L=2,LA
        K = KBT(L)
        IF( SEDB(L,K,NS) < 0. .OR. HBED(L,K) < 0. )THEN
          IF( IFLAG == 0 )THEN  
            OPEN(1,FILE=OUTDIR//'NEGSEDSND.OUT',POSITION='APPEND')  
            IFLAG=1  
          ENDIF  
          WRITE(1,108) TIME,NS,IL(L),JL(L),K,SEDB(L,K,NS),SEDF(L,0,NS),HBED(L,K)
          IF( SEDB(L,K,NS) < 0. ) SEDB(L,K,NS) = 0.0    ! *** Continue with warning
          IF( HBED(L,K) < 0. )    HBED(L,K) = 0.0       ! *** Continue with warning
        ENDIF  
      ENDDO  
    ENDDO  

    IF( IFLAG == 1 )CLOSE(1) 
  ENDIF
  
  107 FORMAT(' Warning: WC  SED < 0: TIME, NS, I, J, K, NEGSED = ',F12.4,4I5,4E13.4)       
  108 FORMAT(' Warning: BED SED < 0: TIME, NS, I, J, K, NEGSEDB, SEDF, HBED = ',F12.4,4I5,4E13.4)       

  !**********************************************************************!

  RETURN

END

