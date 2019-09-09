SUBROUTINE CALSND  

  !**********************************************************************!
  ! *** SUBROUTINE CALSND CALCULATES NONCOHESIVER SEDIMENT SETTLING,  
  ! *** DEPOSITION AND RESUSPENSION AND IS CALLED FOR SSEDTOX  
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-03-02        PAUL M. CRAIG    REWRITTEN TO CODE F90
  !                                     REMOVED KC DEPENDENT DUPLCIATE CODE
  !                                     ADDED OMP
  ! 2010-XX-XX        SCOTT JAMES      ADDED MHK
  !
  !**********************************************************************!
  
  USE GLOBAL  

  IMPLICIT NONE

  INTEGER :: NX,NS,K,L,LP,ISGPFLAG,ISEHFLAG,KTOP,NXX,NSS,IFLAG,LE,LN,ND,LF,LL,L1,NSB,LUTMP,LDTMP
  REAL    :: TIME,GRADSED,SIGP,CRNUM,DUM1,DUM3,DUM4
  REAL    :: ZEQMAX,CSHIELDS,TMPVAL
  REAL    :: SHIELDS,TOP,BOT,WSFAC,WESE,FLUXFAC
  REAL    :: WSETMP,WVEL,CLEFT,CRIGHT,SNDBTMP,SEDAVG
  REAL, EXTERNAL :: FSEDMODE,CSNDZEQ,CSNDSET,CSNDEQC

  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: SEDPHI

  ! *** FIRST CALL ALLOCATIONS AND ASSIGNMENTS
  IF(  .NOT. ALLOCATED(SEDPHI) )THEN
    ALLOCATE(SEDPHI(NSTM))

    ! *** GP - CONVERT SEDIMENT DIAMETERS IN M TO MM AND SET PHI SIZE  
    DO NX=1,NSND  
      NS=NSED+NX  
      SEDPHI(NS)=-LOG(1000.*SEDDIA(NS))/LOG(2.)  
    ENDDO  

    ! *** IF ISNDAL == 0 THEN NO ARMORING, CONSTANTS FOR ENTIRE SIMULATION
    DO K=1,KB  
      DO L=2,LA  
        SIGPHI(L,K)=0.  
      ENDDO  
    ENDDO  

    DO NX=1,NSND  
      DO L=2,LA  
        PEXP(L,NX)=1.  
        PHID(L,NX)=1.  
      ENDDO  
    ENDDO  

    ! *** SET MAXIMUM NONCOHESIVE SEDIMENT DIAMETER  
    DO NX=1,NSND  
      NS=NSED+NX  
      SNDDMX=MAX(SNDDMX,SEDDIA(NS))  
    ENDDO  
  ENDIF
  
  !**********************************************************************!  
  IF( ISDYNSTP == 0 )THEN  
    TIME=(DT*FLOAT(N)+TCON*TBEGIN)/TCON  
  ELSE  
    TIME=TIMESEC/TCON  
  ENDIF  

  !**********************************************************************!  
  ! *** SET TIME/SPATIALLY VARYING NONCOHESIVE ARMORING PARAMETERS  
  IF( ISNDAL >= 1 )THEN  
    ! *** 1 ACTIVATE NON-COHESIVE ARMORING EFFECTS (GARCIA & PARKER)
    ! *** 2 SAME AS 1 WITH ACTIVE-PARENT LAYER FORMULATION
    ISGPFLAG=0  
    ISEHFLAG=0  
    DO NX=1,NSND  
      NS=NX+NSED  
      IF( ISNDEQ(NS) == 1 ) ISGPFLAG=1  
      IF( ISBDLD(NS) >= 2 ) ISEHFLAG=1  
    ENDDO  

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,K,NS,NX,KTOP,NXX,NSS)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      ! *** SET SIGPHI FOR GARCIA AND PARKER (1991) EQS 46 AND 47  
      IF( ISGPFLAG == 1 )THEN  

        ! *** GP - SET MEAN PHI FOR TOP LAYER OF BED  
        DO LP=LF,LL
          L=LWET(LP)
          TVAR3W(L)=0.  
          TVAR3E(L)=0.  
          SIGPHI(L,KBT(L))=0.  
        ENDDO  

        DO NX=1,NSND  
          NS=NSED+NX  
          DO LP=LF,LL
            L=LWET(LP)
            KTOP=KBT(L)  
            TVAR3W(L)=TVAR3W(L)+SEDPHI(NS)*VFRBED(L,KTOP,NS)  
            TVAR3E(L)=TVAR3E(L)+VFRBED(L,KTOP,NS)  
          ENDDO  
        ENDDO  

        DO LP=LF,LL
          L=LWET(LP)
          IF( TVAR3E(L) <= 0.) TVAR3E(L)=1.  
          TVAR3W(L)=TVAR3W(L)/TVAR3E(L)  
        ENDDO  

        DO NX=1,NSND  
          NS=NSED+NX  
          DO LP=LF,LL
            L=LWET(LP)
            KTOP=KBT(L)  
            SIGPHI(L,KTOP) = SIGPHI(L,KTOP)+((SEDPHI(NS)-TVAR3W(L))**2)*VFRBED(L,KTOP,NS)/TVAR3E(L)
          ENDDO  
        ENDDO  

        DO LP=LF,LL
          L=LWET(LP)
          KTOP=KBT(L)
          IF( SIGPHI(L,KTOP) < 0. )THEN  
            SIGPHI(L,KTOP)=-SQRT(ABS(SIGPHI(L,KTOP)))   ! *** PMC TEMP PATCH  
          ELSE
            SIGPHI(L,KTOP)=SQRT(SIGPHI(L,KTOP))  
          ENDIF
        ENDDO  

      ENDIF    ! *** END CALCULATION OF SIGPHI FOR GARCIA AND PARKER (1991)  

      ! *** SET EXPOSURE AND HIDING FUNCTIONS FOR ENGULAND-HANSEN AND WU,WANG,  
      IF( ISEHFLAG == 1 .OR. ISGPFLAG == 1 )THEN  

        DO NX=1,NSND  
          DO LP=LF,LL
            L=LWET(LP)
            PEXP(L,NX)=0.0  
            PHID(L,NX)=0.0  
          ENDDO  
        ENDDO  

        DO LP=LF,LL
          L=LWET(LP)
          SNDBT(L,KBT(L))=0.0  
        ENDDO  

        DO NX=1,NSND  
          NS=NSED+NX  
          DO LP=LF,LL
            L=LWET(LP)
            K=KBT(L)  
            SNDBT(L,K)=SNDBT(L,K)+SNDB(L,K,NX)  
            DO NXX=1,NSND  
              NSS=NSED+NXX  
              PEXP(L,NX)=PEXP(L,NX)+SNDB(L,K,NXX)*SEDDIA(NS) /(SEDDIA(NS)+SEDDIA(NSS))  
              PHID(L,NX)=PHID(L,NX)+SNDB(L,K,NXX)*SEDDIA(NSS)/(SEDDIA(NS)+SEDDIA(NSS))  
            ENDDO  
          ENDDO  
        ENDDO  

        DO NX=1,NSND  
          DO LP=LF,LL
            L=LWET(LP)
            K=KBT(L)  
            IF( SNDBT(L,K) > 0.0 )THEN  
              PEXP(L,NX)=PEXP(L,NX)/SNDBT(L,K)  
              PHID(L,NX)=PHID(L,NX)/SNDBT(L,K)  
            ELSE  
              PEXP(L,NX)=1.0  
              PHID(L,NX)=1.0  
            END IF  
          ENDDO  
        ENDDO  
      ENDIF   ! *** END SET EXPOSURE AND HIDING FUNCTIONS FOR ENGULAND-HANSEN  

    ENDDO  ! *** END OF DOMAIN
    !$OMP END PARALLEL DO 

  ENDIF     ! *** END OF ARMORING SECTION

  !**********************************************************************!  
  ! ***  SET CRITICAL SHIELED'S PARAMETER FOR D50  
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,DUM1,DUM3,DUM4)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    DO LP=LF,LL
      L=LWET(LP)
      CALL SETSHLD(DUM1,CSHIELDS50(L),SEDDIA50(L,KBT(L)),SSG(NSED+1),DUM3,DUM4)  
    ENDDO  
  ENDDO
  !$OMP END PARALLEL DO 
  
  !**********************************************************************!  
  ! *** NONCOHESIVE SEDIMENT SCOUR, DEPOSITION AND VERTICAL PROCESSES 
  DO NX=1,NSND  
    NS=NX+NSED  

    DSEDGMM  = 1./(1.E6*SSG(NS))      ! *** SPECIFIC VOLUME (M**3/G)
    DIASED   = SEDDIA(NS)  
    GPDIASED = G*(SSG(NS)-1.)*DIASED  

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,K,LP,L,KTOP) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      !----------------------------------------------------------------------!  
      ! *** SET SETTLING VELOCITIES  
      IF( ISNDVW == 0 )THEN  
        DO K=0,KS  
          DO LP=1,LLWET(K+1,ND)
            L=LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = WSEDO(NS)  
          ENDDO  
        ENDDO  
      ENDIF  

      IF( ISNDVW >= 1 )THEN  
        DO K=0,KS  
          DO LP=1,LLWET(K+1,ND)
            L=LKWET(LP,K+1,ND)  
            WSETA(L,K,NS) = WSEDO(NS)*CSNDSET(SNDT(L,K+1),SDEN(NS),ISNDVW)  
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
        DO LP=1,LLWET(K+1,ND)
          L=LKWET(LP,K+1,ND)  
          WSETA(L,K,NS) = SPB(L)*WSETA(L,K,NS)
        ENDDO
      ENDDO

      
      ! *** SET THE ROUSE PARAMETER
      IF( IROUSE(NX) == 0 )THEN  
        DO LP=LF,LL
          L=LWET(LP)
          IF( USTAR(L) > 1E-12 )THEN  
            ROUSE(L) = WSETA(L,0,NS)/(VKC*USTAR(L))  
          ELSE  
            ROUSE(L) = 250000.*WSETA(L,0,NS)  
          END IF  
        ENDDO  
      ELSE  
        DO LP=LF,LL
          L=LWET(LP)
          IF( USTARSND(L) > 1E-12 )THEN  
            ROUSE(L) = WSETA(L,0,NS)/(VKC*USTARSND(L))  
          ELSE  
            ROUSE(L) = 250000.*WSETA(L,0,NS)  
          END IF  
        ENDDO  
      ENDIF  

    ENDDO  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO 
    
    !----------------------------------------------------------------------!  
    ! *** COMPUTE THE BEDLOAD COMPONENT FOR THE CURRENT SEDIMENT CLASS
    IF( ISBDLDBC > 0 ) CALL BEDLOAD(NX,NS)  

    !----------------------------------------------------------------------!  
    ! *** HORIZONTAL LOOP  
    
    !$OMP PARALLEL DO DEFAULT(NONE) PRIVATE(ND,LF,LL,LP,L,K,KTOP,LE,LN) &
    !$OMP             PRIVATE(WVEL,CLEFT,CRIGHT,WSFAC,WESE,ZEQMAX,SIGP,CSHIELDS,SHIELDS,TMPVAL,SNDBTMP,WSETMP)  &
    !$OMP             PRIVATE(TOP,BOT,FLUXFAC,CRNUM,GRADSED,SEDAVG,NSB,LUTMP,LDTMP,L1)    &
    !$OMP  SHARED(NDM,LDMWET,LAWET,LLWET,LKWET,LWET,KC,SNDF,SND,SND1,SNDS,WSETA,USTAR,USTARSND,RSNDM,ISNDM1,ISNDM2,NX,NS,KS) &
    !$OMP  SHARED(DIASED,DSEDGMM,GPDIASED,DELTI,SNDB,SNDB1,KBT,SNDFBL,TAUR,TAUBSND,SEDDIA50,HP,ISNDEQ,SSG,FACSUSL,CSHIELDS50) &
    !$OMP  SHARED(DZC,ZEQ,ZEQD,ZEQDI,SIGPHI,VDRBED,ISNDAL,ISLTAUC,TCSHIELDS,ISEDEFF,PHID,PEXP,VFRBED,VDRDEPO,HPK,HPKI,DZIC)  &
    !$OMP  SHARED(ROUSE,COEHEFF2,COEHEFF,FRACCOH,SNDEQB,SNDEQ,SNDEQSAV,DTSED,LEC,LNC,S2TL,S3TL,N,IBMECH,SEDVRDT) &
    !$OMP  SHARED(QSBDLDX,QSBDLDY,QSBDLDOT,DXYIP,QSBDLDIN,QSBDTOP,QWBDTOP,ISTOPT,LMASKDRY,LSC,LWC) &
    !$OMP  SHARED(TVAR1S,TVAR2S,TVAR3S,TVAR1N,TVAR2N,TVAR1W,TVAR1E,KSZ,ISBDLDBC,NSBDLDBC,LSBLBCU,LSBLBCD)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      
      ! *** SET FLUX FOR THE BOTTOM OF THE TOP LAYER
      IF( KC > 1 )THEN
        K=KC  
        DO LP=LF,LL
          L=LWET(LP)
          SNDF(L,K,NX)=0.                                !  (G/M2/S)
          WVEL=DTSED*HPKI(L,K)  
          CLEFT=1.+WSETA(L,K-1,NS)*WVEL  
          CRIGHT=MAX(SND(L,K,NX),0.)  
          SND(L,K,NX)=CRIGHT/CLEFT  
          SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)  
        ENDDO  

        ! *** SET FLUX FOR MIDDLE LAYERS
        DO K=KS,2,-1  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            WVEL=DTSED*HPKI(L,K)  
            CLEFT=1.+WSETA(L,K-1,NS)*WVEL  
            CRIGHT=MAX(SND(L,K,NX),0.)-SNDF(L,K,NX)*WVEL  
            SND(L,K,NX)=CRIGHT/CLEFT  
            SNDF(L,K-1,NX)=-WSETA(L,K-1,NS)*SND(L,K,NX)  
          ENDDO  
        ENDDO  
      ENDIF
      
      ! *** FSEDMODE SETS BEDLOAD (IMODE=1) AND SUSPENDED LOAD (IMODE=2) TRANSPORT FRACTIONS
      DO LP=LF,LL
        L=LWET(LP)
        FACSUSL(L) = FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),RSNDM(NX),ISNDM1(NX),ISNDM2(NX),2)  
      ENDDO  

      ! *** UPDATE SEDIMENT BED MASS & BOTTOM LAYER CONCENTRATION
      WSFAC=1.0  ! *** FOR KC>1
      DO LP=LF,LL
        L=LWET(LP)
        WESE=0.
        
        IF( ISNDEQ(NS) == 0 )THEN
          ! *** USER SPECIFIED EQUILIBRIUM CONCENTRATION
          ZEQ(L) = 0.05         ! *** GARCIA & PARKER
          ZEQD(L) = ZEQ(L)*DZIC(L,KSZ(L))  
          ZEQDI(L) = 1./ZEQD(L)  
          SNDEQB(L) = TAUR(NS)
        ELSE
          ! *** SET EQUILIBRUIM CONCENTRATION REFERENCE HEIGHT (DIMENSIONLESS)  
          ZEQ(L) = CSNDZEQ(ISNDEQ(NS),DIASED,GPDIASED,TAUR(NS),TAUBSND(L),SEDDIA50(L,KBT(L)),HP(L),SSG(NS),WSETA(L,0,NS))
          
          ZEQMAX=0.5*DZC(L,KSZ(L))  
          ZEQ(L)=MIN(ZEQ(L),ZEQMAX)  
          ZEQD(L)=ZEQ(L)*DZIC(L,KSZ(L))  
          ZEQDI(L)=1./ZEQD(L)  
          SIGP=SIGPHI(L,KBT(L))  

          ! *** SET EQUILIBRUIM CONCENTRATION  
          SNDEQB(L) = CSNDEQC(ISNDEQ(NS),DIASED,SSG(NS),WSETA(L,0,NS),TAUR(NS),TAUBSND(L),SEDDIA50(L,KBT(L)),SIGP,ZEQ(L),VDRBED(L,KBT(L)),ISNDAL)      !  

          ! *** APPLIED LIMITOR TO GARCIA AND PARKER  
          IF( ISNDEQ(NS) == 1 )THEN  
            IF( ISLTAUC(NS) == 1 )THEN  
              CSHIELDS=TCSHIELDS(NS)  
              IF( ISEDEFF == 2 )THEN  
                TMPVAL=1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
                CSHIELDS=TMPVAL*CSHIELDS  
              ENDIF  
              SHIELDS=TAUBSND(L)/GPDIASED  
              IF( SHIELDS < CSHIELDS) SNDEQB(L)=0.  
            ENDIF  
            IF( ISLTAUC(NS) == 2 )THEN  
              CSHIELDS=SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED  
              IF( ISEDEFF == 2 )THEN  
                TMPVAL=1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
                CSHIELDS=TMPVAL*CSHIELDS  
              ENDIF  
              SHIELDS=TAUBSND(L)/GPDIASED  
              IF( SHIELDS < CSHIELDS) SNDEQB(L)=0.  
            ENDIF  
            IF( ISLTAUC(NS) == 3 )THEN  
              CSHIELDS=0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)  
              IF( ISEDEFF == 2 )THEN  
                TMPVAL=1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )  
                CSHIELDS=TMPVAL*CSHIELDS  
              ENDIF  
              SHIELDS=TAUBSND(L)/GPDIASED  
              IF( SHIELDS < CSHIELDS) SNDEQB(L)=0.  
            ENDIF  
          ENDIF  

          IF( ISEDEFF == 1 ) SNDEQB(L) = SNDEQB(L)*EXP(-COEHEFF*FRACCOH(L,KBT(L)))  
        ENDIF    ! *** ISNDEQ(NS) == 0
        
        IF( ROUSE(L) < 0.999 .OR. ROUSE(L) > 1.001 )THEN  
          TOP = (ZEQD(L)**(ROUSE(L)-1.))-1.  
          BOT = (1.-ROUSE(L))*(ZEQDI(L)-1.)  
          SNDEQ(L) = SNDEQB(L)*TOP/BOT  
          SNDEQ(L) = FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)  
          SNDEQSAV(L,NX) = SNDEQ(L)  
        ELSE  
          TOP=LOG(ZEQDI(L))  
          BOT=(ZEQDI(L)-1.)  
          SNDEQ(L)=SNDEQB(L)*TOP/BOT  
          SNDEQ(L)=FACSUSL(L)*VFRBED(L,KBT(L),NS)*MAX(SNDEQ(L),0.)  
          SNDEQSAV(L,NX)=SNDEQ(L)  
        ENDIF  

        ! *** SET RESUSPENSION FLUX
        IF( KSZ(L) == KC )THEN    ! *** Alberta
          WSFAC=2.*(1.+ROUSE(L))/(2.+ROUSE(L)*(1.-ZEQ(L)))  
        ENDIF
        WESE = WSFAC*WSETA(L,0,NS)*SNDEQ(L)  

        ! *** SET DEPOSITION VELOCITY  
        WSETMP = WSFAC*WSETA(L,0,NS)  
        WVEL   = DTSED*HPKI(L,KSZ(L))
        
        ! *** SET BOTTOM LAYER SND CONCENTRATION AND FLUX TO SUSPENDED LOAD
        CLEFT  = 1. + WSETMP*WVEL  
        CRIGHT = MAX(SND(L,KSZ(L),NX),0.) + ( WESE-SNDF(L,KSZ(L),NX) )*WVEL  
        SND(L,KSZ(L),NX) = CRIGHT/CLEFT  
        SNDF(L,0,NX) = -WSETMP*SND(L,KSZ(L),NX) + WESE  

        ! *** ADD BED LOAD FLUX TO SUSPENDED LOAD FLUX
        SNDBTMP = SNDB(L,KBT(L),NX) - DTSED*SNDF(L,0,NX) - DTSED*SNDFBL(L,NX)  

        ! *** HANDLE CASE IF INSUFFICIENT BED MATERIAL
        IF( SNDBTMP < 0.0 )THEN  
          ! *** ADJUST BEDLOAD FLUX, TRYING TO KEEP SUSPENDED LOAD FLUX CONSTANT
          SNDBTMP = SNDB(L,KBT(L),NX) - DTSED*SNDF(L,0,NX) 
          LE=LEC(L)
          LN=LNC(L)
          IF( SNDBTMP < 0.0 )THEN
            IF( ISBDLDBC > 0 )THEN
              ! *** ZERO BEDLOAD OUTFLUX
              QSBDLDX(L,NX) = MAX(QSBDLDX(L,NX),0.0)
              QSBDLDY(L,NX) = MAX(QSBDLDY(L,NX),0.0)
              QSBDLDX(LE,NX) = MIN(QSBDLDX(LE,NX),0.0)
              QSBDLDY(LN,NX) = MIN(QSBDLDY(LN,NX),0.0)
              IF( NSBDLDBC > 0 )THEN
                DO NSB=1,NSBDLDBC
                  LUTMP = LSBLBCU(NSB)
                  LDTMP = LSBLBCD(NSB)
                  IF( L == LUTMP )THEN
                    IF( LDTMP > 0 )QSBDLDIN(LDTMP,NX) = QSBDLDIN(LDTMP,NX) - QSBDLDOT(LUTMP,NX)
                    QSBDLDOT(LUTMP,NX) = 0.
                    EXIT
                  ENDIF
                ENDDO
              ENDIF
              ! *** UPDATE BEDLOAD FLUX FOR CURRENT CELL
              SNDFBL(L,NX) = DXYIP(L)*( QSBDLDX(LE,NX)-QSBDLDX(L,NX) + QSBDLDY(LN,NX)-QSBDLDY(L,NX) + QSBDLDOT(L,NX)-QSBDLDIN(L,NX) )  
            ENDIF
            !print *,'ZERO BEDLOAD OUTFLUX',n,L
            
            ! *** REDUCE SUSPENDED LOAD FLUX
            SNDF(L,0,NX)     = SNDB(L,KBT(L),NX)/DTSED
            SND(L,KSZ(L),NX) = SNDS(L,KSZ(L),NX)+(SNDF(L,0,NX)-SNDF(L,KSZ(L),NX))*WVEL  
          ELSE
            ! *** REDUCE BEDLOAD FLUX
            SNDBTMP = SNDBTMP/DTSED
            FLUXFAC = SNDBTMP/SNDFBL(L,NX)

            IF( QSBDLDX(L,NX) < 0. ) QSBDLDX(L,NX) = FLUXFAC*QSBDLDX(L,NX)  
            IF( QSBDLDY(L,NX) < 0. ) QSBDLDY(L,NX) = FLUXFAC*QSBDLDY(L,NX)  
            IF( NSBDLDBC > 0 )THEN
              DO NSB=1,NSBDLDBC
                LUTMP = LSBLBCU(NSB)
                LDTMP = LSBLBCD(NSB)
                IF( L == LUTMP )THEN
                  IF( LDTMP > 0 )QSBDLDIN(LDTMP,NX) = QSBDLDIN(LDTMP,NX) - QSBDLDOT(LUTMP,NX)
                  QSBDLDOT(L,NX) = FLUXFAC*QSBDLDOT(L,NX) 
                  EXIT
                ENDIF
              ENDDO
            ENDIF
            !print *,'REDUCE BEDLOAD FLUX',n,L
            ! *** UPDATE BEDLOAD FLUX FOR CURRENT CELL
            SNDFBL(L,NX) = DXYIP(L)*( QSBDLDX(LE,NX)-QSBDLDX(L,NX) + QSBDLDY(LN,NX)-QSBDLDY(L,NX) + QSBDLDOT(L,NX)-QSBDLDIN(L,NX) )  
          ENDIF
          
          ! *** UPDATE BEDLOAD FLUX FOR SURROUNDING CELLS
          IF( ISBDLDBC > 0 )THEN
            L1 = LEC(L)
            LE=LEC(L1)
            LN=LNC(L1)
            SNDFBL(L1,NX) = DXYIP(L1)*( QSBDLDX(LE,NX)-QSBDLDX(L1,NX) + QSBDLDY(LN,NX)-QSBDLDY(L1,NX) + QSBDLDOT(L1,NX)-QSBDLDIN(L1,NX) )  
            L1 = LWC(L)
            LE=LEC(L1)
            LN=LNC(L1)
            SNDFBL(L1,NX) = DXYIP(L1)*( QSBDLDX(LE,NX)-QSBDLDX(L1,NX) + QSBDLDY(LN,NX)-QSBDLDY(L1,NX) + QSBDLDOT(L1,NX)-QSBDLDIN(L1,NX) )  
            L1 = LSC(L)
            LE=LEC(L1)
            LN=LNC(L1)
            SNDFBL(L1,NX) = DXYIP(L1)*( QSBDLDX(LE,NX)-QSBDLDX(L1,NX) + QSBDLDY(LN,NX)-QSBDLDY(L1,NX) + QSBDLDOT(L1,NX)-QSBDLDIN(L1,NX) )  
            L1 = LNC(L)
            LE=LEC(L1)
            LN=LNC(L1)
            SNDFBL(L1,NX) = DXYIP(L1)*( QSBDLDX(LE,NX)-QSBDLDX(L1,NX) + QSBDLDY(LN,NX)-QSBDLDY(L1,NX) + QSBDLDOT(L1,NX)-QSBDLDIN(L1,NX) )  
          ENDIF
          SNDBTMP=0.0
        ENDIF
            
        SNDB1(L,KBT(L),NX) = S3TL*SNDB(L,KBT(L),NX) + S2TL*SNDB1(L,KBT(L),NX)  
        SNDB(L,KBT(L),NX)  = SNDBTMP  
        SNDF(L,0,NX) = SNDF(L,0,NX) + SNDFBL(L,NX)       ! *** BED/WATER INTEFACE SND FLUX COMBINES BEDLOAD AND SUSPENSION
        QSBDTOP(L) = QSBDTOP(L) + DSEDGMM*SNDF(L,0,NX) 
        IF( IBMECH == 1 .AND. SEDVRDT < 0.00001 )THEN
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*VDRBED(L,KBT(L))*SNDF(L,0,NX)  ! *** IF EITHER CONSTANT OR INSTANTLY CONSOLIDATING, USE BED VR  
        ELSE
          QWBDTOP(L) = QWBDTOP(L) + DSEDGMM*( VDRBED(L,KBT(L))*MAX(SNDF(L,0,NX),0.) + VDRDEPO(NS)*MIN(SNDF(L,0,NX),0.) )
        ENDIF
      ENDDO  

      !----------------------------------------------------------------------!  
      ! **  ANTI-DIFFUSION OF NON-COHESIVE SEDIMENT  KC > 1
      IF( ISTOPT(7) == 1 .AND. KC > 1 )THEN
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            CRNUM=1.+DTSED*WSETA(L,K,NS)*HPKI(L,K+1)
            GRADSED=(SND(L,K+1,NX)-SND(L,K,NX))/(DZC(L,K+1)+DZC(L,K))  
            SEDAVG=0.5*(SND(L,K+1,NX)-SND(L,K,NX)+1.E-16)  
            WSETA(L,K,NS)=-CRNUM*DZC(L,K+1)*WSETA(L,K,NS)*GRADSED/SEDAVG  
          ENDDO  
        ENDDO  

        ! *** TVAR1S=LOWER DIAGONAL  
        DO LP=LF,LL
          L=LWET(LP)
          TVAR1S(L,KSZ(L))=0  
        ENDDO  
        DO K=2,KC  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
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
          TVAR1W(L,KSZ(L))=DELTI*HPK(L,KSZ(L))-MIN(WSETA(L,KSZ(L),NS),0.)
          TVAR1W(L,KC)=DELTI*DZC(L,KC)*HP(L)+MAX(WSETA(L,KC-1,NS),0.)  
        ENDDO  
        DO K=2,KS  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            TVAR1W(L,K)=DELTI*DZC(L,KC)*HP(L)+MAX(WSETA(L,K-1,NS),0.)-MIN(WSETA(L,K,NS),0.)  
          ENDDO  
        ENDDO  

        ! *** TVAR1E=RIGHT HAND SIDE  
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            TVAR1E(L,K)=DELTI*DZC(L,KC)*HP(L)*SND(L,K,NX)  
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
          DO L=LF,LL  
            TVAR2N(L,K)=TVAR2N(L,K)-TVAR2S(L,K+1)*TVAR2N(L,K+1)  
          ENDDO  
        ENDDO  
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SND(L,K,NX)=TVAR2N(L,K)  
          ENDDO  
        ENDDO  

      ENDIF   ! *** END OF ANTI DIFFUSION
      
      !----------------------------------------------------------------------!  
      ! *** FINAL FLUX  
      IF( KC > 1 )THEN
        DO LP=1,LLWET(KS,ND)
          L=LKWET(LP,KS,ND)  
          SNDF(L,KS,NX) = DELTI*DZC(L,KC)*HP(L)*(SND(L,KC,NX)-SNDS(L,KC,NX))  
        ENDDO  

        DO K=KS-1,1,-1  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            SNDF(L,K,NX) = DELTI*DZC(L,K+1)*HP(L)*(SND(L,K+1,NX)-SNDS(L,K+1,NX))+SNDF(L,K+1,NX)  
          ENDDO  
        ENDDO  
      ENDIF
            
    ENDDO  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO 
    
  ENDDO  ! *** END OF NON-COHESIVE CLASS LOOP

  !**********************************************************************!  
  IF( ISDTXBUG > 0 )THEN
    IFLAG=0  
    DO NS=1,NSND  
      DO L=2,LA  
        K = KSZ(L)
        IF( SND(L,K,NS) < 0. )THEN  
          IF( IFLAG == 0 )THEN  
            OPEN(1,FILE=OUTDIR//'NEGSEDSND.OUT',POSITION='APPEND')  
            IFLAG=1  
          ENDIF  
          WRITE(1,107) TIME,NS,IL(L),JL(L),K,SND(L,K,NS)  
          SND(L,K,NS) = 0.0    ! *** Continue with warning
        ENDIF  
      ENDDO  
    ENDDO  

    DO NS=1,NSND  
      DO L=2,LA
        K = KBT(L)
        IF( SNDB(L,K,NS) < 0. .OR. HBED(L,K) < 0. .OR. SEDDIA50(L,K) < 0.0 )THEN
          IF( IFLAG == 0 )THEN  
            OPEN(1,FILE=OUTDIR//'NEGSEDSND.OUT',POSITION='APPEND')  
            IFLAG=1  
          ENDIF  
          WRITE(1,108) TIME,NS,IL(L),JL(L),K,SNDB(L,K,NS),SNDF(L,0,NS),HBED(L,K),SEDDIA50(L,K)*1000000.
          IF( SNDB(L,K,NS) < 0. )  SNDB(L,K,NS) = 0.0    ! *** Continue with warning
          IF( HBED(L,K) < 0. )     HBED(L,K) = 0.0       ! *** Continue with warning
          IF( SEDDIA50(L,K) < 0. ) SEDDIA50(L,K) = 0.0   ! *** Continue with warning
        ENDIF  
      ENDDO  
    ENDDO  

    IF( IFLAG == 1 )CLOSE(1) 
  ENDIF
  
  107 FORMAT(' Warning: WC  SND < 0: TIME, NS, I, J, K, NEGSND = ',F12.4,4I5,4E13.4)       
  108 FORMAT(' Warning: BED SND < 0: TIME, NS, I, J, K, NEGSNDB, SNDF, HBED, D50 = ',F12.4,4I5,4E13.4)       

  !**********************************************************************!  
  RETURN  

END  
