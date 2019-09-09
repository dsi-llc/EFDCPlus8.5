SUBROUTINE CALUVW(ISTL_)  
  !**********************************************************************C
  !
  ! *** CALCULATE THE INTERNAL SOLUTION AT TIME LEVEL (N+1)  
  ! *** THE VALUE OF ISTL INDICATES THE NUMBER OF TIME LEVELS IN THE STEP  
  !
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP

  USE GLOBAL  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: ISTL_
  INTEGER :: NMD,ND,LF,LL,L,LS,LN,K,LE,LW,LNN,LP,LHOST,LCHNU,LCHNV
  INTEGER :: ICFL,JCFL,KCFL,IVAL,IDTCFL,LTMP,L1P,IOBC,IFILE
  INTEGER,SAVE :: NKCE,NKCN,NHOLE,NSTEP
  
  REAL   :: Q1,Q2,CMU,CMV,CRU,CRV,EU,EV,RCDZM,RCDZU,RCDZL,HPPTMP,TMPVALN,DTMAXX,HPDEL
  REAL   :: RLAMN,RLAMO,TMPVAL,CFLUUUT,CFLVVVT,CFLWWWT,CFLCACT,DTCFL,UWTMP,UETMP,VSTMP,VNTMP,WBTMP,WTTMP,DTMAXI
  REAL   :: STOKESU,STOKESX,STOKESY,STOKESW,C1,C2,THETA,AMPL,SINHH
  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LKCE
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LKCN
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LHOLE
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:) :: LSTEP

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: UHDYEE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: VHDXEE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RCXX
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: RCYY
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TVARE
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TVARN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TVARW
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TVARS
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: CERRU
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: CERRV
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SINH2
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: SINH4

  IF( .NOT. ALLOCATED(RCXX) )THEN
    ALLOCATE(LKCE(LCM))  
    ALLOCATE(LKCN(LCM))  
    ALLOCATE(LHOLE(LCM))  
    ALLOCATE(LSTEP(LCM))  
    ALLOCATE(UHDYEE(LCM))  
    ALLOCATE(VHDXEE(LCM))  
    ALLOCATE(RCXX(LCM))  
    ALLOCATE(RCYY(LCM))
    ALLOCATE(TVARE(LCM))  
    ALLOCATE(TVARN(LCM))
    ALLOCATE(TVARW(LCM))  
    ALLOCATE(TVARS(LCM))
    ALLOCATE(CERRU(LCM,KCM))
    ALLOCATE(CERRV(LCM,KCM))

    IF( ISWAVE > 0 .AND. ISWVSD >= 1 )THEN 
      ALLOCATE(SINH2(LCM))  
      ALLOCATE(SINH4(LCM))  
      SINH2 = 0.
      SINH4 = 0.
    ENDIF
    LKCE=0
    LKCN=0
    LHOLE=0
    LSTEP=0
    UHDYEE=0.0
    VHDXEE=0.0
    RCXX=0.0
    RCYY=0.0
    TVARE=0.
    TVARN=0.
    TVARW=0.
    TVARS=0.
    CERRU=0.
    CERRV=0.
    
    DO L=2,LA
      ! *** CERRU
      IF (KSZU(L) < KC) THEN
        Q1 = 0.
        DO K=KSZU(L),KC
          Q1 = Q1 + 1./(1.-SGZU(L,K))
        ENDDO
        IF( Q1 > 0. )THEN
          DO K=KSZU(L),KC
            CERRU(L,K) = (1./(1.-SGZU(L,K)))/Q1
            CERRU(L,K) = CERRU(L,K)*REAL(KC-KSZU(L)+1)
          ENDDO
        ENDIF
      ELSE
        CERRU(L,KC) = 1
      ENDIF
      
      ! *** CERRV
      IF (KSZV(L) < KC) THEN
        Q1 = 0.
        DO K=KSZV(L),KC
          Q1 = Q1 + 1./(1.-SGZV(L,K))
        ENDDO
        IF( Q1 > 0. )THEN
          DO K=KSZV(L),KC
            CERRV(L,K) = (1./(1.-SGZV(L,K)))/Q1
            CERRV(L,K) = CERRV(L,K)*REAL(KC-KSZV(L)+1)
          ENDDO
        ENDIF
      ELSE
        CERRV(L,KC) = 1
      ENDIF      
    ENDDO
    
    ! *** HANDLE KSZ E&N = KC CASES
    NKCE = 0
    NKCN = 0
    DO L=2,LA
      IF( KSZ(LWC(L)) < KC .AND. SUBO(L) > 0.5 .AND. KSZ(L) == KC )THEN
        NKCE = NKCE+1
        LKCE(NKCE) = L 
      ENDIF
      IF( KSZ(LSC(L)) < KC .AND. SVBO(L) > 0.5 .AND. KSZ(L) == KC  )THEN
        NKCN = NKCN+1
        LKCN(NKCN) = L 
      ENDIF
    ENDDO  

    ! *** BUILD A LIST OF ALL THE ISOLATED HOLES IN BATHYMETRY
    NHOLE = 0
    DO L=2,LA
      DO K=KSZ(L),KC
        Q1 = SGZU(L,K) + SGZU(LEC(L),K) + SGZV(L,K) + SGZV(LNC(L),K)
        IF( Q1 == 0. )THEN
          NHOLE= NHOLE+1
          LHOLE(NHOLE) = L
          EXIT
        ENDIF
      ENDDO
    ENDDO
    
    ! *** BUILD A LIST OF ALL THE STEPS BETWEEN KC AND KS BOTTOM LAYERS
    NSTEP = 0
    IF( IGRIDV > 0 .AND. KMINV == 1 )THEN
      DO L=2,LA
        IF( KSZ(L) == KS )THEN
          K = KS
          Q1 = SGZU(L,K) + SGZU(LEC(L),K) + SGZV(L,K) + SGZV(LNC(L),K)
          Q2 = 0.
      
          ! *** EXCLUDE HOLES
          IF( Q1 /= 0. )THEN
        
            ! *** ONE OR MORE U FACES ARE ACTIVE
            IF( ( SUBO(LEC(L)) > 0. .AND. KSZ(LEC(L)) == KC ) .OR. ( SUBO(L) > 0. .AND. KSZ(LWC(L)) == KC ) )THEN
              Q2 = Q2 + 1.
            ENDIF
            
            ! *** ONE OR MORE V FACES ARE ACTIVE
            IF( ( SVBO(LNC(L)) > 0. .AND. KSZ(LNC(L)) == KC ) .OR. ( SVBO(L) > 0. .AND. KSZ(LSC(L)) == KC ) )THEN
              Q2 = Q2 + 1.
            ENDIF
            IF( Q2 > 0. )THEN
              NSTEP = NSTEP+1
              LSTEP(NSTEP) = L
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDIF
  
  IF( ISDYNSTP == 0 )THEN  
    DELT=DT2  
    IF( ISTL_ == 2 )THEN  
      DELT=DT  
    ENDIF  
    DELTI=1./DELT  
  ELSE  
    DELT=DTDYN  
    DELTI=1./DELT  
  ENDIF  
  IFILE = -1
  
  ! *** ZERO NEWLY DRY CELLS
  IF( LADRY > 0 )THEN
    DO LP=1,LADRY
      L=LDRY(LP)
      AAU(L)=0.  
      AAV(L)=0.  
      BBU(L)=1.  
      BBV(L)=1.
      RCXX(L)=0.
      RCYY(L)=0.
      TVARE(L)=0.
      TVARN(L)=0.
      TVARW(L)=0.
      TVARS(L)=0.
    ENDDO

    DO K=1,KC
      DO LP=1,LADRY
        L=LDRY(LP) 
        DU(L,K)=0.
        DV(L,K)=0.
        UUU(L,K)=0.
        VVV(L,K)=0.
        UHDY(L,K)=0. 
        VHDX(L,K)=0.
        UHDYF(L,K)=0. 
        VHDXF(L,K)=0.
        U(L,K)=0.
        V(L,K)=0.
        W(L,K)=0.
        UHDY1(L,K)=0. 
        VHDX1(L,K)=0.
        UHDYF1(L,K)=0. 
        VHDXF1(L,K)=0.
        U1(L,K)=0.
        V1(L,K)=0.
        W1(L,K)=0.

        ! *** ZERO THE TRANSPORT FLUXES
        UHDY2(L,K)=0. 
        VHDX2(L,K)=0.
        UHDYF2(L,K)=0. 
        VHDXF2(L,K)=0.
        U2(L,K)=0.
        V2(L,K)=0.
        W2(L,K)=0.
        CU1(L,K)=0.
        CU2(L,K)=0.

        ! *** Adjacent Cells
        LE=LEC(L)
        LN=LNC(L)
        UHDY(LE,K)=0. 
        VHDX(LN,K)=0.
        UHDYF(LE,K)=0. 
        VHDXF(LN,K)=0.
        V(LN,K)=0.
        U(LE,K)=0.
        UHDY1(LE,K)=0. 
        VHDX1(LN,K)=0.
        UHDYF1(LE,K)=0. 
        VHDXF1(LN,K)=0.
        U1(LE,K)=0.
        V1(LN,K)=0.
        UHDY2(LE,K)=0. 
        VHDX2(LN,K)=0.
        UHDYF2(LE,K)=0. 
        VHDXF2(LN,K)=0.
        U2(LE,K)=0.
        V2(LN,K)=0.
      ENDDO
    ENDDO  
  ENDIF
  
  IF( KC == 1 ) GOTO 30  

  IF( IGRIDV > 0 .AND. KMINV == 1 )THEN
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,LASGZ1,LDMSGZ1,LSGZ1,KSZ,UHDYE,UHDYF,UHDY,HUI,DYIU,U,VHDXE,VHDXF,VHDX,HVI,DXIV,V,W) PRIVATE(ND,LF,LL,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMSGZ1+1  
      LL=MIN(LF+LDMSGZ1-1,LASGZ1)
      
      DO LP=LF,LL
        L = LSGZ1(LP)
        UHDYF(L,KSZ(L)) = UHDYE(L)  
        UHDY(L,KSZ(L))  = UHDYE(L)  
        U(L,KSZ(L))     = UHDYE(L)*HUI(L)*DYIU(L)
        VHDXF(L,KSZ(L)) = VHDXE(L)  
        VHDX(L,KSZ(L))  = VHDXE(L)  
        V(L,KSZ(L))     = VHDXE(L)*HVI(L)*DXIV(L)
        W(L,KSZ(L))     = 0.  
      ENDDO  
    ENDDO
    !$OMP END PARALLEL DO
  ENDIF
  
  !$OMP PARALLEL DEFAULT(SHARED)

  ! ***************************************************************************
  ! *** CALCULATE BOTTOM FRICTION COEFFICIENT  
  IF( ISTL_ == 3 )THEN  
    !$OMP DO PRIVATE(ND,LP,L) 
    DO ND=1,NDM  
      DO LP=1,LLWETZ(KC,ND)
        L=LKWETZ(LP,KC,ND)  
        RCXX(L) = AVCON1/H1U(L) + STBX(L)*SQRT( V1U(L)*V1U(L) + U1(L,KSZU(L))*U1(L,KSZU(L)) )  
        RCYY(L) = AVCON1/H1V(L) + STBY(L)*SQRT( U1V(L)*U1V(L) + V1(L,KSZV(L))*V1(L,KSZV(L)) )  
      ENDDO 
    ENDDO
    !$OMP END DO 
  ELSE  
    IF( AVCON1 < 0.00001 )THEN  
      ! *** FOR 2TL U1 & U AND V1 & V ARE THE SAME
      ! *** THESE ARE ONLY DIFFERENCE FOR 3TL ISTL=2 TRAP CORRECTION STEP

      !$OMP DO PRIVATE(ND,LP,L,Q1,Q2) 
      DO ND=1,NDM  
        DO LP=1,LLWETZ(KC,ND)
          L = LKWETZ(LP,KC,ND)  
          Q1      = SQRT( U1(L,KSZU(L))*U1(L,KSZU(L)) + V1U(L)*V1U(L) )  
          Q2      = SQRT( U(L,KSZU(L)) *U(L,KSZU(L))  + VU(L) *VU(L) )  
          RCXX(L) = STBX(L)*SQRT(Q1*Q2)  
          Q1      = SQRT( V1(L,KSZV(L))*V1(L,KSZV(L)) + U1V(L)*U1V(L) )  
          Q2      = SQRT( V(L,KSZV(L)) *V(L,KSZV(L))  + UV(L) *UV(L) )  
          RCYY(L) = STBY(L)*SQRT(Q1*Q2)  
        ENDDO  
      ENDDO
      !$OMP END DO 
    ELSE  
      !$OMP DO PRIVATE(ND,LP,L,Q1,Q2) 
      DO ND=1,NDM  
        DO LP=1,LLWETZ(KC,ND)
          L=LKWETZ(LP,KC,ND)  
          Q1      = SQRT( U1(L,KSZU(L))*U1(L,KSZU(L)) + V1U(L)*V1U(L) )  
          Q2      = SQRT( U(L,KSZU(L)) *U(L,KSZU(L))  + VU(L) *VU(L) )  
          RCXX(L) = AVCON1/SQRT(H1U(L)*HU(L))+STBX(L)*SQRT(Q1*Q2)  
          Q1      = SQRT( V1(L,KSZV(L))*V1(L,KSZV(L)) + U1V(L)*U1V(L) )  
          Q2      = SQRT( V(L, KSZV(L))*V(L, KSZV(L)) + UV(L) *UV(L) )  
          RCYY(L) = AVCON1/SQRT(H1V(L)*HV(L))+STBY(L)*SQRT(Q1*Q2)  
        ENDDO  
      ENDDO
      !$OMP END DO 
    ENDIF  
  ENDIF  
  
  ! ***************************************************************************
  ! *** SPLIT THE DEPTH AVERAGED FLOWS INTO THE LAYER SPECIFIC FLOWS 
  !$OMP DO PRIVATE(ND,K,LP,L,LE,LN) 
  DO ND=1,NDM  
    DO K=1,KS
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)  
        LN=LNC(L)
        LE=LEC(L)
        UHDYEK(L,K) = DZC(L,K)*( UHDYE(L) - UHDYE(LE) )
        VHDXEK(L,K) = DZC(L,K)*( VHDXE(L) - VHDXE(LN) )
      ENDDO
    ENDDO
  ENDDO
  !$OMP END DO

  ! *** MAKE ADJUSTMENTS FOR SIGMA-ZED STEPS AND HOLES
  IF( IGRIDV > 0 )THEN
    !$OMP SINGLE
    DO LP=1,NHOLE
      L = LHOLE(LP) 
      LN=LNC(L)
      LE=LEC(L)
      UHDYEE(L) = UHDYE(L)
      VHDXEE(L) = VHDXE(L)
      UHDYEE(LE) = UHDYE(LE)
      VHDXEE(LN) = VHDXE(LN)
    ENDDO
  
    DO LP=1,NHOLE
      L = LHOLE(LP) 
      LN=LNC(L)
      LE=LEC(L)
      Q2 = 0.
      DO K=KSZ(L),KS
        Q1 = SGZU(L,K)+SGZU(LE,K)+SGZV(L,K)+SGZV(LN,K)
        IF( Q1 == 0. )THEN
          ! *** ISOLATED HOLE IN BATHYMETRY
            
          ! *** DUE TO MACHINE PRECISION, THE HP-H1P APPROACH CAN PRODUCE FLUX IMBALANCES
          !TMPVAL = DZC(L,K)*(HP(L)-H1P(L))*DXYP(L)/DELT - DZC(L,K)*QSUME(L)
          TMPVAL = -0.5*DZC(L,K)*( UHDYE(LE)+UHDY1E(LE) - UHDYE(L)-UHDY1E(L) + VHDXE(LN)+VHDX1E(LN) - VHDXE(L)-VHDX1E(L) - QSUME(L) ) 
          Q2 = Q2 + DZC(L,K)
          
          UWTMP = TMPVAL*DXP(L)/(DXP(L)+DYP(L))
          UHDYEK(L,K) = UWTMP
          UHDY1EK(L,K) = UWTMP
          UHDYEE(L) = UHDYEE(L) - UWTMP
          
          VSTMP = TMPVAL*DYP(L)/(DXP(L)+DYP(L))
          VHDXEK(L,K) = VSTMP
          VHDX1EK(L,K) = VSTMP
          VHDXEE(L) = VHDXEE(L) - VSTMP
        ELSE
          UHDYEK(L,K) = DZC(L,K)/(1.-Q2)*( UHDYEE(L) - UHDYEE(LE) )
          VHDXEK(L,K) = DZC(L,K)/(1.-Q2)*( VHDXEE(L) - VHDXEE(LN) )
        ENDIF
      ENDDO
    ENDDO
    
    ! *** HANDLE KS-KC STEPS
    IF( KMINV == 1 )THEN
      DO LP=1,NSTEP
        L = LSTEP(LP) 
        LN=LNC(L)
        LE=LEC(L)
        
        K = KS
        HPDEL = HP(L)-H1P(L)
        Q1 = DZC(L,K)*HPDEL*DXYP(L)/DELT - DZC(L,K)*QSUME(L)                           ! *** FLUX REQUIRED DUE TO HPK CHANGE
        Q2 = -0.5*DZC(L,K)*( UHDYE(LE)+UHDY1E(LE) - UHDYE(L)-UHDY1E(L) + &
                             VHDXE(LN)+VHDX1E(LN) - VHDXE(L)-VHDX1E(L) - QSUME(L) )    ! *** FLUX DUE TO HORIZONTAL IN/OUT FLOWS
        IF( HPDEL > 0. )THEN
          UHDYEK(L,K) = MAX(Q1,Q2)
        ELSE
          UHDYEK(L,K) = MIN(Q1,Q2)
        ENDIF
        VHDXEK(L,K) = 0.0
      ENDDO
    ENDIF
    !$OMP END SINGLE
  ENDIF

  ! ***************************************************************************
  ! *** CALCULATE THE U AND V SHEARS  
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LE,LW,LS,LN,RCDZL,RCDZM,RCDZU,CMU,EU,CMV,EV,CRU,CRV) 
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    DO K=1,KS  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        DU(L,K) = SUB3D(L,K)*DU(L,K)
        IF( K > KSZU(L) )THEN  
          RCDZL = CDZLU(L,K)
          RCDZM = CDZMU(L,K)*DELTI  
          RCDZU = CDZUU(L,K)  
          CMU = 1.+RCDZM*HU(L)*AVUI(L,K) 
          EU  = 1./(CMU-RCDZL*CU1(L,K-1))  
          CU1(L,K) = RCDZU*EU  
          DU(L,K)  = (DU(L,K)-RCDZL*DU(L,K-1))*EU  
          UUU(L,K) = -RCDZL*UUU(L,K-1)*EU 
        ELSEIF( K == KSZU(L) )THEN
          ! *** K = KSZU(L)
          RCDZL = CDZLU(L,K)  
          RCDZM = CDZMU(L,K)*DELTI  
          RCDZU = CDZUU(L,K)  
          CMU = 1.+RCDZM*HU(L)*AVUI(L,K)
          EU  = 1./CMU  
          CU1(L,K) = RCDZU*EU  
          DU(L,K)  = (DU(L,K)-RCDZL*RCXX(L)*UHE(L)*HUI(L))*EU
          UUU(L,K) = EU
        ELSE
          CU1(L,K)=0.
          DU(L,K)=0.
          UUU(L,K)=0.
        ENDIF

        DV(L,K) = SVB3D(L,K)*DV(L,K)
        IF( K > KSZV(L) )THEN  
          RCDZL = CDZLV(L,K)
          RCDZM = CDZMV(L,K)*DELTI  
          RCDZU = CDZUV(L,K)  
          CMV = 1.+RCDZM*HV(L)*AVVI(L,K) 
          EV  = 1./(CMV-RCDZL*CU2(L,K-1))  
          CU2(L,K) = RCDZU*EV  
          DV(L,K)  = (DV(L,K)-RCDZL*DV(L,K-1))*EV  
          VVV(L,K) = -RCDZL*VVV(L,K-1)*EV  
        ELSEIF( K == KSZV(L) )THEN
          ! *** K = KSZV(L)
          RCDZL = CDZLV(L,K)  
          RCDZM = CDZMV(L,K)*DELTI  
          RCDZU = CDZUV(L,K)  
          CMV = 1.+RCDZM*HV(L)*AVVI(L,K)
          EV  = 1./CMV  
          CU2(L,K) = RCDZU*EV  
          DV(L,K)  = (DV(L,K)-RCDZL*RCYY(L)*VHE(L)*HVI(L))*EV  
          VVV(L,K) = EV 
        ELSE 
          CU2(L,K)=0.  
          DV(L,K)=0.
          VVV(L,K)=0.
        ENDIF
      ENDDO  
    ENDDO
    
    ! ***  BACK SUBSTITUTION
    DO K=KS-1,1,-1  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        DU(L,K) = DU(L,K)-CU1(L,K)*DU(L,K+1)  
        DV(L,K) = DV(L,K)-CU2(L,K)*DV(L,K+1)  
        UUU(L,K) = UUU(L,K)-CU1(L,K)*UUU(L,K+1)  
        VVV(L,K) = VVV(L,K)-CU2(L,K)*VVV(L,K+1)  
      ENDDO  
    ENDDO  

    ! *** SHERMAN-MORRISON BACK SUBSTITUTION
    DO LP=1,LLWETZ(KC,ND)
      L=LKWETZ(LP,KC,ND)  
      AAU(L)=0.  
      AAV(L)=0.  
      BBU(L)=1.  
      BBV(L)=1.  
    ENDDO
    
    DO K=1,KS  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        CRU = CDZRU(L,K)*RCXX(L)*AVUI(L,K)
        AAU(L) = AAU(L) + CRU*DU(L,K)  
        BBU(L) = BBU(L) + CRU*UUU(L,K)  
        CRV = CDZRV(L,K)*RCYY(L)*AVVI(L,K)  
        AAV(L) = AAV(L) + CRV*DV(L,K)  
        BBV(L) = BBV(L) + CRV*VVV(L,K)  
      ENDDO  
    ENDDO  
    DO LP=1,LLWETZ(KC,ND)
      L=LKWETZ(LP,KC,ND)  
      AAU(L) = AAU(L)/BBU(L)  
      AAV(L) = AAV(L)/BBV(L)  
    ENDDO

    DO K=1,KS  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        DU(L,K) = SUB3D(L,K)*DZGU(L,K)*HU(L)*AVUI(L,K)*( DU(L,K) - AAU(L)*UUU(L,K) ) 
        DV(L,K) = SVB3D(L,K)*DZGV(L,K)*HV(L)*AVVI(L,K)*( DV(L,K) - AAV(L)*VVV(L,K) )  
      ENDDO  
    ENDDO  
    
  ENDDO    ! ***  END OF DOMAIN LOOP
  !$OMP END DO 

  ! ***************************************************************************
  ! *** CALCULATE U AND V USING THE INTERNAL SHEARS DU AND DV (DU/DV ARE IN M2/S)
  ! *** DUSUM+UHE=UHE, DVSUM+VHE=VHE  
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L) 
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    ! *** ADJUST FLOWS BASED ON INTERNAL SHEARS
    
    ! *** INTERIM: COMPUTE SUM OF UNIT WIDTH Q PLUS UNIT FLOWS BY LAYER
    DO K=1,KS  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        UHE(L) = UHE(L) + CDZDU(L,K)*DU(L,K)
        VHE(L) = VHE(L) + CDZDV(L,K)*DV(L,K)
      ENDDO  
    ENDDO  

    ! *** LAYER SPECIFIC DEPTH AVERAGED FLOWS (M3/S)
    DO LP=1,LLWETZ(KC,ND)
      L=LKWETZ(LP,KC,ND)  
      UHDYF(L,KC) = UHE(L)*SUB(L) 
      VHDXF(L,KC) = VHE(L)*SVB(L) 
    ENDDO  
    DO K=KS,1,-1  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        UHDYF(L,K) = SUB3D(L,K)*(UHDYF(L,K+1) - DU(L,K))
        VHDXF(L,K) = SVB3D(L,K)*(VHDXF(L,K+1) - DV(L,K))
      ENDDO  
    ENDDO  

    ! *** COMPUTE FLOWS IN M3/S FROM UNIT WIDTH Q
    DO K=1,KC  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        UHDYF(L,K) = UHDYF(L,K)*DYU(L)  
        VHDXF(L,K) = VHDXF(L,K)*DXV(L)
      ENDDO  
    ENDDO  
  ENDDO    ! ***  END OF DOMAIN LOOP
  !$OMP END DO 
  
  ! ***************************************************************************
  ! *** ADD ADJUSTMENT TO 3D HORIZONTAL TRANSPORT  
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K) 
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    DO LP=1,LLWETZ(KC,ND)
      L=LKWETZ(LP,KC,ND)  
      TVARE(L)=0.  
      TVARN(L)=0.  
    ENDDO
    
    ! *** SCALE TO DEPTH AVERAGE FLOW
    DO K=1,KC  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        TVARE(L) = TVARE(L) + UHDYF(L,K)*SGZU(L,K)
        TVARN(L) = TVARN(L) + VHDXF(L,K)*SGZV(L,K)
      ENDDO  
    ENDDO  
    
    ! *** COMPUTE DIFFERENCE FROM EXTERNAL SOLUTION
    DO LP=1,LLWETZ(KC,ND)
      L=LKWETZ(LP,KC,ND)  
      TVARE(L) = TVARE(L)-UHDYE(L)  
      TVARN(L) = TVARN(L)-VHDXE(L)  
    ENDDO  
    
    ! *** CORRECT INTERNAL SOLUTION LAYER SPECIFIC FLOWS
    DO K=1,KC  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        UHDYF(L,K) = UHDYF(L,K)-TVARE(L)*CERRU(L,K)
        UHDY(L,K)  = UHDYF(L,K)*SGZU(L,K)               ! *** LAYER SPECIFIC FLOWS
        VHDXF(L,K) = VHDXF(L,K)-TVARN(L)*CERRV(L,K)
        VHDX(L,K)  = VHDXF(L,K)*SGZV(L,K)               ! *** LAYER SPECIFIC FLOWS
      ENDDO
    ENDDO  
  ENDDO    ! ***  END OF DOMAIN LOOP
  !$OMP END DO 

  ! ***************************************************************************
  ! *** RESET VELOCITIES  
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,K) 
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)
    
    DO LP=1,LLWETZ(KC,ND)
      L=LKWETZ(LP,KC,ND)  
      UHE(L)=0.  
      VHE(L)=0.  
    ENDDO  
    DO K=1,KC  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        UHE(L) = UHE(L) + UHDYF(L,K)*SGZU(L,K)
        VHE(L) = VHE(L) + VHDXF(L,K)*SGZV(L,K)
        U(L,K) = UHDYF(L,K)*HUI(L)  
        V(L,K) = VHDXF(L,K)*HVI(L)  
      ENDDO  
    ENDDO  
    DO K=1,KC  
      DO LP=1,LLWETZ(K,ND)
        L=LKWETZ(LP,K,ND)  
        U(L,K) = U(L,K)*DYIU(L)  
        V(L,K) = V(L,K)*DXIV(L)  
      ENDDO  
    ENDDO  
    DO LP=1,LLWETZ(KC,ND)
      L=LKWETZ(LP,KC,ND)  
      UHE(L) = UHE(L)*DYIU(L)  
      VHE(L) = VHE(L)*DXIV(L)  
    ENDDO  
  ENDDO    ! ***  END OF DOMAIN LOOP
  !$OMP END DO 

  ! ***************************************************************************
  ! *** CALCULATE W  
  IF( ISTL_ == 3 )THEN  
    !$OMP DO PRIVATE(ND,K,LP,L,LN,LE) 
    DO ND=1,NDM  
      DO K=1,KS
        DO LP=1,LLWETZ(K,ND)
          L=LKWETZ(LP,K,ND)  
          LN=LNC(L)
          LE=LEC(L)
          W(L,K) = W(L,K-1) - W2(L,K) + W2(L,K-1)  &
                   - DXYIP(L)*( UHDY(LE,K)-UHDY(L,K) + UHDYEK(L,K) + UHDY2(LE,K)-UHDY2(L,K) + UHDY2EK(L,K)   &
                              + VHDX(LN,K)-VHDX(L,K) + VHDXEK(L,K) + VHDX2(LN,K)-VHDX2(L,K) + VHDX2EK(L,K) ) &
                   + 2.*( QSUM(L,K)-DZC(L,K)*QSUME(L) )*DXYIP(L)
        ENDDO
      ENDDO
    ENDDO    ! ***  END OF DOMAIN LOOP
    !$OMP END DO 
 
  ELSEIF( ISTL_ == 2 )THEN  
    ! *** TWO TIME LEVEL SOLUTION OR 3TL CORRECTOR
    !$OMP DO PRIVATE(ND,K,LP,L,LN,LE) 
    DO ND=1,NDM  
      DO K=1,KS  
        DO LP=1,LLWETZ(K,ND)
          L=LKWETZ(LP,K,ND)  
          LN=LNC(L)
          LE=LEC(L)  
          W(L,K) = W(L,K-1) - 0.5*DXYIP(L)  &
                             *( UHDY(LE,K)-UHDY(L,K) + UHDYEK(L,K) + UHDY1(LE,K)-UHDY1(L,K) + UHDY1EK(L,K)   &
                              + VHDX(LN,K)-VHDX(L,K) + VHDXEK(L,K) + VHDX1(LN,K)-VHDX1(L,K) + VHDX1EK(L,K) ) &
                            + ( QSUM(L,K) - DZC(L,K)*QSUME(L) )*DXYIP(L)
         ENDDO  
      ENDDO  
    ENDDO    ! ***  END OF DOMAIN LOOP
    !$OMP END DO 
  ENDIF  
  !$OMP END PARALLEL
  
  ! *** APPLY OPEN BOUNDARYS 
  DO LL=1,NBCSOP
    L=LOBCS(LL)
    DO K=1,KS  
      W(L,K)=0.0
    ENDDO  
  ENDDO 
  
  ! ***************************************************************************
  ! *** JUMP TO POINT FOR KC=1
  30 CONTINUE   
  
  ! ***************************************************************************
  ! *** CALCULATE U AND V ON OPEN BOUNDARIES  
  DO K=1,KC  
    DO LL=1,NCBS  
      IF( ISPBS(LL) /= 2 )THEN
        L=LCBS(LL)  
        LN=LNC(L)  
        LNN=LNC(LN)  
        IF( LN /= LC .AND. K >= KSZ(L) )THEN  
          VHDXF(LN,K) = VHDXF(LNN,K)-VHDXE(LNN)+VHDXE(LN)  
          VHDX(LN,K)  = VHDXF(LN,K)*SGZV(LN,K)
          V(LN,K)     = VHDXF(LN,K)/(HV(LN)*DXV(LN))
          W(LN,K)     = 0.
        ELSE  
          W(LN,K)     = 0.
        ENDIF
      ENDIF
    ENDDO  
  ENDDO 
 
  DO K=1,KC  
    DO LL=1,NCBW
      IF( ISPBW(LL) /= 2 )THEN
        L=LCBW(LL)  
        LE=LEC(L)  
        L1P=LEC(LE)  
        IF( LE /= LC .AND. K >= KSZ(L) )THEN  
          UHDYF(LE,K) = UHDYF(L1P,K)-UHDYE(L1P)+UHDYE(LE)  
          UHDY(LE,K)  = UHDYF(LE,K)*SGZU(LE,K)
          U(LE,K)     = UHDYF(LE,K)/(HU(LE)*DYU(LE))  
          W(LE,K)     = 0.
        ELSE  
          W(LE,K)     = 0.
        ENDIF  
      ENDIF
    ENDDO  
  ENDDO
  
  DO K=1,KC  
    DO LL=1,NCBE
      L=LCBE(LL)
      LW=LWC(L)
      IF( ISPBE(LL) /= 2 .AND. K >= KSZ(L) )THEN  
        UHDYF(L,K) = UHDYF(LW,K)-UHDYE(LW)+UHDYE(L)  
        UHDY(L,K)  = UHDYF(L,K)*SGZU(L,K)
        U(L,K)     = UHDYF(L,K)/(HU(L)*DYU(L))  
        W(L,K)     = 0.
      ELSE
        W(L,K)     = 0.
      ENDIF
    ENDDO  
  ENDDO
  
  DO K=1,KC  
    DO LL=1,NCBN
      L=LCBN(LL)  
      IF( ISPBN(LL) /= 2 .AND. K >= KSZ(L) )THEN 
        LS=LSC(L)  
        VHDXF(L,K) = VHDXF(LS,K)-VHDXE(LS)+VHDXE(L)  
        VHDX(L,K)  = VHDXF(L,K)*SGZV(L,K)
        V(L,K)     = VHDXF(L,K)/(HV(L)*DXV(L))  
        W(L,K)     = 0.
      ELSE
        W(L,K)     = 0.
      ENDIF
    ENDDO  
  ENDDO  

  !$OMP PARALLEL DEFAULT(SHARED)

  ! ***************************************************************************
  ! *** CALCULATE AVERAGE CELL FACE TRANSPORTS FOR SALT, TEMPERATURE AND  
  ! *** SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2 
  IF( ISTL_ == 2 )THEN  
    ! *** 3TL ISTL=2 OR IS2TIM>0 
    !$OMP DO PRIVATE(ND,K,LP,L,LN) 
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND) 
          UHDYF2(L,K) = 0.5*(UHDYF(L,K)+UHDYF1(L,K))
          VHDXF2(L,K) = 0.5*(VHDXF(L,K)+VHDXF1(L,K))
          UHDY2(L,K) = 0.5*(UHDY(L,K)+UHDY1(L,K))
          VHDX2(L,K) = 0.5*(VHDX(L,K)+VHDX1(L,K))  
          U2(L,K) = 0.5*(U(L,K)+U1(L,K))  
          V2(L,K) = 0.5*(V(L,K)+V1(L,K))  
          W2(L,K) = 0.5*(W(L,K)+W1(L,K))  
        ENDDO  
      ENDDO  
    ENDDO 
    !$OMP END DO

  ELSE  

    ! *** 3TL ISTL=3
    !$OMP DO PRIVATE(ND,K,LP,L) 
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          UHDYF2(L,K) = 0.5*(UHDYF(L,K)+UHDYF2(L,K))
          VHDXF2(L,K) = 0.5*(VHDXF(L,K)+VHDXF2(L,K))
          UHDY2(L,K) = 0.5*(UHDY(L,K)+UHDY2(L,K))  
          VHDX2(L,K) = 0.5*(VHDX(L,K)+VHDX2(L,K))  
          U2(L,K) = 0.5*(U(L,K)+U2(L,K))  
          V2(L,K) = 0.5*(V(L,K)+V2(L,K))  
          W2(L,K) = 0.5*(W(L,K)+W2(L,K))  
        ENDDO  
      ENDDO  
    ENDDO 
    !$OMP END DO
  ENDIF  
  
  ! ***************************************************************************
  ! *** ADDITIONAL 3D CONTINUITY ADJUSTED ADDED BELOW  
  IF( KC > 1 )THEN  
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      
      DO LP=1,LLWETZ(KC,ND)
        L = LKWETZ(LP,KC,ND)  
        TVARE(L)=0.  
        TVARN(L)=0.  
      ENDDO  
      DO K=1,KC  
        DO LP=1,LLWETZ(K,ND)
          L = LKWETZ(LP,K,ND)  
          TVARE(L) = TVARE(L)+UHDY2(L,K)
          TVARN(L) = TVARN(L)+VHDX2(L,K)
        ENDDO  
      ENDDO
    ENDDO
    !$OMP END DO

    ! *** HANDLE KSZ EAST = KC CASES
    !$OMP SINGLE
    DO LP=1,NKCE
      L = LKCE(LP)  
      TVARE(L) = UHDY2(L,KC)
    ENDDO
    
    ! *** HANDLE KSZ NORTH = KC CASES
    DO LP=1,NKCN
      L = LKCN(LP)  
      TVARN(L) = VHDX2(L,KC)
    ENDDO  
    !$OMP END SINGLE
      
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LE,LW,LS,LN,K,IOBC,HPPTMP) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      
      IF( ISTL_ == 3 )THEN  
        ! *** 3TL AND ISTL=3
        DO LP=1,LLWETZ(KC,ND)
          L=LKWETZ(LP,KC,ND)  
          LE=LEC(L)
          LN=LNC(L)
          HPPTMP = H2P(L) + DELT*DXYIP(L)*( QSUME(L) - TVARE(LE)+TVARE(L)-TVARN(LN)+TVARN(L) )  
          IF( ISGWIE >= 1 ) HPPTMP = HPPTMP - DELT*DXYIP(L)*(EVAPSW(L)-QGW(L))

          ! *** NEGATIVE DEPTH CHECK
          IF( HPPTMP <= 0. .AND. HP(L) > 0. )THEN
            ! *** CHECK OPEN BC
            LW = 0
            DO IOBC=1,NBCSOP
              LW = LOBCS(IOBC)
              IF( L == LW )EXIT
            ENDDO  
            IF( L == LW )CYCLE

            LW=LWC(L)
            LS=LSC(L)
            IF( HP(L) <= HDRY .AND. ISDRY > 0 )THEN
              PRINT '(A,I8,I5,3F10.4,F12.4)',' WARNING! NEG DEPTH IN CONTINUITY CHECK: N,L,H2P,HP,HPNEW,TIMEDAY',N,L,H2P(L),HP(L),HPPTMP,TIMEDAY
              IF( IFILE == -1 )THEN
                IFILE = 8
                OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
              ENDIF
              WRITE(8,'(A,I5,3F10.4,F12.4)')' WARNING!  NEG DEPTH IN CONTINUITY CHECK:'
              WRITE(8,'(A,I10,F12.4,I5,3F10.4)')'            N,TIMEDAY,L,H2P,HP,HPNEW',N,TIMEDAY,L,H2P(L),HP(L),HPPTMP
              
              HPPTMP = HP(L) 
              UHDYE(L)=0.
              VHDXE(L)=0.
              DO K=1,KC
                UHDYF(L,K)=0.
                VHDXF(L,K)=0.
                UHDY(L,K)=0.
                VHDX(L,K)=0.
                U(L,K)=0.
                V(L,K)=0.
                W(L,K)=0.

                UHDYF(LE,K)=0.
                VHDXF(LN,K)=0.
                UHDY(LE,K)=0.
                VHDX(LN,K)=0.
                U(LE,K)=0.
                V(LN,K)=0.

                UHDYF2(L,K)=0.
                VHDXF2(L,K)=0.
                UHDY2(L,K)=0.
                VHDX2(L,K)=0.
                U2(L,K)=0.
                V2(L,K)=0.
                W2(L,K)=0.
                UHDY2(LE,K)=0.
                VHDX2(LN,K)=0.
                U2(LE,K)=0.
                V2(LN,K)=0.
              ENDDO
            ELSE
              PRINT '(A,I8,I5,3F10.4,F12.4)',' ERROR! NEG DEPTH IN CONTINUITY CHECK. N,L,H2P,HP,HPNEW,TIMEDAY',N,L,H2P(L),HP(L),HPPTMP,TIMEDAY
              IF( IFILE == -1 )THEN
                IFILE = 8
                OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
              ENDIF
              WRITE(8,'(A)')              ' ERROR!  NEGATIVE DEPTH IN CONTINUITY CHECK.'
              WRITE(8,'(A,I15,3I5,F15.5)')'         N L ISTL NCTBC TIMEDAY   ',N,L,ISTL_,NCTBC,TIMEDAY
              WRITE(8,'(A,4F10.4)')       '         H1P H2P HP HPNEW         ',H1P(L),H2P(L),HP(L),HPPTMP
              WRITE(8,'(A,4F10.4)')       '         HP  WESN                 ',HP(LW),   HP(LE),   HP(LS),    HP(LN)
              WRITE(8,'(A,4F10.4)')       '         H1P WESN                 ',H1P(LW),  H1P(LE),  H1P(LS),   H1P(LN)
              WRITE(8,'(A,4F10.4)')       '         H2P WESN                 ',H2P(LW),  H2P(LE),  H2P(LS),   H2P(LN)
              WRITE(8,'(A,4F10.4)')       '         SUB/SVB                  ',SUB(L),   SUB(LE),   SVB(L),   SVB(LN)
              WRITE(8,'(A,4E12.4)')       '         UHDYE/VHDXE              ',UHDYE(L), UHDYE(LE), VHDXE(L), VHDXE(LN)
              WRITE(8,'(A,4E12.4)')       '         UHDYE1/VHDXE1            ',UHDY1E(L),UHDY1E(LE),VHDX1E(L),VHDX1E(LN)
              WRITE(8,'(A,4E12.4)')       '         UHDYE2/VHDXE2            ',UHDY2E(L),UHDY2E(LE),VHDX2E(L),VHDX2E(LN)
              WRITE(8,'(A,4E12.4,L5)')    '         WESN FLOWS SUM LAYERS    ',TVARE(L), TVARE(LE), TVARN(L), TVARN(LN), LMASKDRY(L)
              DO K=KC,KSZ(L),-1  
                WRITE(8,'(A,I5,4E12.4)')  '         LAYER FLOWS/CURRENT      ',K,UHDY(L,K),UHDY(LE,K),VHDX(L,K),VHDX(LN,K)
              ENDDO
              DO K=KC,KSZ(L),-1  
                WRITE(8,'(A,I5,4E12.4)')  '         LAYER FLOWS/OLD          ',K,UHDY1(L,K),UHDY1(LE,K),VHDX1(L,K),VHDX1(LN,K)
              ENDDO
              DO K=KC,KSZ(L),-1  
                WRITE(8,'(A,I5,4E12.4)')  '         SUB3D/SVB3D              ',K,SUB3D(L,K),SUB3D(LE,K),SVB3D(L,K),SVB3D(LN,K)
              ENDDO
            ENDIF
          ENDIF

          HP(L)  = HPPTMP   !SPB(L)*HPPTMP + (1.-SPB(L))*(GI*P(L)-BELV(L))  
          HPI(L) = 1./HP(L)  
        ENDDO
        
      ELSE
        ! *** 2TL / ISTL=2
        DO LP=1,LLWETZ(KC,ND)
          L=LKWETZ(LP,KC,ND)  
          LE=LEC(L)
          LN=LNC(L)
          HPPTMP = H1P(L) + DELT*DXYIP(L)*( QSUME(L) + TVARE(L)-TVARE(LE) + TVARN(L)-TVARN(LN) )  
          IF( ISGWIE >= 1 ) HPPTMP = HPPTMP - DELT*DXYIP(L)*(EVAPSW(L)-QGW(L))
          
          ! *** NEGATIVE DEPTH CHECK
          IF( HPPTMP <= 0. .AND. HP(L) > 0. )THEN
            ! *** CHECK OPEN BC
            LW = 0
            DO IOBC=1,NBCSOP
              LW = LOBCS(IOBC)
              IF( L == LW )EXIT
            ENDDO  
            IF( L == LW )CYCLE
            
            LW=LWC(L)
            LS=LSC(L)
            IF( HP(L) <= HDRY .AND. ISDRY > 0 )THEN
              PRINT '(A,I8,I5,3F10.4,F12.4)',' WARNING!  NEG DEPTH IN CONTINUITY CHECK. N,L,H1P,HP,HPNEW,TIMEDAY',N,L,H1P(L),HP(L),HPPTMP,TIMEDAY
              IF( IFILE == -1 )THEN
                IFILE = 8
                OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
              ENDIF
              WRITE(8,'(A,I5,3F10.4,F12.4)')    ' WARNING!  NEG DEPTH IN CONTINUITY CHECK:'
              WRITE(8,'(A,I10,F12.4,I5,3F10.4)')'                N,TIMEDAY,L,H1P,HP,HPNEW',N,TIMEDAY,L,H1P(L),HP(L),HPPTMP
              WRITE(8,'(A,4E12.4,L5)')          '                WESN FLOWS              ',TVARE(L),TVARE(LE),TVARN(L),TVARN(LN),LMASKDRY(L)
              
              HPPTMP = HP(L)
              UHDYE(L)=0.
              VHDXE(L)=0.
              DO K=1,KC
                LN=LNC(L)
                LE=LEC(L)
                
                UHDYF(L,K)=0.
                VHDXF(L,K)=0.
                UHDY(L,K)=0.
                VHDX(L,K)=0.
                U(L,K)=0.
                V(L,K)=0.
                W(L,K)=0.

                UHDYF(LE,K)=0.
                VHDXF(LN,K)=0.
                UHDY(LE,K)=0.
                VHDX(LN,K)=0.
                U(LE,K)=0.
                V(LN,K)=0.

                UHDY2(L,K)=0.
                VHDX2(L,K)=0.
                U2(L,K)=0.
                V2(L,K)=0.
                W2(L,K)=0.
                UHDY2(LE,K)=0.
                VHDX2(LN,K)=0.
                U2(LE,K)=0.
                V2(LN,K)=0.
              ENDDO
            ELSE
              PRINT '(A,I8,I5,3F10.4,F12.4)',' ERROR! NEG DEPTH IN CONTINUITY CHECK. N,L,H1P,HP,HPNEW,TIMEDAY',N,L,H1P(L),HP(L),HPPTMP,TIMEDAY
              IF( IFILE == -1 )THEN
                IFILE = 8
                OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')
              ENDIF
              WRITE(8,'(A)')              ' ERROR!  NEGATIVE DEPTH IN CONTINUITY CHECK.'
              WRITE(8,'(A,I15,3I5,F15.5)')'         N L ISTL NCTBC TIMEDAY   ',N,L,ISTL_,NCTBC,TIMEDAY
              WRITE(8,'(A,4F10.4)')       '         H1P H2P HP HPNEW         ',H1P(L),H2P(L),HP(L),HPPTMP
              WRITE(8,'(A,4F10.4)')       '         HP  WESN                 ',HP(LW),   HP(LE),   HP(LS),    HP(LN)
              WRITE(8,'(A,4F10.4)')       '         H1P WESN                 ',H1P(LW),  H1P(LE),  H1P(LS),   H1P(LN)
              IF( IS2TIM /= 0 )WRITE(8,'(A,4F10.4)')       '         H2P WESN                 ',H2P(LW),  H2P(LE),  H2P(LS),   H2P(LN)
              WRITE(8,'(A,4F12.5)')       '         SUB/SVB                  ',SUB(L),   SUB(LE),   SVB(L),   SVB(LN)
              WRITE(8,'(A,4E12.4)')       '         UHDYE/VHDXE              ',UHDYE(L), UHDYE(LE), VHDXE(L), VHDXE(LN)
              WRITE(8,'(A,4E12.4)')       '         UHDYE1/VHDXE1            ',UHDY1E(L),UHDY1E(LE),VHDX1E(L),VHDX1E(LN)
              WRITE(8,'(A,4E12.4)')       '         UHDYE2/VHDXE2            ',UHDY2E(L),UHDY2E(LE),VHDX2E(L),VHDX2E(LN)
              WRITE(8,'(A,4E12.4,L5)')    '         WESN FLOWS SUM LAYERS    ',TVARE(L), TVARE(LE), TVARN(L), TVARN(LN), LMASKDRY(L)
              DO K=KC,KSZ(L),-1  
                WRITE(8,'(A,I5,4E12.4)')  '         LAYER FLOWS/CURRENT      ',K,UHDY(L,K),UHDY(LE,K),VHDX(L,K),VHDX(LN,K)
              ENDDO
              DO K=KC,KSZ(L),-1  
                WRITE(8,'(A,I5,4E12.4)')  '         LAYER FLOWS/OLD          ',K,UHDY1(L,K),UHDY1(LE,K),VHDX1(L,K),VHDX1(LN,K)
              ENDDO
              DO K=KC,KSZ(L),-1  
                WRITE(8,'(A,I5,4E12.4)')  '         SUB3D/SVB3D              ',K,SUB3D(L,K),SUB3D(LE,K),SVB3D(L,K),SVB3D(LN,K)
              ENDDO
            ENDIF
          ENDIF

          HP(L)  = HPPTMP    !SPB(L)*HPPTMP+(1.-SPB(L))*(GI*P(L)-BELV(L))  
          HPI(L) = 1./HP(L)  
        ENDDO  
      ENDIF  

    ENDDO 
    !$OMP END DO
  ENDIF

  ! ***************************************************************************
  ! *** Include Nondiverg Wave Stokes Drift In Mass Transport
  IF( ISWAVE > 0 .AND. ISWVSD >= 1 )THEN  
    ! *** WV(L).HEIGHT  - WAVE HEIGHT (M)
    ! *** WV(L).DIR     - WAVE DIRECTION (RADIANS) COUNTER-CLOCKWISE (CELL-EAST AXIS,WAVE)
    ! *** WV(L).FREQ    - WAVE FREQENCY (SEC)
    ! *** WV(L).PERIOD  - WAVE PERIOD (SEC)
    ! *** WV(L).K       - WAVE NUMBER
    ! *** WV(L).LENGTH  - WAVE LENGTH (M)

    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,AMPL,SINHH,C1,C2,STOKESU,STOKESW,THETA,STOKESX,STOKESY) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
    
      IF( ISSSMMT == 2 )THEN
        DO LP=LF,LL  
          L=LWET(LP)  
          IF( LWVMASK(L) )THEN
            C1 = WV(L).KHP
            IF( C1 < 20. .AND. C1 > 1E-8 )THEN
              AMPL = 0.5*WV(L).HEIGHT
              AMPL = AMPL**2
              SINHH = SINH(C1)**2              
              SINHH = AMPL*WV(L).FREQ*WV(L).K/SINHH
              SINH2(L) = 0.50*SINHH
              IF( ABS(UHDYE(L)) > ABS(VHDXE(L)) )THEN
                SINH4(L) = 0.25*SINHH/WV(L).FREQ*WV(L).K * SUB(L)*(HP(LWC(L))-HP(L))*DXIU(L)
              ELSE
                SINH4(L) = 0.25*SINHH/WV(L).FREQ*WV(L).K * SUB(L)*(HP(LSC(L))-HP(L))*DYIV(L)
              ENDIF
            ELSE
              SINH2(L) = 0.
              SINH4(L) = 0.
            ENDIF
        
          ENDIF
        ENDDO
        
        DO K=1,KC  
          DO LP=LF,LL  
            L=LWET(LP)  
            IF( LWVMASK(L) )THEN
              ! *** COMPUTE STOKES DRIFT U/V
              C1 = 2.*WV(L).KHP*ZZ(L,K)
              STOKESU = SINH2(L)*COSH(C1)
              STOKESW = SINH4(L)*SINH(C1)
              THETA = WV(L).DIR 
              STOKESX = SUB(L)*STOKESU*COS(THETA)
              STOKESY = SVB(L)*STOKESU*SIN(THETA)            
              UHDY2(L,K) = UHDY2(L,K) + SUB(L)*STOKESX*DYU(L)*HU(L)
              VHDX2(L,K) = VHDX2(L,K) + SVB(L)*STOKESY*DXV(L)*HV(L)
              U2(L,K)    = U2(L,K) + STOKESX 
              V2(L,K)    = V2(L,K) + STOKESY
              W2(L,K)    = W2(L,K) + STOKESW 
            ENDIF
          ENDDO  
        ENDDO  
      
      ELSEIF( NTSMMT > 1 )THEN
        DO K=1,KC  
          DO LP=LF,LL  
            L=LWET(LP)  
            IF( LWVMASK(L) )THEN
              UHDY2(L,K) = UHDY2(L,K) + UVPT(L,K)*DYU(L)
              VHDX2(L,K) = VHDX2(L,K) + VVPT(L,K)*DXV(L)
              U2(L,K)    = U2(L,K) + UVPT(L,K)*HUI(L)  
              V2(L,K)    = V2(L,K) + VVPT(L,K)*HVI(L)
              W2(L,K)    = W2(L,K) + WVPT(L,K) 
            ENDIF
          ENDDO  
        ENDDO  
      ENDIF
    ENDDO 
    !$OMP END DO
  ENDIF  

  
  !$OMP SINGLE
  
  ! *** RESET OPEN BC CONCENTRATIONS
  DO IOBC=1,NBCSOP  
    L = LOBCS(IOBC)  
    HP(L)  = GI*P(L)-BELV(L)
    HPI(L) = 1./HP(L)  
  ENDDO  

  IF( MDCHH >= 1 )THEN
    RLAMN=QCHERR  
    RLAMO=1.-RLAMN  
    DO NMD=1,MDCHH  
      LHOST=LMDCHH(NMD)  
      LCHNU=LMDCHU(NMD)  
      LCHNV=LMDCHV(NMD)  
      IF( MDCHTYP(NMD) == 1 )THEN  
        TMPVAL=DELT*(RLAMN*QCHANU(NMD)+RLAMO*QCHANUN(NMD))  
        HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)  
        HP(LCHNU)=HP(LCHNU)-TMPVAL*DXYIP(LCHNU)  
        HPI(LHOST)=1./HP(LHOST)  
        HPI(LCHNU)=1./HP(LCHNU)  
      ENDIF  
      IF( MDCHTYP(NMD) == 2 )THEN  
        TMPVAL=DELT*(RLAMN*QCHANV(NMD)+RLAMO*QCHANVN(NMD))  
        HP(LHOST)=HP(LHOST)+TMPVAL*DXYIP(LHOST)  
        HP(LCHNV)=HP(LCHNV)-TMPVAL*DXYIP(LCHNV)  
        HPI(LHOST)=1./HP(LHOST)  
        HPI(LCHNV)=1./HP(LCHNV)  
      ENDIF  
    ENDDO  
  ENDIF  
  !$OMP END SINGLE

  ! *** COMPUTE LAYER THICKNESSES
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L) 
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    DO K=1,KC 
      DO LP=1,LLWET(K,ND)
        L=LKWET(LP,K,ND)
        HPK(L,K) = HP(L)*DZC(L,K)
        HPKI(L,K) = 1./HPK(L,K)
      ENDDO  
    ENDDO

    IF( IGRIDV > 0 )THEN
      DO LP=LF,LL
        L=LWET(LP)  
        ! *** SET DIRECTIONAL DEPTHS
        HPW(L) = HP(L)+BELV(L) - BELVW(L)
        HPE(L) = HP(L)+BELV(L) - BELVE(L)
        HPS(L) = HP(L)+BELV(L) - BELVS(L)
        HPN(L) = HP(L)+BELV(L) - BELVN(L)
      ENDDO
    ENDIF
    
  ENDDO 
  !$OMP END DO

  !$OMP END PARALLEL

  ! ***************************************************************************
  ! *** ACCUMULTATE MAX COURANT NUMBERS  
  IF( ISINWV == 1 .OR. ISNEGH > 0 )THEN
    DO K=1,KC  
      DO L=2,LA  
        CFLUUUT=DELT*ABS(DXIU(L)*U(L,K))  
        CFLUUU(L,K)=MAX(CFLUUUT,CFLUUU(L,K))  
        CFLVVVT=DELT*ABS(DYIV(L)*V(L,K))  
        CFLVVV(L,K)=MAX(CFLVVVT,CFLVVV(L,K))  
        CFLWWWT=DELT*ABS(HPI(L)*DZIG(L,K)*W(L,K))  
        CFLWWW(L,K)=MAX(CFLWWWT,CFLWWW(L,K))  
        CFLCACT=DELT*ABS(CAC(L,K)*DXYIP(L)*HPI(L))  
        CFLCAC(L,K)=MAX(CFLCACT,CFLCAC(L,K))  
      ENDDO  
    ENDDO  
  ENDIF

  ! ***************************************************************************
  ! ***CALCULATE NONHYDROSTATIC PRESSURE  
  IF( KC > 1 .AND. ISPNHYDS >= 1 ) CALL CALPNHS  

  ! ***************************************************************************
  ! *** WRITE TO DIAGNOSTIC FILE CFL.OUT WITH DIAGNOSTICS OF MAXIMUM TIME STEP  
  ! *** SEDIMENT TRANSPORT AND PLACE IN UHDY2, VHDX2 AND W2  
  IF( ISCFL >= 1 .AND. ISTL_ == 3 .AND. DEBUG )THEN  
    OPEN(1,FILE=OUTDIR//'CFL.OUT',STATUS='UNKNOWN',POSITION='APPEND')  
    IF( ISCFLM >= 1 .AND. NITER == 1 )THEN  
      OPEN(2,FILE=OUTDIR//'CFLMP.OUT',STATUS='UNKNOWN')  
      CLOSE(2,STATUS='DELETE')  
      DO L=1,LC  
        ICFLMP(L)=0  
      ENDDO  
    ENDIF  
    DTCFL=1.E+18  
    K=1  
    DO L=2,LA  
      LE=LEC(L)
      LN=LNC(L)  
      UWTMP=ABS(DXIU(L) *U2(L,K))  
      UETMP=ABS(DXIU(LE)*U2(LE,K))  
      VSTMP=ABS(DYIV(L) *V2(L,K))  
      VNTMP=ABS(DYIV(LN)*U2(LN,K))  
      WBTMP=0.  
      WTTMP=ABS(HPKI(L,K)*W2(L,K))  
      DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  +1.0E-12  
      DTMAXX=0.5/DTMAXI  
      IF( DTMAXX < DTCFL )THEN  
        DTCFL=DTMAXX  
        ICFL=IL(L)  
        JCFL=JL(L)  
        KCFL=K  
      ENDIF  
    ENDDO  
    IF( KC > 1 )THEN  
      K=KC  
      DO L=2,LA  
        LN=LNC(L)  
        UWTMP=ABS(DXIU(L  )*U2(L  ,K))  
        UETMP=ABS(DXIU(LE)*U2(LE,K))  
        VSTMP=ABS(DYIV(L  )*V2(L  ,K))  
        VNTMP=ABS(DYIV(LN )*U2(LN ,K))  
        WTTMP=0.  
        WBTMP=ABS(HPKI(L,K)*W2(L,K-1))  
        DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  +1.0E-12  
        DTMAXX=0.5/DTMAXI  
        IF( DTMAXX < DTCFL )THEN  
          DTCFL=DTMAXX  
          ICFL=IL(L)  
          JCFL=JL(L)  
          KCFL=K  
        ENDIF  
      ENDDO  
    ENDIF  
    IF( KC > 2 )THEN  
      DO K=2,KS  
        DO L=2,LA  
          LN=LNC(L)  
          UWTMP=ABS(DXIU(L) *U2(L,K))  
          UETMP=ABS(DXIU(LE)*U2(LE,K))  
          VSTMP=ABS(DYIV(L) *V2(L,K))  
          VNTMP=ABS(DYIV(LN)*U2(LN,K))  
          WBTMP=ABS(HPKI(L,K)*W2(L,K-1))  
          WTTMP=ABS(HPKI(L,K)*W2(L,K))  
          DTMAXI=MAX(UWTMP,UETMP)+MAX(VSTMP,VNTMP)+MAX(WBTMP,WTTMP)  +1.0E-12  
          DTMAXX=0.5/DTMAXI  
          IF( DTMAXX < DTCFL )THEN  
            DTCFL=DTMAXX  
            ICFL=IL(L)  
            JCFL=JL(L)  
            KCFL=K  
          ENDIF  
        ENDDO  
      ENDDO  
    ENDIF  
    IVAL=MOD(N,ISCFL)  
    IDTCFL=NINT(DTCFL)  
    IF( ISCFL == 1 ) WRITE(1,1212)DTCFL,N,ICFL,JCFL,KCFL  
    IF( ISCFL >= 2 .AND. IVAL == 0 )WRITE(1,1213)IDTCFL  
    IF( ISCFLM >= 1 )THEN  
      LTMP=LIJ(ICFL,JCFL)  
      ICFLMP(LTMP)=ICFLMP(LTMP)+1  
    ENDIF  
    IF( ISCFLM >= 1 .AND. N == NTS )THEN  
      OPEN(2,FILE=OUTDIR//'CFLMP.OUT',STATUS='UNKNOWN')  
      TMPVALN=1./FLOAT(NTS)  
      DO L=2,LA  
        TMPVAL=TMPVALN*FLOAT(ICFLMP(L))  
        WRITE(2,1214)IL(L),JL(L),ICFLMP(L),TMPVAL  
      ENDDO  
      CLOSE(2)  
    ENDIF  
    CLOSE(1)  
    1212 FORMAT(' MAX TIME STEP =',F10.2,' SEC FOR N,I,J,K =',I8,3I5)  
    1213 FORMAT(I4)  
    1214 FORMAT(2I5,I12,F10.2)  
  ENDIF  

  IF( IFILE == 8 )CLOSE(8)
  ! ***************************************************************************

  RETURN  

  END  
  
