SUBROUTINE CALEXP2T  

  !**********************************************************************!
  ! **  SUBROUTINE CALEXP2T CALCULATES EXPLICIT MOMENTUM EQUATION TERMS  
  ! **  USING A TWO TIME LEVEL SCHEME  
  
  ! *** VARIABLES   DESCRIPTION                                 UNITS
  ! *** FUHU,FVHU   GROSS MOMENTUM, U COMPONENTS                M4/S2
  ! *** FVHV,FUHV   GROSS MOMENTUM, V COMPONENTS                M4/S2
  ! *** FX, FY      INTERNAL MODE FORCING BY LAYER              M4/S2
  ! *** FBBX, FBBY  INTERNAL MODE BOUYANCY FORCING BY LAYER     M4/S2
  ! *** FCAX,FCAY   CORIOLIS FORCING BY LAYER                   M4/S2
  ! *** DU, DV      INTERNAL SHEARS BY LAYER                    M2/S2
  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! ** 2018-01     PAUL M. CRAIG      CORRECTED WITHDRAWAL/RETURN MOMENTUM OPTION
  ! **                                  ADDED MOMENTUM OPTION FOR STANDARD FLOW BC
  ! ** 2016-02     PAUL M. CRAIG      UPDATED SIGMA-Z (SGZ) FOR EE8.0 
  ! ** 2015-06     PAUL M. CRAIG      IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  ! ** 2014-09     PAUL M. CRAIG      ADDED THE LWET BYPASS APPROACH
  ! ** 2014-01     Paul M. Craig      Fixed the TOT FXVEGE/FYVEGE when partial vegetation penetration using VEGK
  ! ** 2011-03     John Hamrick/Scott James  
  ! **                                Fixed the Vegetative Resistance from AVG to TOT using FKC
  ! ** 2011-03     Paul M. Craig      Converted to F90, added OMP
  ! ** 2010-10     Scott James        Added MHK
  ! ** 2008-12     SANG YUK/PMC       Corrected The Explicit Internal Buoyancy Forcings
  ! ** 11/07/2001  John Hamrick       Added body forces fbodyfx and fbodyfy to external momentum equations
  ! ** 11/14/2001  John Hamrick       Corrected orientation of momentum fluxes from sinks and source
  ! ** 01/02/2002  John Hamrick       Corrected 2 layer (kc=-2) curvature acceleration correction
  ! ** 01/11/2002  John Hamrick       Added ick2cor,ck2uum,ck2vvm,ck2uvm,ck2uuc,ck2vvc,ck2uvc,ck2fcx,
  ! **                                  ck2fcy to generalize two layer momentum flux and curvature 
  ! **                                  acceleration correction
  ! ** 01/15/2002  John Hamrick       Modified calculation of coriolis-curvature accelerations at tidal
  ! **                                  open boundaries
  ! ** 01/23/2002  John Hamrick       Added virtual momentum sources and sinks for subgrid scale channel
  ! **                                  interactions, including local variables tmpvec1,tmpvec2,qmcsinkx,
  ! **                                  qmcsinky,qmcsourx,qmsoury
  ! ** 03/19/2002  John Hamrick       Added dry cell bypass and consistent initialization of dry values

  USE GLOBAL  

  IMPLICIT NONE
  INTEGER :: L, LP, K, LN, LS, ID, JD, KD, NWR, IU, JU, KU, LU, NS, LNW, LSE, LL, IT
  INTEGER :: LD, NMD, LHOST, LCHNU, LW, LE, LCHNV, ND, LF
  INTEGER,SAVE :: NSTB

  REAL :: TMPANG, WUU, WVV, CACSUMT, CFEFF, VEAST2, VWEST2, FCORE, FCORW
  REAL :: UNORT1, USOUT1, UNORT2, USOUT2, FCORN, FCORS, VTMPATU
  REAL :: UTMPATV, UMAGTMP, VMAGTMP, DZPU, DZPV
  REAL :: TMPVAL, WVFACT, DETH, CI11H, CI12H, CI22H, DETU
  REAL :: CI11V, CI12V, CI21V, CI22V, CI21H, CI12U, CI21U, CI22U, DETV, CI11U
  REAL :: UHC, UHB, VHC, VHB, UHC1, UHB1, VHC1, VHB1, UHC2, UHB2, VHC2, VHB2
  REAL :: UHB1MX, UHB1MN, VHC1MX, VHC1MN, UHC1MX, UHC1MN, VHB1MX
  REAL :: VHB1MN, UHB2MX, UHB2MN, VHC2MX, VHC2MN, UHC2MX, UHC2MN, VHB2MX
  REAL :: VHB2MN, BOTT, QMF, QUMF, VEAST1, VWEST1, QWRABS, VDIR

  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: CACSUM  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DZPC
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TMPVEC1  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: TMPVEC2  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: FUHJ  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: FVHJ  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QMCSINKX  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QMCSINKY  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QMCSOURX  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QMCSOURY  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: BK
  
  IF(  .NOT. ALLOCATED(TMPVEC1) )THEN
    ALLOCATE(CACSUM(NTHREADS))  
    ALLOCATE(FUHJ(LCM,KCM))  
    ALLOCATE(FVHJ(LCM,KCM))  
    ALLOCATE(TMPVEC1(KCM))  
    ALLOCATE(TMPVEC2(KCM))
    ALLOCATE(BK(KCM,NDM))
    FUHJ=0.
    FVHJ=0.
    TMPVEC1=0.
    TMPVEC2=0.
    NSTB=0
    CACSUM=0.
    BK=0.
    IF( MDCHH >= 1 )THEN
      ALLOCATE(QMCSINKX(LCM,KCM))  
      ALLOCATE(QMCSINKY(LCM,KCM))  
      ALLOCATE(QMCSOURX(LCM,KCM))  
      ALLOCATE(QMCSOURY(LCM,KCM))  
      QMCSINKX=0.
      QMCSINKY=0.
      QMCSOURX=0.
      QMCSOURY=0.
    ENDIF
    IF( ISPNHYDS > 0 )THEN
      ALLOCATE(DZPC(LCM,KCM))
      DZPC=0.
    ENDIF
  ENDIF

  IF( ISDYNSTP == 0 )THEN  
    DELT=DT  
  ELSE  
    DELT=DTDYN  
  ENDIF  

  IF( IS2TIM == 2 )THEN  
    DELT=0.5*DT  
  ENDIF  
 
  DELTI=1./DELT  

  IF( N == 1 .AND. DEBUG )THEN  
    OPEN(1,FILE=OUTDIR//'MFLUX.DIA')  
    CLOSE(1,STATUS='DELETE')  
  ENDIF  

  ! *** WAVE RAMPUP FACTOR
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN
    IF( N < NTSWV )THEN  
      TMPVAL = FLOAT(N)/FLOAT(NTSWV)  
      WVFACT = 0.5-0.5*COS(PI*TMPVAL)  
    ELSE  
      WVFACT = 1.0  
    ENDIF  
  ENDIF
  CACSUMT = 0.

  !**********************************************************************!  
  ! *** ZERO NEWLY DRY CELL SHEARS/MOMENTUM
  IF( LADRY > 0 )THEN
    DO LP=1,LADRY
      L=LDRY(LP)  
      FCAXE(L)=0.  
      FCAYE(L)=0.  
      FXE(L)=0.  
      FYE(L)=0.  
    ENDDO

    DO K=1,KC
      DO LP=1,LADRY
        L=LDRY(LP)  
        FUHU(L,K)=0.  
        FVHU(L,K)=0.  
        FUHV(L,K)=0.  
        FVHV(L,K)=0.  
        FWU(L,K)=0.
        FWV(L,K)=0.
        CAC(L,K)=0.0
        FCAX(L,K)=0.
        FCAY(L,K)=0.
        FX(L,K)=0.
        FY(L,K)=0.
        FBBX(L,K)=0.
        FBBY(L,K)=0.
        DU(L,K)=0.0  
        DV(L,K)=0.0  
        
        ! *** TWO LAYER ROTATIONAL EFFECTS OR WITHDRAWAL/RETURN
        FUHJ(L,K)=0.  
        FVHJ(L,K)=0.  
      ENDDO
    ENDDO
    
    IF( ISVEG > 0 )THEN
      DO LP=1,LADRY
        L=LDRY(LP)  
        FXVEGE(L)=0.
        FYVEGE(L)=0.
      ENDDO
      DO K=1,KC
        DO LP=1,LADRY
          L=LDRY(LP)  
          FXVEG(L,K)=0.0
          FYVEG(L,K)=0.0
        ENDDO
      ENDDO
    ENDIF
  ENDIF

  !**********************************************************************!  
  !$OMP PARALLEL DEFAULT(SHARED)

  !**********************************************************************!  
  !  
  ! **  INITIALIZE EXTERNAL CORIOLIS-CURVATURE AND ADVECTIVE FLUX TERMS  
  !  
  !----------------------------------------------------------------------!  
  !$OMP DO PRIVATE(ND,LF,LL,LP,L)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    DO LP=LF,LL
      L=LWET(LP)  
      FCAXE(L)=0.  
      FCAYE(L)=0.  
      FXE(L)=0.  
      FYE(L)=0.  
    ENDDO  
  ENDDO
  !$OMP END DO

  !----------------------------------------------------------------------!  
  IF( IS2LMC /= 1 )THEN 
    ! *** STANDARD CALCULATION. WITHOUT MOMENTUM-CURVATURE CORRECTION  
    !$OMP DO PRIVATE(ND,K,L,LF,LL,LP,LE,LN,LS,LW,UHC,UHB,VHC,VHB)                
    DO ND=1,NDM 
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)  
          LS = LSC(L)
          LW = LWC(L)

          ! *** U COMPONENTS  
          UHB = 0.5*(UHDY(L,K)+UHDY(LE,K))
          VHC = 0.5*(VHDX(L,K)+VHDX(LW,K)) 
          
          ! ***      |-- EAST FLOWING --|  |-- WEST FLOWING --|
          FUHU(L,K) = MAX(UHB,0.)*U(L,K)  + MIN(UHB,0.)*U(LE,K)
          ! ***      |-- NORTH FLOWING -|  |-- SOUTH FLOWING -|
          FVHU(L,K) = MAX(VHC,0.)*U(LS,K) + MIN(VHC,0.)*U(L,K) 
          
          ! *** V COMPONENTS
          VHB = 0.5*(VHDX(L,K)+VHDX(LN,K))
          UHC = 0.5*(UHDY(L,K)+UHDY(LS,K))
          
          ! ***      |-- NORTH FLOWING -|  |-- SOUTH FLOWING -|
          FVHV(L,K) = MAX(VHB,0.)*V(L,K)  + MIN(VHB,0.)*V(LN,K)
          ! ***      |-- EAST FLOWING --|  |-- WEST FLOWING --|
          FUHV(L,K) = MAX(UHC,0.)*V(LW,K) + MIN(UHC,0.)*V(L,K) 

        ENDDO
      ENDDO  
    ENDDO
    !$OMP END DO

  ELSEIF( IS2LMC == 1 .AND. KC == 2 )THEN  
    ! *** CALCULATION FOR MOMENTUM-CURVATURE CORRECTION (TWO LAYER ONLY)
    !$OMP DO PRIVATE(ND,K,L,LP,LF,LL,LE,LN,LS,LW,UHC1,UHB1,VHC1,VHB1) &
    !$OMP    PRIVATE(UHC2,UHB2,VHC2,VHB2,UHB1MX,UHB1MN,VHC1MX)  &
    !$OMP    PRIVATE(VHC1MN,UHC1MX,UHC1MN,VHB1MX,VHB1MN,UHB2MX,UHB2MN)  &
    !$OMP    PRIVATE(VHC2MX,VHC2MN,UHC2MX,UHC2MN,VHB2MX,VHB2MN,BOTT) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L = LWET(LP)
        LE = LEC(L)
        LN = LNC(L)  
        LS = LSC(L)  
        LW = LWC(L)
        UHC1 = 0.5*(UHDYF(L,1)+UHDYF(LS,1))  
        UHB1 = 0.5*(UHDYF(L,1)+UHDYF(LE,1))  
        VHC1 = 0.5*(VHDXF(L,1)+VHDXF(LW,1))  
        VHB1 = 0.5*(VHDXF(L,1)+VHDXF(LN,1))  
        UHC2 = 0.5*(UHDYF(L,2)+UHDYF(LS,2))  
        UHB2 = 0.5*(UHDYF(L,2)+UHDYF(LE,2))  
        VHC2 = 0.5*(VHDXF(L,2)+VHDXF(LW,2))  
        VHB2 = 0.5*(VHDXF(L,2)+VHDXF(LN,2))  
   
        UHB1MX = 0.  
        UHB1MN = 0.  
        VHC1MX = 0.  
        VHC1MN = 0.  
        UHC1MX = 0.  
        UHC1MN = 0.  
        VHB1MX = 0.  
        VHB1MN = 0.  
        UHB2MX = 0.  
        UHB2MN = 0.  
        VHC2MX = 0.  
        VHC2MN = 0.  
        UHC2MX = 0.  
        UHC2MN = 0.  
        VHB2MX = 0.  
        VHB2MN = 0.  
   
        BOTT = ABS(UHB1*U(L,1))  
        IF( BOTT > 0.0)UHB1MX = 1.+CK2UUM*(UHB2-UHB1)*(U(L,2)  -U(L,1))  /UHB1*U(L,1)  
        BOTT = ABS(UHB1*U(LE,1))  
        IF( BOTT > 0.0)UHB1MN = 1.+CK2UUM*(UHB2-UHB1)*(U(LE,2)-U(LE,1))/UHB1*U(LE,1)  
        BOTT = ABS(VHC1*U(LS,1))  
        IF( BOTT > 0.0)VHC1MX = 1.+CK2UVM*(VHC2-VHC1)*(U(LS,2) -U(LS,1)) /VHC1* U(LS,1)  
        BOTT = ABS(VHC1*U(L,1))  
        IF( BOTT > 0.0)VHC1MN = 1.+CK2UVM*(VHC2-VHC1)*(U(L,2)  -U(L,1))  /VHC1*U(L,1)  
        BOTT = ABS(UHC1*V(LW,1))  
        IF( BOTT > 0.0)UHC1MX = 1.+CK2UVM*(UHC2-UHC1)*(V(LW,2)-V(LW,1))/UHC1*V(LW,1)  
        BOTT = ABS(UHC1*V(L,1))  
        IF( BOTT > 0.0)UHC1MN = 1.+CK2UVM*(UHC2-UHC1)*(V(L,2)  -V(L,1))  /UHC1*V(L,1)  
        BOTT = ABS(VHB1*V(L,1))  
        IF( BOTT > 0.0)VHB1MX = 1.+CK2VVM*(VHB2-VHB1)*(V(L,2)  -V(L,1))  /VHB1*V(L,1)  
        BOTT = ABS(VHB1*V(LN,1))  
        IF( BOTT > 0.0)VHB1MN = 1.+CK2VVM*(VHB2-VHB1)*(V(LN,2) -V(LN,1)) /VHB1* V(LN,1)  
        BOTT = ABS(UHB2*U(L,2))  
        IF( BOTT > 0.0)UHB2MX = 1.+CK2UUM*(UHB2-UHB1)*(U(L,2)  -U(L,1))  /UHB2*U(L,2)  
        BOTT = ABS(UHB2*U(LE,2))  
        IF( BOTT > 0.0)UHB2MN = 1.+CK2UUM*(UHB2-UHB1)*(U(LE,2)-U(LE,1))/UHB2*U(LE,2)  
        BOTT = ABS(VHC2*U(LS,2))  
        IF( BOTT > 0.0)VHC2MX = 1.+CK2UVM*(VHC2-VHC1)*(U(LS,2) -U(LS,1)) /VHC2* U(LS,2)  
        BOTT = ABS(VHC2*U(L,2))  
        IF( BOTT > 0.0)VHC2MN = 1.+CK2UVM*(VHC2-VHC1)*(U(L,2)  -U(L,1))  /VHC2*U(L,2)  
        BOTT = ABS(UHC2*V(LW,2))  
        IF( BOTT > 0.0)UHC2MX = 1.+CK2UVM*(UHC2-UHC1)*(V(LW,2)-V(LW,1))/UHC2*V(LW,2)  
        BOTT = ABS(UHC2*V(L,2))  
        IF( BOTT > 0.0)UHC2MN = 1.+CK2UVM*(UHC2-UHC1)*(V(L,2)  -V(L,1))  /UHC2*V(L,2)  
        BOTT = ABS(VHB2*V(L,2))  
        IF( BOTT > 0.0)VHB2MX = 1.+CK2VVM*(VHB2-VHB1)*(V(L,2)  -V(L,1))  /VHB2*V(L,2)  
        BOTT = ABS(VHB2*V(LN,2))  
        IF( BOTT > 0.0)VHB2MN = 1.+CK2VVM*(VHB2-VHB1)*(V(LN,2) -V(LN,1)) /VHB2* V(LN,2)  
   
        FUHU(L,1) = UHB1MX*MAX(UHB1,0.)*U(L,1)  +UHB1MN*MIN(UHB1,0.)*U(LE,1)  
        FVHU(L,1) = VHC1MX*MAX(VHC1,0.)*U(LS,1) +VHC1MN*MIN(VHC1,0.)*U(L,1)  
        FUHV(L,1) = UHC1MX*MAX(UHC1,0.)*V(LW,1)+UHC1MN*MIN(UHC1,0.)*V(L,1)  
        FVHV(L,1) = VHB1MX*MAX(VHB1,0.)*V(L,1)  +VHB1MN*MIN(VHB1,0.)*V(LN,1)  
        FUHJ(L,1) = 0.  
        FVHJ(L,1) = 0.  
        FUHU(L,2) = UHB2MX*MAX(UHB2,0.)*U(L,2)  +UHB2MN*MIN(UHB2,0.)*U(LE,2)  
        FVHU(L,2) = VHC2MX*MAX(VHC2,0.)*U(LS,2) +VHC2MN*MIN(VHC2,0.)*U(L,2)  
        FUHV(L,2) = UHC2MX*MAX(UHC2,0.)*V(LW,2)+UHC2MN*MIN(UHC2,0.)*V(L,2)  
        FVHV(L,2) = VHB2MX*MAX(VHB2,0.)*V(L,2)  +VHB2MN*MIN(VHB2,0.)*V(LN,2)  
        FUHJ(L,2) = 0.  
        FVHJ(L,2) = 0.  
      ENDDO
    ENDDO  ! *** END OF DOMAIN
    !$OMP END DO

  ENDIF
  
  !$OMP SINGLE
  ! *** ADD WITHDRAWAL/RETURN FLOW MOMENTUM FLUXES
  DO NWR=1,NQWR  
    IF( ABS(NQWRMFU(NWR)) > 0 )THEN  
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NS = NQWRSERQ(NWR)  
      IF( QWRSERT(NS) >= 0. )THEN
        ! *** Original Withdrawal/Return
        IU=IQWRU(NWR)  
        JU=JQWRU(NWR)  
        KU=KQWRU(NWR)
        VDIR = 1.
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        IU=IQWRD(NWR)  
        JU=JQWRD(NWR)  
        KU=KQWRD(NWR) 
        VDIR = -1.
      ENDIF
      LU=LIJ(IU,JU)  

      QWRABS = ABS(QWRSERT(NS))
      QMF = QWR(NWR) + QWRABS 
      QUMF = VDIR*QMF*( QMF/( HPK(LU,KU)*BQWRMFU(NWR) ) )   ! *** M4/S2
      IF( NQWRMFU(NWR) ==  1 ) FUHJ(LU     ,KU) = -QUMF  
      IF( NQWRMFU(NWR) ==  2 ) FVHJ(LU     ,KU) = -QUMF  
      IF( NQWRMFU(NWR) ==  3 ) FUHJ(LEC(LU),KU) = -QUMF  
      IF( NQWRMFU(NWR) ==  4 ) FVHJ(LNC(LU),KU) = -QUMF  
      IF( NQWRMFU(NWR) == -1 ) FUHJ(LU     ,KU) = QUMF  
      IF( NQWRMFU(NWR) == -2 ) FVHJ(LU     ,KU) = QUMF  
      IF( NQWRMFU(NWR) == -3 ) FUHJ(LEC(LU),KU) = QUMF  
      IF( NQWRMFU(NWR) == -4 ) FVHJ(LNC(LU),KU) = QUMF  
    ENDIF  
    IF( ABS(NQWRMFD(NWR)) > 0 )THEN  
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      NS = NQWRSERQ(NWR)  
      IF( QWRSERT(NS) >= 0. )THEN
        ! *** Original Withdrawal/Return
        ID=IQWRD(NWR)  
        JD=JQWRD(NWR)  
        KD=KQWRD(NWR)  
        VDIR = 1.
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        ID=IQWRU(NWR)  
        JD=JQWRU(NWR)  
        KD=KQWRU(NWR) 
        VDIR = -1.
      ENDIF
      LD=LIJ(ID,JD)

      QWRABS = ABS(QWRSERT(NS))
      QMF = QWR(NWR) + QWRABS 
      QUMF = VDIR*QMF*( QMF/( HPK(LD,KD)*BQWRMFD(NWR) ) )   ! *** M4/S2
      IF( NQWRMFD(NWR) ==  1 ) FUHJ(LD     ,KD) = -QUMF  
      IF( NQWRMFD(NWR) ==  2 ) FVHJ(LD     ,KD) = -QUMF  
      IF( NQWRMFD(NWR) ==  3 ) FUHJ(LEC(LD),KD) = -QUMF  
      IF( NQWRMFD(NWR) ==  4 ) FVHJ(LNC(LD),KD) = -QUMF  
      IF( NQWRMFD(NWR) == -1 ) FUHJ(LD     ,KD) = QUMF  
      IF( NQWRMFD(NWR) == -2 ) FVHJ(LD     ,KD) = QUMF  
      IF( NQWRMFD(NWR) == -3 ) FUHJ(LEC(LD),KD) = QUMF  
      IF( NQWRMFD(NWR) == -4 ) FVHJ(LNC(LD),KD) = QUMF  
    ENDIF  
  ENDDO  
  
  ! *** ADD QSER MOMENTUM FLUXES
  DO LL=1,NQSIJ
    IF( ABS(NQSMF(LL)) > 0 )THEN  
      L=LQS(LL)
      DO K=KSZ(L),KC
        ! *** Handle reversing flows in/out of domain
        IF( QSERCELL(K,LL) >= 0. )THEN
          VDIR = 1.
        ELSE
          VDIR = -1.
        ENDIF
        !LD = LIJ(IQS(LL),JQS(LL))
        
        QMF = ABS(QSERCELL(K,LL))
        QUMF = VDIR*QMF*( QMF/( HPK(L,K)*QWIDTH(LL) ) )   ! *** M4/S2
        IF( NQSMF(LL) ==  1 ) FUHJ(L     ,K) = -QUMF  
        IF( NQSMF(LL) ==  2 ) FVHJ(L     ,K) = -QUMF  
        IF( NQSMF(LL) ==  3 ) FUHJ(LEC(L),K) = -QUMF  
        IF( NQSMF(LL) ==  4 ) FVHJ(LNC(L),K) = -QUMF  
        IF( NQSMF(LL) == -1 ) FUHJ(L     ,K) = QUMF  
        IF( NQSMF(LL) == -2 ) FVHJ(L     ,K) = QUMF  
        IF( NQSMF(LL) == -3 ) FUHJ(LEC(L),K) = QUMF  
        IF( NQSMF(LL) == -4 ) FVHJ(LNC(L),K) = QUMF  
      ENDDO
    ENDIF  
  ENDDO
  !$OMP END SINGLE

  !----------------------------------------------------------------------!  
  !  
  ! *** COMPUTE VERTICAL ACCELERATIONS
  IF( KC > 1 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L,LW,LS,WUU,WVV)
    DO ND=1,NDM  
      DO K=1,KS  
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LW = LWC(L)
          LS = LSC(L)
          WUU = 0.5*DXYU(L)*(W(L,K)+W(LW,K))  
          WVV = 0.5*DXYV(L)*(W(L,K)+W(LS,K))  
          FWU(L,K) = MAX(WUU,0.)*U(L,K) + MIN(WUU,0.)*U(L,K+1)
          FWV(L,K) = MAX(WVV,0.)*V(L,K) + MIN(WVV,0.)*V(L,K+1)
        ENDDO  
      ENDDO
    ENDDO  
    !$OMP END DO
  
    ! *** APPLY OPEN BOUNDARYS
    !$OMP SINGLE
    DO LL=1,NBCSOP2
      L=LOBCS2(LL)
      DO K=1,KS  
        FWU(L,K) = 0.0
        FWV(L,K) = 0.0
      ENDDO  
    ENDDO 
    !$OMP END SINGLE 
  ENDIF
  
  ! **********************************************************************!  
  ! ** BLOCK MOMENTUM FLUX ON LAND SIDE OF TRIANGULAR CELLS  
  IF( ITRICELL > 0 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          FUHU(L,K) = STCUV(L)*FUHU(L,K)
          FVHV(L,K) = STCUV(L)*FVHV(L,K)
        ENDDO
      ENDDO
    ENDDO   ! ***  END OF DOMAIN
    !$OMP END DO
  ENDIF
  
  !**********************************************************************!  
  ! *** CALCULATE CORIOLIS AND CURVATURE ACCELERATION COEFFICIENTS  
  !$OMP SINGLE
  CACSUM = 0. 
  CFMAX = CF  
  !$OMP END SINGLE
  
  IF( ISCURVATURE )THEN

    IF( ISDCCA == 0 )THEN  
      ! *** STANDARD CALCULATIONS, NO DIAGNOSTICS
      !$OMP DO PRIVATE(ND,K,LP,L,LE,LN)
      DO ND=1,NDM  

        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LE=LEC(L)
            LN=LNC(L)
            CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)
            CACSUM(ND) = CACSUM(ND)+CAC(L,K)
          ENDDO
        ENDDO
      ENDDO   ! ***  END OF DOMAIN
      !$OMP END DO
      
    ELSE  
      ! *** STANDARD CALCULATIONS, WITH DIAGNOSTICS
      !$OMP SINGLE
      IT = 1
      DO K=1,KC  
        DO L=2,LA
          LE = LEC(L)
          LN = LNC(L)  
          CAC(L,K) = ( FCORC(L)*DXYP(L) + 0.5*SNLT*(V(LN,K)+V(L,K))*DYDI(L) - 0.5*SNLT*(U(LE,K)+U(L,K))*DXDJ(L) )*HP(L)  
          CFEFF = ABS(CAC(L,K))*DXYIP(L)*HPI(L)  
          CFMAX = MAX(CFMAX,CFEFF)  
          CACSUM(IT) = CACSUM(IT)+CAC(L,K) 
        ENDDO  
      ENDDO  
      !$OMP END SINGLE
    ENDIF  

    ! *** ENSURE FCAY & FCAX ARE RESET
    !$OMP SINGLE
    CACSUMT = ABS(SUM(CACSUM(:)))
    !$OMP END SINGLE
    
    IF( CACSUMT < 1.E-7 )THEN
      !$OMP DO PRIVATE(ND,K,LP,L)
      DO ND=1,NDM  
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            FCAX(L,K) = 0.
            FCAY(L,K) = 0.
          ENDDO
        ENDDO
      ENDDO  
      !$OMP END DO
    ENDIF
        
  ENDIF

  !**********************************************************************!  
  !  
  ! **  CALCULATE CORIOLIS-CURVATURE AND ADVECTIVE ACCELERATIONS  
  !  
  !----------------------------------------------------------------------!  

  ! **  STANDARD CALCULATION
  IF( CACSUMT > 1.E-7 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW,LNW,LSE)
    DO ND=1,NDM
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LE = LEC(L)
          LN = LNC(L)  
          LS = LSC(L)  
          LW = LWC(L)
          LNW = LNWC(L)  
          LSE = LSEC(L)  
          FCAX(L,K) = 0.25*SCAX(L)*( CAC(L,K)*(V(LN,K)+V(L,K)) + CAC(LW,K)*(V(LNW,K)+V(LW,K)) )  
          FCAY(L,K) = 0.25*SCAY(L)*( CAC(L,K)*(U(LE,K)+U(L,K)) + CAC(LS,K)*(U(LSE,K)+U(LS,K)) )  
        ENDDO
      ENDDO  
    ENDDO  
    !$OMP END DO

  ENDIF

  !----------------------------------------------------------------------!  
  !  
  ! *** CALCULATION FOR MOMENTUM-CURVATURE CORRECTION  
  IF( IS2LMC == 1 .AND. CACSUMT > 1.E-7 )THEN  
    !$OMP SINGLE
    DO LP=1,LAWET
      L = LWET(LP)
      LE = LEC(L)
      LN = LNC(L)  
      LS = LSC(L)  
      LW = LWC(L)
      LNW = LNWC(L)  
      LSE = LSEC(L)  

      VEAST1 = V(LN,1)+V(L,1)  
      VWEST1 = V(LNW,1)+V(LW,1)  
      VEAST2 = V(LN,2)+V(L,2)  
      VWEST2 = V(LNW,2)+V(LW,2)  
        
      FCORE = CK2FCX*(CAC(L,2) -CAC(L,1)) *(VEAST2-VEAST1)  
      FCORW = CK2FCX*(CAC(LW,2)-CAC(LW,1))*(VWEST2-VWEST1)
        
      FCAX(L,1) = 0.25*SCAX(L)*(CAC(L,1)*VEAST1+FCORE +CAC(LW,1)     *VWEST1+FCORW)  
      FCAX(L,2) = 0.25*SCAX(L)*(CAC(L,2)*VEAST2+FCORE +CAC(LWC(LW),2)*VWEST2+FCORW)

      UNORT1 = U(LE,1)+U(L,1)  
      USOUT1 = U(LSE,1)+U(LS,1)  
      UNORT2 = U(LE,2)+U(L,2)  
      USOUT2 = U(LSE,2)+U(LS,2)  
      FCORN = CK2FCY*(CAC(L ,2)-CAC(L ,1))*(UNORT2-UNORT1)  
      FCORS = CK2FCY*(CAC(LS,2)-CAC(LS,1))*(USOUT2-USOUT1)  

      FCAY(L,1) = 0.25*SCAY(L)*(CAC(L,1)*UNORT1+FCORN +CAC(LS,1)*USOUT1+FCORS)  
      FCAY(L,2) = 0.25*SCAY(L)*(CAC(L,2)*UNORT2+FCORN +CAC(LS,2)*USOUT2+FCORS)  
    ENDDO  
    !$OMP END SINGLE
  ENDIF  
  
  !----------------------------------------------------------------------!  
  ! **  MODIFICATION FOR TYPE 2 OPEN BOUNDARIES  
  !$OMP SINGLE
  DO LL=1,NPBW  
    IF( ISPBW(LL) == 2 )THEN  
      L = LPBW(LL)+1  
      LN = LNC(L)  
      DO K=KSZ(L),KC  
        FCAX(L,K) = 0.5*SCAX(L)*CAC(L,K)*(V(LN,K)+V(L,K))  
      ENDDO  
    ENDIF  
  ENDDO  

  DO LL=1,NPBE  
    IF( ISPBE(LL) == 2 )THEN  
      L=LPBE(LL)  
      LNW=LNWC(L)  
      DO K=KSZ(L),KC  
        FCAX(L,K) = 0.5*SCAX(L)*CAC(LWC(L),K)*(V(LNW,K)+V(LWC(L),K))  
      ENDDO  
    ENDIF  
  ENDDO  

  DO LL=1,NPBS  
    IF( ISPBS(LL) == 2 )THEN  
      L=LNC(LPBS(LL))  
      DO K=KSZ(L),KC  
        FCAY(L,K) = 0.5*SCAY(L)*CAC(L,K)*(U(LEC(L),K)+U(L,K))  
      ENDDO  
    ENDIF  
  ENDDO  

  DO LL=1,NPBN  
    IF( ISPBN(LL) == 2 )THEN  
      L=LPBN(LL)  
      LS=LSC(L)  
      LSE=LSEC(L)  
      DO K=KSZ(L),KC  
        FCAY(L,K) = 0.5*SCAY(L)*CAC(LS,K)*(U(LSE,K)+U(LS,K))  
      ENDDO  
    ENDIF  
  ENDDO  
  !$OMP END SINGLE

  !----------------------------------------------------------------------!  
  !  INITIALIZE EXTERNAL FORCINGS WITH GROSS MOMENTUM
  !$OMP DO PRIVATE(ND,K,LP,L,LE,LN,LS,LW)
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        LE = LEC(L)
        LN = LNC(L) 
        LS = LSC(L) 
        LW = LWC(L)
        FX(L,K) = FSGZU(L,K)*( FUHU(L,K)-FUHU(LW,K) + FVHU(LN,K)-FVHU(L,K) ) + FUHJ(L,K)   ! ***  M4/S2
        FY(L,K) = FSGZV(L,K)*( FVHV(L,K)-FVHV(LS,K) + FUHV(LE,K)-FUHV(L,K) ) + FVHJ(L,K)   ! ***  M4/S2
      ENDDO  
    ENDDO
  ENDDO
  !$OMP END DO

  ! *** TREAT BC'S NEAR EDGES
  !$OMP SINGLE
  DO LL=1,NBCS
    ! *** BC CELL
    L = LBCS(LL)
    DO K=KSZ(L),KC
      FX(L,K) = SAAX(L)*FX(L,K) + FUHJ(L,K)   ! ***  M4/S2
      FY(L,K) = SAAY(L)*FY(L,K) + FVHJ(L,K)   ! ***  M4/S2
    ENDDO

    ! *** EAST/WEST ADJACENT CELL
    L = MAX(1,LBERC(LL))
    DO K=KSZ(L),KC
      FX(L,K) = SAAX(L)*FX(L,K) + FUHJ(L,K)   ! ***  M4/S2
    ENDDO

    ! *** NORTH/SOUTH ADJACENT CELL
    L = MAX(1,LBNRC(LL))
    DO K=KSZ(L),KC
      FY(L,K) = SAAY(L)*FY(L,K) + FVHJ(L,K)   ! ***  M4/S2
    ENDDO
    
  ENDDO  
  !$OMP END SINGLE
      
  !**********************************************************************!  
  !  
  ! ***  ADD VEGETATION DRAG TO HORIZONTAL ADVECTIVE ACCELERATIONS  
  !  
  !----------------------------------------------------------------------!  
  IF( ISVEG > 0 )THEN  

    ! *** ADD IN MHK DEVICES, IF NEEDED
    !$OMP SINGLE
    IF( LMHK ) CALL MHKPWRDIS
    !$OMP END SINGLE

    IF( ISDRY > 0 .AND. LADRY > 0 )THEN
      ! *** UPDATE ACTIVE CELL BY LAYER LIST
      !$OMP DO PRIVATE(ND,K,LN,LP,L)
      DO ND=1,NDM  
        DO K=1,KC  
          LN=0
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            IF( LVEG(L) )THEN
              LN = LN+1
              LKVEG(LN,K,ND) = L
            ENDIF
          ENDDO
          LLVEG(K,ND) = LN
        ENDDO
      ENDDO
      !$OMP END DO
    ENDIF
    
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LW,LE,LS,LN,LNW,LSE) &
    !$OMP    PRIVATE(UTMPATV,VTMPATU,UMAGTMP,VMAGTMP)
    DO ND=1,NDM  

      DO LP=1,LLVEG(KC,ND)
        L=LKVEG(LP,KC,ND) 
        FXVEGE(L) = 0.  
        FYVEGE(L) = 0.  
      ENDDO  

      DO K=1,KC  
        DO LP=1,LLVEG(K,ND)
          L=LKVEG(LP,K,ND)
          LW = LWC(L)    !west cell
          LE = LEC(L)    !east cell
          LS = LSC(L)    !south cell
          LN = LNC(L)    !north cell
          LNW = LNWC(L)  !northwest cell
          LSE = LSEC(L)  !southeast cell
          UTMPATV = 0.25*(U(L,K)+U(LE,K)+U(LS,K)+U(LSE,K))       !u-velocity at v face
          VTMPATU = 0.25*(V(L,K)+V(LW,K)+V(LN,K)+V(LNW,K))       !v-velocity at u face
          UMAGTMP = SQRT( U(L,K)*U(L,K) +   VTMPATU*VTMPATU )    !u-face velocity vector magnitude
          VMAGTMP = SQRT( UTMPATV*UTMPATV + V(L,K)*V(L,K) )      !v-face velocity vector magnitude

          !FXVEG/FYVEG come from CALTBXY unitless, but they are really just a form of the drag coefficient with terms accounting for the area density
          !FXVEG/FYVEG only change inasmuch as the water depth changes and are zero in layers not penetrated by vegetation
          !FXVEG/FYVEG are C_d(N/L^2) 
          !FXVEG/FYVEG are now multiplied by the cell area and cell-averaged velocity
          !FXVEG/FYVEG are C_d(N/L^2)A|q|
          FXVEG(L,K) = UMAGTMP*SUB3D(L,K)*FXVEG(L,K)   ![m/s] q_xC_d
          FYVEG(L,K) = VMAGTMP*SVB3D(L,K)*FYVEG(L,K)   ![m/s] q_yC_d
            
          !FXVEG/FXVEGE are multiplied by the local velocity to yield units of [m^4/s^2]
          !FXVEG/FXVEGE are added to the body forces as C_d(N/L^2)A|q|q
          FXVEGE(L) = FXVEGE(L)+FXVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FUHDXE 
          FYVEGE(L) = FYVEGE(L)+FYVEG(L,K)*DZC(L,K)    !Integrating the vegetative resistance [m/s] used in FVHDYE
        ENDDO  
      ENDDO  

      ! *** ADD VEGETATIVE DRAG TO INTERNAL MODE SHEAR AND EXTERNAL MODE MOMENTUM
      DO K=1,KC  
        DO LP=1,LLVEG(K,ND)
          L=LKVEG(LP,K,ND)
          IF( (K-KSZ(L)+1) > INT(VEGK(L)+1.) )CYCLE
          FX(L,K) = FX(L,K) + (FXVEG(L,K)-FXVEGE(L))*U(L,K)*DXYU(L) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) FXVEGE goes into FUHDXE for momentum conservation
          FY(L,K) = FY(L,K) + (FYVEG(L,K)-FYVEGE(L))*V(L,K)*DXYV(L) ![m^4/s^2] adding vegetative resistance to the body force (no net force added) FYVEGE goes into FVHDYE for momentum conservation
        ENDDO  
      ENDDO  

      ! *** CONVERT THE AVG FXVEGE/FYVEGE TO TOTAL FXVEGE/FYVEGE
      DO LP=1,LLVEG(KC,ND)
        L=LKVEG(LP,KC,ND) 
        FXVEGE(L) = FXVEGE(L)*HUI(L)*VEGK(L) !Calculate vegetative dissipation for FUHDYE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
        FYVEGE(L) = FYVEGE(L)*HVI(L)*VEGK(L) !Calculate vegetative dissipation for FVHDXE for momentum conservation in CALPUV (need to have sum of forces, not average provided to CALPUV)
      ENDDO 
      
      IF( LMHK )THEN 
        DO LP=1,LLVEG(KC,ND)
          L=LKVEG(LP,KC,ND) 
          FXVEGE(L) = FXVEGE(L) + FXMHKE(L)   ! Add MHK to vegetative dissipation in FUHDYE for momentum conservation in CALPUV multiply by HUI took place in MHKPWRDIS
          FYVEGE(L) = FYVEGE(L) + FYMHKE(L)   ! Add MHK to vegetative dissipation in FVHDXE for momentum conservation in CALPUV multiply by HVI took place in MHKPWRDIS
          FXVEGE(L) = FXVEGE(L) + FXSUPE(L)   ! Add MHK support to vegetative dissipation in FUHDYE for momentum conservation in CALPUV
          FYVEGE(L) = FYVEGE(L) + FYSUPE(L)   ! Add MHK support to vegetative dissipation in FVHDXE for momentum conservation in CALPUV
        ENDDO
      ENDIF
    ENDDO  ! *** END OF DOMAIN
    !$OMP END DO

  ENDIF  

  !**********************************************************************!  
  !  
  ! **  ADD HORIZONTAL MOMENTUM DIFFUSION TO ADVECTIVE ACCELERATIONS  
  !  
  !----------------------------------------------------------------------!  
  IF( ISHDMF >= 1 )THEN

    !$OMP DO PRIVATE(ND,K,LP,L) 
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          FX(L,K) = FX(L,K) - SDX(L)*( FMDUX(L,K) + FMDUY(L,K) )
          FY(L,K) = FY(L,K) - SDY(L)*( FMDVX(L,K) + FMDVY(L,K) )  
        ENDDO  
      ENDDO
    ENDDO
    !$OMP END DO

  ENDIF  

  !**********************************************************************!  
  !  
  ! *** ADD BODY FORCE TO ADVECTIVE ACCELERATIONS
  ! *** DISTRIBUTE UNIFORMLY OVER ALL LAYERS IF ISBODYF=1
  ! *** DISTRIBUTE OVER SURFACE LAYER IF ISBODYF=2
  !  
  !----------------------------------------------------------------------!  
  IF( ISBODYF == 1 )THEN  
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          FX(L,K) = FX(L,K) - DYU(L)*HU(L)*FBODYFX(L,K)
          FY(L,K) = FY(L,K) - DXV(L)*HV(L)*FBODYFY(L,K)
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF  

  IF( ISBODYF == 2 )THEN  
    !$OMP DO PRIVATE(ND,LF,LL,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      DO LP=LF,LL
        L=LWET(LP)
        FX(L,KC) = FX(L,KC) - DZIC(L,KC)*DYU(L)*HU(L)*FBODYFX(L,KC)
        FY(L,KC) = FY(L,KC) - DZIC(L,KC)*DXV(L)*HV(L)*FBODYFY(L,KC)
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF

  !**********************************************************************!  
  ! ** ADD EXPLICIT NONHYDROSTATIC PRESSURE   
  IF( KC > 1 .AND. ISPNHYDS >= 1 )THEN  
    !$OMP DO PRIVATE(ND,L,LP,K,TMPVAL) 
    DO ND=1,NDM  
      DO LP=1,LLWET(KC,ND)
        L=LKWET(LP,KC,ND)  
        TMPVAL = 2./(DZC(L,KSZ(L))+DZC(L,KSZ(L)+1))  
        DZPC(L,1) = TMPVAL*(PNHYDS(L,2)-PNHYDS(L,1))  
      ENDDO
    ENDDO  
    !$OMP END DO

    !$OMP DO PRIVATE(ND,L,LP,TMPVAL) 
    DO ND=1,NDM  
      DO LP=1,LLWET(KC,ND)
        L=LKWET(LP,KC,ND)  
        TMPVAL = 2./(DZC(L,KC)+DZC(L,KC-1))  
        DZPC(L,KC) = TMPVAL*(PNHYDS(L,KC)-PNHYDS(L,KC-1))  
      ENDDO
    ENDDO  
    !$OMP END DO

    IF( KC >= 3 )THEN  
      !$OMP DO PRIVATE(ND,K,LP,L,TMPVAL) 
      DO ND=1,NDM  
        DO K=2,KS  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            TMPVAL = 2./(DZC(L,K+1)+2.*DZC(L,K)+DZC(L,K-1) )  
            DZPC(L,K) = TMPVAL*(PNHYDS(L,K+1)-PNHYDS(L,K-1))  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO
    ENDIF  

   !$OMP DO PRIVATE(ND,K,LP,L,LS,LW,DZPU,DZPV)
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)  
          LW = LWC(L)
          LS = LSC(L)  
          DZPU = 0.5*(DZPC(L,K)+DZPC(LW,K))  
          DZPV = 0.5*(DZPC(L,K)+DZPC(LS,K))  
          FX(L,K) = FX(L,K) + SUB3D(L,K)*DYU(L)*( HU(L)*(PNHYDS(L,K)-PNHYDS(LW,K) ) - ( BELV(L)-BELV(LW) + ZZ(L ,K)*(HP(L)-HP(LW)) )*DZPU )
          FY(L,K) = FY(L,K) + SVB3D(L,K)*DXV(L)*( HV(L)*(PNHYDS(L,K)-PNHYDS(LS,K) ) - ( BELV(L)-BELV(LS) + ZZ(LS,K)*(HP(L)-HP(LS)) )*DZPV )
        ENDDO
      ENDDO  
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF

  IF( IATMP > 0 )THEN
    ! *** ADD AIR PRESSURE GRADIENT
    !$OMP DO PRIVATE(ND,K,LP,L,LS,LW)
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LW=LWC(L)
          LS=LSC(L)
          !                                M       M         M2/S2        
          FX(L,K) = FX(L,K) + SUB3D(L,K)*DYU(L)*HU(L)*( ATMP(L)-ATMP(LW) )
          FY(L,K) = FY(L,K) + SVB3D(L,K)*DXV(L)*HV(L)*( ATMP(L)-ATMP(LS) )
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF

  !----------------------------------------------------------------------!  
  ! **  ADD NET WAVE REYNOLDS STRESSES TO EXTERNAL ADVECTIVE ACCEL.  
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN

    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L) 
    DO ND=1,NDM  
      LF=(ND-1)*NWVCELLS+1  
      LL=MIN(LF+NWVCELLS-1,NWVCELLS)
      
      DO K=1,KC  
        DO LP=LF,LL
          L = LWVCELL(LP)  
          IF( LKSZ(L,K) )CYCLE  
          FX(L,K) = FX(L,K)+SUB3D(L,K)*WVFACT*SAAX(L)*FXWAVE(L,K)
          FY(L,K) = FY(L,K)+SVB3D(L,K)*WVFACT*SAAY(L)*FYWAVE(L,K)
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO
  ENDIF  
  
  !**********************************************************************!  
  ! **  CALCULATE TOTAL EXTERNAL ACCELERATIONS  
  !$OMP DO PRIVATE(ND,L,K,LP)
  DO ND=1,NDM  
      
    DO K=1,KC  
      DO LP=1,LLWET(K,ND)
        L = LKWET(LP,K,ND)  
        FCAXE(L) = FCAXE(L) + FCAX(L,K)*SGZU(L,K)  
        FCAYE(L) = FCAYE(L) + FCAY(L,K)*SGZV(L,K)  
        FXE(L) = FXE(L) + FX(L,K)*SGZU(L,K)  
        FYE(L) = FYE(L) + FY(L,K)*SGZV(L,K)  
      ENDDO  
    ENDDO
  ENDDO   ! *** END OF DOMAIN
  !$OMP END DO

  !**********************************************************************!  
  ! **  COMPLETE CALCULATION OF INTERNAL MODE ADVECTIVE ACCELERATIONS  
  IF( KC > 1 )THEN

    !**********************************************************************!  
    ! *** ADD VERTICAL MOMENTUM FLUX INTO HORIZONTAL MOMENTUM FLUX
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWETZ(K,ND)
          L = LKWETZ(LP,K,ND)
          FX(L,K) = FX(L,K) + SAAX(L)*( FWU(L,K)-FWU(L,K-1) )*DZIC(L,K)
          FY(L,K) = FY(L,K) + SAAY(L)*( FWV(L,K)-FWV(L,K-1) )*DZIC(L,K)
        ENDDO
      ENDDO  
    ENDDO   ! *** END OF DOMAIN
    !$OMP END DO

    !**********************************************************************!  
    ! **  ADD SUBGRID SCALE CHANNEL VIRTURAL MOMENTUM SOURCES AND SINKS  
    IF( MDCHH >= 1 .AND. ISCHAN == 3 )THEN  
      !$OMP DO PRIVATE(ND,LF,LL,L,LP,K)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)

        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            QMCSOURX(L,K) = 0.  
            QMCSOURY(L,K) = 0.  
            QMCSINKX(L,K) = 0.  
            QMCSINKY(L,K) = 0.  
          ENDDO  
        ENDDO  
      ENDDO  
      !$OMP END DO

      !$OMP SINGLE
      DO NMD=1,MDCHH  

        LHOST = LMDCHH(NMD)  
        LCHNU = LMDCHU(NMD)  
        LCHNV = LMDCHV(NMD)  

        DETH = CUE(LHOST)*CVN(LHOST)-CUN(LHOST)*CVE(LHOST)  
        CI11H = CVN(LHOST)/DETH  
        CI12H = -CUN(LHOST)/DETH  
        CI21H = -CVE(LHOST)/DETH  
        CI22H = CUE(LHOST)/DETH  

        DETU = CUE(LCHNU)*CVN(LCHNU)-CUN(LCHNU)*CVE(LCHNU)  
        CI11U = CVN(LCHNU)/DETU  
        CI12U = -CUN(LCHNU)/DETU  
        CI21U = -CVE(LCHNU)/DETU  
        CI22U = CUE(LCHNU)/DETU  

        DETV = CUE(LCHNV)*CVN(LCHNV)-CUN(LCHNV)*CVE(LCHNV)  
        CI11V = CVN(LCHNV)/DETV  
        CI12V = -CUN(LCHNV)/DETV  
        CI21V = -CVE(LCHNV)/DETV  
        CI22V = CUE(LCHNV)/DETV  

        ! *** X-DIRECTION CHANNEL  
        IF( MDCHTYP(NMD) == 1 )THEN  
          IF( QCHANU(NMD) > 0.0 )THEN  
            DO K=1,KC  
              QMCSINKX(LCHNU,K) = QMCSINKX(LCHNU,K)-0.5*DZC(L,K)*QCHANU(NMD)*(U(LCHNU,K)+U(LCHNU+1,K))  
              QMCSINKY(LCHNU,K) = QMCSINKY(LCHNU,K)-0.5*DZC(L,K)*QCHANU(NMD)*(V(LCHNU,K)+V(LNC(LCHNU),K))  
            ENDDO  
            DO K=1,KC  
              TMPVEC1(K) = CUE(LCHNU)*QMCSINKX(LCHNU,K)+CVE(LCHNU)*QMCSINKY(LCHNU,K)  
              TMPVEC2(K) = CUN(LCHNU)*QMCSINKX(LCHNU,K)+CVN(LCHNU)*QMCSINKY(LCHNU,K)  
            ENDDO  
            DO K=1,KC  
              QMCSOURX(LHOST,K) = QMCSOURX(LHOST,K)+CI11H*TMPVEC1(K)+CI12H*TMPVEC2(K)  
              QMCSOURY(LHOST,K) = QMCSOURY(LHOST,K)+CI21H*TMPVEC1(K)+CI22H*TMPVEC2(K)  
            ENDDO  
          ELSE  
            DO K=1,KC  
              QMCSINKX(LHOST,K) = QMCSINKX(LHOST,K)+0.5*DZC(L,K)*QCHANU(NMD)*(U(LHOST,K)+U(LHOST+1,K))  
              QMCSINKY(LHOST,K) = QMCSINKY(LCHNU,K)+0.5*DZC(L,K)*QCHANU(NMD)*(V(LHOST,K)+V(LNC(LHOST),K))  
            ENDDO  
            DO K=1,KC  
              TMPVEC1(K) = CUE(LHOST)*QMCSINKX(LHOST,K)+CVE(LHOST)*QMCSINKY(LHOST,K)  
              TMPVEC2(K) = CUN(LHOST)*QMCSINKX(LCHNU,K)+CVN(LHOST)*QMCSINKY(LHOST,K)  
            ENDDO  
            DO K=1,KC  
              QMCSOURX(LCHNU,K) = QMCSOURX(LCHNU,K)-CI11U*TMPVEC1(K)-CI12U*TMPVEC2(K)  
              QMCSOURY(LCHNU,K) = QMCSOURY(LCHNU,K)-CI21U*TMPVEC1(K)-CI22U*TMPVEC2(K)  
            ENDDO  
          ENDIF  
        ENDIF  

        ! *** Y-DIRECTION CHANNEL  
        IF( MDCHTYP(NMD) == 2 )THEN  
          IF( QCHANV(NMD) > 0.0 )THEN  
            DO K=1,KC  
              QMCSINKX(LCHNV,K) = QMCSINKX(LCHNV,K)-0.5*DZC(L,K)*QCHANV(NMD)*(U(LCHNV,K)+U(LCHNV+1,K))  
              QMCSINKY(LCHNV,K) = QMCSINKY(LCHNV,K)-0.5*DZC(L,K)*QCHANV(NMD)*(V(LCHNV,K)+V(LNC(LCHNV),K))  
            ENDDO  
            DO K=1,KC  
              TMPVEC1(K) = CUE(LCHNV)*QMCSINKX(LCHNV,K)+CVE(LCHNV)*QMCSINKY(LCHNV,K)  
              TMPVEC2(K) = CUN(LCHNV)*QMCSINKX(LCHNV,K)+CVN(LCHNV)*QMCSINKY(LCHNV,K)  
            ENDDO  
            DO K=1,KC  
              QMCSOURX(LHOST,K) = QMCSOURX(LHOST,K)+CI11H*TMPVEC1(K)+CI12H*TMPVEC2(K)  
              QMCSOURY(LHOST,K) = QMCSOURY(LHOST,K)+CI21H*TMPVEC1(K)+CI22H*TMPVEC2(K)  
            ENDDO  
          ELSE  
            DO K=1,KC  
              QMCSINKX(LHOST,K) = QMCSINKX(LHOST,K)+0.5*DZC(L,K)*QCHANV(NMD)*(U(LHOST,K)+U(LHOST+1,K))  
              QMCSINKY(LHOST,K) = QMCSINKY(LCHNV,K)+0.5*DZC(L,K)*QCHANV(NMD)*(V(LHOST,K)+V(LNC(LHOST),K))  
            ENDDO  
            DO K=1,KC  
              TMPVEC1(K) = CUE(LHOST)*QMCSINKX(LHOST,K)+CVE(LHOST)*QMCSINKY(LHOST,K)  
              TMPVEC2(K) = CUN(LHOST)*QMCSINKX(LCHNU,K)+CVN(LHOST)*QMCSINKY(LHOST,K)  
            ENDDO  
            DO K=1,KC  
              QMCSOURX(LCHNV,K) = QMCSOURX(LCHNV,K)-CI11V*TMPVEC1(K)-CI12V*TMPVEC2(K)  
              QMCSOURY(LCHNV,K) = QMCSOURY(LCHNV,K)-CI21V*TMPVEC1(K)-CI22V*TMPVEC2(K)  
            ENDDO  
          ENDIF  
        ENDIF  

      ENDDO  

      DO K=1,KC
        DO L=2,LA
          LE = LEC(L)
          LN = LNC(L)  
          IF( QMCSOURX(L,K) /= 0.0 )THEN  
            TMPVAL = SUB3D(L,K)+SUB(LE)  
            TMPVAL = MAX(TMPVAL,1.0)  
            FX(L,K)  = FX(L,K)  - SUB3D(L,K)*QMCSOURX(L,K)/TMPVAL  
            FX(LE,K) = FX(LE,K) - SUB3D(LE,K)*QMCSOURX(L,K)/TMPVAL  
          ENDIF  
          IF( QMCSOURY(L,K) /= 0.0 )THEN  
            TMPVAL = SVB(L)+SVB(LN)  
            TMPVAL = MAX(TMPVAL,1.0)  
            FY(L,K)  = FY(L,K)  - SVB3D(L ,K)*QMCSOURX(L,K)/TMPVAL  
            FY(LN,K) = FY(LN,K) - SVB3D(LN,K)*QMCSOURX(L,K)/TMPVAL  
          ENDIF  
          IF( QMCSINKX(L,K) /= 0.0 )THEN  
            TMPVAL = SUB3D(L,K)+SUB(LE)  
            TMPVAL = MAX(TMPVAL,1.0)  
            FX(L,K)  = FX(L,K)  - SUB3D(L,K)*QMCSINKX(L,K)/TMPVAL  
            FX(LE,K) = FX(LE,K) - SUB3D(LE,K)*QMCSINKX(L,K)/TMPVAL  
          ENDIF  
          IF( QMCSINKY(L,K) /= 0.0 )THEN  
            TMPVAL = SVB(L)+SVB(LNC(L))  
            TMPVAL = MAX(TMPVAL,1.0)  
            FY(L,K)  = FY(L,K)  - SVB3D(L,K)*QMCSINKX(L,K)/TMPVAL  
            FY(LN,K) = FY(LN,K) - SVB3D(LN,K)*QMCSINKX(L,K)/TMPVAL  
          ENDIF  
        ENDDO  
      ENDDO  
      !$OMP END SINGLE

    ENDIF  ! *** END OF CHANNEL MODIFIER SECTION
    
  ENDIF    ! *** END OF KC > 1 SECTION
  
  !**********************************************************************!  
  !  
  ! **  CALCULATE EXPLICIT INTERNAL BUOYANCY FORCINGS CENTERED AT N FOR  
  ! **  THREE TIME LEVEL STEP AND AT (N+1/2) FOR TWO TIME LEVEL STEP  
  ! **  SBX=SBX*0.5*DYU & SBY=SBY*0.5*DXV  
  !  
  !----------------------------------------------------------------------!  

  ! *** ORIGINAL  
  IF( BSC > 1.E-6 .AND. KC > 1 )THEN

    IF( IGRIDV == 1 )THEN
      ! *** SIGMA-ZED BOUYANCY SHEARS
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      DO ND=1,NDM  
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            LS=LSC(L) 
            LW=LWC(L)
            FBBX(L,K) = SUB3D(L,K)*SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*SGZU(L,K+1) + (B(L,K)-B(LW,K))*SGZU(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = SVB3D(L,K)*SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*SGZV(L,K+1) + (B(L,K)-B(LS,K))*SGZV(L,K) )                                &
                                                    - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO
      
    ELSEIF( IGRIDV == 2 )THEN
      ! *** SIGMA-ZED BOUYANCY SHEARS
      
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      DO ND=1,NDM
        ! *** ALL LAYERS
        DO K=1,KS
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND) 
            LS = LSC(L) 
            LW = LWC(L)
            FBBX(L,K) = SUB3D(L,K)*SBX(L)*GP*HU(L)*( BW(L,K+1)*HPW(L)*SGZW(L,K+1) - BE(LW,K+1)*HPE(LW)*SGZE(LW,K+1) + BW(L,K)*HPW(L)*SGZW(L,K) - BE(LW,K)*HPE(LW)*SGZE(LW,K)   &
                                                    - (BW(L,K+1)-BW(L,K)+BE(LW,K+1)-BE(LW,K))*( BELVW(L)+ZW(L,K)*HPW(L) - (BELVE(LW)+ZE(LW,K)*HPE(LW)) ) )  
            FBBY(L,K) = SVB3D(L,K)*SBY(L)*GP*HV(L)*( BS(L,K+1)*HPS(L)*SGZS(L,K+1) - BN(LS,K+1)*HPN(LS)*SGZN(LS,K+1) + BS(L,K)*HPS(L)*SGZS(L,K) - BN(LS,K)*HPN(LS)*SGZN(LS,K)   &
                                                    - (BS(L,K+1)-BS(L,K)+BN(LS,K+1)-BN(LS,K))*( BELVS(L)+ZS(L,K)*HPS(L) - (BELVN(LS)+ZN(LS,K)*HPN(LS)) ) )  
          ENDDO  
        ENDDO
      ENDDO
      !$OMP END DO
      
    ELSEIF( IINTPG == 0 )THEN  
      ! *** IINTPG=0  
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L,LS,LW) 
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
      
        DO K=1,KS  
          DO LP=LF,LL
            L=LWET(LP)  
            LS=LSC(L)  
            LW=LWC(L)
            FBBX(L,K) = SBX(L)*GP*HU(L)*(HU(L)*( (B(L,K+1)-B(LW,K+1))*DZCK(K+1) + (B(L,K)-B(LW,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K))*(BELV(L)-BELV(LW)+Z(L,K)*(HP(L)-HP(LW))) )  
            FBBY(L,K) = SBY(L)*GP*HV(L)*(HV(L)*( (B(L,K+1)-B(LS,K+1))*DZCK(K+1) + (B(L,K)-B(LS,K))*DZCK(K) ) - (B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K))*(BELV(L)-BELV(LS)+Z(L,K)*(HP(L)-HP(LS))) )  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO

    ELSEIF( IINTPG == 1. )THEN
      ! *** JACOBIAN
      
      IF( KC <= 2 )THEN
        ! *** KC  <= 2 LAYERS
        
        !$OMP DO PRIVATE(ND,LF,LL,L,LP,K,LS,LW) 
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
      
          DO K=1,KS
            DO LP=1,LLWET(K,ND)
              L = LKWET(LP,K,ND)  
              LS = LSC(L)
              LW = LWC(L)
              FBBX(L,K) = SBX(L)*GP*HU(L)*( HU(L)*( (B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+(B(L,K)-B(LW,K))*DZC(L,K) )-( B(L,K+1)-B(L,K)+B(LW,K+1)-B(LW,K) )*(BELV(L)-BELV(LW)+Z(L,K)*(HP(L)-HP(LW))) )
              FBBY(L,K) = SBY(L)*GP*HV(L)*( HV(L)*( (B(L,K+1)-B(LS,K+1))*DZC(L,K+1)+(B(L,K)-B(LS,K))*DZC(L,K) )-( B(L,K+1)-B(L,K)+B(LS,K+1)-B(LS,K) )*(BELV(L)-BELV(LS)+Z(L,K)*(HP(L)-HP(LS))) )
            ENDDO
          ENDDO

        ENDDO   ! *** END OF DOMAIN
        !$OMP END DO
      
      ELSE
        ! *** KC  > 2 LAYERS
        
        ! *** BOTTOM LAYER
        !$OMP DO PRIVATE(ND,LF,LL,L,LP,K,LS,LW) 
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
      
          DO LP=LF,LL
            L = LWET(LP)  
            LW = LWC(L)
            LS = LSC(L)
            K = KSZ(L)
            FBBX(L,K) = SBX(L)*GP*HU(L)*( 0.5*HU(L)*( 0.5*(B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+2.0*(B(L,K  )-B(LW,K  ))*DZC(L,K  ) )-0.25*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.50*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW))) )
                      
            FBBY(L,K) = SBY(L)*GP*HV(L)*( 0.5*HV(L)*( 0.5*(B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+2.0*(B(L,K  )-B(LS ,K  ))*DZC(L,K  ) )-0.25*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))  &
                        *(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))  -0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))  &
                        *(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS))) )
          ENDDO
        ENDDO   ! *** END OF DOMAIN
        !$OMP END DO
      
        ! *** LAYER AT KC-1
        !$OMP DO PRIVATE(ND,LF,LL,L,LP,K,LS,LW) 
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
      
          K=KS
          DO LP=LF,LL
            L = LWET(LP)  
            LW = LWC(L)
            LS = LSC(L)
            FBBX(L,K) = SBX(L)*GP*HU(L)*( 0.5*HU(L)*( 2.0*(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LW,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) )-0.50*(B(L,K+1)-B(L,K+1)+B(LW,K+1)-B(LW,K+1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW))) - 0.25*(B(L,K)-B(L,K-1)+B(LW,K)-B(LW,K-1))  &
                        *(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
                      
            FBBY(L,K) = SBY(L)*GP*HV(L)*( 0.5*HV(L)*( 2.0*(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LS ,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) )-0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))  &
                        *(BELV(L)-BELV(LS)+Z(L,K  )*(HP(L)-HP(LS)))   - 0.25*(B(L,K)-B(L,K-1)+B(LS,K)-B(LS ,K-1))   &
                        *(BELV(L)-BELV(LS)+Z(L,K-1)*(HP(L)-HP(LS))) )
          ENDDO
        ENDDO   ! *** END OF DOMAIN
        !$OMP END DO

        IF( KC > 3 )THEN
        
          ! *** MIDDLE LAYERS

          !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
          DO ND=1,NDM  
            DO K=2,KS-1
              DO LP=1,LLWET(K-1,ND)
                L = LKWET(LP,K-1,ND) 
                LW = LWC(L)
                LS = LSC(L)
                FBBX(L,K) = SBX(L)*GP*HU(L)*( 0.5*HU(L)*( 0.5*(B(L,K+2)-B(LW,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LW,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LW,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LW,K-1))*DZC(L,K-1) ) - 0.25*(B(L,K+2)-B(L,K+1)+B(LW,K+2)-B(LW,K+1))  &
                           *(BELV(L)-BELV(LW)+Z(L,K+1)*(HP(L)-HP(LW)))-0.50*(B(L,K+1)-B(L,K  )+B(LW,K+1)-B(LW,K  )) & 
                           *(BELV(L)-BELV(LW)+Z(L,K  )*(HP(L)-HP(LW)))-0.25*(B(L,K  )-B(L,K-1)+B(LW,K  )-B(LW,K-1)) &
                           *(BELV(L)-BELV(LW)+Z(L,K-1)*(HP(L)-HP(LW))) )
                         
                FBBY(L,K) = SBY(L)*GP*HV(L)*( 0.5*HV(L)*( 0.5*(B(L,K+2)-B(LS ,K+2))*DZC(L,K+2)+1.5*(B(L,K+1)-B(LS ,K+1))*DZC(L,K+1)+1.5*(B(L,K)-B(LS ,K))*DZC(L,K)+0.5*(B(L,K-1)-B(LS ,K-1))*DZC(L,K-1) )-0.25*(B(L,K+2)-B(L,K+1)+B(LS ,K+2)-B(LS ,K+1))    &
                           *(BELV(L)-BELV(LS)+Z(L,K+1)*(HP(L)-HP(LS)))-0.50*(B(L,K+1)-B(L,K  )+B(LS ,K+1)-B(LS ,K  ))  &
                           *(BELV(L)-BELV(LS)+Z(L,K  )*(HP(L)-HP(LS)))-0.25*(B(L,K  )-B(L,K-1)+B(LS ,K  )-B(LS ,K-1))  &
                           *(BELV(L)-BELV(LS)+Z(L,K-1)*(HP(L)-HP(LS)))  )
              ENDDO
            ENDDO

          ENDDO   ! *** END OF DOMAIN
          !$OMP END DO

        ENDIF     ! *** END OF KC > 3
      ENDIF       ! *** END OF KC > 2
      ! *** END OF JACOBIAN SECTION

    ELSEIF( IINTPG == 2 )THEN  
      ! *** FINITE VOLUME  
      !$OMP DO PRIVATE(ND,K,LP,L,LS,LW) 
      DO ND=1,NDM  
        DO K=1,KS  
          DO LP=1,LLWET(K-1,ND)
            L = LKWET(LP,K-1,ND)  
            LW = LWC(L)
            LS = LSC(L)
              
            FBBX(L,K) = SBX(L)*GP*HU(L)*( ( HP(L)*B(L,K+1)-HP(LW)*B(LW,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LW)*B(LW,K  ) )*DZC(L,K  ) )  & 
                       -SBX(L)*GP*(BELV(L)-BELV(LW))*( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LW)*B(LW,K+1)-HP(LW)*B(LW,K) )                  &
                       -SBX(L)*GP*(HP(L)-HP(LW))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K)+HP(LW)*ZZ(L,K+1)*B(LW,K+1)-HP(LW)*ZZ(L,K)*B(LW,K) )
                     
            FBBY(L,K) = SBY(L)*GP*HV(L)*( ( HP(L)*B(L,K+1)-HP(LS )*B(LS ,K+1) )*DZC(L,K+1)+( HP(L)*B(L,K  )-HP(LS )*B(LS ,K  ) )*DZC(L,K  ) )  &
                       -SBY(L)*GP*(BELV(L)-BELV(LS ))*( HP(L)*B(L,K+1)-HP(L)*B(L,K)+HP(LS)*B(LS ,K+1)-HP(LS)*B(LS ,K) )                    &
                       -SBY(L)*GP*(HP(L)-HP(LS ))*( HP(L)*ZZ(L,K+1)*B(L,K+1)-HP(L)*ZZ(L,K)*B(L,K) +HP(LS)*ZZ(L,K+1)*B(LS ,K+1)-HP(LS)*ZZ(L,K)*B(LS ,K) )
          ENDDO  
        ENDDO  
    
      ENDDO   ! *** END OF DOMAIN
      !$OMP END DO

    ENDIF  
    
  ENDIF  ! *** END OF BOUYANCY 

  !**********************************************************************!
  !
  ! **  CALCULATE EXPLICIT INTERNAL U AND V SHEAR EQUATION TERMS
  !
  !----------------------------------------------------------------------!
  IF( KC > 1 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      
      ! *** COMPUTE THE INTERNAL SHEARS FOR THE LOWER LAYERS
      DO K=1,KS
        DO LP=1,LLWET(K,ND)
          L = LKWET(LP,K,ND)
          DU(L,K) = CDZFU(L,K)*( HU(L)*( U(L,K+1)-U(L,K) )*DELTI + DXYIU(L)*( FCAX(L,K+1)-FCAX(L,K)   + FBBX(L,K) + SNLT*(FX(L,K)-FX(L,K+1)) ) )   ! *** M2/S2
          DV(L,K) = CDZFV(L,K)*( HV(L)*( V(L,K+1)-V(L,K) )*DELTI + DXYIV(L)*( FCAY(L,K)  -FCAY(L,K+1) + FBBY(L,K) + SNLT*(FY(L,K)-FY(L,K+1)) ) )   ! *** M2/S2
        ENDDO
      ENDDO

      ! *** ADD WIND SHEAR TO THE KC/KS INTERFACE
      IF( NWSER > 0 )THEN
        DO LP=1,LLWET(KS,ND)
          L = LKWET(LP,KS,ND) 
          DU(L,KS) = DU(L,KS) - CDZUU(L,KS)*TSX(L)
          DV(L,KS) = DV(L,KS) - CDZUV(L,KS)*TSY(L)
        ENDDO
      ENDIF
    ENDDO    ! *** END OF DOMAIN
    !$OMP END DO

  ENDIF
  !$OMP END PARALLEL

  !**********************************************************************!
  RETURN
END
