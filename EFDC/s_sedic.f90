SUBROUTINE SEDIC
  
  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER :: CORE,I,INCORE,J,L,LL,M,K,NS,VAR_BED,FDIR,NWV,NSC
  INTEGER :: IWV,JWV,NSKIP
  CHARACTER(LEN=80)  :: STR_LINE
  CHARACTER(LEN=120) :: STR_120
  
  !PT- real values are written in DOUBLE PRECISION. 7/16/08
  DOUBLE PRECISION :: STWVHTMP,STWVTTMP,STWVDTMP
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: BDEN      !(INCORE,KB)
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: TAUTEMP   !(KB)
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:,:) :: PNEW      !(INCORE,KB,NSCM)
  DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:,:)   :: TSED0S    !(KB,INCORE)

  ! Reads in Initial Erosion data.
  ! Reads in Erosion Data for Newly deposited 
  ! Sediments.
  ! Calculates bed parameters.
  ! REVISION DATE :  May 24, 2006
  ! Craig Jones and Scott James
  !**************************************************************************
  ! Set Maximum Number of Size Clases to be read in
  ! for newly deposited sediment erosion rates NSICM 
  ! Open Input Files
  ! Craig Jones
  ! *************************************************************************
  ! 2016-12
  ! Rearranged SEDZLJ initialization and added parameters for better toxics 
  ! similations.  Deprecated NSEDFLUME = 99 (i.e. SEDZLJ toxics) 
  ! Toxics are handled by ISTRAN(5)>0.  NEQUIL no longer used
  ! Paul M. Craig
  
  WRITE(*,'(A)')'READING SEDZLJ FILES'
  
  !CALL SEDDATA !calls routine to convert SEDflume data into a form useable by this code
  OPEN(UNIT=10,FILE='erate.sdf')
  OPEN(UNIT=20,FILE='core_field.sdf')
  OPEN(UNIT=30,FILE='bed.sdf')
  
  ! Read in Sediment Transport Parameters
  ! VAR_BED  = 1 for variable sediment bed
  ! NCALC_BL > 0 for bedload  calculation
  READ (30,'(A80)') STR_LINE
  READ (30,*) VAR_BED, KB, NCALC_BL, SEDSTEP, SEDSTART, IHTSTRT, IMORPH, ISWNWAVE, MAXDEPLIMIT, HPMIN
  IF( HPMIN < 0.003 .OR. HPMIN >= 1.0 ) HPMIN = 0.25
  
  IF( SEDSTEP < TIDALP/REAL(NTSPTC) ) SEDSTEP = TIDALP/REAL(NTSPTC)
  DTSED = SEDSTEP
  DTSEDJ = REAL(DTSED,8)
  
  IF( SEDSTART <= TBEGIN ) SEDSTART=TBEGIN
  
  ! *** NSCM  = Maximum number of grainsize classes defined.  NSCM=NSED most times.
  ! *** ITBM  = Number of Sedflume Shear Categories
  ! *** NSICM = Maximum number of redeposited grainsize classes defined.  NSICM=NSCM most times.
  READ (30,'(A80)') STR_LINE 
  READ (30,'(A80)') STR_LINE    ! DATA READ BY SCANSEDZLJ; ITBM and NSICM.
  
  READ (30,'(A80)') STR_LINE
  READ (30,*)  ZBSKIN,TAUCONST,ISSLOPE,BEDLOAD_CUTOFF

  ! *** Median grain size for each size class
  READ (30,'(A80)') STR_LINE
  READ (30,*) (D50(K),K=1,NSCM)  

  ! *** Erosion
  READ (30,'(A80)') STR_LINE
  READ (30,*) (TCRE(K),K=1,NSCM) 

  ! *** Suspension
  READ (30,'(A80)') STR_LINE
  READ (30,*) (TCRSUS(K),K=1,NSCM)

  ! *** Settling Velocities
  READ (30,'(A80)') STR_LINE
  READ (30,*) (DWSIN(K),K=1,NSCM)

  ! Hydrophobic Contaminant Information
  IF( NSEDFLUME == 2 )THEN
    ! PMC - 2016-11-09 - DEPRECATED.  HANDLED BY CALTOX AND CALTOXB
    !READ (30,'(A80)') STR_LINE
    !READ (30,*) (KPART(K),K=1,NSCM)
    !READ (30,'(A80)') STR_LINE
    !READ (30,*) (DIFFCOFF(K),K=1,NSCM)
    !READ (30,'(A80)') STR_LINE
    !DO LL=1,KB
    !   READ(30,*) (PCONTEMP(K,LL),K=1,NSCM)
    !ENDDO
  ENDIF
  
  !**************************************************************************
  !Reading in Erate.sdf starting with the layer's thickness.
  READ (10,'(A80)') STR_LINE
  READ (10,*) TACTM !read in active layer multiplier
  ! Read in Initial Erosion Data
  
  IF( VAR_BED >= 1 )THEN    
    ! Variable Bed *************************************************

    ! Determine File Format
    READ(20,'(A120)')STR_120
    I = 1
    DO WHILE ( i < 117 )
      IF( STR_120(I:I+2) == 'DSI' .OR. STR_120(I:I+2) == 'dsi' )EXIT
      I = I + 1
    END DO
     
    CLOSE(20)
    OPEN(UNIT=20,FILE='core_field.sdf')
    IF( I < 117 )THEN
      ! *** DSI STANDARD
      READ(20,'(A80)') STR_LINE
      READ(20,'(A80)') STR_LINE
      READ(20,'(A80)') STR_LINE
      READ(20,*)INCORE !read the number of cores    
      READ(20,'(A80)') STR_LINE
      READ(20,'(A80)') STR_LINE
      DO L=2,LA
        READ(20,*)I,J,CORE
        NCORENO(I,J) = CORE
      ENDDO
    ELSE
      ! *** SNL STANDARD
      READ(20,*)INCORE !read the number of cores    
      DO J=JC,1,-1
        READ(20,'(120(I1,1X))')(NCORENO(I,J),I=1,IC)
      ENDDO
    ENDIF
     
    ALLOCATE(BDEN(INCORE,KBM))
    ALLOCATE(PNEW(INCORE,KBM,NSCM+1))
    ALLOCATE(SEDDENS(INCORE))
    ALLOCATE(TAUTEMP(INCORE,KBM))
    ALLOCATE(TSED0S(KB,INCORE))
    BDEN=0.0   
    PNEW=0.0
    SEDDENS=0.0 
    TAUTEMP=0.0
    TSED0S = 0.0
     
    !*************************************************************   
    DO CORE=1,INCORE  ! for each core of data      
      READ (10,'(A80)') STR_LINE
      READ(10,*)(TAUTEMP(CORE,K),K=1,KB)  ! read the critical shear stresses of the core
      READ (10,'(A80)') STR_LINE
      READ (10,*) (TSED0S(K,CORE),K=1,KB)
      READ (10,'(A80)') STR_LINE      
      READ(10,*)(BDEN(CORE,K),K=1,KB)     ! read in the bulk density of the core 
      READ (10,'(A80)') STR_LINE      
      READ(10,*) WATERDENS, SEDDENS(CORE) ! read in the water density and sediment solid's density
      READ (10,'(A80)') STR_LINE  
      DO K=1,KB
        READ(10,*)(PNEW(CORE,K,NS),NS=1,NSCM)
      ENDDO       
      READ (10,'(A80)') STR_LINE
      DO M=1,ITBM
        READ(10,*)TAULOC(M) !shear stress used to erode a portion of the core
        READ(10,*)(ERATETEMP(CORE,K,M),K=1,KB) !erosion rate for each layer subject to shear stress TAULOC
      ENDDO
    ENDDO     
    
  ELSE
  
    ! Constant Erosion in Horizontal ********************************  
    INCORE = 1
    ALLOCATE(BDEN(INCORE,KBM))
    ALLOCATE(PNEW(INCORE,KBM,NSCM+1))
    ALLOCATE(SEDDENS(INCORE))
    ALLOCATE(TAUTEMP(INCORE,KBM))
    ALLOCATE(TSED0S(KB,INCORE))
    BDEN=0.0   
    PNEW=0.0
    SEDDENS=0.0 
    TAUTEMP=0.0
    TSED0S = 0.0

    CORE=1    
    READ (10,'(A80)') STR_LINE
    READ(10,*)(TAUTEMP(CORE,LL),LL=1,KB) !read the critical shear stresses of the core
    READ (10,'(A80)') STR_LINE
    READ (10,*) (TSED0S(K,1),K=1,KB)
    READ (10,'(A80)') STR_LINE      
    READ(10,*)(BDEN(CORE,LL),LL=1,KB) !read in the bulk density of the core 
    READ (10,'(A80)') STR_LINE      
    READ(10,*) WATERDENS, SEDDENS(1) !read in the water density and sediment solid's density 
    READ (10,'(A80)') STR_LINE  
    DO LL=1,KB
      READ(10,*)(PNEW(CORE,LL,K),K=1,NSCM)
    ENDDO
     
    READ (10,'(A80)') STR_LINE
    DO K=1,ITBM
      READ(10,*)TAULOC(K) !shear stress used to erode a portion of the core
      READ(10,*)(ERATETEMP(CORE,LL,K),LL=1,KB) !erosion rate for each layer subject to shear stress TAULOC
    ENDDO
     
    DO L=2,LA
      I=IL(L)
      J=JL(L)
      NCORENO(I,J)=1
    ENDDO

  ENDIF

  DO L=2,LA
    I=IL(L)   ! *** I location as a function of L
    J=JL(L)   ! *** J location as a function of L
    
    IF( NCORENO(I,J) > 0 )THEN
      CORE = NCORENO(I,J)
      DO K=1,KB
        TAUCOR(K,L) = TAUTEMP(CORE,K)                 ! *** Critical shear stresses from cores
        DO M=1,ITBM
          ERATE(K,L,M) = ERATETEMP(CORE,K,M)          ! *** Set erosion rate to measured value
        ENDDO
        DO NS=1,NSCM
          PERSED(NS,K,L) = PNEW(CORE,K,NS)/100.0_8    ! *** Set mass fraction to measured value
        ENDDO

        ! *** Removed the BSC dependency since BSC is only for the water column, not sediment bed (PMC)

        ! *** Compute bed porosity and dry bulk density, depedning on whether the SEDZLJ input files have wet or dry density
        IF( .true. )THEN
          ! *** BDEN is dry bulk density in g/cm3
          PORBED(L,K) = 1.-BDEN(CORE,K)/SEDDENS(CORE)
          BULKDENS(K,L) = BDEN(CORE,K)                                             ! *** Dry Bulk Density (BULKDENS)
        ELSE
          ! *** BDEN is wet bulk density in g/cm3
          PORBED(L,K) = ( BDEN(1,K)-SEDDENS(CORE) ) / ( 1.-SEDDENS(CORE) )
          BULKDENS(K,L) = (1.-PORBED(L,K))*BDEN(CORE,K)                            ! *** Dry Bulk Density (BULKDENS)
        ENDIF
        IF( BULKDENS(K,L) <= 0. )THEN
          PRINT '(" INVALID BULK DENSITY FOR CORE",I5," AT LAYER = ",I3)',NCORENO(I,J),K
          STOP
        ENDIF
      ENDDO
    ENDIF
  ENDDO  
    
  ! *** Disable Non-Cohesives.  Handled by SEDZLJ and CALTRAN
  SNDBT = 0.0
  ISTRAN(7) = 0
  
  !**************************************************************************
  ! Set Initial Layer and Thickness Values
  FORALL(L=2:LA)
    WHERE(TSED0S(1:KB,NCORENO(IL(L),JL(L))) > 0.0 )
      LAYERACTIVE(1:KB,L)=1
    ELSEWHERE
      LAYERACTIVE(1:KB,L)=0
    ENDWHERE
    FORALL(K=1:KB)
      TSED(K,L)  = TSED0S(K,NCORENO(IL(L),JL(L)))*BULKDENS(K,L)  ! PT TSED  in units of (g/cm^2).
      TSED0(K,L) = TSED0S(K,NCORENO(IL(L),JL(L)))*BULKDENS(K,L)  ! PT TSED0 in units of (g/cm^2).
    ENDFORALL
  ENDFORALL
  
  ! *** PRE-PROCESS THE DATA TO ENSURE VALID ASSIGNMENTS
  DO L=2,LA
    DO K=1,KB
      IF( TSED0(K,L)/BULKDENS(K,L) < 1.E-8 .OR. K <= 2 )THEN
        HBED(L,K) = 0.0
        TSED(K,L)  = 0.0
        TSED0(K,L) = 0.0
        TAUCOR(K,L) = 1000.
        DO M=1,ITBM
          ERATE(K,L,M) = 0.
        ENDDO
        DO NS=1,NSCM
          PERSED(NS,K,L) = 1./FLOAT(NSCM)
        ENDDO
      ENDIF
    ENDDO
  ENDDO

  ! *** GET HARD BOTTOM ELEVATIONS AND TOTAL SEDIMENT THICKNESS
  DO L=2,LA
    TSET0T(L)=0.0
    HBEDA(L)=0.0
    KBT(L)=-1                          ! *** Initialize KBT to the first layer with mass
    DO K=1,KB                          ! *** Topdown layer loop
      IF( TSED0(K,L) > 0. .AND. KBT(L) == -1 .AND. K > 2 ) KBT(L) = K
      TSET0T(L) = TSET0T(L) + TSED0(K,L)/BULKDENS(K,L)
      HBED(L,K) = 0.01*TSED(K,L)/BULKDENS(K,L)
      HBEDA(L)  = HBEDA(L) + HBED(L,K)
    ENDDO
    IF( KBT(L) == -1 ) KBT(L) = KB
    ZELBEDA(L) = BELV(L) - HBEDA(L)    ! *** Hard Bottom Elevation
  ENDDO
        
  ! *** Back to BED.SDF
  ! *** Read in Newly Deposited Sediments
  ! ***    Erosion Rates (ERATEND) and Critical Shear Stress for Erosion (TAUCRITE).
  READ (30,'(A80)') STR_LINE
  READ (30,*)  (SCND(NSC),NSC=1,NSICM)
  READ (30,'(A80)') STR_LINE
  READ (30,*)  (TAUCRITE(NSC),NSC=1,NSICM)
  READ (30,'(A80)') STR_LINE
  DO NSC=1,NSICM
     READ(30,*)(ERATEND(NSC,M),M=1,ITBM)
  ENDDO

  DO K=1,NSCM
    ! *** Removed the BSC dependency since BSC is only for the water column, not sediment bed (PMC)
    DISTAR(K)=D50(K)/10000.0*(((SEDDENS(1)/WATERDENS)-1.0)*980.0/0.01**2)**(1.0/3.0)

    ! Settling speed (DWS) calculated from D50 (micron) if DWS is < 0
    IF( DWSIN(K) == -1 )THEN
      ! input diameter using Cheng's model (1998).  Settling speed in cm/s
      DWS(K) = 0.01/(D50(K)*0.0001)*(SQRT(25.0+1.2*DISTAR(K)**2)-5.0)**1.5
    !ELSEIF( DWSIN(K) == -2 )THEN
      ! TBD
    ELSE
      DWS(K) = DWSIN(K)
    ENDIF
  ENDDO
  
  IF( IHTSTRT > 0 ) THEN
    WRITE(*,'(A)')'READING SEDBED_HOT.SDF'

    OPEN(514,FILE='SEDBED_HOT.SDF',FORM='FORMATTED',STATUS='old')
    READ(514,34569)((LAYERACTIVE(K,L),K=1,KB),L=2,LA)
    READ(514,34569)(KBT(L),L=2,LA)
    READ (514,34567)(D50AVG(L),L=2,LA)
    READ (514,34568)((BULKDENS(LL,L),LL=1,KB),L=2,LA)
    READ (514,34568)((TSED(LL,L),LL=1,KB),L=2,LA)
    READ (514,34568)(((PERSED(K,LL,L),K=1,NSED),LL=1,KB),L=2,LA)
    CLOSE(514)
    
    DO L=2,LA
      DO LL=1,KB
        TSED0(LL,L)=TSED(LL,L)
      END DO
    END DO
  ENDIF
  
34567 FORMAT(E17.9)
34568 FORMAT(6E17.9)
34569 FORMAT(8I8)

  !**************************************************************************
  ! Contaminant Transport Model
  IF( NSEDFLUME == 2 )THEN
    ! Initialize Contaminant Transport Variables
    STOP 'NSEDFLUME=2 Deprecated.  Handled by EFDC CALTOX and CALTOXB'
  ENDIF
  ! END Contaminant Transport
  !**************************************************************************
  
  ! Read in Wave Fetch or STWAVE Data if Used
  IF( ISWNWAVE == 1 )THEN
    WRITE(*,'(A)')'READING SEDZLJ: FETCH.INP'
    OPEN(UNIT=50,FILE='fetch.inp')
    DO L=2,LA
      READ (50,*) I,J,(FWDIR(LIJ(I,J),FDIR),FDIR=1,8)
    ENDDO
    CLOSE(50)
     
  ELSEIF (ISWNWAVE == 2 )THEN
    WRITE(*,'(A)')'READING SEDZLJ: STWAVE.INP'
    OPEN(UNIT=51,FILE='stwave.inp')
     
    DO NSKIP=1,5  
      READ (51,'(A80)') STR_LINE   
    ENDDO
     
    READ(51,*) STWVNUM,STWVTIM
     
    STWVTIM=STWVTIM*DT/3600.0
     
    DO NWV=1,STWVNUM  
      READ (51,'(A80)') STR_LINE
      !READ(51,*)IWV,JWV,STWVHT(2,NWV),STWVTP(2,NWV),STWVDR(2,NWV)
        
      DO L=2,LA
        READ(51,*)IWV,JWV,STWVHTMP,STWVTTMP,STWVDTMP
        STWVHT(LIJ(IWV,JWV),NWV)=STWVHTMP
        STWVTP(LIJ(IWV,JWV),NWV)=STWVTTMP
        STWVDR(LIJ(IWV,JWV),NWV)=STWVDTMP
        !STWVHT(L,NWV)=STWVHT(2,NWV)
        !STWVTP(L,NWV)=STWVTP(2,NWV)
        !STWVDR(L,NWV)=STWVDR(2,NWV)     
      ENDDO
      ! Incremental counter for which wave data set we are on
      STINC=0
      NWVCOUNT=STWVTIM-1
    ENDDO
     
    CLOSE(51)
  ENDIF
  
  !**************************************************************************
  
  FORALL(NS=1:NSCM) SSGI(NS) = 1.0/(1.0E6*SSG(NS))  !initialize SSGI.  Specific Volume (m**3/g)
  
  WRITE(6,*)'**************************************'
  WRITE(6,*)'Input Sizes (micron)'
  WRITE(6,*)(SCND(K),K=1,NSICM) 
  WRITE(6,*)'**************************************'
  WRITE(6,*)'Input Critical Shear Ero dynes/cm^2'
  WRITE(6,*)(TAUCRITE(K),K=1,NSICM) 
  WRITE(6,*)'**************************************'
  WRITE(6,*)'Sediment Sizes (micron)'
  WRITE(6,*)(D50(K),K=1,NSCM) 
  WRITE(6,*)'**************************************'
  WRITE(6,*)'DISTAR for each size class'
  WRITE(6,*)(DISTAR(K),K=1,NSCM) 
  WRITE(6,*)'**************************************'
  WRITE(6,*)'Critical Shear Sus dynes/cm^2 '
  WRITE(6,*)(TCRSUS(K),K=1,NSCM) 
  WRITE(6,*)'**************************************'
  WRITE(6,*)'Critical Shear Ero dynes/cm^2 '
  WRITE(6,*)(TCRE(K),K=1,NSCM) 
  WRITE(6,*)'**************************************'
  WRITE(6,*)'Settling Speeds in cm/s'
  WRITE(6,*)(DWS(K),K=1,NSCM) 
  WRITE(6,*)'**************************************'
  WRITE(6,*)'Surface Sediment Density  (g/cm^3) '
  WRITE(6,*)BULKDENS(3,3)
  WRITE(6,*)'**************************************'
  
  ! *** EFDCLOG.OUT
  WRITE(8,*)'**************************************'
  WRITE(8,*)'Input Sizes (micron)'
  WRITE(8,*)(SCND(K),K=1,NSICM) 
  WRITE(8,*)'**************************************'
  WRITE(8,*)'Input Critical Shear Ero dynes/cm^2'
  WRITE(8,*)(TAUCRITE(K),K=1,NSICM) 
  WRITE(8,*)'**************************************'
  WRITE(8,*)'Sediment Sizes (micron)'
  WRITE(8,*)(D50(K),K=1,NSCM) 
  WRITE(8,*)'**************************************'
  WRITE(8,*)'DISTAR for each size class'
  WRITE(8,*)(DISTAR(K),K=1,NSCM) 
  WRITE(8,*)'**************************************'
  WRITE(8,*)'Critical Shear Sus dynes/cm^2 '
  WRITE(8,*)(TCRSUS(K),K=1,NSCM) 
  WRITE(8,*)'**************************************'
  WRITE(8,*)'Critical Shear Ero dynes/cm^2 '
  WRITE(8,*)(TCRE(K),K=1,NSCM) 
  WRITE(8,*)'**************************************'
  WRITE(8,*)'Settling Speeds in cm/s'
  WRITE(8,*)(DWS(K),K=1,NSCM) 
  WRITE(8,*)'**************************************'
  WRITE(8,*)'Surface Sediment Density  (g/cm^3) '
  WRITE(8,*)BULKDENS(3,3)
  WRITE(8,*)'**************************************'

  CLOSE(10)
  CLOSE(20)
  CLOSE(30)
  CLOSE(40)
  
  ! *** ADD SMALL OFFSET FOR PRECISION ISSUES
  SCND(1)     = SCND(1)     - 1E-6
  SCND(NSICM) = SCND(NSICM) + 1E-6

  DO L=2,LA
    SH_SCALE(L) = 1.0
    IF( IHTSTRT == 0 ) THEN
      DO K=1,KB
        IF( HBED(L,K) > 0.0 )THEN
          D50AVG(L) = SUM(PERSED(1:NSCM,K,L)*D50(1:NSCM))     ! *** Calculate local d50 at sediment bed surface
          EXIT
        ENDIF
      ENDDO
      D50AVG(L) = MAX(D50AVG(L),D50(1))
    ENDIF
  ENDDO

  ! *** WSEDO - Fixed/Specified settling velocity.  Convert from cm/s (DWS) to m/s (WSEDO)
  WSEDO(1:NSCM) = DWS(1:NSCM)/100.0_8
  
  ! *** Bedload size cutoff for bedload (e.g. > 64) and approach for probability of deposition
  ! *** approach (Gessler > or Krone <)
  IF( BEDLOAD_CUTOFF < 10. ) BEDLOAD_CUTOFF = 64.
  
  ! *** INITIALIZE BED ARRAYS FOR USE BY SSEDTOX AND TOXICS SUBROUTINES CALTOX AND CALTOXB
  ! *** SEDZLJ DOES NOT HAVE CONSOLIDATION.  SO PORBED, PORBED1, VDRBED ARE TEMPORALLY CONSTANT
  DO L=2,LA
    FORALL(K=1:KB) HBED(L,K)   = 0.01*MAX(1.0E-12,TSED0S(K,NCORENO(IL(L),JL(L))))     ! *** Sediment bed layer thickness (m)
    FORALL(K=1:KB) VDRBED(L,K) = PORBED(L,K)/(1.0-PORBED(L,K))                        ! *** Sediment bed void ratio. (dimensionless)  

    FORALL(K=1:KB) SEDBT(L,K) = TSED(K,L)*10000.                                      ! *** Total sediment mass (g/m^2) in a layer, TSED-sediment layer unit mass (g/cm^2)
    FORALL(K=1:KB) SEDB(L,K,1:NSCM)  = SEDBT(L,K)*PERSED(1:NSCM,K,L)                  ! *** Sediment mass (g/m^2) by class in each layer
    FORALL(K=1:KB) SEDDIA50(L,K) = SUM(PERSED(1:NSCM,K,L)*D50(1:NSCM))                ! *** D50 for sediment layer.  

    FORALL(K=1:KB) HBED1(L,K)   = HBED(L,K)
    FORALL(K=1:KB) PORBED1(L,K) = PORBED(L,K)
    FORALL(K=1:KB) VDRBED1(L,K) = VDRBED(L,K)
    FORALL(K=1:KB) SEDB1(L,K,1:NSCM) = SEDB(L,K,1:NSCM)
  ENDDO
  SNDVDRD = BEDPORC/(1.-BEDPORC)                                                     ! *** Non-Cohesive settling void ratio (not used in SEDZLJ)
  
  DEALLOCATE(TAUTEMP,BDEN,PNEW)

  NNONCO = 1
  IF( NCALC_BL >= 1 )THEN
    !DEALLOCATE(QSBDLDOT)
    !ALLOCATE(QSBDLDOT(LCM,NSCM))
    !QSBDLDOT = 0.0
    DO NSC=1,NSCM
      IF( D50(NSC) > BEDLOAD_CUTOFF )THEN
        NNONCO = NSC
        EXIT
      ENDIF
    ENDDO
  ENDIF
  
  RETURN

  END SUBROUTINE SEDIC

