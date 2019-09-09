SUBROUTINE VARZEROSNL

  !C *** THIS SUBROUTINE ZERO'S MANY ARRAYS AFTER ALLOCATION
  
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 

  USE GLOBAL
  IMPLICIT NONE
  
  ! Begin SEDZLJ variables
  ALPHA_PX=0.0    !(LCM)
  ALPHA_PY=0.0    !(LCM)
  ALPHA_RX=0.0    !(LCM,NSCM)
  ALPHA_RY=0.0    !(LCM,NSCM)
  BLFLAG=0.0      !(LCM,NSCM)
  BULKDENS=0.0    !(KB,LCM)
  CBL=0.0         !(LCM,NSCM) 
  D50=0.0         !(NSCM)
  D50AVG=0.0      !(LCM)
  DEPO=0.0        !(LCM)
  DISTAR=0.0      !(NSCM)
  DZBL=0.0        !(LCM,NSCM)
  !DZBL_LAST=0.0   !(LCM,NSCM)
  DWS=0.0         !(NSCM)
  DWSIN=0.0       !(NSCM)
  ERATEND=0.0     !(NSICM,ITBM)
  ERATE=0.0       !(KB,LCM,ITBM)
  ERATETEMP = 0.0 !(INCORE,KB,ITBM)   Initialization only
  ESUS=0.0        !(LCM)
  ETOTO=0.0       !(LCM)
  HPCM=0.0        !(LCM)
  KPART=0.0       !(NTXM)
  PERSED=0.0      !(NSCM,KB,LCM)
  PSUS=0.0        !(LCM,NSCM)
  QBLFLUX=0.0       !(LCM,NSCM)
  SCND=0.0        !(NSICM)
  SH_SCALE=0.0    !(LCM)
  SSGI=0.0        !(NSCM) DSEDGMM array
  STWVHT=0.0      !(LCM,200)
  STWVTP=0.0      !(LCM,200)
  STWVDR=0.0      !(LCM,200)
  TAU=0.0         !(LCM)
  TAUCOR=0.0      !(KB,LCM)
  TAUCRITE=0.0    !(NSICM)
  TAULOC=0.0      !(ITBM)
  TCRE=0.0        !(NSCM)
  TCRSUS=0.0      !(NSCM)
  TRANS=0.0       !(LCM,NSCM)
  TSED=0.0        !(KB,LCM)
  TSED0=0.0       !(KB,LCM)
  TSEDT=0.0       !(LCM)
  TSET0T=0.0      !(LCM)
  USW=0.0         !(LCM,NSCM)
  IF( NCALC_BL > 0 )THEN
    BLVEL=0.0       !(LCM,NSCM)
    DBL=0.0         !(LCM,NSCM)
    EBL=0.0         !(LCM,NSCM)
    UBL=0.0         !(LCM,NSCM)
    VBL=0.0         !(LCM,NSCM)
    IF( ISTRAN(5) > 0 )THEN
      CBLTOX=0.0
      CBLTXCON=0.0
    ENDIF
  ENDIF
  
  IF( ISTRAN(5) > 0 .AND. NCALC_BL > 0 )THEN
    CBLTOX=0.     !(LCM,NTXM)
  ENDIF

  ! Begin SEDZLJ integer variables
  LAYERACTIVE=0.0 !(KB,LCM)
  NCORENO = 0

  ! *** SCALARS
  HPMIN=0.25      ! *** Minimum depth to compute shears
  MAXDEPLIMIT=0.0 ! *** The maximum fraction of mass from bottom water column layer that can be deposited on active bed layer.
  RHO=1000.0      ! *** Density of water in kg/m^3.
  TACTM=0.0       ! *** Active layer multiplier
  WATERDENS=1.0   ! *** Density of water in g/cm^3

  ! End SEDZLJ variables

  IF ( ISTRAN(8)==1 ) THEN
    ! Begin dissolved carbon dioxide variables
    CDOSATIDX=0.0
    WQCDOS=0.0
    WQITOP=0.0
    WQKRCDOS=0.0
    WQP22=0.0
    ! End dissolved carbon dioxide variables
  ENDIF
  

  ! Begin MHK variables SCJ
  CTMHK=0.0
  CDSUP=0.0
  DENMHK=0.0   ! *** DENSITY OF MHK DEVICES (#/CELL)
  EMHK=0.0     ! *** zero the allocatable arrays for MHK energy generation
  ESUP=0.0     ! *** zero the allocatable arrays for support energy dissipatio
  FXMHK=0.0    ! *** X"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK DEVICE
  FXMHKE=0.0   ! *** COLUMN X"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK DEVICE
  FYMHK=0.0    ! *** Y"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK DEVICE
  FYMHKE=0.0   ! *** COLUMN Y"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK DEVICE
  FXSUP=0.0    ! *** X"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK SUPPORT
  FXSUPE=0.0   ! *** COLUMN X"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK SUPPORT
  FYSUP=0.0    ! *** Y"FORCE" FOR HORIZONTAL ADVECTION CALCULATIONS FROM MHK SUPPORT
  FYSUPE=0.0   ! *** COLUMN Y"FORCE" FOR MOMENTUM CALCULATIONS FROM MHK DEVICE
  IJLTURB=0 
  PMHK=0.0     ! *** array that accumulates layer-wise MHK power in x direction
  PSUP=0.0     ! *** array that accumulates layer-wise MHK support power in x direction 
  VMAXCUT=0.0
  VMINCUT=0.0
  WIDTHMHK=0.0
  WIDTHSUP=0.0
  ZMAXMHK=0.0
  ZMAXSUP=0.0
  ZMINMHK=0.0
  ZMINSUP=0.0
  ! End MHK variables SCJ
  
END SUBROUTINE VARZEROSNL
