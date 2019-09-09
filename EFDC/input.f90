SUBROUTINE INPUT(TITLE)
                                                                  
  ! *** SUBROUTINE INPUT READS ALL INPUT DATA EXCEPT DATA IN LXLY.INP,                                                    
  ! *** MASK.INP AND RESTART.INP
  
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 

  USE GLOBAL
  USE RESTART_MODULE,ONLY:GENRSTINP
  USE INFOMOD,ONLY:SKIPCOM,READSTR
  USE HYDSTRUCMOD,ONLY: HYDSTRUCT_CHECK
  USE IFPORT
  
  IMPLICIT NONE

  CHARACTER*80 TEXT,TITLE,STR*200
  CHARACTER*10 CDUM
  CHARACTER*3  NCARD
  CHARACTER    ADUMMY*5,RESTARTF*50,STRC*650

  REAL :: SEEPRATE(1000)
  REAL,ALLOCATABLE,DIMENSION(:) :: RMULADS
  REAL,ALLOCATABLE,DIMENSION(:) :: ADDADS
  REAL,ALLOCATABLE,DIMENSION(:) :: BOFFMHK,BOFFSUP,TOFFMHK,TOFFSUP
  REAL,ALLOCATABLE,DIMENSION(:,:) :: PFX2
  REAL,ALLOCATABLE,DIMENSION(:,:) :: QSERSM,RULES

  INTEGER,ALLOCATABLE,DIMENSION(:) :: IPARTSP   ! *** NEW BEDFORD
  INTEGER,ALLOCATABLE,DIMENSION(:) :: IDX
  
  INTEGER :: ITURB=0,IPMC,IS,NS,NX,NT,NW,NA,NI
  INTEGER :: IDUM,NDUM,LCM2T,I,J,K,M,L,KDUM,IMDXDY,NPFORN
  INTEGER :: ISO,ITIDASM,NPFOR,NPFORS,NPFORW,NPFORE,ITSSS
  INTEGER :: NTXSRJP,NSDSRJP,NSNSRJP,NDUM1,NDUM2,NTMP,ITYPE
  INTEGER :: MTMP,NTOXSRQ,NSEDSRQ,NSNDSRQ,MMAX,MS,MMIN
  INTEGER :: NTOXSRC,NSEDSRC,NSNDSRC,JSFDCH,IDUMMY,NPP,MD,MU
  INTEGER :: MTSSS,IACROSS,JCTMP,JACROSS,JT,JF,JLAST,NMD
  INTEGER :: IFIRST,ILAST,IT,NP,LT,MVEGIJT,NMDXDY
  INTEGER :: ITMP,JTMP,LTMP,LL,ITMPU,JTMPU,ITMPD,JTMPD,NSEEPCLASSES
  INTEGER :: IZONE,LDUM,JDUM,IVAL,ISALTYP,IREAD,KBINPUT,IISTMP,ITXINTT
  INTEGER :: ISTYP,ISMOOTH,ICHGQS,NCTMP,NFSED,NTT,NFSND,NFTOX
  INTEGER :: LD,ID,JD,NFLAGPWR,ITMPPMX,NPMXZ,NPMXPTS,PMIXSF,NZ
  INTEGER :: NPBPH,IASERVER
  INTEGER :: ISQCTRL(6)
  
  REAL   :: ADMAX,ADMIN,AHMAX,AHMIN
  REAL   :: DXIJ,DYIJ,HIJ,BELVIJ,ZBRIJ,RVALUE,PSERTMP,DIAMMHK
  REAL   :: DXYCVT,HADADJ,RAD,AMP,T1,T2,TMPAMP,TMPPHS,BOTTOM
  REAL   :: TOP1,TOP2,QSSE,SEDVDRT,DSTR,USTR,DUM,QSERTMP
  REAL   :: RMDX,RMDY,CVTFACX,CVTFACY,FBODY1,FBODY2,BDLTMP
  REAL   :: RMULADJ,ADDADJ,RMULADJS,ADDADJS,PSERTMPS,CSERTMP
  REAL   :: QCTLTMP,SOLRCVT,CLDCVT,WINDSCT,RMULADJCOV,RMULADJTHK
  REAL,EXTERNAL :: SETSTVEL,PARSE_REAL

  ALLOCATE(RMULADS(NSTM))
  ALLOCATE(ADDADS(NSTM))
  ALLOCATE(QSERSM(NDQSER,KCM))
  ALLOCATE(PFX2(NPFORM,MTM))
  ALLOCATE(IPARTSP(NTXM))

  RMULADS=0.
  ADDADS=0.
  QSERSM=0.
                                                                  
  G=9.81
  PI=3.1415926535898
  PI2=2.*PI
  2 FORMAT(A80)                                                                                                       
                                                                  
  ! *** READ MAIN INPUT FILE EFDC.INP                                                                                     
                                                                  
  WRITE(*,'(A)')'READING THE MAIN EFDC CONTROL FILE: EFDC.INP' 
  OPEN(1,FILE='efdc.inp',STATUS='UNKNOWN')
                                                                  
  !1**  READ TITLE CARD                                                                                                   
  NCARD='1'
  CALL SEEK('C1')
  READ(1,2) TITLE
  WRITE(7,1002)NCARD
  WRITE(7,2) TITLE
  
  !C1A**  READ MODE OPTIONS                                                                                                
  NCARD='1A'
  CALL SEEK('C1A')
  READ(1,*,IOSTAT=ISO)IS2TIM,IGRIDH,IGRIDV,KMINV,SGZHPDELTA
  WRITE(7,1002)NCARD
  WRITE(7,*) IS2TIM,IGRIDH,IGRIDV,KMINV,SGZHPDELTA
  IF( ISO > 0 ) GOTO 100
                                                              
  !C2**  READ RESTART AND DIAGNOSTIC SWITCHES                                                                              
  NCARD='2'
  CALL SEEK('C2')
  READ(1,*,IOSTAT=ISO) ISRESTI,ISRESTO,ISRESTR,ISGREGOR,ISLOG,ISDIVEX,ISNEGH,ISMMC,ISBAL,ICONTINUE,ISHOW
  
  ! *** HANDLE BATHYMETRY ADJUSTMENTS
  ISRESTIOPT = 0
  IF( ISRESTI == -1 )THEN
    ISRESTIOPT = 1
    ISRESTI = 1
  ENDIF
  
  IF( ISRESTI == 1 .AND. ICONTINUE == 1 )THEN
    ! ** FOR CONTINUATION MODE:
    CALL SEEK('C2A')
    READ(1,*,IOSTAT=ISO) RESTARTF   
  ENDIF
  
  WRITE(7,1002)NCARD
  WRITE(7,*) ISRESTI,ISRESTO,ISRESTR,ISGREGOR,ISLOG,ISDIVEX,ISNEGH,ISMMC,ISBAL,ICONTINUE,ISHOW
  IF( ISMMC<0 )THEN
    DEBUG=.TRUE.
    ISMMC=0
    WRITE(*,'(A)')' DEBUG ON'
  ELSE
    DEBUG=.FALSE.
    WRITE(*,'(A)')' DEBUG OFF'
  ENDIF
  IF( ISO > 0 ) GOTO 100
  
  !C3**  READ RELAXATION PARAMETERS                                                                                        
  NCARD='3'
  CALL SEEK('C3')
  READ(1,*,IOSTAT=ISO) RP,RSQM,ITERM,IRVEC,IATMP,IWDRAG,NRAMPUP,ITERHPM,IDRYCK,ISDSOLV,FILT3TL
  WRITE(7,1002)NCARD
  WRITE(7,*) RP,RSQM,ITERM,IRVEC,IATMP,IWDRAG,NRAMPUP,ITERHPM,IDRYCK,ISDSOLV,FILT3TL
  ! *** Not Used: IWDRAG
  IF( ISO > 0 ) GOTO 100
  IF( IRVEC/= 0 .AND. IRVEC/= 9 .AND. IRVEC/=99 .AND. IRVEC/=9999 )STOP 'INVALID IRVEC'
                                                                  
  !C4**  READ LONGTERM MASS TRANSPORT INTEGRATION ONLY SWITCHES                                                            
  NCARD='4'
  CALL SEEK('C4')
  READ(1,*,IOSTAT=ISO) ISLTMT,ISSSMMT,ISLTMTS,ISIA,RPIA,RSQMIA,ITRMIA,ISAVEC
  WRITE(7,1002)NCARD
  WRITE(7,*) ISLTMT,ISSSMMT,ISLTMTS,ISIA,RPIA,RSQMIA,ITRMIA,ISAVEC
  IF( ISO > 0 ) GOTO 100
                                                                  
  !C5**  READ MOMENTUM ADVECTION AND DIFFUSION SWITCHES AND MISC                                                           
  NCARD='5'
  CALL SEEK('C5')
  READ(1,*,IOSTAT=ISO) ISCDMA,ISHDMF,ISDISP,ISWASP,ISDRY,ISQQ,ISRLID,ISVEG,ISVEGL,ISITB,IHMDSUB,IINTPG
  WRITE(7,1002)NCARD
  WRITE(7,*) ISCDMA,ISHDMF,ISDISP,ISWASP,ISDRY,ISQQ,ISRLID,ISVEG,ISVEGL,ISITB,IHMDSUB,IINTPG
  ! *** Not Used: IDRYCK
  IF( ISO > 0 ) GOTO 100
  IDRYTBP=0
  IF( ISDRY < 0 )THEN
    ISDRY=ABS(ISDRY)
    IDRYTBP=1
  ENDIF
  IF( ISVEG < 1 ) ISITB=0
  IF( ISWASP == 99 ) ISICM=1
  IF( ISRLID == 1  ) ISDRY=-1
  IF( ISWASP == 10 ) ISRCA=1
  JSWAVE=0
  ! PMC      IS1DCHAN=0                                                                                                   
  ! PMC  IF( ISCDMA == 10) IS1DCHAN=1                              Y                                                        
  ISCOSMIC=0
                                                                  
  !C6**  DISSOLVED AND SUSPENDED CONSTITUENT TRANSPORT SWITCHES                                                         
  NCARD='6'
  CALL SEEK('C6')
  DO NS=0,8
    READ(1,*,IOSTAT=ISO) ISTRAN(NS),ISTOPT(NS),ISCDCA(NS),ISADAC(NS),ISFCT(NS),ISPLIT(NS),ISADAH(NS),ISADAV(NS),ISCI(NS),ISCO(NS)
    IF( ISCDCA(NS) >= 4) ISCOSMIC=1  
    WRITE(7,1002)NCARD
    WRITE(7,*) ISTRAN(NS),ISTOPT(NS),ISCDCA(NS),ISADAC(NS),ISFCT(NS),ISPLIT(NS),ISADAH(NS),ISADAV(NS),ISCI(NS),ISCO(NS)
    IF( ISO > 0 ) GOTO 100
  ENDDO
  
  IF( ISTRAN(8) >= 1 .AND. ISTRAN(2) == 0 )THEN
    PRINT *,'*** WARNING: TEMPERATURE SHOULD BE ACTIVATED FOR WQ CALCULATIONS!'
    CALL SLEEPQQ(5000)
  ENDIF
  
  IF( ISRESTI == 1 .AND. ICONTINUE == 1 )THEN
    ! ** FOR CONTINUATION MODE:
    CLOSE(1)
    CALL GENRSTINP(RESTARTF)           ! GENERATE RESTART FILES, OTHERWISE EFDC LOADS EXISTING ONES
    OPEN(1,FILE='efdc.inp',STATUS='UNKNOWN')
  ENDIF
  
  !  *** DEACTIVATE ANY UNUSED OPTIONS
  IF( ISTRAN(2) < 1 ) ISTOPT(2) = 0
  
  ! *** SET TRANSPORT FLAG
  ISTRANACTIVE=0
  DO I=1,8
    IF( ISTRAN(I) > 0 )ISTRANACTIVE=1
  ENDDO
                                                                  
  !C7**  READ TIME-RELATED INTEGER PARAMETERS                                                                              
  NCARD='7'
  CALL SEEK('C7')
  READ(1,*,IOSTAT=ISO) NTC,NTSPTC,NLTC,NTTC,NTCPP,NTSTBC,NTCNB,NTCVB,NTSMMT,NFLTMT,NDRYSTP,NRAMPUP,NUPSTEP
  IF( NRAMPUP < 1 ) NRAMPUP=1
  IF( NUPSTEP < 2 ) NUPSTEP=2

  WRITE(7,1002)NCARD
  WRITE(7,*) NTC,NTSPTC,NLTC,NTTC,NTCPP,NTSTBC,NTCNB,NTCVB,NTSMMT,NFLTMT,NDRYSTP,NRAMPUP,NUPSTEP
  ! *** Not Used:  NTCPP
  ! *** Not Used:  NTCNB
  IF( ISO > 0 ) GOTO 100
                                                                  
  !C8**  READ TIME-RELATED REAL PARAMETERS                                                                                 
  NCARD='8'
  CALL SEEK('C8')
  READ(1,*,IOSTAT=ISO) TCON,TBEGIN,TIDALP,CF,ISCORV,ISDCCA,ISCFL,ISCFLM,DTSSFAC,DTSSDHDT,DTMAX
  
  IF( ISRESTI == 1 .AND. ICONTINUE == 1 )THEN
    ! ** FOR CONTINUATION MODE (TBEGINC IS SET FROM THE RESTART FILE)
    NTC = NTC + INT(TBEGIN - TBEGINC)
    TBEGIN = TBEGINC
  ENDIF
  
  WRITE(7,1002)NCARD
  WRITE(7,*) TCON,TBEGIN,TIDALP,CF,ISCORV,ISDCCA,ISCFL,ISCFLM,DTSSFAC,DTSSDHDT,DTMAX
  IF( ISO > 0 ) GOTO 100

  IF( DTSSFAC > 0.0 )THEN
    ISDYNSTP = 1
    DT = TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)
    IF( DTMAX <= DT*2. ) DTMAX = 3600.
  ELSE
    ISDYNSTP = 0
  ENDIF
  IF( IS2TIM == 0 ) ISDYNSTP = 0
                                                                  
  !C9**  READ SPACE RELATED AND SMOOTHING PARAMETERS                                                                       
  NCARD='9'
  CALL SEEK('C9')
  READ(1,*,IOSTAT=ISO) IC,JC,LC,LVC,ISCLO,NDM,LDM,ISMASK,ISCONNECT,NSHMAX,NSBMAX,WSMH,WSMB

  WRITE(7,1002)NCARD
  WRITE(7,*) IC,JC,LC,LVC,ISCLO,NDM,LDM,ISMASK,ISCONNECT,NSHMAX,NSBMAX,WSMH,WSMB
  IF( ISO > 0 ) GOTO 100
  IS2LMC=0
  IF( KC<0 )THEN
    KC=-KC
    IS2LMC=1
  ENDIF

  ! *** DOMAIN DECOMPOSITION CHECKS FOR HORIZONTAL LOOPS
  ! *** MAKE CONSISTENT WITH NEW OMP APPROACH
  LCM2T=LC-2
  NDM=NTHREADS
  LDM=INT(FLOAT(LCM2T)/FLOAT(NTHREADS))+1

  IF( KC >= 2) ISITB=0
                                                                  
  !C9A**  READ VERTICAL SPACE RELATED  PARAMETERS                                                                          
  NCARD='9A'
  CALL SEEK('C9A')
  READ(1,*,IOSTAT=ISO)KC,KSIG,ISETGVC,SELVREF,BELVREF,ISGVCCK
  WRITE(7,1002)NCARD
  WRITE(7,*) KC,KSIG,ISETGVC,SELVREF,BELVREF,ISGVCCK
  ! *** Not Used: KSIG
  IF( ISO > 0 ) GOTO 100
                                                                  
  !C10*  READ LAYER THICKNESS IN VERTICAL                                                                                  
  NCARD='10'
  CALL SEEK('C10')
  DO K=1,KC
    READ(1,*,IOSTAT=ISO)KDUM,DZCK(K)
    WRITE(7,1002)NCARD
    WRITE(7,*)KDUM,DZCK(K)
    IF( ISO > 0 ) GOTO 100
  ENDDO
                                                                  
  !C11*  READ GRID, ROUGHNESS, MASKING AND DEPTH PARAMETERS                                                                
  NCARD='11'
  CALL SEEK('C11')
  READ(1,*,IOSTAT=ISO) DX,DY,DXYCVT,IMDXDY,ZBRADJ,ZBRCVRT,HMIN,HADADJ,HCVRT,HDRY,HWET,BELADJ,BELCVRT

  WRITE(7,1002)NCARD
  WRITE(7,*) DX,DY,DXYCVT,IMDXDY,ZBRADJ,ZBRCVRT,HMIN,HADADJ,HCVRT,HDRY,HWET,BELADJ,BELCVRT
  IF( ISO > 0 ) GOTO 100
  HDRYICE = 0.91*HDRY
  HDRYWAV = 1.2*HWET
  
  !C11A* READ TWO-LAYER MOMENTUM FLUX AND CURVATURE ACCELERATION                                                           
  !     CORRECTION FACTORS                                                                                                
  NCARD='11A'
  CALL SEEK('C11A')
  READ(1,*,IOSTAT=ISO) ICK2COR,CK2UUM,CK2VVM,CK2UVM,CK2UUC,CK2VVC,CK2UVC,CK2FCX,CK2FCY
  WRITE(7,1002)NCARD
  WRITE(7,*) ICK2COR,CK2UUM,CK2VVM,CK2UVM,CK2UUC,CK2VVC,CK2UVC,CK2FCX,CK2FCY
  IF( ISO > 0 ) GOTO 100
  IF( ICK2COR >= 1 )THEN
    IS2LMC=ICK2COR
  END IF
                                                                  
  !C11B* READ CORNER CELL BOTTOM STRESS CORRECTION OPTIONS                                                                 
  NCARD='11B'
  CALL SEEK('C11B')
  READ(1,*,IOSTAT=ISO)ISCORTBC,ISCORTBCD,FSCORTBC
  WRITE(7,1002)NCARD
  WRITE(7,*) ISCORTBC,ISCORTBCD,FSCORTBC
  IF( ISO > 0 ) GOTO 100
                                                                  
  !C12*  READ TURBULENT DIFFUSION PARAMETERS                                                                               
  NCARD='12'
  CALL SEEK('C12')
  READ(1,*,IOSTAT=ISO) AHO,AHD,AVO,ABO,AVMX,ABMX,VISMUD,AVCON,ZBRWALL
  WRITE(7,1002)NCARD
  WRITE(7,*) AHO,AHD,AVO,ABO,AVMX,ABMX,VISMUD,AVCON,ZBRWALL
  IF( ISO > 0 ) GOTO 100
                                                                  
  !C12A*  READ TURBULENCE CLOSURE OPTIONS                                                                                  
  NCARD='12A'
  CALL SEEK('C12A')
  READ(1,*,IOSTAT=ISO)ISTOPT(0),ISSQL,ISAVBMX,ISFAVB,ISINWV,ISLLIM,IFPROX,XYRATIO,BC_EDGEFACTOR 

  WRITE(7,1002)NCARD
  WRITE(7,*)ISTOPT(0),ISSQL,ISAVBMX,ISFAVB,ISINWV,IFPROX,XYRATIO,BC_EDGEFACTOR
  IF( ISO > 0 ) GOTO 100
  IF( BC_EDGEFACTOR < 0 )BC_EDGEFACTOR=0.0
  IF( BC_EDGEFACTOR > 1)BC_EDGEFACTOR=1.0
                                                                  
  !C13*  READ TURBULENCE CLOSURE PARAMETERS                                                                                
  NCARD='13'
  CALL SEEK('C13')
  ! *** PMC - CTE2 NOT USED
  READ(1,*,IOSTAT=ISO) VKC,CTURB,CTURB2B,CTE1,CTE2,CTE3,CTE4,CTE5,RIQMAX,QQMIN,QQLMIN,DMLMIN

  WRITE(7,1002)NCARD
  WRITE(7,*) VKC,CTURB,CTURB2B,CTE1,CTE2,CTE3,CTE4,CTE5,RIQMAX,QQMIN,QQLMIN,DMLMIN
  IF( ISO > 0 ) GOTO 100
                                                                  
  !C14*  READ TIDAL & ATMOSPHERIC FORCING, GROUND WATER AND SUBGRID CHANNEL PARAMETERS                                                                                    
  NCARD='14'
  CALL SEEK('C14')
  READ(1,*,IOSTAT=ISO) MTIDE,NWSER,NASER,ISGWIT,ISCHAN,ISWAVE,ITIDASM,ISPERC,ISBODYF,ISPNHYDS

  WRITE(7,1002)NCARD
  WRITE(7,*) MTIDE,NWSER,NASER,ISGWIT,ISCHAN,ISWAVE,ITIDASM,ISPERC,ISBODYF,ISPNHYDS
  ISWCBL=0
  ISWVSD=0
  IF( ISO > 0 ) GOTO 100
  IF( ISPERC > 0 ) ISGWIT = 3   ! *** DS-INTL

  !C14A* READ SAND GRAIN NIKURADSE ROUGHNESS                                                                               
  KSW = 0.00001 * 2.5  ! *** DEFAULT IS 10 MICRON
  IF( ISWAVE >= 1 )THEN
    NCARD='14A'
    CALL SEEK('C14A')
    READ(1,*,IOSTAT=ISO) KSW,IUSEWVCELLS,IFWAVE,SWANGRP,ISSTEAD
    IF( ISO > 0 ) GOTO 100
    IF( ISWAVE == 1 .OR. ISWAVE == 2 .OR. ISWAVE == 4 )THEN
      NCARD='14B'
      CALL SEEK('C14B')
      READ(1,*,IOSTAT=ISO)ISWRSR,ISWRSI,WVDISV,WVLSH,WVLSX,ISWVSD,WVLCAL,NTSWV,ISWCBL,ISDZBR
      IF( ISO > 0 ) GOTO 100
      NTSWV = MAX(NTSWV,1)
    ENDIF
  ENDIF
                                                                  
  IF( MTIDE > 0 )THEN
    !C15*  READ PERIODIC FORCING (TIDAL) CONSTITUENT SYMBOLS AND PERIODS                                                     
    NCARD='15'
    CALL SEEK('C15')
    DO M=1,MTIDE
      READ(1,*,IOSTAT=ISO) SYMBOL(M),TCP(M)

      WRITE(7,1002)NCARD
      WRITE(7,*) SYMBOL(M),TCP(M)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
                                                                  
  !C16*  READ SURFACE ELEVATION OR PRESSURE BOUNDARY CONDITION PARAMETERS                                                  
  NCARD='16'
  CALL SEEK('C16')
  READ(1,*,IOSTAT=ISO) NPBS,NPBW,NPBE,NPBN,NPFOR,NPFORT,NPSER,PDGINIT

  WRITE(7,1002)NCARD
  WRITE(7,*) NPBS,NPBW,NPBE,NPBN,NPFOR,NPFORT,NPSER,PDGINIT
  IF( ISO > 0 ) GOTO 100
                                                                  
  IF( NPFOR > 0 )THEN
    !C17*  READ PERIODIC FORCING (TIDAL) SURFACE ELEVATION OR                                                                
    NCARD='17'
    CALL SEEK('C17')
    DO NP=1,NPFOR
      DO M=1,MTIDE
        IF( NPFORT == 0 )THEN
          READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
          WRITE(7,1002)NCARD
          WRITE(7,*)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
          IF( ISO > 0 ) GOTO 100
        ELSEIF( NPFORT >= 1 )THEN
          READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
          RAD=PI2*PFPH(NP,M)/TCP(M)
          CPFAM0(NP,M)=PFAM(NP,M)*COS(RAD)
          SPFAM0(NP,M)=PFAM(NP,M)*SIN(RAD)
          WRITE(7,1002)NCARD
          WRITE(7,*)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),CPFAM0(NP,M),SPFAM0(NP,M)
          IF( ISO > 0 ) GOTO 100
          READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M)
          RAD=PI2*PFPH(NP,M)/TCP(M)
          CPFAM1(NP,M)=PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
          SPFAM1(NP,M)=PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
          WRITE(7,1002)NCARD
          WRITE(7,*)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),CPFAM1(NP,M),SPFAM1(NP,M)
          CPFAM2(NP,M)=0.0
          SPFAM2(NP,M)=0.0
        ELSEIF( NPFORT == 2 )THEN
          READ(1,*,IOSTAT=ISO)NDUM,CDUM,PFAM(NP,M),PFPH(NP,M),PFX2(NP,M)
          RAD=PI2*PFPH(NP,M)/TCP(M)
          IF( PFX2(NP,M)>0.0 )THEN
            CPFAM2(NP,M)=PFAM(NP,M)*COS(RAD)-CPFAM0(NP,M)
            SPFAM2(NP,M)=PFAM(NP,M)*SIN(RAD)-SPFAM0(NP,M)
          ELSE
            CPFAM2(NP,M)=0.
            SPFAM2(NP,M)=0.
          ENDIF
         ENDIF
                                                                  
      ENDDO
    ENDDO
  ENDIF
                                                                  
  IF( NPBS > 0 )THEN
    !C18*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS                                                        
    !     ON SOUTH OPEN BOUNDARIES                                                                                          
    NCARD='18'
    CALL SEEK('C18')
    IF( NPFORT == 0 )THEN
      DO L=1,NPBS
        READ(1,*,IOSTAT=ISO)IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          IF( NPFORS == 0) EXIT
          RAD=PI2*PFPH(NPFORS,M)/TCP(M)
          AMP=G*PFAM(NPFORS,M)
          PCBS(L,M)=AMP*COS(RAD)
          PSBS(L,M)=AMP*SIN(RAD)
        ENDDO
      ENDDO
   ELSEIF( NPFORT == 1 )THEN
      DO L=1,NPBS
        READ(1,*,IOSTAT=ISO) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L),NPSERS1(L),TPCOORDS(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L),NPSERS1(L),TPCOORDS(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          PCBS(L,M)=CPFAM0(NPFORS,M)+TPCOORDS(L)*CPFAM1(NPFORS,M) + TPCOORDS(L)*TPCOORDS(L)*CPFAM2(NPFORS,M)
          PSBS(L,M)=SPFAM0(NPFORS,M)+TPCOORDS(L)*SPFAM1(NPFORS,M) + TPCOORDS(L)*TPCOORDS(L)*SPFAM2(NPFORS,M)
          TMPAMP=SQRT(PCBS(L,M)*PCBS(L,M)+PSBS(L,M)*PSBS(L,M))
          TMPPHS=ATAN2(PSBS(L,M),PCBS(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBS(L,M)=G*PCBS(L,M)
          PSBS(L,M)=G*PSBS(L,M)
        ENDDO
      ENDDO
   ELSEIF( NPFORT == 2 )THEN
      DO L=1,NPBS
        READ(1,*,IOSTAT=ISO) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L),NPSERS1(L),TPCOORDS(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBS(L),JPBS(L),ISPBS(L),NPFORS,NPSERS(L),NPSERS1(L),TPCOORDS(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          BOTTOM=PFX2(NPFORS,M)*(1.0-PFX2(NPFORS,M))
          TOP1=TPCOORDS(L)*PFX2(NPFORS,M)*(TPCOORDS(L)-PFX2(NPFORS,M))
          TOP2=TPCOORDS(L)*(1.0-TPCOORDS(L))
          IF( BOTTOM == 0.0 )THEN
            TOP1=TPCOORDS(L)
            TOP2=TPCOORDS(L)*TPCOORDS(L)
          ELSE
            TOP1=TOP1/BOTTOM
            TOP2=TOP2/BOTTOM
          ENDIF
          PCBS(L,M)=CPFAM0(NPFORS,M)+TOP1*CPFAM1(NPFORS,M)+TOP2*CPFAM2(NPFORS,M)
          PSBS(L,M)=SPFAM0(NPFORS,M)+TOP1*SPFAM1(NPFORS,M)+TOP2*SPFAM2(NPFORS,M)
          TMPAMP=SQRT(PCBS(L,M)*PCBS(L,M)+PSBS(L,M)*PSBS(L,M))
          TMPPHS=ATAN2(PSBS(L,M),PCBS(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBS(L,M)=G*PCBS(L,M)
          PSBS(L,M)=G*PSBS(L,M)
        ENDDO
      ENDDO
    ENDIF
   2068 FORMAT(I4,3X,A2,5X,E14.4,3E14.5,5X,2I5)                                                                           
   2069 FORMAT(I4,3X,A2,5X,2E14.4,5X,2I5)                                                                                 
  ENDIF
                                                                  
  IF( NPBW > 0 )THEN
    !C19*  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS                                                        
    !     ON WEST OPEN BOUNDARIES                                                                                           
    NCARD='19'
    CALL SEEK('C19')
    IF( NPFORT == 0 )THEN
      DO L=1,NPBW
        READ(1,*,IOSTAT=ISO)IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L)
        WRITE(7,1002)NCARD
        WRITE(7,*) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          IF( NPFORW == 0) EXIT
          RAD=PI2*PFPH(NPFORW,M)/TCP(M)
          AMP=G*PFAM(NPFORW,M)
          PCBW(L,M)=AMP*COS(RAD)
          PSBW(L,M)=AMP*SIN(RAD)
        ENDDO
      ENDDO
    ELSEIF( NPFORT == 1 )THEN
      DO L=1,NPBW
        READ(1,*,IOSTAT=ISO) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L),NPSERW1(L),TPCOORDW(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L),NPSERW1(L),TPCOORDW(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          PCBW(L,M)=CPFAM0(NPFORW,M)+TPCOORDW(L)*CPFAM1(NPFORW,M)+TPCOORDW(L)*TPCOORDW(L)*CPFAM2(NPFORW,M)
          PSBW(L,M)=SPFAM0(NPFORW,M)+TPCOORDW(L)*SPFAM1(NPFORW,M)+TPCOORDW(L)*TPCOORDW(L)*SPFAM2(NPFORW,M)
          TMPAMP=SQRT(PCBW(L,M)*PCBW(L,M)+PSBW(L,M)*PSBW(L,M))
          TMPPHS=0.0
          IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBW(L,M),PCBW(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBW(L,M)=G*PCBW(L,M)
          PSBW(L,M)=G*PSBW(L,M)
        ENDDO
      ENDDO
    ELSEIF( NPFORT == 2 )THEN
      DO L=1,NPBW
        READ(1,*,IOSTAT=ISO) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L),NPSERW1(L),TPCOORDW(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBW(L),JPBW(L),ISPBW(L),NPFORW,NPSERW(L),NPSERW1(L),TPCOORDW(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          BOTTOM=PFX2(NPFORW,M)*(1.0-PFX2(NPFORW,M))
          TOP1=TPCOORDW(L)*PFX2(NPFORW,M)*(TPCOORDW(L)-PFX2(NPFORW,M))
          TOP2=TPCOORDW(L)*(1.0-TPCOORDW(L))
          IF( BOTTOM == 0.0 )THEN
            TOP1=TPCOORDW(L)
            TOP2=TPCOORDW(L)*TPCOORDW(L)
          ELSE
            TOP1=TOP1/BOTTOM
            TOP2=TOP2/BOTTOM
          ENDIF
          PCBW(L,M)=CPFAM0(NPFORW,M)+TOP1*CPFAM1(NPFORW,M)+TOP2*CPFAM2(NPFORW,M)
          PSBW(L,M)=SPFAM0(NPFORW,M)+TOP1*SPFAM1(NPFORW,M)+TOP2*SPFAM2(NPFORW,M)
          TMPAMP=SQRT(PCBW(L,M)*PCBW(L,M)+PSBW(L,M)*PSBW(L,M))
          TMPPHS=0.0
          IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBW(L,M),PCBW(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBW(L,M)=G*PCBW(L,M)
          PSBW(L,M)=G*PSBW(L,M)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
                                                                  
  IF( NPBE > 0 )THEN
    !C20*  READ PERIODIC FORCING (TIDAL)ELEVATION BOUNDARY CONDTIONS                                                         
    !     ON EAST OPEN BOUNDARIES                                                                                           
    NCARD='20'
    CALL SEEK('C20')
    IF( NPFORT == 0 )THEN
      DO L=1,NPBE
        READ(1,*,IOSTAT=ISO)IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          IF( NPFORE == 0) EXIT
          RAD=PI2*PFPH(NPFORE,M)/TCP(M)
          AMP=G*PFAM(NPFORE,M)
          PCBE(L,M)=AMP*COS(RAD)
          PSBE(L,M)=AMP*SIN(RAD)
        ENDDO
      ENDDO
    ELSEIF( NPFORT == 1 )THEN
      DO L=1,NPBE
        READ(1,*,IOSTAT=ISO) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L),NPSERE1(L),TPCOORDE(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L),NPSERE1(L),TPCOORDE(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          PCBE(L,M)=CPFAM0(NPFORE,M)+TPCOORDE(L)*CPFAM1(NPFORE,M)+TPCOORDE(L)*TPCOORDE(L)*CPFAM2(NPFORE,M)
          PSBE(L,M)=SPFAM0(NPFORE,M)+TPCOORDE(L)*SPFAM1(NPFORE,M)+TPCOORDE(L)*TPCOORDE(L)*SPFAM2(NPFORE,M)
          TMPAMP=SQRT(PCBE(L,M)*PCBE(L,M)+PSBE(L,M)*PSBE(L,M))
          TMPPHS=0.0
          IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBE(L,M),PCBE(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBE(L,M)=G*PCBE(L,M)
          PSBE(L,M)=G*PSBE(L,M)
        ENDDO
      ENDDO
   ELSEIF( NPFORT == 2 )THEN
      DO L=1,NPBE
        READ(1,*,IOSTAT=ISO) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L),NPSERE1(L),TPCOORDE(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBE(L),JPBE(L),ISPBE(L),NPFORE,NPSERE(L),NPSERE1(L),TPCOORDE(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          BOTTOM=PFX2(NPFORE,M)*(1.0-PFX2(NPFORE,M))
          TOP1=TPCOORDE(L)*PFX2(NPFORE,M)*(TPCOORDE(L)-PFX2(NPFORE,M))
          TOP2=TPCOORDE(L)*(1.0-TPCOORDE(L))
          IF( BOTTOM == 0.0 )THEN
            TOP1=TPCOORDE(L)
            TOP2=TPCOORDE(L)*TPCOORDE(L)
          ELSE
            TOP1=TOP1/BOTTOM
            TOP2=TOP2/BOTTOM
          ENDIF
          PCBE(L,M)=CPFAM0(NPFORE,M)+TOP1*CPFAM1(NPFORE,M)+TOP2*CPFAM2(NPFORE,M)
          PSBE(L,M)=SPFAM0(NPFORE,M)+TOP1*SPFAM1(NPFORE,M)+TOP2*SPFAM2(NPFORE,M)
          TMPAMP=SQRT(PCBE(L,M)*PCBE(L,M)+PSBE(L,M)*PSBE(L,M))
          TMPPHS=0.0
          IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBE(L,M),PCBE(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBE(L,M)=G*PCBE(L,M)
          PSBE(L,M)=G*PSBE(L,M)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
                                                                  
  IF( NPBN > 0 )THEN
    !
    ! **  READ PERIODIC FORCING (TIDAL) ELEVATION BOUNDARY CONDTIONS                                                        
    !     ON NORTH OPEN BOUNDARIES                                                                                          
    NCARD='21'
    CALL SEEK('C21')
    IF( NPFORT == 0 )THEN
      DO L=1,NPBN
        READ(1,*,IOSTAT=ISO)IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          IF( NPFORN == 0) EXIT
          RAD=PI2*PFPH(NPFORN,M)/TCP(M)
          AMP=G*PFAM(NPFORN,M)
          PCBN(L,M)=AMP*COS(RAD)
          PSBN(L,M)=AMP*SIN(RAD)
        ENDDO
      ENDDO
    ELSEIF( NPFORT >= 1 )THEN
      DO L=1,NPBN
        READ(1,*,IOSTAT=ISO) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L),NPSERN1(L),TPCOORDN(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L),NPSERN1(L),TPCOORDN(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          PCBN(L,M)=CPFAM0(NPFORN,M)+TPCOORDN(L)*CPFAM1(NPFORN,M)+TPCOORDN(L)*TPCOORDN(L)*CPFAM2(NPFORN,M)
          PSBN(L,M)=SPFAM0(NPFORN,M)+TPCOORDN(L)*SPFAM1(NPFORN,M)+TPCOORDN(L)*TPCOORDN(L)*SPFAM2(NPFORN,M)
          TMPAMP=SQRT(PCBN(L,M)*PCBN(L,M)+PSBN(L,M)*PSBN(L,M))
          TMPPHS=0.0
          IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBN(L,M),PCBN(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBN(L,M)=G*PCBN(L,M)
          PSBN(L,M)=G*PSBN(L,M)
       ENDDO
      ENDDO
    ELSEIF( NPFORT == 2 )THEN
      DO L=1,NPBN
        READ(1,*,IOSTAT=ISO) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L),NPSERN1(L),TPCOORDN(L)

        WRITE(7,1002)NCARD
        WRITE(7,*) IPBN(L),JPBN(L),ISPBN(L),NPFORN,NPSERN(L),NPSERN1(L),TPCOORDN(L)
        IF( ISO > 0 ) GOTO 100
        DO M=1,MTIDE
          BOTTOM=PFX2(NPFORN,M)*(1.0-PFX2(NPFORN,M))
          TOP1=TPCOORDN(L)*PFX2(NPFORN,M)*(TPCOORDN(L)-PFX2(NPFORN,M))
          TOP2=TPCOORDN(L)*(1.0-TPCOORDN(L))
          IF( BOTTOM == 0.0 )THEN
            TOP1=TPCOORDN(L)
            TOP2=TPCOORDN(L)*TPCOORDN(L)
          ELSE
            TOP1=TOP1/BOTTOM
            TOP2=TOP2/BOTTOM
          ENDIF
          PCBN(L,M)=CPFAM0(NPFORN,M)+TOP1*CPFAM1(NPFORN,M)+TOP2*CPFAM2(NPFORN,M)
          PSBN(L,M)=SPFAM0(NPFORN,M)+TOP1*SPFAM1(NPFORN,M)+TOP2*SPFAM2(NPFORN,M)
          TMPAMP=SQRT(PCBN(L,M)*PCBN(L,M)+PSBN(L,M)*PSBN(L,M))
          TMPPHS=0.0
          IF( TMPAMP>0.0) TMPPHS=ATAN2(PSBN(L,M),PCBN(L,M))
          TMPPHS=TMPPHS*TCP(M)/PI2
          IF( TMPPHS<0.0)TMPPHS=TMPPHS+TCP(M)
          PCBN(L,M)=G*PCBN(L,M)
          PSBN(L,M)=G*PSBN(L,M)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
                                                                  
  IF( NPFORT >= 1 )THEN
    CLOSE(2)
  ENDIF
                                                                  
  !22*  READ NUM OF SEDIMENT AMD TOXICS AND NUM OF CONCENTRATION TIME SERIES                                              
  NCARD='22'
  CALL SEEK('C22')
  READ(1,*,IOSTAT=ISO) NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),NCSER(4),NTOXSER,NSEDSER,NSNDSER,ISSBAL
  
  ! *** REMOVE UNUSED SETTINGS TO ALLOW FOR SELECTIVE ALLOCATIONS
  IF ( ISTRAN(6) < 1)  NSED = 0   
  IF ( ISTRAN(7) < 1 ) NSND = 0
  IF ( NSED == 0 .AND. NSND == 0 ) ISTRAN(5) = 0
  IF ( ISTRAN(5) < 1 ) NTOX = 0
  IF( NTOX == 0 ) NTOXSER = 0
  IF( NSED == 0 ) NSEDSER = 0
  IF( NSND == 0 ) NSNDSER = 0
  
  WRITE(7,1002)NCARD
  WRITE(7,*) NTOX,NSED,NSND,NCSER(1),NCSER(2),NCSER(3),NCSER(4),NTOXSER,NSEDSER,NSNDSER,ISSBAL
  IF( ISO > 0 ) GOTO 100
  
  MTMP=4
  DO NS=1,NTOX
    MTMP=MTMP+1
    MSVTOX(NS)=MTMP
  ENDDO
  DO NS=1,NSED
    MTMP=MTMP+1
    MSVSED(NS)=MTMP
  ENDDO
  DO NS=1,NSND
    MTMP=MTMP+1
    MSVSND(NS)=MTMP
  ENDDO
  
  ! *** NCSER
  DO NS=1,NTOX
    M=MSVTOX(NS)
    NCSER(M)=NTOXSER
  ENDDO
  DO NS=1,NSED
    M=MSVSED(NS)
    NCSER(M)=NSEDSER
  ENDDO
  DO NS=1,NSND
    M=MSVSND(NS)
    NCSER(M)=NSNDSER
  ENDDO
  IF( ISTRAN(6) == 0 .AND. ISTRAN(7) == 0 )THEN
    ISSBAL=0  ! *** PMC SINGLE LINE
  ENDIF
                                                                  
  !C23*  READ VELOCITY, VOL SOUR/SINK, FLOW CONTROL, & WITHDRAW/RETURN DATA                                                
  NCARD='23'
  CALL SEEK('C23')
  READ(1,*,IOSTAT=ISO) NQSIJ,NQJPIJ,NQSER,NQCTL,NQCTLT,NHYDST,NQWR,NQWRSR,ISDIQ,NQCTLSER,NQCRULES
  
  WRITE(7,1002)NCARD
  WRITE(7,*) NQSIJ,NQJPIJ,NQSER,NQCTL,NQCTLT,NHYDST,NQWR,NQWRSR,ISDIQ,NQCTLSER,NQCRULES

  !IF( ISO > 0 ) GOTO 100
                                                                  
  IF( NQSIJ > 0 )THEN
    !C24*  READ VOLUME SOURCE/SINK LOCATIONS, MAGNITUDES, & VOL & CONC SERIES                                                
    NCARD='24'
    CALL SEEK('C24')
    DO L=1,NQSIJ
      READ(1,*,IOSTAT=ISO)IQS(L),JQS(L),QSSE,NQSMUL(L),NQSMF(L),NQSERQ(L),NCSERQ(L,1),NCSERQ(L,2),NCSERQ(L,3),NCSERQ(L,4),NTOXSRQ,NSEDSRQ,NSNDSRQ,QWIDTH(L),QFACTOR(L)

      WRITE(7,1002)NCARD
      WRITE(7,*)IQS(L),JQS(L),QSSE,NQSMUL(L),NQSMF(L),NQSERQ(L),NCSERQ(L,1),NCSERQ(L,2),NCSERQ(L,3),NCSERQ(L,4),NTOXSRQ,NSEDSRQ,NSNDSRQ,QWIDTH(L),QFACTOR(L)
      IF( ISO > 0 ) GOTO 100
      DO K=1,KC
        QSS(K,L) = QSSE*DZCK(K)
      ENDDO
      DO NS=1,NTOX
        M=MSVTOX(NS)
        NCSERQ(L,M)=NTOXSRQ
      ENDDO
      DO NS=1,NSED
        M=MSVSED(NS)
        NCSERQ(L,M)=NSEDSRQ
      ENDDO
      DO NS=1,NSND
        M=MSVSND(NS)
        NCSERQ(L,M)=NSNDSRQ
      ENDDO
    ENDDO
                                                                  
    !C25*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS                                                        
    !     SAL,TEM,DYE,SFL,TOX(1 TO NOTX)                                                                                    
    NCARD='25'
    CALL SEEK('C25')
    MMAX=4+NTOX
    DO L=1,NQSIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
      DO MS=1,MMAX
        DO K=1,KC
          CQS(K,L,MS)=CQSE(MS)
        ENDDO
      ENDDO
    ENDDO
                                                                  
    !C26*  READ TIME CONSTANT VOLUMETRIC SOURCE INFLOW CONCENTRATIONS                                                        
    !     SED(1 TO NSED),SND(1 TO NSND)                                                                                     
    NCARD='26'
    CALL SEEK('C26')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NQSIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
      DO MS=MMIN,MMAX
        DO K=1,KC
          CQS(K,L,MS)=CQSE(MS)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
                                                                  
  IF( NQJPIJ > 0 )THEN
    !C27*  READ JET/PLUME SOURCE LOCATIONS AND PARAMETERS                                                                    
    NCARD='27'
    CALL SEEK('C27')
    DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) IDUM,ICALJP(L),IQJP(L),JQJP(L),KQJP(L),NPORTJP(L),XJETL(L),YJETL(L),ZJET(L),PHJET(L),THJET(L),DJET(L),CFRD(L),DJPER(L)

      WRITE(7,1002)NCARD
      WRITE(7,*)IDUM,ICALJP(L),IQJP(L),JQJP(L),KQJP(L),NPORTJP(L),XJETL(L),YJETL(L),ZJET(L),PHJET(L),THJET(L),DJET(L),CFRD(L),DJPER(L)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C28*  READ JET/PLUME SOURCE LOCATIONS AND PARAMETERS                                                                    
    NCARD='28'
    CALL SEEK('C28')
    DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) IDUM,NJEL(L),NJPMX(L),ISENT(L),ISTJP(L),NUDJP(L),IOUTJP(L),NZPRJP(L),ISDJP(L),IUPCJP(L),JUPCJP(L),KUPCJP(L)

      WRITE(7,1002)NCARD
      WRITE(7,*) IDUM,NJEL(L),NJPMX(L),ISENT(L),ISTJP(L),NUDJP(L),IOUTJP(L),NZPRJP(L),ISDJP(L),IUPCJP(L),JUPCJP(L),KUPCJP(L)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C29*  READ ADDITIONAL JET/PLUME PARAMETERS                                                                              
    NCARD='29'
    CALL SEEK('C29')
    DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) IDUM,QQCJP(L),NQSERJP(L),NQWRSERJP(L),NCSERJP(L,1),NCSERJP(L,2),NCSERJP(L,3),NCSERJP(L,4),NTXSRJP,NSDSRJP,NSNSRJP

      WRITE(7,1002)NCARD
      WRITE(7,*) IDUM,QQCJP(L),NQSERJP(L),NQWRSERJP(L),NCSERJP(L,1),NCSERJP(L,2),NCSERJP(L,3),NCSERJP(L,4),NTXSRJP,NSDSRJP,NSNSRJP
      NUDJPC(L)=1
      IF( ISO > 0 ) GOTO 100
      DO NS=1,NTOX
        M=MSVTOX(NS)
        NCSERJP(L,M)=NTXSRJP
      ENDDO
      DO NS=1,NSED
        M=MSVSED(NS)
        NCSERJP(L,M)=NSDSRJP
      ENDDO
      DO NS=1,NSND
        M=MSVSND(NS)
        NCSERJP(L,M)=NSNSRJP
      ENDDO
      IF( ICALJP(L) == 2 )THEN
        QWRCJP(L)=QQCJP(L)
        QQCJP(L)=0.
      ELSE
        QWRCJP(L)=0.
      ENDIF
    ENDDO
    IF( NQJPIJ > 1 )THEN
      DO L=2,NQJPIJ
        NUDJP(L)=NUDJP(1)
      ENDDO
    ENDIF
                                                                  
    !C30*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT                                                        
    !     JET/PLUME SOURCES SAL,TEM,DYE,SFL,TOX(1 TO NOTX)                                                                  
    NCARD='30'
    CALL SEEK('C30')
    MMAX=4+NTOX
    DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
      IF( ICALJP(L) == 1 )THEN
        DO MS=1,MMAX
          CWRCJP(L,MS)=0.0
          DO K=1,KC
            CQCJP(K,L,MS)=CQSE(MS)
          ENDDO
        ENDDO
      ELSE
        DO MS=1,MMAX
          CWRCJP(L,MS)=CQSE(MS)
          DO K=1,KC
            CQCJP(K,L,MS)=0.0
          ENDDO
        ENDDO
      ENDIF
    ENDDO
                                                                  
    !C31*  READ TIME CONSTANT INFLOW CONCENTRATIONS FOR TIME CONSTANT                                                        
    !     JET/PLUME SOURCES SED(1 TO NSED),SND(1 TO NSND)                                                                   
    NCARD='31'
    CALL SEEK('C31')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NQJPIJ
      READ(1,*,IOSTAT=ISO) (CQSE(M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQSE(M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
      IF( ICALJP(L) == 1 )THEN
        DO MS=MMIN,MMAX
          CWRCJP(L,MS)=0.
          DO K=1,KC
            CQCJP(K,L,MS)=CQSE(MS)
          ENDDO
        ENDDO
      ELSE
        DO MS=MMIN,MMAX
          CWRCJP(L,MS)=CQSE(MS)
          DO K=1,KC
            CQCJP(K,L,MS)=0.
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
                                                                  
  IF( NQCTL > 0 )THEN
    !C32*  READ SURF ELEV OR PRESS DEPENDENT FLOW CONTROL STRUCTURE INFO                                                     
    NCARD='32'
    CALL SEEK('C32')
    DO L=1,NQCTL
      READ(1,*,IOSTAT=ISO)IQCTLU(L),JQCTLU(L),IQCTLD(L),JQCTLD(L),NQCTYP(L),NQCTLQ(L),NQCMUL(L),HQCTLU(L),HQCTLD(L),QCTLMU(L),QCTLGRP(L),BQCLCE(L),NQCMINS(L),HS_FACTOR(L),HS_NTIMES(L),HS_TRANSITION(L)

      WRITE(7,1002)NCARD
      WRITE(7,*)IQCTLU(L),JQCTLU(L),IQCTLD(L),JQCTLD(L),NQCTYP(L),NQCTLQ(L),NQCMUL(L),HQCTLU(L),HQCTLD(L),QCTLMU(L),QCTLGRP(L),BQCLCE(L),NQCMINS(L),HS_FACTOR(L),HS_NTIMES(L),HS_TRANSITION(L)
      IF( ISO > 0 ) GOTO 100
      DO K=1,KC
        QCTLTO(K,L)=0.
        QCTLT(K,L,1)=0.
        QCTLT(K,L,2)=0.
      ENDDO
      IF( NQCTYP(L) > 4 )THEN
        ! 2018-03-26, NTL: TODO: FOLLOWING CHECK SHOULD BE REMOVED TO ALLOW A LARGER NUMBER OF LOOKUP TABLES
        IF( NQCTLQ(L) > NHYDST )THEN
          PRINT*,' *** BAD HYDRAULIC STRUCTURE, L,NQCTLQ = ',L,NQCTLQ(L)
          STOP
        ENDIF
        IF (HS_FACTOR(L) <= 0 ) THEN
          PRINT*,' *** BAD HS_FACTOR: L,HS_FACTOR = ',L,HS_FACTOR(L)
          STOP
        ENDIF
      ENDIF

    ENDDO

    !C32A*  READ THE EQUATION PARAMETERS FOR EACH OF THE HYDRAULIC STRUCTURE EQUATIONS
    IF( NHYDST > 0 )THEN
      NCARD='32A'
      CALL SEEK('C32A')
      DO L=1,NHYDST
        ! *** COMPUTE FLOWS USING HYDRAULIC STRUCTURE EQUATIONS
        READ(1,*,IOSTAT=ISO)NS,NX,HS_REVERSE(L),HS_XSTYPE(L),HS_WIDTH(L),HS_HEIGHT(L),HS_LENGTH(L),HS_MANN(L),HS_ANGLE(L),HS_USELEV(L),HS_DSELEV(L),HS_COEFF(L,1),HS_COEFF(L,2),HS_COEFF(L,3),HS_COEFF(L,4)

        CALL HYDSTRUCT_CHECK(NX,L)  ! *** CHECK STRUCTURE DEFINITIONS

        ! *** USE THE SAME INVERT ELEVATION FOR US/DS FOR SLUICE GATES AND WEIRS
        IF( NX > 5 )THEN
          HS_DSELEV(L) = HS_USELEV(L)
        ENDIF
        
        WRITE(7,1002)NCARD
        WRITE(7,*)NS,NX,HS_REVERSE(L),HS_XSTYPE(L),HS_WIDTH(L),HS_HEIGHT(L),HS_LENGTH(L),HS_MANN(L),HS_ANGLE(L),HS_USELEV(L),HS_DSELEV(L),HS_COEFF(L,1),HS_COEFF(L,2),HS_COEFF(L,3),HS_COEFF(L,4)
        HS_ANGLE(L) = HS_ANGLE(L)*PI/180.0D0
      ENDDO
    ENDIF
  
    !C32B*  READ THE CONTROL INFO FOR ALL HYDRAULIC STRUCTURES
    IF( NQCTLSER > 0 .OR. NQCRULES > 0 )THEN
      NCARD='32B'
      CALL SEEK('C32B')
      DO L=1,NQCTL
        HSCTL(L).ITYPE = 0
        HSCTL(L).ID = 0
        !HSCTL(L).SUBID = 0
        HSCTL(L).IQCTL1 = IQCTLU(L)
        HSCTL(L).JQCTL1 = JQCTLU(L)
        HSCTL(L).IQCTL2 = IQCTLD(L)
        HSCTL(L).JQCTL2 = JQCTLD(L)
        READ(1,*,IOSTAT=ISO) HSCTL(L).ITYPE,HSCTL(L).ID,ITMPU,JTMPU,ITMPD,JTMPD, &
          HSCTL(L).CUR.STATE, HSCTL(L).CUR.HEIGHT, HSCTL(L).CUR.WIDTH, &
          HSCTL(L).CUR.SILL, HSCTL(L).CUR.UNITS, HSCTL(L).CUR.FLOW
        !HSCTL(L).CUR.ID, 
        IF (ITMPU >= IMN .AND. JTMPU >= JMN) THEN
          HSCTL(L).IQCTL1 = ITMPU
          HSCTL(L).JQCTL1 = JTMPU
        ENDIF
        IF (ITMPD >= IMN .AND. JTMPD >= JMN) THEN
          HSCTL(L).IQCTL2 = ITMPD
          HSCTL(L).JQCTL2 = JTMPD
        ENDIF
        
        WRITE(7,1002)NCARD
        WRITE(7,*) HSCTL(L).ITYPE,HSCTL(L).ID,HSCTL(L).IQCTL1,HSCTL(L).JQCTL1,HSCTL(L).IQCTL2,HSCTL(L).JQCTL2, &
          HSCTL(L).CUR.STATE, HSCTL(L).CUR.HEIGHT, HSCTL(L).CUR.WIDTH, &
          HSCTL(L).CUR.SILL, HSCTL(L).CUR.UNITS, HSCTL(L).CUR.FLOW  !HSCTL(L).CUR.ID, 
      ENDDO
    ENDIF  
  ENDIF
                                                                  
  IF( NQWR > 0 )THEN
    !C33*  READ FLOW WITHDRAWAL, HEAT OR MATERIAL ADDITION, FLOW RETURN DATA                                                 
    NCARD='33'
    CALL SEEK('C33')
    DO L=1,NQWR
      READ(1,*,IOSTAT=ISO)IQWRU(L),JQWRU(L),KQWRU(L),IQWRD(L),JQWRD(L),KQWRD(L),QWR(L),NQWRSERQ(L),NQWRMFU(L),NQWRMFD(L),BQWRMFU(L),BQWRMFD(L),ANGWRMFD(L)

      WRITE(7,1002)NCARD
      WRITE(7,*)IQWRU(L),JQWRU(L),KQWRU(L),IQWRD(L),JQWRD(L),KQWRD(L),QWR(L),NQWRSERQ(L),NQWRMFU(L),NQWRMFD(L),BQWRMFU(L),BQWRMFD(L),ANGWRMFD(L)
      IF( ISO > 0 ) GOTO 100
    ENDDO
    
    !C33A*  READ FLOW WITHDRAWAL/RETURN CONTROL
    NCARD='33A'
    CALL SEEK('C33A')
    DO L=1,NQWR
      WRCTL(L).ITYPE = 0
      WRCTL(L).ID = NQWRSERQ(L)
      !WRCTL(L).SUBID = 0
      WRCTL(L).IQCTL1 = IQWRU(L)
      WRCTL(L).JQCTL1 = JQWRU(L)
      WRCTL(L).IQCTL2 = IQWRD(L)
      WRCTL(L).JQCTL2 = JQWRD(L)
      READ(1,*,IOSTAT=ISO) WRCTL(L).ITYPE,WRCTL(L).ID,ITMPU,JTMPU,ITMPD,JTMPD, &
        WRCTL(L).CUR.STATE, WRCTL(L).CUR.FLOW
      IF (ITMPU >= IMN .AND. JTMPU >= JMN) THEN
        WRCTL(L).IQCTL1 = ITMPU
        WRCTL(L).JQCTL1 = JTMPU
      ENDIF
      IF (ITMPD >= IMN .AND. JTMPD >= JMN) THEN
        WRCTL(L).IQCTL2 = ITMPD
        WRCTL(L).JQCTL2 = JTMPD
      ENDIF
        
      WRITE(7,1002)NCARD
      WRITE(7,*) WRCTL(L).ITYPE,WRCTL(L).ID,WRCTL(L).IQCTL1,WRCTL(L).JQCTL1,WRCTL(L).IQCTL2,WRCTL(L).JQCTL2, &
        WRCTL(L).CUR.STATE, WRCTL(L).CUR.FLOW
    ENDDO
    
    !C34*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES                                                 
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='34'
    CALL SEEK('C34')
    MMAX=4+NTOX
    DO L=1,NQWR
      READ(1,*,IOSTAT=ISO) (CQWR(L,MS),MS=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQWR(L,MS),MS=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C35*  READ TIME CONSTANT WITHDRAWAL,ADD,RETURN, CONCENTRATION INCREASES                                                 
    !     SED(1 TO NSED),SND(1 TO NSND)                                                                                     
    NCARD='35'
    CALL SEEK('C35')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NQWR
      READ(1,*,IOSTAT=ISO) (CQWR(L,MS),MS=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CQWR(L,MS),MS=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
                                                                  
  IF( NSED > 0 .OR. NSND > 0 )THEN
    !C36*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATqION OPTIONS                                               
    NCARD='36'
    CALL SEEK('C36')
    READ(1,*,IOSTAT=ISO)ISEDINT,ISEDBINT,ISEDWC,ISMUD,ISNDWC,ISEDVW,ISNDVW,KB,ISDTXBUG

    WRITE(7,1002)NCARD
    WRITE(7,*)ISEDINT,ISEDBINT,ISEDWC,ISMUD,ISNDWC,ISEDVW,ISNDVW,KB,ISDTXBUG
    IF( ISO > 0 ) GOTO 100
                                                                  
    !C36A*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS                                              
    NCARD='36A'
    CALL SEEK('C36A')
    READ(1,*,IOSTAT=ISO)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,ISBKERO

    WRITE(7,1002)NCARD
    WRITE(7,*)ISBEDSTR,ISBSDIAM,ISBSDFUF,COEFTSBL,VISMUDST,ISBKERO
    IF( ISO > 0 ) GOTO 100
                                                                  
    !C36B*  SEDIMENT INITIALIZATION AND WATER COLUMN/BED REPRESENTATION OPTIONS                                              
    NCARD='36B'
    CALL SEEK('C36B')
    IF( NSED > 0 .OR. NSND > 0 )THEN
      READ(1,*,IOSTAT=ISO)ISEDAL,ISNDAL,IALTYP,IALSTUP,ISEDEFF,HBEDAL,COEHEFF,COEHEFF2

      WRITE(7,1002)NCARD
      WRITE(7,*)ISEDAL,ISNDAL,IALTYP,IALSTUP,HBEDAL,COEHEFF,COEHEFF2
      IF( ISO > 0 ) GOTO 100
    ENDIF
    ! *** FORCE BED ARMORING OPTION "INITIALIZATION AT STARTUP" TO BE OFF
    IF( ISRESTI > 0 ) IALSTUP = 0      
    
    !C37*  BED MECHANICAL PROPERTIES PARAMETER SET 1                                                                         
    NCARD='37'
    CALL SEEK('C37')
    IF( NSED>0 .OR. NSND > 0 )THEN
      READ(1,*,IOSTAT=ISO) SEDSTEP,SEDSTART,IBMECH,IMORPH,HBEDMAX,BEDPORC,SEDMDMX,SEDMDMN,SEDVDRD,SEDVDRM,SEDVRDT

      WRITE(7,1002)NCARD
      WRITE(7,*) SEDSTEP,SEDSTART,IBMECH,IMORPH,HBEDMAX,BEDPORC,SEDMDMX,SEDMDMN,SEDVDRD,SEDVDRM,SEDVRDT
      IF( ISO > 0 ) GOTO 100
      IF( SEDSTART <= TBEGIN ) SEDSTART=TBEGIN
    ELSE
      BEDPORC=0.4
      SEDVDRD=1.*(1-BEDPORC)
      SEDVDRM=SEDVDRD
      SEDVDRT=0
    ENDIF

    IF( IBMECH == 0 )THEN
      SNDVDRD=BEDPORC/(1.-BEDPORC)
      SEDVDRM=SEDVDRD
    END IF

    SNDVDRD=BEDPORC/(1.-BEDPORC)
    DO NS=1,NSED
      VDRDEPO(NS)=SEDVDRD
    ENDDO
    DO NS=1,NSND
      NX=NS+NSED
      VDRDEPO(NX)=SNDVDRD
    ENDDO
                                                                  
    !C38*  BED MECHANICAL PROPERTIES PARAMETER SET 2                                                                         
    NCARD='38'
    CALL SEEK('C38')
    READ(1,*,IOSTAT=ISO)IBMECHK,BMECH1,BMECH2,BMECH3,BMECH4,BMECH5,BMECH6
    WRITE(7,1002)NCARD
    WRITE(7,*)IBMECHK,BMECH1,BMECH2,BMECH3,BMECH4,BMECH5,BMECH6
    IF( ISO > 0 ) GOTO 100
  ENDIF
                                                                  
  IF( NSED > 0 )THEN
    !C39*  READ COHESIVE SEDIMENT PARAMETER SET 1 REPEAT DATA LINE NSED TIMES                                                
    NCARD='39'
    CALL SEEK('C39')
    HADJ=0.0
    IF( NSED > 0 )THEN
      DO NS=1,NSED
        READ(1,*,IOSTAT=ISO)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISEDSCOR(NS),ISPROBDEP(NS)

        WRITE(7,1002)NCARD
        WRITE(7,*)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISPROBDEP(NS)
        IF( ISO > 0 ) GOTO 100
        SEDDIA(NS)=0.
        HADJ=SEDN(1)
      ENDDO
      IF( HADJ<HWET)HADJ=HWET  ! *** PMC-PROVIDE MORE CONTROL FOR !MORPH CHANGE LIMITS
    ENDIF
                                                                  
    !C40*  READ COHESIVE SEDIMENT PARAMETER SET 2 REPEAT DATA LINE NSED TIMES                                                
    NCARD='40'
    CALL SEEK('C40')
    DO NS=1,NSED
      READ(1,*,IOSTAT=ISO)IWRSP(NS),IWRSPB(NS),WRSPO(NS),TAUR(NS),TAUN(NS),TEXP(NS),VDRRSPO(NS),COSEDHID(NS)
      
      WRITE(7,1002)NCARD
      WRITE(7,*)IWRSP(NS),IWRSPB(NS),WRSPO(NS),TAUR(NS),TAUN(NS),TEXP(NS),VDRRSPO(NS),COSEDHID(NS)
      IF( ISO > 0 ) GOTO 100
                                                                  
      IF( NS == 1 .AND. IWRSP(NS) == 999 )THEN
        WRITE(*,'(A)')'READING TAU_CRIT_COH.INP'
        OPEN(1001,FILE='tau_crit_coh.inp',STATUS='OLD')
        DO L = 2, 4393
          READ(1001,*,IOSTAT=ISO) (TAUCRCOH(L,K),K=1,10)
        ENDDO
        CLOSE(1001)
      ENDIF
      IF( ISO > 0 ) GOTO 100
      ISNDEQ(NS)=0
      
    ENDDO
  ENDIF

  TAUCMIN=1000.
  ISBDLDBC=0
  SSG=2.65         ! *** DEFAULT GRAIN DENSITY USING QUARTZ
  IF( NSND > 0 )THEN
    !C41*  READ NONCOHESIVE SEDIMENT PARAMETER SET 1 REPEAT DATA LINE NSND TIMES                                             
    NCARD='41'
    CALL SEEK('C41')
    DO NX=1,NSND
      NS=NX+NSED
      READ(1,*,IOSTAT=ISO)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),SEDDIA(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISEDSCOR(NS)

      ! *** IF SETTLING VELOCITY IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA
      IF( WSEDO(NS) < 0.0 )THEN
        WSEDO(NS) = SETSTVEL(SEDDIA(NS),SSG(NS))
      ENDIF

      WRITE(7,1002)NCARD
      WRITE(7,*)SEDO(NS),SEDBO(NS),SDEN(NS),SSG(NS),SEDDIA(NS),WSEDO(NS),SEDN(NS),SEXP(NS),TAUD(NS),ISEDSCOR(NS)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C42*  READ NONCOHESIVE SEDIMENT PARAMETER SET 2 REPEAT DATA LINE NSND TIMES                                             
    NCARD='42'
    CALL SEEK('C42')
    DO NX=1,NSND
      NS=NX+NSED
      READ(1,*,IOSTAT=ISO)ISNDEQ(NS),ISBDLD(NS),TAUR(NS),TAUN(NS),TCSHIELDS(NS),ISLTAUC(NS),IBLTAUC(NS),IROUSE(NX),ISNDM1(NX),ISNDM2(NX),RSNDM(NX)
                                                                  
      ! IF TAUR(NS) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA                                                              
      !       TAUR:     CRITICAL SHIELDS STRESS IN (M/S)**2   (ISNDEQ=2)                                                      
      !       TAUN:     EQUAL TO TAUR FOR NONCHOESIVE SED TRANS  (ISNDEQ=2)                                                   
      !       TEXP:     CRITICAL SHIELDS PARAMETER  (ISNDEQ=2)                                                                
      DSTR=0.0
      USTR=0.0

      ! *** IF TAUR(NS) IS NEGATIVE, COMPUTE USING VAN RIJN'S FORMULA                                                              
      IF( TAUR(NS) < 0.0 )THEN
        CALL SETSHLD(TAUR(NS),TCSHIELDS(NS),SEDDIA(NS),SSG(NS),DSTR,USTR)
        TAUN(NS)=TAUR(NS)
      ENDIF
      TAUCMIN=MIN(TAUCMIN,TAUR(NS))
                                                                  
      WRITE(7,1002)NCARD
      WRITE(7,*)ISNDEQ(NS),TAUR(NS),TAUN(NS),TCSHIELDS(NS),SEDDIA(NS),SSG(NS),DSTR,USTR
      IF( ISO > 0 ) GOTO 100
      IWRSP(NS)=0
      WRSPO(NS)=0
    ENDDO
                                                                  
    !C42A*  READ NONCOHESIVE SEDIMENT BED LOAD PARAMETERS                                                                    
    NCARD='42A'
    CALL SEEK('C42A')
    DO NS=1,NSND
      READ(1,*,IOSTAT=ISO)ISBDLDBC,SBDLDA(NS),SBDLDB(NS),SBDLDG1(NS),SBDLDG2(NS),SBDLDG3(NS),SBDLDG4(NS),SBDLDP(NS),ISBLFUC,BLBSNT

      WRITE(7,1002)NCARD
      WRITE(7,*)ISBDLDBC,SBDLDA(NS),SBDLDB(NS),SBDLDG1(NS),SBDLDG2(NS),SBDLDG3(NS),SBDLDG4(NS),SBDLDP(NS),ISBLFUC,BLBSNT
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
                                                                  
  IF( NTOX > 0 )THEN
    !C43A*  READ TOXIC CONTAMINANT INITIAL CONDITIONS         
    NCARD='43A'
    CALL SEEK('C43A')
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO)NDUM,ITXINT(NT),ITXBDUT(NT),TOXINTW(NT),TOXINTB(NT)
      WRITE(7,1002)NCARD
      WRITE(7,*)NDUM,ITXINT(NT),ITXBDUT(NT),TOXINTW(NT),TOXINTB(NT)
      IF( ISO > 0 ) GOTO 100
    ENDDO

    !C43B*  READ TOXIC KINETIC OPTION FLAGS
    NCARD='43B'
    CALL SEEK('C43B')
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO) NDUM,ITOXKIN(1,NT),ITOXKIN(2,NT),ITOXKIN(3,NT),ITOXKIN(4,NT),ITOXKIN(5,NT)
      WRITE(7,1002)NCARD
      WRITE(7,*) NDUM,ITOXKIN(1,NT),ITOXKIN(2,NT),ITOXKIN(3,NT),ITOXKIN(4,NT),ITOXKIN(5,NT)
      IF( ISO > 0 ) GOTO 100
    ENDDO

    !C43C*  READ TOXIC TIMING AND VOLATILIZATION SWITCHES
    NCARD='43C'
    CALL SEEK('C43C')
    READ(1,*,IOSTAT=ISO) TOXSTEPW, TOXSTEPB, TOX_VEL_MAX, TOX_DEP_MAX,ITOXTEMP,TOXTEMP
    WRITE(7,1002)NCARD
    WRITE(7,*) TOXSTEPW, TOXSTEPB, TOX_VEL_MAX, TOX_DEP_MAX,ITOXTEMP,TOXTEMP
    IF( ISO > 0 ) GOTO 100
    IF( TOXSTEPW < SEDSTEP ) TOXSTEPW = SEDSTEP    ! *** TOXIC KINETICS CANNOT OPERATE ON A FINER TIMESCALE THAN SEDIMENTS
    IF( TOXSTEPB < SEDSTEP ) TOXSTEPB = SEDSTEP    ! *** TOXIC BED MIXING AND DIFFUSION CANNOT OPERATE ON A FINER TIMESCALE THAN SEDIMENTS
    IF( ISTRAN(5) > 0 .AND. ISTRAN(2) == 0 .AND. ( ITOXTEMP < 0 .OR. ITOXTEMP-1 > NCSER(2) ) ) STOP 'ITOXTEMP IS OUT OF RANGE'
    
    !C43D*  READ TOXIC BULK DECAY AND BIODEGRADATION PARAMETERS
    NCARD='43D'
    CALL SEEK('C43D')
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO) NDUM,TOX_BLK_KW(NT), TOX_BLK_KB(NT), TOX_BLK_MXD(NT), TOX_BIO_KW(NT), TOX_BIO_KB(NT), TOX_BIO_MXD(NT), &
                                TOX_BIO_Q10W(NT), TOX_BIO_Q10B(NT), TOX_BIO_TW(NT), TOX_BIO_TB(NT)
      WRITE(7,1002)NCARD
      WRITE(7,*)  NDUM,TOX_BLK_KW(NT), TOX_BLK_KB(NT), TOX_BIO_MXD(NT), TOX_BIO_KW(NT), TOX_BIO_KB(NT), TOX_BIO_MXD(NT), &
                       TOX_BIO_Q10W(NT), TOX_BIO_Q10B(NT), TOX_BIO_TW(NT), TOX_BIO_TB(NT)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                          
    !C43D*  READ TOXIC VOLATILIZATION PARAMETERS
    NCARD='43E'
    CALL SEEK('C43E')
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO) NDUM,TOX_MW(NT), TOX_HE(NT), TOX_KV_TCOEFF(NT), TOX_ATM(NT), TOX_ADJ(NT)  
      WRITE(7,1002)NCARD
      WRITE(7,*) NDUM,TOX_MW(NT), TOX_HE(NT), TOX_KV_TCOEFF(NT), TOX_ATM(NT), TOX_ADJ(NT)
      IF( ISO > 0 ) GOTO 100
    ENDDO

    ! *** 
    
    !C44*  READ TOXIC CONTAMINANT PARAMETERS: SORBTION                                                                                    
    NCARD='44'
    CALL SEEK('C44')
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO)NDUM,ISTOC(NT),DIFTOX(NT),DIFTOXS(NT),PDIFTOX(NT),DPDIFTOX(NT)

      WRITE(7,1002)NCARD
      WRITE(7,*)NDUM,ISTOC(NT),DIFTOX(NT),DIFTOXS(NT),PDIFTOX(NT),DPDIFTOX(NT)
      IF( ISO > 0 ) GOTO 100
      ISPMXZ(NT)=0
      IF( PDIFTOX(NT) < 0.0 ) ISPMXZ(NT)=1
      ISDIFBW(NT)=0  
      IF( DIFTOXS(NT)<0.0 )THEN
        DIFTOXS(NT)=ABS(DIFTOXS(NT))
        ISDIFBW(NT)=1
      ENDIF
    ENDDO
                                                                  
    !C45*  READ TOXIC CONTAMINANT-SEDIMENT INTERACTION PARAMETERS                                                            
    NCARD='45'
    CALL SEEK('C45')
    DO NT=1,NTOX
      IF( NSED > 0 )THEN
        DO NS=1,NSED
          READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
          WRITE(7,1002)NCARD
          WRITE(7,*)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
          IF( ISO > 0 ) GOTO 100
        ENDDO
      ENDIF
      IF( NSND > 0 )THEN
        DO NX=1,NSND
          NS=NX+NSED
          READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
          WRITE(7,1002)NCARD
          WRITE(7,*)NDUM1,NDUM2,ITXPARW(NS,NT),TOXPARW(NS,NT),CONPARW(NS,NT),ITXPARB(NS,NT),TOXPARB(NS,NT),CONPARB(NS,NT)
          IF( ISO > 0 ) GOTO 100
        ENDDO
      ENDIF
    ENDDO
                                                                  
    !C45A*  READ TOXIC CONTAMINANT ORGANIC CARBON PARAMETERS                                                                 
    NCARD='45A'
    CALL SEEK('C45A')
    IF( NTOX > 0 )THEN
      READ(1,*,IOSTAT=ISO)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,STDOCWC,STPOCWC,STDOCBC,STPOCBC
      IF( ISTRAN(5) == 0 )THEN
        ISTDOCW=0
        ISTPOCW=0
        ISTDOCB=0
        ISTPOCB=0
      ENDIF
      WRITE(7,1002)NCARD
      WRITE(7,*)ISTDOCW,ISTPOCW,ISTDOCB,ISTPOCB,STDOCWC,STPOCWC,STDOCBC,STPOCBC
      IF( ISO > 0 ) GOTO 100
    ENDIF
                                                                  
    !C45B* READ TOXIC CONTAMINANT-ORGANIC CARBON INTERACTION PARAMETERS                                                      
    NCARD='45B'
    CALL SEEK('C45B')
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARWC(1,NT),TOXPARWC(1,NT),CONPARWC(1,NT),ITXPARBC(1,NT),TOXPARBC(1,NT),CONPARBC(1,NT)
      WRITE(7,1002)NCARD
      WRITE(7,*)NDUM1,NDUM2,ITXPARWC(1,NT),TOXPARWC(1,NT),CONPARWC(1,NT),ITXPARBC(1,NT),TOXPARBC(1,NT),CONPARBC(1,NT)
      IF( ISO > 0 ) GOTO 100
      READ(1,*,IOSTAT=ISO)NDUM1,NDUM2,ITXPARWC(2,NT),TOXPARWC(2,NT),CONPARWC(2,NT),ITXPARBC(2,NT),TOXPARBC(2,NT),CONPARBC(2,NT)
      WRITE(7,1002)NCARD
      WRITE(7,*)NDUM1,NDUM2,ITXPARWC(2,NT),TOXPARWC(2,NT),CONPARWC(2,NT),ITXPARBC(2,NT),TOXPARBC(2,NT),CONPARBC(2,NT)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C45C* READ TOXIC CONTAMINANT-ORGANIC CARBON WATER COLUMN POC FRACTIONS                                                  
    NCARD='45C'
    CALL SEEK('C45C')
    WRITE(7,1002)NCARD
    NTMP=NSED+NSND
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO)NDUM1,(FPOCWST(NS,NT),NS=1,NTMP)
      IF( ISO > 0 ) GOTO 100
      WRITE(7,*)NDUM1,(FPOCWST(NS,NT),NS=1,NTMP)
    ENDDO
                                                                  
    ! RESET INORGANIC SEDIMENT PARTITION COEFFICIENTS BASED ON                                                              
    ! FRACTION OF POC ASSIGNED TO EACH SEDIMENT CLASS IN WATER COLUMN                                                       
    DO NT=1,NTOX
      IF( ISTOC(NT) >= 2 )THEN
        IF( NSED > 0 )THEN
          DO NS=1,NSED
            ITXPARW(NS,NT)=0
            TOXPARW(NS,NT)=TOXPARWC(2,NT)
            CONPARW(NS,NT)=0.
          ENDDO
        ENDIF
        IF( NSND > 0 )THEN
          DO NX=1,NSND
            NS=NSED+NX
            ITXPARW(NS,NT)=0
            TOXPARW(NS,NT)=TOXPARWC(2,NT)
            CONPARW(NS,NT)=0.
          ENDDO
        ENDIF
      ENDIF
    ENDDO
                                                                  
    !C45D* READ TOXIC CONTAMINANT-ORGANIC CARBON SED BED POC FRACTIONS                                                       
    NCARD='45D'
    CALL SEEK('C45D')
    WRITE(7,1002)NCARD
    NTMP=NSED+NSND
    DO NT=1,NTOX
      READ(1,*,IOSTAT=ISO)NDUM1,(FPOCBST(NS,NT),NS=1,NTMP)
      IF( ISO > 0 ) GOTO 100
      WRITE(7,*)NDUM1,(FPOCBST(NS,NT),NS=1,NTMP)
    ENDDO
                                                                  
    ! RESET INORGANIC SEDIMENT PARTITION COEFFICIENTS BASED ON                                                              
    ! FRACTION OF POC ASSIGNED TO EACH SEDIMENT CLASS IN SEDIMENT BED                                                       
    DO NT=1,NTOX
      IF( ISTOC(NT) >= 2 )THEN
        IF( NSED > 0 )THEN
          DO NS=1,NSED
            ITXPARB(NS,NT)=0
            TOXPARB(NS,NT)=TOXPARBC(2,NT)
            CONPARB(NS,NT)=0
          ENDDO
        ENDIF
        IF( NSND > 0 )THEN
          DO NX=1,NSND
            NS=NSED+NX
            ITXPARB(NS,NT)=0
            TOXPARB(NS,NT)=TOXPARBC(2,NT)
            CONPARB(NS,NT)=0
          ENDDO
        ENDIF
      ENDIF
    ENDDO

  ENDIF  ! *** END OF NTOX>0 BLOCK
                                                                  
  !C46*  READ BUOYANCY, TEMPERATURE, DYE DATA AND CONCENTRATION BC DATA                                                    
  NCARD='46'
  CALL SEEK('C46')
  READ(1,*,IOSTAT=ISO)BSC,TEMO,HEQT,ITMP,KBH,RKDYE,NCBS,NCBW,NCBE,NCBN

  WRITE(7,1002)NCARD
  WRITE(7,*)BSC,TEMO,HEQT,KBH,RKDYE,NCBS,NCBW,NCBE,NCBN
  IF( ISO > 0 ) GOTO 100
  IF( BSC == 2. )THEN
    BSC=1.
    IBSC=1
  ELSE
    IBSC=0
  ENDIF
  ! *** DISABLE BOUYANCY IF ALL CONSTITUENTS IMPACTING DENSITY ARE OFF (7.2)
  IF( BSC > 0. .AND. (ISTRAN(1) < 1 .AND. ISTRAN(2) < 1 .AND.  .NOT. (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0)) )THEN
    BSC=0.
    IBSC=0
  ENDIF
  IF( TEMO < 0.0 )THEN
    TEMO=ABS(TEMO)
    INITTEMP=1
  ELSE
    INITTEMP=0
  ENDIF

  IF( ISICE > 0 )THEN
    !C46A*   READ ICE EFFECTS                                                                                                
    NCARD='46A'
    CALL SEEK('C46A') 
    READ(1,*,IOSTAT=ISO)ISICE,NISER,TEMPICE,CDICE,ICETHMX,RICETHK0
    WRITE(7,1002)NCARD
    WRITE(7,*)ISICE,NISER,TEMPICE,CDICE,ICETHMX,RICETHK0
    IF( ISO > 0 ) GOTO 100
  ENDIF
  
  RHOI = 917.
  RISEVEL = 0.01
  IF( ISICE > 2 )THEN
    !C46B*   READ ICE MODULE PARAMETERS       
    NCARD='46B'
    CALL SEEK('C46B')
    READ(1,*,IOSTAT=ISO) HWI,ICEK,ALBEDOI,BETAI,GAMMAI,MINICETHICK,RHOI,ISRHEVAP,RISEVEL,MELTFACTOR,AFWI,BFWI,CFWI
    WRITE(7,1002) NCARD
    WRITE(7,*) HWI,ICEK,ALBEDOI,BETAI,GAMMAI,MINICETHICK,ICETHMX,RHOI,ISRHEVAP,RISEVEL,MELTFACTOR,AFWI,BFWI,CFWI
    IF( ISO > 0 ) GOTO 100
    
    IF( CFWI <= 0. )THEN
      AFWI = 9.2    ! *** WIND FUNCTION A, Edinger, et. al. (1974)
      BFWI = 0.46   ! *** WIND FUNCTION B, Edinger, et. al. (1974)
      CFWI = 2.0    ! *** WIND FUNCTION C, Edinger, et. al. (1974)
    ELSE
      ! *** Units Conversion = 1 mmhg = 1.33322 millibars
      AFWI = AFWI*1.33322
      BFWI = BFWI*1.33322
    ENDIF

    ! *** FRAZIL ICE TRANSPORT DEFAULTS
    ISCDCA(10) = 0  ! *** UPWIND DIFFERENCE (3TL ONLY)
    ISFCT(10)  = 0 
  ENDIF
  IF( ISICE == 2) NISER = 1
  
  IF( NASER > 0 )THEN
    ! ** READING TWO NEW CARDS FOR EFDC_073 INSTAED OF ASER.INP FOR 072
    DS_LAT=0.0
    DS_LONG=0.0
    COMPUTESOLRAD=.FALSE.
    
    NCARD='46C'
    CALL SEEK('C46C')
    READ(1,*,IOSTAT=ISO) DS_LONG,DS_LAT,COMPUTESOLRAD,USESHADE,IEVAP,WINDFA,WINDFB,WINDFC
    WRITE(7,1002) NCARD
    WRITE(7,*) DS_LONG,DS_LAT,COMPUTESOLRAD,USESHADE,IEVAP,WINDFA,WINDFB,WINDFC
    IF( ISO > 0 ) GOTO 100
    
    !IF( IEVAP > 2 .AND. ISTRAN(2) > 0 .AND. ISTOPT(2) == 4 )THEN
    !  ! *** USE EVAPORATION WIND FUNCTION COEFFICIENTS FOR EQUILIBRIUM TEMPERATURE CALCULATIONS
    !  ! *** Units Conversion = 1 mmhg = 1.33322 millibars
    !  AFW = WINDFA*1.33322
    !  IF( WINDFC <= 0.0 )THEN
    !    ! *** WIND SQUARED TERM NOT USED
    !    BFW = WINDFB*1.33322
    !    CFW = 1.0
    !  ELSE
    !    ! *** WIND SQUARED TERM NOT USED
    !    BFW = WINDFC*1.33322
    !    CFW = 2.0
    !  ENDIF
    !ELSE
      AFW = 9.2        ! *** Edinger, et. al. (1974)  W/m2/mmHg
      BFW = 0.46       ! *** Edinger, et. al. (1974)  W/m2/mmHg
      CFW = 2.0        ! *** Edinger, et. al. (1974)  W/m2/mmHg
    !ENDIF
    
    ! *** CONVERT WIND FACTOR COEFFICIENTS FROM W/M2/MILLIBAR TO M/S/MILLIBAR
    ! *** Latent Heat of Evaporation = 2259 KJ/KG
    ! *** Density of Water           = 1000 kg/m3
    WINDFA = WINDFA*4.42674E-10
    WINDFB = WINDFB*4.42674E-10
    WINDFC = WINDFC*4.42674E-10
    
    NCARD='46D'
    CALL SEEK('C46D')
    READ(1,*,IOSTAT=ISO) IASWRAD,REVC,RCHC,ISVHEAT,SWRATNF,SWRATNS,FSWRATF,DABEDT,TBEDIT,HTBED1,HTBED2,WQKEB(1),WQKETSS
    WRITE(7,1002) NCARD
    WRITE(7,*) IASWRAD,REVC,RCHC,ISVHEAT,SWRATNF,SWRATNS,FSWRATF,DABEDT,TBEDIT,HTBED1,HTBED2,WQKEB(1),WQKETSS
    IF( ISO > 0 ) GOTO 100
    
  ENDIF
  
  ! ****************************************************************************C                                         
  IF( NCBS > 0 )THEN

    !C47*  READ LOCATIONS OF CONC BC'S ON SOUTH BOUNDARIES                                                                   
    NCARD='47'
    CALL SEEK('C47')
    DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) ICBS(L),JCBS(L),NTSCRS(L),NCSERS(L,1),NCSERS(L,2),NCSERS(L,3),NCSERS(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      WRITE(7,1002)NCARD
      WRITE(7,*) ICBS(L),JCBS(L),NTSCRS(L),NCSERS(L,1),NCSERS(L,2),NCSERS(L,3),NCSERS(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      IF( ISO > 0 ) GOTO 100
      DO NS=1,NTOX
        M=MSVTOX(NS)
        NCSERS(L,M)=NTOXSRC
      ENDDO
      DO NS=1,NSED
        M=MSVSED(NS)
        NCSERS(L,M)=NSEDSRC
      ENDDO
      DO NS=1,NSND
        M=MSVSND(NS)
        NCSERS(L,M)=NSNDSRC
      ENDDO
    ENDDO
                                                                  
    !C48*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES                                                       
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='48'
    CALL SEEK('C48')
    MMAX=4+NTOX
    DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,1,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C49*  READ CONSTANT BOTTOM CONCENTRATION ON SOUTH CONC BOUNDARIES                                                       
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='49'
    CALL SEEK('C49')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,1,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C50*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES                                                      
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='50'
    CALL SEEK('C50')
    MMAX=4+NTOX
    DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,2,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C51*  READ CONSTANT SURFACE CONCENTRATION ON SOUTH CONC BOUNDARIES                                                      
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='51'
    CALL SEEK('C51')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBS
      READ(1,*,IOSTAT=ISO) (CBS(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBS(L,2,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
                                                                  
  IF( NCBW > 0 )THEN
    !C52*  READ LOCATIONS OF CONC BC'S ON WEST BOUNDARIES                                                                    
    NCARD='52'
    CALL SEEK('C52')
    DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) ICBW(L),JCBW(L),NTSCRW(L),NCSERW(L,1),NCSERW(L,2),NCSERW(L,3),NCSERW(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      WRITE(7,1002)NCARD
      WRITE(7,*) ICBW(L),JCBW(L),NTSCRW(L),NCSERW(L,1),NCSERW(L,2),NCSERW(L,3),NCSERW(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      IF( ISO > 0 ) GOTO 100
      DO NS=1,NTOX
        M=MSVTOX(NS)
        NCSERW(L,M)=NTOXSRC
      ENDDO
      DO NS=1,NSED
        M=MSVSED(NS)
        NCSERW(L,M)=NSEDSRC
      ENDDO
      DO NS=1,NSND
        M=MSVSND(NS)
        NCSERW(L,M)=NSNDSRC
      ENDDO
    ENDDO
                                                                  
    !C53*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES                                                        
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='53'
    CALL SEEK('C53')
    MMAX=4+NTOX
    DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,1,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C54*  READ CONSTANT BOTTOM CONCENTRATION ON WEST CONC BOUNDARIES                                                        
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='54'
    CALL SEEK('C54')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,1,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C55*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES                                                       
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='55'
    CALL SEEK('C55')
    MMAX=4+NTOX
    DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,2,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C56*  READ CONSTANT SURFACE CONCENTRATION ON WEST CONC BOUNDARIES                                                       
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='56'
    CALL SEEK('C56')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBW
      READ(1,*,IOSTAT=ISO) (CBW(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBW(L,2,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
                                                                  
  IF( NCBE > 0 )THEN
    !C57*  READ LOCATIONS OF CONC BC'S ON EAST BOUNDARIES                                                                    
    NCARD='57'
    CALL SEEK('C57')
    DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) ICBE(L),JCBE(L),NTSCRE(L),NCSERE(L,1),NCSERE(L,2),NCSERE(L,3),NCSERE(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      WRITE(7,1002)NCARD
      WRITE(7,*) ICBE(L),JCBE(L),NTSCRE(L),NCSERE(L,1),NCSERE(L,2),NCSERE(L,3),NCSERE(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      IF( ISO > 0 ) GOTO 100
      DO NS=1,NTOX
        M=MSVTOX(NS)
        NCSERE(L,M)=NTOXSRC
      ENDDO
      DO NS=1,NSED
        M=MSVSED(NS)
        NCSERE(L,M)=NSEDSRC
      ENDDO
      DO NS=1,NSND
        M=MSVSND(NS)
        NCSERE(L,M)=NSNDSRC
      ENDDO
    ENDDO
                                                                  
    !C58*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES                                                        
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='58'
    CALL SEEK('C58')
    MMAX=4+NTOX
    DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,1,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C59*  READ CONSTANT BOTTOM CONCENTRATION ON EAST CONC BOUNDARIES                                                        
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='59'
    CALL SEEK('C59')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,1,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C60*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES                                                       
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='60'
    CALL SEEK('C60')
    MMAX=4+NTOX
    DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,2,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C61*  READ CONSTANT SURFACE CONCENTRATION ON EAST CONC BOUNDARIES                                                       
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='61'
    CALL SEEK('C61')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBE
      READ(1,*,IOSTAT=ISO) (CBE(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBE(L,2,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF    ! *** NCBE>0
                                                                  
  IF( NCBN > 0 )THEN
    !C62*  READ LOCATIONS OF CONC BC'S ON NORTH BOUNDARIES                                                                   
    NCARD='62'
    CALL SEEK('C62')
    DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) ICBN(L),JCBN(L),NTSCRN(L),NCSERN(L,1),NCSERN(L,2),NCSERN(L,3),NCSERN(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      WRITE(7,1002)NCARD
      WRITE(7,*) ICBN(L),JCBN(L),NTSCRN(L),NCSERN(L,1),NCSERN(L,2),NCSERN(L,3),NCSERN(L,4),NTOXSRC,NSEDSRC,NSNDSRC
      IF( ISO > 0 ) GOTO 100
      DO NS=1,NTOX
        M=MSVTOX(NS)
        NCSERN(L,M)=NTOXSRC
      ENDDO
      DO NS=1,NSED
        M=MSVSED(NS)
        NCSERN(L,M)=NSEDSRC
      ENDDO
      DO NS=1,NSND
        M=MSVSND(NS)
        NCSERN(L,M)=NSNDSRC
      ENDDO
    ENDDO
                                                                  
    !C63*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES                                                       
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='63'
    CALL SEEK('C63')
    MMAX=4+NTOX
    DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,1,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,1,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C64*  READ CONSTANT BOTTOM CONCENTRATION ON NORTH CONC BOUNDARIES                                                       
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='64'
    CALL SEEK('C64')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,1,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,1,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C65*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES                                                      
    !     SAL,TEM,DYE,SFL,TOX(1 TO NTOX)                                                                                    
    NCARD='65'
    CALL SEEK('C65')
    MMAX=4+NTOX
    DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,2,M),M=1,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,2,M),M=1,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                  
    !C66*  READ CONSTANT SURFACE CONCENTRATION ON NORTH CONC BOUNDARIES                                                      
    !     SED(1 TO NSED),SND(1,NSND)                                                                                        
    NCARD='66'
    CALL SEEK('C66')
    MMIN=MMAX+1
    MMAX=MMAX+NSED+NSND
    DO L=1,NCBN
      READ(1,*,IOSTAT=ISO) (CBN(L,2,M),M=MMIN,MMAX)
      WRITE(7,1002)NCARD
      WRITE(7,*) (CBN(L,2,M),M=MMIN,MMAX)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF    ! *** NCBN>0
                                                                  
  !C66A*  READ CONCENTRATION DATA ASSIMILATION PARAMETERS                                                                  
  NCARD='66A'
  CALL SEEK('C66A')
  READ(1,*,IOSTAT=ISO) NLCDA,TSCDA,(ISCDA(K),K=1,7)
  WRITE(7,1002)NCARD
  WRITE(7,*) NLCDA,TSCDA,(ISCDA(K),K=1,7)
  IF( ISO > 0 ) GOTO 100
                                                                  
  IF( NLCDA > 0 )THEN
    !C66B*  READ CONCENTRATION DATA ASSIMILATION LOCATIONS AND                                                               
    !      SERIES IDENTIFIERS                                                                                               
    NCARD='66B'
    CALL SEEK('C66B')
    WRITE(7,1002)NCARD
    DO L=1,NLCDA
     READ(1,*,IOSTAT=ISO) ITPCDA(L),ICDA(L),JCDA(L),ICCDA(L),JCCDA(L),(NCSERA(L,K),K=1,7)
     WRITE(7,*)           ITPCDA(L),ICDA(L),JCDA(L),ICCDA(L),JCCDA(L),(NCSERA(L,K),K=1,7)
     IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
                                                                  
  !C67*  READ NEUTRALLY BUOYANT PARTICLE DRIFTER DATA                                                                      
  NCARD='67'
  CALL SEEK('C67')
  READ(1,*,IOSTAT=ISO) ISPD,NPD,NPDRT,NWPD,ISLRPD,ILRPD1,ILRPD2,JLRPD1, JLRPD2, MLRPDRT,IPLRPD
  NPD = 0
  WRITE(7,1002)NCARD
  WRITE(7,*) ISPD,NPD,NPDRT,NWPD,ISLRPD,ILRPD1,ILRPD2,JLRPD1, JLRPD2, MLRPDRT,IPLRPD
  ! *** Not Used:  NPDRT to IPLRPD
  IF( ISO > 0 ) GOTO 100
                                                                  
  IF( NPD > 0 )THEN
    !C68*  READ NEUTRALLY BUOYANT PARTICLE INITIAL POSITIONS                                                                 
    NCARD='68'
    CALL SEEK('C68')
    DO NP=1,NPD
      READ(1,*,IOSTAT=ISO) RI(NP),RJ(NP),RK(NP)
      WRITE(7,1002)NCARD
      WRITE(7,*) RI(NP),RJ(NP),RK(NP)
    ENDDO
  ENDIF
                                                                  
  !C69*  CONSTANTS FOR LONGITUDE AND LATITUDE OF CELL CENTERS                                                              
  NCARD='69'
  CALL SEEK('C69')
  READ(1,*,IOSTAT=ISO) CDLON1,CDLON2,CDLON3,CDLAT1,CDLAT2,CDLAT3
  WRITE(7,1002)NCARD
  WRITE(7,*) CDLON1,CDLON2,CDLON3,CDLAT1,CDLAT2,CDLAT3
  IF( ISO > 0 ) GOTO 100
                                                                  
  !C70*  CONTROLS FOR WRITING ASCII OR BINARY DUMP FILES  (NOT USED)                                                                  
  NCARD='70'
  CALL SEEK('C70')
  READ(1,*,IOSTAT=ISO)ISDUMP,ISADMP,NSDUMP,TSDUMP,TEDUMP,ISDMPP,ISDMPU,ISDMPW,ISDMPT,IADJDMP

  WRITE(7,1002)NCARD
  WRITE(7,*)ISDUMP,ISADMP,NSDUMP,TSDUMP,TEDUMP,ISDMPP,ISDMPU,ISDMPW,ISDMPT,IADJDMP
  IF( ISO > 0 ) GOTO 100
  JSDUMP=1
  NCDUMP=1
                                                                  
  !C71*  CONTROLS FOR HORIZONTAL PLANE SCALAR FIELD CONTOURING                                                             
  NCARD='71'
  CALL SEEK('C71')
  DO NS=1,7
    READ(1,*,IOSTAT=ISO) ISSPH(NS),NPSPH(NS),ISRSPH(NS),ISPHXY(NS)
    WRITE(7,1002)NCARD
    WRITE(7,*) ISSPH(NS),NPSPH(NS),ISRSPH(NS),ISPHXY(NS)
  ENDDO
  IF( ISO > 0 ) GOTO 100

  ! *** SET WATER COLUMN LINKAGE FLAG IF ANY CONSTITUENT OUTPUT IS ENABLED
  DO NS=1,7
    IF( ISTRAN(NS) >= 1 )ISSPH(8)=1
  ENDDO
  IF( ISSPH(4) >= 1 )ISSPH(8)=1
                                                                                    
  !C71A*  CONTROLS FOR HORIZONTAL PLANE SEDIMENT BED PROPERTIES                                                        
  NCARD='71A'
  CALL SEEK('C71A')
  READ(1,*,IOSTAT=ISO) ITMP,ISBEXP,NPBPH

  WRITE(7,1002)NCARD
  WRITE(7,*) ISBEXP,NPBPH
  IF( ISO > 0 ) GOTO 100

  JSBPH=1
  JSBPHA=1

  !C71B*  CONTROLS FOR FOOD CHAIN MODEL OUTPUT                                                                             
  NCARD='71B'
  CALL SEEK('C71B')
  READ(1,*,IOSTAT=ISO) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
  WRITE(7,1002)NCARD
  WRITE(7,*) ISFDCH,NFDCHZ,HBFDCH,TFCAVG
  IF( ISO > 0 ) GOTO 100

  !C72*  CONTROLS FOR EFDC_EXPLORER LINKAGE AND SURFACE ELEVATION RESIDUAL OUTPUT                                              
  JSFDCH=1
  NCARD='72'
  CALL SEEK('C72')
  READ(1,*,IOSTAT=ISO) ISPPH,NPPPH,ISRPPH,IPPHXY

  WRITE(7,1002)NCARD
  WRITE(7,*) ISPPH,NPPPH,ISRPPH,IPPHXY
  IF( ISO > 0 ) GOTO 100
  IF( ISPPH < 0 ) ISPPH = 1
                                                                                                                       
  !C73*  CONTROLS FOR HORIZONTAL PLANE VELOCITY PLOTTING                                                                   
  NCARD='73'
  CALL SEEK('C73')
  READ(1,*,IOSTAT=ISO) ISVPH,NPVPH,ISRVPH,IVPHXY

  WRITE(7,1002)NCARD
  WRITE(7,*) ISVPH,NPVPH,ISRVPH,IVPHXY
  IF( ISO > 0 ) GOTO 100
                                                                                    
  !C74*  CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING                                                               
  NCARD='74'
  CALL SEEK('C74')
  READ(1,*,IOSTAT=ISO) ISECSPV,NPSPV(1),ISSPV(1),ISRSPV(1),ISHPLTV(1)

  WRITE(7,1002)NCARD
  WRITE(7,*) ISECSPV,NPSPV(1),ISSPV(1),ISRSPV(1),ISHPLTV(1)
  SHPLTV(1)=FLOAT(ISHPLTV(1))
  SBPLTV(1)=1.0-SHPLTV(1)
  DO NS=2,7
    READ(1,*,IOSTAT=ISO) IDUMMY,NPSPV(NS),ISSPV(NS),ISRSPV(NS),ISHPLTV(NS)

    WRITE(7,1002)NCARD
    WRITE(7,*) IDUMMY,NPSPV(NS),ISSPV(NS),ISRSPV(NS),ISHPLTV(NS)
    SHPLTV(NS)=FLOAT(ISHPLTV(NS))
    SBPLTV(NS)=1.0-SHPLTV(NS)
  ENDDO
  IF( ISO > 0 ) GOTO 100
                                                                                    
  IF( ISECSPV > 0 )THEN
    !C75*  MORE CONTROLS FOR VERTICAL PLANE SCALAR FIELD CONTOURING                                                          
    NCARD='75'
    CALL SEEK('C75')
    DO IS=1,ISECSPV
      READ(1,*,IOSTAT=ISO) DUM,NIJSPV(IS),CCTITLE(10+IS)

      WRITE(7,1002)NCARD
      WRITE(7,*) DUM,NIJSPV(IS),CCTITLE(10+IS)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                                    
    !C76*  I,J LOCATIONS DEFINING VERTICAL PLANE FOR CONTOURING                                                              
    NCARD='76'
    CALL SEEK('C76')
    DO IS=1,ISECSPV
      DO NPP=1,NIJSPV(IS)
        READ(1,*,IOSTAT=ISO) DUM,ISPV(NPP,IS),JSPV(NPP,IS)

        WRITE(7,1002)NCARD
        WRITE(7,*) DUM,ISPV(NPP,IS),JSPV(NPP,IS)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    ENDDO
  ENDIF
                                                                                    
  !C77*  CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
  NCARD='77'
  CALL SEEK('C77')
  READ(1,*,IOSTAT=ISO) ISECVPV,NPVPV,ISVPV,ISRVPV
  WRITE(7,1002)NCARD
  WRITE(7,*) ISECVPV,NPVPV,ISVPV,ISRVPV
  IF( ISO > 0 ) GOTO 100
                                                                                    
  IF( ISECVPV > 0 )THEN
    !C78*   MORE CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
    NCARD='78'
    CALL SEEK('C78')
    DO IS=1,ISECVPV
      READ(1,*,IOSTAT=ISO) DUM,NIJVPV(IS),ANGVPV(IS),CVTITLE(10+IS)
      WRITE(7,1002)NCARD
      WRITE(7,*) DUM,NIJVPV(IS),ANGVPV(IS),CVTITLE(10+IS)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                                    
    !C79*   MORE CONTROLS FOR VERTICAL PLANE VELOCITY VECTOR PLOTTING
    NCARD='79'
    CALL SEEK('C79')
    DO IS=1,ISECVPV
      DO NPP=1,NIJVPV(IS)
        READ(1,*,IOSTAT=ISO) DUM,IVPV(NPP,IS),JVPV(NPP,IS)
        WRITE(7,1002)NCARD
        WRITE(7,*) DUM,IVPV(NPP,IS),JVPV(NPP,IS)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    ENDDO
  ENDIF
                                                                                    
  !C80*  CONTROLS FOR 3D FIELD OUTPUT
  NCARD='80'
  CALL SEEK('C80')
  READ(1,*,IOSTAT=ISO)IS3DO,ISR3DO,NP3DO,KPC,NWGG,I3DMIN,I3DMAX,J3DMIN,J3DMAX,I3DRW,SELVMAX,BELVMIN
  WRITE(7,1002)NCARD
  WRITE(7,*)IS3DO,ISR3DO,NP3DO,KPC,NWGG,I3DMIN,I3DMAX,J3DMIN,J3DMAX,I3DRW,SELVMAX,BELVMIN
  IF( ISO > 0 ) GOTO 100
  NCALL3D=0
  NRCAL3D=0
  
  !C81* OUTPUT ACTIVATION AND SCALES FOR 3D FIELD OUTPUT                                                                        
  NCARD='81'
  CALL SEEK('C81')
  READ(1,*,IOSTAT=ISO)CDUM,IS3DUUU,JS3DUUU,UUU3DMA,UUU3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DUUU,JS3DUUU,UUU3DMA,UUU3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DVVV,JS3DVVV,VVV3DMA,VVV3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DVVV,JS3DVVV,VVV3DMA,VVV3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DWWW,JS3DWWW,WWW3DMA,WWW3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DWWW,JS3DWWW,WWW3DMA,WWW3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DSAL,JS3DSAL,SAL3DMA,SAL3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DSAL,JS3DSAL,SAL3DMA,SAL3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DTEM,JS3DTEM,TEM3DMA,TEM3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DTEM,JS3DTEM,TEM3DMA,TEM3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DDYE,JS3DDYE,DYE3DMA,DYE3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DDYE,JS3DDYE,DYE3DMA,DYE3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DSED,JS3DSED,SED3DMA,SED3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DSED,JS3DSED,SED3DMA,SED3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DSND,JS3DSND,SND3DMA,SND3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DSND,JS3DSND,SND3DMA,SND3DMI
  IF( ISO > 0 ) GOTO 100
  READ(1,*,IOSTAT=ISO)CDUM,IS3DTOX,JS3DTOX,TOX3DMA,TOX3DMI
  WRITE(7,1002)NCARD
  WRITE(7,*)CDUM,IS3DTOX,JS3DTOX,TOX3DMA,TOX3DMI
  IF( ISO > 0 ) GOTO 100

  IF( ISTRAN(1) < 1 )THEN
    IS3DSAL = 0
    JS3DSAL = 0
  ENDIF
  IF( ISTRAN(2) < 1 )THEN
    IS3DTEM = 0
    JS3DTEM = 0
  ENDIF
  IF( ISTRAN(3) < 1 )THEN
    IS3DDYE = 0
    JS3DDYE = 0
  ENDIF
  IF( ISTRAN(5) < 1 )THEN
    IS3DTOX = 0
    JS3DTOX = 0
  ENDIF
  IF( ISTRAN(6) < 1 )THEN
    IS3DSED = 0
    JS3DSED = 0
  ENDIF
  IF( ISTRAN(7) < 1 )THEN
    IS3DSND = 0
    JS3DSND = 0
  ENDIF
  
  !C82* INPLACE HARMONIC ANALYSIS PARAMETERS
  NCARD='82'
  CALL SEEK('C82')
  READ(1,*,IOSTAT=ISO) ISLSHA,MLLSHA,NTCLSHA,ISLSTR,ISHTA
  WRITE(7,1002)NCARD
  WRITE(7,*) ISLSHA,MLLSHA,NTCLSHA,ISLSTR,ISHTA
  IF( ISO > 0 ) GOTO 100
                                                                                    
  IF( MLLSHA > 0 )THEN
    !C83* HARMONIC ANALYSIS LOCATIONS AND SWITCHES
    NCARD='83'
    CALL SEEK('C83')
    DO M=1,MLLSHA
      READ(1,*,IOSTAT=ISO) ILLSHA(M),JLLSHA(M),LSHAP(M),LSHAB(M),LSHAUE(M),LSHAU(M),CLSL(M)
      WRITE(7,1002)NCARD
      WRITE(7,*) ILLSHA(M),JLLSHA(M),LSHAP(M),LSHAB(M),LSHAUE(M),LSHAU(M),CLSL(M)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
  
  !C84* CONTROLS FOR WRITING TO TIME SERIES FILES                                                                                  
  NCARD='84'
  CALL SEEK('C84')
  READ(1,*,IOSTAT=ISO)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,TCTMSR

  WRITE(7,1002)NCARD
  WRITE(7,*)ISTMSR,MLTMSR,NBTMSR,NSTMSR,NWTMSR,NTSSTSP,TCTMSR
  IF( ISO > 0 ) GOTO 100

  JSTMSR=1
  NCTMSR=1
  JSHYDOUT=1
  NCHYDOUT=1

  IF( NTSSTSP > 0 )THEN
    NCARD='85'
    CALL SEEK('C85')
    DO ITSSS=1,NTSSTSP
      READ(1,*,IOSTAT=ISO)IDUM,MTSSTSP(ITSSS)
      WRITE(7,1002)NCARD
      WRITE(7,*)IDUM,MTSSTSP(ITSSS)
      IF( ISO > 0 ) GOTO 100
    ENDDO
                                                                                    
    NCARD='86'
    CALL SEEK('C86')
    DO ITSSS=1,NTSSTSP
      DO MTSSS=1,MTSSTSP(ITSSS)
        READ(1,*,IOSTAT=ISO)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),TSSTOP(MTSSS,ITSSS)

        WRITE(7,1002)NCARD
        WRITE(7,*)IDUM,IDUM,TSSTRT(MTSSS,ITSSS),TSSTOP(MTSSS,ITSSS)
        IF( ISO > 0 ) GOTO 100
      ENDDO
    ENDDO
  ENDIF
                                                                                    
  IF( MLTMSR > 0 )THEN
    NCARD='87'
    CALL SEEK('C87')
    DO M=1,MLTMSR
      READ(1,*,IOSTAT=ISO)ILTMSR(M),JLTMSR(M),NTSSSS(M),MTMSRP(M),MTMSRC(M),MTMSRA(M),MTMSRUE(M),MTMSRUT(M),MTMSRU(M),MTMSRQE(M),MTMSRQ(M),CLTMSR(M)
      WRITE(7,1002)NCARD
      WRITE(7,*)ILTMSR(M),JLTMSR(M),NTSSSS(M),MTMSRP(M),MTMSRC(M),MTMSRA(M),MTMSRUE(M),MTMSRUT(M),MTMSRU(M),MTMSRQE(M),MTMSRQ(M),CLTMSR(M)
      IF( ISO > 0 ) GOTO 100
    ENDDO
  ENDIF
                                                                                    
  !C88 ** HIGH FREQUENCY OUTPUT
  NCARD='88'
  CALL SEEK('C88')
  READ(1,*,IOSTAT=ISO) HFREOUT
  WRITE(7,1002)NCARD
  WRITE(7,*) HFREOUT
  IF( ISO > 0 ) GOTO 100

  !NCARD='88'
  !CALL SEEK('C88')
  !READ(1,*,IOSTAT=ISO)ISVSFP,MDVSFP,MLVSFP,TMVSFP,TAVSFP
  !WRITE(7,1002)NCARD
  !WRITE(7,*)ISVSFP,MDVSFP,MLVSFP,TMVSFP,TAVSFP
  !IF(ISO > 0 ) GOTO 100
  JSVSFP=1
  
  ! *** NETCDF GENERATION CONTROLS
  NCARD='91'                                                                                                       
  CALL SEEK('C91')                                                                                                 
  READ(1,'(A)',IOSTAT=ISO) STR     
  !                                       NOT      NOT
  READ(STR,*,ERR=100) NCDFOUT,DEFLEV,ROTA,BLK,UTMZ,HREST,BASEDATE,BASETIME,PROJ
  IF (UTMZ < 0) THEN
    UTMZ = ABS(UTMZ)
    HEMI = 2
  ELSE
    HEMI = 1
  ENDIF
  
  WRITE(7,1002)NCARD                                                                                               
  WRITE(7,*) NCDFOUT,DEFLEV,ROTA,BLK,HEMI,UTMZ,HREST,BASEDATE,BASETIME,PROJ                                                              

  IF( NCDFOUT > 0 )THEN
    NCARD='91A'                                                                                                       
    CALL SEEK('C91A')                                                                                                 
    READ(1,*,ERR=100) ISSGLFIL,TBEGNCDF,TENDNCDF                                                                                  
    
    NCARD='91B'                                                                                                       
    CALL SEEK('C91B')                                                                                                 
    READ(1,*,ERR=100) ISNCDF                                                                                     
    DO NS=1,8
      IF( ISTRAN(NS) == 0 ) ISNCDF(NS) = 0
    ENDDO
    
    IF(ISPD < 2 )   ISNCDF(9)  = 0
    IF(NWSER == 0 ) ISNCDF(11) = 0
    IF(ISWAVE < 3 ) ISNCDF(12) = 0
    
    WRITE(7,1002)NCARD                                                                                               
    WRITE(7,*) ISNCDF(1:12)
  ENDIF
  
  NTS=INT8(NTC)*NTSPTC
  NBVSFP=NTC*NTSPTC
  NSVSFP=0
                                                                                    
  DT = TIDALP*FLOAT(NFLTMT)/FLOAT(NTSPTC)
                                                                                    
  GOTO 2000
                                                                                    
  ! *** WRITE INPUT ERROR MESSAGES AND TERMINATE RUN                                                                      
  100 WRITE(6,1001)NCARD                                                                                                
  WRITE(8,1001)NCARD
  WRITE(7,1001)NCARD

  STOP
  
2000 CONTINUE                                                                                                          
                                                                                    
  ! *** NOW REWIND UNIT 1 & READ IN AS CHARACTER TO WRITE TO UNIT 7                                                                                    
  REWIND (1)
  21 READ(1,22,END=24) TEXT                                                                                            
  WRITE (7,23) TEXT
  GOTO 21
  24 CONTINUE                                                                                                          
  CLOSE(1)
  22 FORMAT (A80)                                                                                                      
  23 FORMAT (1X,A80)                                                                                                   
                                                                                    
  ! *** READ CELL TYPES FROM FILES CELL.INP                                                       
  WRITE(*,'(A)')'READING CELL.INP'
  OPEN(1,FILE='cell.inp',STATUS='UNKNOWN')
                                                                                    
  ! *** SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT                                                    
  STRC=READSTR(1)
  READ(STRC,*)JCTMP

  IF( JCTMP/=JC )THEN                                                                                    
    
    ! ***   READ OLD FILE FORMAT                                                                                                                                                                               
    JACROSS=JC
    IF( JC>640)JACROSS=640
    DO JT=1,JC,JACROSS
      JF=JT
      JLAST=JT+JACROSS-1
      IF( JLAST>JC) JLAST=JC
      WRITE (7,8)JF,JLAST
      DO I=1,IC
        READ(1,6,IOSTAT=ISO) (IJCT(I,J),J=JF,JLAST)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CELL.INP'
        WRITE (7,16) (IJCT(I,J),J=JF,JLAST)
      ENDDO
      WRITE(7,15)
    ENDDO
  ELSE
                                                                                    
    IF( IC>640 )THEN
      IACROSS=640
      DO IT=1,IC,IACROSS
        IFIRST=IT
        ILAST=IT+IACROSS-1
        IF( ILAST>IC) ILAST=IC
        WRITE (7,88)IFIRST,ILAST
        DO J=JC,1,-1
          READ(1,66,IOSTAT=ISO)ADUMMY,(IJCT(I,J),I=IFIRST,ILAST)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CELL.INP'
          WRITE (7,166)ADUMMY,(IJCT(I,J),I=IFIRST,ILAST)
        ENDDO
        WRITE(7,15)
      ENDDO
    ELSE
      IFIRST=1
      ILAST=IC
      WRITE (7,88)IFIRST,ILAST
      DO J=JC,1,-1
        READ(1,66,IOSTAT=ISO)ADUMMY,(IJCT(I,J),I=IFIRST,ILAST)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CELL.INP'
        WRITE (7,166)ADUMMY,(IJCT(I,J),I=IFIRST,ILAST)
      ENDDO
      WRITE(7,15)
    ENDIF
  ENDIF
  CLOSE(1)
                                                                                    
  8 FORMAT ('   CELL TYPE ARRAY,J=',I5,2X,'TO J=',I5,//)                                                              
                                                                                    
  !----------------------------------------------------------------------C                                                
                                                                                    
  WRITE(*,'(A)')'READING CELLLT.INP'
  OPEN(1,FILE='celllt.inp',STATUS='UNKNOWN')
                                                                                    
  ! *** SKIP OVER TITLE AND AND HEADER LINES AND DETERMINE FILE FORMAT                                                    
  STRC=READSTR(1)
  READ(STRC,*)JCTMP
                                                                                                                                                                        
  IF( JCTMP/=JC )THEN
    ! ***   READ OLD FILE FORMAT                                                                                            
    JACROSS=JC
    IF( JC>640)JACROSS=640
    DO JT=1,JC,JACROSS
      JF=JT
      JLAST=JT+JACROSS-1
      IF( JLAST>JC) JLAST=JC
      WRITE (7,8)JF,JLAST
      DO I=1,IC
        READ(1,6,IOSTAT=ISO) (IJCTLT(I,J),J=JF,JLAST)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CELLLT.INP'
        WRITE (7,16) (IJCTLT(I,J),J=JF,JLAST)
      ENDDO
      WRITE(7,15)
    ENDDO

  ELSE
    ! ***   READ NEW FILE FORMAT                                                                                            
    IF( IC>640 )THEN
      IACROSS=640
      DO IT=1,IC,IACROSS
        IFIRST=IT
        ILAST=IT+IACROSS-1
        IF( ILAST>IC) ILAST=IC
        WRITE (7,88)IFIRST,ILAST
        DO J=JC,1,-1
          READ(1,66,IOSTAT=ISO)ADUMMY,(IJCTLT(I,J),I=IFIRST,ILAST)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CELLLT.INP'
          WRITE (7,166)ADUMMY,(IJCTLT(I,J),I=IFIRST,ILAST)
        ENDDO
        WRITE(7,15)
      ENDDO
    ELSE
      IFIRST=1
      ILAST=IC
      WRITE (7,88)IFIRST,ILAST
      DO J=JC,1,-1
        READ(1,66,IOSTAT=ISO)ADUMMY,(IJCTLT(I,J),I=IFIRST,ILAST)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CELLLT.INP'
        WRITE (7,166)ADUMMY,(IJCTLT(I,J),I=IFIRST,ILAST)
      ENDDO
      WRITE(7,15)
    ENDIF
  ENDIF
                                                                                    
  CLOSE(1)
                                                                                    
  88 FORMAT ('   CELLLT TYPE ARRAY,I=',I5,2X,'TO I=',I5,//)                                                            
                                                                                    
  ! *** IF ISCONNECT GE 2, READ IN EAST-WEST BOUNDARY CELLS FROM                                                           
  ! *** FILE MAPPGEW.INP TO SPECIFY A PERIODIC DOMAIN IN THE EAST-WEST DIRECTION                                                                                                         
  IF( ISCONNECT >= 2 )THEN
    WRITE(*,'(A)')'READING MAPPGEW.INP'
    OPEN(1,FILE='mappgew.inp',STATUS='UNKNOWN')
                                                                                    
  ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    READ(1,*,IOSTAT=ISO) NPEWBP
    IF( ISO > 0 ) STOP '  READ ERROR FOR FILE MAPPGEW.INP'
    DO NP=1,NPEWBP
      READ(1,*,IOSTAT=ISO) IWPEW(NP),JWPEW(NP),IEPEW(NP),JEPEW(NP)
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE MAPPGEW.INP'
    ENDDO
    CLOSE(1)
  ENDIF
                                                                                    
  ! *** IF ISCONNECT GE 1, READ IN NORTH-SOUTH BOUNDARY CELLS FROM                                                           
  ! *** FILE MAPPGNS.INP TO SPECIFY A PERIODIC DOMAIN IN THE NORTH-SOUTH DIRECTION                                                                                                         
  IF( ISCONNECT == 1 .OR. ISCONNECT == 3 )THEN
    WRITE(*,'(A)')'READING MAPPGNS.INP'
    OPEN(1,FILE='mappgns.inp',STATUS='UNKNOWN')
                                                                                    
  ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    READ(1,*,IOSTAT=ISO) NPNSBP
    IF( ISO > 0 ) STOP '  READ ERROR FOR FILE MAPPGNS.INP'
    DO NP=1,NPNSBP
      READ(1,*,IOSTAT=ISO) ISPNS(NP),JSPNS(NP),INPNS(NP),JNPNS(NP)
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE MAPPGNS.INP'
    ENDDO
    CLOSE(1)
  ENDIF
                                                                                    
  ! *** GENERATE CELL MAPPINGS                                                                                            
                                                                                    
  CALL CELLMAP

  ! *** FORMAT STATEMENTS FOR EFDC.INP
     15 FORMAT (/)                                                                                                        
      6 FORMAT (640I1)                                                                                                    
     66 FORMAT (A5,640I1)                                                                                                 
     16 FORMAT (1X,640I1)                                                                                                 
    166 FORMAT (1X,A5,120I1)                                                                                              
                                                                                    
  ! *** READ CURVILINEAR-ORTHOGONAL OR VARIABLE CELL DATA FROM FILE DXDY.INP                                                            

  ! *** INITIALIZE CELL DIMENSIONS TO CONSTANT CARTESIAN OR DUMMY VALUES                                                  
  DO L=1,LC
    DXP(L)=DX*DXYCVT
    DYP(L)=DY*DXYCVT
    ZBR(L)=ZBRADJ
  ENDDO
                                                                                    
  ! *** READ IN DX, DY, DEPTH AND BOTTOM ELEVATION AT CELL CENTERS OF VARIABLE CELLS                                                     
  LMHK=.FALSE.
  IF( LVC > 0 )THEN
    WRITE(*,'(A)')'READING DXDY.INP'
    OPEN(1,FILE='dxdy.inp',STATUS='UNKNOWN')

    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    IMN = 1000000
    JMN = 1000000
    IMX = 0
    JMX = 0
    IF( ISVEG == 0 )THEN
      DO LT=1,LVC
        READ(1,*,IOSTAT=ISO)I,J,DXIJ,DYIJ,HIJ,BELVIJ,ZBRIJ
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DXDY.INP'
        L=LIJ(I,J)
        IMN = MIN(IMN,I)
        JMN = MIN(JMN,J)
        IMX = MAX(IMX,I)
        JMX = MAX(JMX,J)
        DXP(L)=DXYCVT*DXIJ
        DYP(L)=DXYCVT*DYIJ
        HMP(L)=HADADJ + HCVRT*HIJ
        HMP(L)=MAX(HMP(L),HMIN)
        BELV(L)=BELADJ + BELCVRT*BELVIJ
        BELV1(L)=BELADJ + BELCVRT*BELVIJ
        ZBR(L)=ZBRADJ + ZBRCVRT*ZBRIJ
      ENDDO
    ELSE
      DO LT=1,LVC
        READ(1,*,IOSTAT=ISO)I,J,DXIJ,DYIJ,HIJ,BELVIJ,ZBRIJ,MVEGIJT  ! !SCJ adding MHK devices when MVEGIJT>90
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DXDY.INP'
        L=LIJ(I,J)
        IMN = MIN(IMN,I)
        JMN = MIN(JMN,J)
        IMX = MAX(IMX,I)
        JMX = MAX(JMX,J)
        DXP(L)=DXYCVT*DXIJ
        DYP(L)=DXYCVT*DYIJ
        HMP(L)=HADADJ + HCVRT*HIJ
        HMP(L)=MAX(HMP(L),HMIN)
        BELV(L)=BELADJ + BELCVRT*BELVIJ
        BELV1(L)=BELADJ + BELCVRT*BELVIJ
        ZBR(L)=ZBRADJ + ZBRCVRT*ZBRIJ
        MVEGL(L)=MVEGIJT
        !*** MHK
        IF( MVEGIJT > 90 )THEN
          LMHK=.TRUE.
          ITURB=ITURB+1
          IJLTURB(ITURB,1)=I
          IJLTURB(ITURB,2)=J
          IJLTURB(ITURB,3)=LIJ(I,J)
        ENDIF
      ENDDO
    ENDIF
    CLOSE(1)
  ENDIF

  ! *** OPEN FILE MODDXDY.INP TO MODIFY INPUT VALUES OF DX AND DY                                                         
  IF( IMDXDY > 0 )THEN
    WRITE(*,'(A)')'READING MODDXDY.INP'
    OPEN(1,FILE='moddxdy.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND HEADER LINES                                                                                  
    STR=READSTR(1)   
    READ(1,*) NMDXDY
    IF( NMDXDY >= 1 )THEN
      DO NMD=1,NMDXDY
        READ(1,*)ITMP,JTMP,RMDX,RMDY
        LTMP=LIJ(ITMP,JTMP)
        DXP(LTMP)=RMDX*DXP(LTMP)
        DYP(LTMP)=RMDY*DYP(LTMP)
      ENDDO
    ENDIF
    CLOSE(1)
  ENDIF
                                                                  
  ! *** OPEN FILE MODCHAN.INP TO INSERT SUBGRID CHANNELS INTO HOST CELLS                                                                                                        
  MDCHH=0
  IF( ISCHAN > 0 )THEN
    WRITE(*,'(A)')'READING MODCHAN.INP'
    OPEN(1,FILE='modchan.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    IF( ISCHAN == 1 )THEN
      READ(1,*) MDCHH,MDCHHD,MDCHHD2
      READ(1,*) MDCHITM,MDCHHQ,QCHERR
      IF( MDCHH >= 1 )THEN
        DO NMD=1,MDCHH
          READ(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),IMDCHV(NMD),JMDCHV(NMD)
          QCHANU(NMD)=0.
          QCHANUN(NMD)=0.
          QCHANV(NMD)=0.
          QCHANVN(NMD)=0.
        ENDDO
      ENDIF
    ENDIF
    IF( ISCHAN == 2 )THEN
      READ(1,*) MDCHH,MDCHHD,MDCHHD2
      READ(1,*) MDCHITM,MDCHHQ,QCHERR
      IF( MDCHH >= 1 )THEN
        DO NMD=1,MDCHH
          READ(1,*)MDCHTYP(NMD),IMDCHH(NMD),JMDCHH(NMD),IMDCHU(NMD),JMDCHU(NMD),IMDCHV(NMD),JMDCHV(NMD),CHANLEN(NMD),PMDCH(NMD)
          QCHANU(NMD)=0.
          QCHANUN(NMD)=0.
          QCHANV(NMD)=0.
          QCHANVN(NMD)=0.
        ENDDO
      ENDIF
    ENDIF
    CLOSE(1)
    IF( MDCHH >= 1 )THEN
      DO NMD=1,MDCHH
        LMDCHH(NMD)=LIJ(IMDCHH(NMD),JMDCHH(NMD))
        IF( IMDCHU(NMD) == 1 .AND. JMDCHU(NMD) == 1 )THEN
          LMDCHU(NMD)=1
        ELSE
          LMDCHU(NMD)=LIJ(IMDCHU(NMD),JMDCHU(NMD))
        ENDIF
        IF( IMDCHV(NMD) == 1 .AND. JMDCHV(NMD) == 1 )THEN
          LMDCHV(NMD)=1
        ELSE
          LMDCHV(NMD)=LIJ(IMDCHV(NMD),JMDCHV(NMD))
        ENDIF
      ENDDO
    ENDIF
  ENDIF
                                                                  
  ! *** ENABLE SPATIALLY VARIABLE BACKGROUND AHO
  IF( AHO < 0. )THEN
    AHMAX=ABS(AHO)
    AHMIN=1.0E32
    DO L=2,LA
      AHOXY(L) = ABS(AHO)*DXP(L)*DYP(L)
      AHMAX = MAX(AHOXY(L),AHMAX)
      AHMIN = MIN(AHOXY(L),AHMIN)
    ENDDO
    PRINT '(A,2E16.5)','VARIABLE AHO USED (MIN,MAX): ',AHMIN,AHMAX
  ELSE
    AHOXY = AHO
  ENDIF
  AHO = ABS(AHO)
  
  ! *** ENABLE SPATIALLY VARIABLE SMAGORINSKY AND BACKGROUND DIFFUSIVITY
  IF( AHD < 0. )THEN
    AHMAX=-1.0E32
    AHMIN=1.0E32
    ADMAX=-1.0E32
    ADMIN=1.0E32

    WRITE(*,'(A)')'READING AHMAP.INP'
    OPEN(1,FILE='ahmap.inp')
    STR=READSTR(1)

    DO WHILE(1)
      READ(1,*,END=200) LD,ID,JD,T1,T2
      L = LIJ(ID,JD)
      IF( T1 < 0. )THEN
        AHOXY(L) = ABS(T1)*DXP(L)*DYP(L)
      ELSE
        AHOXY(L) = T1
      ENDIF
      AHDXY(L) = T2

      AHMAX = MAX(AHOXY(L),AHMAX)
      AHMIN = MIN(AHOXY(L),AHMIN)
      ADMAX = MAX(AHDXY(L),ADMAX)
      ADMIN = MIN(AHDXY(L),ADMIN)
    ENDDO
200 CLOSE(1)
    AHD = ABS(AHD)
    PRINT '(A,4E12.4)','VARIABLE AHO & AHD USED (MIN,MAX): ',AHMIN,AHMAX,ADMIN,ADMAX
  ELSE
    AHDXY = AHD
  ENDIF
  
  ! *** ENABLE SPATIALLY VARIABLE VERTICAL VISCOSITY AND DIFFUSIVITY
  IF( AVO < 0. )THEN
    AHMAX=-1.0E32
    AHMIN=1.0E32
    ADMAX=-1.0E32
    ADMIN=1.0E32

    AVO = ABS(AVO)
    ABO = ABS(ABO)

    WRITE(*,'(A)')'READING AVMAP.INP'
    OPEN(1,FILE='avmap.inp')
    STR=READSTR(1)

    DO WHILE(1)
      READ(1,*,END=300) LD,ID,JD,T1,T2
      L = LIJ(ID,JD)
      AVOXY(L) = T1
      AVBXY(L) = T2

      AHMAX = MAX(AVOXY(L),AHMAX)
      AHMIN = MIN(AVOXY(L),AHMIN)
      ADMAX = MAX(AVBXY(L),ADMAX)
      ADMIN = MIN(AVBXY(L),ADMIN)
    ENDDO
300 CLOSE(1)
    PRINT '(A,4E12.4)','VARIABLE AVO & ABO USED (MIN,MAX): ',AHMIN,AHMAX,ADMIN,ADMAX
  ELSE
    AVOXY = AVO
    AVBXY = ABO
  ENDIF

  ! *** OPEN FILE CHANSEC.INP FOR 1-D CHANNEL CROSS SECTION DATA                                                          
  ! *** REMOVED 2004-09-19  PMC                                                                                           
                                                                  
  ! *** OPEN FILE GWATER.INP TO SPECIFY GROUNDWATER INTERACTION                                                           
  ! *** BY INFILTRATION AND EVAPOTRANSPIRATION                                                                            
  ISGWIE = 0
  ! NDL ADDED 2018-03-21 TO FIXED ISSUE DIVIDE BY REZO RNPOR IN LINE 903 OF CALPUV2C.F90
  RNPOR=1.E-12			
  IF( ISGWIT == 1 )THEN
    WRITE(*,'(A)')'READING GWATER.INP'
    OPEN(1,FILE='gwater.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    READ(1,*) ISGWIE
    IF( ISGWIE >= 1 )THEN
      READ(1,*) DAGWZ,RNPOR,RIFTRM
    ELSE
      DAGWZ=0.0
      RNPOR=1.E-12
      RIFTRM=0.0
    ENDIF
    CLOSE(1)
  ENDIF
    339 FORMAT(2I5,6F14.5)                                                                                                
                                                                  
  ! *** OPEN FILE FBODY.INP TO READ IN SPATIALLY VARYING BODY FORCES                                                      
  IF( ISBODYF >= 1 )THEN
    WRITE(*,'(A)')'READING FBODY.INP'
    OPEN(1,FILE='fbody.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    READ(1,*)CVTFACX,CVTFACY
    DO LL=2,LA
      READ(1,*)ITMP,JTMP,FBODY1,FBODY2
      L=LIJ(ITMP,JTMP)
      IF( ISBODYF == 1 )THEN
        DO K=1,KC                                                                                                        
          FBODYFX(L,K)=CVTFACX*FBODY1
          FBODYFY(L,K)=CVTFACY*FBODY2
        ENDDO
      ENDIF
      IF( ISBODYF == 2 )THEN
        DO K=1,KC-1                                                                                                      
          FBODYFX(L,K)=0.0
          FBODYFY(L,K)=0.0
        ENDDO
        FBODYFX(L,KC)=CVTFACX*FBODY1
        FBODYFY(L,KC)=CVTFACY*FBODY2
      ENDIF
    END DO
    DO K=1,KC                                                                                                            
      FBODYFX(1,K)=0.0
      FBODYFY(1,K)=0.0
      FBODYFX(LC,K)=0.0
      FBODYFY(LC,K)=0.0
    END DO
                                                                  
    CLOSE(1)
  ENDIF
                                                                  
  !**********************************************************************C                                                
  ! *** OPEN FILE SEDBLBC.INP TO READ IN SEDIMENT BEDLOAD OUTFLOW                                                         
  ! *** OR RECIRCULATION BOUNDARY CONDITIONS                                                                              
  NSBDLDBC=0
  IF( (ISBDLDBC == 1 .AND. .NOT. LSEDZLJ) .OR. (NCALC_BL >= 2 .AND. LSEDZLJ) )THEN
    WRITE(*,'(A)')'READING SEDBLBC.INP'
    OPEN(1,FILE='sedblbc.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    READ(1,*)NSBDLDBC
    DO NS=1,NSBDLDBC
      READ(1,*)ITMPU,JTMPU,ITMPD,JTMPD,ISDBLDIR(NS)
      LSBLBCU(NS) = LIJ(ITMPU,JTMPU)
      IF( ITMPD > 0 .AND. JTMPD > 0 )THEN
        LSBLBCD(NS) = LIJ(ITMPD,JTMPD)
      ELSE
        LSBLBCD(NS) = 0
      ENDIF

      ISDBLDIR(NS) = ABS(ISDBLDIR(NS))
      IF( ISDBLDIR(NS) == 1 )THEN
        ! *** EAST-WEST
        IF( SUB(LSBLBCU(NS)) < 0.5 .AND. SUB(LEC(LSBLBCU(NS))) > 0.5 )ISDBLDIR(NS) = -1
      ENDIF
      IF( ISDBLDIR(NS) == 2 )THEN
        ! *** NORTH-SOUTH
        IF( SVB(LSBLBCU(NS)) < 0.5 .AND. SVB(LEC(LSBLBCU(NS))) > 0.5 )ISDBLDIR(NS) = -2
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** OPEN FILE GWMAP.INP TO SPECIFY GROUNDWATER INTERACTION BY AMBIENT GROUNDWATER FLOW                                                       
  IF( ISGWIT > 1 )THEN
    !IF( ISGWIT == 2 )THEN
    !  ISGWIE=2
    !ELSE
      ISGWIE = 0
    !ENDIF

    IF( ISGWIT == 3 )THEN
      WRITE(*,'(A)')'READING GWSEEP.INP'
      OPEN(1,FILE='gwseep.inp',STATUS='UNKNOWN')

      STR=READSTR(1)     
      READ(1,*)NSEEPCLASSES
      DO M=1,NSEEPCLASSES
        READ(1,*)IDUM,SEEPRATE(M)
      ENDDO
      CLOSE(1)
    ENDIF
    
    WRITE(*,'(A)')'READING GWMAP.INP'
    OPEN(1,FILE='gwmap.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    DO LL=2,LA
      READ(1,*) ITMP,JTMP,IZONE,RVALUE
      L=LIJ(ITMP,JTMP)
      IF( ISGWIT == 3 )THEN
        IF( IZONE > NSEEPCLASSES )THEN
          WRITE(6,*)'BAD SEEPAGE CLASS AT I,J=',ITMP,JTMP
          STOP
        ENDIF
        QGW(L) = SEEPRATE(IZONE)*RVALUE
      ELSE
        NGWSL(L) = MAX(IZONE,1)
        GWFAC(L) = RVALUE
      ENDIF
    END DO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ IN SPATIALLY VARYING SEDIMENT ROUGHNESS HEIGHT FOR                                                           
  ! *** DETERMINING GRAIN STRESS                                                                                          
  IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN
    IF( ISBEDSTR == 3 )THEN
      WRITE(*,'(A)')'READING SEDROUGH.INP'
      OPEN(1,FILE='sedrough.inp')
      STR=READSTR(1)      
      DO L=2,LC-1
        READ(1,*) LDUM,IDUM,JDUM,ZBRSED(L)
        IF( ZBRSED(L) <= 0.0 )THEN
          STOP ' BAD SEDIMENT ROUGHNESS IN SEDROUGH.INP'
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
                                                                  
  ! *** OPEN FILE DOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT                                                      
  ! *** DISSOLVED ORGANIC CARBON IN WATER COLUMN                                                                          
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1 .OR. ISTOC(NT) == 2)IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTDOCW == 1 )THEN
      WRITE(*,'(A)')'READING DOCW.INP'
      OPEN(1,FILE='docw.inp',STATUS='UNKNOWN')
      STR=READSTR(1)
      READ(1,*)ISALTYP
      IF( ISALTYP == 0 )THEN
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) (STDOCW(L,K),K=1,KC)
          IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCW.INP'
        ENDDO
      ELSE
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STDOCW(L,K),K=1,KC)
          IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCW.INP'
        ENDDO
      ENDIF
    ENDIF
  ENDIF
                                                                  
  ! *** OPEN FILE POCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT                                                      
  ! *** PARTICULATE ORGANIC CARBON IN WATER COLUMN                                                                        
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1)IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTPOCW == 1 )THEN
      WRITE(*,'(A)')'READING POCW.INP'
      OPEN(1,FILE='pocw.inp',STATUS='UNKNOWN')
      STR=READSTR(1)
      READ(1,*)ISALTYP
      IF( ISALTYP == 0 )THEN
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) (STPOCW(L,K),K=1,KC)
          IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCW.INP'
        ENDDO
      ELSE
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STPOCW(L,K),K=1,KC)
          IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCW.INP'
        ENDDO
      ENDIF
    ENDIF
  ENDIF
                                                                  
  ! *** OPEN FILE FPOCW.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT                                                     
  ! *** PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS                                                       
  ! *** IN WATER COLUMN                                                                                                   
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 2 .OR. ISTOC(NT) == 3)IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTPOCW == 3 )THEN
      WRITE(*,'(A)')'READING FPOCW.INP'
      OPEN(1,FILE='fpocw.inp',STATUS='UNKNOWN')
      DO NS=1,NSED+NSND
        DO IS=1,8
          READ(1,*)
        ENDDO
        READ(1,*)ISALTYP
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (STFPOCW(L,K,NS),K=1,KC)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCW.INP'
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STFPOCW(L,K,NS),K=1,KC)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCW.INP'
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  ! *** OPEN FILE DOCB.INP TO SPECIFY SPATIAL VARYING, TIME CONSTANT                                                      
  ! *** DISSOLVED ORGANIC CARBON IN SEDIMENT BED                                                                          
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1 .OR. ISTOC(NT) == 2)IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTDOCB == 1 )THEN
      WRITE(*,'(A)')'READING DOCB.INP'
      OPEN(1,FILE='docb.inp',STATUS='UNKNOWN')
      STR=READSTR(1)
      READ(1,*)ISALTYP,IREAD,KBINPUT
      DO K=1,KB
        DO L=2,LC-1
          STDOCB(L,K)=0.0
        ENDDO
      ENDDO
      IF( IREAD == 0 )THEN
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) STDOCB(L,1)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCB.INP'
            DO K=2,KB
              STDOCB(L,K)=STDOCB(L,1)
            ENDDO
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STDOCB(L,1)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCB.INP'
            DO K=2,KB
              STDOCB(L,K)=STDOCB(L,1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      IF( IREAD == 1 )THEN
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (STDOCB(L,K),K=1,KB)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCB.INP'
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STDOCB(L,K),K=1,KB)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCB.INP'
          ENDDO
        ENDIF
      ENDIF
      IF( IREAD == 2 )THEN
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (STDOCB(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCB.INP'
            DO K=KBINPUT,KB
              STDOCB(L,K)=STDOCB(L,KBINPUT)
            ENDDO
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STDOCB(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE DOCB.INP'
            DO K=KBINPUT,KB
              STDOCB(L,K)=STDOCB(L,KBINPUT)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** OPEN FILE POCB.INP TO READ SPATIALY VARYING, TIME CONSTANT                                                        
  ! *** PARTICULATE ORGANIC CARBON IN BED                                                                                 
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 1)IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTPOCB == 1 )THEN
      WRITE(*,'(A)')'READING POCB.INP'
      OPEN(1,FILE='pocb.inp',STATUS='UNKNOWN')
      STR=READSTR(1)
      READ(1,*)ISALTYP,IREAD,KBINPUT
      DO K=1,KB
        DO L=2,LC-1
          STPOCB(L,K)=0.0
        ENDDO
      ENDDO
      IF( IREAD == 0 )THEN
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) STPOCB(L,1)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
            DO K=2,KB
              STPOCB(L,K)=STPOCB(L,1)
            ENDDO
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STPOCB(L,1)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
            DO K=2,KB
              STPOCB(L,K)=STPOCB(L,1)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      IF( IREAD == 1 )THEN
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (STPOCB(L,K),K=1,KB)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STPOCB(L,K),K=1,KB)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
          ENDDO
        ENDIF
      ENDIF
      IF( IREAD == 2 )THEN
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (STPOCB(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
            DO K=KBINPUT,KB
              STPOCB(L,K)=STPOCB(L,KBINPUT)
            ENDDO
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STPOCB(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP 'READ ERROR FOR FILE POCB.INP'
            DO K=KBINPUT,KB
              STPOCB(L,K)=STPOCB(L,KBINPUT)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** OPEN FILE FPOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT                                                       
  ! *** PARTICULATE ORGANIC CARBON FRACTION FOR EACH SEDIMENT CLASS                                                       
  ! *** IN BED                                                                                                            
  IVAL=0
  DO NT=1,NTOX
    IF( ISTOC(NT) == 2 .OR. ISTOC(NT) == 3)IVAL=1
  ENDDO
  IF( IVAL == 1 )THEN
    IF( ISTPOCB == 3 )THEN
      WRITE(*,'(A)')'READING FPOCB.INP'
      OPEN(1,FILE='fpocb.inp',STATUS='UNKNOWN')
      DO NS=1,NSED+NSND
        DO IS=1,8
          READ(1,*)
        ENDDO
        READ(1,*)ISALTYP,IREAD,KBINPUT
                                                                  
        IF( IREAD == 0 )THEN
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) STFPOCB(L,1,NS)
              IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCB.INP'
              DO K=2,KB
                STFPOCB(L,K,NS)=STFPOCB(L,1,NS)
              ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,STFPOCB(L,1,NS)
              IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCB.INP'
              DO K=2,KB
                STFPOCB(L,K,NS)=STFPOCB(L,1,NS)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
        IF( IREAD == 1 )THEN
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STFPOCB(L,K,NS),K=1,KB)
              IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCB.INP'
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STFPOCB(L,K,NS),K=1,KB)
              IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCB.INP'
            ENDDO
          ENDIF
        ENDIF
        IF( IREAD == 2 )THEN
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (STFPOCB(L,K,NS),K=1,KBINPUT)
              IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCB.INP'
              DO K=KBINPUT,KB
                STFPOCB(L,K,NS)=STFPOCB(L,KBINPUT,NS)
              ENDDO
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(STFPOCB(L,K,NS),K=1,KBINPUT,NS)
              IF( ISO > 0 ) STOP 'READ ERROR FOR FILE FPOCB.INP'
              DO K=KBINPUT,KB
                STFPOCB(L,K,NS)=STFPOCB(L,KBINPUT,NS)
              ENDDO
            ENDDO
          ENDIF
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF

  !**********************************************************************C                                                
  !###########################################################################                                            
  ! HQI Change to include sptially varying, but time constant bulk foc                                                    
  ! FPOCB  - Bulk foc from data                                                                                           
  ! PFPOCB - Pseudo foc from data, to be used for all partitioning calculations                                           
  ! RM, 02/29/04                                                                                                          
  !**********************************************************************C                                                
                                                                  
  ! *** OPEN FILE FOCB.INP TO READ SPATIALY VARYING, TIME CONSTANT                                                        
  ! *** PARTICULATE ORGANIC CARBON IN BED AND PSEUDO-POC IN BED                                                           
  IF( ISTPOCB == 4 )THEN
                                                                  
    WRITE(*,'(A)')'READING FOCB.INP'
    OPEN(1,FILE='focb.inp',STATUS='UNKNOWN')
    STR=READSTR(1)
    READ(1,*)ISALTYP,IREAD,KBINPUT
    DO K=1,KB
      DO L=2,LC-1
        FPOCB(L,K)=0.0
      ENDDO
    ENDDO
    DO L=2,LC-1
       READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(FPOCB(L,K),K=1,KBINPUT)
       DO K=1,KBINPUT
          FPOCB(L,K) = FPOCB(L,K)/1000000.
       ENDDO
       DO K=KBINPUT+1,KB
          FPOCB(L,K)=FPOCB(L,KBINPUT)
       ENDDO
       IF( ISO > 0 ) STOP '  READ ERROR FOR FILE FOCB.INP'
    ENDDO
    CLOSE(1)
                                                                  
    WRITE(*,'(A)')'READING PSEUDO_FOCB.INP'
    OPEN(1,FILE='pseudo_focb.inp',STATUS='UNKNOWN')
    STR=READSTR(1)
    READ(1,*)ISALTYP,IREAD,KBINPUT
    DO L=2,LC-1
       READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(PFPOCB(L,K),K=1,KBINPUT)
       DO K=1,KBINPUT
          PFPOCB(L,K) = PFPOCB(L,K)/1000000.
       ENDDO
       DO K=KBINPUT+1,KB
          PFPOCB(L,K)=PFPOCB(L,KBINPUT)
       ENDDO
       IF( ISO > 0 ) STOP '  READ ERROR FOR FILE PSEUDO_FOCB.INP'
    ENDDO
    CLOSE(1)
  ENDIF
                        
  IF( ISVHEAT > 0 .AND. ISTOPT(2) /= 5 ) ISVHEAT = 0
  DO L=2,LA  
    IF( ISTOPT(2) == 1 )THEN
      SVKEBACK(L) = SWRATNF
    ELSE
      SVKEBACK(L) = WQKEB(1)
    ENDIF
  ENDDO
    
  IF( ISTRAN(2) > 0 .AND. ISTOPT(2) == 5 .AND. ISVHEAT > 0 )THEN
    WRITE(*,'(A)')'READING SVHTFACT.INP'
    OPEN(1,FILE='svhtfact.inp',STATUS='UNKNOWN')

    ! *** DEFAULT VALUES
    DO L=2,LA  
      SVREVC(L)=0.001*ABS(REVC)
      SVRCHC(L)=0.001*ABS(RCHC)
      IF( REVC < 0. ) LSVHTWINDE(L) = .TRUE.
      IF( RCHC < 0. ) LSVHTWINDC(L) = .TRUE.
    ENDDO 

    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    
    ! *** LOOP OVER THE DOMAIN AND ONLY SET THE CELLS THAT ARE IN THE FILE
    DO NS=2,LA
      READ(1,*,END=101) LL,ITMP,JTMP,RMDX,RMDY,RAD
      L=LIJ(ITMP,JTMP)
      LSVHTWINDE(L) = RMDX < 0.
      LSVHTWINDC(L) = RMDY < 0.
      SVREVC(L) = 0.001*ABS(RMDX)
      SVRCHC(L) = 0.001*ABS(RMDY)
      SVKEBACK(L) = RAD
    ENDDO
    101 CLOSE(1)
  ENDIF

  ! *** READ IN INITIAL SALINITY, TEMPERATURE, DYE, SED, SND, TOX                                                         
  ! *** FOR COLD STARTS FORM FILE XXXX.INP                                                                                
  ! *** SALINITY                                                                                                          
  DO K=1,KC
    DO L=2,LA
      SALINIT(L,K)=0.
    ENDDO
  ENDDO

  IF( ISTRAN(1) >= 1 .AND. (ISRESTI == 0                  .OR. &
                        (ISRESTI >= 1 .AND. ISCI(1) == 0) .OR.  &
                        (ISTOPT(1)>1)) )THEN  ! *** PMC SINGLE LINE - FORCE IC
    IF( ISTOPT(1) >= 1 )THEN
      WRITE(*,'(A)')'READING SALT.INP'
      OPEN(1,FILE='salt.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)ISALTYP
      IF( ISALTYP == 0 )THEN
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) (SALINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SALT.INP'
        ENDDO
      ELSE
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SALINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SALT.INP'
        ENDDO
      ENDIF
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** TEMPERATURE                                                                                                       
  DO K=1,KC
    DO L=2,LA
      TEMINIT(L,K)=TEMO
    ENDDO
  ENDDO

  IF( ISTRAN(2) >= 1 .AND. (ISRESTI == 0 .OR. (ISRESTI >= 1 .AND. ISCI(2) == 0) .OR. (ISTOPT(2)>9) ) )THEN 
    IF( ISTOPT(2) >= 1 .OR. INITTEMP > 0 )THEN
      WRITE(*,'(A)')'READING TEMP.INP'
      OPEN(1,FILE='temp.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)ISALTYP
      IF( ISALTYP == 0 )THEN
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) (TEMINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TEMP.INP'
        ENDDO
      ELSE
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(TEMINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TEMP.INP'
        ENDDO
      ENDIF
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** DYE                                                                                                               
  DO K=1,KC
    DO L=2,LA
      DYEINIT(L,K)=0.
    ENDDO
  ENDDO
  IF( ISTRAN(3) >= 1 .AND. (ISRESTI == 0 .OR. (ISRESTI >= 1 .AND. ISCI(3) == 0) .OR. (ISTOPT(3)>1)) )THEN 
    IF( ISTOPT(3) >= 1 )THEN
      WRITE(*,'(A)')'READING DYE.INP'
      OPEN(1,FILE='dye.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)ISALTYP
      IF( ISALTYP == 0 )THEN
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) (DYEINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DYE.INP'
        ENDDO
      ELSE
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(DYEINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DYE.INP'
        ENDDO
      ENDIF
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** SFL                                                                                                               
  DO K=1,KC
    DO L=2,LA
      SFLINIT(L,K)=0.
    ENDDO
  ENDDO
  IF( ISRESTI == 0 .AND. ISTRAN(4) >= 1 )THEN
    IF( ISTOPT(4) >= 1 )THEN
      WRITE(*,'(A)')'READING SFL.INP'
      OPEN(1,FILE='sfl.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)ISALTYP
      IF( ISALTYP == 0 )THEN
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) (SFLINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFL.INP'
        ENDDO
      ELSE
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SFLINIT(L,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFL.INP'
        ENDDO
      ENDIF
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** TOXICS                                                                                                            
  IF( ISTRAN(5) >= 1 )THEN
    DO NT=1,NTOX
      DO K=1,KC
        DO L=2,LA
          TOXINIT(L,K,NT)=TOXINTW(NT)
        ENDDO
      ENDDO
    ENDDO
    DO NT=1,NTOX
      DO K=1,KB
        DO L=2,LA
          TOXBINIT(L,K,NT)=TOXINTB(NT)
        ENDDO
      ENDDO
    ENDDO
    IISTMP=1
    IF( ISRESTI == 0 ) IISTMP=0
    IF( ISRESTI >= 1 .AND. ISCI(5) == 0 ) IISTMP=0
    IF( IISTMP == 0 .AND. ISTRAN(5) >= 1 )THEN
      IF( ISLTMT == 0 )THEN
        WRITE(*,'(A)')'READING TOXW.INP'
        OPEN(1,FILE='toxw.inp',STATUS='UNKNOWN')
        IF( ITXINT(1) == 1 .OR. ITXINT(1) == 3 )THEN
          DO NT=1,NTOX
            ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
            DO IS=1,8
              READ(1,*)
            ENDDO
            READ(1,*)ISALTYP,ITOXWU(NT)
            IF( ISALTYP == 0 )THEN
              DO L=2,LC-1
                READ(1,*,IOSTAT=ISO) (TOXINIT(L,K,NT),K=1,KC)
                IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TOXW.INP'
              ENDDO
            ELSE
              DO L=2,LC-1
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(TOXINIT(L,K,NT),K=1,KC)
                IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TOXW.INP'
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
    ENDIF
    IISTMP=1
    IF( ISRESTI == 0 ) IISTMP=0
    IF( ISRESTI >= 1 .AND. ISCI(5) == 0 ) IISTMP=0
    IF( IISTMP == 0 .AND. ISTRAN(5) >= 1 )THEN
      IF( ISLTMT == 0. )THEN
        WRITE(*,'(A)')'READING TOXB.INP'
        OPEN(1,FILE='toxb.inp',STATUS='UNKNOWN')
        IF( ITXINT(1) == 2 .OR. ITXINT(1) == 3 )THEN
          DO NT=1,NTOX
                                                                  
            ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
            DO IS=1,8
              READ(1,*)
            ENDDO
            ! *** ITOXBU is not used
            READ(1,*)ISALTYP,ITOXBU(NT),KBINPUT
            IF( ISALTYP == 0 )THEN
              DO L=2,LC-1
                READ(1,*,IOSTAT=ISO) (TOXBINIT(L,K,NT),K=1,KBINPUT)
                IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TOXB.INP'
              ENDDO
            ELSE
              DO L=2,LC-1
                READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(TOXBINIT(L,K,NT),K=1,KBINPUT)
                IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TOXB.INP'
              ENDDO
            ENDIF
          ENDDO
        ENDIF
        CLOSE(1)
      ENDIF
    ENDIF
  ENDIF
                                                                  
  ! *** COHESIVE SEDIMENT                                                                                                 
  IF( ISTRAN(6)  >= 1 )THEN
    DO NS=1,NSED
      DO K=1,KC
        DO L=2,LA
          SEDINIT(L,K,NS)=SEDO(NS)
        ENDDO
      ENDDO
    ENDDO
    DO NS=1,NSED
      DO K=1,KB
        DO L=2,LA
          SEDBINIT(L,K,NS)=SEDBO(NS)
        ENDDO
      ENDDO
    ENDDO
    ITXINTT=0
    IF( ISEDINT == 1) ITXINTT=1
    IF( ISEDINT == 3) ITXINTT=1
    IISTMP=1
    IF( ISRESTI == 0 ) IISTMP=0
    IF( ISRESTI >= 1 .AND. ISCI(6) == 0 ) IISTMP=0
    IF( IISTMP == 0 .AND. ISTRAN(6) >= 1 )THEN
      IF( ITXINTT >= 1 )THEN
        WRITE(*,'(A)')'READING SEDW.INP'
        OPEN(1,FILE='sedw.inp',STATUS='UNKNOWN')
        DO NS=1,NSED
          ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
          STR=READSTR(1)
          READ(1,*)ISALTYP
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (SEDINIT(L,K,NS),K=1,KC)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SEDW.INP'
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SEDINIT(L,K,NS),K=1,KC)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SEDW.INP'
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
    ENDIF
    ITXINTT=0
    IF( ISEDINT == 2) ITXINTT=1
    IF( ISEDINT == 3) ITXINTT=1
    IISTMP=1
    IF( ISRESTI == 0 ) IISTMP=0
    IF( ISRESTI >= 1 .AND. ISCI(6) == 0 ) IISTMP=0
    IF( IISTMP == 0 .AND. ISTRAN(6) >= 1 .AND. .NOT. LSEDZLJ )THEN
      IF( ITXINTT >= 1 )THEN
        WRITE(*,'(A)')'READING SEDB.INP'
        OPEN(1,FILE='sedb.inp',STATUS='UNKNOWN')
        DO NS=1,NSED
                                                                  
          ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
          STR=READSTR(1)
          READ(1,*)ISALTYP,ISEDBU(NS),KBINPUT
          DO K=1,KB
            DO L=2,LC-1
              SEDBINIT(L,K,NS)=0.0
            ENDDO
          ENDDO
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (SEDBINIT(L,K,NS),K=1,KBINPUT)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SEDB.INP'
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SEDBINIT(L,K,NS),K=1,KBINPUT)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SEDB.INP'
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
    ENDIF
  ENDIF
                                                                  
  ! *** NON-COHESIVE SEDIMENT                                                                                             
  IF( ISTRAN(7) >= 1 )THEN
    DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KC
        DO L=2,LA
          SNDINIT(L,K,NX)=SEDO(NS)
        ENDDO
      ENDDO
    ENDDO
    DO NX=1,NSND
      NS=NX+NSED
      DO K=1,KB
        DO L=2,LA
          SNDBINIT(L,K,NX)=SEDBO(NS)
        ENDDO
      ENDDO
    ENDDO
    ITXINTT=0
    IF( ISEDINT == 1) ITXINTT=1
    IF( ISEDINT == 3) ITXINTT=1
    IISTMP=1
    IF( ISRESTI == 0) IISTMP=0
    IF( ISRESTI >= 1 .AND. ISCI(7) == 0) IISTMP=0
    IF( IISTMP == 0 .AND. ISTRAN(7) >= 1 )THEN
      IF( ITXINTT >= 1 )THEN
        WRITE(*,'(A)')'READING SNDW.INP'
        OPEN(1,FILE='sndw.inp',STATUS='UNKNOWN')
        DO NX=1,NSND
          ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
          STR=READSTR(1)
          READ(1,*)ISALTYP
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (SNDINIT(L,K,NX),K=1,KC)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNDW.INP'
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SNDINIT(L,K,NX),K=1,KC)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNDW.INP'
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
    ENDIF
    ITXINTT=0
    IF( ISEDINT == 2) ITXINTT=1
    IF( ISEDINT == 3) ITXINTT=1
    IISTMP=1
    IF( ISRESTI == 0) IISTMP=0
    IF( ISRESTI >= 1 .AND. ISCI(7) == 0) IISTMP=0
    IF( IISTMP == 0 .AND. ISTRAN(7) >= 1 )THEN
      IF( ITXINTT >= 1 )THEN
        WRITE(*,'(A)')'READING SNDB.INP'
        OPEN(1,FILE='sndb.inp',STATUS='UNKNOWN')
        DO NX=1,NSND
          ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
          STR=READSTR(1)
          DO K=1,KB
            DO L=2,LC-1
              SNDBINIT(L,K,NX)=0.0
            ENDDO
          ENDDO
          READ(1,*)ISALTYP,ISNDBU(NX),KBINPUT
          IF( ISALTYP == 0 )THEN
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO) (SNDBINIT(L,K,NX),K=1,KBINPUT)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNDB.INP'
            ENDDO
          ELSE
            DO L=2,LC-1
              READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(SNDBINIT(L,K,NX),K=1,KBINPUT)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNDB.INP'
            ENDDO
          ENDIF
        ENDDO
        CLOSE(1)
      ENDIF
    ENDIF
  ENDIF
                                                                  
  !  ** SEDIMENT BED MECHANICAL INITIAL CONDITIONS                                                                        
  IISTMP=1
  IF( ISRESTI == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(6) == 0 ) IISTMP=0
  IF( ISRESTI >= 1 .AND. ISCI(7) == 0 ) IISTMP=0
  IF( (ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 ) .AND. ISEDINT>1 )THEN
                                                                  
    !  ** BED LAYER THICKNESS                                                                                               
    IF( IISTMP == 0 .AND. IBMECH >= 1 )THEN
      WRITE(*,'(A)')'READING BEDLAY.INP'
      OPEN(1,FILE='bedlay.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)IBEDLAYU,ISALTYP,KBINPUT
      IF( IBEDLAYU > 0 )THEN
        DO K=1,KB
          DO L=2,LC-1
            BEDLINIT(L,K)=0.0
          ENDDO
        ENDDO
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (BEDLINIT(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE BEDLAY.INP'
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(BEDLINIT(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE BEDLAY.INP'
          ENDDO
        ENDIF
      ENDIF
      CLOSE(1)
    ENDIF
                                                                  
    !  ** BED LAYER BULK DENSITY                                                                                            
    IF( IISTMP == 0 .AND. IBMECH >= 1 )THEN
      WRITE(*,'(A)')'READING BEDBDN.INP'
      OPEN(1,FILE='bedbdn.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)IBEDBDNU,ISALTYP,KBINPUT
      IF( IBEDBDNU > 0 )THEN
        DO K=1,KB
          DO L=2,LC-1
            BEDBINIT(L,K)=0.0
          ENDDO
        ENDDO
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (BEDBINIT(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE BEDBDN.INP'
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(BEDBINIT(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE BEDBDN.INP'
          ENDDO
        ENDIF
      ENDIF
      CLOSE(1)
    ENDIF
                                                                  
    !  ** BED LAYER DRY DENSITY, POROSITY OR VOID RATIO                                                                     
    IF( IISTMP == 0 .AND. IBMECH >= 1 )THEN
      WRITE(*,'(A)')'READING BEDDDN.INP'
      OPEN(1,FILE='bedddn.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)IBEDDDNU,ISALTYP,KBINPUT
      IF( IBEDDDNU > 0 )THEN
        DO K=1,KB
          DO L=2,LC-1
            BEDDINIT(L,K)=0.0
          ENDDO
        ENDDO
        IF( ISALTYP == 0 )THEN
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO) (BEDDINIT(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE BEDDDN.INP'
          ENDDO
        ELSE
          DO L=2,LC-1
            READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,(BEDDINIT(L,K),K=1,KBINPUT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE BEDDDN.INP'
          ENDDO
        ENDIF
      ENDIF
      CLOSE(1)
    ENDIF
                                                                  
    !  ** CONSOLIDATION MAP                                                                                                 
    IF( IBMECH == 9 )THEN
      WRITE(*,'(A)')'READING CONSOLMAP.INP'
      OPEN(1,FILE='consolmap.inp',STATUS='UNKNOWN')
                                                                  
      ! ***   SKIP OVER TITLE AND AND HEADER LINES                                                                            
      STR=READSTR(1)
      READ(1,*)ISALTYP
      IF( ISALTYP == 0 )THEN
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO) LCONSOL(L)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CONSOLMAP.INP'
        ENDDO
      ELSE
        DO L=2,LC-1
          READ(1,*,IOSTAT=ISO)LDUM,IDUM,JDUM,LCONSOL(L)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE CONSOLMAP.INP'
        ENDDO
      ENDIF
      CLOSE(1)
    ENDIF
  ENDIF
     19 FORMAT (/,' INITIAL BUOYANCY ARRAY:',//)                                                                          
    907 FORMAT(12F6.2)                                                                                                    
                                                                  
  ! *** READ IN OPEN BOUNDARY SURFACE ELEVATION TIME SERIES FROM THE                                                      
  ! *** FILE PSER.INP                                                                                                     
  IF( NPSER >= 1 )THEN
    WRITE(*,'(A)')'READING PSER.INP'
    OPEN(1,FILE='pser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    DO NS=1,NPSER
      READ(1,*,IOSTAT=ISO) ITYPE,MPSER(NS),TCPSER(NS),TAPSER(NS),RMULADJ,ADDADJ,PSERZDF(NS),INTPSER(NS)
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE PSER.INP'
      PSERZDF(NS)=G*PSERZDF(NS)
     IF( ITYPE == 1 )THEN
        READ(1,*,IOSTAT=ISO) RMULADJS,ADDADJS,PSERZDS(NS)
       IF( ISO > 0 ) STOP '  READ ERROR FOR FILE PSER.INP'
        PSERZDS(NS)=G*PSERZDS(NS)
      ELSE
        RMULADJS=0
          ADDADJS=0.
          PSERZDS(NS)=0.
      ENDIF
      IF( ITYPE == 0 )THEN
        DO M=1,MPSER(NS)
          READ(1,*,IOSTAT=ISO)TPSER(M,NS),PSERTMP
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE PSER.INP'
          TPSER(M,NS)=TPSER(M,NS)+TAPSER(NS)
          PSER(M,NS)=G*(PSERTMP+ADDADJ)*RMULADJ
          PSERS(M,NS)=0.0
        ENDDO
      ELSEIF( ITYPE == 1 )THEN
        DO M=1,MPSER(NS)
          READ(1,*,IOSTAT=ISO)TPSER(M,NS),PSERTMP,PSERTMPS
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE PSER.INP'
          TPSER(M,NS)=TPSER(M,NS)+TAPSER(NS)
          PSER(M,NS)=G*(PSERTMP+ADDADJ)*RMULADJ
          PSERS(M,NS)=G*(PSERTMPS+ADDADJS)*RMULADJS
       ENDDO
      ENDIF

    ENDDO
    CLOSE(1)
  ENDIF
   6776 FORMAT(A20)                                                                                                       
                                                                  
  ! *** READ IN VOLUMETRIC SOURCE OR RIVER INFLOW TIME SERIES FROM THE                                                    
  ! *** FILE QSER.INP                                                                                                     
  IF( NQSER >= 1 )THEN
    WRITE(*,'(A)')'READING QSER.INP'
    OPEN(1,FILE='qser.inp',STATUS='UNKNOWN')
    
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    DO NS=1,NQSER
      ISMOOTH=0
      READ(1,*,IOSTAT=ISO)ISTYP, MQSER(NS),TCQSER(NS),TAQSER(NS),RMULADJ,ADDADJ,ICHGQS
      IF( MQSER(NS)<0 )THEN
         ISMOOTH=1
         MQSER(NS)=-MQSER(NS)
      ENDIF
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QSER.INP'

      IF( ISTYP == 1 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QSER.INP'
        DO M=1,MQSER(NS)
          READ(1,*,IOSTAT=ISO)QSER(NS).TIM(M),QSERTMP
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QSER.INP'
          QSER(NS).TIM(M) = QSER(NS).TIM(M)+TAQSER(NS)
          QSERTMP = RMULADJ*(QSERTMP+ADDADJ)
          IF( ICHGQS == 1)  QSERTMP = MAX(QSERTMP,0.0)
          IF( ICHGQS == -1) QSERTMP = MIN(QSERTMP,0.0)
          DO K=1,KC
            QSER(NS).VAL(M,K) = QSERTMP*WKQ(K)
          ENDDO
        ENDDO
      ELSE
        DO M=1,MQSER(NS)
          READ(1,*,IOSTAT=ISO)QSER(NS).TIM(M),(QSER(NS).VAL(M,K), K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QSER.INP'
          QSER(NS).TIM(M) = QSER(NS).TIM(M)+TAQSER(NS)
          DO K=1,KC
            QSER(NS).VAL(M,K) = RMULADJ*(QSER(NS).VAL(M,K)+ADDADJ)
            IF( ICHGQS == 1)  QSER(NS).VAL(M,K) = MAX(QSER(NS).VAL(M,K),0.0)
            IF( ICHGQS == -1) QSER(NS).VAL(M,K) = MIN(QSER(NS).VAL(M,K),0.0)
          ENDDO
        ENDDO
      ENDIF

      IF( ISMOOTH == 1 )THEN
        DO K=1,KC
          DO M=1,MQSER(NS)
            QSERSM(M,K) = QSER(NS).VAL(M,K)
          ENDDO
        ENDDO
        DO K=1,KC
          DO M=2,MQSER(NS)-1
            QSER(NS).VAL(M,K) = 0.25*(QSERSM(M-1,K)+QSERSM(M+1,K))+0.5*QSERSM(M,K)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
   2222 FORMAT(2I5,F12.7,F12.4)                                                                                           
                                                                  
  ! *** READ IN FLOW WITHDRAWL-RETURN FLOW AND CONCENTRATION RISE                                                         
  ! *** TIME SERIES FROM THE FILE QWRS.INP                                                                                
  IF( NQWRSR >= 1 )THEN
    WRITE(*,'(A)')'READING QWRS.INP'
    OPEN(1,FILE='qwrs.inp',STATUS='UNKNOWN')
                                                                  
    NCTMP=4+NSED+NSND+NTOX
    ipmc = NWQV
    IF( ISTRAN(8) > 0 )NCTMP=NCTMP+ipmc   ! *** IF NQWV CHANGES THIS !SHOULD BE UPDATED

    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    DO NS=1,NQWRSR
      READ(1,*,IOSTAT=ISO)ISTYP,MQWRSR(NS),TCQWRSR(NS),TAQWRSR(NS),RMULADJ,ADDADJ
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QWRS.INP'
      IF( ISTYP == 0 )THEN
        ! *** FLOW ONLY.  NO RISE/FALL
        DO NC=1,NCTMP
          DO M=1,MQWRSR(NS)
            CQWRSER(M,NS,NC)=0.
          ENDDO
        ENDDO
        DO M=1,MQWRSR(NS)
          READ(1,*,IOSTAT=ISO)TQWRSER(M,NS),QWRSER(M,NS)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QWRS.INP'
          TQWRSER(M,NS)=TQWRSER(M,NS)+TAQWRSR(NS)
          QWRSER(M,NS)=(RMULADJ*(QWRSER(M,NS)+ADDADJ))
        ENDDO
      ELSE
        ! *** FLOW WITH RISE/FALL
        DO M=1,MQWRSR(NS)
          READ(1,*,IOSTAT=ISO)TQWRSER(M,NS),QWRSER(M,NS),(CQWRSER(M,NS,NC),NC=1,NCTMP)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QWRS.INP'
          TQWRSER(M,NS)=TQWRSER(M,NS)+TAQWRSR(NS)
          QWRSER(M,NS)=(RMULADJ*(QWRSER(M,NS)+ADDADJ))
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ IN GROUNDWATER INFLOW/OUTFLOW AND CONCENTRATION TIME                                                         
  ! *** SERIES FROM THE FILE GWSER.INP                                                                                    
  IF( ISGWIT == 2 )THEN
    WRITE(*,'(A)')'READING GWSER.INP'
    OPEN(1,FILE='gwser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    NCTMP=4+NSED+NSND+NTOX
    STR=READSTR(1)
    READ(1,*)NGWSER
    IF( NGWSER > 0 )THEN
      DO NS=1,NGWSER
        READ(1,*,IOSTAT=ISO) MGWSER(NS),TCGWSER(NS),TAGWSER(NS),RMULADJ,ADDADJ,IGWSER(NS)
        
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE GWSER.INP'
        DO M=1,MGWSER(NS)
          READ(1,*,IOSTAT=ISO)TGWSER(M,NS),GWSER(M,NS),(GWCSER(M,NS,NC),NC=1,NCTMP)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE GWSER.INP'
          TGWSER(M,NS) = TGWSER(M,NS)+TAGWSER(NS)
          GWSER(M,NS)  = RMULADJ*(GWSER(M,NS) + ADDADJ)
        ENDDO
      ENDDO
    ENDIF
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ IN SPATIAL MAPS AND TIME SERIES FOR EXTERNAL SPECIFICATION OF                                                
  ! *** PARTICULATE ORGANIC CARBON FOR USE IN TOXIC CONTAMINANT SORPTION                                                  
  ! *** DISSOLVED ORGANIC CARBON                                                                                          
  ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
  !        READ(1,*)                                                                                                      
  !        READ(1,*)NOCSER                                                                                                
  !            READ(1,*,IOSTAT=ISO)MOCSER(NS),TCOCSER(NS),TAOCSER(NS),                                                    
  !              READ(1,*,IOSTAT=ISO)TOCSER(M,NS),DOCWSER(M,NS),                                                          
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SALINITY TIME SERIES                                                   
  ! *** FROM THE FILE SSER.INP                                                                                            
 IF( NCSER(1) >= 1 )THEN
    WRITE(*,'(A)')'READING SSER.INP'
    OPEN(1,FILE='sser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    NC=1
    DO NS=1,NCSER(NC)
      READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SSER.INP'
      IF( ISTYP == 1 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SSER.INP'
        DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO) TSSAL(NS).TIM(M),CSERTMP
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SSER.INP'
          TSSAL(NS).TIM(M)=TSSAL(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSSAL(NS).VAL(M,K) =(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
          ENDDO
        ENDDO
      ELSE
        DO M=1,MCSER(NS,NC)
         READ(1,*,IOSTAT=ISO) TSSAL(NS).TIM(M),(TSSAL(NS).VAL(M,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SSER.INP'
          TSSAL(NS).TIM(M)=TSSAL(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSSAL(NS).VAL(M,K) =RMULADJ*(TSSAL(NS).VAL(M,K)+ADDADJ)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE TEMPERATURE TIME                                                       
  ! *** SERIES FROM THE FILE TSER.INP                                                                                     
 IF( NCSER(2) >= 1 .OR. ( ISTRAN(5) > 0 .AND. ITOXTEMP > 1 ) )THEN
    WRITE(*,'(A)')'READING TSER.INP'
    OPEN(1,FILE='tser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    NC=2
    DO NS=1,NCSER(NC)
      READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TSER.INP'
      IF( ISTYP == 1 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TSER.INP'
        DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO) TSTEM(NS).TIM(M),CSERTMP
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TSER.INP'
          TSTEM(NS).TIM(M)=TSTEM(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSTEM(NS).VAL(M,K) =(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
          ENDDO
        ENDDO
      ELSE
        DO M=1,MCSER(NS,NC)
         READ(1,*,IOSTAT=ISO) TSTEM(NS).TIM(M),(TSTEM(NS).VAL(M,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TSER.INP'
          TSTEM(NS).TIM(M)=TSTEM(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSTEM(NS).VAL(M,K)=RMULADJ*(TSTEM(NS).VAL(M,K)+ADDADJ)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE DYE CONCENTRATION                                                      
  ! *** TIME SERIES FROM THE FILE DSER.INP                                                                                
 IF( NCSER(3) >= 1 )THEN
    WRITE(*,'(A)')'READING DSER.INP'
    OPEN(1,FILE='dser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    NC=3
    DO NS=1,NCSER(NC)
      READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DSER.INP'
      IF( ISTYP == 1 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DSER.INP'
        DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO) TSDYE(NS).TIM(M),CSERTMP
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DSER.INP'
          TSDYE(NS).TIM(M)=TSDYE(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSDYE(NS).VAL(M,K) =(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
          ENDDO
        ENDDO
      ELSE
        DO M=1,MCSER(NS,NC)
         READ(1,*,IOSTAT=ISO) TSDYE(NS).TIM(M),(TSDYE(NS).VAL(M,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE DSER.INP'
          TSDYE(NS).TIM(M)=TSDYE(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSDYE(NS).VAL(M,K) =RMULADJ*(TSDYE(NS).VAL(M,K)+ADDADJ)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE COHESIVE SEDIMENT                                                      
  ! *** CONCENTRATION TIME SERIES FROM THE FILE SDSER.INP                                                                 
                                                                  
  IF( NSED > 0 )THEN
    NFSED=MSVSED(1)
   IF( NCSER(NFSED) >= 1 )THEN
      WRITE(*,'(A)')'READING SDSER.INP'
      OPEN(1,FILE='sdser.inp',STATUS='UNKNOWN')
                                                                  
      ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
      STR=READSTR(1)
      NC=MSVSED(1)
      DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADS(1),ADDADS(1)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SDSER.INP'
        IF( NSED>1 )THEN
          DO NT=2,NSED
            READ(1,*,IOSTAT=ISO) RMULADS(NT),ADDADS(NT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SDSER.INP'
            NTT=NT-1
            MCSER(NS,NC+NTT)=MCSER(NS,NC)
            TCCSER(NS,NC+NTT)=TCCSER(NS,NC)
            TACSER(NS,NC+NTT)=TACSER(NS,NC)
          ENDDO
        ENDIF
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SDSER.INP'
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSSED(NS,1).TIM(M),CSERTMP
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SDSER.INP'
            TSSED(NS,1).TIM(M)=TSSED(NS,1).TIM(M)+TACSER(NS,NC)
            DO K=1,KC
             TSSED(NS,1).VAL(M,K) =(RMULADS(1)*(CSERTMP+ADDADS(1)))*WKQ(K)
            ENDDO
            DO NT=2,NSED
              !NTT=NT-1
              TSSED(NS,NT).TIM(M)=TSSED(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO) CSERTMP
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SDSER.INP'
              DO K=1,KC
                TSSED(NS,NT).VAL(M,K) =(RMULADS(NT)*(CSERTMP+ADDADS(NT)))*WKQ(K)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSSED(NS,1).TIM(M),(TSSED(NS,1).VAL(M,K),K=1,KC)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SDSER.INP'
            TSSED(NS,1).TIM(M)=TSSED(NS,1).TIM(M)+TACSER(NS,NC)
            DO K=1,KC
              TSSED(NS,1).VAL(M,K) = RMULADS(1)*(TSSED(NS,1).VAL(M,K)+ADDADS(1))
            ENDDO
            DO NT=2,NSED
              !NTT=NT-1
              TSSED(NS,NT).TIM(M)=TSSED(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO)(TSSED(NS,NT).VAL(M,K), K=1,KC)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SDSER.INP'
              DO K=1,KC
                TSSED(NS,NT).VAL(M,K) =RMULADS(NT)*(TSSED(NS,NT).VAL(M,K)+ADDADS(NT))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** CHECK SEDIMENT SERIES                                                                                             
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE NONCOHESIVE SEDIMENT                                                   
  ! *** CONCENTRATION TIME SERIES FROM THE FILE SNSER.INP                                                                 
  IF( NSND > 0 )THEN
    NFSND=MSVSND(1)
   IF( NCSER(NFSND) >= 1 )THEN
      WRITE(*,'(A)')'READING SNSER.INP'
      OPEN(1,FILE='snser.inp',STATUS='UNKNOWN')
                                                                  
      ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
      STR=READSTR(1)
      NC=MSVSND(1)
      DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADS(1),ADDADS(1)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNSER.INP'
        IF( NSND>1 )THEN
          DO NT=2,NSND
            READ(1,*,IOSTAT=ISO)RMULADS(NT),ADDADS(NT)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNSER.INP'
            NTT=NT-1
            MCSER(NS,NC+NTT)=MCSER(NS,NC)
            TCCSER(NS,NC+NTT)=TCCSER(NS,NC)
            TACSER(NS,NC+NTT)=TACSER(NS,NC)
          ENDDO
        ENDIF
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNSER.INP'
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSSND(NS,1).TIM(M),CSERTMP
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNSER.INP'
            TSSND(NS,1).TIM(M)=TSSND(NS,1).TIM(M)+TACSER(NS,NC)
            DO K=1,KC
             TSSND(NS,1).VAL(M,K)=(RMULADS(1)*(CSERTMP+ADDADS(1)))*WKQ(K)
            ENDDO
            DO NT=2,NSND
              !NTT=NT-1
              TSSND(NS,NT).TIM(M)=TSSND(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO)CSERTMP
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNSER.INP'
              DO K=1,KC
                TSSND(NS,NT).VAL(M,K) =(RMULADS(NT)*(CSERTMP+ADDADS(NT)))*WKQ(K)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSSND(NS,1).TIM(M),(TSSND(NS,1).VAL(M,K),K=1,KC)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNSER.INP'
            TSSND(NS,1).TIM(M)=TSSND(NS,1).TIM(M)+TACSER(NS,NC)
            DO K=1,KC
              TSSND(NS,1).VAL(M,K)=RMULADS(1)*(TSSND(NS,1).VAL(M,K)+ADDADS(1))
            ENDDO
            DO NT=2,NSND
              !NTT=NT-1
              TSSND(NS,NT).TIM(M)=TSSND(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO)(TSSND(NS,NT).VAL(M,K), K=1,KC)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SNSER.INP'
              DO K=1,KC
                TSSND(NS,NT).VAL(M,K) =RMULADS(NT)*(TSSND(NS,NT).VAL(M,K)+ADDADS(NT))
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** CHECK SEDIMENT SERIES                                                                                             
                                                                  
   2001 FORMAT(3I5,2F12.5)                                                                                                
                                                                  
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE TOXIC CONTAMINANT                                                      
  ! *** CONCENTRATION TIME SERIES FROM THE FILE TXSER.INP                                                                 
  IF( NTOX > 0 )THEN
    NFTOX=MSVTOX(1)
   IF( NCSER(NFTOX) >= 1 )THEN
      WRITE(*,'(A)')'READING TXSER.INP'
      OPEN(1,FILE='txser.inp',STATUS='UNKNOWN')
                                                                  
      ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
      STR=READSTR(1)
      NC=MSVTOX(1)
      DO NS=1,NCSER(NC)
        READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TXSER.INP'
        DO NT=2,NTOX
          NTT=NT-1
          MCSER(NS,NC+NTT)=MCSER(NS,NC)
          TCCSER(NS,NC+NTT)=TCCSER(NS,NC)
          TACSER(NS,NC+NTT)=TACSER(NS,NC)
        ENDDO
        IF( ISTYP == 1 )THEN
          READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TXSER.INP'
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSTOX(NS,1).TIM(M),CSERTMP
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TXSER.INP'
            TSTOX(NS,1).TIM(M)=TSTOX(NS,1).TIM(M)+TACSER(NS,NC)
            DO K=1,KC
              TSTOX(NS,1).VAL(M,K)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
            ENDDO
            DO NT=2,NTOX
              !NTT=NT-1
              TSTOX(NS,NT).TIM(M)=TSTOX(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO) CSERTMP
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TXSER.INP'
              DO K=1,KC
               TSTOX(NS,NT).VAL(M,K)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
              ENDDO
            ENDDO
          ENDDO
        ELSE
          DO M=1,MCSER(NS,NC)
            READ(1,*,IOSTAT=ISO) TSTOX(NS,1).TIM(M),(TSTOX(NS,1).VAL(M,K), K=1,KC)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TXSER.INP'
            TSTOX(NS,1).TIM(M)=TSTOX(NS,1).TIM(M)+TACSER(NS,NC)
            DO K=1,KC
              TSTOX(NS,1).VAL(M,K)=RMULADJ*(TSTOX(NS,1).VAL(M,K)+ADDADJ)
            ENDDO
            DO NT=2,NTOX
              !NTT=NT-1
              TSTOX(NS,NT).TIM(M)=TSTOX(NS,1).TIM(M)
              READ(1,*,IOSTAT=ISO)(TSTOX(NS,NT).VAL(M,K), K=1,KC)
              IF( ISO > 0 ) STOP '  READ ERROR FOR FILE TXSER.INP'
              DO K=1,KC
                TSTOX(NS,NT).VAL(M,K)=RMULADJ*(TSTOX(NS,NT).VAL(M,K) +ADDADJ)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      CLOSE(1)
    ENDIF
  ENDIF
                                                                  
  ! *** READ IN OPEN BOUNDARY OR VOLUMETRIC SOURCE SHELL FISH LARVAE                                                      
  ! *** TIME SERIES FROM THE FILE SFSER.INP                                                                               
 IF( NCSER(4) >= 1 )THEN
    WRITE(*,'(A)')'READING SFSER.INP'
    OPEN(1,FILE='sfser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    NC=7
    DO NS=1,NCSER(NC)
      READ(1,*,IOSTAT=ISO)ISTYP,MCSER(NS,NC),TCCSER(NS,NC),TACSER(NS,NC),RMULADJ,ADDADJ
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFSER.INP'
      IF( ISTYP == 1 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFSER.INP'
        DO M=1,MCSER(NS,NC)
          READ(1,*,IOSTAT=ISO) TSSFL(NS).TIM(M),CSERTMP
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFSER.INP'
          TSSFL(NS).TIM(M)=TSSFL(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSSFL(NS).VAL(M,K)=(RMULADJ*(CSERTMP+ADDADJ))*WKQ(K)
          ENDDO
        ENDDO
      ELSE
        DO M=1,MCSER(NS,NC)
         READ(1,*,IOSTAT=ISO) TSSFL(NS).TIM(M),(TSSFL(NS).VAL(M,K),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFSER.INP'
          TSSFL(NS).TIM(M)=TSSFL(NS).TIM(M)+TACSER(NS,NC)
          DO K=1,KC
            TSSFL(NS).VAL(M,K)=RMULADJ*(TSSFL(NS).VAL(M,K)+ADDADJ)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ IN FREE SURFACE ELEVATION OR PRESSURE CONTROLLED FLOW                                                        
  ! *** SPECIFICATION FROM THE FILE QCTL.INP                                                                              
  ! *** THE FLOW IS GIVEN BY:                                                                                             
  !         FREE SURFACE                                                                                                  
  !         FREE SURFACE                                                                                                  
  !        FLOW=0                                                                                                         
  !      ELSE                                                                                                             
  !        ENTER QCTL(M,K,NS) VS HDIFCTL(M,NS) TABLE WITH DELH TO GIVE                                                    
  IF( NQCTLT >= 1 )THEN
    WRITE(*,'(A)')'READING QCTL.INP'
    OPEN(1,FILE='qctl.inp',STATUS='UNKNOWN')
    IF( DEBUG )THEN
      OPEN(99,FILE=OUTDIR//'qctlck.inp',STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE=OUTDIR//'qctlck.inp',STATUS='UNKNOWN')
    ENDIF
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    DO NS=1,NQCTLT
      READ(1,*, IOSTAT=ISO)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
      IF( DEBUG )THEN
        WRITE(99,991)NS
        WRITE(99,992)ISTYP,MQCTL(NS),HCTLUA(NS),HCTLUM(NS),HCTLDA(NS),HCTLDM(NS),RMULADJ,ADDADJ,AQCTL(NS)
      ENDIF
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTL.INP'
      IF( ISTYP == 0 )THEN
        DO M=1,MQCTL(NS)
          READ(1,*,IOSTAT=ISO) HDIFCTL(M,NS),(QCTL(M,1,K,NS),K=1,KC)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTL.INP'
          DO K=1,KC
            QCTL(M,1,K,NS)=RMULADJ*(QCTL(M,1,K,NS)+ADDADJ)
          ENDDO
        ENDDO
      ENDIF
      IF( ISTYP == 1 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTL.INP'
        DO M=1,MQCTL(NS)
          READ(1,*,IOSTAT=ISO) HDIFCTL(M,NS),QCTLTMP
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTL.INP'
          DO K=1,KC
            QCTL(M,1,K,NS)=RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
          ENDDO
        ENDDO
      ENDIF
      IF( ISTYP == 2 )THEN
        DO MD=1,MQCTL(NS)
          DO MU=1,MQCTL(NS)
            READ(1,*,IOSTAT=ISO) HDIFCTL(MU,NS),HDIFCTD(MD,NS),(QCTL(MU,MD,K,NS),K=1,KC)
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTL.INP'
            DO K=1,KC
              QCTL(MU,MD,K,NS)=RMULADJ*(QCTL(MU,MD,K,NS)+ADDADJ)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF( ISTYP == 3 )THEN
        READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTL.INP'
        DO MD=1,MQCTL(NS)
          DO MU=1,MQCTL(NS)
           READ(1,*,IOSTAT=ISO)HDIFCTL(MU,NS),HDIFCTD(MD,NS),QCTLTMP
            IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTL.INP'
            DO K=1,KC
              QCTL(MU,MD,K,NS)=RMULADJ*(QCTLTMP+ADDADJ)*WKQ(K)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF( DEBUG )THEN
        IF( ISTYP <= 1 )THEN
          DO M=1,MQCTL(NS)
            WRITE(99,993)M,HDIFCTL(M,NS),(QCTL(M,1,K,NS),K=1,KC)
          ENDDO
        ENDIF
        IF( ISTYP >= 2 )THEN
          DO MD=1,MQCTL(NS)
            DO MU=1,MQCTL(NS)
              WRITE(99,994)MU,MD,HDIFCTL(MU,NS),HDIFCTD(MD,NS),(QCTL(MU,MD,K,NS),K=1,KC)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDDO
    CLOSE(1)
    IF( DEBUG)CLOSE(99)
  ENDIF

  ! *** READ CONTROL TIME-SERIES
  IF( NQCTLSER > 0 )THEN
    ALLOCATE(QCTLSER(NQCTLSER))
    WRITE(*,'(A)')'READING QCTLSER.INP'
    OPEN(1,FILE='qctlser.inp',STATUS='UNKNOWN')
    IF( DEBUG )THEN
      OPEN(99,FILE=OUTDIR//'qctlser.log',STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE=OUTDIR//'qctlser.log',STATUS='UNKNOWN')
    ENDIF
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    DO NS=1,NQCTLSER
      READ(1,*, IOSTAT=ISO) QCTLSER(NS).ITYPE,NX,QCTLSER(NS).TMUL,QCTLSER(NS).TADD,(ISQCTRL(M),M=1,6)
      QCTLSER(NS).PARAM = 0
      IVAL = 1
      DO M=1,6
        IF(ISQCTRL(M) > 0) QCTLSER(NS).PARAM = QCTLSER(NS).PARAM + IVAL
        IVAL = 2*IVAL
      ENDDO
      IF( DEBUG )THEN
        WRITE(99,991) NS
        WRITE(99,992) QCTLSER(NS).ITYPE,NX,QCTLSER(NS).TMUL,QCTLSER(NS).TADD, QCTLSER(NS).PARAM
      ENDIF
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTLSER.INP'
      QCTLSER(NS).COUNT = NX
      ALLOCATE(QCTLSER(NS).TIME(NX))
      ALLOCATE(QCTLSER(NS).ID(NX))
      ALLOCATE(QCTLSER(NS).HEIGHT(NX))
      ALLOCATE(QCTLSER(NS).WIDTH(NX))
      ALLOCATE(QCTLSER(NS).SILL(NX))
      ALLOCATE(QCTLSER(NS).NUM(NX))
      ALLOCATE(QCTLSER(NS).FLOW(NX))
      !ALLOCATE(QCTLSER(NS).RATE(NX))
      DO M=1,NX
          READ(1,*,IOSTAT=ISO) QCTLTMP,QCTLSER(NS).HEIGHT(M),QCTLSER(NS).WIDTH(M),QCTLSER(NS).SILL(M), &
          QCTLSER(NS).NUM(M),QCTLSER(NS).FLOW(M)  !,QCTLSER(NS).RATE(M)
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCTLSER.INP'
          QCTLSER(NS).TIME(M)=QCTLTMP*QCTLSER(NS).TMUL+QCTLSER(NS).TADD
      ENDDO
    ENDDO
    CLOSE(1)
    IF( DEBUG)CLOSE(99)
  ENDIF
  
  ! *** READ CONTROL RULES
  IF( NQCRULES > 0 )THEN
    ALLOCATE(QCRULES(NQCRULES))
    
    WRITE(*,'(A)')'READING QCRULES.INP'
    OPEN(1,FILE='qcrules.inp',STATUS='UNKNOWN')
    IF( DEBUG )THEN
      OPEN(99,FILE=OUTDIR//'qcrules.log',STATUS='UNKNOWN')
      CLOSE(99,STATUS='DELETE')
      OPEN(99,FILE=OUTDIR//'qcrules.log',STATUS='UNKNOWN')
    ENDIF
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    DO NS=1,NQCRULES
      READ(1,*, IOSTAT=ISO) NX,NX,(ISQCTRL(M),M=1,6)  !,QCRULES(NS).PARAM
      QCRULES(NS).PARAM = 0
      IVAL = 1
      DO M=1,6
        IF(ISQCTRL(M) > 0) QCRULES(NS).PARAM = QCRULES(NS).PARAM + IVAL
        IVAL = 2*IVAL
      ENDDO
      IF( DEBUG )THEN
        WRITE(99,*) NS,NX,QCRULES(NS).PARAM
      ENDIF
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCRULES.INP'
            
      QCRULES(NS).NTRGON = 0
      QCRULES(NS).NTROFF = 0
      
      ALLOCATE(RULES(NX,9))
      ALLOCATE(IDX(NX))
      DO M=1,NX
          READ(1,*,IOSTAT=ISO) RULES(M,1),ISTYP,RULES(M,4),RULES(M,5),RULES(M,6),NDUM,RULES(M,8),RULES(M,9),IVAL
          IF( ISO > 0 ) STOP '  READ ERROR FOR FILE QCRULES.INP'
          IDX(M) = M
          RULES(M,2) = ISTYP
          RULES(M,3) = IVAL
          RULES(M,7) = NDUM
          IF (ISTYP == 0) THEN
            QCRULES(NS).NTROFF = QCRULES(NS).NTROFF + 1
          ELSE
            QCRULES(NS).NTRGON = QCRULES(NS).NTRGON + 1
          ENDIF
      ENDDO
      ! SORT TRIGGER LEVELS DESCENDING
      DO M=1,NX-1
        DO L =M+1,NX
          IF(RULES(IDX(L),1) > RULES(IDX(M),1)) THEN
            IVAL = IDX(L)
            IDX(L) = IDX(M)
            IDX(M) = IVAL
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(QCRULES(NS).TRGON(QCRULES(NS).NTRGON))
      ALLOCATE(QCRULES(NS).TROFF(QCRULES(NS).NTROFF))
      L = 0
      DO M=1,NX
        IF (RULES(M,2) == 1.) THEN
          L = L + 1
          QCRULES(NS).TRGON(L).LEVEL = RULES(M,1)             ! TRIGGER LEVEL
          QCRULES(NS).TRGON(L).STATE = INT(RULES(M,2))        ! 0 = OFF, 1 = ON
          QCRULES(NS).TRGON(L).ID = INT(RULES(M,3))           ! INDEX OF RATING CURVE/MASK FOR CONTROL PARAMETERS
          QCRULES(NS).TRGON(L).HEIGHT = RULES(M,4)            ! OPENING HEIGHT (M), FOR UPWARD OPENING
          QCRULES(NS).TRGON(L).WIDTH = RULES(M,5)             ! OPENING WIDTH (M), FOR SIDEWARD OPENING
          QCRULES(NS).TRGON(L).SILL = RULES(M,6)              ! SILL LEVEL CHANGE (M), FOR DOWNWARD OPENING
          QCRULES(NS).TRGON(L).UNITS = INT(RULES(M,7))        ! NUMBER OF GATES, PUMP UNITS
          QCRULES(NS).TRGON(L).FLOW = RULES(M,8)              ! FLOW RATE (M3/S), FOR PUMPS
          QCRULES(NS).TRGON(L).RATE = RULES(M,9)              ! RATE OF OPENING/CLOSING
          IF( DEBUG )THEN
            WRITE(99,*) L,QCRULES(NS).TRGON(L).LEVEL,QCRULES(NS).TRGON(L).STATE, &
              QCRULES(NS).TRGON(L).ID,QCRULES(NS).TRGON(L).HEIGHT,QCRULES(NS).TRGON(L).WIDTH, &
              QCRULES(NS).TRGON(L).SILL,QCRULES(NS).TRGON(L).UNITS,QCRULES(NS).TRGON(L).FLOW, &
              QCRULES(NS).TRGON(L).RATE
          ENDIF          
        ENDIF
      ENDDO
      L = 0
      DO M=NX,1,-1
        IF (RULES(M,2) == 0.) THEN
          L = L + 1
          QCRULES(NS).TROFF(L).LEVEL = RULES(M,1)             ! TRIGGER LEVEL
          QCRULES(NS).TROFF(L).STATE = INT(RULES(M,2))        ! 0 = OFF, 1 = ON
          QCRULES(NS).TROFF(L).ID = INT(RULES(M,3))           ! INDEX OF RATING CURVE/MASK FOR CONTROL PARAMETERS
          QCRULES(NS).TROFF(L).HEIGHT = RULES(M,4)            ! OPENING HEIGHT (M), FOR UPWARD OPENING
          QCRULES(NS).TROFF(L).WIDTH = RULES(M,5)             ! OPENING WIDTH (M), FOR SIDEWARD OPENING
          QCRULES(NS).TROFF(L).SILL = RULES(M,6)              ! SILL LEVEL CHANGE (M), FOR DOWNWARD OPENING
          QCRULES(NS).TROFF(L).UNITS = INT(RULES(M,7))        ! NUMBER OF GATES, PUMP UNITS
          QCRULES(NS).TROFF(L).FLOW = RULES(M,8)              ! FLOW RATE (M3/S), FOR PUMPS
          QCRULES(NS).TROFF(L).RATE = RULES(M,9)              ! RATE OF OPENING/CLOSING
          IF( DEBUG )THEN
            WRITE(99,*) L,QCRULES(NS).TROFF(L).LEVEL,QCRULES(NS).TROFF(L).STATE, &
              QCRULES(NS).TROFF(L).ID,QCRULES(NS).TROFF(L).HEIGHT,QCRULES(NS).TROFF(L).WIDTH, &
              QCRULES(NS).TROFF(L).SILL,QCRULES(NS).TROFF(L).UNITS,QCRULES(NS).TROFF(L).FLOW, &
              QCRULES(NS).TROFF(L).RATE
          ENDIF 
        ENDIF
      ENDDO
      DEALLOCATE(RULES,IDX)
    ENDDO
    CLOSE(1)
    IF( DEBUG)CLOSE(99)
  ENDIF
  
  ! *** CHECK CONTROL DATA
  IF (NQCTL > 0 .AND. (NQCTLSER > 0 .OR. NQCRULES > 0)) THEN
    DO L=1,NQCTL
      IVAL = HSCTL(L).ID
      IF (HSCTL(L).ITYPE == 1) THEN    
        ! *** STRUCTURE IS CONTROLLED BY TIME-SERIES
        IF (IVAL < 1 .OR. IVAL > NQCTLSER) THEN
          STOP 'DATA ERROR: TIME-SERIES INDEX FOR HYDRAULIC STRUCTURE CONTROL IS INVALID!'
        ENDIF
        !IF (QCTLSER(IVAL).TIME(1) > TBEGIN .OR. QCTLSER(IVAL).TIME(QCTLSER(IVAL).COUNT) < TEND) STOP '  DATA ERROR'
      ELSEIF (HSCTL(L).ITYPE == 2 .OR. HSCTL(L).ITYPE == 3) THEN
        ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES 
        IF (IVAL < 1 .OR. IVAL > NQCRULES) THEN
          STOP 'DATA ERROR: RULES INDEX FOR HYDRAULIC STRUCTURE CONTROL IS INVALID!'
        ENDIF
        IF ((QCRULES(IVAL).PARAM .AND. 32) == 32) THEN
          ! *** SET OF LOOKUP TABLE
          DO I=1,QCRULES(IVAL).NTRGON
            ID = QCRULES(IVAL).TRGON(I).ID
            IF (ID < 1 .OR. ID > NQCTLT) STOP '  DATA ERROR: INVALID INDEX FOR LOOKUP TABLE DEFINED IN CONTROL RULES'
          ENDDO
          DO I=1,QCRULES(IVAL).NTROFF 
            ID = QCRULES(IVAL).TROFF(I).ID
            IF (ID < 1 .OR. ID > NQCTLT) STOP '  DATA ERROR: INVALID INDEX FOR LOOKUP TABLE DEFINED IN CONTROL RULES'
          ENDDO
        ELSEIF (QCRULES(IVAL).PARAM > 0) THEN
          ! *** GATE OPENING
        ELSE
          ! *** PUMP FLOWS
        ENDIF
        IF (HSCTL(L).ITYPE == 2) THEN
          ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON SPECIFIC WATER SURFACE ELEVATION
          !LU = LIJ(HSCTL(NCTL).IQCTL1, HSCTL(NCTL).JQCTL1)
        ELSEIF (HSCTL(L).ITYPE == 3) THEN
          ! *** STRUCTURE IS CONTROLLED BY OPERATION RULES DEPENDING ON DIFFERENCE OF WATER SURFACE ELEVATIONS
          !LU = LIJ(HSCTL(NCTL).IQCTL1, HSCTL(NCTL).JQCTL1)
          !LD = LIJ(HSCTL(NCTL).IQCTL2, HSCTL(NCTL).JQCTL2)      
        ENDIF
      ELSE
        ! *** STRUCTURE IS UNCONTROLLED
      ENDIF
    ENDDO
  ENDIF
  IF (NQWR > 0 .AND. NQCRULES > 0) THEN
    DO L=1,NQWR
      IVAL = WRCTL(L).ID
      IF (WRCTL(L).ITYPE == 1 .OR. WRCTL(L).ITYPE == 2) THEN
        ! *** W/R IS CONTROLLED BY OPERATION RULES 
        IF (IVAL < 1 .OR. IVAL > NQCRULES) THEN
          STOP 'DATA ERROR: RULES INDEX FOR WITHDRAWAL/RETURN CONTROL IS INVALID!'
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  
    991 FORMAT(/,'CONTROL TABLE NS =',I5,/)                                                                               
    992 FORMAT(2I5,7F10.4)                                                                                                
    993 FORMAT(I5,11F10.4)                                                                                                
    994 FORMAT(2I5,11F10.4)                                                                                               
   1001 FORMAT(/,'READ ERROR FROM FILE EFDC.INP ON CARD ',A3/)                                                            
   1002 FORMAT(/,'INPUT ECHO NCARD = ',A/)                                                                                
                                                                  
  DO L=1,LC
    PATMT(L)=1000.
    TATMT(L)=0.
    RAINT(L)=0.
    EVAPT(L)=0.
    SOLSWRT(L)=1.  ! *** Address SUNDAY.INP Option
    CLOUDT(L)=0.
    RHAT(L)=0.
    VPAT(L)=0.
    CLEVAP(L)=0.
    CCNHTT(L)=0.
  ENDDO

  ! *** DEACTIVATE THE ATMOSPHERIC FILE IF TEMPERATURE IS NOT SIMUATED
  !IF( ISTRAN(2) == 0 )NASER = 0  2015-06-22 DEPRECATED TO ALLOW EVAP AND RAIN
  LDAYLIGHT = .FALSE.
  IF( NASER > 0 )THEN
    WRITE(*,'(A)')'READING ASER.INP'
    OPEN(1,FILE='aser.inp',STATUS='UNKNOWN')
    
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    CALL SKIPCOM(1,'*')    
    READ(1,'(A)')TEXT
    IASERVER=PARSE_REAL(TEXT)*1000.
    STR=READSTR(1)
    
    IF( IASERVER < 7300 )THEN
      IEVAP=-1
    ENDIF
    
    WRITE(*,'(A,I5)')'  NUMBER OF ATMOSPHERIC SERIES=',NASER
    WRITE(*,'(A,L3)')'  COMPUTESOLRAD=',COMPUTESOLRAD
    WRITE(*,'(A,F10.2)')'  DS_LONG=',DS_LONG
    WRITE(*,'(A,F10.2)')'  DS_LAT=',DS_LAT

    DO NS=1,NASER
      READ(1,*,IOSTAT=ISO) M,TCASER(NS),TAASER(NS),IRELH(NS),RAINCVT,EVAPCVT,SOLRCVT,CLDCVT
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE ASER.INP'
      IF( NS == 1 .AND. IEVAP == -1 )THEN
        IEVAP=1                     ! *** USE ASER DATA
        IF( EVAPCVT < 0 ) IEVAP=2   ! *** LEGACY INPUT TO COMPUTE EVAPORATION USING EFDC ORIGINAL APPROACH
      ELSEIF( IEVAP == 0 )THEN
        EVAPCVT=0.
        RAINCVT=0.
      ENDIF
      WRITE(*,'(A,I5,I10,F6.1)')'  SERIES, NPTS =',NS,TSATM(NS).NREC

      ! *** These parameters are read in for every series but only the last is actually used
      IF( IASERVER<7300 )THEN
        READ(1,*,IOSTAT=ISO) IASWRAD,REVC,RCHC,SWRATNF,SWRATNS,FSWRATF,DABEDT,TBEDIT,HTBED1,HTBED2
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE ASER.INP'
      ENDIF
      
      DO M=1,TSATM(NS).NREC
        READ(1,*,IOSTAT=ISO) TSATM(NS).TIM(M),(TSATM(NS).VAL(M,I),I=1,7)
        
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE ASER.INP'
      ENDDO
      DO M=1,TSATM(NS).NREC
        TSATM(NS).TIM(M)  = TSATM(NS).TIM(M)+TAASER(NS)
        TSATM(NS).VAL(M,4)= RAINCVT*TSATM(NS).VAL(M,4)
        TSATM(NS).VAL(M,5)= EVAPCVT*TSATM(NS).VAL(M,5)
        TSATM(NS).VAL(M,6)= SOLRCVT*TSATM(NS).VAL(M,6)
        TSATM(NS).VAL(M,7)=  CLDCVT*TSATM(NS).VAL(M,7)
      ENDDO
    ENDDO
    CLOSE(1)
    
    ! *** HTBED2 = Bottom Heat Ex Coeff (W / m2 / deg C)
    ! *** RHO    = 1000.0  Density (kg / m^3)
    ! *** CP     = 4179.0  Specific Heat (J / kg / degC)
    ! *** 0.2393E-6 = 1/RHO/CP
    HTBED2 = HTBED2*0.2393E-6    ! *** m/s
  ENDIF
  
  IF( NASER > 1 )THEN 
    WRITE(*,'(A)')'  READING ATMMAP.INP'
    OPEN(1,FILE='atmmap.inp',STATUS='UNKNOWN')
    STR=READSTR(1)    
    READ(1,*) NATMMAP

    DO NA=1,NATMMAP
      READ(1,*)TATMMAPBEG(NA),TATMMAPEND(NA)    
      STR=READSTR(1) 
      DO L=2,LA
        READ(1,*) ID,JD,(ATMWHT(N,L,NA),N=1,NASER)
      ENDDO    
    ENDDO
    CLOSE(1)
  ENDIF
  TEMB(1)=ABS(TBEDIT)
  TEMB(LC)=TEMB(1)
  TEMB1(1)=TEMB(1)
  TEMB1(LC)=TEMB(1)
  IF( ISRESTI == 0 )THEN
    DO L=2,LA
      TEMB(L)=TEMB(1)
      TEMB1(L)=TEMB(1)
    ENDDO
  ENDIF
                                                                  
  ! *** READ IN ABOVE WATER SURFACE WIND TIME SERIES FROM THE                                                             
  ! *** FILE WSER.INP                                                                                                     
  DO L=2,LA
    WINDST(L)=0.
    TSX(L)=0.
    TSY(L)=0.
  ENDDO
  
  IF( NWSER > 0 )THEN
    WRITE(*,'(A)')'READING WSER.INP'
    OPEN(1,FILE='wser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** PERIOD TO TURN OFF WIND SHEAR (DEPRECATED IN EE7.3 SINCE ADDITION OF ISICE)
    WRITE(*,'(A,I5)')  '  NUMBER OF WIND SERIES=',NWSER

    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    CALL SKIPCOM(1,'*')
    
    ! *** LOOP OVER EACH TIME SERIES
    STR=READSTR(1)
    DO NS=1,NWSER
      READ(1,*,IOSTAT=ISO) M,TCWSER(NS),TAWSER(NS),WINDSCT,ISWDINT(NS),WINDH(NS)
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE WSER.INP'
      
      IF( WINDH(NS) <= 0.1)WINDH(NS) = 2.0     ! *** FIXING THE WIND SPEED MEASUREMENT HEIGHT (I.E. USE THE WIND SPEED AS ENTERED)
      WRITE(*,'(A,I5,I10,F6.1)')'  SERIES, NPTS, ANEMOMETER HEIGHT (m)=',NS,TSWND(NS).NREC,WINDH(NS)

      DO M=1,TSWND(NS).NREC
        READ(1,*,IOSTAT=ISO) TSWND(NS).TIM(M),(TSWND(NS).VAL(M,I),I=1,2) !TWSER(M,NS),WINDS(M,NS),WINDD(M,NS)   
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE WSER.INP'
      ENDDO
      
      ! *** APPLY FACTORS AND OFFSETS AND COMPLETE TSWND DATA STRUCTURE
      DO M=1,TSWND(NS).NREC
        TSWND(NS).TIM(M)=TSWND(NS).TIM(M)+TAWSER(NS)
      ENDDO
      IF( ISWDINT(NS) <= 1 )THEN
        DO M=1,TSWND(NS).NREC
          TSWND(NS).VAL(M,1)=WINDSCT*TSWND(NS).VAL(M,1)
        ENDDO
      ENDIF
      IF( ISWDINT(NS) == 1 )THEN
        DO M=1,TSWND(NS).NREC
          IF( TSWND(NS).VAL(M,2) <= 180. )THEN
            TSWND(NS).VAL(M,2)=TSWND(NS).VAL(M,2)+180.
            IF( TSWND(NS).VAL(M,2) == 360.) TSWND(NS).VAL(M,2)=0.
          ELSE
            TSWND(NS).VAL(M,2)=TSWND(NS).VAL(M,2)-180.
            IF( TSWND(NS).VAL(M,2) == 360.) TSWND(NS).VAL(M,2)=0.
          ENDIF
        ENDDO
      ENDIF
      IF( ISWDINT(NS) == 2 )THEN
        DO M=1,TSWND(NS).NREC
          TSWND(NS).VAL(M,1)=WINDSCT*TSWND(NS).VAL(M,1)
          TSWND(NS).VAL(M,2)=WINDSCT*TSWND(NS).VAL(M,2)
        ENDDO
      ENDIF
    ENDDO
    CLOSE(1)
  ENDIF
  
  IF( NWSER > 1 )THEN
    WRITE(*,'(A)')'  READING WNDMAP.INP'
    OPEN(1,FILE='wndmap.inp',STATUS='UNKNOWN')    
    STR=READSTR(1)    
    READ(1,*) NWNDMAP

    DO NW=1,NWNDMAP
      READ(1,*)TWNDMAPBEG(NW),TWNDMAPEND(NW)       
      STR=READSTR(1)
      DO L=2,LA
        READ(1,*) ID,JD,(WNDWHT(N,L,NW),N=1,NWSER)
      ENDDO
    ENDDO
    CLOSE(1)
  ENDIF

  !**********************************************************************C                                                
  ! **  READ IN EXTERNALLY SPECIFIED ICE COVER INFORMATION                                                                
  ! **  FROM THE FILE ISER.INP                                                                                        
  !----------------------------------------------------------------------C                                                
  IF( ISICE == 1 .AND. NISER >= 1 )THEN
    WRITE(*,'(A)')'READING ISER.INP'
    OPEN(1,FILE='iser.inp')
    STR=READSTR(1)
    
    DO NS=1,NISER
      RICECOVT(NS)=0.                                                                                                     
      RICETHKT(NS)=0.                                                                                                     
      MITLAST(NS) =2

      READ(1,*,IOSTAT=ISO) M,TCISER(NS),TAISER(NS),RMULADJCOV,RMULADJTHK
      IF( ISO > 0 ) STOP 'ISER.INP: READING ERROR'

      DO M=1,TSICE(NS).NREC                                                                                                    
        READ(1,*,IOSTAT=ISO) TSICE(NS).TIM(M),TSICE(NS).VAL(M,1) !TISER(M,NS),RICETHKS(M,NS)
        
        IF( TSICE(NS).VAL(M,1) >= MINICETHICK )THEN
          TSICE(NS).VAL(M,2) = 1.0
        ELSE
          TSICE(NS).VAL(M,2) = 0.0
        ENDIF
        IF( ISO > 0 ) STOP 'ISER.INP: READING ERROR'
      ENDDO

      DO M=1,TSICE(NS).NREC                                                                                                    
        TSICE(NS).TIM(M)=TSICE(NS).TIM(M)+TAISER(NS)
        TSICE(NS).VAL(M,2)=RMULADJCOV*TSICE(NS).VAL(M,2)                                                                            
        TSICE(NS).VAL(M,1)=RMULADJTHK*TSICE(NS).VAL(M,1)                                                                         
      ENDDO

    ENDDO

    CLOSE(1)
    
  ELSEIF( ISICE == 2 )THEN
    WRITE(*,'(A)')'READING ISTAT.INP'
    OPEN(1,FILE='istat.inp')
    STR=READSTR(1)
    NS = 1    ! ** ALWAYS NISER = 1
    MITLAST(NS) =2
    READ(1,*,IOSTAT=ISO) M,TCISER(NS),TAISER(NS)
    IF( ISO > 0 ) STOP 'ISTAT.INP: READING ERROR'
    DO M=1,TSICE(NS).NREC
      READ(1,*,IOSTAT=ISO)TSICE(NS).TIM(M),TSICE(NS).VAL(M,2)  ! ** ICE COVER: 0/1
      IF( ISO > 0 ) STOP 'ISTAT.INP: READING ERROR'
      TSICE(NS).TIM(M) = TSICE(NS).TIM(M)+TAISER(NS)
      TSICE(NS).VAL(M,1) = RICETHK0
    ENDDO
  ENDIF                                                                                                                  

  !----------------------------------------------------------------------C
  IF( ISICE == 1 .AND. NISER > 1 )THEN
    WRITE(*,'(A)')'READING ICEMAP.INP'
    OPEN(1,FILE='icemap.inp')
    STR=READSTR(1)
    READ(1,*) NICEMAP

    DO NI=1,NICEMAP
      READ(1,*)TICEMAPBEG(NI),TICEMAPEND(NI)   
      STR=READSTR(1)
      DO L=2,LA
        READ(1,*) ID,JD,(RICEWHT(NI,L,N),N=1,NISER)
      ENDDO
    ENDDO
    CLOSE(1)
  ENDIF

  ! *** End of ICE *******************************                                                                        

  ! *** READ IN SHELL FISH LARAVE BEHAVIOR DATA                                                                           
  ! *** FROM THE FILE SFBSER.INP                                                                                          
  IF( ISTRAN(4) >= 1 )THEN
    WRITE(*,'(A)')'READING SFBSER.INP'
    OPEN(1,FILE='sfbser.inp',STATUS='UNKNOWN')
                                                                  
  ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    
    READ(1,*,IOSTAT=ISO) MSFSER,TCSFSER,TASFSER,TSRSF,TSSSF,ISSFLDN,ISSFLFE,SFLKILL
    IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFBSER.INP'
    DO M=1,MSFSER
      READ(1,*,IOSTAT=ISO) TSFSER(M),RKDSFL(M),WSFLST(M),WSFLSM(M),DSFLMN(M),DSFLMX(M),SFNTBE(M),SFATBT(M)
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE SFBSER.INP'
    ENDDO
    CLOSE(1)
  ENDIF
                                                                  
  ! *** READ VEGETATION DATA FROM VEGE.INP AND VEGSER.INP                                                                 
                                                                  
  IF( ISVEG >= 1 )THEN
    WRITE(*,'(A)')'READING VEGE.INP'
    OPEN(1,FILE='vege.inp',STATUS='UNKNOWN')
    STR=READSTR(1)
    
    READ(1,*)MVEGTYP,MVEGOW,NVEGSER,UVEGSCL
    DO M=1,MVEGTYP
     READ(1,*,ERR=3120)IDUM,NVEGSERV(M),RDLPSQ(M),BPVEG(M),HPVEG(M),ALPVEG(M),BETVEG(M),GAMVEG(M),SCVEG(M)
      BDLTMP=BPVEG(M)*BPVEG(M)*RDLPSQ(M)
      PVEGX(M)=1.-BETVEG(M)*BDLTMP
      PVEGY(M)=1.-BETVEG(M)*BDLTMP
      PVEGZ(M)=1.-ALPVEG(M)*BDLTMP
      BDLPSQ(M)=BPVEG(M)*RDLPSQ(M)
    ENDDO
    CLOSE(1)
    IF( NVEGSER > 0 )THEN
      DO M=1,NVEGSER
        MVEGTLAST(M)=1
      ENDDO
    ENDIF

    !!!BEGIN SCJ BLOCK
    IF( LMHK )THEN !MHK devices
      WRITE(*,'(A)')'READING MHK.INP'
      OPEN(1,FILE='mhk.inp',STATUS='UNKNOWN')

      STR=READSTR(1)
      
      READ(1,*,ERR=3122)MHKTYP,NFLAGPWR,UPSTREAM,OUTPUTFLAG

      ALLOCATE(BOFFMHK(MHKTYP),BOFFSUP(MHKTYP))
      ALLOCATE(TOFFMHK(MHKTYP),TOFFSUP(MHKTYP))

      IF( NFLAGPWR == 1 )THEN
        DO M=1,MHKTYP
          READ(1,*,ERR=3122)WIDTHMHK(M),WIDTHSUP(M), BOFFMHK(M),BOFFSUP(M),TOFFMHK(M),TOFFSUP(M),CTMHK(M),CDSUP(M),VMINCUT(M),VMAXCUT(M),DENMHK(M)
          CTMHK(M)=CTMHK(M)*DENMHK(M)
          CDSUP(M)=CDSUP(M)*DENMHK(M)

          DO L=2,LA
            IF( M+90 == MVEGL(L) )THEN
              ZMINMHK(M,L)=BELV(L)+BOFFMHK(M)
              ZMAXMHK(M,L)=BELV(L)+TOFFMHK(M)
              ZMINSUP(M,L)=BELV(L)+BOFFSUP(M)
              ZMAXSUP(M,L)=BELV(L)+TOFFSUP(M)
              DIAMMHK=ZMAXMHK(M,L)-ZMINMHK(M,L)
              IF( DIAMMHK<0.0 )THEN   !error check
                PRINT*,'MHK ZMIN > ZMAX'
                STOP
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        READ(1,*) !skip the header line
        READ(1,*)BETAMHK_D,BETAMHK_P,CE4MHK,PB_COEF

      ELSEIF( NFLAGPWR == 2 )THEN
        PRINT*,'Not available yet'

      ELSEIF( NFLAGPWR == 3 )THEN !FFP input style
        READ(1,*,ERR=3122)WIDTHMHK(M),WIDTHSUP(M),BOFFMHK(M),HEIGHTMHK(M),HEIGHTSUP(M),REFELEV(M),CTMHK(M),CDSUP(M),VMINCUT(M),VMAXCUT(M),DENMHK(M)
	CTMHK(M)=CTMHK(M)*DENMHK(M)
	CDSUP(M)=CDSUP(M)*DENMHK(M)
	DO L=2,LA
	  IF( M+90 == MVEGL(L) )THEN
	    ZMINMHK(M,L)=BELV(L)+REFELEV(M)
	    ZMAXMHK(M,L)=BELV(L)+REFELEV(M)+HEIGHTMHK(M)
	    ZMINSUP(M,L)=BELV(L)
	    ZMAXSUP(M,L)=BELV(L)+REFELEV(M)+HEIGHTSUP(M)
            DIAMMHK=ZMAXMHK(M,L)-ZMINMHK(M,L)
            IF( DIAMMHK<0.0 )THEN
              PRINT*,'MHK ZMIN > ZMAX'
              STOP
            ENDIF
          ENDIF
        ENDDO

        READ(1,*) !skip the header line
        READ(1,*)BETAMHK_D,BETAMHK_P,CE4MHK,PB_COEF

      ENDIF
      CLOSE(1)

      DO L=2,LA
        IF( MVEGL(L)>90 )THEN
          IF( (DIAMMHK>DXP(L) .OR. DIAMMHK>DYP(L)) .AND. DENMHK(MVEGL(L)-90)>1.0 )THEN
            PRINT*,'MHK DIAMETER EXCEEDS CELL SIZE'
            PRINT*,'AND DENSITY >= 1'
            STOP
          ENDIF
        ENDIF
      ENDDO
    ENDIF
    !!!END SCJ BLOCK

  ENDIF

  GOTO 3124
   3120 WRITE(6,3121)                                                                                                     
   3121 FORMAT('  READ ERROR FOR FILE VEGE.INP ')                                                                         
  STOP

   3122 WRITE(6,3123)                                                                                                     
   3123 FORMAT('  READ ERROR FOR FILE MHK.INP ')                                                                          
  STOP

   3124 CONTINUE                                                                                                          
  IF( NVEGSER >= 1 )THEN
    WRITE(*,'(A)')'READING VEGSER.INP'
    OPEN(1,FILE='vegser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)    
    DO NS=1,NVEGSER
      READ(1,*,IOSTAT=ISO) MVEGSER(NS),TCVEGSER(NS),TAVEGSER(NS)
      IF( ISO > 0 ) STOP '  READ ERROR FOR FILE VEGSER.INP'
      DO M=1,MVEGSER(NS)
        READ(1,*,IOSTAT=ISO)TVEGSER(M,NS),VEGSERR(M,NS),VEGSERB(M,NS),VEGSERH(M,NS)
        IF( ISO > 0 ) STOP '  READ ERROR FOR FILE VEGSER.INP'
        TVEGSER(M,NS)=TVEGSER(M,NS)+TAVEGSER(NS)
      ENDDO
    ENDDO
    CLOSE(1)
                                                                  
    ! *** REINITIALIZE CLASSES HAVING TIME SERIES INFORMATION                                                               
    DO M=1,MVEGTYP
      IF( NVEGSERV(M) > 0 )THEN
        NS=NVEGSERV(M)
        RDLPSQ(M)=VEGSERR(1,NS)
        BPVEG(M)=VEGSERB(1,NS)
        HPVEG(M)=VEGSERH(1,NS)
        BDLTMP=BPVEG(M)*BPVEG(M)*RDLPSQ(M)
        PVEGX(M)=1.-BETVEG(M)*BDLTMP
        PVEGY(M)=1.-BETVEG(M)*BDLTMP
        PVEGZ(M)=1.-ALPVEG(M)*BDLTMP
        BDLPSQ(M)=BPVEG(M)*RDLPSQ(M)
      ENDIF
    ENDDO
  ENDIF
  GOTO 7122
   7120 WRITE(6,7121)                                                                                                     
   7121 FORMAT('  READ ERROR FOR FILE VEGSER.INP ')                                                                       
  STOP
   7122 CONTINUE                                                                                                          
                                                                  
  !**********************************************************************C                                                
  ! *** READ BANK EROSION MAP AND TIME SERIES FILE                                                                        
  IF( (ISTRAN(6)>0 .OR. ISTRAN(7) > 0 ) .AND. ISBKERO >= 1 )THEN
                                                                  
    WRITE(*,'(A)')'READING BEMAP.INP'
    OPEN(1,FILE='bemap.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)                                                                  
    READ(1,*)NBEPAIR,NBESER
                                                                  
    DO NS=1,NBEPAIR
      READ(1,*)IBANKBE(NS),JBANKBE(NS),ICHANBE(NS),JCHANBE(NS),NBESERN(NS),FBESER(NS)
    ENDDO
                                                                  
    CLOSE(1)
                                                                  
  ENDIF
                                                                  
                                                                  
  IF( (ISTRAN(6)>0 .OR. ISTRAN(7) > 0 ) .AND. ISBKERO >= 1 .AND. NBESER > 0 )THEN
                                                                  
    WRITE(*,'(A)')'READING BESER.INP'
    OPEN(1,FILE='beser.inp',STATUS='UNKNOWN')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    
    DO NS=1,NBESER
                                                                  
      READ(1,*)MBESER(NS),TCBESER(NS),TABESER(NS),RMULADJ,ADDADJ
      MBETLAST(NS)=1
                                                                  
      DO M=1,MBESER(NS)
        READ(1,*)TBESER(M,NS),BESER(M,NS),FWCBESER(M,NS)
        TBESER(M,NS)=TBESER(M,NS)+TABESER(NS)
        BESER(M,NS)=RMULADJ*(BESER(M,NS)+ADDADJ)
      ENDDO
                                                                  
    ENDDO
                                                                  
    CLOSE(1)
                                                                  
  ENDIF
                                                                  
  !**********************************************************************C                                                
  ! *** READ ZONALLY VARYING SEDIMENT BED PARTICLE MIXING                                                                 
  !----------------------------------------------------------------------C                                                
  ITMPPMX=0
  DO NT=1,NTOX
    IF( ISPMXZ(NT) == 1 )ITMPPMX=1
  ENDDO
                                                                  
  IF( ISTRAN(5) > 0 .AND. ITMPPMX == 1 )THEN
                                                                  
    WRITE(*,'(A)')'READING PARTMIX.INP'
    OPEN(1,FILE='partmix.inp')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    
    !#######################################################################                                                
    !     HQI change to input multiplication scale factor for particle mixing rate                                          
    !     RM 10/06/05                                                                                                       
    READ(1,*)NPMXZ,NPMXPTS,PMIXSF

    !#######################################################################                                                
    DO NZ=1,NPMXZ
      DO NP=1,NPMXPTS
          READ(1,*)PMXDEPTH(NP,NZ),PMXCOEF(NP,NZ)
          !#######################################################################                                                
          !     HQI change to input multiplication scale factor for particle mixing rate                                          
          !     RM 10/06/05                                                                                                       
          PMXCOEF(NP,NZ) = PMIXSF*PMXCOEF(NP,NZ)
          !#######################################################################                                                
      ENDDO
    ENDDO
                                                                  
    CLOSE(1)
                                                                  
                                                                  
    WRITE(*,'(A)')'READING PMXMAP.INP'
    OPEN(1,FILE='pmxmap.inp')
                                                                  
    ! *** SKIP OVER TITLE AND AND HEADER LINES                                                                              
    STR=READSTR(1)
    
    DO L=2,LC-1
       READ(1,*)LDUM,IDUM,JDUM,LPMXZ(L)
    ENDDO
                                                                  
    CLOSE(1)
                                                                  
  ENDIF
                                                                                                            
END

! ***************************************************************************************
! *** DS-INTL UTILITIES
FUNCTION PARSE_REAL(INLINE)
  USE GLOBAL,ONLY: IK4
  IMPLICIT NONE
  INTEGER(IK4) :: I1,I2,ILEN,IPOS,ILEN2
  REAL :: PARSE_REAL
  CHARACTER*(*) INLINE
  CHARACTER*15  CVAL

  ILEN=LEN_TRIM(INLINE)
  PARSE_REAL=0.
  DO I1=1,ILEN
    IF( INLINE(I1:I1) == ':' )THEN
      DO IPOS=I1+1,ILEN
        IF( INLINE(IPOS:IPOS)/=' ')EXIT
      ENDDO
      IF( IPOS>ILEN)RETURN
      CVAL=INLINE(IPOS:ILEN)
      ILEN2=LEN_TRIM(CVAL)
      DO I2=1,ILEN2
        IF( CVAL(I2:I2) == ' ' .OR. CVAL(I2:I2) == ',' .OR. I2 == ILEN2 )THEN
          READ(CVAL(1:I2),'(F12.1)',ERR=999)PARSE_REAL
          RETURN
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  999 STOP ' ERROR PARSING REAL'                                                                                  
END FUNCTION

FUNCTION PARSE_LOGICAL(INLINE)

  USE GLOBAL,ONLY: IKV
  IMPLICIT NONE

  INTEGER :: ILEN,IC,JC,IPOS,ILEN2
  CHARACTER*(*) INLINE
  CHARACTER*12  CVAL
  LOGICAL       PARSE_LOGICAL

  ILEN=LEN_TRIM(INLINE)
  DO IC=1,ILEN
    IF( INLINE(IC:IC) == ':' )THEN
      DO IPOS=IC+1,ILEN
        IF( INLINE(IPOS:IPOS)/=' ')EXIT
      ENDDO
      IF( IPOS>ILEN)RETURN
      CVAL=INLINE(IPOS:ILEN)
      ILEN2=LEN_TRIM(CVAL)
      DO JC=1,ILEN2
        IF( CVAL(JC:JC) == ' ' .OR. CVAL(JC:JC) == ',' .OR. JC == ILEN2 )THEN
          IF( CVAL(1:1) == 'T' .OR. CVAL(1:1) == 'Y' )THEN
            PARSE_LOGICAL=.TRUE.
            RETURN
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  PARSE_LOGICAL=.FALSE.
  
END FUNCTION
