! ******************************************************************************************
! *** Writes the hyd file for WINASP/WASP-HYDRO water quality model.
Subroutine WASP7HYDRO

  ! *** ----- Previous subroutine name -- SUBROUTINE WASP7EPA
  ! *** ------------------------------------------------------------------------------
  ! ***  PURPOSE:
  ! *** 
  ! ***    Subroutine WASPHYDRO writes the hyd file for WINASP water quality model.
  ! *** 
  ! ***  MODIFICATION HISTORY:
  ! *** 
  ! ***    Date       Author         Comments
  ! ***    ---------- -------------- -------------------------------------------------
  ! ***    2013-05-08 Paul M. Craig  Converted to f90 and cleaned code 
  ! ***     5/23/2006 Hugo Rodriguez Add a limit maximum dispersion specific to certain cells
  ! ***     4/11/2006 Hugo Rodriguez Add qfactor
  ! ***    10/27/2005 Hugo Rodriguez Calculate the initial date of wasp instead of reading it
  ! ***    06/24/2005 Hugo Rodriguez Limit Ab to Abwmax
  ! ***    06/15/2005 Hugo Rodriguez Zero out AH coefficient transferred to wasp IF
  ! ***                              ISHDMF<2
  ! ***    05/04/2005 Hugo Rodriguez Limit maximum dispersion to cell volume/time step
  ! ***    05/02/2005 Hugo Rodriguez Eliminate dispersion for boundary flows.
  ! ***    30/10/2003 Hugo Rodriguez This version uses a dll to transfer the hydrodynamic
  ! ***               Tim Wool       data, including salinity and temperature, from EFDC
  ! ***                              to WASP version 6.1.
  ! *** ------------------------------------------------------------------------------
  ! *** 
  USE GLOBAL
  USE IFPORT

#ifdef WASPOUT
  IMPLICIT NONE

  INTEGER(IK4) :: LT,L,IYEAR,IMON,IDAY,I,IS,LCLTM2,LTYPE,KWASP,LBELOW,LAUX,K,IONE
  INTEGER(IK4) :: LWSPTMP,IZERO,IM1,IERROR,NJUN,IOS,KK,NCHNH,NCHNV,NCHN,NQ,NS,KMUL,LDTM,LUTM
  INTEGER(IK4) :: LSLT,LN,KMUL1,KMUL2,LW,LS,IPTMP,JMTMP,JPTMP,LCELTMP

  INTEGER(IK4), SAVE :: IBEGIN,IDAYS,nsg,nf

  REAL(RK4)    :: SVPT,SCALR,WSS1,WSS2,WSS3,VOLUME,DXYSUM,VELX,VELY,VELZ,VELMAG,ADDLW,ADDLS,FLOWX
  REAL(RK4)    :: UDDXTMP,ADDL,ADDL1,TMPVAL,FLOWY,VDDYTMP,FLOWZ,WDDZTMP,VOLUM,DEPTH

  !INCLUDE 'EFDC_WASPHYDRO.CMN'
  CHARACTER(80)    :: errstring
  LOGICAL(4) RES

  ! *** DECLARE PERSISTENT VARIABLES
  REAL(RK4),    SAVE    :: AD,ADCOEFF,ABWMAX
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: ABWMX  
  CHARACTER*6,SAVE,ALLOCATABLE,DIMENSION(:) :: SSN          !  SSN(lcm)
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: crnu        !  crnu(nf),brintt(nf),flow(nf)  
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: brintt
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: flow
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: SegVolume   !  SegVolume(nsg),SegDepth(nsg),SegVel(nsg),SegSalt(nsg),SegTemp(nsg)
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: SegDepth  
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: SegVel  
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: SegSalt  
  REAL(RK4),    SAVE,ALLOCATABLE,DIMENSION(:) :: SegTemp  

  ! *** BEGIN OF OLD INCLUDE HydroLink_Set.INT
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetdebug(hl_debug)
      !ms$attributes c,dllimport,alias:'__hlsetdebug' :: hlsetdebug
      Integer*4 hl_debug
      !ms$attributes reference :: hl_debug
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlgetlasterror(message)
      !ms$attributes c,dllimport,alias:'__hlgetlasterror' :: hlgetlasterror
      character*(*) message
      !ms$attributes reference :: message
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlopen(FName, hl_mode, hl_handle, ierror)
      !ms$attributes c,dllimport,alias:'__hlopen' :: hlopen
      character*(*) Fname
      Integer*4 hl_mode, hl_handle, ierror
      !ms$attributes reference :: FName, hl_mode, hl_handle, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetlanguage(hl_handle, hl_language, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetlanguage' :: hlsetlanguage
      Integer*4 hl_handle, hl_language, ierror
      !ms$attributes reference :: hl_handle, hl_language, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetcreator(hl_handle, hl_creator, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetcreator' :: hlsetcreator
      Integer*4 hl_handle, hl_creator, ierror
      !ms$attributes reference :: hl_handle, hl_creator, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hladddescription(hl_handle, id, string, ierror)
      !ms$attributes c,dllimport,alias:'__hladddescription' :: hladddescription
      Integer*4 hl_handle, id, ierror
      character *(*) string
      !ms$attributes reference :: hl_handle, id, string, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetseedmoment(hl_handle, month, day, year, hour, minute, second, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetseedmoment' :: hlsetseedmoment
      Integer*4 hl_handle, month, day, year, hour, minute,second,ierror
      !ms$attributes reference :: hl_handle, month, day, year, hour, minute, second, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetnumsegments(hl_handle, numsegs, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetnumsegments' :: hlsetnumsegments
      Integer*4 hl_handle, numsegs, ierror
      !ms$attributes reference :: hl_handle, numsegs,ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetsegname(hl_handle,index, segname, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetsegname' :: hlsetsegname
      Integer*4 hl_handle, index, ierror
      Character *(*) segname
      !ms$attributes reference :: hl_handle, index, segname,ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetnumflowpaths(hl_handle, numfp, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetnumflowpaths' :: hlsetnumflowpaths
      Integer*4 hl_handle, numfp, ierror
      !ms$attributes reference :: hl_handle, numfp, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetnumsegconsts(hl_handle, numsc, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetnumsegconsts' :: hlsetnumsegconsts
      Integer*4 hl_handle, numsc, ierror
      !ms$attributes reference :: hl_handle, numsc, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetnumfpconsts(hl_handle, numfpc, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetnumfpconsts' :: hlsetnumfpconsts
      Integer*4 hl_handle, numfpc, ierror
      !ms$attributes reference :: hl_handle, numfpc, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetsegconsttype(hl_handle, sc_index, sc_type, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetsegconsttype' :: hlsetsegconsttype
      Integer*4 hl_handle, sc_index, sc_type, ierror
      !ms$attributes reference :: hl_handle, sc_index, sc_type, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetfpconsttype(hl_handle, fp_index, fp_type, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetfpconsttype' :: hlsetfpconsttype
      Integer*4 hl_handle, fp_index, fp_type, ierror
      !ms$attributes reference :: hl_handle, fp_index, fp_type, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetvartimestep(hl_handle, vardt, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetvartimestep' :: hlsetvartimestep
      Integer*4 hl_handle, vardt, ierror
      !ms$attributes reference :: hl_handle, vardt, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsethydtimestep(hl_handle, timestep, ierror)
      !ms$attributes c,dllimport,alias:'__hlsethydtimestep' :: hlsethydtimestep
      Integer*4 hl_handle, ierror
      Real*4 timestep
      !ms$attributes reference :: hl_handle, timestep, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetupdateint(hl_handle, updateinterval, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetupdateint' :: hlsetupdateint
      Integer*4 hl_handle, ierror
      Real*4 updateinterval
      !ms$attributes reference :: hl_handle,updateinterval, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsethydtowaspratio(hl_handle, iratio, ierror)
      !ms$attributes c,dllimport,alias:'__hlsethydtowaspratio' :: hlsethydtowaspratio
      Integer*4 hl_handle, iratio, ierror
      !ms$attributes reference :: hl_handle, iratio, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetnumlayers(hl_handle, numlayers, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetnumlayers' :: hlsetnumlayers
      Integer*4 hl_handle, numlayers, ierror
      !ms$attributes reference :: hl_handle, numlayers, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetflowpath(hl_handle, flow_index, from_seg, to_seg, direction, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetflowpath' :: hlsetflowpath
      Integer*4 hl_handle, flow_index,from_seg, to_seg,direction,ierror
      !ms$attributes reference :: hl_handle, flow_index, from_seg, to_seg, direction, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsetflowinfo(hl_handle, index, value, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetflowinfo' :: hlsetflowinfo
      Integer*4 hl_handle, index, ierror
      Real*4 value(31000)
      !ms$attributes reference :: hl_handle, index, value, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
  subroutine hlsetseginfo(hl_handle,index, value, ierror)
      !ms$attributes c,dllimport,alias:'__hlsetseginfo' :: hlsetseginfo
      Integer*4 hl_handle, index, ierror
      Real*4 value(12000)
      !ms$attributes reference :: hl_handle, index, value, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlsettimestep(hl_handle,value,ierror)
      !ms$attributes c,dllimport,alias:'__hlsettimestep' :: hlsettimestep
      Integer*4 hl_handle, ierror
      Real*4 value
      !ms$attributes reference :: hl_handle,value, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlmomentcomplete(hl_handle, ierror)
      !ms$attributes c,dllimport,alias:'__hlmomentcomplete' :: hlmomentcomplete
      Integer*4 hl_handle, ierror
      !ms$attributes reference :: hl_handle, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  interface
    subroutine hlclose(hl_handle, ierror)
      !ms$attributes c,dllimport,alias:'__hlclose' :: hlclose
      Integer*4 hl_handle, ierror
      !ms$attributes reference :: hl_handle, ierror
    END subroutine
  END interface
  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == ==  == =
  ! *** END OF OLD INCLUDE HydroLink_Set.INT

  !PARAMETER(nf=31000,nsg=12000)
  DIMENSION LAUX(ICM,JCM,KCM)

  INTEGER(IK4) LU,LD,IU,ID,JU,JD,NAUX
  INTEGER(IK4) istartyear,istartmonth,istartday,istarthour,istartminute,istartsecond,Ihl_debug,Ihl_mode,inumsegconsts,j
  !INTEGER FLAGWASPBC(NQSERM,KCM)

  CHARACTER*6  sn1
  CHARACTER*23 segname
  CHARACTER*17 SN2
  character*80 DESCRIPTION,MODELERNAME
  character*3 Itext, Jtext, Ktext
  character*80 RECHNR

  REAL(RKD)  AUX,thour,tmin,tsec
  REAL(RK4)  rinterval
  REAL(RK4)  vol1,vol2  

  ! *** 
  ! *** 
  ! *** **********************************************************************! *** 
  ! *** 
  ! ***  READ CONTROL DATA FILE EFDC.WSP
  ! *** 
  ! *** ----------------------------------------------------------------------!
  ! *** 

  ! *** ALLOCATE LOCAL ARRAYS
  IF(  .NOT. ALLOCATED(ABWMX) )THEN
    nsg = lcm*kcm
    nf = 31000
 
    ALLOCATE(ABWMX(nsg))
    ALLOCATE(SSN(LCM))
    ALLOCATE(crnu(nf))
    ALLOCATE(brintt(nf))
    ALLOCATE(flow(nf))
    ALLOCATE(SegVolume(nsg))
    ALLOCATE(SegDepth(nsg))
    ALLOCATE(SegVel(nsg))
    ALLOCATE(SegSalt(nsg))
    ALLOCATE(SegTemp(nsg))
  END IF

  SVPT=1.
  IF( NTSMMT < NTSPTC)SVPT=0.
  ! *** 
  IF( JSWASP == 1 )THEN

    ! ***       for jswasp=1 only first entry
    WRITE(*,'(A)')'READING EFDC.WSP'
    OPEN(1,FILE='EFDC.WSP',STATUS='UNKNOWN')
    WRITE(6,*)'EFDC.WSP opened'

    ! *** C1**  READ CELL VOLUME PARAMETERS (PMC - None of these paramters are being used.  Keep for legacy info)
    READ(1,1)
    READ(1,1)
    READ(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
    WRITE(*,*)'<EFDC.WSP1> ',IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP

    ! *** C2**  READ DIFFUSION PARAMETERS (PMC - The 1st 4 of these paramters are not being used.  Keep for legacy info)
    READ(1,1)
    READ(1,1)
    READ(1,*) NRFLD,SCALR,CONVR,ISNKH,ADCOEFF,ABWMAX        
    WRITE(*,*)'<EFDC.WSP2> ',NRFLD,SCALR,CONVR,ISNKH,ADCOEFF,ABWMAX        
    DO LT=2,LALT                                  
      l=lij(illt(lt),jllt(lt))                    
      ABWMX(l)=ABWMAX                             
    END DO                                        

    ! *** C3**  READ ADVECTION PARAMETERS (PMC - Only HYDFIL and IDAYS are being used.  Keep for legacy info)
    READ(1,1)
    READ(1,1)
    READ(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD,IDAYS                                      
    WRITE(*,*)'<EFDC.WSP3> ', IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD,IDAYS                                      

    ! *** C4**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)  (PMC - Only DEPSED is being used.  Keep for legacy info)
    READ(1,1)
    READ(1,1)
    READ(1,*) DEPSED,TDINTS,SEDIFF, WSS1, WSS2, WSS3
    WRITE(*,*)'<EFDC.WSP4> ',DEPSED,TDINTS,SEDIFF,WSS1,WSS2,WSS3

    ! *** C5**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
    READ(1,1)
    READ(1,1)
    READ(1,*) iyear,imon,iday
    WRITE(*,*)'<EFDC.WSP5> ',iyear,imon,iday

    do i=1,5                               
      read(1,*,err=11) RECHNR                    
      WRITE(*,*)'<EFDC.WSP6> ',RECHNR
    END do                                 
    11 continue       

    DO LT=2,LALT                           
      read(1,*,err=12)i,j,ABWMAX           
      l=lij(i,j)                           
      ABWMX(l)=ABWMAX                      
      WRITE(*,*)'<EFDC.WSP7> ',i,j,ABWMAX           
    END DO                                 
    12 continue                               

    CLOSE(1)
    WRITE(*,*)'EFDC.WSP read succesfully, now WRITE ABMax'

    ! *** CREATE THE WASP FOLDER
    RES = MAKEDIRQQ(OUTDIR//'wasp')

    OPEN(1,FILE=OUTDIR//'wasp\ABmax.txt',STATUS='UNKNOWN') 
    WRITE(1,*)'    I    J     ABmax'        
    DO LT=2,LALT                                  
      l=lij(illt(lt),jllt(lt))                    
      WRITE(1,21)illt(lt),jllt(lt),ABWMX(l)       
    END DO                                        
    close(1)
    21 format(2I5,f10.6)
    WRITE(6,*)'EFDC.WSP read succesfully and ABmax.txt written'
    1 FORMAT (80X)

    ! *** read qser file to check for flows only to some layers  (PMC - Not used for anything.  Keep for Legacy)
    IF( NQSER >= 1 )THEN                                            
      OPEN(1,FILE='QSER.INP',STATUS='UNKNOWN')                    
      ! ***  SKIP OVER TITLE AND AND HEADER LINES                          
      DO IS=1,14                                                  
        READ(1,1)                                                  
      ENDDO                                                       
    END IF                                                        

    CLOSE(1)                                                      
    goto 862                                                      
    860 WRITE(6,861)                                                  
    861 FORMAT('  READ ERROR FOR FILE QSER.INP ')                     
    862 continue                                                      

    ! *** **********************************************************************!

    ! ***  DEFINE EFDC-WASP CORRESPONDENCE AND INITIALIZE FILES

    ! *** **********************************************************************!
    ! ***  begin code inserted by Andy Stoddard 2-26-2008
    ! *** ------ Feb-26-2008 A Stoddard add code to WRITE WASP7SEG_EFDCIJK 
    ! *** ------ file listing of EFDC L,IJK cells and Wasp segment/layer data
    ! *** ----------------------------------------------------------------------!
    ! *** ----------------------------------------------------------------------!
    ! *** ---- 7/1/2005 A Stoddard added code to WRITE WASP7SEG_EFDCIJK.DAT file
    ! *** ---- Segment Name file is linked for import to WASP6 & WASP7
    OPEN(970,FILE=OUTDIR//'wasp\WASP7SEG_EFDCIJK.DAT',STATUS='UNKNOWN')  ! 7/1/2005 AS
    CLOSE(970,STATUS='DELETE')                                      ! 7/1/2005 AS 

    ! ***    OPEN(90,FILE='WASPP.OUT',STATUS='UNKNOWN')
    ! ***    OPEN(93,FILE='WASPC.OUT',STATUS='UNKNOWN')
    ! ***    CLOSE(90,STATUS='DELETE')
    ! ***    CLOSE(93,STATUS='DELETE')
    ! ***    OPEN(90,FILE='WASPP.OUT',STATUS='UNKNOWN')
    ! ***    OPEN(93,FILE='WASPC.OUT',STATUS='UNKNOWN')
    OPEN(970,FILE=OUTDIR//'wasp\WASP7SEG_EFDCIJK.DAT',STATUS='UNKNOWN')  ! 7/1/2005 AS

    WRITE(970,9702)              ! WRITE header records

    LCLTM2=LCLT-2
    LWASP=0
    IF( KC > 1 )THEN
      LTYPE=1
      KWASP=1
      DO LT=2,LALT
        LWASP=LWASP+1
        LBELOW=LWASP+LCLTM2
        I=ILLT(LT)
        J=JLLT(LT)
        L=LIJ(I,J)
        DMULT=HLPF(L)*DZC(L,KC)
        VOLUME=DXYP(L)*HLPF(L)*DZC(L,KC)
        IF( NTSMMT < NTSPTC )THEN
          DMULT=HP(L)*DZC(L,KC)
          VOLUME=DXYP(L)*HP(L)*DZC(L,KC)
        END IF
        LAUX(I,J,KC)=LWASP  

        WRITE(970,9701)L,I,J,KC,LWASP,KWASP         ! ,I,J,L,KC               ! 7/1/2005 AS
      END DO
      LTYPE=2
      DO K=KS,2,-1
        KWASP=KC-K+1
        DO LT=2,LALT
          LWASP=LWASP+1
          LBELOW=LWASP+LCLTM2
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          DMULT=HLPF(L)*DZC(L,K)
          VOLUME=DXYP(L)*HLPF(L)*DZC(L,K)
          IF( NTSMMT < NTSPTC )THEN
            DMULT=HP(L)*DZC(L,KC)
            VOLUME=DXYP(L)*HP(L)*DZC(L,KC)
          END IF
          LAUX(I,J,K)=LWASP 
          WRITE(970,9701)L,I,J,K,LWASP,KWASP         ! ,I,J,L,K               ! 7/1/2005 AS
        END DO
      END DO
    END IF
    LTYPE=2
    IF( KC == 1 ) LTYPE=1
    KWASP=KC
    DO LT=2,LALT
      LWASP=LWASP+1
      LBELOW=LWASP+LCLTM2
      I=ILLT(LT)
      J=JLLT(LT)
      L=LIJ(I,J)
      DMULT=HLPF(L)*DZC(L,KSZ(L))
      VOLUME=DXYP(L)*HLPF(L)*DZC(L,KSZ(L))
      IF( NTSMMT < NTSPTC )THEN
        DMULT=HP(L)*DZC(L,KC)
        VOLUME=DXYP(L)*HP(L)*DZC(L,KC)
      END IF
      IONE=1
      LAUX(I,J,1)=LWASP   
      WRITE(970,9701)L,I,J,IONE,LWASP,KWASP         ! ,I,J,L,IONE               ! 7/1/2005 AS
    END DO
    LTYPE=3
    KWASP=KC+1
    DXYSUM=0.
    LWSPTMP=LWASP+1
    DO LT=2,LALT
      LWSPTMP=LWSPTMP+1
    END DO
    ! ***  The following the lower benthic layer.  All upper benthic layer segments
    ! ***  have this layer immediately below them:
    DO LT=2,LALT
      LWASP=LWASP+1
      LBELOW=LWSPTMP
      I=ILLT(LT)
      J=JLLT(LT)
      L=LIJ(I,J)
      DXYSUM=DXYSUM+DXYP(L)
      VOLUME=DXYP(L)*DEPSED
      IZERO=0
    END DO
    ! ***  Next do the lower benthic layer:
    LTYPE=4
    KWASP=KC+2
    LWASP=LWASP+1
    LBELOW=0
    DMULT=DEPSED
    VOLUME=DXYSUM*DEPSED
    IM1=-1

    CLOSE(970)                           ! 7/1/2005

    ! ***  1001 FORMAT(6I5,2F10.4)
    ! ***  1031 FORMAT(2I5,F10.4,10X,A50)
    ! ***  1032 FORMAT(2F10.4)
    ! ***  1033 FORMAT(3I10,F10.1,4F10.3,'   !',4i5)
       

    ! ***  9701 FORMAT('I=',I3,' J=',I3,' L=',I4 ,' K=',I2,
    ! *** +       ' W',I4,' WK=',I2)

    9702 format('######################################################',/,  &
                '# File    = WASP7SEG_EFDCIJK.DAT',/,  &
                '# Field#1 = EFDC Cell L-Number',/,  &
                '# Field#2 = EFDC Cell I-Index',/,  &
                '# Field#3 = EFDC Cell J-Index',/,  &
                '# Field#4 = EFDC Cell K-Layer  [1=Bottom;KC=Surface]',/,  &
                '# Field#5 = WASP Segment Number',/,  &
                '# Field#6 = WASP Segment Layer [KC=Bottom;1=Surface]',/,  &
                '######################################################')     

  9701 FORMAT(I4,',', &            ! EFDC L-coordinate
              I4,',',     &        ! EFDC I-coordinate
              I5,',',     &        ! EFDC J-cell index 
              I2,',',     &        ! EFDC K-layer number   1=bottom;  KC=surface layer
              I5,',',     &        ! WASP Segment Number
              I2)                  ! WASP Layer Number     KC=bottom; 1= surface layer

    ! ***         WASPSEG.OUT data structure
    ! ***         Wasp7 wired for linkage of A30 character string file
    ! ***         W Wasp Segment Number
    ! ***         WK wasp Layer number 1=sfc   K=KC=bottom
    ! ***         I,J = EFDC I, J index 
    ! ***         L = EFDC grid cell number
    ! ***         K = EFDC layer K=KC sfc    K=1 bottom 

    ! ***  END code inserted by Andy Stoddard 2-26-2008
    ! *** **********************************************************************!

    ! *** ----------------------------------------------------------------------!
    ! *** 
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! ***  hlopen parameters
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    Ihl_handle = 0
    Ihl_debug = 1
    Ihl_mode = 1  !Ihl_mode=0 to READ from dll hyd file, =1 to WRITE to dll hyd file
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the debug flag 0=No debug 1=Debug (LOG.OUT)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetdebug(Ihl_debug)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Open the file
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlopen(HYDFIL,Ihl_mode,Ihl_handle,ierror)
    IF( ierror > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000)ierror,errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the language to FORTRAN
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetlanguage(Ihl_handle,1,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Store a description string
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    DESCRIPTION='   '//char(0)
    CALL hladddescription(Ihl_handle,0,DESCRIPTION,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Store the modeler name
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    MODELERNAME='Created by:  '//char(0)
    CALL hladddescription(Ihl_handle,1,MODELERNAME,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the creator
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetcreator(Ihl_handle, 1, ierror)   !1 for fortran program
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the seed moment (start date of the hyd file)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Calculate istartmonth, istartday,istartyear
    j=1461*(iyear+4800+int((imon-14)/12))/4+367*(imon-2-12*int((imon-14)/12))/12-(3*((iyear+4800+int((imon-14)/12)+100)/100))/4+iday-32075+tbegin+IDAYS
    istartyear=100*(int(4*(j+68569)/146097)-49)+int(4000*(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/   &
               1461001)+int(int(80*int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
    istartmonth=int(80*int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)+2-12*int(int(80*    &
                int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
    istartday=int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)-2447*int(80*int(  &
              int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)/80
    istarthour=0
    istartminute=0
    istartsecond=0

    ! *** Adjust hour, minutes and seconds for the start day IF necessary
    IBEGIN=IDAYS*NTSPTC
    AUX=FLOAT(IBEGIN)/FLOAT(NTSMMT)
    NAUX=INT(AUX)
    IF( AUX > NAUX )THEN
      AUX=(NAUX+1.-AUX)
      tsec=AUX*NTSMMT*TCON/NTSPTC
      istartsecond=INT(tsec)
      IF( tsec > 60 )THEN
        TMIN=tsec/60.0
        istartminute=INT(TMIN)
        istartsecond=INT((TMIN-istartminute)*60.)
        IF( IStartminute > 60 )THEN
          THOUR=FLOAT(istartminute)/60.0
          istarthour=INT(THOUR)
          istartminute=INT((THOUR-istarthour)*60.)
        ENDIF
      ENDIF
    ENDIF

    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetseedmoment(Ihl_handle,istartmonth,istartday,istartyear,istarthour,istartminute,istartsecond,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =

    CALL hlsetnumlayers(Ihl_handle,INT(kc,4),ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the number of segments
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    NJUN=KC*(LCLT-2)
    IF( njun > nsg )THEN
      WRITE(6,500)NJUN,NSG
      STOP
    END IF
  500 FORMAT('THE NUMBER OF WASP SEGMENTS IN YOUR APPLICATION',I6,1x,'IS GREATER THAN THE ARRAY DIMENSION:',I7)

    CALL hlsetnumsegments(Ihl_handle,njun,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! ***           Set WASP Segment Names (defaulft as I,J,K)
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    DO LT=2,LALT
      l=lij(illt(lt),jllt(lt))
      SSN(l)='      '//char(0)
    END DO
    OPEN(94,FILE='segname.inp',iostat=ios,STATUS='old')
    IF( ios == 0 )THEN
      do i=1,4
        read(94,*)
      END do
      do kk=1,la
        READ(94,*,err=111)I,J,SN1
        L=LIJ(I,J)
        SSN(l)=sn1
      END do
    END IF
    111 CLOSE(94)
    I=0
    DO K=KC,1,-1
      DO LT=2,LALT
        I=I+1
        LAUX(ILLT(LT),JLLT(LT),K)=I
        l=lij(illt(lt),jllt(lt))
        IF( il(l) >= 100 )THEN
          WRITE(itext,"(I3)")IL(L)
        ELSEIF( il(l) >= 10 )THEN
          WRITE(itext,"(I2)")IL(L)
        ELSE
          WRITE(itext,"(I1)")IL(L)
        END IF
        IF( jl(l) >= 100 )THEN
          WRITE(jtext,"(I3)")jL(L)
        ELSEIF( jl(l) >= 10 )THEN
          WRITE(jtext,"(I2)")jL(L)
        ELSE
          WRITE(jtext,"(I1)")jL(L)
        END IF
        IF( k >= 100 )THEN
          WRITE(ktext,"(I3)")k
        ELSEIF( k >= 10 )THEN
          WRITE(ktext,"(I2)")k
        ELSE
          WRITE(ktext,"(I1)")k
        END IF
        sn2=' I='//itext//' J='//jtext//' K='//ktext
        segname=sn2//SSN(l)//char(0)
        CALL hlsetsegname(ihl_handle,i,segname,ierror)
      END DO
    END DO
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the number of flow paths
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    NCHNH=0
    NCHNV=0
    DO LT=2,LALT
      I=ILLT(LT)
      J=JLLT(LT)
      L=LIJ(I,J)
      NCHNH=NCHNH+INT(SUBO(L))
      IF( IJCTLT(I+1,J) == 8 )THEN
        IF( SUBO(LEC(L)) == 1.) NCHNH=NCHNH+1
      END IF
      NCHNH=NCHNH+INT(SVBO(L))
      IF( IJCTLT(I,J+1) == 8 )THEN
        IF( SVBO(LNC(L)) == 1.) NCHNH=NCHNH+1
      END IF
      NCHNV=NCHNV+INT(SWB(L))
    END DO
    NCHN=KC*NCHNH+(KC-1)*NCHNV

    NQ=NQSIJ
    DO L=1,NQSIJ
      IF( LIJLT(IQS(L),JQS(L)) == 0 ) NQ=NQ-1
    END DO
    NCHN=NCHN+KC*NQ

    NQ=0                                    
    DO L=1,NQSIJ                            
      I=IQS(L)                           
      J=JQS(L)                           
      IF( LIJLT(I,J) == 0 ) GOTO 1001     
        NS=NQSERQ(L)                     
        DO K=1,KC                             
          !IF(flagwaspbC(ns,k) == 1 ) NQ=NQ+1   
        END DO                                
1001 END DO                                  
    NCHN=NCHN-NQ

    NQ=NQCTL
    DO L=1,NQCTL
      IF( LIJLT(IQCTLU(L),JQCTLU(L)) == 0 )THEN
        IF( LIJLT(IQCTLD(L),JQCTLD(L)) == 0 ) NQ=NQ-1
      END IF
    END DO
    NCHN=NCHN+KC*NQ
    IF( NCHN > NF )THEN
      WRITE(6,600)NCHN,NF
      STOP
    END IF
  600 FORMAT('THE NUMBER OF WASP FLOWS IN YOUR APPLICATION',I6,1X,'IS GREATER THAN THE ARRAY DIMENSION:',I7)

    CALL hlsetnumflowpaths(Ihl_handle, NCHN, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF

    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the flow path and direction
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    LCLTM2=LCLT-2
    LWASP=0
    DO K=KC,1,-1
      KMUL=KC-K
      DO LT=2,LALT
        I=ILLT(LT)
        J=JLLT(LT)
        L=LIJ(I,J)
        IF( SUBO(L) == 1. )THEN
          LWASP=LWASP+1
          LDTM=LT-1+KMUL*LCLTM2
          LUTM=LDTM-1
          IF( IJCTLT(I-1,J) == 8 ) LUTM=0
            CALL hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
            IF( ierror  > 0 )THEN
              CALL hlgetlasterror(errstring)
              WRITE(6,6000) ierror, errstring
              STOP
            END IF
        END IF

        IF( IJCTLT(I+1,J) == 8 )THEN
          IF( SUBO(LEC(L)) == 1. )THEN
            LWASP=LWASP+1
            LDTM=0
            LUTM=LT-1+KMUL*LCLTM2
            CALL hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
            IF( ierror  > 0 )THEN
              CALL hlgetlasterror(errstring)
              WRITE(6,6000) ierror, errstring
              STOP
            END IF
          END IF
        END IF
      END DO

      DO LT=2,LALT
        I=ILLT(LT)
        J=JLLT(LT)
        L=LIJ(I,J)
        IF( SVBO(L) == 1. )THEN
          LWASP=LWASP+1
          LSLT=LSCLT(LT)
          LDTM=LT-1+KMUL*LCLTM2
          LUTM=LSLT-1+KMUL*LCLTM2
          IF( IJCTLT(I,J-1) == 8 ) LUTM=0
            CALL hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,2,ierror)
            IF( ierror  > 0 )THEN
              CALL hlgetlasterror(errstring)
              WRITE(6,6000) ierror, errstring
              STOP
            END IF
        END IF
        IF( IJCTLT(I,J+1) == 8 )THEN
          LN=LNC(L)
          IF( SVBO(LN) == 1. )THEN
            LWASP=LWASP+1
            LSLT=LSCLT(LT)
            LDTM=0
            LUTM=LT-1+KMUL*LCLTM2
            CALL hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,2,ierror)
            IF( ierror  > 0 )THEN
              CALL hlgetlasterror(errstring)
              WRITE(6,6000) ierror, errstring
              STOP
            END IF
          END IF
        END IF
      END DO
    END DO

    DO K=KC,1,-1
      DO LT=1,NQSIJ
        I=IQS(LT)
        J=JQS(LT)
        IF( LIJLT(I,J) == 0 ) GOTO 100
        NS=NQSERQ(Lt)                                      
        !IF(flagwaspbC(ns,k) == 1 ) GOTO 100               
        LWASP=LWASP+1
        LDTM=LAUX(I,J,K)
        LUTM=0
        CALL hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF
  100 END DO
    END DO

    DO K=KC,1,-1
      DO LT=1,NQCTL
        I=IQCTLU(LT)
        J=JQCTLU(LT)
        LUTM=LAUX(I,J,K)
        I=IQCTLD(LT)
        J=JQCTLD(LT)
        LDTM=LAUX(I,J,K)
        IF( LUTM == 0 .AND. LDTM == 0 ) GOTO 200
        LWASP=LWASP+1
        CALL hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,1,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF
  200 END DO
    END DO

    IF( KC > 1 )THEN
      DO K=KS,1,-1
        KMUL1=KS-K
        KMUL2=KMUL1+1
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          IF( SWB(L) == 1. )THEN
            LWASP=LWASP+1
            LUTM=LT-1+KMUL1*LCLTM2
            LDTM=LT-1+KMUL2*LCLTM2
            CALL hlsetflowpath(Ihl_handle,INT(LWASP,4),LUTM,LDTM,3,ierror)
            IF( ierror  > 0 )THEN
              CALL hlgetlasterror(errstring)
              WRITE(6,6000) ierror, errstring
              STOP
            END IF
          END IF
        END DO
      END DO
    ENDIF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the number of segment constituents
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =

    inumsegconsts=3  !volume,depth,velocity
    ! ***   IF( ISTRAN(2) >= 1 ) inumsegconsts=inumsegconsts+1   !temperature modeled in EFDC and transfered
    ! ***   IF( ISTRAN(1) >= 1 ) inumsegconsts=inumsegconsts+1   !salinity modeled in EFDC and transfered
    IF( ISTRAN(2) >= 1 ) inumsegconsts=4   !temperature modeled in EFDC and transfered                   
    IF( ISTRAN(1) >= 1 ) inumsegconsts=5   !salinity modeled in EFDC and transfered                      

    CALL hlsetnumsegconsts(Ihl_handle, inumsegconsts, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set the number of flow path constituents
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! ***  when we add sed transport we need to add more
    CALL hlsetnumfpconsts(Ihl_handle, 3, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Now we will set all the segment constituent types
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetsegconsttype(Ihl_handle, 1, 0, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetsegconsttype(Ihl_handle, 2, 1, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetsegconsttype(Ihl_handle, 3, 2, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! ***  next only IF temperature is transfered
    IF( ISTRAN(2) >= 1 )THEN
    CALL hlsetsegconsttype(Ihl_handle, 4, 3, ierror)
    IF( ierror  > 0 )THEN
        CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! ***  next only IF salinity is transfered
    IF( ISTRAN(1) >= 1 )THEN
    CALL hlsetsegconsttype(Ihl_handle, 5, 4, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** Set all the flow constituent types
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetfpconsttype(Ihl_handle, 1, 0, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    CALL hlsetfpconsttype(Ihl_handle, 1, 1, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    CALL hlsetfpconsttype(Ihl_handle, 1, 2, ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsetvartimestep(Ihl_handle,0,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsethydtimestep(Ihl_handle,REAL(dt,4),ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    rinterval=(dt*NTSMMT)/86400.
    CALL hlsetupdateint(Ihl_handle,rinterval,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    CALL hlsethydtowaspratio(Ihl_handle,INT(NTSMMT,4),ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
    ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
    ! *** 
    ! *** INITIAL CONDITIONS WHEN IDAYS=0
    ! *** 
    IF( IDAYS == 0 )THEN
      IF( ISRESTI == 0 )THEN

        ! *** INITIAL CONDITION FOR A COLD START
        LWASP=0
        DO K=KC,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            LWASP=LWASP+1
            SegVel(LWASP)=0.0
            SegDepth(LWASP)=HP(L)*DZC(L,K)
            SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
            SegSalt(LWASP)=SAL(L,K)
            SegTemp(LWASP)=TEM(L,K)
            IF( NTSMMT < NTSPTC )THEN
              SegDepth(LWASP)=HP(L)*DZC(L,K)
              SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
              SegSalt(LWASP)=SALINIT(L,K)
              SegTemp(LWASP)=TEMINIT(L,K)
            END IF
          END DO
        END DO
        DO I=1,NCHN
          FLOW(I)=0.0
          CRNU(I)=0.0
          BRINTT(I)=0.0
        END DO

        CALL hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetseginfo(Ihl_handle,3,SegVel,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        ! ***     next only IF temperature is transfered
        IF( ISTRAN(2) >= 1 )THEN
          CALL hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
          IF( ierror  > 0 )THEN
            CALL hlgetlasterror(errstring)
            WRITE(6,6000) ierror, errstring
            STOP
          END IF
        END IF
        ! ***     next only IF salinity is transfered
        IF( ISTRAN(1) >= 1 )THEN
          CALL hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
          IF( ierror  > 0 )THEN
            CALL hlgetlasterror(errstring)
            WRITE(6,6000) ierror, errstring
            STOP
          END IF
        END IF
        ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
        CALL hlsetflowinfo(Ihl_handle,1,Flow,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetflowinfo(Ihl_handle,2,crnu,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetflowinfo(Ihl_handle,3,brintt,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlmomentcomplete(Ihl_Handle,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF
        ! *** END OF COLD START

      ELSE    

        ! *** INITIAL CONDITIONS FROM A RESTART FILE
        LWASP=0
        ! *** Segment Properties from the RESTART file
        DO K=KC,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            LN=LNC(L)
            LWASP=LWASP+1
            VELX=0.5*(U(L,K)+U(LEC(L),K))
            VELY=0.5*(V(L,K)+V(LN,K))
            VELZ=0.5*(W(L,K-1)+W(L,K))
            VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
            SegVel(LWASP)=VELMAG
            SegDepth(LWASP)=HP(L)*DZC(L,K)
            SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
            SegSalt(LWASP)=SAL(L,K)
            SegTemp(LWASP)=TEM(L,K)
          END DO
        END DO
        ! ***  Advection and dispersion in the X-direction:
        LWASP=0
        DO K=KC,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            ADDLW=0.0
            IF( SUB(L) == 1. )THEN
              LW=LWC(L)
              ADDLW=DYU(L)*0.5*(AH(L,K)+AH(LWC(L),K))*DZC(L,K)*0.5*(HP(L)+HP(LW))*DXIU (L)
                vol1=DXYP(L)*HP(L)*DZC(L,K)                
                vol2=DXYP(LW)*HP(LW)*DZC(L,K)              
                IF( vol1 < vol2) vol2=vol1                
                AD=vol2/dt*ADCOEFF                                 
                IF( addlw > AD )THEN                      
                  addlw=AD                         
                END IF                                     
                IF( IShdmf < 2 )THEN                       
                  addlw=0.                                 
                END IF                                     

            END IF
            IF( SUBO(L) == 1. )THEN
              LWASP=LWASP+1
              ! ***           IF( IJCTLT(I-1,J) == 8 ) addlw=0.0                  
              FLOW(LWASP)=UHDY(L,K)
              CRNU(LWASP)=2.*UHDY(L,K)*DYIU(L)*DXIU(L)/(HP(L)+HP(LWC(L)))
              BRINTT(LWASP)=ADDLW
            END IF
            IF( IJCTLT(I+1,J) == 8 )THEN
              IF( SUBO(LEC(L)) == 1. )THEN
                LWASP=LWASP+1
                FLOW(LWASP)=UHDY(LEC(L),K)
                CRNU(LWASP)=2.*UHDY(LEC(L),K)*DYIU(LEC(L))*DXIU(LEC(L))/(HP(L)+HP(LEC(L)))
                BRINTT(LWASP)=0.0                          
              END IF
            END IF
          END DO
          ! *** 
          ! ***  Advection and dispersion in the Y-direction:
          ! *** 
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            ADDLS=0.0
            LS=LSC(L)
            IF( SVB(L) == 1. )THEN
              ADDLS=DXV(L)*0.5*(AH(L,K)+AH(LS,K))*DZC(L,K)*0.5*(HP(L) +HP(LS))*DYIV (L)
                vol1=DXYP(L)*HP(L)*DZC(L,K)               
                vol2=DXYP(Ls)*HP(Ls)*DZC(L,K)             
                IF( vol1 < vol2) vol2=vol1               
                AD=vol2/dt*ADCOEFF                                
                IF( addls > AD )THEN                     
                  addls=AD                            
                END IF                                    
                IF( IShdmf < 2 )THEN                       
                  addls=0.                                 
                END IF                                     
            END IF
            IF( SVBO(L) == 1. )THEN
              LWASP=LWASP+1
              ! ***            IF( IJCTLT(I,J-1) == 8 ) addls=0.0                  
              FLOW(LWASP)=VHDX(L,K)
              CRNU(LWASP)=2.*VHDX(L,K)*DYIV(L)*DXIV(L)/(HP(L)+HP(LS))
              BRINTT(LWASP)=ADDLS
            END IF
            IF( IJCTLT(I,J+1) == 8 )THEN
              LN=LNC(L)
              IF( SVBO(LN) == 1. )THEN
                LWASP=LWASP+1
                FLOW(LWASP)=VHDX(LN,K)
                CRNU(LWASP)=2.*VHDX(LN,K)*DYIV(LN)*DXIV(LN)/(HP(L)+HP(LN))
                BRINTT(LWASP)=addls                            
              END IF
            END IF
          END DO
        END DO

    ! ***  Advection and dispersion in input flows
        DO K=KC,1,-1
          DO LT=1,nqsij
            I=Iqs(LT)
            J=Jqs(LT)
            IF( LIJLT(I,J) == 0 ) GOTO 310
            NS=NQSERQ(Lt)
            L=LQS(Lt)
            !IF(flagwaspbC(ns,k) == 1 ) GOTO 310                   
            LWASP=LWASP+1
            FLOW(LWASP)=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
            CRNU(LWASP)=flow(LWASP)/DXp(L)/dyp(l)/(HPK(L,K))
            ! ***          BRINTT(LWASP)=dyp(l)*ah(l,k)*HPK(L,K)/dxp(l)
            BRINTT(LWASP)=0.0
    310   END DO
        END DO
        ! ***   Advection and dispersion in structure flows
        DO K=KC,1,-1
          DO LT=1,nqctl
            Iu=Iqctlu(LT)
            Ju=Jqctlu(LT)
            Lu=Lij(iu,ju)
            Id=Iqctld(LT)
            Jd=Jqctld(LT)
            Ld=Lij(id,jd)
            IF( LU == 0 .AND. LD == 0 ) GOTO 410
            flowx=RQCMUL(Lt)*QCTLT(K,Lt,1)
            UDDXTMP=flowx/DXp(L)/dyp(l)/(HPK(L,K))
            IF( iu == id )THEN
              addls=dxv(lu)*ah(lu,k)*DZC(L,K)*0.5*(hp(lu)+hp(ld))*dyiv(lu)
            ELSE
              addls=dyu(lu)*ah(lu,k)*DZC(L,K)*0.5*(hp(lu)+hp(ld))*dxiu(lu)
            END IF
            IF( IShdmf < 2 )THEN                       
              addls=0.                                 
            END IF                                     
            LWASP=LWASP+1
            FLOW(LWASP)=FLOWX
            CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=ADDLS
    410   END DO
        END DO
  
        ! ***  Advection and dispersion in the Z-direction:
        IF( KC > 1 )THEN
          DO K=KS,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              addl=0.0
              ADDL1=ab(l,k)*hp(l)  !hnr Ev=ab*H
              IF( addl1 > ABWMX(l) )THEN
                addl1=ABWMX(l)     
              END IF
              IF( SPB(L) == 1. )THEN
                ADDL=DXYP(L)*addl1/hp(l)*DZIG(L,K)
                vol1=DXYP(L)*HP(L)*DZC(L,K)             
                vol2=DXYP(L)*HP(L)*DZC(L,K+1)           
                IF( vol1 < vol2) vol2=vol1             
                AD=ADCOEFF*vol2/dt                              
                IF( addl > AD )THEN                    
                  addl=AD                           
                END IF                                  
              END IF
              IF( SWB(L) == 1 )THEN
                LWASP=LWASP+1
                FLOW(LWASP)=-DXYP(L)*W(L,K)
                CRNU(LWASP)=W(L,K)*DZIG(L,K)/HP(L)
                BRINTT(LWASP)=ADDL
              END IF
            END DO
          END DO
        END IF
  
        CALL hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetseginfo(Ihl_handle,3,SegVel,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        ! ***  next only IF temperature is transfered
        IF( ISTRAN(2) >= 1 )THEN
          CALL hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
          IF( ierror  > 0 )THEN
            CALL hlgetlasterror(errstring)
            WRITE(6,6000) ierror, errstring
            STOP
          END IF
        END IF

        ! ***  next only IF salinity is transfered
        IF( ISTRAN(1) >= 1 )THEN
          CALL hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
          IF( ierror  > 0 )THEN
            CALL hlgetlasterror(errstring)
            WRITE(6,6000) ierror, errstring
            STOP
          END IF
        END IF

        ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
        CALL hlsetflowinfo(Ihl_handle,1,Flow,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetflowinfo(Ihl_handle,2,crnu,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlsetflowinfo(Ihl_handle,3,brintt,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF

        CALL hlmomentcomplete(Ihl_Handle,ierror)
        IF( ierror  > 0 )THEN
          CALL hlgetlasterror(errstring)
          WRITE(6,6000) ierror, errstring
          STOP
        END IF
        ! *** END OF HOTSTART

      END IF  ! *** END OF ISRESTI == 0
    
    END IF    ! ***  END INITIAL CONDITIONS WHEN IDAYS=0

    ! *** FINISH INITIALIZATION OF FILES (JSWASP=1)
    GOTO 3000
  
  END IF   ! *** END OF JSWASP=1

  ! ***************************************************************************
  IBEGIN=IDAYS*NTSPTC
  IF( N < IBEGIN) GOTO 3000

  ! *** ----------------------------------------------------------------------!
  ! *** 
  ! ***  WRITE TIME STEP DATA
  ! *** 
  ! ***  Advection and dispersion in the X-direction:
  ! *** 
  LWASP=0
  DO K=KC,1,-1
    DO LT=2,LALT
      I=ILLT(LT)
      J=JLLT(LT)
      L=LIJ(I,J)
      ADDLW=0.0
      IF( SUB(L) == 1. )THEN
        LW=LWC(L)
        ADDLW=DYU(L)*AHULPF(L,K)*DZC(L,K)*0.5*(HLPF(L) +HLPF(LW))*DXIU (L)
        vol1=DXYP(L)*HLPF(L)*DZC(L,K)                       
        vol2=DXYP(LW)*HLPF(LW)*DZC(L,K)                     
        IF( vol1 < vol2) vol2=vol1                       
        AD=ADCOEFF*vol2/dt                                        
        IF( addlw > AD )THEN                             
          addlw=AD                                    
        END IF                                            
        IF( IShdmf < 2 )THEN                       
          addlw=0.                                 
        END IF                                     
      END IF
      IF( SUBO(L) == 1. )THEN
        TMPVAL=UHLPF(L,K)+SVPT*UVPT(L,K)
        FLOWX=DYU(L)*TMPVAL*DZC(L,K)
        UDDXTMP=2.*TMPVAL*DXIU(L)/(HLPF(L)+HLPF(LWC(L)))
        LWASP=LWASP+1
        ! ***        IF( IJCTLT(I-1,J) == 8 ) addlw=0.0                  
        FLOW(LWASP)=FLOWX
        CRNU(LWASP)=UDDXTMP
        BRINTT(LWASP)=ADDLW
      END IF
      IF( IJCTLT(I+1,J) == 8 )THEN
        IF( SUBO(LEC(L)) == 1. )THEN
          TMPVAL=UHLPF(LEC(L),K)+SVPT*UVPT(LEC(L),K)
          FLOWX=DYU(LEC(L))*TMPVAL*DZC(L,K)
          UDDXTMP=2.*TMPVAL*DXIU(LEC(L))/(HLPF(LEC(L))+HLPF(L))
          IPTMP=I+1
          ! ***          addlw=0.0                                     
          LWASP=LWASP+1
          FLOW(LWASP)=FLOWX
          CRNU(LWASP)=UDDXTMP
          BRINTT(LWASP)=ADDLW
        END IF
      END IF
    END DO
    ! *** 
    ! ***  Advection and dispersion in the Y-direction:
    ! *** 
    DO LT=2,LALT
      I=ILLT(LT)
      J=JLLT(LT)
      L=LIJ(I,J)
      ADDLS=0.0
      IF( SVB(L) == 1. )THEN
        LS=LSC(L)
        ADDLS=DXV(L)*AHVLPF(L,K)*DZC(L,K)*0.5*(HLPF(L) +HLPF(LS))*DYIV (L)
        vol1=DXYP(L)*HLPF(L)*DZC(L,K)                      
        vol2=DXYP(LS)*HLPF(LS)*DZC(L,K)                    
        IF( vol1 < vol2) vol2=vol1                      
        AD=ADCOEFF*vol2/dt                                       
        IF( addls > AD )THEN                            
          addls=AD                                   
        END IF                                           
        IF( IShdmf < 2 )THEN                       
          addls=0.                                 
        END IF                                     
      END IF
      IF( SVBO(L) == 1. )THEN
        TMPVAL=VHLPF(L,K)+SVPT*VVPT(L,K)
        FLOWY=DXV(L)*TMPVAL*DZC(L,K)
        ! ***        FLOWY=DXV(L)*vlpf(l,k)*hlpf(l)*DZC(L,K)
        VDDYTMP=2.*TMPVAL*DYIV(L)/(HLPF(L)+HLPF(LSC(L)))
        JMTMP=J-1
        ! ***        IF( IJCTLT(I,J-1) == 8 ) addls=0.0                  
        LWASP=LWASP+1
        FLOW(LWASP)=FLOWY
        CRNU(LWASP)=VDDYTMP
        BRINTT(LWASP)=ADDLS
      END IF
      IF( IJCTLT(I,J+1) == 8 )THEN
        LN=LNC(L)
        IF( SVBO(LN) == 1. )THEN
          TMPVAL=VHLPF(LN,K)+SVPT*VVPT(LN,K)
          FLOWY=DXV(LN)*TMPVAL*DZC(L,K)
          VDDYTMP=2.*TMPVAL*DYIV(LN)/(HLPF(LN)+HLPF(L))
          JPTMP=J+1
          ! ***          addls=0.0                  
          LWASP=LWASP+1
          FLOW(LWASP)=FLOWY
          CRNU(LWASP)=VDDYTMP
          BRINTT(LWASP)=ADDLS
        END IF
      END IF
    END DO    ! *** ENDI OF LT=2,LALT LOOP
  END DO      ! *** END OF KC LOOP

  ! ***  Advection and dispersion in input flows
  DO K=KC,1,-1
    DO LT=1,nqsij
      I=Iqs(LT)
      J=Jqs(LT)
      IF( LIJLT(I,J) == 0 ) GOTO 300
      NS=NQSERQ(Lt)
      L=LQS(Lt)
      !IF(flagwaspbC(ns,k) == 1 ) GOTO 300                   
      flowx=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
      UDDXTMP=flowx/DXp(L)/dyp(l)/(HPK(L,K))
      ! ***      addlw=dyp(l)*ahulpf(l,k)*DZC(L,K)*hlpf(l)/dxp(l)  !hnr
      addlw=0.0                                    
      LWASP=LWASP+1
      FLOW(LWASP)=FLOWX
      CRNU(LWASP)=UDDXTMP
      BRINTT(LWASP)=ADDLW
300 END DO
  END DO
 
  ! ***   Advection and dispersion in structure flows
  DO K=KC,1,-1
    DO LT=1,nqctl
      Iu=Iqctlu(LT)
      Ju=Jqctlu(LT)
      Lu=Lij(iu,ju)
      Id=Iqctld(LT)
      Jd=Jqctld(LT)
      Ld=Lij(id,jd)
      IF( LU == 0 .AND. LD == 0 ) GOTO 400
      flowx=RQCMUL(Lt)*QCTLT(K,Lt,1)
      UDDXTMP=flowx/DXp(L)/dyp(l)*(HPKI(L,K))
      IF( iu == id )THEN
        addls=dxv(lu)*ahvlpf(lu,k)*DZC(L,K)*0.5*(hlpf(lu)+hlpf(ld))*dyiv(lu)
      ELSE
        addls=dyu(lu)*ahulpf(lu,k)*DZC(L,K)*0.5*(hlpf(lu)+hlpf(ld))*dxiu(lu)
      END IF
      IF( IShdmf < 2 )THEN                       
        addls=0.                                 
      END IF                                     
      LWASP=LWASP+1
      FLOW(LWASP)=FLOWX
      CRNU(LWASP)=UDDXTMP
      BRINTT(LWASP)=ADDLs
400 END DO
  END DO
  ! *** 
  ! ***  Advection and dispersion in the Z-direction:
  ! *** 
  IF( KC > 1 )THEN
    DO K=KS,1,-1
      DO LT=2,LALT
        I=ILLT(LT)
        J=JLLT(LT)
        L=LIJ(I,J)
        ADDL=0.0
        addl1=ablpf(l,k)*hlpf(l)  !hnr  Ev=ab*hp
        IF( addl1 > ABWMX(l) )THEN
          addl1=ABWMX(l)   !hnr
        END IF
        IF( SPB(L) == 1. )THEN
          ADDL=DXYP(L)*Addl1/HLPF(L)*DZIG(L,K)
          vol1=DXYP(L)*HLPF(L)*DZC(L,K)       
          vol2=DXYP(L)*HLPF(L)*DZC(L,K+1)     
          IF( vol1 < vol2) vol2=vol1       
          AD=ADCOEFF*vol2/dt                        
          IF( addl > AD )THEN              
            addl=AD                     
          END IF                            
        END IF
        IF( SWB(L) == 1 )THEN
          TMPVAL=WLPF(L,K)+SVPT*WVPT(L,K)
          FLOWZ=-DXYP(L)*TMPVAL
          WDDZTMP=TMPVAL*DZIG(L,K)/HLPF(L)
          LWASP=LWASP+1
          FLOW(LWASP)=FLOWZ
          CRNU(LWASP)=WDDZTMP
          BRINTT(LWASP)=ADDL
        END IF
      END DO
    END DO
  END IF

  ! *** 
  ! ***  Segment Properties:
  ! *** 
  LCELTMP=0
  DO K=KC,1,-1
    DO LT=2,LALT
      LCELTMP=LCELTMP+1
      I=ILLT(LT)
      J=JLLT(LT)
      L=LIJ(I,J)
      LN=LNC(L)
      VOLUM=DXYP(L)*HLPF(L)*DZC(L,K)
      IF( NTSMMT < NTSPTC) VOLUM=DXYP(L)*HP(L)*DZC(L,K)
      DEPTH=HLPF(L)*DZC(L,K)
      VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K) +UHLPF(LEC(L),K)+SVPT*UVPT(LEC(L),K))/HLPF(L)
      VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K) +VHLPF(LN,K)+SVPT*VVPT(LN,K) )/HLPF(L)
      VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1) +WLPF(L,K)+SVPT*WVPT(L,K) )
      VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
      SegSalt(LCELTMP)=SALLPF(L,K)
      SegTemp(LCELTMP)=TEMLPF(L,K)
      SegVolume(LCELTMP)=VOLUM
      SegDepth(LCELTMP)=DEPTH
      SegVel(LCELTMP)=VELMAG
    END DO
  END DO

  CALL hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
  IF( ierror  > 0 )THEN
    CALL hlgetlasterror(errstring)
    WRITE(6,6000) ierror, errstring
    STOP
  END IF

  CALL hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
  IF( ierror  > 0 )THEN
    CALL hlgetlasterror(errstring)
    WRITE(6,6000) ierror, errstring
    STOP
  END IF

  CALL hlsetseginfo(Ihl_handle,3,SegVel,ierror)
  IF( ierror  > 0 )THEN
    CALL hlgetlasterror(errstring)
    WRITE(6,6000) ierror, errstring
    STOP
  END IF

  ! ***  next only IF temperature is transfered
  IF( ISTRAN(2) >= 1 )THEN
    CALL hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
  END IF

  ! ***  next only IF salinity is transfered
  IF( ISTRAN(1) >= 1 )THEN
    CALL hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
    IF( ierror  > 0 )THEN
      CALL hlgetlasterror(errstring)
      WRITE(6,6000) ierror, errstring
      STOP
    END IF
  END IF

  ! == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
  CALL hlsetflowinfo(Ihl_handle,1,Flow,ierror)
  IF( ierror  > 0 )THEN
    CALL hlgetlasterror(errstring)
    WRITE(6,6000) ierror, errstring
    STOP
  END IF

  CALL hlsetflowinfo(Ihl_handle,2,crnu,ierror)
  IF( ierror  > 0 )THEN
    CALL hlgetlasterror(errstring)
    WRITE(6,6000) ierror, errstring
    STOP
  END IF

  CALL hlsetflowinfo(Ihl_handle,3,brintt,ierror)
  IF( ierror  > 0 )THEN
    CALL hlgetlasterror(errstring)
    WRITE(6,6000) ierror, errstring
    STOP
  END IF

  CALL hlmomentcomplete(Ihl_Handle,ierror)
  IF( ierror  > 0 )THEN
    CALL hlgetlasterror(errstring)
    WRITE(6,6000) ierror, errstring
    STOP
  END IF

  3000  JSWASP=0

  6000  format('Error ',I10, ' : ', A)
  9966  FORMAT(f12.4,I6,f12.4)
  9977  FORMAT(3I6,6f12.4)

  ! *** CLOSE LINKAGE FILE
  IF( N >= NTS )THEN
    CALL hlclose(IHL_HANDLE,ierror)
    IHL_HANDLE = 0
  ENDIF
#endif

  ! *** ----------------------------------------------------------------------!
  ! *** END SUBROUTINE WASPHYDROLINK
  ! *** ----------------------------------------------------------------------!
  RETURN

END
