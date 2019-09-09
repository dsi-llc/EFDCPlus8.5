      Subroutine WASP7HYDRO             ! writes the hyd file for WINASP water quality model.

C **  Mar 27 2008 A Stoddard renamed linkage subroutine wasphydrolink
C **  to prevent confusion with other version of wasp7epa.for used in efdc_ds
C **  this version is absolutely needed to generate wasp7 hyd linkage files
C **  for flow, volume, salinity and water temperature
C **  the other version of wasp7epa.for that Paul has in efdc_ds source code
C **  dated Feb 28 2008 does NOT provide salinity and water temperature linkage
C **  calls to wasp7epa changed in efdc source code

C ** C:\Models\EFDC\SourceCode\EFDCWASP-2008-MAR-27\
C ** call to wasp7epa changed in CALMMT.for(873):      IF(ISWASP.EQ.17) CALL WASP7EPA
C ** 1 occurrence(s) have been found.

c-----previous subroutine name -- SUBROUTINE WASP7EPA
C------------------------------------------------------------------------------
C PURPOSE:
C
C   Subroutine WASPHYDROLINK writes the hyd file for WINASP water quality model.
C
C VARIABLE LIST:
C
C MODIFICATION HISTORY:
C
C   Date       Author         Comments
C   ---------- -------------- -------------------------------------------------
C   30/10/2003 Hugo Rodriguez This version uses a dll to transfer the hydrodynamic
C              Tim Wool       data, including salinity and temperature, from EFDC
C                             to WASP version 6.1.
C   05/02/2005 Hugo Rodriguez Eliminate dispersion for boundary flows.
C   05/04/2005 Hugo Rodriguez Limit maximum dispersion to cell volume/time step
C   06/15/2005 Hugo Rodriguez Zero out AH coefficient transferred to wasp if
C                             ISHDMF<2
C   06/24/2005 Hugo Rodriguez Limit Ab to Abwmax
C   10/27/2005 Hugo Rodriguez Calculate the initial date of wasp instead of reading it
C    4/11/2006 Hugo Rodriguez Add qfactor
c    5/23/2006 Hugo Rodriguez Add a limit maximum dispersion specific to certain cells
C------------------------------------------------------------------------------
C
      use global

c **  3/27/2008 A Stoddard changed filename from efdc.cmn to efdc_wasphydro.cmn
C **  INCLUDE 'EFDC.CMN'
      INCLUDE 'WASP\Not Used\EFDC_WASPHYDRO.CMN'
      INCLUDE 'WASP\Not Used\hydrolink_set.INT'
      
      parameter(nf=31000,nsg=12000)
C
      DIMENSION LAUX(ICM,JCM,KCM)
      INTEGER LU,LD,IU,ID,JU,JD,NAUX
      INTEGER istartyear,istartmonth,istartday,istarthour,istartminute
     *,istartsecond,Ihl_debug,Ihl_mode,inumsegconsts,j,
     *FLAGWASPBC(NQSERM,KCM)
      CHARACTER*6  sn(lcm),sn1
      CHARACTER*20 HYDFIL
      CHARACTER*23 segname
      CHARACTER*17 SN2
      character*80 DESCRIPTION,MODELERNAME
      character*3 Itext, Jtext, Ktext
      character*80 RECHNR
	REAL*8  AUX,tmin,thour,tsec
        REAL  rinterval,crnu(nf),brintt(nf),flow(nf)
        REAL  SegVolume(nsg),SegDepth(nsg),SegVel(nsg)
        REAL  SegSalt(nsg),SegTemp(nsg),vol1,vol2,ad,adcoeff,abwmax                !hnr
        real  abwmx(nsg)                                                   !hnr
c        REAL*4  SegSalt(nsg),SegTemp(nsg)                                            !hnr
C
C
C**********************************************************************C
C
C **  READ CONTROL DATA FILE EFDC.WSP
C
C----------------------------------------------------------------------C
C

      


      SVPT=1.
      IF(NTSMMT.LT.NTSPTC)SVPT=0.
C
      IF(JSWASP.EQ.1) THEN

c      for jswasp=1 only first entry
        OPEN(1,FILE='EFDC.WSP',STATUS='UNKNOWN')
        WRITE(6,*)'EFDC.WSP opened'

C
C1**  READ CELL VOLUME PARAMETERS
C
        READ(1,1)
        READ(1,1)
        READ(1,*) IVOPT,IBEDV,SCALV,CONVV,VMULT,VEXP,DMULT,DEXP
      write(*,*)'<EFDC.WSP1> ',IVOPT,IBEDV,SCALV,CONVV,
     +             VMULT,VEXP,DMULT,DEXP

C
C2**  READ DIFFUSION PARAMETERS
C
        READ(1,1)
        READ(1,1)
        READ(1,*) NRFLD,SCALR,CONVR,ISNKH,adcoeff,abwmax					!hnr
c        READ(1,*) NRFLD,SCALR,CONVR,ISNKH                                                       !hnr
      write(*,*)'<EFDC.WSP2> ',NRFLD,SCALR,CONVR,ISNKH,adcoeff,abwmax					!hnr
        DO LT=2,LALT                                      !hnr
          l=lij(illt(lt),jllt(lt))                        !hnr
          abwmx(l)=abwmax                                 !hnr
        END DO                                            !hnr
C
C3**  READ ADVECTION PARAMETERS
C
        READ(1,1)
        READ(1,1)
        READ(1,*) IQOPT,NFIELD,SCALQ,CONVQ,HYDFIL,ISWASPD,ISDHD,IDAYS                                          !hnr
        write(*,*)'<EFDC.WSP3> ', IQOPT,NFIELD,SCALQ,CONVQ,
     +  HYDFIL,ISWASPD,ISDHD,IDAYS                                          !hnr

C
C4**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
C
        READ(1,1)
        READ(1,1)
        READ(1,*) DEPSED,TDINTS,SEDIFF, WSS1, WSS2, WSS3
        write(*,*)'<EFDC.WSP4> ',DEPSED,TDINTS,SEDIFF,WSS1,WSS2,WSS3


C      *  istartyear,istartmonth,istartday
C5**  READ SEDIMENT VOLUME DEPTH AND TDINTS(GROUP C RECORD 1)
C
        READ(1,1)
        READ(1,1)
        READ(1,*) iyear,imon,iday
        write(*,*)'<EFDC.WSP5> ',iyear,imon,iday

        do i=1,5                                   !hnr
          read(1,*,err=11) RECHNR                        !hnr
          write(*,*)'<EFDC.WSP6> ',RECHNR
        end do                                     !hnr
11      continue       
        DO LT=2,LALT                               !hnr
          read(1,*,err=12)i,j,abwmax               !hnr
          l=lij(i,j)                               !hnr
          abwmx(l)=abwmax                          !hnr
          write(*,*)'<EFDC.WSP7> ',i,j,abwmax               !hnr
        END DO                                     !hnr
12      continue                                   !hnr

        CLOSE(1)
        WRITE(*,*)'EFDC.WSP read succesfully..now write ABMax'

        OPEN(1,FILE='ABmax.txt',STATUS='UNKNOWN')   !hnr
        write(1,*)'    I    J     ABmax'            !hnr
        DO LT=2,LALT                                      !hnr
          l=lij(illt(lt),jllt(lt))                        !hnr
          write(1,21)illt(lt),jllt(lt),abwmx(l)           !hnr
        END DO                                            !hnr
        close(1)
21      format(2I5,f10.6)
        WRITE(6,*)'EFDC.WSP read succesfully and ABmax.txt written'
C
    1   FORMAT (80X)
C
C     read qser file to check for flows only to some layers

        IF(NQSER.GE.1)THEN                                                !HNR
          OPEN(1,FILE='QSER.INP',STATUS='UNKNOWN')                        !HNR
C **  SKIP OVER TITLE AND AND HEADER LINES                              !HNR
C                                                                       !HNR
          DO IS=1,14                                                      !HNR
           READ(1,1)                                                      !HNR
          ENDDO                                                           !HNR
C                                                                       !HNR
c          DO NS=1,NQSER                                                   !HNR
c            READ(1,*,IOSTAT=ISO)ISTYP, MQSER(NS),TCQSER(NS),TAQSER(NS),   !HNR
c     &                   RMULADJ,ADDADJ,ICHGQS                          !HNR
c            IF(ISO.GT.0) GOTO 860                                         !HNR
c            IF(ISTYP.EQ.1)THEN                                            !HNR
c              READ(1,*,IOSTAT=ISO) (WKQ(K),K=1,KC)                        !HNR
c              do k=1,kc                                                   !HNR
c                if(wkq(k).lt.1.e-6) flagwaspbC(ns,k)=1                    !HNR
c              end do                                                      !HNR
c              IF(ISO.GT.0) GOTO 860                                       !HNR
c              DO M=1,MQSER(NS)                                            !HNR
c                READ(1,*,IOSTAT=ISO)TQSER(M,NS),QSERTMP                   !HNR
c                IF(ISO.GT.0) GOTO 860                                     !HNR
c              END DO                                                      !HNR
c            ELSE                                                          !HNR
c              DO M=1,MQSER(NS)                                            !HNR
c                READ(1,*,IOSTAT=ISO)TQSER(M,NS),(QSER(M,K,NS), K=1,KC)    !HNR
c                IF(ISO.GT.0) GOTO 860                                     !HNR
c              END DO                                                      !HNR
c            END IF                                                        !HNR
c          END DO                                                          !HNR
        END IF                                                            !HNR
C                                                                       !HNR
        CLOSE(1)                                                          !HNR
        goto 862                                                          !HNR
860     WRITE(6,861)                                                      !HNR
861     FORMAT('  READ ERROR FOR FILE QSER.INP ')                         !HNR
862     continue                                                          !HNR
C
C**********************************************************************C
C
C **  DEFINE EFDC-WASP CORRESPONDENCE AND INITIALIZE FILES
C

C**********************************************************************C
C      begin code inserted by Andy Stoddard 2-26-2008
C------ Feb-26-2008 A Stoddard add code to write WASP7SEG_EFDCIJK 
C------ file listing of EFDC L,IJK cells and Wasp segment/layer data
C----------------------------------------------------------------------C
C----------------------------------------------------------------------C
C---- 7/1/2005 A Stoddard added code to write WASP7SEG_EFDCIJK.DAT file
C---- Segment Name file is linked for import to WASP6 & WASP7
      IF (JSWASP.EQ.1) THEN
        OPEN(970,FILE='WASP7SEG_EFDCIJK.DAT',STATUS='UNKNOWN')  ! 7/1/2005 AS
        CLOSE(970,STATUS='DELETE')                              ! 7/1/2005 AS 

c        OPEN(90,FILE='WASPP.OUT',STATUS='UNKNOWN')
c        OPEN(93,FILE='WASPC.OUT',STATUS='UNKNOWN')
c        CLOSE(90,STATUS='DELETE')
c        CLOSE(93,STATUS='DELETE')
c        OPEN(90,FILE='WASPP.OUT',STATUS='UNKNOWN')
c        OPEN(93,FILE='WASPC.OUT',STATUS='UNKNOWN')
        OPEN(970,FILE='WASP7SEG_EFDCIJK.DAT',STATUS='UNKNOWN')  ! 7/1/2005 AS

        WRITE(970,9702)              ! write header records

c        WRITE(93,1031)IVOPT,IBEDV,TDINTS,TITLEC
c        WRITE(93,1032)SCALV,CONVV
        LCLTM2=LCLT-2
        LWASP=0
        IF (KC.GT.1) THEN
          LTYPE=1
          KWASP=1
          DO LT=2,LALT
            LWASP=LWASP+1
            LBELOW=LWASP+LCLTM2
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            DMULT=HLPF(L)*DZC(KC)
            VOLUME=DXYP(L)*HLPF(L)*DZC(KC)
            IF(NTSMMT.LT.NTSPTC) THEN
              DMULT=HMP(L)*DZC(KC)
              VOLUME=DXYP(L)*HMP(L)*DZC(KC)
            END IF
            LAUX(I,J,KC)=LWASP      !HNR

            WRITE(970,9701)L,I,J,KC,LWASP,KWASP         ! ,I,J,L,KC               ! 7/1/2005 AS
c            WRITE(90,1001)LWASP,KWASP,I,J,L,KC
c            WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, DMULT,
c     +      DEXP,I,J,L,KC
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
              DMULT=HLPF(L)*DZC(K)
              VOLUME=DXYP(L)*HLPF(L)*DZC(K)
              IF(NTSMMT.LT.NTSPTC) THEN
                DMULT=HMP(L)*DZC(KC)
                VOLUME=DXYP(L)*HMP(L)*DZC(KC)
              END IF
              LAUX(I,J,K)=LWASP     !HNR
              WRITE(970,9701)L,I,J,K,LWASP,KWASP         ! ,I,J,L,K               ! 7/1/2005 AS
c              WRITE(90,1001)LWASP,KWASP,I,J,L,K
c              WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, DMULT,
c     +        DEXP,I,J,L,KC
            END DO
          END DO
        END IF
        LTYPE=2
        IF (KC.EQ.1) LTYPE=1
        KWASP=KC
        DO LT=2,LALT
          LWASP=LWASP+1
          LBELOW=LWASP+LCLTM2
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          DMULT=HLPF(L)*DZC(1)
          VOLUME=DXYP(L)*HLPF(L)*DZC(1)
          IF(NTSMMT.LT.NTSPTC) THEN
            DMULT=HMP(L)*DZC(KC)
            VOLUME=DXYP(L)*HMP(L)*DZC(KC)
          END IF
          IONE=1
          LAUX(I,J,1)=LWASP       !HNR
          WRITE(970,9701)L,I,J,IONE,LWASP,KWASP         ! ,I,J,L,IONE               ! 7/1/2005 AS
c          WRITE(90,1001)LWASP,KWASP,I,J,L,IONE
c          WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, DMULT,
c     +    DEXP,I,J,L,IONE
        END DO
        LTYPE=3
        KWASP=KC+1
        DXYSUM=0.
        LWSPTMP=LWASP+1
        DO LT=2,LALT
          LWSPTMP=LWSPTMP+1
        END DO
C The following the lower benthic layer.  All upper benthic layer segments
C have this layer immediately below them:
        DO LT=2,LALT
          LWASP=LWASP+1
          LBELOW=LWSPTMP
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          DXYSUM=DXYSUM+DXYP(L)
          VOLUME=DXYP(L)*DEPSED
          IZERO=0
c          WRITE(90,1001)LWASP,KWASP,I,J,L,IZERO
c          WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, DEPSED,
c     +    DEXP,I,J,L,IZERO
        END DO
C Next do the lower benthic layer:
        LTYPE=4
        KWASP=KC+2
        LWASP=LWASP+1
        LBELOW=0
        DMULT=DEPSED
        VOLUME=DXYSUM*DEPSED
        IM1=-1
c        WRITE(90,1001)LWASP,KWASP,I,J,L,IM1
c        WRITE(93,1033)LWASP,LBELOW,LTYPE,VOLUME,VMULT,VEXP, DMULT,DEXP,
c     +  I,J,L,IM1
        CLOSE(970)                           ! 7/1/2005
c        CLOSE(90)
c        CLOSE(93)
      END IF

c 1001 FORMAT(6I5,2F10.4)
c 1031 FORMAT(2I5,F10.4,10X,A50)
c 1032 FORMAT(2F10.4)
c 1033 FORMAT(3I10,F10.1,4F10.3,'   !',4i5)
       

c 9701 FORMAT('I=',I3,' J=',I3,' L=',I4 ,' K=',I2,
c     +       ' W',I4,' WK=',I2)

 9702 format(
     + '######################################################',/,
     + '# File    = WASP7SEG_EFDCIJK.DAT',/,
     + '# Field#1 = EFDC Cell L-Number',/,
     + '# Field#2 = EFDC Cell I-Index',/,
     + '# Field#3 = EFDC Cell J-Index',/,
     + '# Field#4 = EFDC Cell K-Layer  [1=Bottom;KC=Surface]',/,
     + '# Field#5 = WASP Segment Number',/,
     + '# Field#6 = WASP Segment Layer [KC=Bottom;1=Surface]',/,
     + '######################################################')     

 9701 FORMAT(I4,',',             ! EFDC L-coordinate
     +       I4,',',             ! EFDC I-coordinate
     +       I5,',',             ! EFDC J-cell index 
     +       I2,',',             ! EFDC K-layer number   1=bottom;  KC=surface layer
     +       I5,',',             ! WASP Segment Number
     +       I2)                 ! WASP Layer Number     KC=bottom; 1= surface layer

C             WASPSEG.OUT data structure
C             Wasp7 wired for linkage of A30 character string file
C             W Wasp Segment Number
C             WK wasp Layer number 1=sfc   K=KC=bottom
C             I,J = EFDC I, J index 
C             L = EFDC grid cell number
C             K = EFDC layer K=KC sfc    K=1 bottom 

C      end code inserted by Andy Stoddard 2-26-2008
C**********************************************************************C

C----------------------------------------------------------------------C
C
c=======================================================================
c      hlopen parameters
c=======================================================================
        Ihl_handle = 0
        Ihl_debug = 1
        Ihl_mode = 1	!Ihl_mode=0 to READ from dll hyd file, =1 to WRITE to dll hyd file
c=======================================================================
c     Set the debug flag 0=No debug 1=Debug (LOG.OUT)
c=======================================================================
        call hlsetdebug(Ihl_debug)
c=======================================================================
c     Open the file
c=======================================================================
        call hlopen(HYDFIL,Ihl_mode,Ihl_handle,ierror)
        if(ierror.gt.0) then
          call hlgetlasterror(errstring)
          write(6,6000)ierror,errstring
          stop
        end if
c=======================================================================
c     Set the language to FORTRAN
c=======================================================================
        call hlsetlanguage(Ihl_handle,1,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Store a description string
c=======================================================================
        DESCRIPTION='   '//char(0)
        call hladddescription(Ihl_handle,0,DESCRIPTION,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Store the modeler name
c=======================================================================
        MODELERNAME='Created by:  '//char(0)
        call hladddescription(Ihl_handle,1,MODELERNAME,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Set the creator
c=======================================================================
        call hlsetcreator(Ihl_handle, 1, ierror)   !1 for fortran program
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Set the seed moment (start date of the hyd file)
c=======================================================================
c     Calculate istartmonth, istartday,istartyear
        j=1461*(iyear+4800+int((imon-14)/12))/4+367*(imon-2-12*
     &    int((imon-14)/12))/12-(3*((iyear+4800+int((imon-14)/12)
     &    +100)/100))/4+iday-32075+tbegin+idays
	istartyear=100*(int(4*(j+68569)/146097)-49)+int(4000*
     &   (int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)+1)/
     &   1461001)+int(int(80*int(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*
     &   int(4*int(j+68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartmonth=int(80*int(int(j+68569-(146097*int(4*int(j+68569)/
     &   146097)+3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+
     &   68569)/146097)+3)/4)+1)/1461001)/4+31)/2447)+2-12*int(int(80*
     &   int(int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)
     &   -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)+1)/1461001)/4+31)/2447)/11)
	istartday=int(int(j+68569-(146097*int(4*int(j+68569)/146097)
     &   +3)/4)-1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)
     &   /146097)+3)/4)+1)/1461001)/4+31)-2447*int(80*int(
     &   int(j+68569-(146097*int(4*int(j+68569)/146097)+3)/4)
     &   -1461*int(4000*(int(j+68569-(146097*int(4*int(j+68569)/
     &   146097)+3)/4)+1)/1461001)/4+31)/2447)/80
	istarthour=0
	istartminute=0
	istartsecond=0


	

c     Adjust hour, minutes and seconds for the start day if necessary
	    IBEGIN=IDAYS*NTSPTC
	    AUX=FLOAT(IBEGIN)/FLOAT(NTSMMT)
	    NAUX=INT(AUX)
	      IF(AUX.GT.NAUX) then
	        AUX=(NAUX+1.-AUX)
                tsec=AUX*NTSMMT*TCON/NTSPTC
                istartsecond=INT(tsec)
                IF(tsec.GT.60) THEN
                  TMIN=tsec/60.0
                  istartminute=INT(TMIN)
                  istartsecond=INT((TMIN-istartminute)*60.)
                  IF(istartminute.GT.60) THEN
                    THOUR=FLOAT(istartminute)/60.0
                    istarthour=INT(THOUR)
                    istartminute=INT((THOUR-istarthour)*60.)
                  ENDIF
                ENDIF
              ENDIF
c=======================================================================

        call hlsetseedmoment(Ihl_handle,istartmonth,istartday,
     +     istartyear,istarthour,istartminute,istartsecond,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================

	call hlsetnumlayers(Ihl_handle,kc,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Set the number of segments
c=======================================================================
        NJUN=KC*(LCLT-2)
        if(njun.gt.nsg) then
        write(6,500)NJUN,NSG
        STOP
        END IF
500     FORMAT('THE NUMBER OF WASP SEGMENTS IN YOUR APPLICATION',I6,1x,
     +'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumsegments(Ihl_handle,njun,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c               Set WASP Segment Names (defaulft as I,J,K)
c=======================================================================
        DO LT=2,LALT
          l=lij(illt(lt),jllt(lt))
          sn(l)='      '//char(0)
        END DO
        OPEN(94,FILE='segname.inp',iostat=ios,STATUS='old')
        IF(ios.eq.0) then
          do i=1,4
            read(94,*)
          end do
          do kk=1,la
            READ(94,*,err=111)I,J,SN1
            L=LIJ(I,J)
            sn(l)=sn1
          end do
        end if
111     CLOSE(94)
        I=0
        DO K=KC,1,-1
          DO LT=2,LALT
    	    I=I+1
	    LAUX(ILLT(LT),JLLT(LT),K)=I
	    l=lij(illt(lt),jllt(lt))
             if(il(l).ge.100) then
              WRITE(itext,"(I3)")IL(L)
             elseif(il(l).ge.10) then
              WRITE(itext,"(I2)")IL(L)
             else
              WRITE(itext,"(I1)")IL(L)
             end if
             if(jl(l).ge.100) then
             WRITE(jtext,"(I3)")jL(L)
           elseif(jl(l).ge.10) then
            WRITE(jtext,"(I2)")jL(L)
           else
            WRITE(jtext,"(I1)")jL(L)
           end if
           if(k.ge.100) then
             WRITE(ktext,"(I3)")k
           elseif(k.ge.10) then
            WRITE(ktext,"(I2)")k
            else
            WRITE(ktext,"(I1)")k
            end if
	    sn2=' I='//itext//' J='//jtext//' K='//ktext
    	    segname=sn2//sn(l)//char(0)
	    call hlsetsegname(ihl_handle,i,segname,ierror)
          END DO
        END DO
c=======================================================================
c     Set the number of flow paths
c=======================================================================
        NCHNH=0
        NCHNV=0
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          NCHNH=NCHNH+INT(SUBO(L))
          IF (IJCTLT(I+1,J).EQ.8) THEN
            IF (SUBO(L+1).EQ.1.) NCHNH=NCHNH+1
          END IF
          NCHNH=NCHNH+INT(SVBO(L))
          IF (IJCTLT(I,J+1).EQ.8) THEN
            IF (SVBO(LNC(L)).EQ.1.) NCHNH=NCHNH+1
          END IF
          NCHNV=NCHNV+INT(SWB(L))
        END DO
        NCHN=KC*NCHNH+(KC-1)*NCHNV

	NQ=NQSIJ
	DO L=1,NQSIJ
	  IF(LIJLT(IQS(L),JQS(L)).EQ.0) NQ=NQ-1
	END DO
        NCHN=NCHN+KC*NQ
c
        NQ=0                                        !HNR
	DO L=1,NQSIJ                                !HNR
            I=IQS(L)                               !HNR
            J=JQS(L)                               !HNR
	      IF(LIJLT(I,J).EQ.0) GOTO 1001         !HNR
              NS=NQSERQ(L)                         !HNR
	  DO K=1,KC                                 !HNR
	    IF(flagwaspbC(ns,k).EQ.1) NQ=NQ+1       !HNR
	  END DO                                    !HNR
1001    END DO                                      !HNR
C        WRITE(6,*)NQ
        NCHN=NCHN-NQ
c
	NQ=NQCTL
	DO L=1,NQCTL
	  IF(LIJLT(IQCTLU(L),JQCTLU(L)).EQ.0) THEN
	    IF(LIJLT(IQCTLD(L),JQCTLD(L)).EQ.0) NQ=NQ-1
	  END IF
	END DO
        NCHN=NCHN+KC*NQ
        if(NCHN.gt.NF) then
        write(6,600)NCHN,NF
        STOP
        END IF
600     FORMAT('THE NUMBER OF WASP FLOWS IN YOUR APPLICATION',I6,1X,
     +   'IS GREATER THAN THE ARRAY DIMENSION:',I7)
        call hlsetnumflowpaths(Ihl_handle, NCHN, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Set the flow path and direction
c=======================================================================
        LCLTM2=LCLT-2
        LWASP=0
        DO K=KC,1,-1
          KMUL=KC-K
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (SUBO(L).EQ.1.) THEN
              LWASP=LWASP+1
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LDTM-1
              IF (IJCTLT(I-1,J).EQ.8) LUTM=0
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
            END IF

            IF (IJCTLT(I+1,J).EQ.8) THEN
              IF (SUBO(L+1).EQ.1.) THEN
	        LWASP=LWASP+1
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
              END IF
            END IF
          END DO

          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            IF (SVBO(L).EQ.1.) THEN
              LWASP=LWASP+1
              LSLT=LSCLT(LT)
              LDTM=LT-1+KMUL*LCLTM2
              LUTM=LSLT-1+KMUL*LCLTM2
              IF(IJCTLT(I,J-1).EQ.8) LUTM=0
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
            END IF
            IF (IJCTLT(I,J+1).EQ.8) THEN
              LN=LNC(L)
              IF (SVBO(LN).EQ.1.) THEN
                LWASP=LWASP+1
                LSLT=LSCLT(LT)
                LDTM=0
                LUTM=LT-1+KMUL*LCLTM2
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,2,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
              END IF
            END IF
          END DO
        END DO

        DO K=KC,1,-1
          DO LT=1,NQSIJ
            I=IQS(LT)
            J=JQS(LT)
	      IF(LIJLT(I,J).EQ.0) GOTO 100
              NS=NQSERQ(Lt)                                          !HNR
                IF(flagwaspbC(ns,k).EQ.1) GOTO 100                   !HNR
            LWASP=LWASP+1
            LDTM=LAUX(I,J,K)
            LUTM=0
C              WRITE(6,*)LT,K,NS,FLAGWASPBC(NS,K),LDTM,LUTM
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
100       END DO
        END DO

        DO K=KC,1,-1
          DO LT=1,NQCTL
            I=IQCTLU(LT)
            J=JQCTLU(LT)
            LUTM=LAUX(I,J,K)
            I=IQCTLD(LT)
            J=JQCTLD(LT)
            LDTM=LAUX(I,J,K)
	      IF(LUTM.EQ.0.AND.LDTM.EQ.0) GOTO 200
              LWASP=LWASP+1
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,1,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
200       END DO
        END DO

        IF(KC.GT.1)THEN
          DO K=KS,1,-1
            KMUL1=KS-K
            KMUL2=KMUL1+1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              IF (SWB(L).EQ.1.) THEN
                LWASP=LWASP+1
                LUTM=LT-1+KMUL1*LCLTM2
                LDTM=LT-1+KMUL2*LCLTM2
               call hlsetflowpath(Ihl_handle,LWASP,LUTM,LDTM,3,ierror)
               if(ierror .gt. 0)then
                  call hlgetlasterror(errstring)
                  write(6,6000) ierror, errstring
                  stop
               end if
              END IF
            END DO
          END DO
        ENDIF
c=======================================================================
c     Set the number of segment constituents
c=======================================================================

	inumsegconsts=3  !volume,depth,velocity
c	if(istran(2).GE.1) inumsegconsts=inumsegconsts+1   !temperature modeled in EFDC and transfered
c	if(istran(1).GE.1) inumsegconsts=inumsegconsts+1   !salinity modeled in EFDC and transfered
	if(istran(2).GE.1) inumsegconsts=4   !temperature modeled in EFDC and transfered                       !hnr
	if(istran(1).GE.1) inumsegconsts=5   !salinity modeled in EFDC and transfered                          !hnr

        call hlsetnumsegconsts(Ihl_handle, inumsegconsts, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Set the number of flow path constituents
c=======================================================================
c when we add sed transport we need to add more
        call hlsetnumfpconsts(Ihl_handle, 3, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c     Now we will set all the segment constituent types
c=======================================================================
        call hlsetsegconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
        call hlsetsegconsttype(Ihl_handle, 2, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
        call hlsetsegconsttype(Ihl_handle, 3, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
c next only if temperature is transfered
        if(istran(2).GE.1) then
        call hlsetsegconsttype(Ihl_handle, 4, 3, ierror)
        if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        end if
c=======================================================================
c next only if salinity is transfered
        if(istran(1).GE.1) then
        call hlsetsegconsttype(Ihl_handle, 5, 4, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        end if
c=======================================================================
c     Set all the flow constituent types
c=======================================================================
        call hlsetfpconsttype(Ihl_handle, 1, 0, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 1, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
        call hlsetfpconsttype(Ihl_handle, 1, 2, ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
	call hlsetvartimestep(Ihl_handle,0,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
	call hlsethydtimestep(Ihl_handle,dt,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
	rinterval=(dt*NTSMMT)/86400.
	call hlsetupdateint(Ihl_handle,rinterval,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
	call hlsethydtowaspratio(Ihl_handle,NTSMMT,ierror)
        if(ierror .gt. 0)then
          call hlgetlasterror(errstring)
          write(6,6000) ierror, errstring
          stop
        end if
c=======================================================================
C
C     INITIAL CONDITIONS WHEN IDAYS=0
C
        IF(IDAYS.EQ.0) THEN
          IF(ISRESTI.EQ.0) THEN

C     INITIAL CONDITION FOR A COLD START
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              LWASP=LWASP+1
              SegVel(LWASP)=0.0
              SegDepth(LWASP)=HP(L)*DZC(K)
              SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
              SegSalt(LWASP)=SAL(L,K)
              SegTemp(LWASP)=TEM(L,K)
              IF(NTSMMT.LT.NTSPTC) THEN
                SegDepth(LWASP)=HMP(L)*DZC(K)
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
	  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
c         next only if temperature is transfered
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              stop
            end if
          end if
c         next only if salinity is transfered
          if(istran(1).GE.1) then
	    call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              stop
            end if
          end if
c=======================================================================
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
          ELSE
C     INITIAL CONDITIONS FROM A RESTART FILE
          LWASP=0
C     Segment Properties from the RESTART file
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              LN=LNC(L)
              LWASP=LWASP+1
              VELX=0.5*(U(L,K)+U(L+1,K))
              VELY=0.5*(V(L,K)+V(LN,K))
              VELZ=0.5*(W(L,K-1)+W(L,K))
              VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
              SegVel(LWASP)=VELMAG
              SegDepth(LWASP)=HP(L)*DZC(K)
              SegVolume(LWASP)=SegDepth(LWASP)*DXYP(L)
              SegSalt(LWASP)=SAL(L,K)
              SegTemp(LWASP)=TEM(L,K)
            END DO
          END DO
C Advection and dispersion in the X-direction:
          LWASP=0
          DO K=KC,1,-1
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLW=0.0
              IF (SUB(L).EQ.1.) THEN
                LW=L-1
                ADDLW=DYU(L)*0.5*(AH(L,K)+AH(L-1,K))*DZC(K)*
     +          0.5*(HP(L)+HP(LW))*DXIU (L)
                  vol1=DXYP(L)*HP(L)*DZC(K)                    !hnr
                  vol2=DXYP(LW)*HP(LW)*DZC(K)                  !hnr
                  if (vol1.LT.vol2) vol2=vol1                    !hnr
                  ad=vol2/dt*adcoeff                                     !hnr
                  if (addlw.GT.ad) then                          !hnr
                    addlw=ad                             !hnr
                  end if                                         !hnr
                  if(ishdmf.lt.2) then                           !hnr
                    addlw=0.                                     !hnr
                  end if                                         !hnr

              END IF
              IF (SUBO(L).EQ.1.) THEN
                LWASP=LWASP+1
c               IF (IJCTLT(I-1,J).EQ.8) addlw=0.0                      !hnr
                FLOW(LWASP)=UHDY(L,K)*DZC(K)
	        CRNU(LWASP)=2.*UHDY(L,K)*DYIU(L)*DXIU(L)/(HP(L)+HP(L-1))
  	        BRINTT(LWASP)=ADDLW
              END IF
              IF (IJCTLT(I+1,J).EQ.8) THEN
                IF (SUBO(L+1).EQ.1.) THEN
                  LWASP=LWASP+1
	          FLOW(LWASP)=UHDY(L+1,K)*DZC(K)
	          CRNU(LWASP)=2.*UHDY(L+1,K)*DYIU(L+1)*DXIU(L+1)
     +            /(HP(L)+HP(L+1))
                  BRINTT(LWASP)=0.0                              !hnr
                END IF
              END IF
            END DO
C
C Advection and dispersion in the Y-direction:
C
            DO LT=2,LALT
              I=ILLT(LT)
              J=JLLT(LT)
              L=LIJ(I,J)
              ADDLS=0.0
              LS=LSC(L)
              IF (SVB(L).EQ.1.) THEN
                ADDLS=DXV(L)*0.5*(AH(L,K)+AH(LS,K))*DZC(K)*
     +          0.5*(HP(L) +HP(LS))*DYIV (L)
                  vol1=DXYP(L)*HP(L)*DZC(K)                   !hnr
                  vol2=DXYP(Ls)*HP(Ls)*DZC(K)                 !hnr
                  if (vol1.LT.vol2) vol2=vol1                   !hnr
                  ad=vol2/dt*adcoeff                                    !hnr
                  if (addls.GT.ad) then                         !hnr
                    addls=ad                                !hnr
                  end if                                        !hnr
                  if(ishdmf.lt.2) then                           !hnr
                    addls=0.                                     !hnr
                  end if                                         !hnr
              END IF
              IF (SVBO(L).EQ.1.) THEN
	        LWASP=LWASP+1
c                IF (IJCTLT(I,J-1).EQ.8) addls=0.0                      !hnr
                FLOW(LWASP)=VHDX(L,K)*DZC(K)
	        CRNU(LWASP)=2.*VHDX(L,K)*DYIV(L)*DXIV(L)/(HP(L)+HP(LS))
                BRINTT(LWASP)=ADDLS
              END IF
              IF (IJCTLT(I,J+1).EQ.8) THEN
                LN=LNC(L)
                IF (SVBO(LN).EQ.1.) THEN
  	          LWASP=LWASP+1
	          FLOW(LWASP)=VHDX(LN,K)*DZC(K)
                  CRNU(LWASP)=2.*VHDX(LN,K)*DYIV(LN)*DXIV(LN)
     +                        /(HP(L)+HP(LN))
	          BRINTT(LWASP)=addls                                !hnr
                END IF
              END IF
            END DO
          END DO

c Advection and dispersion in input flows
          DO K=KC,1,-1
            DO LT=1,nqsij
              I=Iqs(LT)
              J=Jqs(LT)
              IF(LIJLT(I,J).EQ.0) GOTO 310
              NS=NQSERQ(Lt)
              L=LQS(Lt)
                IF(flagwaspbC(ns,k).EQ.1) GOTO 310                       !HNR
              LWASP=LWASP+1
              FLOW(LWASP)=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
              CRNU(LWASP)=flow(LWASP)/DXp(L)/dyp(l)/(dzc(k)*HP(L))
c              BRINTT(LWASP)=dyp(l)*ah(l,k)*dzc(k)*hp(l)/dxp(l)
              BRINTT(LWASP)=0.0
310         END DO
          END DO
c  Advection and dispersion in structure flows
          DO K=KC,1,-1
            DO LT=1,nqctl
              Iu=Iqctlu(LT)
              Ju=Jqctlu(LT)
              Lu=Lij(iu,ju)
              Id=Iqctld(LT)
              Jd=Jqctld(LT)
              Ld=Lij(id,jd)
              IF(LU.EQ.0.AND.LD.EQ.0) GOTO 410
              flowx=RQCMUL(Lt)*QCTLT(K,Lt)
              UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*HP(L))
              IF(iu.eq.id) THEN
                addls=dxv(lu)*ah(lu,k)*dzc(k)*0.5*(hp(lu)
     $              +hp(ld))*dyiv(lu)
              ELSE
                addls=dyu(lu)*ah(lu,k)*dzc(k)*0.5*(hp(lu)
     $           +hp(ld))*dxiu(lu)
              END IF
                  if(ishdmf.lt.2) then                           !hnr
                    addls=0.                                     !hnr
                  end if                                         !hnr
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWX
	      CRNU(LWASP)=UDDXTMP
              BRINTT(LWASP)=ADDLS
410         END DO
          END DO
C
C Advection and dispersion in the Z-direction:
C
          IF (KC.GT.1) THEN
            DO K=KS,1,-1
              DO LT=2,LALT
                I=ILLT(LT)
                J=JLLT(LT)
                L=LIJ(I,J)
                addl=0.0
                ADDL1=ab(l,k)*hp(l)  !hnr Ev=ab*H
                if (addl1.gt.abwmx(l)) then
                  addl1=abwmx(l)         !hnr
                end if
                IF (SPB(L).EQ.1.) THEN
                  ADDL=DXYP(L)*addl1/hp(l)*DZIG(K)
                  vol1=DXYP(L)*HP(L)*DZC(K)                 !hnr
                  vol2=DXYP(L)*HP(L)*DZC(K+1)               !hnr
                  if (vol1.LT.vol2) vol2=vol1                 !hnr
                  ad=adcoeff*vol2/dt                                  !hnr
                  if (addl.GT.ad) then                        !hnr
                    addl=ad                               !hnr
                  end if                                      !hnr
                END IF
                IF (SWB(L).EQ.1) THEN
                  LWASP=LWASP+1
	          FLOW(LWASP)=-DXYP(L)*W(L,K)
                  CRNU(LWASP)=W(L,K)*DZIG(K)/HP(L)
	          BRINTT(LWASP)=ADDL
                END IF
              END DO
            END DO
          END IF
	  call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
c next only if temperature is transfered
          if(istran(2).GE.1) then
	    call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              stop
            end if
          end if
c next only if salinity is transfered
          if(istran(1).GE.1) then
	    call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	    if(ierror .gt. 0)then
              call hlgetlasterror(errstring)
              write(6,6000) ierror, errstring
              stop
            end if
          end if
c=======================================================================
	  call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	  if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
	  call hlmomentcomplete(Ihl_Handle,ierror)
          if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
          end if
        END IF
C      END INITIAL CONDITIONS WHEN IDAYS=0
       END IF
C     FINISH INITIALIZATION OF FILES (JSWASP=1)
      GOTO 3000
      END IF
	IBEGIN=IDAYS*NTSPTC
      IF(N.Lt.IBEGIN) GOTO 3000
C
C----------------------------------------------------------------------C
C
C **  WRITE TIME STEP DATA
C
C Advection and dispersion in the X-direction:
C
      LWASP=0
      DO K=KC,1,-1
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          ADDLW=0.0
          IF (SUB(L).EQ.1.) THEN
            LW=L-1
            ADDLW=DYU(L)*AHULPF(L,K)*DZC(K)*0.5*(HLPF(L) +HLPF(LW))
     +      *DXIU (L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)                           !hnr
            vol2=DXYP(LW)*HLPF(LW)*DZC(K)                         !hnr
            if (vol1.LT.vol2) vol2=vol1                           !hnr
            ad=adcoeff*vol2/dt                                            !hnr
            if (addlw.GT.ad) then                                 !hnr
              addlw=ad                                        !hnr
            end if                                                !hnr
                  if(ishdmf.lt.2) then                           !hnr
                    addlw=0.                                     !hnr
                  end if                                         !hnr
          END IF
          IF (SUBO(L).EQ.1.) THEN
            TMPVAL=UHLPF(L,K)+SVPT*UVPT(L,K)
            FLOWX=DYU(L)*TMPVAL*DZC(K)
            UDDXTMP=2.*TMPVAL*DXIU(L)/(HLPF(L)+HLPF(L-1))
            LWASP=LWASP+1
c            IF (IJCTLT(I-1,J).EQ.8) addlw=0.0                      !hnr
            FLOW(LWASP)=FLOWX
            CRNU(LWASP)=UDDXTMP
            BRINTT(LWASP)=ADDLW
          END IF
          IF (IJCTLT(I+1,J).EQ.8) THEN
            IF (SUBO(L+1).EQ.1.) THEN
              TMPVAL=UHLPF(L+1,K)+SVPT*UVPT(L+1,K)
              FLOWX=DYU(L+1)*TMPVAL*DZC(K)
              UDDXTMP=2.*TMPVAL*DXIU(L+1)/(HLPF(L+1)+HLPF(L))
              IPTMP=I+1
c              addlw=0.0                                         !hnr
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
	    BRINTT(LWASP)=ADDLW
            END IF
          END IF
        END DO
C
C Advection and dispersion in the Y-direction:
C
        DO LT=2,LALT
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          ADDLS=0.0
          IF (SVB(L).EQ.1.) THEN
            LS=LSC(L)
            ADDLS=DXV(L)*AHVLPF(L,K)*DZC(K)*0.5*(HLPF(L) +HLPF(LS))
     +      *DYIV (L)
            vol1=DXYP(L)*HLPF(L)*DZC(K)                          !hnr
            vol2=DXYP(LS)*HLPF(LS)*DZC(K)                        !hnr
            if (vol1.LT.vol2) vol2=vol1                          !hnr
            ad=adcoeff*vol2/dt                                           !hnr
            if (addls.GT.ad) then                                !hnr
              addls=ad                                       !hnr
            end if                                               !hnr
                  if(ishdmf.lt.2) then                           !hnr
                    addls=0.                                     !hnr
                  end if                                         !hnr
          END IF
          IF (SVBO(L).EQ.1.) THEN
            TMPVAL=VHLPF(L,K)+SVPT*VVPT(L,K)
            FLOWY=DXV(L)*TMPVAL*DZC(K)
c            FLOWY=DXV(L)*vlpf(l,k)*hlpf(l)*DZC(K)
            VDDYTMP=2.*TMPVAL*DYIV(L)/(HLPF(L)+HLPF(LSC(L)))
            JMTMP=J-1
c            IF(IJCTLT(I,J-1).EQ.8) addls=0.0                      !hnr
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWY
	    CRNU(LWASP)=VDDYTMP
	    BRINTT(LWASP)=ADDLS
          END IF
          IF (IJCTLT(I,J+1).EQ.8) THEN
            LN=LNC(L)
            IF (SVBO(LN).EQ.1.) THEN
              TMPVAL=VHLPF(LN,K)+SVPT*VVPT(LN,K)
              FLOWY=DXV(LN)*TMPVAL*DZC(K)
              VDDYTMP=2.*TMPVAL*DYIV(LN)/(HLPF(LN)+HLPF(L))
              JPTMP=J+1
c              addls=0.0                      !hnr
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWY
	    CRNU(LWASP)=VDDYTMP
	    BRINTT(LWASP)=ADDLS
            END IF
          END IF
        END DO
      END DO

c Advection and dispersion in input flows
      DO K=KC,1,-1
        DO LT=1,nqsij
          I=Iqs(LT)
          J=Jqs(LT)
          IF(LIJLT(I,J).EQ.0) GOTO 300
          NS=NQSERQ(Lt)
          L=LQS(Lt)
            IF(flagwaspbC(ns,k).EQ.1) GOTO 300                       !HNR
          flowx=RQSMUL(Lt)*(QSS(K,Lt)+qfactor(lt)*QSERT(K,NS))
          UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*HmP(L))
c          addlw=dyp(l)*ahulpf(l,k)*dzc(k)*hlpf(l)/dxp(l)  !hnr
          addlw=0.0                                        !hnr
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
	    BRINTT(LWASP)=ADDLW
300     END DO
      END DO
c  Advection and dispersion in structure flows
      DO K=KC,1,-1
        DO LT=1,nqctl
          Iu=Iqctlu(LT)
          Ju=Jqctlu(LT)
          Lu=Lij(iu,ju)
          Id=Iqctld(LT)
          Jd=Jqctld(LT)
          Ld=Lij(id,jd)
          IF(LU.EQ.0.AND.LD.EQ.0) GOTO 400
          flowx=RQCMUL(Lt)*QCTLT(K,Lt)
          UDDXTMP=flowx/DXp(L)/dyp(l)/(dzc(k)*HmP(L))
          IF(iu.eq.id) THEN
            addls=dxv(lu)*ahvlpf(lu,k)*dzc(k)*0.5*(hlpf(lu)
     $            +hlpf(ld))*dyiv(lu)
          ELSE
            addls=dyu(lu)*ahulpf(lu,k)*dzc(k)*0.5*(hlpf(lu)
     $           +hlpf(ld))*dxiu(lu)
          END IF
                  if(ishdmf.lt.2) then                           !hnr
                    addls=0.                                     !hnr
                  end if                                         !hnr
	    LWASP=LWASP+1
	    FLOW(LWASP)=FLOWX
	    CRNU(LWASP)=UDDXTMP
	    BRINTT(LWASP)=ADDLs
400     END DO
      END DO
C
C Advection and dispersion in the Z-direction:
C
      IF (KC.GT.1) THEN
        DO K=KS,1,-1
          DO LT=2,LALT
            I=ILLT(LT)
            J=JLLT(LT)
            L=LIJ(I,J)
            ADDL=0.0
            addl1=ablpf(l,k)*hlpf(l)  !hnr  Ev=ab*hp
            if(addl1.gt.abwmx(l)) then
              addl1=abwmx(l)   !hnr
            end if
            IF (SPB(L).EQ.1.) THEN
              ADDL=DXYP(L)*Addl1/HLPF(L)*DZIG(K)
              vol1=DXYP(L)*HLPF(L)*DZC(K)           !hnr
              vol2=DXYP(L)*HLPF(L)*DZC(K+1)         !hnr
              if (vol1.LT.vol2) vol2=vol1           !hnr
              ad=adcoeff*vol2/dt                            !hnr
              if (addl.GT.ad) then                  !hnr
                addl=ad                         !hnr
              end if                                !hnr
            END IF
            IF (SWB(L).EQ.1) THEN
              TMPVAL=WLPF(L,K)+SVPT*WVPT(L,K)
              FLOWZ=-DXYP(L)*TMPVAL
              WDDZTMP=TMPVAL*DZIG(K)/HLPF(L)
              LWASP=LWASP+1
	      FLOW(LWASP)=FLOWZ
	      CRNU(LWASP)=WDDZTMP
	      BRINTT(LWASP)=ADDL
            END IF
          END DO
        END DO
      END IF




C
C Segment Properties:
C
      LCELTMP=0
      DO K=KC,1,-1
        DO LT=2,LALT
          LCELTMP=LCELTMP+1
          I=ILLT(LT)
          J=JLLT(LT)
          L=LIJ(I,J)
          LN=LNC(L)
          VOLUM=DXYP(L)*HLPF(L)*DZC(K)
          IF(NTSMMT.LT.NTSPTC) VOLUM=DXYP(L)*HP(L)*DZC(K)
          DEPTH=HLPF(L)*DZC(K)
          VELX=0.5*(UHLPF(L,K)+SVPT*UVPT(L,K) +UHLPF(L+1,K)+SVPT*UVPT
     +    (L+1,K))/HLPF(L)
          VELY=0.5*(VHLPF(L,K)+SVPT*VVPT(L,K) +VHLPF(LN,K)+SVPT*VVPT
     +    (LN,K) )/HLPF(L)
          VELZ=0.5*(WLPF(L,K-1)+SVPT*WVPT(L,K-1) +WLPF(L,K)+SVPT*WVPT
     +    (L,K) )
          VELMAG=SQRT(VELX*VELX+VELY*VELY+VELZ*VELZ)
            SegSalt(LCELTMP)=SALLPF(L,K)
            SegTemp(LCELTMP)=TEMLPF(L,K)
	    SegVolume(LCELTMP)=VOLUM
	    SegDepth(LCELTMP)=DEPTH
	    SegVel(LCELTMP)=VELMAG
        END DO
      END DO
C
	   call hlsetseginfo(Ihl_handle,1,SegVolume,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetseginfo(Ihl_handle,2,SegDepth,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetseginfo(Ihl_handle,3,SegVel,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
c next only if temperature is transfered
       if(istran(2).GE.1) then
	   call hlsetseginfo(Ihl_handle,4,SegTemp,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
       end if
c next only if salinity is transfered
       if(istran(1).GE.1) then
	   call hlsetseginfo(Ihl_handle,5,SegSalt,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
       end if
c=======================================================================
	   call hlsetflowinfo(Ihl_handle,1,Flow,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,2,crnu,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlsetflowinfo(Ihl_handle,3,brintt,ierror)
	   if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if
	   call hlmomentcomplete(Ihl_Handle,ierror)
         if(ierror .gt. 0)then
            call hlgetlasterror(errstring)
            write(6,6000) ierror, errstring
            stop
         end if

C
3000  JSWASP=0
6000  format('Error ',I10, ' : ', A)
9966  FORMAT(f12.4,I6,f12.4)
9977  FORMAT(3I6,6f12.4)
C----------------------------------------------------------------------c
C ** end SUBROUTINE WASPHYDROLINK
C----------------------------------------------------------------------c
      RETURN
      END