SUBROUTINE CALSTEPD
 
  ! CHANGE RECORD                                                                                                         
  ! **  SUBROUTINE CALSTEP ESTIMATE THE CURRENT MAXIMUM TIME STEP SIZE                                                    
  ! **  FORM LINEAR STABILITY CRITERIA AND A FACTOR OF SAFETY   
  !                                                          
  !----------------------------------------------------------------------C  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2011-03       Paul M. Craig     Rewritten to F90 and added OMP
 
  USE GLOBAL
  
  IMPLICIT NONE

  INTEGER :: ND,L,K,LF,LL,LE,LN,LP,KM,LLOC,ITRNTMP,NX
  INTEGER :: NMD,LMDCHHT,LMDCHUT,LMDCHVT,ITMPR,LLOCOLD
  !INTEGER :: NTMP,NTCTMP,LSE,LNW
  INTEGER,SAVE :: NUP
  
  REAL   :: QUKTMP,QVKTMP,RTMPR,TESTTEMP,DTMAXX
  REAL   :: TMPUUU,DTTMP,TMPVVV,TMPVAL
  !REAL   :: UEAST,UWEST,VSOUTH,VNORTH,VATUUU,UATVVV
  REAL   :: TOP,QXPLUS,QYPLUS,QZPLUS,QXMINS,QYMINS,QZMINS,QTOTAL,QSRC,BOT
  REAL   :: DTCOMP,DTWARN,DTDYNP
  !REAL   :: TAUBC,UCTR,VCTR,UHMAG,FRIFRE,FRIFRE2,ACACTMP
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL2
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL3
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTL4
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QSUBINN
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: QSUBOUT
  REAL,SAVE,ALLOCATABLE,DIMENSION(:) :: DTHISTORY

  IF(  .NOT. ALLOCATED(DTL1) )THEN
    ALLOCATE(DTL1(LCM))
    ALLOCATE(DTL2(LCM))
    ALLOCATE(DTL3(LCM))
    ALLOCATE(DTL4(LCM))
    ALLOCATE(QSUBINN(LCM,KCM))
    ALLOCATE(QSUBOUT(LCM,KCM))
    ALLOCATE(DTHISTORY(10))
    
    DTL1=0.0
    DTL2=0.0
    DTL3=0.0
    DTL4=0.0
    QSUBINN=0.0
    QSUBOUT=0.0
    DTHISTORY=DT
    NUP=0
  ENDIF

  !**********************************************************************C                                                
  ITRNTMP=0
  DO NX=1,7
    ITRNTMP=ITRNTMP+ISTRAN(NX)
  ENDDO

  IF( N <= 1 ) DTDYN = DT

  DTMIN = DT
  DTMAXX = TIDALP

  ! **  DETERMINE SOURCE/SINKS FOR SUBGRID SCALE CHANNEL EXCHANGES                                                        
  IF( MDCHH >= 1 )THEN
    DO K=1,KC
      DO L=2,LA
        QSUBOUT(L,K)=0.0
        QSUBINN(L,K)=0.0
      ENDDO
    ENDDO

    DO K=1,KC
      DO NMD=1,MDCHH
        IF( LKSZ(L,K) )CYCLE
        LMDCHHT=LMDCHH(NMD)
        LMDCHUT=LMDCHU(NMD)
        LMDCHVT=LMDCHV(NMD)
        IF( MDCHTYP(NMD) == 1 )THEN
          QUKTMP=QCHANU(NMD)*DZC(LMDCHUT,K)
          QVKTMP=0.
        ENDIF
        IF( MDCHTYP(NMD) == 2 )THEN
          QVKTMP=QCHANV(NMD)*DZC(LMDCHVT,K)
          QUKTMP=0.
        ENDIF
        IF( MDCHTYP(NMD) == 3 )THEN
          QUKTMP=QCHANU(NMD)*DZC(LMDCHUT,K)
          QVKTMP=QCHANV(NMD)*DZC(LMDCHVT,K)
        ENDIF
        QSUBOUT(LMDCHHT,K)=QSUBOUT(LMDCHHT,K) +MIN(QUKTMP,0.) +MIN(QVKTMP,0.)
        QSUBINN(LMDCHHT,K)=QSUBINN(LMDCHHT,K) +MAX(QUKTMP,0.) +MAX(QVKTMP,0.)
        QSUBOUT(LMDCHUT,K)=QSUBOUT(LMDCHUT,K) -MAX(QUKTMP,0.)
        QSUBINN(LMDCHUT,K)=QSUBINN(LMDCHUT,K) -MIN(QUKTMP,0.)
        QSUBOUT(LMDCHVT,K)=QSUBOUT(LMDCHVT,K) -MAX(QVKTMP,0.)
        QSUBINN(LMDCHVT,K)=QSUBINN(LMDCHVT,K) -MIN(QVKTMP,0.)
      ENDDO
    ENDDO
  ENDIF

  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,L,K,LF,LL,LP,LE,LN,KM) &
  !$OMP  PRIVATE(TMPUUU,DTTMP,TMPVVV,TMPVAL,TESTTEMP) &
  !$OMP  PRIVATE(TOP,QXPLUS,QYPLUS,QZPLUS,QXMINS,QYMINS,QZMINS,QTOTAL,QSRC,BOT)  SCHEDULE(STATIC,1)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)

    DO L=LF,LL
      DTL1(L) = DTMAXX
      DTL2(L) = DTMAXX
      !DTL3(L) = DTMAXX
    ENDDO

    IF( DTSSDHDT > 0. )THEN
      DO L=LF,LL
        DTL4(L) = DTMAXX
      ENDDO
    ENDIF
    
    ! **  METHOD 1: COURANT–FRIEDRICHS–LEWY
    DO LP=1,LLWET(KC,ND)
      L=LKWET(LP,KC,ND)  
      TMPVAL = 1.E-16
      TMPUUU = ABS(UHE(L))*HUI(L)*DXIU(L)
      TMPVVV = ABS(VHE(L))*HVI(L)*DYIV(L)
      TMPVAL = TMPUUU + TMPVVV
      IF( TMPVAL > 0.  ) DTL1(L)  = 1./TMPVAL
    ENDDO

    ! **  METHOD 2: POSITIVITY OF ADVECTED MATERIAL, DTL2                                                                   
    IF( ITRNTMP >= 1 )THEN
      DO K=1,KC
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          LE=LEC(L)
          LN=LNC(L)
          KM=K-1
          TOP    = H1P(L)*DXYP(L)*DZC(L,K)
          QXPLUS = UHDY2(LE,K)
          QXPLUS = MAX(QXPLUS,0.0)
          QYPLUS = VHDX2(LN,K)
          QYPLUS = MAX(QYPLUS,0.0)
          QZPLUS = W2(L,K)*DXYP(L)
          QZPLUS = MAX(QZPLUS,0.0)
          QXMINS = UHDY2(L,K)
          QXMINS = -MIN(QXMINS,0.0)
          QYMINS = VHDX2(L,K)
          QYMINS = -MIN(QYMINS,0.0)
          QZMINS = W2(L,KM)*DXYP(L)
          QZMINS = -MIN(QZMINS,0.0)
          QTOTAL = QSUM(L,K) + QSUBOUT(L,K) + QSUBINN(L,K)
          QSRC   = -MIN(QTOTAL,0.0)
          BOT = QXPLUS+QYPLUS+QZPLUS+QXMINS+QYMINS+QZMINS+QSRC
          IF( BOT > 1.E-12 )THEN
            DTTMP = TOP/BOT
            DTL2(L) = MIN(DTL2(L),DTTMP)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    ! *** METHOD 3: IMPLICIT BOTTOM FRICTION AND ROTATIONAL ACCELERATION DAMPING                                            
    !DO L=LF,LL
    !  IF( LACTIVE(L) )THEN
    !    TMPVAL=SUB(L)+SUB(LEC(L))+SVB(L)+SVB(LNC(L))
    !    IF( TMPVAL > 0.5 )THEN
    !      LN=LNC(L)
    !      TAUBC=QQ(L,0)/CTURB2
    !      UCTR=0.5*(U(L,1)+U(LEC(L),1))
    !      VCTR=0.5*(V(L,1)+V(LN,1))
    !      UHMAG=HP(L)*SQRT(UCTR*UCTR+VCTR*VCTR)
    !      IF( UHMAG > 0.0 )THEN
    !        FRIFRE=TAUBC/UHMAG
    !        FRIFRE2=FRIFRE*FRIFRE
    !        ACACTMP=(CAC(L,KC)*HPI(L)*DXYIP(L))**2
    !        IF( ACACTMP > FRIFRE2 )THEN
    !          DTTMP=2.*FRIFRE/(ACACTMP-FRIFRE2)
    !          DTL3(L)=MIN(DTL3(L),DTTMP)
    !        ENDIF
    !      ENDIF
    !    ENDIF
    !  ENDIF
    !ENDDO

    ! ***  METHOD 4: LIMIT RATE OF DEPTH CHANGE
    IF( DTSSDHDT > 0. )THEN
      DO LP=1,LAWET
        L=LWET(LP) 
        IF( HP(L) < HDRY/10.)  HP(L)=HDRY/10.
        IF( H1P(L) < HDRY/10.) H1P(L)=HDRY/10.
        TESTTEMP = MAX(ABS(HP(L)-H1P(L)),1.E-06)
        TMPVAL   = DTDYN*HP(L)/TESTTEMP
        DTL4(L)  = DTSSDHDT*TMPVAL
      ENDDO
    ENDIF
  ENDDO   ! *** END OF DOMAIN
  !$OMP END PARALLEL DO
  
  ! **  CHOOSE THE MINIMUM OF THE THREE METHODS                                                                           
  DTL1MN = 2.*DTMAXX
  DO L=2,LA
    IF( DTL1MN > DTL1(L) )THEN
      DTL1MN = DTL1(L)
      L1LOC = L
    ENDIF
  ENDDO
  
  DTL2MN = 2.*DTMAXX
  IF( ITRNTMP >= 1 )THEN
    DO L=2,LA
      IF( DTL2MN > DTL2(L) )THEN
        DTL2MN = DTL2(L)
        L2LOC = L
      ENDIF
    ENDDO
  ENDIF
  
  !DTL3MN = 2.*DTMAXX
  !DO L=2,LA
  !  IF( DTL3MN > DTL3(L) )THEN
  !    DTL3MN = DTL3(L)
  !    L3LOC = L
  !  ENDIF
  !ENDDO
  
  DTL4MN = 2.*DTMAXX
  IF( DTSSDHDT > 0. )THEN
    DO L=2,LA
      IF( DTL4MN > DTL4(L) )THEN
        DTL4MN = DTL4(L)
        L4LOC = L
      ENDIF
    ENDDO
  ENDIF

  ! **  FIND MINIMUM & APPLY A SAFETY FACTOR                                                                              
  DTL1MN = DTL1MN*DTSSFAC
  DTTMP  = 2.*DTMAXX
  IF( DTTMP > DTL1MN )THEN
    DTTMP  = DTL1MN
    DTCOMP = DTTMP/DTSSFAC
    LLOC   = L1LOC
  ENDIF
  
  DTL2MN = DTL2MN*0.5     ! *** FIXED SAFETY FACTOR FOR ADVECTION OF 0.5
  IF( DTTMP > DTL2MN )THEN
    DTTMP  = DTL2MN
    DTCOMP = 2.*DTTMP
    LLOC   = L2LOC
  ENDIF
  
  !DTL3MN = DTL3MN*DTSSFAC
  !IF( DTTMP > DTL3MN )THEN
  !  DTTMP  = DTL3MN
  !  DTCOMP = DTTMP/DTSSFAC
  !  LLOC   = L3LOC
  !ENDIF
  
  IF( DTSSDHDT > 0. )THEN
    DTL4MN = DTL4MN*0.5     ! *** FIXED SAFETY FACTOR FOR ADVECTION OF 0.5
    IF( DTTMP > DTL4MN )THEN
      DTTMP  = DTL4MN
      DTCOMP = 2.*DTTMP
      LLOC   = L4LOC
    ENDIF
  ENDIF
  LLOCOLD = LMINSTEP
  LMINSTEP = LLOC

  IF( DTCOMP < DTMIN )THEN   ! *** DS-INTL SINGLE LINE
    OPEN(8,FILE=OUTDIR//'EFDCLOG.OUT',POSITION='APPEND')  
    WRITE(8,800)TIMEDAY,DTTMP,DTMIN,IL(LLOC),JL(LLOC),HP(LLOC)
    WRITE(6,800)TIMEDAY,DTTMP,DTMIN,IL(LLOC),JL(LLOC),HP(LLOC)
    WRITE(8,801)IL(L1LOC),JL(L1LOC),DTL1MN,HP(L1LOC)
    WRITE(6,801)IL(L1LOC),JL(L1LOC),DTL1MN,HP(L1LOC)
    WRITE(8,802)IL(L2LOC),JL(L2LOC),DTL2MN,HP(L2LOC)
    WRITE(6,802)IL(L2LOC),JL(L2LOC),DTL2MN,HP(L2LOC)
    !WRITE(8,803)IL(L3LOC),JL(L3LOC),DTL3MN,HP(L3LOC)
    !WRITE(6,803)IL(L3LOC),JL(L3LOC),DTL3MN,HP(L3LOC)
    IF( DTSSDHDT > 0. )THEN
      WRITE(8,804)IL(L4LOC),JL(L4LOC),DTL4MN
      WRITE(6,804)IL(L4LOC),JL(L4LOC),DTL4MN
    ENDIF
    CLOSE(8)
  
    DTTMP=DTMIN
  
  ELSEIF( DTTMP < DTMIN )THEN
    DTWARN=DTTMP
    DTTMP=DTMIN
  ELSE
    TMPVAL=DTTMP/DTMIN
    ITMPR=NINT(TMPVAL)
    RTMPR=FLOAT(ITMPR)
    IF( RTMPR < TMPVAL )THEN
      DTTMP=RTMPR*DTMIN
    ELSE
      DTTMP=(RTMPR-1.)*DTMIN
    ENDIF
  ENDIF

  ! **  SET TO MINIMUM TIME STEP ON STARTUP                                                                               
  IF( N < 2 ) DTTMP = DTMIN

  ! *** RESTRICT INCREASE IN TIME STEP TO DTMIN                                                                           
  !DTDYNP = DTHISTORY(10) + DTMIN
  DTDYNP = DTDYN + DTMIN
  IF( DTTMP > DTDYNP )THEN
    NUP = NUP + 1
    IF( NUP >= NUPSTEP )THEN
      DTTMP = DTDYNP
      NUP = 0
    ELSE
      DTTMP = DTDYN
    ENDIF
  ELSEIF( DTTMP < DTDYN )THEN
    IF( NUPSTEP >= 25 .AND. TIMEDAY > TBEGIN+1. )THEN
      ! *** REPORT WHEN THERE IS A DECREASE
      OPEN(9,FILE=OUTDIR//'TIME.LOG',POSITION='APPEND')  
      WRITE(9,'(A,I10,F14.5,2(A,F10.4,I7,F8.3))') 'DYNAMIC STEP DOWN',NITER,TIMEDAY,'     [DT,L,HP]   OLD: ',DTDYN,LLOCOLD,H1P(LLOCOLD),'     NEW: ',DTTMP,LMINSTEP,HP(LMINSTEP)
      CLOSE(9)
    ENDIF
    
    NUP = 0
  ELSE
    DTTMP = DTDYN
  ENDIF
  DTDYN = MIN(DTTMP,DTMAX)

  ! **  SET INCREMENTAL INCREASE IN OUTPUT COUNTER                                                                        
  NINCRMT = NINT(DTDYN/DTMIN)
  DTDYN   = FLOAT(NINCRMT)*DTMIN
  
    100 FORMAT(5I5,5F12.5,E13.5)                                                                                          
    101 FORMAT(3I5,E13.5)                                                                                                 
    800 FORMAT('  TIME,DTDYN,DTMIN,I,J,HP = ',F12.5,2E12.4,2I7,F10.3)
    801 FORMAT('  MOM  ADV,I,J,DT,HP = ',2I5,E13.4,F10.3)                                                                         
    802 FORMAT('  MASS ADV,I,J,DT,HP = ',2I5,E13.4,F10.3)                                                                         
    803 FORMAT('  CURV ACC,I,J,DT,HP = ',2I5,E13.4,F10.3)                                                                         
    804 FORMAT('  LIM DHDT,I,J,DT,HP = ',2I5,E13.4,F10.3) 
    880 FORMAT(3I5,8E13.4)                                                                                                
   8899 FORMAT(' DT3 ERROR ',2I5,6E13.5)                                                                                  

  !DO K=10,2,-1
  !  DTHISTORY(K) = DTHISTORY(K-1)
  !ENDDO
  !DTHISTORY(1) = DTDYN

  !**********************************************************************C                                                
  RETURN
END
