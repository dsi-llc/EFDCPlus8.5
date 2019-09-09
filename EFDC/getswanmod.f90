MODULE GETSWANMOD
  ! ** PURPOSE:
  ! ** DIRECTLY READ SWAN OUTPUT CONTAINED IN
  ! ** FRM/GRP FILE OR TBL FILE CORRESPONDING TO LOC FILE
  ! ** SWN FILE: SET INRHOG 1 TO GET DISSIPATION IN W/M2
  ! ** IFWAVE =1   :C14A
  ! ** SWANGRP=1/0 :C14A
  ! ** OTHER PARAM :C14B
  ! ** AUTHOR: DANG HUU CHUNG

  USE GLOBAL 
  USE INFOMOD
  USE WAVELENGTH
  USE XYIJCONV

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE GETSWAN_GRP
  ! ** READ OUTPUT BY SWAN FOR THE WHOLE EFDC GRID
  ! ** AND THE TERMS IN SWN FILE DETERMINED BY:
  ! $************ OUTPUT REQUESTS *************************
  ! ** GROUP  'GRP1' SUBGrid 0 IC-1  0 JC-1    
  ! ** TABLE  'GRP1' HEADer 'EXAMPLE.grp' HS PDIR RTP WLEN DEPTH DISSIP
  !
  ! ** OUTPUT:
  ! ** GENERATE WAVE PARAMETERS FOR EVERY CELL L=2:LA OF EFDC

  CHARACTER(200) :: STR,TERM
  INTEGER :: NL,I,J,L,ISO,NW
  REAL(RKD)   :: VAL(20),WDIR,WPRD,EDIS
  REAL(RKD)   :: WVDX,WVDY,WVCX,WVCY

  IF( JSWAVE == 0 )THEN
  !** FIRST CALL
    WRITE(*,'(A)')'READING SWAN_GRP.INP'
    OPEN(UGRP,FILE='swan_grp.inp',ACTION='READ')  !FRMFILE
    DO NL=1,5
      READ(UGRP,'(A)',IOSTAT=ISO) STR
      IF( ISO > 0 ) STOP 'SWAN_GRP.INP: READING ERROR'  
    ENDDO
    STR  = ADJUSTL(STR)
    TERM = ADJUSTL(STR(2:))
    STR=READSTR(UGRP)
    NCOL = NUMCOL(STR)

    NHS = FINDSTR(TERM,'Hsi',NCOL)   ! *** INCIDENT OFFSHORE SIGNIFICANT WAVE HEIGHT (m)
    NPK = FINDSTR(TERM,'PkD',NCOL)   ! *** WAVE ANGLE (DEG)
    NRT = FINDSTR(TERM,'RTp',NCOL)   ! *** TIME TO PEAK (SEC)
    NWL = FINDSTR(TERM,'Wle',NCOL)   ! *** WAVE LENGTH (M)
    NHP = FINDSTR(TERM,'Dep',NCOL)   ! *** WATER DEPTH (M)
    NDI = FINDSTR(TERM,'Dis',NCOL)   ! *** WAVE DISSIPATION (M)

    ! ** READ BUFFER DATA 
    DO NW=1,IWVCOUNT-1
      DO J=1,JC
        DO I=1,IC
          READ(UGRP,*,IOSTAT=ISO) VAL(1:NCOL)
          IF( ISO > 0 ) STOP 'SWAN_GRP.INP: READING ERROR' 
        ENDDO
      ENDDO
    ENDDO
  ENDIF

  NWVCELLS = 0
  LOOPJ: DO J=1,JC
    DO I=1,IC
      READ(UGRP,*,IOSTAT=ISO) VAL(1:NCOL)   
      IF( ISO > 0 ) STOP 'SWAN_GRP.INP: READING ERROR'  
  
      L = LIJ(I,J)    
      IF( L <= 1 .OR. L > LA )CYCLE
    
      NWVCELLS = NWVCELLS+1
      LWVCELL(NWVCELLS) = L   ! *** PMC - 2015-12 - This was pointing to L-1 but LWVCELL is used as an L array reference.  This has to be L not L-1.
      LWVMASK(L) = .TRUE.
      
      ! *** WAVE HEIGHT (M), THEN ZERO WAVE IMPACTS FOR VEGETATION
      WV(L).HEISIG = VAL(NHS)
      IF( MVEGL(L) /= MVEGOW ) THEN
        WV(L).HEISIG = 0.
      ELSE
        IF( (MVEGL(LWC(L)) /= MVEGOW) .OR. &
            (MVEGL(LEC(L)) /= MVEGOW) .OR. &
            (MVEGL(LSC(L)) /= MVEGOW) .OR. &
            (MVEGL(LNC(L)) /= MVEGOW)      )  WV(L).HEISIG = 0.5*WV(L).HEISIG
      ENDIF
      
      WDIR    = VAL(NPK)                                ! *** WAVE ANGLE = (TRUE EAST,WAVE) ANTI-CLOCKWISE [0,360] [D]
      WPRD    = VAL(NRT)                                ! *** WAVE PERIOD [S] 
      WV(L).LENGTH  = VAL(NWL)                          ! *** WAVE LENGTH [M] 
      EDIS    = VAL(NDI)                                ! *** WAVE DISSIPATION  [ W/m2=Kg/s3]
      !WV(L).DISSIPA(KC) = EDIS/RHO                     ! *** [M3/S3]   
      
      IF( WV(L).HEISIG >= WHMI .AND. WVPRD > 0. )THEN
        WVDX= COS(WDIR*PI/180)
        WVDY= SIN(WDIR*PI/180)
        WVCX =  CVN(L)*WVDX - CVE(L)*WVDY               ! *** TO LOCAL CURVI-LINEAR SYSTEM
        WVCY = -CUN(L)*WVDX + CUE(L)*WVDY
        WV(L).DIR  = ATAN2(WVCY,WVCX)                   ! *** CELL-EAST,WAVE (COUNTER-CLOCKWISE [-pi,pi])
        WV(L).FREQ = 2.*PI/WPRD
        WV(L).DISSIPA(KC) = 0.25*g*WV(L).HEISIG**2/WPRD ! *** Energy dissipation due to breaking wave  
      ELSE
        WV(L).DIR  = 0.
        WV(L).LENGTH     = 0.
        WV(L).FREQ  = 1.
        WV(L).HEISIG   = 0.
        WV(L).DISSIPA(KC)=0
      ENDIF
    ENDDO
  ENDDO LOOPJ
  WV(2:LA).HEIGHT = WV(2:LA).HEISIG  ! *** ASSIGN INCIDENT SIGNIFICANT WAVE HEIGHT TO APPLIED WAVE HEIGHT
  
  IF( IWVCOUNT == NWVTIM )THEN
    CLOSE(UGRP)
    WRITE(*,'(A)') '** WAVE: END OF WAVE RECORD'
  ENDIF 
  END SUBROUTINE

  SUBROUTINE GETSWAN_LOC
  ! ** READ OUTPUT BY SWAN BASED ON A LOC FILE
  ! ** AND THE TERMS IN SWN FILE DETERMINED BY:
  ! 
  ! ** POINTS 'loc' FILE   'EXAMPLE.loc'
  ! ** TABLE  'loc' HEAD   'EXAMPLE.tbl'  HS PDIR RTP WLEN DEPTH DISSIP 
  !
  ! ** OUTPUT:
  ! ** WAVE PARAMETERS FOR THE CELLS OF EFDC

  CHARACTER(200) :: STR,TERM
  INTEGER :: N,L,ISO,NP,NW
  REAL(RKD)   :: VAL(20),WDIR,WPRD,EDIS
  REAL(RKD)   :: WVDX,WVDY,WVCX,WVCY

  IF( JSWAVE == 0 )THEN
    !** FIRST CALL
    WRITE(*,'(A)')'READING SWAN_LOC.INP'
    OPEN(ULOC,FILE='swan_loc.inp',ACTION='READ')
    STR=READSTR(ULOC)  
    NP=0
    DO WHILE(1)
      READ(ULOC,'(A)',END=250) STR
      NP=NP+1
    ENDDO
    250 REWIND(ULOC) 
    NLOC = NP
    ALLOCATE(SWNLOC%ICEL(NLOC),SWNLOC%JCEL(NLOC))
    ALLOCATE(SWNLOC%XCEL(NLOC),SWNLOC%YCEL(NLOC))

    STR=READSTR(ULOC)  
    DO N=1,NLOC
      READ(ULOC,*) SWNLOC%XCEL(N),SWNLOC%YCEL(N)
    ENDDO
    CLOSE(ULOC)
    CALL XY2IJ(SWNLOC)  

    WRITE(*,'(A)')'READING SWAN_TBL.INP'
    OPEN(UTBL,FILE='swan_tbl.inp',ACTION='READ')
    DO N=1,5
      READ(UTBL,'(A)',IOSTAT=ISO) STR
      IF( ISO > 0 )THEN
        STOP 'SWAN_TBL.INP: READING ERROR'  
      ELSEIF( ISO < 0 )THEN
        PRINT*,'END OF SWAN_TBL.INP'
      ENDIF
    ENDDO
    STR  = ADJUSTL(STR)
    TERM = ADJUSTL(STR(2:))
    STR=READSTR(UTBL)   !NUMBER
    NCOL = NUMCOL(STR)
    
    !%       Hsig          PkDir         RTpeak        Wlen          Depth         Dissip        Genera        Redist        Radstr   
    !%       [m]           [degr]        [sec]         [m]           [m]           [m2/s]        [m2/s]        [m2/s]        [m2/s]   
    
    NHS = FINDSTR(TERM,'Hsi',NCOL)   ! *** INCIDENT OFFSHORE SIGNIFICANT WAVE HEIGHT (m)
    NPK = FINDSTR(TERM,'PkD',NCOL)   ! *** WAVE ANGLE (DEG)
    NRT = FINDSTR(TERM,'RTp',NCOL)   ! *** TIME TO PEAK (SEC)
    NWL = FINDSTR(TERM,'Wle',NCOL)   ! *** WAVE LENGTH (M)
    NHP = FINDSTR(TERM,'Dep',NCOL)   ! *** WATER DEPTH (M)
    NDI = FINDSTR(TERM,'Dis',NCOL)   ! *** WAVE DISSIPATION (M)
    
    ! ** READ BUFFER DATA 
    DO NW=1,IWVCOUNT-1
      DO N=1,NLOC
        READ(UTBL,*,IOSTAT=ISO) VAL(1:NCOL)
        IF( ISO > 0 ) STOP 'SWAN_TBL.INP: READING ERROR' 
      ENDDO
    ENDDO
  ENDIF

  NWVCELLS = 0
  DO N=1,NLOC
    READ(UTBL,*,IOSTAT=ISO) VAL(1:NCOL) 
    IF( ISO > 0 ) STOP 'SWAN_TBL.INP: READING ERROR'  
  
    L= LIJ(SWNLOC%ICEL(N),SWNLOC%JCEL(N))
    IF( L <= 1 .OR. L>LA) CYCLE
 
    NWVCELLS = NWVCELLS+1
    LWVCELL(NWVCELLS) = L   ! *** PMC - 2015-12 - This was pointing to L-1 but LWVCELL is used as an L array reference.  This has to be L not L-1.
    LWVMASK(L) = .TRUE.
    WV(L).HEISIG = VAL(NHS)                          ! *** SIGNIFICANT WAVE HEIGHT (M)
    
    ! *** VEGETATION EFFECT
    IF( ISVEG > 0 )THEN
      IF( MVEGL(L) /= MVEGOW ) THEN
        WV(L).HEIGHT = 0.
      ELSE
        IF( (MVEGL(LWC(L)) /= MVEGOW) .OR. &
            (MVEGL(LEC(L)) /= MVEGOW) .OR. &
            (MVEGL(LSC(L)) /= MVEGOW) .OR. &
            (MVEGL(LNC(L)) /= MVEGOW)      )  WV(L).HEIGHT=0.5*WV(L).HEIGHT
      ENDIF
    ENDIF

    WDIR    = VAL(NPK)                          ![D] WAVE ANGLE = (TRUE EAST,WAVE) ANTI-CLOCKWISE [0,360]
    WPRD    = VAL(NRT)                          ![S] PERIOD
    WV(L).LENGTH  = VAL(NWL)                          ![M] W LENGTH
    EDIS    = VAL(NDI)                          ! W/m2=Kg/s3: WAVE DISSIPATION
    !WV(L).DISSIPA(KC) = EDIS/RHO                    ![M3/S3] 
    IF( WV(L).HEISIG >= WHMI .AND. WVPRD > 0 )THEN 
      WVDX= COS(WDIR*PI/180)
      WVDY= SIN(WDIR*PI/180)
      WVCX =  CVN(L)*WVDX - CVE(L)*WVDY         ! TO LOCAL CURVI-LINEAR SYSTEM
      WVCY = -CUN(L)*WVDX + CUE(L)*WVDY
      WV(L).DIR= ATAN2(WVCY,WVCX)               ! CELL-EAST,WAVE (COUNTER-CLOCKWISE [-pi,pi])
      WV(L).FREQ=2.*PI/WPRD
      WV(L).DISSIPA(KC) = 0.25*g*WV(L).HEISIG**2/WPRD    !Energy dissipation due to breaking wave 
    ELSE
      WV(L).DIR  = 0.
      WV(L).LENGTH     = 0.
      WV(L).FREQ  = 1.
      WV(L).HEISIG   = 0.
      WV(L).DISSIPA(KC)=0
    ENDIF    
  ENDDO
  
  WV(2:LA).HEIGHT = WV(2:LA).HEISIG  ! *** ASSIGN INCIDENT SIGNIFICANT WAVE HEIGHT TO APPLIED WAVE HEIGHT
  
  IF( IWVCOUNT == NWVTIM )THEN
    CLOSE(UTBL)
    WRITE(*,'(A)') '** WAVE: END OF WAVE RECORD'
  ENDIF 
  END SUBROUTINE

  SUBROUTINE STRPROC(STR,TERM,N)
  CHARACTER(*),INTENT(IN) :: STR
  CHARACTER(*),INTENT(OUT) :: TERM(:)
  INTEGER(IK4),INTENT(OUT) :: N
  INTEGER(IK4) :: STRL,I
  CHARACTER(120) :: SSTR

  SSTR = ADJUSTL(STR)
  STRL = LEN_TRIM(SSTR)
  SSTR = SSTR(2:STRL)
  N=0
  DO WHILE(1)
    SSTR = ADJUSTL(SSTR)
    STRL = LEN_TRIM(SSTR)
    IF( STRL > 0 )THEN
      N=N+1
      DO I=1,STRL
        IF( ICHAR(SSTR(I:I)) == 32) EXIT
      ENDDO
      TERM(N)=SSTR(1:I-1)
      SSTR = SSTR(I:STRL)
    ELSE
      EXIT
    ENDIF
  ENDDO

  END SUBROUTINE

  SUBROUTINE GETWAVEDAY
  ! ** READ WAVETIME.INP 
  ! ** OUTPUT:
  ! ** WAVEDAY(1:NWVTIM)
  CHARACTER(200) :: STR
  INTEGER :: NP,IOS

  WRITE(*,'(A)')'READING WAVETIME.INP'
  OPEN(1,FILE='wavetime.inp',ACTION='READ')
  STR=READSTR(1)  
  NP=0
  DO WHILE(1)
    READ(1,'(A)',END=100,IOSTAT=IOS) STR
    IF( IOS > 0 ) STOP 'WAVETIME.INP: READING ERROR'  
    NP=NP+1
  ENDDO
  100 REWIND(1) 
  NWVTIM = NP
  ALLOCATE(WAVEDAY(NWVTIM+1))
  STR=READSTR(1)  
  DO NP=1,NWVTIM
    READ(1,*) WAVEDAY(NP)
  ENDDO
  CLOSE(1)
  END SUBROUTINE
    
END MODULE
