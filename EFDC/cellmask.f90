SUBROUTINE CELLMASK

  ! **  SUBROUTINE CELLMASK CONVERTS LAND CELLS TO WATER CELLS BY
  ! **  MASKING VARIABLES.  DEPTHS IN THE MASKED CELLS SHOULD BE INPUT AT
  ! **  THE END OF THE DXDY.INP FILE.

  ! CHANGE RECORD

  USE GLOBAL
  USE INFOMOD,ONLY:READSTR

  IMPLICIT NONE

  INTEGER :: I,J,K,L,LN,M,MMASK,MTYPE
  CHARACTER(200) :: STR
  
  WRITE(*,'(A)')'READING MASK.INP'
  OPEN(1,FILE='mask.inp',STATUS='UNKNOWN')
  STR = READSTR(1)
  READ(1,*,ERR=1000) MMASK
  
  DO M=1,MMASK
    READ(1,*,ERR=1000) I,J,MTYPE
    L=LIJ(I,J)

    IF( MTYPE == 1 )THEN
      SUB(L)=0.
      SUBO(L)=0.
      UHDYE(L)=0.
      UHDY1E(L)=0.
      UHDY2E(L)=0.
      UMASK(L)=1
      DO K=1,KC
        U(L,K)=0.
        U1(L,K)=0.
        U2(L,K)=0.
        UHDY(L,K)=0.
        UHDY1(L,K)=0.
        UHDY2(L,K)=0.
        UHDYF(L,K)=0.
        UHDYF1(L,K)=0.
        UHDYF2(L,K)=0.
      ENDDO
    ENDIF

    IF( MTYPE == 2 )THEN
      SVB(L)=0.
      SVBO(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      VHDX2E(L)=0.
      VMASK(L) = 1
      DO K=1,KC
        V(L,K)=0.
        V1(L,K)=0.
        V2(L,K)=0.
        VHDX(L,K)=0.
        VHDX1(L,K)=0.
        VHDX2(L,K)=0.
        VHDXF(L,K)=0.
        VHDXF1(L,K)=0.
        VHDXF2(L,K)=0.
      ENDDO
    ENDIF

    IF( MTYPE == 3 )THEN   ! *** PMC
      ! *** U Face
      SUB(L)=0.
      SUBO(L)=0.
      UHDYE(L)=0.
      UHDY1E(L)=0.
      UHDY2E(L)=0.
      DO K=1,KC
        U(L,K)=0.
        U1(L,K)=0.
        U2(L,K)=0.
        UHDY(L,K)=0.
        UHDY1(L,K)=0.
        UHDY2(L,K)=0.
        UHDYF(L,K)=0.
        UHDYF1(L,K)=0.
        UHDYF2(L,K)=0.
      ENDDO
      ! *** V Face
      SVB(L)=0.
      SVBO(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      VHDX2E(L)=0.
      DO K=1,KC
        V(L,K)=0.
        V1(L,K)=0.
        V2(L,K)=0.
        VHDX(L,K)=0.
        VHDX1(L,K)=0.
        VHDX2(L,K)=0.
        VHDXF(L,K)=0.
        VHDXF1(L,K)=0.
        VHDXF2(L,K)=0.
      ENDDO
    ENDIF

    IF( MTYPE == 4 )THEN   ! *** PMC - Change to MTYPE 4 for isolated !waters
      LN=LNC(L)
      SUB(L)=0.
      SUBO(L)=0.
      UHDYE(L)=0.
      UHDY1E(L)=0.
      UHDY2E(L)=0.
      SUB(LEC(L))=0.
      SUBO(LEC(L))=0.
      UHDYE(LEC(L))=0.
      UHDY1E(LEC(L))=0.
      UHDY2E(LEC(L))=0.
      SVB(L)=0.
      SVBO(L)=0.
      VHDXE(L)=0.
      VHDX1E(L)=0.
      VHDX2E(L)=0.
      SVB(LN)=0.
      SVBO(LN)=0.
      VHDXE(LN)=0.
      VHDX1E(LN)=0.
      VHDX2E(LN)=0.
      P(L)=0.
      P1(L)=0.
      DO K=1,KC
        B(L,K)=0.
        B1(L,K)=0.
        SAL(L,K)=0.
        SAL1(L,K)=0.
        TEM(L,K)=TEMO
        TEM1(L,K)=TEMO
        DYE(L,K)=0.
        DYE1(L,K)=0.
        SED(L,K,1)=0.
        SED1(L,K,1)=0.
        QQ(L,K)=0.
        QQ1(L,K)=0.
        QQL(L,K)=0.
        QQL1(L,K)=0.
        U(L,K)=0.
        U1(L,K)=0.
        U2(L,K)=0.
        UHDY(L,K)=0.
        UHDY1(L,K)=0.
        UHDY2(L,K)=0.
        UHDYF(L,K)=0.
        UHDYF1(L,K)=0.
        UHDYF2(L,K)=0.
        U(LEC(L),K)=0.
        U1(LEC(L),K)=0.
        U2(LEC(L),K)=0.
        UHDY(LEC(L),K)=0.
        UHDY1(LEC(L),K)=0.
        UHDY2(LEC(L),K)=0.
        V(L,K)=0.
        V1(L,K)=0.
        V2(L,K)=0.
        VHDX(L,K)=0.
        VHDX1(L,K)=0.
        VHDX2(L,K)=0.
        VHDXF(L,K)=0.
        VHDXF1(L,K)=0.
        VHDXF2(L,K)=0.
        V(LN,K)=0.
        V1(LN,K)=0.
        V2(LN,K)=0.
        VHDX(LN,K)=0.
        VHDX1(LN,K)=0.
        VHDX2(LN,K)=0.
      ENDDO
    ENDIF
  ENDDO
  CLOSE(1)

  ! **  WRITE READ ERRORS ON CELLMASK
  GOTO 1002
  
  1000 WRITE(6,1001)
  1001 FORMAT('  READ ERROR ON FILE MASK.INP ')
  STOP

  1002 CONTINUE
  
  RETURN
END

