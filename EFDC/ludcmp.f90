SUBROUTINE LUDCMP(A,N,NP,INDX,D)

  ! CHANGE RECORD

  IMPLICIT NONE

  INTEGER :: N,NP,INDX(N),I,J,K,IMAX
  REAL :: A(NP,NP),TINY,D,AAMAX,DUM,SUM
  
  PARAMETER (TINY=1.0E-20)
  REAL,ALLOCATABLE,DIMENSION(:) :: VV

  ALLOCATE(VV(N))
  !
  D=1.
  DO 12 I=1,N
    AAMAX=0.
    DO 11 J=1,N
      IF( ABS(A(I,J)) > AAMAX) AAMAX=ABS(A(I,J))
     11   CONTINUE
    IF( AAMAX == 0. ) PAUSE 'SINGULAR MATRIX.'
    VV(I)=1./AAMAX
     12 CONTINUE
  DO 19 J=1,N
    IF( J > 1 )THEN
      DO 14 I=1,J-1
        SUM=A(I,J)
        IF( I > 1 )THEN
          DO 13 K=1,I-1
            SUM=SUM-A(I,K)*A(K,J)
     13         CONTINUE
          A(I,J)=SUM
        ENDIF
     14     CONTINUE
    ENDIF
    AAMAX=0.
    DO 16 I=J,N
      SUM=A(I,J)
      IF( J > 1 )THEN
        DO 15 K=1,J-1
          SUM=SUM-A(I,K)*A(K,J)
     15       CONTINUE
        A(I,J)=SUM
      ENDIF
      DUM=VV(I)*ABS(SUM)
      IF( DUM >= AAMAX )THEN
        IMAX=I
        AAMAX=DUM
      ENDIF
     16   CONTINUE
    IF( J /= IMAX )THEN
      DO 17 K=1,N
        DUM=A(IMAX,K)
        A(IMAX,K)=A(J,K)
        A(J,K)=DUM
     17     CONTINUE
      D=-D
      VV(IMAX)=VV(J)
    ENDIF
    INDX(J)=IMAX
    IF( J /= N )THEN
      IF( A(J,J) == 0. )A(J,J)=TINY
      DUM=1./A(J,J)
      DO 18 I=J+1,N
        A(I,J)=A(I,J)*DUM
     18     CONTINUE
    ENDIF
     19 CONTINUE
  IF( A(N,N) == 0. )A(N,N)=TINY
  DEALLOCATE(VV)
  RETURN
END

