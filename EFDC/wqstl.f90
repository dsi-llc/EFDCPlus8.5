SUBROUTINE RWQSTL(IWQTSTL)

  ! READ IN SPATIALLY AND/OR TEMPORALLY VARYING PARAMETERS FOR SETTLING
  ! VELOCITIES OF ALGAE, RPOM, LPOM & PARTICULATE METAL (UNIT INWQSTL).
  ! ALSO SPATIALLY/TEMPORALLY VARYING REAERATION ADJUSTMENT FACTOR.
  !
  ! CHANGE RECORD
  !  2011-08   PMC - THIS NEVER WORKED FOR TIME VARYING SETTLING RATES.
  !                  NEEDS TO BE FIXED!


  ! ***    WQWSC = Settling velocity for cyanobacteria (m/day)
  ! ***    WQWSD = Settling velocity for algae diatoms (m/day)
  ! ***    WQWSG = Settling velocity for algae greens  (m/day)
  ! ***   WQWSRP = Settling velocity for refractory POM (m/day)
  ! ***   WQWSLP = Settling velocity for labile POM (m/day)
  ! ***    WQWSS = Settling velocity for particles sorbed to TAM (m/day)
  ! ***    WQWSM = Settling velocity for macroalgae (m/day = 0.0)
  ! ***    WQWSM = Reaeration adjustment factor (NOT SAVED)

  USE GLOBAL
  IMPLICIT NONE
  INTEGER :: IWQTSTL
  INTEGER :: I,M,MM
  CHARACTER TITLE(3)*79, STLCONT*3
  !
  WRITE(*,'(A)')' WQ: '//STLFN
  OPEN(1,FILE=STLFN,STATUS='UNKNOWN')
  OPEN(2,FILE=OUTDIR//'WQ3D.OUT',STATUS='UNKNOWN',POSITION='APPEND')
  IF( IWQTSTL == 0 )THEN
    READ(1,50) (TITLE(M),M=1,3)
    WRITE(2,999)
    WRITE(2,50) (TITLE(M),M=1,3)
  ENDIF
  WRITE(2,60)'* SETTLING VELOCITY AT  ', IWQTSTL, ' TH DAY FROM MODEL START'
  READ(1,999)
  READ(1,50) TITLE(1)
  WRITE(2,50) TITLE(1)
  DO I=1,IWQZ
    READ(1,*)   MM,WQWSC(I),WQWSD(I),WQWSG(I),WQWSRP(I),WQWSLP(I),WQWSS(I), WQWSM
    WRITE(2,51) MM,WQWSC(I),WQWSD(I),WQWSG(I),WQWSRP(I),WQWSLP(I),WQWSS(I), WQWSM
  ENDDO
  READ(1,52) IWQTSTL, STLCONT
  WRITE(2,52) IWQTSTL, STLCONT
  IF( STLCONT == 'END' )THEN
    CLOSE(1)
    IWQSTL = 0
  ENDIF
  CLOSE(2)
    999 FORMAT(1X)
     50 FORMAT(A79)
     51 FORMAT(I3, 10F8.3)
     52 FORMAT(I7, 1X, A3)
     60 FORMAT(/, A24, I5, A24)
  RETURN
END

