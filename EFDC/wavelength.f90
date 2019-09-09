MODULE WAVELENGTH
  !Author: Dang Huu Chung

  USE GLOBAL,ONLY:RKD
  IMPLICIT NONE

  REAL(RKD),PRIVATE,PARAMETER :: G=9.81
  REAL(RKD),PRIVATE,PARAMETER :: PI=3.14159265358979

  CONTAINS

  FUNCTION DISRELATION(RLS,TP,HD,U,PHI) RESULT(FWL)
    !Dispersion Relation:          FWL=0
    !RLS: Wave length              [m]
    !TP : Wave period              [s]
    !HD : Water depth              [m]
    !U  : Depth-average velocity   [m/s]
    !PHI: Angle of (wave,current)  [Rad]
    REAL(RKD) :: FWL,RLS,TP,HD,U,PHI
    FWL=(RLS/TP-U*COS(PHI))-SQRT(G*RLS/2._8/PI*TANH(2._8*PI*HD/RLS))
  END FUNCTION

  RECURSIVE SUBROUTINE BISEC(FUN,A0,B0,TOL,TP,HD,U,PHI,X)
    REAL(RKD),INTENT(IN)       :: A0,B0,TOL,TP,HD,U,PHI
    REAL(RKD),EXTERNAL         :: FUN
    REAL(RKD),INTENT(OUT)      :: X
    INTEGER                    :: IST
    REAL(RKD)                  :: A,B,FA,FB,FX
    A=A0
    B=B0
    IST=1
    FA=FUN(A,TP,HD,U,PHI)
    FB=FUN(B,TP,HD,U,PHI)
    IF( FA*FB < 0 )THEN
      X=0.5*(A+B)
      FX=FUN(X,TP,HD,U,PHI)
      IF( FA*FX < 0 )THEN
        B=X
      ELSE
        A=X
      ENDIF
      IF( ABS(A-B) <= TOL )THEN
        RETURN
      ELSE
        CALL BISEC(FUN,A,B,TOL,TP,HD,U,PHI,X)
      ENDIF
    ELSEIF (FA == 0 )THEN
      X=A
    ELSEIF (FB == 0 )THEN
      X=B
    ELSE
      IST=-1
      STOP 'DISPERSION RELATION: FA.FB>0'       
    ENDIF
  END SUBROUTINE BISEC
  
  FUNCTION RTBIS(FUNC,X1,X2,XACC,TP,HD,U,PHI)
    REAL(RKD) :: RTBIS,X1,X2,XACC,TP,HD,U,PHI
    REAL(RKD),EXTERNAL :: FUNC
    INTEGER   :: J,JMAX
    REAL(RKD) :: DX,F,FMID,XMID
    PARAMETER (JMAX=40)
  
    FMID=FUNC(X2,TP,HD,U,PHI)
    F=FUNC(X1,TP,HD,U,PHI)
    IF( F*FMID >= 0.) PAUSE 'ROOT MUST BE BRACKETED IN RTBIS'
    IF( F<0. )THEN
      RTBIS=X1
      DX=X2-X1
    ELSE
      RTBIS=X2
      DX=X1-X2
    ENDIF
    DO J=1,JMAX
      DX=DX*.5
      XMID=RTBIS+DX
      FMID=FUNC(XMID,TP,HD,U,PHI)
      IF( FMID <= 0.)RTBIS=XMID
      IF( ABS(DX)<XACC .OR. FMID == 0. ) RETURN
    ENDDO
    PAUSE 'TOO MANY BISECTIONS IN RTBIS'
  END FUNCTION

END MODULE


   
