FUNCTION ZBRENT(L, ISMERR, SMCH4S, SMK1CH4, SMO2JC,                                &
                     SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                     SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                     SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                     RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                     CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, RSMSS)

  ! USING BRENT'S METHOD (ZBRENT), FIND THE ROOT OF A FUNC SEDFLUX KNOWN TO LIE
  !   BETWEEN RMIN & RMAX WITHIN AN ACCURACY OF TOL (P. 253 IN NUMERICAL RECIPE).
  
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2013-03           Paul M. Craig    Restructed to F90, updated variables for OMP
  
  IMPLICIT NONE
  
  INTEGER, INTENT(IN)  :: L
  INTEGER, INTENT(OUT) :: ISMERR
  REAL, INTENT(IN)     :: SMCH4S, SMK1CH4, SMO2JC
  REAL, INTENT(IN)     :: SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM  
  REAL, INTENT(IN)     :: SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM
  REAL, INTENT(IN)     :: SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM
  REAL, INTENT(INOUT)  :: RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3
  REAL, INTENT(INOUT)  :: CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, RSMSS
  
  INTEGER :: IZMAX,II
  REAL    :: EPS,TOL,RMIN,RMAX,ZBRENT
  REAL    :: A,B,C,D,E,FA,FB,FC,TOL1,XM
  REAL    :: P,Q,R,S

  EXTERNAL SEDFLUXNEW

  PARAMETER (IZMAX=100,EPS=3.0E-8,TOL=1.0E-5,RMIN=1.0E-4,RMAX=100.0)

  ISMERR = 0
  A = RMIN
  B = RMAX

  !FA = SEDFLUX(A,      RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3, RSMSS)
  CALL SEDFLUXNEW(L, A, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                              &
                  SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                  SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                  SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                  RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                  CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, FA)

  !FB = SEDFLUX(B,      RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3, RSMSS)
  CALL SEDFLUXNEW(L, B, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                              &
                  SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                  SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                  SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                  RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                  CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, FB)

  ZBRENT = 0
  IF( FA*FB > 0.0 )THEN
    ISMERR = 1
    RETURN
  ENDIF
  FC = FB
  DO II=1,IZMAX
    IF( FB*FC > 0.0 )THEN
      C = A
      FC = FA
      D = B-A
      E = D
    ENDIF
    IF( ABS(FC) < ABS(FB) )THEN
      A = B
      B = C
      C = A
      FA = FB
      FB = FC
      FC = FA
    ENDIF
    TOL1 = 2.0*EPS*ABS(B) + 0.5*TOL
    XM = 0.5 * (C-B)
    IF( ABS(XM) <= TOL1 .OR. FB == 0.0 )THEN
      ZBRENT = B
      RETURN
    ENDIF
    IF( ABS(E) >= TOL1 .AND. ABS(FA) > ABS(FB) )THEN
      S = FB / FA
      IF( A == C )THEN
        P = 2.0 * XM * S
        Q = 1.0 - S
      ELSE
        Q = FA / FC
        R = FB / FC
        P = S * (2.0*XM*Q*(Q-R) - (B-A)*(R-1.0))
        Q = (Q-1.0) * (R-1.0) * (S-1.0)
      ENDIF
      IF( P > 0.0) Q = -Q
      P = ABS(P)
      IF( 2.0*P  <  MIN(3.0*XM*Q-ABS(TOL1*Q), ABS(E*Q)) )THEN
        E = D
        D = P / Q
      ELSE
        D = XM
        E = D
      ENDIF
    ELSE
      D = XM
      E = D
    ENDIF
    A = B
    FA = FB
    IF( ABS(D) > TOL1 )THEN
      B = B + D
    ELSE
      B = B + SIGN(TOL1,XM)
    ENDIF
    !FB = SEDFLUX(B, RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3, RSMSS)
    CALL SEDFLUXNEW(L, B, SMCH4S, SMK1CH4, SMO2JC, RSMSS,                              &
                    SK1NH4SM, A1NH4SM, A2NH4SM,           A22NH4SM, B1NH4SM, B2NH4SM,  &
                    SK1NO3SM, A1NO3SM, A2NO3SM, RK2NO3SM, A22NO3SM, B1NO3SM, B2NO3SM,  &
                    SK1H2SSM, A1H2SSM, A2H2SSM,           A22H2SSM, B1H2SSM, B2H2SSM,  &
                    RSM1H2S, RSM1NH4, RSM1NO3, RSM2H2S, RSM2NH4, RSM2NO3,              &
                    CSODSM, RNSODSM, RJNITSM, RJDENSM, AQJH2SSM, AQJCH4SM, GJCH4SM, FB)

  ENDDO
  ISMERR = 2
  ZBRENT = B

  RETURN

END

