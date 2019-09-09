SUBROUTINE CALDIFF (CON1,IT)
                                                                                                                         
  ! **  SUBROUTINE CALDIFF CALCULATES THE HORIZONTAL DIFFUSIVE                                                            
  ! **  TRANSPORT OF DISSOLVED OR SUSPENDED CONSITITUENT M LEADING TO                                                     
  ! **  A REVISEDED VALUE AT TIME LEVEL (N+1). THE VALUE OF ISTL                                                          
  ! **  INDICATES THE NUMBER OF TIME LEVELS IN THE STEP                                                                   

  ! CHANGE RECORD      
  ! DATE MODIFIED     BY               DESCRIPTION
  !------------------------------------------------------------------------------------!
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH

  USE GLOBAL
  IMPLICIT NONE   
  INTEGER,INTENT(IN) :: IT                                                                                                                                                                                                       
  REAL,   INTENT(IN) :: CON1(LCM,KCM)                                                                                                          
  INTEGER            :: ND,K,LP,L,LS,LW

  ! **  HORIZONTAL DIFFUSIVE FLUX CALCULATION                                                                             
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLHDMF(K,ND)
        L = LKHDMF(LP,K,ND)
        LW=LWC(L)
        ! ***                G/S                        M       M            (      M^2/S     ) (        G/M^3       )  1/M
        FUHUD(L,K,IT) = FUHUD(L,K,IT) + 0.5*SUB3D(L,K)*DYU(L)*HU(L)*DZC(L,K)*(AH(L,K)+AH(LW,K))*(CON1(LW,K)-CON1(L,K))*DXIU(L)
        LS=LSC(L)
        FVHUD(L,K,IT) = FVHUD(L,K,IT) + 0.5*SVB3D(L,K)*DXV(L)*HV(L)*DZC(L,K)*(AH(L,K)+AH(LS,K))*(CON1(LS,K)-CON1(L,K))*DYIV(L)
      ENDDO
    ENDDO
  ENDDO

  RETURN
END SUBROUTINE

