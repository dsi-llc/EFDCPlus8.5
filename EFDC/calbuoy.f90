SUBROUTINE CALBUOY(UPDATE)  

  !**********************************************************************!
  ! *** CALBUOY CALCULATES THE BUOYANCY USING MELLOR'S APPROXIMATION  
  ! *** TO THE UNESCO EQUATION OF STATE   
  ! *** MELLOR, G.L., J. ATM AND OCEAN TECH, VOL 8, P 609  
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! **  2015-06     PAUL M. CRAIG      IMPLEMENTED SIGMA-Z (SGZ) IN EE7.3 
  ! **  2014-09     PAUL M. CRAIG      ADDED THE LWET BYPASS APPROACH
  ! **  2014-08     D H CHUNG          SET EXPLICIT PRECISIONS OF INTEGER & REAL
  ! **  2011-03     Paul M. Craig      Converted to F90, added OMP
  !  

  USE GLOBAL
  IMPLICIT NONE

  LOGICAL,INTENT(IN) :: UPDATE
  INTEGER            :: NS,K,L,ND,LP,NN
  REAL(RKD), SAVE    :: RHOO,ONED,RHO1
  REAL(RKD)          :: SSTMP,TTMP,RHTMP,TEM0
  
  !************************************************************************************
  IF( IBSC == 1 )THEN
    ! *** DENSITY AS A LINEAR FUNCTION OF SALINITY ONLY.  FOR DIAGNOSTIC PURPOSES ONLY  
    DO K=1,KC  
      DO L=2,LA  
        B(L,K)=0.00075_8*SAL(L,K)  
      ENDDO  
    ENDDO  
    RETURN  
  ENDIF
  
  !************************************************************************************
  ! *** DENSITY RHOO AT P=0, S=0, AND T=TEMO.  ONLY COMPUTE ONCE, SINCE TEMO IS A CONSTANT 
  IF( N <= 5 )THEN  
    ONED = 1._8
    TEM0 = ABS(TEMO)
    RHOO = 999.842594 + 6.793952D-2*TEM0 - 9.095290D-3*TEM0*TEM0 + 1.001685D-4*TEM0*TEM0*TEM0 - 1.120083D-6*TEM0*TEM0*TEM0*TEM0 + 6.536332D-9*TEM0*TEM0*TEM0*TEM0*TEM0  
  ENDIF
  
  !$OMP PARALLEL DEFAULT(SHARED)
  
  !************************************************************************************
  ! *** SAVE THE CURRENT BUOYANCY BEFORE UPDATING
  IF( UPDATE )THEN  
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          B1(L,K)=B(L,K)  
        ENDDO
      ENDDO  
    ENDDO
    !$OMP END DO
  ENDIF
  
  !************************************************************************************
  ! *** DENSITY CALCULATIONS
  
  !$OMP DO PRIVATE(ND,K,L,LP,NN,NS,SSTMP,TTMP,RHTMP,TEM0,RHO1)
  DO ND=1,NDM 
      
    ! *** CASE: NO SALINITY AND NO TEMPERATURE
    IF( ISTRAN(1) == 0 .AND. ISTRAN(2) == 0 )THEN  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  

          ! *** REPLACE DENSITY B(L,K) WITH BUOYANCY B(L,K)  
          B(L,K) = 0.  
        ENDDO  
      ENDDO  

    ! *** CASE: SALINITY AND NO TEMPERATURE
    ELSEIF( ISTRAN(1) >= 1 .AND. ISTRAN(2) == 0 )THEN  
      TEM0=ABS(TEMO)
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          SSTMP=MAX(SAL(L,K),0.)  
          RHO1 = RHOO+SSTMP*(0.824493 - 4.0899D-3*TEM0 + 7.6438D-5*TEM0*TEM0 - 8.2467D-7*TEM0*TEM0*TEM0 + 5.3875D-9*TEM0*TEM0*TEM0*TEM0)  &
                + SQRT(SSTMP)*SSTMP*(-5.72466D-3 + 1.0227D-4*TEM0 - 1.6546D-6*TEM0*TEM0) + 4.8314D-4*SSTMP*SSTMP  
          RHOW(L,K) = RHO1
          
          ! *** REPLACE DENSITY B(L,K) WITH BUOYANCY B(L,K)  
          B(L,K) = (RHO1/RHOO)-1._8  
        ENDDO  
      ENDDO  
    
    ! *** CASE: NO SALINITY BUT WITH TEMPERATURE
    ELSEIF( ISTRAN(1) == 0 .AND. ISTRAN(2) >= 1 )THEN  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          TTMP=TEM(L,K)  
          RHO1 = 999.842594+6.793952D-2*TTMP-9.095290D-3*TTMP*TTMP + 1.001685D-4*TTMP*TTMP*TTMP - 1.120083D-6*TTMP*TTMP*TTMP*TTMP + 6.536332D-9*TTMP*TTMP*TTMP*TTMP*TTMP  
          RHOW(L,K) = RHO1
          
          ! *** REPLACE DENSITY B(L,K) WITH BUOYANCY B(L,K)  
          B(L,K) = (RHO1/RHOO)-1._8  
        ENDDO  
      ENDDO  
    
    ! *** CASE: BOTH SALINITY AND TEMPERATURE
    ELSEIF( ISTRAN(1) >= 1 .AND. ISTRAN(2) >= 1 )THEN  
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          SSTMP=MAX(SAL(L,K),0.)  
          TTMP=TEM(L,K)  
          RHTMP = 999.842594 + 6.793952D-2*TTMP - 9.095290D-3*TTMP*TTMP + 1.001685D-4*TTMP*TTMP*TTMP - 1.120083D-6*TTMP*TTMP*TTMP*TTMP + 6.536332D-9*TTMP*TTMP*TTMP*TTMP*TTMP 
          RHO1= RHTMP+SSTMP*(0.824493 - 4.0899D-3*TTMP + 7.6438D-5*TTMP*TTMP - 8.2467D-7*TTMP*TTMP*TTMP + 5.3875D-9*TTMP*TTMP*TTMP*TTMP) &
                      +SQRT(SSTMP)*SSTMP*(-5.72466D-3+1.0227D-4*TTMP - 1.6546D-6*TTMP*TTMP) + 4.8314D-4*SSTMP*SSTMP  
          RHOW(L,K) = RHO1

          ! *** REPLACE DENSITY B(L,K) WITH BUOYANCY B(L,K)  
          B(L,K) = (RHO1/RHOO)-1._8  
        ENDDO  
      ENDDO  
    ENDIF  

    !------------------------------------------------------------------------------------
    ! *** APPLY LOW SEDIMENT CONCENTRATION CORRECTION TO BUOYANCY  
    IF( ISTRAN(6) >= 1 .OR. ISTRAN(7) >= 1 )THEN  
      
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          TVAR1S(L,K)=0.  
          TVAR1W(L,K)=0.  
        ENDDO  
      ENDDO  

      IF( ISTRAN(6) >= 1 )THEN  
        DO NS=1,NSED  
          DO K=1,KC  
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              TVAR1S(L,K)=TVAR1S(L,K)+SDEN(NS)*SED(L,K,NS)  
              TVAR1W(L,K)=TVAR1W(L,K)+(SSG(NS)-1.)*SDEN(NS)*SED(L,K,NS)  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  
      IF( ISTRAN(7) >= 1 )THEN  
        DO NN=1,NSND  
          NS=NN+NSED  
          DO K=1,KC  
            DO LP=1,LLWET(K,ND)
              L=LKWET(LP,K,ND)  
              TVAR1S(L,K)=TVAR1S(L,K)+SDEN(NS)*SND(L,K,NN)  
              TVAR1W(L,K)=TVAR1W(L,K)+(SSG(NS)-1.)*SDEN(NS)*SND(L,K,NN)  
            ENDDO  
          ENDDO  
        ENDDO  
      ENDIF  

      IF( ISTRAN(1) == 0 .AND. ISTRAN(2) == 0 )THEN
        ! *** RESET DENSITY 
        DO K=1,KC  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            RHOW(L,K) = RHOO         
          ENDDO  
        ENDDO  
      ENDIF
      
      DO K=1,KC  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          B(L,K) = B(L,K)*(1.-TVAR1S(L,K)) + TVAR1W(L,K) 

          ! **  CORRECTION FOR SEDIMENT
          RHOW(L,K) = RHOW(L,K)*( 1. - TVAR1S(L,K) + TVAR1W(L,K) )          
        ENDDO  
      ENDDO  
      
    ENDIF
  ENDDO    ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP END PARALLEL
  
  !************************************************************************************
  RETURN  

END  

