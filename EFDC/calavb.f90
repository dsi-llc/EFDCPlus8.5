SUBROUTINE CALAVB ()  

  ! ****************************************************************************
  ! *** SUBROUTINE CALAV CALCULATES VERTICAL VISCOSITY AND DIFFUSIVITY  
  ! *** USING GLAPERIN ET AL'S MODIFICATION OF THE MELLOR-YAMADA MODEL  
  ! *** (NOTE AV, AB, AND AQ ARE ACTUALLY DIVIDED BY H)  
  ! *** IF ISGA=1 VALUES ARE GEOMETRI!AVERAGES WITH THE PREVIOUS VALUES  
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !    2018-10       PAUL M. CRAIG     IMPLEMENTED SPATIALLY VARYING AVO AND ABO
  !    2015-06       PAUL M. CRAIG     IMPLEMENTED SIGMA-Z (SGZ) IN EE8.0
  !    2014-09       PAUL M. CRAIG     ADDED THE LWET BYPASS APPROACH
  !    2011-03-11    PAUL M. CRAIG     MERGED LATEST CODES AND RESTRUCTURED, ADDED OMP
  !                  JOHN HAMRICK      ADDED DRYCELL BYPASS AND CONSISTENT INITIALIZATION OF DRY VALUES  

  USE GLOBAL  
  IMPLICIT NONE

  INTEGER :: L,K,ND,LF,LL,LP
  REAL    :: QQIMAX,RIQMIN,RIQ,SFAV,SFAB,ABTMP,AVTMP
  REAL    :: SFAV0,SFAV1,SFAV2,SFAV3,SFAB0,SFAB1,RIQMAXX

  ! ***  GALPERIN (DEFAULT)
  SFAV0= 0.392010
  SFAV1= 7.760050
  SFAV2=34.676440
  SFAV3= 6.127200
  SFAB0= 0.493928
  SFAB1=34.676440
  RIQMIN=-0.999/SFAB1

  ! *** KANTHA AND CLAYSON (1994)
  IF( ISTOPT(0) == 2 )THEN
    SFAV0= 0.392010
    SFAV1= 8.679790
    SFAV2=30.192000
    SFAV3= 6.127200
    SFAB0= 0.493928
    SFAB1=30.192000
    RIQMIN=-0.999/SFAB1
  ENDIF

  ! ***  KANTHA (2003)
  IF( ISTOPT(0) == 3 )THEN
    SFAV0= 0.392010
    SFAV1=14.509100
    SFAV2=24.388300
    SFAV3= 3.236400
    SFAB0= 0.490025
    SFAB1=24.388300
    RIQMIN=-0.999/SFAB1
  ENDIF

  RIQMAXX = 1.E32
  IF( ISLLIM >= 1 ) RIQMAXX = RIQMAX
  
  QQIMAX = 1./QQMIN  

  ! ****************************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)
  
  IF( BSC <= 1.E-6 )THEN
    ! *** VERTICAL DIFFUSVITY TIME FILTERING - NONE
    !$OMP DO PRIVATE(ND,K,LP,L,ABTMP,AVTMP) 
    DO ND=1,NDM
      IF( ISFAVB == 0 .OR. N == 1 )THEN
        ! *** VERTICAL DIFFUSVITY TIME FILTERING - NONE
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            AB(L,K) = AVCON*SFAB0*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AV(L,K) = AVCON*SFAV0*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = AV(L,K)*HPI(L) 
            AB(L,K) = AB(L,K)*HPI(L)
          ENDDO  
        ENDDO  

      ELSEIF( ISFAVB == 1 .AND. N > 1 )THEN
        ! *** VERTICAL DIFFUSVITY TIME FILTERING - AVERAGE
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ABTMP = AVCON*SFAB0*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = AVCON*SFAV0*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = 0.5*(AV(L,K) + AVTMP*HPI(L))  
            AB(L,K) = 0.5*(AB(L,K)+ABTMP*HPI(L))  
          ENDDO  
        ENDDO  

      ELSEIF( ISFAVB == 2 .AND. N > 1 )THEN
        ! *** VERTICAL DIFFUSVITY TIME FILTERING - SQRT
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            ABTMP = AVCON*SFAB0*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = AVCON*SFAV0*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = SQRT(AV(L,K)*AVTMP*HPI(L))  
            AB(L,K) = SQRT(AB(L,K)*ABTMP*HPI(L))  
          ENDDO  
        ENDDO  
      ENDIF
    ENDDO
    !$OMP END DO
  
  ELSE
    
    !----------------------------------------------------------------------C
    IF( ISFAVB == 0 .OR. N == 1 )THEN
      ! *** VERTICAL DIFFUSVITY TIME FILTERING - NONE
      !$OMP DO PRIVATE(ND,K,LP,L,RIQ,SFAV,SFAB) 
      DO ND=1,NDM  
    
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            RIQ = -GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(L,K)*(B(L,K+1)-B(L,K))/QQ(L,K)  
            RIQ = MAX(RIQ,RIQMIN)  
            RIQ = MIN(RIQ,RIQMAXX)
            SFAV = SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB = SFAB0/(1.+SFAB1*RIQ)
            AB(L,K) = AVCON*SFAB*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)
            AV(L,K) = AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)
            AV(L,K) = AV(L,K)*HPI(L)  
            AB(L,K) = AB(L,K)*HPI(L)
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO

    !----------------------------------------------------------------------C
    ELSEIF( ISFAVB == 1 .AND. N > 1 )THEN
      ! *** VERTICAL DIFFUSVITY TIME FILTERING - AVERAGE
      !$OMP DO PRIVATE(ND,K,LP,L,RIQ,SFAV,SFAB,AVTMP,ABTMP) 
      DO ND=1,NDM  
    
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L=LKWET(LP,K,ND)  
            RIQ = -GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(L,K)*(B(L,K+1)-B(L,K))/QQ(L,K)  
            RIQ = MAX(RIQ,RIQMIN)  
            RIQ = MIN(RIQ,RIQMAXX)  
            SFAV = SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB = SFAB0/(1.+SFAB1*RIQ)
            ABTMP = AVCON*SFAB*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = 0.5*(AV(L,K) + AVTMP*HPI(L))  
            AB(L,K) = 0.5*(AB(L,K)+ABTMP*HPI(L))  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO

    !----------------------------------------------------------------------C
    ELSEIF( ISFAVB == 2 .AND. N > 1 )THEN
      ! *** VERTICAL DIFFUSVITY TIME FILTERING - SQRT
      !$OMP DO PRIVATE(ND,K,LP,L,RIQ,SFAV,SFAB,AVTMP,ABTMP) 
      DO ND=1,NDM  
    
        DO K=1,KS  
          DO LP=1,LLWET(K,ND)
            L = LKWET(LP,K,ND)  
            RIQ = -GP*HP(L)*DML(L,K)*DML(L,K)*DZIG(L,K)*(B(L,K+1)-B(L,K))/QQ(L,K)  
            RIQ = MAX(RIQ,RIQMIN)  
            RIQ = MIN(RIQ,RIQMAXX)  
            SFAV = SFAV0*(1.+SFAV1*RIQ)/((1.+SFAV2*RIQ)*(1.+SFAV3*RIQ))
            SFAB = SFAB0/(1.+SFAB1*RIQ)
            ABTMP = AVCON*SFAB*DML(L,K)*HP(L)*QQSQR(L,K) + AVBXY(L)  
            AVTMP = AVCON*SFAV*DML(L,K)*HP(L)*QQSQR(L,K) + AVOXY(L)  
            AV(L,K) = SQRT(AV(L,K)*AVTMP*HPI(L))  
            AB(L,K) = SQRT(AB(L,K)*ABTMP*HPI(L))  
          ENDDO  
        ENDDO  
      ENDDO
      !$OMP END DO
    ENDIF
    
  ENDIF
  
  ! *** APPLY OPEN BOUNDARYS 
  !$OMP SINGLE
  DO LL=1,NBCSOP
    L=LOBCS(LL)
    DO K=1,KS  
      AB(L,K)=0.0
    ENDDO  
  ENDDO 
  !$OMP END SINGLE

  ! ****************************************************************************
  IF( ISLOG >= 1 )THEN  
    !$OMP SINGLE
    AVMAX=AVO  
    ABMAX=ABO  
    AVMIN=10.  
    ABMIN=10.  
    DO K=1,KS  
      DO L=2,LA  
        IF( LKSZ(L,K) )CYCLE 
        AVMAX=MAX(AVMAX,AV(L,K))  
        ABMAX=MAX(ABMAX,AB(L,K))  
        AVMIN=MIN(AVMIN,AV(L,K))  
        ABMIN=MIN(ABMIN,AB(L,K))  
      ENDDO
    ENDDO
    !$OMP END SINGLE
  ENDIF

  ! ****************************************************************************
  ! *** CHECK FOR DEPTHS LESS THAN ZBR
  IF( ISDRY > 0 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,K) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
    
      DO LP=LF,LL  
        L=LWET(LP)  
        IF( HPK(L,KSZ(L)) < ZBR(L) )THEN
          DO K=KS-1,KSZ(L),-1
            IF( HP(L)*Z(L,K-1) > ZBR(L) )CYCLE
            AV(L,K) = AV(L,K+1)
            AB(L,K) = AB(L,K+1)
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    !$OMP END DO
  ENDIF
        
  ! ****************************************************************************
  ! *** NOW APPLY MAXIMUM, IF REQURIED
  IF( ISAVBMX >= 1 )THEN  
    !$OMP DO PRIVATE(ND,K,LP,L,AVTMP,ABTMP)
    DO ND=1,NDM  
    
      DO K=1,KS  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          AVTMP = AVMX*HPI(L)  
          ABTMP = ABMX*HPI(L)  
          AV(L,K) = MIN(AV(L,K),AVTMP)  
          AB(L,K) = MIN(AB(L,K),ABTMP)  
        ENDDO  
      ENDDO  
    ENDDO
    !$OMP END DO
  ENDIF  

  ! ****************************************************************************
  ! *** COMPUTE THE INVERSE OF THE AVERAGE AV AT THE U AND V INTERFACES
  IF( IGRIDV == 0 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KS  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          AVUI(L,K) = 2.0/( AV(L,K) + AV(LWC(L),K) )
          AVVI(L,K) = 2.0/( AV(L,K) + AV(LSC(L),K) )
        ENDDO  
      ENDDO 
    ENDDO
    !$OMP END DO
  ELSE
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=1,KS  
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          AVUI(L,K) = SUB(L)*(1. + SUB3D(L,K))/( AV(L,K)*HP(L) + SUB3D(L,K)*AV(LWC(L),K)*HP(LWC(L) ) )*HU(L)
          AVVI(L,K) = SVB(L)*(1. + SVB3D(L,K))/( AV(L,K)*HP(L) + SVB3D(L,K)*AV(LSC(L),K)*HP(LSC(L) ) )*HV(L)
        ENDDO  
      ENDDO 
    ENDDO
    !$OMP END DO
  ENDIF
  
  ! ****************************************************************************
  ! *** ISSQL:   0 SETS QQ AND QQL STABILITY FUNCTIONS PROPORTIONAL TO
  ! ***            MOMENTUM STABILITY FUNCTIONS
  IF( ISTOPT(0) <= 2 .OR. ISSQL == 0 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L)
    DO ND=1,NDM  
      DO K=2,KS  
        DO LP=1,LLWET(K-1,ND)
          L = LKWET(LP,K-1,ND) 
          AQ(L,K) = 0.205*(AV(L,K-1) + AV(L,K))  
        ENDDO  
      ENDDO  
      DO LP=1,LLWET(KS,ND) 
        L = LKWET(LP,KS,ND)
        AQ(L,KSZ(L)) = 0.205*AV(L,KSZ(L))  
        AQ(L,KC)     = 0.205*AV(L,KS)  
      ENDDO  
    ENDDO
    !$OMP END DO
    
  ENDIF

  ! ****************************************************************************
  ! *** ISSQL:     1  SETS QQ AND QQL STABILITY FUNCTIONS TO CONSTANTS
  ! *** ISTOPT(0)  3  KANTHA (2003)
  IF( ISTOPT(0) == 3 .OR. ISSQL == 1 )THEN
    IF( ISFAVB == 1 .AND. N > 1 )THEN
      ! *** AVERAGE FILTER
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
    
        DO K=2,KS  
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            AQ(L,K)=0.314*( DML(L,K-1)*QQSQR(L,K-1)+DML(L,K)*QQSQR(L,K) ) + AVOXY(L)*HPI(L)
          ENDDO
        ENDDO

        DO LP=LF,LL  
          L=LWET(LP)  
          AQ(L,KSZ(L))=0.314*DML(L,KSZ(L))**QQSQR(L,KSZ(L)) + AVOXY(L)*HPI(L)
          AQ(L,KC)    =0.314*DML(L,KC    )* QQSQR(L,KC    ) + AVOXY(L)*HPI(L)
        ENDDO
      ENDDO
      !$OMP END DO

    ELSE
      ! *** SQRT FILTER & NO FILTER
      !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
      DO ND=1,NDM  
        LF=(ND-1)*LDMWET+1  
        LL=MIN(LF+LDMWET-1,LAWET)
        
        DO K=2,KS
          DO LP=1,LLWET(K-1,ND)
            L=LKWET(LP,K-1,ND) 
            AQ(L,K)=0.5*( 0.314*( DML(L,K-1)*QQSQR(L,K-1)+DML(L,K  )*QQSQR(L,K) ) + AVOXY(L)*HPI(L) ) + 0.5*AQ(L,K)
          ENDDO
        ENDDO

        DO LP=LF,LL  
          L=LWET(LP)  
          AQ(L,KSZ(L))=0.5*( 0.314*DML(L,KSZ(L))**QQSQR(L,KSZ(L)) + AVOXY(L)*HPI(L) ) + 0.5*AQ(L,KSZ(L))
          AQ(L,KC)    =0.5*( 0.314*DML(L,KC    )* QQSQR(L,KC    ) + AVOXY(L)*HPI(L) ) + 0.5*AQ(L,KC)
        ENDDO
      ENDDO
      !$OMP END DO

    ENDIF
  ENDIF

  !$OMP END PARALLEL
  
  ! ****************************************************************************
  RETURN  

END  
