SUBROUTINE CALHDMF3(ISTL_)

  ! **  SUBROUTINE CALDMF CALCULATES THE HORIZONTAL VISCOSITY AND
  ! **  DIFFUSIVE MOMENTUM FLUXES. THE VISCOSITY, AH IS CALCULATED USING
  ! **  SMAGORINSKY'S SUBGRID SCALE FORMULATION PLUS A CONSTANT AHO
  !
  !**********************************************************************!
  ! *** SUBROUTINE CALEXP CALCULATES EXPLICIT MOMENTUM EQUATION TERMS

  ! CHANGE RECORD
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! ** 2015-12     PAUL M. CRAIG      ADDED OMP AND UPDATED TO F90
  ! ** 2015-12     PAUL M. CRAIG      ADOPTED AQEA ISHDMF>0 FOR 3TL

  !**********************************************************************!
  USE GLOBAL  
  !  USE OMP_LIB

  IMPLICIT NONE

  !**********************************************************************!
  !
  ! **VARIABLE DEFINITION 
  INTEGER, INTENT(IN) :: ISTL_
  INTEGER :: L,LP,LW,K,ND,LN,LS,LE,LNW,LNE,LSW,LSE
  REAL    :: TMPVAL,WVFACT,DTMPH,DTMPX,AHWVX,DTMP,FMDUX0,FMDUY0,FMDVY0,FMDVX0
  
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DYU1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DYV1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DXU1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:,:) :: DXV1
  REAL,SAVE,ALLOCATABLE,DIMENSION(:)   :: HMC

  IF( .NOT. ALLOCATED(DYU1) )THEN
    ALLOCATE(DYU1(LCM,KCM))  
    ALLOCATE(DYV1(LCM,KCM))  
    ALLOCATE(DXU1(LCM,KCM))  
    ALLOCATE(DXV1(LCM,KCM))  
    ALLOCATE(HMC(LCM))  
    
    DYU1=0.
    DYV1=0.
    DXU1=0.
    DXV1=0.
    HMC=0.
    FMDUX = 1E-8
    FMDUY = 1E-8
    FMDVY = 1E-8
    FMDVX = 1E-8
  ENDIF
  
  IF( ISDRY > 0 )THEN
    IF( LADRY > 0 )THEN
      DO K=1,KC
        DO LP=1,LADRY
          L=LDRY(LP)
          AH(L,K)  = AHOXY(L)
          AHC(L,K) = AHOXY(L)
          DXU1(L,K) = 0.0
          DXV1(L,K) = 0.0
          DYU1(L,K) = 0.0
          DYV1(L,K) = 0.0
          FMDUX(L,K) = 1E-16
          FMDUY(L,K) = 1E-16
          FMDVY(L,K) = 1E-16
          FMDVX(L,K) = 1E-16
        ENDDO  
      ENDDO 
    ENDIF
  ENDIF
  IF( ISTL_ /=3 )RETURN
  
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN
    IF( WVLSH > 0.0 .OR. WVLSX > 0.0 )THEN
      IF (N < NTSWV) THEN
        TMPVAL=FLOAT(N)/FLOAT(NTSWV)
        WVFACT=0.5-0.5*COS(PI*TMPVAL)
      ELSE
        WVFACT=1.0
      END IF
      AHWVX=WVFACT*WVLSX*WVPRD*WVPRD
    ENDIF
  ENDIF
  
  ! ****************************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)
  
  IF( ISDRY > 0 .AND. LADRY > 0 )THEN
    !$OMP DO PRIVATE(ND,K,LN,LP,L)
    DO ND=1,NDM  
      DO K=1,KC  
        LN=0
        DO LP=1,LLWET(K,ND)
          L=LKWET(LP,K,ND)  
          IF( LHDMF(L,K) )THEN
            LN = LN+1
            LKHDMF(LN,K,ND) = L
          ENDIF
        ENDDO
        LLHDMF(K,ND)=LN    ! *** NUMBER OF WET HDMF CELLS FOR THE CURRENT LAYER
      ENDDO
    ENDDO
    !$OMP END DO
  ENDIF
  
  !----------------------------------------------------------------------!
  ! **  CALCULATE HORIZONTAL VELOCITY SHEARS
  !$OMP DO PRIVATE(ND,K,LP,L,LS,LN,LE,LW) 
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLHDMF(K,ND)
        L=LKHDMF(LP,K,ND)  
        LS=LSC(L)
        LN=LNC(L)

        LE=LEC(L)
        LW=LWC(L)

        ! *** SHEAR ACROSS CELL CENTERS
        DXU1(L,K) = (U1(LE,K)-U1(L,K))/DXP(L)
        DYV1(L,K) = (V1(LN,K)-V1(L,K))/DYP(L)
        
        ! *** HORIZONTAL X COMPONENT SHEAR AT SW CORNER
        DYU1(L,K) = 2.*SUB(LS)*(U1(L,K)-U1(LS,K))/(DYU(L)+DYU(LS))
        IF( SVB(L) < 0.5 )DYU1(L,K) = 2.0*U1(L,K)/DYU(L)

        ! *** HORIZONTAL Y COMPONENT SHEAR AT SW CORNER
        DXV1(L,K) = 2.*SVB(LW)*(V1(L,K)-V1(LW,K))/(DXV(L)+DXV(LW))
        IF( SUB(L) < 0.5 )DXV1(L,K) = 2.0*V1(L,K)/DXV(L)

      ENDDO
    ENDDO
  ENDDO   ! *** END OF DOMAIN 
  !$OMP END DO
  
  !----------------------------------------------------------------------!
  ! **  CALCULATE HORIZONTAL VISCOSITY
  IF( AHD > 0.0 )THEN
    !$OMP DO PRIVATE(ND,K,LP,L,LN,LS,LE,LW,LNW,LNE,LSW,LSE) 
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLHDMF(K,ND)
          L=LKHDMF(LP,K,ND)  

          LN=LNC(L)
          LS=LSC(L)
          LE=LEC(L)
          LW=LWC(L)
          LNW=LNWC(L)
          LNE=LNEC(L)
          LSW=LSWC(L)
          LSE=LSEC(L)

          ! *** CELL CENTROID
          AH(L,K)  = AHOXY(L) + AHDXY(L)*DXP(L)*DYP(L)                                                     &
                                  *SQRT( 2.*DXU1(L,K)*DXU1(L,K)  + 2.*DYV1(L,K)*DYV1(L,K)                  &
                                        + 0.0625*(DYU1(L,K) + DYU1(LN,K) + DYU1(LE,K) + DYU1(LNE,K)        &
                                        +         DXV1(L,K) + DXV1(LN,K) + DXV1(LE,K) + DXV1(LNE,K))**2. )
          ! *** SW CORNER
          AHC(L,K) = AHOXY(L) + AHDXY(L)*0.0625*( (DXP(L) + DXP(LW) + DXP(LS) + DXP(LSW))**2. )            &
                                  *SQRT(  0.125*(DXU1(L,K) + DXU1(LW,K) + DXU1(LS,K) + DXU1(LSW,K))**2.    &
                                        + 0.125*(DYV1(L,K) + DYV1(LW,K) + DYV1(LS,K) + DYV1(LSW,K))**2.    &
                                        + DYU1(L,K)*DYU1(L,K) + 2.*DYU1(L,K)*DXV1(L,K) + DXV1(L,K)*DXV1(L,K) )
        ENDDO
      ENDDO
    ENDDO  ! *** END OF DOMAIN  
    !$OMP END DO
    
  ELSEIF( N < 10 .OR. ISWAVE == 2 .OR. ISWAVE == 4 )THEN
    ! *** ONLY NEED TO ASSIGN INITIALLY
    !$OMP SINGLE
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          AH(L,K)  = AHOXY(L)
          AHC(L,K) = AHOXY(L)
        ENDDO
      ENDDO
    ENDDO
    !$OMP END SINGLE
  END IF

  !----------------------------------------------------------------------!
  ! **  CALCULATE HORIZONTAL DIFFUSION DUE TO WAVE BREAKING
  IF( ISWAVE == 2 .OR. ISWAVE == 4 )THEN
    IF( WVLSH > 0.0 .OR. WVLSX > 0.0 )THEN
      !$OMP DO PRIVATE(ND,K,LP,L,LN,LS,LW,LSW,LSE,DTMPH,DTMPX,DTMP) 
      DO ND=1,NDM
        DO K=1,KC
          DO LP=1,LLHDMF(K,ND)
            L=LKHDMF(LP,K,ND)  
            DTMPH=(WVFACT*WV(L).DISSIPA(K))**0.3333
            DTMPX=WV(L).DISSIPA(K)/HP(L)
            AH(L,K) = AH(L,K) + WVLSH*DTMPH*HP(L) + AHWVX*DTMPX
          ENDDO
        ENDDO

        DO K=1,KC
          DO LP=1,LLHDMF(K,ND)
            L=LKHDMF(LP,K,ND)  
            LS=LSC(L)
            LW=LWC(L) 
            LSW=LSWC(L)

            DTMP = 0.25*( WV(L).DISSIPA(K) + WV(LW).DISSIPA(K) + WV(LS).DISSIPA(K) + WV(LSW).DISSIPA(K) )

            DTMPH = (WVFACT*DTMP)**0.3333
            DTMPX = DTMP/HMC(L)

            AHC(L,K) = AHC(L,K) + WVLSH*DTMPH*HMC(L) + AHWVX*DTMPX
          ENDDO
        ENDDO
      ENDDO  ! *** END OF DOMAIN  
      !$OMP END DO
    END IF
  END IF

  ! *** SW CORNER AVERAGE DEPTHS (H1C)
  !$OMP DO PRIVATE(ND,LP,L,LS,LW,LSW) 
  DO ND=1,NDM
    DO LP=1,LLHDMF(KC,ND)
      L=LKHDMF(LP,KC,ND)  
      LS=LSC(L)
      LW=LWC(L)
      LSW=LSWC(L)

      !ICNT   = 1
      H1C(L) = H1P(L)

      ! *** WEST
      IF( LW == LC )THEN
        H1C(L) = H1C(L) +  HMIN	       
      ELSE
        H1C(L) = H1C(L) +  H1P(LW)
      END IF

      ! *** SOUTH
      IF( LS == LC )THEN
        H1C(L) = H1C(L) +  HMIN	       
      ELSE 
        H1C(L) = H1C(L) +  H1P(LS)
      END IF

      ! *** SOUTHWEST
      IF( LSW == LC )THEN
        H1C(L) = H1C(L) +  HMIN	       
      ELSE
        H1C(L) = H1C(L) +  H1P(LSW)
      END IF

      H1C(L) = 0.25*H1C(L)
    ENDDO
  ENDDO  ! *** END OF DOMAIN  
  !$OMP END DO
  
  !----------------------------------------------------------------------!
  ! **  CALCULATE DIFFUSIVE MOMENTUM FLUXES
  !$OMP DO PRIVATE(ND,K,LP,L,LN,LS,LE,LW,FMDUX0,FMDUY0,FMDVY0,FMDVX0) 
  DO ND=1,NDM
    DO K=1,KC
      DO LP=1,LLHDMF(K,ND)
        L=LKHDMF(LP,K,ND)  

        LN=LNC(L)
        LS=LSC(L)
        LE=LEC(L)
        LW=LWC(L)

        FMDUX0 =  2.* DYP(L)*H1P(L)*AH(L,K)*DXU1(L,K)
        FMDUY0 = 0.5*( DXU(L) + DXU(LS) )*H1C(L)*AHC(L,K)*( DYU1(L,K) + DXV1(L,K) )
        FMDVY0 =  2.* DXP(L)*H1P(L)*AH(L,K)*DYV1(L,K)
        FMDVX0 = 0.5*( DYV(L) + DYV(LW) )*H1C(L)*AHC(L,K)*( DYU1(L,K) + DXV1(L,K) )
        FMDUX(L,K) = FMDUX0 ! SQRT(FMDUX(L,K)*MAX(FMDUX0,1E-8))
        FMDUY(L,K) = FMDUY0 ! SQRT(FMDUY(L,K)*MAX(FMDUY0,1E-8))
        FMDVY(L,K) = FMDVY0 ! SQRT(FMDVY(L,K)*MAX(FMDVY0,1E-8))
        FMDVX(L,K) = FMDVX0 ! SQRT(FMDVX(L,K)*MAX(FMDVX0,1E-8))
      ENDDO
    ENDDO
  ENDDO  ! *** END OF DOMAIN  
  !$OMP END DO
  !$OMP END PARALLEL
  
  IF( ISDRY > 0 .AND. NASPECT > 0 )THEN
    ! *** ZERO XY COMPONENT FOR CELLS WITH HIGH ASPECT RATIOS
    !$OMP PARALLEL DO DEFAULT(NONE) SHARED(NDM,KC,LLHDMF,LKHDMF,LASPECT,FMDUY,FMDVX,DXP,DYP) PRIVATE(ND,K,LP,L)
    DO ND=1,NDM
      DO K=1,KC
        DO LP=1,LLHDMF(K,ND)
          L = LKHDMF(LP,K,ND)
          IF( LASPECT(L) )THEN
            FMDUY(L,K) = 1E-8
            FMDVX(L,K) = 1E-8
          ENDIF
        ENDDO
      ENDDO
    ENDDO  
    !$OMP END PARALLEL DO
  ENDIF

  ! *** ZERO BOUNDARY CELL MOMENTUM DIFFUSION
  DO LP=1,NBCSOP
    L=LOBCS(LP)
    DO K=1,KC
      AHC(L,K)=0.0
      FMDUX(L,K)=1E-16
      FMDUY(L,K)=1E-16
      FMDVY(L,K)=1E-16
      FMDVX(L,K)=1E-16
    ENDDO
  ENDDO

  RETURN
  
END SUBROUTINE
