SUBROUTINE SEDZLJ_SHEAR
  
  ! CALCULATES Wave and Current Shear Stress Based on Log-Law of Cristofferson Jonsson
  ! 
  ! ORIGINAL:  May 24, 2007
  !  Craig Jones and Scott James
  ! REVISED: SIGMA-ZED AND OMP - 2016-11-07
  !  Paul M. Craig

  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER :: L,ND,LF,LL,LP
  INTEGER :: M1,M2
  INTEGER :: FZONE
  
  REAL(RKD)  :: MMW,SIGMAWV,JJW,KN,FWVTPP
  REAL(RKD)  :: VELMAG,VELANG,DELW,APROUGH
  REAL(RKD)  :: UTMP,VTMP
  REAL(RKD)  :: WVLENGTH,WVANGLE,WFTIM,EXCURSION
  REAL(RKD)  :: FC1,FC2,FWINDSQ,FC,FWW,FWVHT,SHEAR,GROWTH
  REAL(RKD)  :: FWINDS,FWINDD
  REAL(RKD)  :: TDIFF,WTM1,WTM2,AVGDEPTH
  
  REAL(RKD) ,SAVE,ALLOCATABLE,DIMENSION(:) :: WVFREQ
  REAL(RKD) ,SAVE,ALLOCATABLE,DIMENSION(:) :: WVORBIT
  REAL(RKD) ,SAVE,ALLOCATABLE,DIMENSION(:) :: ZBTEMP

  IF( .NOT. ALLOCATED(WVFREQ) )THEN
    ALLOCATE(WVFREQ(LCM))
    ALLOCATE(WVORBIT(LCM))
    ALLOCATE(ZBTEMP(LCM))
    WVFREQ = 0.0
    WVORBIT=0.0
    ZBTEMP=0.0
  ENDIF
  GROWTH = 0.1_8
  
  !**************************************************************************

  IF( TAUCONST == 0 )THEN
    ! *** Constant Tau not used.  Compute spatially variable Tau

    !**************************************************************************
    IF( ISWNWAVE == 1 )THEN
      ! Wind Wave Fetch
      
      ! Convert wind input into current wind info for wind-driven wave calcs
      IF( ISDYNSTP == 0 )THEN  
          WFTIM=DTSEDJ*FLOAT(N)/TCWSER(1)+TBEGIN*(TCON/TCWSER(1))  
      ELSE  
          WFTIM=TIMESEC/TCWSER(1)  
      ENDIF
      M2 = MWTLAST(1)
      DO WHILE (WFTIM > TSWND(1).TIM(M2))
        M2=M2+1
        IF( M2 > TSWND(1).NREC )THEN
          M2=TSWND(1).NREC
          EXIT
        ENDIF
      END DO
      MWTLAST(1) = M2 
      M1 = M2-1
      TDIFF=TSWND(1).TIM(M2)-TSWND(1).TIM(M1)  
      WTM1=(TSWND(1).TIM(M2)-WFTIM)/TDIFF  
      WTM2=(WFTIM-TSWND(1).TIM(M1))/TDIFF 
      FWINDS=WTM1*TSWND(1).VAL(M1,1)+WTM2*TSWND(1).VAL(M2,1)

      IF( FWINDS > 1.0 )THEN
        ! *** ONLY COMPUTE WIND WAVE IF WIDN SPEED > 1 M/S 
        IF( ABS(TSWND(1).VAL(M1,2)-TSWND(1).VAL(M2,2))<180.0 )THEN
          FWINDD=WTM1*TSWND(1).VAL(M1,2)+WTM2*TSWND(1).VAL(M2,2)
        ELSE
          IF( TSWND(1).VAL(M1,2) > TSWND(1).VAL(M2,2) )THEN
              FWINDD=WTM1*TSWND(1).VAL(M1,2)+WTM2*(TSWND(1).VAL(M2,2)+360.0)
          ELSE
              FWINDD=WTM1*(TSWND(1).VAL(M1,2)+360)+WTM2*TSWND(1).VAL(M2,2)
          ENDIF
          IF( FWINDD >= 360.0)FWINDD=FWINDD-360.0 
        ENDIF
        ! Convert wind into direction it is blowing "from"
        IF( FWINDD <= 180.0 )THEN  
          FWINDD=FWINDD+180.0  
          IF( FWINDD == 360.0)FWINDD=0.0         
        ELSE  
          FWINDD=FWINDD-180.0  
          IF( FWINDD == 360.0)FWINDD=0.0 
        ENDIF
        ! Calculate which of the 8 wind zones (FWZONE) the wind is coming from
        ! Also the Waveangle CCW from East.  Waveangle is center of sector.
        IF( FWINDD >= 337.5 .OR. FWINDD<22.5 )THEN
          FZONE=1
          WVANGLE=4.712
        ELSEIF( FWINDD >= 22.5 .AND. FWINDD<67.5 )THEN
          FZONE=2
          WVANGLE=3.927
        ELSEIF( FWINDD >= 67.5 .AND. FWINDD<112.5 )THEN
          FZONE=3
          WVANGLE=3.142
        ELSEIF( FWINDD >= 112.5 .AND. FWINDD<157.5 )THEN
          FZONE=4
          WVANGLE=2.356
        ELSEIF( FWINDD >= 157.5 .AND. FWINDD<202.5 )THEN
          FZONE=5
          WVANGLE=1.571
        ELSEIF( FWINDD >= 202.5 .AND. FWINDD<247.5 )THEN
          FZONE=6
          WVANGLE=0.7854
        ELSEIF( FWINDD >= 247.5 .AND. FWINDD<292.5 )THEN
          FZONE=7
          WVANGLE=0.
        ELSEIF( FWINDD >= 292.5 .AND. FWINDD<337.5 )THEN
          FZONE=8
          WVANGLE=5.4978
        ENDIF
        FWINDSQ = FWINDS*FWINDS
        !Calculate Domain Average Depth
        !Needs to be calculated along each fetch
        !This is sufficient for small systems
        AVGDEPTH = SUM(HP(2:LA))/FLOAT(LA-1)
        
        !Calculate wave height, period, orbital velocity, and length
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,FWVHT,FC1,FC2,FWVTPP,WVLENGTH)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)
            FC1=(FWINDSQ/9.8)*0.283*TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75)
            FC2=TANH(0.0125*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.42/TANH(0.530*(9.8*AVGDEPTH/FWINDSQ)**0.75))   
            FWVHT=MIN(HP(L),FC1*FC2)
              
            FC1=(FWINDS/9.8)*7.54*TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375)
            FC2=TANH(0.077*(9.8*FWDIR(L,FZONE)/FWINDSQ)**0.25/TANH(0.833*(9.8*AVGDEPTH/FWINDSQ)**0.375))   
            FWVTPP=FC1*FC2
            WVFREQ(L)=2.0*PI/FWVTPP
              
            WVLENGTH=FWVTPP*SQRT(9.8*HP(L)/FC2)
            WVLENGTH=MAX(1.0,WVLENGTH)
            WVORBIT(L)=MAX(0.01,PI*FWVHT/(FWVTPP*SINH(HP(L)*2.0*PI/WVLENGTH)))
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ELSE
        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)
            WVFREQ(L)=0.0
            WVORBIT(L)=0.0
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
      ENDIF
      
    ELSEIF( ISWNWAVE == 2 )THEN
      !Read in EFDC STWAVE Wave Field
      NWVCOUNT=NWVCOUNT+DTSEDJ/3600
      IF( NWVCOUNT == STWVTIM )THEN
        NWVCOUNT=0
        STINC=STINC+1

        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,FWVHT,FC1,FC2,WVLENGTH,EXCURSION)
        DO ND=1,NDM  
          LF=(ND-1)*LDMWET+1  
          LL=MIN(LF+LDMWET-1,LAWET)
          DO LP=LF,LL
            L=LWET(LP)
            IF( STINC > STWVNUM ) EXIT
            IF( STWVTP(L,STINC) > 0.0 )THEN
              WVFREQ(L)=2.0*PI/STWVTP(L,STINC)
              FWVHT=MIN(HP(L),STWVHT(L,STINC))
              FC1=9.8*STWVTP(L,STINC)**2
              FC2=2.*PI*SQRT(TANH(4.*PI**2*HP(L)/(STWVTP(L,STINC)**2*9.8)))
              WVLENGTH=FC1/FC2
              EXCURSION=FWVHT/(2.*SINH((2.*PI/WVLENGTH)*HP(L)))
              WVORBIT(L)=EXCURSION*WVFREQ(L)
              WVORBIT(L)=MAX(0.01,WVORBIT(L))
              WVANG(L)=STWVDR(L,STINC)
            ELSE
              WVFREQ(L)=0.0
              WVORBIT(L)=0.0
              WVANG(L)=0.0
            ENDIF
          ENDDO
        ENDDO
        !$OMP END PARALLEL DO
        
      ENDIF
    ENDIF    

    !*************************************************************************
    !Begin shear stress calculations

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,LF,LL,LP,L,KN,UTMP,VTMP,VELMAG,FC,VELANG,FWW,SIGMAWV,MMW,JJW,DELW,APROUGH,SHEAR)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      ! *** Set up roughness, velocities, and angles 
      IF( ZBSKIN <= 0. )THEN
        ! *** VARIABLE BOTTOM SKIN FRICTION.  ZBTEMP(L) = ZBR(L)  [meters]
        DO LP=LF,LL
          L=LWET(LP)
          ZBTEMP(L) = MAX(D50AVG(L)*1.E-6,1.E-8)
        ENDDO
      ELSE
        DO LP=LF,LL
          L=LWET(LP)
          ZBTEMP(L) = ZBSKIN*1.E-6
        ENDDO
      ENDIF
      
      DO LP=LF,LL
        L=LWET(LP)
        IF( HP(L) < HPMIN )THEN
          !TAU(L) = 0.0
          !USTAR(L) = 0.0
          CYCLE
        ENDIF
        
        IF( ISWAVE == 2 ) ZBTEMP(L)=0.0001
        KN = 30.0*ZBTEMP(L)
        
        ! Calculate Cell Centroid Depth Averaged Velocity Magnitude in cm/s
        UTMP = 100.0*STCUV(L)*( UHE(LEC(L))+UHE(L) )/( SUB(LEC(L))*HU(LEC(L))+HU(L) ) + 1.0E-12
        VTMP = 100.0*STCUV(L)*( VHE(LNC(L))+VHE(L) )/( SVB(LNC(L))*HV(LNC(L))+HV(L) )
        VELMAG = SQRT(UTMP**2+VTMP**2)
          
        ! Calculate Initial Friction Factors
        FC = ( 0.42/LOG(HP(L)/(2.0*ZBTEMP(L))) )**2  
          
        ! Current only Shear Stress
        IF( ISWNWAVE == 0 .AND. UWVSQ(L) == 0.0 .OR. ISWAVE ==  0 )THEN
          ! *** Compute Bed Shear but Limit Growth to 10% of Previous Value
          SHEAR = FC*VELMAG**2
          IF( TAU(L) <= 0.0 )THEN
            TAU(L) = GROWTH*SHEAR
          ELSEIF( (1.+GROWTH)*SHEAR > TAU(L) )THEN
            TAU(L) = TAU(L) + GROWTH*SHEAR
          ELSE
            TAU(L) = SHEAR
          ENDIF
        ELSE
          ! Calculate Combined Wave Friction Factor
          ! Calculate Current Angle CCW From X axis
          IF( UTMP>0.0 .AND. VTMP>0.0 )THEN
            VELANG=ATAN(VTMP/UTMP)
          ELSEIF( UTMP<0.0 )THEN
            VELANG=ATAN(VTMP/UTMP)+PI
          ELSEIF( UTMP>0.0 .AND. VTMP<0.0 )THEN
            VELANG=2*PI+ATAN(VTMP/UTMP)
          ELSEIF( UTMP == 0.0 )THEN
            VELANG=SIGN(0.5*PI,VTMP)
          ENDIF
              
          ! Set Orbital velocity in m/s and waveangle and frequency
          IF( ISWAVE > 0 )THEN
              WVFREQ(L)=WVFRQ
              WVORBIT(L)=SQRT(UWVSQ(L))
              WVANG(L)=WV(L).DIR    
          ELSEIF( ISWNWAVE == 1 )THEN
            WVANG(L)=WVANGLE
          ENDIF
              
          ! Calculate wave friction factor
          FWW=2.0*(0.0747*(KN*WVFREQ(L)/WVORBIT(L)))**0.66
          SIGMAWV = FC/FWW*(VELMAG/(WVORBIT(L)*100.0))**2
          MMW=SQRT(1.0+SIGMAWV**2+2.0*SIGMAWV*ABS(COS(VELANG-WVANG(L))))
          JJW=WVORBIT(L)/(KN*WVFREQ(L))*SQRT(MMW*FWW/2.0)
          FWW=MMW*0.15/JJW
              
          !Calculate wave boundary layer info
          DELW=KN*0.273*SQRT(JJW)
          APROUGH=30.0*DELW*EXP(-5.62*DELW/KN*SQRT(SIGMAWV/MMW))
          !Calculate new current friction factor
          FC = 2.0*(1.0/(2.38*LOG(30.0*HP(L)/(2.718*KN))-2.38*LOG(APROUGH/KN)))**2
              
          !Iterate once more
          SIGMAWV = FC/FWW*(VELMAG/(WVORBIT(L)*100.0))**2
          MMW = SQRT(1.0+SIGMAWV**2+2.0*SIGMAWV*ABS(COS(VELANG-WVANG(L))))
          JJW = WVORBIT(L)/(KN*WVFREQ(L))*SQRT(MMW*FWW/2.0)
          FWW = MMW*0.15/JJW
              
          ! *** Calculate total wave and current shear stress (dynes/cm^2)
          TAU(L) = 0.5*FWW*WVORBIT(L)**2*10000.0*MMW
          TAUB(L) = 0.1*TAU(L)             ! *** Conversion from dynes/cm^2 to Pa
          USTAR(L) = SQRT(TAUB(L)/1000.0)  ! *** USTAR=SQRT(Tau/RhoH2O) and RhoH2O is 1000 kg/cm^3. 
        ENDIF
      ENDDO
    ENDDO
    !$OMP END PARALLEL DO
    
    !*****************************
     
  ELSE
    ! *** Set constant tau is TAUCONST (dynes/cm^2) is greater than 0
    DO L=2,LA
      TAU(L) = TAUCONST
      TAUB(L) = 0.1*TAU(L)
      USTAR(L) = SQRT(TAUB(L)/1000.)
    ENDDO
  ENDIF
  
  !************************
  RETURN
  
END SUBROUTINE SEDZLJ_SHEAR
