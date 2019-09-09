SUBROUTINE SEDZLJ_SLOPE
  
  USE GLOBAL
  IMPLICIT NONE
  
  ! CALCULATES slope change on shear
  ! 
  ! REVISION DATE :  Nov 19, 2009
  ! Craig Jones and Scott James
  !**************************************************************************
  ! Check to see if we've set a constant Tau
  INTEGER :: L
  REAL :: DZBETA,DZTHETA,COSB,SINB,COST,TANT,MINTAU
  REAL :: ALPHA,AA,BB,FCO,FB,FG,TEMP1,TUNEP,TUNER
  
  DOUBLE PRECISION ::  XDIST,YDIST,DELZX,DELZY,VELMAG,DELX,DELY

  AA=1.0             !Geometry factor (Lick, 2009, Figure 3.25)
  BB=1.0             !Geometry factor (Lick, 2009, Figure 3.25)
  ALPHA=ATAN(AA/BB)  !Geometry factor (Lick, 2009, p.96)
  FCO=0.0            !Cohesive force = 7 N/m^2 * D50(K)^2, but zero for non-cohesives (Lick, 2009, p.89)
  FB=0.0             !Binding force, but zero for non-cohesives (Lick, 2009, p.96)
  FG=1.0             !Gravitational force (Lick, 2008, Eqn. 3.8), arbitrarily set to 1 because it is irrelevant for non-cohesives
  TUNEP=1.0          !Bedload pitch tuning factor
  TUNER=1.0          !Bedload roll  tuning factor
  MINTAU=10000.
  DO L=NNONCO,NSCM
    MINTAU = MIN(TCRE(L),MINTAU)
  ENDDO
  
  DO L=2,LA
    IF( TAU(L) > MINTAU )THEN
      ! *** For Beta Angle "Pitch" and Theta Angle "Roll"
      IF( SUB(LEC(L)) /= 0. .AND. SUB(L) /= 0. .AND. SUB(LWC(L)) /= 0. )THEN      ! *** Central differencing
        DELZX = BELV(LEC(L))-BELV(LWC(L))                      
        XDIST = 0.5*DXP(LEC(L))+DXP(L)+0.5*DXP(LWC(L))                            ! *** Across 3 cells
      ELSEIF( SUB(LEC(L)) == 0. .AND. SUB(L) /= 0. .AND. SUB(LWC(L)) /= 0. )THEN  ! *** Forward differencing
        DELZX = BELV(L)-BELV(LWC(L))
        XDIST = 0.5*(DXP(L)+DXP(LWC(L)))                                          ! *** Across 2 cells
      ELSEIF( SUB(LEC(L)) /= 0. .AND. SUB(L) /= 0. .AND. SUB(LWC(L)) == 0. )THEN
        DELZX=BELV(LEC(L))-BELV(L)
        XDIST=0.5*(DXP(LEC(L))+DXP(L))                                            ! *** Across 2 cells
      ELSE
        DELZX = 0.0                                                               ! *** No active cells on either side
        XDIST = 1.0
      ENDIF
    
      IF( SVB(LNC(L))/=0. .AND. SVB(L) /= 0. .AND. SVB(LSC(L)) /= 0. )THEN        ! *** Central differencing
        DELZY=BELV(LNC(L))-BELV(LSC(L))
        YDIST=0.5*DYP(LNC(L))+DYP(L)+0.5*DYP(LSC(L))                              ! *** Across 3 cells
      ELSEIF( SVB(LNC(L)) == 0. .AND. SVB(L) /= 0. .AND. SVB(LSC(L)) /= 0. )THEN  ! *** Forward differencing 
        DELZY=BELV(L)-BELV(LSC(L))
        YDIST=0.5*(DYP(L)+DYP(LSC(L)))                                            ! *** Across 2 cells
      ELSEIF( SVB(LNC(L)) /= 0. .AND. SVB(L) /= 0. .AND. SVB(LSC(L)) == 0. )THEN  ! *** Forward differencing
        DELZY=BELV(LNC(L))-BELV(L)
        YDIST=0.5*(DYP(LNC(L))+DYP(L))                                            ! *** Across 2 cells
      ELSE
        DELZY = 0.0                                                               ! *** No active cells on either side
        YDIST = 1.0
      ENDIF
    
      VELMAG = SQRT(U(L,KSZU(L))**2 + V(L,KSZV(L))**2)  ! *** Local flow speed
      IF( VELMAG > 0.0 )THEN                            ! *** Velocity directional vectors
        DELX = U(L,KSZU(L))/VELMAG
        DELY = V(L,KSZV(L))/VELMAG  
      ELSE           
        DELX = 0.0
        DELY = 0.0
      ENDIF  
    
      DZBETA  = DELX/XDIST*DELZX + DELY/YDIST*DELZY  ! *** Change in bed elevation for the velocity pitch angle
      DZTHETA = DELY/XDIST*DELZX + DELX/YDIST*DELZY  ! *** Change in bed elevation fot the velocity roll  angle
      COSB = 1.0/SQRT(1.0+DZBETA**2)                 ! *** Cosine of pitch angle
      SINB = DZBETA/SQRT(1.0+DZBETA**2)              ! *** Sine of pitch angle
      COST = 1.0/SQRT(1.0+DZTHETA**2)                ! *** Cosine of roll angle
      TANT = DZTHETA                                 ! *** Tangent of roll angle
      TEMP1 = (1.0+(FCO+FB)/(FG*COSB*COST))**2-(TANT/TAN(ALPHA))**2 
    
      ! *** Slope erosion from Willy Lick Work
      IF( TEMP1 > 0.0 )THEN                          ! *** Calculate shear stress and erosion scaling factor
        SH_SCALE(L)=MAX(0.1,BB/AA*SINB+COSB*COST*SQRT(TEMP1)) ! *** Lick, 2009, Eqn.3.36
      ELSE
        SH_SCALE(L)=0.1                              ! *** Minimum value
      ENDIF
    
      ! *** Slope bedload factors from Lesser et al. summary
      ALPHA_PX(L) = 1.0-TUNEP*(TAN(ALPHA)/(COS(ATAN(DELZX/XDIST))*(TAN(ALPHA)-DELZX/XDIST))-1.0) ! *** Bedload velocity pitch angle x-correction
      ALPHA_RX(L,1:NSCM) = -TUNER*SQRT(TCRE(1:NSCM)/TAU(L))*DELZY/YDIST                          ! *** Bedload velocity roll  angle x-correction
      ALPHA_PY(L) = 1.0-TUNEP*(TAN(ALPHA)/(COS(ATAN(DELZY/YDIST))*(TAN(ALPHA)-DELZY/YDIST))-1.0) ! *** Bedload velocity pitch angle y-correction
      ALPHA_RY(L,1:NSCM) = -TUNER*SQRT(TCRE(1:NSCM)/TAU(L))*DELZX/XDIST                          ! *** Bedload velocity roll  angle y-correction
    ELSE
      SH_SCALE(L) = 1.0
      ALPHA_PX(L) = 1.0
      ALPHA_RX(L,1:NSCM) = 0.
      ALPHA_PY(L) = 1.0
      ALPHA_RY(L,1:NSCM) = 0.
    ENDIF
  ENDDO
  
  RETURN

END SUBROUTINE SEDZLJ_SLOPE
