SUBROUTINE CALPSER

  ! CHANGE RECORD
  ! ** SUBROUTINE CALPSER UPDATES TIME VARIABLE SURFACE ELEVATION
  ! ** BOUNDARY CONDITIONS

  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER   :: NS,M2,M1
  REAL(RKD) :: TIME,TDIFF
  REAL      :: WTM1,WTM2
  REAL      :: Y1,Y2,Y3,Y4
  
  PSERT(0)=0.
  
  DO NS=1,NPSER
    TIME = TIMESEC/DBLE(TCPSER(NS))

    M2 = MPTLAST(NS)
    DO WHILE (TIME > TPSER(M2,NS))
      M2=M2+1
      IF( M2 > MPSER(NS) )THEN
        M2=MPSER(NS)
        EXIT
      ENDIF
    END DO
    MPTLAST(NS) = M2  
    M1 = M2-1
    
    TDIFF = TPSER(M2,NS) - TPSER(M1,NS)
    
    IF( INTPSER(NS) == 0 .OR. M1 < 2 .OR. M2 > MPSER(NS)-1 )THEN
      ! *** LINEAR INTERPOLATION
      WTM1 = (TPSER(M2,NS) - TIME)/TDIFF 
      WTM2 = (TIME-TPSER(M1,NS))/TDIFF
      PSERT(NS)  = WTM1*PSER(M1,NS)  + WTM2*PSER(M2,NS)  + PDGINIT  ! *** ADD OFFSET
      PSERST(NS) = WTM1*PSERS(M1,NS) + WTM2*PSERS(M2,NS) + PDGINIT  ! *** ADD OFFSET
    ELSE
      ! *** CATMULL–ROM SPLINE
      WTM1 = (TIME - TPSER(M1,NS))/TDIFF
      Y1 = PSER(M1-1,NS)
      Y2 = PSER(M1,NS)
      Y3 = PSER(M2,NS)
      Y4 = PSER(M2+1,NS)
      PSERT(NS) = 0.5*( (2.*Y2) + (-Y1 + Y3)*WTM1 + (2.*Y1 - 5.*Y2 + 4.*Y3 - Y4)*WTM1**2 + (-Y1 + 3.*Y2 - 3.*Y3 + Y4)*WTM1**3 )
      
      ! *** DSI Note - Uncomment out the following for QC tests of the spline interpolation 
      !if(ns==3) WRITE(100+NS,'(F12.6,F10.3,2I10,F8.5,2F10.3)') TIME,PSERT(NS)/g, M1,M2,WTM1, PSER(M1,NS)/g, PSER(M2,NS)/g 
      !IF( MOD(NITER,1000)==0 .AND. NS==3 ) WRITE(*,'(F12.6,F10.3,2I10,F8.5,2F10.3)') TIME,PSERT(NS)/g, M1,M2,WTM1, PSER(M1,NS)/g, PSER(M2,NS)/g
      
      Y1 = PSERS(M1-1,NS)
      Y2 = PSERS(M1,NS)
      Y3 = PSERS(M2,NS)
      Y4 = PSERS(M2+1,NS)
      PSERST(NS) = 0.5*( (2.*Y2) + (-Y1 + Y3)*WTM1 + (2.*Y1 - 5.*Y2 + 4.*Y3 - Y4)*WTM1**2 + (-Y1 + 3.*Y2 - 3.*Y3 + Y4)*WTM1**3 )

    ENDIF
  ENDDO
  
  RETURN
  
END

