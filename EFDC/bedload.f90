SUBROUTINE BEDLOAD(NX,NS)

  !****************************************************************************
  ! *** SUBROUTINE CALSND CALCULATES NONCOHESIVER SEDIMENT SETTLING,                                                      
  ! *** DEPOSITION AND RESUSPENSION AND IS CALLED FOR SSEDTOX                                                             
  ! *** NOT USED BY SEDZLJ
  !
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  !   2011-03       Paul M. Craig      Converted to F90, added OMP

  USE GLOBAL
  IMPLICIT NONE                                                                                                          
  
  INTEGER, INTENT(IN) :: NX,NS
  INTEGER :: L,NSB,LUTMP,LDTMP,LE,LS,LN,LW,ND,LF,LL,LP
  REAL    :: SLOPE,ASNDFBL,SNDBTMP,SNDFBLM, FLUXFAC
  REAL    :: UCELLCTRM,VCELLCTRM,SHIELDS,BDLDTMPB,CSHIELDSC                                                        
  REAL    :: BDLDTMPA,CSHIELDS,TMPVAL,BDLDTMPP,BDLDTMP,XOUT,YOUT                                                           
  REAL, EXTERNAL :: FSEDMODE
  REAL, EXTERNAL :: FSBDLD
  
  ! *** ZERO BOUNDARY BEDLOAD FLUXES
  IF( NSBDLDBC > 0 )THEN
    DO NSB=1,NSBDLDBC
      LUTMP = LSBLBCU(NSB)
      QSBDLDOT(LUTMP,NX) = 0.       ! *** BC OUTFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
      QSBDLDIN(LUTMP,NX) = 0.       ! *** BC INFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
      LDTMP = LSBLBCD(NSB)
      IF( LDTMP > 0 )THEN
        QSBDLDOT(LDTMP,NX) = 0.     ! *** BC OUTFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
        QSBDLDIN(LDTMP,NX) = 0.     ! *** BC INFLOW OR RECIRCULATION BOUNDARY CONDITION (G/S)
      ENDIF
    ENDDO
  ENDIF
  
  ! *** BED LOAD TRANSPORT HORIZONTAL LOOP                                                                                

  !$OMP PARALLEL DEFAULT(SHARED)
  
  !$OMP DO PRIVATE(ND,LF,LL,L)
  DO ND=1,NDM  
    LF=2+(ND-1)*LDM  
    LL=MIN(LF+LDM-1,LA)

    ! *** ZERO BED LOAD TRANSPORTS                                                                                    
    DO L=LF,LL
      QSBDLDP(L) = 0.         ! *** CELL CENTER BED LOAD TRANSPORT RATE (G/M/S)
      QSBDLDX(L,NX) = 0.      ! *** U FACE SND FLUX DUE TO BEDLOAD (G/S)
      QSBDLDY(L,NX) = 0.      ! *** V FACE SND FLUX DUE TO BEDLOAD (G/S) 
      SNDFBL(L,NX) = 0.       ! *** BED/WATER INTERFACE SND FLUX DUE TO BEDLOAD (G/M2/S)
    ENDDO
  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  !*********************************************************************************
  !*** COMPUTE THE CELL CENTERED BEDLOAD TRANSPORT RATES USING THE SPECIFIED OPTION
      
  IF( ISBDLD(NS) == 0 )THEN
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING GENERIC BED LOAD EQUATION

    !$OMP DO PRIVATE(ND,LF,LL,L,LS,LN) &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA,BDLDTMPB) 
    DO ND=1,NDM  
      LF=2+(ND-1)*LDM  
      LL=MIN(LF+LDM-1,LA)

      DO L=LF,LL
        IF( LMASKDRY(L) )THEN
          CSHIELDS = TCSHIELDS(NS)
          IF( ISEDEFF == 2 )THEN
            TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
            CSHIELDS = TMPVAL*CSHIELDS
          ENDIF
          BDLDTMPP = SBDLDP(NX)
          BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM
          IF( BDLDTMPP > 0.0 )THEN
            FACBEDL(L) = FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
            SHIELDS = TAUBSND(L)/GPDIASED
            IF( SHIELDS > CSHIELDS )THEN
              IF( SBDLDA(NX) > 0.0 )THEN
                BDLDTMPA = (SBDLDG1(NX)*SHIELDS-SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
              ELSE
                BDLDTMPA = 1.0
              ENDIF
              IF( SBDLDB(NX) > 0.0 )THEN
                BDLDTMPB = (SBDLDG3(NX)*SQRT(SHIELDS)-SBDLDG4(NX)*SQRT(CSHIELDS))**SBDLDB(NX)
              ELSE
                BDLDTMPB = 1.0
              ENDIF
              QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA*BDLDTMPB     ! *** (G/M/S)
              IF( ISEDEFF == 1 )THEN
                TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
                QSBDLDP(L) = TMPVAL*QSBDLDP(L)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    
    ENDDO   ! *** END OF DOMAIN LOOP FOR GENERIC BED LOAD EQUATION
    !$OMP END DO

  ELSEIF( ISBDLD(NS) == 1 )THEN
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING VAN RIJN BED LOAD EQUATION

    !$OMP DO PRIVATE(ND,LF,LL,LP,L)   &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,CSHIELDSC,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        CSHIELDS = TCSHIELDS(NS)
        IF( ISEDEFF == 2 )THEN
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        ENDIF
        IF( ISNDAL == 1 )THEN
          TMPVAL = LOG10(19.*DIASED/SEDDIA50(L,KBT(L)))
          TMPVAL = 1.66667/(TMPVAL**2)
          CSHIELDSC = CSHIELDS50(L)*TMPVAL
        ELSE
          CSHIELDSC = TCSHIELDS(NS)
        ENDIF
        IF( ISEDEFF == 2 )THEN
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDSC = TMPVAL*CSHIELDS
        ENDIF
        BDLDTMPP = FSBDLD(DIASED,GPDIASED,SEDDIA50(L,KBT(L)),HP(L),PEXP(L,NX),PHID(L,NX),CSHIELDS,SBDLDP(NX),ISBDLD(NS))
        IF( ISNDAL == 1 )THEN
          BDLDTMPP = ((DIASED/SEDDIA50(L,KBT(L)))**0.3)*BDLDTMPP
        ENDIF
        BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM
        FACBEDL(L) = FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
        SHIELDS = TAUBSND(L)/GPDIASED
        IF( SHIELDS > CSHIELDSC .OR. ISBDLD(NS) > 1 )THEN
          BDLDTMPA = (SBDLDG1(NX)*SHIELDS -SBDLDG2(NX)*CSHIELDSC)**SBDLDA(NX)
          QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA
          IF( ISEDEFF == 1 )THEN
            TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
            QSBDLDP(L) = TMPVAL*QSBDLDP(L)
          ENDIF
        ENDIF
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP FOR VAN RIJN BED LOAD EQUATION
    !$OMP END DO

  ELSEIF( ISBDLD(NS) == 2 )THEN
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING ENGELUND-HANSEN

    !$OMP DO PRIVATE(ND,LF,LL,LP,L)   &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      CSHIELDS = 0.
      IF( IBLTAUC(NS) == 1 ) CSHIELDS = TCSHIELDS(NS)
      DO LP=LF,LL
        L=LWET(LP)  

        IF( IBLTAUC(NS) == 2 ) CSHIELDS = SEDDIA50(L,KBT(L))*CSHIELDS50(L)/DIASED
        IF( IBLTAUC(NS) == 3 ) CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
        IF( ISEDEFF == 2 )THEN
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        ENDIF
        SHIELDS = TAUBSND(L)/GPDIASED
        IF( SHIELDS >= CSHIELDS )THEN
          BDLDTMPP = FSBDLD(DIASED,GPDIASED,SEDDIA50(L,KBT(L)),HP(L),PEXP(L,NX),PHID(L,NX),CSHIELDS,SBDLDP(NX),ISBDLD(NS))
          IF( HGDH(L) > 0.0) BDLDTMPP = BDLDTMPP/(HGDH(L)**0.333)
          BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM
          FACBEDL(L) = FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
          BDLDTMPA = (SBDLDG1(NX)*SHIELDS)**SBDLDA(NX)
          QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA
          IF( ISEDEFF == 1 )THEN
            TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
            QSBDLDP(L) = TMPVAL*QSBDLDP(L)
          ENDIF
        ENDIF
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP FOR ENGELUND-HANSEN SECTION
    !$OMP END DO

  ELSEIF( ISBDLD(NS) == 3 )THEN
    ! *** CALCULATE CELL CENTER TRANSPORT RATES USING WU, WANG, AND JIA

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LS,LN) &
    !$OMP    PRIVATE(CSHIELDS,TMPVAL,BDLDTMPP,BDLDTMP,SHIELDS,BDLDTMPA) 
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        CSHIELDS = TCSHIELDS(NS)
        BDLDTMPP = FSBDLD(DIASED,GPDIASED,SEDDIA50(L,KBT(L)),HP(L),PEXP(L,NX),PHID(L,NX),CSHIELDS,SBDLDP(NX),ISBDLD(NS))
        CSHIELDS = 0.03*((PHID(L,NX)/PEXP(L,NX))**0.6)
        IF( ISEDEFF == 2 )THEN
          TMPVAL = 1.+(COEHEFF2-1.)*( 1.-EXP(-COEHEFF*FRACCOH(L,KBT(L))) )
          CSHIELDS = TMPVAL*CSHIELDS
        ENDIF
        BDLDTMP = SQRT(GPDIASED)*DIASED/DSEDGMM
        FACBEDL(L) = FSEDMODE(WSETA(L,0,NS),USTAR(L),USTARSND(L),RSNDM(NX),ISNDM1(NX),ISNDM2(NX),1)
        SHIELDS = TAUBSND(L)/GPDIASED
        IF( SHIELDS > CSHIELDS )THEN
          BDLDTMPA = (SBDLDG1(NX)*SHIELDS-SBDLDG2(NX)*CSHIELDS)**SBDLDA(NX)
          QSBDLDP(L) = FACBEDL(L)*VFRBED(L,KBT(L),NS)*BDLDTMP*BDLDTMPP*BDLDTMPA
          IF( ISEDEFF == 1 )THEN
            TMPVAL = EXP(-COEHEFF*FRACCOH(L,KBT(L)))
            QSBDLDP(L) = TMPVAL*QSBDLDP(L)
          ENDIF
        ENDIF
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP FOR WU, WANG, AND JIA SECTION
    !$OMP END DO
    
  ENDIF  ! *** END OF CELL CENTERED BEDLOAD TRANSPORT CALCULATIONS
  
  ! ********************************************************************************
  ! *** CALCULATE CELL FACE TRANSPORT RATES (G/M/S)

  IF( ISBLFUC == 0 )THEN
    ! *** CALCULATE CELL FACE TRANSPORT RATES BY DOWN WIND PROJECTION

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LW,LS)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        LW = LWC(L)
        LS = LSC(L)
        
        IF( UCELLCTR(LW) > 0.0 )THEN
          QSBDLDX(L,NX) = SUB(L)*QSBDLDP(LW)*UCELLCTR(LW)
        ELSEIF( UCELLCTR(L) < 0.0 )THEN
          QSBDLDX(L,NX) = SUB(L)*QSBDLDP(L)*UCELLCTR(L)
        ENDIF
        
        IF( VCELLCTR(LS) > 0.0 )THEN
          QSBDLDY(L,NX) = SVB(L)*QSBDLDP(LS)*VCELLCTR(LS)
        ELSEIF( VCELLCTR(L) < 0.0 )THEN
          QSBDLDY(L,NX) = SVB(L)*QSBDLDP(L)*VCELLCTR(L)
        ENDIF

      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  ELSEIF( ISBLFUC == 1 )THEN
    ! *** CALCULATE CELL FACE TRANSPORT RATES BY DOWN WIND PROJECTION WITH CORNER EFFECTS CORRECTION                                                                                    

    !$OMP DO PRIVATE(ND,LF,LL,L,LE,LN)       &
    !$OMP    PRIVATE(UCELLCTRM,VCELLCTRM)
    DO ND=1,NDM  
      LF=2+(ND-1)*LDM  
      LL=MIN(LF+LDM-1,LA)
      
      DO L=LF,LL
        IF( LMASKDRY(L) )THEN
          LE = LEC(L)
          LN = LNC(L)
          IF( UCELLCTR(L) >= 0.0 .AND. VCELLCTR(L) >= 0.0 )THEN
            UCELLCTRM = SUB(LE)*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(LN)*ABS(VCELLCTR(L))
            QSBDLDX(LE,NX) = SUB(LE)*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(LN,NX) = SVB(LN)*QSBDLDP(L)*VCELLCTR(L)
            IF( UCELLCTRM < 1.0E-9 )THEN
              QSBDLDY(LN ,NX) = SVB(LN )*SIGN(QSBDLDP(L),VCELLCTR(L))
            ENDIF
            IF( VCELLCTRM < 1.0E-9 )THEN
              QSBDLDX(LE,NX) = SUB(LE)*SIGN(QSBDLDP(L),UCELLCTR(L))
            ENDIF
          ENDIF
          IF( UCELLCTR(L) >= 0.0 .AND. VCELLCTR(L) < 0.0 )THEN
            UCELLCTRM = SUB(LE)*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(L )*ABS(VCELLCTR(L))
            QSBDLDX(LE,NX) = SUB(LE)*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(L ,NX) = SVB(L )*QSBDLDP(L)*VCELLCTR(L)
            IF( UCELLCTRM < 1.0E-9 )THEN
              QSBDLDY(L  ,NX) = SVB(L  )*SIGN(QSBDLDP(L),VCELLCTR(L))
            ENDIF
            IF( VCELLCTRM < 1.0E-9 )THEN
              QSBDLDX(LE,NX) = SUB(LE)*SIGN(QSBDLDP(L),UCELLCTR(L))
            ENDIF
          ENDIF
          IF( UCELLCTR(L) < 0.0 .AND. VCELLCTR(L) >= 0.0 )THEN
            UCELLCTRM = SUB(L  )*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(LN )*ABS(VCELLCTR(L))
            QSBDLDX(L  ,NX) = SUB(L  )*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(LN ,NX) = SVB(LN )*QSBDLDP(L)*VCELLCTR(L)
            IF( UCELLCTRM < 1.0E-9 )THEN
              QSBDLDY(LN ,NX) = SVB(LN )*SIGN(QSBDLDP(L),VCELLCTR(L))
            ENDIF
            IF( VCELLCTRM < 1.0E-9 )THEN
              QSBDLDX(L  ,NX) = SUB(L  )*SIGN(QSBDLDP(L),UCELLCTR(L))
            ENDIF
          ENDIF
          IF( UCELLCTR(L) < 0.0 .AND. VCELLCTR(L) < 0.0 )THEN
            UCELLCTRM = SUB(L)*ABS(UCELLCTR(L))
            VCELLCTRM = SVB(L)*ABS(VCELLCTR(L))
            QSBDLDX(L,NX) = SUB(L)*QSBDLDP(L)*UCELLCTR(L)
            QSBDLDY(L,NX) = SVB(L)*QSBDLDP(L)*VCELLCTR(L)
            IF( UCELLCTRM < 1.0E-9 )THEN
              QSBDLDY(L,NX) = SVB(L)*SIGN(QSBDLDP(L),VCELLCTR(L))
            ENDIF
            IF( VCELLCTRM < 1.0E-9 )THEN
              QSBDLDX(L,NX) = SUB(L)*SIGN(QSBDLDP(L),UCELLCTR(L))
            ENDIF
          ENDIF
        ENDIF
      ENDDO

    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  ELSEIF( ISBLFUC == 2 )THEN
    ! *** CALCULATE CELL FACE TRANSPORT RATES BY AVERAGING VECTOR COMPONENTS FROM CELL CENTERS TO FACES                                                                      
    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LS,LW)   
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        LS = LSC(L)
        LW = LWC(L)
        QSBDLDX(L,NX) = 0.5*SUB(L)*(QSBDLDP(L)*UCELLCTR(L) + QSBDLDP(LW)*UCELLCTR(LW))
        QSBDLDY(L,NX) = 0.5*SVB(L)*(QSBDLDP(L)*VCELLCTR(L) + QSBDLDP(LS)*UCELLCTR(LS))
      ENDDO

    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  ENDIF    ! *** END OF CELL FACE TRANSPORT RATES SECTION
    
    
  ! ********************************************************************************
  ! *** CONVERT TRANSPORT VECTORS TO FACE VECTORS  (G/S)                                                                        
  !$OMP DO PRIVATE(ND,LF,LL,LP,L)   
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    DO LP=LF,LL
      L=LWET(LP)  
      QSBDLDX(L,NX) = SUB(L)*DYU(L)*QSBDLDX(L,NX)
      QSBDLDY(L,NX) = SVB(L)*DXV(L)*QSBDLDY(L,NX)
    ENDDO

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO

  ! ********************************************************************************
  ! *** ELIMINATE BEDLOAD TRANSPORT UP ADVERSE SLOPES IN DIRECTION OF FLOW                                                
  IF( BLBSNT > 0.0 )THEN

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,SLOPE)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        IF( QSBDLDX(L,NX) > 0.0 )THEN
          SLOPE = (BELV(L)-BELV(LWC(L)))*DXIU(L)
          IF( SLOPE > BLBSNT ) QSBDLDX(L,NX) = 0.0
        ENDIF
        IF( QSBDLDX(L,NX) < 0.0 )THEN
          SLOPE = (BELV(LWC(L))-BELV(L))*DXIU(L)
          IF( SLOPE > BLBSNT ) QSBDLDX(L,NX) = 0.0
        ENDIF
        IF( QSBDLDY(L,NX) > 0.0 )THEN
          SLOPE = (BELV(L)-BELV(LSC(L)))*DYIV(L)
          IF( SLOPE > BLBSNT ) QSBDLDY(L,NX) = 0.0
        ENDIF
        IF( QSBDLDY(L,NX) < 0.0 )THEN
          SLOPE = (BELV(LSC(L))-BELV(L))*DYIV(L)
          IF( SLOPE > BLBSNT ) QSBDLDY(L,NX) = 0.0
        ENDIF
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  ENDIF

  ! ********************************************************************************
  ! *** LIMIT OUTGOING FLUXES IN EACH CELL                                                                                
  ! *** PMC - DELETED FOR 8.4.4.  SUFFICIENT SEDIMENT CHECK IN CALSND
  
  ! ********************************************************************************
  !$OMP DO PRIVATE(ND,LF,LL,LP,L,LE,LN)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    ! *** CALCULATE MASS PER UNIT AREA CHANGE IN BED CONCENTRATION DUE TO TO NET BED LOAD (G/M2/S)                                                                                                  
    DO LP=LF,LL
      L  = LWET(LP)  
      LE = LEC(L)
      LN = LNC(L)
      ! ***                              U COMPONENT                    V COMPONENT
      SNDFBL(L,NX) = DXYIP(L)*( QSBDLDX(LE,NX)-QSBDLDX(L,NX) + QSBDLDY(LN,NX)-QSBDLDY(L,NX) )
    ENDDO

  ENDDO   ! *** END OF DOMAIN LOOP
  !$OMP END DO
  !$OMP END PARALLEL
  
  ! ******************************************************************************
  ! *** COMPUTE CELL CENTERED BEDLOAD FLUXES FOR OUTFLOW OR RECIRCULATION BOUNDARY                                                                                               
  IF( NSBDLDBC > 0 )THEN
    DO NSB=1,NSBDLDBC
      LUTMP = LSBLBCU(NSB)
      LDTMP = LSBLBCD(NSB)
      
      QSBDLDOT(LUTMP,NX) = SNDFBL(LUTMP,NX)*DXYP(LUTMP)
      IF( LDTMP > 0 )THEN
        QSBDLDIN(LDTMP,NX) = QSBDLDIN(LDTMP,NX) + QSBDLDOT(LUTMP,NX)
        SNDFBL(LDTMP,NX) = SNDFBL(LDTMP,NX) + QSBDLDOT(LUTMP,NX)*DXYIP(LDTMP)
      ENDIF
      SNDFBL(LUTMP,NX) = 0.0
    ENDDO
    
  ENDIF

  ! ********************************************************************************
  RETURN

END

