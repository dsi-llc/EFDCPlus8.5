MODULE WQ_RPEM_MODULE

! ** EFDC_DSI  ROOTED PLANT AND EPIPHYTE MODEL (RPEM)
! ** MODULE: WQ_RPEM_MODULE

! CHANGE RECORD 
! DATE MODIFIED     BY               DESCRIPTION        
!-- ------------------------------------------------------------------
! 2012-02-25        Paul M. Craig    Fixed Temperature Dependency Table
! 2011-11-29        PAUL M. CRAIG &  OMP'd MODEL AND ADDED COMPUTATIONAL BYPASSES
!                   DANG CHUNG       RESTRUCTURED TO F90 AND USE OF MODULES 

USE GLOBAL

IMPLICIT NONE

REAL,ALLOCATABLE,DIMENSION(:) :: WQTDTEMP

CONTAINS

SUBROUTINE CAL_RPEM
  ! **********************************************************************C                                                
  ! **  SUBROUTINE WQRPEM SIMULATES ROOTED PLANTS (SHOOTS AND ROOTS),                                                     
  ! **  EPIPHYTES GROWING ON SHOOTS, AND SHOOT ORGANIC DETRITUS                                                           
 
  ! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION                                                                     
 
  ! **********************************************************************C                                                
  INTEGER :: L,IWQTRPEM,K,LL,LF,LP,ND,IOBC,LN                         
  REAL    :: RATION,TOP,TMPNIT,RATIOP,TMPPO4,WQAVGIO,TMP1,RATIOHP,HDRY2    
  REAL    :: RKESSAVG,ALPHATOP,ALPHABOT,TMPEXP,SOURSINK,FACIMP,DTDHWQ 
  REAL    :: TMPWAT,TMPBED,TOP1,TOP2,BOT1,BOT2,BOT3,TMPNH4S,TMPNO3S
  REAL    :: TMPNH4E,TMPNO3E,WQKESS                                                            
  REAL    :: RKESSTOP(LCM),RKESSBOT(LCM)
  REAL    :: WQBCV(NBCSOP,NWQV)
  
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:)   :: NLRPEM   ! *** NUMBER OF RPEM ACTIVE CELLS FOR EACH LAYER BY DOMAIN
  INTEGER,SAVE,ALLOCATABLE,DIMENSION(:,:) :: LLRPEM   ! *** L INDEX FOR THE RPEM CELLS
  INTEGER,SAVE :: ICOUNT
  
  ! ** INPUT AND INITIALIZATION                                                                                           
  IF( JSRPEM == 1 )THEN
    NSRPEMSPFR=0                                                                                                         
    NCRPEMRST=0                                                                                                          
    
    ALLOCATE(NLRPEM(NDM))
    ALLOCATE(LLRPEM(LCM,NDM))
    NLRPEM = 0
    LLRPEM = 0
    ICOUNT = 0
  ENDIF
  ICOUNT = ICOUNT+1
  
  ! *** RPEM ONLY INTERACTS WITH THE BOTTOM LAYER K=KSZ(L)
  
  ! *** SAVE VALUES AT OPEN BOUNDARIES
  IF( IRPEMWC == 0 )THEN
    DO IOBC=1,NBCSOP  
      L=LOBCS(IOBC)  
      WQBCV(IOBC,1:NWQV)  = WQV(L,KSZ(L),1:NWQV)
    ENDDO  
  ENDIF
  HDRY2 = 2.*HDRY
  
  !$OMP PARALLEL DEFAULT(SHARED)
  
  ! *** OBTAIN THE ACTIVE CELL LIST
  IF( (ISDRY > 0 .AND. LADRY > 0) .OR. JSRPEM == 1 )THEN
    !$OMP DO PRIVATE(ND,LF,LL,LN,LP,L)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)
      LN = 0
      DO LP=LF,LL
        L=LWET(LP)
        IF( LMASKRPEM(L) )THEN
          LN = LN+1
          LLRPEM(LN,ND) = L
          IF( ISICE > 2 )THEN
            IF( ICECELL(L) .AND. HP(L) < 3.*HDRY )THEN
              ! *** DEACTIVATE RPEM KINETICS IF THE CELL HAS ICE COVER AND VERY SHALLOW DEPTHS
              LN = LN-1
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      NLRPEM(ND) = LN
      IF( JSRPEM == 1 ) WRITE(6,*) 'RPEM DOMAIN LIST',N,ND,LN  
    ENDDO
    !$OMP END DO
  ENDIF
  
  !$OMP DO PRIVATE(ND,LP,L,IWQTRPEM,K) &
  !$OMP PRIVATE(RATION,TOP,TMPNIT,RATIOP,TMPPO4,WQAVGIO,TMP1)  &
  !$OMP PRIVATE(RKESSAVG,ALPHATOP,ALPHABOT,TMPEXP,SOURSINK,FACIMP,DTDHWQ) & 
  !$OMP PRIVATE(TMPWAT,TMPBED,TOP1,TOP2,BOT1,BOT2,BOT3,TMPNH4S,TMPNO3S)   &
  !$OMP PRIVATE(TMPNH4E,TMPNO3E,WQKESS,RKESSTOP,RKESSBOT,RATIOHP)
  DO ND=1,NDM  
    ! **********************************************************************C                                                
    ! ** SET TEMPERATURE DEPENDENCY FOR GROWTH AND RESPIRATION 
    ! ** F3(T), EQ.(13), LOOK-UP TABLE
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      IWQTRPEM = NINT((TWQ(L)-WQTDMIN)/WQTDINC)+1
      IF( IWQTRPEM < 1) IWQTRPEM=1
      IF( IWQTRPEM > NWQTD) IWQTRPEM=NWQTD
      XLIMTPRPS(L) = RPEMTPRPS(IWQTRPEM)  ! Shoot Production (Growth)                           
      XLIMTPRPE(L) = RPEMTPRPE(IWQTRPEM)  ! Epiphyte Production (Growth)
      XLIMTRRPS(L) = RPEMTRRPS(IWQTRPEM)  ! Shoot Respiration
      XLIMTRRPE(L) = RPEMTRRPE(IWQTRPEM)  ! Epiphyte Respiration
      XLIMTRRPR(L) = RPEMTRRPR(IWQTRPEM)  ! Root Respiration
    ENDDO                                                                                                                  
   
    ! **********************************************************************C                                                
    ! ** SET NUTRIENT LIMITATION FOR PLANT SHOOT GROWTH                                                                     
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      K = KSZ(L)
      RATION = RKHNRPS/RKHNRPR                                       ! Ratio of N half saturation constants for water/bed  
      TOP = WQV(L,K,14) + WQV(L,K,15) + RATION*(SM2NH4(L)+SM2NO3(L))                                                           
      TMPNIT = TOP/(TOP+RKHNRPS)
      RATIOP = RKHPRPS/RKHPRPR
      TOP = WQPO4D(L,K)+RATIOP*SM2PO4(L)                                                                                   
      TMPPO4 = TOP/(TOP+RKHPRPS)
      XLIMNRPS(L) = MIN(TMPNIT,TMPPO4)                                ! F1(N), EQ. (6), Minimum of N or P limit
    ENDDO                                                                                                                  
   
    ! **********************************************************************C                                                
    ! ** SET NUTRIENT LIMITATION FOR EPIPHYTE GROWTH                                                                        
    IF( IRPEME > 0 )THEN
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TOP = WQV(L,K,14)+WQV(L,K,15)                                                                                        
        TMPNIT = TOP/(TOP+RKHNRPE+1.E-18)                                                                                           
        TMPPO4=WQPO4D(L,K)/(WQPO4D(L,K)+RKHPRPE+1.E-18)
        XLIMNRPE(L) = MIN(TMPNIT,TMPPO4)                                ! *** Minimum of N or P limits.  EQ. (20)
      ENDDO                                                                                                                  
    ENDIF
     
    ! **********************************************************************C                                                
    ! ** SET LIGHT LIMITATIONS                                                                                              
   
    ! ** NOTE THAT WQI0,WQI1,AND WQI2 ARE PHOTOSYNTHETIC SSW RADIATION                                                      
    !    FROM CURRENT AND PREVIOUS TWO TIME INCREMENTS IN LANGLEY/DAY                                                       
    !    AND ARE PROVIDED BY THE WATER QUALITY MODEL                                                                        
    WQAVGIO = WQCIA*WQI0+WQCIB*WQI1+WQCIC*WQI2                          ! EQ. (12)
   
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      ! *** COMPUTE TOTAL EXTINCTION COEFFICIENT

      ! *** Bottom PLUS 1 Layer
      IF( KSZ(L) /= KC )THEN
        RKESSTOP(L) = RADKE(L,KSZ(L)+1)   ! EQ. (24)
      ELSE
        RKESSTOP(L)=0.                    ! EQ. (24)
      ENDIF
              
      ! *** Bottom Layer
      WQKESS = RADKE(L,KSZ(L))
      RKESSBOT(L) = WQKESS
    ENDDO                                                                                                                  
   
    ! *** LIGHT LIMITATION SECTION
    IF( IRPEME > 0 )THEN
      ! ** LIGHT LIMITATION FOR SHOOTS                                                                                        
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RKESSBOT(L) = RKESSBOT(L) + RKERPE*Z(L,KSZ(L))*HWQ(L)*WQRPE(L)/CCHLRPE   ! EQ. (10), INCLUDE EPIPHYTES SHADING
      ENDDO                                                                                                                  
    ENDIF
   
    IF( IWQSUN == 2 )THEN
      ! *** ASER TIME INTERVAL VARIABLE LIGHT
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSO(L)=RISSOM                                                ! Maximum value for optimum SSW for rooted plant growth (Langleys/day)
        RKESSAVG=0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        WQAVGIO = PARADJ*2.065*RADTOP(L,KC)                            ! RADTOP INCLUDES SHADING AND ICE COVER
        ALPHATOP=-(WQAVGIO/RISSO(L)) *EXP(-RKESSTOP(L)*TMP1*HWQ(L))    ! EQ. (9), ASSUME HRPS IS WITHIN THE BOTTOM LA
        ALPHABOT=-(WQAVGIO/RISSO(L)) *EXP(-RKESSBOT(L)*HWQ(L))         ! EQ. (8)
        TMPEXP=EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPS(L)=2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP         ! F2(I), EQ. (7)
        XLIMLRPS(L) = MIN(XLIMLRPS(L),1.0)
      ENDDO                                                                                                                  
    ELSE
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSO(L)=RISSOM                                                ! Maximum value for optimum SSW for rooted plant growth (Langleys/day)
        RKESSAVG=0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        ALPHATOP=-(WQAVGIO/RISSO(L)) *EXP(-RKESSTOP(L)*TMP1*HWQ(L))    ! EQ. (9), ASSUME HRPS IS WITHIN THE BOTTOM LA
        ALPHABOT=-(WQAVGIO/RISSO(L)) *EXP(-RKESSBOT(L)*HWQ(L))         ! EQ. (8)
        TMPEXP=EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPS(L)=2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP         ! F2(I), EQ. (7), BUT PRACTICALLY, IT DOES NOT
        XLIMLRPS(L) = MIN(XLIMLRPS(L),1.0)
      ENDDO                                                                                                                  
    ENDIF
    
    IF( IRPEME > 0 )THEN
      ! ** LIGHT LIMITATION FOR EPIPHYTES                                                                                     
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RISSOE(L)=RISSOEM                                                                                                  
        RKESSAVG=0.5*(RKESSTOP(L)+RKESSBOT(L))
        TMP1 = 1.-Z(L,KSZ(L))
        ALPHATOP=-(WQAVGIO/RISSOE(L))*EXP(-RKESSTOP(L)*TMP1*HWQ(L))   ! EQ. (21), ASSUME HRPS IS WITHIN THE BOTTOM L
        ALPHABOT=-(WQAVGIO/RISSOE(L))*EXP(-RKESSBOT(L)*HWQ(L))        ! EQ. (22)
        TMPEXP=EXP(ALPHABOT)-EXP(ALPHATOP)
        XLIMLRPE(L)=2.718/(RKESSAVG*Z(L,KSZ(L))*HWQ(L))*TMPEXP        ! EQ. (21)
        XLIMLRPE(L) = MIN(XLIMLRPE(L),1.0)
      ENDDO                                                                                                                  
    ENDIF
  
    IF( IJRPRS == 2 )THEN
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        TMP1 = 1.-Z(L,KSZ(L))
        RISS(L) = WQAVGIO*EXP(-RKESSTOP(L)*TMP1*HWQ(L))     ! EQ. (11), PART OF IT, USED IF IJRPRS=2 IN CAL_RPEM.INP C5
                                                            ! HOPT IS NOT USED HERE!, CAL_RPEM.INP C9
      ENDDO                                                                                                                  
    ENDIF
     
    ! **********************************************************************C                                                
    ! ** UPDATE GROWTH AND RESPIRATION RATES                                                                                
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      
      ! *** Limitors     TEM          N/P         Light         Self Shading
      PRPS(L) = PMRPS*XLIMTPRPS(L)*XLIMNRPS(L)*XLIMLRPS(L)*EXP(-RKSH*WQRPS(L))  ! EQ. (5)   PRPS - Nutrient, light and temperature limited growth rate (/day)
                                                                                !                  Include self-shading of shoots, JI, 7/4/04
      RRPS(L) = RMRPS*XLIMTRRPS(L)                                              ! EQ. (15)  RRPS - Temperature limited shoot respiration (/day)
    ENDDO                                                                                                                  
   
    ! *** Limit RPEM with depth
    IF( ISDRY > 0 )THEN
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RATIOHP = MIN((MAX(HP(L)-HDRY,0.))/HDRY2,1.0)
        PRPS(L) = PRPS(L)*RATIOHP
        RRPS(L) = RRPS(L)*RATIOHP
      ENDDO    
    ENDIF
    
    ! ** ROOT RESPIRATION                                                                                                   
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      RRPR(L)=RMRPR*XLIMTRRPR(L)                                    ! EQ. (18)
    ENDDO                                                                                                                  
   
    ! ** EPIPHYTE GROWTH AND RESPIRATION                                                                                    
    IF( IRPEME > 0 )THEN
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        PRPE(L) = PMRPE*XLIMTPRPE(L)*XLIMNRPE(L)*XLIMLRPE(L)            ! EQ. (19), F3T*F1N*F2I
        RRPE(L) = RMRPE*XLIMTRRPE(L)                                                                                         
      ENDDO                                                                                                                  
    ENDIF
     
    ! **********************************************************************C                                                
   
    ! ** UPDATE ROOT TO SHOOT FLUX                                                                                          
    IF( IJRPRS == 0 .AND. N < 5 )THEN                      ! *** IJRPRS = 0, CONSTANT TRANSPORT RATE
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RJRPRSC
      ENDDO                                                                                                                  
    ENDIF                                                                                                                  
   
    IF( IJRPRS == 1 )THEN
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RKRPORS*(ROSR*WQRPR(L)-WQRPS(L))                                                                         
      ENDDO                                                                                                                  
    ENDIF                                                                                                                  
   
    IF( IJRPRS == 2 )THEN
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        RJRPRS(L) = RKRPRS*RISS(L)/(RISS(L)+RISSS+1E-18)                                                                           
      ENDDO                                                                                                                  
    ENDIF                                                                                                                  
   
    ! **********************************************************************C                                                
    ! **  UPDATE SHOOT, ROOT, EPIPHYTE AND DETRITUS STATE VARIABLES                                                         
   
    ! ** UPDATE SHOOTS BIOMASS                                                                                                      
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      ! ***            Growth     Respiration   Other
      SOURSINK = (1.-FPRPR)*PRPS(L) - RRPS(L) - RLRPS                     ! EQ. (1)
      FACIMP = 1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPS(L) = FACIMP*(WQRPS(L) + DTWQ*RJRPRS(L))
      WQRPS(L) = MAX(WQRPS(L),0.2)   ! KEEP THE "SEED", SINCE IF RPS!=0, IT WILL NEVER GROW AGAIN, ACCORDING TO EQ.
    ENDDO                                                                                                                  
   
    ! ** UPDATE ROOTS BIOMASS                                                                                                      
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      SOURSINK=-RRPR(L)-RLRPR                                       ! EQ. (2)
      FACIMP=1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPR(L)=FACIMP*(WQRPR(L)+DTWQ*FPRPR*PRPS(L)*WQRPS(L)-DTWQ*RJRPRS(L))
    ENDDO                                                                                                                  

    !------------------------------                                                                                         
    !GO TO 904  ! SKIP EPHIPHYTES AND DETRITUS
    
    ! ** UPDATE EPIPHYTES                                                                                                   
    IF( IRPEME > 0 )THEN
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        SOURSINK=PRPE(L)-RRPE(L)-RLRPE                              ! EQ. (3)
        FACIMP=1./(1.-DTWQ*SOURSINK)                                                                                       
        WQRPE(L)=FACIMP*WQRPE(L)
      ENDDO                                                                                                                  
    ENDIF
     
    ! ** UPDATE SHOOT DETRITUS IN WATER COLUMN                                                                              
    SOURSINK=-RLRPD                                                    ! EQ. (4)  RLRPD-Loss rate for plant detritus at bottom of water column (/day)
    DO LP=1,NLRPEM(ND)
      L = LLRPEM(LP,ND)
      FACIMP=1./(1.-DTWQ*SOURSINK)                                                                                       
      WQRPD(L)=FACIMP*(WQRPD(L)+DTWQ*FRPSD*RLRPS*WQRPS(L))
    ENDDO                                                                                                                  
    904   CONTINUE                                                                                                          
   
    ! **********************************************************************C                                                
    ! ** CALCULATE WATER COLUMN SHOOT, EPIPHYTE AND DETRITUS                                                                
    ! ** RESPIRATION AND NON-RESPIRATION LOSSES TO WATER QUALITY ORGANIC                                                    
    ! ** MATTER STATE VARIABLES AND TO PHOSPHATE AND AMMONIA                                                                
   
    !     GO TO 901   ! SKIP COUPLING WITH NUTRIENTS AND DO IN WATER COLUMN                                    
    IF( IRPEMWC == 0 )THEN
      !WQTSNAME(1)  = 'CHC'  
      !WQTSNAME(2)  = 'CHD'  
      !WQTSNAME(3)  = 'CHG'  
      !WQTSNAME(4)  = 'ROC'  
      !WQTSNAME(5)  = 'LOC'  
      !WQTSNAME(6)  = 'DOC'  
      !WQTSNAME(7)  = 'ROP'  
      !WQTSNAME(8)  = 'LOP'  
      !WQTSNAME(9)  = 'DOP'  
      !WQTSNAME(10) = 'P4D'
      !WQTSNAME(11) = 'RON'  
      !WQTSNAME(12) = 'LON'  
      !WQTSNAME(13) = 'DON'  
      !WQTSNAME(14) = 'NHX'  
      !WQTSNAME(15) = 'NOX'  
      !WQTSNAME(16) = 'SUU'  
      !WQTSNAME(17) = 'SAA'  
      !WQTSNAME(18) = 'COD'  
      !WQTSNAME(19) = 'DOX'  
      !WQTSNAME(20) = 'TAM'  
      !WQTSNAME(21) = 'FCB'  

      ! ** UPDATE SHOOT, ROOT AND DETRITUS
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        ! ** LOSSES TO ORGANIC CARBON                                                                                           
        WQRPSR(L)=RRPS(L)*WQRPS(L)                                    ! EQ. (31A)   WQRPSR - Shoot biomass loss due to respiration
        WQRPSL(L)=(1.-FRPSD)*RLRPS*WQRPS(L)                           !             WQRPSL - Shoot biomass loss due to non-respiration processes
        WQRPDL(L)=RLRPD*WQRPD(L)                                      !             WQRPDL - Detritus biomass loss rate

        ! ** LOSSES TO ORGANIC PHOSPHOROUS                                                                                      
        WQRPSRP(L)=RPSPC*WQRPSR(L)                                    ! EQ. (37A)
        WQRPSLP(L)=RPSPC*WQRPSL(L)                                                                                         
        WQRPDLP(L)=RPSPC*WQRPDL(L)                                                                                         

        ! ** LOSSES TO ORGANIC NITROGEN                                                                                         
        WQRPSRN(L)=RPSNC*WQRPSR(L)                                    ! EQ. (41A)
        WQRPSLN(L)=RPSNC*WQRPSL(L)                                                                                         
        WQRPDLN(L)=RPSNC*WQRPDL(L)                                                                                         
      ENDDO
                                                                                      
      ! *** UPDATE EPIPHYTES
      IF( IRPEME > 0 )THEN
        DO LP=1,NLRPEM(ND)
          L = LLRPEM(LP,ND)
          
          ! ** LOSSES TO ORGANIC CARBON                                                                                           
          WQRPER(L)=RRPE(L)*WQRPE(L)                                                                                         
          WQRPEL(L)=RLRPE*WQRPE(L)                                                                                           

          ! ** LOSSES TO ORGANIC PHOSPHOROUS                                                                                      
          WQRPERP(L)=RPEPC*WQRPER(L)                                                                                         
          WQRPELP(L)=RPEPC*WQRPEL(L)                                                                                         

          ! ** LOSSES TO ORGANIC NITROGEN                                                                                         
          WQRPERN(L)=RPENC*WQRPER(L)                                                                                         
          WQRPELN(L)=RPENC*WQRPEL(L)                                                                                         
        ENDDO                                                                                                                  
      ENDIF
     
      ! ** UPDATE WATER COLUMN ORGANIC CARBON, PHOSPHOROUS AND NITROGEN                                                       
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ=DTWQ/(DZC(L,K)*HWQ(L))                                                                                        
     
        ! ** WATER COLUMN REFRACTORY PARTICULATE ORGANIC CARBON                                                                 
        WQV(L,K,4)=WQV(L,K,4)+DTDHWQ*(                                  &! ! RPOC
                    FCRRPS*WQRPSR(L)+FCRLRPS*WQRPSL(L)                  &! ! EQ. (31A) ??
                  +FCRRPE*WQRPER(L)+FCRLRPE*WQRPEL(L) &
                                    +FCRLRPD*WQRPDL(L) )
     
        ! ** WATER COLUMN LABILE PARTICULATE ORGANIC CARBON                                                                     
        WQV(L,K,5)=WQV(L,K,5)+DTDHWQ*(                                  &! ! LPOC
                    FCLRPS*WQRPSR(L)+FCLLRPS*WQRPSL(L)                  &! ! EQ. (32A)
                  +FCLRPE*WQRPER(L)+FCLLRPE*WQRPEL(L) &
                                    +FCLLRPD*WQRPDL(L) )
     
        ! ** WATER COLUMN DISSOLVED ORGANIC CARBON                                                                              
        WQV(L,K,6)=WQV(L,K,6)+DTDHWQ*(                                  &! ! DOC
                    FCDRPS*WQRPSR(L)+FCDLRPS*WQRPSL(L)                  &! ! EQ. (33A)
                  +FCDRPE*WQRPER(L)+FCDLRPE*WQRPEL(L) &
                                    +FCDLRPD*WQRPDL(L) )
     
        ! ** WATER COLUMN REFRACTORY PARTICULATE ORGANIC PHOSPHOROUS                                                            
        WQV(L,K,7)=WQV(L,K,7)+DTDHWQ*(                                  &! ! RPOP
                    FPRRPS*WQRPSRP(L)+FPRLRPS*WQRPSLP(L)                &! ! EQ. (34A)
                  +FPRRPE*WQRPERP(L)+FPRLRPE*WQRPELP(L) &
                                    +FPRLRPD*WQRPDLP(L) )
     
      ! ** WATER COLUMN LABILE PARTICULATE ORGANIC PHOSPHOROUS                                                                
        WQV(L,K,8)=WQV(L,K,8)+DTDHWQ*(                                 &! ! LPOP
                    FPLRPS*WQRPSRP(L)+FPLLRPS*WQRPSLP(L) &
                  +FPLRPE*WQRPERP(L)+FPLLRPE*WQRPELP(L) &
                                    +FPLLRPD*WQRPDLP(L) )
     
        ! ** WATER COLUMN DISSOLVED ORGANIC PHOSPHOROUS                                                                         
        WQV(L,K,9)=WQV(L,K,9)+DTDHWQ*(                                 &! ! DOP
                    FPDRPS*WQRPSRP(L)+FPDLRPS*WQRPSLP(L) &
                  +FPDRPE*WQRPERP(L)+FPDLRPE*WQRPELP(L) &
                                    +FPDLRPD*WQRPDLP(L) )
     
        ! ** WATER COLUMN TOTAL PHOSPHATE                                                                                       
        WQV(L,K,10)=WQV(L,K,10)+DTDHWQ*(                                &! ! PO4T, EQ. (37A), INCOMPLETE
                    FPIRPS*WQRPSRP(L)+FPILRPS*WQRPSLP(L) &
                  +FPIRPE*WQRPERP(L)+FPILRPE*WQRPELP(L) &
                                    +FPILRPD*WQRPDLP(L) )
     
        ! ** WATER COLUMN REFRACTORY PARTICULATE ORGANIC NITROGEN                                                               
        WQV(L,K,11)=WQV(L,K,11)+DTDHWQ*( &
                    FNRRPS*WQRPSRN(L)+FNRLRPS*WQRPSLN(L)                &! ! RPON
                  +FNRRPE*WQRPERN(L)+FNRLRPE*WQRPELN(L) &
                                    +FNRLRPD*WQRPDLN(L) )
     
        ! ** WATER COLUMN LABILE PARTICULATE ORGANIC NITROGEN                                                                   
        WQV(L,K,12)=WQV(L,K,12)+DTDHWQ*(                               &! ! LPON
                    FNLRPS*WQRPSRN(L)+FNLLRPS*WQRPSLN(L) &
                  +FNLRPE*WQRPERN(L)+FNLLRPE*WQRPELN(L) &
                                    +FNLLRPD*WQRPDLN(L) )
     
        ! ** WATER COLUMN DISSOLVED ORGANIC NITROGEN                                                                            
        WQV(L,K,13)=WQV(L,K,13)+DTDHWQ*(                               &! ! DON
                    FNDRPS*WQRPSRN(L)+FNDLRPS*WQRPSLN(L) &
                  +FNDRPE*WQRPERN(L)+FNDLRPE*WQRPELN(L) &
                                    +FNDLRPD*WQRPDLN(L) )
     
        ! ** WATER COLUMN AMMONIA NITROGEN                                                                                      
        WQV(L,K,14)=WQV(L,K,14)+DTDHWQ*(                               &! ! NH4, EQ. (42A), INCOMPLETE
                    FNIRPS*WQRPSRN(L)+FNILRPS*WQRPSLN(L) &
                  +FNIRPE*WQRPERN(L)+FNILRPE*WQRPELN(L) &
                                    +FNILRPD*WQRPDLN(L) )
     
      ENDDO                                                                                                                  
     
      ! **********************************************************************C                                                
      ! ** CALCULATE WATER COLUMN SHOOT AND EPIPHYTE UPTAKE OF                                                                
      ! ** DISSOLVED PHOSPHATE, AMMONIA, AND NO3                                                                              
     
      ! ** CALCULATE THE FRACTION OF DISSOLVED PHOSPHATE UPTAKE FROM WATER                                                    
      ! ** COLUMN (FRPSPW) BY PLANT SHOOTS                                                                                    
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TMPWAT=RKHPRPR*WQV(L,K,10)                                    ! ??, BUT (38) MEANS PO4DW, DISSOLVED ONLY, NO
        TMPBED=RKHPRPS*SM2PO4(L)                                                                                           
        FRPSPW(L)=TMPWAT/(TMPWAT+TMPBED+1E-18)                        ! EQ. (38)
      ENDDO                                                                                                                  
   
      ! ** UPDATE TOTAL PHOSPHATE IN WATER COLUMN                                                                             
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ=DTWQ/(DZC(L,K)*HWQ(L))                                                                                           
        WQV(L,K,10) = WQV(L,K,10) - DTDHWQ*( RPEPC*PRPE(L)*WQRPE(L) + RPSPC*FRPSPW(L)*PRPS(L)*WQRPS(L) )  ! PO4T, EQ. (37A), COMPLETED
        WQV(L,K,10) = MAX(WQV(L,K,10),0.)
      ENDDO                                                                                                                  
     
      ! ** CALCULATE THE FRACTION OF AMMONIA AND NO3 UPTAKE FROM WATER                                                        
      ! ** COLUMN (FRPSNW) BY PLANT SHOOTS                                                                                    
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TMPWAT=RKHNRPR*(WQV(L,K,14)+WQV(L,K,15))
        TMPBED=RKHNRPS*(SM2NH4(L)+SM2NO3(L))                                                                               
        FRPSNW(L)=TMPWAT/(TMPWAT+TMPBED+1E-18)                              ! EQ. (45)
      ENDDO                                                                                                                  
     
      ! ** CALCULATE THE AMMONIA PREFERENCE FOR SHOOTS                                                                        
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        TOP1=(WQV(L,K,14)+SM2NH4(L))*(WQV(L,K,15)+SM2NO3(L))                ! WHY USE WATER COLUMN AND BED, BOTH OF THEM?
        TOP2=RKHNPRPS*(WQV(L,K,14)+SM2NH4(L))                                                                              
        BOT1=RKHNPRPS+(WQV(L,K,14)+SM2NH4(L))                                                                              
        BOT2=RKHNPRPS+(WQV(L,K,15)+SM2NO3(L))                                                                              
        BOT3=WQV(L,K,14)+SM2NH4(L)+(WQV(L,K,15)+SM2NO3(L))     
        BOT1=MAX(BOT1,1E-12)                                                            
        BOT2=MAX(BOT2,1E-12)                                                            
        BOT3=MAX(BOT3,1E-12)                                                            
        PNRPS(L)=TOP1/(BOT1*BOT2)+TOP2/(BOT2*BOT3)                           ! EQ. (44A)
      ENDDO                                                                                                                  
   
      ! ** CALCULATE THE AMMONIA PREFERENCE FOR EPIPHYTES                                                                     
      IF( IRPEME > 0 )THEN
        DO LP=1,NLRPEM(ND)
          L = LLRPEM(LP,ND)
          K = KSZ(L)
          ! ***    NH3        NO3                                                                                                 
          TOP1=WQV(L,K,14)*WQV(L,K,15)                                                                                       
          TOP2=RKHNPRPE*WQV(L,K,14)                                                                                          
          BOT1=RKHNPRPE+WQV(L,K,14)                                                                                          
          BOT2=RKHNPRPE+WQV(L,K,15)                                                                                          
          BOT3=WQV(L,K,14)+WQV(L,K,15)     
          BOT1=MAX(BOT1,1E-12)                                                               
          BOT2=MAX(BOT2,1E-12)                                                               
          BOT3=MAX(BOT3,1E-12)                                                               
          PNRPE(L)=TOP1/(BOT1*BOT2)+TOP2/(BOT2*BOT3)                    ! EQ. (44B)                                                                    
        ENDDO                                                                                                                  
      ENDIF
     
      ! ** UPDATE WATER COLUMN AMMONIA AND NO3 NITROGEN                                                                       
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ=DTWQ/(DZC(L,K)*HWQ(L))                                                                                           
        TMPNH4S=PNRPS(L)                                                                                                   
        TMPNO3S=1.-PNRPS(L)                                                                                                
        TMPNH4E=PNRPE(L)                                                                                                   
        TMPNO3E=1.-PNRPE(L)                                                                                                
        WQV(L,K,14) = WQV(L,K,14) - DTDHWQ*(RPENC*TMPNH4E*PRPE(L)*WQRPE(L) + RPSNC*TMPNH4S*FRPSNW(L)*PRPS(L)*WQRPS(L))  ! NH4, EQ. (42A), COMPLETED
        WQV(L,K,15) = WQV(L,K,15) - DTDHWQ*(RPENC*TMPNO3E*PRPE(L)*WQRPE(L) + RPSNC*TMPNO3S*FRPSNW(L)*PRPS(L)*WQRPS(L))  ! NO3, EQ. (42B), COMPLETED
        WQV(L,K,14) = MAX(WQV(L,K,14),0.)
        WQV(L,K,15) = MAX(WQV(L,K,15),0.)
      ENDDO                                                                                                                  
     
      ! **********************************************************************C                                                
      ! ** CALCULATE WATER COLUMN SOURCES OF DISSOLVED OXYGEN DUE                                                             
      ! ** TO SHOOT AND EPIPHYTE GROWTH                                                                                       
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ = DTWQ/(DZC(L,K)*HWQ(L))       
        ! *** RPSOC - Plant shoot oxygen to carbon ratio (gO2/gC)
        ! *** PRPS  - Net shoot production rate (/day)
        ! *** WQRPS - Shoot mass (gC/m2)
        ! *** RPEOC - Epiphyte oxygen to carbon ratio (gO2/gC)
        ! *** PRPS  - Net shoot production rate (/day)
        ! *** WQRPS - Shoot mass (gC/m2)
        ! *** WQV(L,K,19) - Dissolved  Oxygen (g/m3)
        WQV(L,K,19) = WQV(L,K,19) + DTDHWQ*(RPSOC*WQRPS(L)*(PRPS(L)-RRPS(L)) + RPEOC*WQRPE(L)*(PRPE(L)-RRPE(L)))  ! DO, EQ. (33)
        WQV(L,K,19) = MAX(WQV(L,K,19),0.0)
      ENDDO                                                                                                                  
    ENDIF    ! *** END OF BLOCK FOR SKIPPING WATER COLUMN LINKAGE
   
    !901   CONTINUE                                                                                                          
   
    ! **********************************************************************C                                                
    ! ** CALCULATE SEDIMENT BED ROOT RESPIRATION AND NON-RESPIRATION                                                        
    ! ** LOSSES TO SEDIMENT DIAGENESIS ORGANIC MATTER STATE VARIABLES                                                       
    ! ** AND TO PHOSPHATE AND AMMONIA                                                                                       
   
    !     GO TO 902   ! SKIP COUPLING WITH THE SEDIMENT DIAGENESIS MODEL.  FULL DIAGENESIS MUST BE ON TO COUPLE
    
    IF( IRPEMSD == 0 .AND. IWQBEN == 1 )THEN   
      ! ** LOSSES TO ORGANIC CARBON                                                                                           
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRR(L)=RRPR(L)*WQRPR(L)                                    ! FROM EQ. (32B)
        WQRPRL(L)=RLRPR*WQRPR(L)                                                                                           
      ENDDO                                                                                                                  
     
      ! ** LOSSES TO ORGANIC PHOSPHOROUS                                                                                      
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRRP(L)=RPRPC*WQRPRR(L)                                    !  ! FROM EQ. (34B)
        WQRPRLP(L)=RPRPC*WQRPRL(L)                                                                                         
      ENDDO                                                                                                                  
     
      ! ** LOSSES TO ORGANIC NITROGEN                                                                                         
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        WQRPRRN(L)=RPRNC*WQRPRR(L)                                    ! FROM EQ. (39B)
        WQRPRLN(L)=RPRNC*WQRPRL(L)                                                                                         
      ENDDO                                                                                                                  
     
      ! ** UPDATE SEDIMENT BED ORGANIC CARBON, PHOSPHOROUS AND NITROGEN                                                       
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        DTDHWQ=DTWQ/SMHSED(1)
     
        ! ** SEDIMENT BED REFRACTORY PARTICULATE ORGANIC CARBON FRPRRPG1  FRPRRPG2   FRPRRPG3                                    
     
        SMPOC(L,1)=SMPOC(L,1)+DTDHWQ* FRPRRPG1*(FCRRPR*WQRPRR(L)+FCRLRPR*WQRPRL(L) )  ! EQ. (32B), REFRACTORY
        SMPOC(L,2)=SMPOC(L,2)+DTDHWQ* FRPRRPG2*(FCRRPR*WQRPRR(L)+FCRLRPR*WQRPRL(L) )
        SMPOC(L,3)=SMPOC(L,3)+DTDHWQ* FRPRRPG3*(FCRRPR*WQRPRR(L)+FCRLRPR*WQRPRL(L) )
     
        ! ** SEDIMENT BED LABILE PARTICULATE ORGANIC CARBON                                                                     
        SMPOC(L,1)=SMPOC(L,1)+DTDHWQ*( FCLRPR*WQRPRR(L)+FCLLRPR*WQRPRL(L) )           ! EQ. (32C), LABILE
     
        ! ** SEDIMENT BED DISSOLVED ORGANIC CARBON                                                                              
        SMPOC(L,1)=SMPOC(L,1)+DTDHWQ*(FCDRPR*WQRPRR(L)+FCDLRPR*WQRPRL(L) )            ! EQ. (33B), DISSOLVED IN WATER COLUMN IS AS
     
        ! ** SEDIMENT BED REFRACTORY PARTICULATE ORGANIC PHOSPHOROUS                                                            
        SMPOP(L,1)=SMPOP(L,1)+DTDHWQ* FRPRRPG1*(FPRRPR*WQRPRRP(L)+FPRLRPR*WQRPRLP(L) )
        SMPOP(L,2)=SMPOP(L,2)+DTDHWQ* FRPRRPG2*(FPRRPR*WQRPRRP(L)+FPRLRPR*WQRPRLP(L) )
        SMPOP(L,3)=SMPOP(L,3)+DTDHWQ* FRPRRPG3*(FPRRPR*WQRPRRP(L)+FPRLRPR*WQRPRLP(L) )
     
        ! ** SEDIMENT BED LABILE PARTICULATE ORGANIC PHOSPHOROUS                                                                
        SMPOP(L,1)=SMPOP(L,1)+DTDHWQ*(FPLRPR*WQRPRRP(L)+FPLLRPR*WQRPRLP(L) )
     
        ! ** SEDIMENT BED DISSOLVED ORGANIC PHOSPHOROUS                                                                         
        SMPOP(L,1)=SMPOP(L,1)+DTDHWQ*(FPDRPR*WQRPRRP(L)+FPDLRPR*WQRPRLP(L) )
     
        ! ** SEDIMENT BED TOTAL PHOSPHATE                                                                                       
        SM2PO4(L)=SM2PO4(L)+DTDHWQ*(FPIRPR*WQRPRRP(L)+FPILRPR*WQRPRLP(L) )             ! EQ. (37B), INCOMPLETE
     
        ! ** SEDIMENT BED REFRACTORY PARTICULATE ORGANIC NITROGEN                                                               
        SMPON(L,1)=SMPON(L,1)+DTDHWQ* FRPRRPG1*(FNRRPR*WQRPRRN(L)+FNRLRPR*WQRPRLN(L) )
        SMPON(L,2)=SMPON(L,2)+DTDHWQ* FRPRRPG2*(FNRRPR*WQRPRRN(L)+FNRLRPR*WQRPRLN(L) )
        SMPON(L,3)=SMPON(L,3)+DTDHWQ* FRPRRPG3*(FNRRPR*WQRPRRN(L)+FNRLRPR*WQRPRLN(L) )
     
        ! ** SEDIMENT BED LABILE PARTICULATE ORGANIC NITROGEN                                                                   
        SMPON(L,1)=SMPON(L,1)+DTDHWQ*(FNLRPR*WQRPRRN(L)+FNLLRPR*WQRPRLN(L) )
     
        ! ** SEDIMENT BED DISSOLVED ORGANIC NITROGEN                                                                            
        SMPON(L,1)=SMPON(L,1)+DTDHWQ*(FNDRPR*WQRPRRN(L)+FNDLRPR*WQRPRLN(L) )
     
        ! ** SEDIMENT BED AMMONIA NITROGEN                                                                                      
        SM2NH4(L)=SM2NH4(L)+DTDHWQ*(FNIRPR*WQRPRRN(L)+FNILRPR*WQRPRLN(L) )              ! EQ. (42B), INCOMPLETE
      ENDDO                                                                                                                  
     
      ! **********************************************************************C                                                
      ! ** CALCULATE SEDIMENT BED ROOT UPTAKE OF DISSOLVED PHOSPHATE,                                                         
      ! ** AMMONIA AND NO3                                                                                                    
     
      ! ** UPDATE TOTAL PHOSPHATE IN SEDIMENT BED                                                                             
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ=DTWQ/(DZC(L,K)*HWQ(L))                                       ! ACCORDING TO (37B), IT'S NOT SMHSED(1)!
        SM2PO4(L)=SM2PO4(L)-DTDHWQ*RPRPC*(1.-FRPSPW(L))*PRPS(L)*WQRPS(L)    ! EQ. (37B), COMPLETED
      ENDDO                                                                                                                  
     
      ! ** UPDATE SEDIMENT BED AMMONIA AND NO3 NITROGEN                                                                       
      DO LP=1,NLRPEM(ND)
        L = LLRPEM(LP,ND)
        K = KSZ(L)
        DTDHWQ=DTWQ/(DZC(L,K)*HWQ(L))
        TMPNH4S=PNRPS(L)
        TMPNO3S=1.-PNRPS(L)                                                                                                
        TMPBED=1.-FRPSNW(L)                                                                                                
        SM2NH4(L)=SM2NH4(L)-DTDHWQ*RPSNC*TMPNH4S*TMPBED*PRPS(L)*WQRPS(L)   ! EQ. (42B), COMPLETED
        SM2NO3(L)=SM2NO3(L)-DTDHWQ*RPSNC*TMPNO3S*TMPBED*PRPS(L)*WQRPS(L)   ! EQ. (43B)
      ENDDO                                                                                                                  
    ENDIF   ! *** END OF SEDIMENT BED LINKAGE BYPASS
  
    902   CONTINUE                                                                                                          

  ENDDO    ! *** END OF DOMAIN LOOP 
  !$OMP END DO
   
  !$OMP END PARALLEL

  ! **********************************************************************C 
  ! *** RESTORE OPEN BOUNDARY CONCENTRATIONS                                               
  IF( IRPEMWC == 0 )THEN
    DO IOBC=1,NBCSOP  
      L=LOBCS(IOBC)
      K = KSZ(L)
      WQV(L,K,1:NWQV) = WQBCV(IOBC,1:NWQV)
    ENDDO  
  ENDIF
  
   
  ! **********************************************************************C                                                
  JSRPEM=0                                                                                                             

  RETURN

END SUBROUTINE                                                                                                                    
                                                                                                                      
!**********************************************************************C                                                

SUBROUTINE RPEMINP
  
  !**********************************************************************C                                                
  ! **  SUBROUTINE RPEMINPUT READS ALL INPUT FOR ROOTED PLANTS AND                                                        
  ! **  EPIPHYTES GROWING ON ROOTED PLANTS AND INITIALIZES VARIABLES                                                      
  ! **  INCLUDING TEMPERATURE GROWTH RELATIONSHIPS                                                                        
  ! **
  ! **  THIS SUBROUTINE IS PART OF  EFDC-FULL VERSION                                                                     
  ! **
  ! **  LAST MODIFIED BY JOHN HAMRICK ON 23 JUNE 2004                                                                     
  ! **
  !----------------------------------------------------------------------C                                                
  
  ! CHANGE RECORD                                                                                                         
  ! DATE MODIFIED     BY                 DATE APPROVED    BY                                                              
  
  !**********************************************************************C                                                
  
  INTEGER :: NS,IS,I,J,L,LDATA,LL,NT 
  REAL    :: TP2RPE,RPST,RPRT,RPET,RPDT,TIMETMP,WTEMP                                                                                                    
  REAL    :: SMNH4,SMNO3,SMPO4
  !**********************************************************************C                                                
  
  ! **  READ INPUT PARAMETER FILE                                                                                         
  WRITE(*,'(A)')' WQ: WQRPEM.INP - RPEM CONTROL FILE'
  OPEN(1,FILE='wqrpem.inp',STATUS='UNKNOWN')
  
  ! C0 title                                                                                                                
  CALL SEEK('C0')
  READ(1,1)  ! *** Skip Title
  
  ! C1 options and initial conditions                                                                                       
  CALL SEEK('C1')
  READ(1,*)INITRPEM,IRPEMWC,IRPEMSD,HRPEMIC,RPSO,RPRO,RPEO,RPDO,SMNH4,SMNO3,SMPO4                                                          

  ! C2 stoichemetric ratios                                                                                                 
  CALL SEEK('C2')
  READ(1,*)RPSOC,RPEOC,RPSNC,RPRNC,RPENC,RPSPC,RPRPC,RPEPC                                                               
  
  ! C3 plant shoot growth respiration and loss parameters                                                                   
  CALL SEEK('C3')
  READ(1,*)PMRPS,FPRPR,RMRPS,rLRPS,FRPSD,rKHNPRPS,rKSH,rKHI     ! added self-shading of shoots, Ji, 7/4/04
  ! *** Unused  rKHI
  
  ! C4 plant root respiration and loss parameters and detritus loss parameter                                               
  CALL SEEK('C4')
  READ(1,*)RMRPR,rLRPR,rLRPD
  
  ! C5 root to shoot carbon transport options and parameters                                                                
  CALL SEEK('C5')
  READ(1,*)iJRPRS,rJRPRSC,rKRPORS,ROSR,rKRPRS,rISSS                                                                      
  
  ! C6 epiphyte growth respiration and loss parameters                                                                      
  CALL SEEK('C6')
  READ(1,*)IRPEME,PMRPE,RMRPE,rLRPE,rKHNPRPE                                                                               
  
  ! C7 nutrient limitation effects on shoot and epiphyte growth                                                             
  CALL SEEK('C7')
  READ(1,*)rKHNRPS,rKHNRPR,rKHNRPE,rKHPRPS,rKHPRPR,rKHPRPE                                                               
  
  ! C8 temperature effects on shoot and epiphyte growth                                                                     
  CALL SEEK('C8')
  READ(1,*)TP1RPS,TP2RPS,rKTP1RPS,rKTP2RPS,TP1RPE,TP2RPE,rKTP1RPE,rKTP2RPE
  
  ! C9 shoot and epiphyte light limitation and extinction parameters                                                        
  CALL SEEK('C9')
  READ(1,*)HRPS,HOPT,rKeRPE,CChlRPE,rISSOM,rISSOEM                                                                       
  
  ! C10 salinity toxicity parameters for fresh water plants and epiphytes (Not Used)                                                  
  CALL SEEK('C10')
  READ(1,*)iSTOXRPE,STOXS,STOXE                                                                                          
  
  ! C11 temperature effects on shoot and root respiration                                                                   
  CALL SEEK('C11')
  READ(1,*)TR1RPS,TR2RPS,rKTR1RPS,rKTR2RPS,TR1RPR,TR2RPR,rKTR1RPR,rKTR2RPR
  
  ! C12 temperature effects on epiphyte respiration                                                                         
  CALL SEEK('C12')
  READ(1,*)TR1RPE, TR2RPE, rKTR1RPE, rKTR2RPE                                                                            
  
  ! C13 effects of shoot respiration and loss on organic carbon                                                             
  CALL SEEK('C13')
  READ(1,*)FCRrps,FCLrps,FCDrps,FCRLrps,FCLLrps,FCDLrps                                                                  
  
  ! C14 effects of root respiration and loss on organic carbon                                                              
  CALL SEEK('C14')
  READ(1,*)FCRrpr,FCLrpr,FCDrpr,FCRLrpr,FCLLrpr,FCDLrpr                                                                  
  
  ! C15 effects of epiphyte respiration and loss on organic carbon                                                          
  CALL SEEK('C15')
  READ(1,*)FCRrpe,FCLrpe,FCDrpe,FCRLrpe,FCLLrpe,FCDLrpe                                                                  
  
  ! C16 effects of detritus loss on organic carbon                                                                          
  CALL SEEK('C16')
  READ(1,*)FCRLrpd,FCLLrpd,FCDLrpd                                                                                       
  
  ! C17 effects of shoot respiration and loss on organic and inorganic phosphorous                                          
  CALL SEEK('C17')
  READ(1,*)FPRrps,FPLrps,FPDrps,FPIrps,FPRLrps,FPLLrps,FPDLrps,FPILrps
  
  ! C18 effects of root respiration and loss on organic and inorganic phosphorous                                           
  CALL SEEK('C18')
  READ(1,*)FPRrpr,FPLrpr,FPDrpr,FPIrpr,FPRLrpr,FPLLrpr,FPDLrpr,FPILrpr
  
  ! C19 effects of epiphyte respiration and loss on organic and inorganic phosphorous                                       
  CALL SEEK('C19')
  READ(1,*)FPRrpe,FPLrpe,FPDrpe,FPIrpe,FPRLrpe,FPLLrpe,FPDLrpe,FPILrpe
  
  ! C20 effects of detritus loss on organic and inorganic phosphorous                                                       
  CALL SEEK('C20')
  READ(1,*)FPRLrpd,FPLLrpd,FPDLrpd,FPILrpd                                                                               
  
  ! C21 effects of shoot respiration and loss on organic and inorganic nitrogen                                             
  CALL SEEK('C21')
  READ(1,*)FNRrps,FNLrps,FNDrps,FNIrps,FNRLrps,FNLLrps,FNDLrps,FNILrps
  
  ! C22 effects of root respiration and loss on organic and inorganic nitrogen                                              
  
  CALL SEEK('C22')
  READ(1,*)FNRrpr,FNLrpr,FNDrpr,FNIrpr,FNRLrpr,FNLLrpr,FNDLrpr,FNILrpr
  
  ! C23 effects of epiphyte respiration and loss on organic and inorganic nitrogen                                          
  CALL SEEK('C23')
  READ(1,*)FNRrpe,FNLrpe,FNDrpe,FNIrpe,FNRLrpe,FNLLrpe,FNDLrpe,FNILrpe
  
  ! C24 effects of detritus loss on organic and inorganic nitrogen                                                          
  CALL SEEK('C24')
  READ(1,*)FNRLrpd,FNLLrpd,FNDLrpd,FNILrpd                                                                               
  
  ! C25 diagensis transfer of refractory, liabile and dissolve organic matter                                               
  CALL SEEK('C25')
  READ(1,*)FRPRRPG1,FRPRRPG2,FRPRRPG3                                                                                    
  
  ! C26 output options                                                                                                      
  CALL SEEK('C26')
  READ(1,*)ISRPEMSPAC,ISRPEMSPFR,ISRPEMTIME,ISRPEMTIFR,ISRPEMTILC,NRPEMEE                                                        
  
  ! C27 time series output locations                                                                                        
  CALL SEEK('C27')
  IF( ISRPEMTILC > 0 )THEN                                                                                                
    DO IS=1,ISRPEMTILC                                                                                                   
      READ(1,*)I,J                                                                                                       
      IRPEMTS(IS)=I                                                                                                      
      JRPEMTS(IS)=J                                                                                                      
      LRPEMTS(IS)=LIJ(I,J)                                                                                               
    ENDDO                                                                                                                
  ENDIF                                                                                                                  
  
  CLOSE(1)
  1 FORMAT(1X)                                                                                                        
  
  !**********************************************************************C                                                
  ! ** SET SPATIALLY CONSTANT INITIAL CONDITIONS                                                                          
  DO L=1,LC
    WQRPS(L)=0.0
    WQRPR(L)=0.0                                                                                                         
    WQRPE(L)=0.0                                                                                                         
    WQRPD(L)=0.0                                                                                                         
    LMASKRPEM(L)=.FALSE.                                                                                                 
  ENDDO                                                                                                                  
  
  NRPEM = 0
  IF( INITRPEM == 0 )THEN
    DO L=1,LC
      IF( HP(L) <= HRPEMIC )THEN  ! uniform values at depth <= HRPEMIC                                                    
        WQRPS(L)=RPSO   ! *** Initial condition for shoot mass    ( gC/m^2 )
        WQRPR(L)=RPRO   ! *** Initial condition for root mass     ( gC/m^2 )                                                 
        WQRPE(L)=RPEO   ! *** Initial condition for epiphyte mass ( gC/m^2 )
        WQRPD(L)=RPDO   ! *** Initial condition for detritus mass ( gC/m^2 )
        LMASKRPEM(L)=.TRUE.                                                                                              
      ENDIF                                                                                                              
    ENDDO
    NRPEM = LA-1                                                                                                                
  ENDIF
  
  !**********************************************************************C                                                
  ! ** SET SPATIALLY VARIABLE INITIAL CONDITIONS                                                                          
  IF( INITRPEM == 1 )THEN
    WRITE(*,'(A)')' WQ: WQRPEMSIC.INP'
    OPEN(1,FILE='wqrpemsic.inp')
  
    DO NS=1,6
      READ(1,1)                                                                                                          
    ENDDO
  
    READ(1,*)LDATA
    DO LL=1,LDATA                                                                                                        
      READ(1,*)I,J,RPST,RPRT,RPET,RPDT                                                                                   
      L=LIJ(I,J)                                                                                                         
      WQRPS(L)=RPST
      WQRPR(L)=RPRT                                                                                                      
      WQRPE(L)=RPET                                                                                                      
      WQRPD(L)=RPDT                                                                                                      
      IF( HP(L) <= 10. ) LMASKRPEM(L)=.TRUE.   ! *** Only include the cell if not too deep
      NRPEM = NRPEM+1 
    ENDDO
  
    CLOSE(1)
  ENDIF                                                                                                                  
  
  !**********************************************************************C                                                
  ! ** READ RESTART CONDITIONS                                                                                            
  IF( INITRPEM == 2 )THEN
    WRITE(*,'(A)')' WQ: WQRPEMRST.INP'
    OPEN(1,FILE='wqrpemrst.inp')
  
    READ(1,*)TIMETMP

    100   CONTINUE                                                                                                        
    READ(1,*,END=200)L,RPST,RPRT,RPET,RPDT                                                                             
    WQRPS(L)=RPST
    WQRPR(L)=RPRT                                                                                                      
    WQRPE(L)=RPET                                                                                                      
    WQRPD(L)=RPDT                                                                                                      
    LMASKRPEM(L)=.TRUE.                                                                                                
    NRPEM = NRPEM+1 
    GOTO 100                                                                                                           
  200   CONTINUE                                                                                                        
  
    CLOSE(1)
  ENDIF                                                                                                                  
  
  !**************************************************************************************                                               
  ! ** GENERATE TEMPERATURE DEPENDENCY FOR GROWTH AND RESPIRATION OVER WQTDMIN TO WQTDMAX                                                        
  DO NT=1,NWQTD  
    WTEMP=WQTDTEMP(NT)
  
    ! *** Shoot Production (Growth)
    RPEMTPrps(NT)=1.
    IF( WTEMP < TP1RPS )THEN
      RPEMTPrps(NT)=EXP(-rKTP1RPS*(WTEMP-TP1RPS)*(WTEMP-TP1RPS) )
    ENDIF
    IF( WTEMP > TP2RPS )THEN
      RPEMTPrps(NT)=EXP(-rKTP2RPS*(WTEMP-TP2RPS)*(WTEMP-TP2RPS) )
    ENDIF
  
    ! *** Epiphyte Production (Growth)
    RPEMTPrpe(NT)=1.
    IF( WTEMP < TP1RPE )THEN
      RPEMTPrpe(NT)=EXP(-rKTP1RPE*(WTEMP-TP1RPE)*(WTEMP-TP1RPE) )
    ENDIF
    IF( WTEMP > TP2RPE )THEN
      RPEMTPrpe(NT)=EXP(-rKTP2RPE*(WTEMP-TP2RPE)*(WTEMP-TP2RPE) )
    ENDIF
  
    ! *** Shoot Respiration
    RPEMTRrps(NT)=1.
    IF( WTEMP < TR1RPS )THEN
      RPEMTRrps(NT)=EXP(-rKTR1RPS*(WTEMP-TR1RPS)*(WTEMP-TR1RPS) )
    ENDIF
    IF( WTEMP > TR2RPS )THEN
      RPEMTRrps(NT)=EXP(-rKTR2RPS*(WTEMP-TR2RPS)*(WTEMP-TR2RPS) )
    ENDIF
  
    ! *** Epiphyte Respiration
    RPEMTRrpe(NT)=1.
    IF( WTEMP < TR1RPE )THEN
      RPEMTRrpe(NT)=EXP(-rKTR1RPE*(WTEMP-TR1RPE)*(WTEMP-TR1RPE) )
    ENDIF
    IF( WTEMP > TR2RPE )THEN
      RPEMTRrpe(NT)=EXP(-rKTR2RPE*(WTEMP-TR2RPE)*(WTEMP-TR2RPE) )
    ENDIF
  
    ! *** Root Respiration
    RPEMTRrpr(NT)=1.
    IF( WTEMP < TR1RPR )THEN
      RPEMTRrpr(NT)=EXP(-rKTR1RPR*(WTEMP-TR1RPR)*(WTEMP-TR1RPR) )
    ENDIF
    IF( WTEMP > TR2RPR )THEN
      RPEMTRrpr(NT)=EXP(-rKTR2RPR*(WTEMP-TR2RPR)*(WTEMP-TR2RPR) )
    ENDIF
  
  ENDDO                                                                                                                  

  ! *** ASSIGN BED POREWATER CONCENTRATIONS IF SEDIMENT DIAGENESIS NOT SIMULATED
  IF( IWQBEN /= 1 )THEN
    DO L=2,LA
      SM2NH4(L) = SMNH4
      SM2NO3(L) = SMNO3
      SM2PO4(L) = SMPO4
    ENDDO
  ENDIF
  
  !**********************************************************************C                                                
  RETURN

END SUBROUTINE          

SUBROUTINE INIT_RPEMVARS

  ! *** ALLOCATE
  ALLOCATE(WQRPS(LCM)       , WQRPR(LCM)       , WQRPE(LCM)       , WQRPD(LCM)    ,              &
           PRPS(LCM)        , RRPS(LCM)        , RRPR(LCM)        , PRPE(LCM)     , RRPE(LCM)  , &
           WQRPSR(LCM)      , WQRPSL(LCM)      , WQRPER(LCM)      , WQRPEL(LCM)   , WQRPDL(LCM), &
           WQRPSRP(LCM)     , WQRPSLP(LCM)     , WQRPERP(LCM)     , WQRPELP(LCM)  ,              &
           WQRPDLP(LCM)     , WQRPSRN(LCM)     , WQRPSLN(LCM)     , WQRPERN(LCM)  ,              &
           WQRPELN(LCM)     , WQRPDLN(LCM)     , WQRPRR(LCM)      , WQRPRL(LCM)   ,              &
           WQRPRRP(LCM)     , WQRPRLP(LCM)     , WQRPRRN(LCM)     , WQRPRLN(LCM)  ,              &
           FRPSPW(LCM)      , FRPSNW(LCM)      , PNRPS(LCM)       , PNRPE(LCM)    ,RJRPRS(LCM),  &
           XLIMTPRPS(LCM)   , XLIMTPRPE(LCM)   , XLIMTRRPS(LCM)   , XLIMTRRPE(LCM),              &
           XLIMTRRPR(LCM)   , XLIMNRPS(LCM)    , XLIMNRPE(LCM)    , XLIMLRPS(LCM) ,              &
           XLIMLRPE(LCM)    , RISS(LCM)        , RISSO(LCM)       , RISSOE(LCM)   ,              &
           RPEMTPRPS(NWQTDM), RPEMTPRPE(NWQTDM), RPEMTRRPS(NWQTDM), &
           RPEMTRRPE(NWQTDM), RPEMTRRPR(NWQTDM)  )
  
     
  ALLOCATE(IRPEMTS(LCM), JRPEMTS(LCM), LRPEMTS(LCM))
  ALLOCATE(LMASKRPEM(LCM))

  WQRPS =0
  WQRPR =0
  WQRPE =0
  WQRPD =0
  PRPS =0
  RRPS =0
  RRPR =0
  PRPE =0
  RRPE =0
  WQRPSR =0
  WQRPSL =0
  WQRPER =0
  WQRPEL =0
  WQRPDL =0
  WQRPSRP =0
  WQRPSLP =0
  WQRPERP =0
  WQRPELP =0
  WQRPDLP =0
  WQRPSRN =0
  WQRPSLN =0
  WQRPERN =0
  WQRPELN =0
  WQRPDLN =0
  WQRPRR =0
  WQRPRL =0
  WQRPRRP =0
  WQRPRLP =0
  WQRPRRN =0
  WQRPRLN =0
  FRPSPW =0
  FRPSNW =0
  PNRPS =0
  PNRPE =0
  RJRPRS =0
  XLIMTPRPS =0
  XLIMTPRPE =0
  XLIMTRRPS =0
  XLIMTRRPE =0
  XLIMTRRPR =0
  XLIMNRPS =0
  XLIMNRPE =0
  XLIMLRPS =0
  XLIMLRPE =0
  RISS =0
  RISSO =0
  RISSOE =0
  RPEMTPRPS=0
  RPEMTPRPE=0
  RPEMTRRPS=0
  RPEMTRRPE=0
  RPEMTRRPR=0

  JSRPEM=1
  
END SUBROUTINE

END MODULE
