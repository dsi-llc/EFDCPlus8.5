! *** SUBROUTINES BAL2T CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, AND ENERGY BALANCES

! *** TOXICS UNITS
! ***   TOX(L,K)  - MG/M3 = UG/L = PPB
! ***   TOXB(L,K) - MG/M2
! ***   THEREFORE, THE MASS BALANCE FOR TOXICS WILL BE IN MG INSTEAD OF GRAMS FOR THE OTHER PARAMETERS

  SUBROUTINE BAL2T1

  ! *** SUBROUTINES BAL2T CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, AND ENERGY BALANCES
  ! *** BAL2T1 - INITIALIZES VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC ENERGY
  ! ***          AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-02-01        Paul M. Craig    REFORMATTED AND CLEANED UP CODE

  USE GLOBAL

  IMPLICIT NONE                                                                                                          
  INTEGER :: NS,L,K,NT,LN,LE                                                                                            

  IF( NBAL > 1 ) RETURN

  VOLBEG=0.
  VOLBEG2T=0.
  WVOLBEG2T=0.
  BVOLBEG2T=0.
  SALBEG=0.
  DYEBEG=0.
  UMOBEG=0.
  VMOBEG=0.
  UUEBEG=0.
  VVEBEG=0.
  PPEBEG=0.
  BBEBEG=0.
  DYEBEG2T=0.0
  DO NS=1,NSED
    SEDBEG2T(NS)=0.0
    SEDBEG2TW(NS)=0.0
    SEDBEG2TB(NS)=0.0
  ENDDO
  DO NS=1,NSND
    SNDBEG2T(NS)=0.0
    SNDBEG2TW(NS)=0.0
    SNDBEG2TB(NS)=0.0
  ENDDO
  DO NT=1,NTOX
    TOXBEG2T(NT)=0.0     ! *** MASS OF TOTAL TOXICS IN WATER AND SEDIMENT AT END OF MASS BALANCE TIME STEP
    TOXBEG2TW(NT)=0.0    ! *** MASS OF TOTAL TOXICS IN WATER AT END OF MASS BALANCE TIME STEP
    TOXBEG2TB(NT)=0.0    ! *** MASS OF TOTAL TOXICS IN SEDIMENT AT END OF MASS BALANCE TIME STEP
  ENDDO
  VOLOUT=0.
  WVOLOUT=0.
  BVOLOUT=0.
  SALOUT=0.
  DYEOUT=0.
  UMOOUT=0.
  VMOOUT=0.
  UUEOUT=0.
  VVEOUT=0.
  PPEOUT=0.
  BBEOUT=0.
  DYEOUT2T=0.0
  VOLMORPH2T=0.0
  DO NS=1,NSED
    SEDOUT2T(NS)=0.0
    SEDFLUX2T(NS)=0.0
  ENDDO
  DO NS=1,NSND
    SNDOUT2T(NS)=0.0
    SNDFLUX2T(NS)=0.0
    SBLOUT2T(NS)=0.0
    SNDFBL2T(NS)=0.0
  ENDDO
  DO NT=1,NTOX
    TOXOUT2T(NT)=0.0
    TOXFLUXW2T(NT)=0.0
    TOXFLUXB2T(NT)=0.0
    TADFLUX2T(NT)=0.0
    TOXBLB2T(NT)=0.0
    TOXFBL2T(NT)=0.0
  ENDDO
  
  DO L=2,LA
    LE=LEC(L)
    LN=LNC(L)
    VOLBEG    = VOLBEG   + SPB(L)*DXYP(L)*HP(L)
    VOLBEG2T  = VOLBEG2T + SPB(L)*DXYP(L)*HP(L)
    WVOLBEG2T = WVOLBEG2T+ SPB(L)*DXYP(L)*HP(L)
    UMOBEG    = UMOBEG   + SPB(L)*0.5*DXYP(L)*HP(L)*(DYIU(L)*HUI(L)*UHDYE(L)+DYIU(LE)*HUI(LE)*UHDYE(LE))
    VMOBEG    = VMOBEG   + SPB(L)*0.5*DXYP(L)*HP(L)*(DXIV(L)*HVI(L)*VHDXE(L)+DXIV(LN)*HVI(LN)*VHDXE(LN))
    PPEBEG    = PPEBEG   + SPB(L)*0.5*DXYP(L)*(GI*P(L)*P(L)-G*BELV(L)*BELV(L))
  ENDDO
  AMOBEG=SQRT(UMOBEG*UMOBEG+VMOBEG*VMOBEG)
  
  DO K=1,KC
    DO L=2,LA
      LE=LEC(L)
      LN=LNC(L)
      DYEBEG   = DYEBEG   + SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(L,K)
      DYEBEG2T = DYEBEG2T + SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(L,K)
      SALBEG   = SALBEG   + SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(L,K)
      UUEBEG   = UUEBEG   + SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K) *( (U(L,K)+U(LE,K))*(U(L,K)+U(LE,K)) )
      VVEBEG   = VVEBEG   + SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K)*( (V(L,K)+V(LN,K))*(V(L,K)+V(LN,K)) )
      BBEBEG   = BBEBEG   + SPB(L)*GP*DXYP(L)*HP(L)*DZC(L,K)*( BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B(L,K)
    ENDDO
  ENDDO
  
  ! *** MASS - WATER COLUMN
  DO NS=1,NSED
    DO K=1,KC
      DO L=2,LA
        SEDBEG2T(NS)  = SEDBEG2T(NS)  + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SED(L,K,NS)
        SEDBEG2TW(NS) = SEDBEG2TW(NS) + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SED(L,K,NS)
      ENDDO
    ENDDO
  ENDDO
  
  DO NS=1,NSND
    DO K=1,KC
      DO L=2,LA
        SNDBEG2T(NS)  = SNDBEG2T(NS)  + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SND(L,K,NS)
        SNDBEG2TW(NS) = SNDBEG2TW(NS) + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SND(L,K,NS)
      ENDDO
    ENDDO
  ENDDO
  
  DO NT=1,NTOX
    DO K=1,KC
      DO L=2,LA
        TOXBEG2T(NT)  = TOXBEG2T(NT)  + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*TOX(L,K,NT)
        TOXBEG2TW(NT) = TOXBEG2TW(NT) + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*TOX(L,K,NT)
      ENDDO
    ENDDO
  ENDDO
  
  ! *** MASS - SEDIMENT BED
  DO NS=1,NSED
    DO L=2,LA
      DO K=1,KBT(L)
        SEDBEG2T(NS)  = SEDBEG2T(NS)  + SCB(L)*DXYP(L)*SEDB(L,K,NS)
        SEDBEG2TB(NS) = SEDBEG2TB(NS) + SCB(L)*DXYP(L)*SEDB(L,K,NS)
      ENDDO
    ENDDO
  ENDDO
  
  DO NS=1,NSND
    DO L=2,LA
      DO K=1,KBT(L)
        SNDBEG2T(NS)  = SNDBEG2T(NS)  + SCB(L)*DXYP(L)*SNDB(L,K,NS)
        SNDBEG2TB(NS) = SNDBEG2TB(NS) + SCB(L)*DXYP(L)*SNDB(L,K,NS)
      ENDDO
    ENDDO
  ENDDO
  
  DO NT=1,NTOX
    DO L=2,LA
      DO K=1,KBT(L)
        TOXBEG2T(NT)  = TOXBEG2T(NT)  + SCB(L)*DXYP(L)*TOXB(L,K,NT)
        TOXBEG2TB(NT) = TOXBEG2TB(NT) + SCB(L)*DXYP(L)*TOXB(L,K,NT)
      ENDDO
    ENDDO
  ENDDO
  
  DO L=2,LA
    DO K=1,KBT(L)
      BVOLBEG2T = BVOLBEG2T + SPB(L)*DXYP(L)*HBED(L,K)
      VOLBEG2T  = VOLBEG2T  + SPB(L)*DXYP(L)*HBED(L,K)
    ENDDO
  ENDDO

  RETURN

END


SUBROUTINE BAL2T2

  ! *** SUBROUTINES BAL2T CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, AND ENERGY BALANCES
  ! *** BAL2T2 - ACCUMULATES THE OPEN BC MASS AND ENERGY FLUXES

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-02-01        Paul M. Craig    REFORMATTED AND CLEANED UP CODE


  USE GLOBAL

  IMPLICIT NONE                                                                                                         
  INTEGER :: LS,NS,NT,LL,K,L,LN                                                                                            
  REAL    :: VOLOPEN
  
  IF( ISDYNSTP == 0 )THEN
    DELT=DT
  ELSE
    DELT=DTDYN
  END IF
  VOLOPEN = 0
  
  ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES
  DO K=1,KC
    DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      VOLOPEN = VOLOPEN - DELT*VHDX2(LN,K)
      VOLOUT  = VOLOUT  - DELT*VHDX2(LN,K)
      WVOLOUT = WVOLOUT - DELT*VHDX2(LN,K)
      SALOUT  = SALOUT  - DELT*MIN(VHDX2(LN,K),0.)*SAL(LN,K) - DELT*MAX(VHDX2(LN,K),0.)*SAL(L,K)
      DYEOUT  = DYEOUT  - DELT*MIN(VHDX2(LN,K),0.)*DYE(LN,K) - DELT*MAX(VHDX2(LN,K),0.)*DYE(L,K)
      DO NS=1,NSED
        SEDOUT2T(NS) = SEDOUT2T(NS) - DELT*MIN(VHDX2(LN,K),0.)*SED(LN,K,NS) - DELT*MAX(VHDX2(LN,K),0.)*SED(L,K,NS)
      ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS) = SNDOUT2T(NS) - DELT*MIN(VHDX2(LN,K),0.)*SND(LN,K,NS) - DELT*MAX(VHDX2(LN,K),0.)*SND(L,K,NS)
      ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT) = TOXOUT2T(NT) - DELT*MIN(VHDX2(LN,K),0.)*TOX(LN,K,NT) - DELT*MAX(VHDX2(LN,K),0.)*TOX(L,K,NT)
      ENDDO
      PPEOUT = PPEOUT - DELT*VHDX2(LN,K)*G*(0.5*(BELV(L)+BELV(LN)) + 0.125*(HP(L)+H1P(L)+HP(LN)+H1P(LN))*(Z(L,K)+Z(L,K-1)) )
      BBEOUT = BBEOUT - DELT*MIN(VHDX2(LN,K),0.)*GP*( BELV(LN) + 0.5*HP(LN)*(Z(L,K)+Z(L,K-1)) )*B(LN,K) &
                      - DELT*MAX(VHDX2(LN,K),0.)*GP*( BELV(L ) + 0.5*HP(L )*(Z(L,K)+Z(L,K-1)) )*B(L,K)
      ! ***
      ! ***
      ! ***
    ENDDO
  ENDDO
  
  DO K=1,KC
    DO LL=1,NCBW
      L=LCBW(LL)
      VOLOPEN = VOLOPEN - DELT*UHDY2(LEC(L),K)
      VOLOUT  = VOLOUT  - DELT*UHDY2(LEC(L),K)
      WVOLOUT = WVOLOUT - DELT*UHDY2(LEC(L),K)
      SALOUT  = SALOUT  - DELT*MIN(UHDY2(LEC(L),K),0.)*SAL(LEC(L),K) - DELT*MAX(UHDY2(LEC(L),K),0.)*SAL(L,K)
      DYEOUT  = DYEOUT  - DELT*MIN(UHDY2(LEC(L),K),0.)*DYE(LEC(L),K) - DELT*MAX(UHDY2(LEC(L),K),0.)*DYE(L,K)
      DO NS=1,NSED
        SEDOUT2T(NS) = SEDOUT2T(NS) - DELT*MIN(UHDY2(LEC(L),K),0.)*SED(LEC(L),K,NS) - DELT*MAX(UHDY2(LEC(L),K),0.)*SED(L,K,NS)
      ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS) = SNDOUT2T(NS) - DELT*MIN(UHDY2(LEC(L),K),0.)*SND(LEC(L),K,NS) - DELT*MAX(UHDY2(LEC(L),K),0.)*SND(L,K,NS)
      ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT) = TOXOUT2T(NT) - DELT*MIN(UHDY2(LEC(L),K),0.)*TOX(LEC(L),K,NT) - DELT*MAX(UHDY2(LEC(L),K),0.)*TOX(L,K,NT)
      ENDDO
      PPEOUT = PPEOUT - DELT*UHDY2(LEC(L),K)*G*(0.5*(BELV(L)+ BELV(LEC(L)))+0.125*(HP(L)+H1P(L)+HP(LEC(L))+H1P(LEC(L)))*(Z(L,K)+Z(L,K-1)))
      BBEOUT = BBEOUT - DELT*MIN(UHDY2(LEC(L),K),0.)*GP*( BELV(LEC(L))+0.5*HP(LEC(L))*(Z(L,K)+Z(L,K-1)) )*B(LEC(L),K) &
                      - DELT*MAX(UHDY2(LEC(L),K),0.)*GP*( BELV(L  )+0.5*HP(L  )*(Z(L,K)+Z(L,K-1)) )*B(L,K)
      ! ***
      ! ***
      ! ***
    ENDDO
  ENDDO
  
  DO K=1,KC
    DO LL=1,NCBE
      L=LCBE(LL)
      VOLOPEN = VOLOPEN + DELT*UHDY2(L,K)
      VOLOUT  = VOLOUT  + DELT*UHDY2(L,K)
      WVOLOUT = WVOLOUT + DELT*UHDY2(L,K)
      SALOUT  = SALOUT  + DELT*MIN(UHDY2(L,K),0.)*SAL(L,K)+DELT*MAX(UHDY2(L,K),0.)*SAL(LWC(L),K)
      DYEOUT  = DYEOUT  + DELT*MIN(UHDY2(L,K),0.)*DYE(L,K)+DELT*MAX(UHDY2(L,K),0.)*DYE(LWC(L),K)
      DO NS=1,NSED
        SEDOUT2T(NS) = SEDOUT2T(NS) + DELT*MIN(UHDY2(L,K),0.)*SED(L,K,NS) + DELT*MAX(UHDY2(L,K),0.)*SED(LWC(L),K,NS)
      ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS) = SNDOUT2T(NS) + DELT*MIN(UHDY2(L,K),0.)*SND(L,K,NS) + DELT*MAX(UHDY2(L,K),0.)*SND(LWC(L),K,NS)
      ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT) = TOXOUT2T(NT) + DELT*MIN(UHDY2(L,K),0.)*TOX(L,K,NT) + DELT*MAX(UHDY2(L,K),0.)*TOX(LWC(L),K,NT)
      ENDDO
      PPEOUT = PPEOUT + DELT*UHDY2(L,K)*G*(0.5*(BELV(L)+BELV(LWC(L)))+0.125*(HP(L)+H1P(L)+HP(LWC(L))+H1P(LWC(L)))*(Z(L,K)+Z(L,K-1)) )
      BBEOUT = BBEOUT + DELT*MIN(UHDY2(L,K),0.)*GP*(BELV(L  )+0.5*HP(L  )*(Z(L,K)+Z(L,K-1)) )*B(L,K) &
                      + DELT*MAX(UHDY2(L,K),0.)*GP*(BELV(LWC(L))+0.5*HP(LWC(L))*(Z(L,K)+Z(L,K-1)) )*B(LWC(L),K)
      ! ***
      ! ***
      ! ***
    ENDDO
  ENDDO
  
  DO K=1,KC
    DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      VOLOPEN = VOLOPEN + DELT*VHDX2(L,K)
      VOLOUT  = VOLOUT  + DELT*VHDX2(L,K)
      WVOLOUT = WVOLOUT + DELT*VHDX2(L,K)
      SALOUT = SALOUT   + DELT*MIN(VHDX2(L,K),0.)*SAL(L,K)+DELT*MAX(VHDX2(L,K),0.)*SAL(LS,K)
      DYEOUT = DYEOUT   + DELT*MIN(VHDX2(L,K),0.)*DYE(L,K)+DELT*MAX(VHDX2(L,K),0.)*DYE(LS,K)
      DO NS=1,NSED
        SEDOUT2T(NS) = SEDOUT2T(NS) + DELT*MIN(VHDX2(L,K),0.)*SED(L,K,NS) + DELT*MAX(VHDX2(L,K),0.)*SED(LS,K,NS)
      ENDDO
      DO NS=1,NSND
        SNDOUT2T(NS) = SNDOUT2T(NS) + DELT*MIN(VHDX2(L,K),0.)*SND(L,K,NS) + DELT*MAX(VHDX2(L,K),0.)*SND(LS,K,NS)
      ENDDO
      DO NT=1,NTOX
        TOXOUT2T(NT) = TOXOUT2T(NT) + DELT*MIN(VHDX2(L,K),0.)*TOX(L,K,NT) + DELT*MAX(VHDX2(L,K),0.)*TOX(LS,K,NT)
      ENDDO
      PPEOUT = PPEOUT + DELT*VHDX2(L,K)*G*(0.5*(BELV(L)+BELV(LS))+0.125*(HP(L)+H1P(L)+HP(LS)+H1P(LS))*(Z(L,K)+Z(L,K-1)) )
      BBEOUT = BBEOUT + DELT*MIN(VHDX2(L,K),0.)*GP*( BELV(L )+0.5*HP(L )*(Z(L,K)+Z(L,K-1)) )*B(L,K) &
                      + DELT*MAX(VHDX2(L,K),0.)*GP*( BELV(LS)+0.5*HP(LS)*(Z(L,K)+Z(L,K-1)) )*B(LS,K)
      ! ***
      ! ***
      ! ***
    ENDDO
  ENDDO
  
  ! ***
  
  RETURN
  
END

SUBROUTINE BAL2T3A

  ! *** SUBROUTINES BAL2T CALCULATE GLOBAL CONSTITUENT MASS, MOMENTUM, AND ENERGY BALANCES
  ! *** BAL2T3A - ACCUMULATES BOUNDARY SOURCES AND SINKS AND POTENTIAL ENERGY
  ! ***           AND ASSOCIATED FLUXES, EXCLUDING OPEN BC'S

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-02-01        Paul M. Craig    REFORMATTED AND CLEANED UP CODE
  ! 2002-XX-XX        JOHN HAMRICK     MODIFIED SND MASS BALANCE WITH RESPECT TO BED LOAD OUTFLOW
  !                                    ADDED QDWASTE TO WATER MASS BALANCE

  USE GLOBAL

  IMPLICIT NONE
  
  INTEGER :: LD,K,L,NSX,NS,NWR,NCTL,ID,JD,KU,NT,M,JU,LU,KD,LL,NQSTMP                                                       
  INTEGER :: IU,NCSTMP                                                                                                     
  
  REAL :: RQWD,QWRABS,VOLBCS

  IF( ISDYNSTP == 0 )THEN
    DELT=DT
  ELSE
    DELT=DTDYN
  END IF
  VOLBCS = 0
  
  ! **  ACCUMULATE INTERNAL SOURCES AND SINKS
  DO L=2,LA
    VOLBCS = VOLBCS   + DELT*(QSUME(L)-QDWASTE(L)) 
    VOLOUT  = VOLOUT  - DELT*(QSUME(L)-QDWASTE(L))
    WVOLOUT = WVOLOUT - DELT*(QSUME(L)-QDWASTE(L))
  ENDDO
  ! ***
  
  DO K=1,KC
    DO LL=1,NQSIJ
      L=LQS(LL)
      PPEOUT = PPEOUT - DELT*QSS(K,LL)*G*( 0.5*(BELV(L)+BELV(LWC(L)))+0.125*(HP(L)+H1P(L)+HP(LWC(L))+H1P(LWC(L)))*(Z(L,K)+Z(L,K-1)) )
    ENDDO
  ENDDO
  
  IF( ISTRAN(1) >= 1 )THEN
    DO K=1,KC
      DO L=2,LC
        CONT(L,K)=SAL(L,K)
      ENDDO
    ENDDO
    
    DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,1)
      DO K=1,KC
        ! ***             |        CONSTANT INFLOW         |   |        CONSTANT OUTFLOW      |  |              VARIABLE INFLOW               |   |          VARIABLE OUTFLOW         |
        SALOUT = SALOUT - DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,1) - DELT*MIN(QSS(K,NS),0.)*SAL(L,K) - DELT*MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,1) - DELT*MIN(QSERCELL(K,NS),0.)*SAL(L,K)
      ENDDO
    ENDDO
    
    DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF( ID == 0 .AND. JD == 0 )THEN
        DO K=1,KC
          SALOUT=SALOUT+DELT*QCTLT(K,NCTL,1)*CONT(LU,K)
        ENDDO
      ENDIF
    ENDDO
    
    DO NWR=1,NQWR
      NQSTMP=NQWRSERQ(NWR)
      NCSTMP=NQWRSERQ(NWR)
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      IF( QWRSERT(NQSTMP) >= 0. )THEN
        ! *** Original Withdrawal/Return
        IU=IQWRU(NWR)
        JU=JQWRU(NWR)
        KU=KQWRU(NWR)
        ID=IQWRD(NWR)
        JD=JQWRD(NWR)
        KD=KQWRD(NWR)
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        ID=IQWRU(NWR)
        JD=JQWRU(NWR)
        KD=KQWRU(NWR)
        IU=IQWRD(NWR)
        JU=JQWRD(NWR)
        KU=KQWRD(NWR)
      ENDIF
      QWRABS = ABS(QWRSERT(NQSTMP))
      LU=LIJ(IU,JU)
      LD=LIJ(ID,JD)

      SALOUT = SALOUT + DELT*( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
      IF( LD /= 1 .OR. LD /= LC )THEN
        SALOUT=SALOUT - DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,1))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,1)) )
      ENDIF
    ENDDO
  ENDIF
  
  IF( ISTRAN(3) >= 1 )THEN
    DO K=1,KC
      DO L=2,LC
        CONT(L,K)=DYE(L,K)
      ENDDO
    ENDDO
    
    DO NS=1,NQSIJ
      L=LQS(NS)
      NQSTMP=NQSERQ(NS)
      NCSTMP=NCSERQ(NS,3)
      DO K=1,KC
        DYEOUT   = DYEOUT   - DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,3) - DELT*MIN(QSS(K,NS),0.)*DYE(L,K) - DELT*MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,3) - DELT*MIN(QSERCELL(K,NS),0.)*DYE(L,K)
        DYEOUT2T = DYEOUT2T - DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,3) - DELT*MIN(QSS(K,NS),0.)*DYE(L,K) - DELT*MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,3) - DELT*MIN(QSERCELL(K,NS),0.)*DYE(L,K) 
      ENDDO
    ENDDO
    
    DO NCTL=1,NQCTL
      RQWD=1.
      IU=IQCTLU(NCTL)
      JU=JQCTLU(NCTL)
      LU=LIJ(IU,JU)
      ID=IQCTLD(NCTL)
      JD=JQCTLD(NCTL)
      IF( ID == 0 .AND. JD == 0 )THEN
        DO K=1,KC
          DYEOUT   = DYEOUT   + DELT*QCTLT(K,NCTL,1)*CONT(LU,K)
          DYEOUT2T = DYEOUT2T + DELT*QCTLT(K,NCTL,1)*CONT(LU,K)
        ENDDO
      ENDIF
    ENDDO
    
    DO NWR=1,NQWR
      NQSTMP=NQWRSERQ(NWR)
      NCSTMP=NQWRSERQ(NWR)
      ! *** Handle +/- Flows for Withdrawal/Return Structures
      IF( QWRSERT(NQSTMP) >= 0. )THEN
        ! *** Original Withdrawal/Return
        IU=IQWRU(NWR)
        JU=JQWRU(NWR)
        KU=KQWRU(NWR)
        ID=IQWRD(NWR)
        JD=JQWRD(NWR)
        KD=KQWRD(NWR)
      ELSE
        ! *** Reverse Flow Withdrawal/Return
        ID=IQWRU(NWR)
        JD=JQWRU(NWR)
        KD=KQWRU(NWR)
        IU=IQWRD(NWR)
        JU=JQWRD(NWR)
        KU=KQWRD(NWR)
      ENDIF
      QWRABS = ABS(QWRSERT(NQSTMP))
      LU=LIJ(IU,JU)
      LD=LIJ(ID,JD)

      DYEOUT = DYEOUT   + DELT*( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
      DYEOUT2T=DYEOUT2T + DELT*( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
      IF( LD /= 1 .OR. LD /= LC )THEN
        DYEOUT   = DYEOUT   - DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,3))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,3)) )
        DYEOUT2T = DYEOUT2T - DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,3))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,3)) )
      ENDIF
    ENDDO
  ENDIF
  
  IF( ISTRAN(5) >= 1 )THEN
    DO NT=1,NTOX
      M=MSVTOX(NT)
      DO K=1,KC
        DO L=2,LC
          CONT(L,K)=TOX(L,K,NT)
        ENDDO
      ENDDO
      
      !  TOXOUT2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE
      !  TO WATER COLUMN VOLUME SOURCES AND SINKS
      DO NS=1,NQSIJ
        L=LQS(NS)
        NQSTMP=NQSERQ(NS)
        NCSTMP=NCSERQ(NS,M)
        DO K=1,KC
          TOXOUT2T(NT) = TOXOUT2T(NT) - DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,M) - DELT*MIN(QSS(K,NS),0.)*TOX(L,K,NT) - DELT*MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,M) - DELT*MIN(QSERCELL(K,NS),0.)*TOX(L,K,NT)
        ENDDO 
      ENDDO
      
      DO NCTL=1,NQCTL
        RQWD=1.
        IU=IQCTLU(NCTL)
        JU=JQCTLU(NCTL)
        LU=LIJ(IU,JU)
        ID=IQCTLD(NCTL)
        JD=JQCTLD(NCTL)
        IF( ID == 0 .AND. JD == 0 )THEN
          DO K=1,KC
            TOXOUT2T(NT) = TOXOUT2T(NT) + DELT*QCTLT(K,NCTL,1)*TOX(LU,K,NT)
          ENDDO
        ENDIF
      ENDDO
      
      DO NWR=1,NQWR
        NQSTMP=NQWRSERQ(NWR)
        NCSTMP=NQWRSERQ(NWR)
        ! *** Handle +/- Flows for Withdrawal/Return Structures
        IF( QWRSERT(NQSTMP) >= 0. )THEN
          ! *** Original Withdrawal/Return
          IU=IQWRU(NWR)
          JU=JQWRU(NWR)
          KU=KQWRU(NWR)
          ID=IQWRD(NWR)
          JD=JQWRD(NWR)
          KD=KQWRD(NWR)
        ELSE
          ! *** Reverse Flow Withdrawal/Return
          ID=IQWRU(NWR)
          JD=JQWRU(NWR)
          KD=KQWRU(NWR)
          IU=IQWRD(NWR)
          JU=JQWRD(NWR)
          KU=KQWRD(NWR)
        ENDIF
        QWRABS = ABS(QWRSERT(NQSTMP))
        LU=LIJ(IU,JU)
        LD=LIJ(ID,JD)

        TOXOUT2T(NT) = TOXOUT2T(NT) + DELT*( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
        IF( LD /= 1 .OR. LD /= LC )THEN
          TOXOUT2T(NT) = TOXOUT2T(NT) - DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  
  IF( ISTRAN(6) >= 1 )THEN
    DO NSX=1,NSED
      M=MSVSED(NSX)
      DO K=1,KC
        DO L=2,LC
          CONT(L,K)=SED(L,K,NSX)
        ENDDO
      ENDDO

      ! SEDOUT2T(NSX) IS IS NET COHESIVE MASS GOING OUT OF DOMAIN DUE
      !   TO WATER COLUMN VOLUME SOURCES AND SINKS
      DO NS=1,NQSIJ
        L=LQS(NS)
        NQSTMP=NQSERQ(NS)
        NCSTMP=NCSERQ(NS,M)
        DO K=1,KC
          SEDOUT2T(NSX) = SEDOUT2T(NSX) - DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,M) - DELT*MIN(QSS(K,NS),0.)*SED(L,K,NSX) - DELT*MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,M) - DELT*MIN(QSERCELL(K,NS),0.)*SED(L,K,NSX)
        ENDDO
      ENDDO
      
      DO NCTL=1,NQCTL
        RQWD=1.
        IU=IQCTLU(NCTL)
        JU=JQCTLU(NCTL)
        LU=LIJ(IU,JU)
        ID=IQCTLD(NCTL)
        JD=JQCTLD(NCTL)
        IF( ID == 0 .AND. JD == 0 )THEN
          DO K=1,KC
           SEDOUT2T(NSX) = SEDOUT2T(NSX) + DELT*QCTLT(K,NCTL,1)*CONT(LU,K)
          ENDDO
        ENDIF
      ENDDO
      
      DO NWR=1,NQWR
        NQSTMP=NQWRSERQ(NWR)
        NCSTMP=NQWRSERQ(NWR)
        ! *** Handle +/- Flows for Withdrawal/Return Structures
        IF( QWRSERT(NQSTMP) >= 0. )THEN
          ! *** Original Withdrawal/Return
          IU=IQWRU(NWR)
          JU=JQWRU(NWR)
          KU=KQWRU(NWR)
          ID=IQWRD(NWR)
          JD=JQWRD(NWR)
          KD=KQWRD(NWR)
        ELSE
          ! *** Reverse Flow Withdrawal/Return
          ID=IQWRU(NWR)
          JD=JQWRU(NWR)
          KD=KQWRU(NWR)
          IU=IQWRD(NWR)
          JU=JQWRD(NWR)
          KU=KQWRD(NWR)
        ENDIF
        QWRABS = ABS(QWRSERT(NQSTMP))
        LU=LIJ(IU,JU)
        LD=LIJ(ID,JD)

        SEDOUT2T(NSX) = SEDOUT2T(NSX) + DELT*( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
        IF( LD /= 1 .OR. LD /= LC )THEN
          SEDOUT2T(NSX) = SEDOUT2T(NSX) - DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  
  IF( ISTRAN(7) >= 1 )THEN
    DO NSX=1,NSND
      M=MSVSND(NSX)
      DO K=1,KC
        DO L=2,LC
          CONT(L,K)=SND(L,K,NSX)
        ENDDO
      ENDDO

      !  SNDOUT2T(NSX) IS NET NONCOHESIVE MASS GOING OUT OF DOMAIN DUE
      !  TO WATER COLUMN VOLUME SOURCES AND SINKS
      DO NS=1,NQSIJ
        L=LQS(NS)
        NQSTMP=NQSERQ(NS)
        NCSTMP=NCSERQ(NS,M)
        DO K=1,KC
          SNDOUT2T(NSX) = SNDOUT2T(NSX) - DELT*MAX(QSS(K,NS),0.)*CQS(K,NS,M) - DELT*MIN(QSS(K,NS),0.)*SND(L,K,NSX) - DELT*MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,M) - DELT*MIN(QSERCELL(K,NS),0.)*SND(L,K,NSX)
        ENDDO
      ENDDO
      
      DO NCTL=1,NQCTL
        RQWD=1.
        IU=IQCTLU(NCTL)
        JU=JQCTLU(NCTL)
        LU=LIJ(IU,JU)
        ID=IQCTLD(NCTL)
        JD=JQCTLD(NCTL)
        IF( ID == 0 .AND. JD == 0 )THEN
          DO K=1,KC
           SNDOUT2T(NSX) = SNDOUT2T(NSX)+DELT*QCTLT(K,NCTL,1)*CONT(LU,K)
          ENDDO
        ENDIF
      ENDDO
      
      DO NWR=1,NQWR
        NQSTMP=NQWRSERQ(NWR)
        NCSTMP=NQWRSERQ(NWR)
        ! *** Handle +/- Flows for Withdrawal/Return Structures
        IF( QWRSERT(NQSTMP) >= 0. )THEN
          ! *** Original Withdrawal/Return
          IU=IQWRU(NWR)
          JU=JQWRU(NWR)
          KU=KQWRU(NWR)
          ID=IQWRD(NWR)
          JD=JQWRD(NWR)
          KD=KQWRD(NWR)
        ELSE
          ! *** Reverse Flow Withdrawal/Return
          ID=IQWRU(NWR)
          JD=JQWRU(NWR)
          KD=KQWRU(NWR)
          IU=IQWRD(NWR)
          JU=JQWRD(NWR)
          KU=KQWRD(NWR)
        ENDIF
        QWRABS = ABS(QWRSERT(NQSTMP))
        LU=LIJ(IU,JU)
        LD=LIJ(ID,JD)

        SNDOUT2T(NSX) = SNDOUT2T(NSX) + DELT*( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
        IF( LD /= 1 .OR. LD /= LC )THEN
          SNDOUT2T(NSX) = SNDOUT2T(NSX) - DELT*( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,M))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,M)) )
        ENDIF
      ENDDO
    ENDDO
  ENDIF
    800 FORMAT('N,NS,SNDFBL2T,DEL',2I5,2E14.5)
  RETURN
END

SUBROUTINE BAL2T3B(IBALSTDT)

  ! *** SUBROUTINES BAL2T CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, AND ENERGY BALANCES
  ! *** BAL2T3B - ACCUMULATES INTERNAL FLUXES THAT RESULT IN MASS/ENERGY IN/OUT OF DOMAIN

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-02-01        Paul M. Craig    REFORMATTED AND CLEANED UP CODE

  USE GLOBAL
  
  IMPLICIT NONE
  
  INTEGER :: LUTMP,LDTMP,L,NSX,NSB,IBALSTDT,NT,M,LBANK,LCHAN,NP
	
  IF( ISDYNSTP == 0 )THEN
    DELT=DT
  ELSE
    DELT=DTDYN
  END IF

  ! **  ACCUMULATE INTERNAL SOURCES AND SINKS
  IF( IBALSTDT == 1 )THEN
    DO L=2,LA
      WVOLOUT=WVOLOUT-DTSED*QMORPH(L)
      BVOLOUT=BVOLOUT+DTSED*QMORPH(L)
      VOLMORPH2T=VOLMORPH2T+DTSED*QMORPH(L)
    ENDDO
  ENDIF
  
  IF( ISTRAN(5) >= 1 )THEN
    DO NT=1,NTOX
      M=MSVTOX(NT)

      ! *** TOXBLB2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
      IF( IBALSTDT == 1 )THEN
        IF( NSBDLDBC > 0 )THEN
          TOXBLB2T(NT)=TOXBLB2T(NT)+DTSED*TOXBLB(NT)
        ENDIF

        ! *** TOXFLUXW2T(NT) IS WATER COLUMN SIDE TOXIC FLUX DUE TO SUSPENDED LOAD (POSITIVE INTO WATER COLUMN)
        ! *** TOXFLUXB2T(NT) IS BED SIDE TOXIC FLUX DUE TO SUSPENDED LOAD          (POSITIVE INTO WATER COLUMN)
        ! *** TADFLUX2T(NT)  IS PORE WATER ADVECTION+DIFFUSION FLUX                (POSITIVE INTO WATER COLUMN)
        ! *** TOXFBL2T(NT)   IS NET TOXIC FLUX FROM BED ASSOCIATED WITH BED LOAD TRANPORT
        ! ***                (SHOULD EQUAL TOXBLB2T(NT)
        DO L=2,LA
          TOXFLUXW2T(NT)=TOXFLUXW2T(NT)+DTSED*DXYP(L)*TOXF(L,0,NT)
          TOXFLUXB2T(NT)=TOXFLUXB2T(NT)+DTSED*DXYP(L)*TOXFB(L,KBT(L),NT)
          TADFLUX2T(NT)=TADFLUX2T(NT)+DTSED*DXYP(L)*TADFLUX(L,NT)
        ENDDO

        TOXFBL2T(NT)=TOXFBL2T(NT)+DTSED*TOXFBLT(NT)
        IF( ISBKERO >= 1 )THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            TOXFLUXB2T(NT)=TOXFLUXB2T(NT) + DTSED*DXYP(LBANK)*TOXFBEBKB(LBANK,NT)
            TOXFLUXB2T(NT)=TOXFLUXB2T(NT) + DTSED*DXYP(LCHAN)*TOXFBECHB(LCHAN,NT)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF

  IF( ISTRAN(6) >= 1 )THEN
    DO NSX=1,NSED
      M=MSVSED(NSX)

      ! SEDFLUX2T(NSX) IS IS NET COHESIVE MASS FLUX POSITIVE FROM BED TO WATER COLUMN
      IF( IBALSTDT == 1 )THEN
        DO L=2,LA
          SEDFLUX2T(NSX)=SEDFLUX2T(NSX)+DTSED*DXYP(L)*SEDF(L,0,NSX)
        ENDDO
        IF( ISBKERO >= 1 )THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            SEDFLUX2T(NSX) = SEDFLUX2T(NSX) + DTSED*DXYP(LBANK)*SEDFBEBKB(LBANK,NSX)
            SEDFLUX2T(NSX) = SEDFLUX2T(NSX) + DTSED*DXYP(LCHAN)*SEDFBECHB(LCHAN,NSX)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  
  IF( ISTRAN(7) >= 1 )THEN
    DO NSX=1,NSND
      M=MSVSND(NSX)

      !  *** SBLOUT2T(NSX) IS NET NONCOHESIVE SEDIMENT MASS GOING OUT OF DOMAIN
      !  ***               DUE TO BED LOAD TRANSPORT OUT OF DOMAIN
      IF( IBALSTDT == 1 )THEN
        IF( NSBDLDBC > 0 )THEN
          DO NSB=1,NSBDLDBC
            LUTMP=LSBLBCU(NSB)
            LDTMP=LSBLBCD(NSB)
            IF( LDTMP == 0 )THEN
              SBLOUT2T(NSX)=SBLOUT2T(NSX) + DTSED*QSBDLDOT(LUTMP,NSX)
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      !  SNDFLUX2T(NSX) IS NET NONCOHESIVE SEDIMENT FLUX DUE TO SUSPENDED LOAD (POSITIVE INTO WATER COLUMN)
      !  SNDFBL2T(NSX) IS NET NONCOHESIVE SEDIMENT FLUX FROM BED ASSOCIATED WITH BED LOAD TRANSPORT (SHOULD EQUAL SBLOUT2T(NSX))
      IF( IBALSTDT == 1 )THEN
        DO L=2,LA
          SNDFLUX2T(NSX)=SNDFLUX2T(NSX)+DTSED*DXYP(L)*(SNDF(L,0,NSX) - SNDFBL(L,NSX))
          SNDFBL2T(NSX)=SNDFBL2T(NSX)+DTSED*DXYP(L)*SNDFBL(L,NSX)
        ENDDO
        IF( ISBKERO >= 1 )THEN
          DO NP=1,NBEPAIR
            LBANK=LIJ(IBANKBE(NP),JBANKBE(NP))
            LCHAN=LIJ(ICHANBE(NP),JCHANBE(NP))
            SNDFLUX2T(NSX)=SNDFLUX2T(NSX) + DTSED*DXYP(LBANK)*SNDFBEBKB(LBANK,NSX)
            SNDFLUX2T(NSX)=SNDFLUX2T(NSX) + DTSED*DXYP(LCHAN)*SNDFBECHB(LCHAN,NSX)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
    800 FORMAT('N,NS,SNDFBL2T,DEL',2I5,2E14.5)
  RETURN
END


SUBROUTINE BAL2T4

  ! *** SUBROUTINES BAL2T CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, AND ENERGY BALANCES
  ! *** BAL2T4 - ACCUMULATES MOMENTUM AND ENERGY DISSIPATION

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-02-01        Paul M. Craig    REFORMATTED AND CLEANED UP CODE

  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER :: L,LN,K
  REAL :: DUTMP,DVTMP

  IF( ISDYNSTP == 0 )THEN
    DELT=DT
  ELSE
    DELT=DTDYN
  END IF

  ! **  CALCULATE MOMENTUM AND ENERGY DISSIPATION
  DO L=2,LA
    LN=LNC(L)
    UUEOUT = UUEOUT + 0.5*DELT*SPB(L)*DXYP(L)*(U(L,KSZ(L))*TBX(L)+U(LEC(L),KSZ(L))*TBX(LEC(L))-U(L,KC)*TSX(L)-U(LEC(L),KC)*TSX(LEC(L)))
    VVEOUT = VVEOUT + 0.5*DELT*SPB(L)*DXYP(L)*(V(L,KSZ(L))*TBY(L)+V(LN ,KSZ(L))*TBX(LN )-V(L,KC)*TSY(L)-V(LN ,KC)*TSX(LN))
  ENDDO
  
  DO K=1,KS
    DO L=2,LA
      LN=LNC(L)
      IF( LKSZ(L,K) )CYCLE
      DUTMP=0.5*( U(L,K+1)+U(LEC(L),K+1)-U(L,K)-U(LEC(L),K) )
      DVTMP=0.5*( V(L,K+1)+V(LN,K+1)-V(L,K)-V(LN,K) )
      UUEOUT = UUEOUT + DELT*SPB(L)*2.0*DXYP(L)*AV(L,K)*( DUTMP*DUTMP )/(DZC(L,K+1)+DZC(L,K))
      VVEOUT = VVEOUT + DELT*SPB(L)*2.0*DXYP(L)*AV(L,K)*( DVTMP*DVTMP )/(DZC(L,K+1)+DZC(L,K))
      BBEOUT = BBEOUT + DELT*SCB(L)*DXYP(L)*HP(L)*GP*AB(L,K)*(B(L,K+1)-B(L,K))
      ! ***
      ! ***
      ! ***
    ENDDO
  ENDDO
  RETURN
END

SUBROUTINE BAL2T5

  ! *** SUBROUTINES BAL2T CALCULATE GLOBAL VOLUME, MASS, MOMENTUM, AND ENERGY BALANCES
  ! *** BAL2T5 - FINALIZES THE MASS/ENERGY BALANCE PERIOD AND WRITES THE PERIOD SUMMARIES

  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-02-01        Paul M. Craig    REFORMATTED AND CLEANED UP CODE


  USE GLOBAL
  IMPLICIT NONE
  
  INTEGER :: K,L,NS,NT,LN,LE
  REAL :: VOLBMO2T,BVOLBMO2T,WVOLBMO2T,SALBMO,DYEBMO,UMOBMO
  REAL :: VOLERR,VOLERR2T,BVOLERR2T,SALERR,DYEERR,WVOLERR2T
  REAL :: DYEERR2T,UMOERR,VMOERR,ENEERR,RVERDE,REERDE,RVERDO
  REAL :: RVERDO2T,RWVERDO2T,RSERDO,RDERDO,RDERDO2T,RUERDO,REERDO
  REAL :: RDERDE,RBVERDE2T,RDERDE2T,RUMERDE,RUMERDO
  REAL :: UUEBMO,RVMERDO,VVEBMO,PPEBMO,TMPVAL,BBEBMO,SBLOUT2TT
  REAL :: RVERDE2T,RWVERDE2T,RSERDE,RVMERDE,VOLBMO,VMOBMO,ENEBMO
  REAL :: ENEBEG,ENEEND,ENEOUT,AMOEND,VVEEND,PPEEND,BBEEND
  REAL :: UUEEND,VMOEND,SALEND,TIME,VOLEND,VOLEND2T,BVOLEND2T
  REAL :: WVOLEND2T,DYEEND,UMOEND

  !**********************************************************************C
  ! **  INCREMENT COUNTER
  IF( ISDYNSTP == 0 )THEN
    NBAL=NBAL+1
    TIME=DT*FLOAT(N)+TCON*TBEGIN
    TIME=TIME/TCON
  ELSE
    NBAL=NBAL+NINCRMT
    TIME=TIMESEC/TCON
  ENDIF

  ! **  CHECK FOR END OF BALANCE PERIOD
  IF( NBAL < NTSMMT ) RETURN

  !**********************************************************************C
  !
  ! **  CALCULATE ENDING VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC
  ! **  ENERGY AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES
  !
  !----------------------------------------------------------------------C
  VOLEND=0.
  VOLEND2T=0.
  BVOLEND2T=0.
  WVOLEND2T=0.
  SALEND=0.
  DYEEND=0.
  UMOEND=0.
  VMOEND=0.
  UUEEND=0.
  VVEEND=0.
  PPEEND=0.
  BBEEND=0.

  DYEEND2T=0.0
  DO NS=1,NSED
    SEDEND2T(NS)=0.0
    SEDEND2TW(NS)=0.0
    SEDEND2TB(NS)=0.0
  ENDDO
  DO NS=1,NSND
    SNDEND2T(NS)=0.0
    SNDEND2TW(NS)=0.0
    SNDEND2TB(NS)=0.0
  ENDDO
  DO NT=1,NTOX
    TOXEND2T(NT)=0.0
    TOXEND2TW(NT)=0.0
    TOXEND2TB(NT)=0.0
  ENDDO

  DO L=2,LA
    LE=LEC(L)
    LN=LNC(L)
    VOLEND    = VOLEND    + SPB(L)*DXYP(L)*HP(L)
    VOLEND2T  = VOLEND2T  + SPB(L)*DXYP(L)*HP(L)
    WVOLEND2T = WVOLEND2T + SPB(L)*DXYP(L)*HP(L)
    UMOEND    = UMOEND    + SPB(L)*0.5*DXYP(L)*HP(L)*(DYIU(L)*HUI(L)*UHDYE(L)+DYIU(LE)*HUI(LE)*UHDYE(LE))
    VMOEND    = VMOEND    + SPB(L)*0.5*DXYP(L)*HP(L)*(DXIV(L)*HVI(L)*VHDXE(L)+DXIV(LN )*HVI(LN )*VHDXE(LN ))
    PPEEND    = PPEEND    + SPB(L)*0.5*DXYP(L)*(GI*P(L)*P(L)-G*BELV(L)*BELV(L))
  ENDDO
  
  AMOEND=SQRT(UMOEND*UMOEND+VMOEND*VMOEND)

  DO K=1,KC
    DO L=2,LA
      LE=LEC(L)
      LN=LNC(L)
      SALEND   = SALEND   + SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(L,K)
      DYEEND   = DYEEND   + SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(L,K)
      DYEEND2T = DYEEND2T + SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(L,K)
      UUEEND   = UUEEND   + SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K)*( (U(L,K)+U(LE,K))*(U(L,K)+U(LE,K)) )
      VVEEND   = VVEEND   + SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K)*( (V(L,K)+V(LN,K))*(V(L,K)+V(LN,K)) )
      BBEEND   = BBEEND   + SPB(L)*GP*DXYP(L)*HP(L)*DZC(L,K)*( BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B(L,K)
    ENDDO
  ENDDO

  ! *** WATER COLUMN (SED,SND,TOX)
  DO NS=1,NSED
    DO K=1,KC
      DO L=2,LA
	    SEDEND2T(NS)  = SEDEND2T(NS)  + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SED(L,K,NS)
	    SEDEND2TW(NS) = SEDEND2TW(NS) + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SED(L,K,NS)
      ENDDO
    ENDDO
  ENDDO
  
  DO NS=1,NSND
    DO K=1,KC
      DO L=2,LA
	    SNDEND2T(NS)  = SNDEND2T(NS)  + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SND(L,K,NS)
	    SNDEND2TW(NS) = SNDEND2TW(NS) + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*SND(L,K,NS)
      ENDDO
    ENDDO
  ENDDO

  DO NT=1,NTOX
    DO K=1,KC
      DO L=2,LA
	    TOXEND2T(NT)  = TOXEND2T(NT)  + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*TOX(L,K,NT)
	    TOXEND2TW(NT) = TOXEND2TW(NT) + SCB(L)*DXYP(L)*HP(L)*DZC(L,K)*TOX(L,K,NT)
      ENDDO
    ENDDO
  ENDDO
  
  ! *** SEDIMENT BED (SED,SND,TOX)
  DO NS=1,NSED
    DO L=2,LA
      DO K=1,KBT(L)
	    SEDEND2T(NS)  = SEDEND2T(NS)  + SCB(L)*DXYP(L)*SEDB(L,K,NS)
	    SEDEND2TB(NS) = SEDEND2TB(NS) + SCB(L)*DXYP(L)*SEDB(L,K,NS)
      ENDDO
    ENDDO
  ENDDO
  
  DO NS=1,NSND
    DO L=2,LA
      DO K=1,KBT(L)
	    SNDEND2T(NS)  = SNDEND2T(NS)  + SCB(L)*DXYP(L)*SNDB(L,K,NS)
	    SNDEND2TB(NS) = SNDEND2TB(NS) + SCB(L)*DXYP(L)*SNDB(L,K,NS)
      ENDDO
    ENDDO
  ENDDO
  
  DO NT=1,NTOX
    DO L=2,LA
      DO K=1,KBT(L)
	    TOXEND2T(NT)  = TOXEND2T(NT)  + SCB(L)*DXYP(L)*TOXB(L,K,NT)
	    TOXEND2TB(NT) = TOXEND2TB(NT) + SCB(L)*DXYP(L)*TOXB(L,K,NT)
      ENDDO
    ENDDO
  ENDDO

  ! *** 
  DO L=2,LA
    DO K=1,KBT(L)
      VOLEND2T  = VOLEND2T  + SPB(L)*DXYP(L)*HBED(L,K)
      BVOLEND2T = BVOLEND2T + SPB(L)*DXYP(L)*HBED(L,K)
    ENDDO
  ENDDO

  ! *** ENERGY
  ! ***    U MOM + V MOM+ HEAD + DENSITY
  ENEBEG = UUEBEG+VVEBEG+PPEBEG+BBEBEG
  ENEEND = UUEEND+VVEEND+PPEEND+BBEEND
  ENEOUT = UUEOUT+VVEOUT+PPEOUT+BBEOUT
  
  VOLBMO=VOLBEG-VOLOUT
  VOLBMO2T=VOLBEG2T-VOLOUT
  BVOLBMO2T=BVOLBEG2T-BVOLOUT
  WVOLBMO2T=WVOLBEG2T-WVOLOUT
  SALBMO=SALBEG-SALOUT
  DYEBMO=DYEBEG-DYEOUT
  DYEBMO2T=DYEBEG2T-DYEOUT2T
  UMOBMO=UMOBEG-DYEOUT
  VMOBMO=VMOBEG-DYEOUT
  ENEBMO=ENEBEG-ENEOUT
                                                                                                                         
  !  TOXOUT2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE TO WATER COLUMN VOLUME SOURCES AND SINKS                                                                           
  !  TOXBLB2T(NT) IS NET TOXIC MASS GOING OUT OF DOMAIN DUE TO DUE TO BED LOAD TRANSPORT OUT OF DOMAIN                                                                            
  !  TOXFLUXW2T(NT) IS WATER COLUMN SIDE TOXIC FLUX DUE TO SUSPENDED LOAD (POSITIVE INTO WATER COLUMN)                                                                                       
  !  TOXFLUXB2T(NT) IS BED SIDE TOXIC FLUX DUE TO SUSPENDED LOAD          (POSITIVE INTO WATER COLUMN)                             
  !  TADFLUX2T(NT) IS PORE WATER ADVECTION+DIFFUSION FLUX                 (POSITIVE INTO WATER COLUMN)                                    
  !  TOXFBL2T(NT) IS NET TOXIC FLUX FROM BED ASSOCIATED WITH BED LOAD TRANSPORT (SHOULD EQUAL TOXBLB2T(NT)                                                                                         
  DO NT=1,NTOX                                                                                                           
    TOXOUT2TW(NT) = TOXOUT2T(NT)   - TOXFLUXW2T(NT) - TOXFLUXB2T(NT) - TADFLUX2T(NT)
    TOXOUT2TB(NT) = TOXFLUXW2T(NT) + TOXFLUXB2T(NT) + TADFLUX2T(NT)  + TOXFBL2T(NT)
    ! MODIFY TOXOUT2T TO INCLUDE BED LOAD BOUNDARY OUT                                                                      
    TOXOUT2T(NT)  = TOXOUT2T(NT)  + TOXBLB2T(NT)                                                                              
    TOXBMO2T(NT)  = TOXBEG2T(NT)  - TOXOUT2T(NT)                                                                               
    TOXBMO2TW(NT) = TOXBEG2TW(NT) - TOXOUT2TW(NT)                                                                            
    TOXBMO2TB(NT) = TOXBEG2TB(NT) - TOXOUT2TB(NT)                                                                            
  ENDDO                                                                                                                  

  ! SEDOUT2T(NS) IS IS NET COHESIVE MASS GOING OUT OF DOMAIN DUE TO WATER COLUMN VOLUME SOURCES AND SINKS                                                                            
  ! SEDFLUX2T(NS) IS IS NET COHESIVE MASS FLUX POSITIVE FROM BED TO WATER COLUMN                                                                                                 
  DO NS=1,NSED
    SEDOUT2TW(NS)=SEDOUT2T(NS)-SEDFLUX2T(NS)                                                                             
    SEDOUT2TB(NS)=SEDFLUX2T(NS)                                                                                          
    SEDBMO2T(NS)=SEDBEG2T(NS)-SEDOUT2T(NS)                                                                               
    SEDBMO2TW(NS)=SEDBEG2TW(NS)-SEDOUT2TW(NS)                                                                            
    SEDBMO2TB(NS)=SEDBEG2TB(NS)-SEDOUT2TB(NS)                                                                            
  ENDDO                                                                                                                  

  !  SNDOUT2T(NS) IS NET NONCOHESIVE MASS GOING OUT OF DOMAIN DUE TO WATER COLUMN VOLUME SOURCES AND SINKS
  !  SBLOUT2T(NS) IS NET NONCOHESIVE SEDIMENT MASS GOING OUT OF DOMAIN DUE TO DUE TO BED LOAD TRANSPORT OUT OF DOMAIN                                                                            
  !  SNDFLUX2T(NS) IS NET NONCOHESIVE SEDIMENT FLUX DUE TO SUSPENDED LOAD   (POSITIVE INTO WATER COLUMN)                                                                                       
  !  SNDFBL2T(NS) IS NET NONCOHESIVE SEDIMENT FLUX FROM BED ASSOCIATED WITH BED LOAD TRANSPORT (SHOULD EQUAL SBLOUT2T(NSX))                                                                    
  DO NS=1,NSND                                                                                                           
    SNDOUT2TW(NS)=SNDOUT2T(NS)-SNDFLUX2T(NS)                                                                             
    SNDOUT2TB(NS)=SNDFLUX2T(NS)+SNDFBL2T(NS)                                                                             
    !  MODIFY SNDOUT2T TO INCLUDE BED LOAD TRANSPORT OUT OF DOMAIN                                                          
    SNDOUT2T(NS)=SNDOUT2T(NS)+SBLOUT2T(NS)
    SNDBMO2T(NS)=SNDBEG2T(NS)-SNDOUT2T(NS)                                                                               
    SNDBMO2TW(NS)=SNDBEG2TW(NS)-SNDOUT2TW(NS)                                                                            
    SNDBMO2TB(NS)=SNDBEG2TB(NS)-SNDOUT2TB(NS)                                                                            
  ENDDO                                                                                                                  
                                                                                                                        
  VOLERR = VOLEND-VOLBMO
  VOLERR2T = VOLEND2T-VOLBMO2T
  BVOLERR2T = BVOLEND2T-BVOLBMO2T
  WVOLERR2T = WVOLEND2T-WVOLBMO2T
  SALERR = SALEND-SALBMO
  DYEERR = DYEEND-DYEBMO
  DYEERR2T = DYEEND2T-DYEBMO2T
  UMOERR = UMOEND-UMOBMO
  VMOERR = VMOEND-VMOBMO
  ENEERR = ENEEND-ENEBMO

  ! *** ABSOLUTE MASS ERROR TERM INITIALIZATION
  DO NS=1,NSED
    SEDERR2T(NS)  = SEDEND2T(NS)-SEDBMO2T(NS)                                                                               
    SEDERR2TW(NS) = SEDEND2TW(NS)-SEDBMO2TW(NS)                                                                            
    SEDERR2TB(NS) = SEDEND2TB(NS)-SEDBMO2TB(NS)                                                                            
  ENDDO                                                                                                                  
  DO NS=1,NSND                                                                                                           
    SNDERR2T(NS)=SNDEND2T(NS)-SNDBMO2T(NS)                                                                               
    SNDERR2TW(NS)=SNDEND2TW(NS)-SNDBMO2TW(NS)                                                                            
    SNDERR2TB(NS)=SNDEND2TB(NS)-SNDBMO2TB(NS)                                                                            
  ENDDO                                                                                                                  
  DO NT=1,NTOX                                                                                                           
    TOXERR2T(NT)=TOXEND2T(NT)-TOXBMO2T(NT)                                                                               
    TOXERR2TW(NT)=TOXEND2TW(NT)-TOXBMO2TW(NT)                                                                            
    TOXERR2TB(NT)=TOXEND2TB(NT)-TOXBMO2TB(NT)                                                                            
  ENDDO                                                                                                                  

  ! *** RELATIVE ERROR TERM INITIALIZATION
  DO NS=1,NSED
    RSEDERE2T(NS)=-9999.                                                                                                 
    RSEDERE2TW(NS)=-9999.                                                                                                
    RSEDERE2TB(NS)=-9999.                                                                                                
  ENDDO                                                                                                                  
  DO NS=1,NSND                                                                                                           
    RSNDERE2T(NS)=-9999.                                                                                                 
    RSNDERE2TW(NS)=-9999.                                                                                                
    RSNDERE2TB(NS)=-9999.                                                                                                
  ENDDO
  DO NT=1,NTOX                                                                                                           
    RTOXERE2T(NT)=-9999.                                                                                                 
    RTOXERE2TW(NT)=-9999.                                                                                                
    RTOXERE2TB(NT)=-9999.                                                                                                
  ENDDO                                                                                                                  

  RVERDO=-9999.
  RVERDO2T=-9999.
  RWVERDO2T=-9999.
  RSERDO=-9999.
  RDERDO=-9999.
  RDERDO2T=-9999.
  RUERDO=-9999.
  RVERDO=-9999.
  REERDO=-9999.
  RUMERDO=-9999.  ! PMC
  RVMERDO=-9999.  ! PMC

  DO NS=1,NSED
    RSEDERO2T(NS)=-9999.                                                                                                 
    RSEDERO2TW(NS)=-9999.                                                                                                
    RSEDERO2TB(NS)=-9999.                                                                                                
  ENDDO                                                                                                                  
  DO NS=1,NSND                                                                                                           
    RSNDERO2T(NS)=-9999.                                                                                                 
    RSNDERO2TW(NS)=-9999.                                                                                                
    RSNDERO2TB(NS)=-9999.                                                                                                
  ENDDO                                                                                                                  
  DO NT=1,NTOX                                                                                                           
    RTOXERO2T(NT)=-9999.                                                                                                 
    RTOXERO2TW(NT)=-9999.                                                                                                
    RTOXERO2TB(NT)=-9999.                                                                                                
  ENDDO                                                                                                                  
                                                                                                                       
  RVERDE=-9999.
  REERDE=-9999.
  RVERDE2T=-9999.
  RBVERDE2T=-9999.
  RWVERDE2T=-9999.
  RSERDE=-9999.
  RDERDE=-9999.
  RDERDE2T=-9999.
  RUMERDE=-9999.
  RVMERDE=-9999.

  IF( VOLEND /= 0. )    RVERDE    = VOLERR/VOLEND
  IF( VOLEND2T /= 0. )  RVERDE2T  = VOLERR2T/VOLEND2T
  IF( BVOLEND2T /= 0. ) RBVERDE2T = BVOLERR2T/BVOLEND2T
  IF( WVOLEND2T /= 0. ) RWVERDE2T = WVOLERR2T/WVOLEND2T
  IF( SALEND /= 0. )    RSERDE   = SALERR/SALEND
  IF( DYEEND /= 0. )    RDERDE   = DYEERR/DYEEND
  IF( DYEEND2T /= 0. )  RDERDE2T = DYEERR2T/DYEEND2T
  IF( UMOEND /= 0. )    RUMERDE  = UMOERR/UMOEND
  IF( VMOEND /= 0. )    RVMERDE  = VMOERR/VMOEND
  IF( ENEEND /= 0. )    REERDE   = ENEERR/ENEEND

  DO NS=1,NSED
    TMPVAL=0.5*(SEDEND2T(NS)+SEDBEG2T(NS))
    IF( TMPVAL /= 0. ) RSEDERE2T(NS)=SEDERR2T(NS)/TMPVAL
    
    TMPVAL=0.5*(SEDEND2TW(NS)+SEDBEG2TW(NS))
    IF( TMPVAL /= 0. ) RSEDERE2TW(NS)=SEDERR2TW(NS)/TMPVAL
    
    TMPVAL=0.5*(SEDEND2TB(NS)+SEDBEG2TB(NS))
    IF( TMPVAL /= 0. ) RSEDERE2TB(NS)=SEDERR2TB(NS)/TMPVAL
  ENDDO
  
  DO NS=1,NSND
    TMPVAL=0.5*(SNDEND2T(NS)+SNDBEG2T(NS))
    IF( TMPVAL /= 0. ) RSNDERE2T(NS)=SNDERR2T(NS)/TMPVAL
    
    TMPVAL=0.5*(SNDEND2TW(NS)+SNDBEG2TW(NS))
    IF( TMPVAL /= 0. ) RSNDERE2TW(NS)=SNDERR2TW(NS)/TMPVAL
    
    TMPVAL=0.5*(SNDEND2TB(NS)+SNDBEG2TB(NS))
    IF( TMPVAL /= 0. ) RSNDERE2TB(NS)=SNDERR2TB(NS)/TMPVAL
  ENDDO
  
  DO NT=1,NTOX
    TMPVAL=0.5*(TOXEND2T(NT)+TOXBEG2T(NT))
    IF( TMPVAL /= 0. ) RTOXERE2T(NT)=TOXERR2T(NT)/TMPVAL

    TMPVAL=0.5*(TOXEND2TW(NT)+TOXBEG2TW(NT))
    IF( TMPVAL /= 0. ) RTOXERE2TW(NT)=TOXERR2TW(NT)/TMPVAL
    
    TMPVAL=0.5*(TOXEND2TB(NT)+TOXBEG2TB(NT))
    IF( TMPVAL /= 0. ) RTOXERE2TB(NT)=TOXERR2TB(NT)/TMPVAL
  ENDDO
  
  IF( VOLOUT /= 0. ) RVERDO=VOLERR/VOLOUT
  IF( VOLOUT /= 0. ) RVERDO2T=VOLERR2T/VOLOUT
  IF( VOLOUT /= 0. ) RWVERDO2T=WVOLERR2T/VOLOUT
  IF( SALOUT /= 0. ) RSERDO=SALERR/SALOUT
  IF( DYEOUT /= 0. ) RDERDO=DYEERR/DYEOUT
  IF( DYEOUT2T /= 0. ) RDERDO2T=DYEERR2T/DYEOUT2T
  IF( UMOOUT /= 0. ) RUMERDO=UMOERR/UMOOUT
  IF( VMOOUT /= 0. ) RVMERDO=VMOERR/VMOOUT
  IF( ENEOUT /= 0. ) REERDO=ENEERR/ENEOUT

  DO NS=1,NSED
    IF( SEDOUT2T(NS)  >= 1E-10) RSEDERO2T(NS)  = SEDERR2T(NS)/SEDOUT2T(NS)
    IF( SEDOUT2TW(NS) >= 1E-10) RSEDERO2TW(NS) = SEDERR2TW(NS)/SEDOUT2TW(NS)
    IF( SEDOUT2TB(NS) >= 1E-10) RSEDERO2TB(NS) = SEDERR2TB(NS)/SEDOUT2TB(NS)
  ENDDO
  DO NS=1,NSND
    IF( SNDOUT2T(NS) >= 1E-10) RSNDERO2T(NS)  = SNDERR2T(NS)/SNDOUT2T(NS)
    IF( SNDOUT2TW(NS) >=  1E-10) RSNDERO2TW(NS) = SNDERR2TW(NS)/SNDOUT2TW(NS)
    IF( SNDOUT2TB(NS) >=  1E-10) RSNDERO2TB(NS) = SNDERR2TB(NS)/SNDOUT2TB(NS)
  ENDDO
  
  DO NT=1,NTOX
    IF( TOXOUT2T(NT) >= 1E-10)  RTOXERO2T(NT)  = TOXERR2T(NT)/TOXOUT2T(NT)
    IF( TOXOUT2TW(NT) >=  1E-10)  RTOXERO2TW(NT) = TOXERR2TW(NT)/TOXOUT2TW(NT)
    IF( TOXOUT2TB(NT) >=  1E-10)  RTOXERO2TB(NT) = TOXERR2TB(NT)/TOXOUT2TB(NT)
  ENDDO

  !**********************************************************************C
  !
  ! **  OUTPUT BALANCE RESULTS TO FILE BAL2T.OUT
  !
  !----------------------------------------------------------------------C
  IF( JSBAL == 1 )THEN
    OPEN(89,FILE=OUTDIR//'BAL2T.OUT')
    CLOSE(89,STATUS='DELETE')
    OPEN(89,FILE=OUTDIR//'BAL2T.OUT')
    OPEN(81,FILE=OUTDIR//'BAL2TERSTT.OUT')
    CLOSE(81,STATUS='DELETE')
    OPEN(81,FILE=OUTDIR//'BAL2TERSTT.OUT')
    OPEN(82,FILE=OUTDIR//'BAL2TERSTW.OUT')
    CLOSE(82,STATUS='DELETE')
    OPEN(82,FILE=OUTDIR//'BAL2TERSTW.OUT')
    OPEN(83,FILE=OUTDIR//'BAL2TERSTB.OUT')
    CLOSE(83,STATUS='DELETE')
    OPEN(83,FILE=OUTDIR//'BAL2TERSTB.OUT')
    OPEN(84,FILE=OUTDIR//'BAL2TERVWT.OUT')
    CLOSE(84,STATUS='DELETE')
    OPEN(84,FILE=OUTDIR//'BAL2TERVWT.OUT')
    JSBAL=0
  ELSE
    OPEN(89,FILE=OUTDIR//'BAL2T.OUT',POSITION='APPEND')
    OPEN(81,FILE=OUTDIR//'BAL2TERSTT.OUT',POSITION='APPEND')
    OPEN(82,FILE=OUTDIR//'BAL2TERSTW.OUT',POSITION='APPEND')
    OPEN(83,FILE=OUTDIR//'BAL2TERSTB.OUT',POSITION='APPEND')
    OPEN(84,FILE=OUTDIR//'BAL2TERVWT.OUT',POSITION='APPEND')
  ENDIF

  WRITE(89,890)NTSMMT,TIME

  WRITE(89,892)'INITIAL     (BEG)',VOLBEG,SALBEG,DYEBEG,ENEBEG,UMOBEG,VMOBEG,AMOBEG
  WRITE(89,892)'TOTAL OUT   (OUT)',VOLOUT,SALOUT,DYEOUT,ENEOUT,UMOOUT,VMOOUT
  WRITE(89,892)'INITIAL-OUT (BMO)',VOLBMO,SALBMO,DYEBMO,ENEBMO,UMOBMO,VMOBMO
  WRITE(89,892)'FINAL       (END)',VOLEND,SALEND,DYEEND,ENEEND,UMOEND,VMOEND,AMOEND
  WRITE(89,892)'TOTAL ERROR (ERR)',VOLERR,SALERR,DYEERR,ENEERR,UMOERR,VMOERR
  WRITE(89,892)'REL ERROR:ERR/END',RVERDE,RSERDE,RDERDE,REERDE,RUMERDE,RVMERDE
  WRITE(89,892)'REL ERROR:ERR/OUT',RVERDO,RSERDO,RDERDO,REERDO,RUMERDO,RVMERDO

  UUEBMO=UUEBEG-UUEOUT
  VVEBMO=VVEBEG-VVEOUT
  PPEBMO=PPEBEG-PPEOUT
  BBEBMO=BBEBEG-BBEOUT
  
  WRITE(89,900)
  WRITE(89,901)UUEBEG
  WRITE(89,902)UUEOUT
  WRITE(89,903)UUEBMO
  WRITE(89,904)UUEEND
  
  WRITE(89,900)
  WRITE(89,905)VVEBEG
  WRITE(89,906)VVEOUT
  WRITE(89,907)VVEBMO
  WRITE(89,908)VVEEND
  
  WRITE(89,900)
  WRITE(89,909)PPEBEG
  WRITE(89,910)PPEOUT
  WRITE(89,911)PPEBMO
  WRITE(89,912)PPEEND
  
  WRITE(89,900)
  WRITE(89,913)BBEBEG
  WRITE(89,914)BBEOUT
  WRITE(89,915)BBEBMO
  WRITE(89,916)BBEEND
  
  ! *****************************
  WRITE(81,8888)TIME,(RSEDERE2T(NS),NS=1,NSED),(RSNDERE2T(NS),NS=1,NSND),(RTOXERE2T(NT),NT=1,NTOX)
  WRITE(82,8888)TIME,(RSEDERE2TW(NS),NS=1,NSED),(RSNDERE2TW(NS),NS=1,NSND),(RTOXERE2TW(NT),NT=1,NTOX)
  WRITE(83,8888)TIME,(RSEDERE2TB(NS),NS=1,NSED),(RSNDERE2TB(NS),NS=1,NSND),(RTOXERE2TB(NT),NT=1,NTOX)
  WRITE(84,8888)TIME,RVERDE2T,RWVERDE2T,RBVERDE2T

  ! *****************************
  WRITE(89,900)
  WRITE(89,*)' NEW SEDIMENT-TOXIC VOLUME BALANCE: [WATER_COL + BED]'
  WRITE(89,*)' _BEG_ = BEGINNING GLOBAL VOLUME'
  WRITE(89,*)' _OUT_ = NET VOLUME GOING OUT OF DOMAIN'
  WRITE(89,*)' _BMO_ = BEGINNING - OUT'
  WRITE(89,*)' _END_ = ENDING GLOBAL VOLUME'
  WRITE(89,*)' _ERR_ = END - BMO'
  WRITE(89,*)' _R_E_ = RELATIVE ERROR = ERR/END'
  WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'

  WRITE(89,900)
  WRITE(89,954)
  WRITE(89,960)VOLBEG2T,VOLOUT,VOLBMO2T,VOLEND2T,VOLERR2T,RVERDE2T,RVERDO2T

  WRITE(89,900)
  WRITE(89,*)' NEW SEDIMENT-TOXIC VOLUME BALANCE: [WATER_COL]'
  WRITE(89,*)' _BEG_ = BEGINNING GLOBAL W_COL VOLUME'
  WRITE(89,*)' _OUT_ = NET W_COL VOLUME GOING OUT OF DOMAIN'
  WRITE(89,*)' _BMO_ = BEGINNING - OUT'
  WRITE(89,*)' _END_ = ENDING GLOBAL W_COL VOLUME'
  WRITE(89,*)' _ERR_ = END - BMO'
  WRITE(89,*)' _R_E_ = RELATIVE ERROR = ERR/END'
  WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'

  WRITE(89,900)
  WRITE(89,955)
  WRITE(89,960)WVOLBEG2T,WVOLOUT,WVOLBMO2T,WVOLEND2T,WVOLERR2T,RWVERDE2T,RWVERDO2T,VOLMORPH2T

  ! *****************************
  WRITE(89,900)

  WRITE(89,*)' NEW SEDIMENT-TOXIC MASS BALANCE: [WATER_COL + BED]'
  WRITE(89,*)' _BEG_ = BEGINNING GLOBAL MASS'
  WRITE(89,*)' _OUT_ = NET MASS GOING OUT OF DOMAIN'
  WRITE(89,*)' _BMO_ = BEGINNING - OUT'
  WRITE(89,*)' _END_ = ENDING GLOBAL MASS'
  WRITE(89,*)' _ERR_ = END - BMO'
  WRITE(89,*)' _R_E_ = RELATIVE ERROR = 2.*ERR/(END+BEG)'
  WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'
  WRITE(89,900)
  WRITE(89,949)

  NS=0
  WRITE(89,950)NS,DYEBEG2T,DYEOUT2T,DYEBMO2T,DYEEND2T,DYEERR2T,RDERDE2T,RDERDO2T
  
  SBLOUT2TT=0.0

  WRITE(89,900)
  WRITE(89,951)
  DO NS=1,NSED
    WRITE(89,950)NS,SEDBEG2T(NS),SEDOUT2T(NS),SEDBMO2T(NS),SEDEND2T(NS),SEDERR2T(NS),RSEDERE2T(NS),RSEDERO2T(NS),SEDFLUX2T(NS)
  ENDDO
  
  WRITE(89,900)
  WRITE(89,952)
  DO NS=1,NSND
	WRITE(89,950)NS,SNDBEG2T(NS),SNDOUT2T(NS),SNDBMO2T(NS),SNDEND2T(NS),SNDERR2T(NS),RSNDERE2T(NS),RSNDERO2T(NS),SNDFLUX2T(NS),SBLOUT2T(NS)
    DSEDGMM=1./(1.E6*SSG(NS+NSED))
	SBLOUT2TT=SBLOUT2TT+DSEDGMM*SBLOUT2T(NS)
  ENDDO
  WRITE(89,'(5X,8(17X),E17.9)')SBLOUT2TT
  
  WRITE(89,900)
  WRITE(89,953)
  DO NT=1,NTOX
    WRITE(89,950)NT,TOXBEG2T(NT),TOXOUT2T(NT),TOXBMO2T(NT),TOXEND2T(NT),TOXERR2T(NT),RTOXERE2T(NT),RTOXERO2T(NT),TOXFLUXW2T(NT),TOXFLUXB2T(NT),TADFLUX2T(NT)
  ENDDO
  WRITE(89,900)

  WRITE(89,*)' NEW SEDIMENT-TOXIC MASS BALANCE: [WATER_COLUMN ONLY]'
  WRITE(89,*)' _BEG_ = BEGINNING GLOBAL MASS'
  WRITE(89,*)' _OUT_ = NET MASS GOING OUT OF DOMAIN'
  WRITE(89,*)' _BMO_ = BEGINNING - OUT'
  WRITE(89,*)' _END_ = ENDING GLOBAL MASS'
  WRITE(89,*)' _ERR_ = END - BMO'
  WRITE(89,*)' _R_E_ = RELATIVE ERROR = 2.*ERR/(END+BEG)'
  WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'
  
  WRITE(89,900)
  WRITE(89,951)
  DO NS=1,NSED
    WRITE(89,950)NS,SEDBEG2TW(NS),SEDOUT2TW(NS),SEDBMO2TW(NS),SEDEND2TW(NS),SEDERR2TW(NS),RSEDERE2TW(NS),RSEDERO2TW(NS),SEDFLUX2T(NS)
  ENDDO
  
  WRITE(89,900)
  WRITE(89,952)
  DO NS=1,NSND
	WRITE(89,950)NS,SNDBEG2TW(NS),SNDOUT2TW(NS),SNDBMO2TW(NS),SNDEND2TW(NS),SNDERR2TW(NS),RSNDERE2TW(NS),RSNDERO2TW(NS),SNDFLUX2T(NS),SBLOUT2T(NS)
  ENDDO
  
  WRITE(89,900)
  WRITE(89,953)
  DO NT=1,NTOX
    WRITE(89,950)NT,TOXBEG2TW(NT),TOXOUT2TW(NT),TOXBMO2TW(NT),TOXEND2TW(NT),TOXERR2TW(NT),RTOXERE2TW(NT),RTOXERO2TW(NT),TOXFLUXW2T(NT),TOXFLUXB2T(NT),TADFLUX2T(NT)
  ENDDO
  
  WRITE(89,900)

  WRITE(89,*)' NEW SEDIMENT-TOXIC MASS BALANCE: [BED ONLY]'
  WRITE(89,*)' _BEG_ = BEGINNING GLOBAL MASS'
  WRITE(89,*)' _OUT_ = NET MASS GOING OUT OF DOMAIN'
  WRITE(89,*)' _BMO_ = BEGINNING - OUT'
  WRITE(89,*)' _END_ = ENDING GLOBAL MASS'
  WRITE(89,*)' _ERR_ = END - BMO'
  WRITE(89,*)' _R_E_ = RELATIVE ERROR = 2.*ERR/(END+BEG)'
  WRITE(89,*)' _R_O_ = RELATIVE ERROR = ERR/OUT'

  WRITE(89,900)
  WRITE(89,951)
  DO NS=1,NSED
    WRITE(89,950)NS,SEDBEG2TB(NS),SEDOUT2TB(NS),SEDBMO2TB(NS),SEDEND2TB(NS),SEDERR2TB(NS),RSEDERE2TB(NS),RSEDERO2TB(NS),SEDFLUX2T(NS)
  ENDDO
  
  WRITE(89,900)
  WRITE(89,952)
  DO NS=1,NSND
    WRITE(89,950)NS,SNDBEG2TB(NS),SNDOUT2TB(NS),SNDBMO2TB(NS),SNDEND2TB(NS),SNDERR2TB(NS),RSNDERE2TB(NS),RSNDERO2TB(NS),SNDFLUX2T(NS),SBLOUT2T(NS)
  ENDDO
  WRITE(89,'(5X,8(17X),E17.9)')SBLOUT2TT
  
  WRITE(89,900)
  WRITE(89,953)
  DO NT=1,NTOX
    WRITE(89,950)NT,TOXBEG2TB(NT),TOXOUT2TB(NT),TOXBMO2TB(NT),TOXEND2TB(NT),TOXERR2TB(NT),RTOXERE2TB(NT),RTOXERO2TB(NT),TOXFLUXW2T(NT),TOXFLUXB2T(NT),TADFLUX2T(NT)
  ENDDO
  
  WRITE(89,900)
  WRITE(89,*)' SEE NOTES IN BAL2T5 FOR INTERPRETATION OF FOLLOWING'

  WRITE(89,900)
  DO NS=1,NSED
    WRITE(89,900)
    WRITE(89,8899)'SEDBEG2T(NS)    = ',SEDBEG2T(NS)
    WRITE(89,8899)'SEDBEG2TB(NS)   = ',SEDBEG2TB(NS)
    WRITE(89,8899)'SEDBEG2TW(NS)   = ',SEDBEG2TW(NS)
    WRITE(89,8899)'SEDEND2T(NS)    = ',SEDEND2T(NS)
    WRITE(89,8899)'SEDEND2TB(NS)   = ',SEDEND2TB(NS)
    WRITE(89,8899)'SEDEND2TW(NS)   = ',SEDEND2TW(NS)
    WRITE(89,8899)'SEDOUT2T(NS)    = ',SEDOUT2T(NS)
    WRITE(89,8899)'SEDFLUX2T(NS)   = ',SEDFLUX2T(NS)
  ENDDO

  DO NS=1,NSND
    !  MODIFY SNDOUT2T BACK TO ORGINAL DEFINITION
    SNDOUT2T(NS)=SNDOUT2T(NS)-SBLOUT2T(NS)
    WRITE(89,900)
    WRITE(89,8899)'SNDBEG2T(NS)    = ',SNDBEG2T(NS)
    WRITE(89,8899)'SNDBEG2TB(NS)   = ',SNDBEG2TB(NS)
    WRITE(89,8899)'SNDBEG2TW(NS)   = ',SNDBEG2TW(NS)
    WRITE(89,8899)'SNDEND2T(NS)    = ',SNDEND2T(NS)
    WRITE(89,8899)'SNDEND2TB(NS)   = ',SNDEND2TB(NS)
    WRITE(89,8899)'SNDEND2TW(NS)   = ',SNDEND2TW(NS)
    WRITE(89,8899)'SNDOUT2T(NS)    = ',SNDOUT2T(NS)
    WRITE(89,8899)'SNDFLUX2T(NS)   = ',SNDFLUX2T(NS)
    WRITE(89,8899)'SNDFBL2T(NS)    = ',SNDFBL2T(NS)
    WRITE(89,8899)'SBLOUT2T(NS)    = ',SBLOUT2T(NS)
  ENDDO

  DO NT=1,NTOX
    ! MODIFY TOXOUT2T BACK TO ORIGINAL DEFINITION
    TOXOUT2T(NT)=TOXOUT2T(NT)-TOXBLB2T(NT)
    WRITE(89,900)
    WRITE(89,8899)'TOXBEG2T(NT)    = ',TOXBEG2T(NT)
    WRITE(89,8899)'TOXBEG2TB(NT)   = ',TOXBEG2TB(NT)
    WRITE(89,8899)'TOXBEG2TW(NT)   = ',TOXBEG2TW(NT)
    WRITE(89,8899)'TOXEND2T(NT)    = ',TOXEND2T(NT)
    WRITE(89,8899)'TOXEND2TB(NT)   = ',TOXEND2TB(NT)
    WRITE(89,8899)'TOXEND2TW(NT)   = ',TOXEND2TW(NT)
    WRITE(89,8899)'TOXOUT2T(NT)    = ',TOXOUT2T(NT)
    WRITE(89,8899)'TOXFLUXW2T(NT)  = ',TOXFLUXW2T(NT)
    WRITE(89,8899)'TOXFLUXB2T(NT)  = ',TOXFLUXB2T(NT)
    WRITE(89,8899)'TOXFBL2T(NT)    = ',TOXFBL2T(NT)
    WRITE(89,8899)'TOXBLB2T(NT)    = ',TOXBLB2T(NT)
    WRITE(89,8899)'TADFLUX2T(NT)   = ',TADFLUX2T(NT)
  ENDDO

  CLOSE(89)
  CLOSE(81)
  CLOSE(82)
  CLOSE(83)
  CLOSE(84)

    8899 FORMAT(A18,E15.7)
    950 FORMAT(I5,12E17.9)
    960 FORMAT(5X,12E17.9)
    949 FORMAT('  NDYE DYEBEG2T         DYEOUT2T         DYEBMO2T         DYEEND2T         DYEERR2T         RDYEERE2T        RDYEERO2T                                                           [UNITS G]',/)
    951 FORMAT('  NSED SEDBEG2T(NS)     SEDOUT2T(NS)     SEDBMO2T(NS)     SEDEND2T(NS)     SEDERR2T(NS)     RSEDERE2T(NS)    RSEDERO2T(NS)    SEDFLUX2T(NS)                                      [UNITS G]',/)
    952 FORMAT('  NSND SNDBEG2T(NS)     SNDOUT2T(NS)     SNDBMO2T(NS)     SNDEND2T(NS)     SNDERR2T(NS)     RSNDERE2T(NS)    RSNDERO2T(NS)    SNDFLUX2T(NS)    SBLOUT2T(NS)                      [UNITS G]',/)
    953 FORMAT('   NT  TOXBEG2T(NT)     TOXOUT2T(NT)     TOXBMO2T(NT)     TOXEND2T(NT)     TOXERR2T(NT)     RTOXERE2T(NT)    RTOXERO2T(NT)    TOXFLUXW2T(NT)    TOXFLUXB2T(NT)   TADFLUX2T(NT)   [UNITS MG]',/)
    954 FORMAT('       VOLBEG2T         VOLOUT           VOLBMO2T         VOLEND2T         VOLERR2T         RVERDE2T         RVERDO2T         VOLBLOUT                                           [UNITS M3]',/)
    955 FORMAT('       WVOLBEG2T        WVOLOUT          WVOLBMO2T         WVOLEND2T       WVOLERR2T        RWVERDE2T        RWVERDO2T        VOLMORPH2T                                         [UNITS M3]',/)
    890 FORMAT (///,130('*'),/,' VOLUME, MASS, AND ENERGY BALANCE OVER',I12,' TIME STEPS ENDING AT TIME ',F12.4,//,                  &
                19X,'    VOLUME(M3)        SALT(KG)          DYE(G)   ENERGY(M5/S2)   U MOM (M5/S2)   V MOM (M5/S2)  MAG MOM (M5/S2)')
    891 FORMAT (' INITIAL VOLUME    INITIAL SALT    INITIAL DYE     INITIAL ENER    INITIAL UMO     INITIAL VMO     INITIAL AMO',/)
    892 FORMAT (1X,A17,1X,7(E14.6,2X))
    900 FORMAT (//)
    901 FORMAT(' UUEBEG =  ',E14.6)
    902 FORMAT(' UUEOUT =  ',E14.6)
    903 FORMAT(' UUEBMO =  ',E14.6)
    904 FORMAT(' UUEEND =  ',E14.6)
    905 FORMAT(' VVEBEG =  ',E14.6)
    906 FORMAT(' VVEOUT =  ',E14.6)
    907 FORMAT(' VVEBMO =  ',E14.6)
    908 FORMAT(' VVEEND =  ',E14.6)
    909 FORMAT(' PPEBEG =  ',E14.6)
    910 FORMAT(' PPEOUT =  ',E14.6)
    911 FORMAT(' PPEBMO =  ',E14.6)
    912 FORMAT(' PPEEND =  ',E14.6)
    913 FORMAT(' BBEBEG =  ',E14.6)
    914 FORMAT(' BBEOUT =  ',E14.6)
    915 FORMAT(' BBEBMO =  ',E14.6)
    916 FORMAT(' BBEEND =  ',E14.6)
    8888 FORMAT(F10.2,10E14.6)

  !**********************************************************************C
  !
  !     RESET COUNTER
  !
  NBAL=0

  !**********************************************************************C
  RETURN

END

