SUBROUTINE CALBAL1

  ! CHANGE RECORD
  ! **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,
  ! **  AND ENERGY BALANCES

  USE GLOBAL
  IMPLICIT NONE                                                                                                          
  INTEGER :: L,LN,K                                                                                                        
  IF( NBAL>1) RETURN

  ! **  INITIALIZE VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC ENERGY                                                  
  ! **  AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES                                                                       
  VOLBEG=0.
  SALBEG=0.
  DYEBEG=0.
  UMOBEG=0.
  VMOBEG=0.
  UUEBEG=0.
  VVEBEG=0.
  PPEBEG=0.
  BBEBEG=0.
  VOLOUT=0.
  SALOUT=0.
  DYEOUT=0.
  UMOOUT=0.
  VMOOUT=0.
  UUEOUT=0.
  VVEOUT=0.
  PPEOUT=0.
  BBEOUT=0.
  DO L=2,LA
    LN=LNC(L)
    VOLBEG=VOLBEG+SPB(L)*DXYP(L)*HP(L)
    UMOBEG=UMOBEG+SPB(L)*0.5*DXYP(L)*HP(L)*(DYIU(L)*HUI(L)*UHDYE(L) + DYIU(LEC(L))*HUI(LEC(L))*UHDYE(LEC(L)))
    VMOBEG=VMOBEG+SPB(L)*0.5*DXYP(L)*HP(L)*(DXIV(L)*HVI(L)*VHDXE(L) + DXIV(LN)*HVI(LN)*VHDXE(LN))
    PPEBEG=PPEBEG+SPB(L)*0.5*DXYP(L)*(GI*P(L)*P(L)-G*BELV(L)*BELV(L))
  ENDDO
  AMOBEG=SQRT(UMOBEG*UMOBEG+VMOBEG*VMOBEG)
  DO K=1,KC
    DO L=2,LA
      LN=LNC(L)
      SALBEG=SALBEG+SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(L,K)
      DYEBEG=DYEBEG+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(L,K)
      UUEBEG=UUEBEG+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K)*( (U(L,K)+U(LEC(L),K))*(U(L,K)+U(LEC(L),K)) )
      VVEBEG=VVEBEG+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K)*( (V(L,K)+V(LN,K))*(V(L,K)+V(LN,K)) )
      BBEBEG=BBEBEG+SPB(L)*GP*DXYP(L)*HP(L)*DZC(L,K)*( BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B(L,K)
    ENDDO
  ENDDO
  RETURN
END

SUBROUTINE CALBAL2

  ! CHANGE RECORD
  ! **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,
  ! **  AND ENERGY BALANCES

  USE GLOBAL
  IMPLICIT NONE                                                                                                          
  INTEGER :: LL,K,LS,L,LN                                                                                                  

  ! **  ACCUMULATE FLUXES ACROSS OPEN BOUNDARIES                                                                          
  DO K=1,KC
    DO LL=1,NCBS
      L=LCBS(LL)
      LN=LNC(L)
      VOLOUT=VOLOUT-VHDX2(LN,K)
      SALOUT=SALOUT-MIN(VHDX2(LN,K),0.)*SAL1(LN,K)-MAX(VHDX2(LN,K),0.)*SAL1(L,K)
      DYEOUT=DYEOUT-MIN(VHDX2(LN,K),0.)*DYE1(LN,K)-MAX(VHDX2(LN,K),0.)*DYE1(L,K)
      PPEOUT=PPEOUT-VHDX2(LN,K)*G*( 0.5*(BELV(L)+BELV(LN))+0.125*(HP(L)+H2P(L)+HP(LN)+H2P(LN))*(Z(L,K)+Z(L,K-1)) )
      BBEOUT=BBEOUT-MIN(VHDX2(LN,K),0.)*GP*( BELV(LN)+0.5*HP(LN)*(Z(L,K)+Z(L,K-1)) )*B1(LN,K)-MAX(VHDX2(LN,K),0.)*GP*( BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B1(L,K)
    ENDDO
  ENDDO
  DO K=1,KC
    DO LL=1,NCBW
      L=LCBW(LL)
      VOLOUT=VOLOUT-UHDY2(LEC(L),K)
      SALOUT=SALOUT-MIN(UHDY2(LEC(L),K),0.)*SAL1(LEC(L),K)-MAX(UHDY2(LEC(L),K),0.)*SAL1(L,K)
      DYEOUT=DYEOUT-MIN(UHDY2(LEC(L),K),0.)*DYE1(LEC(L),K)-MAX(UHDY2(LEC(L),K),0.)*DYE1(L,K)
      PPEOUT=PPEOUT-UHDY2(LEC(L),K)*G*( 0.5*(BELV(L)+BELV(LEC(L)))+0.125*(HP(L)+H2P(L)+HP(LEC(L))+H2P(LEC(L)))*(Z(L,K)+Z(L,K-1)) )
      BBEOUT=BBEOUT-MIN(UHDY2(LEC(L),K),0.)*GP*( BELV(LEC(L))+0.5*HP(LEC(L))*(Z(L,K)+Z(L,K-1)) )*B1(LEC(L),K)-MAX(UHDY2(LEC(L),K),0.)*GP*( BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B1(L,K)
    ENDDO
  ENDDO
  DO K=1,KC
    DO LL=1,NCBE
      L=LCBE(LL)
      VOLOUT=VOLOUT+UHDY2(L,K)
      SALOUT=SALOUT+MIN(UHDY2(L,K),0.)*SAL1(L,K)+MAX(UHDY2(L,K),0.)*SAL1(LWC(L),K)
      DYEOUT=DYEOUT+MIN(UHDY2(L,K),0.)*DYE1(L,K)+MAX(UHDY2(L,K),0.)*DYE1(LWC(L),K)
      PPEOUT=PPEOUT+UHDY2(L,K)*G*( 0.5*(BELV(L)+BELV(LWC(L)))+0.125*(HP(L)+H2P(L)+HP(LWC(L))+H2P(LWC(L)))*(Z(L,K)+Z(L,K-1)) )
      BBEOUT=BBEOUT+MIN(UHDY2(L,K),0.)*GP*(BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B1(L,K)+MAX(UHDY2(L,K),0.)*GP*(BELV(LWC(L))+0.5*HP(LWC(L))*(Z(L,K)+Z(L,K-1)) )*B1(LWC(L),K)
    ENDDO
  ENDDO
  DO K=1,KC
    DO LL=1,NCBN
      L=LCBN(LL)
      LS=LSC(L)
      VOLOUT=VOLOUT+VHDX2(L,K)
      SALOUT=SALOUT+MIN(VHDX2(L,K),0.)*SAL1(L,K)+MAX(VHDX2(L,K),0.)*SAL1(LS,K)
      DYEOUT=DYEOUT+MIN(VHDX2(L,K),0.)*DYE1(L,K)+MAX(VHDX2(L,K),0.)*DYE1(LS,K)
      PPEOUT=PPEOUT+VHDX2(L,K)*G*( 0.5*(BELV(L)+BELV(LS))+0.125*(HP(L)+H2P(L)+HP(LS)+H2P(LS))*(Z(L,K)+Z(L,K-1)) )
      BBEOUT=BBEOUT+MIN(VHDX2(L,K),0.)*GP*( BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B1(L,K)+MAX(VHDX2(L,K),0.)*GP*( BELV(LS)+0.5*HP(LS)*(Z(L,K)+Z(L,K-1)) )*B1(LS,K)
    ENDDO
  ENDDO
  
  RETURN
END



SUBROUTINE CALBAL3

  ! CHANGE RECORD
  ! **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,
  ! **  AND ENERGY BALANCES

  USE GLOBAL

  IMPLICIT NONE                                                                                                          
  INTEGER :: L,K,LL,NS,NQSTMP,NCSTMP,NCTL,IU,JU,LU                                                                         
  INTEGER :: ID,JD,LD,NWR,KU,KD                                                                                            
  REAL :: RQWD,QWRABS                                                                                                      

  ! **  ACCUMULATE INTERNAL SOURCES AND SINKS                                                                             
  DO L=2,LA
    VOLOUT=VOLOUT-QSUME(L)
  ENDDO
  DO K=1,KC
    DO LL=1,NQSIJ
      L=LQS(LL)
      PPEOUT=PPEOUT-QSS(K,LL)*G*( 0.5*(BELV(L)+BELV(LWC(L)))+0.125*(HP(L)+H2P(L)+HP(LWC(L))+H2P(LWC(L)))*(Z(L,K)+Z(L,K-1)) )
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
        SALOUT=SALOUT-MAX(QSS(K,NS),0.)*CQS(K,NS,1)-MIN(QSS(K,NS),0.)*SAL(L,K)-MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,1)-MIN(QSERCELL(K,NS),0.)*SAL(L,K)
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
        LD=LC
        RQWD=0.
      ELSE
        LD=LIJ(ID,JD)
      ENDIF
      DO K=1,KC
        SALOUT=SALOUT+QCTLT(K,NCTL,1)*CONT(LU,K)-RQWD*QCTLT(K,NCTL,1)*CONT(LU,K)
      ENDDO
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

      SALOUT=SALOUT+( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
      IF( LD /= 1 .OR. LD /= LC )THEN
        SALOUT=SALOUT-( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,1))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,1)) )
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
        DYEOUT=DYEOUT-MAX(QSS(K,NS),0.)*CQS(K,NS,3)-MIN(QSS(K,NS),0.)*DYE(L,K)-MAX(QSERCELL(K,NS),0.)*CSERT(K,NCSTMP,3)-MIN(QSERCELL(K,NS),0.)*DYE(L,K)
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
        LD=LC
        RQWD=0.
      ELSE
        LD=LIJ(ID,JD)
      ENDIF
      DO K=1,KC
        DYEOUT=DYEOUT+QCTLT(K,NCTL,1)*CONT(LU,K)-RQWD*QCTLT(K,NCTL,1)*CONT(LU,K)
      ENDDO
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

      DYEOUT=DYEOUT+( (QWR(NWR)+QWRABS)*CONT(LU,KU) )
      IF( LD /= 1 .OR. LD /= LC )THEN
        DYEOUT=DYEOUT-( QWR(NWR)*(CONT(LU,KU)+CQWR(NWR,3))+QWRABS*(CONT(LU,KU)+CQWRSERT(NCSTMP,3)) )
      ENDIF
    ENDDO
  ENDIF
  RETURN

END

SUBROUTINE CALBAL4

  ! CHANGE RECORD
  ! **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,
  ! **  AND ENERGY BALANCES

  USE GLOBAL
  IMPLICIT NONE                                                                                                          
  INTEGER :: L,LN,K                                                                                                        
  REAL :: DUTMP,DVTMP                                                                                                      

  ! **  CALCULATE MOMENTUM AND ENERGY DISSIPATION                                                                         
  DO L=2,LA
    LN=LNC(L)
    UUEOUT=UUEOUT+0.5*SPB(L)*DXYP(L)*(U(L,KSZ(L))*TBX(L)+U(LEC(L),KSZ(L))*TBX(LEC(L))-U(L,KC)*TSX(L)-U(LEC(L),KC)*TSX(LEC(L)))
    VVEOUT=VVEOUT+0.5*SPB(L)*DXYP(L)*(V(L,KSZ(L))*TBY(L)+V(LN,KSZ(L))*TBX(LN)-V(L,KC)*TSY(L)-V(LN,KC)*TSX(LN))
  ENDDO
  DO K=1,KS
    DO L=2,LA
      LN=LNC(L)
      DUTMP=0.5*( U(L,K+1)+U(LEC(L),K+1)-U(L,K)-U(LEC(L),K) )
      DVTMP=0.5*( V(L,K+1)+V(LN,K+1)-V(L,K)-V(LN,K) )
      UUEOUT=UUEOUT+SPB(L)*2.0*DXYP(L)*AV(L,K)*( DUTMP*DUTMP )/(DZC(L,K+1)+DZC(L,K))
      VVEOUT=VVEOUT+SPB(L)*2.0*DXYP(L)*AV(L,K)*( DVTMP*DVTMP )/(DZC(L,K+1)+DZC(L,K))
      BBEOUT=BBEOUT+SCB(L)*DXYP(L)*HP(L)*GP*AB(L,K)*(B(L,K+1)-B(L,K))
    ENDDO
  ENDDO
  RETURN
END


SUBROUTINE CALBAL5
  !
  ! CHANGE RECORD
  ! **  SUBROUTINES CALBAL CALCULATE GLOBAL VOLUME, MASS, MOMENTUM,
  ! **  AND ENERGY BALANCES
  !
  USE GLOBAL
  IMPLICIT NONE                                                                                                          
  REAL :: ENEEND,ENEOUT,VOLBMO,SALBMO,DYEBMO,UMOBMO,VMOBMO,ENEBMO                                                          
  REAL :: VOLERR,SALERR,DYEERR,UMOERR,VMOERR,ENEERR,RVERDE,RSERDE                                                          
  REAL :: RDERDE,RUERDE,REERDE,RVERDO,RSERDO,RUERDO,ENEBEG,RDERDO                                                          
  REAL :: REERDO,RUMERDE,RVMERDE,RUMERDO,UUEBMO,VVEBMO,PPEBMO                                                              
  REAL :: BBEBMO,RVMERDO,PPEEND,BBEEND,AMOEND,VOLEND,SALEND,DYEEND                                                         
  REAL :: VMOEND,UUEEND,VVEEND,UMOEND                                                                                      
  INTEGER :: K,L,LN                                                                                                        

  ! **  CHECK FOR END OF BALANCE PERIOD                                                                                   
  IF( NBAL == NTSMMT )THEN

    ! **  CALCULATE ENDING VOLUME, SALT MASS, DYE MASS, MOMENTUM, KINETIC                                                   
    ! **  ENERGY AND POTENTIAL ENERGY, AND ASSOCIATED FLUXES                                                                
    VOLEND=0.
    SALEND=0.
    DYEEND=0.
    UMOEND=0.
    VMOEND=0.
    UUEEND=0.
    VVEEND=0.
    PPEEND=0.
    BBEEND=0.
    DO L=2,LA
      LN=LNC(L)
      VOLEND=VOLEND+SPB(L)*DXYP(L)*HP(L)
     UMOEND=UMOEND+SPB(L)*0.5*DXYP(L)*HP(L)*(DYIU(L)*HUI(L)*UHDYE(L)+DYIU(LEC(L))*HUI(LEC(L))*UHDYE(LEC(L)))
     VMOEND=VMOEND+SPB(L)*0.5*DXYP(L)*HP(L)*(DXIV(L)*HVI(L)*VHDXE(L)+DXIV(LN)*HVI(LN)*VHDXE(LN))
      PPEEND=PPEEND+SPB(L)*0.5*DXYP(L)*(GI*P(L)*P(L)-G*BELV(L)*BELV(L))
    ENDDO
    AMOEND=SQRT(UMOEND*UMOEND+VMOEND*VMOEND)
    DO K=1,KC
      DO L=2,LA
        LN=LNC(L)
        SALEND=SALEND+SCB(L)*DXYP(L)*HP(L)*SAL(L,K)*DZC(L,K)
        DYEEND=DYEEND+SCB(L)*DXYP(L)*HP(L)*DYE(L,K)*DZC(L,K)
        UUEEND=UUEEND+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K)*( (U(L,K)+U(LEC(L),K))*(U(L,K)+U(LEC(L),K)) )
        VVEEND=VVEEND+SPB(L)*0.125*DXYP(L)*HP(L)*DZC(L,K)*( (V(L,K)+V(LN,K))*(V(L,K)+V(LN,K)) )
        BBEEND=BBEEND+SPB(L)*GP*DXYP(L)*HP(L)*DZC(L,K)*( BELV(L)+0.5*HP(L)*(Z(L,K)+Z(L,K-1)) )*B(L,K)
      ENDDO
    ENDDO
    UUEOUT=DT*UUEOUT
    VVEOUT=DT*VVEOUT
    PPEOUT=DT*PPEOUT
    BBEOUT=DT*BBEOUT
    VOLOUT=DT*VOLOUT
    SALOUT=DT*SALOUT
    DYEOUT=DT*DYEOUT
    UMOOUT=DT*UMOOUT
    VMOOUT=DT*VMOOUT
    ENEBEG=UUEBEG+VVEBEG+PPEBEG+BBEBEG
    ENEEND=UUEEND+VVEEND+PPEEND+BBEEND
    ENEOUT=UUEOUT+VVEOUT+PPEOUT+BBEOUT
    VOLBMO=VOLBEG-VOLOUT
    SALBMO=SALBEG-SALOUT
    DYEBMO=DYEBEG-DYEOUT
    UMOBMO=UMOBEG-DYEOUT
    VMOBMO=VMOBEG-DYEOUT
    ENEBMO=ENEBEG-ENEOUT
    VOLERR=VOLEND-VOLBMO
    SALERR=SALEND-SALBMO
    DYEERR=DYEEND-DYEBMO
    UMOERR=UMOEND-UMOBMO
    VMOERR=VMOEND-VMOBMO
    ENEERR=ENEEND-ENEBMO
    RVERDE=-9999.
    RSERDE=-9999.
    RDERDE=-9999.
    RUERDE=-9999.
    RVERDE=-9999.
    REERDE=-9999.
    RVERDO=-9999.
    RSERDO=-9999.
    RDERDO=-9999.
    RUERDO=-9999.
    RVERDO=-9999.
    REERDO=-9999.
    IF( VOLEND /= 0. ) RVERDE=VOLERR/VOLEND
    IF( SALEND /= 0. ) RSERDE=SALERR/SALEND
    IF( DYEEND /= 0. ) RDERDE=DYEERR/DYEEND
    IF( UMOEND /= 0. ) RUMERDE=UMOERR/UMOEND
    IF( VMOEND /= 0. ) RVMERDE=VMOERR/VMOEND
    IF( ENEEND /= 0. ) REERDE=ENEERR/ENEEND
    IF( VOLOUT /= 0. ) RVERDO=VOLERR/VOLOUT
    IF( SALOUT /= 0. ) RSERDO=SALERR/SALOUT
    IF( DYEOUT /= 0. ) RDERDO=DYEERR/DYEOUT
    IF( UMOOUT /= 0. ) RUMERDO=UMOERR/UMOOUT
    IF( VMOOUT /= 0. ) RVMERDO=VMOERR/VMOOUT
    IF( ENEOUT /= 0. ) REERDO=ENEERR/ENEOUT

    ! **  OUTPUT BALANCE RESULTS TO FILE BAL.OUT
    IF( JSBAL == 1 )THEN
      OPEN(89,FILE=OUTDIR//'BAL.OUT',STATUS='UNKNOWN')
      CLOSE(89,STATUS='DELETE')
      OPEN(89,FILE=OUTDIR//'BAL.OUT',STATUS='UNKNOWN')
      JSBAL=0
    ELSE
      OPEN(89,FILE=OUTDIR//'BAL.OUT',POSITION='APPEND',STATUS='UNKNOWN')
    ENDIF
    WRITE(89,890)NTSMMT,N
    WRITE(89,891)
    WRITE(89,892)VOLBEG,SALBEG,DYEBEG,ENEBEG,UMOBEG,VMOBEG,AMOBEG
    WRITE(89,900)
    WRITE(89,893)
    WRITE(89,892)VOLOUT,SALOUT,DYEOUT,ENEOUT,UMOOUT,VMOOUT
    WRITE(89,900)
    WRITE(89,894)
    WRITE(89,892)VOLBMO,SALBMO,DYEBMO,ENEBMO,UMOBMO,VMOBMO
    WRITE(89,900)
    WRITE(89,895)
    WRITE(89,892)VOLEND,SALEND,DYEEND,ENEEND,UMOEND,VMOEND,AMOEND
    WRITE(89,900)
    WRITE(89,896)
    WRITE(89,892)VOLERR,SALERR,DYEERR,ENEERR,UMOERR,VMOERR
    WRITE(89,900)
    WRITE(89,897)
    WRITE(89,892)RVERDE,RSERDE,RDERDE,REERDE,RUMERDE,RVMERDE
    WRITE(89,900)
    WRITE(89,898)
    WRITE(89,892)RVERDO,RSERDO,RDERDO,REERDO,RUMERDO,RVMERDO
    WRITE(89,899)
    UUEBMO=UUEBEG-UUEOUT
    VVEBMO=VVEBEG-VVEOUT
    PPEBMO=PPEBEG-PPEOUT
    BBEBMO=BBEBEG-BBEOUT
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
    WRITE(89,900)
    WRITE(89,899)
    CLOSE(89)
    890 FORMAT (' VOLUME, MASS, AND ENERGY BALANCE OVER',I5,' TIME STEPS',' ENDING AT TIME STEP',I5,//)
    891 FORMAT (' INITIAL VOLUME    INITIAL SALT    INITIAL DYE     INITIAL ENER    INITIAL UMO     INITIAL VMO     INITIAL AMO',/)
    892 FORMAT (1X,7(E14.6,2X))
    893 FORMAT (' VOLUME OUT        SALT OUT        DYE OUT         ENERGY OUT      UMO OUT         VMO OUT',/)
    894 FORMAT (' INITIAL-OUT VOL   INIT-OUT SALT   INIT-OUT DYE    INIT-OUT ENER   INIT-OUT UMO    INIT-OUT VMO',/)
    895 FORMAT (' FINAL VOLUME      FINAL SALT      FINAL DYE       FINAL ENERGY    FINAL UMO       FINAL VMO       FINAL AMO',/)
    896 FORMAT (' VOLUME ERR        SALT ERR        DYE ERR         ENERGY ERR      UMO ERR         VMO ERR',/)
    897 FORMAT (' R VOL/END ER      R SAL/END ER    R DYE/END ER    R ENE/END ER    R UMO/END ER    R VMO/END ER',/)
    898 FORMAT (' R VOL/OUT ER      R SAL/OUT ER    R DYE/OUT ER    R ENE/OUT ER    R UMO/OUT ER    R VMO/OUT ER',/)
    899 FORMAT (////)
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
    NBAL=0
  ENDIF
  NBAL=NBAL+1
  RETURN
  
END

