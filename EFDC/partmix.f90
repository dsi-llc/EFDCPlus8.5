SUBROUTINE PARTMIX(NT)

  !**********************************************************************C
  !
  ! **  SUBROUTINE PARTMIX CALCULATES PARTICLE MIXING OF TOXICS IN THE 
  ! **  TOP LAYER(S) [PMXDEPTH] OF THE SEDIMENT BED
  !
  !**********************************************************************C

  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2015-12           Paul M. Craig      ADDED OMP AND CLEANED UP CODE

  USE GLOBAL

  IMPLICIT NONE
  
  INTEGER :: ND,LF,LL,LP,K,L,NT,LZ,KM,NPMXPTS,NDP,KK,NP
  REAL    :: DEPINBED,WT1,WT2,TMPVAL,DELHBED,TERM1,TERM2,TERM3
  REAL    :: PARTMIXAVG(KBM)

  !**********************************************************************
  !$OMP PARALLEL DEFAULT(SHARED)
  
  !$OMP DO PRIVATE(ND,LF,LL,K,LP,L)
  DO ND=1,NDM  
    LF=(ND-1)*LDMWET+1  
    LL=MIN(LF+LDMWET-1,LAWET)

    DO K=1,KB
      DO LP=LF,LL
        L=LWET(LP)  
        PARTMIXZ(L,K)=0.0
      ENDDO
    ENDDO
  ENDDO   ! *** END OF DOMAIN 
  !$OMP END DO
  
  !----------------------------------------------------------------------
  !  OLD PARTICLE MIXING IS NOW OPTOIN ISPMXZ=2
  IF( ISPMXZ(NT) == 2 )THEN

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LZ,K,KM,NP,NDP,DEPINBED,WT1,WT2,TMPVAL)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        DEPINBED=0.
        LZ=LPMXZ(L)
        
        IF( KBT(L) > 2 )THEN
          DO K=KBT(L),2,-1
            KM = K-1
            DEPINBED=DEPINBED+HBED(L,K)
            DO NP=1,NPMXPTS-1
              NDP=NP+1
              IF( DEPINBED >= PMXDEPTH(NP,LZ) .AND. DEPINBED < PMXDEPTH(NDP,LZ) )THEN
                WT1 = DEPINBED-PMXDEPTH(NP,LZ)
                WT2 = PMXDEPTH(NDP,LZ)-DEPINBED
                TMPVAL = PMXDEPTH(NP+1,LZ)-PMXDEPTH(NP,LZ)
                PARTMIXZ(L,KM) = (WT2*PMXCOEF(NP,LZ) + WT1*PMXCOEF(NDP,LZ))/TMPVAL
              ENDIF
            ENDDO
          ENDDO
        ELSE
          K  = KBT(L)
          KM = K-1
          DEPINBED = DEPINBED + HBED(L,K)
          DO NP=1,NPMXPTS-1
            NDP=NP+1
            IF( DEPINBED >= PMXDEPTH(NP,LZ) .AND. DEPINBED < PMXDEPTH(NDP,LZ) )THEN
              WT1 = DEPINBED-PMXDEPTH(NP,LZ)
              WT2 = PMXDEPTH(NDP,LZ)-DEPINBED
              TMPVAL = PMXDEPTH(NP+1,LZ)-PMXDEPTH(NP,LZ)
              PARTMIXZ(L,KM) = (WT2*PMXCOEF(NP,LZ) + WT1*PMXCOEF(NDP,LZ))/TMPVAL
            ENDIF
          ENDDO
        ENDIF
      ENDDO

      DO LP=LF,LL
        L=LWET(LP)  
        DO K=KBT(L),2,-1
          KM = K-1
          PARTMIXZ(L,KM) = 2.*PARTMIXZ(L,KM)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT)+1.E-12)
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO
  
  ENDIF

  !----------------------------------------------------------------------
  !  NEW PARTICLE MIXING IS NOW OPTOIN ISPMXZ=1
  IF( ISPMXZ(NT) == 1 )THEN

    !$OMP DO PRIVATE(ND,LF,LL,LP,L,LZ,K,KK,KM,NP,NDP,DEPINBED,DELHBED,WT1,WT2,TMPVAL,TERM1,TERM2,TERM3)
    DO ND=1,NDM  
      LF=(ND-1)*LDMWET+1  
      LL=MIN(LF+LDMWET-1,LAWET)

      DO LP=LF,LL
        L=LWET(LP)  
        LZ=LPMXZ(L)

        DO K=1,KB
          PARTMIXAVG(K)=0.0
          PARTMIXZ(L,K)=0.0
        ENDDO

        DEPINBED=0.
        DO K=KBT(L),1,-1
          DELHBED=HBED(L,K)/10.
          DO KK=1,10
            DEPINBED=DEPINBED+DELHBED
            DO NP=1,NPMXPTS-1
              NDP=NP+1
              IF( DEPINBED >= PMXDEPTH(NP,LZ) .AND. DEPINBED < PMXDEPTH(NDP,LZ) )THEN
                WT1 = DEPINBED-PMXDEPTH(NP,LZ)
                WT2 = PMXDEPTH(NDP,LZ)-DEPINBED
                TMPVAL = PMXDEPTH(NP+1,LZ)-PMXDEPTH(NP,LZ)
                PARTMIXAVG(K)=PARTMIXAVG(K) + (WT2*PMXCOEF(NP,LZ)+WT1*PMXCOEF(NDP,LZ))/TMPVAL
              ENDIF
            ENDDO
          ENDDO
          PARTMIXAVG(K)=PARTMIXAVG(K)/10.
        ENDDO

        DO K=1,KBT(L)-1
          TERM1 = HBED(L,K)/PARTMIXAVG(K)
          TERM2 = HBED(L,K+1)/PARTMIXAVG(K+1)
          TERM3 = HBED(L,K)+HBED(L,K+1)
          PARTMIXZ(L,K) = TERM3/(TERM1+TERM2)
        ENDDO

      ENDDO

      DO LP=LF,LL
        L=LWET(LP)  
        DO K=KBT(L),2,-1
          KM = K-1
          PARTMIXZ(L,KM) = 2.*PARTMIXZ(L,KM)/(2.0-TOXPFTB(L,K,NT)-TOXPFTB(L,KM,NT)+1.E-12)
        ENDDO
      ENDDO
    ENDDO   ! *** END OF DOMAIN LOOP
    !$OMP END DO

  ENDIF

  !**********************************************************************C\
  !$OMP END PARALLEL
  
  RETURN

END

