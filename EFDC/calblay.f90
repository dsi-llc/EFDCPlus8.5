SUBROUTINE CALBLAY  

  !**********************************************************************!
  ! *** SUBROUTINE CALBLAY REMOVES OR ADDS LAYERS TO THE SEDIMENT BED
  ! *** NOT USED FOR SEDZLJ
  !  
  !**********************************************************************!
  !
  !----------------------------------------------------------------------!
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2011-03           Paul M. Craig      Rewritten to F90 and added OMP
  !----------------------------------------------------------------------!

  USE GLOBAL  
  IMPLICIT NONE

  INTEGER :: K,NS,L,NT,NX,ND,LF,LL,IFOUND
  REAL   :: TMPBOT2,TMPTOP1,TMPTOP2,TMPVAL,HBEDMXT,HOLDTOP,FKBTP
  REAL   :: SEDBOLD,TOXBOLD,TMPBOT1,FKBT,SNDBOLD,HTMP1,HTMP2,HBEDMAX1

  ! *** FOR TRANSPORT OF COHESIVE SEDIMENT ONLY SET HBEDMIN TO FRACTION
  ! *** OF HBEDMAX
  IF( ISTRAN(7) == 0 )THEN
      HBEDMIN = 0.1*HBEDMAX
  ELSE
    ! *** NOTE: BELOW IS USED BY HQI FOR HOUSATONIC
    HBEDMIN = 0.006
  ENDIF

  ! *** WHEN NONCOHESIVE TRANSPORT IS ACTIVE, WITHOUT ACTIVE-PARENT LAYER
  ! *** FORMULATION, SET HBEDMIN PROPORTIONAL TO MAXIMUM GRAIN DIAMETER
  IF( ISTRAN(7) >= 1 .AND. ISNDAL <= 1 )THEN
    ! *** NOTE HQI CHANGED ORIGINAL FORMULATION FROM
    ! *** TMPVAL = 2.*SNDDMX
    ! *** TO
    TMPVAL = SNDDMX
    HBEDMIN = MAX(HBEDMIN,TMPVAL)
  END IF

  ! *** WHEN NONCOHESIVE TRANSPORT IS ACTIVE, WITH ACTIVE-PARENT LAYER
  ! *** FORMULATION, SET HBEDMIN PROPORTIONAL ACTIVE LAYER THICKNESS
  IF( ISNDAL == 2 )THEN
    HBEDMIN = 2.*HBEDAL  
    HBEDMAX1 = 1.1*HBEDAL
  END IF

  HBEDMXT = HBEDMAX+HBEDMIN  

  ! *** ADD OR REMOVE TOP LAYER (NO ACTIVE LAYER ARMORING OPTION)  
  IF( ISNDAL <= 1 )THEN
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,L,LF,LL,NS,NX,NT)              &
    !$OMP             PRIVATE(TMPBOT1,TMPBOT2,TMPTOP1,TMPTOP2,TMPVAL,HOLDTOP) &
    !$OMP             PRIVATE(SEDBOLD,SNDBOLD,TOXBOLD,FKBT,FKBTP,HTMP1,HTMP2)
    DO ND=1,NDM  
      LF = 2+(ND-1)*LDM  
      LL = MIN(LF+LDM-1,LA)

      ! *** ADD NEW TOP LAYER  
      DO L=LF,LL  
        IF( HBED(L,KBT(L)) > HBEDMXT )THEN  
          IF( KBT(L) < KB )THEN  
            HOLDTOP = HBED(L,KBT(L))
            HTMP1 = HOLDTOP-HBEDMAX
            HTMP2 = HBEDMAX
            HBED(L,KBT(L)+1) = MIN(HTMP1,HTMP2)   ! *** Thinner layer on top
            HBED(L,KBT(L)) = MAX(HTMP1,HTMP2)
            IF( IBMECH == 1 .AND. SEDVRDT < 0.0 )THEN   ! *** These control maintianing initial void ratio k profile 
              VDRBED2(L,KBT(L)+1) = VDRBED2(L,KBT(L))     ! *** Reset void ratio profile
              SDENAVG(L,KBT(L)+1) = SDENAVG(L,KBT(L))     ! *** Reset avg sediment solids density
            ENDIF                                    
            VDRBED(L,KBT(L)+1) = VDRBED(L,KBT(L))  
            PORBED(L,KBT(L)+1) = PORBED(L,KBT(L))  
            STDOCB(L,KBT(L)+1) = STDOCB(L,KBT(L))  
            STPOCB(L,KBT(L)+1) = STPOCB(L,KBT(L))  
            FKBTP = HBED(L,KBT(L)+1)/HOLDTOP  
            FKBT = HBED(L,KBT(L))/HOLDTOP
              
            IF( ISTRAN(6) >= 1 )THEN  
              SEDBT(L,KBT(L)+1) = 0.  
              SEDBT(L,KBT(L)) = 0.  
              DO NS=1,NSED  
                SEDBOLD = SEDB(L,KBT(L),NS)  
                SEDB(L,KBT(L)+1,NS) = FKBTP*SEDBOLD  
                SEDB(L,KBT(L),NS) = FKBT*SEDBOLD  
                SEDBT(L,KBT(L)+1) = SEDBT(L,KBT(L)+1)+SEDB(L,KBT(L)+1,NS)  
                SEDBT(L,KBT(L)) = SEDBT(L,KBT(L))+SEDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NS) = STFPOCB(L,KBT(L),NS)  
                VFRBED(L,KBT(L)+1,NS) = VFRBED(L,KBT(L),NS)  
              ENDDO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              SNDBT(L,KBT(L)+1) = 0.  
              SNDBT(L,KBT(L)) = 0.  
              DO NS=1,NSND  
                NX = NS+NSED  
                SNDBOLD = SNDB(L,KBT(L),NS)  
                SNDB(L,KBT(L)+1,NS) = FKBTP*SNDBOLD  
                SNDB(L,KBT(L),NS) = FKBT*SNDBOLD  
                SNDBT(L,KBT(L)+1) = SNDBT(L,KBT(L)+1)+SNDB(L,KBT(L)+1,NS)  
                SNDBT(L,KBT(L)) = SNDBT(L,KBT(L))+SNDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NX) = STFPOCB(L,KBT(L),NX)  
                VFRBED(L,KBT(L)+1,NX) = VFRBED(L,KBT(L),NX)  
              ENDDO  
            ENDIF  
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TOXBOLD = TOXB(L,KBT(L),NT)  
                TOXB(L,KBT(L)+1,NT) = FKBTP*TOXBOLD  
                TOXB(L,KBT(L),NT) = FKBT*TOXBOLD  
              ENDDO  
            ENDIF  
            KBT(L) = KBT(L)+1  
          ENDIF  
        ENDIF  
      ENDDO  

      ! *** REZONE WITH NEW TOP LAYER ADDED NEXT TIME STEP  
      DO L=LF,LL  
        IF( HBED(L,KBT(L)) > HBEDMXT )THEN  
          IF( KBT(L) == KB .AND. KB > 1 )THEN  
            TMPBOT1 = HBED(L,1)/(1.+VDRBED(L,1))  
            TMPBOT2 = HBED(L,2)/(1.+VDRBED(L,2))  
            TMPTOP1 = VDRBED(L,1)*TMPBOT1  
            TMPTOP2 = VDRBED(L,2)*TMPBOT2  
            VDRBED(L,1) = (TMPTOP1+TMPTOP2)/(TMPBOT1+TMPBOT2)  
            PORBED(L,1) = VDRBED(L,1)/(1.+VDRBED(L,1))  
            HBED(L,1) = HBED(L,1)+HBED(L,2)  
            IF( KB == 2 )THEN  
              HBED(L,2) = 0  
              VDRBED(L,2) = 0.0  
              PORBED(L,2) = 0.0  
              STDOCB(L,2) = 0.0  
              STPOCB(L,2) = 0.0  
            ENDIF  
            IF( KB > 2 )THEN  
              DO K=2,KBT(L)-1  
                HBED(L,K) = HBED(L,K+1)  
                VDRBED(L,K) = VDRBED(L,K+1)  
                PORBED(L,K) = PORBED(L,K+1)  
                STDOCB(L,K) = STDOCB(L,K+1)  
                STPOCB(L,K) = STPOCB(L,K+1)  
              ENDDO  
              HBED(L,KBT(L)) = 0  
              VDRBED(L,KBT(L)) = 0.0  
              PORBED(L,KBT(L)) = 0.0  
              STDOCB(L,KBT(L)) = 0.0  
              STPOCB(L,KBT(L)) = 0.0  
            ENDIF  
            IF( ISTRAN(6) >= 1 )THEN  
              DO NS=1,NSED  
                SEDB(L,1,NS) = SEDB(L,1,NS)+SEDB(L,2,NS)  
                IF( KB == 2 )THEN  
                  SEDB(L,2,NS) = 0.0  
                  STFPOCB(L,2,NS) = 0.0  
                  VFRBED(L,2,NS) = 0.0  
                ENDIF  
                IF( KB > 2 )THEN  
                  DO K=2,KBT(L)-1  
                    SEDB(L,K,NS) = SEDB(L,K+1,NS)  
                    STFPOCB(L,K,NS) = STFPOCB(L,K+1,NS)  
                    VFRBED(L,K,NS) = VFRBED(L,K+1,NS)  
                  ENDDO  
                  SEDB(L,KBT(L),NS) = 0  
                  STFPOCB(L,KBT(L),NS) = 0.0  
                  VFRBED(L,KBT(L),NS) = 0.0  
                ENDIF  
              ENDDO  
              DO K=1,KB  
                SEDBT(L,K) = 0.0  
              END DO  
              DO NS=1,NSED  
                DO K=1,KB  
                  SEDBT(L,K) = SEDBT(L,K)+SEDB(L,K,NS)  
                END DO  
              END DO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              DO NS=1,NSND  
                NX = NS+NSED  
                SNDB(L,1,NS) = SNDB(L,1,NS)+SNDB(L,2,NS)  
                IF( KB == 2 )THEN  
                  SNDB(L,2,NS) = 0.0  
                  STFPOCB(L,2,NX) = 0.0  
                  VFRBED(L,2,NX) = 0.0  
                ENDIF  
                IF( KB > 2 )THEN  
                  DO K=2,KBT(L)-1  
                    SNDB(L,K,NS) = SNDB(L,K+1,NS)  
                    STFPOCB(L,K,NX) = STFPOCB(L,K+1,NX)  
                    VFRBED(L,K,NX) = VFRBED(L,K+1,NX)  
                  ENDDO  
                  SNDB(L,KBT(L),NS) = 0  
                  STFPOCB(L,KBT(L),NX) = 0.0  
                  VFRBED(L,KBT(L),NX) = 0.0  
                ENDIF  
              ENDDO  
              DO K=1,KB  
                SNDBT(L,K) = 0.0  
              END DO  
              DO NS=1,NSND  
                DO K=1,KB  
                  SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NS)  
                END DO  
              END DO  
            ENDIF  
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TOXB(L,1,NT) = TOXB(L,1,NT)+TOXB(L,2,NT)  
                IF( KB == 2 )THEN  
                  TOXB(L,2,NT) = 0  
                ENDIF  
                IF( KB > 2 )THEN  
                  DO K=2,KBT(L)-1  
                    TOXB(L,K,NT) = TOXB(L,K+1,NT)  
                  ENDDO  
                  TOXB(L,KBT(L),NT) = 0  
                ENDIF  
              ENDDO  
            ENDIF  
            IF( KB == 2 )THEN  
              KBT(L) = 1  
            ENDIF  
            IF( KB > 2 )THEN  
              KBT(L) = KBT(L)-1  
            ENDIF  
          ENDIF  
        ENDIF  
      ENDDO  

      ! *** REMOVE TOP LAYER  
      DO L=LF,LL  
        IF( HBED(L,KBT(L)) < HBEDMIN )THEN  
          IF( KBT(L) > 1 )THEN  
            TMPBOT1 = HBED(L,KBT(L)-1)/(1.+VDRBED(L,KBT(L)-1))  
            TMPBOT2 = HBED(L,KBT(L))/(1.+VDRBED(L,KBT(L)))  
            TMPTOP1 = VDRBED(L,KBT(L)-1)*TMPBOT1  
            TMPTOP2 = VDRBED(L,KBT(L))*TMPBOT2  
            VDRBED(L,KBT(L)-1) = (TMPTOP1+TMPTOP2)/(TMPBOT1+TMPBOT2)  
            PORBED(L,KBT(L)-1) = VDRBED(L,KBT(L)-1)/(1.+VDRBED(L,KBT(L)-1))  
            HBED(L,KBT(L)-1) = HBED(L,KBT(L)-1)+HBED(L,KBT(L))  
            HBED(L,KBT(L)) = 0.0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  
            BDENBED(L,KBT(L)) = 0.  ! PMC
            IF( ISTRAN(6) >= 1 )THEN  
              SEDBT(L,KBT(L)-1) = 0.  
              SEDBT(L,KBT(L)) = 0.  
              DO NS=1,NSED  
                SEDB(L,KBT(L)-1,NS) = SEDB(L,KBT(L)-1,NS)+SEDB(L,KBT(L),NS)  
                SEDBT(L,KBT(L)-1) = SEDBT(L,KBT(L)-1)    +SEDB(L,KBT(L)-1,NS)  
                SEDB(L,KBT(L),NS) = 0.0  
                STFPOCB(L,KBT(L),NS) = 0.0  
                VFRBED(L,KBT(L),NS) = 0.0  
              ENDDO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              SNDBT(L,KBT(L)-1) = 0.  
              SNDBT(L,KBT(L)) = 0.  
              DO NS=1,NSND  
                SNDB(L,KBT(L)-1,NS) = SNDB(L,KBT(L)-1,NS)+SNDB(L,KBT(L),NS)  
                SNDBT(L,KBT(L)-1) = SNDBT(L,KBT(L)-1) +SNDB(L,KBT(L)-1,NS)  
                SNDB(L,KBT(L),NS) = 0.0  
                STFPOCB(L,KBT(L),NS+NSED) = 0.0  
                VFRBED(L,KBT(L),NS+NSED) = 0.0  
              ENDDO  
            ENDIF  
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TOXB(L,KBT(L)-1,NT) = TOXB(L,KBT(L)-1,NT)+TOXB(L,KBT(L),NT)  
                TOXB(L,KBT(L),NT) = 0.0  
              ENDDO  
            ENDIF  
            KBT(L) = KBT(L)-1  

          ! *** PMC BEGIN BLOCK  
          ELSEIF( HBED(L,KBT(L)) < 0.0 )THEN  
            ! *** ZERO NEGATIVE THICKNESSES
            HBED(L,KBT(L)) = 0.0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  
          ENDIF  
          ! *** PMC END BLOCK  
        ENDIF  
      ENDDO  

      ! *** UPDATE BULK DENSITY  
      DO L=LF,LL  
        ! ***REMOVED KB LOOP, ONLY COMPUTE THE TOP LAYER.  CONSOLIDATION IS HANDLED LATER IN CALBED
        K = KBT(L)  
        IF( HBED(L,K) > 0. )THEN  
          ! *** UPDATE TOTAL/WET DENSITY
          BDENBED(L,K) = 1000.*PORBED(L,K)+0.001*(SEDBT(L,K)+SNDBT(L,K))/HBED(L,K)  
        ELSE  
          BDENBED(L,K) = 0.  
        ENDIF  
      ENDDO  

    ENDDO  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO
         
    ! ***************************************************************************        
  ELSE  
    ! ***************************************************************************        
    ! *** ADD OR REMOVE PARENT LAYER WHEN ARMORING IS ACTIVE  (ISNDAL == 2 )
  
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ND,K,L,LF,LL,NS,NX,NT,IFOUND)    &
    !$OMP             PRIVATE(TMPBOT1,TMPBOT2,TMPTOP1,TMPTOP2,TMPVAL,HOLDTOP)  &
    !$OMP             PRIVATE(SEDBOLD,TOXBOLD,FKBT,FKBTP,SNDBOLD,HTMP1,HTMP2)
    DO ND=1,NDM  
      LF = 2+(ND-1)*LDM  
      LL = MIN(LF+LDM-1,LA)

      ! *** CHECK TO SEE IF NEED TO ADD NEW LAYER BELOW THE ACTIVE LAYER 
      DO L=LF,LL
        K = KBT(L)-1
        IF( K == 0 )CYCLE             ! *** CYCLE IF NO PARENT LAYER
        IF( HBED(L,K) > HBEDMXT )THEN  
          IF( KBT(L) < KB )THEN  
            ! *** MOVE ACTIVE LAYER UP  
            HBED(L,KBT(L)+1) = HBED(L,KBT(L))  
            VDRBED(L,KBT(L)+1) = VDRBED(L,KBT(L))  
            PORBED(L,KBT(L)+1) = PORBED(L,KBT(L))  
            STDOCB(L,KBT(L)+1) = STDOCB(L,KBT(L))  
            STPOCB(L,KBT(L)+1) = STPOCB(L,KBT(L))  
            
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TOXB(L,KBT(L)+1,NT) = TOXB(L,KBT(L),NT)  
              ENDDO  
            ENDIF  
            IF( ISTRAN(6) >= 1 )THEN  
              SEDBT(L,KBT(L)+1) = SEDBT(L,KBT(L))  
              DO NS=1,NSED  
                SEDB(L,KBT(L)+1,NS) = SEDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NS) = STFPOCB(L,KBT(L),NS)  
                VFRBED(L,KBT(L)+1,NS) = VFRBED(L,KBT(L),NS)  
              ENDDO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              SNDBT(L,KBT(L)+1) = SNDBT(L,KBT(L))  
              DO NS=1,NSND  
                NX = NS+NSED  
                SNDB(L,KBT(L)+1,NS) = SNDB(L,KBT(L),NS)  
                STFPOCB(L,KBT(L)+1,NX) = STFPOCB(L,KBT(L),NX)  
                VFRBED(L,KBT(L)+1,NX) = VFRBED(L,KBT(L),NX)  
              ENDDO  
            ENDIF  

            ! *** SPLIT PARENT INTO TWO LAYERS  
            HOLDTOP = HBED(L,KBT(L)-1)  
            HBED(L,KBT(L)) = HOLDTOP-HBEDMAX  
            HBED(L,KBT(L)-1) = HBEDMAX  
            VDRBED(L,KBT(L)) = VDRBED(L,KBT(L)-1)  
            PORBED(L,KBT(L)) = PORBED(L,KBT(L)-1)  
            STDOCB(L,KBT(L)) = STDOCB(L,KBT(L)-1)  
            STPOCB(L,KBT(L)) = STPOCB(L,KBT(L)-1)  
            FKBTP = HBED(L,KBT(L))/HOLDTOP  
            FKBT = HBED(L,KBT(L)-1)/HOLDTOP  
            
            IF( ISTRAN(6) >= 1 )THEN  
              SEDBT(L,KBT(L)) = 0.  
              SEDBT(L,KBT(L)-1) = 0.  
              DO NS=1,NSED  
                SEDBOLD = SEDB(L,KBT(L)-1,NS)  
                SEDB(L,KBT(L),NS) = FKBTP*SEDBOLD  
                SEDB(L,KBT(L)-1,NS) = FKBT*SEDBOLD  
                SEDBT(L,KBT(L)) = SEDBT(L,KBT(L))  +SEDB(L,KBT(L),NS)  
                SEDBT(L,KBT(L)-1) = SEDBT(L,KBT(L)-1)  +SEDB(L,KBT(L)-1,NS)  
                STFPOCB(L,KBT(L),NS) = STFPOCB(L,KBT(L)-1,NS)  
                VFRBED(L,KBT(L),NS) = VFRBED(L,KBT(L)-1,NS)  
              ENDDO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              SNDBT(L,KBT(L)) = 0.  
              SNDBT(L,KBT(L)-1) = 0.  
              DO NS=1,NSND  
                NX = NS+NSED  
                SNDBOLD = SNDB(L,KBT(L)-1,NS)  
                SNDB(L,KBT(L),NS) = FKBTP*SNDBOLD  
                SNDB(L,KBT(L)-1,NS) = FKBT*SNDBOLD  
                SNDBT(L,KBT(L)) = SNDBT(L,KBT(L)) +SNDB(L,KBT(L),NS)  
                SNDBT(L,KBT(L)-1) = SNDBT(L,KBT(L)-1) +SNDB(L,KBT(L)-1,NS)  
                STFPOCB(L,KBT(L),NX) = STFPOCB(L,KBT(L)-1,NX)  
                VFRBED(L,KBT(L),NX) = VFRBED(L,KBT(L)-1,NX)  
              ENDDO  
            ENDIF  
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TOXBOLD = TOXB(L,KBT(L)-1,NT)  
                TOXB(L,KBT(L),NT) = FKBTP*TOXBOLD  
                TOXB(L,KBT(L)-1,NT) = FKBT*TOXBOLD  
              ENDDO  
            ENDIF  
            KBT(L) = KBT(L)+1
          
          ELSEIF( KB > 2 )THEN
            ! *** KBT(L) = KB,  REZONE LAYERS BELOW ACTIVE LAYER
            
            ! *** COMBINE THE BOTTOM TWO LAYERS
            TMPBOT1 = HBED(L,1)/(1.+VDRBED(L,1))  
            TMPBOT2 = HBED(L,2)/(1.+VDRBED(L,2))  
            TMPTOP1 = VDRBED(L,1)*TMPBOT1  
            TMPTOP2 = VDRBED(L,2)*TMPBOT2  
            VDRBED(L,1) = (TMPTOP1+TMPTOP2)/(TMPBOT1+TMPBOT2)  
            PORBED(L,1) = VDRBED(L,1)/(1.+VDRBED(L,1))  
            HBED(L,1) = HBED(L,1) + HBED(L,2)
            
            ! *** MOVE ALL SEDIMENT DOWN ONE LAYER
            DO K=2,KBT(L)-1  
              HBED(L,K) = HBED(L,K+1)  
              VDRBED(L,K) = VDRBED(L,K+1)  
              PORBED(L,K) = PORBED(L,K+1)  
              STDOCB(L,K) = STDOCB(L,K+1)  
              STPOCB(L,K) = STPOCB(L,K+1)  
            ENDDO  
            
            ! *** CREATE A NEW EMPTY KB LAYER
            HBED(L,KBT(L)) = 0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  

            ! *** REASSIGN MASSES
            IF( ISTRAN(6) >= 1 )THEN  
              DO NS=1,NSED  
                SEDB(L,1,NS) = SEDB(L,1,NS)+SEDB(L,2,NS)  
                DO K=2,KBT(L)-1  
                  SEDB(L,K,NS) = SEDB(L,K+1,NS)  
                  STFPOCB(L,K,NS) = STFPOCB(L,K+1,NS)  
                  VFRBED(L,K,NS) = VFRBED(L,K+1,NS)  
                ENDDO  
                SEDB(L,KBT(L),NS) = 0.0
                STFPOCB(L,KBT(L),NS) = 0.0  
                VFRBED(L,KBT(L),NS) = 0.0  
              ENDDO  
              DO K=1,KB  
                SEDBT(L,K) = 0.0  
              END DO  
              DO NS=1,NSED  
                DO K=1,KB  
                  SEDBT(L,K) = SEDBT(L,K)+SEDB(L,K,NS)  
                END DO  
              END DO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              DO NS=1,NSND  
                NX = NS+NSED  
                SNDB(L,1,NS) = SNDB(L,1,NS)+SNDB(L,2,NS)  
                DO K=2,KBT(L)-1  
                  SNDB(L,K,NS) = SNDB(L,K+1,NS)  
                  STFPOCB(L,K,NX) = STFPOCB(L,K+1,NX)  
                  VFRBED(L,K,NX) = VFRBED(L,K+1,NX)  
                ENDDO  
                SNDB(L,KBT(L),NS) = 0.0
                STFPOCB(L,KBT(L),NX) = 0.0  
                VFRBED(L,KBT(L),NX) = 0.0  
              ENDDO  
              DO K=1,KB  
                SNDBT(L,K) = 0.0  
              END DO  
              DO NS=1,NSND  
                DO K=1,KB  
                  SNDBT(L,K) = SNDBT(L,K)+SNDB(L,K,NS)  
                END DO  
              END DO  
            ENDIF  
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TOXB(L,1,NT) = TOXB(L,1,NT)+TOXB(L,2,NT)  
                DO K=2,KBT(L)-1  
                  TOXB(L,K,NT) = TOXB(L,K+1,NT)  
                ENDDO  
                TOXB(L,KBT(L),NT) = 0  
              ENDDO  
            ENDIF  
            KBT(L) = KBT(L)-1  
          ENDIF  
        ENDIF   
      ENDDO  
      
      ! *** LIMIT LAYER THICKNESSES FOR ACTIVE LAYERS AFTER SCOUR EVENTS
      DO L=LF,LL  
        K = KBT(L)
        IF( K == 1 )THEN
          IF( HBED(L,K) > HBEDMAX1 )THEN  
            TMPBOT1 = HBEDAL/HBED(L,K)
            TMPBOT2 = 1.0 - TMPBOT1
          
            TMPTOP1 = VDRBED(L,KBT(L)-2)*TMPBOT1  
            TMPTOP2 = VDRBED(L,KBT(L)-1)*TMPBOT2
          
            VDRBED(L,K+1) = VDRBED(L,K)  
            PORBED(L,K+1) = PORBED(L,K)  
            STDOCB(L,K+1) = STDOCB(L,K)   
            STPOCB(L,K+1) = STPOCB(L,K) 
            BDENBED(L,K+1) = BDENBED(L,K) 
            HBED(L,K+1) = HBEDAL
            HBED(L,K) = HBED(L,K) - HBEDAL

            ! *** SPLIT MASS
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TMPVAL = TOXB(L,K,NT)
                TOXB(L,K+1,NT) = TMPVAL*TMPBOT1
                TOXB(L,K,NT)   = TMPVAL*TMPBOT2
              ENDDO  
            ENDIF  
            IF( ISTRAN(6) >= 1 )THEN  
              TMPVAL = SEDBT(L,K)
              SEDBT(L,K+1) = TMPVAL*TMPBOT1
              SEDBT(L,K)   = TMPVAL*TMPBOT2
              DO NS=1,NSED  
                TMPVAL = SEDB(L,K,NS)
                SEDB(L,K+1,NS) = TMPVAL*TMPBOT1
                SEDB(L,K,NS)   = TMPVAL*TMPBOT2
                STFPOCB(L,K+1,NS) = STFPOCB(L,K,NS)  
                VFRBED(L,K+1,NS) = VFRBED(L,K,NS)  
              ENDDO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              TMPVAL = SNDBT(L,K)
              SNDBT(L,K+1) = TMPVAL*TMPBOT1
              SNDBT(L,K)   = TMPVAL*TMPBOT2
              DO NS=1,NSND  
                NX = NS+NSED  
                TMPVAL = SNDB(L,K,NS)
                SNDB(L,K+1,NS) = TMPVAL*TMPBOT1
                SNDB(L,K,NS)   = TMPVAL*TMPBOT2
                STFPOCB(L,K+1,NX) = STFPOCB(L,K,NS)  
                VFRBED(L,K+1,NX) = VFRBED(L,K,NS)  
              ENDDO  
            ENDIF
            KBT(L) = 2
          ENDIF
        ENDIF  
      ENDDO  
  
      ! *** COMBINE THIN PARENT LAYERS WITH LAYER BELOW TO FORM NEW PARENT  
      DO L=LF,LL  
        K = KBT(L)-1
        IF( K == 0 )CYCLE
        IF( HBED(L,K) < HBEDMIN )THEN  ! *** LAYER BELOW TOP/ACTIVE LAYER
          IF( KBT(L) > 2 )THEN  
            TMPBOT1 = HBED(L,KBT(L)-2)/(1.+VDRBED(L,KBT(L)-2))    ! *** 2 LAYERS BELOW TOP/ACTIVE LAYER
            TMPBOT2 = HBED(L,KBT(L)-1)/(1.+VDRBED(L,KBT(L)-1))    ! *** 1 LAYERS BELOW TOP/ACTIVE LAYER
            TMPTOP1 = VDRBED(L,KBT(L)-2)*TMPBOT1  
            TMPTOP2 = VDRBED(L,KBT(L)-1)*TMPBOT2  
            VDRBED(L,KBT(L)-2) = (TMPTOP1+TMPTOP2)/(TMPBOT1+TMPBOT2)  
            PORBED(L,KBT(L)-2) = VDRBED(L,KBT(L)-2)/(1.+VDRBED(L,KBT(L)-2))  
            HBED(L,KBT(L)-2) = HBED(L,KBT(L)-2) + HBED(L,KBT(L)-1)  
            HBED(L,KBT(L)-1) = HBED(L,KBT(L))  
            VDRBED(L,KBT(L)-1) = VDRBED(L,KBT(L))  
            PORBED(L,KBT(L)-1) = PORBED(L,KBT(L))  
            
            HBED(L,KBT(L)) = 0.  
            VDRBED(L,KBT(L)) = 0.  
            PORBED(L,KBT(L)) = 0.  
            STDOCB(L,KBT(L)) = 0.  
            STPOCB(L,KBT(L)) = 0.  
            BDENBED(L,KBT(L)) = 0.
            
            IF( ISTRAN(6) >= 1 )THEN  
              SEDBT(L,KBT(L)-2) = 0.  
              SEDBT(L,KBT(L)-1) = 0.  
              DO NS=1,NSED  
                SEDB(L,KBT(L)-2,NS) = SEDB(L,KBT(L)-2,NS)+SEDB(L,KBT(L)-1,NS)  
                SEDBT(L,KBT(L)-2) = SEDBT(L,KBT(L)-2)    +SEDB(L,KBT(L)-2,NS)  
                SEDB(L,KBT(L)-1,NS) = SEDB(L,KBT(L),NS)  
                SEDBT(L,KBT(L)-1) = SEDBT(L,KBT(L)-1)    +SEDB(L,KBT(L)-1,NS)  
                SEDB(L,KBT(L),NS) = 0.0  
                SEDBT(L,KBT(L)) = 0.0  
                STFPOCB(L,KBT(L),NS) = 0.0  
                VFRBED(L,KBT(L),NS) = 0.0  
              ENDDO  
            ENDIF  
            IF( ISTRAN(7) >= 1 )THEN  
              SNDBT(L,KBT(L)-2) = 0.  
              SNDBT(L,KBT(L)-1) = 0.  
              DO NS=1,NSND  
                SNDB(L,KBT(L)-2,NS) = SNDB(L,KBT(L)-2,NS) + SNDB(L,KBT(L)-1,NS)  
                SNDBT(L,KBT(L)-2)   = SNDBT(L,KBT(L)-2)   + SNDB(L,KBT(L)-2,NS)  
                SNDB(L,KBT(L)-1,NS) = SNDB(L,KBT(L),NS)  
                SNDBT(L,KBT(L)-1)   = SNDBT(L,KBT(L)-1)   + SNDB(L,KBT(L)-1,NS)  
                SNDB(L,KBT(L),NS)   = 0.0  
                SNDBT(L,KBT(L))     = 0.0  
                STFPOCB(L,KBT(L),NS+NSED) = 0.0  
                VFRBED(L,KBT(L),NS+NSED) = 0.0  
              ENDDO  
            ENDIF  
            IF( ISTRAN(5) >= 1 )THEN  
              DO NT=1,NTOX  
                TOXB(L,KBT(L)-2,NT) = TOXB(L,KBT(L)-2,NT) + TOXB(L,KBT(L)-1,NT)  
                TOXB(L,KBT(L)-1,NT) = TOXB(L,KBT(L),NT)  
                TOXB(L,KBT(L),NT)   = 0.0  
              ENDDO  
            ENDIF  
            KBT(L) = KBT(L)-1  

          ELSEIF( K == 1 )THEN  ! *** LAYER BELOW TOP/ACTIVE LAYER IS GETTING THIN
            
            TMPVAL = HBED(L,2) + HBED(L,1)
            IF( TMPVAL < HBEDAL )THEN
              ! *** COLLAPSE THE ACTIVE LAYER BY REMOVING THE PARENT LAYER
              HBED(L,1) = TMPVAL
              HBED(L,2) = 0.
              VDRBED(L,1) = VDRBED(L,2)  
              PORBED(L,1) = PORBED(L,2)  
              IF( ISTRAN(6) >= 1 )THEN  
                SEDBT(L,1) = 0.  
                SEDBT(L,2) = 0.  
                DO NS=1,NSED  
                  SEDB(L,1,NS) = SEDB(L,1,NS)+ SEDB(L,2,NS)
                  SEDBT(L,1)   = SEDBT(L,1)  + SEDB(L,1,NS)  
                  SEDB(L,2,NS) = 0.
                ENDDO  
              ENDIF  
              IF( ISTRAN(7) >= 1 )THEN  
                SNDBT(L,1) = 0.  
                SNDBT(L,2) = 0.  
                DO NS=1,NSND  
                  SNDB(L,1,NS) = SNDB(L,1,NS) + SNDB(L,2,NS)  
                  SNDBT(L,1)   = SNDBT(L,1)   + SNDB(L,1,NS)  
                  SNDB(L,2,NS) = 0.
                ENDDO  
              ENDIF  
              IF( ISTRAN(5) >= 1 )THEN  
                DO NT=1,NTOX  
                  TOXB(L,1,NT) = TOXB(L,1,NT) + TOXB(L,2,NT)  
                  TOXB(L,2,NT) = 0.
                ENDDO  
              ENDIF  

              ! *** ZERO OLD TOP LAYER
              VDRBED(L,KBT(L)) = 0.0  
              PORBED(L,KBT(L)) = 0.0  
              STDOCB(L,KBT(L)) = 0.0  
              STPOCB(L,KBT(L)) = 0.0  

              KBT(L) = KBT(L)-1  
            ENDIF
            
          ELSEIF( HBED(L,KBT(L)) < 0.0 )THEN  
            ! *** ZERO NEGATIVE THICKNESSES
            HBED(L,KBT(L)) = 0.0  
            VDRBED(L,KBT(L)) = 0.0  
            PORBED(L,KBT(L)) = 0.0  
            STDOCB(L,KBT(L)) = 0.0  
            STPOCB(L,KBT(L)) = 0.0  
          ENDIF  
          ! *** PMC END BLOCK  
        ENDIF  
      ENDDO  

      ! *** UPDATE BULK DENSITY  
      DO L=LF,LL  
        ! ***REMOVED KB LOOP, ONLY COMPUTE THE TOP LAYER.  CONSOLIDATION IS HANDLED LATER IN CALBED
        K = KBT(L)
        IF( HBED(L,K) > 0. )THEN  
          ! *** UPDATE TOTAL/WET DENSITY
          BDENBED(L,K) = 1000.*PORBED(L,K) + 0.001*(SEDBT(L,K) + SNDBT(L,K))/HBED(L,K)  
        ELSE  
          BDENBED(L,K) = 0.  
        ENDIF  
      ENDDO  

    ENDDO  ! *** END OF DOMAIN LOOP
    !$OMP END PARALLEL DO

  ENDIF    ! *** END OF ARMORING OPTION SECTION

  ! ***************************************************************************
  RETURN  

END  

