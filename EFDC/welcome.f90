
SUBROUTINE WELCOME  
    
  ! *** CHANGE RECORD  
  ! ***   
  ! DATE MODIFIED     BY           
  ! 06/25/2006        Paul M. Craig
  !                   Updated Code to Fortran 90

  WRITE(6,1)  

1 FORMAT('***********************************************************************'  &
      ,/,'*                                                                     *'  &  
      ,/,'*                                                                     *'  &
      ,/,'*         EEEEEEEEE    FFFFFFFFF    DDDDDDDD       CCCCCCCC           *'  &
      ,/,'*        EEE          FFF          DDD     DD    CCC      CC   +      *'  &
      ,/,'*       EEE          FFF          DDD     DD    CCC            +      *'  &
      ,/,'*      EEEEEEEE     FFFFFFFF     DDD     DD    CCC         +++++++++  *'  &
      ,/,'*     EEE          FFF          DDD     DD    CCC              +      *'  &  
      ,/,'*    EEE          FFF          DDD     DD    CCC      CC       +      *'  &  
      ,/,'*   EEEEEEEEE    FFF          DDDDDDDDDD      CCCCCCCCC               *'  &  
      ,/,'*                                                                     *'  &  
      ,/,'*                ENVIRONMENTAL FLUID DYNAMICS CODE (PLUS)             *'  &  
      ,/,'*                ORIGINALLY DEVELOPED BY JOHN M. HAMRICK              *'  &  
      ,/,'*                                                                     *'  &  
      ,/,'*   EFDC+   BY DSI, LLC,  EDMONDS WA, USA                             *'  &  
      ,/,'*              VERTICAL LAYERING WITH SIGMA-STRETCHED OR SIGMA-ZED    *'  &  
      ,/,'*              MULTI-THREADED VERSION USING OpenMP                    *'  &  
      ,/,'*              WITH SEDZLJ, HYDROKINETIC DEVICES (SNL),               *'  &  
      ,/,'*              LAGRANGIAN PARTICLE TRACKING AND RPEM SUB-MODELS       *'  &  
      ,/,'*                                                                     *'  &  
      ,/,'*   EFDC+ HAS BEEN CUSTOMIZED TO WORK WITH DSI''S EFDC_EXPLORER 8.5.0  *'  &  
      ,/,'*                                                                     *'  &  
      ,/,'*                        EFDC+ 8.5 EPA RELEASE                        *'  &  
      ,/,'*                      VERSION DATE: 19 AUG 2019                      *'  &  
      ,/,'*                                                                     *'  &  
      ,/,'***********************************************************************') 
  RETURN  
  
END  
