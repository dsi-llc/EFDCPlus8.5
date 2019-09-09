SUBROUTINE SOLVSMBE(SMV1,SMV2,SMA11,SMA22,SMA1,SMA2,SMB11,SMB22)  
  
  ! SOLVE 2X2 MATRIX  
  !----------------------------------------------------------------------!  
  ! CHANGE RECORD  
  ! DATE MODIFIED     BY               DESCRIPTION
  !----------------------------------------------------------------------!
  ! 2013-03           Paul M. Craig    Restructed to F90, updated variables for OMP
  
  IMPLICIT NONE

  REAL, INTENT(IN)  :: SMA11, SMA22, SMA1, SMA2, SMB11, SMB22
  REAL, INTENT(OUT) :: SMV1, SMV2
  REAL :: SMA12, SMA21, SMDET

  SMA12 = -SMA2  
  SMA21 = -SMA1  
  SMDET = SMA11*SMA22 - SMA12*SMA21  
  IF( SMDET == 0.0 )THEN  
    PRINT*, 'SINGULAR MATRIX: A11, A12, A21, A22, B11, B22'  
    PRINT*, SMA11,SMA12,SMA21,SMA22,SMB11,SMB22  
    STOP  
  ENDIF  
  SMDET = 1.0 / SMDET  
  SMV1 = (SMB11*SMA22 - SMB22*SMA12) * SMDET  
  SMV2 = (SMB22*SMA11 - SMB11*SMA21) * SMDET  

END SUBROUTINE

