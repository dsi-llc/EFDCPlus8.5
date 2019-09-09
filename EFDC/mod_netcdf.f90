MODULE MOD_NETCDF
#ifdef NCOUT

USE GLOBAL
USE JULIANMOD
USE CONVERTWGS84
USE NETCDF
!USE DRIFTER,ONLY:XYZSCL,MOC

IMPLICIT NONE

TYPE NC_VARIABLE
  INTEGER(2)     :: ID = 0
  CHARACTER(15)  :: NAME,UNITS,POSITIVE,GRID_MAPPING,COORDINATES
  INTEGER(2)     :: DIM_TYPE
  REAL(4)        :: FILLVALUE
  CHARACTER(100) :: STANDARD_NAME,LONG_NAME
END TYPE

INTEGER(4), PARAMETER :: CNRCNT=4, TIMECNT=NF90_UNLIMITED, NC_VAR_CNT = 41
REAL(4),    PARAMETER :: MISSING_VALUE=-999.
INTEGER(4) :: NC_ID, NTI, NC_VAR_CNTM
INTEGER(4) :: TIME_DIM, ROW_DIM, COL_DIM, LYR_DIM, CNR_DIM, SED_DIM, SND_DIM, TOX_DIM, CWQ_DIM, LPT_DIM
INTEGER(4) :: ROWCNT, COLCNT, LYRCNT, NCDFTYPE
INTEGER(4) :: NC_XC, NC_YC, NC_XB, NC_YB, NC_CRS, NC_SIG, NC_KC, NC_TIM
INTEGER(4) :: NC_TRK_XC,NC_TRK_YC,NC_TRK_CRS,NC_TRK_SIG,NC_TRK_VOL

CHARACTER(10)                  :: NC_DATESTAMP
INTEGER(4), ALLOCATABLE        :: NC_IDX(:)
TYPE(NC_VARIABLE), ALLOCATABLE :: NC_VAR(:)
REAL(4),ALLOCATABLE            :: TRK_LON(:),TRK_LAT(:),TRK_LEV(:),TRK_VOL(:)
DATA NCDFTYPE /0/

CONTAINS

  SUBROUTINE NC_DEFINE_VARS(NC_ID,IDX)
  
    INTEGER :: NC_ID,IDX,ID,STAT
    
    IF (NCDFTYPE==0) THEN
      NCDFTYPE = 1
      ALLOCATE(NC_VAR(NC_VAR_CNT),NC_IDX(NC_VAR_CNT))
      ALLOCATE(TRK_LON(NPD),TRK_LAT(NPD),TRK_LEV(NPD),TRK_VOL(NPD))

      NC_VAR(1:40) = (/&
        NC_VARIABLE( 1,'Bottom','m',  'down',  'crs','lat lng',2,MISSING_VALUE,'','Bottom Bathymetry'),&
        NC_VARIABLE( 2,'WSEL',  'm',  'up',  'crs','lat lng',3,MISSING_VALUE,'','Water Surface Elevation'),&
        NC_VARIABLE( 3,'Vx',    'm/s','',    'crs','lat lng',4,MISSING_VALUE,'','Eastward Water Velocity'),&
        NC_VARIABLE( 4,'Vy',    'm/s','',    'crs','lat lng',4,MISSING_VALUE,'','Northward Water Velocity'),&
        NC_VARIABLE( 5,'Vz',    'm/s','',  'crs','lat lng',4,MISSING_VALUE,'','Upward Water Velocity'),&
        NC_VARIABLE( 6,'TauB',  'N/m2','',   'crs','lat lng',3,MISSING_VALUE,'','Bottom Shear Stress'),&
        NC_VARIABLE( 7,'Wx',    'm/s','',    'crs','lat lng',3,MISSING_VALUE,'','Eastward Wind Velocity'),&
        NC_VARIABLE( 8,'Wy',    'm/s','',    'crs','lat lng',3,MISSING_VALUE,'','Northward Wind Velocity'),&
        NC_VARIABLE( 9,'Hs',    'm','',     'crs','lat lng',3,MISSING_VALUE,'','Significant Wave Height'),&
        NC_VARIABLE(10,'Dp',    'degree','','crs','lat lng',3,MISSING_VALUE,'','Wave Direction'),&
        NC_VARIABLE(11,'Tp',    's','',     'crs','lat lng',3,MISSING_VALUE,'','Wave Period'),&      
        NC_VARIABLE(12,'Salt',  'ppt','',    'crs','lat lng',4,MISSING_VALUE,'','Salinity'),&
        NC_VARIABLE(13,'Temp',  'C','', 'crs','lat lng',4,MISSING_VALUE,'','Temperature'),&
        NC_VARIABLE(14,'Dye',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dye'),&
        NC_VARIABLE(15,'SFL',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Shellfish Larvae'),&
        NC_VARIABLE(16,'Toxic', 'mg/l','',    'crs','lat lng',7,MISSING_VALUE,'','Toxics'),&
        NC_VARIABLE(17,'SED',   'mg/l','',    'crs','lat lng',5,MISSING_VALUE,'','Cohesive Sediment'),&
        NC_VARIABLE(18,'SND',   'mg/l','',    'crs','lat lng',6,MISSING_VALUE,'','Non-cohesive Sediment'),&
        NC_VARIABLE(19,'CHC',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Cyanobacteria'),&
        NC_VARIABLE(20,'CHD',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Diatom Algae'),&
        NC_VARIABLE(21,'CHG',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Green Algae'),&
        NC_VARIABLE(22,'ROC',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Refractory Particulate Organic Carbon'),&
        NC_VARIABLE(23,'LOC',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Labile Particulate Organic Carbon'),&
        NC_VARIABLE(24,'DOC',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Organic Carbon'),&
        NC_VARIABLE(25,'ROP',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Refractory Particulate Organic Phosphorus'),&
        NC_VARIABLE(26,'LOP',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Labile Particulate Organic Phosphorus'),&
        NC_VARIABLE(27,'DOP',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Organic Phosphorus'),&
        NC_VARIABLE(28,'P4D',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Total Phosphate'),&
        NC_VARIABLE(29,'RON',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Refractory Particulate Organic Nitrogen'),&
        NC_VARIABLE(30,'LON',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Labile Particulate Organic Nitrogen'),&
        NC_VARIABLE(31,'DON',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Organic Nitrogen'),&
        NC_VARIABLE(32,'NHX',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Ammonia Nitrogen'),&
        NC_VARIABLE(33,'NOX',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Nitrate Nitrogen'),&
        NC_VARIABLE(34,'SUU',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Particulate Biogenic Silica'),&
        NC_VARIABLE(35,'SAA',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Available Silica'),&
        NC_VARIABLE(36,'COD',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Chemical Oxygen Demand'),&
        NC_VARIABLE(37,'DOX',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Dissolved Oxygen'),&
        NC_VARIABLE(38,'TAM',   'mg/l','',    'crs','lat lng',4,MISSING_VALUE,'','Total Active Metal'),&
        NC_VARIABLE(39,'FCB',   'mpn/100ml','','crs','lat lng',4,MISSING_VALUE,'','Fecal Coliform Bacteria'),&
        NC_VARIABLE(40,'MAC',   'ug/l','',    'crs','lat lng',4,MISSING_VALUE,'','Macroalgae') /)
      
      IF ( ISPD >= 2 ) THEN
        IF (ANY(ISOILSPI == 1)) NC_VAR(41) = NC_VARIABLE(41,'MOC', 'kg','',  'crs','lat lng',3,MISSING_VALUE,'','Oil Mass')
      ENDIF
      
    ENDIF
  
    IF (NC_VAR(IDX).ID > 0) THEN
      SELECT CASE (NC_VAR(IDX).DIM_TYPE)
        CASE (1)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/TIME_DIM/),ID)
        CASE (2)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/COL_DIM,ROW_DIM/),ID)
        CASE (3)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/COL_DIM,ROW_DIM,TIME_DIM/),ID)
        CASE (4)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/LYR_DIM,COL_DIM,ROW_DIM,TIME_DIM/),ID)
        CASE (5)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/SED_DIM,LYR_DIM,COL_DIM,ROW_DIM,TIME_DIM/),ID)
        CASE (6)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/SND_DIM,LYR_DIM,COL_DIM,ROW_DIM,TIME_DIM/),ID)
        CASE (7)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/TOX_DIM,LYR_DIM,COL_DIM,ROW_DIM,TIME_DIM/),ID)
        CASE (8)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/CWQ_DIM,LYR_DIM,COL_DIM,ROW_DIM,TIME_DIM/),ID) 
        CASE (9)
          STAT=NF90_DEF_VAR(NC_ID,NC_VAR(IDX).NAME,NF90_FLOAT,(/LPT_DIM,TIME_DIM/),ID)         
      END SELECT           
            
      STAT=CHECK_ERR(STAT,'def_'//TRIM(NC_VAR(IDX).NAME))
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,ID,1,1,DEFLEV)
      IF (LEN_TRIM(NC_VAR(IDX).STANDARD_NAME) > 0) THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'standard_name',TRIM(NC_VAR(IDX).STANDARD_NAME))
      ENDIF
      IF (LEN_TRIM(NC_VAR(IDX).LONG_NAME) > 0) THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'long_name',TRIM(NC_VAR(IDX).LONG_NAME))
      ENDIF
      IF (LEN_TRIM(NC_VAR(IDX).UNITS) > 0) THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'units',TRIM(NC_VAR(IDX).UNITS))
      ENDIF
      IF (LEN_TRIM(NC_VAR(IDX).POSITIVE) > 0) THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'positive',TRIM(NC_VAR(IDX).POSITIVE))
        STAT=NF90_PUT_ATT(NC_ID,ID,'_coordinatezispositive',TRIM(NC_VAR(IDX).POSITIVE))
      ENDIF
      IF (LEN_TRIM(NC_VAR(IDX).GRID_MAPPING) > 0) THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'grid_mapping',TRIM(NC_VAR(IDX).GRID_MAPPING))
      ENDIF
      IF (LEN_TRIM(NC_VAR(IDX).COORDINATES) > 0) THEN
        STAT=NF90_PUT_ATT(NC_ID,ID,'coordinates',TRIM(NC_VAR(IDX).COORDINATES))
      ENDIF
      STAT=NF90_PUT_ATT(NC_ID,ID,'_fillvalue',NC_VAR(IDX).FILLVALUE)
      NC_VAR(IDX).ID = ID
      NC_IDX(IDX) = ID
    ENDIF
    
  END SUBROUTINE

  SUBROUTINE NC_SET_VARS(NC_ID)
    INTEGER :: NC_ID,K
  
    ! ** ALWAYS EXPORT: ZBOT,WSEL,VELX,VELY,VELZ
    DO K=1,5
      CALL NC_DEFINE_VARS(NC_ID,K)                    ! K=1:5 
    ENDDO
    
    ! ** OPTIONALLY EXPORT
    IF (ISNCDF(10) >= 1) CALL NC_DEFINE_VARS(NC_ID,6) ! 6:TAUB
    
    IF (ISNCDF(11) >= 1) THEN
      CALL NC_DEFINE_VARS(NC_ID,7)                    ! 7:WINX
      CALL NC_DEFINE_VARS(NC_ID,8)                    ! 8:WINY
    ENDIF 

    IF (ISNCDF(12) >= 1) THEN
      DO K=9,11
        CALL NC_DEFINE_VARS(NC_ID,K)                  ! WINDWAVE: 9:11 WHEI,WANG,WPER
      ENDDO
    ENDIF
    
    DO K=1,7
      IF (ISNCDF(K) >= 1) CALL NC_DEFINE_VARS(NC_ID,  11+K)  ! 12-18: SAL,TEM,...SND
    ENDDO
    
    DO K=1,NWQV
      IF(ISTRWQ(K) > 0 .AND. ISNCDF(8) >= 1) CALL NC_DEFINE_VARS(NC_ID, 18 + K)   ! 19-40: WQV
    ENDDO
    
    ! ** OIL LPT
    IF( ISNCDF(9) >= 1) THEN  
      IF (ANY(ISOILSPI == 1)) CALL NC_DEFINE_VARS(NC_ID, 41)  
    ENDIF
    
  END SUBROUTINE

  INTEGER FUNCTION NC_CREATE_FILE(NC_ID, FILENAME)
    ! ** VARIABLE DECLARATIONS
    CHARACTER(*),INTENT(IN) :: FILENAME
    INTEGER,INTENT(OUT) :: NC_ID
    CHARACTER(20) :: BDATE,ZONESTR*3
    INTEGER :: STAT
    
    ! ** ENTER DEFINE MODE
    STAT=CHECK_ERR(NF90_CREATE(FILENAME, NF90_CLOBBER+NF90_NETCDF4, NC_ID),'nc_create');!NF90_HDF5
    IF(STAT /= NF90_NOERR) THEN
        NC_CREATE_FILE = STAT
        RETURN
    ENDIF

    ROWCNT = IMX - IMN +1 
    COLCNT = JMX - JMN +1 
    LYRCNT = KC
    BDATE = TRIM(BASEDATE)//' '//TRIM(BASETIME)

    CALL UTMPARS

    ! ** DEFINE DIMENSIONS
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'time', TIMECNT, TIME_DIM),'time_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'imax', ROWCNT, ROW_DIM),'row_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'jmax', COLCNT, COL_DIM),'col_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'kmax', LYRCNT, LYR_DIM),'lyr_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'cnr', CNRCNT, CNR_DIM),'cnr_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'ntox', NTOX, TOX_DIM),'tox_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'nsed', NSED, SED_DIM),'sed_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'nsnd', NSND, SND_DIM),'snd_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'nwqv', NWQV, CWQ_DIM),'cwq_dim');
    STAT=CHECK_ERR(NF90_DEF_DIM(NC_ID,'npd', NPD, LPT_DIM),'lpt_dim');   
    
    ! ** DEFINE VARIABLES
    STAT=CHECK_ERR(NF90_DEF_VAR(NCID=NC_ID,NAME='crs',XTYPE=NF90_CHAR,VARID=NC_CRS),'def_crs') !DIMIDS=0,
    STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'grid_mapping_name','latitude_longitude')
    STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'longitude_of_prime_meridian',0.)
    STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'semi_major_axis',R_MJR)
    STAT=NF90_PUT_ATT(NC_ID,NC_CRS,'inverse_flattening',1._8/FLA)

    STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lng',NF90_FLOAT,(/COL_DIM,ROW_DIM/),NC_XC),'def_lng')
    STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_XC,1,1,DEFLEV)
    STAT=NF90_PUT_ATT(NC_ID,NC_XC,'standard_name','longitude')
    STAT=NF90_PUT_ATT(NC_ID,NC_XC,'long_name','longitude at grid cell center')
    STAT=NF90_PUT_ATT(NC_ID,NC_XC,'units','degrees_east')
    STAT=NF90_PUT_ATT(NC_ID,NC_XC,'bounds','lng_bnds')
    STAT=NF90_PUT_ATT(NC_ID,NC_XC,'projection','geographic')
    STAT=NF90_PUT_ATT(NC_ID,NC_XC,'_FillValue',MISSING_VALUE)

    STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lat',NF90_FLOAT,(/COL_DIM,ROW_DIM/),NC_YC),'def_lat')
    STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_YC,1,1,DEFLEV)
    STAT=NF90_PUT_ATT(NC_ID,NC_YC,'standard_name','latitude')
    STAT=NF90_PUT_ATT(NC_ID,NC_YC,'long_name','latitude at grid cell center')
    STAT=NF90_PUT_ATT(NC_ID,NC_YC,'units','degrees_north')
    STAT=NF90_PUT_ATT(NC_ID,NC_YC,'bounds','lat_bnds')
    STAT=NF90_PUT_ATT(NC_ID,NC_YC,'projection','geographic')
    STAT=NF90_PUT_ATT(NC_ID,NC_YC,'_FillValue',MISSING_VALUE)

    STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lng_bnds',NF90_FLOAT,(/CNR_DIM,COL_DIM,ROW_DIM/),NC_XB),'def_xb')
    STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_XB,1,1,DEFLEV)
    STAT=NF90_PUT_ATT(NC_ID,NC_XB,'_FillValue',MISSING_VALUE)

    STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'lat_bnds',NF90_FLOAT,(/CNR_DIM,COL_DIM,ROW_DIM/),NC_YB),'def_yb')
    STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_YB,1,1,DEFLEV)
    STAT=NF90_PUT_ATT(NC_ID,NC_YB,'_FillValue',MISSING_VALUE)

    IF(IGRIDV > 0) THEN
      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'layers',NF90_BYTE,(/COL_DIM,ROW_DIM/),NC_KC),'def_lat')
      STAT=NF90_PUT_ATT(NC_ID,NC_KC,'long_name','number of vertical layers')
      
      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'sigma',NF90_FLOAT,(/LYR_DIM,COL_DIM,ROW_DIM/), NC_SIG),'def_sigma')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'standard_name','ocean_sigma_coordinate')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'long_name','sigma at layer midpoints')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'units','sigma_level')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'positive','up')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'formula_terms','sigma: sigma eta: WSEL depth: Bottom')
    ELSE
      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'sigma',NF90_FLOAT,(/LYR_DIM/), NC_SIG),'def_sigma')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'standard_name','ocean_sigma_coordinate')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'long_name','sigma at layer midpoints')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'units','sigma_level')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'positive','up')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'formula_terms','sigma: sigma eta: WSEL depth: Bottom')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateZisPositive','up')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateTransformType','Vertical')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateAxisType','GeoZ')
      STAT=NF90_PUT_ATT(NC_ID,NC_SIG,'_CoordinateAxes','sigma')
    ENDIF
    
    STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'time',NF90_DOUBLE,(/TIME_DIM/), NC_TIM),'def_time')
    STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'standard_name','time')
    STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'long_name','time')
    STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'units','days since '//TRIM(BDATE))
    STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'calendar','julian')
    STAT=NF90_PUT_ATT(NC_ID,NC_TIM,'axis','T')
    
    ! *** LPT
    IF( ISNCDF(9) >= 1) THEN 
      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_LAT',NF90_FLOAT,(/LPT_DIM,TIME_DIM/),NC_TRK_YC),'def_trk_lat')
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_TRK_YC,1,1,DEFLEV)
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'standard_name','')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'long_name','Drifter Latitude')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'units','degrees_north')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'valid_min','-90')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'valid_max','90')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'grid_mapping','crs')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_YC,'_FillValue',MISSING_VALUE)
    
      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_LON',NF90_FLOAT,(/LPT_DIM,TIME_DIM/),NC_TRK_XC),'def_trk_lon')
      STAT=NF90_DEF_VAR_DEFLATE(NC_ID,NC_TRK_XC,1,1,DEFLEV)
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'standard_name','')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'long_name','Drifter Longitude')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'units','degrees_east')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'valid_min','-180')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'valid_max','180')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'grid_mapping','crs')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_XC,'_FillValue',MISSING_VALUE)

      STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_LEV',NF90_FLOAT,(/LPT_DIM,TIME_DIM/), NC_TRK_SIG),'def_trk_sigma')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'standard_name','')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'long_name','Drifter Elevation')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'units','m')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'grid_mapping','crs')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateZisPositive','up')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateTransformType','Vertical')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateAxisType','GeoZ')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_CoordinateAxes','lev')
      STAT=NF90_PUT_ATT(NC_ID,NC_TRK_SIG,'_FillValue',MISSING_VALUE)

      IF (ANY(ISOILSPI == 1)) THEN
        STAT=CHECK_ERR(NF90_DEF_VAR(NC_ID,'TRK_VOL',NF90_FLOAT,(/LPT_DIM,TIME_DIM/), NC_TRK_VOL),'def_trk_vol')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'standard_name','')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'long_name','Drifter Oil Volume')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'units','m3')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'grid_mapping','crs')
        STAT=NF90_PUT_ATT(NC_ID,NC_TRK_VOL,'_FillValue',MISSING_VALUE)
      ENDIF
    ENDIF
    ! *** END LPT

    CALL NC_SET_VARS(NC_ID)

    ! ** ASSIGN GLOBAL ATTRIBUTES

    WRITE(ZONESTR,'(I3)') ABS(UTMZ)
    STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'Conventions','CF-1.4'),'att_cf')
    STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'Base_date',BDATE),'att_base')
    STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'Project',TRIM(PROJ)),'att_prj')
    IF (UTMZ > 0) THEN
        STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'utm_zone','UTM Zone '//TRIM(ZONESTR)//' Northern Hemisphere'),'')
    ELSE
        STAT=CHECK_ERR(NF90_PUT_ATT(NC_ID,NF90_GLOBAL,'utm_zone','UTM Zone '//TRIM(ZONESTR)//' Southern Hemisphere'),'')
    ENDIF

    ! ** LEAVE DEFINE MODE
    STAT=CHECK_ERR(NF90_ENDDEF(NC_ID),'end_def');
    NC_CREATE_FILE = STAT
  END FUNCTION
  
  SUBROUTINE NC_CLOSE_FILE(NC_ID)
    INTEGER NC_ID,STATUS
    IF (NCDFOUT > 1) THEN
        STATUS=CHECK_ERR(NF90_CLOSE(NC_ID),'close')
        IF(STATUS /= NF90_NOERR) THEN
            PRINT * ,'Cannot close nc file!'
            RETURN
        ENDIF
    ENDIF
  END SUBROUTINE

  INTEGER FUNCTION CHECK_ERR(STATUS,MSG)
    INTEGER :: STATUS
    CHARACTER(*) :: MSG
    CHECK_ERR = STATUS
    IF (STATUS .NE. NF90_NOERR) THEN
        PRINT *, MSG//': '//TRIM(NF90_STRERROR(STATUS))
        STOP
    END IF
  END FUNCTION 

  INTEGER FUNCTION NC_WRITE_SCALAR1D(NC_ID,ID,SNAPSHOT,FIELD1D)
    ! ** WRITE 1D SCALAR FIELD TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT
    REAL(4),INTENT(IN) :: FIELD1D(:)
    INTEGER :: STAT,STR1D(2),CNT1D(2)

    STR1D = (/1, SNAPSHOT/)
    CNT1D = (/NPD,1/)
    STAT=NF90_PUT_VAR(NC_ID, ID, FIELD1D, STR1D, CNT1D)
    NC_WRITE_SCALAR1D=STAT
  END FUNCTION
  
  INTEGER FUNCTION NC_WRITE_SCALAR2D(NC_ID,ID,FIELD2D,SNAPSHOT)
    ! ** WRITE 2D SCALAR FIELD TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT
    REAL(4),INTENT(IN) :: FIELD2D(:,:)
    INTEGER :: STAT,STR2D(3),CNT2D(3)

    STR2D = (/1, 1, SNAPSHOT/)
    CNT2D = (/COLCNT,ROWCNT,1/)
    STAT=NF90_PUT_VAR(NC_ID, ID, FIELD2D, STR2D, CNT2D)
    NC_WRITE_SCALAR2D=STAT
  END FUNCTION

 INTEGER FUNCTION NC_WRITE_SCALAR3D(NC_ID,ID,FIELD3D,SNAPSHOT)
    ! ** WRITE 3D SCALAR FIELD TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT
    REAL(4),INTENT(IN) :: FIELD3D(:,:,:)
    INTEGER :: STAT,STR3D(4),CNT3D(4)
    STR3D = (/1, 1, 1, SNAPSHOT/)
    CNT3D = (/LYRCNT,COLCNT,ROWCNT,1/)
    STAT=NF90_PUT_VAR(NC_ID, ID, FIELD3D, STR3D, CNT3D)
    NC_WRITE_SCALAR3D=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_CLASS3D(NC_ID,ID,FIELD3D,SNAPSHOT,CLSCNT)
    ! ** WRITE 3D SCALAR FIELD TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,ID,SNAPSHOT,CLSCNT
    REAL(4),INTENT(IN) :: FIELD3D(:,:,:,:)
    INTEGER :: STAT,STR3D(5),CNT3D(5)
    STR3D = (/1, 1, 1, 1, SNAPSHOT/)
    CNT3D = (/CLSCNT,LYRCNT,COLCNT,ROWCNT,1/)
    STAT=NF90_PUT_VAR(NC_ID, ID, FIELD3D, STR3D, CNT3D)
    NC_WRITE_CLASS3D=STAT
  END FUNCTION

  SUBROUTINE READCORN
    ! ** USE GLOBAL ,ONLY:XCR,YCR,UCOR,IC,JC
    USE INFOMOD,ONLY:READSTR
    IMPLICIT NONE
    INTEGER(4)::I,J,N
    CHARACTER*80 :: STR*200
    
    OPEN(UCOR,FILE='corners.inp',ACTION='READ')
    STR = READSTR(UCOR)
    ALLOCATE(XCR(4,JC,IC),YCR(4,JC,IC))
    XCR = MISSING_VALUE
    YCR = MISSING_VALUE
    DO WHILE(1)
        READ(UCOR,*,END=100,ERR=998) I,J,(XCR(N,J,I),YCR(N,J,I),N=1,4)
    ENDDO
    100 CLOSE(UCOR)
    
    RETURN
    998 CLOSE(UCOR)
    STOP 'CORNERS.INP READING ERROR!'
  END SUBROUTINE

  INTEGER FUNCTION NC_WRITE_GRID(NC_ID,GRD_X,GRD_Y)
    ! ** WRITE HORIZONTAL GRID TO NETCDF FILE
    REAL,INTENT(IN) :: GRD_X(:,:,:),GRD_Y(:,:,:)
    REAL(4) :: LON(4,COLCNT,ROWCNT), LAT(4,COLCNT,ROWCNT)
    REAL(4) :: XC(COLCNT,ROWCNT),YC(COLCNT,ROWCNT)
    REAL(8)   :: XM(1),YM(1),XLL(1),YLL(1)
    INTEGER,INTENT(IN) :: NC_ID
    INTEGER :: STAT,I,J,N

    LON = MISSING_VALUE
    LAT = MISSING_VALUE
    XC  = MISSING_VALUE
    YC  = MISSING_VALUE

    DO I=1,ROWCNT
        DO J=1,COLCNT
            IF (GRD_X(1,J,I) > 0) THEN
                XC(J,I)=0.
                YC(J,I)=0.
                DO N=1,CNRCNT
                    XM(1)=GRD_X(N,J,I)
                    YM(1)=GRD_Y(N,J,I)
                    CALL UTMR_WGS84(XM,YM,XLL,YLL)
                    LON(N,J,I) = XLL(1)
                    LAT(N,J,I) = YLL(1)
                    XC(J,I)=XC(J,I)+LON(N,J,I)
                    YC(J,I)=YC(J,I)+LAT(N,J,I)
                ENDDO
                XC(J,I)=XC(J,I)/CNRCNT
                YC(J,I)=YC(J,I)/CNRCNT
            ENDIF
        ENDDO
    ENDDO

    STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_XB,LON),'xb_put') 
    STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_YB,LAT),'yb_put') 
    STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_XC,XC),'xc_put') 
    STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_YC,YC),'yc_put') 
    NC_WRITE_GRID=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_ZBOT(NC_ID)
    ! ** WRITE BOTTOM DEPTH AND SIGMA LEVEL TO A NETCDF FILE
    ! 2017-09-07, NTL: Updated NetCDF layers from top=1 to bottom = KC as the
    !                  ocean_sigma_coordinate ranges from 0 at surface to -1 at bottom 
    INTEGER NC_ID,STAT,I,J,K,L,N
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: BATHY
    REAL(4) :: SIGMA(LYRCNT), SUMDZC
    BYTE, ALLOCATABLE, DIMENSION(:,:) :: KCSGZ
    REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: DZSGZ

    IF(IGRIDV > 0) THEN ! Sigma-Zed level
        ALLOCATE(KCSGZ(COLCNT,ROWCNT),DZSGZ(LYRCNT,COLCNT,ROWCNT)) 
        KCSGZ = KMINV
        DZSGZ = MISSING_VALUE
        DO L = 2,LA  
            I = IL(L) - IMN + 1
            J = JL(L) - JMN + 1
            KCSGZ(J,I) = KC - KSZ(L) + 1    ! Number of vertical layers
            SUMDZC = 0.0
            DO N=1,LYRCNT
                K = LYRCNT - N + 1
                DZSGZ(N,J,I) = -(SUMDZC + 0.5*DZC(L,K))
                SUMDZC = SUMDZC + DZC(L,K)
            ENDDO
        ENDDO     
        STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_KC,KCSGZ),'put_kc')
        STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_SIG,DZSGZ),'put_dz')
        DEALLOCATE(KCSGZ,DZSGZ)
    ELSE  ! ** Standard sigma stretched level
        DO N=1,LYRCNT
            K = LYRCNT - N + 1
            SIGMA(N) = SUM(DZCK(1:K)) - 0.5*DZCK(K) - 1.0
        ENDDO
        STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_SIG,SIGMA),'put_sig')
    ENDIF

    ! ** Write bottom bathymetry
    ALLOCATE(BATHY(COLCNT,ROWCNT))
    BATHY = MISSING_VALUE
    DO L = 2,LA  
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        BATHY(J,I) = -BELV(L)
    ENDDO     
    STAT=CHECK_ERR(NF90_PUT_VAR(NC_ID,NC_IDX(1),BATHY),'put_zb')
    DEALLOCATE(BATHY)
    NC_WRITE_ZBOT=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_WSEL(NC_ID, SNAPSHOT)
    ! ** WRITE WATER SURFACE ELEVATION
    INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: WSEL
    INTEGER :: STAT, I,J,L
  
    ALLOCATE(WSEL(COLCNT,ROWCNT))
    WSEL = MISSING_VALUE
    DO L = 2,LA
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        IF (HP(L) > HDRY) THEN  
            WSEL(J,I) = BELV(L) + HP(L)
        ENDIF
    ENDDO
    STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(2),WSEL,SNAPSHOT)
    DEALLOCATE(WSEL)
    NC_WRITE_WSEL=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_VEL(NC_ID, SNAPSHOT)
    ! ** WRITE VELOCITIES TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
    REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: VX, VY, VZ
    REAL(4) :: UTMP,VTMP
    INTEGER :: I,J,L,K,N,LN,STAT

    ALLOCATE(VX(LYRCNT,COLCNT,ROWCNT),VY(LYRCNT,COLCNT,ROWCNT),VZ(LYRCNT,COLCNT,ROWCNT))
    VX = MISSING_VALUE
    VY = MISSING_VALUE
    VZ = MISSING_VALUE
    DO L = 2,LA
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        IF (HP(L) > HDRY) THEN  
            IF (ROTA == 1) THEN
                LN = LNC(L)  
                DO N=1,KC
                    K = KC - N + 1
                    UTMP = 0.5*(RSSBCE(L)*U(L+1,K)+RSSBCW(L)*U(L,K))  
                    VTMP = 0.5*(RSSBCN(L)*V(LN ,K)+RSSBCS(L)*V(L,K))  
                    VX(N,J,I) = CUE(L)*UTMP+CVE(L)*VTMP  
                    VY(N,J,I) = CUN(L)*UTMP+CVN(L)*VTMP  
                    VZ(N,J,I) = W(L,K)
                ENDDO
            ELSE
                DO N=1,KC
                    K = KC - N + 1
                    VX(N,J,I) = U(L,K)
                    VY(N,J,I) = V(L,K)
                    VZ(N,J,I) = W(L,K)
                ENDDO
            ENDIF
        ENDIF
    ENDDO
    STAT=NC_WRITE_SCALAR3D(NC_ID, NC_IDX(3), VX, SNAPSHOT)
    STAT=NC_WRITE_SCALAR3D(NC_ID, NC_IDX(4), VY, SNAPSHOT)
    STAT=NC_WRITE_SCALAR3D(NC_ID, NC_IDX(5), VZ, SNAPSHOT)
    DEALLOCATE(VX, VY, VZ)
    NC_WRITE_VEL=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_SHEAR(NC_ID, SNAPSHOT)
    ! ** WRITE BOTTOM SHEAR STRESSES TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: SHEAR
    INTEGER :: I,J,L,STAT
  
    ALLOCATE(SHEAR(COLCNT,ROWCNT))
    SHEAR = MISSING_VALUE
    DO L = 2,LA
      I = IL(L) - IMN + 1
      J = JL(L) - JMN + 1
      IF (HP(L) > HDRY) THEN  
        IF( ISTRAN(6) > 0 .OR. ISTRAN(7) > 0 )THEN
          IF( LSEDZLJ )THEN
            SHEAR(J,I) = TAU(L)*0.1             ! DYN/CM2 TO N/M2
          ELSEIF( ISBEDSTR >= 1 )THEN
            SHEAR(J,I) = 1000.*TAUBSED(L)       ! M2/S2 TO N/M2
            IF( ISBEDSTR == 1 )THEN
              SHEAR(J,I) = 1000.*TAUBSND(L)
            ENDIF
          ELSE
            SHEAR(J,I) = 1000.*TAUB(L)
          ENDIF
        ELSE
          SHEAR(J,I) = 1000.*MAX(QQ(L,0),QQMIN)/CTURB2
        ENDIF
      ENDIF
    ENDDO
    STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(6),SHEAR,SNAPSHOT)
    DEALLOCATE(SHEAR)
    NC_WRITE_SHEAR=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_WIND(NC_ID, SNAPSHOT)
    ! ** WRITE WIND VELOCITIES TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: WX, WY
    INTEGER :: I,J,L,STAT

    ALLOCATE(WX(COLCNT,ROWCNT),WY(COLCNT,ROWCNT))
    WX = MISSING_VALUE
    WY = MISSING_VALUE
    DO L = 2,LA
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        WX(J,I) = WNDVELE(L)  
        WY(J,I) = WNDVELN(L)
    ENDDO
    STAT=NC_WRITE_SCALAR2D(NC_ID, NC_IDX(7), WX, SNAPSHOT)
    STAT=NC_WRITE_SCALAR2D(NC_ID, NC_IDX(8), WY, SNAPSHOT)
    DEALLOCATE(WX, WY)
    NC_WRITE_WIND=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_WAVE(NC_ID, SNAPSHOT)
    ! ** WRITE WATER SURFACE ELEVATION
    INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: HS,TP,DIR
    INTEGER :: STAT, I,J,L
  
    ALLOCATE(HS(COLCNT,ROWCNT),TP(COLCNT,ROWCNT),DIR(COLCNT,ROWCNT))
    HS = MISSING_VALUE
    TP = MISSING_VALUE
    DIR = MISSING_VALUE
    DO L = 2,LA
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        IF (HP(L) > HDRY) THEN  
            HS(J,I) = WV(L).HEIGHT
            TP(J,I) = WV(L).FREQ
            DIR(J,I) = WV(L).DIR
        ENDIF
    ENDDO
    STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(9),HS,SNAPSHOT)
    STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(10),DIR,SNAPSHOT)
    STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(11),TP,SNAPSHOT)

    DEALLOCATE(HS,TP,DIR)
    NC_WRITE_WAVE=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_CONS(NC_ID, SNAPSHOT, IDX, VALUES)
    ! ** WRITE SAL/TEM/DYE TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX
    REAL,INTENT(IN) :: VALUES(:,:)
    REAL(4), ALLOCATABLE, DIMENSION(:,:,:) :: CONC
    INTEGER :: I,J,L,K,N,STAT

    ALLOCATE(CONC(LYRCNT,COLCNT,ROWCNT))
    CONC = MISSING_VALUE
    DO L = 2,LA
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        IF (HP(L) > HDRY) THEN  
            DO N=1,KC
                K = KC - N + 1
                CONC(N,J,I) = VALUES(L,K)
            ENDDO
        ENDIF
    ENDDO
    STAT=NC_WRITE_SCALAR3D(NC_ID,NC_IDX(IDX),CONC,SNAPSHOT)
    DEALLOCATE(CONC)
    NC_WRITE_CONS=STAT
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_CONM(NC_ID, SNAPSHOT, IDX, VALUES, CLSCNT)
    ! ** WRITE TOX/SED/SND TO A NETCDF FILE
    INTEGER,INTENT(IN) :: NC_ID,SNAPSHOT,IDX, CLSCNT
    REAL,INTENT(IN) :: VALUES(:,:,:)
    REAL(4), ALLOCATABLE, DIMENSION(:,:,:,:) :: CONC
    INTEGER :: I,J,L,K,N,NS,STAT

    ALLOCATE(CONC(CLSCNT,LYRCNT,COLCNT,ROWCNT))
    CONC = MISSING_VALUE
    DO L = 2,LA
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        IF (HP(L) > HDRY) THEN  
            DO N=1,KC
                K = KC - N + 1
                DO NS=1,CLSCNT
                    CONC(NS,N,J,I) = VALUES(L,K,NS)
                ENDDO
            ENDDO
        ENDIF
    ENDDO
    STAT=NC_WRITE_CLASS3D(NC_ID,NC_IDX(IDX),CONC,SNAPSHOT,CLSCNT)
    DEALLOCATE(CONC)
    NC_WRITE_CONM=STAT
  END FUNCTION
  
  INTEGER FUNCTION NC_WRITE_LPT(NC_ID, SNAPSHOT) 
    ! ** WRITE FOR TRK VARIABLES
    INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
    INTEGER :: NP,STAT
    REAL(8)   :: XM(1),YM(1),XLL(1),YLL(1)
  
    TRK_LON = MISSING_VALUE
    TRK_LAT = MISSING_VALUE
    TRK_LEV = MISSING_VALUE
    TRK_VOL = MISSING_VALUE
 
    DO NP=1,NPD
      IF (LLA(NP) < 2 .OR. LLA(NP) > LA) CYCLE
      XM(1)=XLA(NP)
      YM(1)=YLA(NP)
      CALL UTMR_WGS84(XM,YM,XLL,YLL)
      TRK_LON(NP) = XLL(1)
      TRK_LAT(NP) = YLL(1)
      TRK_LEV(NP) = ZLA(NP)
      IF (ANY(ISOILSPI == 1)) TRK_VOL(NP) = DVOL(NP)
    ENDDO
    
    STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_XC,SNAPSHOT,TRK_LON) 
    STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_YC,SNAPSHOT,TRK_LAT)
    STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_SIG,SNAPSHOT,TRK_LEV)
    IF (ANY(ISOILSPI == 1)) STAT=NC_WRITE_SCALAR1D(NC_ID,NC_TRK_VOL,SNAPSHOT,TRK_VOL)   
    NC_WRITE_LPT=STAT
    
  END FUNCTION

  INTEGER FUNCTION NC_WRITE_MOC(NC_ID, SNAPSHOT)
    ! ** WRITE CONCENTRATION FOR OIL
    INTEGER,INTENT(IN) :: NC_ID, SNAPSHOT
    REAL(4), ALLOCATABLE, DIMENSION(:,:) :: C_OIL
    INTEGER :: STAT, I,J,L
  
    ALLOCATE(C_OIL(COLCNT,ROWCNT))
    C_OIL = MISSING_VALUE
    DO L = 2,LA
        I = IL(L) - IMN + 1
        J = JL(L) - JMN + 1
        IF (HP(L) > HDRY) THEN  
            C_OIL(J,I) = MOC(L)
        ENDIF
    ENDDO
    STAT=NC_WRITE_SCALAR2D(NC_ID,NC_IDX(41),C_OIL,SNAPSHOT)
    DEALLOCATE(C_OIL)
    NC_WRITE_MOC=STAT
  END FUNCTION
  
  SUBROUTINE NC_NEW_FILE(NC_ID,FILENAME)
    INTEGER(4) :: NC_ID
    CHARACTER(*),INTENT(IN) :: FILENAME
    INTEGER(4) :: STATUS

    STATUS = NC_CREATE_FILE(NC_ID,FILENAME) 
    STATUS = NC_WRITE_GRID(NC_ID,XCR(1:4,JMN:JMX,IMN:IMX),YCR(1:4,JMN:JMX,IMN:IMX))
    STATUS = NC_WRITE_ZBOT(NC_ID)     
    NTI = 0
  END SUBROUTINE 

  SUBROUTINE NC_WRITE_SNAPSHOT(NC_ID, TIMESEC)

    DOUBLE PRECISION,INTENT(IN) :: TIMESEC
    INTEGER(4),INTENT(IN) :: NC_ID
    INTEGER :: STATUS, N
    REAL(8) :: TI(1)
  
    NTI = NTI + 1
    TI = TIMEDAY
    WRITE(*,'(A22,F15.3)') 'NETCDF OUTPUT @ HOURS: ',TI(1)*24
    STATUS = NF90_PUT_VAR(NC_ID, NC_TIM, TI, (/NTI/), (/1/))
    STATUS = NC_WRITE_WSEL(NC_ID, NTI)
    STATUS = NC_WRITE_VEL(NC_ID, NTI)
    IF (ISNCDF(10) == 1) STATUS = NC_WRITE_SHEAR(NC_ID, NTI)
    IF (ISNCDF(11) == 1) STATUS = NC_WRITE_WIND(NC_ID, NTI)
    IF (ISNCDF(12) == 1) STATUS = NC_WRITE_WAVE(NC_ID, NTI)
	   
    IF (ISNCDF(1) == 1) STATUS = NC_WRITE_CONS(NC_ID, NTI, 12, SAL)
    IF (ISNCDF(2) == 1) STATUS = NC_WRITE_CONS(NC_ID, NTI, 13, TEM)
    IF (ISNCDF(3) == 1) STATUS = NC_WRITE_CONS(NC_ID, NTI, 14, DYE)
    IF (ISNCDF(4) == 1) STATUS = NC_WRITE_CONS(NC_ID, NTI, 15, SFL)	  
    IF (ISNCDF(5) == 1) STATUS = NC_WRITE_CONM(NC_ID, NTI, 16, TOX, NTOX)
    IF (ISNCDF(6) == 1) STATUS = NC_WRITE_CONM(NC_ID, NTI, 17, SED, NSED)
    IF (ISNCDF(7) == 1) STATUS = NC_WRITE_CONM(NC_ID, NTI, 18, SND, NSND)
    
    IF (ISNCDF(8) == 1) THEN
      DO N=1,NWQV
        IF(ISTRWQ(N) > 0) STATUS = NC_WRITE_CONS(NC_ID, NTI, 18 + N, WQV(:,:,N))
      ENDDO
    ENDIF
  
    ! ** LPT
    IF(ISNCDF(9) == 1) THEN     
      IF (ANY(ISOILSPI == 1)) STATUS = NC_WRITE_MOC(NC_ID, NTI)  
      STATUS = NC_WRITE_LPT(NC_ID, NTI) 
    ENDIF
            
  END SUBROUTINE 
  
  SUBROUTINE NETCDF_WRITE(NC_ID)
  
    INTEGER(4),INTENT(IN) :: NC_ID

    REAL(8) :: TIMEHOUR
    
    CHARACTER(10) :: SDATE,DATESTAMP
    CHARACTER(80) :: FILENAME
    INTEGER, SAVE :: IFRST = 1
    
    IF (ISSGLFIL == 1) THEN
      IF (IFRST == 1) THEN
        FILENAME = './/#output//DSI_'//TRIM(PROJ)//'.nc'
        CALL NC_NEW_FILE(NC_ID, FILENAME)
        IFRST = 0
      ENDIF    
      
    ELSEIF (ISSGLFIL == 0) THEN
      TIMEHOUR = TIMEDAY*24
      CALL TOGREGOR(SDATE,TIMEDAY)
      DATESTAMP = SDATE(1:4)//'_'//SDATE(5:6)//'_'//SDATE(7:8)
      IF((NC_DATESTAMP /= DATESTAMP)) THEN
        IF (IFRST == 0) THEN
          CALL NC_WRITE_SNAPSHOT(NC_ID, TIMESEC)
          CALL NC_CLOSE_FILE(NC_ID)  
          IF (TIMEDAY >= TENDNCDF) RETURN
        ENDIF
        FILENAME = './/#output//DSI_'//TRIM(PROJ)//'_'//DATESTAMP//'.nc'
        CALL NC_NEW_FILE(NC_ID, FILENAME)
        IFRST = 0
      ENDIF
      NC_DATESTAMP = DATESTAMP
    ENDIF 
    
    CALL NC_WRITE_SNAPSHOT(NC_ID, TIMESEC)
    
  END SUBROUTINE 
#endif
END MODULE
