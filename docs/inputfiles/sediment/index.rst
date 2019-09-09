.. _sediment:

========
Sediment 
========


Original Sediment Module
========================

================ ====================================================================================================
Input File       Description
================ ====================================================================================================
SEDW.INP         Water column initial conditions for cohesive sediments
SEDB.INP         Sediment bed initial conditions for cohesive sediments
SDSER.INP        Time series file for cohesive boundary conditions
SNDW.INP         Water column initial conditions for non-cohesive sediments
SNDB.INP         Sediment bed initial conditions for non-cohesive sediments
SNSER.INP        Time series file for non-cohesive boundary conditions
BEDBDN.INP       Sediment bed initial conditions for bulk density
BEDDDN.INP       Sediment bed initial conditions for dry density, porosity or void ratio
BEDLAY.INP       Sediment bed initial conditions for layer thickness
SEDBLBC.INP      Non-cohesive bedload outflow or recirculation boundary conditions
SEDROUGH.INP     Spatially varying grain roughness height for determining grain stress, ISBEDSTR = 3
CONSOLMAP.INP    Spatially varying consolidation approach when IBMECH = 9
SSCOHSEDPMAP.INP Spatially varying cohesive critical bed shear stress and surface erosion rate, IWRSP(1) > 98
BEDMAP.INP       Spatially varying flag for hard-bottom bypass of erosion/deposition calculations
BEMAP.INP        Bank erosion cell map
BESER.INP        Bank erosion time series
================ ====================================================================================================


SEDZLJ Module
=============

===============  =============================================================================================
Input File       Description
===============  =============================================================================================
BED.SDF          SEDZLJ control file with active and deposited erosion parameters
ERATE.SDF        SEDFlume core properties for existing sediment bed
CORE_FIELD.SDF   Sptially varying assignment of core ID's from ERATE.SDF
SEDW.INP         Water column initial conditions for cohesive sediments
SDSER.INP        Time series file for cohesive boundary conditions
SNSER.INP        Time series file for non-cohesive boundary conditions
SNDW.INP         Water column initial conditions for non-cohesive sediments
SEDB.INP         Sediment bed initial conditions for cohesive sediments NSEDFLUME = 2 (10.0)
SNDB.INP         Sediment bed initial conditions for non-cohesive sediments NSEDFLUME = 2 (10.0)
BEDBDN.INP       Sediment bed initial conditions for bulk density NSEDFLUME = 2 (10.0)
BEDDDN.INP       Sediment bed initial conditions for dry density, porosity or void ratio NSEDFLUME = 2 (10.0)
BEDLAY.INP       Sediment bed initial conditions for layer thickness NSEDFLUME = 2 (10.0)
SEDBLBC.INP      Non-cohesive bedload outflow or recirculation boundary conditions
SEDBED_HOT.SDF   Restart/hot start file of SEDZLJ sediment conditions
===============  =============================================================================================
