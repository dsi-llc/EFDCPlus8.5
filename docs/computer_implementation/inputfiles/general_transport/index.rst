.. _generaltransport:

=================
General Transport
=================


Hydrodynamic Parameter Files
============================

=============    =====================================================================================
Input File       Description
=============    =====================================================================================
AHMAP.INP        Spatially varying Smagorinsky (AHD) and background eddy viscosity (AHO) if AHD < 0.0
AVMAP.INP        Spatially varying AVO/ABO if AVO < 0.0
MAPHMD.INP       List of cells to compute horizontal momentum diffusion if IHMDSUB > 0
VEGE.INP         Vegetation class definitions
VEGSER.INP       Vegetation class time series
WSER.INP         Time series file for wind speed and direction
WNDMAP.INP       Cell weightings file for WSER series when NWSER > 1
SUBSET.INP       List of cells and timing for high frequency time series output 
SNAPSHOTS.INP    List of additional times to write the EE_*.OUT linkage files
RESTART.INP      Primary restart/hot start file for hydrodynamics and most other modules
RSTWD.INP        Restart file for wetting & drying parameters (ISDRY > 0)
=============    =====================================================================================
        


Volumetric and Level Boundary Conditions
========================================

============     ================================================================================
Input File       Description
============     ================================================================================
QSER.INP         Time series file for flow type boundary conditions
PSER.INP         Time series file for pressure type open boundary conditions
QWRS.INP         Time series file for withdrawl-return flows and concentration rise/fall
QCTL.INP         Lookup tables for free surface elevation or pressure controlled flow
QCTLSER.INP      Time series of equation based parameters control time-series
QCRULES.INP      Operation rules for hydraulic structure control
GWATER.INP       Groundwater interaction by infiltration and evapotranspiration
GWSEEP.INP       Groundwater interaction by ambient groundwater flow
GWSER.INP        Groundwater inflow/outflow and concentration
GWMAP.INP        Spatially varying map of GSWER series ID
============     ================================================================================

.. toctree::
    :maxdepth: 1

    salinity
    lagrangian
    shellfish
    dyes
