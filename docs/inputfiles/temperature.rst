.. _temperature:

====================
Temperature Module
====================

============     ================================================================================
Input File       Description
============     ================================================================================
TEMP.INP         Water column initial conditions for temperature
TSER.INP         Time series file for temperature boundary conditions
ASER.INP         Time series file for atmospheric parameters
ATMMAP.INP       Cell weightings file for ASER series when NASER > 1
PSHADE.INP       Spatially varying solar radiation shading 
SVHTFACT.INP     Spatially varying surface heat exchange parameters for DSI full heat balance if ISVHEAT > 0
TEMB.INP         Spatially varying initial bed temperature and bed thermal thickness
============     ================================================================================



Ice Sub-Module
==============

============     ================================================================================
Input File       Description
============     ================================================================================
ISER.INP         Time series of user specified ice cover for ISICE = 1
ICEMAP.INP       Cell weightings file for ISER series when NISER > 1 for ISICE = 1
ISTAT.INP        Time series of user specified ice for whole domain when ISICE = 2
ICE.INP          Initial conditions for ice cover when using heat coupled ice (ISICE > 2)
============     ================================================================================






