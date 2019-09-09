.. _dxdy_rst:

DXDY Input File
---------------

The dxdy.inp file specifies many of the physical properties of each cell.  Every cell described in the cell.inp file must be described in this file.  

============ =====================================
Variable     Description
============ =====================================
I            Array index in x direction 
J            Array index in y direction 
DX           Cell dimension in x direction, meters 
DY           Cell dimension in y direction, meters 
DEPTH        Initial water depth, meters 
BOTTOM ELEV  Bottom bed elevation, meters 
ZROUGH       Log law roughness height, zo, meters 
VEG TYPE     Vegetation type class, integer value 
============ =====================================

Below is part of a sample input file that specifies part of a single column of the geometry.

.. literalinclude :: dxdy.inp