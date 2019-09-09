.. _lxly_rst:

LXLY Input File
---------------

The lxly.inp file specifies the cell centered location and the rotation of each cell.

============ ====================================================
Variable     Description
============ ====================================================
I            Array index in x direction 
J            Array index in y direction 
X            x cell center coordinate, longitude, meters, or km         
Y            y cell center coordinate, longitude, meters,  or km
CUE          Rotation matrix component, i 
CVE          Rotation matrix component 
CUN          Rotation matrix component 
CVN          Rotation matrix component 
Wind Shelter 
============ ====================================================

Below is part of a sample input file that specifies part of a single column of the geometry.

.. literalinclude :: lxly.inp