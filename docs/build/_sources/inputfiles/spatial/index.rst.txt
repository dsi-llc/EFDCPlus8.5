.. _spatial:


======================
Required Spatial Files
======================

============     ================================================================================
Input File       Description
============     ================================================================================
cell.inp         Describes the cell mapping and which type of cell goes where.
celllt.inp       Auxiliary cell type file
dxdy.inp         Horizontal cell dimensions, depth, bottom elevation, roughness, vegetation class
lxly.inp         Horizontal cell center coordinates and cell orientation
corners.inp      Provides x,y coordinates corners for Lagrangian Particle Tracking module
============     ================================================================================


Optional Spatial Files
======================

==============   ===================================================================================
Input File       Description
==============   ===================================================================================
mask.inp         Specifies thin barriers if NMASK > 0
layermask.inp    Specifies thin barriers for layer faces if NBLOCKED > 0  (for EFDC+ 10.1 and later)
mappgns.inp      Specifies north-south (J/V face) direction grid connections
mappgew.inp      Specifies east-west (I/U face) direction grid connections
moddxdy.inp      Modifies cell dimensions originally specified in dxdy.inp
sgzlayer.inp     Specifies the bottom active water layer if IGRIDV=1
==============   ===================================================================================


The primary input files that specify the geometry of the problem are given in greater detail below.
 
.. toctree:: 
    :maxdepth: 1

    cell_rst
    dxdy_rst
    lxly_rst



