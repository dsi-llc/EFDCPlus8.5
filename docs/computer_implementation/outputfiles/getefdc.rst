.. _getefdc:

=======
GetEFDC
=======

`GetEFDC` is a Fortran utility to read the binary output files produced by EFDC+. `GetEFDC` is a  **starting point** for the user to create their own analysis of EFDC+ output.  It is expected that a user utilizing this tool has knowledge of Fortran and can modify `GetEFDC` to meet their specific needs. 


Source Code
-----------

The source code for `GetEFDC` is listed below.  It written all in Fortran and is straightforward to compile with a modern compiler (e.g. gfortran or Intel).

Main program:   

* getefdc.f90 
 
and 8 Modules: 

* infomod.f90 

* efdcpromod.f90

* tecmod.f90

* geteeoutmod.f90

* xyijconv.f90

* gethfreqout.f90

* globalvars.f90  
 
 
If IGRIDV > 0 then the output is based on the vertical layer defined in sgzlayer.inp.   


Build GetEFDC
-------------

A ``makefile`` is located under the /GetEFDC/src directory.  This can be used to compile on a Linux machine. Alternatively, this program can be compiled on Windows using Visual Studio.


Running GetEFDC
---------------
The syntax for running the utility is as follows:

.. code-block :: none

    GetEFDC.exe getefdc.inp


Input File
----------

``getefdc.inp`` is the master file that stores all the information about the parameters of
interest which the user is trying to extract. This file must be edited for every change to input
parameters. A sample of the master file is included in the GetEFDC folder. The input
parameters in this file are as follows:

* The full path of the folder containing ``efdc.inp`` file

* LAYK The number of layer in the vertical to get data for 2DH display (>0)

	- = 0 Extract the depth-averaged data

	- >0 Extract the data at layer of k
	- =-1 Extract High Frequency output
	- =-2 Extract data for time series (TS) at a height above bed (m)
	- =-3 Read TMP.DAT file and write an array data file for TECPLOT

* ZOPT This parameter is used in case of LAYK=-2

    - =1 Extract TS data at the depth under water surface

    - =2 Extract TS data at the height above bottom	

* JULTIME Julian time point for a selected layer, if > MAXTIME then JULTIME=MAXTIME
  JULTIME = 0 Extract data for all snapshots
* NLOC Number of locations (cells) to extract data. The location can be given
  as Index (I,J) or UTM coordinates (X,Y) via the parameter INDEX.
* ROTA The option for rotation of velocity components (U,V)

  - = 0  Extract (U,V) without rotation

  - = 1  (U,V) components are rotated to the true east and true north directions

* INDEX 

    - = 0 UTM (X,Y) of cells are used

    - = 1 Indices (I,J) of cells are used

* VPROF The option to extract data for vertical profile, 0 (No)/ 1(Yes)
* TECPLOT The option to extract data for 2DH Tecplot, 0 (No)/ 1(Yes)
* NDRIFTER: A successive set of number of particles to extract data for (X,Y,Z)
* I/X I Indices or X abscissa of cells to extract data
* J/Y J Indices or Y coordinates of cells to extract data

* ZINT Height under water surface or above bed to extract data in case LAYK=-2

Please note that the lines which start with `**` in the ``getefdc.inp`` file are comments and
will be ignored.

Sample GetEFDC Input File
-------------------------

.. literalinclude :: ./getefdc.inp

Additional sample GetEFDC input files are available for each of the sample models provided.  Each of these input files can be found under:

* ``EFDCPlus/SampleModels/Lake_T_HYD_WQ/GetEFDC/``

* ``EFDCPlus/SampleModels/Lake_Washington/GetEFDC/``

* ``EFDCPlus/SampleModels/Ohio_River_4/GetEFDC/``


Output Files
------------

After running GetEFDC a sub-folder ``RESULT`` is generated in the folder ``#output`` of the
working model. The extracted files are ASCII with the following conventions for the file
names:

* First characters group shows the constituent, such as SAL for salinity

* Second character group is `TSK_` which is the time series of the layer K, such as
  `TSK_4` is time series for the layer K=4

* The last character group is _DOM for the domain or CEL for the selected cells

* The vertical profiles for the constituents at the selected cells use the group _PROF
  in the file names, such as ``SAL_PROF.DAT``
