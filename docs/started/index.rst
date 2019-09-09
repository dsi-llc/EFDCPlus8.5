.. _started :


===============
Getting Started
===============

@todo put directions on how to obtain the EFDC+ source code, since we do not know exactly how this will work we will have to go back and fill this in.
Likely, we can have this all on a git repository, thus we can provide directions here on how to clone the git repo.


After cloning the EFDC+ repository the following folders will be available under the root directory.

``EFDC`` - Contains source code to build EFDC+, sample executables for different build options.

``NetCDFLib`` - Necessary library files for building EFDC+ so it can write NetCDF files out.

``NetCDFDLL4.4.0`` - @todo are these for 64 bit compilation? 

``GridGenerator`` - Contains the executable for the simple Grid Generator for EFDC+

``GetEFDC`` - Contains source code for building utility that helps extract EFDC+ formatted binary time series data.

``WASP`` - Provides some files necessary for linkage with the WASP code (advanced used feature).

``SampleModels/`` - Contains several sample EFDC+ models.

``docs`` - contains documentation for EFDC+.

The contents of each folder is described below:

.. code-block:: none

	EFDC/ - 
		|-- DebugDP64/ 
			`-- *.dll
		|-- DebugSP/
			`-- *.dll
		|-- DebugSP64/
			`-- EFDC.exe
			`-- *.dll
		|-- ReleaseDP64/
			`--EFDC.exe
			`--*.dll
		|-- ReleaseSP/
			`-- EFDC.exe
			`-- *.dll
		|-- ReleaseSP64/
			`-- EFDC.exe
			`-- *.dll
		`-- *.f90

	NetCDFLib/
		|--	include/
			`-- C++ header files
		|--	lib/
			`-- *.lib 

	NetCDFDLL4.4.0 
		`--	*.dll
		`--	*.exe

	GridGenerator
		`--GridGenerator.exe
		`--	*.dll

	GetEFDC/
		|-- src/
			`-- *.f90

	WASP/
	
	SampleModels/
		|-- Ohio_River_4
		|-- Lake_T_HYR_WQ
		|-- Lake_Washington

	docs/ 
		|-- build/
			`--*.rst
		|-- gridgen/
			`--*.rst
		|-- images/
			`--*.rst
		|-- inputfiles/
			`--*.rst
		|-- outputfiles/
			`--*.rst
		|-- samplemodels/
			`--*.rst
		|-- started/
			`--*.rst
		|-- _images/
		`-- conf.py - Sets up configuration for Sphinx documentation builder
		`-- index.rst - Root RST file
		`-- license.rst

.. toctree::
    :numbered:
    :maxdepth: 1

    build
    running