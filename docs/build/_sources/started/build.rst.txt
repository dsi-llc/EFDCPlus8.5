.. _build:


Build Instructions
==================

Building EFDC+ from the source code is most easily carried out with Visual Studio (VS).  The Intel Fortran compiler is the preferred compiler and the most tested for building EFDC+.  The simplest build method is to open the VS solution file that comes with the source code.  The solution file is located in the root directory of the ``EFDC`` folder and will have the '.sln' extension.  The following discussion assumes the user is running Windows 7 (or greater) and has access to Visual Studio 2015 (or greater).  These tools have been most throughly tested.  However, it is likely that EFDC+ can be built with other versions of the Intel compiler and Visual Studio. 

In addition to the source code, pre-built executables are available under each of the following folders:

* EFDC/DebugSP64/

* EFDC/ReleaseDP64/

* EFDC/ReleaseSP/

* EFDC/ReleaseSP64/

Next, the differences in each of the build configurations will be explained.

Build Configurations
--------------------

The different build configurations are managed by VS.  Visual studio provides a convenient way for maintaining different build configurations for the same project.  The build configurations can be inspected by clicking the `Project` tab in the top of VS and selecting `properties`.  Each of the build configurations is listed below:


* DEBUG SP 

* DEBUG SP 64 

* DEBUG DP 

* DEBUG DP 64

* Release SP 

* Release SP 64

* Release DP 

* Release DP 64


For each of the bulleted build configurations listed above, if `64` is not specified, the executables is assumed to be compiled for a  32 bit system.  The table below explains the shorthand used to signify the differences in build configurations.

===  ==========================
SP   Single Precision
DP   Double Precision
64   64 bit compilation
===  ==========================

OpenMP Compilation
^^^^^^^^^^^^^^^^^^

Compilation of EFDC+ with OpenMP allows multithreading, which typical results in a reduction of the total calculation time.  The build configuration requires specifying several things in the VS `Properties` page.  These settings are already configured in the builds provided.  However, the details are given below in case a user wants to make modifications or use an OpenMP library besides Intel's. 

**Under: Fortran\\Preprocessor**

=============================== ===============================
OpenMP Conditional Compilation  Yes
=============================== ===============================

**Under: Fortran\\Preprocessor**

=============================== =================================
Process OpenMP Directives		Generate Parallel Code (/Qopenmp)
=============================== =================================

**Under: Fortran\\Libraries**

=============================== ======================================
Runtime Library					Multithreaded DLL (/libs.dll /threads)
=============================== ======================================

**Under: Linker\\Input**

=============================== ==============
Additional Dependencies			libiomp5md.lib
=============================== ==============

Below is a summary of the Intel compiler suite and Visual Studio versions that are known to work.

Intel Compiler Versions Tested
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Intel 15 

* Intel 19.3

* Intel 19.4


Visual Studio Versions Tested
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
* 2015
* 2019 (Preview 4)
