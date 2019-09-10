.. _started :


===============
Getting Started
===============

The EFDC+ source code and associated utilities can be access by cloning the repository from the `EFDC+ GitHub repository <https://github.com/dsi-llc/EFDCPlus/>`_.

From the command line execute the following: 

::

    git clone https://github.com/dsi-llc/EFDCPlus

Alternatively, the repository may be downloaded by using the green `Clone or Download` button on the EFDC+ GitHub page.

After cloning the EFDC+ repository the folders listed below will be available under the root directory.

``EFDC`` - Contains source code to build EFDC+, sample executables for different build options.

``NetCDFLib`` - Necessary library files for building EFDC+ so it can write NetCDF files out.

``GridGenerator`` - Contains the executable for the simple Grid Generator for EFDC+

``GetEFDC`` - Contains source code for building utility that helps extract EFDC+ formatted binary time series data.

``WASP`` - Provides some files necessary for linkage with the WASP code (advanced user feature).

``SampleModels/`` - Contains several sample EFDC+ models.

``docs`` - Contains the computer implementation guide and the theory documentation for EFDC+.


.. toctree::
    :numbered:
    :maxdepth: 1

    build
    running