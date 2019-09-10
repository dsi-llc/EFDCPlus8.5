.. _output:

============
Output Files 
============

These output files are written out by EFDC+ in a binary format. The easiest way to view the results is using EE Modeling System (EEMS).  A demo of EEMS is available and can be accessed by going to the `EEMS website <https://www.eemodelingsystem.com/buy/demo-version>`_.  Alternatively, a rudimentary postprocessing tool is available, referred to as `GetEFDC`.  `GetEFDC` is a Fortran 90 program that can read the binary formats and can be modified to output into another format like a text file.  A detailed description of `GetEFDC` is found on the next page.

In the table below each of the binary output files is listed and described.   

.. tip:

    All output is written to the folder ``#output``.

================ =========================================================
Output file name Description
================ =========================================================
EE_WS.OUT        Water depth
EE_WC.OUT        Water column and top layer of sediments
EE_BC.OUT        Computed boundary flows
EE_BED.OUT       Sediment bed layer information
EE_WQ.OUT        Water quality information for the water column
EE_SD.OUT        Sediment diagensis information
EE_RPEM.OUT      Rooted plant and epiphyte model
EE_SEDZLJ.OUT    Sediment bed data for sedzlj sub-model
EE_HYD.OUT       Water depth and velocity
================ =========================================================

Some of these outputs are optional based on the options set in the efdc.inp file.
