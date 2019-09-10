.. _runcontrol:

===================
Primary Run Control 
===================

The run control files contain options to specify calculation types, time step sizes, output options, and other related controls.  The most important of these files is the efdc.inp file and will be explained in detail later on.

============     ================================================================================
Input File       Description
============     ================================================================================
efdc.inp         Master EFDC+ control file
show.inp         Model run time reporting options
efdcwin.inp      Simplified control (deprecated)
============     ================================================================================


Restart  Related Files
^^^^^^^^^^^^^^^^^^^^^^

The restart files allow an EFDC+ calculation to start up from a specified time step part way through a calculation.  Use of these files may be necessary if a calculation ended prematurely and a user wishes to restart from the last saved time step.  

=============   =================================
Input File      File Description
=============   =================================
restart.inp     hydrodynamic restart
rstwd.inp	    wetting  & drying
temp.rst	    bed temperature restart
wqwcrst.inp	    water quality restart
wqsdrst.inp	    sediment diagenesis restart file
wqrpemrst.inp	rooted plant & epiphyte
=============   =================================

EFDC.INP
^^^^^^^^

The efdc.inp file is extensive and specifies all options for running a calculation. Historically, the options were organized by card types.  As such, each card description and input parameter is given below. Each of these cards is placed in a single efdc.inp file and read in by EFDC+ at run time.

Note, when creating and efdc.inp file a line starting with * or - will be ignored and interpretted as a comment.

.. toctree::
    :maxdepth: 1

    card1
    card2 
    card3
    card4 
    card5
    card6 
    card7   
    card8 
    card9 
    card10 
    card11 
    card11a 
    card11b 
    card12 
    card12a 
    card13 
    card14 
    card14a 
    card15 
    card16 
    card17 
    card18 
    card19 
    card20 
    card21 
    card22 
    card22b 
    card23 
    card24 
    card25 
    card26 
    card27 
    card28 
    card29 
    card30 
    card31 
    card32 
    card33 
    card34 
    card35 
    card36 
    card36a 
    card36b 
    card37 
    card38 
    card39 
    card40 
    card41 
    card42 
    card42a 
    card43a 
    card43b 
    card43c 
    card43d 
    card43e 
    card44 
    card45 
    card45a 
    card45b 
    card45c 
    card45d 
    card46 
    card46a 
    card46c 
    card46e 
    card47 
    card48 
    card49 
    card50 
    card51 
    card52 
    card53 
    card54 
    card55 
    card56 
    card57 
    card58 
    card59 
    card60 
    card61 
    card62 
    card63 
    card64 
    card65 
    card66 
    card66a 
    card66b 
    card67 
    card68 
    card69  
    card70 
    card71 
    card71a 
    card71b 
    card72 
    card73 
    card74 
    card75 
    card76 
    card77 
    card78 
    card79 
    card80 
    card81 
    card82 
    card83 
    card84 
    card85 
    card86 
    card87 
    card88 
    card89 
    card90 
    card91 
    card91a 
    card91b 