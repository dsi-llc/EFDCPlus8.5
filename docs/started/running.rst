.. _running :


Running 
=======

Running EFDC+ on a local computer can be done in several different ways.  By far the simplest way is to use the EFDC Explorer.  The EFDC Explorer is a Graphical User Interface (GUI) that allows EFDC+ models to be visualized and run with ease. Details on the EFDC Explorer can be found `here <https://www.eemodelingsystem.com/>`_.   Alternatively, the EFDC executable can be run through the command prompt or using a batch script. 

For this discussion it is assumed the user is not using the EFDC Explorer.

Execution Options
^^^^^^^^^^^^^^^^^

* Serial: EFDC+ can be run simply on a single core of a user's desktop.  No special compilation is required for this option. 
* Multithreaded: If the executable was compiled with OpenMP libraries, multiple cores on a single machine may be utilized.


Multithreading with OpenMP
--------------------------

With the advent of Intel processors with multiple cores, various technologies have been employed to take advantage of this increased computational power.  With DSI’s knowledge of the EFDC+ code and its structure, the approach selected to apply multi-threading to EFDC was OpenMP.

Additionally, Intel has implemented Hyper-Threading for their multi-core implementation.  This allows a program to utilize two threads for each physical core.

.. tip ::

    Since two threads share a single core, for computationally intensive applications, only using the number of threads equal to the number of cores provides the best throughput


OpenMP Performance
~~~~~~~~~~~~~~~~~~

The performance of the different routines within EFDC+ using OpenMP is highlighted in the figure below.

|omp_perf|

It is DSI's commitment that the EFDC_DSI_OMP models produce exactly the same results regardless of how many threads are used.  Model comparisons demonstrate model differences are equal to zero  i.e. model results are exactly the same, within model precision.

Executing and EFDC+ Run
^^^^^^^^^^^^^^^^^^^^^^^

Once an EFDC+ executable is available, it can be run directly through the command prompt or through simple batch script.  A sample of batch script is given below:

.. literalinclude :: ./sample_script.txt

Each line the script is described in greater detail below. 

* SET KMP_AFFINITY=granularity=fine,compact,1,0 - Specifies an environment variable that binds OpenMP threads to physical processing units. This generally gives the best performance.  For additional information on what this environment variable is doing, go to `this article <https://software.intel.com/en-us/articles/using-kmp-affinity-to-create-openmp-thread-mapping-to-os-proc-ids>`_, written by Intel.

* TITLE: A title is optional but can be helpful if you are running multiple calculations at the same time.  This title will show up at the top of the command prompt.

* CD: This command precedes the location of your `Working Directory`, which is the folder containing your EFDC+ inputs. 

* The path to the location of the executable must be known and specified in the script.

* -NT2: This command line argument after the executable specifies the number of threads to use.  In this case only 2 are requested, but one could easily specify 3,4,5, etc. to run more threads.  

.. important ::
	
	Do not specify a number of threads greater than the number of logical cores available on your computer.


Output Screen
^^^^^^^^^^^^^

The EFDC+ end of run screen contains the CPU usage for the run.  The CPU usage is reported as the Total CPU/# Cores. This statistic is best for interpreting impacts to run times.  A sample of an output screen is given below.


|eor_screen|


.. |omp_perf| image:: images/omp_performance.png
.. |eor_screen| image:: images/eor_screen.png