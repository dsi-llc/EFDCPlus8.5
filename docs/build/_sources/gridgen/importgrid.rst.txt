.. _importgrid:

Import Grids from Files
-----------------------

This option allows user to import an existing grid file. Grid file
formats that are supported include: 

- CVLGrid: DSI's curvilinear orthogonal grid generator - `CVL Grid information <https://www.eemodelingsystem.com/ee-modeling-system/cvlgrid/overview>`_

- RGFGrid: Deltares grid generator

- Grid95: @todo check on EEMS website

- DXDY/LYLY: EFDC+ grid descriptors

- ECOMSED: 

- SEAGRID: @todo add link to download page

- CH3D: Army Corp of Engineers model 

- Corners: file containing coordinates of four corners of each cell

The user should click the *Import Grids* radial button, and the *Import
Grid* form will be displayed. From here the user should select the grid
type from the drop-down list of *Grid types* as shown in Figure 12, then
click the *Browse* button to browse to the grid file, and click *OK*.

In the case that there are a number of sub-grids for a water body, the
*Multiple grid files* option needs to be checked, then the user may
browse to the folder containing the grid files. To select multiple grid
files at same time, hold the `Ctrl` key and select the grid files then
click the *OK* button to load grid files (see Figure 13 and Figure 14).

Save the EFDC model in the same way described in section :ref:`uniformgrid`.

|image11|

**Figure 12** Import from a grid file.

|image12|

**Figure 13.** Import from multiple grid files – file selection.

|image13|

**Figure 14.** Import from multiple grid files – result display.

Note that the user should update the \ *UTM Zone* to the correct value
before clicking the *OK* button to generate a new model. The *UTM Zone*
is not used for model computation but it is important for coordinate
conversion, exporting to GIS formats, and writing NetCDF outputs.

.. |image11| image:: media/image13.png
   :width: 6.00000in
   :height: 3.83093in
.. |image12| image:: media/image14.png
   :width: 6.00000in
   :height: 3.82051in
.. |image13| image:: media/image15.png
   :width: 6.00000in
   :height: 3.83093in