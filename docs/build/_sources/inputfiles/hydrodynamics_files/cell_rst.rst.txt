.. _cellinp:

Cell Input File
---------------

The cell.inp file is a 2x2 matrix with length  in the `i` or `x` direction equal to `IC` and a length in the `j` or `y` direction of `JC`.  IC and JC are specified on :ref:`card9` of EFDC.INP

 In the table below each cell type is described.  These numbers are what are inputted into the ICxJC matrix in the cell.inp file.

+------------+-------------------------------------------------------------+
|Cell Number | Description                                                 |
+------------+-------------------------------------------------------------+
|0           |dry land cell not bordering a water cell on a side or corner.|
+------------+-------------------------------------------------------------+
|1           |triangular water cell with land to the northeast             |
+------------+-------------------------------------------------------------+
|2           |triangular water cell with land to the southeast             |
+------------+-------------------------------------------------------------+
|3           |triangular water cell with land to the southwest             |
+------------+-------------------------------------------------------------+
|4           |triangular water cell with land to the northwest             |
+------------+-------------------------------------------------------------+
|5           |quadrilateral water cell                                     |
+------------+-------------------------------------------------------------+
|9           |dry land cell bordering a water cell on a side or corner or a|
|            |fictitious dry land cell bordering an open boundary water    |
|            |cell on a side or a corner                                   |
+------------+-------------------------------------------------------------+

In the example file listed below IC=10 and JC=6.
Note, the first 4 rows are comments as well as the first 4 rows.  The cell mapping begins from the bottom left corner. 

.. literalinclude :: cell.inp
