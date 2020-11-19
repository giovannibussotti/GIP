#######
ternary
#######


Options
-------

+-----------------------+--------------------------------------------------------------+----------------+
|Option                 |Description                                                   |Argument        |
+=======================+==============================================================+================+
|\-\-samples            |Sample names. It determines the plotting order [**required**] |[char ...]      |     
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-gipOut             |GIP output directory.                                         |[char]          |
|                       |                                                              |                |
|                       |If "NA" the directory "./gipOut" is used [default NA]         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-outName            |Output name [default NA]                                      |[char]          |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-chrs               |Chromosomes to use. If "NA" it uses the same chromsomes as GIP|[char ...]      |
|                       |                                                              |                |
|                       |[default NA]                                                  |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-MAPQ               |Label genes with MAPQ < --MAPQ [default 0]                    |[int]           |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-pseudocount        |Normalized mean coverage                                      |[double]        |
|                       |                                                              |                |
|                       |pseudocount value (for plots only)  [default 0.1]             |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-colorByDelta       | Color genes by increasing delta                              |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-highLowDeltaColor  | Colors for high and low delta.                               |[char char]     |
|                       |                                                              |                |
|                       | DEPENDENCY \-\-colorByDelta  [default darkgreen black]       |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-showDensity        | Show log 10 scaled density area                              |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-highLowDensityColor| Colors for high and low density                              |[char char]     |
|                       |                                                              |                |
|                       | DEPENDENCY \-\-showDensity. [default black deepskyblue]      |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-showQuantile       |In density plot show genes/clusters over this                 |[double]        | 
|                       |                                                              |                |
|                       |quantile cut-off. DEPENDENCY --showDensity [default 0.99]     |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-densityN           |Number of grid points in each direction                       |[int]           | 
|                       |                                                              |                |
|                       |DEPENDENCY --showDensity [default 200]                        |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-densityBins        |Number bins over which to calculate intervals                 |[int]           | 
|                       |                                                              |                |
|                       |DEPENDENCY --showDensity [default 8000]                       |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-plot3dMaxCOV       | 3-D scatteplot visualization threshold.                      |[double]        |
|                       |                                                              |                |
|                       | Gene/cluster coverage values greather than this threshold    |                |
|                       |                                                              |                |
|                       | are shown as \-\-plot3dMaxCOV  [default 3]                   |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-plot3dMaxFC        | 3-D scatteplot visualization threshold. Gene/cluster coverage|[double]        |
|                       |                                                              |                |
|                       | fold change values greather than this threshold              |                |
|                       |                                                              |                |
|                       | are shown as \-\-plot3dMaxFC [default 3]                     |                |
+-----------------------+--------------------------------------------------------------+----------------+  
|\-\-debug              |Dump session and quit                                         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-h, \-\-help          |Show help message                                             |                |
+-----------------------+--------------------------------------------------------------+----------------+



Description
-----------
| The ``ternary`` module aims at comparing the gene sequencing coverage of 3 samples to identify gene CNVs.
| The module loads for the three samples the GIP files with the gene sequencing coverage values (.covPerGe.gz files) and generates a ternary diagram of the normalized coverage values. In this representation, the values coverage values in the 3 samples sum to a constant represented for convenience as 100%. Additionally, the module generates a 3D-scatterplot demonstrating the normalized gene coverage in the 3 samples.
 

Output
------





Example
-------
