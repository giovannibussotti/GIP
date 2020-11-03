#######
ternary
#######

#####
geCNV
#####


Purpose
-------
The ``ternary`` module aims at comparing the gene sequencing coverage of 3 samples to identify gene CNVs.


Algorithm
---------

The module loads for the three samples the GIP files with the gene sequencing coverage values (.covPerGe.gz files) and generates a ternary diagram of the normalized coverage values. In this representation, the values coverage values in the 3 samples sum to a constant represented for convenience as 100%. Additionally, the module generates a 3D-scatterplot demonstrating the normalized gene coverage in the 3 samples.

Options
-------

+-------------------+------------------------------------------------------------------+----------------+
|Option             |Description                                                       |Argument        |
+===================+==================================================================+================+
|\-\-samples        |Sample names. It determines the plotting order [**required**]     |[char ...]      | |                   |                                                                  |                |    
+-------------------+------------------------------------------------------------------+----------------+
|\-\-gipOut         |GIP output directory [**required**]                               |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-outName        |Output name [default NA]                                          |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-chrs           |Chromosomes to use. If "NA" it uses the same chromsomes as GIP    |[char ...]      |
|                   |                                                                  |                |
|                   |[default NA]                                                      |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-MAPQ           |Label genes with MAPQ < --MAPQ [default 0]                        |[int]           |
+-------------------+------------------------------------------------------------------+----------------+



+-------------------+------------------------------------------------------------------+----------------+  
|\-\-debug          |Dump session and quit                                             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-h, \-\-help      |Show help message                                                 |                |
+-------------------+------------------------------------------------------------------+----------------+


TO ADD
+-------------------+------------------------------------------------------------------+----------------+
--colorByDelta" , action="store_true" , help="Color genes by increasing delta [default %(default)s]", default=FALSE)
+-------------------+------------------------------------------------------------------+----------------+
--highLowDeltaColor" , nargs="+", help="Colors for high and low delta. DEPENDENCY --colorByDelta [default %(default)s]", default=c("black","green"))
+-------------------+------------------------------------------------------------------+----------------+
--showDensity"  , action="store_true" , help="Show density area [default %(default)s]",--highLowDensityColor" , nargs="+", help="Colors for high and low density. DEPENDENCY --showDensity. [default %(default)s]", default=c("black","gray"))
+-------------------+------------------------------------------------------------------+----------------+
--plot3dMaxCOV" , type="double" , help="3-D scatteplot visualization threshold. Gene/cluster coverage values greather than this threshold are shown as --plot3dMaxCOV [default %(default)s]" , default=3)
+-------------------+------------------------------------------------------------------+----------------+
--plot3dMaxFC"  , type="double" , help="3-D scatteplot visualization threshold. Gene/cluster coverage fold change values greather than this threshold are shown as --plot3dMaxFC [default %(default)s]" , default=3)




Output
------





Example
-------
