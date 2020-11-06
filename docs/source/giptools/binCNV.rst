######
binCNV
######

Purpose
-------
The ``binCNV`` module aims at comparing the bin sequencing coverage of 2 samples. This module is useful to identify intra-chromosomal CNV regions between 2 isolates.


Algorithm
---------

For each sample the module loads the GIP files with the bin sequencing coverage (.covPerBin.gz files) and calculates the ratio of the normalized coverage value between corresponding bins. 


Options
-------

+-------------------+------------------------------------------------------------------+----------------+
|Option             |Description                                                       |Argument        |
+===================+==================================================================+================+
|\-\-samples        |Sample names. It determines the plotting order [**required**]     |[char ...]      |
+-------------------+------------------------------------------------------------------+----------------+ 
|\-\-gipOut         |GIP output directory [**required**]                               |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-outName        |Output name [default NA]                                          |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-chrs           |Chromosomes to use. If "NA" it uses the same chromsomes as GIP    |[char ...]      |
|                   |                                                                  |                |
|                   |[default NA]                                                      |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-MAPQ           |Label bins with MAPQ < --MAPQ [default 0]                         |[int]           |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-ylim           |Plot visualization threshold. Bin ratio values greather than this |[double]        |
|                   |                                                                  |                |   
|                   |threshold are shown as --ylim [default 3]                         |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-highLowRatio   | Provide 2 numbers. Bins with ratio scores > num1                 |[double,double] |
|                   |                                                                  |                |
|                   | or < num2 will be colored differently [default (1.5 , 0.5)]      |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-pseudocount    | Normalized mean coverage pseudo count value preventing           |[double]        |
|                   |                                                                  |                |
|                   |  infinite (1/0) and NaN (0/0) ratio values [default 0.1]         |                |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-debug          |Dump session and quit                                             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-h, \-\-help      |Show help message                                                 |                |
+-------------------+------------------------------------------------------------------+----------------+



Output
------





Example
-------
