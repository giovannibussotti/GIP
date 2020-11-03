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

The module loads for the three samples the GIP files with the gene sequencing coverage values (.covPerGe.gz files) and generates a ternary diagram of the normalized coverage values. 

In this representation, the values coverage values in the 3 samples sum to a constant represented for convenience as 100%.

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
|\-\-highLowRatio   |Provide 2 numbers. Genes with ratio scores > num1                 |[double,double] |
|                   |                                                                  |                |
|                   | or < num2 will be colored differently [default (1.5 , 0.5)]      |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-pseudocount    |Normalized mean coverage pseudo count value preventing            |[double]        |
|                   |                                                                  |                |
|                   |  infinite (1/0) and NaN (0/0) ratio values [default 0.1]         |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-plot1_ylim     |Plot1 visualization threshold. Gene ratio values greather         |[double]        |
|                   |                                                                  |                | 
|                   | than this threshold are shown as --ylim  [default 5]             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-plot3_min      |Plot3 visualization threshold. Min normalized gene coverage       |[double]        |
|                   |                                                                  |                |
|                   |DEPENDENCY:--scaleFree no [default 0]                             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-plot3_max      |Plot3 visualization threshold. Max normalized gene coverage       |[double]        |
|                   |                                                                  |                |
|                   |DEPENDENCY:--scaleFree no [default 100]                           |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-plot24_min     |Plots 2 and 4 visualization threshold. Min normalized gene        |[double]        |
|                   |                                                                  |                |
|                   |coverage (log10 scale). DEPENDENCY:--scaleFree no [default -1]    |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-plot24_max     |Plots 2 and 4 visualization threshold. Max normalized gene        |[double]        |
|                   |                                                                  |                |
|                   |coverage (log10 scale). DEPENDENCY:--scaleFree no [default 3]     |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-scaleFree      | Graphical parameter plots 3 and 4.                               |[yes|no]        |
|                   |                                                                  |                |
|                   | Set scale free axes [default yes]                                |                |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-debug          |Dump session and quit                                             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-h, \-\-help      |Show help message                                                 |                |
+-------------------+------------------------------------------------------------------+----------------+




Output
------





Example
-------
