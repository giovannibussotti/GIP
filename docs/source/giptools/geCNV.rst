#####
geCNV
#####


Purpose
-------
The ``geCNV`` module aims at comparing the gene sequencing coverage of 2 samples to identify gene CNVs.


Algorithm
---------

The module loads for the two samples the GIP files with the gene sequencing coverage values (.covPerGe.gz files) and calculates for each gene the ratio of the normalized coverage values. 


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
|\-\-pseudocount    |Normalized mean coverage pseudo count value preventing           |[double]        |
|                   |                                                                  |                |
|                   |  infinite (1/0) and NaN (0/0) ratio values [default 0.1]         |                |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-debug          |Dump session and quit                                             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-h, \-\-help      |Show help message                                                 |                |
+-------------------+------------------------------------------------------------------+----------------+


TO ADD:
+-------------------+------------------------------------------------------------------+----------------+
|\-\-plot1_ylim     |Plot1 visualization threshold. Gene ratio values greather         |[double]        |
|                   |                                                                  |                | 
|                   | than this threshold are shown as --ylim  [default 5]             |                |
+-------------------+------------------------------------------------------------------+----------------+
\-\-plot3_min"  , type="double" , help="Plot3 visualization threshold. Min normalized gene coverage. DEPENDENCY:--scaleFree no [default %(default)s]" , default=0)
+-------------------+------------------------------------------------------------------+----------------+
\-\-plot3_max"  , type="double" , help="Plot3 visualization threshold. Max normalized gene coverage. DEPENDENCY:--scaleFree no [default %(default)s]" , default=100)
+-------------------+------------------------------------------------------------------+----------------+
\-\-plot24_min" , type="double" , help="Plots 2 and 4 visualization threshold. Min normalized gene coverage (log10 scale). DEPENDENCY:--scaleFree no [default %(default)s]" , default=-1)
+-------------------+------------------------------------------------------------------+----------------+
\-\-plot24_max" , type="double" , help="Plots 2 and 4 visualization threshold. Max normalized gene coverage (log10 scale). DEPENDENCY:--scaleFree no [default %(default)s]" , default=3)
+-------------------+------------------------------------------------------------------+----------------+
\-\-scaleFree"  , help="Graphical parameter plots 3 and 4. Set scale free axes [yes|no] [default %(default)s]" , default="yes")
+-------------------+------------------------------------------------------------------+----------------+

Output
------





Example
-------
