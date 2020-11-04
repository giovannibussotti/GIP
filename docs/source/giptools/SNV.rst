###
SNV
###

Purpose
-------
The ``SNV`` module aims at comparing the SNV in terms of number, position and frequency for a samples set.


Algorithm
---------

For each sample the module loads the GIP files with the filtered SNVs data (singleVariants.df.gz files) and generates multiple plots overlaying in different colors the SNVs sets. 


Options
-------

+-----------------------+--------------------------------------------------------------+----------------+
|Option                 |Description                                                   |Argument        |
+=======================+==============================================================+================+
|\-\-samples            |Sample names. It determines the plotting order [**required**] |[char ...]      |     
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-gipOut             |GIP output directory [**required**]                           |[char]          |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-outName            |Output name [default NA]                                      |[char]          |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-chrs               |Chromosomes to use. If "NA" it uses the same chromsomes as GIP|[char ...]      |
|                       |                                                              |                |
|                       |[default NA]                                                  |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-coordCartesianYlim |Max y-axis value for density plots                            |[double]        |
|                       |                                                              |                |
|                       |If \"NA\" the value is automatically attributed [default NA]  |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-minVRF             |Discard SNVs with frequency < --minVRF [default 0]            |[double]        |
+-----------------------+--------------------------------------------------------------+----------------+  
|\-\-debug              |Dump session and quit                                         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-h, \-\-help          |Show help message                                             |                |
+-----------------------+--------------------------------------------------------------+----------------+



Output
------





Example
-------



