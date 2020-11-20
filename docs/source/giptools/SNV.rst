###
SNV
###


Options
-------

+-----------------------+--------------------------------------------------------------+----------------+
|Option                 |Description                                                   |Argument        |
+=======================+==============================================================+================+
|\-\-samples            |Sample names (max 7). It determines the plotting order.       |[char ...]      |
|                       |                                                              |                |
|                       |If "NA" all samples are used [default NA]                     |                |
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
|\-\-coordCartesianYlim |Max y-axis value for density plots                            |[double]        |
|                       |                                                              |                |
|                       |If \"NA\" the value is automatically attributed [default NA]  |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-minVRF             |Discard SNVs with frequency < --minVRF [default 0]            |[double]        |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-showCoverage       |Overlay sequencing coverage                                   |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-smoothingFactor    |Merge and average the normalized coverage of this number      |[double]        |
|                       |                                                              |                |
|                       |of adjacent bins. DEPENDENCY --showCoverage [default 30]      |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-maxCov             |Normalized coverage values > \-\-maxCov are shown             |[double]        |
|                       |                                                              |                |
|                       |as \-\-maxCov. DEPENDENCY --showCoverage [default 5]          |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-debug              |Dump session and quit                                         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-h, \-\-help          |Show help message                                             |                |
+-----------------------+--------------------------------------------------------------+----------------+

Description
-----------
| The ``SNV`` module aims at comparing the SNV in terms of number, position and frequency for a samples set.
| For each sample the module loads the GIP files with the filtered SNVs data (singleVariants.df.gz files) and generates multiple plots overlaying in different colors the SNVs sets. 




Output
------





Example
-------



