##########
binDensity
########## 

Options
-------

+-----------------------+--------------------------------------------------------------+----------------+
|Option                 |Description                                                   |Argument        |
+=======================+==============================================================+================+
|\-\-samplesList        |File listing in one column the samples to use.                |[char ...]      |
|                       |                                                              |                |
|                       |If "NA" all samples are used. [default NA]                    |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-gipOut             |GIP output directory [default gipOut]                         |[char]          |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-outName            |Output name [default gipOut/sampleComparison/binDensity]      |[char]          |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-chrs               |Chromosomes to use. If "NA" it uses the same chromsomes as GIP|[char ...]      |
|                       |                                                              |                |
|                       |[default NA]                                                  |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-minMAPQ            |Remove bins with MAPQ < --MAPQ [default 0]                    |[int]           |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-pseudocount        |Normalized mean coverage pseudocount value preventing minus   |[double]        |
|                       |                                                              |                |
|                       |infinite values in log10 transformation [default 0.1]         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-highLowCovThresh   |Provide two numbers. Bins with normalized coverage            |[double double] |
|                       |                                                              |                |
|                       |values > num1 or < num2 will be labeled. [default 1.5 0.5]    |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-bandwidth          |Smoothing bandwidth value.                                    |[double double] |
|                       |                                                              |                |
|                       |Provide two numbers to enforce different bandwidths           |                |
|                       |                                                              |                |
|                       |on the x and y axes respectively [default 10000 0.01]         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-nbin               |Number of equally spaced grid points for                      |[int int]       |
|                       |                                                              |                |
|                       |the density estimation. Provide two numbers to use different  |                |
|                       |                                                              |                |
|                       |numbers for the x and y axes respectively [default 1000]      |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-showSubset         |Show a random subset of genomic bins with normalized coverage |[int]           |
|                       |                                                              |                |
|                       |values above or below --highLowCovThresh. This random subset  |                |
|                       |                                                              |                |
|                       |does not affect the density estimation [default 50000]        |                |
+-----------------------+--------------------------------------------------------------+----------------+  
|\-\-debug              |Dump session and quit                                         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-h, \-\-help          |Show help message                                             |                |
+-----------------------+--------------------------------------------------------------+----------------+

Description
-----------
| The ``binDensity`` module aims at visualizing the normalized bin sequencing coverage of multiple samples using a smoothed color density representation.
| The module loads for all samples the GIP files with the bin sequencing coverage values (.covPerBin.gz files) and generates a smoothed color density scatterplot showing the genomic position (x-axis) and the log10 normalized coverage values (y-axis).


Output
------





Example
-------

