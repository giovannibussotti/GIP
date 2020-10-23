#########
karyotype
#########

Purpose
-------

The karyotype modules aims at comparing the chromosome sequencing coverage distributions of multiple samples. This module is useful when trying to detect chromosome ploidy differences in different isolates.

Algorithm
---------

For each sample the module loads the GIP files with the nucleotide sequencing coverage normalised by the median genomic coverage (.covPerNt.gz files). The mean of the coverage values is then computed for contiguous genomic intervals (binning). The bin scores are then converted to *somy scores* which are then used for producing plots and statistics.

#use --pcMapqFiles to remove the bases where the percent of mapping reads with good MAPQ is < 50%


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
|\-\-window         |Split chromosomes in intervals of --window base pairs             |[int]           |
|                   |                                                                  |                |
|                   |[default 2500]                                                    |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-chrs           |Chromosomes to use. If "NA" use the same chromsomes as GIP        |[char ...]      |
|                   |                                                                  |                |
|                   |[default NA]                                                      |                |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-ylim           |Min and max plot y-axis limits [default (0, 10)]                  |[int] [int]     |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-makeQqplots    |Computes Q-Q plots for all chromosomes in all samples combinations|                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-selectQuantiles|Remove bins with coverage <10 and >90 quantiles                   |[yes|no]        |
|                   |                                                                  |                |
|                   |Requires --window >1 [default no]                                 |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-disomicChr     |Normalize by this chromosome [default NA]                         |[char]          |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-customColors   |Tab-separated file where the first 2 columns are:                 |[char]          |
|                   |                                                                  |                |
|                   |  * SAMPLE: samples names                                         |                |
|                   |  * COLOR:  associated colors                                     |                |
|                   |                                                                  |                |
|                   |[default NA]                                                      |                |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-geom           |Select boxplot or violin [default boxplot]                        |[boxplot|violin]|
+-------------------+------------------------------------------------------------------+----------------+
|\-\-pooled         |Pool all samples together (i.e. one box per chromosome            |                |
|                   |                                                                  |                |
|                   |representing the coverage values of all samples)                  |                |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-debug          |Dump session and quit                                             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-h, \-\-help      |Show help message                                                 |                |
+-------------------+------------------------------------------------------------------+----------------+


The ``--disomicChr`` option is useful to recenter somy scores on a user defined disomic chromosome.
Under the assumption that most of the genome is disomic, the somy score is simply calculated multiplying by two the bin scores.
However, depending on the biological system under investigation, many/most chromosomes may show aneuploidy. 
In this case the somy score distribution of a given disomic chromosome may not be exactly centered on 2, and different samples may give slightly shifted readout for the same disomic chromosome.
To address this problem the user can perform a first run of ``giptools karyotype``  to identify a chromosome whose median coverage is as close as possible to a value of 2, and that it is stable across the samples set. In a second run of ``giptools karyotype`` the user can then specify the chromsome name with the ``--disomicChr`` option.
By doing that the somy scores will be calculated by dividing the bins coverage scores by the median score of the chromosome deemed to be disomic and then multiplied by 2.

#As an option (--selectQuantiles) out of each distribution of bins it considers just those between the 10th and the 90th quantile, so to strip out outliers bins. This is because outliers, hotpost loci can also artificially return into statistically significant comparisons, while in fact the genomic coverage distributions are very similar (except for these outliers)


#To compare coverage distribution shapes the script generates qqplots (explained here: http://www.r-bloggers.com/exploratory-data-analysis-quantile-quantile-plots-for-new-yorks-ozone-pollution-data/), comparing the quantiles of the two downsampled distributions. If the dots (quantiles) are on the diagonal, then they distributions have the same shape.

#WARNING: since the coverage vectors of a chromosome in two conditions have the same length, in fact the qqplots is equivalent to sort the two vectors and plot one against the other. So the coordinates x,y of the points are the coverage score in the two conditions (you have as many quantiles as nucleotides, not a fixed quantile interval, like the 25%)


Output
------



#To test whether the chromosome coverage varies in different conditions in theory one could use a wilcoxon (default), ks or aov test. 


#the script runs Wilcoxon, Kolmogorov-Smirnov and AOV tests on the downsampled data for each comparison and returns one table for each test.
#it also return the difference between the median downsampled somy score for each chromosome and for each comparison





Example
-------

