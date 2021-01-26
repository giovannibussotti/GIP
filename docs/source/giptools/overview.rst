########
overview
########

Options
-------

+--------------------+------------------------------------------------------------------+---------------+
|Option              |Description                                                       |Argument       |
+====================+==================================================================+===============+
|\-\-samples         |Sample names. It determines the plotting order                    |[char ...]     |
|                    |                                                                  |               |
|                    |If "NA" all samples are used [default NA]                         |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-gipOut          |GIP output directory [default gipOut]                             |[char]         |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-outName         |Output name [default gipOut/sampleComparison/overview]            |[char]         |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-chrs            |Chromosomes to use. If "NA" it uses the same chromsomes as GIP    |[char ...]     |
|                    |                                                                  |               |
|                    |[default NA]                                                      |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-MAPQ            |Label bins with MAPQ < --MAPQ [default 0]                         |[int]          |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-ploidy          |Genome ploidy. After normalization multiply                       |[int]          |
|                    |                                                                  |               |
|                    |coverage values by --ploidy [default 2]                           |               |
+--------------------+------------------------------------------------------------------+---------------+ 
|\-\-highLowCovThresh|Provide two numbers. Bins with normalized coverage                |[double]       |
|                    |                                                                  |               |
|                    |values > num1 or < num2 will be labeled [default 1.5 0.5]         |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-minDelta        |Remove genes with normalized coverage delta < --minDelta          |[double]       |
|                    |                                                                  |               |
|                    |[default 0]                                                       |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-ylim            |Plot visualization threshold. Bin or gene normalized coverage     |[double]       |
|                    |                                                                  |               |
|                    |values > --ylim are shown as --ylim [default 5]                   |               |  
+--------------------+------------------------------------------------------------------+---------------+
|\-\-ylimInt         |Bin scatterplot y-axis values interval. If \"NA\" it is           |[double]       |
|                    |                                                                  |               |
|                    |automatically assigned [default NA]                               |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-maxGe           |Plot visualization threshold. Max number of genes to show         |[int]          |
|                    |                                                                  |               |
|                    |(ordered by delta coverage) [default 10000]                       |               |  
+--------------------+------------------------------------------------------------------+---------------+  
|\-\-binPlotDim      |Bin coverage plot height and width values [default 10 20]         |[double double]| 
+--------------------+------------------------------------------------------------------+---------------+
|\-\-debug           |Dump session and quit                                             |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-h, \-\-help       |Show help message                                                 |               |
+--------------------+------------------------------------------------------------------+---------------+


Description
-----------

| The ``overview`` module aims at comparing the samples in terms of chromosomes, genomic bins and genes sequencing coverage. Unlike other modules like ``binCNV`` or ``geCNV`` where normalization accounts for chromosome copy number, data in ``overview`` is normalized by median genome coverage only. ``overview`` is suited to display CNV variation in multiple samples with respect to the reference genome. ``overview`` does not perform sequencing coverage ratios between samples. The normalized scores do not represent somy scores. The normalization procedure is such that the coverage of most genomic bins will be centered on 1. The ``--ploidy`` parameter can be used to indicate the organism ploidy level, which is 2 by default (diploid). ``overview`` multiplies the normalized coverage scores of bins, genes and chromosomes by the ploidy value.     


Example
-------

| From the GIP worked example folder execute

| ``giptools overview``

| This will generate the overview output files in the **gipOut/sampleComparison** folder. 

| The output consists in four files: 

* The .chrCov.pdf file represents the normalized chromosome coverage
* The .binCov.pdf file represents the normalized genomic bin coverage
* The .geCov.pdf file represents the normalized gene coverage. The first heatmap reports scaled values. The second heatmap shows the actual normalized gene coverage, but values greather than --ylim are reported as --ylim. 
* The geCov.xlsx file is an excel table reporting the normalized gene coverage with the associated function (if available) 

| In the genomic bins plot is possible to center the coverage to 1, limit the y-axis to 2.5 and remove bin coloring by adding the options  ``--ploidy 1 --ylim 2.5 --highLowCovThresh 100 -1`` to the command. This gives the following plot:

.. figure:: ../_static/overview.binCov.png
      :width: 100 %
 
| The options ``--highLowCovThresh 1.25 0.5 --MAPQ 50`` can be used to color the genomic bins with normalized coverage above 1.25 and to label low MAPQ bins:

 .. figure:: ../_static/overview.binCov2.png
      :width: 100 %








