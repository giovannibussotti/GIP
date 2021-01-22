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
|\-\-highLowCovThresh|Provide two numbers. Bins with normalized coverage                |[double]       |
|                    |                                                                  |               |
|                    |values > num1 or < num2 will be labeled [default 1.5 0.5]         |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-ylim            |Plot visualization threshold. Bin or gene normalized coverage     |[double]       |
|                    |                                                                  |               |
|                    |values > --ylim are shown as --ylim [default 5]                   |               |  
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

| The ``overview`` module aims at comparing the samples in terms of chromosomes, genomic bins and genes sequencing coverage. Unlike other modules like ``binCNV`` or ``geCNV`` where normalization accounts for chromosome copy number, data in ``overview`` is normalized by median genome coverage only. ``overview`` is suited to display CNV variation in multiple samples with respect to the reference genome. ``overview`` does not perform sequencing coverage ratios between samples. The normalized scores do not represent somy scores. The normalization procedure is such that the coverage of most genomic bins will be centered on 1.   


Example
-------

| From the GIP worked example folder execute

| ``giptools overview``

| This will generate the overview output files in the **gipOut/sampleComparison** folder. 

| The output consists in four files: 

* The .chrCov.pdf file represents the chromosome coverage:

.. figure:: ../_static/karyotype1.png
      :width: 100 %

* The .binCov.pdf file represents the genomic bin coverage:

.. figure:: ../_static/karyotype1.png
      :width: 100 %

* The .geCov.pdf file represents the gene coverage:

.. figure:: ../_static/karyotype1.png
      :width: 100 %

* The geCov.xlsx file is an excel table reporting the normalized gene coverage with the associated function (if available) 


