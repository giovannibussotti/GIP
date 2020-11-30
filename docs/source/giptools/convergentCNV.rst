#############
convergentCNV
#############

Options
-------

+-------------------+------------------------------------------------------------------+----------------+
|Option             |Description                                                       |Argument        |
+===================+==================================================================+================+
|\-\-gipOut         |GIP output directory  [default gipOut]                            |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-outName        |Output name [default gipOut/sampleComparison/convergentCNV]       |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-newickTree     |Newick tree file produced with the phylogeny module               |[char]          |
|                   |                                                                  |                |
|                   |[default gipOut/sampleComparison/phylogeny.treefile]              |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-geCNV          |Gene CNV excel file produced with the geInteraction module        |[char]          |
|                   |                                                                  |                |
|                   |[default gipOut/sampleComparison/geInteraction.CNV.xlsx]          |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-ampThresh      |Gene amplification threshold [default 1.5]                        |[double]        |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-notAmpThresh   |Threshold to consider a gene not amplified [default 1]            |[double]        |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-covSaturation  |Gene normalized coverage saturation value [default 2]             |[double]        |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-treeTip        |Color tree tip by this feature.                                   |[char]          |
|                   |                                                                  |                |
|                   |If \"NA\" tree tips are not colored [default NA]                  |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-tipLab         |Show sample names on the tree tips                                |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-hexpand        |Avoid truncating tip labels expanding the plot panel              |[double]        |
|                   |                                                                  |                |
|                   |by a ratio of the x axis range [default 0]                        |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-branchLen      |Branch length ggtree parameter [default none]                     |[branch.length  |
|                   |                                                                  |                |
|                   |                                                                  | | none]        |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-layout         |Layout ggtree parameter [default rectangular]                     |[rectangular|   |
|                   |                                                                  |                |
|                   |                                                                  |slanted |       |
|                   |                                                                  |                |
|                   |                                                                  |circular|       |      
|                   |                                                                  |                |
|                   |                                                                  |fan |           |
|                   |                                                                  |                |
|                   |                                                                  |daylight]       |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-hideHeatmap    |Do not display the heatmap                                        |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-heatCols       |Heatmap gradient colors [default black blue white red3 red]       |[char ...]      |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-heatMin        |Heatmap minimum normalized coverage value [default 0]             |[double]        |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-hideXlabs      |Do not display x-axis text                                        |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-plotDim        |File height and width values [default 6 7]                        |[double ...]    |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-debug          |Dump session and quit                                             |                |
+-------------------+------------------------------------------------------------------+----------------+


Description
-----------

Detect convergent CNV gene amplifications

Output
------





Example
-------
