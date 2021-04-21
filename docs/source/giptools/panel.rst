#####
panel
#####

Options
-------

+--------------------+------------------------------------------------------------------+---------------+
|Option              |Description                                                       |Argument       |
+====================+==================================================================+===============+
|\-\-samplesList     |File with a column named "sampleId" listing samples names,        |[char]         |
|                    |                                                                  |               |
|                    |plus two TSV columns specifying respectively sample sets          |               |               
|                    |                                                                  |               |
|                    |and colors [required]                                             |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-gipOut          |GIP output directory [default gipOut]                             |[char]         |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-outName         |Output name [default gipOut/sampleComparison/panel]               |[char]         |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-panel           |File with three columns named "gene_id", "set" and "color"        |[char]         |
|                    |                                                                  |               |
|                    |listing respectively the gene ID, the panel group, and the color  |               |
|                    |                                                                  |               |
|                    |[required]                                                        |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-addBackground   |Number of random non-panel genes to be added and shown in         |[int|          |
|                    |                                                                  |               |
|                    |scatterplots. If "all" all non-panel genes are added.             |"all"|"none"]  |
|                    |                                                                  |               |
|                    |If "none" no genes are added [default none]                       |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-covNorm         |Plot gene coverage normalized by                                  |["chromosome"| |
|                    |                                                                  |               |
|                    |[chromosome|genome] median coverage[default genome]               |"genome"]      |
+--------------------+------------------------------------------------------------------+---------------+  
|\-\-covPlotDim      |Gene coverage plot height and width values [default 10 20]        |[double double]| 
+--------------------+------------------------------------------------------------------+---------------+
|\-\-varPlotDim      |Gene variants plot height and width values [default 10 20]        |[double double]| 
+--------------------+------------------------------------------------------------------+---------------+
|\-\-contrast        |Compare samples from these two sets.                              |[char char]    | 
|                    |                                                                  |               |
|                    |If "NA" all samples are considered from the same set [defaul NA]  |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-\-debug           |Dump session and quit                                             |               |
+--------------------+------------------------------------------------------------------+---------------+
|\-h, \-\-help       |Show help message                                                 |               |
+--------------------+------------------------------------------------------------------+---------------+


Description
-----------

| The ``panel`` module aims at extracting and summarizing the normalized coverage and synonimous (S) and non-synonimous (N) SNV information for one or more gene panels in one or more sample sets. The gene panels are specified with the ``--panel`` parameter. The sample sets are specified with the ``--sampleList`` parameter. Optionally the use can compare N, S and coverage statistics in two sample sets of interest by specifying the set names with the ``--contrast`` parameter.   


Example
-------

| From the GIP worked example folder execute

| ``giptools panel --sampleList sampleList.tsv --panel panel.tsv``

| This will generate the panel output files in the **gipOut/sampleComparison** folder. 


