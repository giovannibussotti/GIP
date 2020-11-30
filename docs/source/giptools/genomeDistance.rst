##############
genomeDistance
##############


Options
-------

+-------------------+---------------------------------------------------+----------------+
|Option             |Description                                        |Argument        |
+===================+===================================================+================+
|\-\-samples        |Sample names. It determines the plotting order     |[char ...]      |
|                   |                                                   |                |
|                   |If "NA" all samples are used [default NA]          |                |
+-------------------+---------------------------------------------------+----------------+
|\-\-gipOut         |GIP output director [default gipOut]               |[char]          |
+-------------------+---------------------------------------------------+----------------+
|\-\-outName        |Output name                                        |[char]          |
|                   |                                                   |                |
|                   |[default gipOut/sampleComparison/genomeDistance]   |                |
+-------------------+---------------------------------------------------+----------------+
|\-\-contextGenomes |Two columns TSV file listing additional genomes.   |[char]          |                
|                   |                                                   |                |
|                   |Syntax: name<Tab>filePath [default NA]             |                |
+-------------------+---------------------------------------------------+----------------+  
|\-\-debug          |Dump session and quit                              |                |
+-------------------+---------------------------------------------------+----------------+
|\-h, \-\-help      |Show help message                                  |                |
+-------------------+---------------------------------------------------+----------------+

Description
-----------
| The ``genomeDistance`` module is meant to compare the genomic distance of the samples of interest. 
| For each sample the module loads the genome file embedding the predicted filtered SNVs (pseudoReference.fa.gz files), then measures the ANI (average nucleotide identity) scores for all sample pairs. The resulting distance is then used to run a principal component analysis (PCA) and produce a dendrogram plot demonstrating the measured distances.  



Output
------




Example
-------


