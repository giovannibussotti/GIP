#########
phylogeny
#########

Options
-------

+-------------------+------------------------------------------------------------------+----------------+
|Option             |Description                                                       |Argument        |
+===================+==================================================================+================+
|\-\-samples        |Sample names. If \"NA\" all samples are used                      |[char ...]      |
+-------------------+------------------------------------------------------------------+----------------+ 
|\-\-gipOut         |GIP output directory [default gipOut]                             |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-outName        |Output name [default gipOut/sampleComparison/phylogeny]           |[char]          |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-VRFcutoff      |Provide 2 numbers. For each sample, heterozygous SNVs where       |[double double] |
|                   |                                                                  |                |
|                   |NUM2 <= VRF < NUM are labeled with the IUPAC ambiguous notation.  |                | 
|                   |                                                                  |                |
|                   |SNVs with VRF < NUM2 are removed [default 0.70 0.10]              |                |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-regionsGtf     |Select just the SNVs present in regions specified in this GTF file|[char]          |
|                   |                                                                  |                |
|                   |(e.g. genes). If \"NA\" nor region filter is applied [default NA] |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-iqtreeOpts     |IQ-tree options [default "-alrt 1000 -bb 1000 -nt 1 -quiet"]      |[char]          |
+-------------------+------------------------------------------------------------------+----------------+  
|\-\-debug          |Dump session and quit                                             |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-h, \-\-help      |Show help message                                                 |                |
+-------------------+------------------------------------------------------------------+----------------+

Description
-----------
| The ``phylogeny`` module computes the phylogenetic tree by maximum likelihood. In a first step the module extracts the union of all filtered SNVs (i.e. the ensemble of all SNVs found in all samples) and returns a multi-FASTA file (\-\-outName file).
| To reduce the potentially negative impact of neutrally selected variants on the tree inference the user can limit the analyses to the subset of SNVs present in specific genomic regions (e.g. genetic regions, \-\-regionsGtf option).
| To account for heterozygous SNV, the user can specify a frequency range in which the variants are labeled with the ambiguous IUPAC notation (i.e. representing both the reference and the alternate allele).
| The module executes IQ-tree2 to compute the tree (and optionally to choose the model and perform boostrap support). The user can pass all the desired iqtree options via the \-\-iqtreeOpts parameter. Eventually distance matrices based on predicted trees are produced.
| Caveat. The predicted phylogeny is based on SNVs only and excludes indels. Additionally, for the tree prediction the algorithm considers just the SNV positions but ignores the conserved positions. As a consequence, despite the predicted tree branches order is reliable, the braches length may be inaccurate. 


Example
-------
| From the GIP worked example folder execute

| ``giptools phylogeny``

| This will generate the phylogeny output files in the **gipOut/sampleComparison** folder. 
| The **phylogeny** file is the multi-FASTA file including the SNVs union (see above). The ">reference" reports the ensemble of the genome reference alleles at the SNV positions.  
| please refer to the IQ-tree documentation for more details about the output it produces.




