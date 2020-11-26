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
|\-\-gipOut         |GIP output directory. If "NA" the directory "./gipOut" is used    |[char]          |
|                   |                                                                  |                |
|                   |[default NA]                                                      |                |
+-------------------+------------------------------------------------------------------+----------------+
|\-\-outName        |Output name [default NA]                                          |[char]          |
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
| The ``phylogeny`` module to compute the phylogenetic tree by maximum likelihood. In a first step the module extracts the union of all SNVs (i.e. the ensemble of all SNVs found in all samples) and returns a multi-FASTA file (\-\-outName file).
| To reduce the impact of variants following neutral selection on the tree the user can limit the SNVs to a subset og genomic regions (e.g. genetic regions, \-\-regionsGtf option).
| To account for heterozygous SNV, the user can specify a frequency range in which the variants are labeled with the ambiguous IUPAC notation (i.e. representing both the reference and the alternate allele).
| The module executes IQ-tree2 to compute the tree (and optionally to chose the model and perform boostrap support). The user can pass all the desired iqtree options via the \-\-iqtreeOpts parameter. Eventually distance matrices based on predicted trees are produced.

Output
------





Example
-------
