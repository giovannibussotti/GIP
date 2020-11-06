#############
geInteraction
#############

Purpose
-------
The ``geInteraction`` module aims at detecting CNV genes across multiple samples and identifying gene interactions using a correlation-based network approach.


Algorithm
---------

The module loads the GIP files with the gene sequencing coverage values (.covPerGe.gz files) of all samples, then selects CNV genes. These are defined as the genes with a normalized coverage variation within the sample set greather than --minDelta. Next it builds a network and evaluates clusters based on the correlation computed between all CNV gene pairs. 
  


Options
-------

+-----------------------+--------------------------------------------------------------+----------------+
|Option                 |Description                                                   |Argument        |
+=======================+==============================================================+================+
|\-\-samplesList        |File with a column named \"sample\" listing samples names.    |[char]          |
|                       |                                                              |                |
|                       |Additional TSV columns will be used to annotate the output    |                |
|                       |                                                              |                |
|                       |figures. \"field\"_COLOR columns are used to map colors       |                |
|                       |                                                              |                |
|                       |to the additional fields [**required**]                       |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-gipOut             |GIP output directory [**required**]                           |[char]          |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-outName            |Output name [default NA]                                      |[char]          |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-chrs               |Chromosomes to use. If "NA" it uses the same chromsomes as GIP|[char ...]      |
|                       |                                                              |                |
|                       |[default NA]                                                  |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-minMAPQ            |Remove genes with MAPQ < --MAPQ [default 0]                   |[int]           |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-minDelta           |Min normalized coverage delta between samples [default 1]     |[int]           |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-minMaxCov          |Use only genes with normalized coverage >Value1 or <Value2    |[num num]       |
|                       |                                                              |                |
|                       |in at least one sample.                                       |                |
|                       |                                                              |                |
|                       |If \"NA\" no filter is applied [default NA]                   |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-heatmapType        |Gene normalized coverage value transformation used            |[scaled | log10 |
|                       |                                                              |                |
|                       |for the CNV vs samples heatmap.  [default scaled]             |saturated |     |
|                       |                                                              |                |
|                       |                                                              |flatten]        |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-covSaturation      |Gene normalized coverage saturation value. DEPENDENCY         |[int]           |
|                       |                                                              |                |
|                       |\-\-heatmapType \"saturated\" or \"flatten\" [default 3]      |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-quantileSaturation |Provide two numbers. Saturate the colors of the gene CNV      |[double double] |
|                       |                                                              |                |
|                       |vs samples heatmap for quantiles < num1 or > num2             |                |
|                       |                                                              |                |
|                       |DEPENDENCY \-\-heatmapType \"scaled\" or \"log10\""           |                |
|                       |                                                              |                |
|                       |[default 0 1]                                                 |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-doNotClusterSamples|Do not cluster heatmap columns.                               |                |
|                       |                                                              |                |
|                       |Show the samples in the same order as in \-\-samplesList      |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-clusteringMethod   |Heatmaps clustering method [default complete]                 |[ward.D2|ward   |
|                       |                                                              |                |
|                       |                                                              |single|complete | 
|                       |                                                              |                |
|                       |                                                              |average|mcquitty|
|                       |                                                              |                |
|                       |                                                              |median|centroid]|
+-----------------------+--------------------------------------------------------------+----------------+
|\-\cutree_cnv          |Based on the hierarchical clustering,                         |[int]           |
|                       |                                                              |                |
|                       |divide the genes in this number of clusters [default 1]       |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-cutree_samp        |Based on the hierarchical clustering, divide the samples      |[int]           |
|                       |                                                              |                |
|                       |in this number of clusters [default 1]                        |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-show_geneNames     |Show gene names in the heatmaps                               |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-show_sampNames     |Show sample names in the heatmaps                             |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-kmeansClusters     |NETWORK. Use this number of k-means clusters for              |[int]           |
|                       |                                                              |                |
|                       |network clustering. If \"NA\" use mclust [default NA]         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-MCLinflation       |NETWORK. Use this inflation MCL value for network clustering. |[int]           |
|                       |                                                              |                |
|                       |Higher inflation values result in increased                   |                |
|                       |                                                              |                |
|                       | cluster granularity. If \"NA\" use mclust  [default NA]      |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-MCLexpansion       |NETWORK. MCL expansion value.                                 |[int]           |
|                       |                                                              |                |
|                       |DEPENDENCY \-\-MCLinflation not \"NA\" [default 2]            |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-clMaxSDdist        |NETWORK. Gene CNVs with distance from the cluster             |[double]        | 
|                       |                                                              |                |
|                       |centroid > \-\-clMaxSDdist standard deviations from the       |                |
|                       |                                                              |                |
|                       |mean distance are removed from the cluster. High values make  |                |
|                       |                                                              |                |
|                       |this filter unffective. [default Inf]                         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-clMinSize"         |NETWORK. Min number of members in a cluster [default 2]       |[int]           |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-edgesMeanCorFilter |NETWORK. Remove edges representing CNV correlation scores     |                |
|                       |                                                              |                |
|                       |lower than the mean absolute CNV correlation                  |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-\-edgesPvalueFilter  |NETWORK. Remove edges with adjusted pvalue                    |[double]        |
|                       |                                                              |                |
|                       |below this threshold  [default 0.1]                           |                |
+-----------------------+--------------------------------------------------------------+----------------+  
|\-\-debug              |Dump session and quit                                         |                |
+-----------------------+--------------------------------------------------------------+----------------+
|\-h, \-\-help          |Show help message                                             |                |
+-----------------------+--------------------------------------------------------------+----------------+



heatmapType
If \"scaled\" values are first centered                                                                           
subtracting the mean gene normalized coverage across samples, then scaled dividing by the standard deviation. If \"log10\"  
values are log10 transformed. If \"saturated\" values are     
saturated at \-\-covSaturation. If \"flatten\" values are     
first subracted by the min gene normalized coverage across    
samples, then saturated at \-\-covSaturation




Output
------





Example
-------


#given a folder with multiple covPerBin.gz or covPerGe.gz or chrStartEndScore.gz this script:
#1) selects the bins showing high delta coverage (> --minDelta) (and MAPQ > --minMAPQ) (for covPerBin and covPerGe)
#2) when possible it merges together adjacent bins (with cov > --minDelta) averaging the coverage scores, generating a "CNV" dataset (for covPerBin or chrStartEndScore). CNVs can be filtered by --minCNVLength
#3) generates several heatmaps: 
  #1 Scaled
  #2 log10 
  #3 for each CNV, the values are subtracted by the minimum coverage and then saturated. The latter is useful to focus on coverage variation. This is valuable because it shows you the coverage folds variation much better in situations where a peak (or gene) is highly amplified in all samples (say normalized coverage of 10) and it is hard to appreciate the variation of just one unit (e.g. 10, 11, 9, 10) because the color is saturated 
  #4 saturated scores and using just a four colors palette
  #5 sort columns (samples) by in a specific order defined in sampleSelection. exclude the other samples. (Optional)
  #6 correlation scores (all CNVs vs all CNVs) 
#4) a lollipop plot sorted like the all CNVs vs all CNVs correlation heatmap 
#5) PCA analysis on the CNVs 
#6) hist of entropy and SD of both the selected CNVs and the entire unfiltered set (coverage saturated) 
#7) hierachical clustering on the samples eucledian distance estimated on the peaks   

#The second part of the script is about NETWORKS
 #-given the all vs all CNV correlation matrix (cmr)
 #-take the absolute value of the correlation to consider equally negative and positive correlations
 #-compute mclust clusters 
 #-remove small clusters and the element from the cluster that are far away from the centroid. To do that, for each cluster it measures the centroid (multi dimentional vector) and measure the mean euclidian distance and the standard deviation. Members with distance > clMaxSDdist standard deviations from the mean are removed
 #-write in a folder the filtered clusters
 #-make a network plot (see https://rstudio-pubs-static.s3.amazonaws.com/337696_c6b008e0766e46bebf1401bea67f7b10.html)

#The third part of the script regard tries to turn the igraph network into an interactive network with D3
#example: http://kateto.net/network-visualization
 #The inputs are the standard edges and a nodes data frames, but with a few little twists. 
 #The node IDs in the edges data frame must be integers, and they also have to start from 0. An easy was to get there is to sort the IDs, then transform the character IDs to a factor variable, then transform that to integers (and make sure it starts from zero by subtracting 1).
 #WARNING!!! http://kateto.net/network-visualization is wrong because it converts the source and the target node IDs to integer separatelly. The correct way to do this is implemented in this script. Briefly, 1) sort the edge data frame by IDs in "source"  2) append "source" and "target" together, and assign integer IDs 3) sort the nodes in the nodes dataframe following the same order defined by the node IDS integers

#Rscript  binCoverage2cnvs.R --DIR ../../pipeOut/brazilDeletion/lsdOut/ --minMAPQ 50 --minDelta 1 --outName bin2peakDelta --inFormat covPerBin --filePattern .covPerBin.gz --geBedFile /Volumes/BioIT/Giovanni/datasets/projects/p2p5/Linf.ge.bed 
