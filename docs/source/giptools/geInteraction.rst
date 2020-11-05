#############
geInteraction
#############

purpose
-------


example commandline

example figure

Other options


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
