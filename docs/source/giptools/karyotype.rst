#########
karyotype
#########

purpose
-------

#First generate genome coverage files (.covPerNt) where the sequencing depth of each nucleotide is normalised by the median genomic coverage 
#Then given a list of gzipped .covPerNt files (--file), this script downsample the covPerNt scores (e.g. by binning), then convert the bin scores to somy scores which will be used for the boxplots and statistics 

#The somy score is defined as the downsampled score (e.g. bin score) divided by the median downsampled score of a chromosome deemed to be disomic and multiplied by 2
#If no disomic chromosome is deemed, the somy score is just the downsampled score multiplied by 2
#To assess what chromosome is disomic you can do a dry run of the script and pick a chromsome that has median cenetered on 2 and does not vary much across the various samples

#To test whether the chromosome coverage varies in different conditions in theory one could use a wilcoxon (default), ks or aov test. 
#In practice, since there are many many observations (nucleotide) the distributions of the coverages will look always statistically significant (even if from the boxplot sometimes you can see they are very very similar). 
#Since the coverage vectors are big the statistical test have an hard time (wilcoxon or ks tests return mostly pvalues of 0 and 1). 
#that is why it is needed to adopt a downsampling approach. This script explores two alternative strategies.

#One is to use the resample function (--fft yes) which transforms the coverage profile to the Fourier Space and resize the signal before performing the inverse transform. It is explained here: https://support.bioconductor.org/p/66313/
#The size of the coverage distribution after shrinking is regulated by --maxShrinkedLength. The script will automatically estimate a shrinking factor to have a final chromosome coverage vector shorter than --maxShrinkedLength.
#It works OK for short chromosomes but it takes ages for big chromosomes. Also the fourier transform can introduce negative values (that do not make sense since it is a coverage), and the reported pvalues are not stable at all when changing --maxShrinkedLength

#The other, which is recommended, simply bins the chromosomes and estimate the median coverage of each bin (--window).
#It works much better than the fourier approach.
#As an option (--selectQuantiles) out of each distribution of bins it considers just those between the 10th and the 90th quantile, so to strip out outliers bins. This is because outliers, hotpost loci can also artificially return into statistically significant comparisons, while in fact the genomic coverage distributions are very similar (except for these outliers)

#the script runs Wilcoxon, Kolmogorov-Smirnov and AOV tests on the downsampled data for each comparison and returns one table for each test.
#it also return the difference between the median downsampled somy score for each chromosome and for each comparison

#To compare coverage distribution shapes the script generates qqplots (explained here: http://www.r-bloggers.com/exploratory-data-analysis-quantile-quantile-plots-for-new-yorks-ozone-pollution-data/), comparing the quantiles of the two downsampled distributions. If the dots (quantiles) are on the diagonal, then they distributions have the same shape.

#use --pcMapqFiles to remove the bases where the percent of mapping reads with good MAPQ is < 50%


#WARNING: since the coverage vectors of a chromosome in two conditions have the same length, in fact the qqplots is equivalent to sort the two vectors and plot one against the other. So the coordinates x,y of the points are the coverage score in the two conditions (you have as many quantiles as nucleotides, not a fixed quantile interval, like the 25%)
#WARNING: Converting covPerNt scores to somy scores first, and then downsample (e.g. binning) would results in slightly differet results (medians not exactly centered on 2)



example commandline

example figure

Other options

