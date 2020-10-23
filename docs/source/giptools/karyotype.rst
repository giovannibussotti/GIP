#########
karyotype
#########

Purpose
-------

The karyotype modules aims at comparing the chromosome sequencing coverage distributions of multiple samples. This module is useful when trying to detect chromosome ploidy differences in different isolates.

Algorithm
---------

For each sample the module loads the GIP files with the nucleotide sequencing coverage normalised by the median genomic coverage (.covPerNt.gz files). The mean of the coverage values is then computed for contiguous genomic intervals (binning). The bin scores are then converted to *somy scores* which are then used for producing plots and statistics.



Options
-------


--minDV          - min. num. of read pairs supporting the variant [int] 
--minPercentDVDR* - min. percent of read pairs supporting the variant [num] 
--PRECISE*        - delly "PRECISE" attribute [yes|no] 
--maxBanSeq*      - discard SVs where the percent of overlapping repeats or gap sequences is > --maxBanSeq [num]
--chrEndFilter*   - num. of bases spanning from the chromosome ends inwards. SVs overlapping such telomeric or sub-telomeric regions are discarded [int]


+------------------+---------------------------------------------------------------+
* --gipOut         | output directory [default NA]
+------------------+---------------------------------------------------------------+
* --samples        | Sample names. It determines the plotting order [default None]
+------------------+---------------------------------------------------------------+
* --outName        | Output name [default NA]
+------------------+---------------------------------------------------------------+
* --ylim           | min and max ylim used in the plot. [default 0 10]
+------------------+---------------------------------------------------------------+
* --chrs           | chromosomes to use. If not defined, it considers the
                   |    chromosomes selected at gip runtime [default NA]
+------------------+---------------------------------------------------------------+
* --makeQqplots    | for all samples combinations, and for each chromosome
                        it computes qQplots [default False]
+------------------+---------------------------------------------------------------+
  --window WINDOW       split each chr in chunks of this size, then take mean
                        coverage out of each window to compare ditributions
                        [default 2500]
+------------------+---------------------------------------------------------------+
  --selectQuantiles SELECTQUANTILES
                        strip out bins <10 and >90 quantiles before comparing
                        coverage distributions (DEPENDENCY: --window) [default
                        no]
+------------------+---------------------------------------------------------------+
  --disomicChr DISOMICCHR
                        normalize by this chromosome [default NA]
+------------------+---------------------------------------------------------------+
  --customColors CUSTOMCOLORS
                        provide a file with header "SAMPLE HEX" as first two
                        columns, specifying the color for each sample [default
                        NA]
+------------------+---------------------------------------------------------------+
  --geom GEOM           select boxplot or violin [default boxplot]
+------------------+---------------------------------------------------------------+
  --pooled              flag. pool together all the samples (i.e. one box per
                        chromosome representing the coverage values of all
                        samples) [default False]
+------------------+---------------------------------------------------------------+
* *-h, --help*       - help message
  --debug               dump session and quit [default False]





Under the assumption that most of the genome is disomic, the somy score is simply calculated multiplying by two the bin scores.
More complex biological system where not just aneuploidy but also partial aneuploidy

The somy score is defined as the downsampled score (e.g. bin score) divided by the median downsampled score of a chromosome deemed to be disomic and multiplied by 2
#If no disomic chromosome is deemed, the somy score is just the downsampled score multiplied by 2
#To assess what chromosome is disomic you can do a dry run of the script and pick a chromsome that has median cenetered on 2 and does not vary much across the various samples



Output
------



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

