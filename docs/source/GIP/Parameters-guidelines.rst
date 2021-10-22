#####################
Parameters guidelines
#####################

GIP offers several non-architecture-related parameters that can be adjusted and already described in the GIP steps section.
This page provides some tips on how to set these parameters for the analysis of your genomic sequencing experiment.
These are recommendations based on our work experience. The best parametrization that suits your purpose highly depends on your data and your organism of interest. It is therefore recommendable to explore and test different parametrizations.


* **MAPQ**: Use high values (e.g. 30-50) for a stringent parametrization. Set to 0 to consider all genomic positions to call genomic variants. GIP relies on the BWA mapper which, in case of multiple equally scoring mapping positions, randomly assigns the read to one of the possible locations.    
* **delDup**: It is recommendable to set "true" to remove duplicated reads which represent technical artifacts
* **CGcorrect**: It is recommendable to set "true" if interested in stressing copy number differences between genes or genomic bins of the same sample, especially if the library preparation included PCR amplification steps. This is because the measured differences could be explained by different CG content (rather than a genuine biological signal) causing DNA amplification biases. 
* **chrs**: For unfinished genome assemblies it is recommended to define the list of chromosomes identifiers of interest. While sequencing reads will be mapped to the entire genome, downstream analyses will be limited to the genomic regions of interest (while discarding scaffolds or contigs with potentially sub-optimal annotations). This will improve both computing performance and visualization. If the genome assembly is of good quality just set this parameter to "all". 
* **chrPlotYlim**: These are the y-axis minimum and maximum values for the chromosome ploidy plot. The default, "0 8" is large and inclusive. For better, more zoomed visualizations this value can be tuned and adjusted depending on the target species. For instance, it is rare to see more than 6 chromosome copies in Leishmania.  
* **binPlotYlim**: The y-axis minimum and maximum limits of the genomic bins sequencing coverage plots. The default "0 3" provides a compact visualization where bin amplifications with normalized coverage vlue >3 are shown as 3.  
* **binOverviewSize**: Graphical parameter controlling the heights and the widths of the genomic bin coverage visualizations (default "400 1000").
* **customCoverageLimits**: This parameter can be used to enforce additional custom sequencing coverage thresholds. Significant CNV genomic bins and genes will be retained just if presenting a normalized coverage above the first number, or below the second (default "1.5 0.5"). Assuming diploidy, 1.5 indicates that one of the two chromosomes amplified a copy of the gene (bin). Similarly, 0.5 indicates that one of the two chromosomes lost a copy. This parameter also determines the thresholds above which genes/bins are colored in orange (amplification) and below which are colored in blue (depletion). It is often convenient to use this threshold to limit the number of predictions and focus just on the most relevant ones.
* **binSize**: This parameter governs the resolution of genomic bin analyses. The smaller its value, the higher the number of genomic bins. GIP evaluates sequencing coverage at single nucleotide level (i.e. each nucleotide of each bin). The binSize value can be set smaller than the read size without causing measurement problems. It is recommended to adjust this value depending on the reference genome size. The default, 300,  was successfully applied for the analysis of Leishmania, Candida and Plasmodium genomes. For the analysis of the human genome we used a binSize of 50000 nucleotides.
* **covPerBinSigOPT**: For the identification of significant genomic CNV segments (i.e. "collapsed bins") it is recommended to set the "--padjust BY" to execute the Benjamini & Yekutieli (BY) multiple testing correction (as in the default). BY accounts for the variable dependence that arises from adjacent genomic bins that are part of the same CNV event (i.e. physically connected), thus their coverage estimates are not independent.
* **covPerGeSigOPT**: For the identification of significant CNV genes it is possible to control false discovery rate using the Benjamini & Hochberg (BH) correction, which is valid assuming that genes are separate and their sequencing coverage estimates can be measured independently.   
* **covPerGeRepeatRange**: This parameter defines the maximum distance (in nucleotides) from each gene CNVs in which repeats are labelled. The default is 1000 nucleotides.

* **freebayesOPT**: These are the options directly passed to Freebayes to call SNVs......TO CONTINUE   = "--read-indel-limit 1 --read-mismatch-limit 3 --read-snp-limit 3 --min-alternate-fraction 0.05 --min-base-quality 5 --min-alternate-count 2 --pooled-continuous"

    filterFreebayesOPT = "--minFreq 0.1 --maxFreq 1 --minAO 2 --minAOhomopolymer 20 --contextSpan 5 --homopolymerFreq 0.4 --minMQMR 20 --minMQM 20 --MADrange 4 --randomSNVtoShow 50000"

    filterDellyOPT="--rmLowQual --chrEndFilter 100 --minMAPQ 50 --topHqPercentBnd 150 --topHqPercentIns 150 --topHqPercentDel 150 --topHqPercentDup 150 --topHqPercentInv 150"

    binSizeCircos = 25000

    bigWigOPT = "--normalizeUsing RPKM --ignoreDuplicates --binSize 10 --smoothLength 30" 

