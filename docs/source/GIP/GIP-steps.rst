#########
GIP steps
#########

GIP accepts the following mandatory input parameters:

+----------------+-----------------------------------+    
| \-\-genome     | multi-FASTA genome reference file |
+----------------+-----------------------------------+
| \-\-annotation | gene coordinates file (GTF format)|
+----------------+-----------------------------------+
| \-\-index      | list of sequencing data files     |
+----------------+-----------------------------------+
| \-c            | nextflow configuration file       |
+----------------+-----------------------------------+

| The index file must comply with the following syntax rules:

1. tsv format (i.e. <Tab> separated), 
2. first header row with the labels: sampleId   read1    read2
3. all the following rows must indicate the sample identifier, and the position of the first and second pair-end sequencing data files in fastq.gz format

| Example:   
| sampleId        read1    read2  
| sample1 /home/user/data/s1.r1.fastq.gz  /home/user/data/s1.r2.fastq.gz  
| sample2 /home/user/data/s1.r1.fastq.gz  /home/user/data/s1.r2.fastq.gz  

| GIP results are cached in the **work** directory inside subfolders named with the hexadecimal numbers identifying the executed processes.       
| GIP results are organized and accessible from the **gipOut** directory, hosting symbolic links to the data in **work**.
| The ``--resultDir`` parameter can be used to set an alternative name for the result directory.


In the following we provide a description of GIP steps operated by the Nextflow processes.

Prepare genome and annotation
-----------------------------
GIP prepares the genome assembly and annotation files in the first two processes: *processGeneFunction* and *prepareGenome*
The output is stored in the **gipOut/genome** directory and includes:

+-----------------+----------------------------+
| db/             | BWA database               |
+-----------------+----------------------------+
| geneFunction.tsv| gene functions list        |
+-----------------+----------------------------+
| genome.chrSize  | chromosome size (bp)       |
+-----------------+----------------------------+
| genome.dict     | genome sequence dictionary |
+-----------------+----------------------------+
| genome.fa       | genome sequence            |
+-----------------+----------------------------+
| genome.fa.fai   | genome sequence index      |
+-----------------+----------------------------+
| genome.gaps.gz  | genome gaps                |
+-----------------+----------------------------+
| repeats/        | repeats coordinates        |
+-----------------+----------------------------+
| snpEff/         | snpEff database            |
+-----------------+----------------------------+

| In many non-model organisms the gene function may be not know. If gene function information is available (or at least available for some genes) this should be provided to GIP with the ``--geneFunction`` parameter. This file must:

* include all the genes listed in the ``--annotation`` file
* be a <Tab> separated list with the gene identifier in the first column, and the function in the second:   

| e.g.
| LinJ.01.0010	Protein of unknown function (DUF2946)
| LinJ.01.0020	Endonuclease/Exonuclease/phosphatase family
| LinJ.01.0030	Kinesin-13
| LinJ.01.0040	hypothetical protein - conserved

| If ``--geneFunction`` is not specified by default the **geneFunction.tsv** file reports the gene list with not available (NA) functions.
| The **genome.dict** and the **genome.fa.fai** are generated respectivelly with *picard CreateSequenceDictionary* and *samtools faidx* and are required by downstream analysis tools (e.g. GATK or IGV). 
| The **genome.fa** file is a copy of the input ``--genome`` file where the chromosome identifiers containing white spaces are changed into underscores, and the repetitive positions have been lowercased.
| The **repeats/** directory stores the coordinates of the repetitive elements in a .gff formatted file.
| By default repetitive elements are detected with Red.
| Red features include:

* the ability to detect both transposable elements and simple repeats
* fast execution
* the accuracy when dealing with genomes with unusual nucleotide composition (e.g. Leishmania)

| If the repetitive elements in the species of interest are known the user detect the repeats using RepeatMasker instead.
| To enble RepeatMasker search the user must provide the repeat sequences as multi-FASTA file using the parameter ``--repeatLibrary``



Map reads and collect alignment statistics
------------------------------------------

| Genomic reads are aligned to the reference genome in the *map* process. 
| Tools used at this step include:

* Samotools modules including fixmate, sort and index to reformat the aliment files
* GATK RealignerTargetCreator and IndelRealigner to homogenize the indels
* Picad MarkDuplicates with option VALIDATION_STRINGENCY=LENIENT to label potential optical or PCR read duplicates

| The files generated at this step are placed in the **gipOut/samples/sampleId** folder, and include:

+-----------------------------+-----------------------------------------------+
| sampleId.bam                | alignment file                                |
+-----------------------------+-----------------------------------------------+
| sampleId.bam.bai            | alignment index file                          |
+-----------------------------+-----------------------------------------------+
| sampleId.MarkDup.table      | picard MarkDuplicates tabular output          |
+-----------------------------+-----------------------------------------------+
| sampleId.MarkDup.histData   | picard MarkDuplicates histogram data output   |
+-----------------------------+-----------------------------------------------+
| sampleId.MarkDup.hist.png   | plot of picard MarkDuplicates histogram data  |
+-----------------------------+-----------------------------------------------+
| chrCoverageMedians_sampleId | chromosome coverage table                     |
+-----------------------------+-----------------------------------------------+

| Alignment statistics are produced in the *mappingStats* process using the picard CollectAlignmentSummaryMetrics and CollectInsertSizeMetrics (option "MINIMUM_PCT=0") tools, and include the following files:

+---------------------------------+--------------------------------------------------------+
| sampleId.alignmentMetrics.table | picard CollectInsertSizeMetrics tabular output         |
+---------------------------------+--------------------------------------------------------+
| sampleId.insertSize.histData    | picard CollectInsertSizeMetrics histogram data output  |
+---------------------------------+--------------------------------------------------------+
| sampleId.insertSize.hist.png    | plot of picard CollectInsertSizeMetrics histogram data |
+---------------------------------+--------------------------------------------------------+
| sampleId.insertSize.table       | picard CollectInsertSizeMetrics tabular output         |
+---------------------------------+--------------------------------------------------------+



Evaluate chromosome coverage
----------------------------

| Alignment files are used to evaluate the chromosome sequencing coverage in the *covPerChr* process.
| At this step the  **chrCoveraMedians_sampleId** table is generated in the **gipOut/samples/sampleId** folder.
| This table is used by GIP for downstream normalization steps, and reports the following fields:

+--------------------+---------------------------------------------+
| CHR	             | chromosome identifier                       |
+--------------------+---------------------------------------------+
| MEDIANCOV	     | median chromosome sequencing coverage       |
+--------------------+---------------------------------------------+
| MEDIANCOVminus2MAD | MEDIANCOV plus 2 median absolute deviation  |	
+--------------------+---------------------------------------------+
| MEDIANCOVplus2MAD  | MEDIANCOV minus 2 median absolute deviation |
+--------------------+---------------------------------------------+

| While reads are mapped in the previous step against the entire genome, the user may want to instruct GIP to consider for this step and all the downstream analyses just a sub-set of chromosomes. 
| This GIP feature is useful when dealing with unfinished genome assemblies, containing large amounts of unplaced contigs with very poor annotation available.
| For this purpose, the user can set the parameter ``--chromosomes``, listing the identifiers of the chromosomes of interest.
| By default this parameter reports the 36 *Leishmania* chromosome identifiers.


Measure nucleotide coverage
---------------------------

| Mapped reads are used to measure the sequencing coverage of each nucleotide in the *covPerNt* process.
| Tools used at this step include Samtools view and Bedtools genomecov (options "-d -split").
| The reads mapping with the bitflag (see `SAM format specifications <https://samtools.github.io/hts-specs/SAMv1.pdf>`_) value given by the ``--BITFLAG`` parameter (default 1028) are excluded.
| This parameter applies with the same function also to downstream processes, namely: *covPerBin*, *covPerGe* and *delly*.
| To account for differences in sequencing library size and enable comparisons between samples, the nucleotide sequencing coverage is normalized by the median genomic coverage.
| The files generated at this step are placed in the **gipOut/samples/sampleId** folder, and include:

+----------------------------------------+-------------------------------------+
| sampleId.covPerNt.gz                   | nucelotide coverage                 |
+----------------------------------------+-------------------------------------+
| sampleId.covPerNt.medianGenomeCoverage | median genome coverage              |
+----------------------------------------+-------------------------------------+
| sampleId.pcMapqPerNt.gz                | % of high MAPQ reads per nucleotide |
+----------------------------------------+-------------------------------------+

| The syntax of the **sampleId.covPerNt.gz** file is: chromosome<Tab>position<Tab>normalized sequencing coverage
| The **sampleId.pcMapqPerNt.gz** file reports the percent of reads with MAPQ greater or equal to the ``--MAPQ`` value.
| The file syntax is: chromosome<Tab>position<Tab>%reads
| These files are used to evaluate the chromosomes somy score distritributions and generate additional results providing a karyotype overview: 

+----------------------------------+----------------------------------+
| sampleId.covPerNt.allMedians.tsv |  chromosomes median somy scores  | 
+----------------------------------+----------------------------------+
| sampleId.covPerNt.boxplot.png    |  somy scores boxplot             |
+----------------------------------+----------------------------------+
| sampleId.covPerNt.ridges.png     |  somy scores ridge plot          |
+----------------------------------+----------------------------------+

| To reduce noise, CPU and memory requirements GIP downsamples the **sampleId.covPerNt.gz** nucleotide coverage scores by binnig the genome into 2500 nucleotide long windows. 
| Then for each window the somy score is computed measuring the mean nucleotide coverage scores and multiplying by 2.
| The chromosome median somy score reflects the chromosome copy number under the assuption that most nucleotides in the genome are present in two copies (e.g. disomic chromosomes).


Measure genomic bin sequencing coverage
---------------------------------------

| Mapped reads are used to measure the sequencing coverage of genomic bins in the *covPerBin* process.
| The ``--binSize`` parameter (default 300) controls the bin size (i.e. the number of nucleotides for each bin).
| The sequencing coverage of each bin normalized by 

| GIP At this step:

1. Computes the sequencing depth of each nucleotide without normalizing 
2. Divides the genome in contiguous genomic bins whose size is determined by the ``--binSize`` parameter (default 300bp)
3. Computes mean and median sequencing coveage scores for each bin, and normalize them by median chromosome sequencing coverage
4. Estimates the mean MAPQ score for each bin  

| Please note that it is possible to obtain genomic bins with 0 mean or median coverage, but MAPQ greather than 0. This is the case in genomic depletions where very few reads map to the bin with a certain MAPQ score greather than 0. 
| Bin coverage scores are then corrected for GC content to limit potential sequencing biases during DNA amplification. Given the distribution of bin mean coverage scores and GC-content, GIP fits a loess regression using using a 5 folds cross validation to explore the loess *span* parameter (which relates with the fraction of points used to fit the local regressions, and influence the model smoothness).
| Then GIP corrects the original bin coverage by subtracting the values on the loess model, and adding back the difference between the median coverage of all bin before and after subtraction (i.e. recentering the median bin coverage to 1). Genomic bins that after correction have negative coverage are reported with a 0 value.


| The resulting bin are evaluated for significant copy number variation (CNV) with respect to the reference genome. Often, the CNV span regions larger than the bin size. In order to match the size of the CNV region (at a bin size resolution), GIP collapses adjacent significant CNV bins of the same type (i.e. adjacent bins composing a depletion, or adjacent bins composing an amplification), then averages their coverage score. We refer to these sets of collapsed bins as **segments**.

| For the statistical test GIP derives the single nucleotide coverage distribution after binning (SNCDab) where the coverage of each nucleotide is approximated with the mean bin coverage.  
| For the central limit theorem (CLT):

* Regardless the shape of SNCDab, the sampling distribution of the sample means (SDSM) is gaussian
* The mean (mu) and the standard error (se) of SNCDab correspond to the mean (mu) and the standard deviation (sd) of SDSM with sample size equal n
              
| For each bin the null-hypothesis is that it is possible to observe its sequencing coverage just by chance under a normal (i.e. non-CNV) condition due to coverage fluctuations intruduced by the sequencing technology. The competing hypothesis is that the oberved coverage is the readout of a genuine CNV region.
| Based on the CLT, GIP computes the P-value of each bin by measuring how many se away each bin score is from the SNCDab mu.

| The ``--covPerBinSigPeaksOPT`` parameter accepts a string of 3 parameters, and can be used to customize the detection of bin and segments of interest.

* *--minLen*  - minimum segment length (bp) [int]
* *--pThresh* - adjusted p-value threshold [num]
* *--padjust* - multiple-testing correction method [num]

| The ``--covPerBinSigPeaksOPT`` default is ``"--minLen 0 --pThresh 0.001 --padjust BY"``. The available methods for multiple testing corrections are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Please refer to documentation of the `p.adjust <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust>`_ R function for more details.

| The ``--customCoverageLimits`` parameter can be used to enforce an additional custom coverage cut-offs on the statistically significant bins and segments (and genes, see below). This parameter accepts two numbers: N1, N2 (default 1.5 0.5). Significant CNV bins and segments are selected to have a coverage > N1 (for amplifications) or < N2 (for depletions). 

| The *covPerBin* process returns the following files in the **gipOut/samples/sampleId** folder


+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.gz                          | genomic bin coverage                           |
+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.plot.all.png                | bin coverage genome overview                   |
+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.plot.byChr.pdf              | bin coverage chromosome overview (slides)      |
+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.plot.faceting.png           | bin coverage chromosome overview (multi-panel) |      
+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.plot.tsv.gz                 | bin coverage plots data                        |
+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.significant.bins.tsv.gz     | significant bins                               |
+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.significant.segments.tsv.gz | significant segments                           |
+------------------------------------------------+------------------------------------------------+
| sampleId.covPerBin.significant.stats           | statistical test info                          |
+------------------------------------------------+------------------------------------------------+

| In all three plots, the bins with mean MAPQ lower than ``--MAPQ`` are shown in gray. The statistically significant bins corresponding to amplifications and depletions are shown respectivelly in orange and blu. The y-axis minimum and maximum limits can be set with the parameter ``--binPlotYlim`` (default ``"0 3"``). Depending on the genome size the overview plots may result too small and unreadable. The parameter ``--binOverviewSize`` accepts two integers controlling respectivelly the plots heights and the widths (default ``"400 1000"``). The values specified with the ``--customCoverageLimits`` parameter will be highligthed with red dashed lines.



Measure gene sequencing coverage
--------------------------------

| Mapped reads are used to measure the mean sequencing coverage of annotated genes in the *covPerGe* process.  
| To estimate the mean coverage the N bases are not considered. GIP normalizes the coverage scores by the chromosome median coverage. correct for potential GC-content biases at gene level GIP utilizes the same approach described for genomic bins (see above).To detect statistically significant CNV genes GIP fits a gaussian mixture distribution with 2 components. One distribution accounting for the vast majority of observations fitting the coverage of non-CNV genes (central distribution), and another distribution fitting the CNV genes (outliers distribution). The cental distributions represents the-null hypothesis under which a given coverage value is merely caused by artefact fluctuations in sequencing depth, rather than a genuine, biologically meaningful gene amplification or depletion. To test CNV significance GIP uses the mean and the standard deviation of the central distribution and assigns a z-score and a p-value to all genes. Significant genes with a mean MAPQ score lower than ``--MAPQ`` are discarded. In the same way as for genomic bins, the parameter ``--customCoverageLimits``can be used to enforce custom coverage threshold on significant genes. The parameter ``--covPerGeSigPeaksOPT`` accepts  a string of 3 parameters and can be used to control the statical test.

* *--pThresh* - adjusted p-value threshold [num] 
* *--padjust* - method for multiple testing correction [num]
* *--minLen*  - minimum gene size (bp) [int]

| The default is ``covPerGeSigPeaksOPT="--pThresh 0.001 --padjust BH --minLen 0"``. As for genomic bins, the available methods for multiple testing corrections are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Please refer to documentation of the `p.adjust <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust>`_ R function for more details.

| The *covPerGe* process returns the following files in the **gipOut/samples/sampleId** folder


+--------------------------------------+-----------------------------+
| sampleId.covPerGe.gz                 | gene sequencing coverage    |
+--------------------------------------+-----------------------------+
| sampleId.covPerGe.significant.tsv    | significant gene CNVs       |
+--------------------------------------+-----------------------------+
| sampleId.covPerGe.significant.stats  | statistical test info       |
+--------------------------------------+-----------------------------+
| sampleId.covPerGeKaryoplot/          | folder with CNV genes plots |
+--------------------------------------+-----------------------------+

The **sampleId.covPerGeKaryoplot/** folder includes plot generated with the `karyoploteR <https://www.bioconductor.org/packages/release/bioc/html/karyoploteR.html>`_ package. Only chromosomes hosting significant gene CNVs are shown. Amplified genes are shown in orange, whereas depleted genes are shown in blue. If any, the repetitive elements located in proximity of gene CNVs are marked in the bottom part of the plots. The ``--repeatRange`` parameter can be used to set the maximum distance (in nucleotides) from each gene CNVs in which repeats are labelled.


Detect, annotate and filter single nucleotide variants
------------------------------------------------------

| The single nucleotide variants (SNVs) are detected in the *freebayes* process using the `freebayes <https://arxiv.org/abs/1207.3907>`_ program, and their effects are predicted in the *snpEff* process running `snpEff <https://pcingola.github.io/SnpEff/se_introduction/>_` with option "-ud 0".
| Reads with MAPQ score < than ``--MAPQ`` are not used for detecti on. The user can specify freebayes options through the ``--freebayesOPT`` parameter. Its default is:

.. code-block:: bash

 --freebayesOPT="--read-indel-limit 1 --read-mismatch-limit 3 --read-snp-limit 3 \
 --hwe-priors-off --binomial-obs-priors-off --allele-balance-priors-off \
 --min-alternate-fraction 0.05 --min-base-quality 5 --min-alternate-count 2 --pooled-continuous"


Please refer to the `freebayes manual <https://github.com/ekg/freebayes>`_ for more details.
| GIP returns the following outputs in the **gipOut/samples/sampleId/** folder:

+--------------------------------------+---------------------------------------------+
| sampleId.vcf.gz                      | SNVs (gzip compressed vcf file)             |
+--------------------------------------+---------------------------------------------+
| sampleId.vcf.gz.tbi                  | tabix vcf index                             |
+--------------------------------------+---------------------------------------------+
| snpEff_summary_sampleId.genes.txt.gz | SNVs per gene, snpEff summary table         |
+--------------------------------------+---------------------------------------------+
| snpEff_summary_sampleId.html         | snpEff summary (html)                       |
+--------------------------------------+---------------------------------------------+

| SNV mapping to predicted repetitive elements, or mapping inside low-complexity regions (homopolymer) are at higher risk to be sequencing artefacts. 
| To diminish the number of false positives and short-list high quality SNVs GIP operates additional filters. 
| GIP discards all SNVs mapping inside repetitive elements, removes the variant positions with multiple alternate alleles, evaluates the nucleotide composition complexity of the genomic context of each SNV (i.e. the neighbour bases) and allows the user to apply different, more stringent, filterering criteria for variants detected inside homopolymers.  
| For this purpose the ``--filterFreebayesOPT`` parameter can be used to set the following variables:

* *--minFreq*          - Min. variant read frequency (VRF) [num]
* *--maxFreq*          - Max. VRF [num]
* *--minAO*            - Min. number of reads supporting the alternate allele [int] 
* *--minMQMR*          - Min. mean mapping quality of observed reference alleles [num]
* *--minMQM*           - Min. mean mapping quality of observed alternate alleles [num]
* *--MADrange*         - Discard SNVs whose sequencing depth is > or < *MADrange* MADs from the chromosome median coverage [num]
* *--minAOhomopolymer* - Min. number of reads supporting the alternate allele mapping inside an homopolymer [int]
* *--contextSpan*      - Size on each side of SNV genomic context (bp) [int]
* *--homopolymerFreq*  - Base frequency cut-off to consider a genomic context a homopolymer [num]


| The parameter default is:

.. code-block:: bash

   filterFreebayesOPT="--minFreq 0.1 --maxFreq 1 --minAO 2 --minAOhomopolymer 20 \ 
   --contextSpan 5 --homopolymerFreq 0.4 --minMQMR 20 --minMQM 20 --MADrange 4"

| The results relative to the filtered SNVs are stored in the **gipOut/samples/sampleId/sampleId_freebayesFiltered/** folder including:


+-------------------------------------------------+------------------------------------------------------------+
| singleVariants.df.gz                            | SNVs (table)                                               |
+-------------------------------------------------+------------------------------------------------------------+
| singleVariants.vcf.gz                           | SNVs (gzip compressed vcf file)                            |
+-------------------------------------------------+------------------------------------------------------------+
| singleVariants.vcf.gz.tbi                       | tabix vcf index                                            |
+-------------------------------------------------+------------------------------------------------------------+
| single_allDensities.png                         | VRF density plot                                           |
+-------------------------------------------------+------------------------------------------------------------+
| single_allHists.png                             | VRF histogram plot                                         |
+-------------------------------------------------+------------------------------------------------------------+
| single_allHistsSqrt.png                         | VRF histogram plot (sqrt scale)                            |
+-------------------------------------------------+------------------------------------------------------------+
| single_combinedDotPlotAndDistribution.pdf       | position/VRF plot with marginal distribution               |
+-------------------------------------------------+------------------------------------------------------------+
| single_depthVsVRF.png                           | VRF/depth plot                                             |
+-------------------------------------------------+------------------------------------------------------------+
| single_depthVsVRFletters.png                    | VRF/depth plot                                             |
|                                                 |                                                            |
|                                                 | SNV chromosomes are mapped to different colors and letters |
+-------------------------------------------------+------------------------------------------------------------+
| single_onePlotPerChr.pdf                        | position/VRF and density plots per chromosome              |
+-------------------------------------------------+------------------------------------------------------------+
| single_onePlotPerChr_colouredByVariantType.pdf  | position/VRF colored by SNV type                           |
+-------------------------------------------------+------------------------------------------------------------+
| single_totVarPerChr.png                         | num. SNVs per chromsome kb                                 |
+-------------------------------------------------+------------------------------------------------------------+
| single_variantType.png                          | occurrence of different SNV types                          |
+-------------------------------------------------+------------------------------------------------------------+
| single_variantTypeCombined.png                  | occurrence of different SNV types                          |
|                                                 |                                                            |
|                                                 | equivalent variants combined                               |
+-------------------------------------------------+------------------------------------------------------------+
| single_VRFvsAO.png                              | VRF/alternate allele read support                          |
+-------------------------------------------------+------------------------------------------------------------+
| single_VRFvsAOletters.png                       | VRF/alternate allele read support                          |
|                                                 |                                                            |
|                                                 | SNV chromosomes are mapped to different colors and letters |
+-------------------------------------------------+------------------------------------------------------------+
| single_VRFvsPosFaceting.png                     | position/VRF plot with different chromosomes               |
|                                                 |                                                            |
|                                                 | in different panels                                        |
+-------------------------------------------------+------------------------------------------------------------+
| snpEff_summary_sampleId.genes.txt.gz            | SNVs per gene, snpEff summary table                        |
+-------------------------------------------------+------------------------------------------------------------+
| snpEff_summary_sampleId.html                    | snpEff summary (html)                                      |
+-------------------------------------------------+------------------------------------------------------------+
| dNdS.stats                                      | dNdS analysis statistics                                   |
+-------------------------------------------------+------------------------------------------------------------+
| dNdStable.tsv.gz                                | dNdS analysis per gene                                     |
+-------------------------------------------------+------------------------------------------------------------+
| pseudoReference.fa.gz                           | genome sequence incorporating alternate alleles            |
+-------------------------------------------------+------------------------------------------------------------+
| context/                                        | folder containing the nucleotide frequency logo plots of   |
|                                                 |                                                            |
|                                                 | the genomic contexts of different SNV types                |
+-------------------------------------------------+------------------------------------------------------------+


| For the dNdS analysis the snpEff effects counting as synonimous substitutions are:

* SYNONYMOUS_CODING
* SYNONYMOUS_STOP

| The snpEff effects counting as non-synonimous substitutions are:

* NON_SYNONYMOUS_CODING
* NON_SYNONYMOUS_START 
* START_LOST
* STOP_GAINED
* STOP_LOST



Detect structural variants
--------------------------
 
| The genomic structural variants (SVs) are detected in the *delly* process using the `delly <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436805/>_` program. The SVs are predicted based on pair-end mapping orientation and split-read information, and include unbalanced reaffangements (i.e. CNV deletions or amplifications), as well as balanced rearrangements (inversions and translocations). delly is used to predict the four SV types using just the reads passing both the ``--MAPQ`` and ``--BITFLAG`` filters. The output is the .vcf gzip compressed file  **gipOut/samples/sampleId/sampleId.delly.vcf.gz**.
|GIP allows to apply custom quality filters and select a short-list of SV predictions using the ``--filterDellyOPT`` parameter, and setting the following variables:

* *--minDV*          - min. num. of read pairs supporting the variant [int] 
* *--minPercentDVDR* - min. percent of read pairs supporting the variant [num] 
* *--PRECISE*        - delly "PRECISE" attribute [yes|no] 
* *--maxBanSeq*      - discard SVs where the percent of overlapping repeats or gap sequences is > --maxBanSeq [num]
* *--chrEndFilter*   - num. of bases spanning from the chromosome ends inwards. SVs overlapping such telomeric or sub-telomeric regions are discarded [int]

| The parameter default is:

.. code-block:: bash

   filterDellyOPT="--minDV 2 --minPercentDVDR 5 --PRECISE no \
   --maxBanSeq 90 --chrEndFilter 100"


| The results relative to the filtered SVs are stored in the **gipOut/samples/sampleId/sampleId_dellyFiltered/** folder including:


+-------------------------------------+----------------------------------+
| sampleId.delly.DEL.filter           | deletions table                  |
+-------------------------------------+----------------------------------+
| sampleId.delly.DEL.filter.circosBed | deletions coordinates            |
+-------------------------------------+----------------------------------+
| sampleId.delly.DUP.filter           | tandem duplications table        |
+-------------------------------------+----------------------------------+
| sampleId.delly.DUP.filter.circosBed | tandem duplications coordinates  |
+-------------------------------------+----------------------------------+
| sampleId.delly.INV.filter           | inversions table                 |
+-------------------------------------+----------------------------------+
| sampleId.delly.INV.filter.circosBed | inversions coordinates           |
+-------------------------------------+----------------------------------+
| sampleId.delly.TRA.filter           | translocations table             |
+-------------------------------------+----------------------------------+
| sampleId.delly.TRA.filter.circosBed | translocations coordinates       |
+-------------------------------------+----------------------------------+
| sampleId_circosData/                | data for circos plot             |
+-------------------------------------+----------------------------------+
| sampleId.SV.circos.png              | circos plot                      |
+-------------------------------------+----------------------------------+

All coordinates files are in bed format, except for **sampleId.delly.TRA.filter.circosBed**, where the six fields correspond to the coordinates (chromosome<Tab>start<Tab>end) of the two translocation break points.














