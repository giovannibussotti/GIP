#########
GIP steps
#########

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
| The **genome.fa** file is a copy of the input ``--genome`` file where the repetitive positions have been lowercased.
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
* Picad MarkDuplicates with option VALIDATION_STRINGENCY=LENIENT to detect and remove optical or PCR read duplicates

| The ``delDup`` parameter can be set to false to just label instead of removing reads duplicates. The files generated at this step are placed in the **gipOut/samples/sampleId** folder, and include:

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

| The genome sequencing coverage density is available in bigWig format and computed in the *bigWigGenomeCov* process.
| The bigWig file format is compatible with genome browsers such as `IGV <http://software.broadinstitute.org/software/igv/>`_. A description of the bigWig format is available from `here <https://genome.ucsc.edu/goldenPath/help/bigWig.html>`_. GIP generates the bigWig output file **gipOut/samples/sampleId/sampleId.bw** by applying the bamCoverage module of `deepTools2 <https://academic.oup.com/nar/article/44/W1/W160/2499308>`_. The coverage values are generated ignoring duplicated reads and applying an RPKM normalization on separate chromosomes (bamCoverage options "--normalizeUsingRPKM --ignoreDuplicates"). GIP approach makes the coverage density estimates comparable between chromosomes that may have different copy numbers. The user can provide additional options to bamCoverage with the ``--bigWigOPT`` parameter. The default is ``bigWigOPT="--binSize 10 --smoothLength 30"``, where the two options control the sizes of the bigWig bins (bp) and the size of the window to average the number of reads. Please refer to the bamCoverage `documentation <http://gensoft.pasteur.fr/docs/deepTools/2.4.2/content/tools/bamCoverage.html>`_ for more details.



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
| For this purpose, the user can set the parameter ``--chrs``, listing the identifiers of the chromosomes of interest. Please note that the chromosome identifiers must be separated by white spaces and provided as a text string enclosed in backslash escaped quotation marks. For instance if the user what to limit the analysis to chromosomes chr1, chr2 and chr3 then parameter should be ``chrs="\"1 2 3\""``.
| Otherwise, by default this parameter is set to "all", indicating that all chromosomes present in the input fasta genome will be considered.


Measure genomic bin sequencing coverage
---------------------------------------

| Mapped reads are used to measure the sequencing coverage of genomic bins in the *covPerBin* process.
| The ``--binSize`` parameter (default 300) controls the bin size (i.e. the number of nucleotides for each bin).
| The sequencing coverage of each bin normalized by 

| GIP At this step:

1. Computes the sequencing depth of each nucleotide without normalizing 
2. Divides the genome in contiguous genomic bins whose size is determined by the ``--binSize`` parameter (default 300bp)
3. Computes mean sequencing coveage scores for each bin
4. Normalizes the mean bin coverage by median chromosome sequencing coverage
5. Applies a GC-content correction on the normalized mean bin coverage (optional)
6. Estimates the mean MAPQ score for each bin  

| Please note that it is possible to obtain genomic bins with 0 mean coverage, but MAPQ greather than 0. This is the case in genomic depletions where very few reads map to the bin with a certain MAPQ score greather than 0. 
| The GC-content correction is enabled setting the parameter ``CGcorrect = true`` and is meant to limit potential sequencing biases during DNA amplification. Given the distribution of the normalized bin mean coverage scores and their GC-content, GIP fits a loess regression using using a 5 folds cross validation to explore the loess *span* parameter (which relates with the fraction of points used to fit the local regressions, and influence the model smoothness).
| Then GIP corrects the original bin coverage by subtracting the values on the loess model, and adding back the difference between the median coverage of all bin before and after subtraction (i.e. recentering the median of the bin coverage scores to 1). Genomic bins that after correction have negative coverage are reported with a 0 value.


| The resulting bin are evaluated for significant copy number variation (CNV) with respect to the reference genome. Often, the CNV span regions larger than the bin size. In order to match the size of the CNV region (at a bin size resolution), GIP collapses adjacent significant CNV bins of the same type (i.e. adjacent bins composing a depletion, or adjacent bins composing an amplification), then averages their coverage score. We refer to these sets of collapsed bins as **segments**.
    
| For each bin the null-hypothesis is that it is possible to observe its sequencing coverage just by chance under a normal (i.e. non-CNV) condition due to coverage fluctuations intruduced by the sequencing technology. The competing hypothesis is that the oberved coverage is the readout of a genuine CNV region.
| The cental limit teorem (CLT) states that the distribution of the samples means approximates to a normal distribution. As a consequence, the distribution of the genomic bin mean coverage approximates to a gaussian as well. GIP computes the P-value of each bin by measuring the number of standard deviations from the mean. 

| The ``--covPerBinSigOPT`` parameter accepts a string of 3 parameters, and can be used to customize the detection of bin and segments of interest.

* *--minLen*  - minimum segment length (bp) [int]
* *--pThresh* - adjusted p-value threshold [num]
* *--padjust* - multiple-testing correction method [num]

| The ``--covPerBinSigOPT`` default is ``"--minLen 0 --pThresh 0.001 --padjust BY"``. The available methods for multiple testing corrections are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Please refer to documentation of the `p.adjust <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust>`_ R function for more details.

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
| sampleId.bed                                   | mapped reads in bed format                     |
+------------------------------------------------+------------------------------------------------+

| In all three plots, the bins with mean MAPQ lower than ``--MAPQ`` are shown in gray. The statistically significant bins corresponding to amplifications and depletions are shown respectivelly in orange and blu. The y-axis minimum and maximum limits can be set with the parameter ``--binPlotYlim`` (default ``"0 3"``). Depending on the genome size the overview plots may result too small and unreadable. The parameter ``--binOverviewSize`` accepts two integers controlling respectivelly the plots heights and the widths (default ``"400 1000"``). The values specified with the ``--customCoverageLimits`` parameter will be highligthed with red dashed lines. The **sampleId.bed** file is an intermediate file used by GIP from the quantification of genomic intervals. It is not automatically removed by GIP because it allows the user to re-execute the pipeline with the ``-resume`` option. However, if the user is not planning on re-executing GIP he/she can simply delete this file from the **work/** directory to save disk space.     


| Genomic bin sequencing coverage values are also used to compute the chromosome somy score distritributions and evaluate the chromosome copy number. Bins whith mean MAPQ score lower than the ``--MAPQ`` value are not considered.
| To account for differences in sequencing library size and enable comparisons between samples, the mean bin sequencing coverage is normalized by the median of all genomic bins. 
| Then for each window the somy score is computed measuring the mean nucleotide coverage scores and multiplying by 2.
| The chromosome median somy score reflects the chromosome copy number under the assuption that most nucleotides in the genome are present in two copies (e.g. disomic chromosomes).
| The files produced at this step provide an overview of the sample karyotype and include:

+----------------------------------------+--------------------------------+
| sampleId.karyotype.medianCoverage      | median coverage of all bins    |
+----------------------------------------+--------------------------------+
| sampleId.karyotype.allMedians.tsv      | chromosomes median somy scores | 
+----------------------------------------+--------------------------------+
| sampleId.karyotype.boxplot.png         | somy scores boxplot            |
+----------------------------------------+--------------------------------+
| sampleId.karyotype.ridges.png          | somy scores ridge plot         |
+----------------------------------------+--------------------------------+
 


Measure gene sequencing coverage
--------------------------------

| Mapped reads are used to measure the mean sequencing coverage of annotated genes in the *covPerGe* process.  
| GIP normalizes the coverage scores by the chromosome median coverage. To correct for potential GC-content biases at gene level GIP utilizes the same approach described for genomic bins (option enabled by ``CGcorrect = true``, see above).To detect statistically significant CNV genes GIP fits a gaussian mixture distribution with 2 components. One distribution accounting for the vast majority of observations fitting the coverage of non-CNV genes (central distribution), and another distribution fitting the CNV genes (outliers distribution). The cental distributions represents the-null hypothesis under which a given coverage value is merely caused by artefact fluctuations in sequencing depth, rather than a genuine, biologically meaningful gene amplification or depletion. To test CNV significance GIP uses the mean and the standard deviation of the central distribution and assigns a z-score and a p-value to all genes. Significant genes with a mean MAPQ score lower than ``--MAPQ`` are discarded. In the same way as for genomic bins, the parameter ``--customCoverageLimits`` can be used to enforce custom coverage threshold on significant genes. The parameter ``--covPerGeSigOPT`` accepts  a string of 3 parameters and can be used to control the statical test.

* *--pThresh* - adjusted p-value threshold [num] 
* *--padjust* - method for multiple testing correction [num]
* *--minLen*  - minimum gene size (bp) [int]

| The default is ``covPerGeSigOPT="--pThresh 0.001 --padjust BH --minLen 0"``. As for genomic bins, the available methods for multiple testing corrections are: "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none". Please refer to documentation of the `p.adjust <https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/p.adjust>`_ R function for more details.

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

The **sampleId.covPerGeKaryoplot/** folder includes plot generated with the `karyoploteR <https://www.bioconductor.org/packages/release/bioc/html/karyoploteR.html>`_ package. Only chromosomes hosting significant gene CNVs are shown. Amplified genes are shown in orange, whereas depleted genes are shown in blue. If any, the repetitive elements located in proximity of gene CNVs are marked in the bottom part of the plots. The ``--repeatRange`` parameter can be used to set the maximum distance (in nucleotides) from each gene CNVs in which repeats are labelled. To put the gene CNVs in context of possible larger CNV regions the figure also reports a gray slope indicating the normalized bin coverage scores. In most cases the normalized coverage values of genes and bins are very close. However, for certain genes much shorter than the bin size, the plots may show a discrepancy between bin and gene readouts.  


Detect, annotate and filter single nucleotide variants
------------------------------------------------------

| The single nucleotide variants (SNVs) are detected in the *freebayes* process using the `freebayes <https://arxiv.org/abs/1207.3907>`_ program, and their effects are predicted in the *snpEff* process running `snpEff <https://pcingola.github.io/SnpEff/se_introduction/>_` with option "-ud 0".
| Reads with MAPQ score < than ``--MAPQ`` are not used for detecti on. The user can specify freebayes options through the ``--freebayesOPT`` parameter. Its default is:

.. code-block:: bash

 --freebayesOPT="--read-indel-limit 1 --read-mismatch-limit 3 --read-snp-limit 3 \
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

* synonymous_variant
* stop_retained_variant
* start_retained

| The snpEff effects counting as non-synonimous substitutions are:

* missense_variant
* start_lost
* stop_gained
* stop_lost
* coding_sequence_variant


Detect and filter structural variants
-------------------------------------
 
| The genomic structural variants (SVs) are detected in the *delly* process using the `delly <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3436805/>_` program. The SVs are predicted based on pair-end mapping orientation and split-read information, and include unbalanced reaffangements (i.e. CNV deletions, amplifications and insertions), as well as balanced rearrangements (inversions and break ends translocations). delly is used to predict the five SV types using just the reads passing the ``--MAPQ`` filter. The outputs are the .vcf bgzip compressed file  **gipOut/samples/sampleId/sampleId.delly.vcf.gz** and its tabix index with .tbi extension.
| GIP allows to apply custom quality filters and select a short-list of SV predictions using the ``--filterDellyOPT`` parameter, and setting the following variables:

* *--minMAPQ*      - min median mapping quality of paired-ends supporting the SV [int]
* *--chrEndFilter* - num. of bases spanning from the chromosome ends inwards. SVs overlapping such telomeric or sub-telomeric regions are discarded [int]
* *--rmLowQual*    - Remove delly predictions labelled as LowQual
* *--rmImprecise*  - Keep just delly predictions labelled as PRECISE
* *--topRc[Bnd|Ins|Del|Dup|Inv]*        - Select top SVs based on RC score [int]
* *--topHqCount[Bnd|Ins|Del|Dup|Inv]*   - Select top SVs based on DV+RV score [int]
* *--topHqPercent[Bnd|Ins|Del|Dup|Inv]* - Select top SVs based on (DV+RV/DV+RV+DR+RR)*100 score [int]

| The --topRc* --topHqCount* and --topHqPercent* filters are applied sequentially and consist in selecting the best predicted SVs based on 3 different quality metrics: the RC score, the DV+RV score and the (DV+RV/DV+RV+DR+RR)*100 score.
| Unless specified in ``--filterDellyOPT``, none of these filter is used. To use these filters it is needed to specify the filter type with a suffix indicating the SV of interest: "Bnd" (break ends), "Ins" (insertions), "Del" (deletions), "Dup" (duplications) and "Inv" (inversions).
| for instance ``--topHqCountInv 50`` would select the 50 predicted inversions with the best DV+RV score.
| The vcf description of the RC, DV, RV, DR, RR scores is the following: 

* RC: Raw high-quality read counts or base counts for the SV
* DV: # high-quality variant pairs
* RV: # high-quality variant junction reads
* DR: # high-quality reference pairs
* RR: # high-quality reference junction reads

| The parameter default is:

.. code-block:: bash
   
   filterDellyOPT="--rmLowQual --chrEndFilter 100 --minMAPQ 50 --topHqPercentBnd 150 \ 
   --topHqPercentIns 150 --topHqPercentDel 150 --topHqPercentDup 150 --topHqPercentInv 150"


| The duplication and deletion analysis performed by delly is complementary to the analysis performed considering the sequencing coverage only. Genuine deletions or depletions may not always show a variation in sequencing coverage. Whole genome sequencing data obtained from cell populations is such that a given locus under evolutive pressure can be amplified in a sub-population, and deleted in in another sub-population. Moreover, in biological systems with high DNA plasiticity such as the human pathogen *Leishmania*, a genomic region can undergo multiple, complex genomic rearrangements and shuffling whose presence may be revealed by read pair mapping orientation or split-read information, but not necessarily by sequencing coverage variations.            
| The results relative to the filtered SVs are stored in the **gipOut/samples/sampleId/sampleId_dellyFiltered/** folder including:

+------------------------+----------------------------------+
| output.vcf.gz          | compressed vcf file              |
+------------------------+----------------------------------+
| output.vcf.gz.tbi      | tabix index                      |
+------------------------+----------------------------------+
| DEL.bed                | deletions coordinates            |
+------------------------+----------------------------------+
| DUP.bed                | tandem duplications coordinates  |
+------------------------+----------------------------------+
| INV.bed                | inversions coordinates           |
+------------------------+----------------------------------+
| BND.bed                | break end  coordinates           |
+------------------------+----------------------------------+
| INS.bed                | insertions coordinates           |
+------------------------+----------------------------------+
| sampleId_circosData/   | data for circos plot             |
+------------------------+----------------------------------+
| sampleId.SV.circos.png | circos plot                      |
+------------------------+----------------------------------+

For circos plot representation the chromosomes of interest are binned in into genomic intervals whose size (bp) is regulated by ``--binSizeCircos`` (default 25000). In the the inner part of circos plot the predicted break ends translocations events are shown as black lines. The karyotype color reflects the mean reads MAPQ score calculated for each genomic bin. Black indicates a MAPQ < 2, gray indicates a MAPQ ≥ 2 and < 20 and white indicates a MAPQ ≥ 20. Ticks positions and ticks labels are automatically assigned by GIP depending on genome size. If any, the position of insetions is indicated by red stripes on the karyotype. 

Moving outwards the circos plot shows the tracks relative to predicted duplications (orange), deletions (blue) and inversions (green). The outmost track shows the genomic bin sequencing coverage (light blue bars) normalized by chromosome median coverage and ranging from 0 to 3. To ease visualization, amplifications with normalized coverage greather than 3 are shown with a value of 3.      



Define and quantify gene clusters
--------------------------------- 

Depending on the sequencing technology and the experimental design, annotated genes presenting very high levels of sequence similarity may be difficoult to quantify.The length of the genomic reads and the fragment size influence the read MAPQ scores, thus the unicity of the read alignment.Instead of quantifying individual genes, GIP allows to quantify homologous genes as clusters. Given the set of gene coverage (.covPerGe.gz) files generated for each sample, GIP:
 
1. Selects genes that cannot be directly quantified, i.e. have a mean MAPQ lower than the ``--MAPQ`` value in all samples 
2. Runs `cd-hit-est <http://weizhongli-lab.org/cd-hit/>`_ with option "-g 1" to cluster these genes by sequence similarity 
3. Evaluates the sequencing coverage of the genes belonging to clusters
4. Computes mean sequencing coverage for each gene cluster

The gene clusters analysis is run in the *covPerClstr* process, and the results are stored in the **gipOut/covPerClstr** folder.

+-------------------------+--------------------------------------------------------+
| clstrAnn.tsv            | predicted gene clusters (list format)                  |
+-------------------------+--------------------------------------------------------+
| clstrAnnFormat2.tsv     | predicted gene clusters (table format)                 |
+-------------------------+--------------------------------------------------------+
| sampleId.covPerClstr.gz | mean sequencing gene cluster coverage (gzip compressed)|
+-------------------------+--------------------------------------------------------+
| lowMapq.clstr/          | folder storing the gene cluster sequences              |
+-------------------------+--------------------------------------------------------+

Genes with low mean MAPQ in all samples but not clustering by sequence similarity are kept and part of the output. Normally these genes get a low MAPQ score either because they present internal repetitive sequences, or because their gene or pseudogene homologue is not annotated.

















