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
| This GIP feature is useful when dealing with unfinished genome assemblies, containing large amounts of unplaced contigs with very poor annotation available
| For this purpose, the user can set the parameter ``--chromosomes``, listing the identifiers of the chromosomes of interest.
| By default this parameter reports the 36 *Leishmania* chromosome identifiers.


Measure nucleotide coverage
---------------------------

| Mapped reads are used to measure the sequencing coverage of each nucleotide in the *covPerNt* process.
| Tools used at this step include Samtools view and Bedtools genomecov (options "-d -split").
| The reads mapping with the bitflag (see `SAM format specifications <https://samtools.github.io/hts-specs/SAMv1.pdf>`_) value given by the ``--BITFLAG`` parameter (default 1028) are excluded.
| To account for differences in sequencing library size and enable comparisons between samples, the nucleotide sequencing coverage is normalized by the median genomic coverage.
| The files generated at this step are placed in the **gipOut/samples/sampleId** folder, and include:

+----------------------------------------+-------------------------------------+
| sampleId.covPerNt.gz                   | nucelotide coverage                 |
+----------------------------------------+-------------------------------------+
| sampleId.covPerNt.medianGenomeCoverage | median genome coverage              |
+----------------------------------------+-------------------------------------+
| sampleId.pcMapqPerNt.gz                | % of high MAPQ reads per nucleotide |
+----------------------------------------+-------------------------------------+

| The syntax of the **sampleId.covPerNt.gz** file is: chromosome        position        normalized sequencing coverage
| The **sampleId.pcMapqPerNt.gz** file reports the percent of reads with MAPQ greater or equal to the ``--MAPQ`` value.
| The file syntax is: chromosome	position	%reads
| These files are used to evaluate the chromosomes somy score distritributions and generate these additional results: 

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









