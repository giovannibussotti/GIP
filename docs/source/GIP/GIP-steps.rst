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

prepare genome and annotation
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
| The **genome.fa** file is a copy of the input ``--genome`` file where repetitive positions have been lowercased.
| The **repeats/** directory stores the coordinates of the repetitive elements in a .gff formatted file.
| By default repetitive elements are detected with Red.
| Red features include:

* the ability to detect both transposable elements and simple repeats
* fast execution
* the accuracy when dealing with genomes with unusual nucleotide composition (e.g. Leishmania)

| If the repetitive elements in the species of interest are known the user detect the repeats using RepeatMasker instead.
| To enble RepeatMasker search the user must provide the repeat sequences as multi-FASTA file using the parameter ``--repeatLibrary``



Read mapping
------------

Genomic reads are mapped using BWA-mem
say here about:
The process name
GATK steps
markdups
possible BITFLAG or MAPQ options
the output files and their position in gipOut/samples/

