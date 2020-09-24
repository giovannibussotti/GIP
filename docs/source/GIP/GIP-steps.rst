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

* db/
  BWA database 
* geneFunction.tsv
  File provided with the `geneFunction` parameter . The default is a list of gene 
* genome.chrSize

* genome.dict

* genome.fa
* genome.fa.fai
* genome.gaps.gz
* repeatMasker
* repeats
* snpEff

`geneFunction`


