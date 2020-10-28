###################
Introduction to GIP
###################

Distribution
------------
GIP is a tool for scientific investigation, compatible with Linux and OS X systems and distributed as a self-contained package.

GIP consists in 3 files:

* gip: The Nextflow pipeline
* gip.config: The configuration file
* giptools: The Singularity image  

All the (many) required software dependencies already are already embedded in the Singularity container.

Thanks to the Nextflow implementation GIP can be executed on a local machine (default), on cluster resource manager or the cloud.


Why using GIP?
--------------
GIP is a fully automated pipeline requiring minimal configuration.

While GIP can be used for batch computation of large WGS data sets, the minimum required input is just (i) a paired-end whole genome sequencing data set (.fastq files) and (ii) a reference genome assembly (FASTA file).

GIP is flexible by design, which means that it does not include any built-in hardcoded parametrization (e.g. number/names of chromosomes, centromer position..) limiting its use to higher eukaryotes (e.g. human, mouse) or model organisms in general. 

GIP is particularily adapted to the genome analysis of non-model organisms such as *Leishmaina*. 

GIP can be used to explore the genome instability of biological systems exploiting genome instability for adaptation (Leishmania, Candida, Cancer) through frequent DNA dosage variations (e.g. chromosome aneuploidy, gene CNVs) or SNVs.


GIP analyses overview
---------------------

* Prepare the genome 

* Map the reads

* Evaluate chromosome copy number

* Evaluate gene and bin copy number

* Identify and visualize copy number variation wrt the reference genome

* Identify and measure gene clusters

* Identify SNVs and structural variants

* Generate a report file providing summary statistics, tables in exel, figures




