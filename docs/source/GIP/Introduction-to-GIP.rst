###################
Introduction to GIP
###################

Distribution
------------
GIP is a tool for scientific investigation, compatible with Linux and OS X systems and distributed as a self-contained package.

GIP consists in 3 files:

* gip.nf: The Nextflow pipeline
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

* prepare the genome (repeats, gaps...)

* map the reads (mark duplicates, homogenizing indels)

* Evaluate chromosome copy number (median genome coverage normalization)

* Evaluate gene and bin copy number (chromosome coverage normalization)

* Identify and visualize copy number variation wrt the reference (add some gallery image, e.g. the mbio paper for ge CNV and the covPerBin overview)

* provide bigWig file with a normalization accounting for chromosome copy number to be readily used to visualize and compare ina meaningful way regions of interest (e.g. genes) generating snapshots with genome browsers like IGV (add a gallery link to the MBIO paper IGV snapshots of CNV regions)

* Identify and measure gene clusters

* Identify SNVs and structural variants

* Generate report file providing summary statistics, tables in exel, figures




