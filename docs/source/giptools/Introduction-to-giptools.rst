########################
Introduction to giptools
########################

To run GIP is required to specify in the configuration file the path to the **giptools**. This is a singularity container embedding all GIP scripts and software dependencies. But giptools is more than that. GIP focus is on each individual sample, and all results (e.g. the detection of a gene amplification or a SNV) are measured with respect to the genome reference. giptools can be used to address specific questions regarding the comparisons of samples subsets (e.g. compare the gene copy number between two samples in particular). 

The giptools file can be directly executed and used for the downstream comparison of samples processed by GIP. 
Typing ``giptools`` with no options (or with -h or --help) will list the available comparison modules.
To use any of these modules the syntax is simply ``giptools moduleName``. Each module has its own help documentation accessible with -h or --help. By default all comparison modules generate their output in the **gipOut/sampleComparison** folder.   

gitools modules:

+----------------+--------------------------------------------------------------------------+
| karyotype      | Compare the chromosome sequencing coverage distributions                 |
+----------------+--------------------------------------------------------------------------+
| binCNV         | Compare bin sequencing coverage in 2 samples                             |
+----------------+--------------------------------------------------------------------------+
| geCNV          | Compare gene sequencing coverage in 2 samples                            |
+----------------+--------------------------------------------------------------------------+
| ternary        | Compare gene sequencing coverage in 3 samples                            |
+----------------+--------------------------------------------------------------------------+
| ternaryBin     | Compare bin sequencing coverage in 3 samples                             |
+----------------+--------------------------------------------------------------------------+
| SNV            | Compare SNVs in multiple samples                                         |
+----------------+--------------------------------------------------------------------------+
| binDensity     | Density plot of bin sequencing coverage of many samples                  |
+----------------+--------------------------------------------------------------------------+
| geInteraction  | Detect CNV genes in many samples and produce correlation-based networks  |
+----------------+--------------------------------------------------------------------------+
| genomeDistance | Compare samples genomic distance                                         |
+----------------+--------------------------------------------------------------------------+

All modules require as an input the ``--gipOut`` parameter specifying the GIP output folder.
A description of each of the giptools modules and the available options is provided in the following pages of this documentation.

Caveat. Depending on user's operative system, and on the singularity intallation (e.g. if it is configured to bind automatically specific directories or not) the user may encounter a "Read-only file system" error. To avoid this error it is sufficient to specify the host binding directory as it was done to run GIP. In this case the syntax to execute giptools modules becomes the more verbose: ``singularity run -B /bindDir $giptools``. */bindDir* is any of the upstream directory including the **gipOut** directory. If the **gipOut** directory is in your local directory you can simply use ``-B $PWD``. *$giptools* is the giptools container.



