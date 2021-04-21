########################
Introduction to giptools
########################

To run GIP is required to specify in the configuration file the path to the **giptools**. This is a singularity container embedding all GIP scripts and software dependencies. But giptools is more than that. GIP focus is on each individual sample, and all results (e.g. the detection of a gene amplification or a SNV) are measured with respect to the genome reference. giptools can be used to address specific questions regarding the comparisons of samples subsets (e.g. compare the gene copy number between two samples in particular). 

The giptools file can be directly executed and used for the downstream comparison of samples processed by GIP. 
Typing ``giptools`` with no options (or with -h or --help) will list the available comparison modules.
To use any of these modules the syntax is simply ``giptools moduleName``. Each module has its own help documentation accessible with -h or --help. By default all comparison modules generate their output in the **gipOut/sampleComparison** folder.   

gipools modules:

+----------------+--------------------------------------------------------------------------+
| karyotype      |Compare the chromosome sequencing coverage distributions                  |
+----------------+--------------------------------------------------------------------------+
| binCNV         |Compare bin sequencing coverage in 2 samples                              |
+----------------+--------------------------------------------------------------------------+
| geCNV          |Compare gene sequencing coverage in 2 samples                             |
+----------------+--------------------------------------------------------------------------+
| ternary        |Compare gene sequencing coverage in 3 samples                             |
+----------------+--------------------------------------------------------------------------+
| ternaryBin     |Compare bin sequencing coverage in 3 samples                              |
+----------------+--------------------------------------------------------------------------+
| SNV            |Compare SNVs in multiple samples                                          |
+----------------+--------------------------------------------------------------------------+
| binDensity     |Density plot of bin sequencing coverage of many samples                   |
+----------------+--------------------------------------------------------------------------+
| geInteraction  |Detect CNV genes in many samples and produce correlation-based networks   |
+----------------+--------------------------------------------------------------------------+
| genomeDistance |Compare samples genomic distance                                          |
+----------------+--------------------------------------------------------------------------+
| phylogeny      |Extract the SNVs union and infer the phylogenetic tree                    |
+----------------+--------------------------------------------------------------------------+
| convergentCNV  |Detect convergent CNV gene amplifications                                 |
+----------------+--------------------------------------------------------------------------+
| overview       |Overview of the sequencing coverage of chromosomes, genomic bins and genes|
+----------------+--------------------------------------------------------------------------+
| panel          |Extract genomic information of a gene panel                               |
+----------------+--------------------------------------------------------------------------+

All modules require as an input the ``--gipOut`` parameter specifying the GIP output folder.
A description of each of the giptools modules and the available options is provided in the following pages of this documentation.


Troubleshooting
--------------- 

1. *Read-only file system* error: Depending on user's operative system and on the singularity intallation (e.g. if it is configured to bind automatically specific directories) the user may encounter a this error. A simple workaround is to specify the host binding directory as it was done to run GIP. In this case the syntax to execute giptools modules becomes the more verbose: ``singularity run -B /bindDir $giptools``. */bindDir* is any of the upstream directory including the **gipOut** directory. If the **gipOut** directory is in your local directory you can simply use ``-B $PWD``. *$giptools* is the giptools container.

2.  *caught segfault address (nil), cause memory not mapped*: Depending on the user's dataset the number of detected variants can be big and the default temporary directory can run out of space. Should this problem occur the user can manually specify a temporary directory. For instance  ``mkdir -p path/to/tmpDir ; export TMPDIR=/path/to/tmpDir``, then execute again the giptools module (e.g. ``singularity run -B /bindDir $giptools SNV``)
 




