########################
Introduction to giptools
########################

To run GIP is required to specify in the configuration file the path to the giptools. This is a singularity container embedding all GIP scripts and software dependencies. But giptools is more than that. GIP focus is on each individual sample, and all results (e.g. the detection of a gene amplification or a SNV) are measured with respect to the genome reference. giptools can be used to address specific questions regarding the comparisons of samples subsets (e.g. compare the gene copy number between two samples in particular). 

The giptools file can be directly executed and used for the downstream comparison of samples processed by GIP. 
Typing ``giptools`` with no options (or with -h or --help) will list the available comparison modules.
To use any of these modules the syntax is simply ``giptools moduleName``. Each module has its own help documentation accessible with -h or --help. By default all comparison modules generate their output in the **gipOut/sampleComparison** folder.   





