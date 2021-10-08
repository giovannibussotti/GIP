#############################
Input and output organization
#############################

Input
-----

GIP requires the following mandatory input parameters:

+----------------+-----------------------------------+
| \-\-genome     | multi-FASTA genome reference file |
+----------------+-----------------------------------+
| \-\-annotation | gene coordinates file (GTF format)|
+----------------+-----------------------------------+
| \-\-index      | list of sequencing data files     |
+----------------+-----------------------------------+
| \-c            | nextflow configuration file       |
+----------------+-----------------------------------+

| The input **genome** file can be a normal text file or gzip compressed file. The chromosome identifiers can contain white spaces and extra information other than the chromosome identifier itself (e.g. supercontig identifier). However GIP will consider as chromosome identifier just the characters string before the first white space. The user must ensure that the chromosome identifiers are the same between the **genome** and the **annotation** files, considering that GIP does not consider any characters coming after the first white space (if any).
| The **annotation** file must be in standard GTF format, reporting *gene* or *exon* features (or both) in the third field. If available the GTF file can provide *CDS* entries.
| All additional GIP parameters can be passed with the command line execution, or set in the **gip.config** configuration file.
| The configuration file hosts the default values of all parameters under the ``params{}`` scope.
| The ``process{}`` scope can be used to customize the configuration of all GIP processes, including the allocation of memory or CPUs.
| Other important process parameters include:

* ``executor``       - indicates whether the processes must be executed on the local machine (default) or on a computing cluster (e.g. 'slurm');
* ``clusterOptions`` - provides optional cluster configurations, like the nodes partition where to allocate jobs (e.g. '-p hubbioit --qos hubbioit');
* ``container``      - sets the absolute path of the giptools singularity image to use to execute the processes.

|  The ``singularity{}`` scope contains the configuration options to interface the nextflow pipeline with the `singularity container <https://www.nextflow.io/docs/latest/singularity.html>`_. Singularity allows to mount (a.k.a. bind) host input data at specific container locations defined by the user. GIP container (i.e. the giptools file) comes with a set of built-in folders that can be used as access points for data:

* /fq
* /genome
* /annotation
* /repLib
* /geneFunction
* /gipOut
* /mnt


| The ``runOptions`` parameter can be used to specify bind points with the \-\-bind (or -B) option and the following syntax:
| *'-B host_directory:container_binding_point'*.
| If the container binding point is omitted, this will be considered the same as the host directory.
| Then, the user must have the caution to specify all the input parameters not relative to the host system, but relative to where the data is visible in the container.
| For instance, the user can mount the host folder containing the genome file (e.g. /home/user/data/assemblies/reference.fa) to the /genome container folder by specifying ``runOtions='-B /home/user/data/assemblies:/genome'`` in the configuration file.
| Then, when executing GIP, the user can simply pass the input genome command line with ``--genome /genome/reference.fa``.
| Multiple host directories can be mounted with additional -B directives. For instance, to mount also the working directory (e.g. /home/projects) and the directory containing the sequencing data:

.. code-block:: guess

  singularity {
    enabled    = true
    runOptions = '-B /home/projects -B /home/user/data/assemblies:/genome -B /home/user/sequencingData:/fq'
  }

| Alternatively, it can be convenient to set ``autoMounts=true`` and bind just a top-level folder of all data folders.

.. code-block:: guess

  singularity {
    enabled    = true
    autoMounts = true
    runOptions = '--bind /pasteur'
  }

| By doing that the file paths inside the container and in the host will be identical, and the user can provide all the input files with the normal host paths.
| Please refer to the `Nextflow documentation <https://www.nextflow.io/docs/latest/config.html>`_ for the pipeline configuration file to discover all available options.


| The index file must comply with the following syntax rules:

1. Tsv format (i.e. <Tab> separated)
2. First header row with the labels: sampleId   read1    read2
3. All the following rows must indicate the sample identifier, and the file names first and second pair-end sequencing data files in fastq.gz format
4. To combine the reads originating from multiple technical replicates the fastq.gz files must be comma-separated and in the same order between read1 and read2

| Example:
| sampleId        read1    read2
| sample1 /fq/s1.r1.fastq.gz  /fq/s1.r2.fastq.gz
| sample2 /fq/s2.RUN1.r1.fastq.gz,/fq/s2.RUN2.r1.fastq.gz  /fq/s2.RUN1.r2.fastq.gz,/fq/s2.RUN2.r2.fastq.gz



Output
------

| GIP results are accessible from the **gipOut/** output directory which contains the following subfolders:

+------------------+-----------------------------+
| **genome/**      | reference genome data       |
+------------------+-----------------------------+
| **samples/**     | individual samples results  |
+------------------+-----------------------------+
| **covPerClstr/** | gene cluster quantification |
+------------------+-----------------------------+
| **reports/**     | report files                |
+------------------+-----------------------------+

| The *report* process executed at the end of the pipeline returns .html files in the **reports/** subfolder, summarizing main results and figures for each sample, like :download:`this example <../_static/LIPA83.html>`.
| All the other files in the **gipOut/** directory are symbolic links to the data cached in the **work/** directory, which in turn is organized in subfolders named with the hexadecimal numbers identifying the executed processes.
| Thanks to the Nextflow implementation the user can easily test different GIP parameterization without the need to re-execute the entire pipeline. Just by adding ``-resume`` to the command line, GIP will re-run only the processes that are affected by the parameter change, and use the cached results of all the other processes.

| The ``--resultDir`` parameter can be used to set a name alternative to "gipOut" for the result directory.


In the following part we provide a description of GIP steps operated by the Nextflow processes and all result files.

