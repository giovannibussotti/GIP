.. GIP documentation master file, created by
   sphinx-quickstart on Mon Sep 14 10:59:47 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.
   
   image:: _static/logo.jpg
    :width: 120px
    :alt: GIP logo
    :align: center


Genome Instability Pipeline
===========================

| The Genome Instability Pipeline (GIP) will solve your problem of analyzing whole genome sequencing (WGS) data making bioinformatics screenings faster and easier.

| GIP allows the batch processing of multiple WGS experiments, including the mapping of short reads, the quantification of chromosomes, genes and genomic bins.

| The results of GIP are summarized in a *report* page providing sequencing experiment and mapping statistics, graphical representations of genomic features quantifications, and excel tables.

GIP also enables **giptools**, a tool-suite allowing to  **Detect**, **Compare** and **Visualize** 

* Aneuploidy changes
* Gene copy number variants (CNVs)
* Single nucleotide variants (SNVs)
* Structural variants (SVs)

Look how easy it is to use.
Run GIP on a large set of input WGS experiments (e.g. 200 samples)

  ``gip --genome fasta.file --annotation gtf.file --index fastqs.tsv -c gip.config``

  
Compare gene CNVs in 3 specific isolates

  ``giptools ternary --gipOut gipOutDir --samples Lmj_A445 Lmj_1948 Ltr_16``


.. figure:: _static/triangle1.png
    :align: center
    :width: 500px



.. `Gallery karyotype figure <https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click on image to zoom&p=PMC3&id=6222132_mbo0051841230003.jpg>`_



Features
--------

* Scalable workflow implemented in Nextflow
* Reproducible thanks to the Singularity support
* Easy-to-use and easy-to-install
* Support for multi-copy gene clusters
* Custom comparison of samples sub-sets
* Publication quality figures and excel tables

Requirements
------------

* `Singularity`_ 3.5.2+
* `Nextflow`_ 20.04.1.5335+

.. _Singularity: https://sylabs.io/guides/3.5/user-guide/quick_start.html#quick-installation-steps
.. _Nextflow: https://www.nextflow.io

These are the specific versions that were tested. Use other versions at your own risk.


Installation
------------

1) Download GIP

    ``github download command...``

2) Pull the container

    ``singularity pull command...``

3) Move **giptools** to /usr/local/bin/. If the user does not have permission on the folder she/he can keep giptools in any other location, and just update the "container" parameter in the gip.config: ``container='/Path/To/giptools'``

| Done!
| To make GIP available anywhere, the user can place **gip** in a directory that is in PATH.
  
    
Contribute
----------

- Issue Tracker: github.com/GIP/issues
- Source Code: github.com/GIP

Support
-------

If you are having issues, please let us know.
We have a mailing list located at: GIP@google-groups.com

License
-------

GIP and giptools modules are licensed under a Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International (CC BY-NC-ND 4.0) License. These software are provided without any warranty. Third parties software and scripts installed in giptools have their own copyright and license. A list is provided in the next page of this documentation. Users of giptools must take notice and agree with the terms and conditions of all third party software.


.. toctree::
   :maxdepth: 2
   :hidden:

   software/index

.. toctree::
   :maxdepth: 2
   :hidden:

   GIP/index


.. toctree::
   :maxdepth: 2
   :hidden:

   giptools/index



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


