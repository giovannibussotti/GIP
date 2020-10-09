##############
Worked example
##############


Data acquisition
----------------

GIP user may have access to private WGS data sets, or use publicly avaiablable data. In this tutorial we will consider 7 *Leishmania infantum* WGS samples available from the SRA database under the `PRJNA607007 <https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA607007>`_ accession.
An easy way to retrieve the accessions of the individual samples click on "Send to", then "File", then select "Accession list" format.
The resulting file will contain the following accessions:

| SRR11098642
| SRR11098643
| SRR11098644
| SRR11098645
| SRR11098646
| SRR11098647
| SRR11098648

The following code snippet requires the `SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc>`_ and can be copy/pasted in a bash terminal to recursively download each of the 7 experiments:

.. code-block:: bash

 accessions=(SRR11098642 SRR11098643 SRR11098644 SRR11098645 SRR11098646 SRR11098647 SRR11098648)
 mkdir data
 cd data
 for X in "${accessions[@]}" ; do
   #download from accession
   prefetch $X
   #convert to fastq
   fastq-dump --split-files $X
   #compress
   gzip ${X}*fastq
   #remove .sra downloaded file
   rm ~/ncbi/public/sra/${X}.sra
 done
 

Genome sequence (FASTA format) and gene annotations (GTF format) for *Leishmania infantum* are avaialable from multiple providers.
It is possible to downlod them from ENSEMBL protist with wget:

.. code-block:: bash

 #download
 wget ftp://ftp.ensemblgenomes.org/pub/release-29/protists/gtf/protists_euglenozoa1_collection/leishmania_infantum_jpcm5/Leishmania_infantum_jpcm5.GCA_000002875.2.29.gtf.gz
 wget ftp://ftp.ensemblgenomes.org/pub/release-29/protists/fasta/protists_euglenozoa1_collection/leishmania_infantum_jpcm5/dna/Leishmania_infantum_jpcm5.GCA_000002875.2.29.dna.toplevel.fa.gz
 #decompress
 gunzip Leishmania_infantum_jpcm5.GCA_000002875.2.29.gtf.gz
 gunzip Leishmania_infantum_jpcm5.GCA_000002875.2.29.dna.toplevel.fa.gz



GIP configuration
-----------------


  - prepare the configuration file
  - run 
  - modify parameter and re-run


