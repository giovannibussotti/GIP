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

The following code snippet requires the `SRA Toolkit <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc>`_ and can be copy/pasted in a bash terminal to recursively download each of the 7 WGS experiment:

.. code-block:: bash

 accessions=(SRR11098642 SRR11098643 SRR11098644 SRR11098645 SRR11098646 SRR11098647 SRR11098648)
 mkdir fastqs
 cd fastqs
 for X in "${accessions[@]}" ; do
   #download from accession
   prefetch $X
   #convert to fastq
   fastq-dump --split-files $X
   #compress
   gzip ${X}*fastq
   #remove .sra file
   rm -rf $X
 done
 cd ..
 

Depending on your species of interest the user may need to refer to different specialized genomic data banks to retrieve the genome sequence and the gene coordinates, which are both required GIP inputs (respectively --genome and --annotation).
Additionally, the user may want to specify gene function annotations (--geneFunction parameter). Depending on the species, gene function annotation may be not completely available, and manual curation or data integration from different repositories may be needed.
In the code snipped below we use the ENSEMBL protists FTP server to download the latest *Leishmaina infantum* genome sequence and annotations in .gff3 format, then reformat the latter to obtain both the gene coordinates file (.gtf format) and the gene function files.

.. code-block:: bash

 mkdir data
 cd data
 #download genome
 wget ftp://ftp.ensemblgenomes.org/pub/release-48/protists/fasta/protists_euglenozoa1_collection/leishmania_infantum_gca_900500625/dna/Leishmania_infantum_gca_900500625.LINF.dna.toplevel.fa.gz  
 #donwload gene annotations
 wget ftp://ftp.ensemblgenomes.org/pub/release-48/protists/gff3/protists_euglenozoa1_collection/leishmania_infantum_gca_900500625/Leishmania_infantum_gca_900500625.LINF.48.gff3.gz
 #convert
 gff3=Leishmania_infantum_gca_900500625.LINF.48.gff3.gz
 perl -e '
 open (F,">geneFunction.tsv") or die "cannot open geneFunction.tsv: $!";
 open (G,">annotation.gtf") or die "cannot open annotation.gtf: $!";
 open(IN, "gunzip -c '$gff3' |") or die "gunzip '$gff3': $!";
 while(<IN>){
  if($_=~/^(.*)ID=gene:([^;]+).*description=([^;]+)/){
   my $a   = $1;
   my $id  = $2;
   my $des = $3;
   print G "${a}gene_id \"$id\"; transcript_id \"$id\";\n";
   print F "$id\t$des\n";

  }
 }
 close F;
 close G;' 
 cd ..


All the required data is now available. 
 


GIP configuration
-----------------

The user should prepare the index file indicating the sample names and the respective sequencing data files.
The index is a tab separated file with the following heading row: sampleId	read1	read2 
In this example we will use the sample names as reported in `PRJNA607007 <https://www.ncbi.nlm.nih.gov/sra/?term=PRJNA607007>`_ and the data file paths as they are in the host system.  
So if the fastq files are stored in the ``/pasteur/data/fq`` the end the index file should look like :download:`this <../_static/sampleIndexExample.pdf>`.


Next, the user must edit the GIP configuration file (i.e. **gip.config**) according to the available computing resources, and most importantly, binding an up-level directory containing all the data paths. In this example the ``/pasteur`` would be a good choice. 
Assuming that the user copied giptools locally as ``/pasteur/p2p5/gip/giptools``, and that he/she wants to execute GIP on a slurm cluster with special partition and quality of service options ``-p aTeam --qos fast`` while keeping the default for all the rest, the parameters that need to be updated are:

* ``executor='slurm'``
* ``container='/pasteur/p2p5/gip/giptools'`` 
* ``clusterOptions='-p aTeam --qos fast'``
* ``runOptions = '--bind /pasteur'``

The edited configuration file should look like :download:`this <../_static/gip.config.workedExample.pdf>`.
If instead the user cannot take advantage of a computing cluster, then he/she can run GIP locally by simply specifying ``executor='local'``.


GIP execution
-------------

To run GIP:

.. code-block:: bash

 nextflow gip --genome /mnt/Leishmania_infantum_gca_900500625.LINF.dna.toplevel.fa.gz \
              --annotation /mnt/annotation.gtf \
              --geneFunctions /mnt/geneFunction.tsv \
              --index index.tsv \
              -c gip.config

Assuming that the user now wants to test a stricter parametrization for the SNV filtering he/che can edit the filterFreebayesOPT parameter and increase the --minFreq value to 0.3. Then, when re-exeute GIP, it is possible to add the ``-resume`` option to the command line to compute just the processes affected by the parameter change. GIP will re-use the cached results of all the other processes.






