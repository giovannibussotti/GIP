#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                   L-GERT : Leishmania-Genome Reporting Tool 
========================================================================================
 Quantify compare and visualize Leishmania genomic variants . Started March 2019.
 #### Homepage / Documentation
 https://github.com/giovannibussotti/L-GERT
 #### Authors
 Giovanni Bussotti <giovanni.bussotti@pasteur.fr>
----------------------------------------------------------------------------------------
*/

version = 1.1
def helpMessage() {
    log.info"""
    =========================================
     L-GERT : Leishmania-Genome Reporting Tool v${version}
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow lgert.nf --genome ../inputData/dataset/Linf_test.fa --annotation ../inputData/dataset/Linf_test.ge.gtf --index index.tsv -c lgert.config -resume
    Mandatory arguments:
      -params-file                   sample metadata yaml file 
      -c                             nextflow configuration file
    General Options:
      -resultDir                     result directory 
      -MAPQ                          read MAPQ cut-off  
      -BITFLAG                       SAM bitflag filter
      -CHRSj                         List of chromosome identifiers to consider (included in quotes)
    Karyotype Options:
      -plotCovPerNtOPT               karyptype boxplot plotting options
    Bin Coverage Options:  
      -STEP                          bin size
      -PLOTcovPerBinOPT              coverage plotting options 
      -PLOTcovPerBinRegressionOPT    coverage regression plotting options
      -covPerBinSigPeaksOPT          identify statistically significant CNV wrt the reference  
    Gene Coverage Options:      
      -covPerGeMAPQoperation         Measure average gene MAPQ  
      -plotCovPerGeOPT               Gene coverage plotting options
      -covPerGeSigPeaksOPT           identify statistically significant gene CNV wrt the reference
    SNV Options:  
      -freebayesOPT                  Freebayes options
      -filterFreebayesOPT            Downstream SNV quality filters
    SV Options: 
      -filterDellyOPT                Downastream, Delly filtering options
      -minNormCovForDUP              min normalized sequencing coverage for Delly duplications
      -maxNormCovForDEL              max normalized sequencing coverage for Delly deletions
    BigWig Options:  
      -bigWigOpt                     Options to generate the coverage bigWig file
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


// Configurable variables
params.dummyIndex1   = 'dummyIndex1.tsv'
params.dummyIndex2   = 'dummyIndex2.tsv'


Channel
    .fromPath(params.dummyIndex1)
    .splitCsv(header:true , sep:'\t')
    .map{ row-> tuple(row.sampleId, row.read1, row.read2) }
    .set { ch1  }

Channel
    .fromPath(params.dummyIndex2)
    .splitCsv(header:true , sep:'\t')
    .map{ row-> tuple(row.sampleId, row.read1, row.read2) }
    .set { ch2  }



process map {
  input:
  set sampleId , read1 , read2 from ch1    

  output:
  set val(sampleId) , file("${sampleId}*") into (map)

  script:
  """
  echo $sampleId > ${sampleId}.bam
  echo $sampleId > ${sampleId}.bam.bai
  """
}

process covPerGe {
  input:
  set sampleId , read1 , read2 from ch2

  output:
  set val(sampleId) , file("${sampleId}*") into (covPerGe)
  
  script:
  """
  echo $sampleId > ${sampleId}.covPerGe
  echo $sampleId > ${sampleId}.covPerGe.pdf
  """
}

map
  .flatten()
  .set{mapF}

covPerGe
  .flatten()
  .set{covPerGeF}


process joinTheRightSamples {
  input:
  file ('*') from mapF.join(covPerGeF)
  
  output:
  file ('*') into (lalaland)

  script:
  """ 
  ls > out
  """
}


