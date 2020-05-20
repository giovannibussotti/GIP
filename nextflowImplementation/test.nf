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
    nextflow run lgert.nf -params-file sampInfo.yaml -c nextflow_configs/lgert.config -resume
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
params.resultDir = "lgertOut"
params.MAPQ = 5
params.BITFLAG = 1028
params.CHRSj = "1 2" 
params.plotCovPerNtOPT="--ylim 8"
params.STEP=300
params.PLOTcovPerBinOPT="--segmentation no --highRatio 2 --lowRatio 0.5 --ylim 3 --minWindowNormMean 0.1 --minWindowNormMedian 0.1 --divideByOne"
params.PLOTcovPerBinRegressionOPT="--resizeValue 1000 --filterChrEnds 0 --maxCov 5 --cv"
params.covPerBinSigPeaksOPT="--verbose --inFormat covPerBin --arbitraryHighLowPeakThresh 2 0.5 --minLen 0 --pThresh 0.001 --padjust BY"
params.covPerGeMAPQoperation="MEAN"
params.plotCovPerGeOPT="--divideByOne --locusFieldName locus --scoreFieldName normalizedMeanCoverage --plot3_min 2 --plot3_max 5 --plot_highRatio 2 --plot_lowRatio 0.5 --scoreLabel coverage --scaleFree yes --dumpFilteredGenes"
params.covPerGeSigPeaksOPT="--pFilter --pThresh 0.001 --padjust BH --inFormat covPerGe --arbitraryHighLowPeakThresh 2 0.5 --minLen 0"
params.freebayesOPT="--read-indel-limit 1 --read-mismatch-limit 3 --read-snp-limit 3 --hwe-priors-off --binomial-obs-priors-off --allele-balance-priors-off  --min-alternate-fraction 0.05 --min-base-quality 5 --min-alternate-count 2 --pooled-continuous"
params.filterFreebayesOPT="--variants single --snvOnly yes --minFreq 0.1 --maxFreq 1.1 --minAO 2 --minAOhomopolymer 20 --contextSpan 5 --homopolymerFreq 0.4 --minMQMR 20 --minMQM 20 --useENDfield no --howManyMads 4 --generatePseudoReference yes"
//params.filterDellyOPT="--filterSVatTheTelomericEnds 10000 --minPE 20 --minPercentDVDRtest 20 --maxPercentDVDRreference 99999 --minDVtest 20 --maxDVreference 99999 --PRECISE no --maxBadSeqPercent 50"
params.filterDellyOPT="--filterSVatTheTelomericEnds 1 --minPE 0 --minPercentDVDRtest 0 --maxPercentDVDRreference 99999 --minDVtest 0 --maxDVreference 99999 --PRECISE no --maxBadSeqPercent 100"
params.minNormCovForDUP = 2 
params.maxNormCovForDEL = 0.5
params.bigWigOpt="--binSize 10 --smoothLength 30"


//create variables you can use on as many processed as you want
summaryDir=params.resultDir
resultDir = file(params.resultDir + "/files")
resultDir.with{mkdirs()}
MAPQ = params.MAPQ
BITFLAG     = params.BITFLAG
CHRSj = params.CHRSj
//CHRS  = CHRSj.tokenize(' ')  //string to list, but this is not then interpreted as a bash array
plotCovPerNtOPT=params.plotCovPerNtOPT
STEP=params.STEP
PLOTcovPerBinOPT=params.PLOTcovPerBinOPT
PLOTcovPerBinRegressionOPT=params.PLOTcovPerBinRegressionOPT
covPerBinSigPeaksOPT=params.covPerBinSigPeaksOPT
covPerGeMAPQoperation=params.covPerGeMAPQoperation
plotCovPerGeOPT=params.plotCovPerGeOPT
covPerGeSigPeaksOPT=params.covPerGeSigPeaksOPT
freebayesOPT=params.freebayesOPT
filterFreebayesOPT=params.filterFreebayesOPT
filterDellyOPT=params.filterDellyOPT
minNormCovForDUP=params.minNormCovForDUP
maxNormCovForDEL=params.maxNormCovForDEL
bigWigOpt=params.bigWigOpt




// Configurable variables
// nextflow lgert.nf --genome ../inputData/dataset/Linf_test.fa --annotation ../inputData/dataset/Linf_test.ge.gtf -params-file sampInfo.yaml -c lgert.config

params.genome        = "genome.fa"
params.annotation    = "annotations.gtf"
params.repeatLibrary = "default"
genome     = file(params.genome)
annotation = file(params.annotation)

ch = Channel.from(params.sampleList)
ch.into { ch1; ch2; ch3; ch4; ch5; ch6; ch7; ch8; ch9; ch10}

process prepareGenome {
  publishDir resultDir

  output:

  file ("bo") into bwaDb_ch

  """
  realpath $annotation > bo
  touch /home/gbussott/Desktop/L-GERT_githubRepo/inputData/dataset/testaBlu
  """
}



