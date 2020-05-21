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

  file ("db") into bwaDb_ch1
  file("genome.chrSize") into chrSize_ch1
  set file("genome.fa") , file("genome.fa.fai") , file("genome.dict") , file("genome.chrSize") into (genome_ch1 , genome_ch2)
  file("genome.gaps.gz") into gaps_ch1
  file("repeatMasker") into repeatMasker_ch1
  file("snpEff")  into snpEffDb_ch1

  script:
  if( params.repeatLibrary == 'default' )
  """
  A-prepareGenome.sh -f $genome -x $annotation -c $task.cpus
  """

  else
  """
  A-prepareGenome.sh -f $genome -x $annotation -c $task.cpus -l params.repeatLibrary
  """
}


/*
process testa {
  publishDir resultDir
  
  input:
  //set file(fa) , file(fai) , file(dict) , file(size) from genome_ch
  set fa , fai , dict , size from genome_ch1
  
  output:

  file oraaaa into letters

  """
  cat $size > oraaaa
  """
}
*/


process map {
  publishDir resultDir 
    
  input:
  val SAMPLE from ch1
  set fa , fai , dict , size from genome_ch1  
  val db from bwaDb_ch1

  output:
  set file ("${SAMPLE.ID}.bam") , file ("${SAMPLE.ID}.bam.bai") into (map1 , map2 , map3 , map4 , map5, map6 , map7 , map8)
  file ("${SAMPLE.ID}.MarkDup.log") into (mapDump1)
      
  """ 
  bash mapSample.sh $SAMPLE.ID $task.cpus $SAMPLE.FQ_DIR $fa $db/bwa/genome/ $SAMPLE.R1_FQID $SAMPLE.R2_FQID $SAMPLE.MULTIRUN_R1_FQIDS $SAMPLE.MULTIRUN_R2_FQIDS .    
  """
}

process covPerChr {
  publishDir resultDir

  input:
  set file(bam) , file(bai) from map1
  val SAMPLE from ch2
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch2
  file gaps from gaps_ch1

  output:
  file ("chrCoverageMedians_${SAMPLE.SAMPLE_ID}") into (covPerChr1 , covPerChr2 , covPerChr3 , covPerChr4)

  """
  IFS=' ' read -r -a CHRS <<< "$CHRSj"  
  chrMedianCoverages $SAMPLE.ID $size $MAPQ $BITFLAG \$CHRS . $gaps
  """
}

/*
process covPerNt {
  publishDir resultDir

  input:
  set file(bam) , file(bai) from map2
  val SAMPLE from ch3
  
  output:
  file ("${SAMPLE.SAMPLE_ID}.covPerNt.gz") into (covPerNt1 , covPerNt2)
  file ("${SAMPLE.SAMPLE_ID}.pcMapqPerNt.gz") into (pcMapqPerNt1)
  file ("${SAMPLE.SAMPLE_ID}.covPerNt.boxPlot.pdf") into (covPerNtDump1)

  script:
  """
  covPerNt $bam /mnt/data/$SAMPLE.CHRSIZE ${SAMPLE.SAMPLE_ID}.covPerNt MEDIAN $MAPQ $BITFLAG 
  gzip ${SAMPLE.SAMPLE_ID}.covPerNt
  pcMapqPerNt $bam $MAPQ ${SAMPLE.SAMPLE_ID}.pcMapqPerNt
  Rscript /bin/plotGenomeCoverage_V2.R --files ${SAMPLE.SAMPLE_ID}.covPerNt.gz --NAMES ${SAMPLE.SAMPLE_ID} --DIR . --outName ${SAMPLE.SAMPLE_ID}.covPerNt.boxPlot --pcMapqFiles ${SAMPLE.SAMPLE_ID}.pcMapqPerNt.gz --chr $CHRSj $plotCovPerNtOPT
  """
}

process covPerBin {
  publishDir resultDir

  input:
  set file(bam) , file(bai) from map3
  file(chrCoverageMedians) from covPerChr1
  val SAMPLE from ch4
  
  output:
  file ("${SAMPLE.SAMPLE_ID}.covPerBin.gz") into (covPerBin1)
  set file ("${SAMPLE.SAMPLE_ID}.gcLnorm.covPerBin.pdf") , file ("${SAMPLE.SAMPLE_ID}.PLOTcovPerBin_all.pdf") , file ("${SAMPLE.SAMPLE_ID}.PLOTcovPerBin_byChr.pdf") , file ("${SAMPLE.SAMPLE_ID}.PLOTcovPerBin.df.gz") , file ("${SAMPLE.SAMPLE_ID}.PLOTcovPerBin.extremeRatio.bed.gz") , file ("${SAMPLE.SAMPLE_ID}.PLOTcovPerBin_faceting.pdf") into (covPerBinDump1)

  script:
  """ 
  covPerBin $bam /mnt/data/$SAMPLE.CHRSIZE /bin/gencov2intervals.pl $STEP $MAPQ $BITFLAG . $chrCoverageMedians 

  Rscript /bin/covPerBin2loessGCnormalization_v2.R --ASSEMBLY /mnt/data/$SAMPLE.ASSEMBLY --DIR . --SAMPLE ${SAMPLE.SAMPLE_ID} --outName ${SAMPLE.SAMPLE_ID} 
  mv ${SAMPLE.SAMPLE_ID}.gcLnorm.covPerBin.gz ${SAMPLE.SAMPLE_ID}.covPerBin.gz

  Rscript /bin/compareGenomicCoverageBins.R --referenceName _fake_ --testName ${SAMPLE.SAMPLE_ID} --referenceFile ${SAMPLE.SAMPLE_ID}.covPerBin.gz --testFile ${SAMPLE.SAMPLE_ID}.covPerBin.gz --outName ${SAMPLE.SAMPLE_ID}.PLOTcovPerBin --chrs $CHRSj --minMAPQ $MAPQ $PLOTcovPerBinOPT 

  Rscript /bin/covPerBin2regression.R $PLOTcovPerBinRegressionOPT --filterChrsNames $CHRSj --filterMAPQ $MAPQ --DIR . --SAMPLES ${SAMPLE.SAMPLE_ID}.covPerBin.gz  --NAMES ${SAMPLE.SAMPLE_ID} --outName ${SAMPLE.SAMPLE_ID}.PLOTcovPerBinRegression 

  Rscript /bin/sigPeaks_CLT.R --input ${SAMPLE.SAMPLE_ID}.covPerBin.gz --outName ${SAMPLE.SAMPLE_ID}.covPerBin.sigPeaks --minMAPQ $MAPQ $covPerBinSigPeaksOPT 
  """
}

process mappingStats {
  publishDir resultDir

  input:
  set file(bam) , file(bai) from map4
  val SAMPLE from ch5
  
  output:
  file ("${SAMPLE.SAMPLE_ID}.stats") into (mappingStatsDump1)

  script:
  """ 
  Rscript /bin/mappingStats.R --bams $bam --dir . --assembly /mnt/data/$SAMPLE.ASSEMBLY --outName NA --tmpDir ./_tmpDirCollectAlignmentSummaryMetrics_$SAMPLE.ASSEMBLY  --CollectAlignmentSummaryMetrics "java -jar /bin/picard.jar CollectAlignmentSummaryMetrics"
  """
}

process covPerGe {
  publishDir resultDir

  input:
  set file(bam) , file(bai) from map5
  file(chrCoverageMedians) from covPerChr2
  val SAMPLE from ch6
  
  output:
  set file("${SAMPLE.SAMPLE_ID}.covPerGe.gz") , file("${SAMPLE.SAMPLE_ID}.covPerGe.sigPeaks.pdf") , file("${SAMPLE.SAMPLE_ID}.covPerGe.sigPeaks.sigPeaks.tsv") , file("${SAMPLE.SAMPLE_ID}.covPerGe.sigPeaks.stats") , file("${SAMPLE.SAMPLE_ID}.covPerGe.stats.df.gz") , file("${SAMPLE.SAMPLE_ID}.covPerGe.stats.filtered.df.gz") , file("${SAMPLE.SAMPLE_ID}.covPerGe.stats.pdf") , file("${SAMPLE.SAMPLE_ID}.gcLnorm.covPerGe.pdf")  into (covPerGeDump1)


  script:
  """ 
  grep -v "^#" /mnt/data/$SAMPLE.REPS | cut -f 1,4,5 > ${SAMPLE.SAMPLE_ID}_tmp_reps
  
  covPerGe $bam ${SAMPLE.SAMPLE_ID}.covPerGe /mnt/data/$SAMPLE.GENES $chrCoverageMedians $MAPQ $BITFLAG $covPerGeMAPQoperation /mnt/data/$SAMPLE.ASSEMBLY

  Rscript /bin/covPerGe2loessGCnormalization_v2.R --ASSEMBLY /mnt/data/$SAMPLE.ASSEMBLY --DIR . --SAMPLE ${SAMPLE.SAMPLE_ID} --outName ${SAMPLE.SAMPLE_ID}
  mv ${SAMPLE.SAMPLE_ID}.gcLnorm.covPerGe.gz ${SAMPLE.SAMPLE_ID}.covPerGe.gz

  Rscript /bin/compareGeneCoverage.R --NAMES fake ${SAMPLE.SAMPLE_ID} --samples ${SAMPLE.SAMPLE_ID}.covPerGe.gz ${SAMPLE.SAMPLE_ID}.covPerGe.gz --outName ${SAMPLE.SAMPLE_ID}.covPerGe.stats --repeats ${SAMPLE.SAMPLE_ID}_tmp_reps --gaps /mnt/data/$SAMPLE.GAPS --chrs $CHRSj --minMAPQ $MAPQ $plotCovPerGeOPT 

  Rscript /bin/sigPeaks_mixture.R --input ${SAMPLE.SAMPLE_ID}.covPerGe.gz --outName ${SAMPLE.SAMPLE_ID}.covPerGe.sigPeaks --minMAPQ $MAPQ $covPerGeSigPeaksOPT 

  rm -rf ${SAMPLE.SAMPLE_ID}_tmp_reps
  """
}

process freebayes {
 publishDir resultDir

  input:
  set file(bam) , file(bai) from map6
  file(chrCoverageMedians) from covPerChr3
  val SAMPLE from ch7
  
  output:
  set file("${SAMPLE.SAMPLE_ID}.vcf.gz") , file("${SAMPLE.SAMPLE_ID}.vcf.gz.tbi") , file("${SAMPLE.SAMPLE_ID}_freebayesFiltered")   into (freebayes1)

  script:
  """ 
  freebayes -f /mnt/data/$SAMPLE.ASSEMBLY --min-mapping-quality $MAPQ $freebayesOPT --vcf ${SAMPLE.SAMPLE_ID}.vcf $bam 

  bgzip ${SAMPLE.SAMPLE_ID}.vcf

  tabix -p vcf -f ${SAMPLE.SAMPLE_ID}.vcf.gz

  Rscript /bin/vcf2variantsFrequency_V4.R --selectedChrs $CHRSj --vcfFile ${SAMPLE.SAMPLE_ID}.vcf.gz --chrCoverageMediansFile $chrCoverageMedians --chrSizeFile /mnt/data/$SAMPLE.CHRSIZE --outdir ./${SAMPLE.SAMPLE_ID}_freebayesFiltered --reference /mnt/data/$SAMPLE.ASSEMBLY --discardGtfRegions /mnt/data/$SAMPLE.REPS $filterFreebayesOPT
  """
}

process snpEff {
 publishDir resultDir

  input:
  set file(vcf) , file(tbi) , file(freebayesFilteredDir) from freebayes1
  val SAMPLE from ch8
  
  output:
    set file("${SAMPLE.SAMPLE_ID}.vcf.gz") , file("${SAMPLE.SAMPLE_ID}.vcf.gz.tbi") , file("${SAMPLE.SAMPLE_ID}_freebayesFiltered") , file("snpEff_summary_${SAMPLE.SAMPLE_ID}.genes.txt.gz") , file("snpEff_summary_${SAMPLE.SAMPLE_ID}.html") into (snpEffDump1)

  script:
  """ 
  #run snpEff on Freebayes output
  gunzip -c $vcf > tmpSnpEff_${SAMPLE.SAMPLE_ID}.vcf
  java -jar /opt/snpEff/snpEff.jar -o gatk $SAMPLE.SNPEFF_DB tmpSnpEff_${SAMPLE.SAMPLE_ID}.vcf -noLog -c /mnt/data/$SAMPLE.snpEffConfig -stats snpEff_summary_${SAMPLE.SAMPLE_ID} > ${SAMPLE.SAMPLE_ID}.snpEff.vcf
  bgzip ${SAMPLE.SAMPLE_ID}.snpEff.vcf
  gzip snpEff_summary_${SAMPLE.SAMPLE_ID}.genes.txt
  mv snpEff_summary_${SAMPLE.SAMPLE_ID} snpEff_summary_${SAMPLE.SAMPLE_ID}.html
  mv ${SAMPLE.SAMPLE_ID}.snpEff.vcf.gz ${SAMPLE.SAMPLE_ID}.vcf.gz
  rm -rf $tbi
  tabix -p vcf -f ${SAMPLE.SAMPLE_ID}.vcf.gz

  #run snpEff on the filtered SNVs
  gunzip -c ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.vcf.gz > tmpSnpEff_${SAMPLE.SAMPLE_ID}.vcf
  java -jar /opt/snpEff/snpEff.jar -o gatk $SAMPLE.SNPEFF_DB tmpSnpEff_${SAMPLE.SAMPLE_ID}.vcf -noLog -c /mnt/data/$SAMPLE.snpEffConfig -stats snpEff_summary_${SAMPLE.SAMPLE_ID} > ${SAMPLE.SAMPLE_ID}.snpEff.vcf
  bgzip ${SAMPLE.SAMPLE_ID}.snpEff.vcf
  gzip -c snpEff_summary_${SAMPLE.SAMPLE_ID}.genes.txt > ${SAMPLE.SAMPLE_ID}_freebayesFiltered/snpEff_summary_${SAMPLE.SAMPLE_ID}.genes.txt.gz
  mv snpEff_summary_${SAMPLE.SAMPLE_ID} ${SAMPLE.SAMPLE_ID}_freebayesFiltered/snpEff_summary_${SAMPLE.SAMPLE_ID}.html
  mv ${SAMPLE.SAMPLE_ID}.snpEff.vcf.gz ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.vcf.gz
  rm -rf ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.vcf.gz.tbi
  tabix -p vcf -f ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.vcf.gz

  #extract the EFF field
  Rscript /bin/snpEffVcf2table.R --vcfFile ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.vcf.gz --out ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df.EFF

  #update the singleVariants.df table with the EFF field
  gunzip ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df.gz 
  gunzip ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df.EFF
  awk '{print \$NF}' ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df.EFF > ${SAMPLE.SAMPLE_ID}_freebayesFiltered/EFF
  paste ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df ${SAMPLE.SAMPLE_ID}_freebayesFiltered/EFF > ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df2
  mv ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df2 ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df
  gzip ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df

  rm -rf tmpSnpEff_${SAMPLE.SAMPLE_ID}.vcf snpEff_summary_${SAMPLE.SAMPLE_ID}.genes.txt ${SAMPLE.SAMPLE_ID}_freebayesFiltered/singleVariants.df.EFF ${SAMPLE.SAMPLE_ID}_freebayesFiltered/EFF
  """
}

process dellySVref {
 publishDir resultDir

  input:
  set file(bam) , file(bai) from map7
  file(chrCoverageMedians) from covPerChr4
  val SAMPLE from ch9
  
  output:
  set file("${SAMPLE.SAMPLE_ID}.delly.DEL.filter") , file("${SAMPLE.SAMPLE_ID}.delly.DEL.filter.circosBed") , file("${SAMPLE.SAMPLE_ID}.delly.DUP.filter") , file("${SAMPLE.SAMPLE_ID}.delly.DUP.filter.circosBed") , file("${SAMPLE.SAMPLE_ID}.delly.INV.filter") , file("${SAMPLE.SAMPLE_ID}.delly.INV.filter.circosBed") , file("${SAMPLE.SAMPLE_ID}.delly.TRA.filter") , file("${SAMPLE.SAMPLE_ID}.delly.TRA.filter.circosBed") , file("${SAMPLE.SAMPLE_ID}.delly.vcf.gz") , file("${SAMPLE.SAMPLE_ID}.delly.vcf.gz.tbi")  into (dellySVrefDump1)

  script:
  """
  bash /bin/dellySVref.sh ${SAMPLE.SAMPLE_ID} /mnt/data/$SAMPLE.REPS /mnt/data/$SAMPLE.ASSEMBLY . /mnt/data/$SAMPLE.CHRSIZE $MAPQ $BITFLAG /mnt/data/$SAMPLE.GAPS '$CHRSj' '$filterDellyOPT' $covPerGeMAPQoperation $minNormCovForDUP $maxNormCovForDEL /mnt/data/$SAMPLE.GENES
  """
}

process bigWigGenomeCov {
 publishDir resultDir

  input:
  set file(bam) , file(bai) from map8
  val SAMPLE from ch10
  
  output:
  file("${SAMPLE.SAMPLE_ID}.bw") into (bigWigGenomeCovDump1)

  script:
  """
  #generate a bedGraph per chr
  TMP=tmpbw
  mkdir -p \$TMP
  for CHR in `samtools idxstats $bam | cut -f1 | head -n -2`; do
    samtools view  -b $bam \$CHR > \$TMP/\${CHR}.bam    
    samtools index \$TMP/\${CHR}.bam
    bamCoverage -b \$TMP/\${CHR}.bam --outFileName \$TMP/\${CHR}.bg --numberOfProcessors $task.cpus --outFileFormat bedgraph --normalizeUsing RPKM --ignoreDuplicates -r \$CHR $bigWigOpt
  done

  #combine chrs and turn to bigWig
  mkdir -p \${TMP}/sort
  cat \${TMP}/*.bg | sort -k1,1 -k2,2n -T \${TMP}/sort > \${TMP}/${SAMPLE.SAMPLE_ID}.bg
  samtools idxstats ${SAMPLE.SAMPLE_ID}.bam | cut -f 1,2 | head -n -2 > \${TMP}/chrSize
  /bin/bedGraphToBigWig \${TMP}/${SAMPLE.SAMPLE_ID}.bg \${TMP}/chrSize \${TMP}/${SAMPLE.SAMPLE_ID}.bw
  mv \${TMP}/${SAMPLE.SAMPLE_ID}.bw .
  rm -rf \$TMP
  """
}


*/



