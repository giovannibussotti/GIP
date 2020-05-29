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
    nextflow lgert.nf --genome ../inputData/dataset/Linf_test.fa --annotation ../inputData/dataset/Linf_test.ge.gtf -params-file sampInfo.yaml -c lgert.config
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
params.resultDir     = "lgertOut"
params.genome        = "genome.fa"
params.annotation    = "annotations.gtf"
params.repeatLibrary = "default"
genome               = file(params.genome)
annotation           = file(params.annotation)

ch = Channel.from(params.sampleList)
ch.into { ch1; ch2; ch3; ch4; ch5; ch6; ch7; ch8; ch9; ch10; ch11}

//config file variables
resultDir  = file(params.resultDir + "/files")
resultDir.with{mkdirs()}
MAPQ    = params.MAPQ
BITFLAG = params.BITFLAG
CHRSj   = params.CHRSj
plotCovPerNtOPT  = params.plotCovPerNtOPT
STEP             = params.STEP
PLOTcovPerBinOPT = params.PLOTcovPerBinOPT
PLOTcovPerBinRegressionOPT= params.PLOTcovPerBinRegressionOPT
covPerBinSigPeaksOPT      = params.covPerBinSigPeaksOPT
covPerGeMAPQoperation     = params.covPerGeMAPQoperation
plotCovPerGeOPT           = params.plotCovPerGeOPT
covPerGeSigPeaksOPT       = params.covPerGeSigPeaksOPT
freebayesOPT              = params.freebayesOPT
filterFreebayesOPT        = params.filterFreebayesOPT
filterDellyOPT   = params.filterDellyOPT
minNormCovForDUP = params.minNormCovForDUP
maxNormCovForDEL = params.maxNormCovForDEL
bigWigOpt        = params.bigWigOpt


process prepareGenome {
  publishDir resultDir

  output:
  file ("db") into bwaDb_ch1
  file("genome.chrSize") into chrSize_ch1
  set file("genome.fa") , file("genome.fa.fai") , file("genome.dict") , file("genome.chrSize") into (genome_ch1 , genome_ch2 , genome_ch3 , genome_ch4 , genome_ch5 , genome_ch6 , genome_ch7 , genome_ch8)
  file("genome.gaps.gz") into (gaps_ch1 , gaps_ch2 , gaps_ch3)
  file("repeatMasker") into (repeatMasker_ch1 , repeatMasker_ch2 , repeatMasker_ch3)
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


process map {
  publishDir resultDir 
  tag { "${SAMPLE.ID}" }
 
  input:
  val SAMPLE from ch1
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch1  
  val db from bwaDb_ch1

  output:
  set file ("${SAMPLE.ID}.bam") , file ("${SAMPLE.ID}.bam.bai") into (map1 , map2 , map3 , map4 , map5, map6 , map7 , map8 , map9)
  file ("${SAMPLE.ID}.MarkDup.log") into (mapDump1)
      
  """ 
  bash mapSample.sh $SAMPLE.ID $task.cpus $SAMPLE.FQ_DIR $fa $db/bwa/genome/ $SAMPLE.R1_FQID $SAMPLE.R2_FQID $SAMPLE.MULTIRUN_R1_FQIDS $SAMPLE.MULTIRUN_R2_FQIDS .    
  """
}

process covPerChr {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(bam) , file(bai) from map1
  val SAMPLE from ch2
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch2
  file gaps from gaps_ch1

  output:
  file ("chrCoverageMedians_${SAMPLE.ID}") into (covPerChr1 , covPerChr2 , covPerChr3 , covPerChr4, covPerChr5)

  """
  IFS=' ' read -r -a CHRS <<< "$CHRSj"  
  chrMedianCoverages $SAMPLE.ID $size $MAPQ $BITFLAG \$CHRS . $gaps
  """
}

process covPerNt {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(bam) , file(bai) from map2
  val SAMPLE from ch3
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch3
  
  output:
  file ("${SAMPLE.ID}.covPerNt.gz") into (covPerNt1 , covPerNt2)
  file ("${SAMPLE.ID}.pcMapqPerNt.gz") into (pcMapqPerNt1)
  set file ("${SAMPLE.ID}.covPerNt.boxplot.png") , file ("${SAMPLE.ID}.covPerNt.ridges.png") , file ("${SAMPLE.ID}.covPerNt.allMedians.tsv") , file ("${SAMPLE.ID}.covPerNt.medianGenomeCoverage") into (covPerNtDump1)

  script:
  """
  covPerNt $bam $size ${SAMPLE.ID}.covPerNt MEDIAN $MAPQ $BITFLAG 
  gzip ${SAMPLE.ID}.covPerNt
  pcMapqPerNt $bam $MAPQ ${SAMPLE.ID}.pcMapqPerNt
  Rscript /bin/plotGenomeCoverage_V3.R --files ${SAMPLE.ID}.covPerNt.gz --NAMES ${SAMPLE.ID} --DIR . --outName ${SAMPLE.ID}.covPerNt --pcMapqFiles ${SAMPLE.ID}.pcMapqPerNt.gz --chr $CHRSj $plotCovPerNtOPT
  """
}

process covPerBin {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(bam) , file(bai) from map3
  file(chrCoverageMedians) from covPerChr1
  val SAMPLE from ch4
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch4  

  output:
  file ("${SAMPLE.ID}.covPerBin.gz") into (covPerBin1)
  set file ("${SAMPLE.ID}.covPerBin.all.png") , file ("${SAMPLE.ID}.covPerBin_byChr.pdf") , file ("${SAMPLE.ID}.covPerBin.df.gz") , file ("${SAMPLE.ID}.covPerBin.extremeRatio.bed.gz") , file ("${SAMPLE.ID}.covPerBin.faceting.png") , file("${SAMPLE.ID}.covPerBin.sigPeaks.tsv.gz") , file("${SAMPLE.ID}.covPerBin.sigPeaks.stats") into (covPerBinDump1)
  //file ("${SAMPLE.ID}.gcLnorm.covPerBin.pdf") 

  """ 
  covPerBin $bam $size /bin/gencov2intervals.pl $STEP $MAPQ $BITFLAG . $chrCoverageMedians 

  Rscript /bin/covPerBin2loessGCnormalization_v2.R --ASSEMBLY $fa --DIR . --SAMPLE ${SAMPLE.ID} --outName ${SAMPLE.ID} 
  mv ${SAMPLE.ID}.gcLnorm.covPerBin.gz ${SAMPLE.ID}.covPerBin.gz

  Rscript /bin/compareGenomicCoverageBins.R --referenceName _fake_ --testName ${SAMPLE.ID} --referenceFile ${SAMPLE.ID}.covPerBin.gz --testFile ${SAMPLE.ID}.covPerBin.gz --outName ${SAMPLE.ID}.covPerBin --chrs $CHRSj --minMAPQ $MAPQ $PLOTcovPerBinOPT 

  #Rscript /bin/covPerBin2regression.R $PLOTcovPerBinRegressionOPT --filterChrsNames $CHRSj --filterMAPQ $MAPQ --DIR . --SAMPLES ${SAMPLE.ID}.covPerBin.gz  --NAMES ${SAMPLE.ID} --outName ${SAMPLE.ID}.PLOTcovPerBinRegression 

  Rscript /bin/sigPeaks_CLT.R --input ${SAMPLE.ID}.covPerBin.gz --outName ${SAMPLE.ID}.covPerBin.sigPeaks --minMAPQ $MAPQ $covPerBinSigPeaksOPT 
  """
}

process mappingStats {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(bam) , file(bai) from map4
  val SAMPLE from ch5
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch5

  output:
  set file ("${SAMPLE.ID}.alignmentMetrics.table") , file ("${SAMPLE.ID}.insertSize.histData") , file ("${SAMPLE.ID}.insertSize.hist.png") , file ("${SAMPLE.ID}.insertSize.table") into (mappingStatsDump1)

  """ 
  mkdir \$PWD/tmpDir
  #mapping stats
  java -jar /bin/picard.jar CollectAlignmentSummaryMetrics R=$fa I=$bam O=${SAMPLE.ID}.alignmentMetrics TMP_DIR=\$PWD/tmpDir
  #insert size
  java -jar /bin/picard.jar CollectInsertSizeMetrics I=$bam O=${SAMPLE.ID}.insertSize.metrics H=${SAMPLE.ID}.insertSize.pdf REFERENCE_SEQUENCE=$fa TMP_DIR=\$PWD/tmpDir
  rm -rf \$PWD/tmpDir ${SAMPLE.ID}.insertSize.pdf
  reformatMapStats.sh ${SAMPLE.ID}
  """
}

process covPerGe {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(bam) , file(bai) from map5
  file(chrCoverageMedians) from covPerChr2
  val SAMPLE from ch6
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch6
  file gaps from gaps_ch2
  file repeatMasker from repeatMasker_ch1
  file (annotation)

  output:
  file("${SAMPLE.ID}.covPerGe.gz") into covPerGe1 
  set file("${SAMPLE.ID}.covPerGe.sigPeaks.tsv") , file("${SAMPLE.ID}.covPerGe.sigPeaks.stats") , file("${SAMPLE.ID}.covPerGe.stats.df.gz") , file("${SAMPLE.ID}.covPerGe.stats.filtered.df.gz") , file("${SAMPLE.ID}.covPerGe.stats.cloud.png") into (covPerGeDump1)
  //, file("${SAMPLE.ID}.gcLnorm.covPerGe.pdf") ,  file("${SAMPLE.ID}.covPerGe.sigPeaks.pdf") , file("${SAMPLE.ID}.covPerGe.stats.pdf")

  """ 
  grep -v "^#" $repeatMasker/genome.out.gff | cut -f 1,4,5 > ${SAMPLE.ID}_tmp_reps
  
  covPerGe $bam ${SAMPLE.ID}.covPerGe $annotation $chrCoverageMedians $MAPQ $BITFLAG $covPerGeMAPQoperation $fa

  Rscript /bin/covPerGe2loessGCnormalization_v2.R --ASSEMBLY $fa --DIR . --SAMPLE ${SAMPLE.ID} --outName ${SAMPLE.ID}
  mv ${SAMPLE.ID}.gcLnorm.covPerGe.gz ${SAMPLE.ID}.covPerGe.gz

  Rscript /bin/compareGeneCoverage.R --NAMES fake ${SAMPLE.ID} --samples ${SAMPLE.ID}.covPerGe.gz ${SAMPLE.ID}.covPerGe.gz --outName ${SAMPLE.ID}.covPerGe.stats --repeats ${SAMPLE.ID}_tmp_reps --gaps $gaps --chrs $CHRSj --minMAPQ $MAPQ $plotCovPerGeOPT 

  Rscript /bin/sigPeaks_mixture.R --input ${SAMPLE.ID}.covPerGe.gz --outName ${SAMPLE.ID}.covPerGe.sigPeaks --minMAPQ $MAPQ $covPerGeSigPeaksOPT 

  rm -rf ${SAMPLE.ID}_tmp_reps
  """
}

process freebayes {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }
 
  input:
  set file(bam) , file(bai) from map6
  file(chrCoverageMedians) from covPerChr3
  val SAMPLE from ch7
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch7
  file repeatMasker from repeatMasker_ch2

  output:
  set file("${SAMPLE.ID}.vcf.gz") , file("${SAMPLE.ID}.vcf.gz.tbi") , file("${SAMPLE.ID}_freebayesFiltered") into (freebayes1)

  """ 
  freebayes -f $fa --min-mapping-quality $MAPQ $freebayesOPT --vcf ${SAMPLE.ID}.vcf $bam 
  bgzip ${SAMPLE.ID}.vcf
  tabix -p vcf -f ${SAMPLE.ID}.vcf.gz
  Rscript /bin/vcf2variantsFrequency_V4.R --selectedChrs $CHRSj --vcfFile ${SAMPLE.ID}.vcf.gz --chrCoverageMediansFile $chrCoverageMedians --chrSizeFile $size --outdir ./${SAMPLE.ID}_freebayesFiltered --reference $fa --discardGtfRegions $repeatMasker/genome.out.gff $filterFreebayesOPT
  """
}

process snpEff {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(vcf) , file(tbi) , file(freebayesFilteredDir) from freebayes1
  val SAMPLE from ch8
  file(snpEff) from snpEffDb_ch1
  
  output:
    set file("${SAMPLE.ID}.vcf.gz") , file("${SAMPLE.ID}.vcf.gz.tbi") , file("${SAMPLE.ID}_freebayesFiltered") , file("snpEff_summary_${SAMPLE.ID}.genes.txt.gz") , file("snpEff_summary_${SAMPLE.ID}.html") into (snpEffDump1)

  """ 
  #run snpEff on Freebayes output
  gunzip -c $vcf > tmpSnpEff_${SAMPLE.ID}.vcf
  java -jar /opt/snpEff/snpEff.jar -o gatk genome tmpSnpEff_${SAMPLE.ID}.vcf -noLog -c $snpEff/snpEff.config -stats snpEff_summary_${SAMPLE.ID} > ${SAMPLE.ID}.snpEff.vcf
  bgzip ${SAMPLE.ID}.snpEff.vcf
  gzip snpEff_summary_${SAMPLE.ID}.genes.txt
  mv snpEff_summary_${SAMPLE.ID} snpEff_summary_${SAMPLE.ID}.html
  mv ${SAMPLE.ID}.snpEff.vcf.gz ${SAMPLE.ID}.vcf.gz
  rm -rf $tbi
  tabix -p vcf -f ${SAMPLE.ID}.vcf.gz

  #run snpEff on the filtered SNVs
  gunzip -c ${SAMPLE.ID}_freebayesFiltered/singleVariants.vcf.gz > tmpSnpEff_${SAMPLE.ID}.vcf
  java -jar /opt/snpEff/snpEff.jar -o gatk genome tmpSnpEff_${SAMPLE.ID}.vcf -noLog -c $snpEff/snpEff.config -stats snpEff_summary_${SAMPLE.ID} > ${SAMPLE.ID}.snpEff.vcf
  bgzip ${SAMPLE.ID}.snpEff.vcf
  gzip -c snpEff_summary_${SAMPLE.ID}.genes.txt > ${SAMPLE.ID}_freebayesFiltered/snpEff_summary_${SAMPLE.ID}.genes.txt.gz
  mv snpEff_summary_${SAMPLE.ID} ${SAMPLE.ID}_freebayesFiltered/snpEff_summary_${SAMPLE.ID}.html
  mv ${SAMPLE.ID}.snpEff.vcf.gz ${SAMPLE.ID}_freebayesFiltered/singleVariants.vcf.gz
  rm -rf ${SAMPLE.ID}_freebayesFiltered/singleVariants.vcf.gz.tbi
  tabix -p vcf -f ${SAMPLE.ID}_freebayesFiltered/singleVariants.vcf.gz

  #extract the EFF field
  Rscript /bin/snpEffVcf2table.R --vcfFile ${SAMPLE.ID}_freebayesFiltered/singleVariants.vcf.gz --out ${SAMPLE.ID}_freebayesFiltered/singleVariants.df.EFF

  #update the singleVariants.df table with the EFF field
  gunzip ${SAMPLE.ID}_freebayesFiltered/singleVariants.df.gz 
  gunzip ${SAMPLE.ID}_freebayesFiltered/singleVariants.df.EFF
  awk '{print \$NF}' ${SAMPLE.ID}_freebayesFiltered/singleVariants.df.EFF > ${SAMPLE.ID}_freebayesFiltered/EFF
  paste ${SAMPLE.ID}_freebayesFiltered/singleVariants.df ${SAMPLE.ID}_freebayesFiltered/EFF > ${SAMPLE.ID}_freebayesFiltered/singleVariants.df2
  mv ${SAMPLE.ID}_freebayesFiltered/singleVariants.df2 ${SAMPLE.ID}_freebayesFiltered/singleVariants.df
  gzip ${SAMPLE.ID}_freebayesFiltered/singleVariants.df

  rm -rf tmpSnpEff_${SAMPLE.ID}.vcf snpEff_summary_${SAMPLE.ID}.genes.txt ${SAMPLE.ID}_freebayesFiltered/singleVariants.df.EFF ${SAMPLE.ID}_freebayesFiltered/EFF
  """
}

process dellySVref {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(bam) , file(bai) from map7
  file(chrCoverageMedians) from covPerChr4
  val SAMPLE from ch9
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch8
  file repeatMasker from repeatMasker_ch3
  file gaps from gaps_ch3

  output:
  file("*") into dump

  """
  #####
  #RUN#
  #####
  #prepare tmp bam with good HQ reads
  mkdir -p _tmpHQbam 
  samtools view -b -q $MAPQ -F $BITFLAG $bam > _tmpHQbam/${SAMPLE.ID}.bam
  samtools index _tmpHQbam/${SAMPLE.ID}.bam

  TYPES=(DEL DUP INV TRA)
  for T in "\${TYPES[@]}"; do
    delly -q $MAPQ -t \$T -g $fa -o ${SAMPLE.ID}_\${T}.vcf _tmpHQbam/${SAMPLE.ID}.bam
    bgzip ${SAMPLE.ID}_\${T}.vcf
    tabix -p vcf -f ${SAMPLE.ID}_\${T}.vcf.gz
  done
  rm -rf _tmpHQbam
  #concatenate SVs if they exist
  CONC=""
  for T in "\${TYPES[@]}"; do
    if [ -e ${SAMPLE.ID}_\${T}.vcf.gz ]; then
      CONC="\$CONC ${SAMPLE.ID}_\${T}.vcf.gz"
    fi
  done
  if [ -z \$CONC ]; then 
    touch ${SAMPLE.ID}.delly.vcf
  else 
    vcf-concat \$CONC > ${SAMPLE.ID}.delly.vcf
  fi

  rm -rf ${SAMPLE.ID}_DEL.vcf.gz ${SAMPLE.ID}_DUP.vcf.gz ${SAMPLE.ID}_INV.vcf.gz ${SAMPLE.ID}_TRA.vcf.gz ${SAMPLE.ID}_DEL.vcf.gz.tbi ${SAMPLE.ID}_DUP.vcf.gz.tbi ${SAMPLE.ID}_INV.vcf.gz.tbi ${SAMPLE.ID}_TRA.vcf.gz.tbi
  #sort and compress again
  cat ${SAMPLE.ID}.delly.vcf | vcf-sort > ${SAMPLE.ID}.delly.vcf.tmp
  mv ${SAMPLE.ID}.delly.vcf.tmp ${SAMPLE.ID}.delly.vcf 
  bgzip ${SAMPLE.ID}.delly.vcf 
  tabix -p vcf -f ${SAMPLE.ID}.delly.vcf.gz
  
  ########
  #Filter#
  ########
  DF=${SAMPLE.ID}.delly.filter
  mkdir -p \$DF
  grep -v "^#" $repeatMasker/genome.out.gff | cut -f 1,4,5 > ${SAMPLE.ID}.tabu.bed
  zcat $gaps >> ${SAMPLE.ID}.tabu.bed
  SVTYPES=(DEL DUP INV)
  for SV in "\${SVTYPES[@]}"; do
    Rscript /bin/filterDelly.R --SVTYPE \$SV --test $SAMPLE.ID --reference $SAMPLE.ID --vcfFile ${SAMPLE.ID}.delly.vcf.gz --chrSizeFile $size --badSequencesBed ${SAMPLE.ID}.tabu.bed --useENDfield yes --outName \$DF/${SAMPLE.ID}.delly.\$SV
  done
  Rscript /bin/filterDelly.R --SVTYPE TRA --test $SAMPLE.ID --reference $SAMPLE.ID --vcfFile ${SAMPLE.ID}.delly.vcf.gz --chrSizeFile $size --badSequencesBed ${SAMPLE.ID}.tabu.bed --useENDfield no --outName \$DF/${SAMPLE.ID}.delly.TRA
  rm ${SAMPLE.ID}.tabu.bed
  #add coverage
  mkdir -p \${DF}/coverages
  for X in `ls \$DF | grep ".bed\$" | sed -e 's/.bed\$//'`; do
    cat \$DF/\${X}.bed | perl -ne 'if(\$_=~/^(\\S+)\\s+(\\S+)\\s+(\\S+)/){\$chr=\$1;\$start=\$2;\$end=\$3;\$ge="\${chr}_\${start}_\${end}"; print \"\${chr}\\tdelly\\tSV\\t\${start}\\t\${end}\\t.\\t.\\t.\\tgene_id \\"\$ge\\"; transcript_id \\"\$ge\\";\\n\"} ' > \$DF/_tmp.gtf
    covPerGe $bam \$DF/coverages/\${X}.cov \$DF/_tmp.gtf chrCoverageMedians_$SAMPLE.ID $MAPQ $BITFLAG $covPerGeMAPQoperation $fa
    rm \$DF/_tmp.gtf
    #filter by coverage
    filterByCov \$X \$DF/coverages/\${X}.fbc $minNormCovForDUP $maxNormCovForDEL $SAMPLE.ID \$DF
    #ov with genes
    ovWithGenes \$X \$DF
    mv \$DF/coverages/\${X}.ov \${X}.filter
    bedForCircos \${X}.filter \${X}.filter.circosBed
  done
  rm -rf \$DF
  """
}

process bigWigGenomeCov {
  publishDir resultDir
  tag { "${SAMPLE.ID}" }

  input:
  set file(bam) , file(bai) from map8
  val SAMPLE from ch10
  
  output:
  file("${SAMPLE.ID}.bw") into (bigWigGenomeCovDump1)

  """
  #generate a bedGraph per chr
  TMP=tmpbw
  mkdir -p \$TMP
  for CHR in `samtools idxstats $bam | cut -f1 | head -n -2`; do
    samtools view  -b $bam \$CHR > \$TMP/\${CHR}.bam    
    samtools index \$TMP/\${CHR}.bam
    bamCoverage -b \$TMP/\${CHR}.bam --outFileName \$TMP/\${CHR}.bg --numberOfProcessors $task.cpus --outFileFormat bedgraph --normalizeUsingRPKM --ignoreDuplicates -r \$CHR $bigWigOpt
  done

  #combine chrs and turn to bigWig
  mkdir -p \$TMP/sort
  cat \$TMP/*.bg | sort -k1,1 -k2,2n -T \$TMP/sort > \$TMP/${SAMPLE.ID}.bg
  samtools idxstats $bam | cut -f 1,2 | head -n -2 > \$TMP/chrSize
  /bin/bedGraphToBigWig \$TMP/${SAMPLE.ID}.bg \$TMP/chrSize \$TMP/${SAMPLE.ID}.bw
  mv \$TMP/${SAMPLE.ID}.bw .
  rm -rf \$TMP
  """
}

/*
process covPerClstr {
  publishDir resultDir

  input:
  file('*') from map9.collect()
  //val('*') SAMPLE from ch11
  file('*') from covPerChr5.collect()
  file('*') from covPerGe1.collect()

  output:
  file ('*') into lalaland

  """
  ls > booo
  #bash /bin/covPerClstr.sh -f "$covPerGe" -n "$NAMES" -b "$BAMS" -c "$CHRM" -m $minMAPQ -g $GEGTF -a $ASSEMBLY -l 0.9 -s 0.9 -o covPerClstr
  """
}
*/



