#!/usr/bin/env nextflow
/*
vim: syntax=groovy
-*- mode: groovy;-*-
========================================================================================
                   GIP : GENOME INSTABILITY PIPELINE
========================================================================================
 Quantify compare and visualize Leishmania genomic variants . Started March 2019.
 #### Homepage / Documentation
 https://github.com/giovannibussotti/GIP
 #### Authors
 Giovanni Bussotti <giovanni.bussotti@pasteur.fr>
----------------------------------------------------------------------------------------
*/

version = 1.1
def helpMessage() {
    log.info"""
    =========================================
     GIP : Genome Instability Pipeline v${version}
    =========================================
    
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow gip.nf --genome ../inputData/dataset/Linf_test.fa --annotation ../inputData/dataset/Linf_test.ge.gtf --index index.tsv -c gip.config -resume
    
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
      -covPerGeRepeatRange           Visualize repeats within this distance from significant genes
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

// check mandatory options
//if (!params.genomeIndex && !params.genome) {
//    exit 1, "Reference genome not specified"
//}

// Configurable variables
params.resultDir     = "lgertOut"
params.genome        = "genome.fa"
params.annotation    = "annotations.gtf"
params.repeatLibrary = "default"
params.index         = 'index.tsv'
params.geneFunction  = 'NA'
genome               = file(params.genome)
annotation           = file(params.annotation)
geneFunction         = file(params.geneFunction)
repeatLibrary        = file(params.repeatLibrary)

//ch = Channel.from(params.sampleList)
//ch.into { ch1 }

Channel
    .fromPath(params.index)
    .splitCsv(header:true , sep:'\t')
    .map{ row-> tuple(row.sampleId, row.read1, row.read2) }
    .set { ch1  }

//config file variables
resultDir  = file(params.resultDir + "/samples")
resultDir.with{mkdirs()}
MAPQ    = params.MAPQ
BITFLAG = params.BITFLAG
CHRSj   = params.CHRSj
plotCovPerNtOPT  = params.plotCovPerNtOPT
STEP             = params.STEP
PLOTcovPerBinOPT = params.PLOTcovPerBinOPT
PLOTcovPerBinRegressionOPT= params.PLOTcovPerBinRegressionOPT
covPerBinSigPeaksOPT      = params.covPerBinSigPeaksOPT
covPerGeRepeatRange       = params.covPerGeRepeatRange
covPerGeSigPeaksOPT       = params.covPerGeSigPeaksOPT
freebayesOPT              = params.freebayesOPT
filterFreebayesOPT        = params.filterFreebayesOPT
filterDellyOPT   = params.filterDellyOPT
minNormCovForDUP = params.minNormCovForDUP
maxNormCovForDEL = params.maxNormCovForDEL
bigWigOpt        = params.bigWigOpt


process processGeneFunction {
  output:
  file 'geneFunction.tsv' into (geFun_ch , geFun_ch1 , geFun_ch2)

  script:
  if( params.geneFunction == 'NA'  )
  """
  cat $annotation | perl -ne 'if(\$_=~/gene_id \"([^\"]+)\"/){print "\$1\tNA\n";}' > geneFunction.tsv
  """

  else 
  """
  cp $geneFunction geneFunction.tsv
  """
}

process prepareGenome {
  publishDir "$params.resultDir/genome"

  input:
  file(annotation)
  file(genome)

  output:
  file ("db") into bwaDb_ch1
  file("genome.chrSize") into chrSize_ch1
  set file("genome.fa") , file("genome.fa.fai") , file("genome.dict") , file("genome.chrSize") into (genome_ch1 , genome_ch2 , genome_ch3 , genome_ch4 , genome_ch5 , genome_ch6 , genome_ch7 , genome_ch8, genome_ch9)
  file("genome.gaps.gz") into (gaps_ch1 , gaps_ch2)
  file("repeatMasker") into (repeatMasker_ch1 , repeatMasker_ch2 , repeatMasker_ch3)
  file("snpEff")  into snpEffDb_ch1
  
  script:
  if( params.repeatLibrary == 'default' )
  """
  A-prepareGenome.sh -f $genome -x $annotation -c $task.cpus
  """

  else
  """
  A-prepareGenome.sh -f $genome -x $annotation -c $task.cpus -l $repeatLibrary
  """
}


process map { 
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }
 
  input:
  set sampleId , read1 , read2 from ch1
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch1  
  file(db) from bwaDb_ch1

  output:
  set val(sampleId) , file ("${sampleId}.bam") , file ("${sampleId}.bam.bai") into (map1 , map2 , map3)
  set val(sampleId) , file ("${sampleId}.MarkDup.table") , file ("${sampleId}.MarkDup.histData") , file ("${sampleId}.MarkDup.hist.png") into (mapDump1)
      
  """ 
  mapSample.sh $sampleId $task.cpus $fa $db/bwa/genome/ $read1 $read2 .    
  reformatMarkDup.sh $sampleId
  """
}


process covPerChr {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(bam) , file(bai) from map1
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch2
  file gaps from gaps_ch1

  output:
  set val(sampleId), file(bam) , file(bai), file("chrCoverageMedians_$sampleId") into (covPerChr1 , covPerChr2 , covPerChr3, covPerChr4 )

  """
  IFS=' ' read -r -a CHRS <<< "$CHRSj"  
  chrMedianCoverages $sampleId $size $MAPQ $BITFLAG \$CHRS . $gaps
  """
}

process covPerNt {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) from covPerChr1
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch3
  
  output:
  set val(sampleId), file ("${sampleId}.covPerNt.gz") , file ("${sampleId}.pcMapqPerNt.gz") , file ("${sampleId}.covPerNt.boxplot.png") , file ("${sampleId}.covPerNt.ridges.png") , file ("${sampleId}.covPerNt.allMedians.tsv") , file ("${sampleId}.covPerNt.medianGenomeCoverage") into (covPerNt)
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) , file ("${sampleId}.covPerNt.gz") into (covPerNt4delly)

  script:
  """
  covPerNt $bam $size ${sampleId}.covPerNt MEDIAN 0 $BITFLAG 
  gzip -f ${sampleId}.covPerNt
  pcMapqPerNt $bam $MAPQ ${sampleId}.pcMapqPerNt
  Rscript /bin/plotGenomeCoverage_V3.R --files ${sampleId}.covPerNt.gz --NAMES ${sampleId} --DIR . --outName ${sampleId}.covPerNt --pcMapqFiles ${sampleId}.pcMapqPerNt.gz --chr $CHRSj $plotCovPerNtOPT
  """
}

process covPerBin {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) from covPerChr2
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch4  

  output:
  set val(sampleId), file ("${sampleId}.covPerBin.gz"), file ("${sampleId}.covPerBin.plot.all.png") , file ("${sampleId}.covPerBin.plot.byChr.pdf") , file ("${sampleId}.covPerBin.plot.tsv.gz") , file ("${sampleId}.covPerBin.plot.faceting.png") , file("${sampleId}.covPerBin.significant.bins.tsv.gz") , file("${sampleId}.covPerBin.significant.segments.tsv.gz") , file("${sampleId}.covPerBin.significant.stats") into (covPerBin)
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) , file ("${sampleId}.covPerBin.gz") into (covPerBin4covPerGe) 
  
  //file ("${SAMPLE.ID}.gcLnorm.covPerBin.pdf") 

  """ 
  covPerBin $bam $size /bin/gencov2intervals.pl $STEP 0 $BITFLAG . $covPerChr 

  Rscript /bin/covPerBin2loessGCnormalization_v2.R --ASSEMBLY $fa --DIR . --SAMPLE ${sampleId} --outName ${sampleId} 
  mv ${sampleId}.gcLnorm.covPerBin.gz ${sampleId}.covPerBin.gz

  Rscript /bin/sigPeaks_CLT.R --input ${sampleId}.covPerBin.gz --outName ${sampleId}.covPerBin.significant --minMAPQ $MAPQ $covPerBinSigPeaksOPT

  Rscript /bin/plotCovPerBin.R --covPerBin ${sampleId}.covPerBin.gz --outName ${sampleId}.covPerBin.plot --chrs $CHRSj --significant ${sampleId}.covPerBin.significant.bins.tsv.gz --chrSizeFile genome.chrSize --minMAPQ $MAPQ
  """
}

process mappingStats {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(bam) , file(bai) from map2
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch5

  output:
  set val(sampleId), file ("${sampleId}.alignmentMetrics.table") , file ("${sampleId}.insertSize.histData") , file ("${sampleId}.insertSize.hist.png") , file ("${sampleId}.insertSize.table") into (mappingStats)

  """ 
  mkdir \$PWD/tmpDir
  #mapping stats
  java -jar /bin/picard.jar CollectAlignmentSummaryMetrics R=$fa I=$bam O=${sampleId}.alignmentMetrics TMP_DIR=\$PWD/tmpDir
  #insert size
  java -jar /bin/picard.jar CollectInsertSizeMetrics I=$bam O=${sampleId}.insertSize.metrics H=${sampleId}.insertSize.pdf REFERENCE_SEQUENCE=$fa TMP_DIR=\$PWD/tmpDir
  rm -rf \$PWD/tmpDir ${sampleId}.insertSize.pdf
  reformatMapStats.sh ${sampleId}
  """
}

process covPerGe {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) , file(covPerBin) from covPerBin4covPerGe
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch6
  file repeatMasker from repeatMasker_ch1
  file (annotation)
  file geFun from geFun_ch1

  output:
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) , file("${sampleId}.covPerGe.gz") into covPerGe1 
  set val(sampleId) , file("${sampleId}.covPerGe.significant.tsv") , file("${sampleId}.covPerGe.significant.stats") , file("${sampleId}.covPerGeKaryoplot") into (covPerGeDump1)
  //, file("${sampleId}.gcLnorm.covPerGe.pdf") ,  file("${sampleId}.covPerGe.significant.pdf") , file("${sampleId}.covPerGe.stats.pdf")

  """ 
  covPerGe $bam ${sampleId}.covPerGe $annotation $covPerChr 0 $BITFLAG MEAN $fa

  Rscript /bin/covPerGe2loessGCnormalization_v2.R --ASSEMBLY $fa --DIR . --SAMPLE ${sampleId} --outName ${sampleId}
  mv ${sampleId}.gcLnorm.covPerGe.gz ${sampleId}.covPerGe.gz

  Rscript /bin/sigPeaks_mixture.R --input ${sampleId}.covPerGe.gz --outName ${sampleId}.covPerGe.significant --minMAPQ $MAPQ $covPerGeSigPeaksOPT

  Rscript /bin/karyoplotCovPerGe.R --covPerGe ${sampleId}.covPerGe.gz --covPerBin $covPerBin --chrSize $size --CHRS $CHRSj --REPS $repeatMasker/genome.out.gff --significant ${sampleId}.covPerGe.significant.tsv --outDir ${sampleId}.covPerGeKaryoplot --repeatRange $covPerGeRepeatRange --minMAPQ $MAPQ --geneFunction $geFun
  """
}

process freebayes {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }
 
  input:
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) from covPerChr3
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch7
  file repeatMasker from repeatMasker_ch2

  output:
  set val(sampleId) , file("${sampleId}.vcf.gz") , file("${sampleId}.vcf.gz.tbi") , file("${sampleId}_freebayesFiltered") into (freebayes1)

  """ 
  freebayes -f $fa --min-mapping-quality $MAPQ $freebayesOPT --vcf ${sampleId}.vcf $bam 
  bgzip ${sampleId}.vcf
  tabix -p vcf -f ${sampleId}.vcf.gz
  Rscript /bin/vcf2variantsFrequency_V4.R --selectedChrs $CHRSj --vcfFile ${sampleId}.vcf.gz --chrCoverageMediansFile $covPerChr --chrSizeFile $size --outdir ./${sampleId}_freebayesFiltered --reference $fa --discardGtfRegions $repeatMasker/genome.out.gff $filterFreebayesOPT
  """
}

process snpEff {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(vcf) , file(tbi) , file(freebayesFilteredDir) from freebayes1
  file(snpEff) from snpEffDb_ch1
  file geFun from geFun_ch2 
 
  output:
    set val(sampleId) , file("${sampleId}.vcf.gz") , file("${sampleId}.vcf.gz.tbi") , file("${sampleId}_freebayesFiltered") , file("snpEff_summary_${sampleId}.genes.txt.gz") , file("snpEff_summary_${sampleId}.html") into (snpEff)

  """ 
  #run snpEff on Freebayes output
  gunzip -c $vcf > tmpSnpEff_${sampleId}.vcf
  java -jar /opt/snpEff/snpEff.jar -ud 0 -o gatk genome tmpSnpEff_${sampleId}.vcf -noLog -c $snpEff/snpEff.config -stats snpEff_summary_${sampleId} > ${sampleId}.snpEff.vcf
  bgzip ${sampleId}.snpEff.vcf
  gzip -f snpEff_summary_${sampleId}.genes.txt
  mv snpEff_summary_${sampleId} snpEff_summary_${sampleId}.html
  mv ${sampleId}.snpEff.vcf.gz ${sampleId}.vcf.gz
  rm -rf $tbi
  tabix -p vcf -f ${sampleId}.vcf.gz

  #run snpEff on the filtered SNVs
  gunzip -c ${sampleId}_freebayesFiltered/singleVariants.vcf.gz > tmpSnpEff_${sampleId}.vcf
  java -jar /opt/snpEff/snpEff.jar -ud 0 -o gatk genome tmpSnpEff_${sampleId}.vcf -noLog -c $snpEff/snpEff.config -stats snpEff_summary_${sampleId} > ${sampleId}.snpEff.vcf
  bgzip ${sampleId}.snpEff.vcf
  gzip -c snpEff_summary_${sampleId}.genes.txt > ${sampleId}_freebayesFiltered/snpEff_summary_${sampleId}.genes.txt.gz
  mv snpEff_summary_${sampleId} ${sampleId}_freebayesFiltered/snpEff_summary_${sampleId}.html
  mv ${sampleId}.snpEff.vcf.gz ${sampleId}_freebayesFiltered/singleVariants.vcf.gz
  rm -rf ${sampleId}_freebayesFiltered/singleVariants.vcf.gz.tbi
  tabix -p vcf -f ${sampleId}_freebayesFiltered/singleVariants.vcf.gz

  #extract the EFF field
  Rscript /bin/snpEffVcf2table.R --vcfFile ${sampleId}_freebayesFiltered/singleVariants.vcf.gz --out ${sampleId}_freebayesFiltered/singleVariants.df.EFF

  #update the singleVariants.df table with the EFF field
  gunzip ${sampleId}_freebayesFiltered/singleVariants.df.gz 
  gunzip ${sampleId}_freebayesFiltered/singleVariants.df.EFF
  awk '{print \$NF}' ${sampleId}_freebayesFiltered/singleVariants.df.EFF > ${sampleId}_freebayesFiltered/EFF
  paste ${sampleId}_freebayesFiltered/singleVariants.df ${sampleId}_freebayesFiltered/EFF > ${sampleId}_freebayesFiltered/singleVariants.df2
  mv ${sampleId}_freebayesFiltered/singleVariants.df2 ${sampleId}_freebayesFiltered/singleVariants.df

  #compute dnds
  Rscript /bin/dndsRatio.R --ALLGEANN $geFun --SNPEFFDF ${sampleId}_freebayesFiltered/singleVariants.df --ID ${sampleId}
  mv ${sampleId}_cleanEFF.tsv ${sampleId}_freebayesFiltered/singleVariants.df
  gzip -f ${sampleId}_freebayesFiltered/singleVariants.df  
  mv ${sampleId}_dNdStables.tsv  ${sampleId}_freebayesFiltered/dNdStable.tsv
  gzip -f ${sampleId}_freebayesFiltered/dNdStable.tsv
  mv ${sampleId}_cleanEFF.stats ${sampleId}_freebayesFiltered/dNdS.stats  

  rm -rf tmpSnpEff_${sampleId}.vcf snpEff_summary_${sampleId}.genes.txt ${sampleId}_freebayesFiltered/singleVariants.df.EFF ${sampleId}_freebayesFiltered/EFF
  """
}

process dellySVref {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(bam) , file(bai) , file(covPerChr) , file (covPerNt) from covPerNt4delly
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch8
  file repeatMasker from repeatMasker_ch3
  file gaps from gaps_ch2
  file (annotation)

  output:
  set val(sampleId), file("${sampleId}.delly.vcf.gz") , file("${sampleId}.delly.vcf.gz.tbi") , file("${sampleId}_dellyFiltered")  into delly

  """
  #prepare tmp bam with good HQ reads 
  mkdir -p _tmpHQbam 
  samtools view -b -q $MAPQ -F $BITFLAG $bam > _tmpHQbam/${sampleId}.bam
  samtools index _tmpHQbam/${sampleId}.bam

  #call SVs
  TYPES=(DEL DUP INV TRA)
  for T in "\${TYPES[@]}"; do
    delly -q $MAPQ -t \$T -g $fa -o ${sampleId}_\${T}.vcf _tmpHQbam/${sampleId}.bam
    bgzip ${sampleId}_\${T}.vcf
    tabix -p vcf -f ${sampleId}_\${T}.vcf.gz
  done

  #concatenate SVs if they exist
  CONC=""
  for T in "\${TYPES[@]}"; do
    if [ -e ${sampleId}_\${T}.vcf.gz ]; then
      CONC="\$CONC ${sampleId}_\${T}.vcf.gz"
    fi
  done
  if [ -z "\$CONC" ]; then 
    touch ${sampleId}.delly.vcf
  else 
    vcf-concat \$CONC > ${sampleId}.delly.vcf
  fi

  rm -rf ${sampleId}_DEL.vcf.gz ${sampleId}_DUP.vcf.gz ${sampleId}_INV.vcf.gz ${sampleId}_TRA.vcf.gz ${sampleId}_DEL.vcf.gz.tbi ${sampleId}_DUP.vcf.gz.tbi ${sampleId}_INV.vcf.gz.tbi ${sampleId}_TRA.vcf.gz.tbi _tmpHQbam
  
  #sort and compress again
  cat ${sampleId}.delly.vcf | vcf-sort > ${sampleId}.delly.vcf.tmp
  mv ${sampleId}.delly.vcf.tmp ${sampleId}.delly.vcf 
  bgzip ${sampleId}.delly.vcf 
  tabix -p vcf -f ${sampleId}.delly.vcf.gz
  
  ########
  #Filter#
  ########
  #run filterDelly
  DF=${sampleId}.delly.filter
  mkdir -p \$DF
  grep -v "^#" $repeatMasker/genome.out.gff | cut -f 1,4,5 > ${sampleId}.tabu.bed
  zcat $gaps >> ${sampleId}.tabu.bed
  SVTYPES=(DEL DUP INV)
  for SV in "\${SVTYPES[@]}"; do
    Rscript /bin/filterDelly.R --SVTYPE \$SV --test $sampleId --reference $sampleId --vcfFile ${sampleId}.delly.vcf.gz --chrSizeFile $size --badSequencesBed ${sampleId}.tabu.bed --useENDfield yes --outName \$DF/${sampleId}.delly.\$SV
  done
  Rscript /bin/filterDelly.R --SVTYPE TRA --test $sampleId --reference $sampleId --vcfFile ${sampleId}.delly.vcf.gz --chrSizeFile $size --badSequencesBed ${sampleId}.tabu.bed --useENDfield no --outName \$DF/${sampleId}.delly.TRA
  rm ${sampleId}.tabu.bed
  #convert bed to gtf
  mkdir -p \${DF}/coverages
  for X in `ls \$DF | grep ".bed\$" | sed -e 's/.bed\$//'`; do
    cat \$DF/\${X}.bed | perl -ne 'if(\$_=~/^(\\S+)\\s+(\\S+)\\s+(\\S+)/){\$chr=\$1;\$start=\$2;\$end=\$3;\$ge="\${chr}_\${start}_\${end}"; print \"\${chr}\\tdelly\\tSV\\t\${start}\\t\${end}\\t.\\t.\\t.\\tgene_id \\"\$ge\\"; transcript_id \\"\$ge\\";\\n\"} ' > \$DF/_tmp.gtf
    #evaluate coverage
    covPerGe $bam \$DF/coverages/\${X}.cov \$DF/_tmp.gtf $covPerChr $MAPQ $BITFLAG MEAN $fa
    rm \$DF/_tmp.gtf
    #filter by coverage
    filterByCov \$X \$DF/coverages/\${X}.fbc $minNormCovForDUP $maxNormCovForDEL $sampleId \$DF
    #ov with genes
    ovWithGenes \$X \$DF
    mv \$DF/coverages/\${X}.ov \${X}.filter
    bedForCircos \${X}.filter \${X}.filter.circosBed
  done
  rm -rf \$DF
  
  #circos
  cp /bin/ideogram.conf /bin/ticks.conf .
  prepare4Circos.sh -i $sampleId -m $bam -s $size -g $annotation -c "$CHRSj" -n $covPerNt -b 25000 -t /bin/templateCircos.conf -z ${sampleId}.delly.TRA.filter.circosBed -j ${sampleId}.delly.INV.filter.circosBed -k ${sampleId}.delly.DUP.filter.circosBed -w ${sampleId}.delly.DEL.filter.circosBed

  #reorganize output
  mkdir -p ${sampleId}_dellyFiltered/
  mv ${sampleId}.delly.*filter* ${sampleId}_dellyFiltered/
  mv ${sampleId}.SV.circos.png ${sampleId}_dellyFiltered/
  mv ${sampleId}_circosData ${sampleId}_dellyFiltered/
  """
}

process bigWigGenomeCov {
  publishDir "$params.resultDir/samples/$sampleId"
  tag { "${sampleId}" }

  input:
  set val(sampleId), file(bam) , file(bai) from map3
  
  output:
  set val(sampleId), file("${sampleId}.bw") into (bigWigGenomeCov)

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
  cat \$TMP/*.bg | sort -k1,1 -k2,2n -T \$TMP/sort > \$TMP/${sampleId}.bg
  samtools idxstats $bam | cut -f 1,2 | head -n -2 > \$TMP/chrSize
  /bin/bedGraphToBigWig \$TMP/${sampleId}.bg \$TMP/chrSize \$TMP/${sampleId}.bw
  mv \$TMP/${sampleId}.bw .
  rm -rf \$TMP
  """
}


process covPerClstr {
  publishDir "$params.resultDir/covPerClstr"

  input:
  file('*') from covPerGe1.collect()  
  set file(fa) , file(fai) , file(dict) , file(size) from genome_ch9
  file (annotation)
  file geFun from geFun_ch

  output:
  file ('*') into covPerClstrDump

  """
  #all input in the same order
  COVPERGE=''
  NAMES=''
  CHRM=''
  BAMS='' 
  for X in `ls *covPerGe.gz`; do 
   N=\${X%%.covPerGe.gz}
   NAMES="\$NAMES \$N"
   BAMS="\$BAMS \${N}.bam"
   CHRM="\$CHRM chrCoverageMedians_\$N"
   COVPERGE="\$COVPERGE \$X";
  done
  #remove leading ithe space
  NAMES="\$(echo -e "\${NAMES}" | sed -e 's/^[[:space:]]*//')"
  BAMS="\$(echo -e "\${BAMS}" | sed -e 's/^[[:space:]]*//')"
  CHRM="\$(echo -e "\${CHRM}" | sed -e 's/^[[:space:]]*//')"
  COVPERGE="\$(echo -e "\${COVPERGE}" | sed -e 's/^[[:space:]]*//')"

  #find and quantify gene clstrs
  bash /bin/covPerClstr.sh -f "\$COVPERGE" -n "\$NAMES" -b "\$BAMS" -c "\$CHRM" -m $MAPQ -g $annotation -a $fa -l 0.9 -s 0.9 -o covPerClstr
  #annotate
  for CL in `ls covPerClstr/lowMapq.clstr/`; do 
   echo \$CL >> covPerClstr/clstrAnn.tsv
   for X in `grep ">" covPerClstr/lowMapq.clstr/\$CL | cut -f2 -d">"`; do grep -m1 -P "^\$X\\t" $geFun ; done >> covPerClstr/clstrAnn.tsv
   echo >> covPerClstr/clstrAnn.tsv 
   for X in `grep ">" covPerClstr/lowMapq.clstr/\$CL | cut -f2 -d">"`; do grep -m1 -P "^\$X\\t" $geFun ; done | awk '{printf("%s\\t%s\\n", \$0, "'\$CL'") }' >> covPerClstr/clstrAnnFormat2.tsv
  done
  
  """
}

process report {
  input:
  file ('*') from covPerChr4.join(covPerNt).join(covPerBin).join(mappingStats).join(covPerGeDump1).join(snpEff).join(delly).join(mapDump1).join(bigWigGenomeCov)

  output:
  file("report_*.html") into end  

  script:
  """
  sampleId=`cat input.1`
  reportFileName=report_\${sampleId}.html
  cp /bin/buildReport.Rmd .
  R -e "sample='\$sampleId'; version='$version'; rmarkdown::render('buildReport.Rmd' , output_file = '\$reportFileName')"
  """
}

