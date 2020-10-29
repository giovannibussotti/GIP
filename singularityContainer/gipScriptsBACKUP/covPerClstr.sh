#!/bin/bash
#given a set of covPerGe files it 1) selects genes that are low mapq across all covPerGes 2) cluster them 3) select just the genes belonging to clusters 4) evaluate their coverage without MAPQ filters 5) report average coverage for each cluster

#NOTE: you do not have to remove clusters with members on different chromosomes because differences in ploidy are taken into account by covPerGe_allReads normalizing each gene by chrMedianCov

#dev
#covPerGeFilesString="CL3_7145.covPerGe.gz CL8_7149.covPerGe.gz"
#covPerGeNamesString="CL3 CL8"
#chrCoverageMediansFilesString="chrCoverageMedians_CL3_7145 chrCoverageMedians_CL8_7149"
#bamFilesString="../../../pipeOut/ld1sSet/lsdOut/CL3_7145.bam ../../../pipeOut/ld1sSet/lsdOut/CL8_7149.bam"
#minMapq=10
#geGtf=Ld1S.ge.gtf
#genome=Leishmania_donovani_sudanese.fa
#cdHitLenDiffCutoff=0.9
#cdHitSeqIdCutoff=0.9

#ex
#bash covPerFam.sh -f "CL3_7145.covPerGe.gz CL8_7149.covPerGe.gz" -n "CL3 CL8" -b "../../../pipeOut/ld1sSet/lsdOut/CL3_7145.bam ../../../pipeOut/ld1sSet/lsdOut/CL8_7149.bam" -c "chrCoverageMedians_CL3_7145 chrCoverageMedians_CL8_7149" -m 10 -g Ld1S.ge.gtf -a Leishmania_donovani_sudanese.fa -l 0.9 -s 0.9

#default
minMapq=50
outDir=outDirCovPerClstr

while getopts ':f:n:c:b:m:g:a:l:s:o:h' OPTION ; do
case $OPTION in
  f)  covPerGeFilesString=$OPTARG;;
  n)  covPerGeNamesString=$OPTARG;;
  c)  chrCoverageMediansFilesString=$OPTARG;;
  b)  bamFilesString=$OPTARG;;
  m)  minMapq=$OPTARG;;
  g)  geGtf=$OPTARG;;
  a)  genome=$OPTARG;;
  l)  cdHitLenDiffCutoff=$OPTARG;;
  s)  cdHitSeqIdCutoff=$OPTARG;;
  o)  outDir=$OPTARG;;
  h)  echo "USAGE of this program:"
  echo "   -f    space sparated covPerGe.gz files"
  echo "   -n    space sparated sample names"
  echo "   -c    space sparated chrCoverageMedians files relative to the input covPerGe"
  echo "   -b    space sparated bam files relative to the input covPerGe"
  echo "   -m    minMapq (default 50)"
  echo "   -g    gene gtf file"
  echo "   -a    genome assmbly multi-fasta file"
  echo "   -l    cd-hit length difference cutoff. if set to 0.9, the shorter sequences need to be at least 90% length of the representative of the cluster"
  echo "   -s    cd-hit sequence identity threshold, this is the global sequence identity calculated as: number of identical amino acids in alignment divided by the full length of the shorter sequence"
  echo "   -o    outDir"
  exit 0;;
  \?) echo "Unknown argument \"-$OPTARG\"."
  echo $HELP ; exit 1;;
  :)  echo "Option \"-$OPTARG\" needs an argument"
  echo $HELP ; exit 1;;
  *)  echo "Sorry, I made a mistake when programming this "\"$OPTION\";;
esac
done

#clean
mkdir -p ${outDir}
rm -rf ${outDir}/lowMapq.clstr/ ${outDir}/.lowMapq ${outDir}/.lowMapqCount ${outDir}/.lowMapqSelected ${outDir}/.lowMapqSelected.clstr  ${outDir}/.lowMapqSelected.fa  ${outDir}/.tmpGeBed

#string to array
IFS=' ' read -r -a covPerGeFiles <<< "$covPerGeFilesString"
IFS=' ' read -r -a covPerGeNames <<< "$covPerGeNamesString"
IFS=' ' read -r -a chrCoverageMediansFiles <<< "$chrCoverageMediansFilesString"
IFS=' ' read -r -a bamFiles <<< "$bamFilesString"

#select low mapq
TOTfiles=0
for F in "${covPerGeFiles[@]}"; do
  gunzip -c $F | tail -n +2 | awk '{if($5 < '$minMapq') print $1}' >> ${outDir}/.lowMapq
  TOTfiles=$((TOTfiles+1))
done

#count lowMAPQ genes
cat ${outDir}/.lowMapq | sort | uniq -c > ${outDir}/.lowMapqCount

#remove genes that are not low MAPQ in all samples
awk '{if($1 == '$TOTfiles') print $2}' ${outDir}/.lowMapqCount > ${outDir}/.lowMapqSelected

#extract stranded gene sequences
for X in `cat ${outDir}/.lowMapqSelected`; do
  grep "gene_id \"$X\"" $geGtf
done |  perl -ne 'if($_=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)\s+\S+\s+(\S+)/){$c=$1;$s=$2;$e=$3;$st=$4;if($_=~/gene_id \"([^\"]+)\"/){print "$c\t$s\t$e\t$1\t0\t$st\n";}} ' > ${outDir}/.tmpGeBed
bedtools getfasta -name -s -fi $genome -bed ${outDir}/.tmpGeBed -fo ${outDir}/.lowMapqSelected.fa

#clustering
/opt/cdhit/cd-hit-est -T 1 -s $cdHitLenDiffCutoff -c $cdHitSeqIdCutoff -r 0 -d 0 -g 1 -i ${outDir}/.lowMapqSelected.fa -o ${outDir}/.lowMapqSelected #-uL 0.05 -uS 0.05
/opt/cdhit/make_multi_seq.pl ${outDir}/.lowMapqSelected.fa ${outDir}/.lowMapqSelected.clstr ${outDir}/lowMapq.clstr 1
for X in `grep ">" ${outDir}/lowMapq.clstr/* | cut -f2 -d ">"`; do grep -m1 "gene_id \"$X\"" $geGtf  ; done > ${outDir}/lowMapq.clstr.ge.gtf
#rename clusters
readarray oldClNames < <(ls ${outDir}/lowMapq.clstr | sort -n)
for (( c=1; c<=${#oldClNames[@]}; c++ )); do 
  oldCl=`echo ${oldClNames[${c}-1]} | tr -d "\n"`
  mv ${outDir}/lowMapq.clstr/$oldCl ${outDir}/lowMapq.clstr/clstr${c}
done


#compute mean clstrs coverage
for i in "${!covPerGeNames[@]}"; do 
  N=${covPerGeNames[$i]}
  BAM=${bamFiles[$i]}
  CHRM=${chrCoverageMediansFiles[$i]}
  COVPERGE=${covPerGeFiles[$i]}
  echo $N

  #average cluster coverage
  echo -e "gene_id\tlocus\tmeanCoverage\tnormalizedMeanCoverage\tMAPQ" > ${outDir}/${N}.covPerClstr
  for clstr in `ls ${outDir}/lowMapq.clstr`; do
    perl -e '
    #read clstr ids
    open(F,"<'${outDir}'/lowMapq.clstr/'$clstr'"); while(<F>){if($_=~/^>(.+)/){$clIds{$1}=1;}} close F; 
    #read covPerLowMapqClstrGe_allReads file
    open(IN, "gunzip -c '$COVPERGE' |") or die "gunzip '$COVPERGE': $!";
    while(<IN>){
      if($_=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){ $id=$1; $pos=$2; $cov=$3; $nCov=$4;
        if($clIds{$id}){ $clPos=$pos; $clCov += $cov;$clNCov += $nCov; }
      }
    }
    $clCov  = $clCov / scalar(keys %clIds); 
    $clNCov = $clNCov / scalar(keys %clIds); 
    print "'$clstr'\t$clPos\t$clCov\t$clNCov\t99\n";
    close IN;' >> ${outDir}/${N}.covPerClstr
  done
  gzip ${outDir}/${N}.covPerClstr
done


#clean
rm -rf ${outDir}/lowMapq.clstr.ge.gtf ${outDir}/.lowMapq  ${outDir}/.lowMapqCount  ${outDir}/.lowMapqSelected  ${outDir}/.lowMapqSelected.clstr  ${outDir}/.lowMapqSelected.fa  ${outDir}/.tmpGeBed

