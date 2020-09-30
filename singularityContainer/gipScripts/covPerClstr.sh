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
/opt/cdhit/make_multi_seq.pl ${outDir}/.lowMapqSelected.fa ${outDir}/.lowMapqSelected.clstr ${outDir}/lowMapq.clstr 2
for X in `grep ">" ${outDir}/lowMapq.clstr/* | cut -f2 -d ">"`; do grep -m1 "gene_id \"$X\"" $geGtf  ; done > ${outDir}/lowMapq.clstr.ge.gtf
#rename clusters
readarray oldClNames < <(ls ${outDir}/lowMapq.clstr | sort -n)
for (( c=1; c<=${#oldClNames[@]}; c++ )); do 
  oldCl=`echo ${oldClNames[${c}-1]} | tr -d "\n"`
  mv ${outDir}/lowMapq.clstr/$oldCl ${outDir}/lowMapq.clstr/clstr${c}
done



function covPerGe_allReads {
########
#OUTPUT#
########
#.covPerGe file. This format includes both the gene mean coverage and the gene mean coverage normalized by the chromosome median coverage. To estimate the mean coverage the N bases are not considered
#This special version of covPerGe uses all mapping reads (i.e. no MAPQ filter) and does not estimate the average MAPQ to speed up
#NOTE: genes not on the chromosomes of the chrCoverageMedians file are removed from the output
#covPerGe fields:
#gene_id: gene identifier
#locus: gene chr:start-end
#meanCoverage: Mean nucleotide sequencing coverage of the bases belonging to the gene (N bases aren't counted)
#normalizedMeanCoverage: #meanCoverage (N bases aren't counted) / chromosome median coverage
#MAPQ: a symbolic value of 99 for all genes

##WARNING: You can have genes coverage, but still with a MAPQ >0 because you can have very few reads mapping in the gene (with a certain MAPQ). Similarly, you can have genes with 0 coverage and 0 MAPQ (so no mapping reads at all)

#######
#INPUT#
#######
#BAM: A genomic-sequencing BAM
#OUT: covPerGe output name
#GENES: A gene gtf annotation file
#CHRM: chromosome median coverage file (chrCoverageMedians file, syntax: chrXX median .* )
#BITFLAG: reads with this flag won't contribute to the sequencing depth nor to the gene average MAPQ
#ASSEMBLY: genome assembly multifasta file

local BAM=$1
local OUT=$2
local GENES=$3
local CHRM=$4
local BITFLAG=$5
local ASSEMBLY=$6

perl -e '
open(C,"<'$CHRM'");
while(<C>){
if($_=~/^(\S+)\s+(\S+)/){$chr=$1;$m=$2; $mCovs{$chr} = $m;}
}
close C;

open(F,"<'$GENES'");
open(O,">'$OUT'");
print O "gene_id\tlocus\tmeanCoverage\tnormalizedMeanCoverage\tMAPQ\n";
while(<F>){
 open(TMP,">'$OUT'_currentGene.bed");
 if($_=~/gene_id \"([^\"]+)\"/){$ge=$1;}
 if($_=~/^(\S+)\s+\S+\s+\S+\s+(\S+)\s+(\S+)/){$chr=$1;$st=$2;$en=$3;}
 if(! defined $mCovs{$chr}){close(TMP); next;}
 print O "$ge\t${chr}:${st}-${en}\t";
 print TMP "${chr}\t${st}\t${en}\n";

 #subtract from the gene length the number of Ns
 $cmd1="samtools faidx '$ASSEMBLY' ${chr}:${st}-${en}  | tail -n +2  | awk \x27{print toupper(\$0)}\x27 | tr -cd  N | wc -c | tr -d \x27\\n\x27";
 $Ns=`$cmd1`; 
 die ("covPerGe faidx did not work for command\n$cmd1\n$!") if ($?);
 my $realLength = ($en - $st) - $Ns;

 #mean gene coverage
 $cmd2="samtools view -b -F " . "'$BITFLAG'" . " " . "'$BAM'" . " ${chr}:${st}-${en} | bedtools coverage -a " . "'$OUT'" . "_currentGene.bed -b stdin -d -split | awk  \x27{x+=\$5;next}END{print x/$realLength}\x27"  ;
 $meanGeneCov=`$cmd2`; if($?){die "error with $cmds2\n";}
 chomp $meanGeneCov;

 #normalize
 $normalizedMeanGeneCov = $meanGeneCov / $mCovs{$chr};
 print O "$meanGeneCov\t$normalizedMeanGeneCov\t99\n";
 close(TMP);
}
close O;
close F; '
rm -rf ${OUT}_currentGene.bed
gzip $OUT
}


for i in "${!covPerGeNames[@]}"; do 
  N=${covPerGeNames[$i]}
  BAM=${bamFiles[$i]}
  CHRM=${chrCoverageMediansFiles[$i]}
  OUT="${outDir}/${N}.covPerLowMapqClstrGe_allReads"
  echo $N
  #covPerGe_allReads for all samples
  covPerGe_allReads $BAM $OUT ${outDir}/lowMapq.clstr.ge.gtf $CHRM 1028 $genome
  #average cluster coverage
  echo -e "gene_id\tlocus\tmeanCoverage\tnormalizedMeanCoverage\tMAPQ" > ${outDir}/${N}.covPerClstr
  for clstr in `ls ${outDir}/lowMapq.clstr`; do
    perl -e '
    #read clstr ids
    open(F,"<'${outDir}'/lowMapq.clstr/'$clstr'"); while(<F>){if($_=~/^>(.+)/){$clIds{$1}=1;}} close F; 
    #read covPerLowMapqClstrGe_allReads file
    open(IN, "gunzip -c '$OUT'.gz |") or die "gunzip '$OUT': $!";
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
  #clean
  rm -rf ${OUT}.gz
done

#clean
rm -rf ${outDir}/lowMapq.clstr.ge.gtf ${outDir}/.lowMapq  ${outDir}/.lowMapqCount  ${outDir}/.lowMapqSelected  ${outDir}/.lowMapqSelected.clstr  ${outDir}/.lowMapqSelected.fa  ${outDir}/.tmpGeBed





