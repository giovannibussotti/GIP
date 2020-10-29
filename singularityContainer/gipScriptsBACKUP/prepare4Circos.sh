#!/bin/bash
BINSIZE=25000

while getopts ':i:m:s:g:c:b:x:t:z:j:k:w:h' OPTION ; do
case $OPTION in
  i)  SAMPLEID=$OPTARG;;
  m)  BAM=$OPTARG;;
  s)  CHRSIZE=$OPTARG;;
  g)  GEGTF=$OPTARG;;
  c)  CHRtoUSE=$OPTARG;;
  b)  BINSIZE=$OPTARG;;
  x)  covPerChr=$OPTARG;;
  t)  CONFTEMPLATE=$OPTARG;;
  z)  TRADATA=$OPTARG;;
  j)  INVDATA=$OPTARG;;
  k)  DUPDATA=$OPTARG;;
  w)  DELDATA=$OPTARG;;
	h)	echo "USAGE of this program:"
      echo "   -i    sample id"
      echo "   -m    bam file. The same folder must contain the bai file"
      echo "   -s    chr size file"
      echo "   -g    gene gtf file"
      echo "   -c    space separated list of chromosomes to use in quotation marks"
      echo "   -b    bin size"
      echo "   -x    coverage per chromosome file"
      echo "   -t    configuration template"
      echo "   -z    traslocation data"
      echo "   -j    inversion data"
      echo "   -k    duplication data"
      echo "   -w    deletion data"
      echo "example: ./prepare4Circos.sh -x chrCoverageMedians_Ldo_CH33 -i Ldo_CH33 -m Ldo_CH33_EP.bam -s /pasteur/projets/policy01/BioIT/Giovanni/datasets/projects/p2p5/LdBPKv2.chrSize -g /pasteur/projets/policy01/BioIT/Giovanni/datasets/projects/p2p5/LdBPKv2.ge.gtf -c \"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36\" -b 25000 -t templateCircos.conf -z Ldo_CH33_EP+3.delly.TRA.filter.circosBed -j Ldo_CH33_EP+3.delly.INV.filter.circosBed -k Ldo_CH33_EP+3.delly.DUP.filter.circosBed -w Ldo_CH33_EP+3.delly.DEL.filter.circosBed"
	    exit 0;;
	\?)	echo "Unknown argument \"-$OPTARG\"."
	echo $HELP ; exit 1;;
	:)	echo "Option \"-$OPTARG\" needs an argument"
	echo $HELP ; exit 1;;
	*)	echo "Sorry, I made a mistake when programming this "\"$OPTION\";;
esac
done

######################
#chrSize to Karyotype#
######################
D=${SAMPLEID}_circosData
mkdir -p $D
perl -e '
  #hash with chr to use
  $chrString="'"${CHRtoUSE}"'"; @chrs=split(/\s+/,$chrString); foreach $c(@chrs){$chrToUse{$c}=1;}
  #print chr size to karyotype file
  open(F,"<'$CHRSIZE'"); while(<F>){if($_=~/(\S+)\s+(\S+)/){$chr=$1; $len=$2; if($chrToUse{$chr}){print "chr - chr$chr $chr 0 $len black\n";}} }   close F;
' > $D/karyotype.txt
awk '{print $3 "\t" 1 "\t" $6 }' $D/karyotype.txt > $D/chrSize.bed

#######
#genes#
#######
perl -e '
  #hash with chr to use
  $chrString="'"${CHRtoUSE}"'"; @chrs=split(/\s+/,$chrString); foreach $c(@chrs){$chrToUse{$c}=1;}
  #print genes on + and - strands
  open(P,">'$D'/plus.ge.bed");
  open(M,">'$D'/minus.ge.bed");
  open(F,"<'$GEGTF'"); while(<F>){ 
  	@fi=split(/\t/,$_); 
  	if(! $chrToUse{$fi[0]}){next;} 
    $out="chr$fi[0]\t$fi[3]\t$fi[4]\n"; 
    if($fi[6] eq "+"){
    	print P $out;
    } elsif ($fi[6] eq "-"){
    	print M $out;
    } else {print "$fi[6]\n";}
  }'


#######################################################
#bin the genome coverage for circos to be able to read#
#######################################################
#clean
rm -rf ${SAMPLEID}.covPerBin.gz
#covPerBin
covPerBin $BAM $BINSIZE $CHRSIZE $covPerChr _tmp
#select chrs to use
perl -e '
  #hash with chr to use
  $chrString="'"${CHRtoUSE}"'"; @chrs=split(/\s+/,$chrString); foreach $c(@chrs){$chrToUse{$c}=1;}
  open(F, "gunzip -c '${SAMPLEID}'.covPerBin.gz |") or die "cannot open gzip file $!\n"; 
  while(<F>){ 
   if($_=~/(\S+)/){
    if ($chrToUse{$1}){print;}
   }
  }' > _${SAMPLEID}.covPerBin ; 
mv _${SAMPLEID}.covPerBin ${SAMPLEID}.covPerBin
tail -n +2 ${SAMPLEID}.covPerBin | awk '{print "chr"$1 "\t" $2 "\t" $3 "\t" $5}' > $D/sample.covPerBin
tail -n +2 ${SAMPLEID}.covPerBin | awk '{print "chr"$1 "\t" $2 "\t" $3 "\t" $6}' > $D/sample.mapqPerBin
rm -rf ${SAMPLEID}.covPerBin


###########################
#PREPARE CIRCOS .CONF FILE#
###########################
cp $CONFTEMPLATE tmp.conf
chrList=`cut -f1 ${SAMPLEID}_circosData/chrSize.bed | tr "\n" ";"`
chrList=${chrList%%;}
sed -i 's/__DATADIR__/'$D'/' tmp.conf
sed -i 's/__CHRS__/'$chrList'/' tmp.conf
sed -i 's/__TRADATA__/'$TRADATA'/' tmp.conf
sed -i 's/__INVDATA__/'$INVDATA'/' tmp.conf
sed -i 's/__DUPDATA__/'$DUPDATA'/' tmp.conf
sed -i 's/__DELDATA__/'$DELDATA'/' tmp.conf


#run circos
perl /opt/circos/bin/circos -con tmp.conf 
rm -rf circos.svg
mv circos.png ${SAMPLEID}.SV.circos.png


