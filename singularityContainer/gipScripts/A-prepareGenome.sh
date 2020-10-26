#!/bin/bash
ENGINE=rmblast
CPU=1
OUT=genome
LIB=default
while getopts ':f:x:l:c:o:e:h' OPTION ; do
case $OPTION in
        f)      FA=$OPTARG;;
        x)      EX=$OPTARG;;
        l)      LIB=$OPTARG;;
        c)      CPU=$OPTARG;;
        o)      OUT=$OPTARG;;
        e)      ENGINE=$OPTARG;;
        h)      echo "USAGE of this program:"
        echo "   -f    genome multi-fasta file"
        echo "   -x    gene gtf file"
        echo "   -l    repeatMasker library"
        echo "   -c    cpus"
        echo "   -o    out prefix"
        echo "   -e    repeatMasker engine"
        exit 0;;
        \?)     echo "Unknown argument \"-$OPTARG\"."
        echo $HELP ; exit 1;;
        :)      echo "Option \"-$OPTARG\" needs an argument"
        echo $HELP ; exit 1;;
        *)      echo "Sorry, I made a mistake when programming this "\"$OPTION\";;
esac
done

##################
#PREPARE ASSEMBLY#
##################
#clean chr ids#
#zipped input fasta are also accepted#
referencedFA=`readlink -f $FA`
if file $referencedFA | grep -q compressed  ; then
    gunzip -c $FA | perl -ne 'if($_=~/^>(\S+)/){$c=$1; print ">$c\n";}else{print uc($_);} ' > ${OUT}   
else
    cat $FA | perl -ne 'if($_=~/^>(\S+)/){$c=$1;  print ">$c\n";}else{print uc($_);} ' > ${OUT}
fi

#REPEATS
N=`basename $OUT`
D=`dirname $OUT`
mkdir -p $D/repeats
if [ $LIB == "default" ]; then 
  #Red
  mkdir -p _tmp
  mv $OUT _tmp/${OUT}.fa
  Red -gnm _tmp -msk repeats -rpt repeats -frm 2
  perl -e '
  open(F,"<repeats/'$N'.rpt") or die "cannot open Red output"; 
  $i=1; 
  while(<F>){
    $_=~s/^>//; 
    if ($_ =~/(\S+)\s+(\d+)\s+(\d+)/){
      $chr   = $1;
      $start = $2 + 1;
      $end   = $3;
      print "$chr\tRed\tML\t$start\t$end\t\.\t\*\t\.\t\.\trepeat_$i\n";
      $i++;
    }  
  } ' > repeats/${N}.out.gff
  mv repeats/${N}.msk ${N}.fa
  rm -rf _tmp/ repeats/${N}.rpt

  else
  #RepeatMasker
  /opt/RepeatMasker/RepeatMasker -e $ENGINE -pa $CPU -gff -xsmall -lib $LIB -dir repeats -s ${OUT}
  gzip repeats/${N}.out
  gzip repeats/${N}.ori.out
  mv repeats/${N}.masked ${OUT}
  mv ${OUT} ${OUT}.fa

fi

#chrSize#
awk '/^>/ {if (seqlen) print seqlen;printf ($1 "\t") ; seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' ${OUT}.fa | sed -e 's/^>//' | sort -nk1 > ${OUT}.chrSize
#Indexes#
samtools faidx ${OUT}.fa
#BWA DB#
mkdir -p db/bwa/$N
bwa index -p db/bwa/${N}/db ${OUT}.fa
#.DICT#
java -jar /bin/picard.jar CreateSequenceDictionary R=${OUT}.fa O=${OUT}.dict
#GAP COORDINATES#
perl -ne 'chomp;if( />(.*)/){ if($s==1){print "\t$i\n";$s=0;} $head = $1; $i=0; next};@a=split("",uc($_)); foreach(@a){$i++;if($_ eq "N" && $s ==0 ){print "$head\t$i"; $s =1}elsif($s==1 && $_ ne "N"){print "\t$i\n";$s=0}}' ${OUT}.fa | gzip > ${OUT}.gaps.gz
########
#snpEff#
########
mkdir -p snpEff/cache
mkdir -p snpEff/data/$N
FAp=`readlink -e ${OUT}.fa`
cat $EX | sed -e 's/gene/exon/' > snpEff/data/${N}/genes.gtf
ln -fs $FAp snpEff/data/${N}/sequences.fa
cp /opt/snpEff.configTemplate snpEff/snpEff.config
sed -i 's/GENOME/'$N'/g' snpEff/snpEff.config
java -jar /opt/snpEff/snpEff.jar build -gtf22 -v $N -noLog -c snpEff/snpEff.config

