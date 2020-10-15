#the idea is to speed up covPerBin and covPerGe steps
#the covPerBin generated with covPerBin2 does not have the median field so you need to update downstream script to accept a different input, that is:
#1) the script that does GC correction
#2) the script that plots covPerBin coverage
#3) the script that does kayoplote that also uses covPerBin as input

#the sample.bed file generated at this step should be passed to the covPerGe process and reused there to generate covPerGe file

function covPerBin_normalizedByChrMedianCov2 {
########
#OUTPUT#
########
#.covPerBin file normalized by chromosome median coverage    
#######
#INPUT#
#######
#A chromosome median coverage file (chrCoverageMedians file, syntax: chrXX median .* ) and a covPerBin file (syntax: chr start end meanCoverage medianCoverage .*)
#NOTE1: the input covPerBin file need to be not normalized before this function, so you need to generate it with covPerNt NOTNORMALIZED (followed by gencov2intervals.pl for the binning)
#NOTE2: covPerBin input windows that are not on the chromosomes in the chrCoverageMedians file are removed from the output
    local CHRM=$1
    local GCOVBIN=$2
    local OUT=$3
perl -e '
        open(C,"<'$CHRM'");
        while(<C>){
                if($_=~/^(\S+)\s+(\S+)/){$chr=$1;$m=$2; $mCovs{$chr} = $m;}
        }
        close C;

        open(O,">'$OUT'");
        open(F,"<'$GCOVBIN'");
        $header = <F>;
        print O "$header";
        while(<F>){
                if ($_=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(.*)$/){
                        $chr=$1;
                        $start=$2;
                        $end=$3;
                        $mean=$4;
                        $rest=$5;
                        next if(! defined $mCovs{$chr});
                        $normMean   = $mean / $mCovs{$chr};
                        print O "$chr\t$start\t$end\t${normMean}${rest}\n";      
                }
        }
        close F;
        close O;
'
}


function covPerBin2 {
  local BAM=$1
  local OUT=$2
  local BINSIZE=$3
  local CHRSIZE=$4
  local chrCoverageMedians=$5
  local TMP=$6

  mkdir -p $TMP
  N=${BAM/%.bam}
  BED=${BAM/%.bam/.bed} 
  bedtools bamtobed -split -i $BAM | sort -k1,1 -k2,2n > $BED
  bedtools makewindows -g $CHRSIZE -w $BINSIZE | sort -k1,1 -k2,2n > $TMP/windows
  bedtools map -a $TMP/windows -b $BED -c 5 -o mean -null 0 > $TMP/meanMapq
  bedtools coverage -d -a $TMP/windows -b $BED > $TMP/winCov
  echo -e "chromosome\tstart\tend\tmean" > $TMP/covPerBinNotNorm
  awk '{i=$1"\t"1+$2"\t"$3; count[i]++; sum[i]+=$5;} END{ for(i in count) {m = sum[i]/count[i]; print i, m}}  ' $TMP/winCov  | sort -k1,1 -k2,2n >> $TMP/covPerBinNotNorm
  
  covPerBin_normalizedByChrMedianCov2 $chrCoverageMedians $TMP/covPerBinNotNorm $TMP/covPerBin
  perl -e '
  #read mapq
  open(F,"<'$TMP/meanMapq'") or die "covPerBin2 cannot open meanMapq";
  my %h;
  while(<F>){ 
    if ($_=~/^(\S+)\s+(\S+)\s+\S+\s+(\S+)/){
        $h{"${1}_$2"} = $3;
    }
  } 
  close F;
  
  open(O,">'$OUT'") or die "covPerBin2 cannot open out";
  open(F,"<'$TMP/covPerBin'") or die "covPerBin2 cannot open covPerBin";
  $header = <F>;
  chomp $header;
  print O $header . "\tMAPQ\n";
  while(<F>){ 
    if ($_=~/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
        $chr   = $1;
        $start = $2;
        $end   = $3;
        $mean  = $4;
        $s     = $start -1;
        $mapq = $h{"${chr}_$s"};
        print O "$chr\t$start\t$end\t$mean\t$mapq\n";

    }
  } 
  close F;
  close O;
  '
  gzip $OUT
  rm -rf $TMP 
}


#usage example
BAM=Dog3.bam
OUT=testDog3.covPerBin
BINSIZE=500
CHRSIZE=../../genome/genome.chrSize
chrCoverageMedians=chrCoverageMedians_Dog3
TMP=_tmpDog3
covPerBin2 $BAM $OUT $BINSIZE $CHRSIZE $chrCoverageMedians $TMP




