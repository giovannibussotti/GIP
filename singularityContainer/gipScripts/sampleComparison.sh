#!/bin/bash
ALLIN=$*
IFS=' ' read -r -a IN <<< "$ALLIN"
MODE=${IN[0]}
OPTS=${ALLIN/$MODE/}
#GIPOUT=${IN[1]}
#OPTS=${ALLIN/$MODE[[:space:]]*$GIPOUT/}


#function findArgs {
#find the content of a parameter
#returns the ARGfound array where ${ARGfound[0]} contains the parameter, and ${ARGfound[1]} the rest od the command line
#  local ARG=$1
#  local OPTS=$2
#  local res=`echo $OPTS | perl -ne '$arg='$ARG'; chomp $_; if($_=~/--$arg\s+([^-]+)/){print "$1";}else{die "ERROR sampleComparison findArgs. $arg not found";}'`
#  local rest=${OPTS/--$ARG[[:space:]]*$res}
#  ARGfound=("$res" "$rest")
#}
#findArgs samples "$OPTS"
#SAMPLES=${ARGfound[0]}

if [ $MODE == "karyotype" ]; then
    CMD="Rscript /bin/karyotype $OPTS"
    echo executing $CMD
    $CMD

elif [ $MODE == "binCNV" ]; then
    CMD="Rscript /bin/binCNV $OPTS"
    echo executing $CMD
    $CMD

elif [ $MODE == "geCNV" ]; then
    CMD="Rscript /bin/geCNV $OPTS"
    echo executing $CMD
    $CMD

elif [ $MODE == "ternary" ]; then
    CMD="Rscript /bin/ternary $OPTS"
    echo executing $CMD
    $CMD

elif [ $MODE == "ternaryBin" ]; then
    CMD="Rscript /bin/ternaryBin $OPTS"
    echo executing $CMD
    $CMD

 elif [ $MODE == "SNV" ]; then
    CMD="Rscript /bin/SNV $OPTS"
    echo executing $CMD
    $CMD

 elif [ $MODE == "binDensity" ]; then
    CMD="Rscript /bin/binDensity $OPTS"
    echo executing $CMD
    $CMD

 elif [ $MODE == "geInteraction" ]; then
    CMD="Rscript /bin/geInteraction $OPTS"
    echo executing $CMD
    $CMD

 elif [ $MODE == "genomeDistance" ]; then
    CMD="Rscript /bin/genomeDistance $OPTS"
    echo executing $CMD
    $CMD

else 
  echo "Error. $MODE sample comparison mode not recognized"
fi

