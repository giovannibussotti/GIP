#!/bin/bash
#############################################################################
# giptools                                                                  #
#                                                                           #
# Authors: Giovanni Bussotti                                                #
# Copyright (c) 2021  Institut Pasteur                                      #
#                                                                           #
#                                                                           #
# This program is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as                   #
# published by the Free Software Foundation, either version 3 of the        #
# License, or (at your option) any later version.                           #
#                                                                           #
# This program is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of            #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
# GNU General Public License for more details.                              #
#                                                                           #
# You should have received a copy of the GNU General Public License         #
# along with this program.  If not, see <https://www.gnu.org/licenses/>.    #
#                                                                           #
#############################################################################

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

function printHelp {
  echo ""
  echo "gitools modules:"
  echo ""
  echo -e " karyotype       Compare the chromosome sequencing coverage distributions"
  echo -e " binCNV          Compare bin sequencing coverage in 2 samples"
  echo -e " geCNV           Compare gene sequencing coverage in 2 samples"
  echo -e " ternary         Compare gene sequencing coverage in 3 samples"
  echo -e " ternaryBin      Compare bin sequencing coverage in 3 samples"
  echo -e " SNV             Compare SNVs in multiple samples"
  echo -e " binDensity      Density plot of bin sequencing coverage of many samples"
  echo -e " geInteraction   Detect CNV genes in many samples and produce correlation-based networks"
  echo -e " genomeDistance  Compare samples genomic distance"
  echo -e " phylogeny       Extract the SNVs union and infer the phylogenetic tree"
  echo -e " convergentCNV   Detect convergent CNV gene amplifications"
  echo -e " overview        Overview of the sequencing coverage of chromosomes, genomic bins and genes" 
  echo ""
}

if [ -z $MODE ]; then
  printHelp
  exit 0
fi

if [ $MODE == "karyotype" ]; then
    CMD="Rscript /bin/karyotype $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "binCNV" ]; then
    CMD="Rscript /bin/binCNV $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "geCNV" ]; then
    CMD="Rscript /bin/geCNV $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "ternary" ]; then
    CMD="Rscript /bin/ternary $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "ternaryBin" ]; then
    CMD="Rscript /bin/ternaryBin $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "SNV" ]; then
    CMD="Rscript /bin/SNV $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "binDensity" ]; then
    CMD="Rscript /bin/binDensity $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "geInteraction" ]; then
    CMD="Rscript /bin/geInteraction $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "genomeDistance" ]; then
    CMD="Rscript /bin/genomeDistance $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "phylogeny" ]; then
    OPTS2=`echo $OPTS | perl -ne ' chomp $_; $_=~s/--iqtreeOpts\s+\"([^\"]*)\"//; print "$_";'`
    OPTS3=`echo $OPTS | perl -ne ' chomp $_; if($_=~/--iqtreeOpts\s+\"([^\"]*)\"/){print "$1";} '`
    if [ -z "$OPTS3" ]; then
      Rscript /bin/phylogeny $OPTS2
    else
      Rscript /bin/phylogeny $OPTS2 --iqtreeOpts "$OPTS3"
    fi

elif [ $MODE == "convergentCNV" ]; then
    CMD="Rscript /bin/convergentCNV $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "overview" ]; then
    CMD="Rscript /bin/overview $OPTS"
    #echo executing $CMD
    $CMD

elif [ $MODE == "-h" ] || [ $MODE == "--help" ]; then
    printHelp
 
else 
  echo "Error. $MODE sample comparison mode not recognized"
fi

