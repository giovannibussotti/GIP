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

S=$1
C=$2
ASSEMBLY=$3
INDEX=$4
R1=$5
R2=$6
OUTDIR=$7
delDup=$8

if [[ ! $R1 == *","* ]]; then
	bwaMapSample $S $ASSEMBLY $C $INDEX $OUTDIR $delDup $R1 $R2
	exit 0
fi

#MERGE fastqs
TMP=${OUTDIR}/tmpMerge_$S
rm -rf $TMP
mkdir -p $TMP
IFS=',' read -r -a firstPairs <<< "$R1"
for F in "${firstPairs[@]}"; do
	echo "merging firstPair $F with the other firstPairs"
	zcat ${F} >> ${TMP}/${S}_1.fastq
done
IFS=',' read -r -a secondPairs <<< "$R2"
for F in "${secondPairs[@]}"; do
	echo "merging secondPair $F with the other secondPairs"
	zcat ${F} >> ${TMP}/${S}_2.fastq
done


#MAP
echo "bwaMapSample $S $ASSEMBLY $C $INDEX $OUTDIR $TMP/${S}_1.fastq $TMP/${S}_2.fastq"
bwaMapSample $S $ASSEMBLY $C $INDEX $OUTDIR $delDup $TMP/${S}_1.fastq $TMP/${S}_2.fastq

