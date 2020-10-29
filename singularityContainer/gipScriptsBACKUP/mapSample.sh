#!/bin/bash
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

