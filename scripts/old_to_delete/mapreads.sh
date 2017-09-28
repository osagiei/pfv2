#!/bin/sh

module load apps/bowtie

REF=$2
UNIQUE=$3
LEFTREADS="$1/left.fastq"
RIGHTREADS="$1/right.fastq"

TEMPLEFT="$1/tempL.sam"
TEMPRIGHT="$1/tempR.sam"
OUTLEFT="$1/left.sam"
OUTRIGHT="$1/right.sam"
UNLEFT="$1/unleft.sam"
UNRIGHT="$1/unright.sam"
RELEFT="$1/releft.sam"
RERIGHT="$1/reright.sam"
LOG="$1/bowtie.log"

bowtie -v 1 -m $UNIQUE -q -k 1 --best --sam --sam-nohead $REF $LEFTREADS $TEMPLEFT --un $UNLEFT 2>>$LOG
bowtie -v 1 -m $UNIQUE -q -k 1 --best --sam --sam-nohead $REF $RIGHTREADS $TEMPRIGHT --un $UNRIGHT 2>>$LOG


bowtie -v 1 -m 1 --trim3 5 -q --sam --sam-nohead $REF $UNLEFT $RELEFT 2>>$LOG
bowtie -v 1 -m 1 --trim5 5 -q --sam --sam-nohead $REF $UNRIGHT $RERIGHT 2>>$LOG

cat $TEMPLEFT $RELEFT > $OUTLEFT
cat $TEMPRIGHT $RERIGHT > $OUTRIGHT

rm $TEMPLEFT $TEMPRIGHT $RELEFT $RERIGHT
