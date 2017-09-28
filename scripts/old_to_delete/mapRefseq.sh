#!/bin/sh

module load apps/bowtie2

INREADS="$1"

OUT="$2/refseq.sam"
LOG="$2/bowtie.log"

REF="$3"

bowtie2 -p 8 --very-sensitive --no-hd -x $REF $INREADS -S $OUT 2>> $LOG



