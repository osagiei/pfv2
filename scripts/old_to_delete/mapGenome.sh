#!/bin/sh

module load apps/bowtie2


REF="$3"

INREADS="$1"

OUT="$2/genomic.sam"

LOG="$2/bowtie.log"

bowtie2 -p 8 --very-sensitive --no-hd -x $REF $INREADS -S $OUT 2>> $LOG
