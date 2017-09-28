#!/bin/sh

module load apps/bowtie2

INDEXONE="$1/canonical"
INDEXTWO="$1/ptes"
READS="$2"
OUTONE="$1/canonical.sam"
OUTTWO="$1/ptes.sam"
LOG="$1/bowtie.log"

bowtie2 -p8 --very-sensitive --no-hd --no-unal -x $INDEXONE $READS -S $OUTONE 2>>$LOG 
bowtie2 -p8 --very-sensitive --no-hd --no-unal -x $INDEXTWO $READS -S $OUTTWO 2>>$LOG

