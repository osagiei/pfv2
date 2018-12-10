#!/bin/sh

READS=$1
WORKING_DIR=$2

LOG="${WORKING_DIR}/bowtie-build.log"
ONE="${WORKING_DIR}/Can.fa"
TWO="${WORKING_DIR}/Constructs.fa"

INDEXONE="${WORKING_DIR}/canonical"
INDEXTWO="${WORKING_DIR}/ptes"
OUTONE="${WORKING_DIR}/canonical.sam"
OUTTWO="${WORKING_DIR}/ptes.sam"

bowtie2-build $ONE $INDEXONE 2>> $LOG
bowtie2-build $TWO $INDEXTWO 2>> $LOG

LOG="${WORKING_DIR}/bowtie.log"

bowtie2 -p8 --very-sensitive --no-hd --no-unal -x $INDEXTWO $READS -S $OUTTWO 2>>$LOG
bowtie2 -p8 --very-sensitive --no-hd --no-unal -x $INDEXONE $READS -S $OUTONE 2>>$LOG


#rm $WORKING_DIR/*.bt2 
