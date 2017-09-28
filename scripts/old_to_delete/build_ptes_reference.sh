#!/bin/sh

module load apps/bowtie2


CONSTRUCTS="$1/tempConstructs.fasta"
CAN="$1/tempCan.fasta"
LOG="$1/bowtie-build.log"
ONE="$1/Can.fa"
TWO="$1/Constructs.fa"

TEMPONE="$1/CannonicalJoints.fasta"
TEMPTWO="$1/PTESJoints.fasta"

INDEXONE="$1/canonical"
INDEXTWO="$1/ptes"

bowtie2-build $ONE $INDEXONE 2>> $LOG
bowtie2-build $TWO $INDEXTWO 2>> $LOG
