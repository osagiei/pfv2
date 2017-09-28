#!/bin/sh

module load apps/bedtools
module load apps/bowtie
module load apps/bowtie2

inputBed="$3"
exonsBed="$2/exons.bed"
tempExonsBed="$exonsBed.temp"
tempExonsFASTA="$2/temp_exons.fa"
exonsFASTA="$2/exons.fa"
coords="$3.coords"
transcriptFASTA="$exonsFASTA.mrna";
mapper_transcriptomic="$2/mrna";

LOG="$2/bowtie-build.log"

#extracted exon coordinates from refGene using:
cat $inputBed | bed12ToBed6 -i stdin -n > $exonsBed

awk '{print $1"\t"$2"\t"$3"\t"$4"_"$5"\t"$3-$2"\t"$6}' $exonsBed > $tempExonsBed

#extract exons sequences using exon coordinates file
fastaFromBed -fi $4 -bed $tempExonsBed -fo $tempExonsFASTA -s -name

#transform all to uppercase
tr [a-z] [A-Z] <  $tempExonsFASTA > $exonsFASTA

java -Xmx8G -Xms8G -cp $1/GenomicPTESDiscoverySTAR.jar:$1/commons-lang3-3%2e2%2e1.jar bio.igm.utils.init.MergeUCSCExonsToTranscript $exonsFASTA $2

java -Xmx8G -Xms8G -cp $1/GenomicPTESDiscoverySTAR.jar bio.igm.utils.init.ConvertBedToCoordinates $inputBed $2

bowtie-build $transcriptFASTA $mapper_transcriptomic 2>> $LOG
bowtie2-build $transcriptFASTA $mapper_transcriptomic 2>> $LOG
