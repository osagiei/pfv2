#!/bin/sh

wd=$1
code=$2

join -o 1.1,1.2,1.3,1.4,1.6,2.1,2.2,2.3,2.4,2.6 $wd/left.sam $wd/right.sam | awk '($2 == 0 || $2 == 16) && ($7 == 0 || $7 == 16) && ($2 == $7) && ($3 == $8) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' - > $wd/merged_anchors.tsv

join -o 1.1,1.2,1.3,1.4,1.6,2.1,2.2,2.3,2.4,2.6 $wd/left.sam $wd/right.sam | awk '($2 == 0 || $2 == 16) && ($7 == 0 || $7 == 16) && ($2 != $7) && ($3 == $8) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' - > $wd/putative_sense_antisense_reads.tsv


join -o 1.1,1.2,1.3,1.4,1.6,2.1,2.2,2.3,2.4,2.6 $wd/left.sam $wd/right.sam | awk '($2 == 0 || $2 == 16) && ($7 == 0 || $7 == 16) && ($2 == $7) && ($3 != $8) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' - > $wd/putative_fused_transcripts_reads.tsv


java -cp $code/PTESDiscovery.jar bio.igm.utils.discovery.DetectExonShuffling $wd $wd/merged_anchors.tsv true

rm $wd/merged_anchors.tsv
