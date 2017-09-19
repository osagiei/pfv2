#!/bin/sh

module load apps/bedtools

transcripts_bed=$1
can=$2
path=$3

temp_locus=$path/temp_locus.bed
sum_can=$path/sum_can_juncs.bed
mean_can=$path/mean_can_juncs.bed
ptes=$path/annotated-ptes.bed
temp_ptes=$path/temp_ptes.bed

intersectBed -a $transcripts_bed -b $can -wo | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$6"\t"$17}' - > $temp_locus

groupBy -i $temp_locus -grp 1,2,3,4 -c 6 -o sum > $sum_can

groupBy -i $temp_locus -grp 1,2,3,4 -c 6 -o mean > $mean_can

rm $temp_locus


OLDIFS=$IFS; IFS=$'\n'; for i in `cat $ptes `; do refseq=`echo $i | awk '{print $4}' - | awk -F . '{print $1}' -`; sum=`awk -v x=$refseq '$4==x{print $5}' $sum_can `; mean=`awk -v x=$refseq '$4==x{print $5}' $mean_can `; echo -e "$i\t$sum\t$mean" ; done > $temp_ptes

#OLDIFS=$IFS; IFS=$'\n'; for i in `cat $ptes `; do refseq=`echo $i | awk '{print $4}' - | awk -F . '{print $1}' -`; sum=`awk -v x=$refseq '$4==x{print $0}' $sum_can `; mean=`grep $refseq $mean$


#awk -v x=$refseq '$4==x{print $0}'
#grep '^chr' $temp_ptes > $temp_ptes.tmp
#mv $temp_ptes $ptes
