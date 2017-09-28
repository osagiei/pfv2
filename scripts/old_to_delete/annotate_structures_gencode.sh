#!/bin/sh

module load apps/bedtools

ptes=$1
id=$2
path=$3

bed="~/lustre/OsagieTemp/refs/gtf/gencode_v19_basic_genes.bed"

intersectBed -a $ptes -b $bed-wo -f 0.999999 | sort -k 4,4 - | uniq -f3 -w 24 > ${path}/${id}.tsv 
