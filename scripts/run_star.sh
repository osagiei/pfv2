#!/bin/sh


#module load apps/python
#module load apps/star/2.4.0j/gcc-4.4.6
#module load libs/pysam/0.7.7/gcc-4.4.6+python-2.7.3

module load apps/python27/2.7.8
module load apps/STAR/2.4.0j
module load libs/python/pysam/0.7.5

reads=$1
index=$2  #~/lustre/OsagieTemp/gatk_bundle/refs/STAR
fasta=$3 #fasta should be indexed using samtools
wd=$4

STAR --runThreadN 8 --genomeDir $index --genomeLoad LoadAndKeep \
	--readFilesIn $reads --chimSegmentMin 20 --chimScoreMin 1 \
	--alignIntronMax 100000 --outFilterMultimapNmax 7 \
	--outFilterMismatchNmax 2 --alignTranscriptsPerReadNmax 100000 \
	--outSAMattributes NH HI AS NM MD jM \
	--outFileNamePrefix ${wd}/star_


grep -v ^@ $wd/star_Aligned.out.sam > $wd/genomic.sam
