#!/bin/sh
############################################## check_java_version
_java=java
if [[ "$_java" ]]; then
    # VERSION=$("$_java" -Xms500M -Xmx500M -version 2>&1 | awk -F '"' '/version/ {print $2}')
    VERSION=$($_java -version 2>&1 \
          | head -1 \
          | cut -d'"' -f2 \
          | sed 's/^1\.//' \
          | cut -d'.' -f1)
    echo version "$VERSION"
    if [[ "$VERSION" -gt 15 ]]; then
        echo ">>> Java version check completed -- okay "
    else
        echo -e "Error: Java version check failed, please re-run with version higher than 15 \nOr try re-compiling by running 'sh setup.sh' before re-running -- exiting "
        exit 1
    fi
fi

############################################### usage_info

usage(){
echo -e "\n*** PTESFinder v. 2 ***\nTo run PFv2, ensure STAR, bedtools and bowtie2 are installed on your system and on your path. Also, ensure that your system can execute Java programs; minimum version: 1.6.

Input Files:
        RNASeq data in Illumina FASTQ format
        Genome reference in FASTA format
        Pre-built bowtie2 genome reference index
        Pre-built bowtie2 transcriptome reference index
        Path to pre-built STAR genome reference index files

Running PTESFinder:
        $ sh PFv2.sh <options>

Parameters:

Mandatory:
        -r sequence reads in FASTQ format
        -d working directory
        -i sample_id
        -t transcript reference bowtie index
        -g genome reference in FASTA format
        -b genome reference bowtie2 index
        -c PFv2 code directory
        -l average read length
        -S path to STAR genome reference index files
Optional:
        -p PID -- should be <= 1; ideal values between 0.60 and 0.95, default: 0.85
        -j junction Span --should be an even integer, ideal values between 4 and 14, default: 8
        -G turn off all filters flag and run only genomic and junctional filters
        -T turn off all filters flag and run only transcriptomic and junctional filters
        -C Maximum backsplice genomic span - default: 1000000

#example run command:
sh PFv2.sh \
  -i SRR364679 \
  -r SRR364679.fastq \
  -d SRR364679/ \
  -S STAR/ \
  -t transcriptome-index-bowtie \
  -g genome.fasta \
  -b genome-index-bowtie \
  -l 100 \
  -c pfv2/

Note:
  - Path to pre-built STAR index should contain the following files:
      chrLength.txt
      chrNameLength.txt
      chrName.txt
      chrStart.txt
      Genome
      genomeParameters.txt
      SA
      SAindex

  - Bowtie2 index files with extensions .bt2 required. In the example above, transcriptome-index-bowtie.*.bt2

  - PE reads should be pooled into a single FASTQ file with unique read ids.


#email support: osagie.izuogu@gmail.com";
exit 1;
}
#####################################################################

GFLAG=false
BFLAG=false
SFLAG=false
TFLAG=false
LFLAG=false

JSPAN=8
PID=0.85
MAX_GENOMIC_SPAN=1000000
MIN_OVERHANG=15

while getopts ":r:i:d:t:g:b:l:p:j:c:S:C:h" opt; do
	case $opt in
		r)	FASTQ_READS="$OPTARG";;
		d)	OUTPUT_DIR="$OPTARG"
			if [ ! -d "$OUTPUT_DIR" ]; then
				echo "Error: The path provided is not a directory, using current directory"
				OUTPUT_DIR=`pwd`

			fi
			OUTPUT_DIR=$OUTPUT_DIR/PF;;
	  i)  SAMPLE_ID="$OPTARG";;
		t)	TFLAG=true; TRANSCRIPTOME_INDEX="$OPTARG";;
		g)	GFLAG=true; GENOME_FASTA="$OPTARG";;
		b)	BFLAG=true; GENOME_BOWTIE_INDEX="$OPTARG";;
		S)	SFLAG=true; GENOME_STAR_INDEX="$OPTARG";;
    l)  LFLAG=true; READ_LENGTH="$OPTARG";;
		p)	PID="$OPTARG";;
		j)	JSPAN="$OPTARG";;
    C)  MAX_GENOMIC_SPAN="$OPTARG";;
		c)	if [ ! -d "$OPTARG" ]; then
           echo "Error: The path provided is not a directory, please provide path to directory containing PFv2.jar..exiting"
           exit 1
        fi
			CODEBASE="$OPTARG";;

		h)	usage;;
		\?)
		   echo -e "Invalid option: -$OPTARG\n"
	  		 usage
	  		 ;;
		:)
		   echo "Option -$OPTARG requires an argument!"
			usage
		        exit 1
           		 ;;
	esac
done

if ! $GFLAG || ! $SFLAG || ! $TFLAG || ! $BFLAG
then
	echo "Some mandatory options not provided, please ensure you have supplied a path to genomic reference FASTA, Bowtie2 transcriptome index and genome index paths for Bowtie2 & STAR...exiting!"
	usage
	exit 1
fi

echo -e ">>> Starting PTESFinder v.2 `date`..\n"

#################################### initialization
SEGMENT_SIZE=$(($READ_LENGTH - $MIN_OVERHANG))

WORKING_DIR=${OUTPUT_DIR}/${SAMPLE_ID}/
mkdir -p ${WORKING_DIR}

echo -e "INFO: Started mapping reads to the genome with STAR and Bowtie2 `date` \n"

python ${CODEBASE}/scripts/run_star.py \
  --program STAR \
  -s $SAMPLE_ID  \
  --fastq  $FASTQ_READS \
  --output_dir $OUTPUT_DIR \
  --genome_index $GENOME_STAR_INDEX \
  --genome_fasta $GENOME_FASTA

python ${CODEBASE}/scripts/run_bowtie.py \
  -s $SAMPLE_ID  \
  --fastq  $FASTQ_READS \
  --output_dir $OUTPUT_DIR \
  --reference_index $GENOME_BOWTIE_INDEX \
  --logic_name "genomic"

python ${CODEBASE}/scripts/run_bowtie.py \
  -s $SAMPLE_ID  \
  --fastq  $FASTQ_READS \
  --output_dir $OUTPUT_DIR \
  --reference_index $TRANSCRIPTOME_INDEX \
  --logic_name "transcriptomic"

echo -e "SUCCESS: Finished mapping reads to the genome and transcriptome `date`\n"
echo -e "INFO: Screening mapped reads for putative backsplice junctions\n"

java -Xms20G -Xmx20G -cp ${CODEBASE}/PFv2.jar:${CODEBASE}/apache-lib/commons-lang3-3.2.1.jar \
    bio.igm.utils.discovery.ProcessShuffledCoordinates \
    $WORKING_DIR \
    $MAX_GENOMIC_SPAN \
    $SEGMENT_SIZE

echo -e "SUCCESS: Finished identifying putative backsplice junctions `date`\n"
echo -e "INFO: Generating sequence constructs and evaluating models\n"

java -Xms20G -Xmx20G -cp ${CODEBASE}/PFv2.jar:${CODEBASE}/apache-lib/commons-lang3-3.2.1.jar \
    bio.igm.utils.discovery.GenerateSequenceConstructsGenome \
    $WORKING_DIR \
    $GENOME_FASTA

bowtie2-build ${WORKING_DIR}/Can.fa ${WORKING_DIR}/canonical 2>> ${WORKING_DIR}/bowtie-build.log
bowtie2-build ${WORKING_DIR}/Constructs.fa ${WORKING_DIR}/ptes 2>> ${WORKING_DIR}/bowtie-build.log

python ${CODEBASE}/scripts/run_bowtie.py \
  -s $SAMPLE_ID  \
  --fastq  $FASTQ_READS \
  --output_dir $OUTPUT_DIR \
  --reference_index ${WORKING_DIR}/ptes \
  --logic_name "ptes"

python ${CODEBASE}/scripts/run_bowtie.py \
  -s $SAMPLE_ID  \
  --fastq  $FASTQ_READS \
  --output_dir $OUTPUT_DIR \
  --reference_index ${WORKING_DIR}/canonical \
  --logic_name "canonical"

echo -e "SUCCESS: Finished evaluating putative backsplice models `date` \n"
echo -e "INFO: Filtering potential false positive predictions\n"

java -Xms20G -Xmx20G -cp ${CODEBASE}/PFv2.jar:${CODEBASE}/apache-lib/commons-lang3-3.2.1.jar \
    bio.igm.utils.filter.PipelineFilter \
    $WORKING_DIR \
    $JSPAN \
    $PID 1 0 0

# clean up
if [ -f ${WORKING_DIR}/pf-structures.bed ]
then
  CIRC_READS=$(wc -l ${WORKING_DIR}/pf-supporting-reads.tab | awk '{print $1}')
  CJUNCS_READS=$(awk '{sum += $5}END{print sum}' ${WORKING_DIR}/pf-flanking-canonical-junctions.bed)
  TJR=$(( $CIRC_READS + $CJUNCS_READS ))

  echo -e "INFO: Total number of identified circRNAs: $(wc -l ${WORKING_DIR}/pf-structures.bed | awk '{print $1}')"
  echo -e "INFO: Total number of circRNA supporting reads: $CIRC_READS"
  echo -e "INFO: Total number of canonical reads: $CJUNCS_READS"

  awk -v X=$TJR '{print $1":"$2"-"$3":"$6"\t"$5"\t"($5 / X) * 1000000}' ${WORKING_DIR}/pf-structures.bed > ${WORKING_DIR}/${SAMPLE_ID}_jpms.tsv
  rm -r $WORKING_DIR/*.{ebwt,bt2,fa,sam,bam}
  echo -e "\n>>> Finished run successfully @ `date` ..\n"

else
  echo -e "WARNING: Run finished @ `date` without generating final output, check logs \n"
fi

