'''

Author: Osagie Izuogu

Description: Script for mapping sequenced reads to the genome using STAR
             For STAR-specific parameter details,
             see: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

Date: 07/2020
'''

import subprocess
import os
import argparse
import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s.%(msecs)03d %(name)-4s [%(levelname)-4s] %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o',
        '--output_dir',
        action="store",
        dest="output_dir",
        required=True,
        help="Output directory required")

    parser.add_argument(
        '-g',
        '--genome_index',
        action="store",
        dest="genome_index",
        required=True,
        help="Path to STAR index; should exist for mapping, path is used for creating new index")
    parser.add_argument(
        '-f',
        '--genome_fasta',
        dest="genome_fasta",
        action="store",
        required=True,
        help="Genome fasta containing seq_region sequences")
    parser.add_argument(
        '-s',
        '--sample_id',
        dest="sample_id",
        action="store",
        required=True,
        help="Sample ID")
    parser.add_argument(
        '--fastq',
        nargs='+',
        help='FASTQ input. Format: fastq1 [fastq2]')
    parser.add_argument(
        '--prefix',
        default="star",
        help='Prefix for output file names')
    parser.add_argument(
        '--annotation_gtf',
        default=None,
        help='Annotation in GTF format')
    parser.add_argument('--outFilterMultimapNmax', default='10')
    parser.add_argument('--alignSJoverhangMin', default='8')
    parser.add_argument('--alignSJDBoverhangMin', default='1')
    parser.add_argument('--outFilterMismatchNmax', default='2')
    parser.add_argument('--outFilterMismatchNoverLmax', default='0.1')
    parser.add_argument('--alignIntronMin', default='20')
    parser.add_argument('--alignIntronMax', default='1000000')
    parser.add_argument('--alignMatesGapMax', default='1000000')
    parser.add_argument('--outFilterType', default='BySJout')
    parser.add_argument('--outFilterScoreMinOverLread', default='0.33')
    parser.add_argument('--outFilterMatchNminOverLread', default='0.33')
    parser.add_argument('--limitSjdbInsertNsj', default='1200000')
    parser.add_argument('--outSAMstrandField', default='intronMotif')
    parser.add_argument(
        '--outFilterIntronMotifs',
        default='None',
        help="Use 'RemoveNoncanonical' for Cufflinks compatibility")
    parser.add_argument('--alignSoftClipAtReferenceEnds', default='Yes')
    parser.add_argument(
        '--quantMode',
        default=[
            'TranscriptomeSAM',
            'GeneCounts'],
        nargs='+',
        help='Outputs read counts, and a BAM with reads in transcriptome coordinates')
    parser.add_argument(
        '--outSAMtype',
        default=[
            'BAM',
            'SortedByCoordinate'],
        nargs='+')
    parser.add_argument(
        '--outSAMunmapped',
        default='Within',
        help='Keep unmapped reads in output BAM')
    parser.add_argument(
        '--outSAMattrRGline',
        default=[
            'ID:rg1',
            'SM:sm1'],
        nargs='+',
        help='Adds read group line to BAM header; required by GATK')
    parser.add_argument(
        '--outSAMattributes',
        default=[
            'NH',
            'HI',
            'XS',
            'AS',
            'nM',
            'NM',
            'ch'],
        nargs='+')
    parser.add_argument(
        '--chimSegmentMin',
        default='20',
        help='Minimum chimeric segment length; switches on detection of chimeric (fusion) alignments')
    parser.add_argument(
        '--chimJunctionOverhangMin',
        default='20',
        help='Minimum overhang for a chimeric junction')
    parser.add_argument('--chimScoreMin', default=1)
    parser.add_argument(
        '--chimOutType',
        default=[
            'WithinBAM',
            'SoftClip'],
        nargs='+',
        help='')
    parser.add_argument('--chimMainSegmentMultNmax', default='1', help='')
    parser.add_argument('--genomeLoad', default='NoSharedMemory')
    parser.add_argument(
        '--sjdbFileChrStartEnd',
        default=None,
        help='SJ.out.tab file (e.g., from 1st pass). With this option, only one pass will be run')
    parser.add_argument(
        '--STARlong',
        action='store_true',
        help='Use STARlong instead of STAR')
    parser.add_argument(
        '-t',
        '--threads',
        default='16',
        help='Number of threads')
    parser.add_argument(
        '--build_genome_index',
        dest="build_genome_index",
        action="store_true")
    parser.add_argument('--limitGenomeGenerateRAM', default=36183246720)
    parser.add_argument('--program', default="STAR")

    args = parser.parse_args()

    if args.build_genome_index:
        if not os.path.exists(args.genome_index):
            os.makedirs(args.genome_index)

        COMMAND = "{0} --runThreadN {1} " \
                  "--runMode genomeGenerate " \
                  "--genomeFastaFiles {2} " \
                  "--genomeLoad {3} " \
                  "--genomeDir {4} " \
                  "--limitGenomeGenerateRAM {5}". \
            format(args.program, args.threads, args.genome_fasta, args.genomeLoad,
                   args.genome_index, args.limitGenomeGenerateRAM)

        SIGNAL = subprocess.check_call(
            COMMAND, shell=True, stdout=subprocess.PIPE)

        if SIGNAL == 0:
            logging.info("Finished building STAR index for genome")

    else:
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)

        FASTQ = ' '.join(args.fastq)

        COMMAND = "{0} --runThreadN {1} " \
                  "--genomeDir {2} " \
                  "--genomeLoad {3} " \
                  "--readFilesIn {4} " \
                  "--alignIntronMax {5} " \
                  "--outFilterMultimapNmax {6} " \
                  "--outFilterMismatchNmax {7} " \
                  "--outSAMtype {8} " \
                  "--outSAMattributes {9} " \
                  "--outFileNamePrefix {10}/{11}/{12}_ ". \
            format(args.program, args.threads, args.genome_index, args.genomeLoad, FASTQ, args.alignIntronMax,
                   args.outFilterMultimapNmax, args.outFilterMismatchNmax, ' '.join(args.outSAMtype),
                   ' '.join(args.outSAMattributes), args.output_dir, args.sample_id, args.prefix)

        # identify chimeric reads
        COMMAND += ' --chimSegmentMin {0}' ' --chimScoreMin {1}'.format(
            args.chimSegmentMin, args.chimScoreMin)

        if args.fastq[0].endswith('.gz'):
            COMMAND += ' --readFilesCommand zcat'

        SIGNAL = subprocess.check_call(
            COMMAND, shell=True, stdout=subprocess.PIPE)

        if SIGNAL == 0:
            logging.info("Finished mapping reads to genome with STAR")

            # index BAM

            COMMAND = "samtools index -b {0}/{1}/{2}_Aligned.sortedByCoord.out.bam". \
                format(args.output_dir, args.sample_id, args.prefix)

            SIGNAL = subprocess.check_call(
                COMMAND, shell=True, stdout=subprocess.PIPE)

            if SIGNAL == 0:
                logging.info("Finished indexing BAM file")

                os.system(
                    "readlink -f {0}/{1}/{2}_Aligned.sortedByCoord.out.bam >> {0}/bams.list". format(
                        args.output_dir, args.sample_id, args.prefix))
