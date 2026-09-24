'''
Author: Osagie Izuogu

Description: Maps sequenced reads to the genome with STAR, with chimeric
             detection enabled so that backsplice candidates are reported in
             Chimeric.out.junction.

             For STAR-specific parameter details, see:
             https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

Date: 07/2020
'''

import argparse
import logging
import os
import shutil
import subprocess
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s.%(msecs)03d %(name)-4s [%(levelname)-4s] %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

LOG = logging.getLogger('run_star')

# Outputs the rest of the pipeline reads back.
REQUIRED_OUTPUTS = ('Chimeric.out.junction', 'SJ.out.tab')


def fail(message):
    LOG.error(message)
    sys.exit(1)


def run(command):
    '''Runs a command as an argument list, without a shell.'''
    LOG.info('Running: %s', ' '.join(command))
    try:
        subprocess.run(command, check=True)
    except FileNotFoundError:
        fail('Executable not found: %s' % command[0])
    except subprocess.CalledProcessError as error:
        fail('%s exited with status %d' % (command[0], error.returncode))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description='Map reads to the genome with STAR, detecting chimeric junctions.')
    parser.add_argument(
        '-o', '--output_dir', dest='output_dir', required=True,
        help='Output directory; results are written to <output_dir>/<sample_id>')
    parser.add_argument(
        '-g', '--genome_index', dest='genome_index', required=True,
        help='Path to the STAR index; must exist for mapping, is created by --build_genome_index')
    parser.add_argument(
        '-s', '--sample_id', dest='sample_id', required=True,
        help='Sample ID')
    parser.add_argument(
        '-f', '--genome_fasta', dest='genome_fasta', default=None,
        help='Genome FASTA; required only with --build_genome_index')
    parser.add_argument(
        '--fastq', nargs='+',
        help='FASTQ input. Format: fastq1 [fastq2]')
    parser.add_argument(
        '--prefix', default='star',
        help='Prefix for output file names (default: star)')
    parser.add_argument(
        '-t', '--threads', default='16',
        help='Number of threads (default: 16)')
    parser.add_argument('--program', default='STAR',
                        help='STAR executable to use (default: STAR)')

    # Alignment parameters the pipeline depends on. Changing these changes which
    # backsplice candidates are reported.
    parser.add_argument('--alignIntronMax', default='1000000')
    parser.add_argument('--outFilterMultimapNmax', default='10')
    parser.add_argument('--outFilterMismatchNmax', default='2')
    parser.add_argument('--chimSegmentMin', default='20',
                        help='Minimum chimeric segment length; enables chimeric detection')
    parser.add_argument('--chimScoreMin', default='1',
                        help='Minimum score for a chimeric alignment')
    parser.add_argument('--genomeLoad', default='NoSharedMemory')
    parser.add_argument('--outSAMtype', default=['BAM', 'SortedByCoordinate'], nargs='+')
    parser.add_argument('--outSAMattributes',
                        default=['NH', 'HI', 'XS', 'AS', 'nM', 'NM', 'ch'], nargs='+')
    parser.add_argument('--extra', nargs=argparse.REMAINDER, default=[],
                        help='Additional arguments passed straight through to STAR')

    # Index building.
    parser.add_argument('--build_genome_index', dest='build_genome_index', action='store_true')
    parser.add_argument('--limitGenomeGenerateRAM', default='36183246720')
    parser.add_argument('--skip_bam_index', action='store_true',
                        help='Do not index the sorted BAM; the pipeline does not require it')

    return parser.parse_args(argv)


def build_index(args):
    if not args.genome_fasta:
        fail('--genome_fasta is required with --build_genome_index')
    if not os.path.isfile(args.genome_fasta):
        fail('Genome FASTA not found: %s' % args.genome_fasta)

    os.makedirs(args.genome_index, exist_ok=True)
    run([args.program,
         '--runThreadN', str(args.threads),
         '--runMode', 'genomeGenerate',
         '--genomeFastaFiles', args.genome_fasta,
         '--genomeLoad', args.genomeLoad,
         '--genomeDir', args.genome_index,
         '--limitGenomeGenerateRAM', str(args.limitGenomeGenerateRAM)])
    LOG.info('Finished building the STAR index for the genome')


def map_reads(args):
    if not args.fastq:
        fail('--fastq is required when mapping')
    for fastq in args.fastq:
        if not os.path.isfile(fastq):
            fail('FASTQ not found: %s' % fastq)
    if not os.path.isdir(args.genome_index):
        fail('STAR index directory not found: %s' % args.genome_index)

    sample_dir = os.path.join(args.output_dir, args.sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    command = [args.program,
               '--runThreadN', str(args.threads),
               '--genomeDir', args.genome_index,
               '--genomeLoad', args.genomeLoad,
               '--readFilesIn'] + list(args.fastq) + [
               '--alignIntronMax', str(args.alignIntronMax),
               '--outFilterMultimapNmax', str(args.outFilterMultimapNmax),
               '--outFilterMismatchNmax', str(args.outFilterMismatchNmax),
               '--outSAMtype'] + list(args.outSAMtype) + [
               '--outSAMattributes'] + list(args.outSAMattributes) + [
               '--outFileNamePrefix', os.path.join(sample_dir, args.prefix + '_'),
               # Chimeric detection supplies the backsplice candidates.
               '--chimSegmentMin', str(args.chimSegmentMin),
               '--chimScoreMin', str(args.chimScoreMin),
               # Pinned explicitly: the pipeline reads Chimeric.out.junction, and
               # relying on the STAR default here breaks across releases.
               '--chimOutType', 'Junctions']

    if args.fastq[0].endswith('.gz'):
        command += ['--readFilesCommand', 'zcat']
    command += list(args.extra)

    run(command)
    LOG.info('Finished mapping reads to the genome with STAR')

    for name in REQUIRED_OUTPUTS:
        path = os.path.join(sample_dir, '%s_%s' % (args.prefix, name))
        if not os.path.isfile(path):
            fail('STAR did not produce the expected output: %s' % path)

    index_bam(args, sample_dir)


def index_bam(args, sample_dir):
    bam = os.path.join(sample_dir, '%s_Aligned.sortedByCoord.out.bam' % args.prefix)
    if not os.path.isfile(bam):
        LOG.warning('Sorted BAM not found at %s; skipping indexing', bam)
        return

    with open(os.path.join(args.output_dir, 'bams.list'), 'a') as handle:
        handle.write(os.path.abspath(bam) + '\n')

    if args.skip_bam_index:
        return
    if shutil.which('samtools') is None:
        LOG.warning('samtools not found on PATH; skipping BAM indexing')
        return

    try:
        subprocess.run(['samtools', 'index', '-b', bam], check=True)
        LOG.info('Finished indexing the BAM file')
    except subprocess.CalledProcessError as error:
        # The pipeline reads the junction files, not the BAM, so this is not fatal.
        LOG.warning('samtools index failed with status %d; continuing', error.returncode)


def main(argv=None):
    args = parse_args(argv)
    if args.build_genome_index:
        build_index(args)
    else:
        map_reads(args)


if __name__ == '__main__':
    main()
