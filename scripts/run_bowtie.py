'''
Author: Osagie Izuogu

Description: Maps sequenced reads with Bowtie2, against the genome, the
             transcriptome, or the generated junction constructs.

             For Bowtie-specific parameter details, see:
             http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Date: 07/2020
'''

import argparse
import logging
import os
import subprocess
import sys

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s.%(msecs)03d %(name)-4s [%(levelname)-4s] %(message)s',
    datefmt='%d-%m-%Y %H:%M:%S')

LOG = logging.getLogger('run_bowtie')


def fail(message):
    LOG.error(message)
    sys.exit(1)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description='Map reads with Bowtie2.')
    parser.add_argument(
        '-o', '--output_dir', dest='output_dir', required=True,
        help='Output directory; results are written to <output_dir>/<sample_id>')
    parser.add_argument(
        '-r', '--reference_index', dest='reference_index', required=True,
        help='Path to the Bowtie2 index prefix, without the .bt2 suffix')
    parser.add_argument(
        '-f', '--fastq', dest='fastq', required=True,
        help='Single or merged FASTQ file')
    parser.add_argument(
        '-s', '--sample_id', dest='sample_id', required=True,
        help='Sample ID')
    parser.add_argument(
        '--logic_name', dest='logic_name', required=True,
        help='Name of this alignment, used for the output file names')
    parser.add_argument(
        '-t', '--threads', default='16',
        help='Number of threads (default: 16)')
    parser.add_argument(
        '--program', default='bowtie2',
        help='Bowtie2 executable to use (default: bowtie2)')
    return parser.parse_args(argv)


def check_index(prefix):
    '''Bowtie2 indexes use .bt2 for small references and .bt2l for large ones.'''
    if not (os.path.isfile(prefix + '.1.bt2') or os.path.isfile(prefix + '.1.bt2l')):
        fail('Bowtie2 index not found: expected %s.1.bt2 (or .bt2l)' % prefix)


def main(argv=None):
    args = parse_args(argv)

    if not os.path.isfile(args.fastq):
        fail('FASTQ not found: %s' % args.fastq)
    check_index(args.reference_index)

    sample_dir = os.path.join(args.output_dir, args.sample_id)
    os.makedirs(sample_dir, exist_ok=True)

    sam = os.path.join(sample_dir, '%s.sam' % args.logic_name)
    log_file = os.path.join(sample_dir, '%s_bowtie.log' % args.logic_name)

    command = [args.program,
               '-p', str(args.threads),
               '--very-sensitive',
               '--score-min=C,-15,0',
               '--mm',
               # Unaligned records carry no MD, NM or AS tag, so every downstream filter
               # skips them anyway; dropping them here keeps the SAM files far smaller.
               '--no-unal',
               '-x', args.reference_index,
               '-q', '-U', args.fastq,
               '-S', sam]

    LOG.info('Running: %s', ' '.join(command))
    try:
        with open(log_file, 'w') as handle:
            subprocess.run(command, check=True, stderr=handle)
    except FileNotFoundError:
        fail('Executable not found: %s' % args.program)
    except subprocess.CalledProcessError as error:
        fail('bowtie2 exited with status %d; see %s' % (error.returncode, log_file))

    if not os.path.isfile(sam):
        fail('bowtie2 did not produce %s' % sam)

    LOG.info('Finished mapping reads w/ bowtie2; Analysis: %s', args.logic_name)


if __name__ == '__main__':
    main()
