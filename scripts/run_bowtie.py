'''

Author: Osagie Izuogu

Description: Script for mapping sequenced reads using Bowtie2
             For Bowtie-specific parameter details,
             see: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

Date: 07/2020
'''

import subprocess
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
        '-r',
        '--reference_index',
        action="store",
        dest="reference_index",
        required=True,
        help="Path to Bowtie2 index")
    parser.add_argument(
        '-f',
        '--fastq',
        dest="fastq",
        action="store",
        required=True,
        help="Single or merged FASTQ files")
    parser.add_argument(
        '-s',
        '--sample_id',
        dest="sample_id",
        action="store",
        required=True,
        help="Sample ID")

    parser.add_argument('--logic_name', dest="logic_name", action="store_true")
    parser.add_argument(
        '-t',
        '--threads',
        default='16',
        help='Number of threads')

    args = parser.parse_args()

    COMMAND = "bowtie2 -p{0} " \
              "--very-sensitive --score-min=C,-15,0 --mm -x {1} " \
              "-q -U {2} 2> {3}/{4}/{5}_bowtie.log " \
              "-S {3}/{4}/{5}.sam ".\
        format(args.threads, args.reference_index, args.fastq, args.output_dir,
               args.sample_id, args.logic_name)

    SIGNAL = subprocess.check_call(COMMAND, shell=True, stdout=subprocess.PIPE)

    if SIGNAL == 0:
        logging.info("Finished mapping reads w/ bowtie2; Analysis: %s",
                     args.logic_name)
