# PFv2: a modified computational method to identify post-transcriptional exon shuffling (PTES) and baksplice events

LINKS:

    PTESTinder: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0881-4
    PTESFinder v1: https://sourceforge.net/projects/ptesfinder-v1/files/

CHANGES:

    Dependency on bowtie1 in the identification phase has been removed
    Identification of putative backsplice events using STAR
    Annotation-free backsplice identification

USAGE:

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
            -p PID -- should be between 0 and 1; ideal values 0.60 - 0.95, default: 0.85
            -j junction Span --should be an even integer, ideal values between 4 and 14, default: 8
            -C Maximum backsplice genomic span - default: 1000000

DEPENDENCIES:
    
    Python 3 - including argparse, subprocess, os, logging
    Java > 15
    STAR: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
    Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    Memory requirements are genome size and aligner dependent - see STAR aligner information for details