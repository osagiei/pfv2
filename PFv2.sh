#!/usr/bin/env bash
#
# PTESFinder v2 - annotation-free identification of post-transcriptional exon
# shuffling (PTES) / backsplice junctions from RNA-seq data.
#
# Re-exec under bash when invoked as `sh PFv2.sh`, since this script relies on
# bash arithmetic, [[ ]] and pipefail.
if [ -z "${BASH_VERSION:-}" ]; then
  exec bash "$0" "$@"
fi

set -euo pipefail

readonly VERSION="2.2.1"
readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

########################################################################## defaults
JSPAN=8
PID=0.85
MAX_GENOMIC_SPAN=1000000
MIN_GENOMIC_SPAN=50
MIN_OVERHANG=15
THREADS=16
JAVA_MEM=20G
MIN_JAVA_VERSION=15

CODEBASE="$SCRIPT_DIR"
PYTHON="${PYTHON:-python3}"

FASTQ_READS=""
OUTPUT_DIR=""
SAMPLE_ID=""
TRANSCRIPTOME_INDEX=""
GENOME_FASTA=""
GENOME_BOWTIE_INDEX=""
GENOME_STAR_INDEX=""
READ_LENGTH=""

GENOMIC_ONLY=false
TRANSCRIPTOMIC_ONLY=false
KEEP_INTERMEDIATES=false
LEGACY=false
SKIP_VALIDATION=false
NORMALISE_STRAND=true

########################################################################## helpers
log()  { printf '%s [INFO]  %s\n'  "$(date '+%Y-%m-%d %H:%M:%S')" "$*"; }
warn() { printf '%s [WARN]  %s\n'  "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; }
die()  { printf '%s [ERROR] %s\n'  "$(date '+%Y-%m-%d %H:%M:%S')" "$*" >&2; exit 1; }

on_error() {
  local exit_code=$?
  warn "PTESFinder failed at line ${BASH_LINENO[0]} with exit status ${exit_code}"
  if [[ -n "${WORKING_DIR:-}" && -f "${WORKING_DIR}/run.log" ]]; then
    warn "See ${WORKING_DIR}/run.log for the Java stage logs"
  fi
  exit "$exit_code"
}
trap on_error ERR

usage() {
cat <<'USAGE'

*** PTESFinder v2 ***

Identifies post-transcriptional exon shuffling (backsplice) junctions in RNA-seq
data without relying on an annotation.

Input files:
        RNA-seq reads in Illumina FASTQ format (plain or gzipped)
        Genome reference in FASTA format
        Pre-built STAR genome index
        Pre-built Bowtie2 genome index
        Pre-built Bowtie2 transcriptome index

Usage:
        bash PFv2.sh -i <sample> -r <reads.fastq> -d <dir> -S <star_index> \
                     -t <transcriptome_index> -g <genome.fa> -b <genome_index> -l <read_length>

Mandatory:
        -r  sequence reads in FASTQ format
        -d  working directory
        -i  sample id
        -t  transcriptome reference Bowtie2 index prefix
        -g  genome reference in FASTA format
        -b  genome reference Bowtie2 index prefix
        -l  average read length
        -S  path to the pre-built STAR genome index directory

Optional:
        -c  PFv2 code directory (default: the directory holding this script)
        -p  minimum percent identity per flank, 0-1; ideal 0.60-0.95 (default: 0.85)
        -j  junction span, even integer; ideal 4-14 (default: 8)
        -C  maximum backsplice genomic span in bp (default: 1000000)
        -M  minimum backsplice genomic span in bp (default: 50)
        -n  threads for STAR and Bowtie2 (default: 16)
        -m  Java heap for the PFv2 stages, e.g. 8G (default: 20G)
        -G  run the genomic filter only, skipping the transcriptomic comparison
        -T  run the transcriptomic filter only, skipping the genomic comparison
        -k  keep intermediate SAM/FASTA/index files instead of deleting them
        -L  reproduce the filter semantics of releases up to 2.1.0 (see CHANGELOG)
        -V  skip input and reference validation
        -A  report the aligned strand STAR assigned instead of the strand implied by
            the splice motif; for reproducing pre-2.2.0 output only, see CHANGELOG
        -h  show this message

Example:
        bash PFv2.sh \
          -i SRR364679 \
          -r SRR364679.fastq \
          -d SRR364679/ \
          -S STAR/ \
          -t transcriptome-index-bowtie \
          -g genome.fasta \
          -b genome-index-bowtie \
          -l 100

Notes:
  - The STAR index directory must contain chrLength.txt, chrNameLength.txt,
    chrName.txt, chrStart.txt, Genome, genomeParameters.txt, SA and SAindex.
  - Bowtie2 index prefixes are given without the .bt2 suffix, e.g. pass
    "transcriptome-index-bowtie" for transcriptome-index-bowtie.1.bt2.
  - Sequence names in the genome FASTA must match those used to build the STAR
    index ("chr1" and "1" are not interchangeable).
  - Paired-end reads must be pooled into a single FASTQ with unique read ids.

Dependencies: STAR, Bowtie2, samtools, Python 3, Java 16 or newer.

email support: osagie.izuogu@gmail.com
USAGE
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "Required executable not found on PATH: $1"
}

require_file() {
  [[ -f "$1" ]] || die "$2 not found: $1"
  [[ -r "$1" ]] || die "$2 is not readable: $1"
}

require_bowtie2_index() {
  local prefix="$1" label="$2"
  local found=false f
  for f in "${prefix}".1.bt2 "${prefix}".1.bt2l; do
    [[ -f "$f" ]] && found=true
  done
  $found || die "$label Bowtie2 index not found: expected ${prefix}.1.bt2 (or .bt2l)"
}

require_star_index() {
  local dir="$1" f
  [[ -d "$dir" ]] || die "STAR index directory not found: $dir"
  for f in chrName.txt chrStart.txt Genome SA SAindex genomeParameters.txt; do
    [[ -f "${dir}/${f}" ]] || die "STAR index is incomplete: ${dir}/${f} is missing"
  done
}

is_positive_int() { [[ "$1" =~ ^[0-9]+$ ]] && [[ "$1" -gt 0 ]]; }

check_java_version() {
  require_command java
  local version
  version="$(java -version 2>&1 | head -1 | cut -d'"' -f2 | sed 's/^1\.//' | cut -d'.' -f1)"
  if ! [[ "$version" =~ ^[0-9]+$ ]]; then
    warn "Could not determine the Java version; continuing"
    return 0
  fi
  if [[ "$version" -le "$MIN_JAVA_VERSION" ]]; then
    die "Java $version found, but PFv2 needs $((MIN_JAVA_VERSION + 1)) or newer. \
Install a newer JDK, or recompile for your JDK with 'bash setup.sh'."
  fi
  log "Java $version detected"
}

########################################################################## arguments
while getopts ":r:i:d:t:g:b:l:p:j:c:S:C:M:n:m:GTkLVAh" opt; do
  case $opt in
    r) FASTQ_READS="$OPTARG" ;;
    d) OUTPUT_DIR="$OPTARG" ;;
    i) SAMPLE_ID="$OPTARG" ;;
    t) TRANSCRIPTOME_INDEX="$OPTARG" ;;
    g) GENOME_FASTA="$OPTARG" ;;
    b) GENOME_BOWTIE_INDEX="$OPTARG" ;;
    S) GENOME_STAR_INDEX="$OPTARG" ;;
    l) READ_LENGTH="$OPTARG" ;;
    p) PID="$OPTARG" ;;
    j) JSPAN="$OPTARG" ;;
    C) MAX_GENOMIC_SPAN="$OPTARG" ;;
    M) MIN_GENOMIC_SPAN="$OPTARG" ;;
    n) THREADS="$OPTARG" ;;
    m) JAVA_MEM="$OPTARG" ;;
    c) CODEBASE="$OPTARG" ;;
    G) GENOMIC_ONLY=true ;;
    T) TRANSCRIPTOMIC_ONLY=true ;;
    k) KEEP_INTERMEDIATES=true ;;
    L) LEGACY=true ;;
    V) SKIP_VALIDATION=true ;;
    A) NORMALISE_STRAND=false ;;
    h) usage; exit 0 ;;
    \?) printf 'Invalid option: -%s\n' "$OPTARG" >&2; usage >&2; exit 2 ;;
    :)  printf 'Option -%s requires an argument\n' "$OPTARG" >&2; usage >&2; exit 2 ;;
  esac
done

########################################################################## validation
missing=()
[[ -n "$FASTQ_READS"         ]] || missing+=("-r sequence reads")
[[ -n "$OUTPUT_DIR"          ]] || missing+=("-d working directory")
[[ -n "$SAMPLE_ID"           ]] || missing+=("-i sample id")
[[ -n "$TRANSCRIPTOME_INDEX" ]] || missing+=("-t transcriptome Bowtie2 index")
[[ -n "$GENOME_FASTA"        ]] || missing+=("-g genome FASTA")
[[ -n "$GENOME_BOWTIE_INDEX" ]] || missing+=("-b genome Bowtie2 index")
[[ -n "$GENOME_STAR_INDEX"   ]] || missing+=("-S STAR genome index")
[[ -n "$READ_LENGTH"         ]] || missing+=("-l average read length")

if (( ${#missing[@]} > 0 )); then
  printf 'Missing mandatory option(s):\n' >&2
  printf '  %s\n' "${missing[@]}" >&2
  printf '\n' >&2
  usage >&2
  exit 2
fi

is_positive_int "$READ_LENGTH" || die "-l read length must be a positive integer, got '$READ_LENGTH'"
is_positive_int "$JSPAN"       || die "-j junction span must be a positive integer, got '$JSPAN'"
is_positive_int "$THREADS"     || die "-n threads must be a positive integer, got '$THREADS'"
is_positive_int "$MAX_GENOMIC_SPAN" || die "-C maximum span must be a positive integer, got '$MAX_GENOMIC_SPAN'"
is_positive_int "$MIN_GENOMIC_SPAN" || die "-M minimum span must be a positive integer, got '$MIN_GENOMIC_SPAN'"

(( JSPAN % 2 == 0 )) || die "-j junction span must be an even integer, got $JSPAN"
(( JSPAN >= 2 ))     || die "-j junction span must be at least 2, got $JSPAN"
(( READ_LENGTH > MIN_OVERHANG )) \
  || die "-l read length ($READ_LENGTH) must exceed the minimum overhang ($MIN_OVERHANG)"
(( MAX_GENOMIC_SPAN >= MIN_GENOMIC_SPAN )) \
  || die "-C maximum span ($MAX_GENOMIC_SPAN) is below -M minimum span ($MIN_GENOMIC_SPAN)"

awk -v p="$PID" 'BEGIN { exit !(p > 0 && p <= 1) }' \
  || die "-p percent identity must be in (0, 1], got '$PID'"

[[ "$JAVA_MEM" =~ ^[0-9]+[kKmMgG]?$ ]] || die "-m Java heap must look like 8G or 4096M, got '$JAVA_MEM'"

if $GENOMIC_ONLY && $TRANSCRIPTOMIC_ONLY; then
  warn "-G and -T were both given; running both filters, which is the default"
  GENOMIC_ONLY=false
  TRANSCRIPTOMIC_ONLY=false
fi

########################################################################## preflight
log "PTESFinder v${VERSION} starting"

require_command STAR
require_command bowtie2
require_command bowtie2-build
require_command samtools
require_command "$PYTHON"
check_java_version

"$PYTHON" -c 'import sys; sys.exit(0 if sys.version_info[0] >= 3 else 1)' \
  || die "$PYTHON is not Python 3; set PYTHON=/path/to/python3"

[[ -d "$CODEBASE" ]] || die "-c code directory not found: $CODEBASE"
require_file "${CODEBASE}/PFv2.jar" "PFv2.jar"
require_file "${CODEBASE}/scripts/run_star.py" "run_star.py"
require_file "${CODEBASE}/scripts/run_bowtie.py" "run_bowtie.py"
require_file "${CODEBASE}/scripts/validate_inputs.py" "validate_inputs.py"

COMMONS_LANG="$(find "${CODEBASE}/lib" -name 'commons-lang3-*.jar' -not -name '._*' 2>/dev/null | head -1 || true)"
[[ -n "$COMMONS_LANG" ]] || die "commons-lang3 jar not found in ${CODEBASE}/lib"

require_file "$FASTQ_READS" "FASTQ reads"
require_file "$GENOME_FASTA" "Genome FASTA"
require_star_index "$GENOME_STAR_INDEX"
require_bowtie2_index "$GENOME_BOWTIE_INDEX" "Genome"
require_bowtie2_index "$TRANSCRIPTOME_INDEX" "Transcriptome"

if $SKIP_VALIDATION; then
  warn "Skipping input and reference validation (-V)"
else
  log "Validating inputs and references"
  # Contig naming that differs between the genome FASTA and the STAR index produces an
  # empty result many hours later, so it is worth a few seconds up front.
  "$PYTHON" "${CODEBASE}/scripts/validate_inputs.py" \
    --fastq "$FASTQ_READS" \
    --genome "$GENOME_FASTA" \
    --star-index "$GENOME_STAR_INDEX" \
    --bowtie2-genome "$GENOME_BOWTIE_INDEX" \
    --bowtie2-transcriptome "$TRANSCRIPTOME_INDEX" \
    --read-length "$READ_LENGTH" \
    || die "Validation failed. Fix the errors above, or pass -V to run anyway."
fi

mkdir -p "$OUTPUT_DIR" || die "Cannot create working directory: $OUTPUT_DIR"
OUTPUT_DIR="$(cd "$OUTPUT_DIR" && pwd)/PF"
WORKING_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
mkdir -p "$WORKING_DIR" || die "Cannot create working directory: $WORKING_DIR"

# Each construct arm is the read length less the minimum overhang, so a read can cross
# the seam with at least MIN_OVERHANG bases on its short side. The floor keeps a short
# read library from producing arms too small to host the junction window at all.
SEGMENT_SIZE=$(( READ_LENGTH - MIN_OVERHANG ))
SEGMENT_FLOOR=$(( JSPAN + MIN_OVERHANG ))
if (( SEGMENT_SIZE < SEGMENT_FLOOR )); then
  warn "Construct arm of ${SEGMENT_SIZE} bp is below the floor of ${SEGMENT_FLOOR} bp; using the floor"
  SEGMENT_SIZE=$SEGMENT_FLOOR
fi
CLASSPATH="${CODEBASE}/PFv2.jar:${COMMONS_LANG}"
JAVA_OPTS=(-Xms"${JAVA_MEM}" -Xmx"${JAVA_MEM}")

if $GENOMIC_ONLY; then
  FILTER_ARGS=(0 1 0); FILTER_DESC="genomic only"
elif $TRANSCRIPTOMIC_ONLY; then
  FILTER_ARGS=(0 0 1); FILTER_DESC="transcriptomic only"
else
  FILTER_ARGS=(1 0 0); FILTER_DESC="genomic and transcriptomic"
fi

# Legacy mode reproduces pre-2.2.0 output, which reported the aligned strand.
$LEGACY && NORMALISE_STRAND=false
if $NORMALISE_STRAND; then NORMALISE_ARG=1; else NORMALISE_ARG=0; fi

if $LEGACY; then
  LEGACY_ARG=1
  FILTER_DESC="$FILTER_DESC (legacy semantics)"
else
  LEGACY_ARG=0
fi

log "Sample:        $SAMPLE_ID"
log "Working dir:   $WORKING_DIR"
log "Read length:   $READ_LENGTH (construct arm: $SEGMENT_SIZE bp)"
log "Backsplice span: ${MIN_GENOMIC_SPAN}-${MAX_GENOMIC_SPAN} bp"
log "Junction span:   $JSPAN, percent identity: $PID"
log "Filters:         $FILTER_DESC"
if $NORMALISE_STRAND; then
  log "Strand:          from the splice motif"
else
  warn "Strand:          as STAR aligned it (-A); pre-2.2.0 convention, see CHANGELOG"
fi
$LEGACY && warn "Legacy mode: results reproduce releases up to 2.1.0, not the current defaults"
log "Threads:         $THREADS, Java heap: $JAVA_MEM"

########################################################################## 1. mapping
log "Stage 1/5: mapping reads to the genome with STAR and Bowtie2"

"$PYTHON" "${CODEBASE}/scripts/run_star.py" \
  --sample_id "$SAMPLE_ID" \
  --fastq "$FASTQ_READS" \
  --output_dir "$OUTPUT_DIR" \
  --genome_index "$GENOME_STAR_INDEX" \
  --threads "$THREADS"

for target in "genomic:${GENOME_BOWTIE_INDEX}" "transcriptomic:${TRANSCRIPTOME_INDEX}"; do
  "$PYTHON" "${CODEBASE}/scripts/run_bowtie.py" \
    --sample_id "$SAMPLE_ID" \
    --fastq "$FASTQ_READS" \
    --output_dir "$OUTPUT_DIR" \
    --reference_index "${target#*:}" \
    --logic_name "${target%%:*}" \
    --threads "$THREADS"
done

########################################################################## 2. discovery
log "Stage 2/5: screening mapped reads for putative backsplice junctions"

java "${JAVA_OPTS[@]}" -cp "$CLASSPATH" \
  bio.igm.utils.discovery.ProcessShuffledCoordinates \
  "$WORKING_DIR" "$MAX_GENOMIC_SPAN" "$SEGMENT_SIZE" "$MIN_GENOMIC_SPAN" "$LEGACY_ARG"

if [[ ! -s "${WORKING_DIR}/putative_structures.txt" ]]; then
  warn "No putative backsplice junctions were found; nothing further to do"
  warn "Run finished at $(date) without generating final output"
  exit 0
fi

########################################################################## 3. constructs
log "Stage 3/5: generating sequence constructs and building indexes"

java "${JAVA_OPTS[@]}" -cp "$CLASSPATH" \
  bio.igm.utils.discovery.GenerateSequenceConstructsGenome \
  "$WORKING_DIR" "$GENOME_FASTA" "$LEGACY_ARG" "$NORMALISE_ARG"

[[ -s "${WORKING_DIR}/Constructs.fa" ]] \
  || die "No backsplice constructs were generated. Check that the genome FASTA sequence names match the STAR index."
[[ -s "${WORKING_DIR}/Can.fa" ]] \
  || die "No canonical junction constructs were generated; junction-per-million values cannot be computed."

bowtie2-build --threads "$THREADS" "${WORKING_DIR}/Can.fa" "${WORKING_DIR}/canonical" \
  >> "${WORKING_DIR}/bowtie-build.log" 2>&1
bowtie2-build --threads "$THREADS" "${WORKING_DIR}/Constructs.fa" "${WORKING_DIR}/ptes" \
  >> "${WORKING_DIR}/bowtie-build.log" 2>&1

########################################################################## 4. re-mapping
log "Stage 4/5: re-mapping reads to the candidate junctions"

for target in "ptes:${WORKING_DIR}/ptes" "canonical:${WORKING_DIR}/canonical"; do
  "$PYTHON" "${CODEBASE}/scripts/run_bowtie.py" \
    --sample_id "$SAMPLE_ID" \
    --fastq "$FASTQ_READS" \
    --output_dir "$OUTPUT_DIR" \
    --reference_index "${target#*:}" \
    --logic_name "${target%%:*}" \
    --threads "$THREADS"
done

########################################################################## 5. filtering
log "Stage 5/5: filtering potential false positive predictions"

java "${JAVA_OPTS[@]}" -cp "$CLASSPATH" \
  bio.igm.utils.filter.PipelineFilter \
  "$WORKING_DIR" "$JSPAN" "$PID" "${FILTER_ARGS[@]}" "$LEGACY_ARG"

########################################################################## reporting
if [[ ! -f "${WORKING_DIR}/pf-structures.bed" ]]; then
  die "Filtering finished without producing pf-structures.bed; see ${WORKING_DIR}/run.log"
fi

CIRC_STRUCTURES=$(wc -l < "${WORKING_DIR}/pf-structures.bed" | tr -d ' ')
CIRC_READS=$(wc -l < "${WORKING_DIR}/pf-supporting-reads.tab" | tr -d ' ')
CJUNCS_READS=$(awk '{ sum += $5 } END { print sum + 0 }' "${WORKING_DIR}/pf-flanking-canonical-junctions.bed")
TJR=$(( CIRC_READS + CJUNCS_READS ))

log "Identified circRNAs:        $CIRC_STRUCTURES"
log "circRNA supporting reads:   $CIRC_READS"
log "Canonical junction reads:   $CJUNCS_READS"

if (( TJR > 0 )); then
  awk -v X="$TJR" 'BEGIN { OFS = "\t" } { print $1":"$2"-"$3":"$6, $5, ($5 / X) * 1000000 }' \
    "${WORKING_DIR}/pf-structures.bed" > "${WORKING_DIR}/${SAMPLE_ID}_jpms.tsv"
  log "Junctions per million written to ${WORKING_DIR}/${SAMPLE_ID}_jpms.tsv"
else
  warn "No junction-spanning reads were found; skipping junction-per-million normalisation"
  : > "${WORKING_DIR}/${SAMPLE_ID}_jpms.tsv"
fi

########################################################################## cleanup
if $KEEP_INTERMEDIATES; then
  log "Keeping intermediate files (-k)"
else
  log "Removing intermediate alignment, construct and index files"
  find "$WORKING_DIR" -maxdepth 1 -type f \
    \( -name '*.ebwt' -o -name '*.bt2' -o -name '*.bt2l' \
       -o -name '*.fa' -o -name '*.sam' -o -name '*.bam' -o -name '*.bam.bai' \) \
    -delete
fi

log "Finished successfully at $(date)"
