#!/usr/bin/env bash
#
# Builds the aligner indexes PFv2 needs from a genome FASTA and a transcriptome FASTA.
#
# The published reference record carries the sequences and the Bowtie2 indexes. The STAR
# index is not always published, because it is tens of gigabytes that barely compress and it
# is tied to the STAR genome format version that wrote it — an index built here always
# matches the STAR you actually have.
#
#   bash ops/build_indexes.sh \
#     --genome /data/refs/genome.fa \
#     --transcriptome /data/refs/transcriptome.cdna.fa \
#     --read-length 150 \
#     --out /data/refs/indexes
#
# Budget roughly an hour and 32 GB of RAM for a mammalian STAR index; the Bowtie2 indexes
# take a fraction of that. Use --only to build one of them.
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

GENOME=""; TRANSCRIPTOME=""; OUT=""; READ_LENGTH=150; THREADS="${THREADS:-8}"
ONLY=""; GTF=""; RAM="${RAM:-32000000000}"

while [ $# -gt 0 ]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2 ;;
    --transcriptome) TRANSCRIPTOME="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --read-length) READ_LENGTH="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --ram) RAM="$2"; shift 2 ;;
    --only) ONLY="$2"; shift 2 ;;
    -h|--help) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
step() { printf '>>> %s\n' "$*"; }

[ -n "$OUT" ] || die "--out is required"
[ -n "$GENOME" ] && [ -f "$GENOME" ] || die "--genome must be an existing FASTA"
mkdir -p "$OUT"

want() { [ -z "$ONLY" ] || [ "$ONLY" = "$1" ]; }

if want star; then
  command -v STAR >/dev/null 2>&1 || die "STAR is not on PATH"
  mkdir -p "${OUT}/star"
  # sjdbOverhang is read length minus one, as STAR's own guidance has it; an index built for
  # a very different read length should be rebuilt rather than reused.
  OVERHANG=$((READ_LENGTH - 1))
  # For a small reference the default genomeSAindexNbases wastes memory and time; STAR's
  # formula is min(14, log2(length)/2 - 1).
  BASES=$(grep -v '^>' "$GENOME" | tr -d '\n' | wc -c | tr -d ' ')
  SA_N=$(python3 -c "import math,sys; print(min(14, max(4, int(math.log2(int(sys.argv[1]))/2 - 1))))" "$BASES")
  step "STAR index: ${BASES} bases, sjdbOverhang ${OVERHANG}, genomeSAindexNbases ${SA_N}"
  STAR --runMode genomeGenerate --runThreadN "$THREADS" \
    --genomeDir "${OUT}/star" \
    --genomeFastaFiles "$GENOME" \
    --genomeSAindexNbases "$SA_N" \
    --sjdbOverhang "$OVERHANG" \
    --limitGenomeGenerateRAM "$RAM" \
    ${GTF:+--sjdbGTFfile "$GTF"} \
    --outFileNamePrefix "${OUT}/star_build_" \
    || die "STAR genomeGenerate failed"
  step "STAR index written to ${OUT}/star"
fi

if want bowtie2-genome; then
  command -v bowtie2-build >/dev/null 2>&1 || die "bowtie2-build is not on PATH"
  step "Bowtie2 genome index"
  bowtie2-build --threads "$THREADS" "$GENOME" "${OUT}/genome" > "${OUT}/bowtie2-genome.log" 2>&1 \
    || die "bowtie2-build failed on the genome; see ${OUT}/bowtie2-genome.log"
fi

if want bowtie2-transcriptome; then
  [ -n "$TRANSCRIPTOME" ] && [ -f "$TRANSCRIPTOME" ] \
    || die "--transcriptome must be an existing FASTA to build that index"
  command -v bowtie2-build >/dev/null 2>&1 || die "bowtie2-build is not on PATH"
  step "Bowtie2 transcriptome index"
  bowtie2-build --threads "$THREADS" "$TRANSCRIPTOME" "${OUT}/transcriptome" \
    > "${OUT}/bowtie2-transcriptome.log" 2>&1 \
    || die "bowtie2-build failed on the transcriptome; see ${OUT}/bowtie2-transcriptome.log"
fi

step "Done. Verify before running:"
printf '    ptesfinder validate --genome %s --star-index %s/star --bowtie2-genome %s/genome\n' \
  "$GENOME" "$OUT" "$OUT"
