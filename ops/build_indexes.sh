#!/usr/bin/env bash
#
# Builds the aligner indexes PFv2 needs from a genome FASTA and a transcriptome.
#
# The transcriptome index is not optional: the transcriptomic filter is one of the three
# false-positive filters described in the PTESFinder paper, and PFv2.sh refuses to run
# without it. A read only counts as backsplice evidence if it aligns better to the junction
# construct than to a linear mRNA, so the index must be transcript sequences - exons
# concatenated per transcript - and not the genome.
#
# If you have no cDNA FASTA, give an annotation instead and one is derived here:
#   --gtf             GTF/GFF3, converted with gffread (the usual route; this is the same
#                     annotation STAR takes for --sjdbGTFfile)
#   --transcriptome-bed  UCSC BED12, converted with bedtools and PFv2's own
#                     MergeUCSCExonsToTranscript, which is the route PTESFinder v1 used
#
# The published reference record carries the sequences and the Bowtie2 indexes. The STAR
# index is not always published, because it is tens of gigabytes that barely compress and it
# is tied to the STAR genome format version that wrote it - an index built here always
# matches the STAR you actually have.
#
#   bash ops/build_indexes.sh \
#     --genome /data/refs/genome.fa \
#     --gtf /data/refs/annotation.gtf \
#     --read-length 150 \
#     --out /data/refs/indexes
#
# Budget roughly an hour and 32 GB of RAM for a mammalian STAR index; the Bowtie2 indexes
# take a fraction of that. Use --only to build one of them.
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

GENOME=""; TRANSCRIPTOME=""; OUT=""; READ_LENGTH=150; THREADS="${THREADS:-8}"
ONLY=""; GTF=""; TRANSCRIPTOME_BED=""; RAM="${RAM:-32000000000}"

while [ $# -gt 0 ]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2 ;;
    --transcriptome) TRANSCRIPTOME="$2"; shift 2 ;;
    --transcriptome-bed) TRANSCRIPTOME_BED="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --read-length) READ_LENGTH="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --ram) RAM="$2"; shift 2 ;;
    --only) ONLY="$2"; shift 2 ;;
    -h|--help) sed -n '2,33p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
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
  command -v bowtie2-build >/dev/null 2>&1 || die "bowtie2-build is not on PATH"

  # A transcriptome built over sequence names the genome does not use produces an index that
  # nothing aligns to, which silently turns the transcriptomic filter into a no-op rather
  # than failing. Checked before spending time on the build.
  # $2 is a snippet that reads the annotation on stdin and writes one sequence name per line.
  check_contig_names() {
    local annotation="$1" reader="$2" ann_names genome_names shared
    genome_names="$(grep '^>' "$GENOME" | sed 's/^>//; s/[[:space:]].*//' | sort -u)"
    ann_names="$(eval "$reader" < "$annotation" | sort -u)"
    [ -n "$ann_names" ] || die "no sequence names found in ${annotation}"
    shared="$(comm -12 <(printf '%s\n' "$genome_names") <(printf '%s\n' "$ann_names") | wc -l | tr -d ' ')"
    if [ "$shared" -eq 0 ]; then
      die "$(printf '%s\n' \
        "the annotation and the genome name their sequences differently, so the derived" \
        "transcriptome would be empty and the transcriptomic filter would pass everything." \
        "  annotation: $(printf '%s\n' "$ann_names" | head -3 | paste -sd' ' -)" \
        "  genome:     $(printf '%s\n' "$genome_names" | head -3 | paste -sd' ' -)" \
        "Use an annotation and genome from the same source and assembly release.")"
    fi
  }

  if [ -z "$TRANSCRIPTOME" ] && [ -n "$GTF" ]; then
    # gffread -w concatenates each transcript's exons in transcript orientation, which is
    # what the filter needs; -W would add exon coordinates to the headers and is not wanted.
    command -v gffread >/dev/null 2>&1 \
      || die "gffread is not on PATH; install it (conda install -c bioconda gffread), or pass --transcriptome with a cDNA FASTA"
    check_contig_names "$GTF" "grep -v '^#' | cut -f1"
    TRANSCRIPTOME="${OUT}/transcriptome.fa"
    step "Transcript FASTA from ${GTF} with gffread"
    gffread -w "$TRANSCRIPTOME" -g "$GENOME" "$GTF" > "${OUT}/gffread.log" 2>&1 \
      || die "gffread failed; see ${OUT}/gffread.log"
  elif [ -z "$TRANSCRIPTOME" ] && [ -n "$TRANSCRIPTOME_BED" ]; then
    # The BED12 route PTESFinder v1 used. v1 split BED12 into per-exon records, pulled each
    # exon stranded, then rejoined them with its own MergeUCSCExonsToTranscript, which parses
    # the exon number back out of the FASTA header and so only accepts RefSeq or knownGene
    # style ids. "getfasta -split" does the same join in one step for any transcript id.
    command -v bedtools >/dev/null 2>&1 \
      || die "bedtools is not on PATH; install it (conda install -c bioconda bedtools), or pass --transcriptome with a cDNA FASTA"
    check_contig_names "$TRANSCRIPTOME_BED" 'cut -f1'
    awk 'NF < 12 { print "BED12 is required (12 columns); line " NR " has " NF > "/dev/stderr"; exit 1 }' \
      "$TRANSCRIPTOME_BED" || die "${TRANSCRIPTOME_BED} is not BED12"
    TRANSCRIPTOME="${OUT}/transcriptome.fa"
    step "Transcript FASTA from ${TRANSCRIPTOME_BED} with bedtools"
    bedtools getfasta -fi "$GENOME" -bed "$TRANSCRIPTOME_BED" -split -s -nameOnly \
      -fo "${OUT}/transcriptome.raw.fa" || die "bedtools getfasta failed"
    # Sequence is uppercased because the genome's soft-masked regions would otherwise reach
    # Bowtie2 as lower case, which v1 also guarded against. Headers are left alone apart from
    # the "(+)" that -nameOnly -s appends, so reference names still match the annotation.
    awk '/^>/ { sub(/\([+-]\)$/, ""); print; next } { print toupper($0) }' \
      "${OUT}/transcriptome.raw.fa" > "$TRANSCRIPTOME"
    rm -f "${OUT}/transcriptome.raw.fa"
  fi

  [ -n "$TRANSCRIPTOME" ] && [ -f "$TRANSCRIPTOME" ] || die "$(printf '%s\n' \
    "the transcriptomic filter needs transcript sequences. Give one of:" \
    "  --transcriptome <cdna.fa>     a cDNA/transcript FASTA you already have" \
    "  --gtf <annotation.gtf>        derived here with gffread" \
    "  --transcriptome-bed <ann.bed> derived here the way PTESFinder v1 did")"

  COUNT=$(grep -c '^>' "$TRANSCRIPTOME" || true)
  step "Bowtie2 transcriptome index: ${COUNT} transcript(s)"
  bowtie2-build --threads "$THREADS" "$TRANSCRIPTOME" "${OUT}/transcriptome" \
    > "${OUT}/bowtie2-transcriptome.log" 2>&1 \
    || die "bowtie2-build failed on the transcriptome; see ${OUT}/bowtie2-transcriptome.log"
fi

step "Done. Verify before running:"
printf '    ptesfinder validate --genome %s --star-index %s/star --bowtie2-genome %s/genome --bowtie2-transcriptome %s/transcriptome\n' \
  "$GENOME" "$OUT" "$OUT" "$OUT"
