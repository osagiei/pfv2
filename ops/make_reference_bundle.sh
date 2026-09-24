#!/usr/bin/env bash
#
# Packages a prebuilt PFv2 reference set into upload-ready parts.
#
# Aligner indexes are the expensive part of getting started: a mouse STAR index takes about
# an hour and ~32 GB of RAM to build. Publishing the prebuilt set removes that, and makes
# results comparable between labs because everyone loads byte-identical indexes.
#
# Two things the parts record, because getting either wrong wastes hours:
#
#   * The sequence naming. An Ensembl reference names contigs 1..19,X,Y with no "chr"
#     prefix; a UCSC mm10 reference names them chr1..chr19. They are not interchangeable,
#     and mixing them produces an empty result rather than an error. The manifest states
#     which convention the bundle uses and lists the first few names.
#   * The aligner that built the index. STAR stores a genome format version and refuses an
#     index it cannot read, so an index is only usable within a range of STAR releases. The
#     manifest records the versions, and the README pairs the record with an image tag.
#
# Run it on the machine holding the reference: it packages what is already on disk and
# downloads nothing. The paths below are illustrative -- substitute that machine's own.
#
#   bash ops/make_reference_bundle.sh \
#     --genome  /refs/genome.fa \
#     --cdna    /refs/transcriptome.cdna.fa \
#     --star    /refs/indexes/star/<hash> \
#     --bowtie2-genome /refs/indexes/bowtie2/<hash> \
#     --bowtie2-transcriptome /refs/indexes/bowtie2-tx \
#     --name GRCm38-ensembl102 --out /scratch/refbundle
#
# --out needs free space roughly equal to the reference, since the parts are written
# alongside it: about 27 GB for a mouse set whose STAR index alone is 26 GB. Archives above
# --max-part-size (default 2G) are split, because Zenodo's gateway returns 502 on a large
# single upload; `ptesfinder fetch-references` reassembles them.
#
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

GENOME=""; CDNA=""; STAR_DIR=""; BT2_GENOME=""; BT2_TX=""; GTF=""
NAME="reference"; OUT=""; THREADS="${THREADS:-8}"
# Zenodo's gateway returns 502 on a single upload of tens of gigabytes, so archives above
# this are split and reassembled on fetch. Measured: a 3.37 GB PUT succeeded and a 4 GB one
# returned 502, so the proxy gives up somewhere just under 4 GB; 2G leaves real headroom.
MAX_PART="${MAX_PART:-2G}"

while [ $# -gt 0 ]; do
  case "$1" in
    --genome) GENOME="$2"; shift 2 ;;
    --cdna) CDNA="$2"; shift 2 ;;
    --gtf) GTF="$2"; shift 2 ;;
    --star) STAR_DIR="$2"; shift 2 ;;
    --bowtie2-genome) BT2_GENOME="$2"; shift 2 ;;
    --bowtie2-transcriptome) BT2_TX="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --out) OUT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --max-part-size) MAX_PART="$2"; shift 2 ;;
    -h|--help) sed -n '2,30p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
step() { printf '>>> %s\n' "$*"; }

[ -n "$OUT" ] || die "--out is required"
[ -n "$GENOME" ] || die "--genome is required"
[ -f "$GENOME" ] || die "--genome: no such file: ${GENOME}
       This packages a reference that already exists; it does not download one. Run it on the
       machine holding the reference, with that machine's paths."
[ -n "$STAR_DIR" ] || die "--star is required"
[ -d "$STAR_DIR" ] || die "--star: no such directory: ${STAR_DIR}"

# An incomplete STAR index is the failure this checks for: a half-built one leaves a
# .SA.*.tmp behind and no chrName.txt, and nothing downstream notices until a run fails.
for f in chrLength.txt chrNameLength.txt chrName.txt chrStart.txt Genome SA SAindex genomeParameters.txt; do
  [ -f "${STAR_DIR}/${f}" ] || die "STAR index is incomplete: ${STAR_DIR}/${f} is missing"
done
if find "$STAR_DIR" -maxdepth 1 -name '*tmp*' -print -quit | grep -q .; then
  die "STAR index at ${STAR_DIR} has leftover temporary files; it was not built to completion"
fi

mkdir -p "$OUT"
COMPRESS="gzip"
command -v pigz >/dev/null 2>&1 && COMPRESS="pigz -p ${THREADS}"

pack() {
  local label="$1" src="$2"
  local out="${OUT}/${NAME}-${label}.tar.gz"
  step "Packing ${label}"
  # Archived from the parent so the member paths are relative and predictable.
  tar -C "$(dirname "$src")" -cf - "$(basename "$src")" | $COMPRESS > "$out"
  printf '    %s  %s\n' "$(du -h "$out" | cut -f1)" "$(basename "$out")"
  split_if_large "$out"
}

# Splits an archive that is too large to upload in one request. The whole-file checksum is
# recorded first, so a consumer can verify the reassembled archive and not merely the parts.
split_if_large() {
  local file="$1"
  local limit bytes
  limit=$(numfmt --from=iec "$MAX_PART" 2>/dev/null || echo 4294967296)
  bytes=$(wc -c < "$file" | tr -d ' ')
  [ "$bytes" -gt "$limit" ] || return 0

  step "Splitting $(basename "$file") into ${MAX_PART} parts (too large for one upload)"
  local whole
  whole=$(sha256_of "$file")
  printf '%s  %s\n' "$whole" "$(basename "$file")" >> "${OUT}/WHOLE-SHA256SUMS"
  split -b "$MAX_PART" -d -a 2 "$file" "${file}.part"
  rm -f "$file"
  local n
  n=$(find "$(dirname "$file")" -maxdepth 1 -name "$(basename "$file").part*" | wc -l | tr -d ' ')
  printf '    %s part(s)\n' "$n"
}

sha256_of() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

step "Contig naming check"
FIRST_NAMES=$(head -5 "${STAR_DIR}/chrName.txt" | tr '\n' ' ')
if head -1 "${STAR_DIR}/chrName.txt" | grep -q '^chr'; then
  CONVENTION="ucsc"
else
  CONVENTION="ensembl"
fi
printf '    convention: %s (first names: %s)\n' "$CONVENTION" "$FIRST_NAMES"

STAR_BUILD_VERSION=$(grep -m1 '^versionGenome' "${STAR_DIR}/genomeParameters.txt" | awk '{print $2}' || echo unknown)
SJDB_OVERHANG=$(grep -m1 '^sjdbOverhang' "${STAR_DIR}/genomeParameters.txt" | awk '{print $2}' || echo unknown)
STAR_RUNTIME_VERSION=$(command -v STAR >/dev/null 2>&1 && STAR --version 2>/dev/null || echo "not installed")
BT2_VERSION=$(command -v bowtie2 >/dev/null 2>&1 && bowtie2 --version 2>/dev/null | head -1 | awk '{print $NF}' || echo "not installed")
printf '    STAR genome format %s, sjdbOverhang %s\n' "$STAR_BUILD_VERSION" "$SJDB_OVERHANG"

step "Packing the sequence files"
mkdir -p "${OUT}/seq"
cp "$GENOME" "${OUT}/seq/genome.fa"
[ -n "$CDNA" ] && [ -f "$CDNA" ] && cp "$CDNA" "${OUT}/seq/transcriptome.cdna.fa"
[ -n "$GTF" ] && [ -f "$GTF" ] && cp "$GTF" "${OUT}/seq/annotation.gtf"
pack sequence "${OUT}/seq"
rm -rf "${OUT}/seq"

pack star-index "$STAR_DIR"
[ -n "$BT2_GENOME" ] && [ -d "$BT2_GENOME" ] && pack bowtie2-genome "$BT2_GENOME"
[ -n "$BT2_TX" ] && [ -d "$BT2_TX" ] && pack bowtie2-transcriptome "$BT2_TX"

step "Checksums"
( cd "$OUT" && if command -v sha256sum >/dev/null 2>&1; then
    sha256sum ./*.tar.gz > SHA256SUMS
  else
    shasum -a 256 ./*.tar.gz > SHA256SUMS
  fi )

TOTAL=$(du -sh "$OUT" | cut -f1)
cat > "${OUT}/MANIFEST.json" <<JSON
{
  "bundle": "pfv2-reference",
  "name": "${NAME}",
  "built": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "contig_naming": "${CONVENTION}",
  "first_contigs": "$(echo "$FIRST_NAMES" | sed 's/ *$//')",
  "contigs": $(wc -l < "${STAR_DIR}/chrName.txt" | tr -d ' '),
  "star": {
    "genome_format_version": "${STAR_BUILD_VERSION}",
    "sjdb_overhang": "${SJDB_OVERHANG}",
    "built_with": "${STAR_RUNTIME_VERSION}",
    "note": "STAR refuses an index whose genome format it cannot read, so this index is usable only by STAR releases sharing format ${STAR_BUILD_VERSION}. sjdbOverhang ${SJDB_OVERHANG} suits reads one base longer than that."
  },
  "bowtie2": { "built_with": "${BT2_VERSION}" },
  "total_size": "${TOTAL}"
}
JSON

step "Done: ${OUT} (${TOTAL})"
printf '    %s\n' "$(ls "$OUT" | tr '\n' ' ')"
echo
echo "Upload with:  bash ops/zenodo_upload.sh --dir ${OUT} --title '<title>'"
