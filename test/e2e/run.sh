#!/usr/bin/env bash
#
# End-to-end test against the published test bundle.
#
# Unlike test/smoke, this runs the real aligners: STAR nominates the candidates, Bowtie2
# does four passes, and the result is compared against the junctions the full-depth run
# called on the same chromosome. It is the test that would have caught a broken STAR
# profile or an index built with the wrong parameters.
#
#   bash test/e2e/run.sh                      # the fixture committed alongside this script
#   bash test/e2e/run.sh --zenodo             # the larger chr17 bundle from Zenodo
#   bash test/e2e/run.sh --record 22929133     # a specific record or version DOI
#   bash test/e2e/run.sh --bundle DIR         # a bundle already on disk
#   bash test/e2e/run.sh --keep               # keep the working directory
#
# The committed fixture is chromosome 17 with everything outside the expected junctions
# masked to N. Masking rather than slicing keeps every coordinate identical to GRCm38, so
# the expected output is real chromosome 17 coordinates, while the N runs compress to a few
# megabytes -- small enough to live in the repository and need no download.
#
# Needs STAR, Bowtie2, samtools, Python 3 and Java on PATH, about 8 GB of RAM and 4 GB of
# scratch. Takes a few minutes on eight cores.
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
PYTHON="${PYTHON:-python3}"

# The Zenodo record holding the larger chromosome 17 bundle. The committed fixture is used
# by default and needs no download; this is only consulted when --record is given or
# PFV2_TEST_RECORD is set. Until the record is published its files can still change, so it
# is not the default source for a test.
ZENODO_RECORD="${PFV2_TEST_RECORD:-}"
ZENODO_RECORD_DEFAULT="10.5281/zenodo.22929132"
SAMPLE="bl6-chr17"
THREADS="${THREADS:-8}"
BUNDLE=""
WORK=""
KEEP=false
ZENODO_LATEST=false
FIXTURE="${CODE}/test/e2e/fixture"

while [ $# -gt 0 ]; do
  case "$1" in
    --bundle) BUNDLE="$2"; shift 2 ;;
    --work)   WORK="$2"; shift 2 ;;
    --record) ZENODO_RECORD="${2:-$ZENODO_RECORD_DEFAULT}"; shift 2 ;;
    --zenodo) ZENODO_RECORD="$ZENODO_RECORD_DEFAULT"; ZENODO_LATEST=true; shift ;;
    --threads) THREADS="$2"; shift 2 ;;
    --keep)   KEEP=true; shift ;;
    -h|--help)
      sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
      exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

die() { printf 'E2E FAIL: %s\n' "$*" >&2; exit 1; }
step() { printf '>>> %s\n' "$*"; }

for tool in STAR bowtie2 bowtie2-build samtools java "$PYTHON"; do
  command -v "$tool" >/dev/null 2>&1 || die "$tool is not on PATH"
done
[ -f "${CODE}/PFv2.jar" ] || die "PFv2.jar not found; run 'make build' first"

WORK="${WORK:-$(mktemp -d)}"
mkdir -p "$WORK"
$KEEP || trap 'rm -rf "$WORK"' EXIT

# ------------------------------------------------------------------ bundle
if [ -z "$BUNDLE" ] && [ -z "$ZENODO_RECORD" ] && [ -d "$FIXTURE" ]; then
  step "Using the committed fixture at test/e2e/fixture"
  BUNDLE="$FIXTURE"
fi

if [ -z "$BUNDLE" ]; then
  [ -n "$ZENODO_RECORD" ] || die "no bundle given and no fixture found; \
pass --bundle DIR, or --record ID, or set PFV2_TEST_RECORD"
  step "Fetching the test bundle from Zenodo record ${ZENODO_RECORD}"
  BUNDLE="${WORK}/bundle"
  mkdir -p "$BUNDLE"
  LATEST_FLAG=()
  $ZENODO_LATEST && LATEST_FLAG=(--latest)
  "${CODE}/bin/ptesfinder" fetch-references "$ZENODO_RECORD" --dest "$BUNDLE" \
    "${LATEST_FLAG[@]}" --manifest "${BUNDLE}/zenodo-manifest.json"
  # The bundle may be published as a single archive or as loose files.
  for archive in "${BUNDLE}"/*.tar.gz; do
    [ -f "$archive" ] || continue
    step "Unpacking $(basename "$archive")"
    tar xzf "$archive" -C "$BUNDLE" --strip-components=1
  done
fi

for required in transcriptome.fa "${SAMPLE}.fq.gz" expected/pf-structures.bed; do
  [ -f "${BUNDLE}/${required}" ] || die "bundle is missing ${required} (looked in ${BUNDLE})"
done

if [ -f "${BUNDLE}/SHA256SUMS" ]; then
  step "Verifying bundle checksums"
  ( cd "$BUNDLE" && if command -v sha256sum >/dev/null 2>&1; then
      sha256sum -c --quiet SHA256SUMS
    else
      shasum -a 256 -c --status SHA256SUMS
    fi ) || die "bundle checksums do not match; re-download it"
fi

# PFv2 slices the genome FASTA directly, so a compressed reference is decompressed into the
# working directory rather than read in place.
GENOME="${BUNDLE}/genome.fa"
if [ ! -f "$GENOME" ]; then
  [ -f "${BUNDLE}/genome.fa.gz" ] || die "bundle has neither genome.fa nor genome.fa.gz"
  step "Decompressing the reference"
  GENOME="${WORK}/genome.fa"
  gzip -dc "${BUNDLE}/genome.fa.gz" > "$GENOME"
fi

# ------------------------------------------------------------------ indexes
# Built here rather than shipped: an aligner index is tied to the release that wrote it, and
# a 95 Mb reference indexes in about a minute.
REF="${WORK}/ref"
mkdir -p "$REF/star"
step "Building the STAR index"
# genomeSAindexNbases must be scaled down for a single chromosome; STAR's own guidance is
# min(14, log2(length)/2 - 1), and leaving it at the default wastes memory and time.
BASES=$(grep -v '^>' "$GENOME" | tr -d '\n' | wc -c | tr -d ' ')
SA_N=$("$PYTHON" -c "import math,sys; print(min(14, int(math.log2(int(sys.argv[1]))/2 - 1)))" "$BASES")
STAR --runMode genomeGenerate --runThreadN "$THREADS" \
  --genomeDir "$REF/star" --genomeFastaFiles "$GENOME" \
  --genomeSAindexNbases "$SA_N" --outFileNamePrefix "${WORK}/star_build_" \
  > "${WORK}/star_build.log" 2>&1 || die "STAR genomeGenerate failed; see ${WORK}/star_build.log"

step "Building the Bowtie2 indexes"
bowtie2-build --threads "$THREADS" "$GENOME" "$REF/genome" \
  > "${WORK}/bt2_genome.log" 2>&1 || die "bowtie2-build failed on the genome"
bowtie2-build --threads "$THREADS" "${BUNDLE}/transcriptome.fa" "$REF/transcriptome" \
  > "${WORK}/bt2_tx.log" 2>&1 || die "bowtie2-build failed on the transcriptome"

# ------------------------------------------------------------------ run
step "Running PTESFinder"
READ_LENGTH=$("$PYTHON" - "${BUNDLE}/${SAMPLE}.fq.gz" <<'PY'
import gzip, sys
from collections import Counter
lengths = Counter()
with gzip.open(sys.argv[1], "rt") as fh:
    for i, line in enumerate(fh):
        if i % 4 == 1:
            lengths[len(line.rstrip("\n"))] += 1
        if i > 4000:
            break
print(lengths.most_common(1)[0][0])
PY
)
step "Read length detected as ${READ_LENGTH}"

"${CODE}/bin/ptesfinder" run \
  -i "$SAMPLE" \
  -r "${BUNDLE}/${SAMPLE}.fq.gz" \
  -d "$WORK" \
  -g "$GENOME" \
  -S "$REF/star" \
  -b "$REF/genome" \
  -t "$REF/transcriptome" \
  -l "$READ_LENGTH" \
  -n "$THREADS" \
  -m "${JAVA_HEAP:-8G}" \
  -k

# ------------------------------------------------------------------ compare
RESULT="${WORK}/PF/${SAMPLE}/pf-structures.bed"
[ -f "$RESULT" ] || die "no pf-structures.bed was produced"

step "Comparing against the expected junctions"
"$PYTHON" - "$RESULT" "${BUNDLE}/expected/pf-structures.bed" <<'PY'
import sys

def load(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 6:
                out[(f[0], f[1], f[2], f[5])] = f
    return out

got, want = load(sys.argv[1]), load(sys.argv[2])
missing = sorted(set(want) - set(got))
extra = sorted(set(got) - set(want))
shared = set(want) & set(got)

print(f"  expected {len(want)} junction(s), produced {len(got)}, matched {len(shared)}")
for key in missing:
    print(f"  MISSING  {':'.join(key[:3])} {key[3]}")
for key in extra[:10]:
    print(f"  EXTRA    {':'.join(key[:3])} {key[3]} (reads {got[key][4]})")
if len(extra) > 10:
    print(f"  ... and {len(extra) - 10} more extra")

signal_mismatch = [k for k in shared if len(want[k]) > 6 and len(got[k]) > 6
                   and want[k][6] != got[k][6]]
for key in signal_mismatch:
    print(f"  SIGNAL   {':'.join(key[:3])} expected {want[key][6]} got {got[key][6]}")

# Every expected junction must be recovered, and the splice signal and strand must match.
# Extra junctions are reported but tolerated: the bundle's read set is a subset, so support
# can land differently at the margin.
if missing or signal_mismatch:
    raise SystemExit(1)
PY
rc=$?

if [ "$rc" -eq 0 ]; then
  echo ">>> E2E PASS"
else
  echo ">>> E2E FAIL: the expected junctions were not all recovered" >&2
fi
$KEEP && echo ">>> working directory kept at ${WORK}"
exit "$rc"
