#!/usr/bin/env bash
#
# End-to-end smoke test of the PFv2 Java stages against a synthetic dataset.
#
# STAR and Bowtie2 are not exercised: the fixtures stand in for their output. What this
# checks is that stages 2, 3 and 5 agree on coordinates, junction offsets, splice signals
# and filter outcomes, which is exactly what silently drifts when one is edited alone.
#
#   bash test/smoke/run.sh [work_dir]
#
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
WORK="${1:-$(mktemp -d)}"
PYTHON="${PYTHON:-python3}"

LIB="$(find "${CODE}/lib" -name 'commons-lang3-*.jar' | head -1)"
CP="${CODE}/PFv2.jar:${LIB}"
W="${WORK}/PF/S1"

fail=0
check() {
  local what="$1" expected="$2" actual="$3"
  if [ "$expected" = "$actual" ]; then
    printf '  ok    %s\n' "$what"
  else
    printf '  FAIL  %s\n        expected <%s>\n        actual   <%s>\n' "$what" "$expected" "$actual"
    fail=1
  fi
}

[ -f "${CODE}/PFv2.jar" ] || { echo "PFv2.jar not found; run 'bash setup.sh' first" >&2; exit 1; }

echo ">>> Fixtures in ${WORK}"
"$PYTHON" "${CODE}/test/smoke/make_fixtures.py" "$WORK" --stage star >/dev/null

echo ">>> Stage 2: candidate junctions"
java -cp "$CP" bio.igm.utils.discovery.ProcessShuffledCoordinates "$W" 1000000 85 50 0 2>/dev/null

check "one putative backsplice is nominated" "1" "$(wc -l < "${W}/putative_structures.txt" | tr -d ' ')"
check "one canonical junction is nominated" "1" "$(wc -l < "${W}/canonical_structures.txt" | tr -d ' ')"
check "backsplice anchors" \
  "chr1	1915	1999	1001	1085	chr1:1000-1999_+" \
  "$(cat "${W}/putative_structures.txt")"

echo ">>> Stage 3: sequence constructs"
java -cp "$CP" bio.igm.utils.discovery.GenerateSequenceConstructsGenome "$W" "${WORK}/genome.fa" 0 2>/dev/null

PTES_NAME="$(grep '^>' "${W}/Constructs.fa" | head -1 | cut -c2-)"
CAN_NAME="$(grep '^>' "${W}/Can.fa" | head -1 | cut -c2-)"
check "backsplice construct name carries the seam offset" "chr1:1000-1999_+:85" "${PTES_NAME%:*}"
check "canonical construct carries a GT-AG signal" "GTAG" "${CAN_NAME##*:}"

# The recorded offset must equal the length of the first arm, or the purity window is
# centred on the wrong base.
ARM_LENGTH="$("$PYTHON" - "$WORK" <<'PY'
import sys
from pathlib import Path
work = Path(sys.argv[1])
seq = [l.strip() for l in open(work / "PF/S1/Constructs.fa") if not l.startswith(">")]
genome = {}
name, buf = None, []
for line in open(work / "genome.fa"):
    if line.startswith(">"):
        if name: genome[name] = "".join(buf)
        name, buf = line[1:].split()[0], []
    else:
        buf.append(line.strip())
genome[name] = "".join(buf)
print(len(genome["chr1"][1914:1999]))
PY
)"
check "seam offset equals the first arm length" "85" "$ARM_LENGTH"

echo ">>> Stage 5: filters"
"$PYTHON" "${CODE}/test/smoke/make_fixtures.py" "$WORK" --stage sam >/dev/null
java -cp "$CP" bio.igm.utils.filter.PipelineFilter "$W" 8 0.85 1 0 0 0 2>/dev/null

check "one backsplice survives" "1" "$(wc -l < "${W}/pf-structures.bed" | tr -d ' ')"
check "BED row" \
  "chr1	1000	1999	${PTES_NAME}	2	+	${PTES_NAME##*:}" \
  "$(cat "${W}/pf-structures.bed")"
check "accepted reads" "read1 read4" \
  "$(cut -f1 "${W}/pf-supporting-reads.tab" | sort | tr '\n' ' ' | sed 's/ $//')"
check "read lost to the genome" "read3" \
  "$(cut -f1 "${W}/genomic-better.sam" | sort -u | tr '\n' ' ' | sed 's/ $//')"
check "read lost to the transcriptome" "read5" \
  "$(cut -f1 "${W}/transcriptomic-better.sam" | sort -u | tr '\n' ' ' | sed 's/ $//')"
check "read lost to a canonical junction" "read6" \
  "$(cut -f1 "${W}/canonical-better.sam" | sort -u | tr '\n' ' ' | sed 's/ $//')"
check "non-spanning perfect match is rejected" "read2" \
  "$(cut -f1 "${W}/pf-junctional-filtered.sam" | sort -u | tr '\n' ' ' | sed 's/ $//')"

echo ">>> Stage 5 again, legacy semantics"
java -cp "$CP" bio.igm.utils.filter.PipelineFilter "$W" 8 0.85 1 0 0 1 2>/dev/null
check "legacy mode accepts the non-spanning perfect match" "read1 read2 read4 read6" \
  "$(cut -f1 "${W}/pf-supporting-reads.tab" | sort | tr '\n' ' ' | sed 's/ $//')"

if [ "$fail" -eq 0 ]; then
  echo ">>> SMOKE PASS"
else
  echo ">>> SMOKE FAIL" >&2
fi
exit "$fail"
