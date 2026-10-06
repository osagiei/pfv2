#!/usr/bin/env bash
#
# Checks that STAR, Bowtie2 and samtools actually align, rather than merely being installed.
#
# The Java smoke test uses fixtures that stand in for aligner output, so it passes on a build
# where no aligner can run at all. That is not hypothetical: a binary for the wrong
# architecture, or one using instructions the host does not have, is on PATH and dies on its
# first instruction. `ptesfinder version` catches the dead-on-arrival case; this catches a
# tool that reports a version and then fails to do the work.
#
# Deliberately not part of the default selftest, and so not part of the image build: a
# multi-arch build runs the non-native stage under QEMU, where the aligners cannot be
# expected to work. It is run against the published image on native hardware instead.
#
#   bash test/smoke/aligners.sh [work_dir]
#
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

WORK="${1:-$(mktemp -d)}"
PYTHON="${PYTHON:-python3}"
mkdir -p "$WORK"
cd "$WORK"

fail=0
check() {
  if [ "$2" = "$3" ] || { [ "${4:-eq}" = "ge" ] && [ "$3" -ge "$2" ] 2>/dev/null; }; then
    printf '  ok    %s (%s)\n' "$1" "$3"
  else
    printf '  FAIL  %s\n        wanted %s %s, got %s\n' "$1" "${4:-eq}" "$2" "$3"
    fail=1
  fi
}

# A 20 kb random reference and 10 reads taken straight out of it. Random sequence is used so
# every read has exactly one place to go; a repeat would make the expected counts arguable.
"$PYTHON" - <<'PY'
import random
random.seed(1)
g = "".join(random.choice("ACGT") for _ in range(20000))
with open("genome.fa", "w") as f:
    f.write(">chrT\n")
    for i in range(0, len(g), 60):
        f.write(g[i:i + 60] + "\n")
with open("reads.fq", "w") as f:
    for n, start in enumerate(range(1000, 11000, 1000)):
        s = g[start:start + 100]
        f.write("@r%d\n%s\n+\n%s\n" % (n, s, "I" * len(s)))
PY

echo ">>> Bowtie2: build an index and align"
bowtie2-build -q genome.fa gidx > bowtie2-build.log 2>&1 \
  || { echo "  FAIL  bowtie2-build did not complete"; cat bowtie2-build.log; exit 1; }
bowtie2 -x gidx -U reads.fq -S bt2.sam --no-unal > bowtie2.log 2>&1 \
  || { echo "  FAIL  bowtie2 did not complete"; cat bowtie2.log; exit 1; }
check "bowtie2 aligns all 10 reads" 10 "$(grep -vc '^@' bt2.sam || true)"
check "bowtie2 reports them as unique" 10 \
  "$(awk '!/^@/ && /AS:i:/ && !/XS:i:/' bt2.sam | wc -l | tr -d ' ')"

echo ">>> samtools: read the alignments back"
check "samtools counts the same records" 10 "$(samtools view -c bt2.sam | tr -d ' ')"
samtools sort -o bt2.bam bt2.sam > samtools.log 2>&1 \
  || { echo "  FAIL  samtools sort did not complete"; cat samtools.log; exit 1; }
samtools index bt2.bam && echo "  ok    samtools sort and index"

echo ">>> STAR: generate an index and align"
# genomeSAindexNbases is min(14, log2(20000)/2 - 1) for a reference this small; STAR warns
# and performs badly with the default.
STAR --runMode genomeGenerate --genomeDir star --genomeFastaFiles genome.fa \
  --genomeSAindexNbases 6 --outFileNamePrefix sb_ > star-build.log 2>&1 \
  || { echo "  FAIL  STAR genomeGenerate did not complete"; tail -20 star-build.log; exit 1; }
STAR --genomeDir star --readFilesIn reads.fq --outSAMtype SAM \
  --outFileNamePrefix sa_ > star-align.log 2>&1 \
  || { echo "  FAIL  STAR alignment did not complete"; tail -20 star-align.log; exit 1; }
check "STAR maps all 10 reads uniquely" 10 \
  "$(awk -F'\t' '/Uniquely mapped reads number/ {gsub(/[ \t]/, "", $2); print $2}' sa_Log.final.out)"

# Chimeric detection is the one STAR feature PFv2 depends on, and it is a separate code path
# from ordinary alignment, so it is exercised rather than assumed.
echo ">>> STAR: chimeric detection, which is the feature PFv2 relies on"
STAR --genomeDir star --readFilesIn reads.fq --outSAMtype None \
  --chimSegmentMin 20 --chimOutType Junctions --outFileNamePrefix sc_ \
  > star-chim.log 2>&1 \
  || { echo "  FAIL  STAR did not run with chimeric detection on"; tail -20 star-chim.log; exit 1; }
[ -f sc_Chimeric.out.junction ] \
  && echo "  ok    STAR wrote Chimeric.out.junction" \
  || { echo "  FAIL  STAR produced no Chimeric.out.junction"; fail=1; }

if [ "$fail" -eq 0 ]; then
  echo ">>> ALIGNERS PASS"
else
  echo ">>> ALIGNERS FAIL" >&2
fi
exit "$fail"
