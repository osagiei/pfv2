"""Generate a tiny synthetic dataset that exercises every PFv2 stage.

The fixture is deliberately small and deterministic. It is not a biological test set: it
checks that the stages agree on coordinates, junction offsets, splice signals and filter
outcomes, which is what breaks when one of them is changed in isolation.
"""

from __future__ import annotations

import argparse
import random
from pathlib import Path

SEED = 7
CONTIG_LENGTH = 3000
READ_LENGTH = 100
SEGMENT = READ_LENGTH - 15  # matches PFv2.sh's arm length

# A plus-strand backsplice: donor at 2000 splices back to the acceptor at 1000.
DONOR, ACCEPTOR = 2000, 1000
# A canonical junction whose intron runs 1500..1600 (1-based, inclusive).
INTRON_START, INTRON_END = 1500, 1600


def build_contig() -> str:
    random.seed(SEED)
    bases = [random.choice("ACGT") for _ in range(CONTIG_LENGTH)]
    # Plant a GT..AG motif at the canonical intron edges so the construct passes the
    # splice signal test.
    bases[INTRON_START - 1], bases[INTRON_START] = "G", "T"
    bases[INTRON_END - 2], bases[INTRON_END - 1] = "A", "G"
    return "".join(bases)


def sam(name: str, target: str, pos: int, nm: int, md: str, score: int, length: int = READ_LENGTH) -> str:
    return "\t".join([
        name, "0", target, str(pos), "42", f"{length}M", "*", "0", "0",
        "A" * length, "I" * length, f"AS:i:{score}", "XN:i:0", f"NM:i:{nm}", f"MD:Z:{md}",
    ])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("outdir", type=Path)
    ap.add_argument("--stage", choices=["star", "sam"], required=True,
                    help="star writes the STAR fixtures; sam writes the alignment fixtures")
    args = ap.parse_args()

    work = args.outdir / "PF" / "S1"
    work.mkdir(parents=True, exist_ok=True)

    if args.stage == "star":
        contig = build_contig()
        fasta = args.outdir / "genome.fa"
        with open(fasta, "w") as fh:
            # A descriptive Ensembl-style header on the FIRST record: the pre-2.2.0 parser
            # dropped record one entirely and could not read a header with a description.
            fh.write(">chr1 dna:chromosome chromosome:TEST:1:1:3000:1\n")
            for i in range(0, len(contig), 60):
                fh.write(contig[i:i + 60] + "\n")
            fh.write(">chr2\n" + "A" * 500 + "\n")

        with open(work / "star_Chimeric.out.junction", "w") as fh:
            fh.write(f"chr1\t{DONOR}\t+\tchr1\t{ACCEPTOR}\t+\t1\t0\t0\tread1\t1900\t100M\t950\t100M\n")
            # Rejected: inter-chromosomal.
            fh.write(f"chr1\t{DONOR}\t+\tchr2\t{ACCEPTOR}\t+\t1\t0\t0\tread9\t1900\t100M\t950\t100M\n")
            fh.write("# trailing comment block written by STAR >= 2.7.1\n")

        with open(work / "star_SJ.out.tab", "w") as fh:
            fh.write(f"chr1\t{INTRON_START}\t{INTRON_END}\t1\t1\t0\t9\t0\t40\n")
            # Rejected: non-canonical motif, and no uniquely mapping read.
            fh.write("chr1\t2500\t2600\t1\t0\t0\t9\t0\t40\n")
            fh.write("chr1\t2700\t2800\t1\t1\t0\t0\t9\t40\n")
        print(f"wrote STAR fixtures under {work}")
        return

    # The construct names are read back from the FASTA the previous stage produced, so the
    # fixture cannot drift from the junction offset the pipeline actually assigns.
    ptes = _first_header(work / "Constructs.fa")
    canonical = _first_header(work / "Can.fa")

    with open(work / "ptes.sam", "w") as fh:
        fh.write("@HD\tVN:1.0\tSO:unsorted\n")
        # Spans the seam with generous overhang either side. Releases up to 2.1.0 rejected
        # every read from about POS 41 onwards.
        fh.write(sam("read1", ptes, 41, 1, "70A29", -6) + "\n")
        # The junction falls outside this read: not evidence, however well it matches.
        fh.write(sam("read2", ptes, 95, 0, "100", 0) + "\n")
        # The genome explains this one just as well.
        fh.write(sam("read3", ptes, 1, 0, "100", 0) + "\n")
        fh.write(sam("read4", ptes, 60, 0, "100", 0) + "\n")
        # The transcriptome explains this one better.
        fh.write(sam("read5", ptes, 41, 1, "70A29", -6) + "\n")
        # A real canonical junction explains this one.
        fh.write(sam("read6", ptes, 41, 2, "40A30A28", -12) + "\n")

    with open(work / "genomic.sam", "w") as fh:
        fh.write("@HD\tVN:1.0\n")
        fh.write(sam("read1", "chr1", 1900, 3, "40A30A27", -18) + "\n")
        fh.write(sam("read3", "chr1", 1900, 0, "100", 0) + "\n")

    with open(work / "transcriptomic.sam", "w") as fh:
        fh.write("@HD\tVN:1.0\n")
        fh.write(sam("read5", "NM_000001", 100, 0, "100", 0) + "\n")

    with open(work / "canonical.sam", "w") as fh:
        fh.write("@HD\tVN:1.0\n")
        # read6 spans the canonical seam better than it spans the backsplice seam.
        fh.write(sam("read6", canonical, 41, 0, "100", 0) + "\n")
        fh.write(sam("read7", canonical, 41, 0, "100", 0) + "\n")
        fh.write(sam("read8", canonical, 1, 0, "100", 0) + "\n")
    print(f"wrote alignment fixtures under {work}")


def _first_header(path: Path) -> str:
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                return line[1:].strip()
    raise SystemExit(f"no FASTA record in {path}")


if __name__ == "__main__":
    main()
