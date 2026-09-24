'''
Author: Osagie Izuogu

Description: Validates PFv2 inputs and references before a run starts.

             The check that matters most is contig-name concordance. STAR, Bowtie2 and
             the genome FASTA each carry their own copy of the sequence names, and
             nothing downstream notices when they disagree: PFv2 nominates junctions on
             names STAR reported, then looks those names up in the FASTA, finds nothing,
             and produces an empty result. "chr1" and "1" are not interchangeable, and a
             run that mixes them wastes hours before failing in a way that looks like a
             biological result.

             Exits 0 when everything needed is present and consistent, 1 on any error.
             Warnings never fail the run.

Date: 09/2026
'''

from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
from collections import Counter

# STAR writes these into every index directory; a missing one means the index is
# incomplete rather than absent, which is harder to spot.
STAR_INDEX_FILES = (
    "chrLength.txt", "chrNameLength.txt", "chrName.txt", "chrStart.txt",
    "Genome", "genomeParameters.txt", "SA", "SAindex",
)

FASTQ_SAMPLE_RECORDS = 2000


class Report:
    """Collects findings so a single run reports every problem, not just the first."""

    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.facts: dict[str, object] = {}

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def fact(self, key: str, value) -> None:
        self.facts[key] = value

    def ok(self) -> bool:
        return not self.errors

    def render(self, stream=sys.stdout) -> None:
        for key in sorted(self.facts):
            print(f"  {key:<34} {self.facts[key]}", file=stream)
        for w in self.warnings:
            print(f"  WARN  {w}", file=stream)
        for e in self.errors:
            print(f"  ERROR {e}", file=stream)
        verdict = "PASS" if self.ok() else "FAIL"
        print(f"  {'-' * 60}", file=stream)
        print(f"  validation {verdict}: {len(self.errors)} error(s), "
              f"{len(self.warnings)} warning(s)", file=stream)


def readable_file(path: str, label: str, report: Report) -> bool:
    if not os.path.isfile(path):
        report.error(f"{label} not found: {path}")
        return False
    if not os.access(path, os.R_OK):
        report.error(f"{label} is not readable: {path}")
        return False
    if os.path.getsize(path) == 0:
        report.error(f"{label} is empty: {path}")
        return False
    return True


def check_fastq(paths: list[str], declared_length: int | None, report: Report) -> None:
    """Sanity-check FASTQ structure and compare the read length against --read-length."""
    lengths: Counter[int] = Counter()
    total = 0

    for path in paths:
        if not readable_file(path, "FASTQ", report):
            continue
        opener = gzip.open if path.endswith(".gz") else open
        try:
            with opener(path, "rt") as fh:
                for index, line in enumerate(fh):
                    position = index % 4
                    if position == 0 and not line.startswith("@"):
                        report.error(f"{path}: record {index // 4 + 1} does not start with '@'; "
                                     "the file is not 4-line FASTQ")
                        break
                    if position == 2 and not line.startswith("+"):
                        report.error(f"{path}: record {index // 4 + 1} has no '+' separator")
                        break
                    if position == 1:
                        lengths[len(line.rstrip("\n"))] += 1
                        total += 1
                        if total >= FASTQ_SAMPLE_RECORDS:
                            break
        except (OSError, EOFError) as exc:
            report.error(f"{path}: cannot be read as {'gzip' if path.endswith('.gz') else 'text'}: {exc}")
            continue

    if not lengths:
        if not report.errors:
            report.error("no reads could be read from the FASTQ input")
        return

    modal, count = lengths.most_common(1)[0]
    report.fact("reads sampled", total)
    report.fact("modal read length", modal)
    if len(lengths) > 1:
        report.fact("read length range", f"{min(lengths)}-{max(lengths)}")
        if count / total < 0.5:
            report.warn(f"read lengths are highly variable ({len(lengths)} distinct in "
                        f"{total} reads); PFv2 sizes its constructs from one length")

    if declared_length is not None and declared_length != modal:
        # The construct arm is derived from this number, so a wrong value silently changes
        # how much overhang a read needs on the short side of the junction.
        report.warn(f"-l/--read-length is {declared_length} but the modal read length is "
                    f"{modal}; construct arms will be sized for {declared_length}")


def fasta_names(path: str, limit: int | None = None) -> list[str]:
    """Sequence names, taken as the first whitespace token of each header."""
    names = []
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                names.append(line[1:].split()[0] if len(line) > 1 else "")
                if limit and len(names) >= limit:
                    break
    return names


def check_genome(path: str, report: Report) -> list[str]:
    if not readable_file(path, "Genome FASTA", report):
        return []
    if path.endswith(".gz"):
        # STAR and Bowtie2 accept gzip, but PFv2 slices the FASTA directly.
        report.error(f"Genome FASTA must be uncompressed for PFv2: {path}")
        return []

    names = fasta_names(path)
    if not names:
        report.error(f"Genome FASTA contains no records: {path}")
        return []
    if "" in names:
        report.error(f"Genome FASTA has a record with an empty name: {path}")

    duplicates = [n for n, c in Counter(names).items() if c > 1]
    if duplicates:
        report.error(f"Genome FASTA has duplicate sequence names: {duplicates[:5]}")

    report.fact("genome sequences", len(names))
    if not os.path.isfile(path + ".fai"):
        report.warn(f"no FASTA index at {path}.fai; 'samtools faidx' makes other tools faster "
                    "(PFv2 itself does not need it)")
    return names


def check_star_index(path: str, report: Report) -> list[str]:
    if not os.path.isdir(path):
        report.error(f"STAR index directory not found: {path}")
        return []

    missing = [f for f in STAR_INDEX_FILES if not os.path.isfile(os.path.join(path, f))]
    if missing:
        report.error(f"STAR index at {path} is incomplete; missing: {', '.join(missing)}")

    chrname = os.path.join(path, "chrName.txt")
    if not os.path.isfile(chrname):
        return []
    with open(chrname) as fh:
        names = [line.strip() for line in fh if line.strip()]
    report.fact("STAR index sequences", len(names))
    return names


def check_bowtie2_index(prefix: str, label: str, report: Report) -> list[str]:
    small = prefix + ".1.bt2"
    large = prefix + ".1.bt2l"
    if not (os.path.isfile(small) or os.path.isfile(large)):
        report.error(f"{label} Bowtie2 index not found: expected {small} or {large}")
        return []

    if shutil.which("bowtie2-inspect") is None:
        report.warn(f"bowtie2-inspect not on PATH; cannot verify {label} index sequence names")
        return []

    try:
        result = subprocess.run(["bowtie2-inspect", "-n", prefix],
                                capture_output=True, text=True, timeout=300)
    except (OSError, subprocess.TimeoutExpired) as exc:
        report.warn(f"could not inspect the {label} Bowtie2 index: {exc}")
        return []
    if result.returncode != 0:
        report.error(f"{label} Bowtie2 index at {prefix} could not be read by bowtie2-inspect: "
                     f"{result.stderr.strip().splitlines()[-1] if result.stderr.strip() else 'unknown error'}")
        return []

    names = [line.strip() for line in result.stdout.splitlines() if line.strip()]
    report.fact(f"{label} Bowtie2 index sequences", len(names))
    return names


def check_concordance(genome: list[str], star: list[str], bowtie: list[str], report: Report) -> None:
    """Compare the sequence-name sets each component carries.

    The genome FASTA is the reference point, because PFv2 resolves every junction against
    it. A name STAR can report but the FASTA does not contain is fatal: every junction on
    that contig is silently lost.
    """
    if not genome:
        return

    gset = set(genome)

    if star:
        orphans = sorted(set(star) - gset)
        if orphans:
            hint = ""
            stripped = {n[3:] for n in gset if n.startswith("chr")}
            if set(star) & stripped or {f"chr{n}" for n in star} & gset:
                hint = (" The two look like the same assembly under different naming "
                        "conventions ('chr1' versus '1').")
            report.error(f"{len(orphans)} sequence(s) in the STAR index are absent from the "
                         f"genome FASTA, e.g. {orphans[:5]}. Junctions called on them would be "
                         f"dropped without a usable error.{hint}")
        extra = len(gset - set(star))
        if extra:
            report.warn(f"{extra} sequence(s) in the genome FASTA are not in the STAR index; "
                        "no junction can be called on them")

    if bowtie:
        orphans = sorted(set(bowtie) - gset)
        if orphans:
            report.error(f"{len(orphans)} sequence(s) in the genome Bowtie2 index are absent "
                         f"from the genome FASTA, e.g. {orphans[:5]}; the genomic filter would "
                         "be comparing against a different assembly")

    if star and bowtie and set(star) != set(bowtie):
        only_star = len(set(star) - set(bowtie))
        only_bowtie = len(set(bowtie) - set(star))
        report.warn(f"the STAR and genome Bowtie2 indexes do not cover the same sequences "
                    f"({only_star} only in STAR, {only_bowtie} only in Bowtie2)")


def check_transcriptome_index(prefix: str, report: Report) -> None:
    names = check_bowtie2_index(prefix, "Transcriptome", report)
    if names and all(n.startswith("chr") for n in names[: min(len(names), 25)]):
        # A transcriptome index built over chromosomes makes the transcriptomic filter a
        # duplicate of the genomic one, which silently weakens the filter set.
        report.warn("the transcriptome Bowtie2 index looks like a genome index; the "
                    "transcriptomic filter would duplicate the genomic one")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Validate PFv2 inputs and references before a run.")
    p.add_argument("--fastq", nargs="+", help="sequence reads in FASTQ format")
    p.add_argument("--genome", help="genome reference in FASTA format")
    p.add_argument("--star-index", help="pre-built STAR genome index directory")
    p.add_argument("--bowtie2-genome", help="genome Bowtie2 index prefix")
    p.add_argument("--bowtie2-transcriptome", help="transcriptome Bowtie2 index prefix")
    p.add_argument("--read-length", type=int, default=None,
                   help="the -l value PFv2 will be given, checked against the reads")
    p.add_argument("--json", dest="json_out", default=None,
                   help="also write the report as JSON to this path")
    p.add_argument("--warnings-as-errors", action="store_true",
                   help="exit non-zero on warnings too")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    report = Report()

    print("PFv2 input validation")

    if args.fastq:
        check_fastq(args.fastq, args.read_length, report)

    genome_names = check_genome(args.genome, report) if args.genome else []
    star_names = check_star_index(args.star_index, report) if args.star_index else []
    bowtie_names = (check_bowtie2_index(args.bowtie2_genome, "Genome", report)
                    if args.bowtie2_genome else [])
    if args.bowtie2_transcriptome:
        check_transcriptome_index(args.bowtie2_transcriptome, report)

    check_concordance(genome_names, star_names, bowtie_names, report)

    report.render()

    if args.json_out:
        with open(args.json_out, "w") as fh:
            json.dump({
                "ok": report.ok(),
                "facts": report.facts,
                "warnings": report.warnings,
                "errors": report.errors,
            }, fh, indent=2)
            fh.write("\n")

    if not report.ok():
        return 1
    if args.warnings_as_errors and report.warnings:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
