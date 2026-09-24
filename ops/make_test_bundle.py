'''
Author: Osagie Izuogu

Description: Builds the small end-to-end test bundle from a full reference and a completed
             PFv2 run.

             The bundle is one chromosome, the transcripts on it, and the reads needed to
             rediscover the junctions PFv2 called there. It deliberately ships FASTA rather
             than prebuilt aligner indexes: an index is tied to the aligner release that
             wrote it, so shipping one turns an aligner upgrade into a silent failure,
             while a 95 Mb reference indexes in about a minute.

             Read selection has three parts, and all three matter:

             * every chimeric read STAR reported within the chromosome, so candidate
               nomination and the rejection paths are exercised, not just the calls;
             * the reads that supported the calls, so the expected output is reachable at
               all -- at full depth these are a handful of reads out of 187 million;
             * a background sample of reads mapped to the chromosome, so the genomic and
               transcriptomic filters and the canonical denominator have real work to do.

             Without the background the filters pass everything and the test proves nothing.

Date: 09/2026
'''

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import random
import subprocess
import sys
from datetime import datetime, timezone

CHIM_DONOR_CHR, CHIM_ACC_CHR, CHIM_READ = 0, 3, 9
FASTA_WIDTH = 60


def log(message: str) -> None:
    print(f"  {message}", flush=True)


def fail(message: str) -> None:
    print(f"ERROR {message}", file=sys.stderr)
    raise SystemExit(1)


def opener(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def read_contig(genome: str, contig: str) -> str:
    """Returns one contig's sequence, matching on the first whitespace token of the header."""
    parts: list[str] = []
    with opener(genome) as fh:
        keeping = False
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0] if len(line) > 1 else ""
                if keeping:
                    break
                keeping = name == contig
                continue
            if keeping:
                parts.append(line.strip())
    if not parts:
        fail(f"contig {contig} not found in {genome}")
    return "".join(parts)


def windows_from(structures: str, contig: str, flank: int) -> list[tuple[int, int]]:
    """Merged [start, end) windows around every expected junction on the contig.

    A window has to be wide enough to hold both construct arms plus the reads that align
    across them, and wide enough that the background sample inside it is not trivially
    small. Overlapping windows are merged so the masked sequence has no seams.
    """
    spans = []
    with open(structures) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 3 and f[0] == contig:
                spans.append((max(0, int(f[1]) - flank), int(f[2]) + flank))
    if not spans:
        return []
    spans.sort()
    merged = [list(spans[0])]
    for start, end in spans[1:]:
        if start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(a, b) for a, b in merged]


def write_contig(sequence: str, contig: str, out_path: str,
                 windows: list[tuple[int, int]] | None) -> tuple[int, int]:
    """Writes the contig, optionally masking everything outside the given windows to N.

    Masking rather than slicing keeps every coordinate identical to the source assembly, so
    the expected output is real chromosome 17 coordinates and can be read against an
    annotation. A run of N compresses to almost nothing, which is what makes the masked
    reference small enough to keep in the repository.
    """
    length = len(sequence)
    if windows:
        kept = 0
        buf = ["N"] * length
        for start, end in windows:
            start, end = max(0, start), min(length, end)
            buf[start:end] = sequence[start:end]
            kept += end - start
        sequence = "".join(buf)
    else:
        kept = length

    writer = gzip.open(out_path, "wt", compresslevel=9) if out_path.endswith(".gz") else open(out_path, "w")
    with writer as out:
        out.write(f">{contig}\n")
        for i in range(0, length, FASTA_WIDTH):
            out.write(sequence[i:i + FASTA_WIDTH] + "\n")
    return length, kept


def extract_transcripts(cdna: str, contig: str, out_path: str,
                        windows: list[tuple[int, int]] | None = None) -> int:
    """Keeps Ensembl cDNA records on the contig, and inside the windows when masking.

    Ensembl headers carry the locus as ``chromosome:<assembly>:<name>:<start>:<end>:<strand>``,
    so the coordinates can be read straight off the header without touching an annotation.
    A transcript outside every window has no unmasked sequence behind it, so keeping it would
    only inflate the bundle.
    """
    count = 0
    with opener(cdna) as fh, open(out_path, "w") as out:
        keeping = False
        for line in fh:
            if line.startswith(">"):
                keeping = False
                locus = next((f for f in line.split()
                              if f.startswith(("chromosome:", "scaffold:"))), None)
                if locus:
                    parts = locus.split(":")
                    if len(parts) >= 5 and parts[2] == contig:
                        try:
                            start, end = int(parts[3]), int(parts[4])
                        except ValueError:
                            start = end = None
                        if start is None:
                            keeping = True
                        elif not windows:
                            keeping = True
                        else:
                            keeping = any(start - 1 < w_end and end > w_start
                                          for w_start, w_end in windows)
                if keeping:
                    count += 1
                    out.write(line)
                continue
            if keeping:
                out.write(line)
    return count


def in_windows(position: int, windows: list[tuple[int, int]] | None) -> bool:
    if not windows:
        return True
    return any(start <= position < end for start, end in windows)


def chimeric_reads(path: str, contig: str, windows: list[tuple[int, int]] | None) -> set[str]:
    """Read names of chimeric records with both segments on the contig, inside the windows."""
    names: set[str] = set()
    with opener(path) as fh:
        for line in fh:
            if not line.strip() or line[0] == "#" or line.startswith("chr_donorA"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) <= CHIM_READ:
                continue
            if f[CHIM_DONOR_CHR] != contig or f[CHIM_ACC_CHR] != contig:
                continue
            try:
                donor, acceptor = int(f[1]), int(f[4])
            except ValueError:
                continue
            if in_windows(donor, windows) and in_windows(acceptor, windows):
                names.add(f[CHIM_READ])
    return names


def supporting_reads(path: str, contig: str) -> set[str]:
    """Read names from pf-supporting-reads.tab whose construct is on the contig."""
    names: set[str] = set()
    prefix = f"{contig}:"
    if not os.path.isfile(path):
        return names
    with opener(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > 2 and f[2].startswith(prefix):
                names.add(f[0])
    return names


def background_reads(bam: str, contig: str, count: int, seed: int,
                     windows: list[tuple[int, int]] | None = None) -> set[str]:
    """A reproducible sample of read names mapped to the contig."""
    if not count:
        return set()
    if not os.path.isfile(bam):
        log(f"no BAM at {bam}; skipping the background sample")
        return set()

    # Reservoir sampling over the stream, so the whole name list never has to be held.
    rng = random.Random(seed)
    reservoir: list[str] = []
    seen = 0
    regions = [f"{contig}:{start + 1}-{end}" for start, end in (windows or [])] or [contig]
    cmd = ["samtools", "view", "-F", "0x4", bam, *regions]
    with subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            tab = line.find("\t")
            if tab < 0:
                continue
            name = line[:tab]
            seen += 1
            if len(reservoir) < count:
                reservoir.append(name)
            else:
                j = rng.randrange(seen)
                if j < count:
                    reservoir[j] = name
        if proc.wait() != 0:
            fail(f"samtools view failed on {bam}")
    log(f"sampled {len(reservoir)} of {seen} alignment(s) on {contig} as background")
    return set(reservoir)


def extract_fastq(sources: list[str], names: set[str], out_path: str) -> tuple[int, int]:
    """Streams the source FASTQs once, keeping records whose name is wanted.

    Mates are pooled into one single-end file, which is what PFv2 consumes.
    """
    kept = 0
    total = 0
    with gzip.open(out_path, "wt", compresslevel=6) as out:
        for source in sources:
            with opener(source) as fh:
                while True:
                    header = fh.readline()
                    if not header:
                        break
                    seq, plus, qual = fh.readline(), fh.readline(), fh.readline()
                    if not qual:
                        break
                    total += 1
                    name = header[1:].split()[0] if len(header) > 1 else ""
                    # Mates may carry a /1 or /2 suffix on the name itself.
                    if name.endswith(("/1", "/2")):
                        name = name[:-2]
                    if name in names:
                        out.write(header)
                        out.write(seq)
                        out.write(plus)
                        out.write(qual)
                        kept += 1
            log(f"scanned {source}")
    return kept, total


def sha256(path: str) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genome", required=True, help="full genome FASTA")
    p.add_argument("--cdna", required=True, help="full transcriptome cDNA FASTA")
    p.add_argument("--contig", required=True, help="chromosome to build the bundle from")
    p.add_argument("--fastq", nargs="+", required=True, help="source FASTQ files")
    p.add_argument("--chimeric", required=True, help="STAR Chimeric.out.junction from a full run")
    p.add_argument("--supporting", required=True, help="pf-supporting-reads.tab from a full run")
    p.add_argument("--structures", required=True, help="pf-structures.bed from a full run")
    p.add_argument("--bam", default="", help="coordinate-sorted BAM, for the background sample")
    p.add_argument("--background", type=int, default=250_000,
                   help="background reads to sample (default 250000; 0 disables)")
    p.add_argument("--seed", type=int, default=7, help="sampling seed, for reproducibility")
    p.add_argument("--sample-id", default="testdata", help="name used for the reads file")
    p.add_argument("--flank", type=int, default=0,
                   help="when set, mask the contig to N outside each expected junction "
                        "plus this many bases either side; keeps coordinates real while "
                        "making the bundle small enough to commit (default 0: whole contig)")
    p.add_argument("--gzip-genome", action="store_true",
                   help="write genome.fa.gz instead of genome.fa; worth it for a masked "
                        "reference, whose N runs compress to almost nothing")
    p.add_argument("--out", required=True, help="bundle directory to create")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    out = args.out
    os.makedirs(out, exist_ok=True)
    expected_dir = os.path.join(out, "expected")
    os.makedirs(expected_dir, exist_ok=True)

    print(f"Building the PFv2 test bundle from contig {args.contig}")

    sequence = read_contig(args.genome, args.contig)
    windows = windows_from(args.structures, args.contig, args.flank) if args.flank else None
    if args.flank and not windows:
        fail(f"--flank was given but no expected junctions are on contig {args.contig}")
    if windows:
        log(f"masking to {len(windows)} window(s) around the expected junctions "
            f"(+/- {args.flank:,} bp)")

    genome_name = "genome.fa.gz" if args.gzip_genome else "genome.fa"
    genome_out = os.path.join(out, genome_name)
    length, kept = write_contig(sequence, args.contig, genome_out, windows)
    del sequence
    log(f"{genome_name}: {length:,} positions, {kept:,} unmasked "
        f"({100 * kept / length:.2f}%)")

    cdna_out = os.path.join(out, "transcriptome.fa")
    transcripts = extract_transcripts(args.cdna, args.contig, cdna_out, windows)
    log(f"transcriptome.fa: {transcripts:,} transcript(s)")
    if not transcripts:
        log("WARNING no transcripts matched; the transcriptomic filter will have no reference")

    chim = chimeric_reads(args.chimeric, args.contig, windows)
    sup = supporting_reads(args.supporting, args.contig)
    bg = background_reads(args.bam, args.contig, args.background, args.seed, windows)
    names = chim | sup | bg
    log(f"read names: {len(chim):,} chimeric, {len(sup):,} supporting, "
        f"{len(bg):,} background, {len(names):,} unique")
    if not sup:
        log("WARNING no supporting reads found; the expected output may be unreachable")

    reads_out = os.path.join(out, f"{args.sample_id}.fq.gz")
    kept, scanned = extract_fastq(args.fastq, names, reads_out)
    log(f"{args.sample_id}.fq.gz: {kept:,} record(s) from {scanned:,} scanned")

    expected_bed = os.path.join(expected_dir, "pf-structures.bed")
    rows = []
    with open(args.structures) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f and f[0] == args.contig:
                rows.append(line.rstrip("\n"))
    rows.sort(key=lambda r: (int(r.split("\t")[1]), int(r.split("\t")[2])))
    with open(expected_bed, "w") as fh:
        for r in rows:
            fh.write(r + "\n")
    log(f"expected/pf-structures.bed: {len(rows)} junction(s)")

    files = {}
    for root, _dirs, filenames in os.walk(out):
        for name in filenames:
            path = os.path.join(root, name)
            rel = os.path.relpath(path, out)
            if rel in ("MANIFEST.json", "SHA256SUMS"):
                continue
            files[rel] = {"bytes": os.path.getsize(path), "sha256": sha256(path)}

    manifest = {
        "bundle": "pfv2-test-data",
        "built": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "contig": args.contig,
        "contig_length": length,
        "unmasked_bases": kept,
        "windows": [{"start0": a, "end0": b} for a, b in (windows or [])],
        "flank": args.flank,
        "sample_id": args.sample_id,
        "reads": {"records": kept, "chimeric": len(chim), "supporting": len(sup),
                  "background": len(bg), "seed": args.seed},
        "expected_junctions": len(rows),
        "provenance": {
            "genome": os.path.abspath(args.genome),
            "cdna": os.path.abspath(args.cdna),
            "chimeric": os.path.abspath(args.chimeric),
            "structures": os.path.abspath(args.structures),
            "fastq": [os.path.abspath(f) for f in args.fastq],
        },
        "files": files,
    }
    with open(os.path.join(out, "MANIFEST.json"), "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with open(os.path.join(out, "SHA256SUMS"), "w") as fh:
        for rel in sorted(files):
            fh.write(f"{files[rel]['sha256']}  {rel}\n")

    total = sum(f["bytes"] for f in files.values())
    log(f"bundle total: {total / 1e6:.1f} MB across {len(files)} file(s)")
    print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
