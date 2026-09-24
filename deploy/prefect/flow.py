"""PFv2 as a Prefect flow.

Prefect is the orchestration layer here, not the compute substrate. Each PFv2 stage is a
task so that a failure is attributable to a stage and a retry resumes from it, rather than
re-running a two-hour alignment because the filter stage hit a full disk.

The stage boundaries are the same ones PFv2.sh uses, and for the same reason: each writes
files the next one reads, so a stage is the natural unit of both caching and blame.

    prefect work-pool create pfv2 --type process
    prefect deploy --prefect-file deploy/prefect/prefect.yaml
    prefect deployment run 'pfv2/pfv2-local' \
      --param sample_id=S01 \
      --param fastq=/data/reads/S01.fq.gz \
      --param outdir=/data/results \
      --param genome=/refs/genome.fa \
      --param star_index=/refs/indexes/star \
      --param bowtie2_genome=/refs/indexes/bowtie2/genome \
      --param bowtie2_transcriptome=/refs/indexes/bowtie2/transcriptome \
      --param read_length=150
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
from datetime import timedelta
from pathlib import Path

from prefect import flow, get_run_logger, task
from prefect.tasks import task_input_hash

PFV2_HOME = Path(os.environ.get("PFV2_HOME", Path(__file__).resolve().parents[2]))

# A stage's inputs fully determine its outputs, so re-submitting the same run reuses them.
# Kept short enough that a reference rebuilt in place does not serve a stale cache for long.
CACHE_TTL = timedelta(days=7)


class StageError(RuntimeError):
    pass


def _run(cmd: list[str], *, log, cwd: Path | None = None) -> None:
    log.info("$ %s", " ".join(str(c) for c in cmd))
    try:
        subprocess.run([str(c) for c in cmd], check=True, cwd=cwd)
    except FileNotFoundError as exc:
        raise StageError(f"{cmd[0]} is not installed: {exc}") from exc
    except subprocess.CalledProcessError as exc:
        raise StageError(f"{cmd[0]} exited {exc.returncode}") from exc


@task(retries=2, retry_delay_seconds=30, cache_key_fn=task_input_hash, cache_expiration=CACHE_TTL)
def fetch_reference(record: str, dest: str, *, latest: bool = False,
                    only: list[str] | None = None) -> str:
    """Pull a prebuilt reference from Zenodo, verifying each file's published MD5.

    Retried, because this is the one stage whose failures are usually transient.
    """
    log = get_run_logger()
    cmd = [PFV2_HOME / "scripts" / "fetch_references.py", record, "--dest", dest,
           "--manifest", str(Path(dest) / "zenodo-manifest.json")]
    if latest:
        cmd.append("--latest")
    if only:
        cmd += ["--only", *only]
    _run(["python3", *cmd], log=log)
    return dest


@task
def validate(fastq: list[str], genome: str, star_index: str, bowtie2_genome: str,
             bowtie2_transcriptome: str | None, read_length: int, report: str) -> dict:
    """Check inputs and references before any compute is spent.

    Not retried: a naming mismatch between the genome FASTA and the STAR index will fail
    the same way every time, and failing fast is the whole point of this stage.
    """
    log = get_run_logger()
    cmd = ["python3", PFV2_HOME / "scripts" / "validate_inputs.py",
           "--fastq", *fastq,
           "--genome", genome,
           "--star-index", star_index,
           "--bowtie2-genome", bowtie2_genome,
           "--read-length", str(read_length),
           "--json", report]
    if bowtie2_transcriptome:
        cmd += ["--bowtie2-transcriptome", bowtie2_transcriptome]
    _run(cmd, log=log)
    return json.loads(Path(report).read_text())


@task(retries=1, retry_delay_seconds=60)
def run_pfv2(sample_id: str, fastq: list[str], outdir: str, genome: str, star_index: str,
             bowtie2_genome: str, bowtie2_transcriptome: str, read_length: int,
             junction_span: int, pid: float, max_genomic_span: int, min_genomic_span: int,
             threads: int, java_heap: str, filters: str, keep_intermediates: bool,
             legacy: bool, aligned_strand: bool) -> str:
    """Run the pipeline itself.

    Retried once, because the failures worth retrying at this level are environmental --
    a full scratch disk, a preempted node. A second deterministic failure is a real one.
    """
    log = get_run_logger()
    cmd = [PFV2_HOME / "PFv2.sh",
           "-c", PFV2_HOME,
           "-i", sample_id,
           "-r", fastq[0],
           "-d", outdir,
           "-g", genome,
           "-S", star_index,
           "-b", bowtie2_genome,
           "-t", bowtie2_transcriptome,
           "-l", read_length,
           "-j", junction_span,
           "-p", pid,
           "-C", max_genomic_span,
           "-M", min_genomic_span,
           "-n", threads,
           "-m", java_heap,
           # Validation ran as its own task; repeating it would re-read the reference.
           "-V"]
    if filters == "genomic":
        cmd.append("-G")
    elif filters == "transcriptomic":
        cmd.append("-T")
    if keep_intermediates:
        cmd.append("-k")
    if legacy:
        cmd.append("-L")
    if aligned_strand:
        cmd.append("-A")

    _run(["bash", *cmd], log=log)
    return str(Path(outdir) / "PF" / sample_id)


@task
def collect(result_dir: str, sample_id: str) -> dict:
    """Summarise a finished run.

    A run that called nothing is a valid outcome, and this is the record that distinguishes
    it from a run that did not happen.
    """
    log = get_run_logger()
    work = Path(result_dir)
    bed = work / "pf-structures.bed"
    jpms = work / f"{sample_id}_jpms.tsv"

    junctions = sum(1 for _ in open(bed)) if bed.is_file() else 0
    reads = 0
    supporting = work / "pf-supporting-reads.tab"
    if supporting.is_file():
        reads = sum(1 for _ in open(supporting))

    summary = {
        "sample_id": sample_id,
        "junctions": junctions,
        "supporting_reads": reads,
        "outputs": {
            "structures": str(bed) if bed.is_file() else None,
            "jpms": str(jpms) if jpms.is_file() else None,
        },
        "empty_result": junctions == 0,
    }
    (work / "pfv2-summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    log.info("%s: %d junction(s) from %d supporting read(s)", sample_id, junctions, reads)
    return summary


@flow(name="pfv2", log_prints=True)
def pfv2_flow(
    sample_id: str,
    fastq: list[str],
    outdir: str,
    genome: str,
    star_index: str,
    bowtie2_genome: str,
    bowtie2_transcriptome: str,
    read_length: int,
    junction_span: int = 8,
    pid: float = 0.85,
    max_genomic_span: int = 1_000_000,
    min_genomic_span: int = 50,
    threads: int = 16,
    java_heap: str = "20G",
    filters: str = "both",
    keep_intermediates: bool = False,
    legacy: bool = False,
    aligned_strand: bool = False,
    zenodo_record: str | None = None,
    zenodo_dest: str | None = None,
    zenodo_latest: bool = False,
    skip_validation: bool = False,
) -> dict:
    """Identify backsplice junctions in one sample.

    A cohort is this flow submitted once per sample. PFv2 works on one library at a time and
    the samples are independent, so there is nothing for a cohort-level flow to coordinate
    beyond concurrency, which the work pool already bounds.
    """
    log = get_run_logger()
    Path(outdir).mkdir(parents=True, exist_ok=True)

    if zenodo_record:
        if not zenodo_dest:
            raise ValueError("zenodo_dest is required when zenodo_record is given")
        fetch_reference(zenodo_record, zenodo_dest, latest=zenodo_latest)

    if skip_validation:
        log.warning("validation skipped; a reference naming mismatch will not be caught")
    else:
        report = validate(
            fastq, genome, star_index, bowtie2_genome, bowtie2_transcriptome, read_length,
            str(Path(outdir) / f"{sample_id}-validation.json"),
        )
        for warning in report.get("warnings", []):
            log.warning("%s", warning)

    result_dir = run_pfv2(
        sample_id, fastq, outdir, genome, star_index, bowtie2_genome, bowtie2_transcriptome,
        read_length, junction_span, pid, max_genomic_span, min_genomic_span, threads,
        java_heap, filters, keep_intermediates, legacy, aligned_strand,
    )
    return collect(result_dir, sample_id)


if __name__ == "__main__":
    pfv2_flow()
