# PFv2 on Prefect

Each PFv2 stage is a Prefect task, so a failure is attributable to a stage and a retry
resumes from it rather than repeating a two-hour alignment.

## Setup

```bash
pip install prefect
prefect work-pool create pfv2 --type process
prefect deploy --prefect-file deploy/prefect/prefect.yaml
```

## Run one sample

```bash
prefect deployment run 'pfv2/pfv2-local' \
  --param sample_id=SAMEA120815325 \
  --param 'fastq=["/data/reads/SAMEA120815325/reads_se.fq.gz"]' \
  --param outdir=/data/results \
  --param genome=/refs/genome.fa \
  --param star_index=/refs/indexes/star \
  --param bowtie2_genome=/refs/indexes/bowtie2/genome \
  --param bowtie2_transcriptome=/refs/indexes/bowtie2/transcriptome \
  --param read_length=150
```

## Run a cohort

A cohort is this flow submitted once per sample. PFv2 works on one library at a time and
the samples are independent, so there is nothing for a cohort-level flow to coordinate
beyond concurrency, which the work pool already bounds:

```bash
prefect work-pool set-concurrency-limit pfv2 4
while read -r sample fastq; do
  prefect deployment run 'pfv2/pfv2-local' \
    --param "sample_id=$sample" --param "fastq=[\"$fastq\"]" \
    --param outdir=/data/results --param read_length=150 \
    --param genome=/refs/genome.fa \
    --param star_index=/refs/indexes/star \
    --param bowtie2_genome=/refs/indexes/bowtie2/genome \
    --param bowtie2_transcriptome=/refs/indexes/bowtie2/transcriptome
done < samples.tsv
```

## Notes

- `fetch_reference` is the only retried-by-default network task, and it verifies every
  file against the MD5 Zenodo publishes before accepting it.
- `validate` is deliberately not retried: a reference naming mismatch fails identically
  every time, and failing in seconds is the point.
- Stage caching is keyed on task inputs with a 7-day expiry, so re-submitting an identical
  run reuses the reference fetch rather than repeating it.
- `PFV2_HOME` tells the flow where the PFv2 checkout or install lives; it defaults to the
  repository root two levels above `flow.py`, and the container image sets it to
  `/opt/pfv2`.
