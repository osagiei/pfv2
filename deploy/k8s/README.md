# PFv2 on Kubernetes

One Job per sample. PFv2 works on a single library at a time, so a cohort is a set of Jobs
rather than a workflow engine's concern.

## Install

```bash
helm install pfv2 deploy/k8s -f my-values.yaml
```

A minimal `my-values.yaml`:

```yaml
samples:
  - id: SAMEA120815325
    fastq: /reads/SAMEA120815325/reads_se.fq.gz
    readLength: 150

reference:
  genomeFasta: /refs/genome.fa
  starIndex: /refs/indexes/star
  bowtie2Genome: /refs/indexes/bowtie2/genome
  bowtie2Transcriptome: /refs/indexes/bowtie2/transcriptome

storage:
  reference: { claimName: pfv2-references }
  reads:     { claimName: pfv2-reads }
  output:    { claimName: pfv2-output }
```

Without Helm, `job-example.yaml` is the same thing for one sample.

## What the chart renders

For each sample, one Job with:

1. an optional `fetch-reference` init container that pulls a prebuilt reference from a
   Zenodo record and verifies each file's published MD5;
2. a `validate` init container that checks the reads, the genome FASTA and both index sets,
   and in particular that their sequence names agree — a `chr1` versus `1` mismatch
   otherwise produces an empty result hours later;
3. the `pfv2` container, which runs the pipeline with validation disabled because it has
   already run.

## Sizing

STAR holds the whole genome index resident and the construct generator holds the genome as
a string, so memory is the binding constraint. The defaults (40 GiB requested, 64 GiB
limit, 16 CPU) suit a mammalian genome. Under-requesting memory gets the pod OOM-killed
partway through mapping, which reads like a tool failure rather than a scheduling one.

`storage.scratch` is an `emptyDir` sized for the intermediates: construct indexes and four
SAM files, which for a deep library is tens of gigabytes. Results are written to the output
claim; intermediates never touch it.

## Security

Pods run as UID 1000 with a read-only root filesystem, no privilege escalation and all
capabilities dropped. `/tmp` and the scratch directory are the only writable paths. PFv2
talks to no Kubernetes API, so its ServiceAccount is created with no role and its token is
not mounted.
