# PFv2

PTESFinder v2 — an annotation-free computational method to identify
post-transcriptional exon shuffling (PTES) / backsplice junctions from RNA-seq data.

## Changes from v1

- Candidate backsplices are called from STAR chimeric alignments, so no transcript
  annotation is required.
- The bowtie1 dependency in the identification phase has been removed.

## How it works

| Stage | Step | Output |
|-------|------|--------|
| 1 | Map reads to the genome with STAR (chimeric detection on) and to the genome and transcriptome with Bowtie2 | `star_Chimeric.out.junction`, `star_SJ.out.tab`, `genomic.sam`, `transcriptomic.sam` |
| 2 | Call candidate backsplice and canonical junction coordinates | `putative_structures.txt`, `canonical_structures.txt` |
| 3 | Build junction-spanning sequence constructs and Bowtie2 indexes for them | `Constructs.fa`, `Can.fa` |
| 4 | Re-map the reads to the candidate junctions | `ptes.sam`, `canonical.sam` |
| 5 | Filter false positives and count supporting reads | `pf-structures.bed`, `<sample>_jpms.tsv` |

A read counts as evidence for a backsplice only if it aligns *better* to the
construct than to the genome and to the transcriptome (more aligned bases and
strictly fewer mismatches), and it spans the junction by at least `-j/2`
positions on each side with at least `-p` percent identity in each flank.

## Quick start with Docker

The image carries PFv2 and the three tools it shells out to, for linux/amd64 and
linux/arm64, so a run needs nothing from the host but reads and a reference.

```bash
docker run --rm conidiobolus/pfv2:latest version

docker run --rm \
  -v /data/reads:/reads:ro -v /data/refs:/refs:ro -v /data/results:/work \
  conidiobolus/pfv2:latest \
  -i SAMEA120815325 -r /reads/SAMEA120815325/reads_se.fq.gz -d /work \
  -g /refs/genome.fa -S /refs/indexes/star \
  -b /refs/indexes/bowtie2/genome -t /refs/indexes/bowtie2/transcriptome \
  -l 150 -n 16 -m 20G
```

The image also exposes the helpers as subcommands:

| Command | Does |
|---------|------|
| `validate` | check inputs and references without running anything |
| `fetch-references` | pull a prebuilt reference from a Zenodo record |
| `selftest` | run the end-to-end smoke test inside the image |
| `version` | print the PFv2, STAR, Bowtie2, samtools and Java versions |

Aligner indexes are portable between architectures: the STAR and Bowtie2 formats are
little-endian 64-bit with no architecture-specific packing, so an index built on x86_64
loads unchanged under the aarch64 binaries.

## Running on a cluster

- **Kubernetes** — `deploy/k8s` is a Helm chart rendering one Job per sample, with
  validation as an init container so a bad reference fails in seconds rather than hours.
  See [deploy/k8s/README.md](deploy/k8s/README.md).
- **Prefect** — `deploy/prefect` runs each stage as a Prefect task, so a retry resumes
  from the failed stage instead of repeating the alignment. See
  [deploy/prefect/README.md](deploy/prefect/README.md).

## Validating inputs before a run

A run validates automatically; `-V` skips it. To check a reference set on its own, before
committing a cluster to it:

```bash
ptesfinder validate --fastq reads.fq.gz --genome genome.fa --star-index star/ --bowtie2-genome idx/genome --bowtie2-transcriptome idx/transcriptome --read-length 150
```

It checks FASTQ structure and gzip integrity, the read length against `-l`, that the STAR
index is complete, that both Bowtie2 indexes exist and are readable, and — the check that
matters most — that the sequence names in the genome FASTA, the STAR index and the Bowtie2
index all agree. `chr1` and `1` are not interchangeable: mixing them produces an empty
result many hours later, which looks like a biological finding rather than a
misconfiguration.

## Test data and references

### Quick start: the committed fixture

```bash
make e2e
```

No download. `test/e2e/fixture` is chromosome 17 of GRCm38 with everything outside eight
known backsplice junctions masked to `N`, the 215 transcripts overlapping those windows, and
the reads needed to rediscover them. Masking rather than slicing keeps every coordinate
identical to the source assembly, so the expected output is real chromosome 17 coordinates,
while the `N` runs compress the reference to 1.5 MB.

| File | Size | Contents |
|------|------|----------|
| `genome.fa.gz` | 1.5 MB | chromosome 17, 4.4% unmasked |
| `bl6-chr17.fq.gz` | 6.3 MB | 79,592 reads: chimeric, supporting, and background |
| `transcriptome.fa` | 456 KB | 215 transcripts in the windows |
| `expected/pf-structures.bed` | 492 B | the 8 junctions the full-depth run called here |

8.2 MB in total, and the run takes about two minutes on 16 cores. It recovers all eight
junctions with the same 15 supporting reads the full-depth run found on this chromosome.

It ships FASTA, not indexes, and builds the indexes as part of the test. That is deliberate:
a STAR index records the genome format version that wrote it and STAR refuses one it cannot
read, so a committed index would turn an aligner upgrade into a confusing failure. A 95 Mb
reference indexes in about a minute.

### Prebuilt references

Building a mouse STAR index takes about an hour and ~32 GB of RAM, so the prebuilt set is
published rather than rebuilt per site. Fetch it with the CLI, which verifies every file
against the MD5 Zenodo publishes before putting it in place:

```bash
ptesfinder fetch-references <record> -d /data/refs --manifest /data/refs/provenance.json
```

| DOI | Contents | Size |
|-----|----------|------|
| [10.5281/zenodo.22929132](https://doi.org/10.5281/zenodo.22929132) | test data: chromosome 17, unmasked, with a 250k-read background sample | 76 MB |
| [10.5281/zenodo.22929143](https://doi.org/10.5281/zenodo.22929143) | reference: GRCm38 genome FASTA, Ensembl 102 cDNA, STAR index, both Bowtie2 indexes | 27 GB |

These are **concept DOIs**: each always resolves to the newest version of its record, so the
commands below keep working when a record is superseded. Cite the specific version DOI shown
on the record page in a methods section instead.

They are separate records on purpose: a getting-started run should not be a 27 GB download.
The reference set unpacks from four archives, so the STAR index can be skipped by a site that
builds its own:

```bash
ptesfinder fetch-references 10.5281/zenodo.22929143 --latest -d /data/refs --only star_index
```

The reference set is published as four archives. The STAR index is split into 4 GB parts
because a single upload that large is rejected by Zenodo's gateway; `fetch-references` rejoins
them and checks the result against the whole-file checksum, so the split is invisible unless
a part fails to download.

```bash
bash test/e2e/run.sh --zenodo
```

`fetch-references` accepts a record id, a version DOI, a concept DOI with `--latest`, or a
record URL. A file whose MD5 does not match is deleted rather than used — a truncated genome
otherwise produces a run that looks fine and is not — and a verified copy already on disk is
skipped, so an interrupted fetch can simply be re-run. Files are recognised by role from
their names (`genome_fasta`, `star_index`, `bowtie2_genome`, `bowtie2_transcriptome`,
`annotation_gtf`), and `--only` takes a role or a file name, so the STAR index can be
skipped when a site builds its own.

Pin a version DOI for anything whose results have to be reproducible. A concept DOI with
`--latest` follows the newest version, which is what you want for "give me the current
reference" and not what you want in a methods section.

Zenodo's default quota is 50 GB per record, so the ~33 GB reference fits — but only just, and
adding a second assembly to the same record would not.

Two things to check before using a prebuilt reference, because getting either wrong costs
hours and does not raise an error:

- **Sequence naming.** The published GRCm38 set is Ensembl-named: contigs are `1`…`19`, `X`,
  `Y`, with no `chr` prefix. It is **not** interchangeable with a UCSC mm10 reference, whose
  contigs are `chr1`…`chr19`. Mixing conventions produces an empty result, not an error.
  `ptesfinder validate` compares the names across the genome FASTA, the STAR index and the
  Bowtie2 index and fails on a mismatch.
- **Aligner versions.** A STAR index is only loadable by releases sharing its genome format
  version, and its `sjdbOverhang` is tuned to a read length. Both are recorded in the
  bundle's `MANIFEST.json`, and the container image pins the aligner versions the published
  index was built for.

### Rebuilding the bundles

Both are reproducible rather than opaque archives:

```bash
python3 ops/make_test_bundle.py --help
```

```bash
bash ops/make_reference_bundle.sh --help
```

`ops/zenodo_upload.sh` uploads a bundle directory to a Zenodo **draft** and prints the link.
It never publishes: a published record's files cannot be changed afterwards, only superseded
by a new version, so pressing publish is left to a person.

## Dependencies

- Java 16 or newer
- Python 3
- [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- samtools

Memory use is genome-size and aligner dependent; see the STAR documentation.
The Java stages default to a 20 GB heap, which is tuneable with `-m`.

## Build

```bash
make build
```

This compiles the sources, runs the test suite and packages `PFv2.jar`. The jar is only
written when the tests pass. `bash setup.sh` does the same thing directly, and
`bash setup.sh --skip-tests` bypasses the tests.

Run `make` on its own to list every target:

| Target | Does |
|--------|------|
| `make build` | compile, test and package `PFv2.jar` |
| `make smoke` | end-to-end smoke test over a synthetic dataset |
| `make install` | symlink the `pfv2` CLI into `PREFIX/bin` (default `/usr/local`) |
| `make image-local` | build the container image for this machine and self-test it |
| `make image` | build and push a multi-arch image |
| `make clean` | remove build output |

## The ptesfinder command

`bin/ptesfinder` is a single entry point for everything. It is also the container image's
entrypoint, so a subcommand behaves the same from a checkout, from an install on `PATH`, or
inside the image. It is installed as `pfv2` as well, so either name works.

```bash
make install
```

| Command | Does |
|---------|------|
| `ptesfinder run` | identify backsplice junctions in one sample |
| `ptesfinder validate` | check reads and references without running anything |
| `ptesfinder fetch-references` | download a prebuilt reference set from a Zenodo record |
| `ptesfinder selftest` | run the end-to-end smoke test |
| `ptesfinder version` | print the release and the aligner versions |

Bare pipeline flags work without the `run` subcommand, so existing `PFv2.sh` invocations
keep working unchanged. Without installing, call `./bin/ptesfinder` or `bash PFv2.sh` from
the checkout. `make install PREFIX=~/.local` installs without needing root.

`version` names the release rather than leaving it implied, because "PTESFinder" on its own
does not identify a tool: v1 is a different implementation with different semantics, and a
results table that records only the name cannot be reproduced.

## Input files

- RNA-seq reads in Illumina FASTQ format (plain or gzipped)
- Genome reference in FASTA format
- Pre-built STAR genome index
- Pre-built Bowtie2 genome index
- Pre-built Bowtie2 transcriptome index

Sequence names in the genome FASTA must match those used to build the STAR index —
`chr1` and `1` are not interchangeable. Paired-end reads must be pooled into a
single FASTQ with unique read ids.

## Running

```bash
ptesfinder run -i SRR364679 -r SRR364679.fastq -d SRR364679/ -S STAR/ \
         -t transcriptome-index-bowtie -g genome.fasta \
         -b genome-index-bowtie -l 100
```

`bash PFv2.sh` with the same flags is equivalent. `-c` defaults to the directory holding the
script, so it only needs passing when the code lives somewhere else.

### Parameters

**Mandatory**

| Flag | Meaning |
|------|---------|
| `-r` | sequence reads in FASTQ format |
| `-d` | working directory |
| `-i` | sample id |
| `-t` | transcriptome reference Bowtie2 index prefix |
| `-g` | genome reference in FASTA format |
| `-b` | genome reference Bowtie2 index prefix |
| `-l` | average read length |
| `-S` | path to the pre-built STAR genome index directory |

**Optional**

| Flag | Default | Meaning |
|------|---------|---------|
| `-c` | script directory | PFv2 code directory |
| `-p` | 0.85 | minimum percent identity per flank, 0–1; ideal 0.60–0.95 |
| `-j` | 8 | junction span, even integer; ideal 4–14 |
| `-C` | 1000000 | maximum backsplice genomic span, bp |
| `-M` | 50 | minimum backsplice genomic span, bp |
| `-n` | 16 | threads for STAR and Bowtie2 |
| `-m` | 20G | Java heap for the PFv2 stages |
| `-G` | off | run the genomic filter only |
| `-T` | off | run the transcriptomic filter only |
| `-k` | off | keep intermediate SAM/FASTA/index files |
| `-L` | off | reproduce the filter semantics of releases up to 2.1.0 |
| `-A` | off | report STAR's aligned strand instead of the motif-derived one |
| `-V` | off | skip input and reference validation |
| `-h` | | show usage |

### The strand column

STAR reports the strand a chimeric segment *aligned* to. For an unstranded library a
back-splicing read aligns opposite its host gene about half the time, so that strand says
nothing about the gene. The splice motif does.

Since 2.2.0 a junction whose signal reads as the reverse complement of a canonical motif is
reported on the other strand, with the signal and construct turned to match. This is the
convention the canonical junctions already follow, because STAR derives their strand from
the motif in `SJ.out.tab`, and it resolves the strand-convention discrepancy described in
the circRNA detection benchmarking study of Vromman et al. (2023), listed under
[References](#references).

`-A` restores the previous column for anyone reproducing pre-2.2.0 output. It should not be
used to feed a consensus that keys on strand.

## Output

Results are written to `<-d>/PF/<sample_id>/`:

| File | Contents |
|------|----------|
| `<sample>_jpms.tsv` | junction id, supporting read count, junctions per million |
| `pf-structures.bed` | identified backsplices as BED, with read counts and splice signal |
| `pf-structure-counts.tsv` | supporting read count per construct |
| `pf-supporting-reads.tab` | the accepted alignments |
| `pf-pid.tsv` | per-read edit distance and per-flank percent identity |
| `pf-flanking-canonical-junctions.bed` | canonical junctions used for normalisation |
| `run.log` | stage-by-stage log from the Java components |

Intermediate SAM, FASTA and index files are removed at the end of a successful
run unless `-k` is given.

Junctions per million are computed against the total of backsplice-supporting and
canonical junction reads, so values are comparable between samples only when the
same reference and parameters were used.

## Repository layout

```
bin/ptesfinder   CLI entry point (also installed as pfv2), and the image's entrypoint
PFv2.sh          the pipeline itself
setup.sh         build, test and package
Dockerfile       multi-arch image; docker/ holds its entrypoint and tool pins
Makefile         build, test and image targets
deploy/k8s       Helm chart, one Job per sample
deploy/prefect   Prefect flow and deployments
src/             Java sources
  entities/      SAM records and junction structures
  utils/discovery/  STAR output to candidate junctions and constructs
  utils/filter/     genomic, transcriptomic, junction span and identity filters
  utils/annotate/   standalone post-hoc annotation tools, not part of the pipeline
  utils/init/       standalone reference preparation tools, not part of the pipeline
test/            unit suite, test/smoke (no aligners) and test/e2e (real aligners)
scripts/         STAR and Bowtie2 wrappers, input validation, Zenodo reference fetch
ops/             bundle builders and the Zenodo upload helper
lib/             Apache Commons dependencies
```

Only `utils/discovery` and `utils/filter` are invoked by `PFv2.sh`. The classes
under `utils/annotate` and `utils/init` are standalone utilities carried over
from v1, each with its own `main`.

## Licence

MIT. See [LICENSE](LICENSE).

## References

**PTESFinder (v1), the method this implements**

- Izuogu OG, Alhasan AA, Alafghani HM, Santibanez-Koref M, Elliott DJ, Jackson MS.
  *PTESFinder: a computational method to identify post-transcriptional exon shuffling
  (PTES) events.* BMC Bioinformatics 17:31 (2016).
  [PMC4711006](https://pmc.ncbi.nlm.nih.gov/articles/PMC4711006/) ·
  [doi:10.1186/s12859-016-0881-4](https://doi.org/10.1186/s12859-016-0881-4)
  Source: https://sourceforge.net/projects/ptesfinder-v1/

**Benchmarking**

- Vromman M, Anckaert J, Bortoluzzi S, ... Izuogu O, Jackson MS, Santibanez-Koref M, ...
  Vandesompele J, Volders PJ. *Large-scale benchmarking of circRNA detection tools reveals
  large differences in sensitivity but not in precision.* Nature Methods 20(8):1159-1169
  (2023). [PMID 37443337](https://pubmed.ncbi.nlm.nih.gov/37443337/) ·
  [doi:10.1038/s41592-023-01944-6](https://doi.org/10.1038/s41592-023-01944-6)

**Applications**

- Grassi L, Izuogu OG, Jorge NAN, et al. *Cell type-specific novel long non-coding RNA and
  circular RNA in the BLUEPRINT hematopoietic transcriptomes atlas.* Haematologica
  106(10):2613-2623 (2021).
  [PMC8485671](https://pmc.ncbi.nlm.nih.gov/articles/PMC8485671/) ·
  [doi:10.3324/haematol.2019.238147](https://doi.org/10.3324/haematol.2019.238147)

- Whittle B, Izuogu O, Lowes H, et al. *Early-stage idiopathic Parkinson's disease is
  associated with reduced circular RNA expression.* npj Parkinson's Disease 10:25 (2024).
  [nature.com/articles/s41531-024-00636-y](https://www.nature.com/articles/s41531-024-00636-y) ·
  [doi:10.1038/s41531-024-00636-y](https://doi.org/10.1038/s41531-024-00636-y)

## Releasing

`RELEASING.md` tracks release state and the steps that need a person.

## Support

osagie.izuogu@gmail.com
