# Changelog

## 2.2.0

A correctness pass over the filter and construct stages, plus container and cluster
runnables and an input validation step. Several fixes change results; each is called out
below and `-L` restores the previous behaviour where one exists.

### Fixed - these change results

- **The junction offset sat one base before the construct seam.** The anchor intervals are
  inclusive at both ends, so an arm of *n* bases spans `stop - start + 1`, but the offset
  recorded in the construct name was `stop - start`. The purity window was therefore
  centred one base upstream of the seam, covering 6 bases on one side and 4 on the other
  for the default `-j 8`. Verified against the real seam on a synthetic construct and now
  pinned by a test.

- **Construct arms were one base longer than the configured segment.** `-l 150` with the
  15 bp minimum overhang is meant to build 135 bp arms, and built 136, so the shortest
  overhang a read could have was 14 rather than 15. The anchor arithmetic now produces arms
  of exactly the segment length, which also brings the constructs to parity with the Python
  port.

- **A read no longer has to span the junction to be counted, unless `-L` is given.**
  Releases up to 2.1.0 accepted any read matching a construct perfectly (`NM:i:0`), whether
  or not it crossed the seam. A read lying wholly inside one arm says nothing about the
  junction joining the arms. `-L` restores the old behaviour.

- **Competing alignments are ranked by Bowtie2's `AS` tag rather than by MD and NM.**
  Soft-clipped bases contribute to neither MD nor NM, so a short perfect genomic fragment
  looked flawless to the old comparison and beat a full-length construct alignment carrying
  two mismatches - discarding true positives. `AS` accounts for clipping and gap penalties.
  When either alignment lacks `AS` the comparison falls back to the old tags. `-L` selects
  the old metric.

- **Canonical junctions now compete for reads.** A read spanning a real forward splice
  junction is a linear read, however well it also fits a backsplice construct. The genomic
  filter cannot catch these, because such a read does not align contiguously to the genome
  at all. The canonical constructs, previously built only to supply the normalisation
  denominator, are now also used as a competitor; only canonical alignments that themselves
  pass the span and identity filters count. `-L` skips this comparison.

- **The canonical motif rule matched its own comment.** `SJ.out.tab` motifs above 2 were
  dropped with the comment "non-canonical splice site", which in fact admitted motif 0 (the
  genuinely non-canonical one) and excluded the GC-AG pair. Motifs 1 to 4 are now accepted,
  motif 0 rejected, and a junction with no uniquely mapping read is rejected. On a mouse
  sample this moved the denominator from 198,288 to 192,382 junctions. `-L` restores the
  old rule.

- **`CTAC` is no longer accepted as a canonical splice signal.** The signal is already
  reverse-complemented for minus-strand junctions, so a correctly stranded canonical
  junction always reads `GTAG`; accepting `CTAC` as well admitted junctions whose motif
  contradicts the strand STAR assigned. The accepted set is now GT-AG, GC-AG and AT-AC,
  matching the motifs the candidate filter passes through. `-L` restores `{GTAG, CTAC}`.

- **The best alignment per read is kept, not the last one in the file.** Reads were stored
  in a map keyed by read name, so a second alignment silently replaced the first regardless
  of quality.

- **Constructs spanning an assembly gap are dropped.** A run of ten or more Ns means any
  alignment there is meaningless.

- **Optional SAM tags are located by prefix** (`MD:Z:`, `NM:i:`, `AS:i:`) rather than by
  offset from the end of the line.

### Changed - the strand column

Backsplice junctions are now reported on the strand their splice motif implies, rather than
the strand STAR's chimeric segment aligned to. An unstranded library aligns a back-splicing
read opposite its host gene about half the time, so the aligned strand carries no
information about the gene; the motif does, and it is already how the canonical junctions
are reported because STAR derives their strand from the motif in `SJ.out.tab`.

This resolves the strand-convention discrepancy described in the circRNA detection
benchmarking study of Vromman et al., Nature Methods 20(8):1159-1169 (2023),
[PMID 37443337](https://pubmed.ncbi.nlm.nih.gov/37443337/). See the README's References
section for the full list.

Across a mouse test sample the reported splice signal became canonical for the whole
candidate set, evenly split between strands, matching the profile of the canonical
junctions from the same run.

`-A` restores the previous column, and `-L` implies it.

### Added

- **`-A`, to report STAR's aligned strand instead of the motif-derived one.** Off by
  default: the strand column now comes from the splice motif.

- **A single `ptesfinder` CLI** (`bin/ptesfinder`, also installed as `pfv2`) with `run`,
  `validate`, `fetch-references`,
  `selftest` and `version` subcommands. It is also the container image's entrypoint, so the
  subcommands behave identically from a checkout, from an install on `PATH`, or inside the
  image - previously they existed only inside the image and a native user had to call the
  helper scripts by path. Bare pipeline flags still work without the `run` subcommand, so
  existing `PFv2.sh` invocations are unaffected. `make install` symlinks it into
  `PREFIX/bin` under both names. `version` names the release explicitly, since "PTESFinder"
  without one does not identify a tool.

- **Multi-arch Docker image** (`linux/amd64`, `linux/arm64`), carrying pinned STAR 2.7.11b,
  Bowtie2 2.5.4 and samtools 1.21 from bioconda. The image self-tests at build time, so a
  broken build cannot be pushed. `Makefile` has `image-local` and `image` targets and
  `.github/workflows/image.yaml` publishes on a tag.

- **Kubernetes deployment** under `deploy/k8s`: a Helm chart rendering one Job per sample,
  with reference fetching and validation as init containers, a read-only root filesystem
  and a non-root pod security context.

- **Prefect deployment** under `deploy/prefect`: each stage is a task, so a retry resumes
  from the failed stage rather than repeating the alignment.

- **`scripts/validate_inputs.py`**, run automatically by `PFv2.sh` and skippable with `-V`.
  It checks FASTQ structure and gzip integrity, the read length against `-l`, STAR index
  completeness, Bowtie2 index presence and readability, and that the sequence names in the
  genome FASTA, the STAR index and the Bowtie2 index agree.

- **`scripts/fetch_references.py`**, which pulls prebuilt references from a Zenodo record by
  id, version DOI, concept DOI or URL, verifying each file against its published MD5 and
  writing a provenance manifest.

- **`test/smoke/`**, an end-to-end check of stages 2, 3 and 5 against a synthetic dataset,
  run by CI and inside the image.

- **`test/e2e/`**, a full run through the real aligners against a committed 8.2 MB fixture
  (`make e2e`). The fixture is chromosome 17 of GRCm38 masked to `N` outside eight known
  junctions, so coordinates stay identical to the source assembly while the reference
  compresses to 1.5 MB. It recovers all eight junctions with the same supporting reads the
  full-depth run found there, and unlike the smoke test it exercises candidate nomination and
  all four Bowtie2 passes. FASTA is committed rather than indexes: STAR records the genome
  format version that wrote an index and refuses one it cannot read, so a committed index
  would turn an aligner upgrade into a confusing failure.

- **`ops/make_test_bundle.py`** and **`ops/make_reference_bundle.sh`**, which build the test
  and reference bundles reproducibly rather than as opaque archives, and
  **`ops/zenodo_upload.sh`**, which uploads to a Zenodo draft and never publishes. The
  reference packager refuses a STAR index with leftover temporary files, and records the
  contig naming convention and the STAR genome format version in the bundle manifest --
  the two things that silently waste hours when they are wrong.

- **Published bundles**: test data (Zenodo 22929133, 76 MB) and the GRCm38 reference
  (Zenodo 22929144, 27 GB), both Ensembl-named with no `chr` prefix.

- **`-M` minimum span, `-n` threads, `-m` Java heap, `-k` keep intermediates, `-L` legacy,
  `-V` skip validation, `-N` normalise strand.** `bowtie2` now runs with `--no-unal`.

- A floor on the construct arm, so a short-read library cannot produce arms too small to
  host the junction window.

### Notes

- The construct generator holds the genome in memory, which is why `-m` defaults to 20G;
  a mammalian genome needs roughly 4 GB of that before any structures are buffered.
- Structures whose anchors run past the end of a sequence are dropped rather than
  truncated, so junctions within one construct arm of a contig end are not called.

## 2.1.0

### Fixed - these change results

- **Junction span filter rejected most spanning reads.** The guard in
  `MDFilter.checkJunctionSpan` compared the read-relative junction offset against
  the read's absolute position in the construct, two different coordinate spaces.
  For 100 bp reads against an 85 bp construct arm it discarded every read aligning
  from roughly POS 41 onwards, including reads with generous overhang on both
  sides. The guard now requires the junction window to fit inside the alignment
  string, which is the intended "at least `-j/2` positions either side" rule.
  **Expect more supporting reads and more called backsplices than previous
  releases produced.** Re-run any analysis whose counts are being compared across
  versions.

- **The genome FASTA parser silently dropped sequences.** It discarded the first
  record in the file and never flushed the trailing records, so on a typical
  human FASTA `chr1` and the last chromosomes were absent from the genome map.
  Every structure on those chromosomes failed construct generation and was
  swallowed by a per-structure `catch`. The parser now streams the file and keeps
  every record.

- **Sequence names are taken from the first token of the FASTA header,** so
  Ensembl-style headers such as `>1 dna:chromosome chromosome:GRCh38:1:...` now
  resolve. A chromosome present in the structures but absent from the genome is
  reported as an error rather than silently dropping its structures, and the run
  fails when no construct at all can be generated.

- **`-C` (maximum backsplice span) was overridden by a hardcoded 1 Mb ceiling** in
  `ProcessShuffledCoordinates`, so values above 1 Mb had no effect. The bound now
  follows the flag. The lower bound, previously hardcoded at 50 bp, is exposed as `-M`.

- **Reverse complementation emitted the literal text `null`** for IUPAC ambiguity
  codes, corrupting any construct built over one. All IUPAC codes are now
  complemented properly and unknown symbols become `N`.

- **Construct FASTAs were opened in append mode,** so re-running into an existing
  working directory concatenated new constructs onto the previous run's. The same
  applied to `unique.sam` and the genomic filter outputs. All outputs are now truncated.

- **`-G` and `-T` were documented but not parsed,** so passing either was an
  "invalid option" error. Both now work. The transcriptomic-only path previously
  filtered an empty read set and returned nothing; it now loads `ptes.sam` first.
  The genomic-only path returned only the last batch when more than 5 million
  alignments were processed; it now reads back the complete `unique.sam`.

- **Optional SAM tags are located by prefix** (`MD:Z:`, `NM:i:`) rather than by
  column offset from the end of the line, which was fragile across aligners and
  tag orderings.

### Changed

- Reads accepted by the perfect-match shortcut (`NM:i:0`) without spanning the
  junction are now counted and reported separately in the log. This behaviour is
  unchanged, but it was previously invisible.
- Undefined percent identity and junction sequence columns are written as `NA`
  rather than `null`.
- Java components exit non-zero on failure, print usage on bad arguments and
  validate their inputs up front.
- Logs are mirrored to stderr as well as `run.log`, and report per-stage counts of
  records read, retained and rejected.

### Added

- Test suite (`test/`), run by `setup.sh`; the jar is only packaged when it passes.
- `PFv2.sh` validates every mandatory argument, checks that STAR, Bowtie2,
  samtools, Python 3 and Java are on `PATH`, and verifies that the FASTQ, genome
  FASTA, STAR index and both Bowtie2 indexes exist before doing any work.
- `PFv2.sh` runs under `set -euo pipefail` and stops at the first failing stage
  rather than cascading into confusing downstream errors.
- New flags: `-M` minimum span, `-n` threads, `-m` Java heap, `-k` keep intermediates.
- `-c` now defaults to the directory containing `PFv2.sh`.
- Junction-per-million normalisation is skipped, with a warning, when no
  junction-spanning reads were found, instead of dividing by zero.

### Removed

- `apache-lib/` (duplicate of `lib/`), `dist/`, `build/`, the committed `classes/`
  tree (stale - it predated the STAR rewrite), the duplicate `scripts/PFv2.jar`,
  and two 2017 Dropbox conflicted-copy files under `nbproject/`.

### Notes

- The Python wrappers no longer use `shell=True` and no longer depend on GNU
  `readlink -f`, so they work on macOS as well as Linux.
- `run_star.py` had roughly twenty argparse options that were never passed to
  STAR. They have been removed rather than wired up, so STAR's behaviour is
  unchanged; `--extra` passes arbitrary arguments through. `--chimOutType Junctions`
  is now set explicitly, pinning what was previously STAR's default.
- `PFv2.sh` re-execs itself under bash, so `sh PFv2.sh` still works on systems
  where `/bin/sh` is dash.
