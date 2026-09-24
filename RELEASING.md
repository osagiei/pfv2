# Releasing

Steps that need a person, either because they publish something irreversible or because they
need credentials. Everything else is `make build && make smoke && make e2e`.

## The STAR index is not on Zenodo, by decision

The published reference record carries the genome FASTA, the Ensembl cDNA and both Bowtie2
indexes. It does not carry the STAR index, and that is now the intended shape rather than an
outstanding task.

Repeated upload attempts failed: Zenodo's gateway returned 502 on a 23 GB archive, then on
4 GB parts, then on 2 GB parts, and a throughput probe measured 877 KB/s with a 100 MB PUT
timing out at 90 seconds. At that rate no usefully sized part completes inside the proxy's
window. Draft 22939780 exists as a new version and is clean, if the index is ever worth
adding when throughput recovers; `ops/zenodo_upload.sh` skips files already uploaded, so
resuming costs only what is missing.

Shipping the index was always the weaker option anyway. A STAR index records the genome
format version that wrote it and STAR refuses one it cannot read, so a published index ties
users to a range of STAR releases and silently expires. `ops/build_indexes.sh` builds it from
the published genome in about an hour, and always matches the STAR actually installed.

## Rotate the Zenodo token before release

The token currently in `.env` should be revoked and replaced. During the first reference
upload it was passed to `curl` with `-H`, which puts it in the process argument list where
`ps` exposes it to every other user on a shared machine. The upload ran on a shared host, so
treat it as disclosed. `ops/zenodo_upload.sh` now feeds the header through `curl --config` on
stdin instead, so argv never holds it - but that does not undo the earlier exposure.

Revoke at https://zenodo.org/account/settings/applications/tokens/ and put the replacement in
`.env`, which is gitignored.

## Remaining steps

- [x] **Zenodo records published.**
      - test data: [10.5281/zenodo.22929132](https://doi.org/10.5281/zenodo.22929132) - complete
      - reference: [10.5281/zenodo.22929143](https://doi.org/10.5281/zenodo.22929143) - version 22929144 was
        published without the STAR index, because that upload had failed with a 502 before
        publication. A published record's files cannot be changed, so draft **22939780** is a
        new version carrying the STAR index as 4 GB parts. Publish it once the upload
        completes; the concept DOI then resolves to it and the documented commands need no
        change.

- [x] **Container image published.** `conidiobolus/pfv2:2.2.1` and `:latest`, public,
      `linux/amd64` and `linux/arm64`. Verified by an anonymous pull on a host with no Docker
      credentials, and both variants pass `ptesfinder selftest`.

      Re-publish with:

      ```bash
      docker buildx create --use --name multiarch --driver docker-container
      make image
      docker buildx imagetools inspect conidiobolus/pfv2:2.2.1
      ```

## Release checklist

- [ ] `make build` - 101 unit checks pass and the jar is repackaged
- [ ] `make smoke` - stages 2, 3 and 5 against the synthetic dataset
- [ ] `make e2e` - full run through the real aligners against the committed fixture
- [ ] Commit the rebuilt `PFv2.jar`; CI fails when its classes do not match `src/`
- [ ] `CHANGELOG.md` has an entry for the version in `PFv2.sh`
- [ ] Tag `<version>` (bare semver, matching the existing `2.0.0` tag), which triggers the
      image workflow; a `v`-prefixed tag works too
- [ ] Set `DOCKERHUB_USERNAME` and `DOCKERHUB_TOKEN` as repository secrets, and
      `DOCKERHUB_NAMESPACE` as a variable if it is not `conidiobolus`

## Notes for the upstream pull request

- The filter changes alter results. `CHANGELOG.md` separates them from the rest and `-L`
  restores the previous behaviour, which is what a reviewer comparing against published
  numbers will need.
- The strand column changed meaning. `-A` restores the old one. The rationale and the
  citation are in the README under *The strand column*.
- `src/bio/igm/utils/annotate/` and `src/bio/igm/utils/init/` are v1 utilities that the
  pipeline does not invoke. They still compile and are left in place; whether upstream wants
  them is a decision for that repository, not this branch.

## Things deliberately not automated

- `ops/zenodo_upload.sh` creates drafts and never publishes. A published record's files
  cannot be replaced, only superseded by a new version.
- No workflow pushes an image on a branch build. Only a `v*` tag does.
