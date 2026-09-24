#!/usr/bin/env bash
#
# Sets a deposition's metadata and, with --publish, publishes it.
#
# Zenodo will not publish a record without creators, so this fills in the descriptive
# metadata that ops/zenodo_upload.sh deliberately leaves minimal, then optionally publishes.
#
# Publishing is irreversible: a published record's files can never be changed, only
# superseded by a new version under the same concept DOI. --publish therefore has to be
# asked for explicitly; without it this only updates the draft's metadata.
#
#   bash ops/zenodo_publish.sh --deposition 22929133 --kind testdata
#   bash ops/zenodo_publish.sh --deposition 22929144 --kind reference --publish
#
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

DEPOSITION=""; KIND=""; PUBLISH=false; HOST="https://zenodo.org"
CREATOR_NAME="${ZENODO_CREATOR:-Izuogu, Osagie G.}"
CREATOR_AFFIL="${ZENODO_AFFILIATION:-}"
ORCID="${ZENODO_ORCID:-}"

while [ $# -gt 0 ]; do
  case "$1" in
    --deposition) DEPOSITION="$2"; shift 2 ;;
    --kind) KIND="$2"; shift 2 ;;
    --publish) PUBLISH=true; shift ;;
    --sandbox) HOST="https://sandbox.zenodo.org"; shift ;;
    -h|--help) sed -n '2,15p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
step() { printf '>>> %s\n' "$*"; }

[ -n "$DEPOSITION" ] || die "--deposition is required"
case "$KIND" in testdata|reference) ;; *) die "--kind must be testdata or reference" ;; esac

if [ -z "${ZENODO_TOKEN:-}" ] && [ -z "${ZENODO:-}" ]; then
  ENV_FILE="${ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/.env}"
  [ -f "$ENV_FILE" ] && { set -a; . "$ENV_FILE"; set +a; }
fi
ZENODO_TOKEN="${ZENODO_TOKEN:-${ZENODO:-}}"
[ -n "$ZENODO_TOKEN" ] || die "no token: set ZENODO_TOKEN or ZENODO, or put one in .env"

# Through --config on stdin, never argv, so `ps` cannot show it.
auth_config() { printf 'header = "Authorization: Bearer %s"\n' "$ZENODO_TOKEN"; }

RESPONSE="$(mktemp)"; chmod 600 "$RESPONSE"
trap 'rm -f "$RESPONSE" "$META_FILE"' EXIT INT TERM
META_FILE="$(mktemp)"; chmod 600 "$META_FILE"

VERSION="$(grep -m1 '^readonly VERSION=' "$(dirname "${BASH_SOURCE[0]}")/../PFv2.sh" | cut -d'"' -f2)"

python3 - "$KIND" "$VERSION" "$CREATOR_NAME" "$CREATOR_AFFIL" "$ORCID" > "$META_FILE" <<'PY'
import json, sys
kind, version, name, affiliation, orcid = sys.argv[1:6]

creator = {"name": name}
if affiliation:
    creator["affiliation"] = affiliation
if orcid:
    creator["orcid"] = orcid

common = {
    "upload_type": "dataset",
    "access_right": "open",
    "license": "mit",
    "creators": [creator],
    "version": version,
    "keywords": ["circRNA", "circular RNA", "backsplice junction",
                 "post-transcriptional exon shuffling", "PTES", "PTESFinder",
                 "RNA-seq", "GRCm38", "Mus musculus"],
    "related_identifiers": [
        {"identifier": "10.1186/s12859-016-0881-4",
         "relation": "isSupplementTo", "scheme": "doi"},
        {"identifier": "10.1038/s41592-023-01944-6",
         "relation": "isReferencedBy", "scheme": "doi"},
    ],
}

if kind == "testdata":
    common.update({
        "title": ("PTESFinder v2 test data: GRCm38 chromosome 17 (Ensembl 102) "
                  "with RNA-seq reads"),
        "description": (
            "<p>A small end-to-end test dataset for <strong>PTESFinder v2 (PFv2)</strong>, "
            "an annotation-free method for identifying post-transcriptional exon shuffling "
            "(PTES) / backsplice junctions from RNA-seq data.</p>"
            "<p>Contents:</p><ul>"
            "<li><code>genome.fa</code> - chromosome 17 of GRCm38 (Ensembl release 102)</li>"
            "<li><code>transcriptome.fa</code> - the Ensembl 102 cDNA records on that chromosome</li>"
            "<li><code>bl6-chr17.fq.gz</code> - single-end reads sufficient to rediscover the "
            "backsplice junctions PFv2 calls there, comprising the chimeric reads STAR reported "
            "within the chromosome, the reads supporting the calls, and a background sample</li>"
            "<li><code>expected/pf-structures.bed</code> - the junctions a full-depth run called "
            "on this chromosome, for comparison</li></ul>"
            "<p><strong>Sequence naming:</strong> contigs are Ensembl-named "
            "(<code>1</code>...<code>19</code>, <code>X</code>, <code>Y</code>) with no "
            "<code>chr</code> prefix. This is not interchangeable with a UCSC mm10 reference; "
            "mixing the two conventions produces an empty result rather than an error.</p>"
            "<p>Reads derive from C57BL/6 mouse RNA-seq. See <code>MANIFEST.json</code> for "
            "provenance and checksums.</p>"),
    })
else:
    common.update({
        "title": ("PTESFinder v2 reference: GRCm38 (Ensembl 102) genome, transcriptome, "
                  "STAR and Bowtie2 indexes"),
        "description": (
            "<p>A prebuilt reference set for <strong>PTESFinder v2 (PFv2)</strong>, an "
            "annotation-free method for identifying post-transcriptional exon shuffling "
            "(PTES) / backsplice junctions from RNA-seq data. Building these indexes takes "
            "about an hour and ~32 GB of RAM, so they are published rather than rebuilt per "
            "site; using identical indexes also makes results comparable between groups.</p>"
            "<p>Contents (each a gzipped tar archive):</p><ul>"
            "<li><code>sequence</code> - GRCm38 genome FASTA and Ensembl 102 cDNA</li>"
            "<li><code>star-index</code> - STAR genome index, split into 4 GB parts</li>"
            "<li><code>bowtie2-genome</code> - Bowtie2 genome index</li>"
            "<li><code>bowtie2-transcriptome</code> - Bowtie2 transcriptome index</li></ul>"
            "<p><strong>Reassembling the STAR index:</strong> the archive is published as "
            "<code>.partNN</code> files because a single upload of that size is rejected. Rejoin "
            "with <code>cat GRCm38-ensembl102-star-index.tar.gz.part* &gt; "
            "GRCm38-ensembl102-star-index.tar.gz</code> and check the result against "
            "<code>WHOLE-SHA256SUMS</code>. PFv2's <code>ptesfinder fetch-references</code> does "
            "this automatically.</p>"
            "<p><strong>Sequence naming:</strong> contigs are Ensembl-named "
            "(<code>1</code>...<code>19</code>, <code>X</code>, <code>Y</code>) with no "
            "<code>chr</code> prefix, and are not interchangeable with a UCSC mm10 reference.</p>"
            "<p><strong>Aligner versions:</strong> the STAR index has genome format version "
            "2.7.4a and <code>sjdbOverhang</code> 149, built with STAR 2.7.11b. STAR refuses an "
            "index whose genome format it cannot read, so this index is usable only by STAR "
            "releases sharing that format. Bowtie2 indexes were built with Bowtie2 2.5.4. See "
            "<code>MANIFEST.json</code> for full provenance.</p>"),
    })

print(json.dumps({"metadata": common}))
PY

step "Updating metadata on deposition ${DEPOSITION}"
code=$(auth_config | curl -sS --config - -o "$RESPONSE" -w '%{http_code}' \
  -X PUT -H "Content-Type: application/json" \
  --data-binary "@${META_FILE}" \
  "${HOST}/api/deposit/depositions/${DEPOSITION}")
if [ "$code" != "200" ]; then
  die "metadata update failed with HTTP ${code}: $(head -c 500 "$RESPONSE")"
fi
python3 -c "
import json,sys
d=json.load(open('$RESPONSE'))
m=d.get('metadata',{})
print('    title  :', m.get('title','')[:78])
print('    creators:', ', '.join(c.get('name','') for c in m.get('creators',[])))
print('    licence:', m.get('license'), '| version:', m.get('version'), '| files:', len(d.get('files',[])))
"

if ! $PUBLISH; then
  step "Metadata set. Not published (pass --publish to publish)."
  printf '    %s/deposit/%s\n' "$HOST" "$DEPOSITION"
  exit 0
fi

step "Publishing - this cannot be undone"
code=$(auth_config | curl -sS --config - -o "$RESPONSE" -w '%{http_code}' \
  -X POST "${HOST}/api/deposit/depositions/${DEPOSITION}/actions/publish")
if [ "$code" != "202" ] && [ "$code" != "200" ]; then
  die "publish failed with HTTP ${code}: $(head -c 500 "$RESPONSE")"
fi
python3 -c "
import json
d=json.load(open('$RESPONSE'))
print('    state :', d.get('state'))
print('    DOI   :', d.get('doi'))
print('    record:', (d.get('links') or {}).get('record_html'))
"
