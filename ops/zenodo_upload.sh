#!/usr/bin/env bash
#
# Uploads a bundle directory to a Zenodo deposition.
#
# This script never publishes. It creates or fills a draft and prints the link; reviewing
# the metadata and pressing publish is a deliberate, irreversible act that belongs to a
# person, because a published record's files cannot be changed afterwards -- only
# superseded by a new version.
#
# The token is read from ZENODO_TOKEN or ZENODO, or from a .env file holding either, and is
# never written to disk, a log, or the manifest. Create one at
# https://zenodo.org/account/settings/applications/tokens/new/ with the deposit:write and
# deposit:actions scopes.
#
# Keep the .env out of version control. A token committed to git stays in the history after
# a force-push and has to be revoked, not deleted.
#
#   export ZENODO_TOKEN=...
#   bash ops/zenodo_upload.sh --dir /scratch/refbundle \
#     --title 'PFv2 reference: GRCm38 (Ensembl 102), STAR and Bowtie2 indexes'
#
#   # add files to an existing draft
#   bash ops/zenodo_upload.sh --dir /scratch/more --deposition 1234567
#
#   # add files to an ALREADY PUBLISHED record, which needs a new version: a published
#   # record's files can never be changed, only superseded
#   bash ops/zenodo_upload.sh --dir /scratch/more --new-version-of 1234567
#
#   --sandbox targets sandbox.zenodo.org, which is the right place to rehearse a 30 GB upload.
#
if [ -z "${BASH_VERSION:-}" ]; then exec bash "$0" "$@"; fi
set -euo pipefail

DIR=""; TITLE=""; DEPOSITION=""; HOST="https://zenodo.org"
DESCRIPTION=""; DRY_RUN=false; NEW_VERSION_OF=""

while [ $# -gt 0 ]; do
  case "$1" in
    --dir) DIR="$2"; shift 2 ;;
    --title) TITLE="$2"; shift 2 ;;
    --description) DESCRIPTION="$2"; shift 2 ;;
    --deposition) DEPOSITION="$2"; shift 2 ;;
    --new-version-of) NEW_VERSION_OF="$2"; shift 2 ;;
    --sandbox) HOST="https://sandbox.zenodo.org"; shift ;;
    --dry-run) DRY_RUN=true; shift ;;
    -h|--help) sed -n '2,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 2 ;;
  esac
done

die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }
step() { printf '>>> %s\n' "$*"; }

[ -n "$DIR" ] && [ -d "$DIR" ] || die "--dir must be an existing directory"
command -v curl >/dev/null 2>&1 || die "curl is required"
command -v python3 >/dev/null 2>&1 || die "python3 is required"

# Built with a read loop rather than mapfile, which needs bash 4; macOS ships 3.2.
FILES=()
while IFS= read -r line; do
  FILES+=("$line")
done < <(find "$DIR" -maxdepth 1 -type f ! -name '.*' | sort)
[ "${#FILES[@]}" -gt 0 ] || die "no files to upload in $DIR"

# Zenodo file keys are flat, so a bundle with subdirectories cannot be uploaded as loose
# files without silently losing them. Refuse rather than upload an incomplete record.
SUBDIRS=$(find "$DIR" -mindepth 1 -maxdepth 1 -type d ! -name '.*' | wc -l | tr -d ' ')
if [ "$SUBDIRS" -gt 0 ]; then
  die "${DIR} contains subdirectories, which Zenodo cannot store as loose files.
       Pack the bundle first, e.g.
         tar -C $(dirname "$DIR") -czf ${DIR}.tar.gz $(basename "$DIR")
       and upload the directory holding that archive."
fi

TOTAL=$(du -sh "$DIR" | cut -f1)
step "Uploading ${#FILES[@]} file(s), ${TOTAL} total, to ${HOST}"
for f in "${FILES[@]}"; do printf '    %8s  %s\n' "$(du -h "$f" | cut -f1)" "$(basename "$f")"; done

# Zenodo's default quota is 50 GB per record. A larger set needs a quota increase requested
# from Zenodo support before the upload, not after it fails partway through.
# `du -sb` is GNU only; sum the sizes portably instead.
BYTES=$(find "$DIR" -maxdepth 1 -type f ! -name '.*' -exec wc -c {} \; 2>/dev/null \
        | awk '{s += $1} END {print s + 0}')
if [ "${BYTES:-0}" -gt 50000000000 ]; then
  printf 'WARNING: %s exceeds the 50 GB default record quota; request an increase first\n' "$TOTAL" >&2
fi

$DRY_RUN && { step "Dry run: nothing was sent"; exit 0; }

# Fall back to a .env beside the repository, so the token does not have to be exported by
# hand every time. Sourced in a subshell-safe way with `set -a`, and never echoed.
if [ -z "${ZENODO_TOKEN:-}" ] && [ -z "${ZENODO:-}" ]; then
  ENV_FILE="${ENV_FILE:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/.env}"
  if [ -f "$ENV_FILE" ]; then
    step "Reading the token from $(basename "$ENV_FILE")"
    set -a
    # shellcheck disable=SC1090
    . "$ENV_FILE"
    set +a
  fi
fi
ZENODO_TOKEN="${ZENODO_TOKEN:-${ZENODO:-}}"
[ -n "$ZENODO_TOKEN" ] || die "no token: set ZENODO_TOKEN or ZENODO, or put one in .env"

# Responses can carry account details, so they land in a private temp file that is removed
# however the script exits.
RESPONSE_FILE="$(mktemp)"
chmod 600 "$RESPONSE_FILE"
trap 'rm -f "$RESPONSE_FILE"' EXIT INT TERM

# The token is fed to curl through --config on stdin. Passing it with -H would put it in
# curl's argv, where `ps` exposes it to every user on the machine for the whole upload -- and
# a multi-gigabyte upload is a long window.
auth_config() {
  printf 'header = "Authorization: Bearer %s"\n' "$ZENODO_TOKEN"
}

api() {
  local method="$1" path="$2"; shift 2
  auth_config | curl -sS --config - -X "$method" "$@" "${HOST}/api${path}"
}

if [ -n "$NEW_VERSION_OF" ]; then
  step "Creating a new version of published record ${NEW_VERSION_OF}"
  # A published record is immutable. Adding a file means a new version, which carries its own
  # DOI while the concept DOI keeps resolving to whichever version is newest.
  code=$(auth_config | curl -sS --config - -o "$RESPONSE_FILE" -w '%{http_code}' \
    -X POST "${HOST}/api/deposit/depositions/${NEW_VERSION_OF}/actions/newversion")
  if [ "$code" != "201" ] && [ "$code" != "200" ]; then
    die "could not create a new version of ${NEW_VERSION_OF} (HTTP ${code}): $(head -c 400 "$RESPONSE_FILE")"
  fi
  DEPOSITION=$(python3 -c "
import json,sys
d=json.load(open('$RESPONSE_FILE'))
draft=(d.get('links') or {}).get('latest_draft','')
print(draft.rstrip('/').split('/')[-1] if draft else '')
")
  [ -n "$DEPOSITION" ] || die "new version created but its draft id could not be read"
  step "New draft deposition ${DEPOSITION}"
fi

if [ -z "$DEPOSITION" ]; then
  [ -n "$TITLE" ] || die "--title is required when creating a new deposition"
  step "Creating a draft deposition"
  META=$(python3 - "$TITLE" "$DESCRIPTION" <<'PY'
import json, sys
title, description = sys.argv[1], sys.argv[2] or sys.argv[1]
print(json.dumps({"metadata": {
    "title": title,
    "upload_type": "dataset",
    "description": description,
}}))
PY
)
  RESPONSE=$(api POST /deposit/depositions -H "Content-Type: application/json" -d "$META")
  DEPOSITION=$(printf '%s' "$RESPONSE" | python3 -c "import json,sys; print(json.load(sys.stdin).get('id',''))" 2>/dev/null || true)
  [ -n "$DEPOSITION" ] || die "could not create a deposition: $(printf '%s' "$RESPONSE" | head -c 400)"
  step "Draft deposition ${DEPOSITION}"
fi

BUCKET=$(api GET "/deposit/depositions/${DEPOSITION}" \
  | python3 -c "import json,sys; print((json.load(sys.stdin).get('links') or {}).get('bucket',''))" 2>/dev/null || true)
[ -n "$BUCKET" ] || die "could not read the bucket URL for deposition ${DEPOSITION}"

for f in "${FILES[@]}"; do
  name="$(basename "$f")"
  step "Uploading ${name}"
  # The bucket API streams a single PUT per file, which is what makes multi-gigabyte parts
  # workable; the older /files form buffers and falls over.
  code=$(auth_config | curl -sS --config - -o "$RESPONSE_FILE" -w '%{http_code}' \
    -X PUT --upload-file "$f" \
    "${BUCKET}/${name}")
  if [ "$code" != "200" ] && [ "$code" != "201" ]; then
    printf 'ERROR: upload of %s failed with HTTP %s: %s\n' \
      "$name" "$code" "$(head -c 300 "$RESPONSE_FILE")" >&2
    exit 1
  fi
  : > "$RESPONSE_FILE"
done

step "Done. The draft is NOT published."
printf '    Review and publish at: %s/deposit/%s\n' "$HOST" "$DEPOSITION"
printf '    Once published, wire the record into the README and test/e2e/run.sh:\n'
printf '      PFV2_TEST_RECORD=%s bash test/e2e/run.sh\n' "$DEPOSITION"
