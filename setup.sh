#!/usr/bin/env bash
#
# Compiles PFv2, runs the test suite and packages PFv2.jar.
#
if [ -z "${BASH_VERSION:-}" ]; then
  exec bash "$0" "$@"
fi

set -euo pipefail

readonly CODE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly CLASSES="${CODE}/build/classes"
readonly TEST_CLASSES="${CODE}/build/test-classes"

SKIP_TESTS=false
for arg in "$@"; do
  case "$arg" in
    --skip-tests) SKIP_TESTS=true ;;
    -h|--help)
      echo "Usage: bash setup.sh [--skip-tests]"
      exit 0 ;;
    *)
      echo "Unknown option: $arg" >&2
      exit 2 ;;
  esac
done

log() { printf '>>> %s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

command -v javac >/dev/null 2>&1 || die "javac not found on PATH; install a JDK (16 or newer)"
command -v jar   >/dev/null 2>&1 || die "jar not found on PATH; install a JDK (16 or newer)"

LIBS="$(find "${CODE}/lib" -name '*.jar' -not -name '._*' 2>/dev/null | tr '\n' ':')"
[[ -n "$LIBS" ]] || die "No jars found in ${CODE}/lib"

rm -rf "$CLASSES" "$TEST_CLASSES"
mkdir -p "$CLASSES" "$TEST_CLASSES"

log "Compiling PFv2 sources"
javac -Xlint:all,-serial -encoding UTF-8 \
  -cp "${LIBS}" \
  -d "$CLASSES" \
  "${CODE}"/src/bio/igm/entities/*.java \
  "${CODE}"/src/bio/igm/*/*/*.java

if $SKIP_TESTS; then
  log "Skipping tests (--skip-tests)"
else
  log "Compiling tests"
  javac -encoding UTF-8 \
    -cp "${LIBS}${CLASSES}" \
    -d "$TEST_CLASSES" \
    "${CODE}"/test/bio/igm/*.java \
    "${CODE}"/test/bio/igm/*/*.java \
    "${CODE}"/test/bio/igm/*/*/*.java

  log "Running tests"
  java -cp "${LIBS}${CLASSES}:${TEST_CLASSES}" bio.igm.TestRunner \
    || die "Tests failed; PFv2.jar was not rebuilt"
fi

log "Packaging PFv2.jar"
jar --create --file "${CODE}/PFv2.jar" --manifest "${CODE}/manifest.mf" -C "$CLASSES" .

log "Built ${CODE}/PFv2.jar"
log "Run the pipeline with: bash ${CODE}/PFv2.sh -h"
