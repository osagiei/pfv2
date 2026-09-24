'''
Author: Osagie Izuogu

Description: Fetches prebuilt PFv2 references from a Zenodo record.

             Zenodo is a good fit for this: a record gets a DOI, versions are immutable,
             and every file is published with its size and MD5, so a download can be
             verified rather than assumed. Presence is not completeness -- a half-written
             genome FASTA produces a run that looks fine and is not -- so nothing is
             treated as usable until its checksum matches.

             A concept DOI (the "all versions" DOI) resolves to the newest version, which
             is what you want for "give me the current reference". Pin the version DOI
             instead when a result has to be reproducible.

             Only the standard library is used, so this runs in any Python 3 without a
             virtualenv.

Date: 09/2026
'''

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import urllib.error
import urllib.request

ZENODO_API = "https://zenodo.org/api/records"
CHUNK = 1 << 20
USER_AGENT = "pfv2-fetch-references/1.0"

# Roles PFv2 needs, and how they are recognised in a record's file listing. A record may
# publish more, and anything unrecognised is still downloadable by name.
# Order matters: the first match wins, so the specific index patterns are tested before the
# broad sequence one. A file named "bowtie2-genome.tar.gz" contains an index, not a FASTA, and
# classifying it as genome_fasta made --only genome_fasta fetch the wrong archive.
ROLE_PATTERNS = (
    ("star_index", r"star[-_.]?index|star.*\.tar(\.gz)?(\.part\d+)?$"),
    ("bowtie2_genome", r"bowtie2[-_.]genome"),
    ("bowtie2_transcriptome", r"bowtie2[-_.](transcriptome|cdna)"),
    ("annotation_gtf", r"\.gtf(\.gz)?$"),
    ("genome_fasta", r"(^|[-_.])(sequence|genome)([-_.]|\.fa|$)|\.dna\."),
)


class FetchError(RuntimeError):
    pass


def record_id(reference: str) -> str:
    """Accepts a bare id, a Zenodo URL, or a DOI in either form."""
    reference = reference.strip()
    if reference.isdigit():
        return reference
    match = re.search(r"zenodo[./](\d+)", reference)
    if match:
        return match.group(1)
    raise FetchError(f"cannot read a Zenodo record id from '{reference}'")


def get_json(url: str) -> dict:
    request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.load(response)
    except urllib.error.HTTPError as exc:
        if exc.code in (404, 410):
            raise FetchError(
                f"Zenodo has no public record at {url}.\n"
                "       If you were given this id recently it is probably still an "
                "unpublished draft: the public API only serves published records, so the "
                "record has to be published (or you need the owner's token and the "
                "deposit API) before it can be fetched."
            ) from exc
        if exc.code in (401, 403):
            raise FetchError(
                f"Zenodo refused access to {url} ({exc.code}); the record may be "
                "restricted or embargoed."
            ) from exc
        raise FetchError(f"Zenodo returned {exc.code} for {url}") from exc
    except urllib.error.URLError as exc:
        raise FetchError(f"cannot reach Zenodo: {exc.reason}") from exc


def resolve(reference: str, *, latest: bool) -> dict:
    record = get_json(f"{ZENODO_API}/{record_id(reference)}")
    if latest:
        newest = record.get("links", {}).get("latest")
        if newest and newest.rstrip("/").split("/")[-1] != str(record.get("id")):
            record = get_json(newest)
    return record


def files_of(record: dict) -> list[dict]:
    """Normalises the file listing across the record shapes Zenodo has served."""
    entries = record.get("files")
    if isinstance(entries, dict):  # InvenioRDM style
        entries = list(entries.get("entries", {}).values())
    if not entries:
        return []

    listing = []
    for entry in entries:
        key = entry.get("key") or entry.get("filename")
        link = (entry.get("links", {}) or {}).get("self") or entry.get("link")
        checksum = entry.get("checksum") or ""
        if isinstance(checksum, dict):
            checksum = checksum.get("md5", "")
        if not key or not link:
            continue
        listing.append({
            "key": key,
            "url": link,
            "size": entry.get("size"),
            "md5": checksum.split(":")[-1] if checksum else None,
        })
    return listing


def role_of(key: str) -> str | None:
    lowered = key.lower()
    for role, pattern in ROLE_PATTERNS:
        if re.search(pattern, lowered):
            return role
    return None


def md5_of(path: str) -> str:
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            digest.update(block)
    return digest.hexdigest()


def download(entry: dict, dest_dir: str, *, force: bool) -> tuple[str, str]:
    """Fetches one file, verifying it against the published MD5.

    Returns (status, path) where status is one of skipped, fetched.
    """
    path = os.path.join(dest_dir, entry["key"])
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)

    if os.path.isfile(path) and not force:
        if entry["md5"] is None:
            # Nothing to verify against, so size is the only available signal.
            if entry["size"] is None or os.path.getsize(path) == entry["size"]:
                return "skipped", path
        elif md5_of(path) == entry["md5"]:
            return "skipped", path

    partial = path + ".part"
    request = urllib.request.Request(entry["url"], headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(request, timeout=120) as response, open(partial, "wb") as out:
            digest = hashlib.md5()
            while True:
                block = response.read(CHUNK)
                if not block:
                    break
                out.write(block)
                digest.update(block)
    except (urllib.error.HTTPError, urllib.error.URLError, OSError) as exc:
        if os.path.exists(partial):
            os.remove(partial)
        raise FetchError(f"{entry['key']}: download failed: {exc}") from exc

    actual = digest.hexdigest()
    if entry["md5"] and actual != entry["md5"]:
        os.remove(partial)
        raise FetchError(f"{entry['key']}: MD5 mismatch "
                         f"(expected {entry['md5']}, got {actual}); the file was not kept")
    if entry["size"] is not None and os.path.getsize(partial) != entry["size"]:
        os.remove(partial)
        raise FetchError(f"{entry['key']}: size mismatch "
                         f"(expected {entry['size']}, got {os.path.getsize(partial)})")

    # Only rename once the bytes are verified, so a consumer can never open a partial file.
    os.replace(partial, path)
    return "fetched", path


PART_RE = re.compile(r"^(?P<base>.+)\.part(?P<index>\d+)$")


def reassemble(dest: str, fetched: list[dict]) -> list[str]:
    """Concatenates any `<name>.partNN` group back into `<name>`.

    A file too large for one upload is published as parts. Rejoining them here means the
    caller never has to know a record was split, and the parts are removed only once the
    whole file is in place.
    """
    groups: dict[str, list[tuple[int, str]]] = {}
    for entry in fetched:
        match = PART_RE.match(entry["key"])
        if match:
            groups.setdefault(match.group("base"), []).append(
                (int(match.group("index")), entry["key"]))

    rebuilt = []
    for base, parts in sorted(groups.items()):
        parts.sort()
        target = os.path.join(dest, base)
        print(f"  joining {len(parts)} part(s) into {base}")
        tmp = target + ".joining"
        try:
            with open(tmp, "wb") as out:
                for _index, key in parts:
                    with open(os.path.join(dest, key), "rb") as fh:
                        while True:
                            block = fh.read(CHUNK)
                            if not block:
                                break
                            out.write(block)
        except OSError as exc:
            if os.path.exists(tmp):
                os.remove(tmp)
            print(f"  FAIL    could not join {base}: {exc}", file=sys.stderr)
            continue
        os.replace(tmp, target)
        # Only now that the whole file exists is it safe to drop the parts.
        for _index, key in parts:
            os.remove(os.path.join(dest, key))
        rebuilt.append(base)
    return rebuilt


def verify_whole(dest: str, rebuilt: list[str]) -> int:
    """Checks reassembled files against WHOLE-SHA256SUMS when the record publishes one."""
    sums = os.path.join(dest, "WHOLE-SHA256SUMS")
    if not rebuilt or not os.path.isfile(sums):
        return 0
    expected = {}
    with open(sums) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) == 2:
                expected[parts[1]] = parts[0]

    failures = 0
    for name in rebuilt:
        want = expected.get(name)
        if not want:
            continue
        digest = hashlib.sha256()
        with open(os.path.join(dest, name), "rb") as fh:
            for block in iter(lambda: fh.read(CHUNK), b""):
                digest.update(block)
        if digest.hexdigest() == want:
            print(f"  verified {name}")
        else:
            print(f"  FAIL    {name}: reassembled file does not match its published checksum",
                  file=sys.stderr)
            failures += 1
    return failures


def human(size) -> str:
    if size is None:
        return "?"
    value = float(size)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if value < 1024 or unit == "TiB":
            return f"{value:.1f}{unit}" if unit != "B" else f"{int(value)}B"
        value /= 1024
    return f"{value:.1f}TiB"


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Fetch prebuilt PFv2 references from a Zenodo record.")
    p.add_argument("record", help="Zenodo record id, DOI, or record URL")
    p.add_argument("-d", "--dest", default=".", help="destination directory (default: .)")
    p.add_argument("--list", action="store_true", dest="list_only",
                   help="list the record's files and exit without downloading")
    p.add_argument("--only", nargs="+", default=None,
                   help="fetch just these file names, or just these roles "
                        f"({', '.join(r for r, _ in ROLE_PATTERNS)})")
    p.add_argument("--latest", action="store_true",
                   help="follow a concept DOI to the newest version")
    p.add_argument("--force", action="store_true",
                   help="re-download even when a verified copy is already present")
    p.add_argument("--manifest", default=None,
                   help="write a provenance record of what was fetched to this path")
    return p.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    try:
        record = resolve(args.record, latest=args.latest)
    except FetchError as exc:
        print(f"ERROR {exc}", file=sys.stderr)
        return 1

    title = (record.get("metadata", {}) or {}).get("title", "(untitled)")
    doi = record.get("doi") or (record.get("metadata", {}) or {}).get("doi", "")
    listing = files_of(record)

    print(f"Zenodo record {record.get('id')}: {title}")
    if doi:
        print(f"  DOI {doi}")
    if not listing:
        print("  this record publishes no files (it may be restricted or embargoed)",
              file=sys.stderr)
        return 1

    for entry in listing:
        entry["role"] = role_of(entry["key"])

    if args.list_only:
        for entry in listing:
            print(f"  {entry['key']:<50} {human(entry['size']):>10}  {entry['role'] or ''}")
        return 0

    wanted = listing
    if args.only:
        requested = set(args.only)
        wanted = [e for e in listing if e["key"] in requested or e["role"] in requested]
        unmatched = requested - {e["key"] for e in wanted} - {e["role"] for e in wanted if e["role"]}
        if unmatched:
            print(f"ERROR no file or role in this record matches: {sorted(unmatched)}",
                  file=sys.stderr)
            return 1

    os.makedirs(args.dest, exist_ok=True)
    fetched = []
    failures = 0
    for entry in wanted:
        try:
            status, path = download(entry, args.dest, force=args.force)
        except FetchError as exc:
            print(f"  FAIL    {exc}", file=sys.stderr)
            failures += 1
            continue
        print(f"  {status:<7} {entry['key']} ({human(entry['size'])})")
        fetched.append({
            "key": entry["key"], "role": entry["role"], "path": os.path.abspath(path),
            "md5": entry["md5"], "size": entry["size"], "status": status,
        })

    rebuilt = reassemble(args.dest, fetched)
    failures += verify_whole(args.dest, rebuilt)

    if args.manifest:
        with open(args.manifest, "w") as fh:
            json.dump({
                "record_id": record.get("id"),
                "doi": doi,
                "title": title,
                "files": fetched,
            }, fh, indent=2)
            fh.write("\n")
        print(f"  manifest written to {args.manifest}")

    if failures:
        print(f"{failures} file(s) could not be fetched", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
