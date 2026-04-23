#!/usr/bin/env python3
"""
Merge chunked Illumina FASTQ.gz files into one R1 and one R2 per logical prefix.

Example input (one directory):
  CLP163_GACATGGT_L008_R1_001.fastq.gz ... CLP163_GACATGGT_L008_R1_005.fastq.gz
  CLP163_GACATGGT_L008_R2_001.fastq.gz ... CLP163_GACATGGT_L008_R2_005.fastq.gz

Output in the same directory (by default):
  CLP163_GACATGGT_L008_R1.fastq.gz
  CLP163_GACATGGT_L008_R2.fastq.gz

If index reads are present (I1/I2), they are merged separately, e.g. ``CLP713_I1.fastq.gz``,
not mixed into R1/R2.

Default mode: for each *immediate subdirectory* of ROOT, look for *.fastq.gz / *.fq.gz
and merge. Use --flat to merge files that sit directly under ROOT instead.

``--organize-flat`` only sees FASTQs sitting *directly* under ROOT (not inside sample
subfolders). For trees that are already ``ROOT/SampleID/*.fastq.gz``, skip that flag.

Supported naming patterns include:
- `<prefix>_R1_001.fastq.gz` / `<prefix>_R2_001.fastq.gz`
- `<prefix>_I1_001.fastq.gz` / `<prefix>_I2_001.fastq.gz`
- `<prefix>.1.fq.gz` / `<prefix>.2.fq.gz`
- and single-chunk variants without numeric chunk suffix.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import sys
from collections import defaultdict
from pathlib import Path

# Examples:
# CLP163_GACATGGT_L008_R1_001.fastq.gz -> prefix=CLP163_GACATGGT_L008, read=1, chunk=1
# CLP713_I2_001.fastq.gz               -> prefix=CLP713, read=2, chunk=1
# ABJ133.1.fq.gz                       -> prefix=ABJ133, read=1, chunk=1
PATTERNS = [
    re.compile(
        r"^(?P<prefix>.+)_(?P<tag>[RI])(?P<read>[12])_(?P<chunk>\d+)\.(?P<ext>fastq|fq)\.gz$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?P<prefix>.+)_(?P<tag>[RI])(?P<read>[12])\.(?P<ext>fastq|fq)\.gz$",
        re.IGNORECASE,
    ),
    re.compile(
        r"^(?P<prefix>.+)\.(?P<read>[12])(?:\.(?P<chunk>\d+))?\.(?P<ext>fastq|fq)\.gz$",
        re.IGNORECASE,
    ),
]


def parse_fastq_name(name: str) -> tuple[str, str, int, str] | None:
    """Return (prefix, read, chunk, ext) for supported Illumina naming conventions."""
    for pat in PATTERNS:
        m = pat.match(name)
        if not m:
            continue
        prefix = m.group("prefix")
        read = m.group("read")
        chunk = int(m.group("chunk")) if m.groupdict().get("chunk") else 1
        ext = m.group("ext").lower()
        return prefix, read, chunk, ext
    return None


def iter_fastq_gz(dirpath: Path) -> list[Path]:
    out: list[Path] = []
    for pat in ("*.fastq.gz", "*.fq.gz", "*.FASTQ.GZ", "*.FQ.GZ"):
        out.extend(dirpath.glob(pat))
    return sorted((p for p in set(out) if not p.name.startswith("._")), key=lambda p: p.name)


def group_chunks(paths: list[Path]) -> dict[tuple[str, str], list[tuple[int, Path]]]:
    """Map (prefix, '1'|'2') -> sorted list of (chunk_int, path)."""
    groups: dict[tuple[str, str], list[tuple[int, Path]]] = defaultdict(list)
    for p in paths:
        parsed = parse_fastq_name(p.name)
        if not parsed:
            continue
        prefix, read, chunk, _ = parsed
        groups[(prefix, read)].append((chunk, p))
    for key in groups:
        groups[key].sort(key=lambda t: (t[0], t[1].name))
    return groups


def is_gzip_file(path: Path) -> bool:
    """Check gzip magic bytes (1f 8b)."""
    with path.open("rb") as fh:
        return fh.read(2) == b"\x1f\x8b"


def concat_gzip(inputs: list[Path], output: Path) -> None:
    """
    Write a .fastq.gz output from possibly mixed inputs.

    - If an input is gzip, append raw bytes directly (equivalent to `cat a.gz b.gz > out.gz`).
    - If an input is not gzip, compress it as an additional gzip member on the fly.
    """
    with output.open("wb") as out_raw:
        for path in inputs:
            if is_gzip_file(path):
                with path.open("rb") as in_raw:
                    shutil.copyfileobj(in_raw, out_raw, length=1024 * 1024)
            else:
                print(
                    f"Warning: input is not gzip; compressing on the fly: {path}",
                    file=sys.stderr,
                )
                with path.open("rb") as in_raw:
                    with gzip.GzipFile(fileobj=out_raw, mode="wb", compresslevel=6) as gz_out:
                        shutil.copyfileobj(in_raw, gz_out, length=1024 * 1024)


def organize_flat_directory(
    root: Path,
    *,
    dry_run: bool,
    force: bool,
    verbose: bool,
) -> tuple[int, int]:
    """
    Move chunked FASTQs from a flat directory into per-prefix directories.

    Returns (moved_files, created_dirs).
    """
    paths = iter_fastq_gz(root)
    groups = group_chunks(paths)
    if not groups:
        if verbose:
            print(f"(no chunked R1/R2 files to organize in {root})", file=sys.stderr)
        return 0, 0

    moved = 0
    created_dirs = 0
    created_set: set[Path] = set()
    for (_, _), items in sorted(groups.items()):
        for _, src in items:
            parsed = parse_fastq_name(src.name)
            if not parsed:
                continue
            prefix, _, _, _ = parsed
            dst_dir = root / prefix
            dst = dst_dir / src.name

            if src.parent == dst_dir:
                continue
            if not dst_dir.exists():
                if dry_run:
                    if dst_dir not in created_set:
                        created_set.add(dst_dir)
                        created_dirs += 1
                else:
                    dst_dir.mkdir(parents=True, exist_ok=True)
                    if dst_dir not in created_set:
                        created_set.add(dst_dir)
                        created_dirs += 1

            if dst.exists() and not force:
                print(f"Skip move (exists): {dst}", file=sys.stderr)
                continue

            if verbose or dry_run:
                print(f"{'WOULD MOVE' if dry_run else 'MOVE'} {src} -> {dst}")
            if not dry_run:
                shutil.move(str(src), str(dst))
            moved += 1

    return moved, created_dirs


def merge_directory(
    dirpath: Path,
    *,
    dry_run: bool,
    force: bool,
    verbose: bool,
) -> int:
    """Merge chunks in dirpath. Returns number of output files written (or would write)."""
    paths = iter_fastq_gz(dirpath)
    groups = group_chunks(paths)
    if not groups:
        if verbose:
            print(f"  (no matching R1/R2 or I1/I2 files in {dirpath})", file=sys.stderr)
        return 0

    n_out = 0
    for (prefix, read) in sorted(groups.keys()):
        chunks = groups[(prefix, read)]
        ordered = [p for _, p in chunks]
        first = parse_fastq_name(ordered[0].name)
        ext = first[3] if first else "fastq"
        out_path = dirpath / f"{prefix}_R{read}.{ext}.gz"

        if out_path.exists() and not force:
            print(f"Skip (exists): {out_path}", file=sys.stderr)
            continue

        if verbose or dry_run:
            print(f"{'WOULD WRITE' if dry_run else 'WRITE'} {out_path}")
            for p in ordered:
                print(f"    <- {p.name}")

        if dry_run:
            n_out += 1
            continue

        concat_gzip(ordered, out_path)
        n_out += 1

    return n_out


def collect_target_dirs(root: Path, flat: bool) -> list[Path]:
    if flat:
        return [root.resolve()]
    subs = [p for p in root.iterdir() if p.is_dir()]
    if not subs:
        return [root.resolve()]
    return sorted(subs)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Merge Illumina FASTQ(.gz) chunks into single R1/R2 files per prefix. "
            "Supports R1/R2, I1/I2, and .1/.2 conventions."
        )
    )
    parser.add_argument(
        "root",
        type=Path,
        help="Top directory (e.g. HudsonAlpha_Dec2015_SongSparrow)",
    )
    parser.add_argument(
        "--flat",
        action="store_true",
        help="Merge FASTQs directly under ROOT instead of processing each immediate subdirectory",
    )
    parser.add_argument(
        "--organize-flat",
        action="store_true",
        help="For a flat ROOT, move chunked files into per-sample directories before merge",
    )
    parser.add_argument(
        "-n",
        "--dry-run",
        action="store_true",
        help="Print planned merges only",
    )
    parser.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite existing merged outputs",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print per-directory notes even when nothing matches",
    )
    args = parser.parse_args()

    root = args.root.expanduser().resolve()
    if not root.is_dir():
        print(f"Not a directory: {root}", file=sys.stderr)
        sys.exit(1)

    if args.organize_flat:
        moved, created_dirs = organize_flat_directory(
            root,
            dry_run=args.dry_run,
            force=args.force,
            verbose=args.verbose or args.dry_run,
        )
        if args.dry_run:
            print(
                f"Dry run: would move {moved} file(s) into {created_dirs} directory(ies).",
                file=sys.stderr,
            )
        else:
            print(
                f"Organized flat root: moved {moved} file(s) into {created_dirs} directory(ies).",
                file=sys.stderr,
            )

    targets = collect_target_dirs(root, args.flat)
    total = 0
    for d in targets:
        if args.verbose or not args.flat:
            print(f"== {d}", file=sys.stderr)
        total += merge_directory(
            d,
            dry_run=args.dry_run,
            force=args.force,
            verbose=args.verbose or args.dry_run,
        )

    if args.dry_run:
        print(f"Dry run: would create {total} merged file(s).", file=sys.stderr)
    else:
        print(f"Done: wrote {total} merged file(s).", file=sys.stderr)


if __name__ == "__main__":
    main()
