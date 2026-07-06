#!/usr/bin/env python3
"""Move loose SRR*.fastq.gz files into SRR*/ subdirs and cross-reference Table 1."""

from __future__ import annotations

import re
import shutil
from collections import defaultdict
from pathlib import Path

RAW_DATA = Path("/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/raw_data")
CHECKLIST = Path(__file__).resolve().parent.parent / "mitogenome_assembly_checklist.md"


def parse_table1_srr_map(path: Path) -> dict[str, list[dict]]:
    rows: dict[str, list[dict]] = defaultdict(list)
    in_table = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("## Table 1."):
            in_table = True
            continue
        if in_table and line.startswith("## "):
            break
        if not in_table or not line.startswith("| *"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 7:
            continue
        for srr in re.findall(r"SRR\d+", parts[3]):
            rows[srr].append(
                {
                    "species": re.sub(r"^\*|\*$", "", parts[0]),
                    "catalog": parts[1],
                    "field": parts[2],
                    "old_dir": parts[3],
                    "new_dir": parts[4].strip("`"),
                }
            )
    return dict(rows)


def loose_srr_files(raw_root: Path) -> dict[str, list[Path]]:
    grouped: dict[str, list[Path]] = defaultdict(list)
    for path in sorted(raw_root.glob("SRR*.fastq.gz")):
        if not path.is_file():
            continue
        match = re.match(r"(SRR\d+)", path.name)
        if match:
            grouped[match.group(1)].append(path)
    return dict(grouped)


def organize(raw_root: Path, dry_run: bool = False) -> dict[str, Path]:
    moved: dict[str, Path] = {}
    for srr, files in loose_srr_files(raw_root).items():
        dest = raw_root / srr
        if not dry_run:
            dest.mkdir(exist_ok=True)
            for src in files:
                target = dest / src.name
                if target.exists():
                    raise FileExistsError(f"Refusing to overwrite {target}")
                shutil.move(str(src), str(target))
        moved[srr] = dest
    return moved


def main() -> None:
    if not RAW_DATA.is_dir():
        raise SystemExit(f"Raw data directory not found: {RAW_DATA}")

    srr_map = parse_table1_srr_map(CHECKLIST)
    moved = organize(RAW_DATA, dry_run=False)
    remaining = list(RAW_DATA.glob("SRR*.fastq.gz"))

    print(f"Organized {len(moved)} SRR accessions into subdirectories under {RAW_DATA}")
    if remaining:
        print(f"WARNING: {len(remaining)} loose SRR fastq files remain at top level")

    matched = sorted(set(moved) & set(srr_map))
    unmatched_loose = sorted(set(moved) - set(srr_map))
    unmatched_checklist = sorted(set(srr_map) - set(moved))

    print(f"\nCross-reference with Table 1 Old directory name ({len(srr_map)} SRR refs):")
    print(f"  Matched (moved + in checklist): {len(matched)}")
    print(f"  Moved but not in checklist:   {len(unmatched_loose)}")
    for srr in unmatched_loose:
        print(f"    {srr}/")
    print(f"  In checklist but not moved:   {len(unmatched_checklist)}")

    print("\nMatched SRR -> Table 1 rows:")
    for srr in matched:
        for row in srr_map[srr]:
            print(
                f"  {srr}/ -> {row['species']} | {row['catalog']} | {row['field']} "
                f"| {row['new_dir']}"
            )


if __name__ == "__main__":
    main()
