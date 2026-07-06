#!/usr/bin/env python3
"""Rename/merge raw_data directories to checklist New directory name values."""

from __future__ import annotations

import re
import shutil
from collections import defaultdict
from pathlib import Path

RAW_DATA = Path("/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/raw_data")
CHECKLIST = Path(__file__).resolve().parent.parent / "mitogenome_assembly_checklist.md"


def normalize_dir(value: str) -> str:
    value = value.strip().strip("`").strip("/")
    if "/" in value:
        value = value.rsplit("/", 1)[-1]
    return value


def split_old_dirs(value: str) -> list[str]:
    parts = re.split(r"\s*,\s*", value)
    out: list[str] = []
    for part in parts:
        part = part.strip().strip("`").strip("/")
        if not part:
            continue
        if "/" in part:
            part = part.rsplit("/", 1)[-1]
        out.append(part)
    return out


def parse_checklist_rows(path: Path) -> list[tuple[str, list[str]]]:
    rows: list[tuple[str, list[str]]] = []
    in_directory_table = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if "Old directory name | New directory name" in line:
            in_directory_table = True
            continue
        if line.startswith("## "):
            if "Table 3." in line:
                in_directory_table = False
            continue
        if not in_directory_table or not line.startswith("|"):
            continue
        if line.startswith("| - |"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 7:
            continue
        new_dir = normalize_dir(parts[4])
        if not new_dir:
            continue
        old_dirs = split_old_dirs(parts[3])
        rows.append((new_dir, old_dirs))
    return rows


def list_top_level_dirs(root: Path) -> set[str]:
    return {p.name for p in root.iterdir() if p.is_dir() and not p.name.startswith(".")}


def find_sources(target: str, old_dirs: list[str], disk_dirs: set[str]) -> set[str]:
    sources: set[str] = set()
    old_set = set(old_dirs)
    for name in disk_dirs:
        if name == target or name in old_set or name.startswith(f"{target}_"):
            sources.add(name)
    return sources


def merge_path(src: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    for item in src.iterdir():
        if item.name.startswith("."):
            continue
        target = dest / item.name
        if item.is_dir():
            if target.exists() and target.is_dir():
                merge_path(item, target)
                item.rmdir()
            elif target.exists():
                fallback = dest / f"{src.name}__{item.name}"
                shutil.move(str(item), str(fallback))
            else:
                shutil.move(str(item), str(target))
        elif target.exists():
            fallback = dest / f"{src.name}__{item.name}"
            shutil.move(str(item), str(fallback))
        else:
            shutil.move(str(item), str(target))
    if src.exists() and not any(src.iterdir()):
        src.rmdir()


def consolidate_group(root: Path, target: str, sources: set[str]) -> list[str]:
    actions: list[str] = []
    source_names = sorted(name for name in sources if (root / name).is_dir())
    if not source_names:
        return actions

    dest = root / target

    if dest.is_dir():
        pass
    else:
        primary = root / source_names[0]
        primary.rename(dest)
        actions.append(f"renamed {source_names[0]} -> {target}")
        source_names = source_names[1:]

    for name in source_names:
        if name == target:
            continue
        src = root / name
        if not src.exists():
            continue
        merge_path(src, dest)
        if src.exists():
            shutil.rmtree(src)
        actions.append(f"merged {name} -> {target}")

    return actions


def main() -> None:
    if not RAW_DATA.is_dir():
        raise SystemExit(f"Raw data directory not found: {RAW_DATA}")

    rows = parse_checklist_rows(CHECKLIST)
    grouped: dict[str, list[str]] = defaultdict(list)
    for new_dir, old_dirs in rows:
        grouped[new_dir].extend(old_dirs)

    disk_dirs = list_top_level_dirs(RAW_DATA)
    all_actions: list[str] = []
    unmatched_targets: list[str] = []

    for target in sorted(grouped):
        old_dirs = sorted(set(grouped[target]))
        sources = find_sources(target, old_dirs, disk_dirs)
        if not sources:
            unmatched_targets.append(target)
            continue
        actions = consolidate_group(RAW_DATA, target, sources)
        all_actions.extend(actions)
        disk_dirs = list_top_level_dirs(RAW_DATA)

    print(f"Processed {len(grouped)} checklist new directory names")
    print(f"Actions taken: {len(all_actions)}")
    for action in all_actions[:50]:
        print(f"  {action}")
    if len(all_actions) > 50:
        print(f"  ... and {len(all_actions) - 50} more")
    if unmatched_targets:
        print(f"No matching source directories for {len(unmatched_targets)} targets")
        for name in unmatched_targets[:20]:
            print(f"  {name}")
        if len(unmatched_targets) > 20:
            print(f"  ... and {len(unmatched_targets) - 20} more")


if __name__ == "__main__":
    main()
