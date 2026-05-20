#!/usr/bin/env python3
"""
Historical migration script (completed). Moved MitoZ outputs from
``mitogenome_annotations/`` into ``mitogenomes_output/<assembly_dir>/annotation/<suffix>/mitoz/``.

Example::

    python scripts/migrate_annotations_to_assemblies.py --dry-run
    python scripts/migrate_annotations_to_assemblies.py
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import sys

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from mitogenome_paths import ASM_ROOT, BATCH_LOGS_DIR, parse_sample_label  # noqa: E402

LEGACY_ANNOT_ROOT = Path(__file__).resolve().parent.parent / "mitogenome_annotations"

LOG_SUFFIXES = {".log", ".tsv"}
SKIP_NAMES = {".DS_Store"}


def is_batch_artifact(path: Path) -> bool:
    return path.is_file() and path.suffix in LOG_SUFFIXES


def migrate(*, dry_run: bool) -> None:
    if not LEGACY_ANNOT_ROOT.is_dir():
        print(f"No legacy directory: {LEGACY_ANNOT_ROOT}")
        return

    BATCH_LOGS_DIR.mkdir(parents=True, exist_ok=True)

    moved = 0
    merged = 0
    skipped = 0
    unmapped: list[str] = []
    conflicts: list[str] = []

    for child in sorted(LEGACY_ANNOT_ROOT.iterdir()):
        if child.name.startswith("._") or child.name in SKIP_NAMES:
            continue
        if is_batch_artifact(child):
            dest = BATCH_LOGS_DIR / child.name
            if dest.exists() and not dry_run:
                dest.unlink()
            if dry_run:
                print(f"LOG  {child.name} -> _batch_logs/{child.name}")
            else:
                shutil.move(str(child), str(dest))
            moved += 1
            continue
        if not child.is_dir():
            continue

        label = child.name
        asm_name, suffix = parse_sample_label(label)
        asm_dir = ASM_ROOT / asm_name
        if not asm_dir.is_dir():
            unmapped.append(label)
            continue

        dest = asm_dir / "annotation" / suffix
        if dest.exists():
            # Prefer keeping existing destination if it already has mitoz results.
            dest_mitoz = dest / "mitoz"
            src_mitoz = child / "mitoz"
            if dest_mitoz.is_dir() and any(dest_mitoz.rglob("*.result")):
                conflicts.append(f"{label} (dest already has .result)")
                skipped += 1
                continue
            if dry_run:
                print(f"MERGE {label} -> {dest.relative_to(REPO)} (replace existing)")
            else:
                shutil.rmtree(dest)
                shutil.move(str(child), str(dest))
            merged += 1
            continue

        if dry_run:
            print(f"MOVE {label} -> {dest.relative_to(REPO)}")
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(child), str(dest))
        moved += 1

    # Relocate stray mitoz/ directly under assembly dirs into annotation/<suffix>/
    for asm_dir in sorted(ASM_ROOT.iterdir()):
        if not asm_dir.is_dir() or asm_dir.name.startswith(".") or asm_dir.name == "_batch_logs":
            continue
        stray = asm_dir / "mitoz"
        if not stray.is_dir():
            continue
        # Guess suffix from sole novoplasty_* subdir
        novops = [p for p in asm_dir.iterdir() if p.is_dir() and p.name.startswith("novoplasty_")]
        if len(novops) != 1:
            conflicts.append(f"{asm_dir.name}/mitoz (ambiguous novoplasty subdirs)")
            continue
        suffix = novops[0].name[len("novoplasty_") :]
        dest = asm_dir / "annotation" / suffix
        if dest.exists():
            conflicts.append(f"{asm_dir.name}/mitoz (annotation/{suffix} exists)")
            continue
        if dry_run:
            print(f"STRAY {asm_dir.name}/mitoz -> annotation/{suffix}/mitoz")
        else:
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(stray), str(dest / "mitoz"))
        moved += 1

    print()
    print(f"Annotation dirs moved: {moved}")
    print(f"Merged into existing dest: {merged}")
    print(f"Skipped (conflict): {skipped}")
    print(f"Batch logs moved: (included above)")
    print(f"Unmapped (no assembly dir): {len(unmapped)}")
    for name in unmapped[:20]:
        print(f"  {name}")
    if len(unmapped) > 20:
        print(f"  ... and {len(unmapped) - 20} more")
    if conflicts:
        print(f"Conflicts/manual review: {len(conflicts)}")
        for c in conflicts[:15]:
            print(f"  {c}")

    readme = LEGACY_ANNOT_ROOT / "README.md"
    text = (
        "# Legacy annotation directory\n\n"
        "MitoZ outputs now live next to each assembly under\n"
        "`mitogenomes_output/<sample>/annotation/<suffix>/mitoz/`.\n\n"
        "Batch logs are in `mitogenomes_output/_batch_logs/`.\n"
    )
    if not dry_run:
        readme.write_text(text, encoding="utf-8")


def main() -> None:
    p = argparse.ArgumentParser(description="Migrate mitogenome_annotations into assembly dirs.")
    p.add_argument("-n", "--dry-run", action="store_true")
    args = p.parse_args()
    migrate(dry_run=args.dry_run)


if __name__ == "__main__":
    main()
