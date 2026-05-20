#!/usr/bin/env python3
"""
Finish repo reorganization:

1. Move remaining ``mitogenome_annotations/<sample>/`` trees into matching
   ``assembled_mitogenomes/<sample>/annotation/<suffix>/``.
2. Rename ``assembled_mitogenomes`` -> ``mitogenomes_output``.
3. Remove ``mitogenome_annotations/``.
"""
from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(Path(__file__).resolve().parent))

from mitogenome_paths import ASM_ROOT, LEGACY_ANNOT_ROOT  # noqa: E402

# When mitogenome_annotations name differs from assembly dir on disk.
ASSEMBLY_DIR_ALIASES: dict[str, str] = {
    "Mareca_americana_UAM_11919_KSW3039": "Mareca_americana_UAM_11919_JSW3039",
}


def infer_suffix(asm_dir: Path) -> str:
    novops = sorted(p.name[len("novoplasty_") :] for p in asm_dir.iterdir() if p.is_dir() and p.name.startswith("novoplasty_"))
    if len(novops) == 1:
        return novops[0]
    if len(novops) > 1:
        return novops[0]
    # fallback: last segment of directory name (often field #)
    parts = asm_dir.name.split("_")
    return parts[-1] if parts else "default"


def move_legacy_annotations(*, dry_run: bool) -> list[str]:
    if not LEGACY_ANNOT_ROOT.is_dir():
        return []

    moved: list[str] = []
    for child in sorted(LEGACY_ANNOT_ROOT.iterdir()):
        if not child.is_dir() or child.name.startswith(".") or child.name == "README.md":
            continue

        asm_name = ASSEMBLY_DIR_ALIASES.get(child.name, child.name)
        asm_dir = ASM_ROOT / asm_name
        if not asm_dir.is_dir():
            print(f"SKIP (no assembly dir): {child.name}")
            continue

        suffix = infer_suffix(asm_dir)
        dest_parent = asm_dir / "annotation" / suffix
        dest_mitoz = dest_parent / "mitoz"
        src_mitoz = child / "mitoz"

        if not src_mitoz.is_dir():
            print(f"SKIP (no mitoz/): {child.name}")
            continue

        if dest_mitoz.exists():
            print(f"MERGE {child.name} -> {dest_mitoz.relative_to(REPO)} (replace)")
            if not dry_run:
                shutil.rmtree(dest_mitoz)
        else:
            print(f"MOVE {child.name} -> {dest_mitoz.relative_to(REPO)}")

        if not dry_run:
            dest_parent.mkdir(parents=True, exist_ok=True)
            shutil.move(str(src_mitoz), str(dest_mitoz))
        moved.append(child.name)
    return moved


def rename_output_root(*, dry_run: bool) -> None:
    new_root = REPO / "mitogenomes_output"
    if new_root.exists():
        print(f"Already exists: {new_root}")
        return
    if not ASM_ROOT.is_dir():
        print(f"Missing: {ASM_ROOT}")
        return
    print(f"RENAME {ASM_ROOT.relative_to(REPO)} -> mitogenomes_output/")
    if not dry_run:
        shutil.move(str(ASM_ROOT), str(new_root))


def remove_legacy_root(*, dry_run: bool) -> None:
    if not LEGACY_ANNOT_ROOT.is_dir():
        return
    remaining = [p.name for p in LEGACY_ANNOT_ROOT.iterdir() if not p.name.startswith(".")]
    if remaining:
        print(f"Legacy dir not empty, will not remove: {remaining}")
        return
    print(f"REMOVE {LEGACY_ANNOT_ROOT.relative_to(REPO)}/")
    if not dry_run:
        shutil.rmtree(LEGACY_ANNOT_ROOT)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("-n", "--dry-run", action="store_true")
    args = p.parse_args()
    moved = move_legacy_annotations(dry_run=args.dry_run)
    print(f"\nMoved {len(moved)} legacy annotation tree(s).")
    rename_output_root(dry_run=args.dry_run)
    if not args.dry_run:
        # Refresh path after rename for remove step
        legacy = REPO / "mitogenome_annotations"
        if legacy.is_dir():
            readme = legacy / "README.md"
            if readme.is_file():
                readme.unlink()
            if not any(legacy.iterdir()):
                shutil.rmtree(legacy)
            else:
                print(f"Legacy dir not empty: {list(legacy.iterdir())}")
    else:
        remove_legacy_root(dry_run=True)


if __name__ == "__main__":
    main()
