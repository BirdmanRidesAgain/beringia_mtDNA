#!/usr/bin/env python3
"""Apply Anas/Calidris/Mareca naming conventions across repo and mitogenomes_output."""

from __future__ import annotations

import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
ASM = ROOT / "mitogenomes_output"

# Top-level sample directory renames (old -> new)
DIR_RENAMES = [
    ("Anas_crecca_carolinensis_", "Anas_carolinensis_"),
    ("Anas_crecca_crecca_", "Anas_crecca_"),
    ("Calidria_alpina_pacifica_", "Calidris_alpina_pacifica_"),
]

# String replacements in file contents (order matters)
TEXT_REPLACEMENTS = [
    ("Anas crecca carolinensis", "Anas carolinensis"),
    ("Anas crecca crecca", "Anas crecca"),
    ("Anas_crecca_carolinensis_", "Anas_carolinensis_"),
    ("Anas_crecca_crecca_", "Anas_crecca_"),
    ("Calidria alpina pacifica", "Calidris alpina pacifica"),
    ("Calidria_alpina_pacifica_", "Calidris_alpina_pacifica_"),
    ("Calidria", "Calidris"),  # after longer patterns
    ("Mareca_americana_UAM_11919_JSW3039", "Mareca_americana_UAM_11919_KSW3039"),
    ("UAM11919,JSW3039", "UAM11919,KSW3039"),  # depths tsv
]

TEXT_FILES = [
    ROOT / "mitogenomes_assembled_checklist.md",
    ROOT / "mitogenome_assembly_results.csv",
    ROOT / "table_S1.tsv",
    ROOT / "scripts" / "sync_assembly_results_from_checklist.py",
    ROOT / "scripts" / "build_mitogenome_assembly_results.py",
    ROOT / "scripts" / "finalize_mitogenomes_output.py",
]

TEXT_EXTENSIONS = {".md", ".csv", ".tsv", ".py", ".txt", ".log", ".json", ".html", ".conf", ".njs", ".tbl", ".out"}


def rename_top_level_dirs() -> list[tuple[str, str]]:
    done = []
    if not ASM.is_dir():
        return done
    names = sorted(p.name for p in ASM.iterdir() if p.is_dir() and not p.name.startswith("_"))
    for old_prefix, new_prefix in DIR_RENAMES:
        for name in names:
            if not name.startswith(old_prefix):
                continue
            new_name = new_prefix + name[len(old_prefix) :]
            old_path = ASM / name
            new_path = ASM / new_name
            if new_path.exists():
                print(f"SKIP (target exists): {name} -> {new_name}")
                continue
            old_path.rename(new_path)
            print(f"RENAMED DIR: {name} -> {new_name}")
            done.append((name, new_name))
            names = [new_name if n == name else n for n in names]
    return done


def rename_paths_under_asm() -> None:
    """Rename files/dirs under mitogenomes_output whose names still contain old tokens."""
    tokens = [
        ("Anas_crecca_carolinensis_", "Anas_carolinensis_"),
        ("Anas_crecca_crecca_", "Anas_crecca_"),
        ("Calidria_alpina_pacifica_", "Calidris_alpina_pacifica_"),
        ("Calidria_", "Calidris_"),
        ("Mareca_americana_UAM_11919_JSW3039", "Mareca_americana_UAM_11919_KSW3039"),
    ]
    for _ in range(20):  # multiple passes for nested renames
        changed = False
        paths = sorted(ASM.rglob("*"), key=lambda p: len(p.parts), reverse=True)
        for path in paths:
            name = path.name
            new_name = name
            for old, new in tokens:
                new_name = new_name.replace(old, new)
            if new_name == name:
                continue
            target = path.with_name(new_name)
            if target.exists():
                continue
            path.rename(target)
            changed = True
        if not changed:
            break


def patch_text_file(path: Path) -> bool:
    if not path.is_file():
        return False
    text = path.read_text(encoding="utf-8", errors="replace")
    original = text
    for old, new in TEXT_REPLACEMENTS:
        text = text.replace(old, new)
    if text == original:
        return False
    path.write_text(text, encoding="utf-8")
    return True


def patch_repo_text_files() -> None:
    for path in TEXT_FILES:
        if patch_text_file(path):
            print(f"PATCHED: {path.relative_to(ROOT)}")

    for path in ASM.rglob("*"):
        if path.suffix.lower() not in TEXT_EXTENSIONS:
            continue
        if patch_text_file(path):
            print(f"PATCHED: {path.relative_to(ROOT)}")


def main() -> None:
    rename_top_level_dirs()
    rename_paths_under_asm()
    patch_repo_text_files()
    print("Done.")


if __name__ == "__main__":
    main()
