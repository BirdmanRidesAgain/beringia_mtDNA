"""Shared helpers for 2018 resequencing reassembly/reannotation (93 samples)."""
from __future__ import annotations

import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CHECKLIST = REPO / "mitogenomes_assembled_checklist.md"
FASTQ_ROOT = Path("/Volumes/Seagate Por/KAC_mtDNA_reassembly/raw_data")
OUT_ROOT = REPO / "mitogenomes_output"
REFSEQ = REPO / "test/refseq/Taeniopygia_guttata_NC_007897.1_mitogenome.fa"


def parse_reassembly_samples(checklist: Path = CHECKLIST) -> list[tuple[str, str]]:
    """Return (field_num, new_directory_name) for Table 1 rows with multiple datasets? True."""
    text = checklist.read_text(encoding="utf-8")
    start = text.find("## Table 1.")
    end = text.find("\n## Table 2.")
    if start < 0 or end < 0:
        raise ValueError("Could not locate Table 1 in checklist.")
    block = text[start:end]
    samples: list[tuple[str, str]] = []
    for line in block.splitlines():
        if not line.startswith("|") or line.startswith("| -") or "Field num" in line:
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 7:
            continue
        if parts[5] != "True":
            continue
        field = parts[2]
        new_dir = parts[4].strip("`").strip("/")
        if field and new_dir:
            samples.append((field, new_dir))
    return samples


def clear_assembly_artifacts(out_dir: Path, field: str) -> list[str]:
    """Remove prior fastp/NOVOPlasty outputs so a re-run replaces them."""
    removed: list[str] = []
    for name in (f"novoplasty_{field}", f"trim_{field}"):
        target = out_dir / name
        if target.exists():
            shutil.rmtree(target)
            removed.append(name)
    return removed


def clear_annotation_artifacts(out_dir_name: str, asm_root: Path = OUT_ROOT) -> list[str]:
    """Remove all ``annotation/`` outputs under an assembly directory."""
    from mitogenome_paths import clear_all_annotations

    sample_dir = asm_root / out_dir_name
    return clear_all_annotations(sample_dir)
