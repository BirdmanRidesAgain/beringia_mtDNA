"""Path conventions for assemblies and co-located MitoZ annotations."""
from __future__ import annotations

import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ASM_ROOT = REPO / "mitogenomes_output"
ANNOTATION_FOLDER = "annotation"
BATCH_LOGS_DIR = ASM_ROOT / "_batch_logs"


def parse_sample_label(sample_label: str) -> tuple[str, str]:
    """
    Split a batch annotation label into assembly directory name and suffix.

    ``Motacilla_..._ABJ133__ABJ133`` -> (``Motacilla_..._ABJ133``, ``ABJ133``)
    """
    if "__" in sample_label:
        asm_name, suffix = sample_label.split("__", 1)
        return asm_name, suffix
    return sample_label, "default"


def annotation_output_dir(sample_dir: Path, sample_label: str) -> Path:
    """Directory passed to annotate_mitogenome as ``output_root / sample_name``."""
    _, suffix = parse_sample_label(sample_label)
    return sample_dir / ANNOTATION_FOLDER / suffix


def mitoz_dir(sample_dir: Path, sample_label: str) -> Path:
    return annotation_output_dir(sample_dir, sample_label) / "mitoz"


def annotation_done(sample_dir: Path, sample_label: str) -> bool:
    """True if MitoZ wrote a ``*.result/`` directory for this job."""
    mz = mitoz_dir(sample_dir, sample_label)
    if not mz.is_dir():
        return False
    prefix = f"{sample_label}."
    if (mz / f"{sample_label}.result").is_dir():
        return True
    for child in mz.iterdir():
        if child.is_dir() and child.name.endswith(".result"):
            return True
        if child.is_dir() and child.name.startswith(prefix) and child.name.endswith(".result"):
            return True
    return False


def clear_annotation_output(sample_dir: Path, sample_label: str) -> list[str]:
    """Remove one annotation output tree under ``sample_dir/annotation/``."""
    target = annotation_output_dir(sample_dir, sample_label)
    if not target.exists():
        return []
    shutil.rmtree(target)
    return [str(target.relative_to(sample_dir))]


def clear_all_annotations(sample_dir: Path) -> list[str]:
    """Remove every ``annotation/`` subtree for an assembly directory."""
    ann = sample_dir / ANNOTATION_FOLDER
    if not ann.is_dir():
        return []
    removed = []
    for child in list(ann.iterdir()):
        if child.is_dir() and not child.name.startswith("._"):
            shutil.rmtree(child)
            removed.append(f"{ANNOTATION_FOLDER}/{child.name}")
    if ann.is_dir() and not any(ann.iterdir()):
        ann.rmdir()
    return removed
