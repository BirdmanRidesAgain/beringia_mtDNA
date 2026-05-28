"""Path conventions for assemblies and co-located MitoZ annotations."""
from __future__ import annotations

import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ASM_ROOT = REPO / "mitogenomes_output"
ANNOTATION_FOLDER = "annotation"
BATCH_LOGS_DIR = ASM_ROOT / "_batch_logs"

# NOVOPlasty output tiers (lower = better). Matches batch_annotate_assemblies.py.
TIER_CIRC = 1
TIER_OPTION = 2
TIER_CONTIGS = 3
TIER_TMP = 4
TIER_NONE = 99
TIER_NAMES = {1: "circularized", 2: "option", 3: "contigs", 4: "contigs_tmp", 99: "none"}


def _list_novoplasty_files(novop: Path) -> list[Path]:
    return [p for p in novop.iterdir() if p.is_file() and not p.name.startswith("._")]


def classify_novoplasty_dir(novop: Path) -> int:
    """Best NOVOPlasty tier in one ``novoplasty_*`` directory."""
    files = _list_novoplasty_files(novop)
    if any(p.name.startswith("Circularized_assembly_") and p.suffix == ".fasta" for p in files):
        return TIER_CIRC
    if any(p.name.startswith("Option_") and p.suffix == ".fasta" for p in files):
        return TIER_OPTION
    if any(p.name.startswith("Contigs_") and p.suffix == ".fasta" for p in files):
        return TIER_CONTIGS
    if any(
        p.name.startswith("contigs_tmp_") and p.suffix == ".txt" and p.stat().st_size > 0
        for p in files
    ):
        return TIER_TMP
    return TIER_NONE


def tier_to_assembly_level(tier: int) -> str:
    """
    Map best tier to CSV ``assembly_level``.

    - ``success``: circularized FASTA (tier 1)
    - ``partial_success``: Option or Contigs FASTA (tiers 2–3)
    - ``failed``: only contigs_tmp, empty output, or no novoplasty dir (tiers 4, 99)
    """
    if tier == TIER_CIRC:
        return "success"
    if tier in (TIER_OPTION, TIER_CONTIGS):
        return "partial_success"
    return "failed"


def assembly_level_for_sample_dir(sample_dir: Path) -> str:
    """Best tier across all ``novoplasty_*`` subdirectories under a sample folder."""
    if not sample_dir.is_dir():
        return "failed"
    novops = sorted(
        p for p in sample_dir.iterdir() if p.is_dir() and p.name.startswith("novoplasty_")
    )
    if not novops:
        return "failed"
    best = TIER_NONE
    for novop in novops:
        best = min(best, classify_novoplasty_dir(novop))
    return tier_to_assembly_level(best)


def _annotation_files_under(suffix_dir: Path) -> list[Path]:
    """All regular files under one ``annotation/<suffix>/`` tree."""
    return [
        f
        for f in suffix_dir.rglob("*")
        if f.is_file() and not f.name.startswith("._")
    ]


def annotation_status_for_suffix_dir(suffix_dir: Path) -> str:
    """
    Classify one annotation attempt (``annotation/<suffix>/``).

    - ``not_attempted``: no annotation output files
    - ``failed``: files exist but every file is 0 bytes
    - ``success``: at least one non-zero output file
    """
    files = _annotation_files_under(suffix_dir)
    if not files:
        return "not_attempted"
    if all(f.stat().st_size == 0 for f in files):
        return "failed"
    return "success"


def annotation_success_for_sample_dir(sample_dir: Path) -> str:
    """
    Best annotation outcome for a sample directory.

    Scans every ``annotation/<suffix>/`` subtree. If any suffix run succeeded,
    the sample is ``success``; otherwise ``failed`` when outputs exist but are
    all empty; otherwise ``not_attempted``.
    """
    if not sample_dir.is_dir():
        return "not_attempted"
    ann = sample_dir / ANNOTATION_FOLDER
    if not ann.is_dir():
        return "not_attempted"

    statuses: list[str] = []
    for suffix_dir in sorted(ann.iterdir()):
        if not suffix_dir.is_dir() or suffix_dir.name.startswith("._"):
            continue
        statuses.append(annotation_status_for_suffix_dir(suffix_dir))

    if not statuses:
        return "not_attempted"
    if "success" in statuses:
        return "success"
    if "failed" in statuses:
        return "failed"
    return "not_attempted"


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
