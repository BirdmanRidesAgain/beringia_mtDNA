#!/usr/bin/env python3
"""
Batch-annotate mitogenome assemblies under ``mitogenomes_output/``.

For each ``novoplasty_*`` subdirectory, picks the best NOVOPlasty FASTA per
priority (Circularized > Option > Contigs > contigs_tmp) and runs
``scripts/annotate_mitogenome.py`` (MitoZ) inside the conda env ``mitoz-x64``.

Special handling:

* If a ``novoplasty_*`` subdir contains both ``Option_1_*.fasta`` and
  ``Option_2_*.fasta``, both options are annotated as separate runs (``__option1``
  / ``__option2`` suffixes).
* For sample dirs with multiple ``novoplasty_*`` subdirs (the merged trinomial
  Anas_crecca_* dirs), keep only the higher-tier subdir; on a tier tie keep all.
* Sample dirs with no ``novoplasty_*`` subdir, or with only zero-byte
  ``contigs_tmp_*.txt``, are recorded in ``mitogenomes_output/_batch_logs/skipped.tsv``
  and not annotated.

MitoZ outputs are written under each assembly directory:
``mitogenomes_output/<sample>/annotation/<suffix>/mitoz/``.

Examples::

    python scripts/batch_annotate_assemblies.py --dry-run
    python scripts/batch_annotate_assemblies.py --fresh-log
    python scripts/batch_annotate_assemblies.py --first 50 --last 80
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
from collections import Counter
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ASM_ROOT = REPO / "mitogenomes_output"
ANNOT_SCRIPT = REPO / "scripts/annotate_mitogenome.py"
BATCH_LOGS_DIR = ASM_ROOT / "_batch_logs"
LOG = BATCH_LOGS_DIR / "batch_annotate_assemblies.log"
SKIP_TSV = BATCH_LOGS_DIR / "skipped.tsv"

from mitogenome_paths import (  # noqa: E402
    annotation_done,
    clear_annotation_output,
    parse_sample_label,
)

CONDA_ENV = "mitoz-x64"

TIER_CIRC = 1
TIER_OPTION = 2
TIER_CONTIGS = 3
TIER_TMP = 4
TIER_NONE = 99
TIER_NAMES = {1: "circularized", 2: "option", 3: "contigs", 4: "contigs_tmp"}


def conda_executable() -> Path:
    """Prefer known install paths so a stale ``CONDA_EXE`` doesn't point at Miniforge."""
    for candidate in (
        Path("/opt/anaconda3/bin/conda"),
        Path("/opt/miniconda3/bin/conda"),
        Path.home() / "miniconda3/bin/conda",
        Path.home() / "anaconda3/bin/conda",
    ):
        if candidate.is_file():
            return candidate
    if exe := os.environ.get("CONDA_EXE"):
        p = Path(exe)
        if p.is_file():
            return p
    w = shutil.which("conda")
    if w:
        return Path(w)
    raise FileNotFoundError(
        "Could not find conda. Set CONDA_EXE to your conda binary, or install Anaconda/Miniconda."
    )


def conda_install_base() -> Path:
    """``conda info --base`` for the chosen ``conda`` binary, or ``CONDA_ROOT`` override."""
    if override := os.environ.get("CONDA_ROOT"):
        p = Path(override)
        if p.is_dir():
            return p
    r = subprocess.run(
        [str(conda_executable()), "info", "--base"],
        capture_output=True,
        text=True,
        check=True,
    )
    return Path(r.stdout.strip())


def mitoz_python() -> Path:
    """Interpreter for the ``mitoz-x64`` env (must exist)."""
    p = conda_install_base() / "envs" / CONDA_ENV / "bin" / "python"
    if not p.is_file():
        raise FileNotFoundError(
            f"Expected conda env {CONDA_ENV!r} at {p}. "
            f"Create it (e.g. `conda create -n mitoz-x64 -c bioconda mitoz`) "
            f"or set CONDA_ROOT to the right install location."
        )
    return p


def annotate_child_environ() -> dict[str, str]:
    """PATH prefix so ``mitoz`` resolves the same as ``conda activate mitoz-x64``."""
    env = os.environ.copy()
    bin_dir = mitoz_python().parent
    env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")
    return env


@dataclass
class AnnotationJob:
    sample_dir: Path
    novoplasty_dir: Path
    fasta: Path
    sample_label: str
    tier: int


@dataclass
class SkippedItem:
    path: Path
    reason: str


def _list_files(d: Path) -> list[Path]:
    return [p for p in d.iterdir() if p.is_file() and not p.name.startswith("._")]


def _classify(novop: Path) -> tuple[int, list[Path]]:
    """Return ``(tier, candidate_fastas)`` for one ``novoplasty_*`` subdir."""
    files = _list_files(novop)
    circ = sorted(p for p in files if p.name.startswith("Circularized_assembly_") and p.suffix == ".fasta")
    options = sorted(p for p in files if p.name.startswith("Option_") and p.suffix == ".fasta")
    contigs = sorted(p for p in files if p.name.startswith("Contigs_") and p.suffix == ".fasta")
    contigs_tmp = sorted(p for p in files if p.name.startswith("contigs_tmp_") and p.suffix == ".txt")
    if circ:
        return TIER_CIRC, circ
    if options:
        return TIER_OPTION, options
    if contigs:
        return TIER_CONTIGS, contigs
    if contigs_tmp:
        non_empty = [p for p in contigs_tmp if p.stat().st_size > 0]
        if non_empty:
            return TIER_TMP, non_empty
    return TIER_NONE, []


def _novoplasty_token(novop: Path) -> str:
    return novop.name[len("novoplasty_"):]


def _disambiguator(fasta: Path) -> str:
    """Short tag used to disambiguate multiple FASTAs from one subdir."""
    name = fasta.name
    if name.startswith("Circularized_assembly_"):
        idx = name[len("Circularized_assembly_"):].split("_", 1)[0]
        return f"circ{idx}"
    if name.startswith("Option_"):
        idx = name[len("Option_"):].split("_", 1)[0]
        return f"option{idx}"
    if name.startswith("Contigs_"):
        idx = name[len("Contigs_"):].split("_", 1)[0]
        return f"contigs{idx}"
    if name.startswith("contigs_tmp_"):
        return "tmp"
    return fasta.stem


def discover_jobs() -> tuple[list[AnnotationJob], list[SkippedItem]]:
    jobs: list[AnnotationJob] = []
    skipped: list[SkippedItem] = []

    sample_dirs = sorted(
        p for p in ASM_ROOT.iterdir()
        if p.is_dir() and not p.name.startswith(".") and p.name != "bbsplit_ref"
    )

    for sd in sample_dirs:
        novops = sorted(p for p in sd.iterdir() if p.is_dir() and p.name.startswith("novoplasty_"))
        if not novops:
            skipped.append(SkippedItem(sd, "no novoplasty_* subdirectory"))
            continue

        classified = [(n, *_classify(n)) for n in novops]

        # Record TIER_NONE subdirs (e.g. zero-byte contigs_tmp) as skipped.
        for n, t, _c in classified:
            if t == TIER_NONE:
                skipped.append(SkippedItem(n, "no annotatable FASTA (only zero-byte/non-FASTA outputs)"))

        annotatable = [(n, t, c) for n, t, c in classified if t != TIER_NONE]
        if not annotatable:
            continue

        # When a sample has multiple novoplasty_* subdirs (the merged Anas_crecca dirs),
        # keep the higher-tier (lower number); on tie keep all.
        if len(annotatable) > 1:
            best_tier = min(t for _n, t, _c in annotatable)
            kept = [(n, t, c) for n, t, c in annotatable if t == best_tier]
            for n, t, _c in annotatable:
                if t != best_tier:
                    skipped.append(
                        SkippedItem(
                            n,
                            f"superseded by sibling novoplasty_* subdir at higher tier "
                            f"(this tier={t}, kept tier={best_tier})",
                        )
                    )
        else:
            kept = annotatable

        for novop, tier, cands in kept:
            token = _novoplasty_token(novop)
            multi = len(cands) > 1
            for fasta in cands:
                if multi:
                    sample_label = f"{sd.name}__{token}__{_disambiguator(fasta)}"
                else:
                    sample_label = f"{sd.name}__{token}"
                jobs.append(AnnotationJob(sd, novop, fasta, sample_label, tier))

    jobs.sort(key=lambda j: j.sample_label)
    return jobs, skipped


def output_already_done(job: AnnotationJob) -> bool:
    """Resumable check: True if a MitoZ ``*.result/`` directory exists for this job."""
    return annotation_done(job.sample_dir, job.sample_label)


def annotate_cmd(job: AnnotationJob, *, threads: int, clade: str) -> list[str]:
    _, suffix = parse_sample_label(job.sample_label)
    annot_root = job.sample_dir / "annotation"
    return [
        str(mitoz_python()),
        str(ANNOT_SCRIPT),
        "-i", str(job.fasta),
        "-s", suffix,
        "-o", str(annot_root),
        "--threads", str(threads),
        "--mitoz-clade", clade,
    ]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch MitoZ annotation for assembled mitogenomes.")
    p.add_argument("--first", type=int, default=1, help="1-based first job index (default: 1).")
    p.add_argument("--last", type=int, default=None, help="1-based last job index inclusive (default: last job).")
    p.add_argument("--fresh-log", action="store_true", help="Overwrite the batch log instead of appending.")
    p.add_argument("--dry-run", action="store_true", help="List jobs and skips without running anything.")
    p.add_argument(
        "--threads",
        type=int,
        default=max(1, (os.cpu_count() or 1) // 2),
        help="Threads passed through to MitoZ.",
    )
    p.add_argument("--mitoz-clade", default="Chordata", help="MitoZ clade (default: Chordata).")
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-run a job even when a MitoZ *.result/ directory already exists.",
    )
    return p.parse_args()


def write_skip_tsv(skipped: list[SkippedItem]) -> None:
    BATCH_LOGS_DIR.mkdir(parents=True, exist_ok=True)
    with SKIP_TSV.open("w", encoding="utf-8") as f:
        f.write("path\treason\n")
        for s in skipped:
            try:
                rel = s.path.relative_to(REPO)
            except ValueError:
                rel = s.path
            f.write(f"{rel}\t{s.reason}\n")


def main() -> None:
    args = parse_args()
    jobs, skipped = discover_jobs()
    write_skip_tsv(skipped)

    total = len(jobs)
    if total == 0:
        print("No annotation jobs discovered.")
        print(f"Skip log: {SKIP_TSV}")
        return

    first = max(1, args.first)
    last = total if args.last is None else min(args.last, total)
    if first > last or first > total:
        raise SystemExit(f"Invalid range --first {args.first} --last {args.last} (have {total} jobs).")
    subset = jobs[first - 1 : last]

    if args.dry_run:
        tier_counts = Counter(j.tier for j in jobs)
        print(f"Discovered {total} annotation jobs and {len(skipped)} skipped item(s).")
        print(f"Skip log: {SKIP_TSV}")
        for t in sorted(tier_counts):
            print(f"  tier {t} ({TIER_NAMES.get(t,'?')}): {tier_counts[t]} job(s)")
        print()
        print("First 20 jobs in subset:")
        for offset, j in enumerate(subset[:20]):
            idx = first + offset
            print(f"  [{idx}/{total}] tier={j.tier} {j.sample_label}  <- {j.fasta.relative_to(ASM_ROOT)}")
        if len(subset) > 20:
            print("  ...")
            for offset, j in list(enumerate(subset))[-5:]:
                idx = first + offset
                print(f"  [{idx}/{total}] tier={j.tier} {j.sample_label}  <- {j.fasta.relative_to(ASM_ROOT)}")
        print()
        print(f"Skipped items ({len(skipped)}):")
        for s in skipped:
            try:
                rel = s.path.relative_to(REPO)
            except ValueError:
                rel = s.path
            print(f"  {rel}\t{s.reason}")
        return

    log_mode = "w" if args.fresh_log else "a"
    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(f"Total annotation jobs: {total}\n")
        log.write(f"Skipped items: {len(skipped)} (see {SKIP_TSV.relative_to(REPO)})\n")
        log.write(
            f"\n--- START {datetime.now(timezone.utc).isoformat()}: "
            f"jobs [{first}/{total}]–[{last}/{total}] ---\n"
        )
        log.flush()
        for offset, job in enumerate(subset):
            idx = first + offset
            tname = TIER_NAMES.get(job.tier, "?")
            log.write(f"\n[{idx}/{total}] tier={job.tier}({tname}) {job.sample_label}\n")
            log.write(f"  fasta: {job.fasta.relative_to(ASM_ROOT)}\n")
            log.flush()
            if (not args.force) and output_already_done(job):
                log.write("  Skipped: existing mitoz output detected (use --force to re-run).\n")
                log.flush()
                continue
            t0 = datetime.now(timezone.utc)
            r = subprocess.run(
                annotate_cmd(job, threads=args.threads, clade=args.mitoz_clade),
                env=annotate_child_environ(),
            )
            dt = (datetime.now(timezone.utc) - t0).total_seconds()
            log.write(f"  Exit code: {r.returncode}  duration: {dt:.1f}s\n")
            log.flush()


if __name__ == "__main__":
    main()
