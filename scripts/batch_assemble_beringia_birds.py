#!/usr/bin/env python3
"""
Stage and assemble mitogenomes for Table 1 samples with raw reads in ~/Downloads/beringia_birds.

1. Symlink (or copy) paired FASTQs into ``data/beringia_birds/<directory_name>/`` using
   canonical ``<directory_name>_R{1,2}.fastq.gz`` names (``data/`` is gitignored).
2. Run ``assemble_mitogenome.py`` into ``mitogenomes_output/<directory_name>/``.

Example::

    python scripts/batch_assemble_beringia_birds.py --dry-run
    python scripts/batch_assemble_beringia_birds.py --stage-only
    python scripts/batch_assemble_beringia_birds.py --fresh-log --update-checklist
    python scripts/batch_assemble_beringia_birds.py --field KSW2408 --verbose
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from batch_assemble_beringia2018 import (  # noqa: E402
    CHECKLIST,
    REFSEQ,
    assemble_child_environ,
    assemble_cmd,
    mitogenome_python,
    novoplasty_done,
)

REPO = Path(__file__).resolve().parent.parent
DOWNLOADS_ROOT = Path.home() / "Downloads/beringia_birds"
DATA_ROOT = REPO / "data/beringia_birds"
OUT_ROOT = REPO / "mitogenomes_output"
LOG = OUT_ROOT / "_batch_logs/batch_assemble_beringia_birds.log"

OLD_DIR_PREFIX = "beringia_birds/"
STAGED_DIR_PREFIX = "data/beringia_birds/"


@dataclass(frozen=True)
class BeringiaBirdSample:
    species: str
    catalog_num: str
    field_num: str
    directory_name: str


def parse_beringia_birds_samples(checklist: Path) -> list[BeringiaBirdSample]:
    """Table 1 rows whose old directory points at the Downloads ingest folder."""
    samples: list[BeringiaBirdSample] = []
    in_table = False
    for line in checklist.read_text(encoding="utf-8").splitlines():
        if line.startswith("| Species | Catalog"):
            in_table = True
            continue
        if not in_table:
            continue
        if not line.startswith("|"):
            if samples:
                break
            continue
        if line.startswith("| - |"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 7:
            continue
        old_dir = parts[3].strip("`").strip("/")
        if OLD_DIR_PREFIX not in old_dir and STAGED_DIR_PREFIX not in old_dir:
            continue
        directory_name = parts[4].strip("`").strip("/")
        if not directory_name:
            continue
        samples.append(
            BeringiaBirdSample(
                species=re.sub(r"^\*|\*$", "", parts[0]).strip(),
                catalog_num=parts[1],
                field_num=parts[2],
                directory_name=directory_name,
            )
        )
    return samples


def find_read_pair_in_dir(sample_dir: Path) -> tuple[Path, Path]:
    """Return one R1/R2 pair under *sample_dir* (Illumina _R1/_R2 or _1/_2)."""
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"Sample directory missing: {sample_dir}")

    for r1_pat, r2_pat in (
        ("*_R1.fastq.gz", "*_R2.fastq.gz"),
        ("*_R1.fq.gz", "*_R2.fq.gz"),
        ("*_1.fq.gz", "*_2.fq.gz"),
        ("*_1.fastq.gz", "*_2.fastq.gz"),
    ):
        r1 = sorted(sample_dir.glob(r1_pat))
        r2 = sorted(sample_dir.glob(r2_pat))
        if len(r1) == 1 and len(r2) == 1:
            return r1[0], r2[0]

    raise FileNotFoundError(
        f"Expected one R1 and one R2 under {sample_dir}; "
        f"found: {', '.join(p.name for p in sorted(sample_dir.iterdir()) if p.is_file())}"
    )


def link_or_copy(src: Path, dest: Path, *, copy: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        if dest.resolve() == src.resolve():
            return
        dest.unlink()
    if copy:
        shutil.copy2(src, dest)
        return
    try:
        os.link(src, dest)
    except OSError:
        dest.symlink_to(src.resolve())


def stage_sample_reads(
    sample: BeringiaBirdSample,
    *,
    downloads_root: Path,
    data_root: Path,
    copy: bool,
    dry_run: bool,
) -> tuple[Path, Path]:
    """Stage reads at data/beringia_birds/<directory_name>/<directory_name>_R{1,2}.fastq.gz."""
    source_dir = downloads_root / sample.directory_name
    dest_dir = data_root / sample.directory_name
    r1_src, r2_src = find_read_pair_in_dir(source_dir)
    r1_dest = dest_dir / f"{sample.directory_name}_R1.fastq.gz"
    r2_dest = dest_dir / f"{sample.directory_name}_R2.fastq.gz"

    if dry_run:
        return r1_dest, r2_dest

    link_or_copy(r1_src, r1_dest, copy=copy)
    link_or_copy(r2_src, r2_dest, copy=copy)
    return r1_dest, r2_dest


def update_checklist_staged_paths(checklist: Path) -> int:
    """Rewrite old-directory column from beringia_birds/... to data/beringia_birds/..."""
    text = checklist.read_text(encoding="utf-8")
    updated = 0
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith("|") and OLD_DIR_PREFIX in line:
            new_line = line.replace(f"`{OLD_DIR_PREFIX}", f"`{STAGED_DIR_PREFIX}")
            if new_line != line:
                updated += 1
            lines.append(new_line)
        else:
            lines.append(line)
    if updated:
        checklist.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return updated


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Stage Downloads/beringia_birds reads and batch-assemble mitogenomes."
    )
    p.add_argument(
        "--field",
        action="append",
        dest="fields",
        metavar="FIELD",
        help="Process only this field number (repeatable).",
    )
    p.add_argument("--fresh-log", action="store_true", help="Overwrite log file.")
    p.add_argument("-v", "--verbose", action="store_true", help="Stream tool output.")
    p.add_argument("--dry-run", action="store_true", help="Print actions without staging or assembling.")
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip assembly when output already has a non-empty NOVOPlasty FASTA.",
    )
    p.add_argument(
        "--stage-only",
        action="store_true",
        help="Only stage reads into data/beringia_birds/; do not assemble.",
    )
    p.add_argument(
        "--skip-stage",
        action="store_true",
        help="Assemble from existing staged reads (skip Downloads -> data/).",
    )
    p.add_argument(
        "--copy",
        action="store_true",
        help="Copy FASTQs instead of hardlink/symlink into data/beringia_birds/.",
    )
    p.add_argument(
        "--downloads-root",
        type=Path,
        default=DOWNLOADS_ROOT,
        help=f"Source tree (default: {DOWNLOADS_ROOT}).",
    )
    p.add_argument(
        "--data-root",
        type=Path,
        default=DATA_ROOT,
        help=f"Staged read tree in repo (default: {DATA_ROOT}).",
    )
    p.add_argument(
        "--update-checklist",
        action="store_true",
        help="After staging, set Table 1 old directory to data/beringia_birds/<sample>/.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    downloads_root = args.downloads_root.expanduser().resolve()
    data_root = args.data_root.resolve()

    if not args.skip_stage and not downloads_root.is_dir():
        raise SystemExit(f"Downloads directory missing: {downloads_root}")
    if not REFSEQ.is_file():
        raise SystemExit(f"Reference FASTA missing: {REFSEQ}")

    samples = parse_beringia_birds_samples(CHECKLIST)
    if args.fields:
        wanted = set(args.fields)
        samples = [s for s in samples if s.field_num in wanted]
        missing = wanted - {s.field_num for s in samples}
        if missing:
            raise SystemExit(f"Unknown or non-ingest field(s): {', '.join(sorted(missing))}")

    if not samples:
        raise SystemExit("No beringia_birds samples found in Table 1 of the checklist.")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    log_mode = "w" if args.fresh_log else "a"
    total = len(samples)
    staged_any = False

    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(f"Downloads root: {downloads_root}\n")
        log.write(f"Data root: {data_root}\n")
        log.write(f"Samples: {total}\n")
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")

        for idx, sample in enumerate(samples, start=1):
            out_dir = OUT_ROOT / sample.directory_name
            log.write(
                f"\n[{idx}/{total}] {sample.field_num} -> {sample.directory_name}\n"
            )
            log.flush()

            if not args.skip_stage:
                if args.dry_run:
                    src_dir = downloads_root / sample.directory_name
                    r1_src, r2_src = find_read_pair_in_dir(src_dir)
                    r1_dest = data_root / sample.directory_name / (
                        f"{sample.directory_name}_R1.fastq.gz"
                    )
                    r2_dest = r1_dest.with_name(
                        f"{sample.directory_name}_R2.fastq.gz"
                    )
                    log.write(f"  Stage: {r1_src} -> {r1_dest}\n")
                    log.write(f"  Stage: {r2_src} -> {r2_dest}\n")
                else:
                    r1_dest, r2_dest = stage_sample_reads(
                        sample,
                        downloads_root=downloads_root,
                        data_root=data_root,
                        copy=args.copy,
                        dry_run=False,
                    )
                    staged_any = True
                    log.write(f"  Staged R1: {r1_dest}\n  Staged R2: {r2_dest}\n")
            else:
                dest_dir = data_root / sample.directory_name
                r1_dest, r2_dest = find_read_pair_in_dir(dest_dir)

            if args.stage_only or args.dry_run:
                if args.stage_only and not args.dry_run:
                    log.write("  Stage-only; assembly skipped.\n")
                continue

            cmd = assemble_cmd(
                sample.field_num,
                out_dir,
                r1_dest,
                r2_dest,
                verbose=args.verbose,
            )
            log.write(f"  CMD: {' '.join(cmd)}\n")
            log.flush()

            if args.dry_run:
                continue
            if args.skip_existing and novoplasty_done(out_dir, sample.field_num):
                log.write("  Skipped: existing NOVOPlasty FASTA found.\n")
                log.flush()
                continue

            r = subprocess.run(cmd, env=assemble_child_environ())
            log.write(f"  Exit code: {r.returncode}\n")
            log.flush()
            if r.returncode != 0:
                raise SystemExit(
                    f"Assembly failed for {sample.field_num} (exit {r.returncode}); see {LOG}"
                )

        if args.update_checklist and staged_any and not args.dry_run:
            n = update_checklist_staged_paths(CHECKLIST)
            log.write(f"\nChecklist: updated {n} old-directory path(s) to {STAGED_DIR_PREFIX}...\n")


if __name__ == "__main__":
    main()
