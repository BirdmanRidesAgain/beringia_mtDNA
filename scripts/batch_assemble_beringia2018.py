#!/usr/bin/env python3
"""
Batch-assemble the 18 pending mitogenomes from ~/Downloads/new_beringia_birds.

Reads sample metadata from the pending-assembly block in mitogenomes_assembled_checklist.md
(Table 1). Uses conda env ``mitogenome-tools`` and ``assemble_mitogenome.py`` with
``--skip-bbsplit``, matching batch_assemble_sra.py.

Example::

    python scripts/batch_assemble_beringia2018.py --fresh-log
    python scripts/batch_assemble_beringia2018.py --field JJW905 --verbose
"""
from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CHECKLIST = REPO / "mitogenomes_assembled_checklist.md"
FASTQ_ROOT = Path.home() / "Downloads/new_beringia_birds"
REFSEQ = REPO / "test/refseq/Taeniopygia_guttata_NC_007897.1_mitogenome.fa"
OUT_ROOT = REPO / "assembled_mitogenomes"
ASM = REPO / "scripts/assemble_mitogenome.py"
LOG = OUT_ROOT / "batch_assemble_beringia2018.log"

CONDA_ENV = "mitogenome-tools"
PENDING_MARKER = "<!-- Pending assembly:"


def conda_executable() -> Path:
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
    raise FileNotFoundError("Could not find conda.")


def conda_install_base() -> Path:
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


def mitogenome_python() -> Path:
    p = conda_install_base() / "envs" / CONDA_ENV / "bin" / "python"
    if not p.is_file():
        raise FileNotFoundError(
            f"Expected conda env {CONDA_ENV!r} at {p}. "
            "Create it from resources/environment.yml."
        )
    return p


def assemble_child_environ() -> dict[str, str]:
    env = os.environ.copy()
    bin_dir = mitogenome_python().parent
    env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")
    return env


def parse_pending_samples(checklist: Path) -> list[tuple[str, str]]:
    """Return (field_num, output_dir_name) for each pending row."""
    text = checklist.read_text(encoding="utf-8")
    start = text.find(PENDING_MARKER)
    if start < 0:
        raise ValueError(f"No pending block found in {checklist}")
    rest = text[start:]
    end = rest.find("\n## Table 2.")
    if end < 0:
        raise ValueError("Could not find end of pending block (Table 2 header).")
    block = rest[:end]
    samples: list[tuple[str, str]] = []
    for line in block.splitlines():
        if not line.startswith("|") or line.startswith("| -") or "Field num" in line:
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 5:
            continue
        field = parts[2]
        new_dir = parts[4].strip("`").strip("/")
        if field and new_dir:
            samples.append((field, new_dir))
    return samples


def find_read_pair(field: str, root: Path) -> tuple[Path, Path]:
    """Resolve paired FASTQs (flat layout or one subdirectory per field #)."""
    search_roots = [root]
    sub = root / field
    if sub.is_dir():
        search_roots.insert(0, sub)

    r1: list[Path] = []
    r2: list[Path] = []
    for base in search_roots:
        r1 = sorted(base.glob(f"{field}_*_1.fq.gz"))
        r2 = sorted(base.glob(f"{field}_*_2.fq.gz"))
        if r1 or r2:
            break
        r1 = sorted(base.glob("*_1.fq.gz"))
        r2 = sorted(base.glob("*_2.fq.gz"))
        if r1 or r2:
            break

    if len(r1) != 1 or len(r2) != 1:
        raise FileNotFoundError(
            f"Expected one R1 and one R2 for {field} under {root}; "
            f"found R1={len(r1)}, R2={len(r2)}"
        )
    return r1[0], r2[0]


def assemble_cmd(field: str, out_dir: Path, fq1: Path, fq2: Path, *, verbose: bool) -> list[str]:
    cmd = [
        str(mitogenome_python()),
        str(ASM),
        "-1",
        str(fq1),
        "-2",
        str(fq2),
        "-s",
        field,
        "-o",
        str(out_dir),
        "-r",
        str(REFSEQ),
        "-S",
    ]
    if verbose:
        cmd.append("-v")
    return cmd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch assembly for 2018 Beringia resequencing set.")
    p.add_argument(
        "--field",
        action="append",
        dest="fields",
        metavar="FIELD",
        help="Assemble only this field number (repeatable). Default: all pending samples.",
    )
    p.add_argument("--fresh-log", action="store_true", help="Overwrite log file.")
    p.add_argument("-v", "--verbose", action="store_true", help="Stream tool output.")
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running assembly.",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip samples whose output dir already has a non-empty NOVOPlasty FASTA.",
    )
    return p.parse_args()


def novoplasty_done(out_dir: Path, field: str) -> bool:
    nov_dir = out_dir / f"novoplasty_{field}"
    if not nov_dir.is_dir():
        return False
    for pat in (
        "Circularized_assembly_*.fasta",
        "C_*.fasta",
        "Option_1_*.fasta",
        "Contigs_*.fasta",
    ):
        for f in nov_dir.glob(pat):
            if f.is_file() and f.stat().st_size > 0:
                return True
    return False


def main() -> None:
    args = parse_args()
    if not REFSEQ.is_file():
        raise SystemExit(f"Reference FASTA missing: {REFSEQ}")
    if not FASTQ_ROOT.is_dir():
        raise SystemExit(f"FASTQ directory missing: {FASTQ_ROOT}")

    samples = parse_pending_samples(CHECKLIST)
    if args.fields:
        wanted = set(args.fields)
        samples = [(f, d) for f, d in samples if f in wanted]
        missing = wanted - {f for f, _ in samples}
        if missing:
            raise SystemExit(f"Unknown or non-pending field(s): {', '.join(sorted(missing))}")

    if not samples:
        raise SystemExit("No samples to assemble.")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    log_mode = "w" if args.fresh_log else "a"
    total = len(samples)

    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(f"FASTQ root: {FASTQ_ROOT}\n")
        log.write(f"Samples: {total}\n")
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")
        for idx, (field, out_name) in enumerate(samples, start=1):
            out_dir = OUT_ROOT / out_name
            fq1, fq2 = find_read_pair(field, FASTQ_ROOT)
            cmd = assemble_cmd(field, out_dir, fq1, fq2, verbose=args.verbose)
            log.write(f"\n[{idx}/{total}] {field} -> {out_name}\n")
            log.write(f"  R1: {fq1}\n  R2: {fq2}\n")
            log.write(f"  CMD: {' '.join(cmd)}\n")
            log.flush()
            if args.dry_run:
                continue
            if args.skip_existing and novoplasty_done(out_dir, field):
                log.write("  Skipped: existing NOVOPlasty FASTA found.\n")
                log.flush()
                continue
            r = subprocess.run(cmd, env=assemble_child_environ())
            log.write(f"Exit code: {r.returncode}\n")
            log.flush()
            if r.returncode != 0:
                raise SystemExit(f"Assembly failed for {field} (exit {r.returncode}); see {LOG}")


if __name__ == "__main__":
    main()
