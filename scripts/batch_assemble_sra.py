#!/usr/bin/env python3
"""
Batch-assemble mitogenomes for paired FASTQs under ~/Downloads/sra_downloads/fastq.

Uses conda env ``mitogenome-tools``: runs ``assemble_mitogenome.py`` with that env's
``python`` and prepends its ``bin`` to ``PATH`` (same as ``conda activate``).

Examples::

    # Full 128-sample run, new log file
    python scripts/batch_assemble_sra.py --fresh-log

    # Only samples 50–80 (1-based indices), append to log
    python scripts/batch_assemble_sra.py --first 50 --last 80
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
FASTQ_ROOT = Path.home() / "Downloads/sra_downloads/fastq"
REFSEQ = REPO / "test/refseq/Taeniopygia_guttata_NC_007897.1_mitogenome.fa"
OUT_ROOT = REPO / "mitogenomes_output"
ASM = REPO / "scripts/assemble_mitogenome.py"
LOG = OUT_ROOT / "batch_assemble_sra.log"

CONDA_ENV = "mitogenome-tools"


def conda_executable() -> Path:
    """Prefer known install paths so ``CONDA_EXE`` from another shell does not point at Miniforge."""
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


def mitogenome_python() -> Path:
    """Interpreter for ``mitogenome-tools`` (must exist)."""
    p = conda_install_base() / "envs" / CONDA_ENV / "bin" / "python"
    if not p.is_file():
        raise FileNotFoundError(
            f"Expected conda env {CONDA_ENV!r} at {p}. "
            f"Create it or set CONDA_ROOT if conda uses a different install location."
        )
    return p


def assemble_child_environ() -> dict[str, str]:
    """PATH prefix so ``fastp``, ``NOVOPlasty.pl``, etc. match ``conda activate mitogenome-tools``."""
    env = os.environ.copy()
    bin_dir = mitogenome_python().parent
    env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")
    return env


def assemble_cmd(sample: str, fq1: Path, fq2: Path) -> list[str]:
    sample_out = OUT_ROOT / sample
    return [
        str(mitogenome_python()),
        str(ASM),
        "-1",
        str(fq1),
        "-2",
        str(fq2),
        "-s",
        sample,
        "-o",
        str(sample_out),
        "-r",
        str(REFSEQ),
        "-S",
    ]


def discover_runs() -> list[tuple[str, Path, Path]]:
    pattern = re.compile(r"(.+)_([12])\.fastq(?:\.gz)?$")
    pairs: dict[tuple[str, str], dict[str, Path]] = {}
    for p in FASTQ_ROOT.rglob("*"):
        if not p.is_file():
            continue
        m = pattern.match(p.name)
        if not m:
            continue
        key = (str(p.parent), m.group(1))
        pairs.setdefault(key, {})[m.group(2)] = p

    runs: list[tuple[str, Path, Path]] = []
    for (_, sample), d in sorted(pairs.items()):
        if "1" in d and "2" in d:
            runs.append((sample, d["1"], d["2"]))
    return runs


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch mitogenome assembly for SRA FASTQ pairs.")
    p.add_argument("--first", type=int, default=1, help="1-based first sample index (default: 1).")
    p.add_argument(
        "--last",
        type=int,
        default=None,
        help="1-based last sample index inclusive (default: last pair found).",
    )
    p.add_argument(
        "--fresh-log",
        action="store_true",
        help="Overwrite batch_assemble_sra.log instead of appending.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    runs = discover_runs()
    if not runs:
        raise SystemExit(f"No paired FASTQs found under {FASTQ_ROOT}")

    total = len(runs)
    first = max(1, args.first)
    last = total if args.last is None else min(args.last, total)
    if first > last or first > total:
        raise SystemExit(f"Invalid range --first {args.first} --last {args.last} (have {total} pairs).")

    subset = runs[first - 1 : last]

    log_mode = "w" if args.fresh_log else "a"
    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(f"Total paired samples: {total}\n")
        log.write(
            f"\n--- START {datetime.now(timezone.utc).isoformat()}Z: "
            f"samples [{first}/{total}]–[{last}/{total}] ---\n"
        )
        for offset, (sample, fq1, fq2) in enumerate(subset):
            idx = first + offset
            log.write(f"\n[{idx}/{total}] {sample}\n")
            log.flush()
            r = subprocess.run(assemble_cmd(sample, fq1, fq2), env=assemble_child_environ())
            log.write(f"Exit code: {r.returncode}\n")
            log.flush()


if __name__ == "__main__":
    main()
