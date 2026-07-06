#!/usr/bin/env python3
"""
Trim paired FASTQs (fastp), then assemble and annotate with MitoZ.

Trimmed reads are written to the SSD under ``trimmed/<sample_name>/``.
MitoZ outputs go to ``mitogenomes_output/<sample_name>/mitoz/``.

Requires conda envs ``mitogenome-tools`` (fastp) and ``mitoz-x64`` (mitoz).

On macOS, MEGAHIT 1.2.9 segfaults during assembly; the default assembler is
SPAdes (``--assembler spades``). Use ``--fastp-no-reports`` for large batches
to avoid fastp hanging while writing HTML/JSON reports.

After trimming, samples with more than ``--max-read-pairs`` (default 1.5M)
are randomly subsampled with seqkit before MitoZ, targeting ~30-minute SPAdes
runs on this hardware. Use ``--max-read-pairs 0`` to disable.

Example::

    python scripts/trim_assemble_annotate_mitoz.py \\
        -1 /path/to/sample_R1.fastq.gz \\
        -2 /path/to/sample_R2.fastq.gz \\
        -s Anas_acuta_UAM_22113_REW601
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from assemble_mitogenome import default_threads, run_fastp  # noqa: E402
from batch_annotate_assemblies import annotate_child_environ, mitoz_python  # noqa: E402
from mitogenome_paths import REPO  # noqa: E402
from mitoz_runtime import (  # noqa: E402
    DEFAULT_DOWNSAMPLE_SEED,
    default_assembler,
    default_max_read_pairs,
    default_mitoz_memory_gb,
    prepare_mitoz_reads,
    trim_output_paths,
    trim_outputs_ready,
)

DEFAULT_TRIM_ROOT = Path("/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/trimmed")
DEFAULT_ASM_ROOT = REPO / "mitogenomes_output"
CONDA_TRIM_ENV = "mitogenome-tools"


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
    raise FileNotFoundError(
        "Could not find conda. Set CONDA_EXE or install Anaconda/Miniconda."
    )


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


def trim_python() -> Path:
    p = conda_install_base() / "envs" / CONDA_TRIM_ENV / "bin" / "python"
    if not p.is_file():
        raise FileNotFoundError(
            f"Expected conda env {CONDA_TRIM_ENV!r} at {p}. "
            "Create it with fastp before running this script."
        )
    return p


def trim_child_environ() -> dict[str, str]:
    env = os.environ.copy()
    bin_dir = trim_python().parent
    env["PATH"] = str(bin_dir) + os.pathsep + env.get("PATH", "")
    return env


def mitoz_seqkit() -> Path:
    seqkit = mitoz_python().parent / "seqkit"
    if not seqkit.is_file():
        raise FileNotFoundError(f"seqkit not found in mitoz env: {seqkit}")
    return seqkit


def mitoz_done(workdir: Path, outprefix: str) -> bool:
    if not workdir.is_dir():
        return False
    return any(workdir.rglob(f"{outprefix}*.gbf")) or any(workdir.rglob("*.result"))


def run_mitoz_all(
    *,
    fq1: Path,
    fq2: Path,
    sample_name: str,
    workdir: Path,
    threads: int,
    clade: str,
    requiring_taxa: str,
    genetic_code: str,
    species_name: str | None,
    assembler: str,
    memory_gb: str,
    verbose: bool,
) -> None:
    workdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mitoz",
        "all",
        "--skip_filter",
        "--assembler",
        assembler,
        "--memory",
        memory_gb,
        "--fq1",
        str(fq1.resolve()),
        "--fq2",
        str(fq2.resolve()),
        "--workdir",
        str(workdir.resolve()),
        "--outprefix",
        sample_name,
        "--thread_number",
        str(threads),
        "--clade",
        clade,
        "--requiring_taxa",
        requiring_taxa,
        "--genetic_code",
        genetic_code,
    ]
    if species_name:
        cmd.extend(["--species_name", species_name])

    env = annotate_child_environ()
    if verbose:
        print("Running MitoZ:")
        print(" ".join(cmd))
        subprocess.run(cmd, check=True, cwd=str(workdir), env=env)
        return

    completed = subprocess.run(cmd, cwd=str(workdir), env=env, capture_output=True, text=True)
    if completed.returncode != 0:
        if completed.stderr:
            print(completed.stderr.rstrip(), file=sys.stderr)
        if completed.stdout:
            print(completed.stdout.rstrip(), file=sys.stderr)
        raise subprocess.CalledProcessError(
            completed.returncode, cmd, completed.stdout, completed.stderr
        )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Trim reads (fastp on SSD), then assemble+annotate with MitoZ."
    )
    p.add_argument("-1", "--r1", required=True, help="Input R1 FASTQ(.gz)")
    p.add_argument("-2", "--r2", required=True, help="Input R2 FASTQ(.gz)")
    p.add_argument("-s", "--sample-name", required=True, help="Sample / output directory name")
    p.add_argument(
        "--trim-root",
        type=Path,
        default=DEFAULT_TRIM_ROOT,
        help=f"SSD root for trimmed reads (default: {DEFAULT_TRIM_ROOT})",
    )
    p.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_ASM_ROOT,
        help=f"Repo output root (default: {DEFAULT_ASM_ROOT})",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=default_threads(),
        help="Threads for fastp and MitoZ (default: half available CPUs)",
    )
    p.add_argument("--clade", default="Chordata", help="MitoZ clade (default: Chordata)")
    p.add_argument(
        "--requiring-taxa",
        default="Chordata",
        help="MitoZ --requiring_taxa filter (default: Chordata)",
    )
    p.add_argument(
        "--genetic-code",
        default="2",
        help="Mitochondrial genetic code for MitoZ (default: 2)",
    )
    p.add_argument(
        "--assembler",
        default=default_assembler(),
        choices=("spades", "megahit", "mitoassemble"),
        help=f"MitoZ assembler (default: {default_assembler()} on this OS)",
    )
    p.add_argument(
        "--memory",
        default=default_mitoz_memory_gb(),
        help=f"MitoZ memory limit in GB for SPAdes/MEGAHIT (default: {default_mitoz_memory_gb()})",
    )
    p.add_argument(
        "--fastp-no-reports",
        action="store_true",
        help="Skip fastp HTML/JSON reports (recommended for large batch runs)",
    )
    p.add_argument(
        "--max-read-pairs",
        type=int,
        default=default_max_read_pairs(),
        metavar="N",
        help=(
            "Downsample trimmed reads above N pairs before MitoZ "
            f"(default: {default_max_read_pairs():,}). Use 0 to disable."
        ),
    )
    p.add_argument(
        "--downsample-seed",
        type=int,
        default=DEFAULT_DOWNSAMPLE_SEED,
        help=f"Random seed for seqkit subsampling (default: {DEFAULT_DOWNSAMPLE_SEED})",
    )
    p.add_argument(
        "--species-name",
        help="Species name for MitoZ GenBank output (optional)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Re-run even if trimmed reads and/or MitoZ outputs already exist",
    )
    p.add_argument(
        "--skip-mitoz",
        action="store_true",
        help="Only run fastp trimming (for testing)",
    )
    p.add_argument("-v", "--verbose", action="store_true", help="Print commands and stream output")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    r1 = Path(args.r1).expanduser().resolve()
    r2 = Path(args.r2).expanduser().resolve()
    if not r1.is_file():
        raise FileNotFoundError(f"R1 not found: {r1}")
    if not r2.is_file():
        raise FileNotFoundError(f"R2 not found: {r2}")

    sample_name = args.sample_name.strip()
    if not sample_name:
        raise ValueError("sample-name must not be empty")

    trim_python()
    mitoz_python()

    trim_root = args.trim_root.expanduser().resolve()
    trim_root.mkdir(parents=True, exist_ok=True)

    trim_dir = trim_root / f"trim_{sample_name}"
    trimmed_r1, trimmed_r2 = trim_output_paths(trim_root, sample_name)

    if trim_outputs_ready(trimmed_r1, trimmed_r2) and not args.force:
        print(f"Using existing trimmed reads in {trim_dir}")
    else:
        if trim_dir.exists() and args.force:
            shutil.rmtree(trim_dir)
        env = trim_child_environ()
        prev_path = os.environ.get("PATH")
        os.environ["PATH"] = env["PATH"]
        try:
            trimmed_r1, trimmed_r2 = run_fastp(
                fq1=str(r1),
                fq2=str(r2),
                sample_name=sample_name,
                output=str(trim_root),
                threads=args.threads,
                verbose=args.verbose,
                write_reports=not args.fastp_no_reports,
            )
        finally:
            if prev_path is None:
                os.environ.pop("PATH", None)
            else:
                os.environ["PATH"] = prev_path
        print(f"Trimmed R1: {trimmed_r1}")
        print(f"Trimmed R2: {trimmed_r2}")

    if args.skip_mitoz:
        return

    mitoz_workdir = args.output_root.expanduser().resolve() / sample_name / "mitoz"
    if mitoz_done(mitoz_workdir, sample_name) and not args.force:
        print(f"MitoZ output already present in {mitoz_workdir} (use --force to re-run)")
        return

    if mitoz_workdir.exists() and args.force:
        shutil.rmtree(mitoz_workdir)

    mitoz_r1, mitoz_r2, raw_pairs, subsampled_pairs = prepare_mitoz_reads(
        trimmed_r1=Path(trimmed_r1),
        trimmed_r2=Path(trimmed_r2),
        trim_root=trim_root,
        sample_name=sample_name,
        max_read_pairs=args.max_read_pairs,
        seed=args.downsample_seed,
        threads=args.threads,
        seqkit=mitoz_seqkit(),
        force=args.force,
        verbose=args.verbose,
    )
    if subsampled_pairs is not None:
        print(
            f"Downsampled {raw_pairs:,} trimmed pairs to {subsampled_pairs:,} "
            f"for MitoZ (cap {args.max_read_pairs:,})"
        )
    else:
        print(f"MitoZ input: {raw_pairs:,} trimmed read pairs (no downsampling)")

    run_mitoz_all(
        fq1=mitoz_r1,
        fq2=mitoz_r2,
        sample_name=sample_name,
        workdir=mitoz_workdir,
        threads=args.threads,
        clade=args.clade,
        requiring_taxa=args.requiring_taxa,
        genetic_code=args.genetic_code,
        species_name=args.species_name,
        assembler=args.assembler,
        memory_gb=args.memory,
        verbose=args.verbose,
    )
    print(f"MitoZ output: {mitoz_workdir.resolve()}")


if __name__ == "__main__":
    main()
