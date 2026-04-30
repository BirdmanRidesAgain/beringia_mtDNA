#!/usr/bin/env python3
"""
Run MitoZ annotation (and optionally MITOS2) for a mitochondrial assembly FASTA.

This wrapper creates a sample-specific output directory and runs:
1) MitoZ annotate (default)
2) MITOS2 (optional; enabled with --run-mitos2)

Both commands are configurable via template strings so the script remains usable
across different local MITOS2/MitoZ installations.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


def default_threads() -> int:
    cpu_count = os.cpu_count() or 1
    return max(1, cpu_count // 2)


def run_external(cmd: list[str], *, verbose: bool, label: str) -> None:
    if verbose:
        print(f"Running {label}:")
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)
        return

    completed = subprocess.run(cmd, capture_output=True, text=True)
    if completed.returncode != 0:
        if completed.stderr:
            print(completed.stderr.rstrip(), file=sys.stderr)
        if completed.stdout:
            print(completed.stdout.rstrip(), file=sys.stderr)
        raise subprocess.CalledProcessError(
            completed.returncode, cmd, completed.stdout, completed.stderr
        )


def build_command(template: str, values: dict[str, str]) -> list[str]:
    rendered = template.format(**values)
    return shlex.split(rendered)


def require_executable(executable: str, *, install_hint: str) -> None:
    if shutil.which(executable) is None:
        raise FileNotFoundError(
            f"Required executable '{executable}' not found on PATH.\n"
            f"Install hint: {install_hint}"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run MitoZ annotation and optionally MITOS2 for a mitogenome FASTA."
    )
    parser.add_argument(
        "-i",
        "--input-fasta",
        required=True,
        help="Path to input mitogenome FASTA (circularized or noncircularized).",
    )
    parser.add_argument(
        "-s",
        "--sample-name",
        help="Sample name for output folder/prefix. Defaults to FASTA stem.",
    )
    parser.add_argument(
        "-o",
        "--output-root",
        default="annotation_results",
        help="Root output directory (default: annotation_results).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=default_threads(),
        help="Thread count for tools that support it.",
    )
    parser.add_argument(
        "--genetic-code",
        default="2",
        help="Mitochondrial genetic code table (default: 2 for vertebrate mtDNA).",
    )
    parser.add_argument(
        "--mitoz-clade",
        default="Chordata",
        help="MitoZ clade value (default: Chordata).",
    )
    parser.add_argument(
        "--run-mitos2",
        action="store_true",
        help="Also run MITOS2 (default: off).",
    )
    parser.add_argument(
        "--mitos2-template",
        default="mitos -i {fasta} -o {mitos_out}",
        help=(
            "MITOS2 command template. Available fields: {fasta}, {sample}, "
            "{sample_dir}, {mitos_out}, {mitoz_out}, {threads}, {genetic_code}, {clade}"
        ),
    )
    parser.add_argument(
        "--mitoz-template",
        default=(
            "mitoz annotate --genetic_code {genetic_code} --clade {clade} "
            "--outprefix {sample} --thread_number {threads} "
            "--fastafile {fasta} --workdir {mitoz_out}"
        ),
        help=(
            "MitoZ command template. Available fields: {fasta}, {sample}, "
            "{sample_dir}, {mitos_out}, {mitoz_out}, {threads}, {genetic_code}, {clade}"
        ),
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print full commands and stream output.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    fasta = Path(args.input_fasta).expanduser().resolve()
    if not fasta.is_file():
        raise FileNotFoundError(f"Input FASTA not found: {fasta}")

    sample_name = args.sample_name or fasta.stem
    sample_dir = Path(args.output_root) / sample_name
    mitos_out = sample_dir / "mitos2"
    mitoz_out = sample_dir / "mitoz"

    mitoz_out.mkdir(parents=True, exist_ok=True)
    if args.run_mitos2:
        mitos_out.mkdir(parents=True, exist_ok=True)

    values = {
        "fasta": str(fasta),
        "sample": sample_name,
        "sample_dir": str(sample_dir.resolve()),
        "mitos_out": str(mitos_out.resolve()),
        "mitoz_out": str(mitoz_out.resolve()),
        "threads": str(args.threads),
        "genetic_code": str(args.genetic_code),
        "clade": str(args.mitoz_clade),
    }

    require_executable(
        "mitoz",
        install_hint=(
            "Activate your annotation conda environment, then install MitoZ. "
            "Example: conda activate mitogenome-tools && "
            "python -m pip install mitoz"
        ),
    )
    if args.run_mitos2:
        require_executable(
            "mitos",
            install_hint=(
                "Install MITOS2 in a dedicated conda env due version constraints. "
                "Example: conda create --strict-channel-priority "
                "-c conda-forge -c bioconda -n mitos2 'mitos>=2' && "
                "conda activate mitos2"
            ),
        )

    mitoz_cmd = build_command(args.mitoz_template, values)

    run_external(mitoz_cmd, verbose=args.verbose, label="MitoZ")
    if args.run_mitos2:
        mitos_cmd = build_command(args.mitos2_template, values)
        run_external(mitos_cmd, verbose=args.verbose, label="MITOS2")

    print(f"Annotation completed for {sample_name}")
    print(f"MitoZ output: {mitoz_out.resolve()}")
    if args.run_mitos2:
        print(f"MITOS2 output: {mitos_out.resolve()}")


if __name__ == "__main__":
    main()
