#!/usr/bin/env python3
"""
Assemble mitogenomes pipeline.

- fastp: trim reads to ./<output>/trim_<sample_name>/
- bbsplit (optional): mtDNA-enriched reads, or skip and feed trimmed reads to NOVOPlasty
- NOVOPlasty: assemble using repo template config; seed is always --refseq
"""

from __future__ import annotations

import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

PACKAGE_ROOT = Path(__file__).resolve().parent
# Bundled defaults and bbsplit index cache (similar to `src/main/resources` / app “resources” dirs).
RESOURCES_DIR = PACKAGE_ROOT / "resources"
NOVOPLASTY_TEMPLATE = RESOURCES_DIR / "NOVOPlasty_config.txt"


def default_threads() -> int:
    """Use half of available CPUs by default (minimum 1)."""
    cpu_count = os.cpu_count() or 1
    return max(1, cpu_count // 2)


def run_external(
    cmd: list[str],
    *,
    verbose: bool,
    cwd: Path | None = None,
    label: str | None = None,
) -> None:
    """Run a subprocess; by default capture output and only print it on failure."""
    cwd_arg = str(cwd) if cwd is not None else None
    if verbose:
        if label:
            print(f"Running {label}:")
        print(" ".join(cmd))
        subprocess.run(cmd, check=True, cwd=cwd_arg)
        return

    completed = subprocess.run(cmd, cwd=cwd_arg, capture_output=True, text=True)
    if completed.returncode != 0:
        err = (completed.stderr or "").rstrip()
        out = (completed.stdout or "").rstrip()
        if err:
            print(err, file=sys.stderr)
        if out:
            print(out, file=sys.stderr)
        raise subprocess.CalledProcessError(
            completed.returncode, cmd, completed.stdout, completed.stderr
        )


def run_fastp(
    fq1: str,
    fq2: str,
    sample_name: str,
    output: str,
    threads: int,
    verbose: bool,
) -> tuple[Path, Path]:
    """Run fastp and return paths to trimmed paired reads."""
    trim_dir = Path(output) / f"trim_{sample_name}"
    trim_dir.mkdir(parents=True, exist_ok=True)

    trimmed_r1 = trim_dir / f"{sample_name}_R1_trimmed.fastq.gz"
    trimmed_r2 = trim_dir / f"{sample_name}_R2_trimmed.fastq.gz"
    report_html = trim_dir / f"{sample_name}_fastp.html"
    report_json = trim_dir / f"{sample_name}_fastp.json"

    cmd = [
        "fastp",
        "-i",
        fq1,
        "-I",
        fq2,
        "-o",
        str(trimmed_r1),
        "-O",
        str(trimmed_r2),
        "-h",
        str(report_html),
        "-j",
        str(report_json),
        "-w",
        str(threads),
    ]

    run_external(cmd, verbose=verbose, label="fastp")

    return trimmed_r1, trimmed_r2


def stage_reference_for_bbsplit(refseq: str) -> Path:
    """
    Copy the user FASTA into resources/bbsplit_ref/<id>/ and use that path as bbsplit `ref=`.

    BBTools writes index files next to the reference; staging keeps those out of the user's directory.
    """
    src = Path(refseq).expanduser().resolve()
    if not src.is_file():
        raise FileNotFoundError(f"Reference FASTA not found: {src}")

    st = src.stat()
    cache_id = hashlib.sha256(f"{src}\0{st.st_mtime_ns}\0{st.st_size}".encode()).hexdigest()[:16]
    cache_dir = RESOURCES_DIR / "bbsplit_ref" / cache_id
    cache_dir.mkdir(parents=True, exist_ok=True)

    staged = cache_dir / src.name
    if staged.exists():
        dst_st = staged.stat()
        if dst_st.st_size == st.st_size and dst_st.st_mtime_ns == st.st_mtime_ns:
            return staged
    shutil.copy2(src, staged)
    return staged


def run_bbsplit(
    trimmed_r1: Path,
    trimmed_r2: Path,
    sample_name: str,
    refseq: Path,
    output: str,
    memory_gb: int,
    threads: int,
    verbose: bool,
) -> tuple[Path, Path]:
    """Run bbsplit on trimmed reads and return mtDNA-enriched read paths."""
    bbsplit_dir = Path(output) / f"bbsplit_{sample_name}"
    bbsplit_dir.mkdir(parents=True, exist_ok=True)

    mt_r1 = bbsplit_dir / f"{sample_name}_R1_trimmed_mtDNA.fastq.gz"
    mt_r2 = bbsplit_dir / f"{sample_name}_R2_trimmed_mtDNA.fastq.gz"

    cmd = [
        "bbsplit.sh",
        f"-Xmx{memory_gb}g",
        f"in1={trimmed_r1.resolve()}",
        f"in2={trimmed_r2.resolve()}",
        f"ref={refseq.resolve()}",
        f"outm1={mt_r1.resolve()}",
        f"outm2={mt_r2.resolve()}",
        f"threads={threads}",
    ]

    run_external(cmd, verbose=verbose, label="bbsplit")

    return mt_r1, mt_r2


def _replace_first_config_value(config_text: str, key: str, value: str) -> str:
    pattern = rf"^({re.escape(key)}\s*=\s*).*$"
    return re.sub(pattern, rf"\g<1>{value}", config_text, count=1, flags=re.MULTILINE)


def run_novoplasty(
    sample_name: str,
    refseq: str,
    reads_r1: Path,
    reads_r2: Path,
    output: str,
    verbose: bool,
) -> Path:
    """Run NOVOPlasty using a sample-specific config derived from a template."""
    novoplasty_dir = Path(output) / f"novoplasty_{sample_name}"
    novoplasty_dir.mkdir(parents=True, exist_ok=True)

    if not NOVOPLASTY_TEMPLATE.is_file():
        raise FileNotFoundError(
            f"NOVOPlasty template missing: {NOVOPLASTY_TEMPLATE} (expected under {RESOURCES_DIR})"
        )
    template_text = NOVOPLASTY_TEMPLATE.read_text(encoding="utf-8")
    updated_text = _replace_first_config_value(template_text, "Project name", sample_name)
    updated_text = _replace_first_config_value(updated_text, "Seed Input", str(Path(refseq).resolve()))
    updated_text = _replace_first_config_value(updated_text, "Forward reads", str(reads_r1.resolve()))
    updated_text = _replace_first_config_value(updated_text, "Reverse reads", str(reads_r2.resolve()))

    run_config = (novoplasty_dir / f"{sample_name}_NOVOPlasty_config.txt").resolve()
    run_config.write_text(updated_text, encoding="utf-8")

    cmd = ["NOVOPlasty.pl", "-c", str(run_config)]
    run_external(cmd, verbose=verbose, cwd=novoplasty_dir, label="NOVOPlasty")

    return run_config


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Mitogenome assembly scaffold")
    parser.add_argument("-1", "--fq1", required=True, help="Input R1 FASTQ(.gz)")
    parser.add_argument("-2", "--fq2", required=True, help="Input R2 FASTQ(.gz)")
    parser.add_argument("-s", "--sample_name", required=True, help="Sample identifier")
    parser.add_argument("-o", "--output", required=True, help="Output directory root")
    parser.add_argument("-r", "--refseq", required=True, help="Reference sequence fasta (also NOVOPlasty seed)")
    parser.add_argument(
        "-S",
        "--skip-bbsplit",
        action="store_true",
        help="Send trimmed reads directly to NOVOPlasty (skip bbsplit)",
    )
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        default=8,
        help="Memory (GB) for bbsplit; ignored when --skip-bbsplit is set; default: 8",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=default_threads(),
        help="Threads for fastp and bbsplit; default: half available CPUs",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Print commands and stream tool output; default is quiet except errors",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    verbose = args.verbose

    trimmed_r1, trimmed_r2 = run_fastp(
        fq1=args.fq1,
        fq2=args.fq2,
        sample_name=args.sample_name,
        output=args.output,
        threads=args.threads,
        verbose=verbose,
    )

    if verbose:
        print(f"Trimmed R1: {trimmed_r1}")
        print(f"Trimmed R2: {trimmed_r2}")

    if args.skip_bbsplit:
        novoplasty_r1, novoplasty_r2 = trimmed_r1, trimmed_r2
        if verbose:
            print("Skipping bbsplit; NOVOPlasty will use trimmed reads")
    else:
        bbsplit_ref = stage_reference_for_bbsplit(args.refseq)
        if verbose:
            print(f"bbsplit ref (staged): {bbsplit_ref.resolve()}")
        novoplasty_r1, novoplasty_r2 = run_bbsplit(
            trimmed_r1=trimmed_r1,
            trimmed_r2=trimmed_r2,
            sample_name=args.sample_name,
            refseq=bbsplit_ref,
            output=args.output,
            memory_gb=args.memory,
            threads=args.threads,
            verbose=verbose,
        )
        if verbose:
            print(f"mtDNA R1: {novoplasty_r1}")
            print(f"mtDNA R2: {novoplasty_r2}")

    run_config = run_novoplasty(
        sample_name=args.sample_name,
        refseq=args.refseq,
        reads_r1=novoplasty_r1,
        reads_r2=novoplasty_r2,
        output=args.output,
        verbose=verbose,
    )
    if verbose:
        print(f"NOVOPlasty config: {run_config}")


if __name__ == "__main__":
    main()
