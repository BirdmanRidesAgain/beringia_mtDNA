"""Runtime helpers for MitoZ batch jobs (SSD checks, trim watchdog)."""

from __future__ import annotations

import os
import platform
import subprocess
import time
from pathlib import Path

DEFAULT_SSD_MOUNT = Path("/Volumes/KC_SSD")
MIN_TRIM_BYTES = 1024
# Batch timings: >2M pairs often exceed 30 min; 4M+ can take many hours on SPAdes.
DEFAULT_MAX_READ_PAIRS = 1_500_000
DEFAULT_DOWNSAMPLE_SEED = 42


def default_assembler() -> str:
    """MEGAHIT 1.2.9 segfaults on macOS during assembly; SPAdes works under Rosetta."""
    if platform.system() == "Darwin":
        return "spades"
    return "megahit"


def default_mitoz_memory_gb() -> str:
    return "16"


def default_max_read_pairs() -> int:
    return DEFAULT_MAX_READ_PAIRS


def trim_output_paths(trim_root: Path, sample_name: str) -> tuple[Path, Path]:
    trim_dir = trim_root / f"trim_{sample_name}"
    return (
        trim_dir / f"{sample_name}_R1_trimmed.fastq.gz",
        trim_dir / f"{sample_name}_R2_trimmed.fastq.gz",
    )


def mitoz_input_paths(trim_root: Path, sample_name: str) -> tuple[Path, Path]:
    """Downsampled FASTQs fed to MitoZ (only created when trimming exceeds the cap)."""
    trim_dir = trim_root / f"trim_{sample_name}"
    return (
        trim_dir / f"{sample_name}_R1_mitoz.fastq.gz",
        trim_dir / f"{sample_name}_R2_mitoz.fastq.gz",
    )


def count_read_pairs(fastq: Path, *, seqkit: Path) -> int:
    """Return the number of reads in a (gzipped) FASTQ via ``seqkit stats``."""
    completed = subprocess.run(
        [str(seqkit), "stats", "-T", str(fastq)],
        capture_output=True,
        text=True,
        check=True,
    )
    line = completed.stdout.strip().splitlines()[-1]
    return int(line.split("\t")[3])


def downsample_read_pairs(
    *,
    r1: Path,
    r2: Path,
    out_r1: Path,
    out_r2: Path,
    max_pairs: int,
    seed: int,
    threads: int,
    seqkit: Path,
    verbose: bool,
) -> None:
    """Randomly subsample paired reads to *max_pairs* using the same seed on R1/R2."""
    out_r1.parent.mkdir(parents=True, exist_ok=True)
    for src, dest in ((r1, out_r1), (r2, out_r2)):
        cmd = [
            str(seqkit),
            "sample",
            "--two-pass",
            "-s",
            str(seed),
            "-j",
            str(threads),
            "-n",
            str(max_pairs),
            str(src),
            "-o",
            str(dest),
        ]
        if verbose:
            print("Running:", " ".join(cmd))
            subprocess.run(cmd, check=True)
            continue
        completed = subprocess.run(cmd, capture_output=True, text=True)
        if completed.returncode != 0:
            err = (completed.stderr or completed.stdout or "").rstrip()
            raise subprocess.CalledProcessError(completed.returncode, cmd, output=err)


def prepare_mitoz_reads(
    *,
    trimmed_r1: Path,
    trimmed_r2: Path,
    trim_root: Path,
    sample_name: str,
    max_read_pairs: int,
    seed: int,
    threads: int,
    seqkit: Path,
    force: bool,
    verbose: bool,
) -> tuple[Path, Path, int, int | None]:
    """
    Choose MitoZ input reads, downsampling when trimmed depth exceeds *max_read_pairs*.

    Returns ``(r1, r2, raw_pairs, subsampled_pairs_or_none)``.
    """
    raw_pairs = count_read_pairs(trimmed_r1, seqkit=seqkit)
    if max_read_pairs <= 0 or raw_pairs <= max_read_pairs:
        return trimmed_r1, trimmed_r2, raw_pairs, None

    mitoz_r1, mitoz_r2 = mitoz_input_paths(trim_root, sample_name)
    if not force and trim_outputs_ready(mitoz_r1, mitoz_r2):
        return mitoz_r1, mitoz_r2, raw_pairs, max_read_pairs

    if force:
        for path in (mitoz_r1, mitoz_r2):
            if path.exists():
                path.unlink()

    downsample_read_pairs(
        r1=trimmed_r1,
        r2=trimmed_r2,
        out_r1=mitoz_r1,
        out_r2=mitoz_r2,
        max_pairs=max_read_pairs,
        seed=seed,
        threads=threads,
        seqkit=seqkit,
        verbose=verbose,
    )
    return mitoz_r1, mitoz_r2, raw_pairs, max_read_pairs


def ssd_is_mounted(mount_point: Path) -> bool:
    """Return True when *mount_point* exists and is a filesystem mount."""
    try:
        return mount_point.is_dir() and os.path.ismount(mount_point)
    except OSError:
        return False


def trim_outputs_ready(r1: Path, r2: Path) -> bool:
    """Both trimmed FASTQs exist and are non-empty."""
    try:
        return (
            r1.is_file()
            and r2.is_file()
            and r1.stat().st_size >= MIN_TRIM_BYTES
            and r2.stat().st_size >= MIN_TRIM_BYTES
        )
    except OSError:
        return False


def trim_byte_total(r1: Path, r2: Path) -> int:
    total = 0
    for path in (r1, r2):
        try:
            if path.is_file():
                total += path.stat().st_size
        except OSError:
            pass
    return total


class TrimHangError(RuntimeError):
    """Raised when trimming stalls or the SSD is unavailable."""


def run_subprocess_with_trim_watchdog(
    cmd: list[str],
    *,
    trim_root: Path,
    sample_name: str,
    ssd_mount: Path,
    trim_hang_minutes: float,
    poll_seconds: float,
    skip_trim_watch: bool,
    log_write,
) -> int:
    """
    Run *cmd* and abort the batch when trimming stalls or the SSD unmounts.

    Trim monitoring stops once both trimmed FASTQs exist with non-zero size.
    """
    trim_r1, trim_r2 = trim_output_paths(trim_root, sample_name)
    proc = subprocess.Popen(cmd)
    last_total = -1
    last_progress = time.monotonic()
    trim_phase = not skip_trim_watch

    try:
        while proc.poll() is None:
            time.sleep(poll_seconds)

            if not ssd_is_mounted(ssd_mount):
                proc.kill()
                proc.wait(timeout=30)
                msg = f"SSD not mounted at {ssd_mount}; aborting batch"
                log_write(msg)
                raise TrimHangError(msg)

            if not trim_phase:
                continue

            if trim_outputs_ready(trim_r1, trim_r2):
                trim_phase = False
                log_write("Trim outputs ready; disabling trim hang watchdog")
                continue

            total = trim_byte_total(trim_r1, trim_r2)
            if total > last_total:
                last_total = total
                last_progress = time.monotonic()
                continue

            stalled_min = (time.monotonic() - last_progress) / 60.0
            if stalled_min >= trim_hang_minutes:
                proc.kill()
                proc.wait(timeout=30)
                msg = (
                    f"Trimming stalled for {stalled_min:.1f} min "
                    f"(no growth in {trim_r1.parent}); aborting batch"
                )
                log_write(msg)
                raise TrimHangError(msg)

        return proc.returncode if proc.returncode is not None else 1
    except TrimHangError:
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=30)
        raise
