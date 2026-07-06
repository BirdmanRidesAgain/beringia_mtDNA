#!/usr/bin/env python3
"""
Batch trim + MitoZ assemble/annotate for paired FASTQs on the SSD.

By default, only individuals with ``assembly_success`` of ``False`` in
``mitogenome_assembly_results.csv`` are processed. Use ``--all-samples`` to run
every directory on the SSD instead (433 samples).

Discovers one sample per subdirectory under
``/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/raw_data/``, then runs
``scripts/trim_assemble_annotate_mitoz.py`` on each sample sequentially.

Designed for long unattended runs (e.g. overnight). Use ``--skip-completed`` to
resume after interruption.

On macOS the pipeline defaults to MitoZ + SPAdes (MEGAHIT segfaults on Apple
Silicon). Trimming uses ``--fastp-no-reports`` by default so fastp does not
hang generating HTML/JSON on large FASTQs.

The batch aborts immediately if the SSD unmounts or if trimming shows no output
growth for ``--trim-hang-minutes`` (likely a detached drive or stuck fastp).

Trimmed samples with more than ``--max-read-pairs`` (default 1.5M) are randomly
subsampled with seqkit before MitoZ to keep SPAdes runs near ~30 minutes.

Examples::

    # Dry run: list samples and planned commands
    python scripts/batch_mitoz_ssd.py --dry-run

    # Reassembly run (default: partial_success + failed from CSV)
    nohup python scripts/batch_mitoz_ssd.py --skip-completed --fresh-log \\
        > mitogenomes_output/_batch_logs/batch_mitoz_ssd.nohup.out 2>&1 &

    # All SSD samples (not just CSV-filtered reassembly targets)
    python scripts/batch_mitoz_ssd.py --all-samples --dry-run

    # Restart all samples from scratch after a fix
    nohup python scripts/batch_mitoz_ssd.py --force --fresh-log \\
        > mitogenomes_output/_batch_logs/batch_mitoz_ssd.nohup.out 2>&1 &

    # All SSD samples without MitoZ pipeline output (novoplasty-only + never attempted)
    # Completed sample dirs are archived to SSD automatically after each run.
    nohup caffeinate -i python scripts/batch_mitoz_ssd.py --pending-mitoz --skip-completed --fresh-log \\
        > mitogenomes_output/_batch_logs/batch_mitoz_ssd_pending.nohup.out 2>&1 &
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from assemble_mitogenome import default_threads  # noqa: E402
from mitogenome_paths import ASM_ROOT, REPO  # noqa: E402
from mitoz_runtime import (  # noqa: E402
    DEFAULT_SSD_MOUNT,
    TrimHangError,
    default_assembler,
    default_max_read_pairs,
    default_mitoz_memory_gb,
    run_subprocess_with_trim_watchdog,
    ssd_is_mounted,
    trim_output_paths,
    trim_outputs_ready,
)
from trim_assemble_annotate_mitoz import (  # noqa: E402
    DEFAULT_TRIM_ROOT,
    mitoz_done,
)
from archive_mitoz_to_ssd import archive_one_sample, local_free_gb  # noqa: E402

RAW_DATA = Path("/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/raw_data")
ARCHIVE_ROOT = Path("/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/assembled_mitogenomes")
CHECKLIST = REPO / "mitogenome_assembly_checklist.md"
RESULTS_CSV = REPO / "mitogenome_assembly_results.csv"
DEFAULT_ASSEMBLY_SUCCESS = (False,)
PIPELINE = SCRIPTS / "trim_assemble_annotate_mitoz.py"
BATCH_LOGS = ASM_ROOT / "_batch_logs"
LOG = BATCH_LOGS / "batch_mitoz_ssd.log"
LOG_PENDING = BATCH_LOGS / "batch_mitoz_ssd_pending.log"
SKIP_TSV = BATCH_LOGS / "batch_mitoz_ssd_skipped.tsv"
FILTERED_TSV = BATCH_LOGS / "batch_mitoz_ssd_filtered_by_csv.tsv"
FAILED_TSV = BATCH_LOGS / "batch_mitoz_ssd_failed.tsv"
FAILED_PENDING_TSV = BATCH_LOGS / "batch_mitoz_ssd_pending_failed.tsv"

PAIR_SUFFIXES = (
    ("_R1.fastq.gz", "_R2.fastq.gz"),
    ("_R1.fq.gz", "_R2.fq.gz"),
    ("_1.fastq.gz", "_2.fastq.gz"),
    ("_1.fq.gz", "_2.fq.gz"),
)


@dataclass(frozen=True)
class SampleRun:
    sample_name: str
    r1: Path
    r2: Path
    species_name: str | None = None
    assembly_success: bool | None = None


def parse_bool_cell(value: str | None) -> bool | None:
    if value is None:
        return None
    text = value.strip().lower()
    if text in ("true", "1", "yes"):
        return True
    if text in ("false", "0", "no"):
        return False
    return None


def load_results_index(csv_path: Path) -> dict[str, dict[str, str]]:
    if not csv_path.is_file():
        raise FileNotFoundError(f"Results CSV not found: {csv_path}")
    index: dict[str, dict[str, str]] = {}
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not reader.fieldnames or "directory_name" not in reader.fieldnames:
            raise ValueError(f"{csv_path} must include a directory_name column")
        if "assembly_success" not in reader.fieldnames:
            raise ValueError(f"{csv_path} must include an assembly_success column")
        for row in reader:
            directory = row["directory_name"].strip()
            if directory:
                index[directory] = row
    return index


def parse_assembly_success_values(spec: str) -> tuple[bool, ...]:
    values: list[bool] = []
    for part in spec.split(","):
        parsed = parse_bool_cell(part)
        if parsed is None:
            raise ValueError(f"Invalid assembly_success value: {part!r} (use true/false)")
        values.append(parsed)
    if not values:
        raise ValueError("At least one assembly_success value is required")
    return tuple(values)


def filter_runs_by_assembly_success(
    runs: list[SampleRun],
    results_index: dict[str, dict[str, str]],
    allowed: set[bool],
) -> tuple[list[SampleRun], list[tuple[str, str]]]:
    filtered: list[SampleRun] = []
    excluded: list[tuple[str, str]] = []
    for run in runs:
        row = results_index.get(run.sample_name)
        if row is None:
            excluded.append((run.sample_name, "not listed in mitogenome_assembly_results.csv"))
            continue
        success = parse_bool_cell(row.get("assembly_success"))
        if success is None or success not in allowed:
            excluded.append(
                (run.sample_name, f"assembly_success={row.get('assembly_success') or '(empty)'}")
            )
            continue
        species = run.species_name or row.get("species", "").strip() or None
        filtered.append(
            SampleRun(
                sample_name=run.sample_name,
                r1=run.r1,
                r2=run.r2,
                species_name=species,
                assembly_success=success,
            )
        )
    return filtered, excluded


def enrich_runs_from_csv(
    runs: list[SampleRun],
    results_index: dict[str, dict[str, str]],
) -> list[SampleRun]:
    """Fill species and assembly_success from the results CSV when missing."""
    enriched: list[SampleRun] = []
    for run in runs:
        row = results_index.get(run.sample_name)
        if row is None:
            enriched.append(run)
            continue
        species = run.species_name or row.get("species", "").strip() or None
        success = run.assembly_success
        if success is None:
            success = parse_bool_cell(row.get("assembly_success"))
        enriched.append(
            SampleRun(
                sample_name=run.sample_name,
                r1=run.r1,
                r2=run.r2,
                species_name=species,
                assembly_success=success,
            )
        )
    return enriched


def mitoz_complete_local_or_archived(sample_name: str, output_root: Path) -> bool:
    """True when MitoZ finished locally or under the SSD archive."""
    local_mitoz = output_root / sample_name / "mitoz"
    if mitoz_done(local_mitoz, sample_name):
        return True
    archived_mitoz = ARCHIVE_ROOT / sample_name / "mitoz"
    return mitoz_done(archived_mitoz, sample_name)


def load_failed_samples(path: Path) -> dict[str, tuple[int, str]]:
    """Load prior pipeline failures from a failed TSV (sample_name -> exit_code, message)."""
    if not path.is_file():
        return {}
    prior: dict[str, tuple[int, str]] = {}
    with path.open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = (row.get("sample_name") or "").strip()
            if not name:
                continue
            try:
                code = int(row.get("exit_code") or 1)
            except ValueError:
                code = 1
            msg = (row.get("message") or "pipeline failed").strip()
            prior[name] = (code, msg)
    return prior


def write_failed_tsv(path: Path, failed: list[tuple[str, int, str]]) -> None:
    """Write failed samples, preserving any prior entries not superseded."""
    merged = load_failed_samples(path)
    for name, code, msg in failed:
        merged[name] = (code, msg)
    if not merged:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        f.write("sample_name\texit_code\tmessage\n")
        for name, (code, msg) in sorted(merged.items()):
            f.write(f"{name}\t{code}\t{msg}\n")


def filter_pending_mitoz(
    runs: list[SampleRun],
    output_root: Path,
    *,
    prior_failures: dict[str, tuple[int, str]] | None = None,
) -> tuple[list[SampleRun], list[tuple[str, str]]]:
    """Keep SSD samples without complete MitoZ output and without a prior pipeline failure."""
    prior_failures = prior_failures or {}
    pending: list[SampleRun] = []
    skipped: list[tuple[str, str]] = []
    for run in runs:
        if mitoz_complete_local_or_archived(run.sample_name, output_root):
            skipped.append((run.sample_name, "mitoz complete (local or SSD archive)"))
            continue
        if run.sample_name in prior_failures:
            skipped.append((run.sample_name, "prior pipeline failure (not retrying)"))
            continue
        pending.append(run)
    return pending, skipped


def parse_species_by_directory(checklist: Path) -> dict[str, str]:
    """Map checklist New directory name -> species (Table 1 only)."""
    mapping: dict[str, str] = {}
    in_table = False
    for line in checklist.read_text(encoding="utf-8").splitlines():
        if line.startswith("| Species | Catalog"):
            in_table = True
            continue
        if in_table and line.startswith("## "):
            break
        if not in_table or not line.startswith("|") or line.startswith("| - |"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 5:
            continue
        species = re.sub(r"^\*|\*$", "", parts[0]).strip()
        directory = parts[4].strip("`").strip("/")
        if directory and species:
            mapping[directory] = species
    return mapping


def list_read_pairs(sample_dir: Path) -> list[tuple[Path, Path]]:
    pairs: list[tuple[Path, Path]] = []
    seen: set[tuple[Path, Path]] = set()
    for s1, s2 in PAIR_SUFFIXES:
        for r1 in sorted(sample_dir.glob(f"*{s1}")):
            if r1.name.startswith("._"):
                continue
            r2 = sample_dir / (r1.name[: -len(s1)] + s2)
            if not r2.is_file():
                continue
            key = (r1.resolve(), r2.resolve())
            if key in seen:
                continue
            seen.add(key)
            pairs.append((r1, r2))
    return pairs


def find_loose_pair(sample_dir: Path) -> tuple[Path, Path] | None:
    """Fallback when R1/R2 basenames differ but only one of each exists."""
    for r1_pat, r2_pat in (("*R1*.fastq.gz", "*R2*.fastq.gz"), ("*R1*.fq.gz", "*R2*.fq.gz")):
        r1 = sorted(
            p for p in sample_dir.glob(r1_pat) if p.is_file() and not p.name.startswith("._")
        )
        r2 = sorted(
            p for p in sample_dir.glob(r2_pat) if p.is_file() and not p.name.startswith("._")
        )
        if len(r1) == 1 and len(r2) == 1:
            return r1[0], r2[0]
    return None


def find_read_pair(sample_dir: Path) -> tuple[Path, Path]:
    """Pick one R1/R2 pair under a raw_data sample directory."""
    if not sample_dir.is_dir():
        raise FileNotFoundError(f"Sample directory missing: {sample_dir}")

    name = sample_dir.name
    for ext in (".fastq.gz", ".fq.gz"):
        canonical_r1 = sample_dir / f"{name}_R1{ext}"
        canonical_r2 = sample_dir / f"{name}_R2{ext}"
        if canonical_r1.is_file() and canonical_r2.is_file():
            return canonical_r1, canonical_r2

    pairs = list_read_pairs(sample_dir)
    if not pairs:
        loose = find_loose_pair(sample_dir)
        if loose:
            return loose
        files = ", ".join(p.name for p in sorted(sample_dir.iterdir()) if p.is_file())
        raise FileNotFoundError(f"No R1/R2 pair under {sample_dir}; files: {files or '(none)'}")
    if len(pairs) == 1:
        return pairs[0]

    prefixed = [p for p in pairs if p[0].name.startswith(name)]
    if len(prefixed) == 1:
        return prefixed[0]
    if prefixed:
        return max(prefixed, key=lambda p: p[0].stat().st_size + p[1].stat().st_size)

    return max(pairs, key=lambda p: p[0].stat().st_size + p[1].stat().st_size)


def discover_samples(
    raw_root: Path,
    *,
    species_by_dir: dict[str, str],
) -> tuple[list[SampleRun], list[tuple[str, str]]]:
    runs: list[SampleRun] = []
    skipped: list[tuple[str, str]] = []

    sample_dirs = sorted(
        p for p in raw_root.iterdir() if p.is_dir() and not p.name.startswith(".")
    )
    for sample_dir in sample_dirs:
        name = sample_dir.name
        try:
            r1, r2 = find_read_pair(sample_dir)
        except FileNotFoundError as exc:
            skipped.append((name, str(exc)))
            continue
        runs.append(
            SampleRun(
                sample_name=name,
                r1=r1,
                r2=r2,
                species_name=species_by_dir.get(name),
            )
        )
    return runs, skipped


def pipeline_cmd(
    run: SampleRun,
    *,
    trim_root: Path,
    output_root: Path,
    threads: int,
    assembler: str,
    memory_gb: str,
    max_read_pairs: int,
    fastp_no_reports: bool,
    verbose: bool,
    force: bool,
) -> list[str]:
    cmd = [
        sys.executable,
        str(PIPELINE),
        "-1",
        str(run.r1),
        "-2",
        str(run.r2),
        "-s",
        run.sample_name,
        "--trim-root",
        str(trim_root),
        "--output-root",
        str(output_root),
        "--threads",
        str(threads),
        "--assembler",
        assembler,
        "--memory",
        memory_gb,
    ]
    cmd.extend(["--max-read-pairs", str(max_read_pairs)])
    if fastp_no_reports:
        cmd.append("--fastp-no-reports")
    if run.species_name:
        cmd.extend(["--species-name", run.species_name])
    if force:
        cmd.append("--force")
    if verbose:
        cmd.append("-v")
    return cmd


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Batch trim + MitoZ for SSD raw_data samples. "
            "By default only assembly_success=False individuals from "
            "mitogenome_assembly_results.csv are processed."
        )
    )
    p.add_argument(
        "--results-csv",
        type=Path,
        default=RESULTS_CSV,
        help=f"Assembly results CSV for level filtering (default: {RESULTS_CSV})",
    )
    p.add_argument(
        "--assembly-success",
        default=",".join("true" if v else "false" for v in DEFAULT_ASSEMBLY_SUCCESS),
        help=(
            "Comma-separated assembly_success booleans to include from the CSV "
            f"(default: {','.join('true' if v else 'false' for v in DEFAULT_ASSEMBLY_SUCCESS)}). "
            "Ignored with --all-samples."
        ),
    )
    p.add_argument(
        "--all-samples",
        action="store_true",
        help="Process every SSD sample with read pairs, ignoring CSV assembly_success filter",
    )
    p.add_argument(
        "--pending-mitoz",
        action="store_true",
        help=(
            "Process every SSD sample with read pairs that lacks sample/mitoz/ output "
            "(novoplasty-only successes, failed reassembly attempts, and SSD dirs not yet run)"
        ),
    )
    p.add_argument(
        "--raw-root",
        type=Path,
        default=RAW_DATA,
        help=f"SSD raw FASTQ root (default: {RAW_DATA})",
    )
    p.add_argument(
        "--trim-root",
        type=Path,
        default=DEFAULT_TRIM_ROOT,
        help=f"SSD trimmed-read root (default: {DEFAULT_TRIM_ROOT})",
    )
    p.add_argument(
        "--output-root",
        type=Path,
        default=ASM_ROOT,
        help=f"Assembly output root (default: {ASM_ROOT})",
    )
    p.add_argument(
        "--sample",
        action="append",
        dest="samples",
        metavar="NAME",
        help="Only this sample directory name (repeatable)",
    )
    p.add_argument("--first", type=int, default=1, help="1-based first sample index (default: 1)")
    p.add_argument(
        "--last",
        type=int,
        default=None,
        help="1-based last sample index inclusive (default: last sample found)",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=default_threads(),
        help="Threads per sample for fastp and MitoZ (default: half available CPUs)",
    )
    p.add_argument(
        "--skip-completed",
        action="store_true",
        help=(
            "Skip samples whose MitoZ output already exists, and samples with a prior "
            "pipeline failure recorded in the failed TSV (do not retry data failures)"
        ),
    )
    p.add_argument(
        "--assembler",
        default=default_assembler(),
        choices=("spades", "megahit", "mitoassemble"),
        help=f"MitoZ assembler passed to trim_assemble_annotate_mitoz.py (default: {default_assembler()})",
    )
    p.add_argument(
        "--memory",
        default=default_mitoz_memory_gb(),
        help=f"MitoZ memory limit in GB (default: {default_mitoz_memory_gb()})",
    )
    p.add_argument(
        "--max-read-pairs",
        type=int,
        default=default_max_read_pairs(),
        metavar="N",
        help=(
            "Downsample trimmed reads above N pairs before MitoZ "
            f"(default: {default_max_read_pairs():,}; ~30 min SPAdes target). Use 0 to disable."
        ),
    )
    p.add_argument(
        "--no-downsample",
        action="store_true",
        help="Disable read downsampling before MitoZ (same as --max-read-pairs 0)",
    )
    p.add_argument(
        "--fastp-reports",
        action="store_true",
        help="Write fastp HTML/JSON reports (off by default for batch runs)",
    )
    p.add_argument(
        "--ssd-mount",
        type=Path,
        default=DEFAULT_SSD_MOUNT,
        help=f"SSD mount point to verify before/during each sample (default: {DEFAULT_SSD_MOUNT})",
    )
    p.add_argument(
        "--trim-hang-minutes",
        type=float,
        default=90.0,
        help="Abort batch if trim outputs show no growth for this many minutes (default: 90)",
    )
    p.add_argument(
        "--poll-seconds",
        type=float,
        default=60.0,
        help="Watchdog polling interval while a sample runs (default: 60)",
    )
    p.add_argument(
        "--force",
        action="store_true",
        help="Pass --force to trim_assemble_annotate_mitoz.py (re-run trim and MitoZ)",
    )
    p.add_argument(
        "--no-archive-after-sample",
        action="store_true",
        help="Do not rsync each completed sample dir to SSD and remove it locally",
    )
    p.add_argument(
        "--archive-root",
        type=Path,
        default=ARCHIVE_ROOT,
        help=f"SSD archive root for --archive-after-sample (default: {ARCHIVE_ROOT})",
    )
    p.add_argument("--fresh-log", action="store_true", help="Overwrite batch log file")
    p.add_argument("-v", "--verbose", action="store_true", help="Stream per-sample tool output")
    p.add_argument("--dry-run", action="store_true", help="Print planned work only")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.pending_mitoz and args.all_samples:
        raise SystemExit("Use only one of --pending-mitoz or --all-samples")

    raw_root = args.raw_root.expanduser().resolve()
    trim_root = args.trim_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()

    ssd_mount = args.ssd_mount.expanduser().resolve()
    fastp_no_reports = not args.fastp_reports
    max_read_pairs = 0 if args.no_downsample else args.max_read_pairs
    log_path = LOG_PENDING if args.pending_mitoz else LOG
    failed_tsv = FAILED_PENDING_TSV if args.pending_mitoz else FAILED_TSV
    prior_failures = load_failed_samples(failed_tsv)

    if not raw_root.is_dir():
        raise SystemExit(f"Raw data root missing: {raw_root}")
    if not ssd_is_mounted(ssd_mount):
        raise SystemExit(f"SSD not mounted at {ssd_mount}")
    if not PIPELINE.is_file():
        raise SystemExit(f"Pipeline script missing: {PIPELINE}")

    species_by_dir = parse_species_by_directory(CHECKLIST) if CHECKLIST.is_file() else {}
    runs, discovery_skipped = discover_samples(raw_root, species_by_dir=species_by_dir)

    results_index = load_results_index(args.results_csv.expanduser().resolve())
    runs = enrich_runs_from_csv(runs, results_index)

    csv_filtered: list[tuple[str, str]] = []
    pending_mitoz_skipped: list[tuple[str, str]] = []
    allowed_success: tuple[bool, ...] = ()
    discovered_total = len(runs)
    if args.pending_mitoz:
        if not args.force:
            runs, pending_mitoz_skipped = filter_pending_mitoz(
                runs,
                output_root,
                prior_failures=prior_failures if args.skip_completed else None,
            )
    elif args.all_samples:
        pass
    else:
        allowed_success = parse_assembly_success_values(args.assembly_success)
        discovered_total = len(runs)
        runs, csv_filtered = filter_runs_by_assembly_success(
            runs,
            results_index,
            set(allowed_success),
        )

    if args.samples:
        wanted = set(args.samples)
        runs = [r for r in runs if r.sample_name in wanted]
        missing = wanted - {r.sample_name for r in runs}
        if missing:
            raise SystemExit(f"Unknown or unreadable sample(s): {', '.join(sorted(missing))}")

    if not runs:
        raise SystemExit(f"No runnable samples found under {raw_root}")

    total = len(runs)
    first = max(1, args.first)
    last = total if args.last is None else min(args.last, total)
    if first > last or first > total:
        raise SystemExit(f"Invalid range --first {args.first} --last {args.last} (have {total} samples)")

    subset = runs[first - 1 : last]

    skipped_completed = 0
    skipped_prior_failures = 0
    if args.skip_completed and not args.force:
        pending: list[SampleRun] = []
        for run in subset:
            if mitoz_complete_local_or_archived(run.sample_name, output_root):
                skipped_completed += 1
            elif run.sample_name in prior_failures:
                skipped_prior_failures += 1
            else:
                pending.append(run)
        subset = pending

    BATCH_LOGS.mkdir(parents=True, exist_ok=True)
    log_mode = "w" if args.fresh_log else "a"
    failed: list[tuple[str, int, str]] = []
    succeeded = 0
    archived = 0
    aborted = False
    archive_after_sample = not args.no_archive_after_sample
    archive_root = args.archive_root.expanduser().resolve()

    with log_path.open(log_mode, encoding="utf-8") as log:
        log.write(f"Raw root: {raw_root}\n")
        log.write(f"Trim root: {trim_root}\n")
        log.write(f"Output root: {output_root}\n")
        log.write(f"Discovered: {discovered_total} sample(s) with read pairs on SSD\n")
        if args.pending_mitoz:
            log.write(
                f"Pending MitoZ filter: {len(runs)} to process, "
                f"{len(pending_mitoz_skipped)} skipped "
                f"(complete or prior failure)\n"
            )
        if prior_failures and args.skip_completed:
            log.write(f"Prior failures on record: {len(prior_failures)}\n")
        elif not args.all_samples:
            log.write(f"Results CSV: {args.results_csv.expanduser().resolve()}\n")
            log.write(
                f"assembly_success filter: "
                f"{', '.join('true' if v else 'false' for v in allowed_success)} "
                f"({len(runs)} matched, {len(csv_filtered)} excluded by CSV filter)\n"
            )
        if discovery_skipped:
            log.write(f"Discovery skipped: {len(discovery_skipped)} dir(s) without pairs\n")
        log.write(
            f"Running index [{first}/{total}]–[{last}/{total}] "
            f"({len(subset)} to process, {skipped_completed} already done, "
            f"{skipped_prior_failures} prior failures not retried)\n"
        )
        log.write(f"Threads per sample: {args.threads}\n")
        log.write(f"Assembler: {args.assembler}  Memory (GB): {args.memory}\n")
        if max_read_pairs > 0:
            log.write(f"Max read pairs (MitoZ input): {max_read_pairs:,}\n")
        else:
            log.write("Max read pairs (MitoZ input): disabled\n")
        log.write(f"SSD mount: {ssd_mount}\n")
        log.write(
            f"Trim hang watchdog: {args.trim_hang_minutes:g} min "
            f"(poll every {args.poll_seconds:g}s)\n"
        )
        log.write(f"fastp reports: {not fastp_no_reports}\n")
        log.write(
            f"archive after sample: {archive_after_sample} -> {archive_root}\n"
            if archive_after_sample
            else "archive after sample: disabled\n"
        )
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")
        log.flush()

        for offset, run in enumerate(subset):
            if not ssd_is_mounted(ssd_mount):
                msg = f"SSD not mounted at {ssd_mount} before sample {run.sample_name}"
                log.write(f"\nABORT: {msg}\n")
                log.flush()
                raise SystemExit(msg)

            idx = first + offset
            mitoz_dir = output_root / run.sample_name / "mitoz"
            trimmed_r1, trimmed_r2 = trim_output_paths(trim_root, run.sample_name)
            skip_trim_watch = (
                not args.force and trim_outputs_ready(trimmed_r1, trimmed_r2)
            )
            cmd = pipeline_cmd(
                run,
                trim_root=trim_root,
                output_root=output_root,
                threads=args.threads,
                assembler=args.assembler,
                memory_gb=args.memory,
                max_read_pairs=max_read_pairs,
                fastp_no_reports=fastp_no_reports,
                verbose=args.verbose,
                force=args.force,
            )
            block = [
                f"\n[{idx}/{total}] {run.sample_name}",
                f"  R1: {run.r1}",
                f"  R2: {run.r2}",
            ]
            if run.assembly_success is not None:
                block.append(f"  assembly_success: {run.assembly_success}")
            if run.species_name:
                block.append(f"  species: {run.species_name}")
            if skip_trim_watch:
                block.append("  trim watchdog: skipped (reusing existing trimmed reads)")
            block.append(f"  CMD: {' '.join(cmd)}")
            log.write("\n".join(block) + "\n")
            log.flush()

            if args.dry_run:
                continue

            t0 = time.monotonic()

            def log_write(msg: str) -> None:
                log.write(f"  WATCHDOG: {msg}\n")
                log.flush()

            try:
                returncode = run_subprocess_with_trim_watchdog(
                    cmd,
                    trim_root=trim_root,
                    sample_name=run.sample_name,
                    ssd_mount=ssd_mount,
                    trim_hang_minutes=args.trim_hang_minutes,
                    poll_seconds=args.poll_seconds,
                    skip_trim_watch=skip_trim_watch,
                    log_write=log_write,
                )
            except TrimHangError as exc:
                log.write(f"  ABORT: {exc}\n")
                log.flush()
                aborted = True
                break

            elapsed = time.monotonic() - t0
            log.write(f"  Exit code: {returncode}\n")
            log.write(f"  Elapsed: {elapsed / 60:.1f} min\n")
            log.write(f"  Output: {mitoz_dir}\n")
            log.flush()

            if returncode == 0:
                succeeded += 1
            else:
                failed.append((run.sample_name, returncode, "pipeline failed"))

            if archive_after_sample and mitoz_dir.is_dir() and mitoz_done(mitoz_dir, run.sample_name):
                if not archive_root.parent.is_dir():
                    log.write("  ARCHIVE: skipped (SSD archive path not available)\n")
                    log.flush()
                else:
                    try:
                        t_arch = time.monotonic()
                        outcome = archive_one_sample(
                            run.sample_name,
                            archive_root=archive_root,
                            output_root=output_root,
                            log=log,
                        )
                        arch_elapsed = time.monotonic() - t_arch
                        if outcome is None:
                            log.write("  ARCHIVE: skipped (no archivable local sample dir)\n")
                        elif outcome == "already_archived":
                            log.write(
                                f"  ARCHIVE: already on SSD; removed local "
                                f"({arch_elapsed / 60:.1f} min); "
                                f"free {local_free_gb():.1f} GB\n"
                            )
                        else:
                            archived += 1
                            log.write(
                                f"  ARCHIVE: rsync+verify OK ({arch_elapsed / 60:.1f} min); "
                                f"free {local_free_gb():.1f} GB\n"
                            )
                    except Exception as exc:  # noqa: BLE001
                        log.write(f"  ARCHIVE ERROR: {exc}\n")
                    log.flush()

        if aborted:
            log.write(
                f"\n--- ABORT {datetime.now(timezone.utc).isoformat()}Z ---\n"
                f"Batch stopped by trim/SSD watchdog after {succeeded} success(es), "
                f"{archived} archived\n"
            )
        else:
            log.write(
                f"\n--- END {datetime.now(timezone.utc).isoformat()}Z ---\n"
                f"Summary: {succeeded} succeeded, {len(failed)} failed, "
                f"{archived} archived to SSD, "
                f"{skipped_completed} skipped (already done), "
                f"{skipped_prior_failures} skipped (prior failure), "
                f"{len(discovery_skipped)} skipped (no read pair)\n"
            )
        for name, code, msg in sorted(failed):
            log.write(f"  failed [{code}] {name}: {msg}\n")

    if discovery_skipped:
        with SKIP_TSV.open("w", encoding="utf-8") as f:
            f.write("sample_name\treason\n")
            for name, reason in sorted(discovery_skipped):
                f.write(f"{name}\t{reason}\n")

    if csv_filtered:
        with FILTERED_TSV.open("w", encoding="utf-8") as f:
            f.write("sample_name\treason\n")
            for name, reason in sorted(csv_filtered):
                f.write(f"{name}\t{reason}\n")

    # Merge this run's failures into the TSV; never drop prior data failures.
    write_failed_tsv(failed_tsv, failed)

    print(f"Discovered {discovered_total} SSD sample(s); running [{first}–{last}] ({len(subset)} to process)")
    if args.pending_mitoz:
        print(
            f"Pending MitoZ: {len(runs)} to process, "
            f"{len(pending_mitoz_skipped)} skipped (complete or prior failure)"
        )
    elif not args.all_samples:
        print(
            f"CSV filter ({', '.join('true' if v else 'false' for v in allowed_success)}): "
            f"{len(runs)} matched, {len(csv_filtered)} excluded"
        )
    if skipped_completed:
        print(f"Skipped {skipped_completed} already-completed sample(s)")
    if skipped_prior_failures:
        print(f"Skipped {skipped_prior_failures} prior pipeline failure(s) (not retrying)")
    if discovery_skipped:
        print(f"Skipped {len(discovery_skipped)} dir(s) without read pairs (see {SKIP_TSV})")
    if csv_filtered:
        print(f"Excluded {len(csv_filtered)} dir(s) by assembly_success filter (see {FILTERED_TSV})")

    if args.dry_run:
        print(f"Dry run complete. Log: {log_path}")
        return

    if aborted:
        print(
            "Batch aborted by trim/SSD watchdog; reconnect SSD and restart. "
            f"See {LOG}",
            file=sys.stderr,
        )
        raise SystemExit(2)

    if failed:
        print(f"Finished with {len(failed)} failure(s); see {FAILED_TSV} and {LOG}", file=sys.stderr)
        raise SystemExit(1)

    print(f"Done: {succeeded}/{len(subset)} succeeded. Log: {LOG}")


if __name__ == "__main__":
    main()
