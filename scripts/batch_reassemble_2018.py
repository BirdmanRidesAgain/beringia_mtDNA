#!/usr/bin/env python3
"""
Re-assemble the 93 mitogenomes with additional 2018 resequencing data (Table 1,
``multiple datasets?`` = True).

Reads merged FASTQs from the Seagate raw_data volume and writes into existing
``mitogenomes_output/<New_directory_name>/`` trees, replacing prior
``trim_*`` and ``novoplasty_*`` outputs.

Parallelism: run multiple samples at once (``--jobs``). Each sample uses
``--threads-per-job`` for fastp (NOVOPlasty is mostly single-threaded), so total
CPU use is roughly ``jobs * threads_per_job``.

Example::

    python scripts/batch_reassemble_2018.py --dry-run
    python scripts/batch_reassemble_2018.py --jobs 4 --threads-per-job 2 --fresh-log
    python scripts/batch_reassemble_2018.py --skip-completed --jobs 4
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from batch_assemble_beringia2018 import (  # noqa: E402
    assemble_child_environ,
    assemble_cmd,
    find_read_pair,
    mitogenome_python,
    novoplasty_done,
)
from reassembly_2018_common import (  # noqa: E402
    CHECKLIST,
    FASTQ_ROOT,
    OUT_ROOT,
    REFSEQ,
    clear_assembly_artifacts,
    parse_reassembly_samples,
)

LOG = OUT_ROOT / "batch_reassemble_2018.log"
FAILED_TSV = OUT_ROOT / "batch_reassemble_2018_failed.tsv"


@dataclass(frozen=True)
class SampleTask:
    idx: int
    total: int
    field: str
    out_name: str
    fastq_root: str
    threads_per_job: int
    verbose: bool
    no_clean: bool
    dry_run: bool


@dataclass
class SampleResult:
    idx: int
    field: str
    out_name: str
    ok: bool
    exit_code: int
    message: str
    log_block: str


def default_jobs(threads_per_job: int) -> int:
    cpus = os.cpu_count() or 4
    return max(1, cpus // max(1, threads_per_job))


def run_sample(task: SampleTask) -> SampleResult:
    """Worker: assemble one sample (runs in a child process when jobs > 1)."""
    lines: list[str] = []
    field = task.field
    out_name = task.out_name
    fastq_root = Path(task.fastq_root)
    out_dir = OUT_ROOT / out_name

    def add(line: str = "") -> None:
        lines.append(line)

    add(f"\n[{task.idx}/{task.total}] {field} -> {out_name}")
    out_dir.mkdir(parents=True, exist_ok=True)

    try:
        fq1, fq2 = find_read_pair(field, fastq_root)
    except FileNotFoundError as exc:
        add(f"  ERROR: {exc}")
        return SampleResult(
            task.idx,
            field,
            out_name,
            ok=False,
            exit_code=-1,
            message=str(exc),
            log_block="\n".join(lines) + "\n",
        )

    add(f"  R1: {fq1}")
    add(f"  R2: {fq2}")
    if not task.no_clean:
        removed = clear_assembly_artifacts(out_dir, field)
        if removed:
            add(f"  Removed: {', '.join(removed)}")

    cmd = assemble_cmd(
        field,
        out_dir,
        fq1,
        fq2,
        verbose=task.verbose,
        threads=task.threads_per_job,
    )
    add(f"  CMD: {' '.join(cmd)}")

    if task.dry_run:
        return SampleResult(
            task.idx,
            field,
            out_name,
            ok=True,
            exit_code=0,
            message="dry-run",
            log_block="\n".join(lines) + "\n",
        )

    r = subprocess.run(cmd, env=assemble_child_environ())
    add(f"  Exit code: {r.returncode}")
    if r.returncode == 0:
        return SampleResult(
            task.idx,
            field,
            out_name,
            ok=True,
            exit_code=0,
            message="ok",
            log_block="\n".join(lines) + "\n",
        )
    return SampleResult(
        task.idx,
        field,
        out_name,
        ok=False,
        exit_code=r.returncode,
        message="assembly failed",
        log_block="\n".join(lines) + "\n",
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Re-assemble 93 Beringia samples from merged 2018 resequencing FASTQs."
    )
    p.add_argument(
        "--field",
        action="append",
        dest="fields",
        metavar="FIELD",
        help="Only this field number (repeatable).",
    )
    p.add_argument("--fresh-log", action="store_true", help="Overwrite log file.")
    p.add_argument("-v", "--verbose", action="store_true", help="Stream tool output.")
    p.add_argument("--dry-run", action="store_true", help="Print planned work only.")
    p.add_argument(
        "--fastq-root",
        type=Path,
        default=FASTQ_ROOT,
        help=f"Directory of per-field FASTQ folders (default: {FASTQ_ROOT}).",
    )
    p.add_argument(
        "--no-clean",
        action="store_true",
        help="Do not remove existing trim_/novoplasty_* before re-run.",
    )
    p.add_argument(
        "--skip-completed",
        action="store_true",
        help="Skip samples that already have a non-empty NOVOPlasty FASTA.",
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Samples to assemble in parallel (default: CPU count // threads-per-job). "
            "NOVOPlasty is mostly single-threaded; use 3–4 jobs if RAM is ample."
        ),
    )
    p.add_argument(
        "--threads-per-job",
        type=int,
        default=2,
        metavar="T",
        help="Threads for fastp per sample (default: 2). Passed as assemble_mitogenome.py -t.",
    )
    return p.parse_args()


def main() -> None:
    args = parse_args()
    fastq_root = args.fastq_root.expanduser().resolve()
    threads_per_job = max(1, args.threads_per_job)
    jobs = args.jobs if args.jobs is not None else default_jobs(threads_per_job)
    jobs = max(1, jobs)

    if not REFSEQ.is_file():
        raise SystemExit(f"Reference FASTA missing: {REFSEQ}")
    if not fastq_root.is_dir():
        raise SystemExit(f"FASTQ root missing: {fastq_root}")

    mitogenome_python()

    samples = parse_reassembly_samples(CHECKLIST)
    if args.fields:
        wanted = set(args.fields)
        samples = [(f, d) for f, d in samples if f in wanted]
        missing = wanted - {f for f, _ in samples}
        if missing:
            raise SystemExit(f"Unknown or non-reassembly field(s): {', '.join(sorted(missing))}")

    if not samples:
        raise SystemExit("No samples to assemble.")

    skipped_done = 0
    if args.skip_completed:
        pending: list[tuple[str, str]] = []
        for field, out_name in samples:
            if novoplasty_done(OUT_ROOT / out_name, field):
                skipped_done += 1
            else:
                pending.append((field, out_name))
        samples = pending

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    log_mode = "w" if args.fresh_log else "a"
    total = len(samples)
    failed: list[tuple[str, int, str]] = []
    succeeded = 0

    tasks = [
        SampleTask(
            idx=idx,
            total=total,
            field=field,
            out_name=out_name,
            fastq_root=str(fastq_root),
            threads_per_job=threads_per_job,
            verbose=args.verbose,
            no_clean=args.no_clean,
            dry_run=args.dry_run,
        )
        for idx, (field, out_name) in enumerate(samples, start=1)
    ]

    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(f"FASTQ root: {fastq_root}\n")
        log.write(f"Samples: {total}  (skipped already done: {skipped_done})\n")
        log.write(f"Parallel jobs: {jobs}  threads per job: {threads_per_job}\n")
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")
        log.flush()

        if jobs == 1 or args.dry_run:
            for task in tasks:
                result = run_sample(task)
                log.write(result.log_block)
                log.flush()
                if result.ok:
                    succeeded += 1
                else:
                    failed.append((result.field, result.exit_code, result.message))
        else:
            with ProcessPoolExecutor(max_workers=jobs) as pool:
                futures = {pool.submit(run_sample, task): task for task in tasks}
                for fut in as_completed(futures):
                    result = fut.result()
                    log.write(result.log_block)
                    log.flush()
                    if result.ok:
                        succeeded += 1
                    else:
                        failed.append((result.field, result.exit_code, result.message))

        log.write(
            f"\n--- END {datetime.now(timezone.utc).isoformat()}Z ---\n"
            f"Summary: {succeeded} succeeded, {len(failed)} failed, "
            f"{skipped_done} skipped (already assembled)\n"
        )
        for field, code, msg in sorted(failed):
            log.write(f"  failed [{code}] {field}: {msg}\n")

    with FAILED_TSV.open("w", encoding="utf-8") as f:
        f.write("field\texit_code\tmessage\n")
        for field, code, msg in sorted(failed):
            f.write(f"{field}\t{code}\t{msg}\n")

    cpus = os.cpu_count() or 4
    print(
        f"Config: {jobs} parallel sample(s) × {threads_per_job} thread(s)/sample "
        f"(~{jobs * threads_per_job} CPUs vs {cpus} available)"
    )

    if args.dry_run:
        print(f"Dry run: {total} sample(s). See {LOG}")
        return

    if failed:
        print(f"Finished with {len(failed)} failure(s); see {FAILED_TSV} and {LOG}", file=sys.stderr)
        raise SystemExit(1)

    print(f"Done: {succeeded}/{total} assemblies. Log: {LOG}")


if __name__ == "__main__":
    main()
