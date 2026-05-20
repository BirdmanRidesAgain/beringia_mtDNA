#!/usr/bin/env python3
"""
Re-annotate the 93 re-assembled mitogenomes (2018 resequencing set).

Removes prior MitoZ outputs under each assembly's ``annotation/`` subtree and
runs MitoZ via ``annotate_mitogenome.py`` (conda env ``mitoz-x64``).

Parallelism: run multiple annotations at once (``--jobs``). Each uses
``--threads-per-job`` for MitoZ; total CPU use is roughly
``jobs * threads_per_job``.

Example::

    python scripts/batch_reannotate_2018.py --dry-run
    python scripts/batch_reannotate_2018.py --jobs 4 --threads-per-job 2 --fresh-log
    python scripts/batch_reannotate_2018.py --skip-completed --jobs 4
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

from batch_annotate_assemblies import (  # noqa: E402
    ASM_ROOT,
    BATCH_LOGS_DIR,
    AnnotationJob,
    SkippedItem,
    TIER_NAMES,
    TIER_NONE,
    _classify,
    _disambiguator,
    _novoplasty_token,
    annotate_child_environ,
    annotate_cmd,
    mitoz_python,
    output_already_done,
)
from batch_assemble_beringia2018 import novoplasty_done  # noqa: E402
from mitogenome_paths import clear_annotation_output  # noqa: E402
from reassembly_2018_common import (  # noqa: E402
    CHECKLIST,
    clear_annotation_artifacts,
    parse_reassembly_samples,
)

LOG = BATCH_LOGS_DIR / "batch_reannotate_2018.log"
SKIP_TSV = BATCH_LOGS_DIR / "batch_reannotate_2018_skipped.tsv"
FAILED_TSV = BATCH_LOGS_DIR / "batch_reannotate_2018_failed.tsv"


@dataclass(frozen=True)
class AnnotTask:
    idx: int
    total: int
    sample_dir: str
    novoplasty_dir: str
    fasta: str
    sample_label: str
    tier: int
    threads: int
    mitoz_clade: str
    no_clean: bool
    dry_run: bool


@dataclass
class AnnotResult:
    idx: int
    sample_label: str
    ok: bool
    exit_code: int
    duration_s: float
    message: str
    log_block: str


def default_jobs(threads_per_job: int) -> int:
    cpus = os.cpu_count() or 4
    return max(1, cpus // max(1, threads_per_job))


def discover_jobs(
    samples: list[tuple[str, str]],
) -> tuple[list[AnnotationJob], list[SkippedItem]]:
    jobs: list[AnnotationJob] = []
    skipped: list[SkippedItem] = []

    for field, out_name in samples:
        sd = ASM_ROOT / out_name
        if not sd.is_dir():
            skipped.append(SkippedItem(sd, "output directory missing"))
            continue
        if not novoplasty_done(sd, field):
            nov_dir = sd / f"novoplasty_{field}"
            skipped.append(SkippedItem(nov_dir, "assembly failed (no NOVOPlasty FASTA)"))
            continue

        novop = sd / f"novoplasty_{field}"
        tier, cands = _classify(novop)
        if tier == TIER_NONE or not cands:
            skipped.append(SkippedItem(novop, "no annotatable FASTA"))
            continue

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


def job_to_task(job: AnnotationJob, *, idx: int, total: int, **kwargs) -> AnnotTask:
    return AnnotTask(
        idx=idx,
        total=total,
        sample_dir=str(job.sample_dir),
        novoplasty_dir=str(job.novoplasty_dir),
        fasta=str(job.fasta),
        sample_label=job.sample_label,
        tier=job.tier,
        **kwargs,
    )


def run_annotate(task: AnnotTask) -> AnnotResult:
    """Worker: run MitoZ for one assembly (child process when jobs > 1)."""
    lines: list[str] = []

    def add(line: str = "") -> None:
        lines.append(line)

    tname = TIER_NAMES.get(task.tier, "?")
    add(f"\n[{task.idx}/{task.total}] tier={task.tier}({tname}) {task.sample_label}")
    add(f"  fasta: {task.fasta}")

    job = AnnotationJob(
        Path(task.sample_dir),
        Path(task.novoplasty_dir),
        Path(task.fasta),
        task.sample_label,
        task.tier,
    )

    if not task.no_clean:
        removed = clear_annotation_output(job.sample_dir, job.sample_label)
        if removed:
            add(f"  Removed: {', '.join(removed)}")

    if task.dry_run:
        return AnnotResult(
            task.idx,
            task.sample_label,
            ok=True,
            exit_code=0,
            duration_s=0.0,
            message="dry-run",
            log_block="\n".join(lines) + "\n",
        )

    t0 = datetime.now(timezone.utc)
    r = subprocess.run(
        annotate_cmd(job, threads=task.threads, clade=task.mitoz_clade),
        env=annotate_child_environ(),
    )
    dt = (datetime.now(timezone.utc) - t0).total_seconds()
    add(f"  Exit code: {r.returncode}  duration: {dt:.1f}s")

    if r.returncode == 0:
        return AnnotResult(
            task.idx,
            task.sample_label,
            ok=True,
            exit_code=0,
            duration_s=dt,
            message="ok",
            log_block="\n".join(lines) + "\n",
        )
    add("  FAILED (continuing with remaining samples).")
    return AnnotResult(
        task.idx,
        task.sample_label,
        ok=False,
        exit_code=r.returncode,
        duration_s=dt,
        message="annotation failed",
        log_block="\n".join(lines) + "\n",
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Re-annotate 93 Beringia 2018 re-assemblies.")
    p.add_argument(
        "--field",
        action="append",
        dest="fields",
        metavar="FIELD",
        help="Only this field number (repeatable).",
    )
    p.add_argument("--fresh-log", action="store_true", help="Overwrite log file.")
    p.add_argument("--dry-run", action="store_true", help="List jobs only.")
    p.add_argument(
        "--no-clean",
        action="store_true",
        help="Keep existing MitoZ output dirs (default: remove before each run).",
    )
    p.add_argument(
        "--skip-completed",
        action="store_true",
        help="Skip samples that already have MitoZ .result output.",
    )
    p.add_argument(
        "--jobs",
        type=int,
        default=None,
        metavar="N",
        help=(
            "Annotations to run in parallel (default: CPU count // threads-per-job). "
            "MitoZ is memory-heavy; try 2–4 if RAM is limited."
        ),
    )
    p.add_argument(
        "--threads-per-job",
        type=int,
        default=2,
        metavar="T",
        help="MitoZ threads per sample (default: 2).",
    )
    p.add_argument(
        "--threads",
        type=int,
        default=None,
        help=argparse.SUPPRESS,
    )
    p.add_argument("--mitoz-clade", default="Chordata")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    mitoz_python()

    threads_per_job = max(1, args.threads if args.threads is not None else args.threads_per_job)
    jobs_parallel = args.jobs if args.jobs is not None else default_jobs(threads_per_job)
    jobs_parallel = max(1, jobs_parallel)

    samples = parse_reassembly_samples(CHECKLIST)
    if args.fields:
        wanted = set(args.fields)
        samples = [(f, d) for f, d in samples if f in wanted]
        missing = wanted - {f for f, _ in samples}
        if missing:
            raise SystemExit(f"Unknown or non-reassembly field(s): {', '.join(sorted(missing))}")

    annot_jobs, skipped = discover_jobs(samples)
    BATCH_LOGS_DIR.mkdir(parents=True, exist_ok=True)

    with SKIP_TSV.open("w", encoding="utf-8") as f:
        f.write("path\treason\n")
        for s in skipped:
            try:
                rel = s.path.relative_to(REPO)
            except ValueError:
                rel = s.path
            f.write(f"{rel}\t{s.reason}\n")

    if args.dry_run:
        print(f"Jobs: {len(annot_jobs)}, skipped: {len(skipped)}")
        print(
            f"Config: {jobs_parallel} parallel × {threads_per_job} threads "
            f"(~{jobs_parallel * threads_per_job} CPUs)"
        )
        for j in annot_jobs:
            print(f"  tier={j.tier} {j.sample_label} <- {j.fasta.name}")
        for s in skipped:
            print(f"  SKIP {s.path.name}: {s.reason}")
        return

    if not annot_jobs:
        raise SystemExit(f"No annotation jobs; see {SKIP_TSV}")

    skipped_done = 0
    pending = annot_jobs
    if args.skip_completed:
        pending = []
        for j in annot_jobs:
            if output_already_done(j):
                skipped_done += 1
            else:
                pending.append(j)

    total = len(pending)
    if total == 0:
        print(f"All {len(annot_jobs)} annotation(s) already done.")
        return

    tasks = [
        job_to_task(
            j,
            idx=idx,
            total=total,
            threads=threads_per_job,
            mitoz_clade=args.mitoz_clade,
            no_clean=args.no_clean,
            dry_run=False,
        )
        for idx, j in enumerate(pending, start=1)
    ]

    log_mode = "w" if args.fresh_log else "a"
    succeeded = 0
    failed_jobs: list[tuple[str, int]] = []

    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(
            f"Reassembly set: {total} job(s) to run, {len(skipped)} skipped (assembly), "
            f"{skipped_done} skipped (already annotated)\n"
        )
        log.write(f"Parallel jobs: {jobs_parallel}  threads per job: {threads_per_job}\n")
        log.write(f"Skip log: {SKIP_TSV.relative_to(REPO)}\n")
        log.write(f"Failed log: {FAILED_TSV.relative_to(REPO)}\n")
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")
        log.flush()

        if jobs_parallel == 1:
            for task in tasks:
                result = run_annotate(task)
                log.write(result.log_block)
                log.flush()
                if result.ok:
                    succeeded += 1
                else:
                    failed_jobs.append((result.sample_label, result.exit_code))
        else:
            with ProcessPoolExecutor(max_workers=jobs_parallel) as pool:
                futures = {pool.submit(run_annotate, task): task for task in tasks}
                for fut in as_completed(futures):
                    result = fut.result()
                    log.write(result.log_block)
                    log.flush()
                    if result.ok:
                        succeeded += 1
                    else:
                        failed_jobs.append((result.sample_label, result.exit_code))

        log.write(
            f"\n--- END {datetime.now(timezone.utc).isoformat()}Z ---\n"
            f"Summary: {succeeded} succeeded, {len(failed_jobs)} failed, "
            f"{skipped_done} skipped (already annotated), "
            f"{len(skipped)} skipped (no assembly FASTA)\n"
        )
        for label, code in sorted(failed_jobs):
            log.write(f"  failed [{code}]: {label}\n")

    with FAILED_TSV.open("w", encoding="utf-8") as f:
        f.write("sample_label\texit_code\n")
        for label, code in sorted(failed_jobs):
            f.write(f"{label}\t{code}\n")

    cpus = os.cpu_count() or 4
    print(
        f"Config: {jobs_parallel} parallel sample(s) × {threads_per_job} thread(s)/sample "
        f"(~{jobs_parallel * threads_per_job} CPUs vs {cpus} available)"
    )

    if failed_jobs:
        print(
            f"Finished with {len(failed_jobs)} failure(s); see {FAILED_TSV} and {LOG}",
            file=sys.stderr,
        )
        raise SystemExit(1)

    print(f"Done: {succeeded}/{total} annotations. Log: {LOG}")


if __name__ == "__main__":
    main()
