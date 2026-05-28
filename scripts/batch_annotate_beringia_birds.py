#!/usr/bin/env python3
"""
Batch-annotate mitogenomes for the beringia_birds ingest set (Table 1, data/beringia_birds/).

Uses MitoZ via ``annotate_mitogenome.py`` and conda env ``mitoz-x64``.
Skips samples with no annotatable NOVOPlasty FASTA (failed assemblies).

Example::

    python scripts/batch_annotate_beringia_birds.py --dry-run
    python scripts/batch_annotate_beringia_birds.py --fresh-log
    python scripts/batch_annotate_beringia_birds.py --field KSW2408
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
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
from batch_assemble_beringia2018 import CHECKLIST, novoplasty_done  # noqa: E402
from batch_assemble_beringia_birds import (  # noqa: E402
    BeringiaBirdSample,
    parse_beringia_birds_samples,
)

LOG = BATCH_LOGS_DIR / "batch_annotate_beringia_birds.log"
SKIP_TSV = BATCH_LOGS_DIR / "batch_annotate_beringia_birds_skipped.tsv"
FAILED_TSV = BATCH_LOGS_DIR / "batch_annotate_beringia_birds_failed.tsv"


def discover_beringia_bird_jobs(
    samples: list[BeringiaBirdSample],
) -> tuple[list[AnnotationJob], list[SkippedItem]]:
    jobs: list[AnnotationJob] = []
    skipped: list[SkippedItem] = []

    for sample in samples:
        sd = ASM_ROOT / sample.directory_name
        if not sd.is_dir():
            skipped.append(SkippedItem(sd, "output directory missing"))
            continue
        if not novoplasty_done(sd, sample.field_num):
            nov_dir = sd / f"novoplasty_{sample.field_num}"
            skipped.append(SkippedItem(nov_dir, "assembly failed (no NOVOPlasty FASTA)"))
            continue

        novop = sd / f"novoplasty_{sample.field_num}"
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


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="MitoZ annotation for beringia_birds assemblies."
    )
    p.add_argument(
        "--field",
        action="append",
        dest="fields",
        metavar="FIELD",
        help="Annotate only this field number (repeatable).",
    )
    p.add_argument("--fresh-log", action="store_true", help="Overwrite log file.")
    p.add_argument("--dry-run", action="store_true", help="List jobs only.")
    p.add_argument("--force", action="store_true", help="Re-run even if MitoZ output exists.")
    p.add_argument(
        "--threads",
        type=int,
        default=max(1, (os.cpu_count() or 1) // 2),
    )
    p.add_argument("--mitoz-clade", default="Chordata")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    mitoz_python()  # fail fast if env missing

    samples = parse_beringia_birds_samples(CHECKLIST)
    if args.fields:
        wanted = set(args.fields)
        samples = [s for s in samples if s.field_num in wanted]
        missing = wanted - {s.field_num for s in samples}
        if missing:
            raise SystemExit(f"Unknown or non-ingest field(s): {', '.join(sorted(missing))}")

    if not samples:
        raise SystemExit("No beringia_birds samples found in Table 1 of the checklist.")

    jobs, skipped = discover_beringia_bird_jobs(samples)
    BATCH_LOGS_DIR.mkdir(parents=True, exist_ok=True)
    with SKIP_TSV.open("w", encoding="utf-8") as f:
        f.write("path\treason\n")
        for s in skipped:
            try:
                rel = s.path.relative_to(REPO)
            except ValueError:
                rel = s.path
            f.write(f"{rel}\t{s.reason}\n")

    total = len(jobs)
    if args.dry_run:
        print(f"Samples: {len(samples)}, jobs: {total}, skipped: {len(skipped)}")
        for j in jobs:
            print(f"  tier={j.tier} {j.sample_label} <- {j.fasta.name}")
        for s in skipped:
            print(f"  SKIP {s.path.name}: {s.reason}")
        return

    if total == 0:
        raise SystemExit(f"No annotation jobs; see {SKIP_TSV}")

    log_mode = "w" if args.fresh_log else "a"
    succeeded = 0
    failed_jobs: list[tuple[str, int]] = []
    skipped_done = 0

    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(
            f"Beringia birds: {len(samples)} sample(s), {total} job(s), "
            f"{len(skipped)} skipped (assembly)\n"
        )
        log.write(f"Skip log: {SKIP_TSV.relative_to(REPO)}\n")
        log.write(f"Failed log: {FAILED_TSV.relative_to(REPO)}\n")
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")
        for idx, job in enumerate(jobs, start=1):
            tname = TIER_NAMES.get(job.tier, "?")
            log.write(f"\n[{idx}/{total}] tier={job.tier}({tname}) {job.sample_label}\n")
            log.write(f"  fasta: {job.fasta}\n")
            log.flush()
            if (not args.force) and output_already_done(job):
                log.write("  Skipped: existing MitoZ output.\n")
                log.flush()
                skipped_done += 1
                continue
            t0 = datetime.now(timezone.utc)
            r = subprocess.run(
                annotate_cmd(job, threads=args.threads, clade=args.mitoz_clade),
                env=annotate_child_environ(),
            )
            dt = (datetime.now(timezone.utc) - t0).total_seconds()
            log.write(f"  Exit code: {r.returncode}  duration: {dt:.1f}s\n")
            log.flush()
            if r.returncode == 0:
                succeeded += 1
            else:
                failed_jobs.append((job.sample_label, r.returncode))
                log.write("  FAILED (continuing with remaining samples).\n")
                log.flush()

        log.write(
            f"\n--- END {datetime.now(timezone.utc).isoformat()}Z ---\n"
            f"Summary: {succeeded} succeeded, {len(failed_jobs)} failed, "
            f"{skipped_done} skipped (already annotated), "
            f"{len(skipped)} skipped (no assembly FASTA)\n"
        )
        for label, code in failed_jobs:
            log.write(f"  failed [{code}]: {label}\n")

    with FAILED_TSV.open("w", encoding="utf-8") as f:
        f.write("sample_label\texit_code\n")
        for label, code in failed_jobs:
            f.write(f"{label}\t{code}\n")

    if failed_jobs:
        print(
            f"Finished with {len(failed_jobs)} failure(s); see {FAILED_TSV} and {LOG}",
            file=sys.stderr,
        )
        raise SystemExit(1)


if __name__ == "__main__":
    main()
