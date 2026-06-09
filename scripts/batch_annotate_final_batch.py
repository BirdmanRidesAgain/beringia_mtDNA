#!/usr/bin/env python3
"""MitoZ annotation for final_batch assemblies (Table 1, data/final_batch/)."""

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
from final_batch_samples import FinalBatchSample, parse_final_batch_samples  # noqa: E402

LOG = BATCH_LOGS_DIR / "batch_annotate_final_batch.log"
SKIP_TSV = BATCH_LOGS_DIR / "batch_annotate_final_batch_skipped.tsv"
FAILED_TSV = BATCH_LOGS_DIR / "batch_annotate_final_batch_failed.tsv"


def discover_jobs(
    samples: list[FinalBatchSample],
) -> tuple[list[AnnotationJob], list[SkippedItem]]:
    jobs: list[AnnotationJob] = []
    skipped: list[SkippedItem] = []
    for sample in samples:
        sd = ASM_ROOT / sample.directory_name
        if not sd.is_dir():
            skipped.append(SkippedItem(sd, "output directory missing"))
            continue
        if not novoplasty_done(sd, sample.field_num):
            skipped.append(
                SkippedItem(
                    sd / f"novoplasty_{sample.field_num}",
                    "assembly failed (no NOVOPlasty FASTA)",
                )
            )
            continue
        novop = sd / f"novoplasty_{sample.field_num}"
        tier, cands = _classify(novop)
        if tier == TIER_NONE or not cands:
            skipped.append(SkippedItem(novop, "no annotatable FASTA"))
            continue
        token = _novoplasty_token(novop)
        multi = len(cands) > 1
        for fasta in cands:
            label = (
                f"{sd.name}__{token}__{_disambiguator(fasta)}"
                if multi
                else f"{sd.name}__{token}"
            )
            jobs.append(AnnotationJob(sd, novop, fasta, label, tier))
    jobs.sort(key=lambda j: j.sample_label)
    return jobs, skipped


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="MitoZ annotation for final_batch assemblies.")
    p.add_argument("--field", action="append", dest="fields", metavar="FIELD")
    p.add_argument("--fresh-log", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 1) // 2))
    p.add_argument("--mitoz-clade", default="Chordata")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    mitoz_python()

    samples = parse_final_batch_samples(CHECKLIST)
    if args.fields:
        wanted = set(args.fields)
        samples = [s for s in samples if s.field_num in wanted]
        missing = wanted - {s.field_num for s in samples}
        if missing:
            raise SystemExit(f"Unknown field(s): {', '.join(sorted(missing))}")
    if not samples:
        raise SystemExit("No final_batch samples in checklist Table 1.")

    jobs, skipped = discover_jobs(samples)
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
        print(f"Samples: {len(samples)}, jobs: {len(jobs)}, skipped: {len(skipped)}")
        for j in jobs:
            print(f"  tier={j.tier} {j.sample_label} <- {j.fasta.name}")
        for s in skipped:
            print(f"  SKIP {s.path.name}: {s.reason}")
        return
    if not jobs:
        raise SystemExit(f"No annotation jobs; see {SKIP_TSV}")

    log_mode = "w" if args.fresh_log else "a"
    succeeded = 0
    failed: list[tuple[str, int]] = []
    skipped_done = 0

    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(f"Final batch: {len(samples)} samples, {len(jobs)} jobs\n")
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")
        for idx, job in enumerate(jobs, start=1):
            log.write(f"\n[{idx}/{len(jobs)}] {job.sample_label}\n  fasta: {job.fasta}\n")
            if (not args.force) and output_already_done(job):
                log.write("  Skipped: existing MitoZ output.\n")
                skipped_done += 1
                continue
            r = subprocess.run(
                annotate_cmd(job, threads=args.threads, clade=args.mitoz_clade),
                env=annotate_child_environ(),
            )
            log.write(f"  Exit code: {r.returncode}\n")
            if r.returncode == 0:
                succeeded += 1
            else:
                failed.append((job.sample_label, r.returncode))

    with FAILED_TSV.open("w", encoding="utf-8") as f:
        f.write("sample_label\texit_code\n")
        for label, code in failed:
            f.write(f"{label}\t{code}\n")
    if failed:
        raise SystemExit(f"{len(failed)} failure(s); see {FAILED_TSV} and {LOG}")


if __name__ == "__main__":
    main()
