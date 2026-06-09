#!/usr/bin/env python3
"""
Stage and assemble mitogenomes for Table 2 final-batch samples in ~/Downloads/final_batch.

Example::

    python scripts/batch_assemble_final_batch.py --dry-run
    python scripts/batch_assemble_final_batch.py --stage-only --update-checklist
    python scripts/batch_assemble_final_batch.py --skip-stage --fresh-log
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from batch_assemble_beringia2018 import (  # noqa: E402
    CHECKLIST,
    REFSEQ,
    assemble_child_environ,
    assemble_cmd,
    novoplasty_done,
)
from batch_assemble_beringia_birds import find_read_pair_in_dir, mitogenome_python  # noqa: E402
from final_batch_samples import (  # noqa: E402
    parse_final_batch_samples,
    stage_final_batch_reads,
)

REPO = Path(__file__).resolve().parent.parent
DOWNLOADS_ROOT = Path.home() / "Downloads/final_batch"
DATA_ROOT = REPO / "data/final_batch"
OUT_ROOT = REPO / "mitogenomes_output"
LOG = OUT_ROOT / "_batch_logs/batch_assemble_final_batch.log"

OLD_DIR_PREFIX = "final_batch/"
STAGED_DIR_PREFIX = "data/final_batch/"


def update_checklist_staged_paths(checklist: Path) -> int:
    text = checklist.read_text(encoding="utf-8")
    updated = 0
    lines: list[str] = []
    for line in text.splitlines():
        if line.startswith("|") and OLD_DIR_PREFIX in line:
            new_line = line.replace(f"`{OLD_DIR_PREFIX}", f"`{STAGED_DIR_PREFIX}")
            if new_line != line:
                updated += 1
            lines.append(new_line)
        else:
            lines.append(line)
    if updated:
        checklist.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return updated


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Stage and assemble final_batch mitogenomes.")
    p.add_argument("--field", action="append", dest="fields", metavar="FIELD")
    p.add_argument("--fresh-log", action="store_true")
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--skip-existing", action="store_true")
    p.add_argument("--stage-only", action="store_true")
    p.add_argument("--skip-stage", action="store_true")
    p.add_argument("--copy", action="store_true")
    p.add_argument("--downloads-root", type=Path, default=DOWNLOADS_ROOT)
    p.add_argument("--data-root", type=Path, default=DATA_ROOT)
    p.add_argument("--update-checklist", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    downloads_root = args.downloads_root.expanduser().resolve()
    data_root = args.data_root.resolve()

    if not args.skip_stage and not downloads_root.is_dir():
        raise SystemExit(f"Downloads directory missing: {downloads_root}")
    if not REFSEQ.is_file():
        raise SystemExit(f"Reference FASTA missing: {REFSEQ}")

    samples = parse_final_batch_samples(CHECKLIST)
    if args.fields:
        wanted = set(args.fields)
        samples = [s for s in samples if s.field_num in wanted]
        missing = wanted - {s.field_num for s in samples}
        if missing:
            raise SystemExit(f"Unknown or non-ingest field(s): {', '.join(sorted(missing))}")

    if not samples:
        raise SystemExit(
            "No final_batch samples in Table 1 (old directory must be final_batch/ or data/final_batch/)."
        )

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    log_mode = "w" if args.fresh_log else "a"
    total = len(samples)
    staged_any = False

    with LOG.open(log_mode, encoding="utf-8") as log:
        log.write(f"Downloads root: {downloads_root}\nData root: {data_root}\nSamples: {total}\n")
        log.write(f"\n--- START {datetime.now(timezone.utc).isoformat()}Z ---\n")

        for idx, sample in enumerate(samples, start=1):
            out_dir = OUT_ROOT / sample.directory_name
            log.write(f"\n[{idx}/{total}] {sample.field_num} -> {sample.directory_name}\n")
            log.flush()

            if not args.skip_stage:
                if args.dry_run:
                    from final_batch_samples import downloads_sample_dir

                    src = downloads_sample_dir(downloads_root, sample.directory_name)
                    r1_src, r2_src = find_read_pair_in_dir(src)
                    r1_dest = data_root / sample.directory_name / (
                        f"{sample.directory_name}_R1.fastq.gz"
                    )
                    r2_dest = r1_dest.with_name(f"{sample.directory_name}_R2.fastq.gz")
                    log.write(f"  Stage: {r1_src} -> {r1_dest}\n  Stage: {r2_src} -> {r2_dest}\n")
                else:
                    r1_dest, r2_dest = stage_final_batch_reads(
                        sample,
                        downloads_root=downloads_root,
                        data_root=data_root,
                        copy=args.copy,
                        dry_run=False,
                    )
                    staged_any = True
                    log.write(f"  Staged R1: {r1_dest}\n  Staged R2: {r2_dest}\n")
            else:
                r1_dest, r2_dest = find_read_pair_in_dir(data_root / sample.directory_name)

            if args.stage_only or args.dry_run:
                if args.stage_only and not args.dry_run:
                    log.write("  Stage-only; assembly skipped.\n")
                continue

            cmd = assemble_cmd(
                sample.field_num, out_dir, r1_dest, r2_dest, verbose=args.verbose
            )
            log.write(f"  CMD: {' '.join(cmd)}\n")
            log.flush()
            if args.dry_run:
                continue
            if args.skip_existing and novoplasty_done(out_dir, sample.field_num):
                log.write("  Skipped: existing NOVOPlasty FASTA found.\n")
                continue
            r = subprocess.run(cmd, env=assemble_child_environ())
            log.write(f"  Exit code: {r.returncode}\n")
            if r.returncode != 0:
                raise SystemExit(
                    f"Assembly failed for {sample.field_num} (exit {r.returncode}); see {LOG}"
                )

        if args.update_checklist and staged_any and not args.dry_run:
            n = update_checklist_staged_paths(CHECKLIST)
            log.write(f"\nChecklist: updated {n} path(s) to {STAGED_DIR_PREFIX}...\n")


if __name__ == "__main__":
    main()
