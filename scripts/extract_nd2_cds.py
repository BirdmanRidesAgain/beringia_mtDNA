#!/usr/bin/env python3
"""Extract ND2 CDS sequences for samples with ND2_present=True in the results CSV."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mitogenome_paths import (  # noqa: E402
    extract_nd2_records_from_gene_fasta,
    find_gene_fasta_for_sample,
    resolve_sample_dir,
)

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_CSV = ROOT / "mitogenome_assembly_results.csv"


def _truthy(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def write_nd2_fasta(
    out_path: Path,
    records: list[tuple[str, str]],
    *,
    sample_prefix: str,
) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as fh:
        for idx, (header, seq) in enumerate(records):
            if len(records) == 1:
                fh.write(f">{sample_prefix}_ND2 {header}\n")
            else:
                fh.write(f">{sample_prefix}_ND2_{idx + 1} {header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def extract_one(
    directory_name: str,
    assembler: str,
    *,
    dry_run: bool,
    overwrite: bool,
) -> tuple[str, str]:
    sample_dir = resolve_sample_dir(directory_name)
    if sample_dir is None:
        return "missing_dir", f"sample directory not found for {directory_name}"

    gene_fasta = find_gene_fasta_for_sample(sample_dir, assembler)
    if gene_fasta is None:
        return "missing_gene_fasta", f"no *.gbf.gene.fasta found under {sample_dir}"

    nd2_records = extract_nd2_records_from_gene_fasta(gene_fasta)
    if not nd2_records:
        return "missing_nd2", f"no ND2 record in {gene_fasta}"

    out_path = sample_dir / f"{sample_dir.name}_ND2_CDS.fa"
    if out_path.exists() and not overwrite:
        return "skipped", f"{out_path} already exists"

    if dry_run:
        return "dry_run", f"would write {len(nd2_records)} record(s) from {gene_fasta} -> {out_path}"

    write_nd2_fasta(out_path, nd2_records, sample_prefix=sample_dir.name)
    return "written", str(out_path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"Results CSV (default: {DEFAULT_CSV.name})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report actions without writing FASTA files",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing *_ND2_CDS.fa files",
    )
    args = parser.parse_args()

    rows = list(csv.DictReader(args.csv.open(newline="")))
    targets = [r for r in rows if _truthy(r.get("ND2_present", ""))]
    if not targets:
        print("No rows with ND2_present=True")
        return

    counts: dict[str, int] = {}
    problems: list[tuple[str, str, str]] = []

    for row in targets:
        directory_name = row["directory_name"]
        status, detail = extract_one(
            directory_name,
            row.get("assembler") or "",
            dry_run=args.dry_run,
            overwrite=args.overwrite,
        )
        counts[status] = counts.get(status, 0) + 1
        if status not in {"written", "dry_run", "skipped"}:
            problems.append((directory_name, status, detail))

    print(f"ND2_present=True rows: {len(targets)}")
    for status in sorted(counts):
        print(f"  {status}: {counts[status]}")
    if problems:
        print("\nProblems:")
        for directory_name, status, detail in problems:
            print(f"  {directory_name}\t{status}\t{detail}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
