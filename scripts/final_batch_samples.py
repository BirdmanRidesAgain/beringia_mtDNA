"""Table 2 final-batch ingest: sample metadata and Downloads folder aliases."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from batch_assemble_beringia_birds import (  # noqa: F401 re-export
    BeringiaBirdSample,
    find_read_pair_in_dir,
    link_or_copy,
    stage_sample_reads,
)

CHECKLIST_INGEST_PREFIXES = ("final_batch/", "data/final_batch/")

# Downloads subdir name when it differs from canonical directory_name (field # typo).
DOWNLOADS_DIR_ALIASES: dict[str, str] = {
    "Lagopus_lagopus_koreni_UAM_28171_KSW5133": (
        "Lagopus_lagopus_koreni_UAM_28171_KSW5113"
    ),
}


@dataclass(frozen=True)
class FinalBatchSample:
    species: str
    catalog_num: str
    field_num: str
    directory_name: str

    def as_ingest(self) -> BeringiaBirdSample:
        return BeringiaBirdSample(
            species=self.species,
            catalog_num=self.catalog_num,
            field_num=self.field_num,
            directory_name=self.directory_name,
        )


def downloads_sample_dir(downloads_root: Path, directory_name: str) -> Path:
    sub = DOWNLOADS_DIR_ALIASES.get(directory_name, directory_name)
    return downloads_root / sub


def parse_final_batch_samples(checklist: Path) -> list[FinalBatchSample]:
    """Table 1 rows whose old directory is under final_batch/ or data/final_batch/."""
    samples: list[FinalBatchSample] = []
    in_table = False
    for line in checklist.read_text(encoding="utf-8").splitlines():
        if line.startswith("| Species | Catalog"):
            in_table = True
            continue
        if not in_table:
            continue
        if not line.startswith("|"):
            if samples:
                break
            continue
        if line.startswith("| - |"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 7:
            continue
        old_dir = parts[3].strip("`").strip("/")
        if not any(p in old_dir for p in CHECKLIST_INGEST_PREFIXES):
            continue
        directory_name = parts[4].strip("`").strip("/")
        if not directory_name:
            continue
        samples.append(
            FinalBatchSample(
                species=re.sub(r"^\*|\*$", "", parts[0]).strip(),
                catalog_num=parts[1],
                field_num=parts[2],
                directory_name=directory_name,
            )
        )
    return samples


def stage_final_batch_reads(
    sample: FinalBatchSample,
    *,
    downloads_root: Path,
    data_root: Path,
    copy: bool,
    dry_run: bool,
) -> tuple[Path, Path]:
    """Stage reads; resolve Downloads folder aliases when needed."""
    ingest = sample.as_ingest()
    source_dir = downloads_sample_dir(downloads_root, sample.directory_name)
    dest_dir = data_root / sample.directory_name
    r1_src, r2_src = find_read_pair_in_dir(source_dir)
    r1_dest = dest_dir / f"{sample.directory_name}_R1.fastq.gz"
    r2_dest = dest_dir / f"{sample.directory_name}_R2.fastq.gz"
    if dry_run:
        return r1_dest, r2_dest
    link_or_copy(r1_src, r1_dest, copy=copy)
    link_or_copy(r2_src, r2_dest, copy=copy)
    return r1_dest, r2_dest
