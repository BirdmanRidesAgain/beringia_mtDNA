#!/usr/bin/env python3
"""Rebuild mitogenome_assembly_results.csv from mitogenomes_assembled_checklist.md Table 1."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mitogenome_paths import (  # noqa: E402
    ASM_ROOT,
    annotation_success_for_sample_dir,
    assembly_level_for_sample_dir,
)

ROOT = Path(__file__).resolve().parent.parent
CHECKLIST = ROOT / "mitogenomes_assembled_checklist.md"
CSV_IN = ROOT / "mitogenome_assembly_results.csv"
CSV_OUT = ROOT / "mitogenome_assembly_results.csv"

GENUS_TO_ORDER = {
    "Anser": "Anseriformes",
    "Spatula": "Anseriformes",
    "Mareca": "Anseriformes",
    "Anas": "Anseriformes",
    "Aythya": "Anseriformes",
    "Somateria": "Anseriformes",
    "Histrionicus": "Anseriformes",
    "Melanitta": "Anseriformes",
    "Clangula": "Anseriformes",
    "Mergus": "Anseriformes",
    "Lagopus": "Galliformes",
    "Pluvialis": "Charadriiformes",
    "Numenius": "Charadriiformes",
    "Arenaria": "Charadriiformes",
    "Calidris": "Charadriiformes",
    "Gallinago": "Charadriiformes",
    "Tringa": "Charadriiformes",
    "Uria": "Charadriiformes",
    "Larus": "Charadriiformes",
    "Gavia": "Gaviiformes",
    "Picoides": "Piciformes",
    "Pica": "Passeriformes",
    "Corvus": "Passeriformes",
    "Phylloscopus": "Passeriformes",
    "Luscinia": "Passeriformes",
    "Motacilla": "Passeriformes",
    "Anthus": "Passeriformes",
    "Pinicola": "Passeriformes",
    "Leucosticte": "Passeriformes",
    "Calcarius": "Passeriformes",
    "Plectrophenax": "Passeriformes",
}

# Legacy CSV directory_name -> canonical checklist directory_name
LEGACY_DIRECTORY_ALIASES = {
    "UAM29835_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Calcarius_lapponicus_alascensis_UAM_29835",
    "DDG1922_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Larus_brachyrhynchus_UAM_14803_DDG1922",
    "UAM10581_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Picoides_tridactylus_UAM_10581",
    "TransB_SVD_422_ACTCCATC-CTCTGGTT__Taeniopygia_guttata_NC_007897.1_mitogenome.fa": (
        "Pinicola_enucleator_kamtschatschensis_UWBM_51642_SVD422"
    ),
    "JJW2478_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Spatula_clypeata_UAM_35604_JJW2478",
    "Mareca_penelope_UAM9759_DDG1703": "Mareca_penelope_UAM_9359_DDG1703",
}

# Checklist directory_name -> on-disk mitogenomes_output folder (when names differ)
DISK_DIRECTORY_ALIASES: dict[str, str] = {}

_DIR_SUFFIX_RE = re.compile(r"_(UAM|UWBM)_(\d+)_([^_]+)$")

OUT_COLUMNS = [
    "directory_name",
    "species",
    "assembly_level",
    "annotation_success",
    "catalog_num",
    "field_num",
    "coverage",
    "order",
    "taxon_pair",
]


def norm_catalog(value: str) -> str:
    return re.sub(r"\s+", "", value.strip())


def norm_field(value: str) -> str:
    return value.strip()


def match_key(catalog: str, field: str) -> tuple[str, str]:
    return norm_catalog(catalog), norm_field(field)


def sanitize_species_for_dir(species: str) -> str:
    value = re.sub(r"\s+", "_", species.strip())
    return re.sub(r"[^A-Za-z0-9_]", "", value)


def build_directory_name(species: str, catalog: str, field: str) -> str:
    cat = norm_catalog(catalog)
    fld = norm_field(field)
    parts = cat.split(maxsplit=1) if " " in catalog else [cat]
    if len(parts) == 1:
        museum_num = parts[0]
        if museum_num.startswith("UAM"):
            museum, num = "UAM", museum_num[3:]
        elif museum_num.startswith("UWBM"):
            museum, num = "UWBM", museum_num[4:]
        else:
            museum, num = museum_num, ""
        prefix = f"{sanitize_species_for_dir(species)}_{museum}_{num}"
    else:
        prefix = f"{sanitize_species_for_dir(species)}_{cat.replace(' ', '_')}"
    if fld:
        return f"{prefix}_{fld}"
    return prefix


def parse_checklist(path: Path) -> list[dict]:
    rows = []
    in_table = False
    for line in path.read_text().splitlines():
        if line.startswith("| Species | Catalog"):
            in_table = True
            continue
        if not in_table:
            continue
        if not line.startswith("|"):
            if rows:
                break
            continue
        if line.startswith("| - |"):
            continue
        parts = [p.strip() for p in line.split("|")[1:-1]]
        if len(parts) < 7:
            continue
        species = re.sub(r"^\*|\*$", "", parts[0]).strip()
        catalog = parts[1]
        field = parts[2]
        directory = parts[4].strip("`").rstrip("/")
        if not directory:
            directory = build_directory_name(species, catalog, field)
        rows.append(
            {
                "species": species,
                "catalog_num": norm_catalog(catalog),
                "field_num": norm_field(field),
                "directory_name": directory,
            }
        )
    return rows


def directory_suffix(directory_name: str) -> str | None:
    match = _DIR_SUFFIX_RE.search(directory_name)
    if not match:
        return None
    museum, number, field = match.groups()
    return f"_{museum}_{number}_{field}"


def resolve_disk_directory(directory_name: str) -> Path | None:
    """Return existing mitogenomes_output sample dir (checklist or disk alias)."""
    if directory_name in DISK_DIRECTORY_ALIASES:
        alias = ASM_ROOT / DISK_DIRECTORY_ALIASES[directory_name]
        if alias.is_dir():
            return alias
    direct = ASM_ROOT / directory_name
    if direct.is_dir():
        return direct
    suffix = directory_suffix(directory_name)
    if suffix:
        matches = [p for p in ASM_ROOT.iterdir() if p.is_dir() and p.name.endswith(suffix)]
        if len(matches) == 1:
            return matches[0]
    return None


def add_order(species: str) -> str:
    genus = species.split()[0] if species else ""
    return GENUS_TO_ORDER.get(genus, "")


def load_csv_index(path: Path) -> tuple[dict[tuple[str, str], dict], dict[str, dict]]:
    by_key: dict[tuple[str, str], dict] = {}
    by_dir: dict[str, dict] = {}
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            canon = LEGACY_DIRECTORY_ALIASES.get(row["directory_name"], row["directory_name"])
            cat = row.get("catalog_guess") or row.get("catalog_num") or ""
            fld = row.get("accession_guess") or row.get("field_num") or ""
            if cat and fld:
                by_key[match_key(cat, fld)] = row
            by_dir[canon] = row
            by_dir[row["directory_name"]] = row
    return by_key, by_dir


def main() -> None:
    checklist = parse_checklist(CHECKLIST)
    by_key, by_dir = load_csv_index(CSV_IN)

    out_rows = []
    level_changes = []
    species_mismatches = []
    dir_missing_on_disk = []

    for chk in checklist:
        key = match_key(chk["catalog_num"], chk["field_num"])
        old = by_dir.get(chk["directory_name"]) or by_key.get(key)
        if old:
            coverage = old.get("coverage") or ""
            taxon_pair = old.get("taxon_pair") or ""
            old_level = old.get("assembly_level") or ""
            old_species = old.get("species") or ""
        else:
            coverage = taxon_pair = old_level = old_species = ""

        disk_path = resolve_disk_directory(chk["directory_name"])
        new_level = assembly_level_for_sample_dir(disk_path) if disk_path else "failed"
        annotation_success = (
            annotation_success_for_sample_dir(disk_path) if disk_path else "not_attempted"
        )
        if disk_path is None:
            dir_missing_on_disk.append(chk["directory_name"])

        if old_level and new_level and old_level != new_level:
            level_changes.append((chk["directory_name"], old_level, new_level))

        if old_species and old_species != chk["species"] and old_species != "Unknown":
            species_mismatches.append((chk["directory_name"], old_species, chk["species"]))

        out_rows.append(
            {
                "directory_name": chk["directory_name"],
                "species": chk["species"],
                "assembly_level": new_level,
                "annotation_success": annotation_success,
                "catalog_num": chk["catalog_num"],
                "field_num": chk["field_num"],
                "coverage": coverage,
                "order": add_order(chk["species"]),
                "taxon_pair": taxon_pair,
            }
        )

    with CSV_OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_COLUMNS)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote {len(out_rows)} rows to {CSV_OUT}")
    print(f"assembly_level changes: {len(level_changes)}")
    for d, old, new in level_changes[:30]:
        print(f"  {d}: {old} -> {new}")
    if len(level_changes) > 30:
        print(f"  ... and {len(level_changes) - 30} more")
    print(f"species mismatches (old csv vs checklist): {len(species_mismatches)}")
    for item in species_mismatches:
        print(f"  {item}")
    print(f"directories missing under mitogenomes_output: {len(dir_missing_on_disk)}")


if __name__ == "__main__":
    main()
