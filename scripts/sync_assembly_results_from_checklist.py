#!/usr/bin/env python3
"""Rebuild mitogenome_assembly_results.csv from mitogenome_assembly_checklist.md Table 1."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from mitogenome_paths import (  # noqa: E402
    assembly_metrics_for_sample_dir,
    resolve_sample_dir,
)

ROOT = Path(__file__).resolve().parent.parent
CHECKLIST = ROOT / "mitogenome_assembly_checklist.md"
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
    "UAM29835_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Calcarius_lapponicus_coloratus_UAM_29835",
    "Calcarius_lapponicus_alascensis_UAM_29835": "Calcarius_lapponicus_coloratus_UAM_29835",
    "DDG1922_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Larus_brachyrhynchus_UAM_14803_DDG1922",
    "Larus_brachyrhynchus_UAM_14803_DDG2122": "Larus_brachyrhynchus_UAM_14803_DDG1922",
    "Lagopus_lagopus_koreni_UAM_28171_KSW5133": "Lagopus_lagopus_koreni_UAM_28171_KSW5113",
    "Melanitta_americana_UAM_35578_JJW951": "Melanitta_americana_UAM_28065_JJW951",
    "Melanitta_americana_UAM_28065_JJW950": "Melanitta_americana_UAM_27720_JJW950",
    "UAM10581_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Picoides_tridactylus_UAM_10581",
    "TransB_SVD_422_ACTCCATC-CTCTGGTT__Taeniopygia_guttata_NC_007897.1_mitogenome.fa": (
        "Pinicola_enucleator_kamtschatschensis_UWBM_51642_SVD422"
    ),
    "JJW2478_Taeniopygia_guttata_NC_007897.1_mitogenome.fa": "Spatula_clypeata_UAM_35604_JJW2478",
    "Mareca_penelope_UAM9759_DDG1703": "Mareca_penelope_UAM_9359_DDG1703",
}

OUT_COLUMNS = [
    "directory_name",
    "species",
    "assembly_success",
    "assembler",
    "length_bp",
    "circularity",
    "num_protein_coding_genes_totally_found",
    "num_trna_coding_genes_totally_found",
    "num_rrna_coding_genes_totally_found",
    "num_genes_totally_found",
    "ND2_present",
    "catalog_num",
    "field_num",
    "coverage",
    "order",
    "taxon_pair",
    "data_origin",
]

# Samples from McLaughlin et al. 2020 (vs unknown origin elsewhere).
MCLAUGHLIN_GENERA = {"Luscinia", "Clangula", "Pinicola", "Pica", "Mareca", "Numenius"}


def data_origin_for_species(species: str) -> str:
    species = (species or "").strip()
    genus = species.split()[0] if species else ""
    if genus in MCLAUGHLIN_GENERA:
        return "mclaughlin_2020"
    if species.startswith("Tringa brevipes") or species.startswith("Tringa incana"):
        return "mclaughlin_2020"
    return "unknown"


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


def parse_legacy_assembly_success(row: dict) -> str | None:
    """Normalize prior success labels from old or new CSV column names."""
    if "assembly_success" in row and row["assembly_success"] != "":
        text = str(row["assembly_success"]).strip().lower()
        if text in ("complete", "partial", "failure"):
            return text
        if text in ("true", "1", "yes"):
            return "complete"
        if text in ("false", "0", "no"):
            return "failure"
        return text
    level = (row.get("assembly_level") or "").strip()
    if level == "success":
        return "complete"
    if level == "partial_success":
        return "partial"
    if level == "failed":
        return "failure"
    return None


def main() -> None:
    checklist = parse_checklist(CHECKLIST)
    by_key, by_dir = load_csv_index(CSV_IN)

    out_rows = []
    success_changes = []
    species_mismatches = []
    dir_missing_on_disk = []

    for chk in checklist:
        key = match_key(chk["catalog_num"], chk["field_num"])
        old = by_dir.get(chk["directory_name"]) or by_key.get(key)
        if old:
            coverage = old.get("coverage") or ""
            taxon_pair = old.get("taxon_pair") or ""
            old_success = parse_legacy_assembly_success(old)
            old_species = old.get("species") or ""
        else:
            coverage = taxon_pair = old_species = ""
            old_success = None

        disk_path = resolve_sample_dir(chk["directory_name"])
        metrics = assembly_metrics_for_sample_dir(disk_path)
        if disk_path is None:
            dir_missing_on_disk.append(chk["directory_name"])

        if old_success is not None and old_success != metrics["assembly_success"]:
            success_changes.append(
                (chk["directory_name"], old_success, metrics["assembly_success"])
            )

        if old_species and old_species != chk["species"] and old_species != "Unknown":
            species_mismatches.append((chk["directory_name"], old_species, chk["species"]))

        out_rows.append(
            {
                "directory_name": chk["directory_name"],
                "species": chk["species"],
                "assembly_success": metrics["assembly_success"],
                "assembler": metrics["assembler"],
                "length_bp": metrics["length_bp"],
                "circularity": metrics["circularity"],
                "num_protein_coding_genes_totally_found": metrics[
                    "num_protein_coding_genes_totally_found"
                ],
                "num_trna_coding_genes_totally_found": metrics[
                    "num_trna_coding_genes_totally_found"
                ],
                "num_rrna_coding_genes_totally_found": metrics[
                    "num_rrna_coding_genes_totally_found"
                ],
                "num_genes_totally_found": metrics["num_genes_totally_found"],
                "ND2_present": metrics["ND2_present"],
                "catalog_num": chk["catalog_num"],
                "field_num": chk["field_num"],
                "coverage": coverage,
                "order": add_order(chk["species"]),
                "taxon_pair": taxon_pair,
                "data_origin": data_origin_for_species(chk["species"]),
            }
        )

    with CSV_OUT.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=OUT_COLUMNS)
        writer.writeheader()
        writer.writerows(out_rows)

    n_complete = sum(1 for r in out_rows if r["assembly_success"] == "complete")
    n_partial = sum(1 for r in out_rows if r["assembly_success"] == "partial")
    n_failure = sum(1 for r in out_rows if r["assembly_success"] == "failure")
    n_mitoz = sum(1 for r in out_rows if r["assembler"] == "mitoz")
    n_novoplasty = sum(1 for r in out_rows if r["assembler"] == "novoplasty")
    n_nd2 = sum(1 for r in out_rows if r["ND2_present"])
    n_with_genes = sum(1 for r in out_rows if r["num_genes_totally_found"] > 0)
    print(f"Wrote {len(out_rows)} rows to {CSV_OUT}")
    print(
        f"assembly_success complete: {n_complete}, partial: {n_partial}, "
        f"failure: {n_failure}"
    )
    print(
        f"assembler mitoz: {n_mitoz}, novoplasty: {n_novoplasty}, "
        f"unset: {len(out_rows) - n_mitoz - n_novoplasty}"
    )
    print(f"rows with gene counts: {n_with_genes}")
    print(f"ND2_present True: {n_nd2}")
    print(f"assembly_success changes: {len(success_changes)}")
    for d, old, new in success_changes[:30]:
        print(f"  {d}: {old} -> {new}")
    if len(success_changes) > 30:
        print(f"  ... and {len(success_changes) - 30} more")
    print(f"species mismatches (old csv vs checklist): {len(species_mismatches)}")
    for item in species_mismatches:
        print(f"  {item}")
    print(f"directories missing under mitogenomes_output: {len(dir_missing_on_disk)}")


if __name__ == "__main__":
    main()
