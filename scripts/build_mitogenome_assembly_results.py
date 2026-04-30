#!/usr/bin/env python3
"""Build canonical mitogenome assembly results table.

This script produces a union table that retains rows from both:
- assembly status table (new_mitogenome_assembly_stats.csv)
- historical coverage table (old_mitogenome_depths.csv)
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent
ASSEMBLY_STATS = ROOT / "new_mitogenome_assembly_stats.csv"
ASSEMBLY_STATS_FALLBACK = ROOT / "mitogenome_assembly_results.csv"
DEPTHS = ROOT / "old_mitogenome_depths.csv"
OUTPUT = ROOT / "mitogenome_assembly_results.csv"


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
    "Calidria": "Charadriiformes",
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


def parse_directory_keys(directory_name: str) -> tuple[str | None, str | None]:
    """Extract catalog/accession guesses from directory names."""
    match_standard = re.search(r"_(UAM|UWBM)_(\d+)_([^_]+)$", directory_name)
    if match_standard:
        museum, number, accession = match_standard.groups()
        return f"{museum}{number}", accession

    match_fmnh = re.search(r"_(FMNH\d+)_([^_]+)$", directory_name)
    if match_fmnh:
        catalog, accession = match_fmnh.groups()
        return catalog, accession

    return None, None


def normalize_old_columns(df_old: pd.DataFrame) -> pd.DataFrame:
    colmap = {}
    for col in df_old.columns:
        lowered = col.strip().lower()
        if lowered in {"mtdna accession", "accession"}:
            colmap[col] = "accession"
        elif lowered == "species":
            colmap[col] = "species_old"
        elif lowered == "catalog":
            colmap[col] = "catalog"
        elif lowered == "coverage":
            colmap[col] = "coverage"
        elif lowered == "taxon pair":
            colmap[col] = "taxon_pair"

    df = df_old.rename(columns=colmap).copy()
    required = {"species_old", "catalog", "accession", "coverage", "taxon_pair"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing expected columns in old depths CSV: {sorted(missing)}")

    df["catalog"] = df["catalog"].astype(str).str.strip()
    df["accession"] = df["accession"].astype(str).str.strip()
    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    return df


def add_order(series_species: pd.Series) -> pd.Series:
    genus = series_species.fillna("").astype(str).str.split().str[0]
    return genus.map(GENUS_TO_ORDER)


def sanitize_species_for_dir(value: str) -> str:
    value = re.sub(r"\s+", "_", value.strip())
    return re.sub(r"[^A-Za-z0-9_]", "", value)


def load_assembly_base() -> pd.DataFrame:
    """Load base assembly table, with fallback to existing merged file."""
    if ASSEMBLY_STATS.exists():
        df = pd.read_csv(ASSEMBLY_STATS)
    elif ASSEMBLY_STATS_FALLBACK.exists():
        df = pd.read_csv(ASSEMBLY_STATS_FALLBACK)
    else:
        raise FileNotFoundError(
            "Could not find new_mitogenome_assembly_stats.csv or mitogenome_assembly_results.csv"
        )

    required = {"directory_name", "species", "assembly_level"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Assembly base missing columns: {sorted(missing)}")
    return df[list(required)].copy()


def main() -> None:
    df_new = load_assembly_base()
    df_old = normalize_old_columns(pd.read_csv(DEPTHS))

    keys = df_new["directory_name"].apply(parse_directory_keys)
    df_new["catalog_guess"] = keys.str[0]
    df_new["accession_guess"] = keys.str[1]

    merged = df_new.merge(
        df_old,
        how="outer",
        left_on=["catalog_guess", "accession_guess"],
        right_on=["catalog", "accession"],
        indicator=True,
    )

    # Fill key fields for old-only rows.
    merged["catalog_guess"] = merged["catalog_guess"].where(
        merged["catalog_guess"].notna(), merged["catalog"]
    )
    merged["accession_guess"] = merged["accession_guess"].where(
        merged["accession_guess"].notna(), merged["accession"]
    )

    merged["species"] = merged["species"].where(merged["species"].notna(), merged["species_old"])
    merged["directory_name"] = merged["directory_name"].where(
        merged["directory_name"].notna(),
        merged.apply(
            lambda r: (
                f"{sanitize_species_for_dir(str(r['species_old']))}_{r['catalog']}_{r['accession']}"
                if pd.notna(r["species_old"]) and pd.notna(r["catalog"]) and pd.notna(r["accession"])
                else pd.NA
            ),
            axis=1,
        ),
    )

    merged["match_status"] = merged["_merge"].map(
        {"both": "matched", "left_only": "unmatched", "right_only": "unmatched"}
    )
    merged["match_confidence"] = merged["match_status"].map(
        {"matched": "high", "unmatched": "low"}
    )

    def notes(row: pd.Series) -> str:
        if row["_merge"] == "both":
            return ""
        if row["_merge"] == "right_only":
            return "no_matching_assembly_status_in_new"
        if pd.isna(row["catalog_guess"]) or pd.isna(row["accession_guess"]):
            return "control_or_unparsed_directory_name"
        return "no_matching_catalog_accession_in_old_depths"

    merged["match_notes"] = merged.apply(notes, axis=1)
    merged["order"] = add_order(merged["species"])

    # Keep old species label as audit field.
    merged["species_old"] = merged["species_old"].where(
        merged["species_old"].notna(), pd.NA
    )

    out_cols = [
        "directory_name",
        "species",
        "assembly_level",
        "catalog_guess",
        "accession_guess",
        "coverage",
        "order",
        "taxon_pair",
        "species_old",
        "match_status",
        "match_confidence",
        "match_notes",
        "_merge",
    ]
    merged[out_cols].to_csv(OUTPUT, index=False)
    print(f"Wrote {OUTPUT} ({len(merged)} rows)")


if __name__ == "__main__":
    main()
