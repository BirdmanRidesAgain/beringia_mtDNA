#!/usr/bin/env python3
"""Move Table 1 rows not in table_S1.tsv to Table 4."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CHECKLISTS = [
    ROOT / "mitogenomes_assembled_checklist.md",
    ROOT / "mitogenome_assembly_checklist.md",
]
S1 = ROOT / "table_S1.tsv"

TABLE_HEADER = (
    "| Species | Catalog num | Field num | Old directory name | New directory name | "
    "multiple datasets? | In Table S1? |"
)
TABLE_SEP = "| - | - | - | - | - | - | - |"

TABLE2_PENDING = """## Table 2. Pending reassembly (zero-length circularized FASTA)

These individuals have a NOVOPlasty circularized FASTA file on disk but `length_bp` is 0; reassembly is planned.

| Species | Catalog num | Field num | Old directory name | New directory name | multiple datasets? | In Table S1? |
| - | - | - | - | - | - | - |
| *Arenaria melanocephala* | UAM 20098 | JMM327 | `JMM327_Taeniopygia_guttata_NC_007897.1_mitogenome.fa/` | `Arenaria_melanocephala_UAM_20098_JMM327/` | False | True |
| *Larus argentatus smithsonianus* | UAM 35755 | JJW2547 | `JJW2547_Taeniopygia_guttata_NC_007897.1_mitogenome.fa/` | `Larus_argentatus_smithsonianus_UAM_35755_JJW2547/` | False | True |
| *Leucosticte tephrocotis umbrina* | UAM 17748 | UAMX2921 | `UAMX2921_Taeniopygia_guttata_NC_007897.1_mitogenome.fa/` | `Leucosticte_tephrocotis_umbrina_UAM_17748_UAMX2921/` | True | True |
"""

TABLE3_MISSING = """## Table 3. Individuals in `table_S1.csv` without an assembled mitogenome

This **1** individual appears in `table_S1.csv` but does not have a matching row in the main assembled table (catalog + field #).

| Species | Catalog num | Field num | Year | Locality | NCBI-SRA |
| - | - | - | - | - | - |
| *Numenius phaeopus variegatus* | UAM 14228 | TB1364 | 2001 | Philippines: Cordova; Mactan Island | OR243053 |
"""

UAM29835_DIRECTORY = "Calcarius_lapponicus_coloratus_UAM_29835"


def norm_catalog(value: str) -> str:
    return re.sub(r"\s+", "", value.strip())


def norm_field(value: str) -> str:
    return value.strip()


def italicize_species(species: str) -> str:
    inner = re.sub(r"^\*|\*$", "", species.strip())
    return f"*{inner}*"


def load_s1_keys(path: Path) -> set[tuple[str, str]]:
    keys: set[tuple[str, str]] = set()
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 4 or parts[0] == "Species":
            continue
        keys.add((norm_catalog(parts[1]), norm_field(parts[3])))
    return keys


def parse_table_row(line: str) -> list[str] | None:
    if not line.startswith("|") or line.startswith("| - |") or "Species | Catalog" in line:
        return None
    return [p.strip() for p in line.split("|")[1:-1]]


def format_table_row(cells: list[str]) -> str:
    return "| " + " | ".join(cells) + " |"


def apply_uam29835_fix(cells: list[str]) -> list[str]:
    if norm_catalog(cells[1]) != "UAM29835":
        return cells
    cells = cells.copy()
    cells[0] = italicize_species("Calcarius lapponicus coloratus")
    cells[2] = ""
    cells[4] = f"`{UAM29835_DIRECTORY}`"
    cells[6] = "False"
    return cells


def transform(content: str, s1_keys: set[tuple[str, str]]) -> str:
    lines = content.splitlines()
    out: list[str] = []
    table1_rows: list[str] = []
    table4_rows: list[str] = []
    in_table1 = False

    for line in lines:
        if line.startswith("## Table 1."):
            in_table1 = True
            out.append(line)
            continue
        if in_table1 and line.startswith("## Table"):
            in_table1 = False
            break
        if in_table1:
            cells = parse_table_row(line)
            if cells and len(cells) >= 7:
                cells[0] = italicize_species(cells[0])
                key = (norm_catalog(cells[1]), norm_field(cells[2]))
                if key in s1_keys:
                    cells[6] = "True"
                    table1_rows.append(format_table_row(cells))
                else:
                    cells = apply_uam29835_fix(cells)
                    cells[6] = "False"
                    table4_rows.append(format_table_row(cells))
                continue
            out.append(line)
            continue
        out.append(line)

    table4_rows.sort(key=lambda row: row.lower())
    table4_block = [
        "",
        "## Table 4. Individuals with assembled genomes not present in `table_S1.csv`",
        "",
        f"This **{len(table4_rows)}** individual"
        + ("s have" if len(table4_rows) != 1 else " has")
        + " assembled mitogenomes but do not appear in `table_S1.csv` (catalog + field #).",
        "",
        TABLE_HEADER,
        TABLE_SEP,
        *table4_rows,
    ]

    out.extend(table1_rows)
    out.append("")
    out.append(TABLE2_PENDING.rstrip())
    out.append("")
    out.append(TABLE3_MISSING.rstrip())
    out.extend(table4_block)
    return "\n".join(out).rstrip() + "\n"


def main() -> None:
    s1_keys = load_s1_keys(S1)
    source = ROOT / "mitogenomes_assembled_checklist.md"
    transformed = transform(source.read_text(encoding="utf-8"), s1_keys)
    for path in CHECKLISTS:
        path.write_text(transformed, encoding="utf-8")
        print(f"Updated {path}")


if __name__ == "__main__":
    main()
