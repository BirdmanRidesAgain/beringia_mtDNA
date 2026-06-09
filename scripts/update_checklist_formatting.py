#!/usr/bin/env python3
"""Reformat mitogenome checklist: italicize species, split pending-reassembly table."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CHECKLISTS = [
    ROOT / "mitogenomes_assembled_checklist.md",
    ROOT / "mitogenome_assembly_checklist.md",
]

PENDING_REASSEMBLY = {
    ("UAM 20098", "JMM327"),
    ("UAM 35755", "JJW2547"),
    ("UAM 17748", "UAMX2921"),
}

TABLE2_MISSING = """## Table 3. Individuals in `table_S1.csv` without an assembled mitogenome

This **1** individual appears in `table_S1.csv` but does not have a matching row in the main assembled table (catalog + field #).

| Species | Catalog num | Field num | Year | Locality | NCBI-SRA |
| - | - | - | - | - | - |
| *Numenius phaeopus variegatus* | UAM 14228 | TB1364 | 2001 | Philippines: Cordova; Mactan Island | OR243053 |
"""


def italicize_species(species: str) -> str:
    inner = re.sub(r"^\*|\*$", "", species.strip())
    return f"*{inner}*"


def parse_table_row(line: str) -> list[str] | None:
    if not line.startswith("|") or line.startswith("| - |") or "Species | Catalog" in line:
        return None
    return [p.strip() for p in line.split("|")[1:-1]]


def format_table_row(cells: list[str]) -> str:
    return "| " + " | ".join(cells) + " |"


def transform(content: str) -> str:
    lines = content.splitlines()
    out: list[str] = []
    in_table1 = False
    pending_rows: list[str] = []
    table1_done = False

    i = 0
    while i < len(lines):
        line = lines[i]

        if line.startswith("## Table 1."):
            in_table1 = True
            out.append(line)
            i += 1
            continue

        if in_table1 and line.startswith("## Table 2."):
            in_table1 = False
            table1_done = True
            # Replace old Table 2 block with new Table 2 + Table 3
            out.append("")
            out.append("## Table 2. Pending reassembly (zero-length circularized FASTA)")
            out.append("")
            out.append(
                "These individuals have a NOVOPlasty circularized FASTA file on disk "
                "but `length_bp` is 0; reassembly is planned."
            )
            out.append("")
            out.append(
                "| Species | Catalog num | Field num | Old directory name | New directory name | multiple datasets? | In Table S1? |"
            )
            out.append("| - | - | - | - | - | - | - |")
            out.extend(pending_rows)
            out.append("")
            out.append(TABLE2_MISSING.rstrip())
            # Skip until end of file (old Table 2 content)
            while i < len(lines):
                i += 1
            continue

        if in_table1:
            cells = parse_table_row(line)
            if cells and len(cells) >= 7:
                catalog, field = cells[1], cells[2]
                cells[0] = italicize_species(cells[0])
                if (catalog, field) in PENDING_REASSEMBLY:
                    pending_rows.append(format_table_row(cells))
                    i += 1
                    continue
                out.append(format_table_row(cells))
                i += 1
                continue

        if not table1_done and line.startswith("## Table 2."):
            # Should not happen if Table 1 had rows; handled above
            pass

        out.append(line)
        i += 1

    if pending_rows and not table1_done:
        raise RuntimeError("Table 2 header not found after Table 1")

    return "\n".join(out).rstrip() + "\n"


def main() -> None:
    source = ROOT / "mitogenomes_assembled_checklist.md"
    transformed = transform(source.read_text(encoding="utf-8"))
    for path in CHECKLISTS:
        path.write_text(transformed, encoding="utf-8")
        print(f"Updated {path}")


if __name__ == "__main__":
    main()
