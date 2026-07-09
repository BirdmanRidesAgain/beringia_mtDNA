# Mitogenome assembly results - key to columns

One row per assembled mitogenome in **Table 1** of `mitogenome_assembly_checklist.md`. Metrics are read from MitoZ `summary.txt` under `mitoz/` when that pipeline produced a circular 15–17 kbp scaffold; otherwise from NOVOPlasty `summary.txt` under `novoplasty_*/annotation/` when present (local `mitogenomes_output/` or SSD `assembled_mitogenomes/mitogenomes_output/<group>/`).

| Column name | Interpretation |
| - | - |
| `directory_name` | Canonical sample folder under `mitogenomes_output/` (species, museum catalog, field #). |
| `species` | Genus, species, and subspecies if applicable (from the checklist / Table S1). |
| `assembly_success` | Boolean. `True` when a full circular mitogenome was recovered via MitoZ (circular scaffold 15–17 kbp), or when a NOVOPlasty run has an annotated `summary.txt` under `novoplasty_*/annotation/`. A circularized NOVOPlasty FASTA alone is not sufficient. `False` otherwise. |
| `assembler` | `mitoz` when `assembly_success` comes from the MitoZ pipeline; `novoplasty` when success comes from a NOVOPlasty annotation `summary.txt`; empty otherwise. |
| `length_bp` | Semicolon-separated scaffold lengths (bp) from the MitoZ scaffold table (one value per `NODE_*` or `Contig*`). For NOVOPlasty successes, from the annotation `summary.txt`. Without either, lengths of best-tier NOVOPlasty assembly file(s). |
| `circularity` | Semicolon-separated `yes`/`no` flags aligned with `length_bp` nodes from the active summary. Empty when no summary. |
| `num_protein_coding_genes_totally_found` | From the active `summary.txt` (`Protein coding genes totally found`). `0` when no summary. |
| `num_trna_coding_genes_totally_found` | From the active `summary.txt` (`tRNA genes totally found`). `0` when no summary. |
| `num_rrna_coding_genes_totally_found` | From the active `summary.txt` (`rRNA genes totally found`). `0` when no summary. |
| `num_genes_totally_found` | From the active `summary.txt` (`Genes totally found`). `0` when no summary. |
| `ND2_present` | Boolean. `True` if any gene row in the active summary has `Gene_name` `ND2`. `False` otherwise. |
| `catalog_num` | Museum catalog number (e.g. `UAM22113`, `UWBM44114`). |
| `field_num` | Field / accession label (e.g. `REW601`, `JJW905`). |
| `coverage` | mtDNA coverage from Collier et al. 2025 (`samtools depth` on earlier assemblies), when available. |
| `order` | Taxonomic order (derived from genus). |
| `taxon_pair` | Taxon pair label from Collier et al. 2025, when available. |
