# Mitogenome assembly results - key to columns

One row per assembled mitogenome in **Table 1** of `mitogenome_assembly_checklist.md`. Metrics are read from the best MitoZ `summary.txt` under `mitoz/` (local `mitogenomes_output/` or SSD archive), falling back to NOVOPlasty for length/success only when no MitoZ summary exists.

| Column name | Interpretation |
| - | - |
| `directory_name` | Canonical sample folder under `mitogenomes_output/` (species, museum catalog, field #). |
| `species` | Genus, species, and subspecies if applicable (from the checklist / Table S1). |
| `assembly_success` | Boolean. `True` when a full circular mitogenome was recovered: MitoZ circular scaffold 15–17 kbp, or (if no MitoZ output) a NOVOPlasty circularized FASTA. `False` otherwise. |
| `length_bp` | Semicolon-separated scaffold lengths (bp) from the MitoZ scaffold table (one value per `NODE_*`). Without MitoZ, lengths of best-tier NOVOPlasty assembly file(s). |
| `circularity` | Semicolon-separated `yes`/`no` flags aligned with `length_bp` nodes from MitoZ. Empty when no MitoZ summary. |
| `num_protein_coding_genes_totally_found` | From MitoZ `Protein coding genes totally found`. `0` when no MitoZ summary (failed / not run). |
| `num_trna_coding_genes_totally_found` | From MitoZ `tRNA genes totally found`. `0` when no MitoZ summary. |
| `num_rrna_coding_genes_totally_found` | From MitoZ `rRNA genes totally found`. `0` when no MitoZ summary. |
| `num_genes_totally_found` | From MitoZ `Genes totally found`. `0` when no MitoZ summary. |
| `ND2_present` | Boolean. `True` if any gene row in the MitoZ summary has `Gene_name` `ND2`. `False` otherwise. |
| `catalog_num` | Museum catalog number (e.g. `UAM22113`, `UWBM44114`). |
| `field_num` | Field / accession label (e.g. `REW601`, `JJW905`). |
| `coverage` | mtDNA coverage from Collier et al. 2025 (`samtools depth` on earlier assemblies), when available. |
| `order` | Taxonomic order (derived from genus). |
| `taxon_pair` | Taxon pair label from Collier et al. 2025, when available. |
