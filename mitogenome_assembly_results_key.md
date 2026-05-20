# Mitogenome assembly results - key to columns

One row per assembled mitogenome in **Table 1** of `mitogenomes_assembled_checklist.md`. Assembly status is read from `mitogenomes_output/<directory_name>/novoplasty_*` outputs.

| Column name | Interpretation |
| - | - |
| `directory_name` | Canonical sample folder under `mitogenomes_output/` (species, museum catalog, field #). |
| `species` | Genus, species, and subspecies if applicable (from the checklist / Table S1). |
| `assembly_level` | Best NOVOPlasty tier under `novoplasty_*`: `success` = circularized FASTA; `partial_success` = `Option_*` or `Contigs_*` FASTA only; `failed` = only non-empty `contigs_tmp_*`, empty output, or no `novoplasty_*` dir. |
| `catalog_num` | Museum catalog number (e.g. `UAM22113`, `UWBM44114`). |
| `field_num` | Field / accession label (e.g. `REW601`, `JJW905`). |
| `coverage` | mtDNA coverage from Collier et al. 2025 (`samtools depth` on earlier assemblies), when available. |
| `order` | Taxonomic order (derived from genus). |
| `taxon_pair` | Taxon pair label from Collier et al. 2025, when available. |
