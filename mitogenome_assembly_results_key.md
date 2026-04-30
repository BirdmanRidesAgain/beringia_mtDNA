# Mitogenome assembly results - key to columns

This dataset is made from the union of a pair of csvs - one documenting depth from the 'old' mitogenomes, and the other documenting the status of the NOVOplasty (reference-free) mitogenomes.

| Column name | Interpretation |
| - | - |
| `directory_name` | Canonical name of the directory the new mitogenome is found in. Composed of taxonomic identity, UAM number and field number. |
| `species` | Genus, species, and subspecies if applicable. Derived from Table_S1. |
| `assembly_level` | NOVOPlasty assembly level. |
| `catalog_guess` | UAM number, derived from `directory_name` |
| `accession_guess` | Field number, derived from `directory_name` |
| `coverage` | Coverage estimates. Taken from Collier et al. 2025, which was Calculated from `samtools depth` estimates of old mitogenomes. |
| `order` | Taxonomic order |
| `taxon_pair` | Taxon pair from Collier et al. 2025. |
| `species_old` | Artifact from union where sheet was derived from. |
| `match_status` | Whether or not there was a match between the two data sheets. |
| `match_confidence` | |
| `match_notes` | |
| `_merge` | Whether this entry was found in both sheets, the 'old' mitogenomes only (`right_only`), or the new mitogenomes only (`left_only`). `right_only` entries usually correspond to McLaughlin's dataset. `left_only` ones usually correspond to UWBM birds, or individuals not included in Collier et al. 2025 the first time around. |
