# Beringian mitogenome reassembly

-----

| Name | Date |
| - | - |
| Keiler Collier | 9 April 2026 |

## Motivation

This is an effort to re-assemble several hundred reference-contaminated avian mitogenomes from Keiler Collier's MS thesis.
As [Sangster & Luksenburg](https://academic.oup.com/gbe/article/13/9/evab210/6368065?login=false) state, this is a major problem.

In many cases, we variant-called low-cov Illumina data from references with extremely high quality thresholds.
This resulted in genuine variants being lost in favor of calling the reference sequence.

## Directory

| Title | Link |
| - | - |
| Bad mitogenome coverage visualization prior to running new assemblies. | [Bad mdDNA analysis](notebooks/bad_mtDNA_analysis.ipynb) |
| Reassembly protocol testing | [Protocol details](notebooks/mtDNA_assembly_protocol_testing.ipynb) |
| Checklist of reassembled mitogenomes | [Completed mitogenome checklist](./mitogenomes_assembled_checklist.md) |
| Reassembled mitogenome summary statistics and ML analysis | [Assembly stats notebook](notebooks/assembly_stats.ipynb) |
| Reassembled mitogenome summary statistics (CSV) | [Assembly results CSV](./mitogenome_assembly_results.csv) |
| Executable scripts | [Scripts directory](./scripts/) |

-----

## Assembly protocol

We assembled mitogenomes directly from trimmed WGS reads using NOVOplasty as an assembler.
In an *Anser albifrons* test case, this was found to be more effective and less taxon-sensitive than splitting out mtDNA with bbsplit.

We found the taxon of the reference sequence to be irrelevant in determining assembly quality.

## Assembly matching

Raw data from UAF was not uniformly labeled and it was not apparent which species any mitogenome belonged to.

We used `table_S1.csv` (and a 2019 UAM flat file I have retained) to reunite taxon with assembly.
~150 assemblies are still missing, many of which are likely individuals involved in [Jessica McLaughlin's project](https://onlinelibrary.wiley.com/doi/10.1111/mec.15574).
These raw reads will need to be referenced from GenBank.

## Assembly stats

We assembled using all the same values.

## Annotation stats

We annotated with mitoZ.

### Data wrangling

After running NOVOplasty, we classified every output as `failed`, `partial_success` or `success`, depending on the size of the final assembly.

I then added a taxonomic order column and united the dataframe with the depth/coverage values from Collier et al. 2025.
