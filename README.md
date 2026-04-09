# Beringian mitogenome reassembly

-----

| Name | Date |
| - | - |
| Keiler Collier | 9 April 2026 |

## Motivation

This is an effort to re-assemble several hundred reference-contaminated avian mitogenomes from Keiler Collier's MS thesis.
In many cases, we variant-called low-cov Illumina data from references with extremely high quality thresholds.
This resulted in genuine variants being lost in favor of calling the reference sequence.

For motivation and a more in-depth analysis, see my [bad mdDNA analysis](bad_mtDNA_analysis.ipynb).

See the checklist of completed mitogenomes [at this sheet](./mitogenomes_assembled_checklist.md).


## Quickstart

Install the dependencies, and run the following code:

```{bash}
FQ1=sample1_R1.fq.gz
FQ2=sample1_R2.fq.gz
SAMPLE_NAME=sample1
OUTPUT_PREFIX=test_prefix
REFSEQ=closest_mitogenome_to_sample.fa

./assemble_mitogenome.py --fq1 $FQ1 --fq2 $FQ2 --sample_name $SAMPLE_NAME --output $OUTPUT_PREFIX --refseq $REFSEQ -m 32

```

To run assembly without mtDNA splitting (trimmed reads go straight to NOVOPlasty), add `-S` / `--skip-bbsplit`.
This takes longer, but uses less memory, and is less sensitive to refseq/seed selection.

### `resources/` (bundled assets)

The repo keeps read-only pipeline inputs in **`resources/`**:

- **`resources/NOVOPlasty_config.txt`** — template copied per run into `<output>/novoplasty_<sample>/`.
- **`resources/bbsplit_ref/`** — cache of reference FASTA copies (and BBTools index files next to them). Your original `--refseq` path is unchanged; `bbsplit` uses the staged copy so indices do not clutter the directory where you keep references.

You can add `resources/bbsplit_ref/` to `.gitignore` if you do not want index caches under version control.

### `assemble_mitogenome.py` options

| Short | Long | Type | Description |
| - | - | - | - |
| `-1` | `--fq1` | path (file) | Input forward (R1) FASTQ, optionally gzip-compressed. |
| `-2` | `--fq2` | path (file) | Input reverse (R2) FASTQ, optionally gzip-compressed. |
| `-s` | `--sample_name` | string | Sample identifier; used for output subdirectories and NOVOPlasty project name. |
| `-o` | `--output` | path (directory) | Root output directory; creates `trim_<sample>`, optionally `bbsplit_<sample>`, and `novoplasty_<sample>` beneath it. |
| `-r` | `--refseq` | path (file) | Mitochondrial reference FASTA used as `bbsplit` bait and as the NOVOPlasty seed (`Seed Input`). |
| `-S` | `--skip-bbsplit` | flag (boolean) | If set, skip `bbsplit`; NOVOPlasty uses the `fastp` trimmed reads. `-m` is not used in this mode. |
| `-m` | `--memory` | integer (gigabytes) | Java heap for `bbsplit.sh` as `-Xmx<m>g`. Default `8`. Ignored when `--skip-bbsplit` is set. |
| `-t` | `--threads` | integer | Worker threads for `fastp` (`-w`) and `bbsplit` (`threads=`). Default: half of `os.cpu_count()`, minimum `1`. |
| `-v` | `--verbose` | flag (boolean) | Print each command and stream `fastp` / `bbsplit` / NOVOPlasty output. Default is quiet (tools run with captured output); on failure, their stderr/stdout is printed before the traceback. |

## Methods

### Overview

We run three tools:

- [fastp](https://github.com/OpenGene/fastp) to trim data
- [bbsplit](https://github.com/bbushnell/BBTools) to split out mtDNA
- [NOVOPlasty](https://github.com/ndierckx/NOVOPlasty) to assemble mitogenomes

### Dependencies

Dependencies are most easily installed with the [conda](https://anaconda.org/channels/anaconda/packages/conda/overview) environment described in `environment.yml`.
It can be created with

```{bash}
conda env create -f environment.yml
conda activate mitogenome-tools
```

### Refseq selection

This pipeline uses a WGS mitochondrial sequence as bait to split out the mtDNA from UCE Illumina reads.
It's important to use as close a reference as possible, because the greedy splitting algorithm that `bbsplit.sh` uses becomes less effective when its less similar to the target:

| Sample | Refseq | mtDNA split yield (R1+R2) |
| - | - | - |
| `CSW_4521` (*Anser albifrons*) | `Anser_anser_OZ124297.1_mitogenome.fa` | $4.21\text{Mb}$ |
| `CSW_4521` (*Anser albifrons*) | `Taeniopygia_guttata_NC_007897.1_mitogenome.fa` | $1.46 \text{Mb}$ |

The congeneric had over 3x more data.