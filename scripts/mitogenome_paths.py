"""Path conventions for assemblies and co-located MitoZ annotations."""
from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
ASM_ROOT = REPO / "mitogenomes_output"
ANNOTATION_FOLDER = "annotation"
BATCH_LOGS_DIR = ASM_ROOT / "_batch_logs"
# Prefer the current SSD layout; fall back to the older keiler_C path if present.
_SSD_CANDIDATES = (
    Path("/Volumes/KC_SSD/mtDNA_reassembly/assembled_mitogenomes"),
    Path("/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/assembled_mitogenomes"),
)
ARCHIVE_ROOT = next((p for p in _SSD_CANDIDATES if p.is_dir()), _SSD_CANDIDATES[0])
SSD_OUTPUT_ROOT = ARCHIVE_ROOT / "mitogenomes_output"

# NOVOPlasty output tiers (lower = better). Matches batch_annotate_assemblies.py.
TIER_CIRC = 1
TIER_OPTION = 2
TIER_CONTIGS = 3
TIER_TMP = 4
TIER_NONE = 99
TIER_NAMES = {1: "circularized", 2: "option", 3: "contigs", 4: "contigs_tmp", 99: "none"}

# MitoZ ``summary.txt``: full circular mitogenome expected in this size range (bp).
MITOZ_SUCCESS_MIN_BP = 15_000
MITOZ_SUCCESS_MAX_BP = 17_000


@dataclass(frozen=True)
class MitozSummary:
    """Parsed MitoZ ``summary.txt`` scaffold table and gene totals."""

    scaffolds: tuple[tuple[str, int, str], ...]  # (seq_id, length_bp, circularity)
    protein_coding_genes: int
    trna_genes: int
    rrna_genes: int
    genes_total: int
    nd2_present: bool

    @property
    def length_bp(self) -> str:
        if not self.scaffolds:
            return "0"
        return ";".join(str(length) for _, length, _ in self.scaffolds)

    @property
    def circularity(self) -> str:
        if not self.scaffolds:
            return ""
        return ";".join(circ for _, _, circ in self.scaffolds)

    @property
    def assembly_success(self) -> bool:
        return any(
            circ == "yes" and MITOZ_SUCCESS_MIN_BP <= length <= MITOZ_SUCCESS_MAX_BP
            for _, length, circ in self.scaffolds
        )


def _list_novoplasty_files(novop: Path) -> list[Path]:
    return [p for p in novop.iterdir() if p.is_file() and not p.name.startswith("._")]


def classify_novoplasty_dir(novop: Path) -> int:
    """Best NOVOPlasty tier in one ``novoplasty_*`` directory."""
    files = _list_novoplasty_files(novop)
    if any(p.name.startswith("Circularized_assembly_") and p.suffix == ".fasta" for p in files):
        return TIER_CIRC
    if any(p.name.startswith("Option_") and p.suffix == ".fasta" for p in files):
        return TIER_OPTION
    if any(p.name.startswith("Contigs_") and p.suffix == ".fasta" for p in files):
        return TIER_CONTIGS
    if any(
        p.name.startswith("contigs_tmp_") and p.suffix == ".txt" and p.stat().st_size > 0
        for p in files
    ):
        return TIER_TMP
    return TIER_NONE


def tier_to_assembly_success(tier: int) -> bool:
    """True only for a circularized NOVOPlasty FASTA (tier 1)."""
    return tier == TIER_CIRC


def tier_to_assembly_level(tier: int) -> str:
    """
    Map best tier to a legacy three-way label (internal use).

    - ``success``: circularized FASTA (tier 1)
    - ``partial_success``: Option or Contigs FASTA (tiers 2–3)
    - ``failed``: only contigs_tmp, empty output, or no novoplasty dir (tiers 4, 99)
    """
    if tier == TIER_CIRC:
        return "success"
    if tier in (TIER_OPTION, TIER_CONTIGS):
        return "partial_success"
    return "failed"


_DIR_SUFFIX_RE = re.compile(r"_(UAM|UWBM)_(\d+)_([^_]+)$")

# Checklist directory_name -> on-disk mitogenomes_output folder (when names differ)
DISK_DIRECTORY_ALIASES: dict[str, str] = {
    "Calcarius_lapponicus_alascensis_UAM_29835": "Calcarius_lapponicus_coloratus_UAM_29835",
    "Larus_brachyrhynchus_UAM_14803_DDG2122": "Larus_brachyrhynchus_UAM_14803_DDG1922",
    "Lagopus_lagopus_koreni_UAM_28171_KSW5133": "Lagopus_lagopus_koreni_UAM_28171_KSW5113",
    "Melanitta_americana_UAM_35578_JJW951": "Melanitta_americana_UAM_28065_JJW951",
    "Melanitta_americana_UAM_28065_JJW950": "Melanitta_americana_UAM_27720_JJW950",
}


def _nested_ssd_sample_dirs(directory_name: str) -> list[Path]:
    """``mitogenomes_output/<species_group>/<directory_name>/`` on SSD."""
    if not SSD_OUTPUT_ROOT.is_dir():
        return []
    return sorted(
        p
        for p in SSD_OUTPUT_ROOT.glob(f"*/{directory_name}")
        if p.is_dir() and not p.name.startswith("._")
    )


def resolve_sample_dir(directory_name: str) -> Path | None:
    """
    Return an existing sample dir for a checklist/CSV directory_name.

    Prefers local ``mitogenomes_output/``, then SSD ``mitogenomes_output/<group>/``,
    then the legacy flat SSD archive under ``assembled_mitogenomes/<sample>/``.
    """
    candidates: list[Path] = []
    if directory_name in DISK_DIRECTORY_ALIASES:
        alias_name = DISK_DIRECTORY_ALIASES[directory_name]
        candidates.extend(
            [
                ASM_ROOT / alias_name,
                ARCHIVE_ROOT / alias_name,
                *_nested_ssd_sample_dirs(alias_name),
            ]
        )
    candidates.extend(
        [
            ASM_ROOT / directory_name,
            ARCHIVE_ROOT / directory_name,
            *_nested_ssd_sample_dirs(directory_name),
        ]
    )
    for path in candidates:
        if path.is_dir():
            return path

    match = _DIR_SUFFIX_RE.search(directory_name)
    if match:
        suffix = f"_{match.group(1)}_{match.group(2)}_{match.group(3)}"
        for root in (ASM_ROOT, ARCHIVE_ROOT, SSD_OUTPUT_ROOT):
            if not root.is_dir():
                continue
            if root == SSD_OUTPUT_ROOT:
                matches = [
                    p
                    for p in root.glob(f"*/{directory_name}")
                    if p.is_dir()
                ]
                if not matches:
                    matches = [
                        p
                        for group in root.iterdir()
                        if group.is_dir() and not group.name.startswith("._")
                        for p in group.iterdir()
                        if p.is_dir() and p.name.endswith(suffix)
                    ]
            else:
                matches = [p for p in root.iterdir() if p.is_dir() and p.name.endswith(suffix)]
            if len(matches) == 1:
                return matches[0]
    return None


def sequence_length_from_fasta(path: Path) -> int:
    """Total bp across all sequences in a FASTA or FASTA-like text file."""
    total = 0
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            total += len(line)
    return total


def candidate_assembly_paths(novop: Path) -> list[Path]:
    """
    Assembly outputs to measure for one ``novoplasty_*`` directory.

    Prefer Circularized > Option > Contigs FASTAs; use non-empty ``contigs_tmp_*``
    only when no FASTA exists.
    """
    files = _list_novoplasty_files(novop)
    circ = sorted(
        p for p in files if p.name.startswith("Circularized_assembly_") and p.suffix == ".fasta"
    )
    if circ:
        return circ
    options = sorted(p for p in files if p.name.startswith("Option_") and p.suffix == ".fasta")
    if options:
        return options
    contigs = sorted(p for p in files if p.name.startswith("Contigs_") and p.suffix == ".fasta")
    if contigs:
        return contigs
    contigs_tmp = sorted(
        p
        for p in files
        if p.name.startswith("contigs_tmp_") and p.suffix == ".txt" and p.stat().st_size > 0
    )
    return contigs_tmp


def assembly_length_bp_for_sample_dir(sample_dir: Path | None) -> str:
    """
    Colon-separated sequence lengths (bp) for the best-tier assembly file(s).

    When multiple ``novoplasty_*`` subdirs tie at the best tier, lengths from all
    kept subdirs are concatenated (sorted by path).
    """
    if sample_dir is None or not sample_dir.is_dir():
        return "0"

    novops = sorted(
        p for p in sample_dir.iterdir() if p.is_dir() and p.name.startswith("novoplasty_")
    )
    if not novops:
        return "0"

    classified: list[tuple[int, list[Path]]] = []
    for novop in novops:
        tier = classify_novoplasty_dir(novop)
        if tier == TIER_NONE:
            continue
        paths = candidate_assembly_paths(novop)
        if paths:
            classified.append((tier, paths))

    if not classified:
        return "0"

    best_tier = min(tier for tier, _ in classified)
    kept_paths = sorted(p for tier, paths in classified if tier == best_tier for p in paths)
    lengths = [str(sequence_length_from_fasta(p)) for p in kept_paths]
    return ";".join(lengths)


def has_mitoz_pipeline_output(sample_dir: Path | None) -> bool:
    """True when a sample has MitoZ output from trim+assemble (``sample/mitoz/``)."""
    if sample_dir is not None and (sample_dir / "mitoz").is_dir():
        return True
    if sample_dir is not None:
        archived = ARCHIVE_ROOT / sample_dir.name / "mitoz"
        if archived.is_dir():
            return True
    return False


def mitoz_dir_for_sample(sample_dir: Path) -> Path | None:
    """Local or SSD ``mitoz/`` directory for a resolved sample path."""
    local = sample_dir / "mitoz"
    if local.is_dir():
        return local
    archived = ARCHIVE_ROOT / sample_dir.name / "mitoz"
    if archived.is_dir():
        return archived
    return None


def assembly_method_for_sample_dir(sample_dir: Path | None) -> str:
    """``mitoZ`` for pipeline output under ``sample/mitoz/``; otherwise ``novoplasty``."""
    if has_mitoz_pipeline_output(sample_dir):
        return "mitoZ"
    return "novoplasty"


def _mitoz_summary_score(path: Path) -> int:
    rel = str(path)
    score = 0
    if "tmp_" in rel:
        score -= 10
    if ".result/" in rel:
        score += 10
    if "metaspades.mitogenome.fa.result" in rel:
        score += 5
    return score


def find_mitoz_summary_txt(sample_dir: Path) -> Path | None:
    """Best ``summary.txt`` under local or SSD-archived ``mitoz/``."""
    mitoz_dir = mitoz_dir_for_sample(sample_dir)
    if mitoz_dir is None:
        return None
    candidates = [
        p for p in mitoz_dir.rglob("summary.txt") if p.is_file() and not p.name.startswith("._")
    ]
    if not candidates:
        return None
    return max(candidates, key=_mitoz_summary_score)


def _parse_int_after_colon(line: str) -> int:
    return int(line.split(":", 1)[1].strip())


def parse_mitoz_summary(path: Path) -> MitozSummary:
    """Parse scaffold table, gene totals, and ND2 presence from a MitoZ ``summary.txt``."""
    scaffolds: list[tuple[str, int, str]] = []
    protein_coding_genes = 0
    trna_genes = 0
    rrna_genes = 0
    genes_total = 0
    nd2_present = False
    in_scaffold_table = False
    in_gene_table = False

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if (
            line.startswith("#Seq_id")
            and "Length(bp)" in line
            and "Circularity" in line
        ):
            in_scaffold_table = True
            in_gene_table = False
            continue
        if (
            line.startswith("#Seq_id")
            and "Gene_name" in line
        ):
            in_scaffold_table = False
            in_gene_table = True
            continue
        if in_scaffold_table:
            if not line.strip():
                in_scaffold_table = False
                continue
            if line.startswith("#") or line.startswith("-"):
                in_scaffold_table = False
                continue
            parts = line.split()
            if len(parts) >= 3 and parts[1].isdigit():
                scaffolds.append((parts[0], int(parts[1]), parts[2].lower()))
            continue
        if in_gene_table:
            if not line.strip() or line.startswith("-"):
                continue
            if line.startswith("#") or line.startswith("Protein coding"):
                in_gene_table = False
                # fall through to gene-total parsers below
            else:
                parts = line.split()
                # Columns: Seq_id Start End Length Direction Type Gene_name ...
                if len(parts) >= 7 and parts[5] in ("CDS", "tRNA", "rRNA"):
                    gene_name = parts[6]
                    if gene_name.upper() == "ND2":
                        nd2_present = True
                continue
        if line.startswith("Protein coding genes totally found:"):
            protein_coding_genes = _parse_int_after_colon(line)
        elif line.startswith("tRNA genes totally found:"):
            trna_genes = _parse_int_after_colon(line)
        elif line.startswith("rRNA genes totally found:"):
            rrna_genes = _parse_int_after_colon(line)
        elif line.startswith("Genes totally found:"):
            genes_total = _parse_int_after_colon(line)

    return MitozSummary(
        scaffolds=tuple(scaffolds),
        protein_coding_genes=protein_coding_genes,
        trna_genes=trna_genes,
        rrna_genes=rrna_genes,
        genes_total=genes_total,
        nd2_present=nd2_present,
    )


def mitoz_summary_for_sample_dir(sample_dir: Path | None) -> MitozSummary | None:
    """Best MitoZ ``summary.txt`` for a sample, or ``None`` if missing."""
    if sample_dir is None or not sample_dir.is_dir():
        return None
    summary_path = find_mitoz_summary_txt(sample_dir)
    if summary_path is None:
        return None
    return parse_mitoz_summary(summary_path)


GENE_FASTA_METASPADES_SUFFIX = "metaspades.mitogenome.fa_mitoscaf.fa.gbf.gene.fasta"


def find_annotation_summary_txt(sample_dir: Path, assembler: str = "") -> Path | None:
    """Best ``summary.txt`` for the active assembler, or either source when unknown."""
    if assembler == "mitoz":
        return find_mitoz_summary_txt(sample_dir)
    if assembler == "novoplasty":
        return find_novoplasty_annotation_summary_txt(sample_dir)
    mitoz = find_mitoz_summary_txt(sample_dir)
    if mitoz is not None:
        return mitoz
    return find_novoplasty_annotation_summary_txt(sample_dir)


def find_gene_fasta_for_sample(sample_dir: Path, assembler: str = "") -> Path | None:
    """
    Best MitoZ ``*.gbf.gene.fasta`` for a sample.

    Prefer the file co-located with the active ``summary.txt``; otherwise fall back
    to the best ``*metaspades.mitogenome.fa_mitoscaf.fa.gbf.gene.fasta`` under the
    sample directory.
    """
    summary = find_annotation_summary_txt(sample_dir, assembler)
    if summary is not None:
        siblings = [
            p
            for p in summary.parent.glob("*.gbf.gene.fasta")
            if p.is_file() and not p.name.startswith("._")
        ]
        if siblings:
            return max(siblings, key=_mitoz_summary_score)

    candidates = [
        p
        for p in sample_dir.rglob("*")
        if p.is_file()
        and p.name.endswith(GENE_FASTA_METASPADES_SUFFIX)
        and not p.name.startswith("._")
    ]
    if not candidates:
        return None
    return max(candidates, key=_mitoz_summary_score)


def parse_fasta_records(path: Path) -> list[tuple[str, str]]:
    """Return ``(header, sequence)`` pairs from a FASTA file."""
    records: list[tuple[str, str]] = []
    header: str | None = None
    seq_parts: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_parts)))
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line)
        if header is not None:
            records.append((header, "".join(seq_parts)))
    return records


def is_nd2_gene_header(header: str) -> bool:
    """True when a MitoZ gene FASTA header names ND2 (e.g. ``NODE_2;ND2;...``)."""
    return any(part.upper() == "ND2" for part in header.split(";"))


def extract_nd2_records_from_gene_fasta(path: Path) -> list[tuple[str, str]]:
    """Extract ND2 CDS records from a MitoZ ``*.gbf.gene.fasta`` file."""
    return [(h, s) for h, s in parse_fasta_records(path) if is_nd2_gene_header(h)]


def find_novoplasty_annotation_summary_txt(sample_dir: Path) -> Path | None:
    """Best ``summary.txt`` under any ``novoplasty_*/annotation/`` tree."""
    candidates: list[Path] = []
    for novop in sorted(sample_dir.iterdir()):
        if not novop.is_dir() or not novop.name.startswith("novoplasty"):
            continue
        ann = novop / ANNOTATION_FOLDER
        if not ann.is_dir():
            continue
        for summary in ann.rglob("summary.txt"):
            if summary.is_file() and not summary.name.startswith("._"):
                candidates.append(summary)
    if not candidates:
        return None
    return max(candidates, key=_mitoz_summary_score)


def novoplasty_annotation_summary_for_sample_dir(
    sample_dir: Path | None,
) -> MitozSummary | None:
    """Parsed annotation ``summary.txt`` from a NOVOPlasty run, if present."""
    if sample_dir is None or not sample_dir.is_dir():
        return None
    summary_path = find_novoplasty_annotation_summary_txt(sample_dir)
    if summary_path is None:
        return None
    return parse_mitoz_summary(summary_path)


def success_from_gene_count(num_genes: int) -> str:
    """Map gene totals to ``failure`` / ``partial`` / ``complete``."""
    if num_genes <= 0:
        return "failure"
    if num_genes <= 36:
        return "partial"
    return "complete"


def _metrics_from_summary(summary: MitozSummary, assembler: str) -> dict:
    genes = summary.genes_total
    return {
        "assembly_success": success_from_gene_count(genes),
        "length_bp": summary.length_bp,
        "circularity": summary.circularity,
        "num_protein_coding_genes_totally_found": summary.protein_coding_genes,
        "num_trna_coding_genes_totally_found": summary.trna_genes,
        "num_rrna_coding_genes_totally_found": summary.rrna_genes,
        "num_genes_totally_found": genes,
        "ND2_present": summary.nd2_present,
        "assembler": assembler if genes > 0 else "",
    }


def mitoz_results_for_sample_dir(sample_dir: Path | None) -> tuple[bool, str]:
    """``(assembly_success, length_bp)`` from MitoZ ``summary.txt``."""
    summary = mitoz_summary_for_sample_dir(sample_dir)
    if summary is None:
        return False, "0"
    return summary.assembly_success, summary.length_bp


def assembly_level_for_sample_dir(sample_dir: Path) -> str:
    """Best tier across all ``novoplasty_*`` subdirectories under a sample folder."""
    if not sample_dir.is_dir():
        return "failed"
    novops = sorted(
        p for p in sample_dir.iterdir() if p.is_dir() and p.name.startswith("novoplasty_")
    )
    if not novops:
        return "failed"
    best = TIER_NONE
    for novop in novops:
        best = min(best, classify_novoplasty_dir(novop))
    return tier_to_assembly_level(best)


def _preferred_metrics(
    mitoz_summary: MitozSummary | None,
    novo_summary: MitozSummary | None,
) -> dict | None:
    """
    Default preference when not overriding by gene count.

    Circular MitoZ (15–17 kbp) > NOVOPlasty annotation summary > incomplete MitoZ.
    """
    if mitoz_summary is not None and mitoz_summary.assembly_success:
        return _metrics_from_summary(mitoz_summary, "mitoz")
    if novo_summary is not None:
        return _metrics_from_summary(novo_summary, "novoplasty")
    if mitoz_summary is not None:
        return _metrics_from_summary(mitoz_summary, "mitoz")
    return None


def assembly_success_and_length_for_sample_dir(
    sample_dir: Path | None,
) -> tuple[str, str]:
    """
    ``(assembly_success, length_bp)`` preferring MitoZ pipeline output when present.

    ``assembly_success`` is ``failure`` / ``partial`` / ``complete`` from gene counts.
    """
    metrics = assembly_metrics_for_sample_dir(sample_dir)
    return metrics["assembly_success"], metrics["length_bp"]


def assembly_metrics_for_sample_dir(sample_dir: Path | None) -> dict:
    """
    Assembly and annotation metrics for the results CSV.

    When both MitoZ and NOVOPlasty annotation summaries exist:
    - if exactly one has 37 genes, use that assembly;
    - else if the higher gene count is ≤37, use the higher-gene assembly;
    - else keep the default preference (circular MitoZ > NOVOPlasty > incomplete MitoZ).

    ``assembler`` is set whenever at least one gene is found. ``assembly_success``
    is ``failure`` (0 genes), ``partial`` (1–36), or ``complete`` (37+).
    """
    empty = {
        "assembly_success": "failure",
        "length_bp": "0",
        "circularity": "",
        "num_protein_coding_genes_totally_found": 0,
        "num_trna_coding_genes_totally_found": 0,
        "num_rrna_coding_genes_totally_found": 0,
        "num_genes_totally_found": 0,
        "ND2_present": False,
        "assembler": "",
    }
    if sample_dir is None or not sample_dir.is_dir():
        return empty

    mitoz_summary = mitoz_summary_for_sample_dir(sample_dir)
    novo_summary = novoplasty_annotation_summary_for_sample_dir(sample_dir)

    if mitoz_summary is not None and novo_summary is not None:
        mitoz_genes = mitoz_summary.genes_total
        novo_genes = novo_summary.genes_total

        if mitoz_genes == 37 and novo_genes != 37:
            return _metrics_from_summary(mitoz_summary, "mitoz")
        if novo_genes == 37 and mitoz_genes != 37:
            return _metrics_from_summary(novo_summary, "novoplasty")

        if mitoz_genes > novo_genes:
            better, better_name = mitoz_summary, "mitoz"
        elif novo_genes > mitoz_genes:
            better, better_name = novo_summary, "novoplasty"
        else:
            preferred = _preferred_metrics(mitoz_summary, novo_summary)
            return preferred if preferred is not None else empty

        if better.genes_total <= 37:
            return _metrics_from_summary(better, better_name)

        preferred = _preferred_metrics(mitoz_summary, novo_summary)
        return preferred if preferred is not None else empty

    preferred = _preferred_metrics(mitoz_summary, novo_summary)
    if preferred is not None:
        return preferred

    length_bp = assembly_length_bp_for_sample_dir(sample_dir)
    n_lengths = length_bp.count(";") + 1 if length_bp and length_bp != "0" else 0
    if n_lengths:
        circularity = ";".join(["no"] * n_lengths)
    else:
        circularity = ""
    return {
        "assembly_success": "failure",
        "length_bp": length_bp,
        "circularity": circularity,
        "num_protein_coding_genes_totally_found": 0,
        "num_trna_coding_genes_totally_found": 0,
        "num_rrna_coding_genes_totally_found": 0,
        "num_genes_totally_found": 0,
        "ND2_present": False,
        "assembler": "",
    }


def _annotation_files_under(suffix_dir: Path) -> list[Path]:
    """All regular files under one ``annotation/<suffix>/`` tree."""
    return [
        f
        for f in suffix_dir.rglob("*")
        if f.is_file() and not f.name.startswith("._")
    ]


def annotation_status_for_suffix_dir(suffix_dir: Path) -> str:
    """
    Classify one annotation attempt (``annotation/<suffix>/``).

    - ``not_attempted``: no annotation output files
    - ``failed``: files exist but every file is 0 bytes
    - ``success``: at least one non-zero output file
    """
    files = _annotation_files_under(suffix_dir)
    if not files:
        return "not_attempted"
    if all(f.stat().st_size == 0 for f in files):
        return "failed"
    return "success"


def annotation_success_for_sample_dir(sample_dir: Path) -> str:
    """
    Best annotation outcome for a sample directory.

    Scans every ``annotation/<suffix>/`` subtree. If any suffix run succeeded,
    the sample is ``success``; otherwise ``failed`` when outputs exist but are
    all empty; otherwise ``not_attempted``.
    """
    if not sample_dir.is_dir():
        return "not_attempted"
    ann = sample_dir / ANNOTATION_FOLDER
    if not ann.is_dir():
        return "not_attempted"

    statuses: list[str] = []
    for suffix_dir in sorted(ann.iterdir()):
        if not suffix_dir.is_dir() or suffix_dir.name.startswith("._"):
            continue
        statuses.append(annotation_status_for_suffix_dir(suffix_dir))

    if not statuses:
        return "not_attempted"
    if "success" in statuses:
        return "success"
    if "failed" in statuses:
        return "failed"
    return "not_attempted"


def parse_sample_label(sample_label: str) -> tuple[str, str]:
    """
    Split a batch annotation label into assembly directory name and suffix.

    ``Motacilla_..._ABJ133__ABJ133`` -> (``Motacilla_..._ABJ133``, ``ABJ133``)
    """
    if "__" in sample_label:
        asm_name, suffix = sample_label.split("__", 1)
        return asm_name, suffix
    return sample_label, "default"


def annotation_output_dir(sample_dir: Path, sample_label: str) -> Path:
    """Directory passed to annotate_mitogenome as ``output_root / sample_name``."""
    _, suffix = parse_sample_label(sample_label)
    return sample_dir / ANNOTATION_FOLDER / suffix


def mitoz_dir(sample_dir: Path, sample_label: str) -> Path:
    return annotation_output_dir(sample_dir, sample_label) / "mitoz"


def annotation_done(sample_dir: Path, sample_label: str) -> bool:
    """True if MitoZ wrote a ``*.result/`` directory for this job."""
    mz = mitoz_dir(sample_dir, sample_label)
    if not mz.is_dir():
        return False
    prefix = f"{sample_label}."
    if (mz / f"{sample_label}.result").is_dir():
        return True
    for child in mz.iterdir():
        if child.is_dir() and child.name.endswith(".result"):
            return True
        if child.is_dir() and child.name.startswith(prefix) and child.name.endswith(".result"):
            return True
    return False


def clear_annotation_output(sample_dir: Path, sample_label: str) -> list[str]:
    """Remove one annotation output tree under ``sample_dir/annotation/``."""
    target = annotation_output_dir(sample_dir, sample_label)
    if not target.exists():
        return []
    shutil.rmtree(target)
    return [str(target.relative_to(sample_dir))]


def clear_all_annotations(sample_dir: Path) -> list[str]:
    """Remove every ``annotation/`` subtree for an assembly directory."""
    ann = sample_dir / ANNOTATION_FOLDER
    if not ann.is_dir():
        return []
    removed = []
    for child in list(ann.iterdir()):
        if child.is_dir() and not child.name.startswith("._"):
            shutil.rmtree(child)
            removed.append(f"{ANNOTATION_FOLDER}/{child.name}")
    if ann.is_dir() and not any(ann.iterdir()):
        ann.rmdir()
    return removed
