#!/usr/bin/env python3
"""
Align complete/partial mitogenome assemblies by Collier et al. taxon_pair.

For each taxon_pair in mitogenome_assembly_results.csv with ≥2 usable
assemblies (assembly_success in {complete, partial}), gather one sequence
per sample, run MUSCLE, and write a MEGA-friendly aligned FASTA (.fas).

Default output root:
  /Volumes/KC_SSD/mtDNA_reassembly/aligned_mitogenomes
"""
from __future__ import annotations

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "scripts"))

from mitogenome_paths import (  # noqa: E402
    find_annotation_summary_txt,
    find_novoplasty_annotation_summary_txt,
    parse_fasta_records,
    resolve_sample_dir,
)

DEFAULT_RESULTS = REPO / "mitogenome_assembly_results.csv"
DEFAULT_OUT = Path("/Volumes/KC_SSD/mtDNA_reassembly/aligned_mitogenomes")
SCAFFOLD_JOINER = "N" * 20

# Short CSV taxon_pair codes -> PeerJ Table 1 style labels (for filenames).
TAXON_PAIR_LABELS: dict[str, str] = {
    "Anas acuta": "Anas_acuta",
    "Anas c": "Anas_carolinensis_A_crecca",
    "Anas p": "Anas_platyrhynchos",
    "Anser": "Anser_albifrons",
    "Anthus": "Anthus_rubescens",
    "Arenaria": "Arenaria_interpres_A_melanocephala",
    "Aythya": "Aythya_marila",
    "Calcarius": "Calcarius_lapponicus",
    "Calidris": "Calidris_alpina",
    "Clangula": "Clangula_hyemalis",
    "Corvus": "Corvus_corax",
    "Gallinago": "Gallinago_delicata_G_gallinago",
    "Gavia": "Gavia_stellata",
    "Histrionicus": "Histrionicus_histrionicus",
    "Lagopus": "Lagopus_lagopus",
    "Larus a": "Larus_argentatus",
    "Larus c": "Larus_brachyrhynchus_L_canus",
    "Leucosticte": "Leucosticte_arctoa_L_tephrocotis",
    "Luscinia": "Luscinia_svecica",
    "Mareca": "Mareca_americana_M_penelope",
    "Melanitta": "Melanitta_americana",
    "Mergus": "Mergus_merganser",
    "Motacilla": "Motacilla_tschutschensis",
    "Numenius": "Numenius_phaeopus",
    "Phylloscopus": "Phylloscopus_borealis_P_examinandus",
    "Pica": "Pica_hudsonia_P_pica",
    "Picoides": "Picoides_dorsalis_P_tridactylus",
    "Pinicola": "Pinicola_enucleator",
    "Plectrophenax": "Plectrophenax_hyperboreus_P_nivalis",
    "Pluvialis df": "Pluvialis_dominica_P_fulva",
    "Somateria m": "Somateria_mollissima",
    "Somateria s": "Somateria_spectabilis",
    "Spatula": "Spatula_clypeata",
    "Tringa bi": "Tringa_brevipes_T_incana",
    "Tringa nm": "Tringa_melanoleuca_T_nebularia",
    "Uria a": "Uria_aalge",
    "Uria l": "Uria_lomvia",
}


def find_muscle(explicit: str | None = None) -> str:
    candidates: list[str] = []
    if explicit:
        candidates.append(explicit)
    which = shutil.which("muscle")
    if which:
        candidates.append(which)
    candidates.extend(
        [
            str(REPO / "bin" / "run_muscle"),  # wrapper sets DYLD_LIBRARY_PATH
            str(REPO / "bin" / "muscle"),
            "/opt/anaconda3/bin/muscle",
            "/opt/homebrew/bin/muscle",
            "/usr/local/bin/muscle",
        ]
    )
    for path in candidates:
        p = Path(path)
        if p.is_file() and os_access_ok(p):
            return str(p)
    raise SystemExit(
        "MUSCLE executable not found. Place it at bin/muscle or pass --muscle."
    )


def os_access_ok(path: Path) -> bool:
    import os

    return os.access(path, os.X_OK)


def mega_safe_id(name: str) -> str:
    """MEGA-friendly sequence ID: alnum/underscore only."""
    cleaned = re.sub(r"[^A-Za-z0-9_]+", "_", name.strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned or "seq"


def output_stem(taxon_pair: str) -> str:
    if taxon_pair in TAXON_PAIR_LABELS:
        return TAXON_PAIR_LABELS[taxon_pair]
    return mega_safe_id(taxon_pair)


def gbf_fasta_beside_summary(summary: Path) -> Path | None:
    siblings = sorted(
        p
        for p in summary.parent.glob("*.gbf.fasta")
        if p.is_file() and not p.name.startswith("._")
    )
    return siblings[0] if siblings else None


def raw_novoplasty_assembly(sample_dir: Path) -> Path | None:
    from mitogenome_paths import candidate_assembly_paths

    best: Path | None = None
    best_tier = 99
    for novop in sorted(sample_dir.glob("novoplasty_*")):
        if not novop.is_dir():
            continue
        from mitogenome_paths import classify_novoplasty_dir

        tier = classify_novoplasty_dir(novop)
        paths = candidate_assembly_paths(novop)
        if paths and tier < best_tier:
            best_tier = tier
            best = paths[0]
    return best


def mitoz_mitogenome_fa(sample_dir: Path, summary: Path) -> Path | None:
    """Prefer the annotated gbf.fasta; else metaspades.mitogenome.fa near the result."""
    gbf = gbf_fasta_beside_summary(summary)
    if gbf is not None:
        return gbf
    for pat in (
        "*.metaspades.mitogenome.fa",
        "*mitogenome.fa",
    ):
        hits = [
            p
            for p in sample_dir.joinpath("mitoz").rglob(pat)
            if p.is_file() and not p.name.startswith("._") and "result" not in p.name
        ]
        if hits:
            return max(hits, key=lambda p: p.stat().st_size)
    hits = list(sample_dir.joinpath("mitoz").rglob("*.metaspades.mitogenome.fa"))
    hits = [p for p in hits if p.is_file() and not p.name.startswith("._")]
    return max(hits, key=lambda p: p.stat().st_size) if hits else None


def assembly_fasta_for_row(sample_dir: Path, assembler: str) -> Path | None:
    assembler = (assembler or "").strip().lower()
    if assembler == "mitoz":
        summary = find_annotation_summary_txt(sample_dir, "mitoz")
        if summary is None:
            return None
        return mitoz_mitogenome_fa(sample_dir, summary)
    if assembler == "novoplasty":
        summary = find_novoplasty_annotation_summary_txt(sample_dir)
        if summary is not None:
            gbf = gbf_fasta_beside_summary(summary)
            if gbf is not None:
                return gbf
        return raw_novoplasty_assembly(sample_dir)
    # Unknown assembler: try mitoz then novoplasty
    for asm in ("mitoz", "novoplasty"):
        path = assembly_fasta_for_row(sample_dir, asm)
        if path is not None:
            return path
    return None


def one_sequence_from_fasta(path: Path) -> tuple[str, int]:
    """
    Return (sequence, n_scaffolds).

    Multi-scaffold assemblies are concatenated with N-spacers so each sample
    contributes exactly one sequence to the taxon alignment (MEGA-friendly).
    """
    records = parse_fasta_records(path)
    if not records:
        raise ValueError(f"No FASTA records in {path}")
    seqs = [seq.upper().replace(" ", "").replace("\t", "") for _, seq in records]
    if len(seqs) == 1:
        return seqs[0], 1
    return SCAFFOLD_JOINER.join(seqs), len(seqs)


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as fh:
        for seq_id, seq in records:
            fh.write(f">{seq_id}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def run_muscle(
    muscle: str,
    infile: Path,
    outfile: Path,
    *,
    algorithm: str = "super5",
) -> None:
    """Run MUSCLE 5 (-super5/-align) or classic MUSCLE (-in/-out).

    Default ``super5`` is much faster for multi-kb mitogenomes with indels;
    use ``align`` for the slower full progressive algorithm.
    """
    help_txt = subprocess.run(
        [muscle, "-h"],
        capture_output=True,
        text=True,
        check=False,
    )
    blob = (help_txt.stdout or "") + (help_txt.stderr or "")
    is_v5 = "-align" in blob or "MUSCLE v5" in blob or "muscle 5" in blob.lower()
    if is_v5:
        mode = "-super5" if algorithm == "super5" else "-align"
        cmd = [muscle, mode, str(infile), "-output", str(outfile)]
    else:
        cmd = [muscle, "-in", str(infile), "-out", str(outfile)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not outfile.is_file() or outfile.stat().st_size == 0:
        raise RuntimeError(
            f"MUSCLE failed ({proc.returncode}): {' '.join(cmd)}\n"
            f"stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}"
        )


def load_rows(results_csv: Path) -> list[dict[str, str]]:
    with results_csv.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results-csv", type=Path, default=DEFAULT_RESULTS)
    p.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    p.add_argument("--muscle", default=None, help="Path to MUSCLE executable")
    p.add_argument(
        "--taxon-pair",
        action="append",
        default=[],
        help="Restrict to one or more taxon_pair codes (repeatable)",
    )
    p.add_argument(
        "--min-seqs",
        type=int,
        default=2,
        help="Skip taxon pairs with fewer than this many sequences (default: 2)",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve inputs and print plans without running MUSCLE",
    )
    p.add_argument(
        "--skip-existing",
        action="store_true",
        help="Skip taxon pairs that already have a non-empty *_aligned.fas",
    )
    p.add_argument(
        "--algorithm",
        choices=("super5", "align"),
        default="super5",
        help="MUSCLE 5 algorithm (default: super5; use align for full progressive)",
    )
    args = p.parse_args()

    muscle = None if args.dry_run else find_muscle(args.muscle)
    rows = load_rows(args.results_csv)
    wanted = set(args.taxon_pair) if args.taxon_pair else None

    by_taxon: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        status = (row.get("assembly_success") or "").strip().lower()
        if status not in ("complete", "partial"):
            continue
        taxon = (row.get("taxon_pair") or "").strip()
        if not taxon:
            continue
        if wanted is not None and taxon not in wanted:
            continue
        by_taxon[taxon].append(row)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.out_dir / "alignment_manifest.tsv"
    log_path = args.out_dir / "alignment_run.log"

    manifest_rows: list[dict[str, str]] = []
    log_lines: list[str] = []

    print(f"Taxon pairs to process: {len(by_taxon)}")
    if muscle:
        print(f"MUSCLE: {muscle}  algorithm={args.algorithm}")
    print(f"Output: {args.out_dir}")

    for taxon in sorted(by_taxon, key=lambda t: output_stem(t).lower()):
        stem = output_stem(taxon)
        out_fas = args.out_dir / f"{stem}_mitogenomes_aligned.fas"
        unaligned = args.out_dir / f"{stem}_mitogenomes_unaligned.fas"

        if (
            args.skip_existing
            and out_fas.is_file()
            and out_fas.stat().st_size > 0
            and not args.dry_run
        ):
            msg = f"SKIP_EXISTING {taxon} ({stem}): {out_fas.name}"
            print(msg)
            log_lines.append(msg)
            continue

        records: list[tuple[str, str]] = []
        notes: list[str] = []
        for row in sorted(by_taxon[taxon], key=lambda r: r["directory_name"]):
            name = row["directory_name"]
            sample_dir = resolve_sample_dir(name)
            if sample_dir is None:
                notes.append(f"MISSING_DIR\t{name}")
                continue
            asm = (row.get("assembler") or "").strip()
            fasta = assembly_fasta_for_row(sample_dir, asm)
            if fasta is None:
                notes.append(f"MISSING_FASTA\t{name}\tassembler={asm}")
                continue
            try:
                seq, n_scaff = one_sequence_from_fasta(fasta)
            except ValueError as exc:
                notes.append(f"BAD_FASTA\t{name}\t{exc}")
                continue
            seq_id = mega_safe_id(name)
            records.append((seq_id, seq))
            manifest_rows.append(
                {
                    "taxon_pair": taxon,
                    "output_stem": stem,
                    "directory_name": name,
                    "seq_id": seq_id,
                    "assembly_success": row.get("assembly_success", ""),
                    "assembler": asm,
                    "source_fasta": str(fasta),
                    "n_scaffolds": str(n_scaff),
                    "ungapped_length": str(len(seq.replace("N", ""))),
                    "sequence_length": str(len(seq)),
                }
            )

        if len(records) < args.min_seqs:
            msg = (
                f"SKIP {taxon} ({stem}): {len(records)} sequence(s) "
                f"< min {args.min_seqs}"
            )
            print(msg)
            log_lines.append(msg)
            for n in notes:
                log_lines.append(f"  {n}")
            continue

        write_fasta(unaligned, records)
        msg = f"ALIGN {taxon} ({stem}): n={len(records)} -> {out_fas.name}"
        print(msg)
        log_lines.append(msg)
        for n in notes:
            log_lines.append(f"  {n}")

        if args.dry_run:
            continue

        with tempfile.TemporaryDirectory(prefix="muscle_") as tmp:
            tmp_in = Path(tmp) / "in.fas"
            tmp_out = Path(tmp) / "out.fas"
            shutil.copy(unaligned, tmp_in)
            try:
                run_muscle(
                    muscle,
                    tmp_in,
                    tmp_out,
                    algorithm=args.algorithm,
                )
            except RuntimeError as exc:
                err = f"FAIL {taxon}: {exc}"
                print(err)
                log_lines.append(err)
                continue
            # Rewrite with consistent line wrapping for MEGA
            aln_records = parse_fasta_records(tmp_out)
            # Preserve input order when possible
            order = {sid: i for i, (sid, _) in enumerate(records)}
            aln_records.sort(key=lambda r: order.get(r[0].split()[0], 10**9))
            cleaned = [(mega_safe_id(h.split()[0]), s) for h, s in aln_records]
            write_fasta(out_fas, cleaned)
            print(f"  wrote {out_fas.name} ({out_fas.stat().st_size} bytes)")
            log_lines.append(f"  wrote {out_fas.name}")

    with manifest_path.open("w", newline="", encoding="utf-8") as fh:
        fields = [
            "taxon_pair",
            "output_stem",
            "directory_name",
            "seq_id",
            "assembly_success",
            "assembler",
            "source_fasta",
            "n_scaffolds",
            "ungapped_length",
            "sequence_length",
        ]
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for row in manifest_rows:
            w.writerow(row)

    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    print(f"Manifest: {manifest_path}")
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()
