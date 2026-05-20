#!/usr/bin/env bash
# Merge multiple R1/R2 FASTQ.gz files per sample subdirectory into one pair.
#
# Uses: zcat file1 file2 ... | pigz -p N > sample_R1.fq.gz
#
# Skips subdirectories that already have exactly one R1 and one R2 (after ignoring
# AppleDouble "._*" files). Writes merged outputs as <subdir>_R1.fq.gz and
# <subdir>_R2.fq.gz, then removes the source files.
#
# Example:
#   ./scripts/merge_raw_read_pairs.sh -n "/Volumes/Seagate Por/KAC_mtDNA_reassembly/raw_data"
#   ./scripts/merge_raw_read_pairs.sh "/Volumes/Seagate Por/KAC_mtDNA_reassembly/raw_data"

set -euo pipefail

RAW_ROOT=""
DRY_RUN=0
PIGZ_THREADS="${PIGZ_THREADS:-4}"
LOG=""

usage() {
  sed -n '2,12p' "$0"
  echo "Options:"
  echo "  -n, --dry-run     Print actions only"
  echo "  --log FILE        Append log (default: <repo>/mitogenomes_output/merge_raw_read_pairs.log)"
  echo "  -j N              pigz threads (default: 4)"
  exit "${1:-0}"
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    -n|--dry-run) DRY_RUN=1; shift ;;
    -j) PIGZ_THREADS="$2"; shift 2 ;;
    --log) LOG="$2"; shift 2 ;;
    -h|--help) usage 0 ;;
    -*) echo "Unknown option: $1" >&2; usage 1 ;;
    *)
      if [[ -z "$RAW_ROOT" ]]; then
        RAW_ROOT="$1"
      else
        echo "Unexpected argument: $1" >&2
        usage 1
      fi
      shift
      ;;
  esac
done

if [[ -z "$RAW_ROOT" ]]; then
  echo "ERROR: RAW_DATA root directory required." >&2
  usage 1
fi

if [[ ! -d "$RAW_ROOT" ]]; then
  echo "ERROR: Not a directory: $RAW_ROOT" >&2
  exit 1
fi

REPO="$(cd "$(dirname "$0")/.." && pwd)"
if [[ -z "$LOG" ]]; then
  LOG="$REPO/mitogenomes_output/merge_raw_read_pairs.log"
fi
mkdir -p "$(dirname "$LOG")"

if ! command -v pigz >/dev/null 2>&1; then
  echo "ERROR: pigz not found in PATH." >&2
  exit 1
fi

log() {
  local msg="[$(date -u +%Y-%m-%dT%H:%M:%SZ)] $*"
  echo "$msg"
  echo "$msg" >>"$LOG"
}

is_fastq_gz() {
  local n="${1,,}"
  [[ "$n" == *.fq.gz || "$n" == *.fastq.gz ]]
}

is_r1() {
  local n="${1,,}"
  is_fastq_gz "$n" || return 1
  [[ "$n" == *_r1.fq.gz || "$n" == *_r1.fastq.gz ]] && return 0
  [[ "$n" =~ _1\.(fq|fastq)\.gz$ ]] && return 0
  return 1
}

is_r2() {
  local n="${1,,}"
  is_fastq_gz "$n" || return 1
  [[ "$n" == *_r2.fq.gz || "$n" == *_r2.fastq.gz ]] && return 0
  [[ "$n" =~ _2\.(fq|fastq)\.gz$ ]] && return 0
  return 1
}

collect_reads() {
  local dir="$1"
  local which="$2" # R1 or R2
  local -a out=()
  local f
  shopt -s nullglob
  for f in "$dir"/*; do
    [[ -f "$f" ]] || continue
    local base
    base="$(basename "$f")"
    [[ "$base" == ._* ]] && continue
    if [[ "$which" == R1 ]] && is_r1 "$base"; then
      out+=("$f")
    elif [[ "$which" == R2 ]] && is_r2 "$base"; then
      out+=("$f")
    fi
  done
  shopt -u nullglob
  if [[ ${#out[@]} -gt 0 ]]; then
    printf '%s\n' "${out[@]}" | LC_ALL=C sort
  fi
}

merge_group() {
  local dir="$1"
  local which="$2" # R1 or R2
  shift 2
  local -a inputs=("$@")
  local sample
  sample="$(basename "$dir")"
  local out_final="$dir/${sample}_${which}.fq.gz"
  local out_tmp="${out_final}.merge_tmp"

  log "  $which: merge ${#inputs[@]} file(s) -> $(basename "$out_final")"
  for f in "${inputs[@]}"; do
    log "      <- $(basename "$f")"
  done

  if [[ "$DRY_RUN" -eq 1 ]]; then
    return 0
  fi

  # Exclude temp/output from inputs if re-running
  local -a to_cat=()
  local f
  for f in "${inputs[@]}"; do
    [[ "$f" == "$out_final" || "$f" == "$out_tmp" ]] && continue
    to_cat+=("$f")
  done

  zcat "${to_cat[@]}" | pigz -p "$PIGZ_THREADS" > "$out_tmp"
  if [[ ! -s "$out_tmp" ]]; then
    echo "ERROR: empty merge output: $out_tmp" >&2
    rm -f "$out_tmp"
    return 1
  fi

  for f in "${to_cat[@]}"; do
    rm -f "$f"
  done
  mv -f "$out_tmp" "$out_final"
}

merge_sample_dir() {
  local sample="$1"
  local dir="$RAW_ROOT/$sample"
  mapfile -t r1s < <(collect_reads "$dir" R1)
  mapfile -t r2s < <(collect_reads "$dir" R2)

  if [[ ${#r1s[@]} -le 1 && ${#r2s[@]} -le 1 ]]; then
    return 0
  fi
  if [[ ${#r1s[@]} -eq 0 || ${#r2s[@]} -eq 0 ]]; then
    log "SKIP $sample: R1=${#r1s[@]} R2=${#r2s[@]} (unpaired)"
    return 1
  fi

  log "MERGE $sample (R1=${#r1s[@]}, R2=${#r2s[@]})"
  merge_group "$dir" R1 "${r1s[@]}"
  merge_group "$dir" R2 "${r2s[@]}"

  # Remove macOS AppleDouble sidecars if present
  if [[ "$DRY_RUN" -eq 0 ]]; then
    rm -f "$dir"/._* 2>/dev/null || true
  fi
}

log "=== merge_raw_read_pairs.sh ==="
log "ROOT=$RAW_ROOT dry_run=$DRY_RUN pigz_threads=$PIGZ_THREADS"

merged=0
skipped=0
failed=0

shopt -s nullglob
for dir in "$RAW_ROOT"/*/; do
  sample="$(basename "$dir")"
  [[ "$sample" == .* ]] && continue
  mapfile -t r1s < <(collect_reads "$dir" R1)
  mapfile -t r2s < <(collect_reads "$dir" R2)
  if [[ ${#r1s[@]} -le 1 && ${#r2s[@]} -le 1 ]]; then
    ((skipped++)) || true
    continue
  fi
  if merge_sample_dir "$sample"; then
    ((merged++)) || true
  else
    ((failed++)) || true
  fi
done
shopt -u nullglob

log "Done: merged=$merged skipped=$skipped failed=$failed"

if [[ "$failed" -gt 0 ]]; then
  exit 1
fi
