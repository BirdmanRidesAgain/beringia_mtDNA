#!/usr/bin/env python3
"""
Move completed sample output trees to the SSD and remove them from local disk.

When a sample has a complete local ``mitoz/`` run, the entire sample directory is
copied (``novoplasty_*``, ``mitoz/``, ``annotation/``, etc.) to::

    /Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/assembled_mitogenomes/<sample>/

Samples whose ``mitoz/`` was archived earlier but still have other local outputs
(e.g. NOVOPlasty) are merged into the existing SSD sample directory (backfill).

Rsync merges into any partial archive left by an interrupted prior run. Local
directories are deleted only after the SSD copy passes verification.

Incomplete local ``mitoz/`` dirs can be removed with ``--remove-incomplete``.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPTS))

from mitogenome_paths import REPO  # noqa: E402
from trim_assemble_annotate_mitoz import mitoz_done  # noqa: E402

DEFAULT_ARCHIVE_ROOT = Path(
    "/Volumes/KC_SSD/keiler_C/KAC_mtDNA_reassembly/assembled_mitogenomes"
)
DEFAULT_OUTPUT_ROOT = REPO / "mitogenomes_output"
LOG_PATH = DEFAULT_OUTPUT_ROOT / "_batch_logs" / "archive_mitoz_to_ssd.log"
LOCAL_DISK_PATH = Path("/System/Volumes/Data")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Archive completed mitogenome sample dirs (full tree) to SSD."
    )
    p.add_argument(
        "--archive-root",
        type=Path,
        default=DEFAULT_ARCHIVE_ROOT,
        help=f"SSD archive root (default: {DEFAULT_ARCHIVE_ROOT})",
    )
    p.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help=f"Local mitogenomes_output root (default: {DEFAULT_OUTPUT_ROOT})",
    )
    p.add_argument(
        "--remove-incomplete",
        action="store_true",
        help="Delete local mitoz/ dirs that are not complete (for rerun)",
    )
    p.add_argument(
        "--min-free-gb",
        type=float,
        default=0.0,
        help="Stop archiving once local free space reaches this many GB (0 = no limit)",
    )
    p.add_argument("--dry-run", action="store_true", help="Print actions only")
    return p.parse_args()


def log_write(log, msg: str) -> None:
    line = msg.rstrip()
    print(line, flush=True)
    if log is not None:
        log.write(line + "\n")
        log.flush()


def dir_size_bytes(path: Path) -> int:
    if not path.is_dir():
        return 0
    try:
        out = subprocess.run(
            ["du", "-sk", str(path)],
            check=True,
            capture_output=True,
            text=True,
        )
        kb = int(out.stdout.split()[0])
        return kb * 1024
    except (OSError, subprocess.CalledProcessError, ValueError):
        total = 0
        for p in path.rglob("*"):
            if p.is_file():
                try:
                    total += p.stat().st_size
                except OSError:
                    pass
        return total


def local_free_gb() -> float:
    path = LOCAL_DISK_PATH if LOCAL_DISK_PATH.is_dir() else Path("/")
    return shutil.disk_usage(path).free / (1024**3)


def log_local_free_space(log) -> None:
    log_write(log, f"local free space (APFS available): {local_free_gb():.1f} GB")
    log_write(
        log,
        "  note: macOS System Settings may show more space if purgeable caches "
        "or Time Machine local snapshots can be reclaimed",
    )


def rsync_tree(src: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["rsync", "-a", "--partial", f"{src}/", f"{dest}/"],
        check=True,
    )


def sample_dir_has_content(sample_dir: Path) -> bool:
    for entry in sample_dir.iterdir():
        if entry.name.startswith("."):
            continue
        return True
    return False


def ssd_mitoz_complete(dest_sample_dir: Path, sample_name: str) -> bool:
    mitoz_dir = dest_sample_dir / "mitoz"
    return mitoz_dir.is_dir() and mitoz_done(mitoz_dir, sample_name)


def archive_sample_complete(
    src: Path,
    dest: Path,
    sample_name: str,
    *,
    src_size: int,
    pre_dest_size: int | None = None,
    min_dest_ratio: float = 0.95,
) -> bool:
    if not ssd_mitoz_complete(dest, sample_name):
        return False
    dest_size = dir_size_bytes(dest)
    if pre_dest_size is not None:
        return dest_size >= pre_dest_size + src_size * min_dest_ratio
    if src_size == 0:
        return dest_size > 0
    return dest_size >= src_size * min_dest_ratio


def archive_sample_tree(
    src: Path,
    dest: Path,
    *,
    sample_name: str,
    src_size: int,
    dry_run: bool,
    log,
) -> str:
    """
    Archive ``src`` sample directory to ``dest``.

    Returns one of: ``archived``, ``already_archived``.
    """
    if dry_run:
        return "archived"

    pre_dest_size = dir_size_bytes(dest) if dest.is_dir() else 0
    if dest.is_dir() and archive_sample_complete(
        src, dest, sample_name, src_size=src_size
    ):
        shutil.rmtree(src)
        return "already_archived"

    log_write(log, "  rsync in progress...")
    rsync_tree(src, dest)

    if not ssd_mitoz_complete(dest, sample_name):
        raise RuntimeError("SSD archive missing complete mitoz/ after rsync (local copy kept)")

    new_dest_size = dir_size_bytes(dest)
    expected_min = pre_dest_size + src_size * 0.95
    if src_size > 0 and new_dest_size + 1024 * 1024 < expected_min:
        raise RuntimeError(
            f"SSD archive incomplete after rsync "
            f"(local {src_size / (1024**3):.2f} GB, SSD {new_dest_size / (1024**3):.2f} GB; "
            "local copy kept)"
        )

    shutil.rmtree(src)
    return "archived"


def archive_one_sample(
    sample_name: str,
    *,
    archive_root: Path = DEFAULT_ARCHIVE_ROOT,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
    log=None,
    dry_run: bool = False,
) -> str | None:
    """
    Archive one local sample directory when ``mitoz/`` is complete.

    Returns ``archived``, ``already_archived``, or ``None`` if there is nothing to archive.
    """
    sample_dir = output_root / sample_name
    if not sample_dir.is_dir() or not sample_dir_has_content(sample_dir):
        return None
    mitoz_dir = sample_dir / "mitoz"
    if not mitoz_dir.is_dir() or not mitoz_done(mitoz_dir, sample_name):
        return None
    dest = archive_root / sample_name
    src_size = dir_size_bytes(sample_dir)
    return archive_sample_tree(
        sample_dir,
        dest,
        sample_name=sample_name,
        src_size=src_size,
        dry_run=dry_run,
        log=log,
    )


def collect_candidates(
    output_root: Path,
    archive_root: Path,
) -> list[tuple[int, Path, Path, str]]:
    candidates: list[tuple[int, Path, Path, str]] = []

    for sample_dir in sorted(output_root.iterdir()):
        if not sample_dir.is_dir() or sample_dir.name.startswith(("_", ".")):
            continue
        if not sample_dir_has_content(sample_dir):
            continue

        dest = archive_root / sample_dir.name
        local_mitoz = sample_dir / "mitoz"
        local_mitoz_complete = local_mitoz.is_dir() and mitoz_done(
            local_mitoz, sample_dir.name
        )

        if local_mitoz_complete:
            size = dir_size_bytes(sample_dir)
            candidates.append((size, sample_dir, dest, "full"))
            continue

        if local_mitoz.is_dir():
            continue

        if dest.is_dir() and ssd_mitoz_complete(dest, sample_dir.name):
            size = dir_size_bytes(sample_dir)
            candidates.append((size, sample_dir, dest, "backfill"))

    candidates.sort(key=lambda item: item[0], reverse=True)
    return candidates


def main() -> None:
    args = parse_args()
    archive_root = args.archive_root.expanduser().resolve()
    output_root = args.output_root.expanduser().resolve()

    if not archive_root.parent.is_dir():
        raise SystemExit(f"SSD path not available: {archive_root.parent}")

    archived = 0
    already_archived = 0
    backfilled = 0
    removed_incomplete = 0
    errors: list[str] = []

    print("Scanning for sample dirs to archive...", flush=True)
    candidates = collect_candidates(output_root, archive_root)
    full_n = sum(1 for *_, kind in candidates if kind == "full")
    backfill_n = sum(1 for *_, kind in candidates if kind == "backfill")

    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8") as log:
        log_write(log, f"\n--- ARCHIVE START {datetime.now(timezone.utc).isoformat()}Z ---")
        log_write(log, f"archive_root: {archive_root}")
        log_write(log, f"output_root: {output_root}")
        log_write(
            log,
            f"sample dirs to archive: {len(candidates)} "
            f"(full={full_n}, backfill={backfill_n})",
        )
        log_local_free_space(log)

        for size_bytes, sample_dir, dest, kind in candidates:
            if args.min_free_gb > 0 and local_free_gb() >= args.min_free_gb:
                log_write(
                    log,
                    f"STOP: local free space {local_free_gb():.1f} GB "
                    f">= --min-free-gb {args.min_free_gb:g}",
                )
                break

            name = sample_dir.name
            size_gb = size_bytes / (1024**3)
            label = "BACKFILL" if kind == "backfill" else "ARCHIVE"
            try:
                log_write(log, f"{label} {name} ({size_gb:.2f} GB) -> {dest}")
                t0 = time.monotonic()
                outcome = archive_sample_tree(
                    sample_dir,
                    dest,
                    sample_name=name,
                    src_size=size_bytes,
                    dry_run=args.dry_run,
                    log=log,
                )
                elapsed = time.monotonic() - t0
                if outcome == "already_archived":
                    already_archived += 1
                    log_write(log, f"  OK already on SSD; removed local ({elapsed / 60:.1f} min)")
                else:
                    archived += 1
                    if kind == "backfill":
                        backfilled += 1
                    log_write(
                        log,
                        f"  OK rsync+verify ({elapsed / 60:.1f} min); "
                        f"local free {local_free_gb():.1f} GB",
                    )
            except Exception as exc:  # noqa: BLE001
                msg = f"ERROR archiving {name}: {exc}"
                errors.append(msg)
                log_write(log, msg)

        if args.remove_incomplete:
            for sample_dir in sorted(output_root.iterdir()):
                if not sample_dir.is_dir() or sample_dir.name.startswith(("_", ".")):
                    continue
                mitoz_dir = sample_dir / "mitoz"
                if not mitoz_dir.is_dir() or mitoz_done(mitoz_dir, sample_dir.name):
                    continue
                try:
                    log_write(log, f"REMOVE incomplete mitoz: {mitoz_dir}")
                    if not args.dry_run:
                        shutil.rmtree(mitoz_dir)
                    removed_incomplete += 1
                except Exception as exc:  # noqa: BLE001
                    msg = f"ERROR removing {sample_dir.name}: {exc}"
                    errors.append(msg)
                    log_write(log, msg)

        log_write(
            log,
            f"--- ARCHIVE END archived={archived} backfilled={backfilled} "
            f"already_archived={already_archived} "
            f"removed_incomplete={removed_incomplete} errors={len(errors)} "
            f"local_free_gb={local_free_gb():.1f} ---",
        )

    if errors:
        print(f"Finished with {len(errors)} error(s); see {LOG_PATH}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
