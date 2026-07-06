#!/usr/bin/env bash
set -uo pipefail
cd "$(dirname "$0")/.."
LOG_DIR=mitogenomes_output/_batch_logs
mkdir -p "$LOG_DIR"

MIN_FREE_GB="${MIN_FREE_GB:-30}"

echo "=== archive start $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
# caffeinate -i keeps the system awake during multi-hour rsync transfers
if ! caffeinate -i python scripts/archive_mitoz_to_ssd.py --remove-incomplete; then
  echo "WARNING: archive finished with errors (see archive_mitoz_to_ssd.log)"
fi

FREE_GB=$(python3 -c "import shutil; print(int(shutil.disk_usage('/System/Volumes/Data').free / 1024**3))")
echo "=== local free space after archive (APFS available): ${FREE_GB} GB ==="
echo "    (System Settings may show more if purgeable snapshots/caches can be reclaimed)"

if [ "${FREE_GB}" -lt "${MIN_FREE_GB}" ]; then
  echo "ERROR: need at least ${MIN_FREE_GB} GB free to start batch (have ${FREE_GB} GB)"
  exit 1
fi

echo "=== pending batch start $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
python scripts/batch_mitoz_ssd.py --pending-mitoz --skip-completed --fresh-log

echo "=== done $(date -u +%Y-%m-%dT%H:%M:%SZ) ==="
