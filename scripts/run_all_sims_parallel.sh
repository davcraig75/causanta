#!/bin/bash
# Re-run every parameter file in causanta/simulate/params/ in parallel.
# With cleanup turned off (the new core.py default), every per-hour snapshot
# is preserved on disk for inspection.
#
# Disk budget: ~30 GB total (mostly the 6mm runs at ~1.5 GB each in
# intermediate snapshots, ~1.5 GB in final + viz).
set -u

cd "$(dirname "$0")/.."
source .venv/bin/activate

LOGDIR=output/rerun_logs
mkdir -p "$LOGDIR"

START_ALL=$(date +%s)
echo "$(date): launching parallel re-run of all parameter files" \
    > "$LOGDIR/rerun_all.log"

# The 33 parameter files: small (1mm) + large (2mm) + xlarge (6mm) + all multiseed
PARAMS=(
    # 1mm robustness (seed 42)
    causanta/simulate/params/robustness_baseline.json
    causanta/simulate/params/robustness_reduced.json
    causanta/simulate/params/robustness_removed.json
    # 2mm large (seed 42)
    causanta/simulate/params/large_baseline.json
    causanta/simulate/params/large_reduced.json
    causanta/simulate/params/large_removed.json
    # 6mm xlarge (seed 42)
    causanta/simulate/params/xlarge_baseline.json
    causanta/simulate/params/xlarge_reduced.json
    causanta/simulate/params/xlarge_removed.json
    # 2mm multiseed baseline (seeds 43-46)
    causanta/simulate/params/multiseed_2mm_seed43.json
    causanta/simulate/params/multiseed_2mm_seed44.json
    causanta/simulate/params/multiseed_2mm_seed45.json
    causanta/simulate/params/multiseed_2mm_seed46.json
    # 2mm multiseed reduced (seeds 43-46)
    causanta/simulate/params/multiseed_2mm_reduced_seed43.json
    causanta/simulate/params/multiseed_2mm_reduced_seed44.json
    causanta/simulate/params/multiseed_2mm_reduced_seed45.json
    causanta/simulate/params/multiseed_2mm_reduced_seed46.json
    # 2mm multiseed removed (seeds 43-46)
    causanta/simulate/params/multiseed_2mm_removed_seed43.json
    causanta/simulate/params/multiseed_2mm_removed_seed44.json
    causanta/simulate/params/multiseed_2mm_removed_seed45.json
    causanta/simulate/params/multiseed_2mm_removed_seed46.json
    # 6mm multiseed baseline (seeds 43-46)
    causanta/simulate/params/multiseed_6mm_baseline_seed43.json
    causanta/simulate/params/multiseed_6mm_baseline_seed44.json
    causanta/simulate/params/multiseed_6mm_baseline_seed45.json
    causanta/simulate/params/multiseed_6mm_baseline_seed46.json
    # 6mm multiseed reduced (seeds 43-46)
    causanta/simulate/params/multiseed_6mm_reduced_seed43.json
    causanta/simulate/params/multiseed_6mm_reduced_seed44.json
    causanta/simulate/params/multiseed_6mm_reduced_seed45.json
    causanta/simulate/params/multiseed_6mm_reduced_seed46.json
    # 6mm multiseed removed (seeds 43-46)
    causanta/simulate/params/multiseed_6mm_removed_seed43.json
    causanta/simulate/params/multiseed_6mm_removed_seed44.json
    causanta/simulate/params/multiseed_6mm_removed_seed45.json
    causanta/simulate/params/multiseed_6mm_removed_seed46.json
)

PIDS=()
for param in "${PARAMS[@]}"; do
    [[ -f "$param" ]] || { echo "MISSING: $param" >> "$LOGDIR/rerun_all.log"; continue; }
    name=$(basename "$param" .json)
    log="$LOGDIR/${name}.log"
    python -m causanta.simulate.core "$param" > "$log" 2>&1 &
    PIDS+=($!)
    echo "  started $name (pid=$!)" >> "$LOGDIR/rerun_all.log"
done

echo "$(date): all ${#PIDS[@]} sims launched; waiting for completion..." >> "$LOGDIR/rerun_all.log"

# Wait for every sim
for pid in "${PIDS[@]}"; do
    wait "$pid"
done

ELAPSED=$(( $(date +%s) - START_ALL ))
echo "$(date): ALL ${#PIDS[@]} sims complete in ${ELAPSED}s" >> "$LOGDIR/rerun_all.log"
date > "$LOGDIR/rerun_all.done"
