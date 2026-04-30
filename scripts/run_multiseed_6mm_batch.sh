#!/bin/bash
# Run 12 multiseed 6mm simulations as 3 sequential batches of 4 (one per scenario).
# Wall time: ~3.5-4 hr (each batch is ~75 min single-threaded wall time;
# 4-way parallel on each batch incurs ~10-25% slowdown from CPU contention).
set -u

cd "$(dirname "$0")/.."
source .venv/bin/activate

LOGDIR=output/multiseed_logs
mkdir -p "$LOGDIR"

START_ALL=$(date +%s)
echo "$(date): starting 6mm multiseed batch (baseline -> reduced -> removed)" \
    > "$LOGDIR/6mm_batch.log"

run_scenario_batch() {
    local scenario="$1"
    local batch_start=$(date +%s)
    echo "$(date): batch start: $scenario (4 sims, parallel)" >> "$LOGDIR/6mm_batch.log"

    local pids=()
    for seed in 43 44 45 46; do
        local param="causanta/simulate/params/multiseed_6mm_${scenario}_seed${seed}.json"
        local log="$LOGDIR/6mm_${scenario}_seed${seed}.log"
        python -m causanta.simulate.core "$param" > "$log" 2>&1 &
        pids+=($!)
        echo "  started $scenario seed=$seed (pid=$!)" >> "$LOGDIR/6mm_batch.log"
    done

    for pid in "${pids[@]}"; do
        wait "$pid"
    done

    local elapsed=$(( $(date +%s) - batch_start ))
    echo "$(date): batch end: $scenario in ${elapsed}s" >> "$LOGDIR/6mm_batch.log"
    date > "$LOGDIR/6mm_${scenario}.done"
}

run_scenario_batch baseline
run_scenario_batch reduced
run_scenario_batch removed

ELAPSED=$(( $(date +%s) - START_ALL ))
echo "$(date): ALL 12 sims complete in ${ELAPSED}s" >> "$LOGDIR/6mm_batch.log"
date > "$LOGDIR/6mm_batch.done"
