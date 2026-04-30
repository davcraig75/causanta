#!/bin/bash
# Run 8 multiseed 2mm simulations in parallel (reduced + removed scenarios)
# Wall time: ~25-40 min when run alongside the 6mm batch.
set -u

cd "$(dirname "$0")/.."
source .venv/bin/activate

LOGDIR=output/multiseed_logs
mkdir -p "$LOGDIR"

START=$(date +%s)
echo "$(date): starting 2mm multiseed batch (reduced + removed × seeds 43-46)" \
    > "$LOGDIR/2mm_batch.log"

PIDS=()
for scenario in reduced removed; do
    for seed in 43 44 45 46; do
        param="causanta/simulate/params/multiseed_2mm_${scenario}_seed${seed}.json"
        log="$LOGDIR/2mm_${scenario}_seed${seed}.log"
        python -m causanta.simulate.core "$param" > "$log" 2>&1 &
        PIDS+=($!)
        echo "  started $scenario seed=$seed (pid=$!)" >> "$LOGDIR/2mm_batch.log"
    done
done

# Wait for all 8 to finish
for pid in "${PIDS[@]}"; do
    wait "$pid"
done

ELAPSED=$(( $(date +%s) - START ))
echo "$(date): all 8 sims complete in ${ELAPSED}s" >> "$LOGDIR/2mm_batch.log"
date > "$LOGDIR/2mm_batch.done"
