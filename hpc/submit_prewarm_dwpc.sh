#!/bin/bash
#
# HPC array job: pre-warm the on-disk DWPC matrix cache for the web tool.
#
# One slurm array task per metapath. Cache lives on shared FS
# (data/dwpc_cache/) so the web app reads them after the array completes.
#
# Usage:
#   bash hpc/submit_prewarm_dwpc.sh                    # G -> BP
#   SOURCE_TYPE=G TARGET_TYPE=D bash hpc/submit_prewarm_dwpc.sh
#

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
cd "$REPO_ROOT"

SOURCE_TYPE="${SOURCE_TYPE:-G}"
TARGET_TYPE="${TARGET_TYPE:-BP}"
ARRAY_SIZE="${ARRAY_SIZE:-300}"  # oversized; extra tasks no-op
MEM="${MEM:-16G}"
TIME_LIMIT="${TIME_LIMIT:-02:00:00}"
PARTITION="${PARTITION:-amilan}"
QOS="${QOS:-normal}"

LOG_DIR="hpc/logs/prewarm_dwpc_${SOURCE_TYPE}_${TARGET_TYPE}"
mkdir -p "$LOG_DIR"

# Discover metapaths once on the submit host so the array tasks don't repeat it.
# Stored as a text file the tasks read with sed.
MP_LIST_FILE="$LOG_DIR/metapaths.txt"
module load anaconda >/dev/null 2>&1 || true
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate multi_dwpc
python3 scripts/prewarm_dwpc_cache.py \
    --source-type "$SOURCE_TYPE" --target-type "$TARGET_TYPE" \
    --list-metapaths > "$MP_LIST_FILE"
TOTAL=$(wc -l < "$MP_LIST_FILE" | tr -d ' ')
conda deactivate

if [[ "$TOTAL" -eq 0 ]]; then
    echo "No metapaths discovered for $SOURCE_TYPE -> $TARGET_TYPE; aborting."
    exit 1
fi

if [[ "$TOTAL" -lt "$ARRAY_SIZE" ]]; then
    ARRAY_RANGE="1-$TOTAL"
else
    ARRAY_RANGE="1-$ARRAY_SIZE"
fi

echo "Submitting $TOTAL prewarm tasks ($SOURCE_TYPE -> $TARGET_TYPE)"
echo "Metapath list: $MP_LIST_FILE"
echo "Log dir:       $LOG_DIR"

JOB_CMD='cd "'"$REPO_ROOT"'" && module load anaconda && source "$(conda info --base)/etc/profile.d/conda.sh" && conda activate multi_dwpc && \
    if [[ $SLURM_ARRAY_TASK_ID -gt '"$TOTAL"' ]]; then \
        echo "[task $SLURM_ARRAY_TASK_ID] no-op (only '"$TOTAL"' metapaths)"; exit 0; \
    fi && \
    METAPATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "'"$MP_LIST_FILE"'") && \
    echo "[task $SLURM_ARRAY_TASK_ID] metapath=$METAPATH" && \
    python3 scripts/prewarm_dwpc_cache.py --single-metapath "$METAPATH"'

sbatch \
    --parsable \
    --export=ALL \
    --job-name="prewarm-${SOURCE_TYPE}-${TARGET_TYPE}" \
    --partition="$PARTITION" --qos="$QOS" \
    --cpus-per-task=1 --mem="$MEM" --time="$TIME_LIMIT" \
    --array="$ARRAY_RANGE" \
    --output="$LOG_DIR/%A_%a.out" --error="$LOG_DIR/%A_%a.err" \
    --wrap="$JOB_CMD"
