#!/bin/bash
#
# Build the laptop null bundle on Alpine: one array task per G->BP metapath
# (computing its DWPC matrix into CACHE_DIR if absent), then a finalize job that
# merges the parts into BUNDLE_DIR once every task succeeds.
#
# Usage (from the repository root on an Alpine login node):
#   bash hpc/submit_null_bundle.sh
#   DATA_DIR=/projects/$USER/repositories/multi-dwpc/data bash hpc/submit_null_bundle.sh
#

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_DIR="${DATA_DIR:-$REPO_ROOT/data}"
WORK_DIR="${WORK_DIR:-/scratch/alpine/$USER/multi-dwpc/null_bundle}"
CACHE_DIR="${CACHE_DIR:-$WORK_DIR/dwpc_cache}"
PARTS_DIR="${PARTS_DIR:-$WORK_DIR/parts}"
BUNDLE_DIR="${BUNDLE_DIR:-$DATA_DIR/null_bundle}"
ACCOUNT="${ACCOUNT:-amc-general}"
PARTITION="${PARTITION:-acpu}"
QOS="${QOS:-cpu-normal}"
CONDA_SH="${CONDA_SH:-/curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh}"
LOG_DIR="$REPO_ROOT/hpc/logs/null_bundle"
N_METAPATHS=52

mkdir -p "$LOG_DIR" "$CACHE_DIR" "$PARTS_DIR"
ACTIVATE="source $CONDA_SH && conda activate multi_dwpc && cd $REPO_ROOT"

BUILD_JOB=$(sbatch --parsable \
    --job-name=null-bundle-build --account="$ACCOUNT" --partition="$PARTITION" --qos="$QOS" \
    --cpus-per-task=4 --mem=32G --time=04:00:00 --array=0-$((N_METAPATHS - 1)) \
    --output="$LOG_DIR/build_%A_%a.out" --error="$LOG_DIR/build_%A_%a.err" \
    --wrap="$ACTIVATE && python scripts/build_null_bundle.py --index \$SLURM_ARRAY_TASK_ID \
        --data-dir $DATA_DIR --cache-dir $CACHE_DIR --parts-dir $PARTS_DIR")

FINALIZE_JOB=$(sbatch --parsable --dependency=afterok:"$BUILD_JOB" \
    --job-name=null-bundle-finalize --account="$ACCOUNT" --partition="$PARTITION" --qos="$QOS" \
    --cpus-per-task=1 --mem=32G --time=02:00:00 \
    --output="$LOG_DIR/finalize_%j.out" --error="$LOG_DIR/finalize_%j.err" \
    --wrap="$ACTIVATE && python scripts/finalize_null_bundle.py --parts-dir $PARTS_DIR \
        --bundle-dir $BUNDLE_DIR --data-dir $DATA_DIR --expected-parts $N_METAPATHS")

echo "build array: $BUILD_JOB  finalize: $FINALIZE_JOB"
echo "parts: $PARTS_DIR  bundle: $BUNDLE_DIR  logs: $LOG_DIR"
