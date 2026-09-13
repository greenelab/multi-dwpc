#!/bin/bash
#
# Validate the null bundle on Alpine: one array task per metapath compares the
# laptop path with the matrix-based adapter, then a summarize job writes
# tables, figures and the verdict into docs/tasks/laptop-null-bundle/.
#
# Usage (from the repository root on an Alpine login node):
#   DATA_DIR=/projects/$USER/repositories/multi-dwpc/data bash hpc/submit_validate_null_bundle.sh
#   AFTER_JOB=<finalize job id> ... bash hpc/submit_validate_null_bundle.sh   # wait for the build
#

set -euo pipefail

REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
DATA_DIR="${DATA_DIR:-$REPO_ROOT/data}"
WORK_DIR="${WORK_DIR:-/scratch/alpine/$USER/multi-dwpc/null_bundle}"
CACHE_DIR="${CACHE_DIR:-$WORK_DIR/dwpc_cache}"
OUT_DIR="${OUT_DIR:-$WORK_DIR/validate}"
BUNDLE_DIR="${BUNDLE_DIR:-$DATA_DIR/null_bundle}"
ACCOUNT="${ACCOUNT:-amc-general}"
PARTITION="${PARTITION:-acpu}"
QOS="${QOS:-cpu-normal}"
CONDA_SH="${CONDA_SH:-/curc/sw/anaconda3/2023.09/etc/profile.d/conda.sh}"
LOG_DIR="$REPO_ROOT/hpc/logs/null_bundle"
N_METAPATHS=52
DEPENDENCY=()
if [[ -n "${AFTER_JOB:-}" ]]; then
    DEPENDENCY=(--dependency=afterok:"$AFTER_JOB")
fi

mkdir -p "$LOG_DIR" "$OUT_DIR"
ACTIVATE="source $CONDA_SH && conda activate multi_dwpc && cd $REPO_ROOT"

VALIDATE_JOB=$(sbatch --parsable "${DEPENDENCY[@]}" \
    --job-name=null-bundle-validate --account="$ACCOUNT" --partition="$PARTITION" --qos="$QOS" \
    --cpus-per-task=2 --mem=24G --time=04:00:00 --array=0-$((N_METAPATHS - 1)) \
    --output="$LOG_DIR/validate_%A_%a.out" --error="$LOG_DIR/validate_%A_%a.err" \
    --wrap="$ACTIVATE && python scripts/validate_null_bundle.py --index \$SLURM_ARRAY_TASK_ID \
        --bundle-dir $BUNDLE_DIR --data-dir $DATA_DIR --cache-dir $CACHE_DIR --out-dir $OUT_DIR")

SUMMARY_JOB=$(sbatch --parsable --dependency=afterok:"$VALIDATE_JOB" \
    --job-name=null-bundle-validate-summary --account="$ACCOUNT" --partition="$PARTITION" --qos="$QOS" \
    --cpus-per-task=1 --mem=16G --time=01:00:00 \
    --output="$LOG_DIR/validate_summary_%j.out" --error="$LOG_DIR/validate_summary_%j.err" \
    --wrap="$ACTIVATE && python scripts/summarize_null_bundle_validation.py --in-dir $OUT_DIR \
        --task-dir docs/tasks/laptop-null-bundle")

echo "validate array: $VALIDATE_JOB  summary: $SUMMARY_JOB"
