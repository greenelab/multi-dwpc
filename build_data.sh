#!/usr/bin/env bash

# this script originates from the instructions at
# https://github.com/greenelab/multi-dwpc#prepare-the-data

set -euo pipefail

time (
    # Run individual data-prep steps (no top-level pipeline runner; chain
    # the steps you need)
    uv run poe load-data
    uv run poe filter-change
    uv run poe go-hierarchy-analysis
    uv run poe filter-jaccard

    # Generate year null datasets
    uv run poe gen-permutation
    uv run poe gen-random

    # Compute DWPC matrices for the prepared data
    uv run poe compute-dwpc-direct
)
