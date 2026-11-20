#!/usr/bin/env bash
# PAACS Mesh
# Project author: Paul LeDewhurst
# Copyright PAACS Software Ltd.
# SPDX-License-Identifier: MPL-2.0

set -euo pipefail

BENCHMARK=${1:-./paacs-v09-50m-bench}
THREADS=${PAACS_V09_THREADS:-0}
CELLS=${PAACS_V09_CELLS:-128}
RUNS=${PAACS_V09_RUNS:-11}
THRESHOLD=${PAACS_V09_THRESHOLD_MCELLS:-50}

exec "$BENCHMARK" \
  --cells "$CELLS" \
  --threads "$THREADS" \
  --runs "$RUNS" \
  --threshold "$THRESHOLD"
