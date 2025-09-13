#!/usr/bin/env bash
# ------------------------------------------------------------------------------
# Singularity-based Dandelion preprocessing (+ optional Immcantation shell)
# Path-agnostic template for GEO: GSE307464
# ------------------------------------------------------------------------------
# Usage:
#   SIF=./sc-dandelion_latest.sif META=meta.csv BIND=$PWD ./scripts/preprocess_and_env.sh
# ------------------------------------------------------------------------------

set -euo pipefail

SIF="${SIF:-./sc-dandelion_latest.sif}"
META="${META:-meta.csv}"
BIND="${BIND:-$PWD}"
DANDELION_CMD="${DANDELION_CMD:-dandelion-preprocess}"

echo "[INFO] Running Dandelion preprocess inside Singularity"
singularity run -B "$BIND" "$SIF" \
  "$DANDELION_CMD" --filter_to_high_confidence --meta "$META"

# Optional: interactive Immcantation shell (no installation steps baked in)
if [[ "${IMMCANTATION:-0}" == "1" ]]; then
  echo "[INFO] Launching Immcantation/suite:4.5.0 (interactive shell)"
  docker run --rm -it -v "$BIND:/work" -w /work immcantation/suite:4.5.0 /bin/bash
fi
