#!/usr/bin/env bash
#
# run_from_zip.sh
#
# Convenience wrapper to (optionally) download CAR_extraction_inputs.zip,
# unzip it, and run the pipeline using scripts/car_extraction_from_bundle.sh.
#
# Usage:
#   scripts/run_from_zip.sh --zip <path_or_URL> [--output <DIR>] [--threads N]
#
# Example (local zip):
#   scripts/run_from_zip.sh --zip ./archives/CAR_extraction_inputs.zip
#
# Example (Zenodo URL):
#   scripts/run_from_zip.sh --zip https://zenodo.org/record/.../files/CAR_extraction_inputs.zip?download=1
#
# -----------------------------------------------------------------------------

set -euo pipefail

ZIP_ARG=""
OUTPUT_DIR="${PWD}/analysis_output"
THREADS="${THREADS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 16)}"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --zip) ZIP_ARG="$2"; shift 2;;
    --output) OUTPUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help)
      echo "Usage: $0 --zip <path_or_URL> [--output DIR] [--threads N]"
      exit 0;;
    *) echo "Unknown arg: $1"; exit 1;;
  esac
done

if [[ -z "${ZIP_ARG}" ]]; then
  echo "ERROR: --zip is required"; exit 1
fi

WORKDIR="${PWD}/_bundle_work"
mkdir -p "$WORKDIR"

ZIP_LOCAL="${WORKDIR}/CAR_extraction_inputs.zip"
if [[ "$ZIP_ARG" =~ ^https?:// ]]; then
  echo "Downloading bundle from URL…"
  if command -v curl >/dev/null 2>&1; then
    curl -L "$ZIP_ARG" -o "$ZIP_LOCAL"
  elif command -v wget >/dev/null 2>&1; then
    wget -O "$ZIP_LOCAL" "$ZIP_ARG"
  else
    echo "ERROR: neither curl nor wget available to download $ZIP_ARG"; exit 1
  fi
else
  ZIP_LOCAL="$ZIP_ARG"
fi

[[ -f "$ZIP_LOCAL" ]] || { echo "ERROR: ZIP not found at $ZIP_LOCAL"; exit 1; }

echo "Unzipping bundle…"
UNZIP_DIR="${WORKDIR}/unzipped"
rm -rf "$UNZIP_DIR"
mkdir -p "$UNZIP_DIR"
unzip -q "$ZIP_LOCAL" -d "$UNZIP_DIR"

BUNDLE_DIR="${UNZIP_DIR}/CAR_extraction_inputs"
[[ -d "$BUNDLE_DIR" ]] || { echo "ERROR: Expected folder not found: $BUNDLE_DIR"; exit 1; }

# Run the main pipeline
SCRIPTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
"${SCRIPTDIR}/car_extraction_from_bundle.sh" \
  --bundle-dir "$BUNDLE_DIR" \
  --output "$OUTPUT_DIR" \
  --threads "$THREADS"
