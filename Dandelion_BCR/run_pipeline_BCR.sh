#!/usr/bin/env bash
# scripts/run_pipeline.sh
# Reproducible BCR pipeline wrapper around the Dandelion Singularity/Apptainer container.
# Supports: preprocess (optional), merge, integrate, changeo (optional), all.
# Tested with dandelion sc-dandelion_latest.sif (see README for how to pull/locate).

set -euo pipefail

# -------- Defaults (override via flags) ---------------------------------------
ROOT="${ROOT:-$(pwd)/dandelion_inputs}"  # path to your extracted archive root
OUTDIR="${OUTDIR:-$ROOT/analysis_ready}"
META="${META:-auto}"                      # CSV describing samples; 'auto' will resolve
GEX_H5AD="${GEX_H5AD:-auto}"             # gene expression .h5ad; 'auto' will resolve
SIF="${SIF:-auto}"                        # sc-dandelion_latest.sif; 'auto' will resolve
ACTIONS="${ACTIONS:-all}"                 # csv list: preprocess,merge,integrate,changeo,all
FILE_PREFIX="${FILE_PREFIX:-all}"         # 'all' (default in dandelion docs) or 'filtered'
MERGED_TSV="${MERGED_TSV:-$OUTDIR/my_merged_contigs.tsv}"
# -----------------------------------------------------------------------------

have_cmd () { command -v "$1" >/dev/null 2>&1; }
die() { echo "[ERROR] $*" >&2; exit 1; }

detect_runtime() {
  if have_cmd apptainer; then echo "apptainer"; return
  elif have_cmd singularity; then echo "singularity"; return
  else die "Neither 'apptainer' nor 'singularity' found on PATH."
  fi
}

RUNTIME="$(detect_runtime)"

# Resolve defaults
if [[ "$SIF" == "auto" ]]; then
  for c in \
    "$ROOT/base_inputs/containers/sc-dandelion_latest.sif" \
    "$ROOT/analysis_ready/sc-dandelion_latest.sif" \
    "./sc-dandelion_latest.sif"
  do
    [[ -f "$c" ]] && SIF="$c" && break
  done
fi
[[ "$SIF" == "auto" ]] && die "Could not find sc-dandelion_latest.sif. See README to pull it."

if [[ "$META" == "auto" ]]; then
  if [[ -f "$ROOT/base_inputs/meta/meta.csv" ]]; then META="$ROOT/base_inputs/meta/meta.csv"
  elif [[ -f "$ROOT/analysis_ready/meta.csv" ]]; then META="$ROOT/analysis_ready/meta.csv"
  else META=""; fi
fi

if [[ "$GEX_H5AD" == "auto" ]]; then
  if [[ -f "$ROOT/base_inputs/adata/processed_adata_NTX_all_samples.h5ad" ]]; then
    GEX_H5AD="$ROOT/base_inputs/adata/processed_adata_NTX_all_samples.h5ad"
  elif [[ -f "$ROOT/analysis_ready/processed_adata_NTX_all_samples.h5ad" ]]; then
    GEX_H5AD="$ROOT/analysis_ready/processed_adata_NTX_all_samples.h5ad"
  else
    GEX_H5AD=""
  fi
fi

mkdir -p "$OUTDIR"

# Helpers to run inside container with a consistent bind
S_RUN() { "$RUNTIME" run -B "$ROOT":"$ROOT" "$SIF" "$@"; }
S_EXEC() { "$RUNTIME" exec -B "$ROOT":"$ROOT" "$SIF" "$@"; }

preprocess_step() {
  echo "==> [preprocess] Checking for Cell Ranger input folders…"
  # Look for per-sample folders that might contain 'all_contig.fasta'
  mapfile -t candidates < <(find "$ROOT" -maxdepth 3 -type f -name "${FILE_PREFIX}_contig.fasta" -printf '%h\n' | sort -u || true)

  if (( ${#candidates[@]} )); then
    echo "==> [preprocess] Found $((${#candidates[@]})) sample folder(s). Running dandelion-preprocess…"
    # If META is set, pass --meta to group prefixes/individuals as recommended in docs
    if [[ -n "${META}" ]]; then
      S_RUN dandelion-preprocess --meta "$META" --file_prefix "$FILE_PREFIX"
    else
      S_RUN dandelion-preprocess --file_prefix "$FILE_PREFIX"
    fi
  else
    echo "==> [preprocess] No ${FILE_PREFIX}_contig.fasta found. Skipping (probably already preprocessed)."
  fi
}

merge_step() {
  echo "==> [merge] Gathering per-sample AIRR TSVs…"
  # Prefer dandelion outputs in typical locations
  mapfile -t tsvs < <( \
    find "$ROOT" -type f \( -name "*_contig_dandelion.tsv" -o -name "all_contig_dandelion.tsv" \) \
      -path "*/dandelion/*" | sort )
  if (( ${#tsvs[@]} == 0 )); then
    # Fallback to provided base_inputs contigs
    mapfile -t tsvs < <(find "$ROOT/base_inputs/contigs" -type f -name "*.tsv" | sort || true)
  fi
  (( ${#tsvs[@]} )) || die "No dandelion AIRR TSVs found to merge."

  echo "==> [merge] Will merge ${#tsvs[@]} file(s) into: $MERGED_TSV"
  tmp="$MERGED_TSV.tmp"
  : > "$tmp"
  first=1
  for f in "${tsvs[@]}"; do
    if (( first )); then
      cat "$f" >> "$tmp"
      first=0
    else
      # append without header
      tail -n +2 "$f" >> "$tmp"
    fi
  done
  mv "$tmp" "$MERGED_TSV"
  echo "==> [merge] Done."
}

integrate_step() {
  [[ -n "$GEX_H5AD" ]] || die "GEX .h5ad not found; set GEX_H5AD env or flag."
  [[ -f "$MERGED_TSV" ]] || die "Merged contigs TSV not found at $MERGED_TSV"

  echo "==> [integrate] Linking VDJ with GEX (writes .h5ddl and integrated .h5ad)…"
  S_EXEC python - <<'PY'
import os
import dandelion as ddl
import scanpy as sc

ROOT = os.environ.get("ROOT_DIR", "")
OUTDIR = os.environ.get("OUTDIR", "")
MERGED_TSV = os.environ.get("MERGED_TSV", "")
GEX_H5AD = os.environ.get("GEX_H5AD", "")

if not OUTDIR:
    OUTDIR = os.path.join(ROOT, "analysis_ready")

# Read GEX and VDJ (AIRR TSV from dandelion)
adata = sc.read(GEX_H5AD)
vdj = ddl.read_airr(MERGED_TSV)  # dandelion's *_contig_dandelion.tsv is AIRR-compliant

# Quality checks & synchronization (filters ambiguous contigs, aligns with GEX)
vdj, adata = ddl.pp.check_contigs(vdj, adata)  # recommended entry point

# Example: call clones (BCR default = V/J identical & <=15% CDR3 mismatches)
ddl.tl.find_clones(vdj)

# Transfer analysis info back onto .obs / .obsp for convenience
ddl.tl.transfer(adata, vdj)

# Persist both objects
vdj_out = os.path.join(OUTDIR, "vdj.h5ddl")
adata_out = os.path.join(OUTDIR, "adata_with_bcr_integration.h5ad")
vdj.write(vdj_out)
adata.write(adata_out)

print(f"[integrate] Wrote: {vdj_out}")
print(f"[integrate] Wrote: {adata_out}")
PY
  echo "==> [integrate] Done."
}

changeo_step() {
  # Optional: run Change-O clonotype calling inside the dandelion container
  local h5ddl="$OUTDIR/vdj.h5ddl"
  [[ -f "$h5ddl" ]] || die "Expected $h5ddl from integrate_step."
  echo "==> [changeo] Running changeo-clonotypes via SIF…"
  S_RUN changeo-clonotypes --h5ddl "$h5ddl"
  echo "==> [changeo] Wrote SHazaM threshold plot and *_changeo.h5ddl (see README)."
}

usage() {
  cat <<EOF
Usage:
  ROOT=./dandelion_inputs \\
  SIF=/path/to/sc-dandelion_latest.sif \\
  META=auto GEX_H5AD=auto OUTDIR=./dandelion_inputs/analysis_ready \\
  FILE_PREFIX=all ACTIONS=all \\
  bash scripts/run_pipeline.sh

Environment variables / flags (export or prefix on the command):
  ROOT         : Archive root (default: ./dandelion_inputs)
  SIF          : Path to sc-dandelion_latest.sif (default: auto-detect)
  META         : Path to meta.csv (default: auto-detect; optional)
  GEX_H5AD     : Path to processed_adata_*.h5ad (default: auto-detect)
  OUTDIR       : Output dir (default: ROOT/analysis_ready)
  FILE_PREFIX  : 'all' or 'filtered' (default: all)
  MERGED_TSV   : Output merged AIRR TSV (default: OUTDIR/my_merged_contigs.tsv)
  ACTIONS      : comma list or 'all' (preprocess,merge,integrate,changeo)
EOF
}

# Parse minimal flags (optional; most use env vars)
for arg in "$@"; do
  case "$arg" in
    -h|--help) usage; exit 0;;
    --root=*) ROOT="${arg#*=}";;
    --sif=*) SIF="${arg#*=}";;
    --meta=*) META="${arg#*=}";;
    --gex=*) GEX_H5AD="${arg#*=}";;
    --outdir=*) OUTDIR="${arg#*=}";;
    --file-prefix=*) FILE_PREFIX="${arg#*=}";;
    --actions=*) ACTIONS="${arg#*=}";;
    --merged-tsv=*) MERGED_TSV="${arg#*=}";;
    *) echo "[WARN] Unknown arg: $arg";;
  esac
done

export ROOT_DIR="$ROOT" OUTDIR="$OUTDIR" MERGED_TSV="$MERGED_TSV" GEX_H5AD="$GEX_H5AD"

IFS=',' read -r -a acts <<< "$ACTIONS"
run_all=true; [[ "$ACTIONS" == "all" ]] || run_all=false

$run_all && preprocess_step || true
$run_all && merge_step || true
$run_all && integrate_step || true
$run_all && changeo_step || true

if ! $run_all; then
  for a in "${acts[@]}"; do
    case "$a" in
      preprocess) preprocess_step;;
      merge) merge_step;;
      integrate) integrate_step;;
      changeo) changeo_step;;
      *) echo "[WARN] Unknown action: $a";;
    esac
  done
fi

echo "==> Pipeline finished."
