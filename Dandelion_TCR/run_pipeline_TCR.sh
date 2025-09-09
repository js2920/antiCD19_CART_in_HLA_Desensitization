#!/usr/bin/env bash
# scripts/run_pipeline_TCR.sh
# Reproducible TCR pipeline wrapper around the Dandelion Singularit container.
# Actions: preprocess (optional), merge, integrate, all.

set -euo pipefail

# -------- Defaults (override via env or --flags) -------------------------------
ROOT="${ROOT:-$(pwd)/dandelion_inputs_TCR}"   # archive root
OUTDIR="${OUTDIR:-$ROOT/analysis_ready}"
META="${META:-auto}"                          # Dandelion meta.csv; 'auto' will resolve
GEX_H5AD="${GEX_H5AD:-auto}"                  # processed_adata_NTX_all_samples.h5ad; 'auto' resolves
SIF="${SIF:-auto}"                            # sc-dandelion_latest.sif; 'auto' resolves
ACTIONS="${ACTIONS:-all}"                     # csv list: preprocess,merge,integrate,all
FILE_PREFIX="${FILE_PREFIX:-all}"             # 'all' (default in Dandelion) or 'filtered'
MERGED_TSV="${MERGED_TSV:-$OUTDIR/my_merged_TCR_contigs.tsv}"
PER_CELL_TSV="${PER_CELL_TSV:-$OUTDIR/per_cell_TCR_table.csv}"
# ------------------------------------------------------------------------------

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
[[ "$SIF" == "auto" ]] && die "Could not find sc-dandelion_latest.sif. See README."

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

# Container helpers (bind ROOT into container)
S_RUN() { "$RUNTIME" run  -B "$ROOT":"$ROOT" "$SIF" "$@"; }
S_EXEC(){ "$RUNTIME" exec -B "$ROOT":"$ROOT" "$SIF" "$@"; }

preprocess_step() {
  echo "==> [preprocess] Checking for raw 10x contig FASTA folders…"
  mapfile -t candidates < <(find "$ROOT" -maxdepth 3 -type f -name "${FILE_PREFIX}_contig.fasta" -printf '%h\n' | sort -u || true)

  if (( ${#candidates[@]} )); then
    echo "==> [preprocess] Found ${#candidates[@]} candidate sample folder(s). Running dandelion-preprocess…"
    if [[ -n "${META}" ]]; then
      S_RUN dandelion-preprocess --meta "$META" --file_prefix "$FILE_PREFIX"
    else
      S_RUN dandelion-preprocess --file_prefix "$FILE_PREFIX"
    fi
  else
    echo "==> [preprocess] No ${FILE_PREFIX}_contig.fasta found. Skipping (likely already preprocessed)."
  fi
}

merge_step() {
  echo "==> [merge] Gathering per-sample AIRR TSVs (TCR)…"
  mapfile -t tsvs < <( \
    find "$ROOT" -type f -name "all_contig_dandelion.tsv" -path "*/dandelion/*" | sort )
  if (( ${#tsvs[@]} == 0 )); then
    # Fallback to base_inputs/contigs if provided
    mapfile -t tsvs < <(find "$ROOT/base_inputs/contigs" -type f -name "*.tsv" | sort || true)
  fi
  (( ${#tsvs[@]} )) || die "No Dandelion AIRR TSVs found to merge."

  echo "==> [merge] Will merge ${#tsvs[@]} file(s) into: $MERGED_TSV"
  tmp="$MERGED_TSV.tmp"
  : > "$tmp"
  first=1
  for f in "${tsvs[@]}"; do
    if (( first )); then
      cat "$f" >> "$tmp"
      first=0
    else
      tail -n +2 "$f" >> "$tmp"  # append without header
    fi
  done
  mv "$tmp" "$MERGED_TSV"
  echo "==> [merge] Done."
}

integrate_step() {
  [[ -n "$GEX_H5AD" ]] || die "GEX .h5ad not found; set GEX_H5AD env or flag."
  [[ -f "$MERGED_TSV" ]] || die "Merged contigs TSV not found at $MERGED_TSV"

  echo "==> [integrate] Integrating TCR (TRA/TRB) with AnnData; writing per-cell table + augmented .h5ad…"
  S_EXEC python - <<'PY'
import os, re
import pandas as pd
import numpy as np
import scanpy as sc
import dandelion as ddl

ROOT     = os.environ.get("ROOT_DIR", "")
OUTDIR   = os.environ.get("OUTDIR", "")
MERGED   = os.environ.get("MERGED_TSV", "")
GEX      = os.environ.get("GEX_H5AD", "")
PER_CELL = os.environ.get("PER_CELL_TSV", "")

os.makedirs(OUTDIR, exist_ok=True)

# ---- Load GEX and contigs (AIRR) --------------------------------------------
adata = sc.read(GEX)

# Read AIRR directly with pandas for robust column handling (then also via ddl)
with open(MERGED) as fh:
    header = fh.readline().rstrip("\n").split("\t")
cols = [c for c in [
    "cell_id","barcode","locus","chain","productive",
    "umi_count","consensus_count","read_count",
    "v_call","j_call","junction_aa","cdr3_aa","cdr3"
] if c in header]
df = pd.read_csv(MERGED, sep="\t", usecols=cols, low_memory=True)

# Harmonize columns
df.rename(columns={
    "cell_id":"barcode",
    "locus":"chain",
    "junction_aa":"cdr3",
    "cdr3_aa":"cdr3",
}, inplace=True)

# Productive flag -> bool
if "productive" in df.columns:
    df["productive"] = (df["productive"].astype(str).str.strip().str.upper()
                        .map({"TRUE":True,"T":True,"1":True,"FALSE":False,"F":False,"0":False})
                        .fillna(False))
else:
    df["productive"] = True

# Keep TRA/TRB productive
df["chain"] = df.get("chain","").astype(str).str.upper().str.extract(r"(TRA|TRB)", expand=False)
df = df[df["chain"].isin(["TRA","TRB"]) & df["productive"]].copy()

# ---- Drop multiplets (cells with >1 of TRA or >1 of TRB) ---------------------
multi = set()
for ch in ("TRA","TRB"):
    vc = df[df["chain"]==ch]["barcode"].value_counts()
    multi.update(vc[vc>1].index)
if multi:
    df = df[~df["barcode"].isin(multi)].copy()
    adata = adata[~adata.obs_names.isin(multi)].copy()

# ---- Keep one contig per (barcode, chain) using best available metric --------
metric = next((m for m in ("umi_count","consensus_count","read_count") if m in df.columns), None)
if metric is None:
    df["__metric"] = 1
    metric = "__metric"
df[metric] = pd.to_numeric(df[metric], errors="coerce").fillna(0)
idx = (df.sort_values(metric, ascending=False)
         .groupby(["barcode","chain"], sort=False).head(1).index)
df = df.loc[idx].copy().reset_index(drop=True)

# ---- Per-cell wide table -----------------------------------------------------
keep = [c for c in ("cdr3","v_call","j_call") if c in df.columns]
base = df[["barcode","chain",*keep]].set_index("barcode")

def pivot(chain):
    sub = base[base["chain"]==chain]
    if sub.empty:
        return pd.DataFrame(columns=[f"{chain}_{c}" for c in keep], dtype=str)
    sub = sub[keep].astype(str)
    sub.columns = [f"{chain}_{c}" for c in keep]
    return sub

per_cell = pivot("TRA").join(pivot("TRB"), how="outer")
per_cell.to_csv(PER_CELL)

# ---- Join into AnnData; simple boolean convenience ---------------------------
adata.obs = adata.obs.join(per_cell, how="left")
adata.obs["has_TCR"] = adata.obs.get("TRA_cdr3","").astype(str).ne("") | \
                       adata.obs.get("TRB_cdr3","").astype(str).ne("")

# ---- Optional: clone calling and transfer via Dandelion ----------------------
vdj = ddl.read_airr(MERGED)
# Align & basic QC (will silently drop barcodes that don't match)
vdj, adata = ddl.pp.check_contigs(vdj, adata)
ddl.tl.find_clones(vdj, key="junction_aa")
ddl.tl.transfer(adata, vdj)

# ---- Save outputs ------------------------------------------------------------
vdj_out   = os.path.join(OUTDIR, "vdj_TCR.h5ddl")
adata_out = os.path.join(OUTDIR, "adata_with_TCR_integration.h5ad")
vdj.write(vdj_out)
adata.write(adata_out)

print(f"[integrate] Wrote per-cell table      → {PER_CELL}")
print(f"[integrate] Wrote vdj_TCR object      → {vdj_out}")
print(f"[integrate] Wrote augmented AnnData   → {adata_out}")
PY
  echo "==> [integrate] Done."
}

usage() {
  cat <<EOF
Usage:
  ROOT=./dandelion_inputs_TCR \\
  SIF=/path/to/sc-dandelion_latest.sif \\
  META=auto GEX_H5AD=auto OUTDIR=./dandelion_inputs_TCR/analysis_ready \\
  FILE_PREFIX=all ACTIONS=all \\
  bash scripts/run_pipeline_TCR.sh

Environment variables / flags:
  ROOT         : Archive root (default: ./dandelion_inputs_TCR)
  SIF          : Path to sc-dandelion_latest.sif (default: auto-detect)
  META         : Path to meta.csv (default: auto-detect; optional)
  GEX_H5AD     : Path to processed_adata_*.h5ad (default: auto-detect)
  OUTDIR       : Output dir (default: ROOT/analysis_ready)
  FILE_PREFIX  : 'all' or 'filtered' (default: all)
  MERGED_TSV   : Output merged AIRR TSV (default: OUTDIR/my_merged_TCR_contigs.tsv)
  PER_CELL_TSV : Per-cell TCR table CSV (default: OUTDIR/per_cell_TCR_table.csv)
  ACTIONS      : comma list or 'all' (preprocess,merge,integrate)
EOF
}

# Minimal flag parser (most config via env vars)
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
    --per-cell=*) PER_CELL_TSV="${arg#*=}";;
    *) echo "[WARN] Unknown arg: $arg";;
  esac
done

export ROOT_DIR="$ROOT" OUTDIR="$OUTDIR" MERGED_TSV="$MERGED_TSV" GEX_H5AD="$GEX_H5AD" PER_CELL_TSV="$PER_CELL_TSV"

IFS=',' read -r -a acts <<< "$ACTIONS"
run_all=true; [[ "$ACTIONS" == "all" ]] || run_all=false

$run_all && preprocess_step || true
$run_all && merge_step || true
$run_all && integrate_step || true

if ! $run_all; then
  for a in "${acts[@]}"; do
    case "$a" in
      preprocess) preprocess_step;;
      merge) merge_step;;
      integrate) integrate_step;;
      *) echo "[WARN] Unknown action: $a";;
    esac
  done
fi

echo "==> TCR pipeline finished."
