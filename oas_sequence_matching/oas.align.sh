#!/usr/bin/env bash
#==============================================================================
# align_cdr3_against_oas.sh — align query CDR‑H3 clonotypes to a prebuilt
# OAS IGHG CDR3 MMseqs2 database and generate (optional) UMAP figures.
#
# Steps:
#   1) Build query FASTA from clonotypes CSV/TSV
#   2) MMseqs2 search (>=95% ID / >=95% coverage by default) + global rescore
#   3) Merge hits with full query CDR3
#   4) (optional) Add CellTypist labels and draw UMAPs
#
# Usage (typical):
#   scripts/align_cdr3_against_oas.sh \
#     --db-dir     /home/jensf/data/oas_ighg_cdr3 \
#     --clonotypes dandelion_inputs/analysis_ready/my_merged_contigs.tsv \
#     --adata      totalvi_celltypist_inputs/processed_adata/processed_adata_NTX_all_samples.h5ad \
#     --outdir     review_align_outputs \
#     --threads    16
#
# Required:
#   --db-dir       directory containing: cdr3_ighg.faa  and  cdr3_ighg_db*
#   --clonotypes   CSV/TSV with cell IDs + CDR‑H3 AA column (many common names auto-detected)
#
# Optional:
#   --adata        AnnData .h5ad with X_umap or X_totalVI_umap + CellTypist label column
#   --mmseqs       path to mmseqs binary (auto-discovered if not provided)
#   --min-id 0.95  minimum sequence identity (default 0.95)
#   --min-cov 0.95 minimum coverage (default 0.95)
#   --cov-mode 0   coverage mode (default 0, see MMseqs2 docs)
#   --rescore-mode 3  global SW-like rescoring (default 3)
#   --pre "NTX_1_,NTX_2_,NTX_3_"    sample prefixes to tag as "pre"
#   --post "NTX_4_,NTX_5_"          sample prefixes to tag as "post"
#   --seq-map path/to/seq_to_category.tsv  (2 columns: CDR3\tCategory)
#   --outdir ./cdr3_align_out
#   --threads 8
#
# Environment:
#   Use env/cdr3-align.yml for a reproducible setup (see below).
#==============================================================================

set -Eeuo pipefail
LC_ALL=C

# --- defaults -----------------------------------------------------------------
OUTDIR="$PWD/cdr3_align_out"
THREADS=8
MIN_ID=0.95
MIN_COV=0.95
COV_MODE=0
RESCORE_MODE=3
DB_DIR=""
CLONOTYPES=""
ADATA=""
SEQ_MAP=""
PRE_PREFIXES="NTX_1_,NTX_2_,NTX_3_"
POST_PREFIXES="NTX_4_,NTX_5_"
MMSEQS_BIN=""

# --- helpers ------------------------------------------------------------------
die(){ printf "[ERROR] %s\n" "$*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing tool: $1"; }
ap(){ python3 - "$@"; }  # run inline python

resolve_mmseqs(){
  local bin="$1"
  if [[ -n "$bin" ]]; then command -v "$bin" >/dev/null 2>&1 && { echo "$bin"; return; }
                          [[ -x "$bin" ]] && { echo "$bin"; return; }
  fi
  if command -v mmseqs >/dev/null 2>&1; then echo "mmseqs"; return; fi
  if command -v conda >/dev/null 2>&1; then
    local base; base="$(conda info --base 2>/dev/null || true)"
    [[ -n "$base" && -x "$base/bin/mmseqs" ]] && { echo "$base/bin/mmseqs"; return; }
    while IFS= read -r -d '' p; do echo "$p"; return; done \
      < <(find "$base/envs" -type f -name mmseqs -perm -111 -print0 2>/dev/null || true)
  fi
  return 1
}

usage(){
cat <<EOF
Usage:
  $(basename "$0") --db-dir <dir> --clonotypes <csv/tsv> [--adata <adata.h5ad>] [options]

Options:
  --db-dir DIR         Directory with cdr3_ighg.faa and cdr3_ighg_db*
  --clonotypes PATH    CSV/TSV containing cell IDs + CDR3 AA column
  --adata PATH         AnnData .h5ad (for optional UMAP figures)
  --seq-map PATH       2-col TSV/CSV: <cdr3aa> <category>
  --pre STR            Comma-separated sample prefixes tagged as "pre"  [$PRE_PREFIXES]
  --post STR           Comma-separated sample prefixes tagged as "post" [$POST_PREFIXES]
  --mmseqs PATH        Path to mmseqs binary (auto if omitted)
  --min-id FLOAT       Min sequence identity (default $MIN_ID)
  --min-cov FLOAT      Min coverage (default $MIN_COV)
  --cov-mode INT       Coverage mode (default $COV_MODE)
  --rescore-mode INT   Rescore mode   (default $RESCORE_MODE)
  --threads INT        Threads for search/convertalis (default $THREADS)
  --outdir DIR         Output directory (default $OUTDIR)
  -h, --help           This help
EOF
}

# --- arg parsing --------------------------------------------------------------
[[ $# -eq 0 ]] && { usage; exit 0; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --db-dir)       DB_DIR="$2"; shift 2;;
    --clonotypes)   CLONOTYPES="$2"; shift 2;;
    --adata)        ADATA="$2"; shift 2;;
    --seq-map)      SEQ_MAP="$2"; shift 2;;
    --pre)          PRE_PREFIXES="$2"; shift 2;;
    --post)         POST_PREFIXES="$2"; shift 2;;
    --mmseqs)       MMSEQS_BIN="$2"; shift 2;;
    --min-id)       MIN_ID="$2"; shift 2;;
    --min-cov)      MIN_COV="$2"; shift 2;;
    --cov-mode)     COV_MODE="$2"; shift 2;;
    --rescore-mode) RESCORE_MODE="$2"; shift 2;;
    --threads)      THREADS="$2"; shift 2;;
    --outdir)       OUTDIR="$2"; shift 2;;
    -h|--help)      usage; exit 0;;
    *) die "Unknown option: $1";;
  esac
done

# --- checks -------------------------------------------------------------------
[[ -n "$DB_DIR"     ]] || die "--db-dir is required"
[[ -n "$CLONOTYPES" ]] || die "--clonotypes is required"
[[ -d "$DB_DIR"     ]] || die "DB dir not found: $DB_DIR"
[[ -r "$CLONOTYPES" ]] || die "Cannot read clonotypes: $CLONOTYPES"

need python3
MMSEQS_BIN="$(resolve_mmseqs "$MMSEQS_BIN")" || die "mmseqs not found; set --mmseqs /abs/path/to/mmseqs"

mkdir -p "$OUTDIR"
LOG="$OUTDIR/LOG.txt"
WORK="$OUTDIR/work"; mkdir -p "$WORK"
echo "# CDR3 alignment — $(date)" > "$LOG"

DB_FASTA="$DB_DIR/cdr3_ighg.faa"
DB_PREFIX="$DB_DIR/cdr3_ighg_db"
[[ -r "$DB_FASTA" ]] || die "Missing FASTA in DB dir: $DB_FASTA"
[[ -e "${DB_PREFIX}.index" || -e "${DB_PREFIX}.dbtype" || -e "${DB_PREFIX}" ]] || die "MMseqs DB not found under: $DB_PREFIX*"

echo "DB dir     : $DB_DIR"     | tee -a "$LOG"
echo "Clonotypes : $CLONOTYPES" | tee -a "$LOG"
[[ -n "$ADATA" ]] && echo "AnnData    : $ADATA" | tee -a "$LOG"
echo "Outdir     : $OUTDIR"     | tee -a "$LOG"
echo "mmseqs     : $("$MMSEQS_BIN" --version 2>/dev/null | tr '\n' ' ')" | tee -a "$LOG"

# --- STEP 1: query FASTA from clonotypes -------------------------------------
PAIR_CSV="$WORK/cdr3_pairs.csv"
QUERY_FAA="$WORK/query.faa"
echo "▸ Building query FASTA → $QUERY_FAA" | tee -a "$LOG"

ap "$CLONOTYPES" "$PAIR_CSV" "$QUERY_FAA" <<'PY'
import sys, pandas as pd, pathlib as P
src, out_csv, out_faa = sys.argv[1], sys.argv[2], sys.argv[3]
df = pd.read_csv(src, sep=None, engine="python")
cand_cell = [c for c in df.columns if c.lower() in ("cell_id","cellid","barcode","obs_names","cell","cellbarcode")]
if not cand_cell:
    raise SystemExit("[ERROR] Could not find a cell ID column (cell_id/cellid/barcode/obs_names/cell)")
cell_col = cand_cell[0]
cands = ["heavy_junction_aa","full_query_cdr3","query_cdr3","cdr3aa",
         "cdr3","cdr3_aa","junction_aa","IR_VDJ_1_junction_aa"]
seq_col = next((c for c in cands if c in df.columns), None)
if seq_col is None:
    raise SystemExit("[ERROR] Could not find a CDR‑H3 column in clonotypes CSV/TSV.")
pairs = df[[cell_col, seq_col]].rename(columns={cell_col:"cell_id", seq_col:"cdr3aa"}).dropna()
pairs["cell_id"] = pairs["cell_id"].astype(str)
pairs["cdr3aa"]  = pairs["cdr3aa"].astype(str)
pairs = pairs[(pairs["cdr3aa"].str.len() > 0)].drop_duplicates(subset=["cell_id"])
pairs.to_csv(out_csv, index=False)
with open(out_faa, "w") as fh:
    for r in pairs.itertuples(index=False):
        fh.write(f">{r.cell_id}\n{r.cdr3aa}\n")
print(f"[INFO] Wrote {len(pairs)} unique cell_id/CDR3 pairs", file=sys.stderr)
PY

# --- STEP 2: mmseqs search + convertalis -------------------------------------
TMPDIR=$(mktemp -d); trap 'rm -rf "$TMPDIR"' EXIT
HITS_TSV="$WORK/hits.tsv"

echo "▸ MMseqs2 search (min-id=$MIN_ID, min-cov=$MIN_COV, cov-mode=$COV_MODE, rescore=$RESCORE_MODE) …" | tee -a "$LOG"
"$MMSEQS_BIN" createdb "$QUERY_FAA" "$WORK/query_db" >/dev/null
"$MMSEQS_BIN" search "$WORK/query_db" "$DB_PREFIX" "$WORK/aln" "$TMPDIR" \
  --min-seq-id "$MIN_ID" \
  --cov-mode "$COV_MODE" -c "$MIN_COV" \
  --rescore-mode "$RESCORE_MODE" \
  --threads "$THREADS"
"$MMSEQS_BIN" convertalis "$WORK/query_db" "$DB_PREFIX" "$WORK/aln" "$HITS_TSV" \
  --format-output query,target,pident,qlen,tlen,qstart,qend,tstart,tend,bits,qseq,tseq \
  --threads "$THREADS"

# --- STEP 3: merge hits with full query CDR3 ---------------------------------
MATCHES_TSV="$OUTDIR/cdr3_matches_local95.tsv"
MATCHES_CSV="$OUTDIR/cdr3_matches_local95.csv"

ap "$HITS_TSV" "$PAIR_CSV" "$MATCHES_TSV" "$MATCHES_CSV" <<'PY'
import sys, pandas as pd
hits_tsv, pairs_csv, out_tsv, out_csv = sys.argv[1:5]
hits = pd.read_csv(hits_tsv, sep='\t', header=None,
                   names=["cell_id","target_id","pident","qlen","tlen",
                          "qstart","qend","tstart","tend","global_bitscore","query_cdr3","target_cdr3"])
pairs = pd.read_csv(pairs_csv).rename(columns={"cdr3aa":"full_query_cdr3"})
m = hits.merge(pairs, how="left", on="cell_id")
m = m[["cell_id","target_id","pident","qlen","tlen","qstart","qend","tstart","tend",
       "global_bitscore","query_cdr3","target_cdr3","full_query_cdr3"]]
m.to_csv(out_tsv, sep="\t", index=False)
m.to_csv(out_csv, index=False)
print(f"[INFO] Matches: {len(m)}  |  Unique query cells (with hits): {m['cell_id'].nunique()}", file=sys.stderr)
PY

echo "▸ Tables:"
echo "   $MATCHES_TSV" | tee -a "$LOG"
echo "   $MATCHES_CSV" | tee -a "$LOG"

# --- STEP 4 (optional): UMAP figures + CellTypist labels ---------------------
if [[ -n "$ADATA" ]]; then
  PLOT_DIR="$OUTDIR/svg_figures"; mkdir -p "$PLOT_DIR"
  echo "▸ UMAP figures (AnnData) → $PLOT_DIR" | tee -a "$LOG"
  ap "$ADATA" "$MATCHES_CSV" "$PLOT_DIR" "$SEQ_MAP" "$PRE_PREFIXES" "$POST_PREFIXES" <<'PY'
import sys, os, warnings, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

adata_path, matches_csv, outdir, seq_map, pre_str, post_str = sys.argv[1:7]
warnings.filterwarnings("ignore")
os.makedirs(outdir, exist_ok=True)

# read adata via scanpy if available; otherwise via anndata
try:
    import scanpy as sc
    adata = sc.read_h5ad(adata_path)
    has_scanpy = True
except Exception:
    import anndata as ad
    adata = ad.read_h5ad(adata_path)
    has_scanpy = False

if "X_umap" not in adata.obsm:
    if "X_totalVI_umap" in adata.obsm: adata.obsm["X_umap"] = adata.obsm["X_totalVI_umap"]
    else: raise SystemExit("[ERROR] No UMAP coords found (need X_umap or X_totalVI_umap).")

raw = pd.read_csv(matches_csv, sep=None, engine="python")
seq_col = "query_cdr3" if "query_cdr3" in raw.columns else ("full_query_cdr3" if "full_query_cdr3" in raw.columns else None)
if seq_col is None: raise SystemExit("[ERROR] Could not find a CDR‑H3 column.")
matches = raw.loc[:, ["cell_id", seq_col]].rename(columns={seq_col: "cdr3aa"}).dropna()

# try common label columns
possible_cols = ["celltypist","celltypist_label","CellTypist_prediction","cell_type","celltype"]
label_col = next((c for c in possible_cols if c in adata.obs.columns), None)
if label_col is None:
    raise SystemExit("[ERROR] Could not locate a CellTypist label column in adata.obs")

matches["celltypist_label"] = matches["cell_id"].map(adata.obs[label_col].astype(str))
enriched_csv = os.path.join(outdir, "cdr3_matches_local95_with_celltypist.csv")
matches.to_csv(enriched_csv, index=False)

# optional custom categories
seq_categories = {}
if seq_map and os.path.isfile(seq_map):
    tmp = pd.read_csv(seq_map, sep=None, engine="python", header=None, names=["cdr3aa","category"])
    seq_categories = dict(zip(tmp["cdr3aa"].astype(str), tmp["category"].astype(str)))

# group by sample prefixes
pre_prefixes  = tuple([p.strip() for p in pre_str.split(",") if p.strip()])
post_prefixes = tuple([p.strip() for p in post_str.split(",") if p.strip()])
obs_names = adata.obs_names.astype(str)
all_ids   = matches["cell_id"].astype(str).unique().tolist()
pre_ids   = [cid for cid in all_ids if any(cid.startswith(p) for p in pre_prefixes)]
post_ids  = [cid for cid in all_ids if any(cid.startswith(p) for p in post_prefixes)]

umap = adata.obsm["X_umap"]
idx_pre  = obs_names.isin(pre_ids).values
idx_post = obs_names.isin(post_ids).values

# masks for custom categories
cat2ids = {cat: matches.loc[matches["cdr3aa"]==cdr3, "cell_id"].astype(str).unique().tolist()
           for cdr3,cat in seq_categories.items()}
idx_cat = {cat: obs_names.isin(ids).values for cat, ids in cat2ids.items()}
idx_any_cat = None
for m in idx_cat.values():
    idx_any_cat = m if idx_any_cat is None else (idx_any_cat | m)

def draw_umap(group_name, group_mask, out_stub):
    if group_mask.sum() == 0:
        print(f"[INFO] No cells for {group_name} – skipped.", file=sys.stderr); return
    plt.figure(figsize=(6,6))
    other = ~(idx_pre | idx_post)
    plt.scatter(umap[other,0], umap[other,1], s=8, c="lightgrey", linewidths=0, rasterized=True)
    plt.scatter(umap[group_mask & ~(idx_any_cat if idx_any_cat is not None else np.zeros_like(group_mask,dtype=bool)),0],
                umap[group_mask & ~(idx_any_cat if idx_any_cat is not None else np.zeros_like(group_mask,dtype=bool)),1],
                s=22, edgecolors="black", linewidths=0.25, label=f"{group_name}")
    for cat, mask in idx_cat.items():
        submask = mask & group_mask
        if not submask.any(): continue
        plt.scatter(umap[submask,0], umap[submask,1], s=32, edgecolors="black", linewidths=0.35, label=f"{cat}")
    plt.axis("off"); plt.title(f"UMAP — {group_name}")
    plt.legend(frameon=False, loc="upper right", fontsize=6)
    plt.savefig(f"{out_stub}.svg", bbox_inches="tight")
    plt.savefig(f"{out_stub}.png", bbox_inches="tight", dpi=300)
    print(f"[INFO] Figure → {out_stub}.svg / .png", file=sys.stderr)

draw_umap("NTX 1/2/3 (pre)",  idx_pre,  os.path.join(outdir, "UMAP_cdr3_matches_local95_NTX123"))
draw_umap("NTX 4/5 (post)",  idx_post, os.path.join(outdir, "UMAP_cdr3_matches_local95_NTX45"))
print(f"[INFO] Unique matched cells: {len(all_ids)}  | pre={idx_pre.sum()}  post={idx_post.sum()}", file=sys.stderr)
PY
fi

# --- manifest + checksums -----------------------------------------------------
{
  echo "# Manifest — $(date)"
  echo "DB_DIR         : $DB_DIR"
  echo "CLONOTYPES     : $CLONOTYPES"
  echo "ADATA          : ${ADATA:-<none>}"
  echo "Parameters     : min_id=$MIN_ID  min_cov=$MIN_COV  cov_mode=$COV_MODE  rescore=$RESCORE_MODE  threads=$THREADS"
  echo
  echo "Outputs:"
  echo "  $MATCHES_TSV"
  echo "  $MATCHES_CSV"
  [[ -n "$ADATA" ]] && {
    echo "  $OUTDIR/svg_figures/UMAP_cdr3_matches_local95_NTX123.svg(.png)"
    echo "  $OUTDIR/svg_figures/UMAP_cdr3_matches_local95_NTX45.svg(.png)"
    echo "  $OUTDIR/svg_figures/cdr3_matches_local95_with_celltypist.csv"
  }
} > "$OUTDIR/MANIFEST.txt"

if command -v sha256sum >/dev/null 2>&1; then
  {
    echo "# SHA256"
    for f in "$MATCHES_TSV" "$MATCHES_CSV" "$OUTDIR"/svg_figures/*.{svg,png,csv}; do
      [[ -e "$f" ]] && sha256sum "$f"
    done
  } > "$OUTDIR/CHECKSUMS.sha256"
fi

echo "✅ Done. See:"
echo "   • $OUTDIR/MANIFEST.txt"
echo "   • $OUTDIR/LOG.txt"
