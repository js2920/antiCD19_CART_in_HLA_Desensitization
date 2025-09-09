#!/usr/bin/env bash
# scripts/build_oas_ighg_cdr3_db.sh
# Build a CDR-H3-only FASTA from OAS IGHG CSVs and (re)create an MMseqs2 DB.
#
# see: https://opig.stats.ox.ac.uk/webapps/oas/downloads/
# and their Publication Observed Antibody Space: A diverse database of cleaned, annotated, and translated unpaired and paired antibody sequences
# Tobias H Olsen, Fergus Boyles, Charlotte M Deane  
# Protein Sci
# 2022 Jan;31(1):141-146.
# doi: 10.1002/pro.4205. Epub 2021 Oct 29.
# for downloading the Observed Antobody Space Dataset on which this script builds
#
# Usage:
#   scripts/build_oas_ighg_cdr3_db.sh \
#     --input /path/to/oas_csv_ighg \
#     --outdir /path/to/outdir \
#     [--column sequence_alignment_aa] [--min 12] [--max 32] [--force]
#
# Example:
#   scripts/build_oas_ighg_cdr3_db.sh -i ~/data/oas_csv_ighg -o ~/data/oas_ighg_cdr3
#
# Requirements provided by the env file:
#   mmseqs2, csvkit (csvcut), gawk, coreutils (realpath, sha256sum), grep, sed

set -Eeuo pipefail

# ---- defaults (can be overridden via CLI) -----------------------------------
OAS_CSV="${OAS_CSV:-$HOME/data/oas_csv_ighg}"
DB_DIR="${DB_DIR:-$HOME/data/oas_ighg_cdr3}"
COL="${COL:-sequence_alignment_aa}"
MIN="${MIN:-12}"
MAX="${MAX:-32}"
FORCE="${FORCE:-0}"
LOG=""

# ---- usage ------------------------------------------------------------------
usage() {
  cat <<EOF
Build CDR-H3 FASTA and MMseqs2 DB from OAS IGHG CSVs.

Options:
  -i, --input DIR        Directory with *.csv (OAS IGHG) [default: $OAS_CSV]
  -o, --outdir DIR       Output dir for FASTA & DB      [default: $DB_DIR]
  -c, --column NAME      CSV column with AA seqs        [default: $COL]
      --min N            Min CDR-H3 length (aa)         [default: $MIN]
      --max N            Max CDR-H3 length (aa)         [default: $MAX]
  -f, --force            Overwrite existing outputs
  -l, --log FILE         Log file (default: outdir/build_*.log)
  -h, --help             Show this help

Notes:
  • Heuristic: match IMGT-like CDR-H3 motif C[A-Z]{min,max}W within the AA column.
  • 'mmseqs createdb' is single-threaded; no --threads option is used here.
EOF
}

# ---- parse args --------------------------------------------------------------
if [[ $# -eq 0 ]]; then usage; exit 0; fi
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input)   OAS_CSV="$2"; shift 2;;
    -o|--outdir)  DB_DIR="$2"; shift 2;;
    -c|--column)  COL="$2"; shift 2;;
    --min)        MIN="$2"; shift 2;;
    --max)        MAX="$2"; shift 2;;
    -f|--force)   FORCE=1; shift 1;;
    -l|--log)     LOG="$2"; shift 2;;
    -h|--help)    usage; exit 0;;
    *) echo "Unknown option: $1"; usage; exit 2;;
  esac
done

# ---- deps check --------------------------------------------------------------
require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "Missing dependency: $1"; exit 127; }; }
for cmd in mmseqs csvcut gawk realpath sha256sum grep sed; do require_cmd "$cmd"; done

# ---- inputs & outputs --------------------------------------------------------
if [[ ! -d "$OAS_CSV" ]]; then
  echo "Input directory not found: $OAS_CSV"; exit 1
fi
mkdir -p "$DB_DIR"

FASTA="$DB_DIR/cdr3_ighg.faa"
DB_PREFIX="$DB_DIR/cdr3_ighg_db"
CHECKSUMS="$DB_DIR/checksums.sha256"
TIMESTAMP="$(date +"%Y%m%d_%H%M%S")"
LOG="${LOG:-$DB_DIR/build_${TIMESTAMP}.log}"

# ---- logging ----------------------------------------------------------------
# tee everything to log
exec > >(tee -a "$LOG") 2>&1
echo "=== $(basename "$0") started: $(date -Is) ==="
echo "Input dir: $OAS_CSV"
echo "Out dir  : $DB_DIR"
echo "Column   : $COL"
echo "Len range: ${MIN}-${MAX}"
echo "Force    : $FORCE"
echo "Log file : $LOG"
echo

# ---- discover CSVs -----------------------------------------------------------
mapfile -d '' CSV_FILES < <(find -L "$OAS_CSV" -maxdepth 1 -type f -name '*.csv' -print0)
if [[ ${#CSV_FILES[@]} -eq 0 ]]; then
  echo "No CSV files found in: $OAS_CSV"; exit 1
fi
echo "Found ${#CSV_FILES[@]} CSV files."

# ---- guard existing outputs --------------------------------------------------
if [[ -e "$FASTA" || -e "${DB_PREFIX}"* ]]; then
  if [[ "$FORCE" -ne 1 ]]; then
    echo "Outputs exist. Use --force to overwrite:"
    ls -lh "$FASTA" 2>/dev/null || true
    ls -lh "${DB_PREFIX}"* 2>/dev/null || true
    exit 3
  fi
  rm -f "$FASTA" "${DB_PREFIX}"* "$CHECKSUMS"
fi

# ---- build FASTA -------------------------------------------------------------
TMP_FASTA="${FASTA}.tmp"
> "$TMP_FASTA"
echo "▸ Building new CDR-H3 FASTA …"

# For reproducible locale behavior
export LC_ALL=C

for csv in "${CSV_FILES[@]}"; do
  # Extract the target AA column; drop header; regex match CDR-H3; write FASTA.
  csvcut -c "$COL" "$csv" | tail -n +2 | \
  gawk -v file="$(realpath "$csv")" -v min="$MIN" -v max="$MAX" '
    BEGIN { OFS="" }
    {
      seq=$1
      if (seq=="" || seq=="NULL") next
      regex=sprintf("C[A-Z]{%d,%d}W", min, max)
      if (match(seq, regex)) {
        cdr3=substr(seq, RSTART, RLENGTH)
        printf ">%s_%d\n%s\n", file, NR+1, cdr3
      }
    }' >> "$TMP_FASTA"
done

mv -f "$TMP_FASTA" "$FASTA"
CNT=$(grep -c '^>' "$FASTA" || echo 0)
echo "   $CNT CDR-H3 sequences written to $(basename "$FASTA")"

# ---- create / overwrite MMseqs2 DB ------------------------------------------
echo "▸ Re-building MMseqs2 DB …"
mmseqs createdb "$FASTA" "$DB_PREFIX"
echo "   DB files:"
ls -lh "${DB_PREFIX}"* 2>/dev/null || true

# ---- checksums for reproducibility ------------------------------------------
echo "▸ Writing SHA256 checksums …"
sha256sum "$FASTA" > "$CHECKSUMS"
# include all DB artifacts
find "$DB_DIR" -maxdepth 1 -type f -name "$(basename "$DB_PREFIX")*" -print0 \
  | xargs -0 -r sha256sum >> "$CHECKSUMS"
echo "   $(wc -l < "$CHECKSUMS" | sed 's/^/   /') entries in $(basename "$CHECKSUMS")"

# ---- summary ----------------------------------------------------------------
echo "✅ Done – database ready: $DB_PREFIX"
du -sh "$DB_DIR"/cdr3_ighg_* 2>/dev/null | awk '{print "   " $2 "  " $1}'
echo "Versions:"
echo "   $(mmseqs --version 2>/dev/null || true)"
echo "   csvcut $(csvcut --version 2>&1 | head -n1)"
echo "   gawk $(gawk --version | head -n1)"
echo "   coreutils $(sha256sum --version | head -n1)"
echo "Log saved to: $LOG"
echo "=== Finished: $(date -Is) ==="
