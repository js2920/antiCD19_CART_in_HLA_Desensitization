#!/usr/bin/env bash
#
# car_extraction_from_bundle.sh
#
# Run the KZV‑101 extraction + alignment + barcode intersection pipeline
# using ONLY the inputs included in CAR_extraction_inputs.zip.
#
# It will:
#   1) Detect samples from base_inputs/fastq/*_R[12]_*fastq.gz
#      (or read an optional sample_map.tsv if present)
#   2) Extract 10x barcodes & UMIs (UMI-tools)
#   3) Align R2 to the KZV‑101 mini-reference (STAR; uses bundled index)
#   4) Count KYV‑101 UMIs per barcode (≥1 read, ≥1 UMI)
#   5) Keep only barcodes present in the bundled CellBender list
#
# Notes
#  - Expects bundle layout as shown by `zipinfo` in the user’s message:
#       CAR_extraction_inputs/
#         base_inputs/
#           fastq/
#           reference/TCR_KYV101_specific.fasta
#           star_index/STAR_index_KZV101/ (prebuilt)
#           cellbender_lists/*.csv
#         sample_map.tsv              (optional; see README for format)
#  - Output goes to:  <OUTPUT_DIR>/<sample>/
#
# Requirements in PATH: umi_tools, STAR, samtools, awk, zcat, unzip (if unzipping)
# Recommended install method: conda env from environment.yml
#
# -----------------------------------------------------------------------------

set -euo pipefail
IFS=$'\n\t'

show_help() {
cat <<'EOF'
Usage:
  car_extraction_from_bundle.sh --bundle-dir <DIR> [--output <DIR>] [--threads N] [--bam-sort-ram BYTES]

Arguments:
  --bundle-dir   Path to extracted CAR_extraction_inputs/ directory (required).
  --output       Output directory (default: ./analysis_output)
  --threads      Threads for STAR/samtools (default: detected via nproc or 16)
  --bam-sort-ram Bytes for STAR sort RAM (default: 128000000000 = 128 GB)

Behavior:
  - Auto-discovers samples from base_inputs/fastq/*_R1_*.fastq.gz and matching R2.
  - Derives sample name by trimming everything from first '_GEX' onward.
    Example: ES_NTX_4_GEX_S4_R1_001.fastq.gz -> sample 'ES_NTX_4'.
  - CellBender file expected at base_inputs/cellbender_lists/<sample>_cellbender_barcodes.csv
    (If not exact match, will try to find a unique CSV that starts with <sample>.)
  - If sample_map.tsv exists at bundle root, it will be read first. See README.

Outputs (per sample):
  KZV101_UmiCounts.tsv
  CellBender_cells.norm
  KZV101_positive_cells.tsv
  KZV101_mapped.bam/.bai
  umi_extract.log
  STAR logs (Aligned.sortedByCoord.out.bam + STAR outputs)
EOF
}

# ---- defaults
OUTPUT_DIR="${PWD}/analysis_output"
THREADS="${THREADS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 16)}"
BAM_SORT_RAM="${BAM_SORT_RAM:-128000000000}"

BUNDLE_DIR=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bundle-dir) BUNDLE_DIR="$2"; shift 2;;
    --output) OUTPUT_DIR="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --bam-sort-ram) BAM_SORT_RAM="$2"; shift 2;;
    -h|--help) show_help; exit 0;;
    *) echo "Unknown arg: $1"; show_help; exit 1;;
  esac
done

if [[ -z "${BUNDLE_DIR}" ]]; then
  echo "ERROR: --bundle-dir is required"; show_help; exit 1
fi
if [[ ! -d "${BUNDLE_DIR}" ]]; then
  echo "ERROR: bundle dir not found: ${BUNDLE_DIR}"; exit 1
fi

# ---- tool sanity check
for t in umi_tools STAR samtools awk zcat; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not in PATH"; exit 1; }
done

# ---- derive key paths from bundle
FASTQ_DIR="${BUNDLE_DIR}/base_inputs/fastq"
CB_DIR="${BUNDLE_DIR}/base_inputs/cellbender_lists"
REF_FASTA="${BUNDLE_DIR}/base_inputs/reference/TCR_KYV101_specific.fasta"
STAR_INDEX="${BUNDLE_DIR}/base_inputs/star_index/STAR_index_KZV101"

mkdir -p "${OUTPUT_DIR}"
exec > >(tee -a "${OUTPUT_DIR}/KZV101_CellBender.log") 2>&1
echo "[$(date)] KZV‑101 pipeline (bundle mode) started."
echo "Bundle: ${BUNDLE_DIR}"
echo "Output: ${OUTPUT_DIR}"
echo "Threads: ${THREADS}, STAR sort RAM: ${BAM_SORT_RAM}"

# ---- STAR index (use bundled if present; else build from bundled FASTA)
if [[ -f "${STAR_INDEX}/SA" ]]; then
  echo "[$(date)] Using bundled STAR index at: ${STAR_INDEX}"
else
  echo "[$(date)] Bundled STAR index not found; building from bundled FASTA…"
  mkdir -p "${STAR_INDEX}"
  STAR  --runThreadN           "${THREADS}"          \
        --runMode              genomeGenerate        \
        --genomeDir            "${STAR_INDEX}"       \
        --genomeFastaFiles     "${REF_FASTA}"        \
        --genomeSAindexNbases  3                     \
        --outFileNamePrefix    "${STAR_INDEX}/"
fi

# ---- helpers

# Turn a possibly relative path into absolute, rooted at FASTQ_DIR if no leading slash.
abs_fastq() {
  local p="$1"
  if [[ "$p" == /* ]]; then
    printf "%s" "$p"
  else
    printf "%s/%s" "$FASTQ_DIR" "$p"
  fi
}

# Find CellBender CSV for sample key (e.g., ES_NTX_4)
find_cellbender_csv() {
  local sample_key="$1"
  local exact="${CB_DIR}/${sample_key}_cellbender_barcodes.csv"
  if [[ -f "$exact" ]]; then
    printf "%s" "$exact"; return 0
  fi
  # fallback: any unique CSV starting with sample_key
  local cand
  mapfile -t cand < <(find "${CB_DIR}" -maxdepth 1 -type f -name "${sample_key}*.csv" | sort)
  if [[ "${#cand[@]}" -eq 1 ]]; then
    printf "%s" "${cand[0]}"; return 0
  fi
  return 1
}

# UMI counting (same logic as your original script)
count_umis() {
  local bam="$1" out_tsv="$2"
  samtools view "$bam" | awk -v OFS='\t' '
    {
      split($1,a,"_");                          # readName …_<BARCODE>_<UMI>
      if (length(a) >= 3) {
        cell=a[length(a)-1]; umi=a[length(a)];
        if (index(cell,"N")==0) {               # skip ambiguous barcodes
          key=cell "_" umi
          if (!(key in seen)) {seen[key]=1; count[cell]++}
        }
      }
    }
    END { for (c in count) print c, count[c] }
  ' > "$out_tsv"
}

# Normalize CellBender list (drop 10x suffix and N-containing barcodes)
normalize_cellbender() {
  local in_csv="$1" out_norm="$2"
  tail -n +2 "$in_csv" | sed -E 's/-[0-9]+$//' | grep -v 'N' > "$out_norm"
}

# ---- sample discovery

declare -a SAMPLES
declare -A R1S R2S CB

SAMPLE_MAP="${BUNDLE_DIR}/sample_map.tsv"
if [[ -f "${SAMPLE_MAP}" ]]; then
  echo "[$(date)] Reading sample_map.tsv from bundle: ${SAMPLE_MAP}"
  # Expected columns (tab-separated). Header optional:
  # sample_id   r1_fastq   r2_fastq   cellbender_csv(optional)
  # r1/r2 can be basenames (relative to base_inputs/fastq) or absolute paths.
  while IFS=$'\t' read -r c1 c2 c3 c4; do
    # skip blank and comment lines
    [[ -z "${c1:-}" || "${c1:0:1}" == "#" ]] && continue
    # skip header if present
    if [[ "${c1,,}" == "sample" || "${c1,,}" == "sample_id" ]]; then
      continue
    fi
    local_sample="${c1}"
    r1_path="$(abs_fastq "${c2}")"
    r2_path="$(abs_fastq "${c3}")"
    if [[ ! -f "$r1_path" || ! -f "$r2_path" ]]; then
      echo "ERROR in sample_map.tsv: FASTQ not found for sample ${local_sample}"
      exit 1
    fi
    cb_csv=""
    if [[ -n "${c4:-}" ]]; then
      cb_csv="${c4}"
      [[ "$cb_csv" != /* ]] && cb_csv="${CB_DIR}/$(basename "$cb_csv")"
      [[ ! -f "$cb_csv" ]] && { echo "ERROR: CellBender CSV not found: $cb_csv"; exit 1; }
    else
      # derive sample key used for CellBender naming (truncate at first _GEX)
      sample_key="${local_sample%%_GEX*}"
      if ! cb_csv="$(find_cellbender_csv "$sample_key")"; then
        echo "ERROR: Could not find CellBender list for sample ${local_sample} (key: ${sample_key})"
        exit 1
      fi
    fi
    SAMPLES+=("$local_sample")
    R1S["$local_sample"]="$r1_path"
    R2S["$local_sample"]="$r2_path"
    CB["$local_sample"]="$cb_csv"
  done < "${SAMPLE_MAP}"
else
  echo "[$(date)] sample_map.tsv not found; auto-discovering samples from FASTQ names…"
  mapfile -t r1s < <(find "${FASTQ_DIR}" -maxdepth 1 -type f -name "*_R1_*.fastq.gz" | sort)
  for r1 in "${r1s[@]}"; do
    base="$(basename "$r1")"
    # derive matching R2
    r2="${r1/_R1_/_R2_}"
    [[ -f "$r2" ]] || { echo "WARN: No matching R2 for $r1; skipping"; continue; }
    # sample name = trim at first '_GEX'
    sample_key="${base%%_GEX*}"          # e.g. ES_NTX_4
    sample="${sample_key}"
    # prefer user-friendly names without _GEX*: if trimming changed nothing, also drop trailing _R1...
    # (but for your files the line above already yields ES_NTX_4 / ES_NTX_5)
    if ! cb_csv="$(find_cellbender_csv "$sample_key")"; then
      echo "ERROR: Could not find CellBender list for ${sample_key} in ${CB_DIR}"
      exit 1
    fi
    SAMPLES+=("$sample")
    R1S["$sample"]="$r1"
    R2S["$sample"]="$r2"
    CB["$sample"]="$cb_csv"
  done
fi

if [[ "${#SAMPLES[@]}" -eq 0 ]]; then
  echo "ERROR: No samples discovered."; exit 1
fi

echo "[$(date)] Discovered samples:"
for s in "${SAMPLES[@]}"; do
  echo "  - ${s}"
  echo "      R1: ${R1S[$s]}"
  echo "      R2: ${R2S[$s]}"
  echo "      CB: ${CB[$s]}"
done

# ---- per-sample processing
process_sample () {
  local SAMPLE="$1"
  local R1="${R1S[$SAMPLE]}"
  local R2="${R2S[$SAMPLE]}"
  local CB_LIST="${CB[$SAMPLE]}"

  echo "[$(date)] ▶  ${SAMPLE}"
  local SAMP_DIR="${OUTPUT_DIR}/${SAMPLE}"; mkdir -p "$SAMP_DIR"

  # A) Barcode + UMI extraction
  umi_tools extract \
    --bc-pattern='CCCCCCCCCCCCCCCCNNNNNNNNNN' \
    --stdin         "$R1"             --stdout "$SAMP_DIR/R1_ex.fq.gz" \
    --read2-in      "$R2"             --read2-out "$SAMP_DIR/R2_ex.fq.gz" \
    --log           "$SAMP_DIR/umi_extract.log"

  # B) STAR alignment (R2 only)
  STAR  --runThreadN             "${THREADS}"                  \
        --genomeDir              "${STAR_INDEX}"               \
        --readFilesIn            "$SAMP_DIR/R2_ex.fq.gz"       \
        --readFilesCommand       zcat                          \
        --outFileNamePrefix      "$SAMP_DIR/"                  \
        --outSAMtype             BAM SortedByCoordinate        \
        --outFilterMatchNmin     60                            \
        --outFilterMismatchNmax  4                             \
        --chimSegmentMin         15                            \
        --chimScoreMin           1                             \
        --outFilterMultimapNmax  1                             \
        --limitBAMsortRAM        "${BAM_SORT_RAM}"

  # Keep mapped reads only
  samtools view -@ "${THREADS}" -b -F 4 \
      "$SAMP_DIR/Aligned.sortedByCoord.out.bam" \
      > "$SAMP_DIR/KZV101_mapped.bam"
  samtools index "$SAMP_DIR/KZV101_mapped.bam"

  # C) UMI counts per cell barcode
  count_umis "$SAMP_DIR/KZV101_mapped.bam" "$SAMP_DIR/KZV101_UmiCounts.tsv"

  # D) Normalize CellBender list
  normalize_cellbender "$CB_LIST" "$SAMP_DIR/CellBender_cells.norm"

  # E) Intersect
  awk 'NR==FNR { ok[$1]=1; next } ($1 in ok) { print }' \
      "$SAMP_DIR/CellBender_cells.norm" \
      "$SAMP_DIR/KZV101_UmiCounts.tsv" \
      > "$SAMP_DIR/KZV101_positive_cells.tsv"

  echo "[$(date)] ▲  ${SAMPLE} done – $(wc -l < "$SAMP_DIR/KZV101_positive_cells.tsv") barcodes kept."

  # Cleanup intermediates
  rm -f "$SAMP_DIR/R1_ex.fq.gz" "$SAMP_DIR/R2_ex.fq.gz"
}

for s in "${SAMPLES[@]}"; do
  process_sample "$s"
done

echo "[$(date)] Pipeline complete. Results in ${OUTPUT_DIR}/<sample>/KZV101_positive_cells.tsv"
