# see: https://opig.stats.ox.ac.uk/webapps/oas/downloads/
and their Publication Observed Antibody Space: A diverse database of cleaned, annotated, and translated unpaired and paired antibody sequences
Tobias H Olsen, Fergus Boyles, Charlotte M Deane   Protein Sci
2022 Jan;31(1):141-146.
doi: 10.1002/pro.4205. Epub 2021 Oct 29.
for downloading the Observed Antobody Space Dataset on which this script builds

**What this sub repository provides**
- Build a CDR‑H3 FASTA + MMseqs2 DB from OAS IGHG CSVs.
- Align clonotype CDR‑H3s to that DB.
- Optional: UMAP overlays using CellTypist labels from your `.h5ad`.

## Quick start

### 1) Build the database
```bash
mamba env create -f env/oas-cdr3-db.yml
conda activate oas-cdr3-db

scripts/build_oas_ighg_cdr3_db.sh \
  --input  /home/data/oas_csv_ighg \
  --outdir /home/oas_ighg_cdr3 \
  --force

mamba env create -f env/cdr3-align.yml
conda activate cdr3-align

scripts/align_cdr3_against_oas.sh \
  --db-dir     /home/jensf/data/oas_ighg_cdr3 \
  --clonotypes dandelion_inputs/analysis_ready/my_merged_contigs.tsv \
  --adata      totalvi_celltypist_inputs/processed_adata/processed_adata_NTX_all_samples.h5ad \
  --outdir     review_align_outputs \
  --threads    16
