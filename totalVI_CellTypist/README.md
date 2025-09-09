# TOTALVI + CellTypist (Reproducible from Zenodo Archive)

This repository reproduces the analysis using only the Zenodo-deposited zip archive
from https://zenodo.org/records/17063869

## Inputs

Download the Zenodo zip archive and keep the structure intact:

totalvi_celltypist_inputs/
processed_adata/processed_adata_NTX_all_samples.h5ad # optional
cellbender_h5/NTX_1/cellbender_output_NTX_1_filtered.h5
cellbender_h5/NTX_2/cellbender_output_NTX_2_filtered.h5
cellbender_h5/NTX_3/cellbender_output_NTX_3_filtered.h5
cellbender_h5/NTX_4/cellbender_NTX_4_filtered.h5
cellbender_h5/NTX_5/cellbender_NTX_5_filtered.h5
adt_reference/TotalSeq_C_Human_Universal_Cocktail_399905_Antibody_reference_UMI_counting.csv
models/Immune_All_Low.pkl
models/Cells_Human_Tonsil.pkl


## Environment

```bash
# Option A: pip (CPU example)
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
# If you have CUDA, install a matching torch build from https://pytorch.org/get-started/locally/

## RUN

# Using a zip file
python scripts/run_totalvi_celltypist_repro.py \
  --zenodo-zip /path/to/totalvi_celltypist_inputs.zip \
  --outdir outputs \
  --max-epochs 150 \
  --seed 0

# Or using an already-extracted directory
python scripts/run_totalvi_celltypist_repro.py \
  --zenodo-dir /path/to/totalvi_celltypist_inputs \
  --outdir outputs

Outputs

outputs/processed_adata_NTX_all_samples.h5ad
outputs/PREvsPOST/plasma_DE_NTX13_vs_NTX45_bonferroni.csv
outputs/figures/ (UMAPs, boxplots, dotplots, MA plot, cell-cycle proportions)

Notes
The script uses only files in the Zenodo archive 
Random seeds are set for reproducibility.
GPU will be used if available; otherwise CPU is used.
