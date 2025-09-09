# KYV-101 CAR Extraction (Single-cell) — Bundle Runner

This repository runs a reproducible KYV‑101 extraction pipeline directly from a data bundle `CAR_extraction_inputs.zip`.

**What it does:**

1. Extracts 10x barcodes & UMIs from R1 / R2 (UMI-tools).
2. Aligns R2 reads to a KYV‑101 mini-reference (STAR) using the bundled index.
3. Counts UMIs per barcode (≥1 read, unique UMI per cell).
4. Intersects barcodes with CellBender's filtered list (bundled).

## Quick start

```bash
# 1) Create environment
conda env create -f environment_car_extraction.yml
conda activate car-extraction

# 2) Run with a local ZIP
scripts/run_from_zip_car_extraction.sh --zip ./archives/CAR_extraction_inputs.zip

#    or download from Zenodo 
scripts/run_from_zip_car_extraction.sh --zip "https://zenodo.org/records/17063869"
