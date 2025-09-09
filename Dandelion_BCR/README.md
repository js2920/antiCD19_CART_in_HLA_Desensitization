# Reproducible single‑cell BCR inputs (Dandelion) — pipeline wrapper

the authors of Dandelion and Immcantation have to be acknowledged and thanked at the beginning

Suo C, Polanski K, Dann E, Lindeboom RGH, Vilarrasa-Blasi R, Vento-Tormo R, Haniffa M, Meyer KB, Dratva LM, Tuong ZK, Clatworthy MR, Teichmann SA. Dandelion uses the single-cell adaptive immune receptor repertoire to explore lymphocyte developmental origins. Nat Biotechnol. 2024;42(1):40-51.

&

Gabernet G, Marquez S, Bjornson R, Peltzer A, Meng H, Aron E, Lee NY, Jensen CG, Ladd D, Polster M, Hanssen F, Heumos S, nf-core c, Yaari G, Kowarik MC, Nahnsen S, Kleinstein SH. nf-core/airrflow: An adaptive immune receptor repertoire analysis workflow employing the Immcantation framework. PLoS Comput Biol. 2024;20(7):e1012265.

This repository contains a lightweight, reproducible wrapper script and documentation to (optionally) **preprocess**, **merge**, and **integrate** single‑cell BCR V(D)J data using the **Dandelion** container, and (optionally) run **Change‑O** clonotype calling, matching the archive layout obtained via the zenodo repository:

```text
dandelion_inputs_TCR/
├── analysis_ready/
│   ├── meta.csv
│   ├── per_cell_TCR_table.csv
│   ├── sc-dandelion_latest.sif
│   ├── adata_with_TCR_integration.h5ad
│   ├── processed_adata_NTX_all_samples.h5ad
│   ├── NTX_1_TCR_fastq/
│   │   └── dandelion/all_contig_dandelion.tsv
│   ├── NTX_2_TCR_fastq/
│   │   └── dandelion/all_contig_dandelion.tsv
│   ├── NTX_3_TCR_fastq/
│   │   └── dandelion/all_contig_dandelion.tsv
│   ├── NTX_4_TCR_fastq/
│   │   └── dandelion/all_contig_dandelion.tsv
│   └── NTX_5_TCR_fastq/
│       └── dandelion/all_contig_dandelion.tsv
├── base_inputs/
│   ├── adata/processed_adata_NTX_all_samples.h5ad
│   ├── contigs/
│   │   ├── NTX_1.all_contig_dandelion.tsv
│   │   ├── NTX_2.all_contig_dandelion.tsv
│   │   ├── NTX_3.all_contig_dandelion.tsv
│   │   ├── NTX_4.all_contig_dandelion.tsv
│   │   └── NTX_5.all_contig_dandelion.tsv
│   ├── containers/sc-dandelion_latest.sif
│   ├── germlines/corrected_germline.fasta
│   ├── igblast_internal_data/
│   └── meta/meta.csv
├── README_dandelion_inputs_TCR.txt
└── manifest_dandelion_TCR.tsv
```


> **What this repo provides**
>
> * `scripts/run_pipeline.sh` — one‑file wrapper around the **Dandelion** container:
>   * optional `dandelion-preprocess` (Cell Ranger outputs → corrected AIRR tables)
>   * merge per‑sample Dandelion AIRR TSVs
>   * integrate VDJ with your `.h5ad` GEX (saves `vdj.h5ddl` + `adata_with_bcr_integration.h5ad`)
>   * optional **Change‑O** clonotype calling (`changeo-clonotypes`)
> * A clear overview of container dependencies (Dandelion SIF; Immcantation Docker as an alternative).

---

## Dependencies

### 1) Dandelion Singularity container

Dandelion ships a ready‑to‑run container that bundles IgBLAST, Change‑O, SHazaM, TIgGER, etc. Pull it once and point the script at the resulting `sc-dandelion_latest.sif`:

```bash
# Pull (creates sc-dandelion_latest.sif in the CWD)
singularity pull library://kt16/default/sc-dandelion:latest

# from repo root:
chmod +x scripts/run_pipeline.sh

If you prefer to run IgBLAST/Change‑O outside the Dandelion SIF, use Immcantation’s Docker images:
# pull the suite (contains Change-O, SHazaM, TIgGER, IgBLAST, etc.)
docker pull immcantation/suite:4.6.0   # or a specific stable tag

