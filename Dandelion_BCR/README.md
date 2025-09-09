# Reproducible single‑cell BCR inputs (Dandelion) — pipeline wrapper

the authors of Dandelion and Immcantation have to be acknowledged and thanked at the beginning

Suo C, Polanski K, Dann E, Lindeboom RGH, Vilarrasa-Blasi R, Vento-Tormo R, Haniffa M, Meyer KB, Dratva LM, Tuong ZK, Clatworthy MR, Teichmann SA. Dandelion uses the single-cell adaptive immune receptor repertoire to explore lymphocyte developmental origins. Nat Biotechnol. 2024;42(1):40-51.

&

Gabernet G, Marquez S, Bjornson R, Peltzer A, Meng H, Aron E, Lee NY, Jensen CG, Ladd D, Polster M, Hanssen F, Heumos S, nf-core c, Yaari G, Kowarik MC, Nahnsen S, Kleinstein SH. nf-core/airrflow: An adaptive immune receptor repertoire analysis workflow employing the Immcantation framework. PLoS Comput Biol. 2024;20(7):e1012265.

This repository contains a lightweight, reproducible wrapper script and documentation to (optionally) **preprocess**, **merge**, and **integrate** single‑cell BCR V(D)J data using the **Dandelion** container, and (optionally) run **Change‑O** clonotype calling, matching the archive layout obtained via the zenodo repository:

dandelion_inputs/
├── manifest_dandelion.tsv
├── README_dandelion_inputs.txt
├── analysis_ready/
│   ├── meta.csv
│   ├── sc-dandelion_latest.sif
│   ├── processed_adata_NTX_all_samples.h5ad
│   ├── adata_with_bcr_integration.h5ad   — often produced; can be re-created
│   ├── my_merged_contigs.tsv             — often produced; can be re-created
│   └── NTX_fastq/
│       └── dandelion/
│           └── all_contig_dandelion.tsv  — per‑sample AIRR TSV from Dandelion
├── base_inputs/
│   ├── containers/
│   │   └── sc-dandelion_latest.sif       — SIF (duplicate of above is OK)
│   ├── meta/
│   │   └── meta.csv                      — Dandelion CLI metadata
│   ├── adata/
│   │   └── processed_adata_NTX_all_samples.h5ad  — GEX
│   ├── contigs/
│   │   └── NTX.all_contig_dandelion.tsv  — optional centralized contig TSVs
│   ├── igblast_internal_data/            — optional custom IgBLAST data
│   └── germlines/
│       └── corrected_germline.fasta      — optional custom germlines
└── env/
    ├── conda_env_export.yml
    ├── pip_freeze.txt
    ├── python_versions.txt
    ├── docker_version.txt
    ├── igblast_version.txt
    └── singularity_version.txt



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

### 1) Dandelion Singularity/Apptainer container (required)

Dandelion ships a ready‑to‑run container that bundles IgBLAST, Change‑O, SHazaM, TIgGER, etc. Pull it once and point the script at the resulting `sc-dandelion_latest.sif`:

```bash
# Pull (creates sc-dandelion_latest.sif in the CWD)
singularity pull library://kt16/default/sc-dandelion:latest
