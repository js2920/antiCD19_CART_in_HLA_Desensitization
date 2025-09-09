#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Reproducible single-cell multi-modal pipeline using only files from a Zenodo archive.

What it does:
  • TOTALVI for RNA+ADT
  • Scrublet for doublets
  • CellTypist (using the bundled model .pkl, no internet)
  • Cell cycle scoring
  • UMAPs (clusters, sample_id, cell cycle)
  • Cell-cycle proportion bar plot by merged_label
  • Saves processed AnnData

Inputs (required):
  - Either --zenodo-zip /path/to/totalvi_celltypist_inputs.zip
    OR     --zenodo-dir /path/to/extracted/totalvi_celltypist_inputs

Assumed Zenodo layout:
  totalvi_celltypist_inputs/
    processed_adata/
      processed_adata_NTX_all_samples.h5ad        (optional; not required)
    cellbender_h5/
      NTX_1/cellbender_output_NTX_1_filtered.h5
      NTX_2/cellbender_output_NTX_2_filtered.h5
      NTX_3/cellbender_output_NTX_3_filtered.h5
      NTX_4/cellbender_NTX_4_filtered.h5
      NTX_5/cellbender_NTX_5_filtered.h5
    adt_reference/
      TotalSeq_C_Human_Universal_Cocktail_399905_Antibody_reference_UMI_counting.csv
    models/
      Cells_Human_Tonsil.pkl
      Immune_All_Low.pkl

Outputs:
  outputs/
    processed_adata_NTX_all_samples.h5ad
    figures/*.png  (UMAPs, celltypist confidence, cell-cycle proportions)

Repro notes:
  - No internet access required.
  - Seeds set for numpy and scvi.
  - Uses GPU if available, else CPU.
"""

import os
import re
import sys
import argparse
import warnings
import zipfile
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
import anndata as ad
import torch
from scipy.sparse import issparse

import scvi
from scvi.model import TOTALVI

import scrublet as scr
import celltypist

warnings.filterwarnings("ignore")

###############################################################################
# Utilities
###############################################################################

def _log_versions():
    pkgs = [
        ("python", f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"),
        ("numpy", np.__version__),
        ("pandas", pd.__version__),
        ("scanpy", sc.__version__),
        ("anndata", ad.__version__),
        ("scvi-tools", scvi.__version__),
        ("torch", torch.__version__),
        ("seaborn", sns.__version__),
        ("scrublet", getattr(scr, "__version__", "NA")),
        ("celltypist", getattr(celltypist, "__version__", "NA")),
    ]
    print("\n[ENV VERSIONS]")
    for k, v in pkgs:
        print(f"  {k:<12s} {v}")
    print("")

def _download_or_use(path_or_url: str, workdir: Path) -> Path:
    """
    If path_or_url is a local path, return it.
    If it looks like a URL, error out (offline script).
    """
    p = str(path_or_url)
    if re.match(r"^https?://", p):
        raise RuntimeError(
            "A URL was provided, but this script is offline by design. "
            "Please download the Zenodo zip manually and pass --zenodo-zip /local/path.zip"
        )
    local = Path(p).expanduser().resolve()
    if not local.exists():
        raise FileNotFoundError(f"Could not find: {local}")
    return local

def _extract_if_zip(zip_or_dir: Path, tmp_root: Path) -> Path:
    """
    If zip_or_dir is a zip, extract into tmp_root and return extracted directory path.
    If it's a directory, return as-is.
    """
    if zip_or_dir.is_dir():
        return zip_or_dir

    if zipfile.is_zipfile(zip_or_dir):
        extract_dir = tmp_root / "zenodo_extract"
        extract_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_or_dir, "r") as zf:
            zf.extractall(extract_dir)
        # Try to locate the top-level folder 'totalvi_celltypist_inputs'
        candidates = list(extract_dir.glob("totalvi_celltypist_inputs"))
        if not candidates:
            matches = list(extract_dir.rglob("totalvi_celltypist_inputs"))
            if matches:
                return matches[0]
            return extract_dir
        return candidates[0]

    raise ValueError(f"Not a directory or zip: {zip_or_dir}")

def _find_cellbender_h5(sample_dir: Path) -> Path:
    """
    Find cellbender h5 file in sample_dir matching known patterns.
    """
    patterns = [
        "cellbender_output_*_filtered.h5",
        "cellbender_*_filtered.h5",
        "*.h5"
    ]
    for pat in patterns:
        hits = sorted(sample_dir.glob(pat))
        if hits:
            return hits[0]
    raise FileNotFoundError(f"No cellbender .h5 found in {sample_dir}")

def _ensure_figdir(base_out: Path) -> Path:
    figs = base_out / "figures"
    base_out.mkdir(parents=True, exist_ok=True)
    figs.mkdir(parents=True, exist_ok=True)
    return figs

###############################################################################
# Plot helpers still used
###############################################################################

def plot_cell_cycle_proportions(adata, outdir: Path):
    adata.obs["NTX_group"] = adata.obs["sample_id"].map(
        lambda x: "NTX_1_2_3" if x in ["NTX_1","NTX_2","NTX_3"] else "NTX_4_5"
    )
    adata.obs["merged_label"] = [
        f"{ct}_{grp}"
        for ct,grp in zip(adata.obs["cell_type"], adata.obs["NTX_group"])
    ]
    if "cell_cycle_stage" not in adata.obs.columns:
        print("[WARN] cell_cycle_stage not found => skip.\n")
        return

    df = (
        adata.obs
        .groupby(["merged_label","cell_cycle_stage"])
        .size()
        .reset_index(name="count")
    )
    total_by_label = df.groupby("merged_label")["count"].transform("sum")
    df["proportion"] = df["count"]/total_by_label

    g = sns.catplot(
        data=df,
        x="merged_label",
        y="proportion",
        hue="cell_cycle_stage",
        kind="bar",
        height=5, aspect=2
    )
    g.set_xticklabels(rotation=90)
    g.set_axis_labels("Merged Label (CellType_NTXgroup)", "Proportion")
    plt.tight_layout()
    fpath = outdir / "cell_cycle_proportions_merged_label.png"
    plt.savefig(fpath, dpi=300)
    plt.close()
    print(f"[INFO] => {fpath}\n")

###############################################################################
# Main pipeline
###############################################################################

def main():
    parser = argparse.ArgumentParser(description="Reproducible TOTALVI + CellTypist pipeline (Zenodo-only inputs).")
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument("--zenodo-zip", type=str, help="Path to totalvi_celltypist_inputs.zip")
    g.add_argument("--zenodo-dir", type=str, help="Path to extracted totalvi_celltypist_inputs directory")
    parser.add_argument("--outdir", type=str, default="outputs", help="Directory to write outputs into")
    parser.add_argument("--max-epochs", type=int, default=150, help="TOTALVI training epochs")
    parser.add_argument("--seed", type=int, default=0, help="Random seed")
    parser.add_argument("--celltypist-model", type=str, default="Immune_All_Low.pkl",
                        help="Which bundled CellTypist model to use (filename in models/)")
    args = parser.parse_args()

    # Seeds
    np.random.seed(args.seed)
    scvi.settings.seed = args.seed
    print(f"Seed set to {scvi.settings.seed}\n")

    _log_versions()

    workdir = Path.cwd()
    outdir = (workdir / args.outdir).resolve()
    figs_dir = _ensure_figdir(outdir)
    sc.settings.figdir = figs_dir.as_posix()
    sc.settings.autoshow = False

    with tempfile.TemporaryDirectory() as tmpd:
        tmp_root = Path(tmpd)
        if args.zenodo_zip:
            zip_path = _download_or_use(args.zenodo_zip, tmp_root)
            zenodo_root = _extract_if_zip(zip_path, tmp_root)
        else:
            zenodo_root = _extract_if_zip(Path(args.zenodo_dir).expanduser().resolve(), tmp_root)

        # Expect this folder
        if zenodo_root.name != "totalvi_celltypist_inputs":
            candidate = list(Path(zenodo_root).rglob("totalvi_celltypist_inputs"))
            if candidate:
                zenodo_root = candidate[0]
        print(f"[INFO] Using Zenodo root: {zenodo_root}")

        # Paths within the archive
        cellbender_root = zenodo_root / "cellbender_h5"
        model_root      = zenodo_root / "models"
        adt_ref_path    = zenodo_root / "adt_reference" / "TotalSeq_C_Human_Universal_Cocktail_399905_Antibody_reference_UMI_counting.csv"

        # Load ADT reference (optional)
        try:
            adt_ref_df = pd.read_csv(adt_ref_path)
            print(f"ADT reference loaded with {adt_ref_df.shape[0]} entries.\n")
        except Exception as e:
            print(f"Warning: Unable to load ADT reference file.\nError: {e}\n")

        # Sample discovery based on NTX_1..NTX_5 folders
        samples = {}
        for s in ["NTX_1","NTX_2","NTX_3","NTX_4","NTX_5"]:
            sdir = cellbender_root / s
            if not sdir.exists():
                print(f"[WARN] Missing folder {sdir}; skipping.")
                continue
            h5 = _find_cellbender_h5(sdir)
            samples[s] = h5

        if not samples:
            raise RuntimeError("No cellbender .h5 files were found in the Zenodo archive.")

        # Build per-sample RNA + ADT
        adata_list = []
        for sample_id, cb_h5_path in samples.items():
            print(f"\n--- Processing {sample_id} ---")
            adata = sc.read_10x_h5(cb_h5_path.as_posix(), gex_only=False)
            adata.var_names_make_unique()

            rna = adata[:, adata.var["feature_types"]=="Gene Expression"].copy()
            adt = adata[:, adata.var["feature_types"]=="Antibody Capture"].copy()
            rna.layers["counts"] = rna.X.copy()
            rna.obs_names = [f"{sample_id}_{x}" for x in rna.obs_names]

            if adt.n_vars > 0:
                adt.obs_names = [f"{sample_id}_{x}" for x in adt.obs_names]
                adt = adt[rna.obs_names, :].copy()
                adt.layers["counts"] = adt.X.copy()
                adt.var_names_make_unique()
                if issparse(adt.X):
                    rna.obsm["protein_expression"] = adt.X.toarray().astype(np.float32)
                else:
                    rna.obsm["protein_expression"] = adt.X.copy().astype(np.float32)
                rna.uns["protein_names"] = adt.var_names.tolist()
            else:
                rna.obsm["protein_expression"] = np.zeros((rna.n_obs,0), dtype=np.float32)
                rna.uns["protein_names"] = []

            rna.obs["sample_id"] = sample_id
            rna.obs["individual_id"] = "H"
            adata_list.append(rna)

        # Concatenate & align proteins
        all_protein_names = list(set().union(*(a.uns["protein_names"] for a in adata_list)))
        for a in adata_list:
            new_expr = np.zeros((a.n_obs, len(all_protein_names)), dtype=np.float32)
            pn = a.uns["protein_names"]
            if pn:
                idxs = [all_protein_names.index(x) for x in pn]
                new_expr[:, idxs] = a.obsm["protein_expression"]
            a.obsm["protein_expression"] = new_expr
            a.uns["protein_names"] = all_protein_names

        adata_combined = sc.concat(
            adata_list,
            join="outer",
            merge="same",
            uns_merge="unique",
            label="sample_id",
            keys=list(samples.keys())
        )
        print(f"\n[INFO] Combined shape: {adata_combined.shape}\n")

        # Preprocessing
        adata_combined.raw = adata_combined.copy()
        adata_combined.var["mt"] = adata_combined.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata_combined, qc_vars=["mt"], log1p=False, inplace=True)
        adata_combined = adata_combined[adata_combined.obs["total_counts"]>500, :]
        adata_combined = adata_combined[adata_combined.obs["pct_counts_mt"]<10, :]
        sc.pp.filter_genes(adata_combined, min_cells=3)

        # Scrublet (use counts if available)
        scrub_counts = adata_combined.layers.get("counts", adata_combined.X.copy())
        if issparse(scrub_counts):
            scrub_counts = scrub_counts.toarray()
        scrub = scr.Scrublet(scrub_counts, expected_doublet_rate=0.06)
        d_scores, d_preds = scrub.scrub_doublets()
        adata_combined.obs["doublet_score"] = d_scores
        adata_combined.obs["predicted_doublet"] = d_preds
        adata_combined = adata_combined[~adata_combined.obs["predicted_doublet"]].copy()

        # Filter Hemoglobin
        hb_genes = ["HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"]
        hb_in_data = [g for g in hb_genes if g in adata_combined.var_names]
        tot_sum = np.asarray(adata_combined.X.sum(axis=1)).flatten()
        hb_sum = np.asarray(adata_combined[:, hb_in_data].X.sum(axis=1)).flatten() if hb_in_data else np.zeros_like(tot_sum)
        adata_combined.obs["pct_hb"] = (hb_sum / np.clip(tot_sum, 1e-12, None))*100
        adata_combined = adata_combined[adata_combined.obs["pct_hb"]<1.0, :]

        # Normalize & log
        sc.pp.normalize_total(adata_combined, target_sum=1e4)
        sc.pp.log1p(adata_combined)
        adata_combined.layers["log1p"] = adata_combined.X.copy()
        adata_combined.X = adata_combined.X.astype(np.float32)

        # Cell Cycle
        S_GENES = [
            "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
            "CCNA2","HELLS","RFC2","RPA2","NASP","RAD51AP1","GMNN","WDR76","SLBP",
            "UBR7","POLD3","MSH2","ATAD2","RAD51","RRM2","CDC45","CDC6","EXO1","TIPIN",
            "DSCC1","BLM","CASP8AP2","USP1","CLSPN","POLA1","CHAF1B","BRIP1","E2F8"
        ]
        G2M_GENES = [
            "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
            "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","CKAP2",
            "AURKB","BUB1","KIF11","ANP32E","TUBB4B","GTSE1","KIF20B","HJURP","CDCA3",
            "CDC20","TTK","CDC25C","KIF2C","RANGAP1","NASP","BUB1B","CENPE","CDCA2",
            "CDCA8","ECT2","KIF23","HMMR","AURKA","PSRC1","ANLN","LBR","CKAP2L","CENPA"
        ]
        valid_s = [g for g in S_GENES if g in adata_combined.var_names]
        valid_g2m = [g for g in G2M_GENES if g in adata_combined.var_names]
        if valid_s and valid_g2m:
            sc.tl.score_genes_cell_cycle(adata_combined, s_genes=valid_s, g2m_genes=valid_g2m)
            adata_combined.obs["cell_cycle_stage"] = adata_combined.obs["phase"]

        # TOTALVI
        print("[INFO] Setting up TOTALVI...\n")
        scvi.model.TOTALVI.setup_anndata(
            adata_combined,
            layer="counts",
            protein_expression_obsm_key="protein_expression",
            batch_key="sample_id"
        )
        model = TOTALVI(
            adata_combined,
            n_layers_encoder=2,
            n_layers_decoder=2,
            n_hidden=384,
            n_latent=30,
            gene_likelihood="nb"
        )
        print("[INFO] Training TOTALVI...\n")
        model.train(
            max_epochs=args.max_epochs, lr=1e-3,
            plan_kwargs={
                "weight_decay":1e-6,
                "reduce_lr_on_plateau":True,
                "lr_factor":0.5,
                "lr_patience":20,
                "lr_scheduler_metric":"elbo_validation"
            },
            early_stopping=True,
            early_stopping_patience=15,
            accelerator="gpu" if torch.cuda.is_available() else "cpu",
            devices=1
        )

        adata_combined.obsm["X_totalVI"] = model.get_latent_representation()
        sc.pp.neighbors(adata_combined, use_rep="X_totalVI")
        sc.tl.umap(adata_combined)
        sc.tl.leiden(adata_combined, resolution=1.0)

        # UMAPs
        sc.pl.umap(
            adata_combined,
            color="leiden",
            legend_loc="on data",
            title="UMAP (TOTALVI) - Clusters on data",
            save="_clusters_onmap.png",
            show=False
        )
        sc.pl.umap(
            adata_combined,
            color="sample_id",
            title="UMAP (TOTALVI) - sample",
            save="_sample.png",
            show=False
        )
        if "cell_cycle_stage" in adata_combined.obs.columns:
            sc.pl.umap(
                adata_combined,
                color="cell_cycle_stage",
                title="UMAP (TOTALVI) by Cell Cycle",
                save="_cell_cycle_stage.png",
                show=False
            )

        # Denoised protein (arcsinh) — computed for completeness (not plotted here)
        all_batches = adata_combined.obs["sample_id"].unique().tolist()
        rna_samps, prot_samps = model.get_normalized_expression(
            adata=adata_combined,
            return_numpy=True, n_samples=10,
            library_size=1e4,
            scale_protein=True,
            transform_batch=all_batches,
            return_mean=False
        )
        denoised_prot_mean = prot_samps.mean(axis=-1)
        denoised_prot_std = prot_samps.std(axis=-1)
        adata_combined.obsm["denoised_protein"] = denoised_prot_mean
        adata_combined.obsm["denoised_protein_std"] = denoised_prot_std
        adata_combined.obsm["denoised_protein_asinh"] = np.arcsinh(denoised_prot_mean/5)

        prot_names = adata_combined.uns.get("protein_names", [])
        for i, p in enumerate(prot_names):
            adata_combined.obs[f"{p}_denoised"] = denoised_prot_mean[:, i]
            adata_combined.obs[f"{p}_denoised_std"] = denoised_prot_std[:, i]
            adata_combined.obs[f"{p}_denoised_arcsin"] = np.arcsinh(denoised_prot_mean[:, i]/5)

        # CellTypist using bundled model (no download)
        print("[INFO] Annotating with CellTypist (bundled model)...\n")
        model_path = (model_root / args.celltypist_model)
        if not model_path.exists():
            raise FileNotFoundError(f"CellTypist model not found at {model_path}.")
        adata_ct = adata_combined.copy()
        adata_ct.X = adata_ct.layers["log1p"]
        predictions = celltypist.annotate(
            adata_ct, model=model_path.as_posix(),
            majority_voting=True, mode="best match"
        )
        adata_ct = predictions.to_adata()
        adata_combined.obs["cell_type"] = adata_ct.obs["majority_voting"]
        adata_combined.obs["predicted_labels"] = adata_ct.obs["predicted_labels"]
        adata_combined.obs["cell_type_confidence"] = adata_ct.obs["conf_score"]

        cluster2type = {}
        for c in adata_combined.obs["leiden"].unique():
            top_type = (
                adata_combined.obs
                .loc[adata_combined.obs["leiden"]==c, "cell_type"]
                .value_counts()
                .idxmax()
            )
            cluster2type[c] = top_type
        adata_combined.obs["cluster_cell_type"] = adata_combined.obs["leiden"].map(cluster2type)

        sc.pl.umap(
            adata_combined,
            color="cluster_cell_type",
            save="_cluster_cell_type.png",
            show=False
        )
        sc.pl.umap(
            adata_combined,
            color="cell_type_confidence",
            color_map="viridis",
            save="_celltypist_confidence.png",
            show=False
        )

        # Save AnnData
        out_h5ad = outdir / "processed_adata_NTX_all_samples.h5ad"
        print("\n[INFO] Writing final AnnData to disk...\n")
        adata_combined.write(out_h5ad.as_posix())

        # Cell cycle proportions
        print("[INFO] Plot cell cycle proportions by merged_label...\n")
        plot_cell_cycle_proportions(adata_combined, figs_dir)

        print("\nDone. Outputs are under:", outdir.as_posix())

if __name__ == "__main__":
    main()
