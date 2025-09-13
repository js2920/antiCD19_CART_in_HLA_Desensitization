#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TOTALVI + CellTypist pipeline (template)
=======================================

Path-agnostic template for the NTX dataset (GEO: GSE307464).

What this script does
---------------------
1. Reads CellBender-filtered 10x HDF5 files (GEX + ADT in the same H5 if present).
2. QC → Scrublet doublet removal → mito / haemoglobin filters.
3. Trains TOTALVI on RNA+ADT; computes UMAP / Leiden.
4. CellTypist annotation + cluster-to-type mapping.
5. Summary figures (UMAPs, cell-cycle proportions).
6. Writes a processed AnnData.

Inputs
------
• A data directory with per-sample subfolders or files containing CellBender-filtered H5:
  expected patterns (anywhere under --data-dir):
    - cellbender_output_<SAMPLE>_filtered.h5
    - cellbender_<SAMPLE>_filtered.h5
• (Optional) TotalSeq C feature reference CSV

Outputs
-------
• <out_dir>/processed_adata_NTX_all_samples.h5ad
• Figures under <out_dir>/figures/
• Logs printed to stdout

Example
-------
python scripts/01_totalvi_celltypist.py \
  --data-dir /path/to/GSE307464_download \
  --samples NTX_1,NTX_2,NTX_3,NTX_4,NTX_5 \
  --out-dir results/totalvi_celltypist

Requirements
------------
scanpy, anndata, numpy, pandas, matplotlib, seaborn, scrublet, scvi-tools, celltypist, torch
"""

# ── Imports ───────────────────────────────────────────────────────────────────
import os
import re
import sys
import glob
import warnings
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")  # non-interactive back-end
import matplotlib.pyplot as plt
import seaborn as sns
import torch
import anndata as ad

from scipy.sparse import issparse
import scrublet as scr

import scvi
from scvi.model import TOTALVI

import celltypist
from celltypist import models

warnings.filterwarnings("ignore")
np.random.seed(0)
scvi.settings.seed = 0


# ── CLI ───────────────────────────────────────────────────────────────────────
def get_args():
    ap = argparse.ArgumentParser(
        description="TOTALVI + CellTypist pipeline template (GEO: GSE307464)"
    )
    ap.add_argument(
        "--data-dir", required=True,
        help="Root folder that contains per-sample files/subfolders (downloaded from GEO)."
    )
    ap.add_argument(
        "--samples", default="NTX_1,NTX_2,NTX_3,NTX_4,NTX_5",
        help="Comma-separated sample IDs to include (default: NTX_1..NTX_5)."
    )
    ap.add_argument(
        "--feature-ref", default=None,
        help="(Optional) Path to TotalSeq C feature reference CSV."
    )
    ap.add_argument(
        "--out-dir", default="results/totalvi_celltypist",
        help="Output directory (H5AD + figures)."
    )
    ap.add_argument("--min-counts", type=int, default=500,
                    help="Min UMI counts per cell (default: 500).")
    ap.add_argument("--max-mt", type=float, default=10.0,
                    help="Max %% MT counts per cell (default: 10).")
    ap.add_argument("--max-hb", type=float, default=1.0,
                    help="Max %% haemoglobin gene counts (default: 1).")
    ap.add_argument("--doublet-rate", type=float, default=0.06,
                    help="Expected doublet rate for Scrublet (default: 0.06).")
    return ap.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────
def find_cellbender_h5(data_dir, sample_id):
    """
    Search recursively under data_dir for a CellBender-filtered H5 for `sample_id`.
    Accepts flexible filename patterns to accommodate deposit differences.
    """
    patterns = [
        f"**/cellbender_output_{sample_id}_filtered.h5",
        f"**/cellbender_{sample_id}_filtered.h5",
        f"**/{sample_id}/cellbender_output_{sample_id}_filtered.h5",
        f"**/{sample_id}/cellbender_{sample_id}_filtered.h5",
    ]
    for pat in patterns:
        hits = glob.glob(os.path.join(data_dir, pat), recursive=True)
        if hits:
            return sorted(hits)[0]
    return None


def arcsinh(x, cofactor=5.0):
    return np.arcsinh(x / cofactor)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    args = get_args()
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    out_dir = args.out_dir
    fig_dir = os.path.join(out_dir, "figures")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    sc.settings.autoshow = False

    print(f"[INFO] Seed = {scvi.settings.seed}")
    print(f"[INFO] Samples: {samples}")
    print(f"[INFO] Data dir: {args.data_dir}")
    print(f"[INFO] Output: {out_dir}")

    # 1) Read optional ADT reference (book-keeping only)
    if args.feature_ref and os.path.exists(args.feature_ref):
        try:
            adt_ref_df = pd.read_csv(args.feature_ref)
            print(f"[INFO] ADT reference → {adt_ref_df.shape[0]} antibodies")
        except Exception as e:
            print(f"[WARN] ADT reference failed to read: {e}")
            adt_ref_df = pd.DataFrame()
    else:
        adt_ref_df = pd.DataFrame()

    # 2) Load & process each sample
    adata_list = []
    for sid in samples:
        print(f"\n── Processing {sid} ──")
        cb_h5 = find_cellbender_h5(args.data_dir, sid)
        if cb_h5 is None:
            sys.exit(f"[ERROR] No CellBender file found for {sid}. "
                     f"Looked under {args.data_dir} with multiple patterns.")
        print(f"[INFO] Using: {cb_h5}")
        adata = sc.read_10x_h5(cb_h5, gex_only=False)
        adata.var_names_make_unique()

        # split modalities
        rna = adata[:, adata.var["feature_types"] == "Gene Expression"].copy()
        adt = adata[:, adata.var["feature_types"] == "Antibody Capture"].copy()

        # store raw counts for RNA
        rna.layers["counts"] = rna.X.copy()

        # give unique, sample‑prefixed cell IDs
        rna.obs_names = [f"{sid}_{x}" for x in rna.obs_names]

        # place proteins into obsm as a float32 dense array
        if adt.n_vars:
            adt.obs_names = rna.obs_names
            mat = adt.X.toarray() if issparse(adt.X) else adt.X.copy()
            rna.obsm["protein_expression"] = mat.astype(np.float32)
            rna.uns["protein_names"] = adt.var_names.tolist()
        else:
            rna.obsm["protein_expression"] = np.zeros((rna.n_obs, 0), dtype=np.float32)
            rna.uns["protein_names"] = []

        rna.obs["sample_id"] = sid
        adata_list.append(rna)

    # 3) Harmonise ADT space & concatenate
    all_prots = list(set().union(*(a.uns["protein_names"] for a in adata_list)))
    all_prots = sorted(all_prots)
    for a in adata_list:
        mat = np.zeros((a.n_obs, len(all_prots)), dtype=np.float32)
        if a.uns["protein_names"]:
            idx = [all_prots.index(p) for p in a.uns["protein_names"]]
            mat[:, idx] = a.obsm["protein_expression"]
        a.obsm["protein_expression"] = mat
        a.uns["protein_names"] = all_prots

    adata = sc.concat(
        adata_list, join="outer", merge="same",
        uns_merge="unique", label="sample_id", keys=samples
    )
    print(f"[INFO] Combined AnnData shape: {adata.shape}")

    # 4) QC + Scrublet + filters
    adata.raw = adata.copy()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], log1p=False, inplace=True)

    # cell-level filters
    adata = adata[adata.obs["total_counts"] > args.min_counts, :]
    adata = adata[adata.obs["pct_counts_mt"] < args.max_mt, :]

    # gene-level filter
    sc.pp.filter_genes(adata, min_cells=3)

    # Scrublet doublets (on counts layer if available)
    scrub_counts = adata.layers.get("counts", adata.X)
    scrub_counts = scrub_counts.toarray() if issparse(scrub_counts) else scrub_counts
    scrub = scr.Scrublet(scrub_counts, expected_doublet_rate=args.doublet_rate)
    scrub_scores, scrub_preds = scrub.scrub_doublets()
    adata.obs["doublet_score"] = scrub_scores
    adata.obs["predicted_doublet"] = scrub_preds
    adata = adata[~adata.obs["predicted_doublet"], :]

    # Haemoglobin filter
    hb_genes = ["HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ"]
    tot = adata.X.sum(axis=1).A1 if issparse(adata.X) else adata.X.sum(1)
    hb = adata[:, [g for g in hb_genes if g in adata.var_names]].X
    hb = hb.sum(axis=1).A1 if issparse(hb) else hb.sum(1)
    adata.obs["pct_hb"] = 100 * hb / np.maximum(tot, 1e-8)
    adata = adata[adata.obs["pct_hb"] < args.max_hb, :]

    # normalise / log1p for visualisation & DE downstream (RNA)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()
    adata.X = adata.X.astype(np.float32)

    # Cell cycle scoring (graceful fallback if lists incomplete)
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
    valid_s = [g for g in S_GENES if g in adata.var_names]
    valid_g2m = [g for g in G2M_GENES if g in adata.var_names]
    if valid_s and valid_g2m:
        sc.tl.score_genes_cell_cycle(adata, s_genes=valid_s, g2m_genes=valid_g2m)
        adata.obs["cell_cycle_stage"] = adata.obs["phase"]
    else:
        adata.obs["cell_cycle_stage"] = "NA"

    # 5) TOTALVI + UMAP/Leiden
    print("[INFO] Training TOTALVI ...")
    scvi.model.TOTALVI.setup_anndata(
        adata,
        layer="counts",
        protein_expression_obsm_key="protein_expression",
        batch_key="sample_id",
    )
    model = TOTALVI(
        adata,
        n_layers_encoder=2,
        n_layers_decoder=2,
        n_hidden=384,
        n_latent=30,
        gene_likelihood="nb",
    )
    model.train(
        max_epochs=150,
        lr=1e-3,
        plan_kwargs=dict(
            weight_decay=1e-6,
            reduce_lr_on_plateau=True,
            lr_factor=0.5,
            lr_patience=20,
            lr_scheduler_metric="elbo_validation",
        ),
        early_stopping=True,
        early_stopping_patience=15,
    )

    adata.obsm["X_totalVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_totalVI")
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=1.0)

    sc.pl.umap(adata, color="leiden", legend_loc="on data",
               save="_clusters.png", show=False)
    sc.pl.umap(adata, color="sample_id", save="_sample.png", show=False)
    if "cell_cycle_stage" in adata.obs:
        sc.pl.umap(adata, color="cell_cycle_stage",
                   save="_cell_cycle_stage.png", show=False)

    # 6) CellTypist annotation
    print("[INFO] CellTypist annotation ...")
    models.download_models(force_update=False)
    ct_pred = celltypist.annotate(
        adata.copy(),
        model="Immune_All_Low.pkl",
        majority_voting=True,
        mode="best match",
    ).to_adata()
    adata.obs["cell_type"] = ct_pred.obs["majority_voting"]
    adata.obs["predicted_labels"] = ct_pred.obs["predicted_labels"]
    adata.obs["cell_type_confidence"] = ct_pred.obs["conf_score"]

    # majority cell-type per Leiden cluster
    cluster2type = (
        adata.obs.groupby("leiden")["cell_type"]
        .agg(lambda x: x.value_counts().idxmax())
        .to_dict()
    )
    adata.obs["cluster_cell_type"] = adata.obs["leiden"].map(cluster2type)

    sc.pl.umap(adata, color="cluster_cell_type",
               save="_cluster_cell_type.png", show=False)
    sc.pl.umap(adata, color="cell_type_confidence", cmap="viridis",
               save="_celltypist_conf.png", show=False)

    # 7) Cell‑cycle proportions (kept)
    def plot_cell_cycle_proportions(adata, outdir):
        if "cell_cycle_stage" not in adata.obs:
            return
        adata = adata.copy()
        adata.obs["NTX_group"] = adata.obs["sample_id"].map(
            lambda x: "NTX_1_2_3" if str(x).startswith(("NTX_1","NTX_2","NTX_3")) else "NTX_4_5"
        )
        adata.obs["merged_label"] = adata.obs["cell_type"] + "_" + adata.obs["NTX_group"]
        df = (adata.obs.groupby(["merged_label","cell_cycle_stage"])
              .size().reset_index(name="n"))
        df["prop"] = df["n"] / df.groupby("merged_label")["n"].transform("sum")
        g = sns.catplot(data=df, x="merged_label", y="prop", hue="cell_cycle_stage",
                        kind="bar", height=5, aspect=2)
        g.set_xticklabels(rotation=90)
        g.set_axis_labels("label", "proportion")
        fp = os.path.join(outdir, "cell_cycle_props.png")
        plt.tight_layout()
        plt.savefig(fp, dpi=300)
        plt.close()
        print("[INFO] →", fp)

    plot_cell_cycle_proportions(adata, fig_dir)

    # 8) Combined dot‑plots (RNA kept; ADT is optional and will no‑op if missing)
    def rename_ct_grp(adata_):
        adata_.obs["NTX_group"] = adata_.obs["sample_id"].map(
            lambda x: "NTX_1_2_3" if str(x).startswith(("NTX_1","NTX_2","NTX_3")) else "NTX_4_5"
        )
        adata_.obs["merged_label"] = adata_.obs["cell_type"] + "_" + adata_.obs["NTX_group"]

    def combined_dotplot_top100(adata_, n=10):
        a = adata_.copy()
        rename_ct_grp(a)
        a.X = a.layers["log1p"]
        a.raw = None
        sc.tl.rank_genes_groups(a, groupby="merged_label", method="wilcoxon")
        r = a.uns["rank_genes_groups"]
        cats = r["names"].dtype.names
        genes = list({g for c in cats for g in r["names"][c][:n]})[:100]
        sc.pl.dotplot(
            a, var_names=genes, groupby="merged_label",
            standard_scale="var", dendrogram=False,
            save="_topGenes.png", show=False
        )

    def combined_dotplot_top100_adt(adata_, n=10):
        # Optional: falls back to raw protein counts if available.
        if "protein_expression" not in adata_.obsm:
            return
        X = adata_.obsm["protein_expression"]
        prots = adata_.uns.get("protein_names", [f"ADT_{i}" for i in range(X.shape[1])])
        adt = ad.AnnData(X=np.asarray(X), obs=adata_.obs.copy(), var=pd.DataFrame(index=prots))
        rename_ct_grp(adt)
        sc.tl.rank_genes_groups(adt, groupby="merged_label", method="wilcoxon")
        r = adt.uns["rank_genes_groups"]
        cats = r["names"].dtype.names
        top = list({p for c in cats for p in r["names"][c][:n]})[:100]
        sc.pl.dotplot(
            adt, var_names=top, groupby="merged_label",
            standard_scale="var", dendrogram=False,
            save="_topADTs_raw.png", show=False
        )

    combined_dotplot_top100(adata)
    combined_dotplot_top100_adt(adata)

    # 9) Save H5AD
    out_h5ad = os.path.join(out_dir, "processed_adata_NTX_all_samples.h5ad")
    print("[INFO] Saving AnnData ...")
    adata.write_h5ad(out_h5ad)
    print(f"[DONE] TOTALVI + CellTypist pipeline finished. → {out_h5ad}")


if __name__ == "__main__":
    main()
