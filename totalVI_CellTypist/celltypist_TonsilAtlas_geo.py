#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Refinement of Plasma + Memory-B with Tonsil CellTypist model (template)
======================================================================

Focus module for re-annotation and visualisation on the NTX dataset (GEO: GSE307464).

What it does
------------
• Reads a processed AnnData (output from 01_totalvi_celltypist.py).
• Subsets to Plasma cells + a selected Leiden cluster (default: "5") of memory B.
• Re-annotates with the Tonsil CellTypist model.
• Generates global and per-sample UMAPs.
• Pre vs Post UMAPs (Memory B vs Plasma).
• Writes a refined H5AD.

"""

import os
import argparse
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import celltypist
from celltypist import models


PALETTE_POOL = (
    list(matplotlib.colormaps["tab20b"].colors) +
    list(matplotlib.colormaps["tab20c"].colors) +
    list(matplotlib.colormaps["Set3"].colors)
)


def get_args():
    ap = argparse.ArgumentParser(
        description="Refine Plasma + Memory-B with Tonsil CellTypist (template)"
    )
    ap.add_argument(
        "--input-h5ad", required=True,
        help="Processed AnnData from script 01 (e.g., results/.../processed_adata_NTX_all_samples.h5ad)."
    )
    ap.add_argument(
        "--out-dir", default="results/tonsil_refinement",
        help="Output directory for figures and refined h5ad."
    )
    ap.add_argument(
        "--memory-leiden", default="5",
        help='Leiden cluster to include as memory-B alongside Plasma cells (default: "5").'
    )
    return ap.parse_args()


def main():
    args = get_args()
    os.makedirs(args.out_dir, exist_ok=True)
    sc.settings.figdir = args.out_dir
    sc.settings.autoshow = False

    # 1) Read & subset
    adata = sc.read_h5ad(args.input_h5ad)
    mask = (adata.obs.get("cell_type", pd.Series(index=adata.obs.index)) == "Plasma cells") | \
           (adata.obs.get("leiden", pd.Series(index=adata.obs.index)) == args.memory_leiden)
    sub = adata[mask].copy()
    if "log1p" in sub.layers:
        sub.X = sub.layers["log1p"]

    # 2) CellTypist (Tonsil)
    models.download_models(force_update=False)
    pred = celltypist.annotate(
        sub, model="Cells_Human_Tonsil.pkl",
        majority_voting=True, mode="best match"
    ).to_adata()

    lbl = "Tonsil_CellTypist_Label"
    pred.obs[lbl] = pred.obs["majority_voting"]
    pred.obs["ct_conf"] = pred.obs["conf_score"]

    # Keep the same UMAP embedding from TOTALVI space for consistency when available
    if "X_umap" in adata.obsm:
        pred.obsm["X_umap"] = adata.obsm["X_umap"][mask]
        pred.uns["umap"] = adata.uns.get("umap", {})
    elif "X_totalVI" in adata.obsm:
        pred.obsm["X_totalVI"] = adata.obsm["X_totalVI"][mask]
        sc.pp.neighbors(pred, use_rep="X_totalVI")
        sc.tl.umap(pred)

    # palette
    cats = pred.obs[lbl].astype("category").cat.categories
    pal = dict(zip(cats, PALETTE_POOL[:len(cats)]))
    pal.update({
        "MBC derived IgG+ PC": "#7b3294",
        "MBC derived early PC precursors": "#d73027",
    })

    # 3) UMAPs (global)
    sc.pl.umap(pred, color="leiden", palette=PALETTE_POOL,
               legend_loc="right margin", save="_sub_leiden.png", show=False)
    sc.pl.umap(pred, color="sample_id", save="_sub_sample.png", show=False)
    sc.pl.umap(pred, color=lbl, palette=pal,
               legend_loc="right margin", save="_sub_labels.png", show=False)
    sc.pl.umap(pred, color="ct_conf", cmap="viridis",
               save="_sub_conf.png", show=False)

    # 4) Per-sample UMAPs
    for sid in sorted(pred.obs["sample_id"].astype(str).unique()):
        sub2 = pred[pred.obs["sample_id"].astype(str) == sid]
        if sub2.n_obs:
            sc.pl.umap(sub2, color=lbl, palette=pal,
                       legend_loc="right margin", title=sid,
                       save=f"_{sid}_labels.png", show=False)

    # 5) Pre/Post UMAPs for Memory B vs Plasma
    sel = pred.obs[lbl].isin(["Memory B cells", "Plasma cells"])
    pre = pred[sel & pred.obs["sample_id"].astype(str).str.startswith(("NTX_1","NTX_2","NTX_3"))]
    post = pred[sel & pred.obs["sample_id"].astype(str).str.startswith(("NTX_4","NTX_5"))]
    if pre.n_obs:
        sc.pl.umap(pre, color=lbl, palette=pal, save="_preCAR_memB_PC.png", show=False)
    if post.n_obs:
        sc.pl.umap(post, color=lbl, palette=pal, save="_postCAR_memB_PC.png", show=False)

    # 6) Save refined object
    out_h5ad = os.path.join(args.out_dir, "refined_plasma_memory_b.h5ad")
    pred.write_h5ad(out_h5ad)
    print(f"[DONE] Refinement script finished. → {out_h5ad}")


if __name__ == "__main__":
    main()
