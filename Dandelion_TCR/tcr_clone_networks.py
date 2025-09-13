#!/usr/bin/env python3
# ========================================================================
#  T‑cell clone networks + diversity (TOTALVI UMAP overlay) – path‑agnostic
# ========================================================================
"""
tcr_clone_networks.py
=====================

• Builds clone networks from productive single‑pair TCRs
• Calculates Gini diversity and renders network / UMAP plots

Run example
-----------
python scripts/tcr_clone_networks.py \
    --contigs  /PATH/TO/NTX_2_TCR_fastq/dandelion/all_contig_dandelion.tsv \
               /PATH/TO/NTX_3_TCR_fastq/dandelion/all_contig_dandelion.tsv \
               /PATH/TO/NTX_4_TCR_fastq/dandelion/all_contig_dandelion.tsv \
               /PATH/TO/NTX_5_TCR_fastq/dandelion/all_contig_dandelion.tsv \
    --adata    results/totalvi_celltypist/processed_adata_NTX_all_samples.h5ad \
    --out-dir  results/tcr_clone_networks
"""
from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import List

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import dandelion as ddl
from dandelion.utilities import read_airr

# ──────────────────────────────────────────────────────────────────────
# Constants (adjust PRE/POST labels if your sample naming differs)
# ──────────────────────────────────────────────────────────────────────
_PREFIX_RE = re.compile(r"(NTX_\d+)")
UMAP_KEYS = ("X_totalVI_umap", "X_umap_totalvi", "X_umap")

PRE = {"NTX_1", "NTX_2", "NTX_3"}
POST = {"NTX_4", "NTX_5"}

MIN_CLONE_SZ = 3
NODE_SCALE = 5
CRIMSON = "crimson"


# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
def detect_prefix(path: Path) -> str:
    m = _PREFIX_RE.search(str(path))
    return m.group(1) if m else "UNK"


def load_vdj(paths: List[Path], logger: logging.Logger) -> ddl.Dandelion:
    parts = [read_airr(str(p), prefix=detect_prefix(p)) for p in paths]
    for p, v in zip(paths, parts):
        logger.info("Loaded %-34s → %5d contigs", p.name, v.n_contigs)
    return ddl.concat(parts)


def ensure_umap(a: ad.AnnData, logger: logging.Logger) -> str:
    for k in UMAP_KEYS:
        if k in a.obsm and a.obsm[k].shape[1] == 2:
            return k
    logger.info("UMAP not found – computing from TOTALVI latent space")
    sc.pp.neighbors(a, use_rep="X_totalVI", n_neighbors=15, key_added="totalVI")
    sc.tl.umap(a, min_dist=0.3, random_state=0)
    a.obsm["X_totalVI_umap"] = a.obsm["X_umap"].copy()
    return "X_totalVI_umap"


def add_clone_size(a: ad.AnnData):
    vc = a.obs["clone_id"].value_counts()
    a.obs["clone_size"] = a.obs["clone_id"].map(vc).fillna(0).astype(int)


def size_legend(ax, scale=NODE_SCALE):
    for lbl, n in (("1 cell", 1), ("5 cells", 5), ("≥10 cells", 10)):
        ax.scatter([], [], s=n * scale, c=CRIMSON, label=lbl)
    ax.legend(
        title="Clone size",
        frameon=False,
        labelspacing=1,
        scatterpoints=1,
        fontsize="small",
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
    )


def guarantee_clone_size(vdj: ddl.Dandelion):
    if "clone_size" not in vdj.metadata:
        if "clone_id_size" not in vdj.metadata:
            ddl.tl.clone_size(vdj)
        vdj.metadata.rename(columns={"clone_id_size": "clone_size"}, inplace=True)


def dandelion_clonality(vdj: ddl.Dandelion) -> float:
    guarantee_clone_size(vdj)
    if vdj.n_contigs == 0:
        return np.nan
    vdj.metadata["__all"] = "all"
    t = ddl.tl.clone_diversity(vdj, groupby="__all", method="gini", metric="clone_size", return_table=True)
    return float(t.iat[0, 0]) if not t.empty else np.nan


def safe(txt: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", txt).strip("_")


def plot_clone_network(adata_sub: ad.AnnData, ax: plt.Axes):
    sizes = adata_sub.obs["clone_size"].fillna(1).astype(float) * NODE_SCALE
    ddl.pl.clone_network(
        adata_sub,
        edges=False,
        size=sizes.values,
        show=False,
        ax=ax,
    )
    for coll in ax.collections[:1]:
        coll.set_facecolor(CRIMSON)
        coll.set_edgecolor("none")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, description="Clone networks & diversity metrics"
    )
    p.add_argument("--contigs", type=Path, nargs="+", required=True, help="all_contig_dandelion.tsv files")
    p.add_argument("--adata", type=Path, required=True, help="Input AnnData (h5ad)")
    p.add_argument("--out-dir", type=Path, required=True, help="Directory for figures & CSVs")
    p.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING"], default="INFO")
    return p


def main() -> None:
    args = build_parser().parse_args()
    out_dir: Path = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "clone_size_gini.csv"

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s │ %(levelname)-8s │ %(message)s",
        datefmt="%H:%M:%S",
    )
    log = logging.getLogger(__name__)
    log.info("▶ Clone‑network pipeline started")

    # 1 ░ contigs ➜ clones
    vdj = load_vdj(args.contigs, log)
    vdj.threshold = 0.85
    ddl.tl.find_clones(vdj, key="junction_aa")
    vdj = vdj[vdj.metadata.chain_status == "Single pair"]
    guarantee_clone_size(vdj)
    ddl.tl.generate_network(vdj, key="junction_aa", min_size=MIN_CLONE_SZ, layout_method="mod_fr")

    # 2 ░ AnnData
    adata = ad.read_h5ad(args.adata)
    umap_key = ensure_umap(adata, log)
    ddl.tl.transfer(adata, vdj)
    add_clone_size(adata)

    # 3 ░ UMAP overlay for single‑pair T cells
    sp = adata.obs["chain_status"] == "Single pair"
    t_labels = [
        l
        for l in adata.obs["cell_type"].unique().tolist()
        if l.lower() == "mait cells" or "t cells" in l.lower()
    ]
    X = adata.obsm[umap_key]
    cmap = plt.colormaps["tab20"]
    plt.figure(figsize=(6, 6), dpi=300)
    plt.scatter(X[:, 0], X[:, 1], s=1, c="lightgrey", alpha=0.2)
    for i, lab in enumerate(sorted(t_labels)):
        sel = sp & (adata.obs["cell_type"] == lab)
        plt.scatter(X[sel, 0], X[sel, 1], s=4, c=[cmap(i)], alpha=0.8, label=lab)
    plt.title("Single‑pair T‑cell subsets")
    plt.xlabel("UMAP‑1")
    plt.ylabel("UMAP‑2")
    plt.legend(markerscale=2, frameon=False, fontsize=6)
    plt.tight_layout()
    plt.savefig(out_dir / "totalvi_tcell_subsets_overlay.png")
    plt.close()

    rows = [["Label", "Phase", "Gini"]]

    # 4 ░ combined
    g_all = dandelion_clonality(vdj)
    rows.append(["Combined", "All", f"{g_all:.4f}" if not np.isnan(g_all) else "NA"])
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300, facecolor="white")
    plot_clone_network(adata, ax)
    size_legend(ax)
    ax.set_title("Combined")
    plt.tight_layout(rect=[0, 0, 0.82, 1.0])
    plt.savefig(out_dir / "clone_network_all_samples.png")
    plt.close(fig)

    # helper
    def do(label: str, phase: str, mask):
        idx = adata.obs_names[mask]
        if idx.empty:
            return
        g = dandelion_clonality(vdj[vdj.metadata.index.isin(idx)].copy())
        rows.append([label, phase, f"{g:.4f}"])
        fig_, ax_ = plt.subplots(figsize=(6, 6), dpi=300, facecolor="white")
        plot_clone_network(adata[idx], ax_)
        size_legend(ax_)
        ax_.set_title(f"{label} ({phase})")
        plt.tight_layout(rect=[0, 0, 0.82, 1.0])
        plt.savefig(out_dir / f"clone_network_{safe(label)}_{phase}.png")
        plt.close(fig_)
        log.info("Saved network – %s (%s)", label, phase)

    # 5 ░ combined T‑cells pre/post
    tmask = adata.obs["cell_type"].isin(t_labels)
    do("Combined T‑cells", "PRE", tmask & adata.obs["sample_id"].isin(PRE))
    do("Combined T‑cells", "POST", tmask & adata.obs["sample_id"].isin(POST))

    # 6 ░ each subset
    for lab in t_labels:
        sel = adata.obs["cell_type"] == lab
        do(lab, "PRE", sel & adata.obs["sample_id"].isin(PRE))
        do(lab, "POST", sel & adata.obs["sample_id"].isin(POST))

    # 7 ░ write CSV
    import csv as _csv
    with csv_path.open("w", newline="") as f:
        _csv.writer(f).writerows(rows)
    log.info("Wrote Gini CSV → %s", csv_path)
    log.info("✓ DONE")


if __name__ == "__main__":
    main()
