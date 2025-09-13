#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full Memory‑B + Plasma clonotyping & reporting workflow (template).

Outputs are written into --out-dir (default: results/clonotyping).
"""
import os
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import dandelion as ddl
import matplotlib.pyplot as plt

# Optional deps (graceful fallbacks)
try:
    from upsetplot import from_contents, UpSet
except Exception:  # pragma: no cover
    from_contents, UpSet = None, None

# ── CLI ───────────────────────────────────────────────────────────────────────
def parse_args():
    ap = argparse.ArgumentParser(description="B/Plasma clonotyping + reports")
    ap.add_argument("--adata", required=True,
                    help="H5AD with BCR already integrated (from merge_bcr_and_integrate.py).")
    ap.add_argument("--contigs", required=True,
                    help="merged_contigs.tsv (from merge_bcr_and_integrate.py).")
    ap.add_argument("--out-dir", default="results/clonotyping",
                    help="Output directory for CSVs and figures.")
    ap.add_argument("--pre-samples", default="NTX_1,NTX_2,NTX_3")
    ap.add_argument("--post-samples", default="NTX_4,NTX_5")
    return ap.parse_args()


# ── Helpers ───────────────────────────────────────────────────────────────────
def rename_barcode(bc: str) -> str:
    return bc.replace("_fastq_", "_") + "-1" if "_fastq_" in bc else bc


def build_heavy_light_table(vdj_obj: ddl.Dandelion) -> pd.DataFrame:
    cols = ["junction_aa", "cell_id", "locus", "c_call", "c_call_10x",
            "fwr1", "cdr1", "fwr2", "cdr2", "fwr3", "cdr3", "fwr4"]
    cols = [c for c in cols if c in vdj_obj.data.columns]
    rec = []
    for r in vdj_obj.data[cols].itertuples(index=False):
        d = r._asdict()
        for k in ("fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4"):
            if k in d: d[k] = str(d[k]).replace(".", "").replace("-", "").upper()
        rec.append(d)
    df = pd.DataFrame(rec)
    heavy = df[df["locus"].str.upper() == "IGH"].drop_duplicates("cell_id")
    light = df[df["locus"].str.upper().isin(["IGK","IGL"])].drop_duplicates("cell_id")
    heavy = heavy.rename(columns={c:f"heavy_{c}" for c in heavy.columns if c not in ["cell_id","locus"]})
    light = light.rename(columns={c:f"light_{c}" for c in light.columns if c not in ["cell_id","locus"]})
    return heavy.merge(light, on="cell_id", how="outer")


def get_dominant_heavy_c_call(cell_ids, df_hl):
    if "heavy_c_call" not in df_hl.columns:
        return ""
    sub = df_hl[df_hl["cell_id"].isin(cell_ids)]
    vc = sub["heavy_c_call"].value_counts()
    return vc.index[0] if len(vc) else ""


# ── Reporting utilities (CSV & plots) ─────────────────────────────────────────
def export_top5000_clonotypes(adata, vdj_obj, outdir, prefix=""):
    if "clone_id" not in adata.obs:
        vdj_obj.threshold = 0.85
        ddl.tl.find_clones(vdj_obj, key="junction_aa")
        ddl.tl.transfer(adata, vdj_obj, overwrite=True)

    clone_counts = adata.obs["clone_id"].value_counts(dropna=True)
    top_clone_ids = clone_counts.index[:5000]

    mask_cells = adata.obs["clone_id"].isin(top_clone_ids)
    top_cells = adata.obs.loc[mask_cells].copy()
    df_hl = build_heavy_light_table(vdj_obj)

    cinfo = top_cells[["clone_id"]].copy()
    for col in ("cell_type", "sample_id"):
        if col in top_cells.columns:
            cinfo[col] = top_cells[col]
    cinfo["cell_id"] = cinfo.index

    df_merged = pd.merge(cinfo.reset_index(drop=True), df_hl, on="cell_id", how="left")

    lead = ["cell_id","clone_id","cell_type","sample_id"]
    df_merged = df_merged[[c for c in lead if c in df_merged] +
                          [c for c in df_merged.columns if c not in lead]]

    out = os.path.join(outdir, f"top5000_full_clonotypes_{prefix}.csv")
    df_merged.to_csv(out, index=False)
    print("[INFO] →", out)


def plot_venn_cloneid_plasmaG1_plasmaG2_memoryBG1(adata, outdir, pre_ids, post_ids, prefix=""):
    try:
        from matplotlib_venn import venn3  # local import to make dependency optional
    except Exception:
        print("[WARN] matplotlib-venn not installed → Venn plot skipped.")
        return

    if "clone_id" not in adata.obs:
        return

    set_p1 = set(adata.obs.loc[(adata.obs["sample_id"].isin(pre_ids)) &
                               (adata.obs["cell_type"]=="Plasma cells"),
                               "clone_id"].dropna())
    set_p2 = set(adata.obs.loc[(adata.obs["sample_id"].isin(post_ids)) &
                               (adata.obs["cell_type"]=="Plasma cells"),
                               "clone_id"].dropna())
    set_mb = set(adata.obs.loc[(adata.obs["sample_id"].isin(pre_ids)) &
                               (adata.obs["cell_type"]=="Memory B cells"),
                               "clone_id"].dropna())

    plt.figure(figsize=(5,5))
    venn3([set_p1, set_p2, set_mb], set_labels=("Plasma G1", "Plasma G2", "Memory B G1"))
    fn = os.path.join(outdir, f"venn_cloneid_{prefix}.png")
    plt.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close()
    print("[INFO] →", fn)


def plot_upset_cloneid_clusters(adata, clusters, outdir, pre_ids, post_ids, prefix=""):
    if from_contents is None or UpSet is None:
        print("[WARN] upsetplot not installed → UpSet skipped.")
        return
    if "clone_id" not in adata.obs:
        return
    sets = {}
    for cl in clusters:
        mask = adata.obs["leiden"] == str(cl)
        if mask.sum() == 0:
            continue
        ad = adata[mask]
        sets[f"Cl{cl}_NTX123"] = set(ad.obs.loc[ad.obs["sample_id"].isin(pre_ids), "clone_id"].dropna())
        sets[f"Cl{cl}_NTX45"]  = set(ad.obs.loc[ad.obs["sample_id"].isin(post_ids), "clone_id"].dropna())

    if not sets:
        return

    from upsetplot import from_contents, UpSet  # ensure available if we reached here
    up = UpSet(from_contents(sets), subset_size="count")

    plt.figure(figsize=(10, 6))
    up.plot()
    fn = os.path.join(outdir, f"upset_cloneid_clusters_{prefix}.png")
    plt.savefig(fn, dpi=150, bbox_inches="tight")
    plt.close()
    print("[INFO] →", fn)


# ── Main pipeline ─────────────────────────────────────────────────────────────
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    pre_ids  = set([s.strip() for s in args.pre_samples.split(",") if s.strip()])
    post_ids = set([s.strip() for s in args.post_samples.split(",") if s.strip()])

    adata = sc.read_h5ad(args.adata)

    contigs = pd.read_csv(args.contigs, sep="\t")
    contigs["cell_id"] = contigs["cell_id"].map(rename_barcode)

    vdj_all = ddl.Dandelion(contigs)
    vdj_all, adata = ddl.pp.check_contigs(vdj_all, adata, productive_only=True)

    # Keep Memory B + Plasma only
    keep = adata.obs["cell_type"].isin(["Memory B cells", "Plasma cells"])
    ad_bpl = adata[keep].copy()
    vdj_bpl = vdj_all[vdj_all.metadata.index.isin(ad_bpl.obs_names)].copy()
    vdj_bpl.simplify()

    # Clones & transfer
    vdj_bpl.threshold = 0.85
    ddl.tl.find_clones(vdj_bpl, key="junction_aa")
    ddl.tl.generate_network(vdj_bpl, min_size=4, layout_method="mod_fr")
    ddl.tl.transfer(ad_bpl, vdj_bpl, overwrite=True)

    df_hl = build_heavy_light_table(vdj_bpl)
    export_top5000_clonotypes(ad_bpl, vdj_bpl, args.out_dir, prefix="bpl_only")
    plot_venn_cloneid_plasmaG1_plasmaG2_memoryBG1(ad_bpl, args.out_dir, pre_ids, post_ids, "bpl_only")
    plot_upset_cloneid_clusters(ad_bpl, [1,3,5,14,16], args.out_dir, pre_ids, post_ids, "bpl_only")

    print("[DONE] B/Plasma clonotyping pipeline →", args.out_dir)


if __name__ == "__main__":
    main()
