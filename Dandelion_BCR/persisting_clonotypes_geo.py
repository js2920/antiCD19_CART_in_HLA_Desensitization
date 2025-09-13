#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
B/Plasma clonotype visualisation pipeline (template)

Outputs → --out-dir (default: results/persisting_clonotypes)

Generates:
• UpSet overlaps for specified pairs / clusters
• totalVI UMAP overlays for key overlap sets
• Sequence exports for Cl1Pre∩Cl1Post and Cl1Pre∩Cl3Pre

"""
import os
import argparse
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import dandelion as ddl
import matplotlib.pyplot as plt
import seaborn as sns

# Optional deps (graceful fallbacks)
try:
    from upsetplot import from_contents, UpSet
except Exception:  # pragma: no cover
    from_contents, UpSet = None, None

try:
    from Bio.Seq import Seq  # optional; used for AA translation of concatenated fragments
except Exception:  # pragma: no cover
    Seq = None

warnings.filterwarnings("ignore", category=FutureWarning)

# --------------------------------------------------------------------- #
# CLI                                                                   #
# --------------------------------------------------------------------- #
def parse_args():
    ap = argparse.ArgumentParser(description="Persisting B/Plasma clonotypes visualisation")
    ap.add_argument("--adata", required=True,
                    help="H5AD with BCR integrated.")
    ap.add_argument("--contigs", required=True,
                    help="merged_contigs.tsv (same used for integration).")
    ap.add_argument("--out-dir", default="results/persisting_clonotypes",
                    help="Output directory.")
    ap.add_argument("--pre-samples", default="NTX_1,NTX_2,NTX_3")
    ap.add_argument("--post-samples", default="NTX_4,NTX_5")
    return ap.parse_args()


# --------------------------------------------------------------------- #
# settings & palettes                                                   #
# --------------------------------------------------------------------- #
PAL = {
    "Cl1PreCl1Post": "#e41a1c",
    "Cl1Cl3Cl5_Pre": "#377eb8",
    "Cl1PreCl3Pre":  "#4daf4a",
    "background":    "lightgrey",
}
UMAP_STYLE = {"bg_marker": "o", "bg_size": 10, "bg_alpha": 0.5,
              "hl_marker": "s", "hl_size": 40, "hl_edge": 0.3}

PAIR_SPECS = [("5", "pre", "3", "pre"),
              ("5", "pre", "3", "post"),
              ("5", "pre", "1", "pre"),
              ("5", "pre", "1", "post")]


# --------------------------------------------------------------------- #
# helpers – barcode, heavy/light, isotype                               #
# --------------------------------------------------------------------- #
def rename_bc(bc: str) -> str:
    return bc.replace("_fastq_", "_") + "-1" if "_fastq_" in bc else bc

def _clean(s: str) -> str:
    return str(s).replace(".", "").replace("-", "").upper()

def build_heavy_light(vdj: ddl.Dandelion) -> pd.DataFrame:
    cols = ["junction_aa","cell_id","locus","c_call","c_call_10x",
            "fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4"]
    cols = [c for c in cols if c in vdj.data.columns]
    rec = []
    for r in vdj.data[cols].itertuples(index=False):
        d = r._asdict()
        for k in ["fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4"]:
            if k in d: d[k] = _clean(d[k])
        rec.append(d)
    df = pd.DataFrame(rec)
    heavy = df[df["locus"].str.upper() == "IGH"].drop_duplicates("cell_id")
    light = df[df["locus"].str.upper().isin(["IGK","IGL"])].drop_duplicates("cell_id")
    heavy = heavy.rename(columns={c:f"heavy_{c}" for c in heavy.columns if c not in ["cell_id","locus"]})
    light = light.rename(columns={c:f"light_{c}" for c in light.columns if c not in ["cell_id","locus"]})
    return heavy.merge(light,on="cell_id",how="outer")

def dominant_isotype(df_hl: pd.DataFrame, cid: str, ad: sc.AnnData) -> str:
    cells = ad.obs.index[ad.obs["clone_id"] == cid]
    sub   = df_hl[df_hl["cell_id"].isin(cells)]
    for col in ["heavy_c_call","heavy_c_call_10x"]:
        if col in sub:
            vc = sub[col].dropna().value_counts()
            if len(vc): return vc.idxmax()
    return "NA"

def annotate_heavy_isotype(ad: sc.AnnData, df_hl: pd.DataFrame):
    if "heavy_iso" in ad.obs:
        return
    iso_col = "heavy_c_call" if "heavy_c_call" in df_hl.columns else "heavy_c_call_10x"
    ad.obs["heavy_iso"] = ad.obs.index.map(df_hl.set_index("cell_id")[iso_col]).fillna("")


# --------------------------------------------------------------------- #
# AA/fragments table & helpers                                          #
# --------------------------------------------------------------------- #
def build_aa_table_for_all_contigs(vdj_obj):
    needed = ["junction_aa","cell_id","locus","fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4"]
    cols = [c for c in needed if c in vdj_obj.data.columns]
    rec = []
    for row in vdj_obj.data[cols].itertuples(index=False):
        d = row._asdict()
        for p in ["fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4"]:
            if p in d: d[p] = _clean(str(d[p]))
        full = "".join(d.get(k,"") for k in ["fwr1","cdr1","fwr2","cdr2","fwr3","cdr3","fwr4"])
        d["full_nuc_seq"] = full
        d["full_aa_seq"]  = str(Seq(full).translate(to_stop=False)) if full and Seq else ""
        rec.append(d)
    return pd.DataFrame(rec)

def add_CDR3_and_celltype(df_in: pd.DataFrame, cell_to_type: pd.Series) -> pd.DataFrame:
    df = df_in.copy()
    if "junction_aa" in df.columns and "CDR3" not in df.columns:
        df.insert(df.columns.get_loc("junction_aa"), "CDR3", df["junction_aa"])
    if "cell_id" in df.columns:
        mapper = cell_to_type.astype("object").to_dict()
        df["celltypist_label"] = df["cell_id"].map(mapper).fillna("")
    return df


# --------------------------------------------------------------------- #
# UpSet helpers                                                         #
# --------------------------------------------------------------------- #
def _label_bars(up):
    ax = None
    if hasattr(up, "intersections"):
        ax = up.intersections
    if not ax:
        return
    for p in ax.patches:
        if p.get_height() > 0:
            ax.text(p.get_x()+p.get_width()/2, p.get_height()+0.3,
                    f"{int(p.get_height())}", ha="center", va="bottom", fontsize=8)

def plot_upset(sets, fn_out):
    if from_contents is None or UpSet is None:
        print("[WARN] upsetplot not installed → UpSet skipped.")
        return
    sets = {k: v for k, v in sets.items() if v}
    if len(sets) < 2:
        return
    plt.figure(figsize=(10, 6))
    up = UpSet(from_contents(sets), subset_size="count", show_counts=True)
    up.plot()
    _label_bars(up)
    plt.savefig(fn_out, dpi=150, bbox_inches="tight")
    plt.close()
    print("[INFO] UpSet ⇒", fn_out)


# --------------------------------------------------------------------- #
# Overlaps & UMAPs                                                      #
# --------------------------------------------------------------------- #
def _cells(ad, cl, stage, pre_set, post_set):
    mask = (ad.obs["leiden"] == str(cl)) & \
           (ad.obs["sample_id"].isin(pre_set if stage == "pre" else post_set))
    return ad.obs[mask]

def export_overlap_isotype(ad, vdj, df_hl, outdir, pre_set, post_set):
    print("[INFO] Overlap & isotype tables …")
    sets = {}
    for A, sa, B, sb in PAIR_SPECS:
        tagA, idxA = f"Cl{A}_{sa.capitalize()}", _cells(ad, A, sa, pre_set, post_set).index
        tagB, idxB = f"Cl{B}_{sb.capitalize()}", _cells(ad, B, sb, pre_set, post_set).index
        setA = set(ad.obs.loc[idxA, "clone_id"].dropna()); sets[tagA] = setA
        setB = set(ad.obs.loc[idxB, "clone_id"].dropna()); sets[tagB] = setB
        ov   = setA & setB
        if not ov:
            continue
        rows = []
        for cid in ov:
            rows.append({
                "clone_id": cid,
                "dominant_isotype": dominant_isotype(df_hl, cid, ad),
                f"{tagA}_cells": (ad.obs.loc[idxA, "clone_id"] == cid).sum(),
                f"{tagB}_cells": (ad.obs.loc[idxB, "clone_id"] == cid).sum()
            })
        fn = os.path.join(outdir, f"overlap_clonotypes_{tagA}_vs_{tagB}.csv")
        pd.DataFrame(rows).to_csv(fn, index=False)
        print("   ↳", fn, f"({len(rows)} clonotypes)")
    return sets

def _totalvi_umap_key(obsm):
    for k in obsm:
        if "totalvi" in k.lower() and "umap" in k.lower():
            return k
    return "X_umap" if "X_umap" in obsm else None

def plot_overlap_umaps(ad, outdir, pre_set, post_set):
    basis = _totalvi_umap_key(ad.obsm)
    if basis is None or "clone_id" not in ad.obs:
        return

    cl1_pre  = set(ad.obs.loc[(ad.obs["leiden"]=="1") & ad.obs["sample_id"].isin(pre_set),  "clone_id"].dropna())
    cl1_post = set(ad.obs.loc[(ad.obs["leiden"]=="1") & ad.obs["sample_id"].isin(post_set), "clone_id"].dropna())
    cl3_pre  = set(ad.obs.loc[(ad.obs["leiden"]=="3") & ad.obs["sample_id"].isin(pre_set),  "clone_id"].dropna())
    cl5_pre  = set(ad.obs.loc[(ad.obs["leiden"]=="5") & ad.obs["sample_id"].isin(pre_set),  "clone_id"].dropna())

    sets = {"Cl1PreCl1Post": cl1_pre & cl1_post,
            "Cl1Cl3Cl5_Pre": cl1_pre & cl3_pre & cl5_pre,
            "Cl1PreCl3Pre":  cl1_pre & cl3_pre}

    sample_ids = sorted(ad.obs["sample_id"].unique())
    for tag, cset in sets.items():
        if not cset:
            continue
        n = len(cset)
        for sid in sample_ids:
            sub = ad[ad.obs["sample_id"] == sid].copy()
            if sub.n_obs == 0:
                continue
            highlight = sub.obs["clone_id"].isin(cset).values
            if not highlight.any():
                continue
            x, y = sub.obsm[basis][:, 0], sub.obsm[basis][:, 1]

            fig, ax = plt.subplots(figsize=(5, 5))
            ax.scatter(x[~highlight], y[~highlight],
                       c="lightgrey", s=UMAP_STYLE["bg_size"],
                       marker=UMAP_STYLE["bg_marker"], alpha=UMAP_STYLE["bg_alpha"],
                       linewidths=0)
            ax.scatter(x[highlight], y[highlight],
                       c=PAL[tag], s=UMAP_STYLE["hl_size"],
                       marker=UMAP_STYLE["hl_marker"], edgecolors="k",
                       linewidths=UMAP_STYLE["hl_edge"],
                       label=f"{tag} ({n} clonotypes)")
            ax.set_title(f"{sid}: {tag}", pad=12)
            ax.set_xlabel("UMAP‑1"); ax.set_ylabel("UMAP‑2")
            ax.set_xticks([]); ax.set_yticks([])
            ax.legend(frameon=False, loc="lower right")
            plt.tight_layout()
            fn = os.path.join(outdir, f"{sid}_{tag}_custom_umap.png")
            plt.savefig(fn, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print("[INFO] UMAP ⇒", fn)
            del sub


# --------------------------------------------------------------------- #
# Sequence exports                                                      #
# --------------------------------------------------------------------- #
def export_cl1_overlap_sequences(ad, vdj, outdir, pre_set, post_set):
    cl1_pre  = set(ad.obs.loc[(ad.obs["leiden"]=="1") & ad.obs["sample_id"].isin(pre_set),  "clone_id"].dropna())
    cl1_post = set(ad.obs.loc[(ad.obs["leiden"]=="1") & ad.obs["sample_id"].isin(post_set), "clone_id"].dropna())
    overlap  = cl1_pre & cl1_post
    if not overlap:
        print("[INFO] No Cl1 PRE/POST overlapping clonotypes."); return

    mask   = ad.obs["clone_id"].isin(overlap)
    cells  = ad.obs[mask].copy()
    cells["cluster"] = "Cl" + cells["leiden"].astype(str)
    cells["stage"]   = np.where(cells["sample_id"].isin(pre_set), "Pre", "Post")

    aa_tbl = build_aa_table_for_all_contigs(vdj)
    aa_tbl = add_CDR3_and_celltype(aa_tbl, ad.obs["cell_type"])
    merged = aa_tbl.merge(cells[["cluster", "stage"]],
                          left_on="cell_id", right_index=True, how="inner")

    fn = os.path.join(outdir, "Cl1PreCl1Post_overlap_sequences.csv")
    merged.to_csv(fn, index=False)
    print("[INFO] Sequence CSV ⇒", fn, f"({merged.shape[0]} contigs)")


def export_cl1_cl3_pre_overlap_sequences(ad, vdj, outdir, pre_set):
    cl1_pre = set(ad.obs.loc[(ad.obs["leiden"]=="1") & ad.obs["sample_id"].isin(pre_set), "clone_id"].dropna())
    cl3_pre = set(ad.obs.loc[(ad.obs["leiden"]=="3") & ad.obs["sample_id"].isin(pre_set), "clone_id"].dropna())
    overlap = cl1_pre & cl3_pre
    if not overlap:
        print("[INFO] No Cl1 PRE / Cl3 PRE overlapping clonotypes."); return

    mask  = ad.obs["leiden"].isin(["1","3"]) & ad.obs["sample_id"].isin(pre_set) & ad.obs["clone_id"].isin(overlap)
    cells = ad.obs[mask].copy()
    cells["cluster"] = cells["leiden"].map({"1":"Cl1","3":"Cl3"})
    cells["stage"]   = "Pre"

    aa_tbl = build_aa_table_for_all_contigs(vdj)
    aa_tbl = add_CDR3_and_celltype(aa_tbl, ad.obs["cell_type"])
    merged = aa_tbl.merge(cells[["cluster","stage"]],
                          left_on="cell_id", right_index=True, how="inner")
    fn = os.path.join(outdir, "Cl1PreCl3Pre_overlap_sequences.csv")
    merged.to_csv(fn, index=False)
    print("[INFO] Sequence CSV ⇒", fn, f"({merged.shape[0]} contigs)")


# --------------------------------------------------------------------- #
# Cluster‑3‑pre overlap dot‑plot of top HVGs                            #
# --------------------------------------------------------------------- #
def _choose_hvg_flavor():
    try:
        from skmisc import loess  # noqa: F401
        return "seurat_v3"
    except Exception:
        print("[WARN] scikit-misc not installed → using 'seurat' HVG method.")
        return "seurat"

def _is_logged(adata): 
    try:
        return float(adata.X.max()) < 30
    except Exception:
        return True

def analyze_cl3_pre_overlap_dotplot(ad, df_hl, outdir, pre_set, top_n=50):
    cl1_pre = set(ad.obs.loc[(ad.obs["leiden"]=="1") & ad.obs["sample_id"].isin(pre_set), "clone_id"].dropna())
    cl3_pre = set(ad.obs.loc[(ad.obs["leiden"]=="3") & ad.obs["sample_id"].isin(pre_set), "clone_id"].dropna())
    overlap = cl1_pre & cl3_pre
    if not overlap:
        print("[INFO] No Cl1‑PRE / Cl3‑PRE overlap → dot‑plot skipped."); return

    mask = ad.obs["clone_id"].isin(overlap) & ad.obs["sample_id"].isin(pre_set) & ad.obs["leiden"].isin(["1","3"])
    sub = ad[mask].copy()
    if sub.n_obs == 0:
        print("[WARN] No cells matched overlap mask."); return

    # annotate heavy isotype (if available)
    if "heavy_iso" not in sub.obs:
        # this requires df_hl to be from the same AnnData subset
        pass

    def _grp(r):
        if r["leiden"]=="3" and r["cell_type"]=="Memory B cells": return "csMBCs"
        if r["leiden"]=="3" and r["cell_type"]=="Plasma cells" and str(r.get("heavy_iso","")).startswith("IGHG"):
            return "IgG+ MBC‑derived PCs"
        if r["leiden"]=="1" and r["cell_type"]=="Plasma cells" and str(r.get("heavy_iso","")).startswith("IGHG"):
            return "IgG+ Mature PCs"
        return np.nan

    sub.obs["dotplot_group"] = sub.obs.apply(_grp, axis=1)
    sub = sub[sub.obs["dotplot_group"].notna()].copy()
    if sub.n_obs == 0:
        print("[WARN] No cells qualify for the dot‑plot groups."); return

    if not _is_logged(sub):
        sc.pp.normalize_total(sub, target_sum=1e4); sc.pp.log1p(sub)

    hvg_flavor = _choose_hvg_flavor()
    sc.pp.highly_variable_genes(sub, flavor=hvg_flavor, n_top_genes=max(2000, top_n))
    hvgs = (sub.var.query("highly_variable").sort_values("highly_variable_rank")
            .index[:top_n].tolist())

    dp = sc.pl.dotplot(sub, var_names=hvgs, groupby="dotplot_group",
                       dendrogram=False, standard_scale="var", show=False)
    fn = os.path.join(outdir, f"cl3pre_overlap_top{top_n}HVG_dotplot.png")
    if hasattr(dp, "savefig"):      dp.savefig(fn, dpi=150, bbox_inches="tight")
    elif hasattr(dp, "figure"):     dp.figure.savefig(fn, dpi=150, bbox_inches="tight")
    else:                           plt.gcf().savefig(fn, dpi=150, bbox_inches="tight")
    plt.close("all")
    print("[INFO] Dot‑plot ⇒", fn)


# --------------------------------------------------------------------- #
# main                                                                  #
# --------------------------------------------------------------------- #
def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True); sc.settings.figdir = args.out_dir

    pre_set  = set([s.strip() for s in args.pre_samples.split(",") if s.strip()])
    post_set = set([s.strip() for s in args.post_samples.split(",") if s.strip()])

    ad = sc.read_h5ad(args.adata)
    df = pd.read_csv(args.contigs, sep="\t")
    df["cell_id"] = df["cell_id"].map(rename_bc)

    vdj_all, ad = ddl.pp.check_contigs(ddl.Dandelion(df), ad, productive_only=True)
    ad = ad[ad.obs["cell_type"].isin(["Memory B cells","Plasma cells"])].copy()
    vdj = vdj_all[vdj_all.metadata.index.isin(ad.obs_names)].copy(); vdj.simplify()

    vdj.threshold = 0.85
    ddl.tl.find_clones(vdj, key="junction_aa")
    ddl.tl.generate_network(vdj, min_size=4, layout_method="mod_fr")
    ddl.tl.transfer(ad, vdj, overwrite=True)

    df_hl = build_heavy_light(vdj)

    sets_overlap = export_overlap_isotype(ad, vdj, df_hl, args.out_dir, pre_set, post_set)
    plot_upset(sets_overlap, os.path.join(args.out_dir, "upset_overlap_specified_pairs_bpl.png"))

    sets_cluster = {}
    for cl in [1,3,5,14]:
        mask = ad.obs["leiden"] == str(cl)
        if mask.sum() == 0: 
            continue
        pre  = set(ad.obs.loc[mask & ad.obs["sample_id"].isin(pre_set),  "clone_id"].dropna())
        post = set(ad.obs.loc[mask & ad.obs["sample_id"].isin(post_set), "clone_id"].dropna())
        if pre:  sets_cluster[f"Cl{cl}_Pre"]  = pre
        if post: sets_cluster[f"Cl{cl}_Post"] = post
    plot_upset(sets_cluster, os.path.join(args.out_dir, "upset_clusters_1_3_5_14_bpl.png"))

    plot_overlap_umaps(ad, args.out_dir, pre_set, post_set)
    export_cl1_overlap_sequences(ad, vdj, args.out_dir, pre_set, post_set)
    export_cl1_cl3_pre_overlap_sequences(ad, vdj, args.out_dir, pre_set)
    analyze_cl3_pre_overlap_dotplot(ad, df_hl, args.out_dir, pre_set, top_n=50)

    print("\n[INFO] All outputs written to", args.out_dir)


if __name__ == "__main__":
    main()
