#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SCENIC on LLP overlap sequences — combined run (LLP_all), 10 iterations.
Aggregation + pre/post meta‑enrichment + downstream dotplots.
"""
from __future__ import annotations

import argparse
import inspect
import logging
import math
import pickle
import re
import warnings
from collections import Counter, defaultdict
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster
from matplotlib.colors import Normalize
from matplotlib.legend import Legend
from pyscenic.aucell import aucell as scenic_aucell
from pyscenic.prune import df2regulons, prune2df
from pyscenic.utils import modules_from_adjacencies

try:
    from scipy.stats import mannwhitneyu, norm
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

log = logging.getLogger("SCENIC_LLP_all")

# ---------- helpers ----------
def _safe_toarray(X): return X.toarray() if hasattr(X, "toarray") else X
def _normalize_adj_columns(df: pd.DataFrame):
    lower = {c.lower(): c for c in df.columns}
    tf_col = lower.get("tf", lower.get("regulator"))
    target_col = lower.get("target", lower.get("gene"))
    imp_col = lower.get("importance", lower.get("weight", lower.get("score")))
    if tf_col is None or target_col is None:
        raise ValueError("Could not find TF/target columns in adjacency DataFrame.")
    return tf_col, target_col, imp_col
def _get_regulon_targets(regulon):
    if hasattr(regulon, "genes"):
        try: return set(regulon.genes)
        except Exception: pass
    if hasattr(regulon, "gene2weight"): return set(getattr(regulon, "gene2weight").keys())
    return set()
def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float); n = p.size
    order = np.argsort(p); ranks = np.arange(1, n + 1, dtype=float)
    q = np.empty(n, dtype=float); q[order] = (p[order] * n) / ranks
    for i in range(n - 2, -1, -1): q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.minimum(q, 1.0)

# ---------- SCENIC multi-iteration ----------
def run_scenic(expr: pd.DataFrame, *, ranking_db: Path, motif_tbl: Path, tf_list: Path,
               out_dir: Path, n_iter: int, n_workers: int, threads_per_worker: int, seeds: list[int]):
    out_dir.mkdir(parents=True, exist_ok=True)
    # TFs present
    tfs = pd.read_csv(tf_list, header=None)[0].astype(str).tolist()
    tfs = [t for t in tfs if t in expr.columns]
    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker,
                           processes=True, dashboard_address=None)
    cli = Client(cluster)
    db = RankingDatabase(ranking_db, name="hg38_mc10")
    try:
        for i in range(1, n_iter + 1):
            iter_out = out_dir / f"iter_{i:02d}"; iter_out.mkdir(exist_ok=True, parents=True)
            seed = seeds[i - 1] if i - 1 < len(seeds) else None
            kwargs = dict(gene_names=expr.columns, tf_names=tfs, client_or_address=cli, verbose=False)
            if "seed" in inspect.signature(grnboost2).parameters and seed is not None:
                kwargs["seed"] = seed
            log.info(f"[{iter_out.name}] GRNBoost2 on {expr.shape[0]}×{expr.shape[1]} (seed={seed})")
            adj = grnboost2(expr, **kwargs)
            adj.to_csv(iter_out / f"adjacencies_LLP_all_iter{i:02d}.csv", index=False)
            mods = [m for m in modules_from_adjacencies(adj, expr)
                    if len(getattr(m, "genes", getattr(m, "gene2weight", []))) >= 10]
            pickle.dump(mods, open(iter_out / f"modules_LLP_all_iter{i:02d}.pkl", "wb"))
            with ProgressBar():
                pruned = prune2df([db], mods, motif_tbl, num_workers=1, client_or_address="custom_multiprocessing")
            pruned.to_csv(iter_out / f"pruned_motifs_LLP_all_iter{i:02d}.csv", index=False)
            regs = df2regulons(pruned)
            pickle.dump(regs, open(iter_out / f"regulons_LLP_all_iter{i:02d}.pkl", "wb"))
            auc = scenic_aucell(expr, regs, num_workers=1)
            auc.to_csv(iter_out / f"auc_mtx_LLP_all_iter{i:02d}.csv")
            pd.DataFrame([(r.name, len(_get_regulon_targets(r))) for r in regs],
                         columns=["Regulon", "NumTargets"]).sort_values("NumTargets", ascending=False)\
              .to_csv(iter_out / f"top_regulons_by_targetcount_LLP_all_iter{i:02d}.csv", index=False)
            (auc.mean(axis=0).rename("MeanAUC").to_frame().reset_index()
                .rename(columns={"index": "Regulon"}).sort_values("MeanAUC", ascending=False)
                .to_csv(iter_out / f"top_regulons_by_meanAUC_LLP_all_iter{i:02d}.csv", index=False))
            log.info(f"[{iter_out.name}] done ✓")
    finally:
        cli.close(); cluster.close()

# ---------- aggregation ----------
def aggregate_LLP_all(base_out: Path):
    agg_dir = base_out / "_aggregate"; agg_dir.mkdir(parents=True, exist_ok=True)
    iter_dirs = [p for p in base_out.iterdir() if p.is_dir() and p.name.startswith("iter_")]
    if not iter_dirs: log.warning("No iterations to aggregate."); return
    auc_series_per_iter = {}; run_metrics_rows = []
    edge_counts = Counter(); edge_imp_values = defaultdict(list)
    for idir in iter_dirs:
        run_id = idir.name
        adj_files = list(idir.glob("adjacencies_LLP_all_iter*.csv"))
        if adj_files:
            adj = pd.read_csv(adj_files[0]); tfc, tgtc, iw = _normalize_adj_columns(adj)
            for row in adj[[tfc, tgtc]].itertuples(index=False):
                edge_counts[(str(row[0]), str(row[1]))] += 1
            if iw is not None:
                for row in adj[[tfc, tgtc, iw]].itertuples(index=False):
                    edge_imp_values[(str(row[0]), str(row[1]))].append(float(row[2]))
            n_edges = len(adj)
        else:
            n_edges = np.nan
        auc_sum_files = list(idir.glob("top_regulons_by_meanAUC_LLP_all_iter*.csv"))
        if auc_sum_files:
            auc_sum = pd.read_csv(auc_sum_files[0])
            if {"Regulon", "MeanAUC"}.issubset(auc_sum.columns):
                auc_series_per_iter[run_id] = pd.Series(auc_sum["MeanAUC"].values,
                                                        index=auc_sum["Regulon"].astype(str))
        n_cells = np.nan
        auc_mtx_files = list(idir.glob("auc_mtx_LLP_all_iter*.csv"))
        if auc_mtx_files:
            try:
                with open(auc_mtx_files[0], "r") as f: n_cells = sum(1 for _ in f) - 1
            except Exception: pass
        run_metrics_rows.append(dict(Iteration=run_id, n_edges=n_edges, n_cells=n_cells))
    if auc_series_per_iter:
        auc_wide = pd.DataFrame(auc_series_per_iter).sort_index()
        auc_wide.to_csv(agg_dir / "AUC_mean_per_regulon_wide_LLP_all.csv")
        stats_reg = pd.DataFrame({
            "n_runs_detected": auc_wide.notna().sum(axis=1),
            "mean_meanAUC": auc_wide.mean(axis=1, skipna=True),
            "std_meanAUC": auc_wide.std(axis=1, ddof=1, skipna=True),
            "min_meanAUC": auc_wide.min(axis=1, skipna=True),
            "max_meanAUC": auc_wide.max(axis=1, skipna=True),
        })
        stats_reg["cv_meanAUC"] = stats_reg["std_meanAUC"] / stats_reg["mean_meanAUC"]
        stats_reg.sort_values("mean_meanAUC", ascending=False)\
          .to_csv(agg_dir / "auc_summary_by_regulon_LLP_all.csv")
        auc_wide.copy().set_axis([_base.split("(")[0].strip() for _base in auc_wide.index], axis=0)\
          .groupby(level=0).max().agg(["count", "mean"]).to_csv(agg_dir / "auc_summary_by_TF_LLP_all.csv")
        auc_wide.corr(min_periods=3).to_csv(agg_dir / "auc_mean_correlation_across_runs_LLP_all.csv")
    # edge stability
    if edge_counts:
        rows = []
        for (tf, tgt), cnt in edge_counts.items():
            vals = edge_imp_values.get((tf, tgt), [])
            rows.append(dict(TF=tf, target=tgt, n_runs=cnt,
                             frac_runs=cnt / len(run_metrics_rows) if run_metrics_rows else np.nan,
                             mean_importance=(float(np.mean(vals)) if vals else np.nan),
                             std_importance=(float(np.std(vals, ddof=1)) if len(vals) >= 2 else np.nan)))
        pd.DataFrame(rows).sort_values(["n_runs", "mean_importance"], ascending=[False, False])\
          .to_csv(agg_dir / "edge_stability_LLP_all.csv", index=False)
    pd.DataFrame(run_metrics_rows).set_index("Iteration").to_csv(agg_dir / "run_metrics_LLP_all.csv")

# ---------- meta‑enrichment pre vs post ----------
def meta_enrichment_pre_vs_post(all_dir: Path, pre_ids: list[str], post_ids: list[str], min_cells=3):
    agg_dir = all_dir / "_aggregate"; per_run_dir = agg_dir / "per_run"
    per_run_dir.mkdir(parents=True, exist_ok=True)
    iter_paths = sorted(all_dir.glob("iter_*/auc_mtx_*.csv"))
    if not iter_paths:
        log.warning("No AUCell matrices; skipping meta-enrichment."); return
    per_reg = defaultdict(lambda: dict(
        z_list=[], w_list=[], p_list=[], delta_list=[], rbc_list=[],
        n1_list=[], n2_list=[], pos_signs=0, neg_signs=0, n_runs_with_reg=0))
    n_total = 0
    for auc_path in iter_paths:
        n_total += 1; run_id = auc_path.parent.name
        auc = pd.read_csv(auc_path, index_col=0)
        pre = [cid for cid in pre_ids if cid in auc.index]
        post = [cid for cid in post_ids if cid in auc.index]
        n1, n2 = len(pre), len(post)
        if n1 < min_cells or n2 < min_cells:
            log.warning(f"{run_id}: insufficient cells (pre={n1}, post={n2})."); continue
        iter_rows = []
        for regulon in auc.columns:
            pre_vals = auc.loc[pre, regulon].astype(float).values
            post_vals = auc.loc[post, regulon].astype(float).values
            pre_mean = float(np.mean(pre_vals)); post_mean = float(np.mean(post_vals))
            delta = pre_mean - post_mean; sign = 1.0 if delta > 0 else (-1.0 if delta < 0 else 0.0)
            if HAS_SCIPY:
                try:
                    U, p = mannwhitneyu(pre_vals, post_vals, alternative="two-sided")
                    rbc = (2.0 * U / (n1 * n2)) - 1.0
                    z = (norm.isf(max(min(p, 1-1e-300), 1e-300)/2.0)) * (1 if sign >= 0 else -1)
                except Exception:
                    U = p = rbc = z = float("nan")
            else:
                U = p = rbc = z = float("nan")
            w = math.sqrt(n1 * n2)
            acc = per_reg[regulon]
            acc["z_list"].append(z); acc["w_list"].append(w); acc["p_list"].append(p)
            acc["delta_list"].append(delta); acc["rbc_list"].append(rbc)
            acc["n1_list"].append(n1); acc["n2_list"].append(n2); acc["n_runs_with_reg"] += 1
            if sign > 0: acc["pos_signs"] += 1
            elif sign < 0: acc["neg_signs"] += 1
            iter_rows.append(dict(Regulon=regulon, pre_mean=pre_mean, post_mean=post_mean,
                                  delta=delta, U_stat=U, p_value=p, r_rank_biserial=rbc,
                                  n_pre=n1, n_post=n2))
        pd.DataFrame(iter_rows).sort_values(["p_value", "delta"], ascending=[True, False])\
          .to_csv(per_run_dir / f"{run_id}_enrichment.csv", index=False)
    meta_rows = []
    for regulon, acc in per_reg.items():
        k = acc["n_runs_with_reg"]; 
        if k == 0: continue
        if HAS_SCIPY and any(np.isfinite(acc["z_list"])):
            z_vals = np.array([z for z in acc["z_list"] if np.isfinite(z)], dtype=float)
            w_vals = np.array([w for z, w in zip(acc["z_list"], acc["w_list"]) if np.isfinite(z)], dtype=float)
            if z_vals.size > 0 and w_vals.size > 0:
                Z_meta = float(np.sum(w_vals * z_vals) / max(np.sqrt(np.sum(w_vals ** 2)), 1e-300))
                p_meta = float(2.0 * norm.sf(abs(Z_meta)))
                enriched_in = "LLP_pre" if Z_meta > 0 else "LLP_post"
            else:
                Z_meta = p_meta = float("nan"); enriched_in = "NA"
        else:
            Z_meta = p_meta = float("nan"); enriched_in = "NA"
        meta_rows.append(dict(Regulon=regulon, k_runs=k, Z_meta=Z_meta, p_meta=p_meta, enriched_in=enriched_in))
    meta_df = pd.DataFrame(meta_rows)
    if not meta_df.empty and HAS_SCIPY and meta_df["p_meta"].notna().any():
        meta_df["FDR_BH"] = _bh_fdr(meta_df["p_meta"].fillna(1.0).values)
    else:
        meta_df["FDR_BH"] = np.nan
    meta_df.sort_values(["FDR_BH", "Z_meta"], ascending=[True, False])\
      .to_csv(all_dir / "_aggregate" / "regulon_enrichment_pre_vs_post_meta.csv", index=False)

# ---------- dotplots ----------
def _dot_nice_levels_pct(vals, k=3):
    v = np.asarray(vals, dtype=float) * 100.0; v = v[np.isfinite(v)]
    if v.size == 0: return [10, 50, 90]
    lo, hi = np.nanpercentile(v, 10), np.nanpercentile(v, 90)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = float(np.min(v)), float(np.max(v) if np.max(v) > np.min(v) else np.min(v) + 1)
    raw = np.linspace(lo, hi, k); return [int(round(x)) for x in raw]
def _dot_map_size(val, vmin, vmax, smin=20, smax=300):
    if not np.isfinite(val): return smin
    if vmin == vmax: return (smin + smax) / 2.0
    frac = np.clip((val - vmin) / (vmax - vmin), 0.0, 1.0)
    return smin + frac * (smax - smin)
def make_dotplot_pre_post(df: pd.DataFrame, selected: list[str], out_png: Path, thr=0.10):
    sns.set_style("white"); plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 300, "font.size": 8,
                                                 "legend.frameon": False, "axes.grid": False})
    sub = df[df["Regulon"].isin(selected)].copy()
    plot_df = pd.DataFrame({
        "Regulon": np.repeat(sub["Regulon"].values, 2),
        "Phase":   np.tile(["Pre", "Post"], len(sub)),
        "MeanAUC": np.concatenate([sub["MeanAUC_pre"].values,  sub["MeanAUC_post"].values]),
        "FracHigh":np.concatenate([sub["FracHigh_pre"].values, sub["FracHigh_post"].values]),
    })
    plot_df["Regulon"] = pd.Categorical(plot_df["Regulon"], categories=list(reversed(selected)), ordered=True)
    vmin = float(np.nanmin(plot_df["MeanAUC"])) if np.isfinite(np.nanmin(plot_df["MeanAUC"])) else 0.0
    vmax = float(np.nanmax(plot_df["MeanAUC"])) if np.isfinite(np.nanmax(plot_df["MeanAUC"])) else 1.0
    if vmin == vmax: vmax = vmin + 1e-6
    norm = Normalize(vmin=vmin, vmax=vmax)
    h = max(3.8, 0.28 * len(selected) + 1.2); fig, ax = plt.subplots(figsize=(4.8, h))
    smin, smax = 20, 300
    sns.scatterplot(data=plot_df, x="Phase", y="Regulon", hue="MeanAUC", palette="viridis", hue_norm=norm,
                    size="FracHigh", sizes=(smin, smax), linewidth=0, edgecolor="none", ax=ax, legend=False)
    ax.set_ylabel(""); ax.set_xlabel("")
    ax.set_title(f"Regulon activity (color=mean AUCell, size=% cells ≥ {thr:g}) (≥8/10 runs)", fontsize=9)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap("viridis"), norm=norm); sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.04); cbar.set_label("Mean AUCell", fontsize=6)
    pct_levels = _dot_nice_levels_pct(plot_df["FracHigh"], k=3)
    vmin_s = float(np.nanmin(plot_df["FracHigh"])) if np.isfinite(np.nanmin(plot_df["FracHigh"])) else 0.0
    vmax_s = float(np.nanmax(plot_df["FracHigh"])) if np.isfinite(np.nanmax(plot_df["FracHigh"])) else 1.0
    handles, labels = [], []
    for p in pct_levels:
        marker_area = _dot_map_size(p / 100.0, vmin_s, vmax_s, smin=smin, smax=smax)
        hndl = ax.scatter([], [], s=marker_area, facecolor="none", edgecolor="k")
        handles.append(hndl); labels.append(f"{p}%")
    Legend(ax, handles, labels, title="% cells ≥ thr", bbox_to_anchor=(1.03, 1.0), loc="upper left",
           borderaxespad=0., fontsize=6, title_fontsize=6)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight"); fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="SCENIC on LLP overlap sequences (LLP_all, 10×) + meta enrichment + dotplots.")
    ap.add_argument("--anndata", required=True, type=Path)
    ap.add_argument("--cell-list-csv", required=True, type=Path)
    ap.add_argument("--out-base", required=True, type=Path)
    ap.add_argument("--ranking-db", required=True, type=Path)
    ap.add_argument("--motif-table", required=True, type=Path)
    ap.add_argument("--tf-list", required=True, type=Path)
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--threads-per-worker", type=int, default=1)
    ap.add_argument("--n-iter", type=int, default=10)
    ap.add_argument("--dot-threshold", type=float, default=0.10)
    ap.add_argument("--dot-min-runs", type=int, default=8)
    ap.add_argument("--top-n", nargs="+", type=int, default=[30, 15])
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()
    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(levelname)-8s %(message)s", datefmt="%H:%M:%S")

    OUT_BASE = args.out_base; ALL_DIR = OUT_BASE / "scenic_LLP_all"
    DOT_OUT_DIR = OUT_BASE / "LLP_pre_vs_post"; DOT_OUT_DIR.mkdir(parents=True, exist_ok=True)

    log.info(f"Reading AnnData → {args.anndata}")
    adata = ad.read_h5ad(args.anndata); log.info(f"AnnData: {adata.n_obs} cells, {adata.n_vars} genes")

    log.info(f"Reading cell list → {args.cell-list-csv}")
    df = pd.read_csv(args.cell_list_csv, usecols=["cell_id"])
    cell_ids = df["cell_id"].astype(str).unique()
    grp_pre  = [cid for cid in cell_ids if re.match(r"^NTX_[123]_", cid)]
    grp_post = [cid for cid in cell_ids if re.match(r"^NTX_[45]_", cid)]
    pre_in_adata  = [cid for cid in grp_pre  if cid in adata.obs_names]
    post_in_adata = [cid for cid in grp_post if cid in adata.obs_names]
    all_in_adata  = [cid for cid in pre_in_adata + post_in_adata if cid in adata.obs_names]
    if len(all_in_adata) < 2: raise SystemExit("LLP_all: <2 cells → aborting.")

    X = adata.raw.to_adata() if adata.raw else adata.copy()
    if "counts" in X.layers: X.X = X.layers["counts"].copy()
    X.var["highly_variable"] = True
    expr = pd.DataFrame(_safe_toarray(X.X), index=X.obs_names, columns=X.var_names).loc[all_in_adata]

    seeds = [1337 + i for i in range(args.n_iter)]
    run_scenic(expr, ranking_db=args.ranking_db, motif_tbl=args.motif_table, tf_list=args.tf_list,
               out_dir=ALL_DIR, n_iter=args.n_iter, n_workers=args.n_workers,
               threads_per_worker=args.threads_per_worker, seeds=seeds)

    aggregate_LLP_all(ALL_DIR)

    if len(pre_in_adata) >= 3 and len(post_in_adata) >= 3:
        if not HAS_SCIPY:
            log.warning("SciPy not found — meta p-values/FDR will be NaN; install scipy for full stats.")
        meta_enrichment_pre_vs_post(ALL_DIR, pre_in_adata, post_in_adata, min_cells=3)

    # Dotplots
    # Load AUCell means and frac>=thr per run, then average across runs
    iter_dirs = sorted([p for p in (ALL_DIR).glob("iter_*") if p.is_dir()],
                       key=lambda p: int(re.search(r'(\d+)', p.name).group(1)) if re.search(r'(\d+)', p.name) else 0)
    rows = []
    for idir in iter_dirs:
        auc_file = next(iter(idir.glob("auc_mtx_*.csv")), None)
        if auc_file is None: continue
        auc = pd.read_csv(auc_file, index_col=0)
        pre = [cid for cid in pre_in_adata if cid in auc.index]
        post = [cid for cid in post_in_adata if cid in auc.index]
        if len(pre) == 0 or len(post) == 0: continue
        rows.append(pd.DataFrame({
            "Regulon": auc.columns.astype(str),
            "MeanAUC_pre":   auc.loc[pre].mean(axis=0).values.astype(float),
            "FracHigh_pre":  (auc.loc[pre]  >= args.dot_threshold).mean(axis=0).values.astype(float),
            "MeanAUC_post":  auc.loc[post].mean(axis=0).values.astype(float),
            "FracHigh_post": (auc.loc[post] >= args.dot_threshold).mean(axis=0).values.astype(float),
            "Iteration": idir.name,
        }))
    if rows:
        df_all = pd.concat(rows, ignore_index=True)
        agg = (df_all.groupby("Regulon")
               .agg(MeanAUC_pre=("MeanAUC_pre", "mean"),
                    FracHigh_pre=("FracHigh_pre", "mean"),
                    MeanAUC_post=("MeanAUC_post", "mean"),
                    FracHigh_post=("FracHigh_post", "mean"),
                    n_runs_detected=("Iteration", "nunique"))
               .reset_index())
        agg = agg[agg["n_runs_detected"] >= args.dot_min_runs].copy()
        src = DOT_OUT_DIR / "dotplot_source_table_LLP_pre_vs_LLP_post.csv"
        agg.to_csv(src, index=False)
        # rank regulons: prefer meta file if exists
        meta_path = ALL_DIR / "_aggregate" / "regulon_enrichment_pre_vs_post_meta.csv"
        if meta_path.exists():
            meta = pd.read_csv(meta_path)
            if {"Regulon", "Z_meta", "FDR_BH"}.issubset(meta.columns):
                meta_sig = meta[meta["FDR_BH"] <= 0.10].copy()
                meta_sig["absZ"] = meta_sig["Z_meta"].abs()
                ranked = meta_sig.sort_values("absZ", ascending=False)["Regulon"].tolist()
            else:
                ranked = (agg.assign(DeltaAUC=(agg["MeanAUC_post"] - agg["MeanAUC_pre"]).abs())
                              .sort_values("DeltaAUC", ascending=False)["Regulon"].tolist())
        else:
            ranked = (agg.assign(DeltaAUC=(agg["MeanAUC_post"] - agg["MeanAUC_pre"]).abs())
                          .sort_values("DeltaAUC", ascending=False)["Regulon"].tolist())
        for n in args.top_n:
            sel = ranked[:n]
            make_dotplot_pre_post(agg, sel, DOT_OUT_DIR / f"dotplot_LLP_pre_vs_LLP_post_top{n}.png",
                                  thr=args.dot_threshold)
    log.info("LLP_all pipeline complete ✔")

if __name__ == "__main__":
    main()
