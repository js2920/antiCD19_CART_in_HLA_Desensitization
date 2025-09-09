#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SCENIC on refined CellTypist PRE-only subsets — union run (10 iterations)
Aggregation + pairwise meta‑enrichment + dotplots split back into 3 groups.

Groups:
  • MBC_IgG_PC_pre
  • Mature_IgG_PC_pre
  • csMBC_pre
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
from itertools import combinations
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

# Optional stats
try:
    from scipy.stats import mannwhitneyu, norm
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False

log = logging.getLogger("SCENIC_CTpre_all")
GROUPS = ["MBC_IgG_PC_pre", "Mature_IgG_PC_pre", "csMBC_pre"]

# ---------- helpers ----------
def _safe_toarray(X): return X.toarray() if hasattr(X, "toarray") else X
def _iter_dirs_for_label(base: Path, label: str):
    p = base / f"scenic_{label}"
    return sorted([d for d in p.iterdir() if d.is_dir() and d.name.startswith("iter_")]) if p.exists() else []
def _extract_tf_from_regulon(regulon) -> str:
    for attr in ("transcription_factor", "tf", "name"):
        if hasattr(regulon, attr):
            val = getattr(regulon, attr)
            if isinstance(val, str) and val:
                return re.split(r"\s*\(", val)[0].strip()
    return re.split(r"\s*\(", str(getattr(regulon, "name", regulon)))[0].strip()
def _get_regulon_targets(regulon):
    if hasattr(regulon, "genes"):
        try: return set(regulon.genes)
        except Exception: pass
    if hasattr(regulon, "gene2weight"): return set(getattr(regulon, "gene2weight").keys())
    return set()
def _normalize_adj_columns(df: pd.DataFrame):
    lower = {c.lower(): c for c in df.columns}
    tf_col = lower.get("tf", lower.get("regulator"))
    target_col = lower.get("target", lower.get("gene"))
    imp_col = lower.get("importance", lower.get("weight", lower.get("score")))
    if tf_col is None or target_col is None:
        raise ValueError("Could not find TF/target columns in adjacency DataFrame.")
    return tf_col, target_col, imp_col
def _base_regulon_name(name: str) -> str:
    return re.split(r"\s*\(", str(name))[0].strip()
def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float); n = p.size
    order = np.argsort(p); ranks = np.arange(1, n + 1, dtype=float)
    q = np.empty(n, dtype=float); q[order] = (p[order] * n) / ranks
    for i in range(n - 2, -1, -1): q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.minimum(q, 1.0)
def _signed_z_from_two_sided_p(p_two_sided: float, sign: float, eps=1e-300) -> float:
    if not HAS_SCIPY: return float("nan")
    p = float(max(min(p_two_sided, 1 - eps), eps)); z_abs = norm.isf(p / 2.0)
    return math.copysign(z_abs, sign if sign != 0 else 1.0)

# ---------- SCENIC multi-iteration ----------
def run_scenic(expr: pd.DataFrame, *, ranking_db: Path, motif_tbl: Path, tf_list: Path,
               out_dir: Path, n_iter: int, n_workers: int, threads_per_worker: int, seeds: list[int]):
    out_dir.mkdir(parents=True, exist_ok=True)
    # TFs present
    tfs = pd.read_csv(tf_list, header=None)[0].astype(str).tolist()
    tfs = [t for t in tfs if t in expr.columns]
    # Dask cluster reused across iterations
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
            adj.to_csv(iter_out / f"adjacencies_CTpre_all_iter{i:02d}.csv", index=False)
            mods = [m for m in modules_from_adjacencies(adj, expr)
                    if len(getattr(m, "genes", getattr(m, "gene2weight", []))) >= 10]
            pickle.dump(mods, open(iter_out / f"modules_CTpre_all_iter{i:02d}.pkl", "wb"))
            with ProgressBar():
                pruned = prune2df([db], mods, motif_tbl, num_workers=1, client_or_address="custom_multiprocessing")
            pruned.to_csv(iter_out / f"pruned_motifs_CTpre_all_iter{i:02d}.csv", index=False)
            regs = df2regulons(pruned)
            pickle.dump(regs, open(iter_out / f"regulons_CTpre_all_iter{i:02d}.pkl", "wb"))
            auc = scenic_aucell(expr, regs, num_workers=1)  # cells × regulons
            auc.to_csv(iter_out / f"auc_mtx_CTpre_all_iter{i:02d}.csv")
            pd.DataFrame([(r.name, len(_get_regulon_targets(r))) for r in regs],
                         columns=["Regulon", "NumTargets"]).sort_values("NumTargets", ascending=False)\
              .to_csv(iter_out / f"top_regulons_by_targetcount_CTpre_all_iter{i:02d}.csv", index=False)
            (auc.mean(axis=0).rename("MeanAUC").to_frame().reset_index()
                .rename(columns={"index": "Regulon"}).sort_values("MeanAUC", ascending=False)
                .to_csv(iter_out / f"top_regulons_by_meanAUC_CTpre_all_iter{i:02d}.csv", index=False))
            log.info(f"[{iter_out.name}] done ✓")
    finally:
        cli.close(); cluster.close()

# ---------- aggregation ----------
def aggregate_CTpre_all(base_out: Path):
    agg_dir = base_out / "_aggregate"; agg_dir.mkdir(parents=True, exist_ok=True)
    iter_dirs = [p for p in base_out.iterdir() if p.is_dir() and p.name.startswith("iter_")]
    if not iter_dirs:
        log.warning("No iterations found to aggregate."); return
    auc_series_per_iter = {}; tf_sets_per_iter = {}; regulon_target_sets = defaultdict(list)
    run_metrics_rows = []; edge_counts = Counter(); edge_imp_values = defaultdict(list)
    for idir in iter_dirs:
        run_id = idir.name
        # adj
        adj_files = list(idir.glob("adjacencies_CTpre_all_iter*.csv"))
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
        # regs
        reg_files = list(idir.glob("regulons_CTpre_all_iter*.pkl")); tf_set = set(); target_sizes = []
        if reg_files:
            try:
                regs = pickle.load(open(reg_files[0], "rb"))
                for r in regs:
                    tf = _extract_tf_from_regulon(r); tf_set.add(tf)
                    targets = _get_regulon_targets(r); target_sizes.append(len(targets))
                    if len(targets) > 0: regulon_target_sets[tf].append(set(targets))
                n_regulons = len(regs); n_unique_tfs = len(tf_set)
                median_targets = float(np.median(target_sizes)) if target_sizes else np.nan
            except Exception:
                n_regulons = n_unique_tfs = median_targets = np.nan
        else:
            n_regulons = n_unique_tfs = median_targets = np.nan
        tf_sets_per_iter[run_id] = tf_set
        # AUC summaries
        auc_sum_files = list(idir.glob("top_regulons_by_meanAUC_CTpre_all_iter*.csv"))
        if auc_sum_files:
            auc_sum = pd.read_csv(auc_sum_files[0])
            if {"Regulon", "MeanAUC"}.issubset(auc_sum.columns):
                auc_series_per_iter[run_id] = pd.Series(auc_sum["MeanAUC"].values,
                                                        index=auc_sum["Regulon"].astype(str))
        # AUC matrix row count
        n_cells = np.nan
        auc_mtx_files = list(idir.glob("auc_mtx_CTpre_all_iter*.csv"))
        if auc_mtx_files:
            try:
                with open(auc_mtx_files[0], "r") as f: n_cells = sum(1 for _ in f) - 1
            except Exception: pass
        run_metrics_rows.append(dict(Iteration=run_id, n_edges=n_edges, n_regulons=n_regulons,
                                     n_unique_TFs=n_unique_tfs, median_targets_per_regulon=median_targets,
                                     n_cells=n_cells))
    run_ids = list(tf_sets_per_iter.keys())
    # AUC wide
    if auc_series_per_iter:
        auc_wide = pd.DataFrame(auc_series_per_iter).sort_index()
        auc_wide.to_csv(agg_dir / "AUC_mean_per_regulon_wide_CTpre_all.csv")
        stats_reg = pd.DataFrame({
            "n_runs_detected": auc_wide.notna().sum(axis=1),
            "mean_meanAUC": auc_wide.mean(axis=1, skipna=True),
            "std_meanAUC": auc_wide.std(axis=1, ddof=1, skipna=True),
            "min_meanAUC": auc_wide.min(axis=1, skipna=True),
            "max_meanAUC": auc_wide.max(axis=1, skipna=True),
        })
        stats_reg["cv_meanAUC"] = stats_reg["std_meanAUC"] / stats_reg["mean_meanAUC"]
        stats_reg.sort_values("mean_meanAUC", ascending=False)\
          .to_csv(agg_dir / "auc_summary_by_regulon_CTpre_all.csv")
        auc_by_tf = auc_wide.copy(); auc_by_tf.index = [_base_regulon_name(i) for i in auc_by_tf.index]
        auc_by_tf = auc_by_tf.groupby(auc_by_tf.index).max()
        stats_tf = pd.DataFrame({
            "n_runs_detected": auc_by_tf.notna().sum(axis=1),
            "mean_meanAUC": auc_by_tf.mean(axis=1, skipna=True),
            "std_meanAUC": auc_by_tf.std(axis=1, ddof=1, skipna=True),
            "min_meanAUC": auc_by_tf.min(axis=1, skipna=True),
            "max_meanAUC": auc_by_tf.max(axis=1, skipna=True),
        })
        stats_tf["cv_meanAUC"] = stats_tf["std_meanAUC"] / stats_tf["mean_meanAUC"]
        stats_tf.sort_values("mean_meanAUC", ascending=False)\
          .to_csv(agg_dir / "auc_summary_by_TF_CTpre_all.csv")
        auc_wide.corr(min_periods=3).to_csv(agg_dir / "auc_mean_correlation_across_runs_CTpre_all.csv")
    # Jaccard TF presence
    if run_ids:
        j = pd.DataFrame(index=run_ids, columns=run_ids, dtype=float)
        for a in run_ids:
            for b in run_ids:
                A, B = tf_sets_per_iter[a], tf_sets_per_iter[b]; u = len(A | B)
                j.loc[a, b] = (len(A & B) / u) if u else np.nan
        j.to_csv(agg_dir / "jaccard_regulon_TF_sets_CTpre_all.csv")
    # edge stability
    if edge_counts:
        rows = []
        for (tf, tgt), cnt in edge_counts.items():
            vals = edge_imp_values.get((tf, tgt), [])
            rows.append(dict(TF=tf, target=tgt, n_runs=cnt,
                             frac_runs=cnt / len(run_ids) if run_ids else np.nan,
                             mean_importance=(float(np.mean(vals)) if vals else np.nan),
                             std_importance=(float(np.std(vals, ddof=1)) if len(vals) >= 2 else np.nan)))
        pd.DataFrame(rows).sort_values(["n_runs", "mean_importance"], ascending=[False, False])\
          .to_csv(agg_dir / "edge_stability_CTpre_all.csv", index=False)
    # run metrics
    pd.DataFrame(run_metrics_rows).set_index("Iteration").to_csv(agg_dir / "run_metrics_CTpre_all.csv")

# ---------- meta‑enrichment across the 3 PRE groups ----------
def meta_enrichment_CTpre_all(all_dir: Path, group2ids: dict[str, list[str]], *,
                              min_cells: int = 3):
    base = all_dir; agg_dir = base / "_aggregate"; agg_dir.mkdir(exist_ok=True, parents=True)
    pairs = [(a, b) for a, b in combinations(GROUPS, 2) if a in group2ids and b in group2ids]
    if not pairs: 
        log.warning("No valid group pairs for enrichment."); return
    per_pair_reg = {(a, b): defaultdict(lambda: dict(
        z_list=[], w_list=[], p_list=[], delta_list=[], rbc_list=[], n1_list=[], n2_list=[], n_runs_with_reg=0
    )) for (a, b) in pairs}
    iter_paths = sorted(base.glob("iter_*/auc_mtx_*.csv"))
    if not iter_paths:
        log.warning("No AUCell matrices found; skipping meta-enrichment."); return
    for auc_path in iter_paths:
        run_id = auc_path.parent.name
        auc = pd.read_csv(auc_path, index_col=0)
        for a, b in pairs:
            a_ids = [cid for cid in group2ids[a] if cid in auc.index]
            b_ids = [cid for cid in group2ids[b] if cid in auc.index]
            n1, n2 = len(a_ids), len(b_ids)
            if n1 < min_cells or n2 < min_cells:
                log.warning(f"{run_id} [{a} vs {b}]: insufficient cells ({n1}/{n2})."); continue
            a_mean = auc.loc[a_ids].mean(axis=0); b_mean = auc.loc[b_ids].mean(axis=0)
            iter_rows = []
            for regulon in auc.columns:
                a_vals = auc.loc[a_ids, regulon].astype(float).values
                b_vals = auc.loc[b_ids, regulon].astype(float).values
                delta = float(a_mean[regulon] - b_mean[regulon])
                sign = 1.0 if delta > 0 else (-1.0 if delta < 0 else 0.0)
                if HAS_SCIPY:
                    try:
                        U, p = mannwhitneyu(a_vals, b_vals, alternative="two-sided")
                        rbc = (2.0 * U / (n1 * n2)) - 1.0
                        z = _signed_z_from_two_sided_p(p, sign)
                    except Exception:
                        U = p = rbc = z = float("nan")
                else:
                    U = p = rbc = z = float("nan")
                w = math.sqrt(n1 * n2)
                acc = per_pair_reg[(a, b)][regulon]
                acc["z_list"].append(z); acc["w_list"].append(w); acc["p_list"].append(p)
                acc["delta_list"].append(delta); acc["rbc_list"].append(rbc)
                acc["n1_list"].append(n1); acc["n2_list"].append(n2); acc["n_runs_with_reg"] += 1
                iter_rows.append(dict(Regulon=regulon, mean_A=a_mean[regulon], mean_B=b_mean[regulon],
                                      delta=delta, U_stat=U, p_value=p, r_rank_biserial=rbc,
                                      n_A=n1, n_B=n2))
            pr_dir = agg_dir / "per_run" / f"{re.sub('[^A-Za-z0-9]+','_',a)}_vs_{re.sub('[^A-Za-z0-9]+','_',b)}"
            pr_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(iter_rows).sort_values(["p_value", "delta"], ascending=[True, False])\
              .to_csv(pr_dir / f"{run_id}_enrichment.csv", index=False)
    # combine
    wide_rows = []
    for a, b in pairs:
        rows = []
        for regulon, acc in per_pair_reg[(a, b)].items():
            k = acc["n_runs_with_reg"]; 
            if k == 0: continue
            if HAS_SCIPY and any(np.isfinite(acc["z_list"])):
                z_vals = np.array([z for z in acc["z_list"] if np.isfinite(z)], dtype=float)
                w_vals = np.array([w for z, w in zip(acc["z_list"], acc["w_list"]) if np.isfinite(z)], dtype=float)
                if z_vals.size > 0 and w_vals.size > 0:
                    Z_meta = float(np.sum(w_vals * z_vals) / max(np.sqrt(np.sum(w_vals ** 2)), 1e-300))
                    p_meta = float(2.0 * norm.sf(abs(Z_meta)))
                    enriched_in = a if Z_meta > 0 else b
                else:
                    Z_meta = p_meta = float("nan"); enriched_in = "NA"
            else:
                Z_meta = p_meta = float("nan"); enriched_in = "NA"
            rows.append(dict(Regulon=regulon, k_runs=k, Z_meta=Z_meta, p_meta=p_meta,
                             enriched_in=enriched_in))
        meta_df = pd.DataFrame(rows)
        if not meta_df.empty and HAS_SCIPY and meta_df["p_meta"].notna().any():
            meta_df["FDR_BH"] = _bh_fdr(meta_df["p_meta"].fillna(1.0).values)
        else:
            meta_df["FDR_BH"] = np.nan
        meta_df = meta_df.sort_values(["FDR_BH", "Z_meta"], ascending=[True, False])
        out_csv = (all_dir / "_aggregate" / f"meta_enrichment_{re.sub('[^A-Za-z0-9]+','_',a)}_vs_{re.sub('[^A-Za-z0-9]+','_',b)}.csv")
        meta_df.to_csv(out_csv, index=False)
        for _, r in meta_df.iterrows():
            wide_rows.append(dict(Regulon=r["Regulon"], pair=f"{a}_vs_{b}",
                                  Z_meta=r["Z_meta"], p_meta=r["p_meta"],
                                  FDR_BH=r["FDR_BH"], enriched_in=r["enriched_in"]))
    if wide_rows:
        pd.DataFrame(wide_rows).to_csv(all_dir / "_aggregate" / "meta_enrichment_all_pairs_wide.csv", index=False)

# ---------- dotplots (from CTpre_all AUCell) ----------
def _dot_iter_dirs(base: Path):
    iters = [p for p in base.glob("iter_*") if p.is_dir()]
    def _key(p): m = re.search(r'(\d+)', p.name); return int(m.group(1)) if m else 0
    return sorted(iters, key=_key)

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

def _dot_save_fig(fig: plt.Figure, path_png: Path):
    fig.savefig(path_png, bbox_inches="tight"); fig.savefig(path_png.with_suffix(".svg"), bbox_inches="tight")
    log.info(f"Saved → {path_png} & {path_png.with_suffix('.svg')}")

def _dot_write_csv(df: pd.DataFrame, path: Path):
    df.to_csv(path, index=False); log.info(f"Wrote → {path}")

def _dot_load_multigroup_from_all(all_dir: Path, group2ids: dict, thr=0.10) -> pd.DataFrame:
    iter_dirs = _dot_iter_dirs(all_dir)
    if not iter_dirs: raise FileNotFoundError(f"No iteration directories found in {all_dir}")
    rows = []; used = 0
    for idir in iter_dirs:
        auc_file = next(iter(idir.glob("auc_mtx_*.csv")), None)
        if auc_file is None: continue
        auc = pd.read_csv(auc_file, index_col=0)
        per_iter = None; ok_iter = True
        for grp in GROUPS:
            ids = [cid for cid in group2ids.get(grp, []) if cid in auc.index]
            if len(ids) == 0:
                log.warning(f"{idir.name}: no cells for {grp}; skipping this iteration."); ok_iter = False; break
            mean_s  = auc.loc[ids].mean(axis=0).rename(f"MeanAUC_{grp}")
            frac_s  = (auc.loc[ids] >= thr).mean(axis=0).rename(f"FracHigh_{grp}")
            df_g = pd.concat([mean_s, frac_s], axis=1).reset_index().rename(columns={"index": "Regulon"})
            per_iter = df_g if per_iter is None else per_iter.merge(df_g, on="Regulon", how="outer")
        if not ok_iter or per_iter is None: continue
        per_iter["Iteration"] = idir.name; rows.append(per_iter); used += 1
    if not rows: raise FileNotFoundError(f"No usable AUCell CSVs found in {all_dir} for all groups.")
    df_all = pd.concat(rows, ignore_index=True)
    agg_ops = {f"MeanAUC_{g}": "mean" for g in GROUPS}
    for g in GROUPS: agg_ops[f"FracHigh_{g}"] = "mean"
    agg = (df_all.groupby("Regulon").agg({**agg_ops, "Iteration": "nunique"})
           .rename(columns={"Iteration": "n_runs_detected"}).reset_index())
    log.info(f"[DOTPLOT] Averaged over {used} iterations → {agg.shape[0]} regulons")
    return agg

def _rank_regulons_multigroup(df_merged: pd.DataFrame, agg_dir: Path, top_n: int, fdr=0.10) -> list[str]:
    meta_files = list((agg_dir).glob("meta_enrichment_*_vs_*.csv")); ranked = []
    if meta_files:
        frames = []
        for mf in meta_files:
            try:
                tmp = pd.read_csv(mf, usecols=["Regulon", "Z_meta", "p_meta", "FDR_BH", "enriched_in"])
                tmp["pair"] = mf.stem.replace("meta_enrichment_", ""); frames.append(tmp)
            except Exception: continue
        if frames:
            meta_all = pd.concat(frames, ignore_index=True)
            meta_all["absZ"] = meta_all["Z_meta"].abs()
            summary = (meta_all.groupby("Regulon", as_index=False)
                       .agg(best_FDR=("FDR_BH", "min"), best_absZ=("absZ", "max")))
            sig = summary[summary["best_FDR"] <= fdr]
            if not sig.empty:
                ranked = sig.sort_values(["best_FDR", "best_absZ"], ascending=[True, False])["Regulon"].tolist()
    if not ranked:
        tmp = df_merged.copy(); means = tmp[[f"MeanAUC_{g}" for g in GROUPS]].fillna(0.0).values
        var = np.var(means, axis=1); tmp["_varAUC"] = var
        ranked = tmp.sort_values("_varAUC", ascending=False)["Regulon"].tolist()
    return ranked[:top_n]

def make_multigroup_dotplot(df: pd.DataFrame, selected: list[str], out_png: Path, thr=0.10):
    sns.set_style("white"); plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 300, "font.size": 8,
                                                 "legend.frameon": False, "axes.grid": False})
    sub = df[df["Regulon"].isin(selected)].copy()
    rows = []
    for _, r in sub.iterrows():
        for g in GROUPS:
            rows.append(dict(Regulon=r["Regulon"], Group=g,
                             MeanAUC=r.get(f"MeanAUC_{g}", np.nan),
                             FracHigh=r.get(f"FracHigh_{g}", np.nan)))
    plot_df = pd.DataFrame(rows)
    plot_df["Regulon"] = pd.Categorical(plot_df["Regulon"], categories=list(reversed(selected)), ordered=True)
    plot_df["Group"]   = pd.Categorical(plot_df["Group"], categories=GROUPS, ordered=True)
    vmin = float(np.nanmin(plot_df["MeanAUC"])) if np.isfinite(np.nanmin(plot_df["MeanAUC"])) else 0.0
    vmax = float(np.nanmax(plot_df["MeanAUC"])) if np.isfinite(np.nanmax(plot_df["MeanAUC"])) else 1.0
    if vmin == vmax: vmax = vmin + 1e-6
    norm = Normalize(vmin=vmin, vmax=vmax)
    h = max(4.0, 0.28 * len(selected) + 1.3); fig, ax = plt.subplots(figsize=(5.0, h))
    smin, smax = 20, 300
    sns.scatterplot(data=plot_df, x="Group", y="Regulon", hue="MeanAUC", palette="viridis", hue_norm=norm,
                    size="FracHigh", sizes=(smin, smax), linewidth=0, edgecolor="none", ax=ax, legend=False)
    ax.set_ylabel(""); ax.set_xlabel("")
    ax.set_title(f"CTpre groups: color=mean AUCell, size=% cells ≥ {thr:g} (≥7/10 runs)", fontsize=9)
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
    size_leg = Legend(ax, handles, labels, title="% cells ≥ thr",
                      bbox_to_anchor=(1.03, 1.0), loc="upper left", borderaxespad=0., fontsize=6, title_fontsize=6)
    ax.add_artist(size_leg); ax.tick_params(axis="x", labelsize=8)
    _dot_save_fig(fig, out_png); plt.close(fig)

# ---------- main ----------
def main():
    ap = argparse.ArgumentParser(description="SCENIC CTpre groups (union run, 10×) with aggregation, meta‑enrichment, and dotplots.")
    ap.add_argument("--in-h5ad", required=True, type=Path)
    ap.add_argument("--out-base", required=True, type=Path)
    ap.add_argument("--ranking-db", required=True, type=Path)
    ap.add_argument("--motif-table", required=True, type=Path)
    ap.add_argument("--tf-list", required=True, type=Path)
    ap.add_argument("--n-workers", type=int, default=12)
    ap.add_argument("--threads-per-worker", type=int, default=1)
    ap.add_argument("--n-iter", type=int, default=10)
    ap.add_argument("--dot-threshold", type=float, default=0.10)
    ap.add_argument("--dot-min-runs", type=int, default=7)
    ap.add_argument("--top-n", nargs="+", type=int, default=[30, 15])
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(level=logging.INFO if args.verbose else logging.WARNING,
                        format="%(asctime)s %(levelname)-8s %(message)s", datefmt="%H:%M:%S")

    OUT_BASE = args.out_base; ALL_DIR = OUT_BASE / "scenic_CTpre_all"
    DOT_OUT_DIR = OUT_BASE / "CTpre_groups"; DOT_OUT_DIR.mkdir(parents=True, exist_ok=True)

    log.info(f"Reading AnnData → {args.in_h5ad}")
    adata = ad.read_h5ad(args.in_h5ad); log.info(f"AnnData: {adata.n_obs} cells, {adata.n_vars} genes")

    # Build annotations
    if "sample_id" in adata.obs:
        def _grp_from_sample(x): return "pre" if str(x) in {"NTX_1", "NTX_2", "NTX_3"} else "post"
        adata.obs["NTX_group"] = adata.obs["sample_id"].map(_grp_from_sample)
    else:
        adata.obs["NTX_group"] = "pre"
    if "Human_Tonsil_Classifier_Reclassification" in adata.obs:
        adata.obs["cell_type"] = adata.obs["Human_Tonsil_Classifier_Reclassification"]
    else:
        adata.obs["cell_type"] = "unknown"

    def merged_label(row) -> str | None:
        ct = str(row["cell_type"]).strip().lower(); grp = row["NTX_group"]
        if ct == "mbc derived igg+ pc" and grp == "pre": return "MBC_IgG_PC_pre"
        if ct == "mature igg+ pc" and grp == "pre": return "Mature_IgG_PC_pre"
        if (grp == "pre") and (("csmbc" in ct) or all(k in ct for k in ("class", "switched", "memory"))):
            return "csMBC_pre"
        return None

    adata.obs["merged_label"] = adata.obs.apply(merged_label, axis=1)
    group2ids = {g: adata.obs_names[adata.obs["merged_label"] == g].astype(str).tolist() for g in GROUPS}
    all_ids = sorted(set().union(*group2ids.values()))
    if len(all_ids) < 2: raise SystemExit("CTpre_all: <2 total cells in union — aborting.")

    # Expression matrix (prefer raw)
    X = adata.raw.to_adata() if adata.raw else adata.copy()
    if "counts" in X.layers: X.X = X.layers["counts"].copy()
    X.var["highly_variable"] = True
    expr = pd.DataFrame(_safe_toarray(X.X), index=X.obs_names, columns=X.var_names)
    expr = expr.loc[all_ids]  # subset to union

    # Run SCENIC
    seeds = [1337 + i for i in range(args.n_iter)]
    run_scenic(expr, ranking_db=args.ranking_db, motif_tbl=args.motif_table, tf_list=args.tf_list,
               out_dir=ALL_DIR, n_iter=args.n_iter, n_workers=args.n_workers,
               threads_per_worker=args.threads_per_worker, seeds=seeds)

    # Aggregate
    aggregate_CTpre_all(ALL_DIR)

    # Meta‑enrichment
    if any(len(group2ids[g]) >= 3 for g in GROUPS):
        if not HAS_SCIPY:
            log.warning("SciPy not found — meta p-values/FDR will be NaN; install scipy for full stats.")
        meta_enrichment_CTpre_all(ALL_DIR, group2ids, min_cells=3)

    # Dotplots
    try:
        df = _dot_load_multigroup_from_all(ALL_DIR, group2ids, thr=args.dot_threshold)
        df = df[df["n_runs_detected"] >= args.dot_min_runs].copy()
        _dot_write_csv(df, DOT_OUT_DIR / "dotplot_source_table_CTpre_groups.csv")
        for n in args.top_n:
            selected = _rank_regulons_multigroup(df, ALL_DIR / "_aggregate", top_n=n, fdr=0.10)
            if not selected: 
                log.warning(f"No regulons selected for top_n={n}. Skipping.")
                continue
            make_multigroup_dotplot(df, selected, DOT_OUT_DIR / f"dotplot_CTpre_groups_top{n}.png",
                                    thr=args.dot_threshold)
        log.info("CTpre groups dotplots complete ✔")
    except Exception as e:
        log.warning(f"[DOTPLOT] Skipped due to error: {e}")

if __name__ == "__main__":
    main()
