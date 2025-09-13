#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
SCENIC end-to-end pipeline (single executable)
==============================================

Stages
------
A) Reference T-cell SCENIC -> AUCell(all) -> CAR subset (CD4/CD8), multi-run
B) LLPCs union SCENIC (multi-run) + aggregation + pre/post meta + dotplots
C) CT pre-only union SCENIC (multi-run) + pairwise meta + dotplots
D) NetworkX consensus & overlaps across SCENIC_CAR runs (CD4)

This template is path-agnostic:  outputs
default under `results/`. It preserves the logic of your original scripts.

Required tools: arboreto, pyscenic, ctxcore, dask, networkx, scanpy/anndata, seaborn.
"""

from __future__ import annotations
import argparse
import inspect
import logging
import math
import os
import pickle
import re
import warnings
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple, Dict

# ---- Core scientific stack
import numpy as np
import pandas as pd
import anndata as ad

# GRN / SCENIC
from arboreto.algo import grnboost2
from ctxcore.rnkdb import FeatherRankingDatabase as RankingDatabase
from pyscenic.utils import modules_from_adjacencies
from pyscenic.prune import prune2df, df2regulons
from pyscenic.aucell import aucell as scenic_aucell

# Parallelism
from dask.distributed import LocalCluster, Client
from dask.diagnostics import ProgressBar

# Plotting (dotplots)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.legend import Legend

# Network analysis
import networkx as nx

# Optional stats
try:
    from scipy.stats import mannwhitneyu, norm
    HAS_SCIPY = True
except Exception:
    HAS_SCIPY = False


# ──────────────────────────────────────────────────────────────────────
# Defaults (generic; override via CLI)
# ──────────────────────────────────────────────────────────────────────
# Output roots
DEF_CAR_ROOT   = Path("results/SCENIC_CAR_runs")
DEF_LLP_BASE   = Path("results/SCENIC_LLP_union")
DEF_CTPRE_BASE = Path("results/SCENIC_CTpre_union")
DEF_NX_DIR     = Path("results/SCENIC_CAR_runs/NetworkX_CD4_vsIterations")

# A) ref->CAR (multi-run)
DEF_ITERS_REFCAR = 20
DEF_WORKERS_A    = 8
DEF_SEED_BASE    = 42000
DEF_MIN_REGULON  = 10
DEF_Q_ACTIVITY   = 0.95
DEF_MIN_FRAC_CD4 = 0.30
DEF_MIN_FRAC_CD8 = 0.30

# Regexes to pick reference CD4/CD8 T-cells from the **reference** atlas
CD4_REF_PATTERNS = [r"effector\s*helper\s*T", r"regulatory\s*T|^t[-\s]*reg"]
CD8_REF_PATTERNS = [
    r"\btrm\b.*cytotoxic", r"\btemra\b.*cytotoxic",
    r"tissue[-\s]*resident.*cytotoxic", r"effector.*memory.*RA.*cytotoxic"
]

# Dataset-specific CAR barcodes (can be overridden with --car-cd4-file / --car-cd8-file)
CAR_CD4_BARCODES = [
    "NTX_5_TGTTCCGAGCCACCTG-1","NTX_5_GAAACTCCACATGGGA-1","NTX_5_TGAGGGATCTTGCAAG-1",
    "NTX_4_CGTTCTGTCACCTCGT-1","NTX_4_GGACAAGTCAAACGGG-1","NTX_4_CGACTTCAGCCGGTAA-1",
    "NTX_5_TGTATTCCATGAACCT-1","NTX_5_TGGTTAGCAGCCTTTC-1","NTX_5_ACTTGTTGTTATGTGC-1",
    "NTX_5_CGCTATCTCACAACGT-1","NTX_5_TTCGAAGTCGTCACGG-1",
]
CAR_CD8_BARCODES = [
    "NTX_5_TTTATGCAGCCGATTT-1","NTX_5_TTTATGCTCTAACTCT-1","NTX_5_ATTACTCAGCGTAATA-1",
]

# B/C) union runs & plots
DEF_ITERS_BC     = 10
DEF_WORKERS_BC   = 12
DOT_AUC_THR      = 0.10
DOT_MIN_RUNS_LLP = 8
DOT_MIN_RUNS_CTP = 7
DOT_META_FDR     = 0.10
MIN_CELLS_GROUP  = 3

CT_GROUPS = ["MBC_IgG_PC_pre", "Mature_IgG_PC_pre", "csMBC_pre"]

# D) NetworkX CD4
PREFER_ADJ_CSV        = True
TOP_EDGES_PER_TF      = 50
TARGET_EDGE_RANGE     = (1500, 20000)
OVERLAP_FRAC          = 0.60
EDGE_OVERLAP_FR       = 0.60
STRICT_NODE_FR        = 0.95
STRICT_AUTO_MIN_EDGES = 50
NODE_OVERLAP_FRACS    = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95, 1.00]
LABEL_TOP             = 50
SPRING_SEED           = 4_572_321


# ──────────────────────────────────────────────────────────────────────
# Logging
# ──────────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("SCENIC-ALL")


# ──────────────────────────────────────────────────────────────────────
# Utilities
# ──────────────────────────────────────────────────────────────────────
def _safe_toarray(X): return X.toarray() if hasattr(X, "toarray") else X

def _require_path(p: Optional[Path], what: str, for_stage: str):
    if p is None:
        raise SystemExit(f"[ERROR] --{what} is required for stage '{for_stage}'.")

def _load_lines(path: Path) -> List[str]:
    """Read a file of barcodes (one per line or CSV/TSV)."""
    txt = Path(path).read_text().strip().splitlines()
    out = []
    for row in txt:
        row = row.strip()
        if not row: 
            continue
        if "," in row or "\t" in row:
            out += [x.strip() for x in re.split(r"[,\t]", row) if x.strip()]
        else:
            out.append(row)
    return out


# ──────────────────────────────────────────────────────────────────────
# Expression prep & helpers
# ──────────────────────────────────────────────────────────────────────
def prepare_expr(adata: ad.AnnData) -> pd.DataFrame:
    X = adata.raw.to_adata() if adata.raw else adata.copy()
    if "counts" in X.layers: 
        X.X = X.layers["counts"].copy()
    X.var["highly_variable"] = True
    mat = _safe_toarray(X.X)
    return pd.DataFrame(mat, index=X.obs_names, columns=X.var_names)

def load_tf_list(path: Path, genes: Sequence[str]) -> List[str]:
    tfs = pd.read_csv(path, header=None)[0].astype(str).tolist()
    return [t for t in tfs if t in genes]

def find_celltype_column(adata: ad.AnnData) -> str:
    for c in ["celltypist","celltypist_prediction","majority_voting",
              "predicted_labels","celltype","cell_type","annotation"]:
        if c in adata.obs: 
            return c
    raise SystemExit("Cell-type column not found; pass --celltype-col to override.")

def match_indices(adata: ad.AnnData, column: str, patterns: Iterable[str]) -> List[str]:
    ser = adata.obs[column].astype(str)
    mask = np.zeros(len(ser), dtype=bool)
    for pat in patterns:
        mask |= ser.str.contains(pat, case=False, regex=True, na=False)
    return ser.index[mask].tolist()

def normalize_tf_name(reg_name: str) -> str:
    return reg_name.split("(")[0]

def _bh_fdr(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    if n == 0: 
        return p
    order = np.argsort(p)
    ranks = np.arange(1, n + 1, dtype=float)
    q = np.empty(n, dtype=float)
    q[order] = (p[order] * n) / ranks
    for i in range(n - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.minimum(q, 1.0)

def _signed_z_from_two_sided_p(p_two_sided: float, sign: float) -> float:
    if not HAS_SCIPY: 
        return np.nan
    p = float(min(max(p_two_sided, 1e-300), 1 - 1e-300))
    z_abs = norm.isf(p / 2.0)
    return math.copysign(z_abs,  1.0 if sign > 0 else (-1.0 if sign < 0 else 1.0))


# ──────────────────────────────────────────────────────────────────────
# A) SCENIC reference->CAR
# ──────────────────────────────────────────────────────────────────────
def scenic_infer(expr: pd.DataFrame, tfs: List[str], out_dir: Path, label: str,
                 gene_db: Path, mc10_tbl: Path, seed: Optional[int],
                 n_workers: int, threads_per_worker: int, min_regulon_size: int,
                 force: bool) -> Tuple[list, pd.DataFrame]:
    out_dir.mkdir(parents=True, exist_ok=True)
    reg_pkl = out_dir / f"regulons_{label}.pkl"
    auc_csv = out_dir / f"auc_mtx_{label}.csv"

    if reg_pkl.exists() and auc_csv.exists() and not force:
        regs = pickle.load(open(reg_pkl, "rb"))
        auc = pd.read_csv(auc_csv, index_col=0)
        return regs, auc

    cluster = LocalCluster(n_workers=n_workers, threads_per_worker=threads_per_worker, dashboard_address=None)
    cli = Client(cluster)
    try:
        log.info(f"[{label}] GRNBoost2 on {expr.shape[0]} cells × {expr.shape[1]} genes (seed={seed})")
        kwargs = dict(gene_names=expr.columns, tf_names=tfs, client_or_address=cli, verbose=False)
        if "seed" in inspect.signature(grnboost2).parameters and seed is not None:
            kwargs["seed"] = seed
        adj = grnboost2(expr, **kwargs)
        adj.to_csv(out_dir / f"adjacencies_{label}.csv", index=False)

        mods = [m for m in modules_from_adjacencies(adj, expr)
                if len(getattr(m, "genes", getattr(m, "gene2weight", []))) >= min_regulon_size]
        pickle.dump(mods, open(out_dir / f"modules_{label}.pkl", "wb"))

        db = RankingDatabase(gene_db, name="hg38_mc10")
        with ProgressBar():
            pruned = prune2df([db], mods, mc10_tbl, num_workers=1, client_or_address="custom_multiprocessing")
        pruned.to_csv(out_dir / f"pruned_motifs_{label}.csv", index=False)

        regs = df2regulons(pruned)
        pickle.dump(regs, open(reg_pkl, "wb"))

        auc = scenic_aucell(expr, regs, num_workers=n_workers)
        auc.to_csv(auc_csv)
    finally:
        cli.close(); cluster.close()
    return regs, auc


def active_tfs_from_auc_all(auc_all: pd.DataFrame, subset_cells: List[str],
                            q: float, min_frac: float) -> Tuple[List[str], pd.Series]:
    subset_cells = [c for c in subset_cells if c in auc_all.index]
    if not subset_cells:
        return [], pd.Series(dtype=float)
    cutoff = auc_all.quantile(q, axis=0)
    is_act = (auc_all.loc[subset_cells].ge(cutoff)).astype(int)
    frac_active = is_act.mean(axis=0)
    active_regs = frac_active[frac_active >= min_frac].index.tolist()
    active_tfs = sorted({normalize_tf_name(r) for r in active_regs})
    return active_tfs, frac_active


def run_ref_to_car(adata_path: Path, car_root: Path, gene_db: Path, mc10_tbl: Path, tf_file: Path,
                   iters: int, seed_base: int, n_workers: int, min_regulon_size: int,
                   q: float, min_frac_cd4: float, min_frac_cd8: float,
                   force: bool, celltype_col: Optional[str],
                   car_cd4_cells: List[str], car_cd8_cells: List[str]) -> None:
    car_root.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(adata_path)
    log.info(f"[REF->CAR] Loaded {adata.n_obs} cells × {adata.n_vars} genes")
    col = celltype_col or find_celltype_column(adata)

    idx_cd4 = match_indices(adata, col, CD4_REF_PATTERNS)
    idx_cd8 = match_indices(adata, col, CD8_REF_PATTERNS)
    if len(idx_cd4) == 0 or len(idx_cd8) == 0:
        raise SystemExit(f"No reference matches. Check column '{col}' and regex patterns.")

    expr_cd4 = prepare_expr(adata[idx_cd4])
    expr_cd8 = prepare_expr(adata[idx_cd8])
    expr_all = prepare_expr(adata)

    tfs_cd4 = load_tf_list(tf_file, expr_cd4.columns)
    tfs_cd8 = load_tf_list(tf_file, expr_cd8.columns)

    for i in range(1, iters + 1):
        run_dir = car_root / f"SCENIC_CAR_{i}"
        out_ref_cd4 = run_dir / "REF_CD4"
        out_ref_cd8 = run_dir / "REF_CD8"
        out_car_cd4 = run_dir / "CAR_CD4"
        out_car_cd8 = run_dir / "CAR_CD8"
        out_car_cd4.mkdir(parents=True, exist_ok=True)
        out_car_cd8.mkdir(parents=True, exist_ok=True)

        # CD4 reference → all cells AUCell → pick CAR CD4 subset
        regs_cd4, _ = scenic_infer(
            expr_cd4, tfs_cd4, out_ref_cd4, "REF_CD4",
            gene_db, mc10_tbl, seed_base + i, n_workers, 1, min_regulon_size, force
        )
        auc_all_cd4 = scenic_aucell(expr_all, regs_cd4, num_workers=n_workers)
        auc_all_cd4.to_csv(out_ref_cd4 / "auc_mtx_ALLcells_REF_CD4.csv")
        act_tfs_cd4, frac_cd4 = active_tfs_from_auc_all(auc_all_cd4, car_cd4_cells, q, min_frac_cd4)
        pickle.dump([r for r in regs_cd4 if normalize_tf_name(getattr(r, "name","")) in set(act_tfs_cd4)],
                    open(out_car_cd4 / "regulons_CAR_CD4.pkl","wb"))
        keep_cells_cd4 = [c for c in car_cd4_cells if c in auc_all_cd4.index]
        auc_all_cd4.loc[keep_cells_cd4].to_csv(out_car_cd4/"auc_mtx_CAR_CD4.csv")
        pd.Series(act_tfs_cd4).to_csv(out_car_cd4/"active_TFs_in_CAR_CD4.txt", index=False, header=False)
        frac_cd4.rename("frac_active_in_subset").to_frame().to_csv(out_car_cd4/"regulon_activity_fraction_CAR_CD4.csv")

        # CD8 reference → all cells AUCell → pick CAR CD8 subset
        regs_cd8, _ = scenic_infer(
            expr_cd8, tfs_cd8, out_ref_cd8, "REF_CD8",
            gene_db, mc10_tbl, seed_base + 10_000 + i, n_workers, 1, min_regulon_size, force
        )
        auc_all_cd8 = scenic_aucell(expr_all, regs_cd8, num_workers=n_workers)
        auc_all_cd8.to_csv(out_ref_cd8 / "auc_mtx_ALLcells_REF_CD8.csv")
        act_tfs_cd8, frac_cd8 = active_tfs_from_auc_all(auc_all_cd8, car_cd8_cells, q, min_frac_cd8)
        pickle.dump([r for r in regs_cd8 if normalize_tf_name(getattr(r, "name","")) in set(act_tfs_cd8)],
                    open(out_car_cd8 / "regulons_CAR_CD8.pkl","wb"))
        keep_cells_cd8 = [c for c in car_cd8_cells if c in auc_all_cd8.index]
        auc_all_cd8.loc[keep_cells_cd8].to_csv(out_car_cd8/"auc_mtx_CAR_CD8.csv")
        pd.Series(act_tfs_cd8).to_csv(out_car_cd8/"active_TFs_in_CAR_CD8.txt", index=False, header=False)
        frac_cd8.rename("frac_active_in_subset").to_frame().to_csv(out_car_cd8/"regulon_activity_fraction_CAR_CD8.csv")

        log.info(f"[REF->CAR] Iteration {i}/{iters} complete.")


# ──────────────────────────────────────────────────────────────────────
# Generic SCENIC multi-iteration (used by LLP_all and CTpre_all)
# ──────────────────────────────────────────────────────────────────────
def scenic_multi_run(adata: ad.AnnData, out_base: Path, label: str, n_iter: int,
                     gene_db: Path, mc10_tbl: Path, tf_file: Path,
                     dask_workers: int, threads_per_worker: int, force: bool) -> None:
    base = out_base / f"scenic_{label}"
    base.mkdir(parents=True, exist_ok=True)

    X = adata.raw.to_adata() if adata.raw else adata.copy()
    if "counts" in X.layers: 
        X.X = X.layers["counts"].copy()
    X.var["highly_variable"] = True
    tfs = [t for t in pd.read_csv(tf_file, header=None)[0].astype(str).tolist() if t in X.var_names]
    expr = pd.DataFrame(_safe_toarray(X.X), index=X.obs_names, columns=X.var_names)

    cluster = LocalCluster(n_workers=dask_workers, threads_per_worker=threads_per_worker, processes=True, dashboard_address=None)
    cli = Client(cluster)
    db = RankingDatabase(gene_db, name="hg38_mc10")
    seeds = [1337 + i for i in range(n_iter)]
    has_seed = "seed" in inspect.signature(grnboost2).parameters
    try:
        for i in range(1, n_iter + 1):
            iter_out = base / f"iter_{i:02d}"
            iter_out.mkdir(parents=True, exist_ok=True)
            reg_pkl = iter_out / f"regulons_{label}_iter{i:02d}.pkl"
            auc_csv = iter_out / f"auc_mtx_{label}_iter{i:02d}.csv"
            if reg_pkl.exists() and auc_csv.exists() and not force:
                log.info(f"[{label}] iter_{i:02d} exists; skipping (use --force to recompute).")
                continue

            kwargs = dict(gene_names=expr.columns, tf_names=tfs, client_or_address=cli, verbose=False)
            if has_seed: 
                kwargs["seed"] = seeds[i-1]
            adj = grnboost2(expr, **kwargs)
            adj.to_csv(iter_out / f"adjacencies_{label}_iter{i:02d}.csv", index=False)

            mods = [m for m in modules_from_adjacencies(adj, expr) if len(getattr(m,"genes",getattr(m,"gene2weight",[])))>=10]
            pickle.dump(mods, open(iter_out / f"modules_{label}_iter{i:02d}.pkl","wb"))

            with ProgressBar():
                pruned = prune2df([db], mods, mc10_tbl, num_workers=1, client_or_address="custom_multiprocessing")
            pruned.to_csv(iter_out / f"pruned_motifs_{label}_iter{i:02d}.csv", index=False)

            regs = df2regulons(pruned)
            pickle.dump(regs, open(reg_pkl,"wb"))

            auc = scenic_aucell(expr, regs, num_workers=1)
            auc.to_csv(auc_csv)

            # Quick summaries
            (pd.DataFrame([(getattr(r,"name",""), len(getattr(r,"gene2weight",{}))) for r in regs],
                          columns=["Regulon","NumTargets"])
             .sort_values("NumTargets", ascending=False)
             .to_csv(iter_out / f"top_regulons_by_targetcount_{label}_iter{i:02d}.csv", index=False))
            (auc.mean(axis=0).rename("MeanAUC").to_frame().reset_index()
             .rename(columns={"index":"Regulon"})
             .sort_values("MeanAUC", ascending=False)
             .to_csv(iter_out / f"top_regulons_by_meanAUC_{label}_iter{i:02d}.csv", index=False))
            log.info(f"[{label}] iter_{i:02d} complete.")
    finally:
        cli.close(); cluster.close()


def iter_dirs(base: Path, label: str) -> List[Path]:
    root = base / f"scenic_{label}"
    if not root.exists(): 
        return []
    iters = [p for p in root.iterdir() if p.is_dir() and p.name.startswith("iter_")]
    def _key(p):
        m = re.search(r"(\d+)", p.name)
        return int(m.group(1)) if m else 0
    return sorted(iters, key=_key)


def aggregate_label(base: Path, label: str) -> None:
    base_label = base / f"scenic_{label}"
    agg_dir = base_label / "_aggregate"
    agg_dir.mkdir(parents=True, exist_ok=True)

    runs = iter_dirs(base, label)
    if not runs:
        log.warning(f"[{label}] no iterations to aggregate.")
        return

    auc_series_per_iter = {}
    tf_sets_per_iter = {}
    regulon_target_sets = defaultdict(list)
    run_metrics_rows = []
    edge_counts = Counter()
    edge_imp_values = defaultdict(list)

    def _norm_cols(df):
        lower = {c.lower(): c for c in df.columns}
        tf_col = lower.get("tf", lower.get("regulator", None))
        target_col = lower.get("target", lower.get("gene", None))
        imp_col = lower.get("importance", lower.get("weight", lower.get("score", None)))
        if tf_col is None or target_col is None:
            raise ValueError("TF/target columns not found in adjacency CSV.")
        return tf_col, target_col, imp_col

    for run in runs:
        rid = run.name
        adj_files = list(run.glob(f"adjacencies_{label}_iter*.csv"))
        if adj_files:
            adj = pd.read_csv(adj_files[0])
            TF, TG, IMP = _norm_cols(adj)
            for row in adj[[TF, TG]].itertuples(index=False):
                edge_counts[(str(row[0]), str(row[1]))] += 1
            if IMP is not None:
                for row in adj[[TF, TG, IMP]].itertuples(index=False):
                    edge_imp_values[(str(row[0]), str(row[1]))].append(float(row[2]))
            n_edges = len(adj)
        else:
            n_edges = np.nan

        mod_files = list(run.glob(f"modules_{label}_iter*.pkl"))
        if mod_files:
            try:
                mods = pickle.load(open(mod_files[0], "rb")); n_modules = len(mods)
            except Exception: n_modules = np.nan
        else: n_modules = np.nan

        reg_files = list(run.glob(f"regulons_{label}_iter*.pkl"))
        tf_set = set(); target_sizes = []
        if reg_files:
            try:
                regs = pickle.load(open(reg_files[0], "rb"))
                for r in regs:
                    tf = re.split(r"\s*\(", getattr(r, "name", "") or str(r))[0].strip()
                    tf_set.add(tf)
                    g2w = getattr(r, "gene2weight", {})
                    target_sizes.append(len(g2w))
                    if len(g2w) > 0:
                        regulon_target_sets[tf].append(set(g2w.keys()))
                n_regulons = len(regs); n_unique_tfs = len(tf_set)
                median_targets = float(np.median(target_sizes)) if target_sizes else np.nan
            except Exception:
                n_regulons = np.nan; n_unique_tfs = np.nan; median_targets = np.nan
        else:
            n_regulons = np.nan; n_unique_tfs = np.nan; median_targets = np.nan
        tf_sets_per_iter[rid] = tf_set

        auc_sum_files = list(run.glob(f"top_regulons_by_meanAUC_{label}_iter*.csv"))
        if auc_sum_files:
            auc_sum = pd.read_csv(auc_sum_files[0])
            if {"Regulon","MeanAUC"}.issubset(auc_sum.columns):
                s = pd.Series(auc_sum["MeanAUC"].values, index=auc_sum["Regulon"].astype(str))
                auc_series_per_iter[rid] = s

        auc_mtx_files = list(run.glob(f"auc_mtx_{label}_iter*.csv"))
        if auc_mtx_files:
            try:
                with open(auc_mtx_files[0], "r") as f:
                    n_cells = sum(1 for _ in f) - 1
            except Exception:
                n_cells = np.nan
        else: 
            n_cells = np.nan

        run_metrics_rows.append(dict(
            Iteration=rid, n_edges=n_edges, n_modules=n_modules,
            n_regulons=n_regulons, n_unique_TFs=n_unique_tfs,
            median_targets_per_regulon=median_targets, n_cells=n_cells))

    # AUC summaries
    if auc_series_per_iter:
        auc_wide = pd.DataFrame(auc_series_per_iter).sort_index()
        auc_wide.to_csv(agg_dir / f"AUC_mean_per_regulon_wide_{label}.csv")
        stats_reg = pd.DataFrame({
            "n_runs_detected": auc_wide.notna().sum(axis=1),
            "mean_meanAUC": auc_wide.mean(axis=1, skipna=True),
            "std_meanAUC": auc_wide.std(axis=1, ddof=1, skipna=True),
            "min_meanAUC": auc_wide.min(axis=1, skipna=True),
            "max_meanAUC": auc_wide.max(axis=1, skipna=True),
        })
        stats_reg["cv_meanAUC"] = stats_reg["std_meanAUC"]/stats_reg["mean_meanAUC"]
        stats_reg.sort_values("mean_meanAUC", ascending=False).to_csv(agg_dir / f"auc_summary_by_regulon_{label}.csv")
        auc_by_tf = auc_wide.copy()
        auc_by_tf.index = [re.split(r"\s*\(", i)[0].strip() for i in auc_by_tf.index]
        auc_by_tf = auc_by_tf.groupby(auc_by_tf.index).max()
        stats_tf = pd.DataFrame({
            "n_runs_detected": auc_by_tf.notna().sum(axis=1),
            "mean_meanAUC": auc_by_tf.mean(axis=1, skipna=True),
            "std_meanAUC": auc_by_tf.std(axis=1, ddof=1, skipna=True),
            "min_meanAUC": auc_by_tf.min(axis=1, skipna=True),
            "max_meanAUC": auc_by_tf.max(axis=1, skipna=True),
        })
        stats_tf["cv_meanAUC"] = stats_tf["std_meanAUC"]/stats_tf["mean_meanAUC"]
        stats_tf.sort_values("mean_meanAUC", ascending=False).to_csv(agg_dir / f"auc_summary_by_TF_{label}.csv")
        auc_wide.corr(min_periods=3).to_csv(agg_dir / f"auc_mean_correlation_across_runs_{label}.csv")

    # Jaccard TF presence
    run_ids = list(tf_sets_per_iter.keys())
    if run_ids:
        jacc = pd.DataFrame(index=run_ids, columns=run_ids, dtype=float)
        for a in run_ids:
            for b in run_ids:
                A, B = tf_sets_per_iter[a], tf_sets_per_iter[b]
                u = len(A|B); jacc.loc[a,b] = (len(A&B)/u) if u else np.nan
        jacc.to_csv(agg_dir / f"jaccard_regulon_TF_sets_{label}.csv")

    # Target set stability
    rows = []
    for tf, sets_list in regulon_target_sets.items():
        if not sets_list: 
            continue
        n_runs_tf = len(sets_list)
        comb = list(combinations(sets_list, 2))
        if comb:
            js = []
            for a,b in comb:
                u = a|b; js.append(len(a&b)/len(u) if len(u) else 1.0)
            mean_j = float(np.mean(js)) if js else np.nan
        else:
            mean_j = 1.0
        counts = Counter()
        for s in sets_list: 
            counts.update(s)
        core50 = sorted([g for g,c in counts.items() if c/n_runs_tf >= 0.5])
        core80 = sorted([g for g,c in counts.items() if c/n_runs_tf >= 0.8])
        rows.append(dict(TF=tf, runs_with_TF=n_runs_tf, mean_pairwise_jaccard=mean_j,
                         union_targets=len(set().union(*sets_list)),
                         core50_count=len(core50), core80_count=len(core80),
                         core50_targets=";".join(core50), core80_targets=";".join(core80)))
    if rows:
        (pd.DataFrame(rows).sort_values(["runs_with_TF","mean_pairwise_jaccard","core80_count"],
                                        ascending=[False,False,False])
         .to_csv(agg_dir / f"target_stability_by_TF_{label}.csv", index=False))

    # Edge stability
    if edge_counts:
        rows = []
        for (tf,tgt), cnt in edge_counts.items():
            vals = edge_imp_values.get((tf,tgt),[])
            rows.append(dict(TF=tf, target=tgt, n_runs=cnt,
                             frac_runs=cnt/len(run_ids) if run_ids else np.nan,
                             mean_importance=float(np.mean(vals)) if vals else np.nan,
                             std_importance=float(np.std(vals, ddof=1)) if len(vals)>=2 else np.nan))
        (pd.DataFrame(rows).sort_values(["n_runs","mean_importance"], ascending=[False,False])
         .to_csv(agg_dir / f"edge_stability_{label}.csv", index=False))

    # Run metrics
    pd.DataFrame(run_metrics_rows).set_index("Iteration").to_csv(agg_dir / f"run_metrics_{label}.csv")
    log.info(f"[{label}] aggregation complete → {agg_dir}")


# ──────────────────────────────────────────────────────────────────────
# B) LLP union: meta-enrichment & dotplots
# ──────────────────────────────────────────────────────────────────────
def llp_meta_enrichment(base: Path, pre_ids: List[str], post_ids: List[str], label: str="LLP_all") -> None:
    base_label = base / f"scenic_{label}"
    agg_dir = base_label / "_aggregate"; (agg_dir / "per_run").mkdir(parents=True, exist_ok=True)
    iter_paths = sorted(base_label.glob("iter_*/auc_mtx_*.csv"))
    if not iter_paths:
        log.warning("[LLP] no AUCell matrices found; skip meta.")
        return

    per_reg = defaultdict(lambda: dict(z=[], w=[], p=[], delta=[], rbc=[], n1=[], n2=[], k=0))
    for auc_path in iter_paths:
        auc = pd.read_csv(auc_path, index_col=0)
        pre  = [cid for cid in pre_ids  if cid in auc.index]
        post = [cid for cid in post_ids if cid in auc.index]
        n1, n2 = len(pre), len(post)
        if n1 < MIN_CELLS_GROUP or n2 < MIN_CELLS_GROUP:
            log.warning(f"[LLP] {auc_path.parent.name}: insufficient cells (pre={n1}, post={n2}).")
            continue
        rows = []
        a_mean = auc.loc[pre].mean(axis=0); b_mean = auc.loc[post].mean(axis=0)
        for regulon in auc.columns:
            pv = np.nan; U = np.nan; rbc = np.nan; z = np.nan
            delta = float(a_mean[regulon] - b_mean[regulon])
            sign = 1.0 if delta>0 else (-1.0 if delta<0 else 0.0)
            if HAS_SCIPY:
                try:
                    U, pv = mannwhitneyu(auc.loc[pre,regulon].values, auc.loc[post,regulon].values, alternative="two-sided")
                    rbc   = (2.0 * U / (n1 * n2)) - 1.0
                    z     = _signed_z_from_two_sided_p(pv, sign)
                except Exception: 
                    pass
            w = math.sqrt(n1*n2)
            acc = per_reg[regulon]
            acc["z"].append(z); acc["w"].append(w); acc["p"].append(pv); acc["delta"].append(delta); acc["rbc"].append(rbc)
            acc["n1"].append(n1); acc["n2"].append(n2); acc["k"] = acc["k"] + 1
            rows.append(dict(Regulon=regulon, pre_mean=a_mean[regulon], post_mean=b_mean[regulon],
                             delta=delta, U_stat=U, p_value=pv, r_rank_biserial=rbc, n_pre=n1, n_post=n2))
        pd.DataFrame(rows).sort_values(["p_value","delta"], ascending=[True,False]).to_csv(
            agg_dir / "per_run" / f"{auc_path.parent.name}_enrichment.csv", index=False)

    meta_rows = []
    for regulon, acc in per_reg.items():
        if acc["k"] == 0: 
            continue
        z_vals = np.array([z for z in acc["z"] if np.isfinite(z)], float)
        w_vals = np.array([w for z, w in zip(acc["z"], acc["w"]) if np.isfinite(z)], float)
        if HAS_SCIPY and z_vals.size>0 and w_vals.size>0:
            Z_meta = float(np.sum(w_vals*z_vals)/max(math.sqrt(np.sum(w_vals**2)), 1e-300))
            p_meta = float(2.0 * norm.sf(abs(Z_meta)))
            enriched_in = "LLP_pre" if Z_meta > 0 else "LLP_post"
        else:
            Z_meta = np.nan; p_meta = np.nan; enriched_in = "NA"
        meta_rows.append(dict(Regulon=regulon, k_runs=acc["k"], Z_meta=Z_meta, p_meta=p_meta,
                              enriched_in=enriched_in,
                              delta_median=float(np.nanmedian(acc["delta"])) if len(acc["delta"]) else np.nan,
                              r_rank_biserial_mean=float(np.nanmean(acc["rbc"])) if len(acc["rbc"]) else np.nan))
    meta_df = pd.DataFrame(meta_rows)
    if not meta_df.empty and HAS_SCIPY and meta_df["p_meta"].notna().any():
        meta_df["FDR_BH"] = _bh_fdr(meta_df["p_meta"].fillna(1.0).values)
    else:
        meta_df["FDR_BH"] = np.nan
    meta_df.sort_values(["FDR_BH","Z_meta","delta_median"], ascending=[True,False,False]).to_csv(
        agg_dir / "regulon_enrichment_pre_vs_post_meta.csv", index=False)
    log.info("[LLP] meta-enrichment written.")


def _nice_pct_levels(vals, k=3):
    v = np.asarray(vals, dtype=float) * 100.0
    v = v[np.isfinite(v)]
    if v.size == 0: 
        return [10,50,90]
    lo, hi = np.nanpercentile(v, 10), np.nanpercentile(v, 90)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo==hi:
        lo, hi = float(np.min(v)), float(np.max(v) if np.max(v)>np.min(v) else np.min(v)+1)
    raw = np.linspace(lo, hi, k)
    return [int(round(x)) for x in raw]


def llp_dotplots(base: Path, pre_ids: List[str], post_ids: List[str], label: str="LLP_all") -> None:
    all_dir = base / f"scenic_{label}"
    out_dir = base / "LLP_pre_vs_post"; out_dir.mkdir(parents=True, exist_ok=True)
    iters = iter_dirs(base, label)
    rows = []
    for idir in iters:
        auc_file = next(iter(idir.glob("auc_mtx_*.csv")), None)
        if auc_file is None: 
            continue
        auc = pd.read_csv(auc_file, index_col=0)
        pre  = [cid for cid in pre_ids  if cid in auc.index]
        post = [cid for cid in post_ids if cid in auc.index]
        if len(pre)==0 or len(post)==0: 
            continue
        it = pd.DataFrame({
            "Regulon": auc.columns.astype(str),
            "MeanAUC_pre": auc.loc[pre].mean(axis=0).values,
            "FracHigh_pre": (auc.loc[pre] >= DOT_AUC_THR).mean(axis=0).values,
            "MeanAUC_post": auc.loc[post].mean(axis=0).values,
            "FracHigh_post": (auc.loc[post] >= DOT_AUC_THR).mean(axis=0).values,
            "Iteration": idir.name,
        })
        rows.append(it)
    if not rows:
        log.warning("[LLP] dotplot skipped; no AUCell CSVs usable.")
        return
    df_all = pd.concat(rows, ignore_index=True)
    agg = (df_all.groupby("Regulon")
           .agg(MeanAUC_pre=("MeanAUC_pre","mean"), FracHigh_pre=("FracHigh_pre","mean"),
                MeanAUC_post=("MeanAUC_post","mean"), FracHigh_post=("FracHigh_post","mean"),
                n_runs=("Iteration","nunique"))).reset_index()
    keep = (agg["n_runs"] >= DOT_MIN_RUNS_LLP)
    df = agg[keep].copy()
    df["DeltaAUC"] = df["MeanAUC_post"] - df["MeanAUC_pre"]
    out_dir.joinpath("dotplot_source_table_LLP_pre_vs_LLP_post.csv").write_text(df.to_csv(index=False))

    meta_path = all_dir / "_aggregate" / "regulon_enrichment_pre_vs_post_meta.csv"
    for top_n in (30, 15):
        if meta_path.exists():
            meta = pd.read_csv(meta_path)
            if {"Regulon","Z_meta","FDR_BH"}.issubset(meta.columns):
                hit = meta[meta["FDR_BH"] <= DOT_META_FDR].copy()
                hit["absZ"] = hit["Z_meta"].abs()
                ranked = [r for r in hit.sort_values(["absZ"], ascending=False)["Regulon"].tolist() if r in set(df["Regulon"])]
            else:
                ranked = df.reindex(df["DeltaAUC"].abs().sort_values(ascending=False).index)["Regulon"].tolist()
        else:
            ranked = df.reindex(df["DeltaAUC"].abs().sort_values(ascending=False).index)["Regulon"].tolist()
        selected = ranked[:top_n]
        sub = df[df["Regulon"].isin(selected)].copy()
        plot_df = pd.DataFrame({
            "Regulon": np.repeat(sub["Regulon"].values, 2),
            "Phase":   np.tile(["Pre","Post"], len(sub)),
            "MeanAUC": np.concatenate([sub["MeanAUC_pre"].values, sub["MeanAUC_post"].values]),
            "FracHigh":np.concatenate([sub["FracHigh_pre"].values, sub["FracHigh_post"].values]),
        })
        yorder = list(reversed(selected))
        plot_df["Regulon"] = pd.Categorical(plot_df["Regulon"], categories=yorder, ordered=True)
        vmin, vmax = float(np.nanmin(plot_df["MeanAUC"])), float(np.nanmax(plot_df["MeanAUC"]))
        if vmin==vmax: 
            vmax = vmin + 1e-6
        norm = Normalize(vmin=vmin, vmax=vmax)
        h = max(3.8, 0.28*len(selected)+1.2)
        fig, ax = plt.subplots(figsize=(4.8, h))
        sns.scatterplot(data=plot_df, x="Phase", y="Regulon", hue="MeanAUC", palette="viridis",
                        hue_norm=norm, size="FracHigh", sizes=(20,300), linewidth=0, edgecolor="none", ax=ax, legend=False)
        ax.set_ylabel(""); ax.set_xlabel("")
        ax.set_title(f"Regulon activity (mean AUC color, % ≥ {DOT_AUC_THR:g} size) — LLP pre vs post", fontsize=9)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap("viridis"), norm=norm); sm.set_array([])
        fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.04).set_label("Mean AUCell", fontsize=6)
        levels = _nice_pct_levels(plot_df["FracHigh"], k=3)
        vmin_s, vmax_s = float(np.nanmin(plot_df["FracHigh"])), float(np.nanmax(plot_df["FracHigh"]))
        handles, labels = [], []
        for p in levels:
            area = 20 + ((p/100.0 - vmin_s) / max(vmax_s - vmin_s, 1e-9)) * (300-20)
            hndl = ax.scatter([], [], s=area, facecolor="none", edgecolor="k")
            handles.append(hndl); labels.append(f"{p}%")
        leg = Legend(ax, handles, labels, title="% cells ≥ thr",
                     bbox_to_anchor=(1.03,1.0), loc="upper left", borderaxespad=0., fontsize=6, title_fontsize=6)
        ax.add_artist(leg); fig.tight_layout()
        fig.savefig(out_dir / f"dotplot_LLP_pre_vs_LLP_post_top{top_n}.png", dpi=300)
        fig.savefig(out_dir / f"dotplot_LLP_pre_vs_LLP_post_top{top_n}.svg"); plt.close(fig)
    log.info("[LLP] dotplots saved.")


# ──────────────────────────────────────────────────────────────────────
# C) CT pre-only union: pairwise meta & dotplots (3 groups)
# ──────────────────────────────────────────────────────────────────────
def ctpre_build_groups(adata: ad.AnnData) -> Dict[str, List[str]]:
    if "sample_id" in adata.obs:
        def _grp(x): return "pre" if str(x) in {"NTX_1","NTX_2","NTX_3"} else "post"
        adata.obs["NTX_group"] = adata.obs["sample_id"].map(_grp)
    else:
        adata.obs["NTX_group"] = "pre"
    if "Human_Tonsil_Classifier_Reclassification" in adata.obs:
        adata.obs["cell_type"] = adata.obs["Human_Tonsil_Classifier_Reclassification"]
    else:
        adata.obs["cell_type"] = "unknown"
    def merged(row):
        ct = str(row["cell_type"]).strip().lower(); grp = row["NTX_group"]
        if ct == "mbc derived igg+ pc" and grp=="pre": return "MBC_IgG_PC_pre"
        if ct == "mature igg+ pc"     and grp=="pre": return "Mature_IgG_PC_pre"
        if (grp=="pre") and (("csmbc" in ct) or all(k in ct for k in ("class","switched","memory"))): return "csMBC_pre"
        return None
    adata.obs["merged_label"] = adata.obs.apply(merged, axis=1)
    g2ids = {}
    for g in CT_GROUPS:
        ids = adata.obs_names[adata.obs["merged_label"] == g].astype(str).tolist()
        g2ids[g] = ids
    return g2ids


def ctpre_meta_pairwise(base: Path, group2ids: Dict[str,List[str]], label: str="CTpre_all") -> None:
    base_label = base / f"scenic_{label}"
    agg_dir = base_label / "_aggregate"
    agg_dir.mkdir(parents=True, exist_ok=True)
    pairs = [(a,b) for a,b in combinations(CT_GROUPS, 2) if a in group2ids and b in group2ids]
    iter_paths = sorted(base_label.glob("iter_*/auc_mtx_*.csv"))
    if not iter_paths or not pairs:
        log.warning("[CTpre] no AUCell matrices or pairs; skip meta.")
        return
    per_pair_reg = { (a,b): defaultdict(lambda: dict(z=[],w=[],p=[],delta=[],rbc=[],k=0)) for (a,b) in pairs }
    for auc_path in iter_paths:
        auc = pd.read_csv(auc_path, index_col=0)
        for a,b in pairs:
            A = [cid for cid in group2ids[a] if cid in auc.index]
            B = [cid for cid in group2ids[b] if cid in auc.index]
            n1, n2 = len(A), len(B)
            if n1 < MIN_CELLS_GROUP or n2 < MIN_CELLS_GROUP:
                log.warning(f"[CTpre] {auc_path.parent.name} {a} vs {b}: insufficient cells.")
                continue
            a_mean = auc.loc[A].mean(axis=0); b_mean = auc.loc[B].mean(axis=0)
            rows = []
            for regulon in auc.columns:
                delta = float(a_mean[regulon] - b_mean[regulon])
                sign = 1.0 if delta>0 else (-1.0 if delta<0 else 0.0)
                U = np.nan; pv = np.nan; rbc = np.nan; z = np.nan
                if HAS_SCIPY:
                    try:
                        U, pv = mannwhitneyu(auc.loc[A,regulon].values, auc.loc[B,regulon].values, alternative="two-sided")
                        rbc   = (2.0 * U / (n1*n2)) - 1.0
                        z = _signed_z_from_two_sided_p(pv, sign)
                    except Exception: 
                        pass
                w = math.sqrt(n1*n2)
                acc = per_pair_reg[(a,b)][regulon]
                acc["z"].append(z); acc["w"].append(w); acc["p"].append(pv); acc["delta"].append(delta); acc["rbc"].append(rbc); acc["k"] += 1
                rows.append(dict(Regulon=regulon, mean_A=a_mean[regulon], mean_B=b_mean[regulon],
                                 delta=delta, U_stat=U, p_value=pv, r_rank_biserial=rbc, n_A=n1, n_B=n2))
            pr_dir = agg_dir / "per_run" / (re.sub(r"[^A-Za-z0-9]+","_", a)+"_vs_"+re.sub(r"[^A-Za-z0-9]+","_", b))
            pr_dir.mkdir(parents=True, exist_ok=True)
            pd.DataFrame(rows).sort_values(["p_value","delta"], ascending=[True,False]).to_csv(
                pr_dir / f"{auc_path.parent.name}_enrichment.csv", index=False)

    wide_rows = []
    for a,b in pairs:
        rows = []
        for regulon, acc in per_pair_reg[(a,b)].items():
            if acc["k"] == 0: 
                continue
            z_vals = np.array([z for z in acc["z"] if np.isfinite(z)], float)
            w_vals = np.array([w for z,w in zip(acc["z"], acc["w"]) if np.isfinite(z)], float)
            if HAS_SCIPY and z_vals.size>0 and w_vals.size>0:
                Z_meta = float(np.sum(w_vals*z_vals)/max(math.sqrt(np.sum(w_vals**2)), 1e-300))
                p_meta = float(2.0 * norm.sf(abs(Z_meta)))
                enriched_in = a if Z_meta > 0 else b
            else:
                Z_meta = np.nan; p_meta=np.nan; enriched_in="NA"
            rows.append(dict(Regulon=regulon, k_runs=acc["k"], Z_meta=Z_meta, p_meta=p_meta,
                             enriched_in=enriched_in,
                             delta_median=float(np.nanmedian(acc["delta"])) if len(acc["delta"]) else np.nan,
                             r_rank_biserial_mean=float(np.nanmean(acc["rbc"])) if len(acc["rbc"]) else np.nan))
        meta_df = pd.DataFrame(rows)
        if not meta_df.empty and HAS_SCIPY and meta_df["p_meta"].notna().any():
            meta_df["FDR_BH"] = _bh_fdr(meta_df["p_meta"].fillna(1.0).values)
        else:
            meta_df["FDR_BH"] = np.nan
        meta_df.sort_values(["FDR_BH","Z_meta","delta_median"], ascending=[True,False,False]).to_csv(
            agg_dir / f"meta_enrichment_{re.sub(r'[^A-Za-z0-9]+','_',a)}_vs_{re.sub(r'[^A-Za-z0-9]+','_',b)}.csv", index=False)
        for _, r in meta_df.iterrows():
            wide_rows.append(dict(Regulon=r["Regulon"], pair=f"{a}_vs_{b}",
                                  Z_meta=r["Z_meta"], p_meta=r["p_meta"], FDR_BH=r["FDR_BH"], enriched_in=r["enriched_in"]))
    if wide_rows:
        pd.DataFrame(wide_rows).to_csv(agg_dir / "meta_enrichment_all_pairs_wide.csv", index=False)
        log.info("[CTpre] meta-enrichment written.")


def ctpre_dotplots(base: Path, group2ids: Dict[str,List[str]], label: str="CTpre_all") -> None:
    all_dir = base / f"scenic_{label}"
    out_dir = base / "CTpre_groups"; out_dir.mkdir(parents=True, exist_ok=True)
    iters = iter_dirs(base, label)
    rows = []
    for idir in iters:
        auc_file = next(iter(idir.glob("auc_mtx_*.csv")), None)
        if auc_file is None: 
            continue
        auc = pd.read_csv(auc_file, index_col=0)
        per_iter = None; ok = True
        for g in CT_GROUPS:
            ids = [cid for cid in group2ids.get(g,[]) if cid in auc.index]
            if len(ids)==0: 
                ok=False; break
            mean_s = auc.loc[ids].mean(axis=0).rename(f"MeanAUC_{g}")
            frac_s = (auc.loc[ids] >= DOT_AUC_THR).mean(axis=0).rename(f"FracHigh_{g}")
            df_g = pd.concat([mean_s, frac_s], axis=1).reset_index().rename(columns={"index":"Regulon"})
            per_iter = df_g if per_iter is None else per_iter.merge(df_g, on="Regulon", how="outer")
        if not ok or per_iter is None: 
            continue
        per_iter["Iteration"] = idir.name
        rows.append(per_iter)
    if not rows:
        log.warning("[CTpre] dotplot skipped; no usable AUCell CSVs.")
        return
    df_all = pd.concat(rows, ignore_index=True)
    agg_ops = {f"MeanAUC_{g}":"mean" for g in CT_GROUPS} | {f"FracHigh_{g}":"mean" for g in CT_GROUPS}
    agg = (df_all.groupby("Regulon").agg(agg_ops | {"Iteration":"nunique"}).rename(columns={"Iteration":"n_runs_detected"}).reset_index())
    df = agg[agg["n_runs_detected"] >= DOT_MIN_RUNS_CTP].copy()
    df.to_csv(out_dir / "dotplot_source_table_CTpre_groups.csv", index=False)

    # Regulator ranking preference: meta if available, else variance across groups
    agg_dir = all_dir / "_aggregate"
    ranked = []
    metas = list(agg_dir.glob("meta_enrichment_*_vs_*.csv"))
    if metas:
        frames = []
        for mf in metas:
            try:
                tmp = pd.read_csv(mf, usecols=["Regulon","Z_meta","p_meta","FDR_BH","enriched_in"])
                tmp["pair"] = mf.stem.replace("meta_enrichment_",""); frames.append(tmp)
            except Exception: 
                continue
        if frames:
            meta_all = pd.concat(frames, ignore_index=True)
            meta_all["absZ"] = meta_all["Z_meta"].abs()
            summary = (meta_all.groupby("Regulon", as_index=False)
                       .agg(best_FDR=("FDR_BH","min"), best_absZ=("absZ","max")))
            sig = summary[summary["best_FDR"] <= DOT_META_FDR]
            if not sig.empty:
                ranked = sig.sort_values(["best_FDR","best_absZ"], ascending=[True,False])["Regulon"].tolist()
    if not ranked:
        means = df[[f"MeanAUC_{g}" for g in CT_GROUPS]].fillna(0.0).values
        var = np.var(means, axis=1)
        df["_varAUC"] = var
        ranked = df.sort_values("_varAUC", ascending=False)["Regulon"].tolist()

    for top_n in (30, 15):
        selected = ranked[:top_n]
        sub = df[df["Regulon"].isin(selected)].copy()
        plot_rows = []
        for _, r in sub.iterrows():
            for g in CT_GROUPS:
                plot_rows.append(dict(Regulon=r["Regulon"], Group=g,
                                      MeanAUC=r.get(f"MeanAUC_{g}", np.nan),
                                      FracHigh=r.get(f"FracHigh_{g}", np.nan)))
        plot_df = pd.DataFrame(plot_rows)
        yorder = list(reversed(selected))
        plot_df["Regulon"] = pd.Categorical(plot_df["Regulon"], categories=yorder, ordered=True)
        plot_df["Group"] = pd.Categorical(plot_df["Group"], categories=CT_GROUPS, ordered=True)
        vmin, vmax = float(np.nanmin(plot_df["MeanAUC"])), float(np.nanmax(plot_df["MeanAUC"]))
        if vmin==vmax: 
            vmax=vmin+1e-6
        norm = Normalize(vmin=vmin, vmax=vmax)
        h = max(4.0, 0.28*len(selected)+1.3)
        fig, ax = plt.subplots(figsize=(5.0, h))
        sns.scatterplot(data=plot_df, x="Group", y="Regulon", hue="MeanAUC", palette="viridis",
                        hue_norm=norm, size="FracHigh", sizes=(20,300), linewidth=0, edgecolor="none", ax=ax, legend=False)
        ax.set_ylabel(""); ax.set_xlabel("")
        ax.set_title(f"CTpre groups: mean AUC color, % ≥ {DOT_AUC_THR:g} size", fontsize=9)
        sm = plt.cm.ScalarMappable(cmap=plt.cm.get_cmap("viridis"), norm=norm); sm.set_array([])
        fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.04).set_label("Mean AUCell", fontsize=6)
        levels = _nice_pct_levels(plot_df["FracHigh"], k=3)
        vmin_s, vmax_s = float(np.nanmin(plot_df["FracHigh"])), float(np.nanmax(plot_df["FracHigh"]))
        handles, labels = [], []
        for p in levels:
            area = 20 + ((p/100.0 - vmin_s) / max(vmax_s - vmin_s, 1e-9)) * (300-20)
            hndl = ax.scatter([], [], s=area, facecolor="none", edgecolor="k")
            handles.append(hndl); labels.append(f"{p}%")
        leg = Legend(ax, handles, labels, title="% cells ≥ thr",
                     bbox_to_anchor=(1.03,1.0), loc="upper left", borderaxespad=0., fontsize=6, title_fontsize=6)
        ax.add_artist(leg); fig.tight_layout()
        fig.savefig(out_dir / f"dotplot_CTpre_groups_top{top_n}.png", dpi=300)
        fig.savefig(out_dir / f"dotplot_CTpre_groups_top{top_n}.svg"); plt.close(fig)
    log.info("[CTpre] dotplots saved.")


# ──────────────────────────────────────────────────────────────────────
# D) NetworkX across SCENIC_CAR runs (CD4)
# ──────────────────────────────────────────────────────────────────────
def nx_cd4_across_runs(base: Path, cell_label: str="CAR_CD4", outdir: Optional[Path]=None) -> None:
    outdir = outdir or base / "NetworkX_CD4_vsIterations"
    outdir.mkdir(parents=True, exist_ok=True)

    def scenic_idx(path_str):
        m = re.search(r"SCENIC_CAR_(\d+)$", Path(path_str).name)
        return int(m.group(1)) if m else 10**9

    runs_all = sorted([p for p in (base.parent if base.name=="SCENIC_CAR" else base).glob("SCENIC_CAR_*") if p.is_dir()], key=lambda p: scenic_idx(str(p)))
    if not runs_all:
        log.error("[NX] No SCENIC_CAR_* folders found."); return

    def _extract_edges_from_pickle(path):
        import ctxcore  # noqa: F401
        edges = []
        regs = pickle.load(open(path, "rb"))
        for reg in regs:
            tf = getattr(reg,"name",None)
            g2w = getattr(reg,"gene2weight",{}) or {}
            if tf:
                for tgt,w in g2w.items():
                    if pd.notna(w): edges.append((str(tf), str(tgt), float(w)))
        return edges

    def _extract_edges_from_csv(path):
        df = pd.read_csv(path)
        cmap = {c.lower(): c for c in df.columns}
        TF  = cmap.get("tf", "TF"); TGT = cmap.get("target","target"); WEI = cmap.get("importance","importance")
        out = []
        for _, r in df.iterrows():
            tf, tgt, w = str(r[TF]), str(r[TGT]), r[WEI]
            if pd.notna(w): out.append((tf, tgt, float(w)))
        return out

    def extract_edges(path_str):
        p = Path(path_str)
        if p.suffix.lower() == ".csv":
            return _extract_edges_from_csv(p)
        try:
            return _extract_edges_from_pickle(p)
        except Exception:
            csv_alt = p.with_name(f"adjacencies_{cell_label}.csv")
            if csv_alt.is_file():
                warnings.warn(f"[INFO] Using CSV fallback: {csv_alt.name}")
                return _extract_edges_from_csv(csv_alt)
            raise

    pkl_or_csv = {}
    for r in runs_all:
        pkl = r / cell_label / f"regulons_{cell_label}.pkl"
        csv = r / cell_label / f"adjacencies_{cell_label}.csv"
        if PREFER_ADJ_CSV and csv.is_file():
            pkl_or_csv[r.name] = str(csv)
        elif pkl.is_file():
            pkl_or_csv[r.name] = str(pkl)
        elif csv.is_file():
            pkl_or_csv[r.name] = str(csv)

    if not pkl_or_csv:
        log.error("[NX] No pkl/CSV for runs."); return

    def build_graph_from_edges(edges, wmin, dmin):
        G = nx.DiGraph()
        for u,v,w in edges:
            if w >= wmin: 
                G.add_edge(u,v,weight=float(w))
        if dmin > 0:
            G.remove_nodes_from([n for n,d in G.degree() if d<dmin])
        return G

    def densify_per_run(raw_edges):
        if len(raw_edges)==0: 
            return nx.DiGraph(), nx.DiGraph(), None, None, 0, 0
        by_tf = defaultdict(list)
        for u,v,w in raw_edges: 
            by_tf[u].append((u,v,float(w)))
        edges = []
        for tf,lst in by_tf.items():
            lst.sort(key=lambda x: x[2], reverse=True)
            edges.extend(lst[:TOP_EDGES_PER_TF])
        weights = np.array([w for _,_,w in edges], float)
        uniq = np.unique(weights); uniq.sort()
        chosen = None
        for wmin in uniq[::-1]:
            G = build_graph_from_edges(edges, wmin, dmin=0)
            e = G.number_of_edges()
            if TARGET_EDGE_RANGE[0] <= e <= TARGET_EDGE_RANGE[1]:
                chosen = (float(wmin), G); break
        if chosen is None:
            wmin = float(uniq[0]); chosen = (wmin, build_graph_from_edges(edges, wmin, dmin=0))
        wmin, G = chosen
        if G.number_of_nodes()>0:
            comp = max(nx.weakly_connected_components(G), key=len)
            G_plot = G.subgraph(comp).copy()
        else:
            G_plot = nx.DiGraph()
        return G, G_plot, wmin, 0, len(raw_edges), G.number_of_edges()

    def bc_and_out(G):
        if G.number_of_nodes()==0: 
            return {}, {}
        Gu = G.to_undirected()
        for _,_,d in Gu.edges(data=True):
            d["dist"] = 1.0 / max(d.get("weight", 1e-12), 1e-12)
        bc = nx.betweenness_centrality(Gu, weight="dist")
        wout = dict(G.out_degree(weight="weight"))
        return bc, wout

    def draw_graph(G, title, hubs, out_png):
        if G.number_of_nodes()==0 or G.number_of_edges()==0:
            log.warning(f"[NX] {title}: empty graph; skip plot"); return
        pos = nx.spring_layout(G, k=0.15, seed=SPRING_SEED, weight="weight")
        fig, ax = plt.subplots(figsize=(20,15))
        nx.draw_networkx_edges(G, pos, edge_color="gainsboro", alpha=0.25, ax=ax)
        Gu = G.to_undirected()
        bc = nx.betweenness_centrality(Gu, weight=None)
        vals = np.array(list(bc.values()))
        if vals.size==0 or vals.ptp()==0: sizes = [300]*G.number_of_nodes()
        else: sizes = list(np.interp(vals, (vals.min(), vals.max()), (300, 6000)))
        nx.draw_networkx_nodes(G, pos, node_size=sizes, node_color="steelblue", linewidths=0, ax=ax, alpha=0.9)
        labels = {n:n for n in hubs[:LABEL_TOP]}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=12, font_color="white", ax=ax)
        ax.set_title(title, fontsize=22, fontweight="bold"); ax.set_axis_off(); fig.tight_layout()
        out_png = Path(out_png); fig.savefig(out_png, dpi=300); fig.savefig(out_png.with_suffix(".svg")); plt.close(fig)

    graphs = {}
    node_presence = Counter()
    for tag in sorted(pkl_or_csv, key=lambda k: int(re.search(r"\d+", k).group(0)) if re.search(r"\d+", k) else 0):
        infile = pkl_or_csv[tag]
        raw_edges = extract_edges(infile)
        G_agg, G_plot, wmin_used, _, n_raw, n_kept = densify_per_run(raw_edges)
        if G_agg.number_of_edges()==0:
            log.warning(f"[NX] {tag}: empty after densification."); 
            continue
        node_presence.update([n for n,deg in G_agg.out_degree() if deg>0])
        bc_plot, _ = bc_and_out(G_plot)
        hubs_plot = [n for n,_ in sorted(bc_plot.items(), key=lambda x: x[1], reverse=True)[:LABEL_TOP]]
        draw_graph(G_plot, tag, hubs_plot, outdir / f"{tag}_network.png")
        graphs[tag] = G_agg
        _, wout_agg = bc_and_out(G_agg)
        pd.Series([n for n,_ in sorted(wout_agg.items(), key=lambda x: x[1], reverse=True)[:50]]).to_csv(outdir/f"{tag}_top50_hubs.txt", header=False, index=False)

    if not graphs:
        log.error("[NX] No usable run graphs."); 
        return

    def aggregate(graphs_dict, keep_fn=lambda u,v: True):
        agg = nx.DiGraph()
        for G in graphs_dict.values():
            for u,v,d in G.edges(data=True):
                if not keep_fn(u,v): 
                    continue
                if agg.has_edge(u,v): 
                    agg[u][v]["_w"].append(d["weight"])
                else: 
                    agg.add_edge(u,v,_w=[d["weight"]])
        for u,v,d in agg.edges(data=True):
            d["weight"] = float(np.mean(d["_w"])); del d["_w"]
        return agg

    consensus = aggregate(graphs)
    nx.write_graphml(consensus, outdir / "consensus_CD4.graphml")
    bc_c, _ = bc_and_out(consensus)
    hubs_c = [n for n,_ in sorted(bc_c.items(), key=lambda x: x[1], reverse=True)[:LABEL_TOP]]
    draw_graph(consensus, "CONSENSUS", hubs_c, outdir / "consensus_network.png")

    n_runs = len(graphs)
    min_runs_nodes = int(math.ceil(OVERLAP_FRAC * n_runs))
    occ_nodes = Counter(); [occ_nodes.update(list(G.nodes())) for G in graphs.values()]
    overlap_nodes = {n for n,c in occ_nodes.items() if c >= min_runs_nodes}
    def keep_if_nodes(u,v): return (u in overlap_nodes) and (v in overlap_nodes)
    overlapN = aggregate(graphs, keep_if_nodes)
    if overlapN.number_of_edges()>0:
        bc_n,_ = bc_and_out(overlapN)
        hubs_n = [n for n,_ in sorted(bc_n.items(), key=lambda x: x[1], reverse=True)[:LABEL_TOP]]
        draw_graph(overlapN, f"OVERLAP-NODES (≥{min_runs_nodes}/{n_runs})", hubs_n, outdir / "overlapNodes_network.png")
        nx.write_graphml(overlapN, outdir / f"overlapNodes_p{int(OVERLAP_FR*100)}_CAR_CD4.graphml")

    min_runs_edges = int(math.ceil(EDGE_OVERLAP_FR * n_runs))
    def keep_if_edge(u,v): return sum(G.has_edge(u,v) for G in graphs.values()) >= min_runs_edges
    overlapE = aggregate(graphs, keep_if_edge)
    if overlapE.number_of_edges()>0:
        bc_e,_ = bc_and_out(overlapE)
        hubs_e = [n for n,_ in sorted(bc_e.items(), key=lambda x: x[1], reverse=True)[:LABEL_TOP]]
        draw_graph(overlapE, f"OVERLAP-EDGES (≥{min_runs_edges}/{n_runs})", hubs_e, outdir / "overlapEdges_network.png")
        nx.write_graphml(overlapE, outdir / f"overlapEdges_p{int(EDGE_OVERLAP_FR*100)}_CAR_CD4.graphml")

    min_runs_nodes_strict = int(math.ceil(STRICT_NODE_FR * n_runs))
    occ_nodes_strict = Counter()
    for G in graphs.values(): 
        occ_nodes_strict.update([n for n,deg in G.out_degree() if deg>0])
    strict_nodes = {n for n,c in occ_nodes_strict.items() if c >= min_runs_nodes_strict}
    edge_count = Counter(); [edge_count.update(list(G.edges())) for G in graphs.values()]
    for frac in [1.0, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70]:
        thr = int(math.ceil(frac * n_runs))
        n_edges = sum(1 for (u,v),c in edge_count.items() if c>=thr and (u in strict_nodes and v in strict_nodes))
        if n_edges >= STRICT_AUTO_MIN_EDGES:
            chosen_frac = frac; break
    else:
        chosen_frac = 0.60
    thr_edges = int(math.ceil(chosen_frac * n_runs))
    def keep_if_strict(u,v): return (u in strict_nodes) and (v in strict_nodes) and (edge_count[(u,v)] >= thr_edges)
    overlapS = aggregate(graphs, keep_if_strict)
    if overlapS.number_of_edges()>0:
        bc_s,_ = bc_and_out(overlapS)
        hubs_s = [n for n,_ in sorted(bc_s.items(), key=lambda x: x[1], reverse=True)[:LABEL_TOP]]
        draw_graph(overlapS, f"OVERLAP-STRICT (nodes≥{min_runs_nodes_strict}/{n_runs}, edges≥{thr_edges}/{n_runs})",
                   hubs_s, outdir / "overlapStrict_network.png")
        nx.write_graphml(overlapS, outdir / "overlapStrict_CAR_CD4.graphml")

    dfp = pd.DataFrame(sorted(node_presence.items()), columns=["node","runs_with_outgoing_edge"])
    dfp["fraction"] = dfp["runs_with_outgoing_edge"]/n_runs
    dfp.sort_values(["runs_with_outgoing_edge","node"], ascending=[False,True]).to_csv(outdir/"node_presence_counts.tsv", sep="\t", index=False)
    rows = [(u,v,c, c/n_runs, (u in strict_nodes and v in strict_nodes)) for (u,v),c in edge_count.items()]
    pd.DataFrame(rows, columns=["u","v","presence_count","presence_fraction","both_nodes_strict"]).sort_values(
        ["presence_count","u","v"], ascending=[False,True,True]).to_csv(outdir/"edge_presence_counts.tsv", sep="\t", index=False)
    log.info("[NX] NetworkX outputs written.")


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description="SCENIC full pipeline (one file).")
    ap.add_argument("stage", choices=["all","refcar","llp","ctpre","nx"], help="Which stage to run.")
    # Inputs (stage-dependent; required are checked at runtime)
    ap.add_argument("--adata-ntx", type=Path, default=None, help="Processed NTX AnnData (for 'refcar'/'llp').")
    ap.add_argument("--adata-ctpre", type=Path, default=None, help="Refined Tonsil AnnData (for 'ctpre').")
    ap.add_argument("--gene-db", type=Path, default=None, help="Feather ranking database (pyscenic mc9/mc10).")
    ap.add_argument("--mc10-tbl", type=Path, default=None, help="Motif-to-TF annotation table (.tbl).")
    ap.add_argument("--tf-list", type=Path, default=None, help="TF list (plain text, one per line).")
    # Optional: CAR barcode lists
    ap.add_argument("--car-cd4-file", type=Path, default=None, help="File with CAR CD4 barcodes (one per line).")
    ap.add_argument("--car-cd8-file", type=Path, default=None, help="File with CAR CD8 barcodes (one per line).")
    # Outputs
    ap.add_argument("--car-root", type=Path, default=DEF_CAR_ROOT)
    ap.add_argument("--llp-base", type=Path, default=DEF_LLP_BASE)
    ap.add_argument("--ctpre-base", type=Path, default=DEF_CTPRE_BASE)
    ap.add_argument("--nx-dir", type=Path, default=DEF_NX_DIR)
    # Controls
    ap.add_argument("--iters", type=int, default=None, help="Iterations for 'refcar' (default: 20).")
    ap.add_argument("--dask-workers", type=int, default=None, help="Workers for 'refcar' (default: 8).")
    ap.add_argument("--celltype-col", type=str, default=None, help="Override cell-type column for 'refcar'.")
    ap.add_argument("--force", action="store_true", help="Recompute even if outputs exist.")
    args = ap.parse_args()

    # Resolve optional barcode overrides
    car_cd4_cells = CAR_CD4_BARCODES[:]
    car_cd8_cells = CAR_CD8_BARCODES[:]
    if args.car_cd4_file and Path(args.car_cd4_file).exists():
        car_cd4_cells = _load_lines(args.car_cd4_file)
        log.info(f"[refcar] Using {len(car_cd4_cells)} CAR CD4 barcodes from file.")
    if args.car_cd8_file and Path(args.car_cd8_file).exists():
        car_cd8_cells = _load_lines(args.car_cd8_file)
        log.info(f"[refcar] Using {len(car_cd8_cells)} CAR CD8 barcodes from file.")

    # Stage A: ref->CAR
    if args.stage in ("all","refcar"):
        for what in ("adata_ntx","gene_db","mc10_tbl","tf_list"):
            _require_path(getattr(args, what), what.replace("_","-"), "refcar")
        run_ref_to_car(
            adata_path=args.adata_ntx,
            car_root=args.car_root,
            gene_db=args.gene_db, mc10_tbl=args.mc10_tbl, tf_file=args.tf_list,
            iters=(args.iters or DEF_ITERS_REFCAR),
            seed_base=DEF_SEED_BASE,
            n_workers=(args.dask_workers or DEF_WORKERS_A),
            min_regulon_size=DEF_MIN_REGULON,
            q=DEF_Q_ACTIVITY, min_frac_cd4=DEF_MIN_FRAC_CD4, min_frac_cd8=DEF_MIN_FRAC_CD8,
            force=args.force, celltype_col=args.celltype_col,
            car_cd4_cells=car_cd4_cells, car_cd8_cells=car_cd8_cells
        )

    # Stage B: LLP_all
    if args.stage in ("all","llp"):
        for what in ("adata_ntx","gene_db","mc10_tbl","tf_list"):
            _require_path(getattr(args, what), what.replace("_","-"), "llp")
        adata = ad.read_h5ad(args.adata_ntx)
        obs_names = adata.obs_names.astype(str).tolist()
        pre_ids  = [cid for cid in obs_names if re.match(r"^NTX_[123]_", cid)]
        post_ids = [cid for cid in obs_names if re.match(r"^NTX_[45]_", cid)]
        all_ids  = [cid for cid in obs_names if re.match(r"^NTX_[1-5]_", cid)]
        subset = adata[all_ids].copy() if len(all_ids)>=2 else adata.copy()
        scenic_multi_run(subset, out_base=args.llp_base, label="LLP_all",
                         n_iter=DEF_ITERS_BC, gene_db=args.gene_db, mc10_tbl=args.mc10_tbl, tf_file=args.tf_list,
                         dask_workers=DEF_WORKERS_BC, threads_per_worker=1, force=args.force)
        aggregate_label(args.llp_base, "LLP_all")
        llp_meta_enrichment(args.llp_base, pre_ids, post_ids, label="LLP_all")
        llp_dotplots(args.llp_base, pre_ids, post_ids, label="LLP_all")

    # Stage C: CTpre_all
    if args.stage in ("all","ctpre"):
        for what in ("adata_ctpre","gene_db","mc10_tbl","tf_list"):
            _require_path(getattr(args, what), what.replace("_","-"), "ctpre")
        adata = ad.read_h5ad(args.adata_ctpre)
        group2ids = ctpre_build_groups(adata)
        all_ids = sorted(set().union(*[set(v) for v in group2ids.values()])) if group2ids else []
        if len(all_ids) < 2:
            log.error("[CTpre] <2 total cells in union — aborting CTpre stage.")
        else:
            subset = adata[all_ids].copy()
            scenic_multi_run(subset, out_base=args.ctpre_base, label="CTpre_all",
                             n_iter=DEF_ITERS_BC, gene_db=args.gene_db, mc10_tbl=args.mc10_tbl, tf_file=args.tf_list,
                             dask_workers=DEF_WORKERS_BC, threads_per_worker=1, force=args.force)
            aggregate_label(args.ctpre_base, "CTpre_all")
            ctpre_meta_pairwise(args.ctpre_base, group2ids, label="CTpre_all")
            ctpre_dotplots(args.ctpre_base, group2ids, label="CTpre_all")

    # Stage D: NetworkX CD4 (across SCENIC_CAR runs)
    if args.stage in ("all","nx"):
        args.nx_dir.mkdir(parents=True, exist_ok=True)
        nx_cd4_across_runs(base=args.car_root, cell_label="CAR_CD4", outdir=args.nx_dir)


if __name__ == "__main__":
    main()
