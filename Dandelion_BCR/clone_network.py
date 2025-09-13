#!/usr/bin/env python3
"""
combined_clonotype_networks.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Generate clone‑network maps for

• every T‑cell subset
• all T‑cells combined
• selected B‑cell / plasma‑cell populations

Each map is stored as <name>.png AND <name>.svg in the chosen output directory.

Parameters
----------
Clone definition   : junction_aa, 85% AA identity (overridable)
Visual filter      : clones with ≥3 cells (overridable)
Colours            : crimson (T lineage) • royalblue (B / plasma lineage)

Configuration
-------------
You can pass inputs via CLI or environment variables:

  CLI (preferred):
    --workdir, --outdir
    --tcr-contigs, --tcr-adata
    --bcr-contigs-tsv, --bcr-adata
    --bcr-obs-column, --pop-...-value (for B-cell labels)
    --pre-ids, --post-ids
    --thresh-id, --min-clone, --node-scale, --layout
    --skip-tcr, --skip-bcr

  Environment (fallbacks):
    CLONENET_WORKDIR, CLONENET_OUTDIR

Defaults mirror the original layout under <workdir>, but you can override all paths.

Circular rendering is robust to NaNs/Infs; it falls back safely and avoids invalid axes.
"""

from __future__ import annotations

import argparse
import logging
import os
import pathlib
import re
from typing import List, Dict, Optional, Set

import numpy as np
import anndata as ad
import pandas as pd
import matplotlib.pyplot as plt
import dandelion as ddl
from dandelion.utilities import read_airr

# ─── Globals populated by configure_paths() ─────────────────────────────────
WORK: pathlib.Path
OUTDIR: pathlib.Path

TCR_CONTIGS: List[pathlib.Path]
TCR_ADATA: pathlib.Path

BCR_CONTIGS_TSV: pathlib.Path
BCR_ADATA: pathlib.Path

# PRE/POST, thresholds & style (overridable via CLI)
PRE_IDS: Set[str] = {"NTX_1", "NTX_2", "NTX_3"}
POST_IDS: Set[str] = {"NTX_4", "NTX_5"}

THRESH_ID: float = 0.85
MIN_CLONE: int = 3
NODE_SCALE: int = 5
LAYOUT: str = "mod_fr"

COLOR_T = "crimson"
COLOR_B = "royalblue"

# ─── Logging ────────────────────────────────────────────────────────────────
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-7s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ─── CLI / environment configuration ────────────────────────────────────────
def _resolve_path(p: str) -> pathlib.Path:
    return pathlib.Path(os.path.expanduser(p)).resolve()

def _pick_dir(cli_value: Optional[str], env_key: str, default: str) -> pathlib.Path:
    if cli_value:
        return _resolve_path(cli_value)
    env_val = os.environ.get(env_key)
    if env_val:
        return _resolve_path(env_val)
    return _resolve_path(default)

def _split_ids(s: Optional[str]) -> Set[str]:
    if not s:
        return set()
    return {x.strip() for x in s.split(",") if x.strip()}

def configure_paths(args: argparse.Namespace) -> None:
    """Populate module‑level path variables according to CLI/env/defaults."""
    global WORK, OUTDIR, TCR_CONTIGS, TCR_ADATA, BCR_CONTIGS_TSV, BCR_ADATA
    global PRE_IDS, POST_IDS, THRESH_ID, MIN_CLONE, NODE_SCALE, LAYOUT

    WORK   = _pick_dir(args.workdir, "CLONENET_WORKDIR", "/work")
    OUTDIR = _pick_dir(args.outdir,  "CLONENET_OUTDIR",  str(WORK / "Clonotype_Networks"))
    OUTDIR.mkdir(parents=True, exist_ok=True)

    # PRE/POST sets
    if args.pre_ids:
        PRE_IDS  = _split_ids(args.pre_ids)
    if args.post_ids:
        POST_IDS = _split_ids(args.post_ids)

    # Thresholds & style
    if args.thresh_id is not None:
        THRESH_ID = float(args.thresh_id)
    if args.min_clone is not None:
        MIN_CLONE = int(args.min_clone)
    if args.node_scale is not None:
        NODE_SCALE = int(args.node_scale)
    if args.layout is not None:
        LAYOUT = str(args.layout)

    # Inputs (either explicit via CLI or defaults relative to WORK)
    if args.tcr_contigs:
        TCR_CONTIGS = [_resolve_path(p) for p in args.tcr_contigs]
    else:
        TCR_CONTIGS = [
            WORK / f"NTX_{n}_TCR_fastq/dandelion/all_contig_dandelion.tsv"
            for n in (2, 3, 4, 5)
        ]
    TCR_ADATA = _resolve_path(args.tcr_adata) if args.tcr_adata else (WORK / "TCR/adata_with_TCR_integration.h5ad")

    BCR_CONTIGS_TSV = _resolve_path(args.bcr_contigs_tsv) if args.bcr_contigs_tsv else (WORK / "my_merged_contigs.tsv")
    BCR_ADATA       = _resolve_path(args.bcr_adata)       if args.bcr_adata       else (WORK / "adata_with_bcr_integration.h5ad")

    log.info("Configured WORK    : %s", WORK)
    log.info("Configured OUTDIR  : %s", OUTDIR)
    log.info("TCR contigs (n=%d) : %s", len(TCR_CONTIGS), ", ".join(map(str, TCR_CONTIGS)))
    log.info("TCR AnnData        : %s", TCR_ADATA)
    log.info("BCR contigs TSV    : %s", BCR_CONTIGS_TSV)
    log.info("BCR AnnData        : %s", BCR_ADATA)

# ─── Helpers ────────────────────────────────────────────────────────────────
PREFIX_RE = re.compile(r"(NTX_\d+)")
SAFE_RE   = re.compile(r"[^A-Za-z0-9]+")

def detect_prefix(p: pathlib.Path) -> str:
    m = PREFIX_RE.search(str(p))
    return m.group(1) if m else "UNK"

def safe_name(txt: str) -> str:
    return SAFE_RE.sub("_", txt).strip("_")

def load_vdj(paths: List[pathlib.Path]) -> ddl.Dandelion:
    existing = [p for p in paths if p.exists()]
    missing  = [p for p in paths if not p.exists()]
    for p in missing:
        log.warning("Missing contig file: %s", p)
    if not existing:
        raise FileNotFoundError("No VDJ contig files found among: " + ", ".join(map(str, paths)))
    parts = [read_airr(str(p), prefix=detect_prefix(p)) for p in existing]
    for p, part in zip(existing, parts):
        log.info("Loaded %-40s → %5d contigs", p.name, part.n_contigs)
    return ddl.concat(parts)

def make_clones(vdj: ddl.Dandelion) -> None:
    vdj.threshold = THRESH_ID
    ddl.tl.find_clones(vdj, key="junction_aa")
    # ensure clone_size available
    if "clone_id_size" not in vdj.metadata.columns:
        ddl.tl.clone_size(vdj)
    vdj.metadata.rename(columns={"clone_id_size": "clone_size"}, inplace=True)

def add_network(vdj: ddl.Dandelion) -> None:
    ddl.tl.generate_network(
        vdj,
        clone_col="clone_id",
        min_size=MIN_CLONE,
        layout_method=LAYOUT,
    )

def add_clone_size_obs(adata: ad.AnnData) -> None:
    """Ensure adata.obs['clone_size'] exists (counts per clone_id)."""
    if "clone_size" not in adata.obs.columns and "clone_id" in adata.obs.columns:
        vc = adata.obs["clone_id"].value_counts(dropna=False)
        adata.obs["clone_size"] = adata.obs["clone_id"].map(vc).fillna(0).astype(int)

def size_legend(ax, color):
    for lbl, n in (("1 cell", 1), ("5 cells", 5), ("≥10 cells", 10)):
        ax.scatter([], [], s=n*NODE_SCALE, c=color, edgecolors="none", label=lbl)
    ax.legend(title="Clone size", frameon=False, scatterpoints=1,
              fontsize="small", loc="center left", bbox_to_anchor=(1.02, 0.5))

# ─── helpers for circular axes ─────────────────────────────────────────────
def _collect_node_points(ax) -> np.ndarray:
    """
    Gather all node (x,y) positions from PathCollections on the axes,
    dropping masked and non‑finite rows. Returns (N,2) array or empty array.
    """
    chunks = []
    for coll in ax.collections:
        if not hasattr(coll, "get_offsets"):
            continue
        offs = coll.get_offsets()
        if offs is None:
            continue
        arr = np.asarray(offs)
        if arr.ndim != 2 or arr.shape[1] != 2 or arr.size == 0:
            continue
        if np.ma.isMaskedArray(arr):
            mask = np.ma.getmaskarray(arr)
            arr = arr[~mask.any(axis=1)]
            arr = np.asarray(arr.filled(np.nan))
        finite = np.isfinite(arr).all(axis=1)
        arr = arr[finite]
        if arr.size:
            chunks.append(arr)
    if chunks:
        return np.vstack(chunks)
    return np.empty((0, 2), dtype=float)

def enforce_circular_axes(ax, pad: float = 0.05) -> None:
    """
    Force the network to appear circular by:
      • setting equal data aspect
      • applying symmetric x/y limits about the data center
    Robust to missing/non‑finite coordinates.
    """
    ax.relim(); ax.autoscale_view()
    x0, x1 = ax.get_xlim(); y0, y1 = ax.get_ylim()
    finite_limits = np.isfinite([x0, x1, y0, y1]).all()
    span_x = (x1 - x0) if finite_limits else np.nan
    span_y = (y1 - y0) if finite_limits else np.nan

    if not finite_limits or span_x <= 0 or span_y <= 0:
        pts = _collect_node_points(ax)
        if pts.size >= 2:
            x0, x1 = float(np.min(pts[:, 0])), float(np.max(pts[:, 0]))
            y0, y1 = float(np.min(pts[:, 1])), float(np.max(pts[:, 1]))
            finite_limits = np.isfinite([x0, x1, y0, y1]).all()
            span_x = (x1 - x0) if finite_limits else np.nan
            span_y = (y1 - y0) if finite_limits else np.nan
        else:
            x0, x1, y0, y1 = 0.0, 1.0, 0.0, 1.0
            span_x, span_y = 1.0, 1.0

    cx = (x0 + x1) / 2.0
    cy = (y0 + y1) / 2.0
    r  = max(span_x, span_y) / 2.0
    if not np.isfinite(r) or r <= 0:
        r = 0.5
    r *= (1.0 + pad)

    try:
        ax.set_box_aspect(1)  # mpl >= 3.3
    except Exception:
        pass
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(cx - r, cx + r)
    ax.set_ylim(cy - r, cy + r)

def plot_network(adata_sub: ad.AnnData, color: str,
                 title: str, stem: pathlib.Path):
    stem.parent.mkdir(parents=True, exist_ok=True)
    if getattr(adata_sub, "n_obs", 0) == 0:
        log.info("No cells to plot for %s; skipping.", title)
        return

    add_clone_size_obs(adata_sub)
    sizes = (adata_sub.obs.get("clone_size", 1)
             .fillna(1).astype(float).values * NODE_SCALE)

    fig, ax = plt.subplots(figsize=(6, 6), dpi=300, facecolor="white")
    ddl.pl.clone_network(
        adata_sub,
        edges=False,
        size=sizes,
        show=False,
        ax=ax,
    )
    if ax.collections:
        try:
            ax.collections[0].set_facecolor(color)
            ax.collections[0].set_edgecolor("none")
        except Exception:
            pass

    try:
        enforce_circular_axes(ax, pad=0.05)
    except Exception as e:
        log.warning("Circular axis adjustment skipped (%s); using autoscale.", e)
        ax.relim(); ax.autoscale_view()
        try:
            ax.set_aspect("equal", adjustable="box")
            ax.set_box_aspect(1)
        except Exception:
            pass

    size_legend(ax, color)
    ax.set_title(title)
    ax.axis("off")
    plt.tight_layout(rect=[0, 0, 0.82, 1.0])

    fig.savefig(stem.with_suffix(".png"))
    fig.savefig(stem.with_suffix(".svg"))
    plt.close(fig)
    try:
        rel = stem.relative_to(OUTDIR)
        log.info("Wrote %s (.png + .svg)", rel)
    except Exception:
        log.info("Wrote %s (.png + .svg)", stem)

# ─── barcode helper for BCR contigs ───────────────────────────────────────
def rename_barcode(bc: str) -> str:
    return bc.replace("_fastq_", "_") + "-1" if "_fastq_" in bc else bc

def coerce_numeric(df: pd.DataFrame):
    for col in ("UMI_count", "UMIs", "umis", "umi_count",
                "consensus_count", "read_count", "reads"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0).astype(int)

# ═══════════════════════════════════════════════════════════════════════════
# 1)  T‑cell workflow
# ═══════════════════════════════════════════════════════════════════════════
def run_tcr():
    log.info("[TCR] loading VDJ contigs …")
    vdj = load_vdj(TCR_CONTIGS)
    vdj = vdj[vdj.metadata.chain_status == "Single pair"].copy()

    log.info("[TCR] clone calling …")
    make_clones(vdj)
    add_network(vdj)

    log.info("[TCR] loading AnnData …")
    adata = ad.read_h5ad(TCR_ADATA)
    ddl.tl.transfer(adata, vdj)

    # pick T subsets (incl. MAIT Cells)
    subsets = sorted(
        lab for lab in adata.obs["cell_type"].unique()
        if lab.lower() == "mait cells" or "t cells" in lab.lower()
    )
    if not subsets:
        log.warning("[TCR] no T‑cell subsets found")
        return

    # per‑subset maps
    for lab in subsets:
        mask = adata.obs["cell_type"] == lab
        for phase, ids in (("Pre", PRE_IDS), ("Post", POST_IDS)):
            sel = mask & adata.obs["sample_id"].isin(ids)
            if not sel.any():
                continue
            stem = OUTDIR / f"Tcell_{safe_name(lab)}_{phase}"
            plot_network(adata[sel], COLOR_T, f"{lab} ({phase})", stem)

    # ALL T‑cells maps
    tmask = adata.obs["cell_type"].isin(subsets)

    if tmask.any():
        stem = OUTDIR / "Tcell_All"
        plot_network(adata[tmask], COLOR_T, "All T‑cells (combined)", stem)

    for phase, ids in (("Pre", PRE_IDS), ("Post", POST_IDS)):
        sel = tmask & adata.obs["sample_id"].isin(ids)
        if sel.any():
            stem = OUTDIR / f"Tcell_All_{phase}"
            plot_network(adata[sel], COLOR_T, f"All T‑cells ({phase})", stem)

# ═══════════════════════════════════════════════════════════════════════════
# 2)  B‑cell / plasma‑cell workflow
# ═══════════════════════════════════════════════════════════════════════════
def find_bcr_label_column(adata: ad.AnnData, prefer: Optional[str]) -> Optional[str]:
    candidates = [prefer] if prefer else []
    candidates += [
        "celltypist_label",
        "Tonsil_CellTypist_Label",
        "Human_Tonsil_Classifier_Reclassification",
        "majority_voting",
        "predicted_labels",
        "cell_type",
    ]
    for c in candidates:
        if c and c in adata.obs.columns:
            return c
    return None

def run_bcr(bcr_obs_col: Optional[str],
            pop_mbc_igg_pc: str,
            pop_mature_igg_pc: str,
            pop_csmbc: str):
    log.info("[BCR] reading contig table …")
    df = pd.read_csv(BCR_CONTIGS_TSV, sep="\t", dtype=str)
    if "cell_id" in df.columns:
        df["cell_id"] = df["cell_id"].map(rename_barcode)
    coerce_numeric(df)

    vdj = ddl.Dandelion(df)
    log.info("[BCR] %d contigs loaded", vdj.n_contigs)

    adata = ad.read_h5ad(BCR_ADATA)
    vdj, adata = ddl.pp.check_contigs(vdj, adata, productive_only=True)
    vdj.simplify()

    make_clones(vdj)
    add_network(vdj)
    ddl.tl.transfer(adata, vdj, overwrite=True)

    adata.obs["CART_phase"] = adata.obs["sample_id"].map(
        lambda x: "Pre" if x in PRE_IDS else "Post" if x in POST_IDS else "Other"
    )

    label_col = find_bcr_label_column(adata, bcr_obs_col)
    if not label_col:
        log.warning("[BCR] No suitable label column found; available obs: %s", ", ".join(adata.obs.columns))
        return
    log.info("[BCR] Using label column: %s", label_col)

    populations: Dict[str, Dict[str, str]] = {
        "MBC_IgG_PC":    {"col": label_col, "value": pop_mbc_igg_pc},
        "Mature_IgG_PC": {"col": label_col, "value": pop_mature_igg_pc},
        "csMBC":         {"col": label_col, "value": pop_csmbc},
    }

    for pop, spec in populations.items():
        m = (adata.obs[spec["col"]].astype(str) == spec["value"])
        if not m.any():
            log.info("[BCR] No cells for population '%s' (= '%s').", pop, spec["value"])
            continue
        for phase in ("Pre", "Post"):
            if pop == "csMBC" and phase == "Post":
                # preserve the original skip rule
                continue
            sel = m & (adata.obs["CART_phase"] == phase)
            if not sel.any():
                continue
            stem = OUTDIR / f"Bcell_{safe_name(pop)}_{phase}"
            plot_network(adata[sel], COLOR_B, f"{pop} ({phase})", stem)

# ═══════════════════════════════════════════════════════════════════════════
# 3)  CLI & main
# ═══════════════════════════════════════════════════════════════════════════
def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Clonotype network figure generator (T & B lineages).")
    # Directories
    p.add_argument("--workdir", type=str, default=None,
                   help="Base working directory for inputs (default: /work or $CLONENET_WORKDIR).")
    p.add_argument("--outdir", type=str, default=None,
                   help="Directory for outputs (default: <workdir>/Clonotype_Networks or $CLONENET_OUTDIR).")
    # Inputs (optional overrides)
    p.add_argument("--tcr-contigs", type=str, nargs="+", default=None,
                   help="Paths to TCR all_contig_dandelion.tsv files.")
    p.add_argument("--tcr-adata", type=str, default=None,
                   help="AnnData with TCR integration (default: <workdir>/TCR/adata_with_TCR_integration.h5ad).")
    p.add_argument("--bcr-contigs-tsv", type=str, default=None,
                   help="Merged BCR contigs TSV (default: <workdir>/my_merged_contigs.tsv).")
    p.add_argument("--bcr-adata", type=str, default=None,
                   help="AnnData with BCR + labels (default: <workdir>/adata_with_bcr_integration.h5ad).")
    # Labels & populations
    p.add_argument("--bcr-obs-column", type=str, default=None,
                   help="obs column to use for B-cell populations (auto-detects if omitted).")
    p.add_argument("--pop-mbc-igg-pc-value", type=str, default="MBC derived IgG+ PC",
                   help="Exact label value for MBC-derived IgG+ PCs.")
    p.add_argument("--pop-mature-igg-pc-value", type=str, default="Mature IgG+ PC",
                   help="Exact label value for Mature IgG+ PCs.")
    p.add_argument("--pop-csmbc-value", type=str, default="csMBC",
                   help="Exact label value for class-switched MBCs.")
    # Cohort sets
    p.add_argument("--pre-ids", type=str, default=None,
                   help="Comma-separated PRE sample IDs (default: NTX_1,NTX_2,NTX_3).")
    p.add_argument("--post-ids", type=str, default=None,
                   help="Comma-separated POST sample IDs (default: NTX_4,NTX_5).")
    # Thresholds & rendering
    p.add_argument("--thresh-id", type=float, default=None, help="AA identity for clones (default: 0.85).")
    p.add_argument("--min-clone", type=int, default=None, help="Minimum clone size to plot (default: 3).")
    p.add_argument("--node-scale", type=int, default=None, help="Scalar for node sizes (default: 5).")
    p.add_argument("--layout", type=str, default=None, help="Dandelion layout method (default: mod_fr).")
    # Stage toggles
    p.add_argument("--skip-tcr", action="store_true", help="Skip T-cell networks.")
    p.add_argument("--skip-bcr", action="store_true", help="Skip B-cell/plasma networks.")
    return p.parse_args()

def main():
    args = parse_args()
    configure_paths(args)

    if not args.skip_tcr:
        run_tcr()
    else:
        log.info("Skipping TCR stage (per --skip-tcr).")

    if not args.skip_bcr:
        run_bcr(
            bcr_obs_col=args.bcr_obs_column,
            pop_mbc_igg_pc=args.pop_mbc_igg_pc_value,
            pop_mature_igg_pc=args.pop_mature_igg_pc_value,
            pop_csmbc=args.pop_csmbc_value,
        )
    else:
        log.info("Skipping BCR stage (per --skip-bcr).")

    log.info("Clone‑network figures written to %s", OUTDIR)

if __name__ == "__main__":
    main()
