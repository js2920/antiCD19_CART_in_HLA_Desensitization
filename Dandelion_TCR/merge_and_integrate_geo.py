#!/usr/bin/env python3
# ========================================================================
#  TCR integration – Dandelion contigs ➜ AnnData.obs  (path‑agnostic)
# ========================================================================
"""
dandelion_integration_TCR.py
============================

• Merges productive TRA/TRB contigs (Dandelion) from multiple samples
• Ensures single‑pair chains per cell, drops ambiguous cells
• Writes per‑cell CSV + augmented AnnData

Run example
-----------
python scripts/dandelion_integration_TCR.py \
    --contigs  /PATH/TO/NTX_2_TCR_fastq/dandelion/all_contig_dandelion.tsv \
               /PATH/TO/NTX_3_TCR_fastq/dandelion/all_contig_dandelion.tsv \
               /PATH/TO/NTX_4_TCR_fastq/dandelion/all_contig_dandelion.tsv \
               /PATH/TO/NTX_5_TCR_fastq/dandelion/all_contig_dandelion.tsv \
    --adata    results/totalvi_celltypist/processed_adata_NTX_all_samples.h5ad \
    --out-dir  results/tcr_integration
"""
from __future__ import annotations

import argparse
import logging
import re
import sys
from pathlib import Path
from typing import List, Set

import anndata as ad
import pandas as pd

# ──────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────
USECOLS = [
    "cell_id",
    "barcode",
    "locus",
    "chain",
    "productive",
    "umi_count",
    "consensus_count",
    "read_count",
    "v_call",
    "j_call",
    "junction_aa",
    "cdr3_aa",
    "cdr3",
]
# Default prefix detector matches 'NTX_<digit>' anywhere in the path
_PREFIX_RE = re.compile(r"(NTX_\d+)")


def detect_prefix(path: Path) -> str:
    """
    Extracts a sample prefix (e.g. 'NTX_3') from the file path.
    Adjust _PREFIX_RE if your naming differs.
    """
    m = _PREFIX_RE.search(str(path))
    if not m:
        raise ValueError(f"Could not detect an ‘NTX_#’ prefix in {path}")
    return m.group(1)


def unify_barcode(raw: str, prefix: str) -> str:
    """
    Harmonise a raw Cell Ranger barcode to '<prefix>_<barcode>-1'.
    Removes lane suffixes (e.g. '-2') and ensures single '-1'.
    """
    raw = str(raw)
    raw = re.sub(r"-\d+$", "", raw)  # strip lane suffix like '-2'
    if raw.endswith("-1"):
        raw = raw[:-2]
    return f"{prefix}_{raw}-1"


def read_contig_file(fp: Path, logger: logging.Logger) -> pd.DataFrame:
    prefix = detect_prefix(fp)
    logger.info("Reading %-32s", fp.name)

    # Read only the columns available in the file
    with fp.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
    cols = [c for c in USECOLS if c in header]
    df = pd.read_csv(fp, sep="\t", usecols=cols, low_memory=True)

    # Harmonise column names
    df.rename(
        columns={
            "cell_id": "barcode",
            "locus": "chain",
            "junction_aa": "cdr3",
            "cdr3_aa": "cdr3",
        },
        inplace=True,
    )

    # Barcode fix
    df["barcode"] = df["barcode"].map(lambda b: unify_barcode(b, prefix))

    # Productive → bool
    if "productive" in df.columns:
        df["productive"] = (
            df["productive"]
            .astype(str)
            .str.upper()
            .map({"TRUE": True, "T": True, "1": True, "FALSE": False, "F": False, "0": False})
            .fillna(False)
        )
    else:
        df["productive"] = True

    # Keep only productive TRA/TRB
    df["chain"] = df["chain"].astype(str).str.upper().str.extract(r"(TRA|TRB)", expand=False)
    df = df[df["chain"].isin(["TRA", "TRB"]) & df["productive"]]

    logger.info("   %6d productive TRA/TRB contigs", len(df))
    return df.reset_index(drop=True)


def load_all_contigs(paths: List[Path], logger: logging.Logger) -> pd.DataFrame:
    return pd.concat([read_contig_file(p, logger) for p in paths], ignore_index=True)


def find_multiplets(df: pd.DataFrame) -> Set[str]:
    """Cells with > 1 productive chain of a given type are removed."""
    multi: Set[str] = set()
    for ch in ("TRA", "TRB"):
        bad = df[df["chain"] == ch]["barcode"].value_counts()
        multi.update(bad[bad > 1].index)
    return multi


def deduplicate(df: pd.DataFrame) -> pd.DataFrame:
    """
    Keep a single chain per (barcode, chain) using a priority metric:
    umi_count > consensus_count > read_count > constant 1
    """
    metric = next((m for m in ("umi_count", "consensus_count", "read_count") if m in df.columns), None)
    if metric is None:
        df["__metric"] = 1
        metric = "__metric"

    df[metric] = pd.to_numeric(df[metric], errors="coerce").fillna(0)
    idx = (
        df.sort_values(metric, ascending=False)
        .groupby(["barcode", "chain"], sort=False)
        .head(1)
        .index
    )
    return df.loc[idx].copy().reset_index(drop=True)


def pivot_per_cell(df: pd.DataFrame) -> pd.DataFrame:
    """Wide per‑cell table: TRA_cdr3, TRB_v_call, …"""
    keep = ["cdr3", "v_call", "j_call"]
    df = df[["barcode", "chain", *keep]].set_index("barcode")

    def _pivot(chain: str) -> pd.DataFrame:
        sub = df[df["chain"] == chain]
        if sub.empty:
            return pd.DataFrame(columns=[f"{chain}_{c}" for c in keep], dtype=str)
        sub = sub[keep].astype(str)
        sub.columns = [f"{chain}_{c}" for c in sub.columns]
        return sub

    return _pivot("TRA").join(_pivot("TRB"), how="outer")


# ──────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Integrate productive Dandelion TCR contigs into AnnData.obs",
    )
    p.add_argument("--contigs", type=Path, nargs="+", required=True, help="all_contig_dandelion.tsv files")
    p.add_argument("--adata", type=Path, required=True, help="Input AnnData (h5ad)")
    p.add_argument("--out-dir", type=Path, required=True, help="Output directory (will be created)")
    p.add_argument("--log-level", choices=["DEBUG", "INFO", "WARNING"], default="INFO")
    return p


def main() -> None:
    args = build_parser().parse_args()

    # Logging
    out_dir: Path = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    log_file = out_dir / "tcr_integration.log"

    logging.basicConfig(
        level=args.log_level,
        format="%(asctime)s │ %(levelname)-8s │ %(message)s",
        datefmt="%H:%M:%S",
        handlers=[logging.FileHandler(log_file, mode="w"), logging.StreamHandler(sys.stdout)],
    )
    log = logging.getLogger(__name__)
    log.info("▶ TCR integration started")

    # Load & filter contigs
    contigs = load_all_contigs(args.contigs, log)
    multiplets = find_multiplets(contigs)
    log.info("Dropping %d multiplet cells", len(multiplets))
    contigs = contigs[~contigs["barcode"].isin(multiplets)]
    contigs = deduplicate(contigs)

    per_cell = pivot_per_cell(contigs)
    per_cell_csv = out_dir / "per_cell_TCR_table.csv"
    per_cell.to_csv(per_cell_csv)
    log.info("Per‑cell table  → %s", per_cell_csv)

    # AnnData augmentation
    log.info("Reading AnnData  → %s", args.adata)
    adata = ad.read_h5ad(args.adata)
    if multiplets:
        adata = adata[~adata.obs_names.isin(multiplets)].copy()
        log.info("AnnData after multiplet removal: %d cells", adata.n_obs)

    adata.obs = adata.obs.join(per_cell, how="left")
    adata.obs["has_TCR"] = adata.obs["TRA_cdr3"].notna() | adata.obs["TRB_cdr3"].notna()
    adata.obs = adata.obs.applymap(lambda x: str(x) if pd.notna(x) else "")

    out_h5ad = out_dir / "adata_with_TCR_integration.h5ad"
    adata.write(out_h5ad)
    log.info("Augmented AnnData → %s", out_h5ad)
    log.info("✓ DONE")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:  # noqa: BLE001
        logging.getLogger(__name__).exception("Script failed: %s", exc)
        sys.exit(1)
