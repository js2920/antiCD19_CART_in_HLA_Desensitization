#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
eplet_similarity.py 

Given just an antibodies CSV, this script will:
  1) Detect HLA antigens in the file.
  2) If needed, synthesize a minimal MFI CSV and call:
       Rscript scripts/buildepletmap_A.R  -> builds eplet_map.csv + Registry snapshot
  3) Build antigen×eplet incidence (Python).
  4) Compute two similarity panels + clustered heatmaps:
       - IDF-weighted Jaccard
       - Szymkiewicz–Simpson (overlap coefficient)

Typical use:
  python scripts/eplet_similarity_for_antibodies_A.py \
    --antibodies antibodies_python.csv \
    --outdir outputs/similarity_antibodies

Optional:
  --mfi /path/to/MFI.csv       (if you already have an MFI-style table with antigen columns)
  --eplet-ref-dir outputs/eplet_ref
  --rscript Rscript            (path to Rscript; default: 'Rscript' on PATH)
  --r-build-script scripts/buildepletmap_A.R
  --verified-only              (passed to the R builder)

Outputs in --outdir:
  EpletIDFWeighted_Jaccard_matrix_antibodies.csv (+ distance, + .png/.svg heatmap)
  EpletOverlap_SzymkiewiczSimpson_matrix_antibodies.csv (+ distance, + .png/.svg heatmap)
  similarity_MANIFEST.txt
"""
from __future__ import annotations
import argparse, re, sys, subprocess, shlex, tempfile
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
from scipy.spatial.distance import squareform

# ----------------- tokenization -----------------
ALLELE_PAT = re.compile(r"\b(?:HLA-)?[A-Z0-9]{1,5}\*\d+(?::\d+){1,3}[A-Z]*\b", re.IGNORECASE)

def normalize_antigen(s: str) -> str:
    s = (s or "").replace("\u00A0", " ").strip()
    s = re.sub(r"\s*\([^)]*\)\s*$", "", s)
    s = s.replace("HLA-", "")
    return s.upper()

def extract_hla_alleles(text: str) -> List[str]:
    toks = ALLELE_PAT.findall((text or "").replace("\u00A0", " "))
    seen, out = set(), []
    for t in toks:
        tt = normalize_antigen(t)
        if tt not in seen:
            seen.add(tt); out.append(tt)
    return out

def detect_antigens_from_df(df: pd.DataFrame) -> List[str]:
    # (1) common single-column headers
    for pat in (r"^antigen$", r"^allele$", r"^hla$", r"^antibody$", r"^specificity$"):
        cols = [c for c in df.columns if re.search(pat, str(c), re.IGNORECASE)]
        if cols:
            vals = pd.Series(df[cols[0]]).dropna().astype(str).tolist()
            out = []
            for v in vals:
                out.extend(extract_hla_alleles(v) or [normalize_antigen(v)])
            return sorted(set(a for a in out if "*" in a))
    # (2) wide tables whose column names are antigens like A*02:01
    wide_cols = [c for c in df.columns if ("*" in str(c)) and not re.search(r"^(?i)(pc_|nc_)", str(c))]
    if wide_cols:
        return sorted(set(normalize_antigen(c) for c in wide_cols))
    # (3) scan cells
    out: List[str] = []
    for col in df.columns:
        for v in pd.Series(df[col]).astype(str).tolist():
            out.extend(extract_hla_alleles(v))
    return sorted(set(out))

# ----------------- eplet similarity -----------------
def idf_weights(inc: pd.DataFrame) -> pd.Series:
    N = float(inc.shape[0])
    df = inc.sum(axis=0).astype(float).clip(lower=1.0)
    w = np.log(N / df).clip(lower=0.0)
    w.name = "idf_weight"
    return w

def weighted_jaccard_matrix(inc: pd.DataFrame, weights: pd.Series) -> pd.DataFrame:
    W = weights.reindex(inc.columns).fillna(0.0).astype(float)
    A = inc.values.astype(float)
    X = A * W.values
    inter = A.dot(X.T)
    sums  = X.sum(axis=1, keepdims=True)
    union = sums + sums.T - inter
    with np.errstate(divide="ignore", invalid="ignore"):
        S = np.where(union > 0, inter / union, 1.0)
    return pd.DataFrame(S, index=inc.index, columns=inc.index)

def overlap_coefficient_matrix(inc: pd.DataFrame) -> pd.DataFrame:
    A = inc.values.astype(int)
    sizes = A.sum(axis=1, keepdims=True)
    inter = A.dot(A.T)
    mins  = np.minimum(sizes, sizes.T)
    with np.errstate(divide="ignore", invalid="ignore"):
        S = np.where(mins > 0, inter / mins, 1.0)
    return pd.DataFrame(S, index=inc.index, columns=inc.index)

def save_heatmap_with_dendrogram(names: List[str], S: np.ndarray, out_png: Path, title: str) -> None:
    from scipy.cluster.hierarchy import linkage, leaves_list, dendrogram
    from scipy.spatial.distance import squareform

    names = list(names)
    S = np.asarray(S, float)
    S = np.clip(S, 0.0, 1.0)
    D = 1.0 - S

    Z = linkage(squareform(D, checks=False), method="average", optimal_ordering=True)
    order = leaves_list(Z)
    names_ord = [names[i] for i in order]
    S_ord = S[np.ix_(order, order)]

    fig = plt.figure(figsize=(min(20, 4 + 0.13*len(names)), min(20, 4 + 0.13*len(names))), dpi=220)
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 6], height_ratios=[1, 6], wspace=0.02, hspace=0.02)
    ax_top  = fig.add_subplot(gs[0, 1])
    ax_left = fig.add_subplot(gs[1, 0])
    ax_main = fig.add_subplot(gs[1, 1])

    dendrogram(Z, no_labels=True, color_threshold=0, ax=ax_top);  ax_top.axis("off")
    dendrogram(Z, orientation="right", no_labels=True, color_threshold=0, ax=ax_left); ax_left.axis("off")

    im = ax_main.imshow(S_ord, vmin=0.0, vmax=1.0, cmap="viridis", aspect="auto")
    ax_main.set_xticks(range(len(names_ord))); ax_main.set_yticks(range(len(names_ord)))
    ax_main.set_xticklabels(names_ord, rotation=90, fontsize=7)
    ax_main.set_yticklabels(names_ord, fontsize=7)

    cb = fig.colorbar(im, ax=ax_main, fraction=0.046, pad=0.02)
    cb.set_label("Similarity", rotation=270, labelpad=12)

    fig.suptitle(title, fontsize=12)
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, bbox_inches="tight")
    fig.savefig(out_png.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

# ----------------- orchestration helpers -----------------
def synthesize_mfi_from_antigens(antigens: List[str], out_csv: Path) -> Path:
    """Create a minimal 'MFI-like' CSV with antigen columns (one empty row)."""
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame([{}], columns=antigens)  # one row of NAs across antigen columns
    df.to_csv(out_csv, index=False)
    return out_csv

def call_r_build_map(rscript: str, r_build_script: Path, mfi_csv: Path, eplet_ref_dir: Path, verified_only: bool) -> Path:
    eplet_ref_dir.mkdir(parents=True, exist_ok=True)
    cmd = [rscript, str(r_build_script), str(mfi_csv), str(eplet_ref_dir)]
    if verified_only:
        cmd.append("--verified-only")
    print("[R] ", " ".join(shlex.quote(c) for c in cmd))
    cp = subprocess.run(cmd, capture_output=True, text=True)
    if cp.returncode != 0:
        sys.stderr.write(cp.stdout + "\n" + cp.stderr + "\n")
        raise SystemExit(f"R builder failed (exit {cp.returncode}).")
    print(cp.stdout)
    print(cp.stderr)
    map_path = eplet_ref_dir / "eplet_map.csv"
    if not map_path.exists():
        raise SystemExit(f"R builder completed but eplet_map.csv not found in {eplet_ref_dir}")
    return map_path

def build_incidence_from_map(map_path: Path, out_csv: Path) -> Path:
    # Import the local builder function
    from importlib import util
    builder_py = Path(__file__).parent / "build_eplet_incidence_A.py"
    spec = util.spec_from_file_location("build_incidence_mod", builder_py)
    mod = util.module_from_spec(spec); spec.loader.exec_module(mod)  # type: ignore
    inc = mod.build_incidence(pd.read_csv(map_path))
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    inc.to_csv(out_csv)
    return out_csv

# ----------------- main -----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--antibodies", required=True, help="Path to antibodies_python.csv")
    ap.add_argument("--outdir", default="outputs/similarity_antibodies", help="Output directory")

    # Optional inputs / builders
    ap.add_argument("--mfi", default=None, help="Optional: existing MFI CSV with antigen columns")
    ap.add_argument("--eplet-ref-dir", default="outputs/eplet_ref", help="Where to write/read eplet_map.csv")
    ap.add_argument("--rscript", default="Rscript", help="Path to Rscript (on PATH or absolute)")
    ap.add_argument("--r-build-script", default=None, help="Path to buildepletmap_A.R (default: scripts/buildepletmap_A.R)")
    ap.add_argument("--verified-only", action="store_true", help="Pass --verified-only to the R builder")

    args = ap.parse_args()
    ab_path = Path(args.antibodies)
    outdir  = Path(args.outdir); outdir.mkdir(parents=True, exist_ok=True)
    eplet_ref_dir = Path(args.eplet_ref_dir)

    # Resolve R builder script
    r_build_script = Path(args.r_build_script) if args.r_build_script else (Path(__file__).parent / "buildepletmap_A.R")
    if not r_build_script.exists():
        raise SystemExit(f"R builder not found at: {r_build_script}")

    # Detect antigens present in antibodies CSV
    if not ab_path.exists():
        raise SystemExit(f"Antibodies CSV not found: {ab_path}")
    df_ab = pd.read_csv(ab_path, encoding="utf-8", keep_default_na=False)
    antigens = [normalize_antigen(a) for a in detect_antigens_from_df(df_ab)]
    if len(antigens) < 2:
        raise SystemExit("Fewer than 2 HLA antigens detected in antibodies CSV; nothing to cluster.")

    # Ensure eplet_map.csv exists (build via R as needed)
    eplet_map = eplet_ref_dir / "eplet_map.csv"
    if eplet_map.exists():
        print(f"[OK] Using existing map: {eplet_map}")
    else:
        if args.mfi:
            mfi_csv = Path(args.mfi)
            if not mfi_csv.exists():
                raise SystemExit(f"--mfi provided but not found: {mfi_csv}")
        else:
            # Synthesize a minimal MFI with antigen columns only
            tmpdir = Path(tempfile.mkdtemp(prefix="mfi_synth_"))
            mfi_csv = synthesize_mfi_from_antigens(antigens, tmpdir / "synthetic_mfi.csv")
            print(f"[AUTO] Synthesized MFI with {len(antigens)} antigen columns at: {mfi_csv}")

        # Call R builder
        eplet_map = call_r_build_map(args.rscript, r_build_script, mfi_csv, eplet_ref_dir, args.verified_only)

    # Build incidence
    inc_path = eplet_ref_dir / "similarity" / "antigen_x_eplet_incidence.csv"
    if inc_path.exists():
        print(f"[OK] Using existing incidence: {inc_path}")
        inc = pd.read_csv(inc_path, index_col=0)
    else:
        inc_path.parent.mkdir(parents=True, exist_ok=True)
        inc_path = build_incidence_from_map(eplet_map, inc_path)
        inc = pd.read_csv(inc_path, index_col=0)
        print(f"[SAVE] {inc_path.resolve()}")

    # Filter to present antigens
    inc.index = [normalize_antigen(i) for i in inc.index]
    present = [a for a in antigens if a in inc.index]
    missing = [a for a in antigens if a not in inc.index]
    if len(present) < 2:
        raise SystemExit("Too few antibodies matched to eplet incidence; check antigen tokens and registry mapping.")
    inc_sub = inc.loc[present].copy()

    # (a) IDF-weighted Jaccard
    w_idf = idf_weights(inc_sub)
    S_jac = weighted_jaccard_matrix(inc_sub, w_idf)
    S_jac.to_csv(outdir / "EpletIDFWeighted_Jaccard_matrix_antibodies.csv")
    pd.DataFrame(1.0 - np.clip(S_jac.values, 0, 1), index=S_jac.index, columns=S_jac.columns)\
      .to_csv(outdir / "EpletIDFWeighted_Jaccard_distance_antibodies.csv")
    save_heatmap_with_dendrogram(
        list(S_jac.index), S_jac.values,
        outdir / "EpletIDFWeighted_Jaccard_heatmap_antibodies.png",
        "Eplet similarity (IDF‑weighted Jaccard) across antibodies"
    )

    # (b) Overlap coefficient
    S_ovlp = overlap_coefficient_matrix(inc_sub)
    S_ovlp.to_csv(outdir / "EpletOverlap_SzymkiewiczSimpson_matrix_antibodies.csv")
    pd.DataFrame(1.0 - np.clip(S_ovlp.values, 0, 1), index=S_ovlp.index, columns=S_ovlp.columns)\
      .to_csv(outdir / "EpletOverlap_SzymkiewiczSimpson_distance_antibodies.csv")
    save_heatmap_with_dendrogram(
        list(S_ovlp.index), S_ovlp.values,
        outdir / "EpletOverlap_SzymkiewiczSimpson_heatmap_antibodies.png",
        "Eplet similarity (Szymkiewicz–Simpson overlap) across antibodies"
    )

    # Manifest
    with open(outdir / "similarity_MANIFEST.txt", "w") as fh:
        fh.write("Antibody eplet similarity build (from scratch)\n")
        fh.write(f"Antibodies CSV: {ab_path.resolve()}\n")
        fh.write(f"Eplet ref dir: {eplet_ref_dir.resolve()}\n")
        fh.write(f"Eplet map:     {eplet_map.resolve()}\n")
        fh.write(f"Incidence:     {inc_path.resolve()}\n")
        fh.write(f"#Antibodies detected: {len(antigens)}\n")
        fh.write(f"#Matched in incidence: {len(present)}\n")
        if missing:
            fh.write("\nMissing (not in incidence):\n")
            for m in missing:
                fh.write(f"  - {m}\n")

    print(f"[DONE] Wrote outputs to {outdir.resolve()}")

if __name__ == "__main__":
    main()
