#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge per-sample Dandelion `all_contig_dandelion.tsv`, build a single Dandelion
object, filter to productive chains, and integrate into a processed AnnData.

Outputs:
  <out_dir>/merged_contigs.tsv
  <out_dir>/adata_with_bcr_integration.h5ad

Template for GEO: GSE307464 (NTX_1..NTX_5).
"""
import os
import re
import argparse
import pandas as pd
import scanpy as sc
import dandelion as ddl


def parse_args():
    ap = argparse.ArgumentParser(description="Merge BCR contigs and integrate into AnnData")
    ap.add_argument("--contigs", required=True,
                    help="Comma-separated TSVs (e.g. '.../NTX_1/...tsv,.../NTX_2/...tsv').")
    ap.add_argument("--samples", default=None,
                    help="Comma-separated sample IDs matching --contigs (e.g. 'NTX_1,NTX_2'). "
                         "If omitted, IDs are inferred from file paths with regex '(NTX_\\d+)'.")
    ap.add_argument("--adata", required=True,
                    help="Processed AnnData H5AD (e.g. results/totalvi_celltypist/processed_adata_NTX_all_samples.h5ad).")
    ap.add_argument("--out-dir", default="results/bcr_integration",
                    help="Output folder (default: results/bcr_integration).")
    return ap.parse_args()


def rename_barcode(bc: str) -> str:
    return bc.replace("_fastq_", "_") + "-1" if "_fastq_" in bc else bc


def infer_id_from_path(path: str) -> str:
    m = re.search(r"(NTX_\d+)", os.path.basename(path)) or re.search(r"(NTX_\d+)", path)
    return m.group(1) if m else "UNKNOWN"


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    contig_paths = [p.strip() for p in args.contigs.split(",") if p.strip()]
    if args.samples:
        sample_ids = [s.strip() for s in args.samples.split(",")]
        if len(sample_ids) != len(contig_paths):
            raise SystemExit("[ERROR] --samples must have same length as --contigs")
    else:
        sample_ids = [infer_id_from_path(p) for p in contig_paths]

    # Load & merge with sample_id
    frames = []
    for p, sid in zip(contig_paths, sample_ids):
        df = pd.read_csv(p, sep="\t")
        df["sample_id"] = sid
        frames.append(df)
    merged = pd.concat(frames, ignore_index=True)
    merged["cell_id"] = merged["cell_id"].map(rename_barcode)
    merged_tsv = os.path.join(args.out_dir, "merged_contigs.tsv")
    merged.to_csv(merged_tsv, sep="\t", index=False)
    print(f"[INFO] merged contigs → {len(merged):,} rows\n      {merged_tsv}")

    # Build Dandelion object & integrate into AnnData
    vdj = ddl.Dandelion(merged)
    adata = sc.read_h5ad(args.adata)

    vdj, adata = ddl.pp.check_contigs(vdj, adata, productive_only=True)
    ddl.tl.transfer(adata, vdj, overwrite=True)

    out_h5ad = os.path.join(args.out_dir, "adata_with_bcr_integration.h5ad")
    adata.write(out_h5ad)
    print(f"[DONE] AnnData with BCR integrated → {out_h5ad}")


if __name__ == "__main__":
    main()
