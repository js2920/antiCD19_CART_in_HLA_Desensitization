# SCENIC & NetworkX reproducibility pack

Great thanks to the (py)SCENIC authors

Van de Sande B, Flerin C, Davie K, De Waegeneer M, Hulselmans G, Aibar S, Seurinck R, Saelens W, Cannoodt R, Rouchon Q, Verbeiren T, De Maeyer D, Reumers J, Saeys Y, Aerts S. A scalable SCENIC workflow for single-cell gene regulatory network analysis. Nat Protoc. 2020;15(7):2247-76.

This repository contains three analysis scripts to reproduce:
1) **SCENIC on memory B and Plasma pre CAR T groups (union run, 10×) + aggregation + pairwise meta‑enrichment + dotplots**  
2) **SCENIC on LLP overlap sequences (union run, 10×) + aggregation + pre/post meta‑enrichment + dotplots**  
3) **NetworkX consensus/overlap graphs across multiple SCENIC_CAR runs (CAR_CD4)**

All scripts are CLI tools with configurable paths and parameters. They are designed to run against the Zenodo reference bundle **https://doi.org/10.5281/zenodo.17085448**

> 🔗 **Data bundle**: 

---

## Quick start

### 1) Create and activate environment

Using **conda** (recommended):

```bash
mamba env create -f env/environment.yml
conda activate scenic-nx

*unpack the Zenodo bundle*

mkdir -p refs
unzip scenic_reference_bundle_XXXXX.zip -d refs
# The archive expands into: refs/scenic_reference_bundle_XXXXX/
export REF_ROOT="$(pwd)/refs/scenic_reference_bundle_XXXXX"

a) analysis for B memory and Plasma cells

python scripts/SCENIC_B_Plasma.py \
  --in-h5ad        "$REF_ROOT/data/refined_plasma_memory_b.h5ad" \
  --out-base       "./results/SCENIC_memB_memdPCs_maturePCs" \
  --ranking-db     "$REF_ROOT/db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  --motif-table    "$REF_ROOT/db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
  --tf-list        "$REF_ROOT/db/allTFs_hg38.txt" \
  --n-workers 12 --n-iter 10

b) analysis for long lived plasma cells

python scripts/SCENIC_LLPCS.py \
  --anndata        "$REF_ROOT/data/processed_adata_NTX_all_samples.h5ad" \
  --cell-list-csv  "$REF_ROOT/metadata/Cl1_14PreCl1_14Post_overlap_sequences.csv" \
  --out-base       "./results/SCENIC_LLPs_combined" \
  --ranking-db     "$REF_ROOT/db/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather" \
  --motif-table    "$REF_ROOT/db/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl" \
  --tf-list        "$REF_ROOT/db/allTFs_hg38.txt" \
  --n-workers 12 --n-iter 10

c) Betweenness Centrality analysis by NetworkX for CD4+ CAR T cells

python scripts/NetworkX_CD4_CAR.py \
  --base    "$REF_ROOT/optional/scenic_car_runs" \
  --cell    "CAR_CD4" \
  --outdir  "./results/NetworkX_CD4_vsIterations" \
  --top-edges-per-tf 50 \
  --target-edge-range 1500 20000 \
  --overlap-fraction 0.6 \
  --strict-node-overlap-fraction 0.95

