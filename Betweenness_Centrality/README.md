# SCENIC & NetworkX reproducibility pack

Great thanks to the (py)SCENIC authors

Van de Sande B, Flerin C, Davie K, De Waegeneer M, Hulselmans G, Aibar S, Seurinck R, Saelens W, Cannoodt R, Rouchon Q, Verbeiren T, De Maeyer D, Reumers J, Saeys Y, Aerts S. A scalable SCENIC workflow for single-cell gene regulatory network analysis. Nat Protoc. 2020;15(7):2247-76.

This repository contains three analysis scripts to reproduce:
1) **SCENIC on CTpre groups (union run, 10×) + aggregation + pairwise meta‑enrichment + dotplots**  
2) **SCENIC on LLP overlap sequences (union run, 10×) + aggregation + pre/post meta‑enrichment + dotplots**  
3) **NetworkX consensus/overlap graphs across multiple SCENIC_CAR runs (CAR_CD4)**

All scripts are CLI tools with configurable paths and parameters. They are designed to run against the Zenodo reference bundle

> 🔗 **Data bundle**: 

---

## Quick start

### 1) Create and activate environment

Using **conda** (recommended):

```bash
mamba env create -f env/environment.yml
conda activate scenic-nx

