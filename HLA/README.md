
# HLA eplet similarity — build from scratch
Eplet definitions were obtained from the International HLA Epitope (Eplet) Registry using the community hlapro R package. It was created under the auspices of the 16th International HLA and Immunogenetics Workshop (IHIW) and is maintained by the EpRegistry team.

Great thanks to

Duquesnoy RJ, Marrari M, Tambur AR, Mulder A, Sousa LC, da Silva AS, do Monte SJ. First report on the antibody verification of HLA-DR, HLA-DQ and HLA-DP epitopes recorded in the HLA Epitope Registry. Hum Immunol. 2014 Nov;75(11):1097-103. doi: 10.1016/j.humimm.2014.09.012. Epub 2014 Oct 13. PMID: 25305456.

&

Bezstarosti S, Bakker KH, Kramer CSM, de Fijter JW, Reinders MEJ, Mulder A, Claas FHJ, Heidt S. A Comprehensive Evaluation of the Antibody-Verified Status of Eplets Listed in the HLA Epitope Registry. Front Immunol. 2022 Jan 28;12:800946. doi: 10.3389/fimmu.2021.800946. PMID: 35154076; PMCID: PMC8831796.



This repo lets a user start **only with** `antibodies_python.csv` and produce:
- `eplet_map.csv` (via R + hlapro)
- `antigen_x_eplet_incidence.csv` (Python expansion)
- Similarity matrices + clustered heatmaps across the antibodies:
  - IDF‑weighted Jaccard
  - Szymkiewicz–Simpson (overlap)

## Setup

```bash
# Python tools
conda env create -f envs/python-hla.yml
conda activate hla-clustering-py

# R builder
conda env create -f envs/r-hlapro.yml
# You can run the R builder via 'conda run' or expose its Rscript on PATH.

# If Rscript is on PATH (e.g., your default R):
python scripts/eplet_similarity_for_antibodies_A.py \
  --antibodies antibodies_python.csv \
  --outdir outputs/similarity_antibodies \
  --eplet-ref-dir outputs/eplet_ref \
  --verified-only

# If using the conda R env (without altering PATH):
python scripts/eplet_similarity_for_antibodies_A.py \
  --antibodies antibodies_python.csv \
  --outdir outputs/similarity_antibodies \
  --eplet-ref-dir outputs/eplet_ref \
  --rscript "$(conda run -n hla-eplet-r which Rscript)" \
  --verified-only

