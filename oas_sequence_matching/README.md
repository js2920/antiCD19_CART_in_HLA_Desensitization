# see: https://opig.stats.ox.ac.uk/webapps/oas/downloads/
and their Publication Observed Antibody Space: A diverse database of cleaned, annotated, and translated unpaired and paired antibody sequences
Tobias H Olsen, Fergus Boyles, Charlotte M Deane   Protein Sci
2022 Jan;31(1):141-146.
doi: 10.1002/pro.4205. Epub 2021 Oct 29.
for downloading the Observed Antobody Space Dataset on which this script builds
envs
  sas-cdr3-db.yml          # pinned environment for database building
  cdr3-align.yml           # pinned environment for alignment + plots
scripts
  build_oas_ighg_cdr3_db.sh     # Step 1: build FASTA + MMseqs2 DB
  align_cdr3_against_oas.sh     # Step 2: align + merge + (optional) UMAPs
