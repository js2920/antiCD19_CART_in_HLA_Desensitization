# single‑cell TCR inputs (Dandelion)

The authors of Dandelion and Immcantation have to be acknowledged and thanked at the beginning

Suo C, Polanski K, Dann E, Lindeboom RGH, Vilarrasa-Blasi R, Vento-Tormo R, Haniffa M, Meyer KB, Dratva LM, Tuong ZK, Clatworthy MR, Teichmann SA. Dandelion uses the single-cell adaptive immune receptor repertoire to explore lymphocyte developmental origins. Nat Biotechnol. 2024;42(1):40-51.

&

Gabernet G, Marquez S, Bjornson R, Peltzer A, Meng H, Aron E, Lee NY, Jensen CG, Ladd D, Polster M, Hanssen F, Heumos S, nf-core c, Yaari G, Kowarik MC, Nahnsen S, Kleinstein SH. nf-core/airrflow: An adaptive immune receptor repertoire analysis workflow employing the Immcantation framework. PLoS Comput Biol. 2024;20(7):e1012265.
TCR single cell processing using the  **https://doi.org/10.5281/zenodo.17084542** repository or the **GEO GSE307464 repository**
Files structure from Zenodo zip archive

```text
dandelion_inputs_TCR/
|-- analysis_ready/
|   |-- meta.csv
|   |-- per_cell_TCR_table.csv
|   |-- sc-dandelion_latest.sif
|   |-- adata_with_TCR_integration.h5ad
|   |-- processed_adata_NTX_all_samples.h5ad
|   |-- NTX_1_TCR_fastq/
|   |   `-- dandelion/all_contig_dandelion.tsv
|   |-- NTX_2_TCR_fastq/
|   |   `-- dandelion/all_contig_dandelion.tsv
|   |-- NTX_3_TCR_fastq/
|   |   `-- dandelion/all_contig_dandelion.tsv
|   |-- NTX_4_TCR_fastq/
|   |   `-- dandelion/all_contig_dandelion.tsv
|   `-- NTX_5_TCR_fastq/
|       `-- dandelion/all_contig_dandelion.tsv
|-- base_inputs/
|   |-- adata/processed_adata_NTX_all_samples.h5ad
|   |-- contigs/
|   |   |-- NTX_1.all_contig_dandelion.tsv
|   |   |-- NTX_2.all_contig_dandelion.tsv
|   |   |-- NTX_3.all_contig_dandelion.tsv
|   |   |-- NTX_4.all_contig_dandelion.tsv
|   |   `-- NTX_5.all_contig_dandelion.tsv
|   |-- containers/sc-dandelion_latest.sif
|   |-- germlines/corrected_germline.fasta
|   |-- igblast_internal_data/
|   `-- meta/meta.csv
|-- README_dandelion_inputs_TCR.txt
`-- manifest_dandelion_TCR.tsv
```

Rest follows BCR examples
