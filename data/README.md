# Data Directory

Analysis-ready inputs for reproducing the manuscript. All files live in
`data/clean/`. Scripts locate them through `here::here("data", "clean", ...)`,
so no path editing is needed when the repository root is the working directory.

## Files in this repository (`data/clean/`)

### Clinical / phenotype metadata
- **`Metadata_20240925.csv` / `.txt`** — participant phenotype table: `Enterotype`
  (DMM-derived Type 1–4), `Sex`, `Age`, `BMI`, `Hyperuricemia` status, serum
  `UA`, and other routine clinical indices. Used for demographic, regression and
  mediation analyses. (Identical content, two delimited formats.)

### Gut metagenome — phyloseq objects
Each `phyloseq` object bundles an abundance table, taxonomy/feature metadata and
the sample metadata (including enterotype and uric-acid fields).
- **`phyloseq_Taxa_with_enterotype_seed10.rds`** — taxonomic composition
  (Phylum→Species). Primary input for most differential-abundance analyses.
- **`phyloseq_Pathway_with_enterotype_20250210.rds`** — MetaCyc pathway abundances.
- **`phyloseq_KO_with_enterotype_20250210.rds`** — KEGG Ortholog abundances.
- **`phyloseq_EC_with_enterotype_20250210.rds`** — Enzyme Commission (EC) abundances.

### Serum metabolome
- **`metabolites_normalized_multiIS_20250904.csv` / `.txt`** — multi-internal-standard
  normalized metabolite feature table.
- **`metabolites_normalized_multiIS_and_metadata_20250904.rds`** — the normalized
  metabolome merged with sample metadata, ready for statistical analysis.
- **`metabolome_20250904.rds`** — packaged metabolome dataset object.

### Stratified functional tables (HUMAnN output)
- **`pathabundance_stratified.tsv`** (~25 MB) — per-species stratified MetaCyc
  pathway abundances; input for the taxonomic-attribution analyses.

### Early baseline dataset
- **`UA_Analysis_20230811.rds`** — 2023 baseline `phyloseq` object (305 samples,
  509 species, 65 clinical indices). **Still required**: it is the root input
  for the DMM enterotype-selection step (`scripts/figure1/panel_A_dmm_selection.R`),
  which derives `phyloseq_Taxa_with_enterotype_seed10.rds`.

## Files hosted externally

| File | Why external | Location |
|---|---|---|
| `ec_stratified.tsv` (~137 MB) | exceeds the GitHub 100 MB file-size limit | Zenodo, DOI `[10.5281/zenodo.XXXXXXX]` |
| Raw shotgun metagenomic reads | sequence-data deposition requirement | NGDC GSA, `[GSA: CRAXXXXXX]` |
| Raw serum metabolomics | metabolomics-data deposition requirement | NGDC OMIX, `[OMIX: OMIXXXXXX]` |

To run the EC-attribution scripts (`scripts/figure3/attribution/`,
`scripts/figure3/panel_GH_convergent_ec.R`), download `ec_stratified.tsv` from
Zenodo and place it in `data/clean/`.
