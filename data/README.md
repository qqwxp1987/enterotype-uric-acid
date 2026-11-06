# Data Inputs

Scripts locate repository files through `here::here()` or paths resolved from the script location; no local absolute-path editing is required.

## Analysis-ready files (`data/clean/`)

### Clinical and phenotype metadata

- `Metadata_20240925.csv` / `.txt`: clinical and phenotype fields used across the analyses. The matched-metabolomics table contains 303 records.
- `results/tables/section_3_3/section3_6_analysis_data.csv`: the 305-participant clinical analysis table used for the regression and QC workflows.

### Gut metagenome — phyloseq objects

Each object bundles feature abundances and sample metadata.

- `phyloseq_Taxa_with_enterotype_seed10.rds`: taxonomic composition and the prespecified seed-10 enterotype assignments.
- `phyloseq_Pathway_with_enterotype_20250210.rds`: MetaCyc pathway abundances.
- `phyloseq_KO_with_enterotype_20250210.rds`: KEGG Ortholog abundances.
- `phyloseq_EC_with_enterotype_20250210.rds`: Enzyme Commission abundances.

### Serum metabolome

- `metabolites_normalized_multiIS_20250904.csv` / `.txt`: normalized metabolite feature table.
- `metabolites_normalized_multiIS_and_metadata_20250904.rds`: normalized metabolome merged with sample metadata.
- `metabolome_20250904.rds`: packaged metabolome dataset object.

### Stratified functional tables and QC

- `pathabundance_stratified.tsv`: per-species stratified MetaCyc pathway abundances used for taxonomic attribution and denominator checks.
- `kneaddata_read_count_table.tsv`: per-sample KneadData read counts. It contains 325 sequenced samples; the QC builder restricts calculations to the 305 project samples and retains the other records for provenance.

### Baseline dataset for DMM

- `UA_Analysis_20230811.rds`: baseline phyloseq object used by `scripts/figure1/panel_A_dmm_selection.R` and the canonical enterotype-robustness workflow.

## Frozen purine curation (`data/curated/purine/`)

- `frozen_inclusion_register.csv`: 474-row machine-readable register, including the 34 included purine-related pathways, module assignments, curation roles, and module-score membership.
- `pathway_composition.csv`: exact pathway-to-EC composition used to construct double-count-free EC-union module scores.

These are analysis inputs, not files to be recreated during an ordinary results rerun.

## Large or raw data

The raw metagenomic and metabolomic datasets are deposited in CNSA under accession **CNP0009610**.

The full `ec_stratified.tsv` (~137 MB) is not redistributed in this repository. The derived files needed for the published EC-attribution results are already included:

- `results/tables/section_3_4/attribution/ec_stratified_subset.rds`
- `results/tables/section_3_4/attribution/prep_data.rds`

Only `scripts/figure3/attribution/step1_prep_subsets.R`, which regenerates those subsets from the original full HUMAnN output, requires `data/clean/ec_stratified.tsv`. All downstream attribution and final-figure steps use the included subsets.
