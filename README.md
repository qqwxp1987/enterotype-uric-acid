# Gut Microbiome Enterotypes and Hyperuricemia — Analysis Code and Key Results

This repository provides the analysis code, analysis-ready inputs, and key numerical results for:

> **Gut Microbiome Enterotypes Associated with Hyperuricemia Exhibit Distinct Multi-Omics Architectures and Concordant Purine-Related Metagenomic Signatures**

The study includes 305 adults with shotgun metagenomics, of whom 303 have matched serum untargeted metabolomics. The public workflow covers enterotype definition, UA/HUA associations and robustness models, multi-omics core features and taxonomic attribution, purine-module analysis, microbiome–metabolome integration, exploratory mediation, and within-Type 4 heterogeneity.

## What is included

- **One authoritative version of each required analysis script**, with portable repository-relative paths and fixed random seeds where applicable.
- **Final main and supplementary figures** as PDFs.
- **Directly queryable final CSVs** for Tables S1–S19 in `results/tables/final/`.
- **The final 29-sheet supplementary workbook**, including the Contents sheet and companion sheets S8a–S8i.
- **Frozen statistical outputs and audit records** needed to query the published results without immediately repeating long bootstrap or clustering runs.

This is a **main-results release**, not an archive of every exploratory analysis. Superseded versions, reviewer-only inspection files, manuscript/response generators, numbering scripts, PDF/image QA utilities, and obsolete preliminary analyses are intentionally excluded. See [ANALYSIS_MAP.md](ANALYSIS_MAP.md) for the inclusion rationale and the figure/table-to-script mapping.

## Repository layout

```text
data/
  clean/                  analysis-ready clinical, metagenomic and metabolomic inputs
  curated/purine/         frozen pathway-composition and inclusion registers
scripts/
  figure1/ ... figure6/   statistical analyses and final figure generation
  utils/                  shared plotting helpers
  verification/           self-consistency and key-result checks
results/
  tables/final/           canonical final-numbered CSVs (S1–S19 and S8a–S8i)
  tables/                 upstream statistical outputs and small required objects
  figures/                final Figures 1–6 (PDF)
  r1_2_r2_1/              enterotype-robustness protocol and audit records
  r1_5_covariate_sensitivity/
supplementary/
  figures/                final Figures S1–S7 (PDF)
  tables/                 final supplementary workbook
  build_supplementary_xlsx.py
```

## Quick result lookup

- Main figures: `results/figures/Figure1.pdf` through `Figure6.pdf`.
- Supplementary figures: `supplementary/figures/`.
- Final tables: `results/tables/final/`.
- Workbook: `supplementary/tables/Supplementary_Tables.xlsx`.
- Exact script mapping: [ANALYSIS_MAP.md](ANALYSIS_MAP.md).
- Validation scope and limitations: [VALIDATION.md](VALIDATION.md).
- File hashes: `MANIFEST.sha256`.

## Software requirements

- **R 4.3.1** was used for the manuscript analyses.
- Core R packages include `here`, `phyloseq`, `microbiome`, `DirichletMultinomial`, `vegan`, `MicrobiomeStat`, `limma`, `nlme`, `logistf`, `emmeans`, `ppcor`, `mediation`, `tidyverse`, `ComplexHeatmap`, `patchwork`, and `magick`.
- **Python 3** with `pandas`, `numpy`, `openpyxl`, and `pypdf` is required for table builders and release validation.

Package versions matter for exact numerical reproduction, particularly for LinDA, GLS optimization, Firth models, and bootstrap procedures. See the `session_info.txt` files under the relevant results directories and [VALIDATION.md](VALIDATION.md) before interpreting a rerun.

## Running selected analyses

Run commands from the repository root. R scripts resolve paths with `here::here()`; Python scripts resolve the repository root from their own location.

```bash
# Nested-model decomposition and Firth stability (Figure 3; Tables S6–S7)
Rscript scripts/figure3/nested_model_stability.R
python scripts/figure3/build_nested_model_tables.py

# Purine-module scores (Figure S6; Tables S3, S9, S10)
Rscript scripts/figure4/purine_module_analysis.R
python scripts/figure4/build_purine_module_tables.py
Rscript scripts/figure4/figureS6_purine_modules.R

# Continuous-gradient and k=3 sensitivity tables (Tables S5 and S19)
python scripts/figure1/build_gradient_sensitivity_tables.py

# QC and missingness tables (Tables S1–S2)
python scripts/figure2/build_qc_missingness_tables.py

# Rebuild the supplementary workbook from the final CSVs
python supplementary/build_supplementary_xlsx.py

# Check release structure and key numerical anchors
python scripts/verification/validate_release.py
```

Long-running steps—5,000-replicate mediation bootstraps, 100 DMM bootstrap refits, and 100 DMM initialization refits—do not need to be repeated merely to inspect the published results. Their frozen outputs are retained. Figure assembly scripts are mapped separately from statistical scripts in [ANALYSIS_MAP.md](ANALYSIS_MAP.md).

## Data availability

The metagenomic and metabolomic datasets generated and analyzed in the study are deposited in the **China National GeneBank Sequence Archive (CNSA), accession CNP0009610**.

Analysis-ready inputs are documented in `data/README.md`. The full `ec_stratified.tsv` file is not redistributed here because of its size. The derived EC-stratified subset required for the published attribution analyses is included under `results/tables/section_3_4/attribution/`; only regenerating that subset from scratch requires the original full HUMAnN output.

## Citation and license

Please cite the associated paper when using the code or data. Code is released under the MIT License; see `LICENSE`.
