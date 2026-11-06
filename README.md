# Gut Microbiome Enterotypes and Hyperuricemia

Analysis code, analysis-ready inputs, figures, and numerical results for:

> **Gut Microbiome Enterotypes Associated with Hyperuricemia Exhibit Distinct Multi-Omics Architectures and Concordant Purine-Related Metagenomic Signatures**

The study includes 305 adults with shotgun metagenomics, of whom 303 have matched serum untargeted metabolomics. Analyses cover enterotype definition, serum uric acid (UA) and hyperuricemia (HUA) associations, multi-omics core features, taxonomic attribution, purine-module scores, microbiome–metabolome integration, exploratory mediation, and within-enterotype heterogeneity.

## Repository structure

```text
data/
  clean/                       analysis-ready clinical, metagenomic and metabolomic inputs
  curated/purine/              pathway-composition and inclusion registers
scripts/
  figure1/ ... figure6/        statistical analyses and figure generation
  utils/                       shared plotting helpers
  verification/                structural and numerical checks
results/
  figures/                     Figures 1–6 (PDF)
  tables/final/                Tables S1–S19 and S8a–S8i (CSV)
  tables/                      upstream statistical outputs and required objects
  enterotype_robustness/       enterotype assignment and clustering robustness
  covariate_sensitivity/       covariate-adjustment sensitivity outputs
  verification/                machine-readable validation report
supplementary/
  figures/                     Figures S1–S7 (PDF)
  tables/                      supplementary workbook
  build_supplementary_xlsx.py
```

## Analysis workflows

| Manuscript section | Main outputs | Principal scripts |
|---|---|---|
| **3.1 Enterotype definition and community structure** | Figure 1; Figures S1–S2; Tables S5 and S19 | `scripts/figure1/` |
| **3.2 Enterotype associations with UA and HUA** | Figure 2; Table 1; Figure S3 | `scripts/figure2/`; `scripts/figure6/type4_detailed_characterization.R` |
| **3.3 Renal and cardiometabolic adjustment** | Figure 3; Tables S6–S7 | `scripts/figure3/` |
| **3.4 Multi-omics architectures and purine-related functions** | Figure 4; Figures S4–S6; Tables S3 and S8–S11 | `scripts/figure3/section3_3_differential_analysis.R`; `scripts/figure3/attribution/`; `scripts/figure4/`; `scripts/figure5/extended_covariate_sensitivity.R` |
| **3.5 Microbiome–metabolome integration** | Figure 5; Figure S7; Table S12 | `scripts/figure5/`; `scripts/figure6/extended_covariate_feature_attribution.R` |
| **3.6 Exploratory mediation and within-Type 4 heterogeneity** | Figure 6; Table S4; Tables S13–S18 | `scripts/figure6/` |
| **Data quality and analysis-sample accounting** | Tables S1–S2 | `scripts/figure2/build_qc_missingness_tables.py` |

Final figures are in `results/figures/` and `supplementary/figures/`. Final tables are in `results/tables/final/`; the corresponding 29-sheet workbook is `supplementary/tables/Supplementary_Tables.xlsx`.

## Software requirements

- **R 4.3.1** was used for the manuscript analyses.
- Core R packages: `here`, `phyloseq`, `microbiome`, `DirichletMultinomial`, `vegan`, `MicrobiomeStat`, `limma`, `nlme`, `logistf`, `emmeans`, `ppcor`, `mediation`, `tidyverse`, `ComplexHeatmap`, `patchwork`, and `magick`.
- **Python 3** with `pandas`, `numpy`, `openpyxl`, and `pypdf` for table builders and structural validation.

Package versions can affect numerical reproduction, particularly for LinDA, GLS optimization, Firth models, and bootstrap procedures. Relevant `session_info.txt` files are retained with the analysis outputs.

## Running analyses

Run commands from the repository root. R scripts resolve paths with `here::here()`; Python scripts resolve the repository root from their own location.

```bash
# Enterotype robustness, continuous-gradient and k=3 sensitivity
Rscript scripts/figure1/enterotype_robustness.R
python scripts/figure1/build_gradient_sensitivity_tables.py

# Clinical covariate-adjustment models and Firth stability
Rscript scripts/figure3/nested_model_stability.R
python scripts/figure3/build_nested_model_tables.py

# Purine-module scores and associated tables
Rscript scripts/figure4/purine_module_analysis.R
python scripts/figure4/build_purine_module_tables.py
Rscript scripts/figure4/figureS6_purine_modules.R

# QC and missingness tables
python scripts/figure2/build_qc_missingness_tables.py

# Rebuild the supplementary workbook from final CSVs
python supplementary/build_supplementary_xlsx.py

# Check the repository structure and key numerical anchors
python scripts/verification/validate_release.py
```

The mediation workflow uses 5,000 bootstrap resamples. Enterotype robustness uses 100 bootstrap refits and 100 initialization refits. Existing statistical outputs are provided so that the reported results can also be inspected without repeating these longer runs.

## Data availability

The metagenomic and metabolomic datasets are deposited in the **China National GeneBank Sequence Archive (CNSA), accession CNP0009610**.

Analysis-ready inputs are documented in `data/README.md`. The full `ec_stratified.tsv` file is not redistributed because of its size. The derived EC-stratified subset required for the published attribution analyses is included under `results/tables/section_3_4/attribution/`; only regenerating that subset from the original full HUMAnN output requires an additional input.

`MANIFEST.sha256` records hashes for the published file set. Code is released under the MIT License; see `LICENSE`.
