# Gut Microbiome Enterotypes and Serum Uric Acid — Analysis Code & Data

This repository contains the data and analysis code required to reproduce the
figures and tables of the study:

> **Gut Microbiome Enterotypes Associate with Serum Uric Acid through Distinct
> Multi-Omics Architectures**

We analysed 305 adults from a health-examination cohort with matched clinical
phenotyping and shotgun metagenomics (303 with serum untargeted metabolomics).
Dirichlet multinomial mixture (DMM) modelling identified four genus-level
enterotypes along a *Bacteroides*–*Prevotella* gradient; the code here covers
enterotype definition, enterotype–uric-acid association and robustness testing,
multi-omics core-feature extraction, microbiome–metabolome integration, and
causal mediation analysis.

---

## Repository structure

```
.
├── data/
│   ├── README.md            # description of every data file
│   └── clean/               # analysis-ready inputs (phyloseq objects, metabolome, metadata)
├── scripts/
│   ├── figure1/ … figure6/  # analysis & plotting code, grouped by manuscript figure
│   └── utils/               # shared helpers (colour palettes)
├── results/
│   ├── tables/              # intermediate statistical result tables
│   ├── figures/             # final composite figures + individual panels
│   └── cache/               # cached intermediate plot objects (.rds)
├── supplementary/
│   ├── figures/             # supplementary figures S1–S3
│   └── tables/              # supplementary tables + builder script
├── README.md
├── LICENSE                  # MIT
└── .gitignore
```

> **Note on folder names.** The `scripts/figureN/` folders reflect the analysis
> history; because figures were renumbered during revision, a folder may feed
> more than one manuscript figure. The **figure → script map** below is the
> authoritative reference.

---

## Data

Analysis-ready inputs are under `data/clean/` (see `data/README.md` for a full
description of each file). Two categories of data are hosted **externally**
because of size or sequence-data deposition requirements:

| Data | Location |
|---|---|
| Raw shotgun metagenomic sequencing reads | NGDC Genome Sequence Archive — accession `[GSA: CRAXXXXXX]` |
| Serum untargeted metabolomics (raw) | NGDC OMIX — accession `[OMIX: OMIXXXXXX]` |
| `ec_stratified.tsv` — per-species stratified EC abundance table (~137 MB) | Zenodo — DOI `[10.5281/zenodo.XXXXXXX]` |

To run the EC-attribution scripts, download `ec_stratified.tsv` from Zenodo and
place it in `data/clean/`.

---

## Software requirements

- **R ≥ 4.2** (analyses were run on Windows; `here::here()` makes paths portable).
- Core R packages:
  - Microbiome / stats: `phyloseq`, `microbiome`, `DirichletMultinomial`,
    `vegan`, `MicrobiomeStat`, `limma`, `nlme`, `car`, `emmeans`, `mediation`
  - Data handling & plotting: `tidyverse` (`dplyr`, `ggplot2`, `readr`,
    `stringr`, `tidyr`), `here`, `ComplexHeatmap`, `circlize`, `magick`
- A small number of scripts call the KEGG REST API for enzyme annotation
  (`scripts/figure3/attribution/fetch_kegg_ecs.py`, Python 3).

Installing the exact package versions used for the paper (recommended: record
your own `sessionInfo()` or use `renv`) is advised for bit-for-bit reproduction.

---

## How to reproduce

The project uses the [`here`](https://here.r-lib.org/) package: open the
RStudio project (`HUA_FCTJ2021.Rproj`) or set the working directory to the
repository root, and all relative paths resolve automatically.

Within each figure folder, run the scripts in this order:

1. `section*_step1_*` / `panel_*` data-preparation and per-panel scripts;
2. `assemble_figure*.R` — stitches the panels into the final composite figure.

Outputs are written to `results/figures/`, `results/tables/` and
`supplementary/`.

### Figure / table → script map

| Manuscript item | Script(s) |
|---|---|
| **Figure 1** (enterotype landscape) | `scripts/figure1/panel_A`–`panel_E` + `assemble_figure1.R` |
| **Figure 2** (UA / HUA association) | `scripts/figure2/panel_A`–`panel_D` + `assemble_figure2.R` |
| **Figure 3** (robustness adjustment) | `scripts/figure3/section3_6_step1`–`step4` |
| **Figure 4** (multi-omics core signatures) | `scripts/figure3/section3_3_differential_analysis.R`, `panel_C_core_pathways_heatmap_rotated.R`, `panel_GH_convergent_ec.R`, `attribution/`; `scripts/figure4/panel_A_core_summary.R` + `assemble_figure4.R` |
| **Figure 5** (microbiome–metabolome integration) | `scripts/figure5/panel_A_dbRDA.R`, `panel_B_top_adjusted.R`, `panel_CD_function_integration.R` + `assemble_figure5.R` |
| **Figure 6** (mediation & heterogeneity) | `scripts/figure6/section3_5_step1`–`step3` + `panel_A`–`panel_D` + `assemble_figure6.R` |
| **Table 1** (cohort demographics) | `scripts/figure2/table1_upgrade.R` |
| **Figure S1** (per-subject composition) | `scripts/figure1/suppl_S1_individual_bar.R` |
| **Figure S2** (additional core-feature panels) | `scripts/figure3/panel_B_radial_taxa.R`, `panel_E_Species_Attribution_Bar.R`; `scripts/figure4/panel_D_type1_core_metabolites.R` + `assemble_figureS2.R` |
| **Figure S3** (full 75-species heatmap) | `scripts/figure5/panel_B_full_adjusted.R` |
| **Tables S1–S6** | `supplementary/build_supplementary_xlsx.R` |

---

## Citation

If you use this code or data, please cite the associated paper (full citation
to be added on publication) and the Zenodo archive of this repository
(DOI `[10.5281/zenodo.XXXXXXX]`).

## License

Code in this repository is released under the MIT License (see `LICENSE`).

## Contact

For questions about the analyses, please contact the corresponding author of
the paper.
