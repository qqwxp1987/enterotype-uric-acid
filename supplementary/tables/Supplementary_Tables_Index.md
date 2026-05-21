# Supplementary Tables — HUA_FCTJ2021

Numbered Supplementary Tables S1–S6, compiled from existing result CSVs under
`results/tables/`. Source CSVs were not modified.

> Note on format: the approved plan (P4) preferred a single multi-sheet
> `Supplementary_Tables.xlsx`. Code execution (Python / pandas / openpyxl) was
> not available in this environment, so the fallback CSV deliverables below were
> produced instead. A ready-to-run builder, `_build_supplementary_tables.py`, is
> included; running it where pandas + openpyxl are installed will regenerate the
> equivalent `Supplementary_Tables.xlsx` (one sheet per table plus a Contents
> sheet).

## Table list

| Table | Title | Description | File |
|---|---|---|---|
| Table S1 | Independent core multi-omics features by enterotype | Counts of enterotype-independent core features per omics layer (taxa, pathway, EC, metabolome), with per-layer E1-only / E4-only / shared core feature lists. | `Table_S1_core_features.csv` (+ companion lists below) |
| Table S2 | Causal mediation analysis — main results (Type 1 and Type 4 models) | Bootstrapped causal mediation estimates (ACME, ADE, total effect, proportion mediated) for each candidate mediator in the Type 1 vs reference and Type 4 vs reference models; BH-adjusted ACME p-values included. | `Table_S2_mediation_main.csv` |
| Table S3 | Mediation sensitivity analysis — binary hyperuricemia outcome | Mediation models refitted with a binary hyperuricemia (HUA) outcome instead of continuous serum urate. | `Table_S3_mediation_binary_sensitivity.csv` |
| Table S4 | Mediation sensitivity analysis — reference-group split (Type 2 vs Type 3) | Mediation analysis splitting the reference group: Type 1 and Type 4 each compared against Type 2 (E2) and Type 3 (E3) separately. | `Table_S4_mediation_refgroup_sensitivity.csv` |
| Table S5 | Mediation sensitivity analysis — sequential-ignorability (medsens) | Sequential-ignorability sensitivity analysis (medsens): sensitivity parameter rho at which the mediated effect changes sign (NA = no sign change within the tested range). | `Table_S5_mediation_medsens.csv` |
| Table S6 | Within-Type 4 heterogeneity — HUA+ vs HUA− comparison | Within-Type 4 comparison of candidate mediators between HUA-positive and HUA-negative individuals (Wilcoxon rank-sum test, Cliff's delta). | `Table_S6_type4_heterogeneity.csv` |

## Table S1 companion per-layer core feature lists

`Table_S1_core_features.csv` gives the per-layer counts. The following files
provide the actual feature lists (DESeq2/limma/LinDA log2 fold changes and
BH-adjusted p-values for the four enterotype contrasts E1vE2, E1vE3, E4vE2,
E4vE3):

| File | Layer | Set | n |
|---|---|---|---|
| `Table_S1a_taxa_E1_only.csv` | Taxa | E1-only core | 72 |
| `Table_S1b_taxa_Shared.csv` | Taxa | Shared core | 1 |
| `Table_S1c_taxa_E4_only.csv` | Taxa | E4-only core | 2 |
| `Table_S1d_pathway_E1_only.csv` | Pathway | E1-only core | 43 |
| `Table_S1e_pathway_E4_only.csv` | Pathway | E4-only core | 52 |
| `Table_S1f_metabolome_E1_only.csv` | Metabolome | E1-only core | 24 |

EC-layer core lists are large (E1-only n=296, E4-only n=371, shared n=36) and
were not re-typed into this folder; they remain available unchanged in the
source directory as `results/tables/section_3_4/ec_E1_only.csv`,
`ec_E4_only.csv` and `ec_Shared.csv` (see source paths below).

## Source CSV paths used

| Compiled table | Source CSV (read-only) |
|---|---|
| Table S1 (summary) | `results/tables/section_3_4/independent_core_summary.csv` |
| Table S1 (taxa lists) | `results/tables/section_3_4/taxa_E1_only.csv`, `taxa_Shared.csv`, `taxa_E4_only.csv` |
| Table S1 (pathway lists) | `results/tables/section_3_4/pathway_E1_only.csv`, `pathway_E4_only.csv` |
| Table S1 (metabolome list) | `results/tables/section_3_4/metabolome_E1_only.csv` |
| Table S1 (EC lists, not re-copied) | `results/tables/section_3_4/ec_E1_only.csv`, `ec_E4_only.csv`, `ec_Shared.csv` |
| Table S2 | `results/tables/section_3_6/mediation_main_results.csv` |
| Table S3 | `results/tables/section_3_6/mediation_binary_sensitivity.csv` |
| Table S4 | `results/tables/section_3_6/mediation_refgroup_sensitivity.csv` |
| Table S5 | `results/tables/section_3_6/mediation_medsens_results.csv` |
| Table S6 | `results/tables/section_3_6/type4_heterogeneity_results.csv` |

### Notes on source selection (files present in more than one location)

The mediation/heterogeneity files exist in both `section_3_5/` and
`section_3_6/`. Both copies were inspected:

- `mediation_main_results.csv` — identical in both directories. Used `section_3_6/`.
- `mediation_binary_sensitivity.csv` — the two copies differ. The `section_3_5/`
  copy has coarse p-values (multiples of 0.02, indicating ~50 bootstrap
  replicates); the `section_3_6/` copy has fine-grained p-values (~2500
  replicates). The more complete/authoritative `section_3_6/` copy was used.
- `mediation_refgroup_sensitivity.csv` — both copies use ~2500 replicates but
  come from different bootstrap runs (slightly different CI/p values).
  `section_3_6/` used for consistency with S2/S3.
- `mediation_medsens_results.csv` — identical in both directories. Used `section_3_6/`.
- `type4_heterogeneity_results.csv` — identical in both directories. Used `section_3_6/`.

No source CSV was missing or empty.
