# Manuscript-to-Code Map

## Scope decision

### Required for the main-results release

1. **All six main figures and Table 1**, including the statistical scripts that produce their reported estimates.
2. **Analyses needed to support the main conclusions**: enterotype robustness, clinical covariate adjustment, multi-omics core features, taxonomic attribution, purine-module scores, cross-omics integration, exploratory mediation, and Type 4 diagnostics.
3. **Second-round additions formally incorporated into the manuscript**: nested-model decomposition, Firth stability, continuous-gradient and k = 3 outcome sensitivity, and purine-module analysis.
4. **Final-numbered supplementary outputs** needed to query those results: Figures S1–S7 and Tables S1–S19.
5. **Final plotting code** where the submitted figure layout differs from an earlier revision.

### Intentionally not included

- Superseded analysis versions, duplicate copies, preliminary 20-iteration bootstrap workflows, the retired GAM-based nonlinearity interpretation, and 1,000-iteration mediation outputs replaced by 5,000-iteration analyses.
- Exploratory analyses that were not incorporated into the final manuscript, including the withdrawn k = 3 heterogeneity argument and contextual external-dataset extraction.
- Reviewer-only inspection bundles and verbose candidate-selection matrices when a concise final table is published.
- Word tracked-change tools, response-letter generators, figure/table renumbering utilities, PDF measurement scripts, contact sheets, and visual-QA scripts.
- Large plot caches, historical PNG/TIFF variants, local R session files, and raw-data downloads.
- The original full `ec_stratified.tsv`; the derived subset needed for all published downstream attribution results is included.

The exclusions reduce clutter; they do not remove any script required to query or recompute a main published result from the provided analysis-ready inputs.

## Main figures and Table 1

| Manuscript item | Statistical / analytical scripts | Final figure or table scripts | Canonical outputs |
|---|---|---|---|
| **Figure 1** — DMM enterotypes, PCoA, composition, P/B ratio, alpha diversity | `scripts/figure1/panel_A_dmm_selection.R` through `panel_E_alpha_diversity.R` | `scripts/figure1/assemble_figure1.R` | `results/figures/Figure1.pdf`; `results/tables/section_3_1/` |
| **Figure 2** — clinical screen, UA/HUA and adjusted regression | `scripts/figure2/panel_A_clinical_screening.R`, `panel_B_ua_violin.R`, `panel_C_hua_prevalence.R`, `panel_D_forest_plot.R` | `scripts/figure2/assemble_figure2.R` | `results/figures/Figure2.pdf`; `results/tables/section_3_2/` |
| **Table 1** — cohort characteristics | `scripts/figure2/table1_upgrade.R` | Same script | `results/tables/section_3_2/Table1_Demographics_v2.csv` |
| **Figure 3** — renal/cardiometabolic robustness and subgroup analyses | `scripts/figure3/section3_6_step1_data_preparation.R`, `section3_6_step2_main_models.R`, `section3_6_step3_sensitivity_analysis.R` | `scripts/figure3/section3_6_step4_visualization_v2.R` | `results/figures/Figure3.pdf`; `results/tables/section_3_3/` |
| **Tables S6–S7** — fixed-population nested models, EPV and Firth stability | `scripts/figure3/nested_model_stability.R` | `scripts/figure3/build_nested_model_tables.py` | `results/tables/final/Table_S6_nested_models.csv`; `Table_S7_model_stability.csv` |
| **Figure 4** — core multi-omics features and convergent EC attribution | `scripts/figure3/section3_3_differential_analysis.R`; `scripts/figure3/attribution/step1_prep_subsets.R`, `step2_pathway_contribution.R`, `step3_ec_contribution.R` | `scripts/figure4/assemble_figure4.R` | `results/figures/Figure4.pdf`; `results/tables/section_3_4/`; `results/tables/final/Table_S8*` |
| **Figure 5** — db-RDA and cross-omics partial correlations | `scripts/figure5/panel_A_dbRDA.R`, `panel_B_top_adjusted.R`, `panel_CD_function_integration.R` | `scripts/figure5/assemble_figure5.R` | `results/figures/Figure5.pdf`; `results/tables/section_3_4/` |
| **Figure 6** — exploratory mediation, sensitivity and within-Type 4 heterogeneity | `scripts/figure6/section3_5_step1_data_preparation.R`, `section3_5_step2_mediation_analysis.R`, `section3_5_step3_type4_heterogeneity.R` | `scripts/figure6/panel_A_mediation_forest.R` through `panel_D_type4_boxplots.R`; `assemble_figure6.R` | `results/figures/Figure6.pdf`; `results/tables/section_3_5/`; final Tables S13–S16 and S18 |

**Note on folder names:** historical `section3_6_*` scripts feed the robustness results now presented in Sections 3.2–3.3 and Figure 3. The table above, rather than the historical filename, is the authoritative mapping.

## Supplementary figures

| Final item | Script(s) | Notes |
|---|---|---|
| **Figure S1** — posterior assignment, bootstrap recovery, component-count sensitivity | `scripts/figure1/r12_r21_enterotype_robustness.R` | Canonical 100-bootstrap / 100-initialization workflow; historical internal output stem was Figure S4. |
| **Figure S2** — individual-level composition | `scripts/figure1/figureS2_individual_composition.R` | Historical internal stem was Figure S1. |
| **Figure S3** — Type 4 influence and leave-one-out diagnostics | `scripts/figure6/type4_detailed_characterization.R` | Historical internal stem was Figure S5. |
| **Figure S4** — core taxa, pathway attribution, metabolites and pathway heatmap | `scripts/figure3/figureS4A_radial_taxa.R`, `figureS4B_species_attribution.R`; `scripts/figure4/figureS4C_type1_core_metabolites.R`, `figureS4D_core_pathways_heatmap.R`; assembly: `scripts/figure4/assemble_figureS4.R` | Historical internal stem was Figure S2. |
| **Figure S5** — KEGG pathway-context schematic | `scripts/figure4/figureS5_purine_context.R` | Annotation schematic, not a statistical test; historical internal stem was Figure S7. |
| **Figure S6** — purine-module scores | `scripts/figure4/purine_module_analysis.R`; plotting: `scripts/figure4/figureS6_purine_modules.R` | Uses the frozen 34-pathway register and EC-union scores. |
| **Figure S7** — all-75-core-species sensitivity heatmap | `scripts/figure5/figureS7_full_species_heatmap.R` | Historical internal stems were Figure S3 and Figure S5. |

## Supplementary tables

| Final table(s) | Source / builder script(s) |
|---|---|
| **S1–S2** — QC, missingness and analysis-specific sample sizes | `scripts/figure2/build_qc_missingness_tables.py` |
| **S3, S9–S10** — purine curation, module statistics and module–UA associations | `scripts/figure4/purine_module_analysis.R`; `scripts/figure4/build_purine_module_tables.py` |
| **S4** — concise 22-candidate traceability | `scripts/figure6/build_mediator_traceability.R` |
| **S5, S19** — continuous P/B gradient and k = 3 outcome sensitivity | Statistics: `scripts/figure1/r12_r21_enterotype_robustness.R`; final tables: `scripts/figure1/build_gradient_sensitivity_tables.py` |
| **S6–S7** — nested models and Firth stability | `scripts/figure3/nested_model_stability.R`; `scripts/figure3/build_nested_model_tables.py` |
| **S8 and S8a–S8i** — independent core features and per-layer lists | `scripts/figure3/section3_3_differential_analysis.R`; attribution uses `scripts/figure3/attribution/` |
| **S11–S12** — extended-covariate multi-omics summary and direct candidate associations | Statistics: `scripts/figure5/extended_covariate_sensitivity.R` and `scripts/figure6/extended_covariate_feature_attribution.R`; final tables: `scripts/figure5/build_extended_covariate_tables.R` |
| **S13–S16** — main mediation, binary-HUA, reference-split and medsens analyses | `scripts/figure6/section3_5_step1_data_preparation.R`; `section3_5_step2_mediation_analysis.R` |
| **S17** — extended-covariate mediation | `scripts/figure6/extended_covariate_feature_attribution.R` |
| **S18** — within-Type 4 heterogeneity | `scripts/figure6/section3_5_step3_type4_heterogeneity.R` |

`supplementary/build_supplementary_xlsx.py` reads only `results/tables/final/` and writes the final 29-sheet workbook. **S8h is the shared EC list; S8i is the Type 4-exclusive EC list.** S16 includes the complete `rho_at_sign_change` column; `NA` means no sign reversal within the explored range.

## Statistical computation versus final layout

- `scripts/figure4/assemble_figure4.R` reads authoritative statistical tables and included attribution objects, recomputes the displayed contribution summaries, and applies the second-round layout.
- `scripts/figure5/assemble_figure5.R` recomputes the statistical panels by default from the analysis-ready inputs. Layout-only reuse of caches is optional and is not required for this release.
- `scripts/figure6/assemble_figure6.R` deliberately reads the frozen mediation and heterogeneity tables; it checks their structure and does not repeat the 5,000-replicate bootstrap analyses.
- The final PDFs under `results/figures/` and `supplementary/figures/` are the submitted versions, not historical intermediate renders.
