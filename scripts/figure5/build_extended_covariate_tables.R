# Build final Table S12 and Table S11 source CSVs for the covariate sensitivity analyses.
# Source analyses are immutable inputs; this script only joins and summarizes.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(here)
})

out_dir <- here("results", "tables", "final")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Table S8: candidate attribution plus direct UA/HUA sensitivity
# -----------------------------------------------------------------------------

hier <- read_csv(
  here("results", "tables", "supplementary", "feature_attribution_hierarchy.csv"),
  show_col_types = FALSE
)
direct <- read_csv(
  here("results", "tables", "supplementary", "direct_feature_UA_association.csv"),
  show_col_types = FALSE
)
hua <- read_csv(
  here("results", "r1_5_covariate_sensitivity", "direct_feature_HUA_sensitivity.csv"),
  show_col_types = FALSE
)

ua_base <- direct %>%
  filter(model == "LM_UA_base") %>%
  transmute(
    feature, UA_n_sex_age = n,
    UA_beta_sex_age = UA_beta, UA_FDR_sex_age = UA_padj
  )
ua_bmi_cr <- direct %>%
  filter(model == "LM_UA_extended") %>%
  transmute(
    feature, UA_n_sex_age_BMI_Cr = n,
    UA_beta_sex_age_BMI_Cr = UA_beta,
    UA_FDR_sex_age_BMI_Cr = UA_padj
  )

table_s8 <- hier %>%
  select(
    feature, label, layer, attribution,
    within_enterotype_FDR_supported = layer3_sig,
    n_enterotypes_FDR = layer3_n_et_sig,
    enterotypes_FDR = layer3_et_sig
  ) %>%
  left_join(ua_base, by = "feature") %>%
  left_join(ua_bmi_cr, by = "feature") %>%
  left_join(
    hua %>% transmute(
      feature,
      HUA_n_sex_age_full = original_n,
      HUA_events_sex_age_full = original_events,
      HUA_OR_sex_age_full = original_OR,
      HUA_FDR_sex_age_full = original_padj,
      HUA_n_matched = base_cc_n,
      HUA_events_matched = base_cc_events,
      HUA_OR_sex_age_matched = base_cc_OR,
      HUA_FDR_sex_age_matched = base_cc_padj,
      HUA_OR_sex_age_BMI_Cr = extended_OR,
      HUA_CI_lo_sex_age_BMI_Cr = extended_OR_CI_lo,
      HUA_CI_hi_sex_age_BMI_Cr = extended_OR_CI_hi,
      HUA_FDR_sex_age_BMI_Cr = extended_padj,
      HUA_direction_retained = direction_retained_matched,
      HUA_significance_transition = significance_transition
    ),
    by = "feature"
  ) %>%
  arrange(factor(layer, levels = c("Metabolite", "Pathway", "EC")), feature)

stopifnot(nrow(table_s8) == 22, !anyDuplicated(table_s8$feature))
write_csv(
  table_s8,
  file.path(out_dir, "Table_S12_candidate_feature_associations.csv"),
  na = "NA"
)

# -----------------------------------------------------------------------------
# Table S10: compact reader-facing summary
# -----------------------------------------------------------------------------

da <- read_csv(
  here("results", "r1_5_covariate_sensitivity", "differential_core_stability_summary.csv"),
  show_col_types = FALSE
)
cross <- read_csv(
  here("results", "r1_5_covariate_sensitivity", "crossomics_correlation_summary.csv"),
  show_col_types = FALSE
) %>% filter(target_scope == "crossomics_only")
dbrda <- read_csv(
  here("results", "r1_5_covariate_sensitivity", "crossomics_dbrda_summary.csv"),
  show_col_types = FALSE
)
bmi_et <- read_csv(
  here("results", "r1_5_covariate_sensitivity", "bmi_missingness_by_enterotype.csv"),
  show_col_types = FALSE
)
bmi_tests <- read_csv(
  here("results", "r1_5_covariate_sensitivity", "bmi_complete_case_selection_tests.csv"),
  show_col_types = FALSE
)

pair_n_note <- function(layer, et) {
  if (layer == "Metabolome" && et == "E1") {
    "Original pairwise n=216/166; matched n=166/123."
  } else if (layer == "Metabolome" && et == "E4") {
    "Original pairwise n=137/87; matched n=101/58."
  } else if (et == "E1") {
    "Original pairwise n=218/168; matched n=167/124."
  } else {
    "Original pairwise n=137/87; matched n=101/58."
  }
}

da_rows <- da %>%
  rowwise() %>%
  transmute(
    analysis_component = "Differential features",
    layer_or_model = layer,
    subgroup_or_scope = enterotype,
    original_n = NA_real_, matched_n = NA_real_,
    items_or_original_core_n = original_core_n,
    matched_sex_age_FDR_or_core_n = base_cc_core_n,
    sex_age_BMI_Cr_FDR_or_core_n = extended_core_retained_n,
    direction_concordance = direction_retained_both_rate,
    coefficient_concordance = NA_real_,
    sex_age_BMI_Cr_effect = NA_real_, sex_age_BMI_Cr_p = NA_real_,
    sample_size_note = pair_n_note(layer, enterotype),
    interpretation = paste0(
      "FDR-level original-core retention=",
      ifelse(is.na(extended_core_retained_rate), "NA",
             sprintf("%.1f%%", 100 * extended_core_retained_rate)),
      "; new adjusted core features=", extended_new_core_n, "."
    )
  ) %>% ungroup()

cross_rows <- cross %>% transmute(
  analysis_component = "Cross-omics partial correlations",
  layer_or_model = family,
  subgroup_or_scope = "12 serum metabolites (UA excluded)",
  original_n = n_original, matched_n = n_matched,
  items_or_original_core_n = pairs_n,
  matched_sex_age_FDR_or_core_n = base_cc_FDR_n,
  sex_age_BMI_Cr_FDR_or_core_n = extended_FDR_n,
  direction_concordance = direction_concordance_base_cc_vs_extended,
  coefficient_concordance = rho_spearman_base_cc_vs_extended,
  sex_age_BMI_Cr_effect = NA_real_, sex_age_BMI_Cr_p = NA_real_,
  sample_size_note = "Matched sex-and-age and sex-age-BMI-creatinine models used the same 224 participants.",
  interpretation = paste0(
    "FDR retained/lost/gained=", FDR_retained_n, "/", FDR_lost_n,
    "/", FDR_gained_n, "."
  )
)

dbrda_ext <- dbrda %>% filter(model == "extended_cc")
dbrda_base_cc <- dbrda %>% filter(model == "base_cc")
dbrda_row <- tibble(
  analysis_component = "Global microbiome-metabolome association",
  layer_or_model = "Partial db-RDA",
  subgroup_or_scope = "Metabolome PCs constraining species Bray-Curtis dissimilarity",
  original_n = 303, matched_n = dbrda_ext$n,
  items_or_original_core_n = dbrda_ext$n_metabolome_PCs,
  matched_sex_age_FDR_or_core_n = NA_real_, sex_age_BMI_Cr_FDR_or_core_n = NA_real_,
  direction_concordance = NA_real_, coefficient_concordance = NA_real_,
  sex_age_BMI_Cr_effect = dbrda_ext$adjusted_R2,
  sex_age_BMI_Cr_p = dbrda_ext$permutation_P,
  sample_size_note = paste0(
    "Matched sex-and-age adjusted R2=", sprintf("%.3f", dbrda_base_cc$adjusted_R2),
    ", P=", sprintf("%.3f", dbrda_base_cc$permutation_P), "."
  ),
  interpretation = "The global association remained detectable after BMI and creatinine adjustment."
)

highlight_row <- tibble(
  analysis_component = "Highlighted shared EC sensitivity",
  layer_or_model = "EC 2.7.6.3; EC 2.6.1.1; EC 2.1.1.13",
  subgroup_or_scope = "Dual-contrast Type 1 and Type 4 core criterion",
  original_n = NA_real_, matched_n = NA_real_,
  items_or_original_core_n = 3, matched_sex_age_FDR_or_core_n = NA_real_,
  sex_age_BMI_Cr_FDR_or_core_n = 2, direction_concordance = 1,
  coefficient_concordance = NA_real_, sex_age_BMI_Cr_effect = NA_real_,
  sex_age_BMI_Cr_p = NA_real_, sample_size_note = "See differential-feature rows above.",
  interpretation = paste(
    "EC 2.6.1.1 and EC 2.7.6.3 retained the dual-contrast core criterion;",
    "EC 2.1.1.13 retained positive directions but not Type 1 FDR-level core status."
  )
)

type4 <- bmi_et %>% filter(Enterotype == 4)
bmi_p <- bmi_tests %>% filter(comparison == "BMI completeness by enterotype")
ua_p <- bmi_tests %>% filter(comparison == "UA by BMI completeness")
bmi_row <- tibble(
  analysis_component = "Complete-case selection",
  layer_or_model = "BMI availability",
  subgroup_or_scope = "All participants",
  original_n = 305, matched_n = sum(bmi_et$n_BMI_complete),
  items_or_original_core_n = 80, matched_sex_age_FDR_or_core_n = NA_real_,
  sex_age_BMI_Cr_FDR_or_core_n = NA_real_, direction_concordance = NA_real_,
  coefficient_concordance = NA_real_,
  sex_age_BMI_Cr_effect = type4$BMI_complete_rate,
  sex_age_BMI_Cr_p = bmi_p$pvalue,
  sample_size_note = paste0(
    "BMI missing for 80/305; Type 4 complete=", type4$n_BMI_complete,
    "/", type4$n_total, "."
  ),
  interpretation = paste0(
    "Completeness differed by enterotype; mean UA was 341.6 in missing vs 318.3 in complete participants (P=",
    sprintf("%.4f", ua_p$pvalue), ")."
  )
)

table_s10 <- bind_rows(da_rows, cross_rows, dbrda_row, highlight_row, bmi_row)
stopifnot(nrow(table_s10) == 13)
write_csv(
  table_s10,
  file.path(out_dir, "Table_S11_multiomics_covariate_sensitivity.csv"),
  na = "NA"
)

cat("Table S12 rows:", nrow(table_s8), "columns:", ncol(table_s8), "\n")
cat("Table S11 rows:", nrow(table_s10), "columns:", ncol(table_s10), "\n")
