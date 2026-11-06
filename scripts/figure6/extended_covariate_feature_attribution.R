# ==============================================================================
# Extended Covariate Sensitivity Analysis & Direct Feature-UA Association
#
# Purpose:
#   Respond to reviewer comments R1.5 (inconsistent covariate adjustment) and
#   R1.3 (feature attribution hierarchy).
#
#   Analysis 1 (R1.5): Re-run mediation with extended covariates
#     (Gender + Age + BMI + Cr) and compare to base (Gender + Age).
#   Analysis 2 (R1.3, Layer 2): Direct feature-UA association across enterotypes
#     (enterotype as covariate).
#   Analysis 3 (R1.3, Layer 3): Within-enterotype feature-UA association.
#   Analysis 4: Feature attribution hierarchy summary.
#
# Input:
#   - results/tables/section_3_5/mediation_analysis_data.rds
#   - data/clean/phyloseq_Taxa_with_enterotype_seed10.rds (for BMI, Cr, eGFR)
#   - results/tables/section_3_5/mediation_main_results.csv (original base results)
#
# Output (results/tables/supplementary/):
#   - mediation_extended_covariate_sensitivity.csv  (Analysis 1)
#   - direct_feature_UA_association.csv             (Analysis 2)
#   - within_enterotype_feature_UA.csv              (Analysis 3)
#   - feature_attribution_hierarchy.csv             (Analysis 4)
# ==============================================================================

library(mediation)
library(tidyverse)
library(here)

set.seed(2024)

cat("================================================================\n")
cat("Extended Covariate Sensitivity & Feature Attribution Hierarchy\n")
cat("================================================================\n\n")

# --- 0. Setup ---
out_dir <- here("results", "tables", "supplementary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Mediator definitions (22 candidates: 9 metabolites, 7 pathways, 6 ECs)
mediator_cols <- c(
  "MeOIA", "pCresol", "pCresolSulfate", "IndoxylSulfate", "IPA", "PAGln",
  "DCA", "TMAO", "Adenosine",
  "PWY5695", "PWY0_1296", "PWY5004", "PWY7219",
  "PWY6609", "ARGORNPROST", "PWY7234",
  "EC_2_6_1_1", "EC_2_3_1_129", "EC_4_2_1_8",
  "EC_2_7_6_3", "EC_3_1_3_15", "EC_2_1_1_13"
)

mediator_labels <- c(
  MeOIA = "5-Methoxyindoleacetate", pCresol = "p-Cresol",
  pCresolSulfate = "p-Cresol sulfate", IndoxylSulfate = "Indoxyl sulfate",
  IPA = "Indole-3-propionic acid", PAGln = "Phenylacetylglutamine",
  DCA = "Deoxycholic acid", TMAO = "TMAO", Adenosine = "Adenosine",
  PWY5695 = "PWY-5695 (urate biosynthesis)",
  PWY0_1296 = "PWY0-1296 (purine ribonucleoside degradation)",
  PWY5004 = "PWY-5004 (citrulline metabolism)",
  PWY7219 = "PWY-7219 (adenosine ribonucleotides de novo)",
  PWY6609 = "PWY-6609 (adenine/adenosine salvage)",
  ARGORNPROST = "ARGORNPROST-PWY (Arg/Orn/Pro)",
  PWY7234 = "PWY-7234 (IMP biosynthesis)",
  EC_2_6_1_1 = "EC 2.6.1.1 (aspartate transaminase)",
  EC_2_3_1_129 = "EC 2.3.1.129 (UDP-GlcNAc O-acyltransferase)",
  EC_4_2_1_8 = "EC 4.2.1.8",
  EC_2_7_6_3 = "EC 2.7.6.3 (diphosphokinase)",
  EC_3_1_3_15 = "EC 3.1.3.15 (histidinol-phosphatase)",
  EC_2_1_1_13 = "EC 2.1.1.13 (Met synthase)"
)

mediator_layer <- c(
  rep("Metabolite", 9), rep("Pathway", 7), rep("EC", 6)
)
names(mediator_layer) <- mediator_cols

SIMS <- 5000  # R1.5 sensitivity: match the primary mediation analysis

# Skip Analysis 1 when output already exists unless an explicit verified rerun is requested.
FORCE_MEDIATION <- identical(Sys.getenv("R15_FORCE_MEDIATION", unset = "0"), "1")
SKIP_MEDIATION <- file.exists(
  here("results", "tables", "supplementary",
       "mediation_extended_covariate_sensitivity.csv")
) && !FORCE_MEDIATION
cat(sprintf("Mediation bootstrap iterations: %d; force rerun: %s\n",
            SIMS, FORCE_MEDIATION))

# ==============================================================================
# 1. Load data and merge extended covariates
# ==============================================================================

cat("--- Loading analysis data ---\n")
df <- readRDS(here("results", "tables", "section_3_5", "mediation_analysis_data.rds"))
cat(sprintf("Mediation data: %d samples x %d variables\n", nrow(df), ncol(df)))

cat("--- Loading phyloseq metadata for BMI, Cr, eGFR ---\n")
library(phyloseq)
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)

ext_cov <- meta %>%
  select(SampleID, BMI, Cr, eGFR) %>%
  mutate(
    BMI  = as.numeric(BMI),
    Cr   = as.numeric(Cr),
    eGFR = as.numeric(eGFR)
  )

cat(sprintf("Extended covariates: BMI NA=%d, Cr NA=%d, eGFR NA=%d (of %d)\n",
            sum(is.na(ext_cov$BMI)), sum(is.na(ext_cov$Cr)),
            sum(is.na(ext_cov$eGFR)), nrow(ext_cov)))

# Merge into analysis data frame
df <- df %>%
  left_join(ext_cov, by = "SampleID")

cat(sprintf("Merged data: %d samples\n", nrow(df)))
cat(sprintf("Complete cases (UA + Gender + Age + BMI + Cr): %d\n",
            sum(complete.cases(df[, c("UA", "Gender", "Age", "BMI", "Cr")]))))

# Ensure Enterotype is character for modeling
df$Enterotype <- as.character(df$Enterotype)
df$Gender <- as.character(df$Gender)

cat("\n")

# ==============================================================================
# ANALYSIS 1: Extended Covariate Mediation Sensitivity (R1.5)
# ==============================================================================

cat("================================================================\n")
cat("ANALYSIS 1: Extended Covariate Mediation Sensitivity\n")
cat("================================================================\n\n")

if (SKIP_MEDIATION) {
  cat(">>> Skipping mediation computation (output exists). Loading saved results. <<<\n")
  sensitivity_compare <- read.csv(file.path(
    out_dir, "mediation_extended_covariate_sensitivity.csv"
  ))
  cat(sprintf("Loaded %d comparisons from saved file.\n", nrow(sensitivity_compare)))
} else {

# Core mediation runner (reuses pattern from section3_5_step2)
run_mediation <- function(data, treat_col, mediator_col, outcome_col,
                          covariates, sims = SIMS, outcome_family = "linear") {
  cov_str <- paste(covariates, collapse = " + ")
  med_formula <- as.formula(paste(mediator_col, "~", treat_col, "+", cov_str))
  out_formula <- as.formula(paste(outcome_col, "~", treat_col, "+",
                                  mediator_col, "+", cov_str))

  med_model <- eval(bquote(lm(.(med_formula), data = data)))

  if (outcome_family == "linear") {
    out_model <- eval(bquote(lm(.(out_formula), data = data)))
  } else {
    out_model <- eval(bquote(glm(.(out_formula), data = data,
                                  family = binomial("logit"))))
  }

  med_result <- tryCatch(
    mediate(med_model, out_model,
      treat = treat_col, mediator = mediator_col,
      boot = TRUE, sims = sims, boot.ci.type = "perc"
    ),
    error = function(e) {
      warning(sprintf("Mediation failed for %s: %s", mediator_col, e$message))
      return(NULL)
    }
  )

  if (is.null(med_result)) return(NULL)

  tibble(
    mediator       = mediator_col,
    ACME_estimate  = med_result$d0,
    ACME_CI_lo     = med_result$d0.ci[1],
    ACME_CI_hi     = med_result$d0.ci[2],
    ACME_p         = med_result$d0.p,
    ADE_estimate   = med_result$z0,
    ADE_p          = med_result$z0.p,
    Total_estimate = med_result$tau.coef,
    Total_p        = med_result$tau.p,
    Prop_mediated  = med_result$n0,
    Prop_p         = med_result$n0.p
  )
}

# Run mediation for one model across all mediators under a covariate set
run_model_set <- function(model_name, treat_col, covariates, data) {
  sub <- data[!is.na(data[[treat_col]]), ]
  cov_str <- paste(covariates, collapse = "+")
  cat(sprintf("--- %s | covariates: %s | n=%d ---\n",
              model_name, cov_str, nrow(sub)))

  results_list <- list()
  for (m in mediator_cols) {
    if (all(is.na(sub[[m]]))) {
      cat(sprintf("  %-18s: SKIPPED (all NA)\n", m))
      next
    }
    res <- run_mediation(sub, treat_col, m, "UA", covariates, sims = SIMS)
    if (!is.null(res)) {
      res$model <- model_name
      res$covariates <- cov_str
      res$label <- mediator_labels[m]
      results_list[[m]] <- res
      cat(sprintf("  %-18s: ACME=%.3f [%.3f, %.3f] p=%.4f\n",
                  m, res$ACME_estimate, res$ACME_CI_lo,
                  res$ACME_CI_hi, res$ACME_p))
    } else {
      cat(sprintf("  %-18s: FAILED\n", m))
    }
  }
  bind_rows(results_list)
}

# Use complete cases for BOTH base and extended so comparison is fair
# (same samples, only covariate set differs)
df_cc <- df[complete.cases(df[, c("UA", "Gender", "Age", "BMI", "Cr")]), ]
cat(sprintf("Complete-case subset for sensitivity: n=%d\n", nrow(df_cc)))
cat(sprintf("  Type1 model: T=%d, Ref=%d\n",
            sum(df_cc$Treat_T1 == 1, na.rm = TRUE),
            sum(df_cc$Treat_T1 == 0, na.rm = TRUE)))
cat(sprintf("  Type4 model: T=%d, Ref=%d\n\n",
            sum(df_cc$Treat_T4 == 1, na.rm = TRUE),
            sum(df_cc$Treat_T4 == 0, na.rm = TRUE)))

base_covs  <- c("Gender", "Age")
ext_covs   <- c("Gender", "Age", "BMI", "Cr")

cat(">>> Running BASE adjustment (Gender + Age) <<<\n")
base_t1  <- run_model_set("Type1_vs_Ref",  "Treat_T1", base_covs, df_cc)
base_t4  <- run_model_set("Type4_vs_Ref",  "Treat_T4", base_covs, df_cc)

cat("\n>>> Running EXTENDED adjustment (Gender + Age + BMI + Cr) <<<\n")
ext_t1   <- run_model_set("Type1_vs_Ref",  "Treat_T1", ext_covs,  df_cc)
ext_t4   <- run_model_set("Type4_vs_Ref",  "Treat_T4", ext_covs,  df_cc)

base_all <- bind_rows(base_t1, base_t4) %>% mutate(adjustment = "Base (Gender+Age)")
ext_all  <- bind_rows(ext_t1,  ext_t4)  %>% mutate(adjustment = "Extended (Gender+Age+BMI+Cr)")

# BH correction within each model x adjustment combination
base_all <- base_all %>%
  group_by(model) %>%
  mutate(ACME_padj = p.adjust(ACME_p, method = "BH")) %>%
  ungroup()

ext_all <- ext_all %>%
  group_by(model) %>%
  mutate(ACME_padj = p.adjust(ACME_p, method = "BH")) %>%
  ungroup()

# Build sensitivity comparison table
sensitivity_compare <- base_all %>%
  select(mediator, label, model,
         base_ACME = ACME_estimate,
         base_CI_lo = ACME_CI_lo, base_CI_hi = ACME_CI_hi,
         base_p = ACME_p, base_padj = ACME_padj) %>%
  left_join(
    ext_all %>%
      select(mediator, model,
             ext_ACME = ACME_estimate,
             ext_CI_lo = ACME_CI_lo, ext_CI_hi = ACME_CI_hi,
             ext_p = ACME_p, ext_padj = ACME_padj),
    by = c("mediator", "model")
  ) %>%
  mutate(
    direction_consistent = sign(base_ACME) == sign(ext_ACME),
    base_sig    = base_padj < 0.05,
    ext_sig     = ext_padj  < 0.05,
    sig_stable  = base_sig == ext_sig,
    ACME_change = ext_ACME - base_ACME,
    ACME_pct_change = ifelse(abs(base_ACME) > 1e-10,
                             100 * (ext_ACME - base_ACME) / abs(base_ACME), NA),
    stability_summary = case_when(
      base_sig & ext_sig  ~ "Stable significant",
      base_sig & !ext_sig ~ "Lost significance",
      !base_sig & ext_sig ~ "Gained significance",
      TRUE                ~ "Stable non-significant"
    )
  )

write.csv(sensitivity_compare,
          file.path(out_dir, "mediation_extended_covariate_sensitivity.csv"),
          row.names = FALSE)
cat("\n>>> Analysis 1 results saved: mediation_extended_covariate_sensitivity.csv\n")

cat("\n--- Sensitivity Summary ---\n")
cat(sprintf("Total comparisons: %d\n", nrow(sensitivity_compare)))
cat(sprintf("Direction consistent: %d / %d\n",
            sum(sensitivity_compare$direction_consistent, na.rm = TRUE),
            nrow(sensitivity_compare)))
cat(sprintf("Significance stable: %d / %d\n",
            sum(sensitivity_compare$sig_stable, na.rm = TRUE),
            nrow(sensitivity_compare)))
cat("\nStability breakdown:\n")
print(table(sensitivity_compare$stability_summary))

cat("\nFeatures significant under BOTH adjustments:\n")
both_sig <- sensitivity_compare %>%
  filter(base_sig & ext_sig) %>%
  select(mediator, label, model, base_ACME, base_padj, ext_ACME, ext_padj)
if (nrow(both_sig) > 0) print(both_sig, n = Inf) else cat("  None\n")

cat("\nFeatures losing significance under extended adjustment:\n")
lost <- sensitivity_compare %>%
  filter(base_sig & !ext_sig) %>%
  select(mediator, label, model, base_ACME, base_padj, ext_ACME, ext_padj)
if (nrow(lost) > 0) print(lost, n = Inf) else cat("  None\n")

# Save combined full results for reference
mediation_sensitivity_full <- bind_rows(base_all, ext_all)
saveRDS(mediation_sensitivity_full,
        file.path(out_dir, "mediation_extended_covariate_full.rds"))

} # end of if (!SKIP_MEDIATION) else block

cat("\n")

# ==============================================================================
# ANALYSIS 2: Direct Feature-UA Association (Layer 2 — R1.3)
# ==============================================================================

cat("================================================================\n")
cat("ANALYSIS 2: Direct Feature-UA Association (Layer 2)\n")
cat("================================================================\n\n")

# Layer 2: feature ~ UA + Enterotype + covariates
# Tests whether UA directly associates with feature abundance,
# independent of enterotype membership.

run_direct_lm <- function(feature, data, covariates, outcome_var = "UA") {
  cov_str <- paste(covariates, collapse = " + ")
  formula_str <- sprintf("%s ~ %s + Enterotype + %s", feature, outcome_var, cov_str)
  f <- as.formula(formula_str)

  sub <- data[complete.cases(data[, c(feature, outcome_var,
                                       "Enterotype", covariates)]), ]
  if (nrow(sub) < 10) return(NULL)

  fit <- tryCatch(lm(f, data = sub), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  # Coefficients for UA and enterotype levels
  coefs <- summary(fit)$coefficients

  # UA coefficient
  ua_row <- grep("^UA$", rownames(coefs))
  if (length(ua_row) == 0) return(NULL)
  ua_row <- ua_row[1]

  # Enterotype overall p (F-test via anova comparing with/without Enterotype)
  fit_no_et <- tryCatch(
    update(fit, as.formula(paste(". ~ . - Enterotype"))),
    error = function(e) NULL
  )
  et_p <- NA
  if (!is.null(fit_no_et)) {
    aov <- anova(fit_no_et, fit)
    et_p <- aov$`Pr(>F)`[2]
  }

  tibble(
    feature       = feature,
    label         = mediator_labels[feature],
    layer         = mediator_layer[feature],
    n             = nrow(sub),
    UA_beta       = coefs[ua_row, "Estimate"],
    UA_SE         = coefs[ua_row, "Std. Error"],
    UA_t          = coefs[ua_row, "t value"],
    UA_p          = coefs[ua_row, "Pr(>|t|)"],
    Enterotype_p  = et_p,
    adjustment    = paste(covariates, collapse = "+")
  )
}

run_direct_logistic <- function(feature, data, covariates) {
  # HUA as binary outcome, feature as predictor (logistic regression)
  cov_str <- paste(covariates, collapse = " + ")
  formula_str <- sprintf("HUA ~ %s + Enterotype + %s", feature, cov_str)
  f <- as.formula(formula_str)

  sub <- data[complete.cases(data[, c(feature, "HUA", "Enterotype", covariates)]), ]
  if (nrow(sub) < 10) return(NULL)

  # Ensure HUA is numeric 0/1 for binomial GLM
  sub$HUA <- as.numeric(sub$HUA)

  fit <- tryCatch(glm(f, data = sub, family = binomial("logit")),
                  error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  coefs <- summary(fit)$coefficients
  # Extract the feature's coefficient (not HUA's)
  feat_row <- grep(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", feature), "$"),
                   rownames(coefs))
  if (length(feat_row) == 0) return(NULL)
  feat_row <- feat_row[1]

  tibble(
    feature       = feature,
    label         = mediator_labels[feature],
    layer         = mediator_layer[feature],
    n             = nrow(sub),
    HUA_beta      = coefs[feat_row, "Estimate"],
    HUA_SE        = coefs[feat_row, "Std. Error"],
    HUA_z         = coefs[feat_row, "z value"],
    HUA_p         = coefs[feat_row, "Pr(>|z|)"],
    adjustment    = paste(covariates, collapse = "+")
  )
}

# 2a. Base adjustment: feature ~ UA + Enterotype + Gender + Age
cat("--- 2a: Direct UA association, base adjustment (Gender + Age) ---\n")
direct_base <- map_dfr(mediator_cols, function(m) {
  if (all(is.na(df[[m]]))) return(NULL)
  run_direct_lm(m, df, c("Gender", "Age"))
})
direct_base$UA_padj <- p.adjust(direct_base$UA_p, method = "BH")
direct_base$ET_padj <- p.adjust(direct_base$Enterotype_p, method = "BH")
direct_base$model <- "LM_UA_base"
cat(sprintf("  Tested %d features\n", nrow(direct_base)))
print(direct_base %>%
        select(feature, label, n, UA_beta, UA_p, UA_padj, Enterotype_p, ET_padj) %>%
        arrange(UA_padj), n = Inf)

# 2b. Extended adjustment: feature ~ UA + Enterotype + Gender + Age + BMI + Cr
cat("\n--- 2b: Direct UA association, extended adjustment (Gender + Age + BMI + Cr) ---\n")
direct_ext <- map_dfr(mediator_cols, function(m) {
  if (all(is.na(df[[m]]))) return(NULL)
  run_direct_lm(m, df, c("Gender", "Age", "BMI", "Cr"))
})
direct_ext$UA_padj <- p.adjust(direct_ext$UA_p, method = "BH")
direct_ext$ET_padj <- p.adjust(direct_ext$Enterotype_p, method = "BH")
direct_ext$model <- "LM_UA_extended"
cat(sprintf("  Tested %d features\n", nrow(direct_ext)))
print(direct_ext %>%
        select(feature, label, n, UA_beta, UA_p, UA_padj, Enterotype_p, ET_padj) %>%
        arrange(UA_padj), n = Inf)

# 2c. HUA binary: logistic regression HUA ~ feature + Enterotype + Gender + Age
cat("\n--- 2c: Direct HUA association (logistic, Gender + Age) ---\n")
direct_hua <- map_dfr(mediator_cols, function(m) {
  if (all(is.na(df[[m]]))) return(NULL)
  run_direct_logistic(m, df, c("Gender", "Age"))
})
if (nrow(direct_hua) > 0) {
  direct_hua$HUA_padj <- p.adjust(direct_hua$HUA_p, method = "BH")
  direct_hua$model <- "LOGISTIC_HUA_base"
}
cat(sprintf("  Tested %d features\n", nrow(direct_hua)))
if (nrow(direct_hua) > 0) {
  print(direct_hua %>%
          select(feature, label, n, HUA_beta, HUA_p, HUA_padj) %>%
          arrange(HUA_padj), n = Inf)
} else {
  cat("  No features tested (all failed)\n")
}

# Combine into one output table
direct_all_parts <- list(
  direct_base %>% select(feature, label, layer, model, n, adjustment,
                         UA_beta, UA_SE, UA_p, UA_padj, Enterotype_p, ET_padj),
  direct_ext  %>% select(feature, label, layer, model, n, adjustment,
                         UA_beta, UA_SE, UA_p, UA_padj, Enterotype_p, ET_padj)
)
if (nrow(direct_hua) > 0) {
  direct_all_parts[[3]] <- direct_hua %>%
    select(feature, label, layer, model, n, adjustment,
           HUA_beta, HUA_SE, HUA_p, HUA_padj) %>%
    mutate(UA_beta = NA, UA_SE = NA, UA_p = NA, UA_padj = NA,
           Enterotype_p = NA, ET_padj = NA)
}
direct_all <- bind_rows(direct_all_parts)

write.csv(direct_all,
          file.path(out_dir, "direct_feature_UA_association.csv"),
          row.names = FALSE)
cat("\n>>> Analysis 2 results saved: direct_feature_UA_association.csv\n")

# Summary
cat("\n--- Layer 2 Summary ---\n")
cat("Base adjustment (UA significant at FDR<0.05):\n")
sig_base <- direct_base %>% filter(UA_padj < 0.05)
if (nrow(sig_base) > 0) {
  print(sig_base %>% select(feature, label, UA_beta, UA_padj))
} else cat("  None\n")

cat("\nExtended adjustment (UA significant at FDR<0.05):\n")
sig_ext <- direct_ext %>% filter(UA_padj < 0.05)
if (nrow(sig_ext) > 0) {
  print(sig_ext %>% select(feature, label, UA_beta, UA_padj))
} else cat("  None\n")

cat("\nHUA logistic (significant at FDR<0.05):\n")
if (nrow(direct_hua) > 0 && "HUA_padj" %in% colnames(direct_hua)) {
  sig_hua <- direct_hua %>% filter(HUA_padj < 0.05)
  if (nrow(sig_hua) > 0) {
    print(sig_hua %>% select(feature, label, HUA_beta, HUA_padj))
  } else cat("  None\n")
} else cat("  None (no HUA results)\n")

cat("\n")

# ==============================================================================
# ANALYSIS 3: Within-Enterotype Feature-UA Association (Layer 3 — R1.3)
# ==============================================================================

cat("================================================================\n")
cat("ANALYSIS 3: Within-Enterotype Feature-UA Association (Layer 3)\n")
cat("================================================================\n\n")

enterotypes <- c("1", "2", "3", "4")

within_results <- list()
for (et in enterotypes) {
  sub_et <- df[df$Enterotype == et, ]
  cat(sprintf("--- Enterotype %s: n=%d ---\n", et, nrow(sub_et)))

  for (m in mediator_cols) {
    if (all(is.na(sub_et[[m]]))) next

    formula_str <- sprintf("%s ~ UA + Gender + Age", m)
    f <- as.formula(formula_str)
    sub_cc <- sub_et[complete.cases(sub_et[, c(m, "UA", "Gender", "Age")]), ]

    if (nrow(sub_cc) < 10) next

    fit <- tryCatch(lm(f, data = sub_cc), error = function(e) NULL)
    if (is.null(fit)) next

    coefs <- summary(fit)$coefficients
    ua_row <- grep("^UA$", rownames(coefs))
    if (length(ua_row) == 0) next

    within_results[[paste(et, m, sep = "_")]] <- tibble(
      enterotype = et,
      feature    = m,
      label      = mediator_labels[m],
      layer      = mediator_layer[m],
      n          = nrow(sub_cc),
      UA_beta    = coefs[ua_row[1], "Estimate"],
      UA_SE      = coefs[ua_row[1], "Std. Error"],
      UA_t       = coefs[ua_row[1], "t value"],
      UA_p       = coefs[ua_row[1], "Pr(>|t|)"]
    )
  }
}

within_df <- bind_rows(within_results)
within_df$UA_padj <- p.adjust(within_df$UA_p, method = "BH")

write.csv(within_df,
          file.path(out_dir, "within_enterotype_feature_UA.csv"),
          row.names = FALSE)
cat("\n>>> Analysis 3 results saved: within_enterotype_feature_UA.csv\n")

# Summary: which features significant in each enterotype
cat("\n--- Layer 3: Within-enterotype UA associations (FDR<0.05) ---\n")
within_sig <- within_df %>% filter(UA_padj < 0.05)
if (nrow(within_sig) > 0) {
  print(within_sig %>% select(enterotype, feature, label, n, UA_beta, UA_p, UA_padj))
} else {
  cat("  None at FDR<0.05\n")
  cat("\n  Nominal (p<0.05):\n")
  within_nom <- within_df %>% filter(UA_p < 0.05)
  if (nrow(within_nom) > 0) {
    print(within_nom %>% select(enterotype, feature, label, UA_beta, UA_p, UA_padj))
  } else cat("  None\n")
}

# Feature x enterotype beta(p) matrix
cat("\n--- Feature x Enterotype beta(p) matrix ---\n")
beta_matrix <- within_df %>%
  select(feature, enterotype, UA_beta, UA_p) %>%
  pivot_wider(names_from = enterotype, values_from = c(UA_beta, UA_p),
              names_sep = "_ET")
beta_matrix <- beta_matrix %>%
  mutate(label = mediator_labels[feature]) %>%
  select(feature, label, starts_with("UA_beta_ET"), starts_with("UA_p_ET"))
print(beta_matrix, n = Inf)

# Count how many enterotypes each feature is significant in (nominal p<0.05)
robust_count <- within_df %>%
  filter(UA_p < 0.05) %>%
  group_by(feature) %>%
  summarize(n_enterotypes_sig = n(), .groups = "drop")

cat("\n--- Robustness: features significant (p<0.05) in multiple enterotypes ---\n")
robust_count <- robust_count %>%
  arrange(desc(n_enterotypes_sig)) %>%
  mutate(label = mediator_labels[feature])
if (nrow(robust_count) > 0) {
  print(robust_count, n = Inf)
} else cat("  None\n")

cat("\n")

# ==============================================================================
# ANALYSIS 4: Feature Attribution Hierarchy Summary
# ==============================================================================

cat("================================================================\n")
cat("ANALYSIS 4: Feature Attribution Hierarchy Summary\n")
cat("================================================================\n\n")

# Classification logic:
#   Layer 2 significant: direct_base UA_padj < 0.05 (feature ~ UA + ET + cov)
#   Layer 3 significant: within_df UA_p < 0.05 in >=1 enterotype
#   Enterotype-associated: direct_base ET_padj < 0.05
#
# Categories:
#   "UA-associated candidate" — Layer 2 sig AND Layer 3 sig (>=1 ET)
#   "UA-associated (Layer 2 only)" — Layer 2 sig, Layer 3 not sig
#   "UA-associated (Layer 3 only)" — Layer 3 sig, Layer 2 not sig
#   "enterotype-associated marker only" — ET sig, Layer 2/3 not sig
#   "both UA & enterotype" — UA-associated (any layer) AND ET sig
#   "neither" — none significant

layer2_sig <- direct_base %>%
  select(feature, layer2_padj = UA_padj, layer2_beta = UA_beta,
         layer2_p = UA_p) %>%
  mutate(layer2_sig = layer2_padj < 0.05)

layer2_ext_sig <- direct_ext %>%
  select(feature, layer2ext_padj = UA_padj, layer2ext_beta = UA_beta) %>%
  mutate(layer2ext_sig = layer2ext_padj < 0.05)

layer3_sig <- within_df %>%
  filter(UA_p < 0.05) %>%
  group_by(feature) %>%
  summarize(
    layer3_sig = TRUE,
    layer3_n_et_sig = n(),
    layer3_et_sig = paste(enterotype, collapse = ","),
    layer3_min_p = min(UA_p),
    .groups = "drop"
  )

et_sig <- direct_base %>%
  select(feature, et_padj = ET_padj, et_p = Enterotype_p) %>%
  mutate(et_sig = et_padj < 0.05)

hua_sig <- if (nrow(direct_hua) > 0 && "HUA_padj" %in% colnames(direct_hua)) {
  direct_hua %>%
    select(feature, hua_padj = HUA_padj, hua_beta = HUA_beta) %>%
    mutate(hua_sig = hua_padj < 0.05)
} else {
  tibble(feature = mediator_cols, hua_padj = NA_real_,
         hua_beta = NA_real_, hua_sig = FALSE)
}

hierarchy <- tibble(feature = mediator_cols) %>%
  mutate(
    label = mediator_labels[feature],
    layer = mediator_layer[feature]
  ) %>%
  left_join(layer2_sig,  by = "feature") %>%
  left_join(layer2_ext_sig, by = "feature") %>%
  left_join(layer3_sig,  by = "feature") %>%
  left_join(et_sig,      by = "feature") %>%
  left_join(hua_sig,     by = "feature") %>%
  mutate(
    layer2_sig = replace_na(layer2_sig, FALSE),
    layer2ext_sig = replace_na(layer2ext_sig, FALSE),
    layer3_sig = replace_na(layer3_sig, FALSE),
    et_sig = replace_na(et_sig, FALSE),
    hua_sig = replace_na(hua_sig, FALSE),
    layer3_n_et_sig = replace_na(layer3_n_et_sig, 0L),
    layer3_et_sig = replace_na(layer3_et_sig, "none"),
    layer3_min_p = replace_na(layer3_min_p, NA_real_),
    ua_associated = layer2_sig | layer3_sig
  ) %>%
  mutate(
    attribution = case_when(
      ua_associated & et_sig  ~ "Both UA-associated & enterotype marker",
      layer2_sig & layer3_sig ~ "UA-associated candidate",
      layer2_sig & !layer3_sig ~ "UA-associated (Layer 2 only)",
      !layer2_sig & layer3_sig ~ "UA-associated (Layer 3 only)",
      et_sig ~ "Enterotype-associated marker only",
      TRUE ~ "Neither"
    )
  ) %>%
  arrange(layer, desc(ua_associated), desc(et_sig), layer2_padj)

write.csv(hierarchy,
          file.path(out_dir, "feature_attribution_hierarchy.csv"),
          row.names = FALSE)
cat(">>> Analysis 4 results saved: feature_attribution_hierarchy.csv\n")

cat("\n--- Attribution Hierarchy Summary ---\n")
cat("Classification counts:\n")
print(table(hierarchy$attribution))

cat("\nFull hierarchy table:\n")
print(hierarchy %>%
        select(feature, label, layer, attribution,
               layer2_sig, layer2_p, layer3_sig, layer3_n_et_sig,
               et_sig, et_p, hua_sig, hua_padj),
      n = Inf)

# Key metrics for reviewer response
cat("\n--- Key Metrics for Reviewer Response ---\n")
n_ua_candidate <- sum(hierarchy$attribution == "UA-associated candidate")
n_ua_both <- sum(hierarchy$attribution == "Both UA-associated & enterotype marker")
n_ua_l2only <- sum(hierarchy$attribution == "UA-associated (Layer 2 only)")
n_ua_l3only <- sum(hierarchy$attribution == "UA-associated (Layer 3 only)")
n_et_only <- sum(hierarchy$attribution == "Enterotype-associated marker only")
n_neither <- sum(hierarchy$attribution == "Neither")

cat(sprintf("Total features: %d\n", nrow(hierarchy)))
cat(sprintf("UA-associated candidate (Layer 2 + Layer 3): %d\n", n_ua_candidate))
cat(sprintf("Both UA & enterotype: %d\n", n_ua_both))
cat(sprintf("UA-associated (Layer 2 only): %d\n", n_ua_l2only))
cat(sprintf("UA-associated (Layer 3 only): %d\n", n_ua_l3only))
cat(sprintf("Enterotype-associated marker only: %d\n", n_et_only))
cat(sprintf("Neither: %d\n", n_neither))

cat(sprintf("\nFeatures upgraded from 'enterotype marker' to 'UA-associated': %d\n",
            n_ua_candidate + n_ua_both + n_ua_l2only + n_ua_l3only))

# Mediation sensitivity summary
cat("\n--- Mediation Sensitivity (Analysis 1) Key Metrics ---\n")
cat(sprintf("Comparisons: %d\n", nrow(sensitivity_compare)))
cat(sprintf("Direction consistent: %d / %d (%.1f%%)\n",
            sum(sensitivity_compare$direction_consistent, na.rm = TRUE),
            nrow(sensitivity_compare),
            100 * sum(sensitivity_compare$direction_consistent, na.rm = TRUE) /
              nrow(sensitivity_compare)))
cat(sprintf("Significance stable: %d / %d (%.1f%%)\n",
            sum(sensitivity_compare$sig_stable, na.rm = TRUE),
            nrow(sensitivity_compare),
            100 * sum(sensitivity_compare$sig_stable, na.rm = TRUE) /
              nrow(sensitivity_compare)))
cat(sprintf("Significant under base: %d\n", sum(sensitivity_compare$base_sig, na.rm = TRUE)))
cat(sprintf("Significant under extended: %d\n", sum(sensitivity_compare$ext_sig, na.rm = TRUE)))
cat(sprintf("Lost significance: %d\n",
            sum(sensitivity_compare$stability_summary == "Lost significance", na.rm = TRUE)))
cat(sprintf("Gained significance: %d\n",
            sum(sensitivity_compare$stability_summary == "Gained significance", na.rm = TRUE)))

cat("\n================================================================\n")
cat("All analyses complete.\n")
cat("================================================================\n")

# Save session info for reproducibility
cat("\n--- Session Info ---\n")
print(sessionInfo())
