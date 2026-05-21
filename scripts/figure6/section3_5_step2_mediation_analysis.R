# ==============================================================================
# Section 3.5 — Step 2: Causal Mediation Analysis (Revised — 22 candidates)
#
# Purpose:
#   Test whether specific metabolites/pathways/ECs mediate the enterotype→UA
#   association, separately for Type 1 and Type 4 models.
#
# Input:
#   - results/tables/section_3_5/mediation_analysis_data.rds (from Step 1)
#
# Output (all in results/tables/section_3_5/):
#   - mediation_main_results.csv         Main analysis (Table S-X)
#   - mediation_binary_sensitivity.csv   HUA binary outcome (Table S-Y)
#   - mediation_medsens_results.csv      Sequential ignorability (Table S-Z)
#   - mediation_refgroup_sensitivity.csv Reference group sensitivity (Table S-V)
#   - mediation_results_all.rds          Full R objects for downstream use
#
# Methods:
#   - Baron-Kenny framework with R `mediation` package (Imai et al. 2010)
#   - Non-parametric bootstrap, sims = 5000
#   - BH correction within each model (Type 1 / Type 4)
# ==============================================================================

library(mediation)
library(tidyverse)
library(here)

cat("========================================\n")
cat("Section 3.5 Step 2: Mediation Analysis (Revised — 22 candidates)\n")
cat("========================================\n\n")

set.seed(2024) # For reproducibility of bootstrap

# --- 0. Setup ---
out_dir <- here("results", "tables", "section_3_5")
df <- readRDS(file.path(out_dir, "mediation_analysis_data.rds"))

mediator_cols <- c(
  # Metabolites (9)
  "MeOIA", "pCresol", "pCresolSulfate", "IndoxylSulfate", "IPA", "PAGln",
  "DCA", "TMAO", "Adenosine",
  # Pathways (7)
  "PWY5695", "PWY0_1296", "PWY5004", "PWY7219",
  "PWY6609", "ARGORNPROST", "PWY7234",
  # Convergent shared ECs (6)
  "EC_2_6_1_1", "EC_2_3_1_129", "EC_4_2_1_8",
  "EC_2_7_6_3", "EC_3_1_3_15", "EC_2_1_1_13"
)

mediator_labels <- c(
  MeOIA          = "5-Methoxyindoleacetate",
  pCresol        = "p-Cresol",
  pCresolSulfate = "p-Cresol sulfate",
  IndoxylSulfate = "Indoxyl sulfate",
  IPA            = "Indole-3-propionic acid",
  PAGln          = "Phenylacetylglutamine",
  DCA            = "Deoxycholic acid",
  TMAO           = "TMAO",
  Adenosine      = "Adenosine",
  PWY5695        = "PWY-5695 (urate biosynthesis)",
  PWY0_1296      = "PWY0-1296 (purine ribonucleoside degradation)",
  PWY5004        = "PWY-5004 (citrulline metabolism)",
  PWY7219        = "PWY-7219 (adenosine ribonucleotides de novo)",
  PWY6609        = "PWY-6609 (adenine/adenosine salvage)",
  ARGORNPROST    = "ARGORNPROST-PWY (Arg/Orn/Pro)",
  PWY7234        = "PWY-7234 (IMP biosynthesis)",
  EC_2_6_1_1     = "EC 2.6.1.1 (aspartate transaminase)",
  EC_2_3_1_129   = "EC 2.3.1.129 (UDP-GlcNAc O-acyltransferase)",
  EC_4_2_1_8     = "EC 4.2.1.8",
  EC_2_7_6_3     = "EC 2.7.6.3 (diphosphokinase)",
  EC_3_1_3_15    = "EC 3.1.3.15 (histidinol-phosphatase)",
  EC_2_1_1_13    = "EC 2.1.1.13 (Met synthase)"
)

cat(sprintf("Testing %d candidate mediators\n\n", length(mediator_cols)))

# ==============================================================================
# 1. Core mediation function
# ==============================================================================

run_mediation <- function(data, treat_col, mediator_col, outcome_col,
                          covariates = c("Gender", "Age"),
                          sims = 5000, outcome_family = "linear") {
  cov_str <- paste(covariates, collapse = " + ")
  med_formula <- as.formula(paste(mediator_col, "~", treat_col, "+", cov_str))
  out_formula <- as.formula(paste(outcome_col, "~", treat_col, "+",
                                  mediator_col, "+", cov_str))

  med_model <- eval(bquote(lm(.(med_formula), data = data)))

  if (outcome_family == "linear") {
    out_model <- eval(bquote(lm(.(out_formula), data = data)))
  } else {
    out_model <- eval(bquote(glm(.(out_formula), data = data, family = binomial("logit"))))
  }

  med_result <- tryCatch(
    {
      mediate(med_model, out_model,
        treat = treat_col, mediator = mediator_col,
        boot = TRUE, sims = sims, boot.ci.type = "perc"
      )
    },
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
    ADE_CI_lo      = med_result$z0.ci[1],
    ADE_CI_hi      = med_result$z0.ci[2],
    ADE_p          = med_result$z0.p,
    Total_estimate = med_result$tau.coef,
    Total_CI_lo    = med_result$tau.ci[1],
    Total_CI_hi    = med_result$tau.ci[2],
    Total_p        = med_result$tau.p,
    Prop_mediated  = med_result$n0,
    Prop_CI_lo     = med_result$n0.ci[1],
    Prop_CI_hi     = med_result$n0.ci[2],
    Prop_p         = med_result$n0.p
  )
}

# ==============================================================================
# 2. Main analysis: continuous UA outcome
# ==============================================================================

cat("=== Main Analysis: Continuous UA ===\n\n")

mediate_objects <- list()

run_model <- function(model_name, treat_col) {
  sub <- df[!is.na(df[[treat_col]]), ]
  cat(sprintf(
    "--- %s (n=%d: Treatment=%d, Reference=%d) ---\n",
    model_name, nrow(sub),
    sum(sub[[treat_col]] == 1), sum(sub[[treat_col]] == 0)
  ))

  results_list <- list()
  objects_list <- list()

  for (m in mediator_cols) {
    cat(sprintf("  %s ... ", mediator_labels[m]))

    if (all(is.na(sub[[m]]))) {
      cat("SKIPPED (all NA)\n")
      next
    }

    res <- run_mediation(sub, treat_col, m, "UA", sims = 5000)
    if (!is.null(res)) {
      res$model <- model_name
      res$label <- mediator_labels[m]
      results_list[[m]] <- res
      cat(sprintf(
        "ACME=%.2f [%.2f, %.2f] P=%.4f\n",
        res$ACME_estimate, res$ACME_CI_lo, res$ACME_CI_hi, res$ACME_p
      ))

      # Store raw mediate object for medsens
      cov_str <- "Gender + Age"
      m_formula <- as.formula(paste(m, "~", treat_col, "+", cov_str))
      out_formula <- as.formula(paste("UA ~", treat_col, "+", m, "+", cov_str))
      med_mod <- eval(bquote(lm(.(m_formula), data = sub)))
      out_mod <- eval(bquote(lm(.(out_formula), data = sub)))
      objects_list[[m]] <- mediate(med_mod, out_mod,
        treat = treat_col, mediator = m,
        boot = TRUE, sims = 50,
        boot.ci.type = "perc"
      )
    } else {
      cat("FAILED\n")
    }
  }

  mediate_objects[[model_name]] <<- objects_list
  bind_rows(results_list)
}

res_t1 <- run_model("Type1_vs_Ref", "Treat_T1")
cat("\n")
res_t4 <- run_model("Type4_vs_Ref", "Treat_T4")

# BH correction within each model
res_t1$ACME_padj <- p.adjust(res_t1$ACME_p, method = "BH")
res_t4$ACME_padj <- p.adjust(res_t4$ACME_p, method = "BH")

main_results <- bind_rows(res_t1, res_t4)

cat("\n=== Main Results Summary ===\n")
main_results %>%
  select(model, label, ACME_estimate, ACME_p, ACME_padj, Prop_mediated) %>%
  arrange(model, ACME_padj) %>%
  print(n = Inf)

write.csv(main_results, file.path(out_dir, "mediation_main_results.csv"),
  row.names = FALSE
)
cat("\nMain results saved.\n")

# ==============================================================================
# 3. Sensitivity Analysis 1: Binary HUA outcome
# ==============================================================================

cat("\n=== Sensitivity Analysis 1: Binary HUA Outcome ===\n\n")

run_binary_model <- function(model_name, treat_col) {
  sub <- df[!is.na(df[[treat_col]]) & !is.na(df$HUA), ]
  cat(sprintf(
    "--- %s (n=%d, HUA+=%d) ---\n",
    model_name, nrow(sub), sum(sub$HUA == 1)
  ))

  results_list <- list()
  for (m in mediator_cols) {
    cat(sprintf("  %s ... ", mediator_labels[m]))
    if (all(is.na(sub[[m]]))) {
      cat("SKIPPED\n")
      next
    }

    res <- run_mediation(sub, treat_col, m, "HUA",
      sims = 5000, outcome_family = "logistic"
    )
    if (!is.null(res)) {
      res$model <- model_name
      res$label <- mediator_labels[m]
      results_list[[m]] <- res
      cat(sprintf("ACME=%.4f P=%.4f\n", res$ACME_estimate, res$ACME_p))
    } else {
      cat("FAILED\n")
    }
  }
  bind_rows(results_list)
}

bin_t1 <- run_binary_model("Type1_vs_Ref_HUA", "Treat_T1")
cat("\n")
bin_t4 <- run_binary_model("Type4_vs_Ref_HUA", "Treat_T4")

bin_t1$ACME_padj <- p.adjust(bin_t1$ACME_p, method = "BH")
bin_t4$ACME_padj <- p.adjust(bin_t4$ACME_p, method = "BH")

binary_results <- bind_rows(bin_t1, bin_t4)

write.csv(binary_results, file.path(out_dir, "mediation_binary_sensitivity.csv"),
  row.names = FALSE
)
cat("\nBinary HUA sensitivity results saved.\n")

# ==============================================================================
# 4. Sensitivity Analysis 2: medsens (sequential ignorability)
# ==============================================================================

cat("\n=== Sensitivity Analysis 2: Sequential Ignorability (medsens) ===\n\n")

sig_mediators <- main_results %>%
  filter(ACME_padj < 0.05) %>%
  select(model, mediator)

if (nrow(sig_mediators) == 0) {
  cat("No significant mediators after BH correction. Skipping medsens.\n")
  medsens_results <- tibble()
} else {
  cat(sprintf(
    "Running medsens for %d significant mediator(s)...\n",
    nrow(sig_mediators)
  ))

  medsens_list <- list()
  for (i in seq_len(nrow(sig_mediators))) {
    mod_name <- sig_mediators$model[i]
    med_name <- sig_mediators$mediator[i]
    cat(sprintf("  %s / %s ... ", mod_name, mediator_labels[med_name]))

    med_obj <- mediate_objects[[mod_name]][[med_name]]
    if (is.null(med_obj)) {
      cat("no mediate object found, SKIP\n")
      next
    }

    sens <- tryCatch(
      {
        medsens(med_obj, rho.by = 0.05, effect.type = "indirect", sims = 50)
      },
      error = function(e) {
        cat(sprintf("ERROR: %s\n", e$message))
        return(NULL)
      }
    )

    if (!is.null(sens)) {
      rho_critical <- sens$rho.at.sign.change
      if (is.null(rho_critical) || length(rho_critical) == 0) {
        rho_critical <- NA_real_
        cat("rho_critical: ACME never crosses zero (very robust)\n")
      } else {
        cat(sprintf("rho_critical = %.3f\n", rho_critical))
      }

      medsens_list[[paste(mod_name, med_name, sep = "_")]] <- tibble(
        model = mod_name,
        mediator = med_name,
        label = mediator_labels[med_name],
        rho_at_sign_change = rho_critical
      )
    }
  }
  medsens_results <- bind_rows(medsens_list)
}

if (nrow(medsens_results) > 0) {
  write.csv(medsens_results, file.path(out_dir, "mediation_medsens_results.csv"),
    row.names = FALSE
  )
  cat("\nmedsens results saved.\n")
} else {
  cat("No medsens results to save.\n")
}

# ==============================================================================
# 5. Sensitivity Analysis 3: Separate reference groups (E2, E3)
# ==============================================================================

cat("\n=== Sensitivity Analysis 3: Separate Reference Groups ===\n\n")

run_refgroup_sensitivity <- function(type_num, ref_num) {
  model_name <- sprintf("Type%s_vs_E%s", type_num, ref_num)

  sub <- df[df$Enterotype %in% c(type_num, ref_num), ]
  sub$Treat <- ifelse(sub$Enterotype == type_num, 1L, 0L)

  cat(sprintf(
    "--- %s (n=%d: Treatment=%d, Reference=%d) ---\n",
    model_name, nrow(sub),
    sum(sub$Treat == 1), sum(sub$Treat == 0)
  ))

  results_list <- list()
  for (m in mediator_cols) {
    if (all(is.na(sub[[m]]))) next

    res <- run_mediation(sub, "Treat", m, "UA", sims = 5000)
    if (!is.null(res)) {
      res$model <- model_name
      res$label <- mediator_labels[m]
      results_list[[m]] <- res
    }
  }

  result <- bind_rows(results_list)
  if (nrow(result) > 0) {
    result$ACME_padj <- p.adjust(result$ACME_p, method = "BH")
  }
  cat(sprintf("  Completed: %d mediators tested\n", nrow(result)))
  result
}

ref_t1_e2 <- run_refgroup_sensitivity("1", "2")
ref_t1_e3 <- run_refgroup_sensitivity("1", "3")
ref_t4_e2 <- run_refgroup_sensitivity("4", "2")
ref_t4_e3 <- run_refgroup_sensitivity("4", "3")

refgroup_results <- bind_rows(ref_t1_e2, ref_t1_e3, ref_t4_e2, ref_t4_e3)

write.csv(refgroup_results, file.path(out_dir, "mediation_refgroup_sensitivity.csv"),
  row.names = FALSE
)
cat("\nReference group sensitivity results saved.\n")

# ==============================================================================
# 6. Save all R objects for downstream use
# ==============================================================================

all_results <- list(
  main = main_results,
  binary = binary_results,
  medsens = medsens_results,
  refgroup = refgroup_results,
  mediate_objects = mediate_objects
)
saveRDS(all_results, file.path(out_dir, "mediation_results_all.rds"))

# ==============================================================================
# 7. Summary report
# ==============================================================================

cat("\n========================================\n")
cat("Step 2 Summary Report\n")
cat("========================================\n\n")

cat("--- Main analysis: significant mediators (ACME padj < 0.05) ---\n")
sig_main <- main_results %>% filter(ACME_padj < 0.05)
if (nrow(sig_main) > 0) {
  sig_main %>%
    select(
      model, label, ACME_estimate, ACME_CI_lo, ACME_CI_hi,
      ACME_p, ACME_padj, Prop_mediated
    ) %>%
    print(n = Inf)
} else {
  cat("None at padj < 0.05.\n")
  cat("\nNominal significant (ACME p < 0.05):\n")
  main_results %>%
    filter(ACME_p < 0.05) %>%
    select(model, label, ACME_estimate, ACME_p, ACME_padj, Prop_mediated) %>%
    print(n = Inf)
}

cat("\n--- Binary HUA sensitivity: significant mediators (ACME padj < 0.05) ---\n")
sig_bin <- binary_results %>% filter(ACME_padj < 0.05)
if (nrow(sig_bin) > 0) {
  sig_bin %>%
    select(model, label, ACME_estimate, ACME_p, ACME_padj) %>%
    print(n = Inf)
} else {
  cat("None at padj < 0.05.\n")
}

cat("\n--- medsens results ---\n")
if (nrow(medsens_results) > 0) {
  print(medsens_results)
} else {
  cat("No medsens results (no significant mediators or skipped).\n")
}

cat("\n--- Reference group sensitivity: consistency check ---\n")
if (nrow(sig_main) > 0) {
  for (i in seq_len(nrow(sig_main))) {
    m <- sig_main$mediator[i]
    type_num <- ifelse(grepl("Type1", sig_main$model[i]), "1", "4")
    cat(sprintf("\n  %s in %s model:\n", mediator_labels[m], sig_main$model[i]))

    ref_sub <- refgroup_results %>%
      filter(mediator == m, grepl(paste0("Type", type_num), model))
    if (nrow(ref_sub) > 0) {
      ref_sub %>%
        select(model, ACME_estimate, ACME_p, ACME_padj) %>%
        print()
    }
  }
} else {
  cat("No significant mediators to check.\n")
}

cat("\nStep 2 complete.\n")
