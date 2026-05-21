# ============================================================================== 
# Section 3.6 — Step 3: Sensitivity Analyses
#
# Purpose:
#   Test whether the enterotype effect direction remains consistent in key host
#   subgroups and restricted samples.
# ============================================================================== 

library(tidyverse)
library(nlme)
library(here)

cat("========================================\n")
cat("Section 3.6 Step 3: Sensitivity Analyses\n")
cat("========================================\n\n")

out_dir <- here("results", "tables", "section_3_6")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

analysis_df <- readRDS(file.path(out_dir, "section3_6_analysis_data.rds"))

subsets <- list(
  Male = analysis_df %>% filter(Gender == "Male"),
  Female = analysis_df %>% filter(Gender == "Female"),
  Non_obese = analysis_df %>% filter(!is.na(obese_cn) & obese_cn == FALSE),
  eGFR_normal = analysis_df %>% filter(!is.na(egfr_normal) & egfr_normal == TRUE)
)

fit_subset_models <- function(df, subset_name) {
  df <- df %>% filter(complete.cases(UA, Hyperuricemia, Enterotype, Age))

  ua_formula <- if (dplyr::n_distinct(df$Gender) > 1) {
    UA ~ Enterotype + Gender + Age
  } else {
    UA ~ Enterotype + Age
  }

  hua_formula <- if (dplyr::n_distinct(df$Gender) > 1) {
    Hyperuricemia ~ Enterotype + Gender + Age
  } else {
    Hyperuricemia ~ Enterotype + Age
  }

  wls_fit <- gls(
    ua_formula,
    data = df,
    weights = varIdent(form = ~ 1 | Enterotype),
    method = "REML",
    na.action = na.omit
  )

  logit_fit <- glm(hua_formula, data = df, family = binomial)
  logit_null <- glm(update(hua_formula, . ~ . - Enterotype), data = df, family = binomial)

  wls_tab <- summary(wls_fit)$tTable
  logit_tab <- summary(logit_fit)$coefficients
  logit_ci <- suppressMessages(confint.default(logit_fit))
  lr <- anova(logit_null, logit_fit, test = "Chisq")

  wls_df <- data.frame(
    subset = subset_name,
    n = nrow(df),
    outcome = "UA",
    term = rownames(wls_tab),
    estimate = wls_tab[, "Value"],
    se = wls_tab[, "Std.Error"],
    p = wls_tab[, "p-value"],
    stringsAsFactors = FALSE
  ) %>%
    filter(grepl("^Enterotype", term)) %>%
    mutate(
      ci_low = estimate - 1.96 * se,
      ci_high = estimate + 1.96 * se,
      overall_enterotype_p = anova(wls_fit, type = "marginal")["Enterotype", "p-value"]
    )

  logit_df <- data.frame(
    subset = subset_name,
    n = nrow(df),
    outcome = "HUA",
    term = rownames(logit_tab),
    estimate = exp(logit_tab[, "Estimate"]),
    se = logit_tab[, "Std. Error"],
    p = logit_tab[, "Pr(>|z|)"],
    ci_low = exp(logit_ci[, 1]),
    ci_high = exp(logit_ci[, 2]),
    stringsAsFactors = FALSE
  ) %>%
    filter(grepl("^Enterotype", term)) %>%
    mutate(overall_enterotype_p = lr$`Pr(>Chi)`[2])

  bind_rows(wls_df, logit_df)
}

subgroup_results <- bind_rows(
  lapply(names(subsets), function(name) fit_subset_models(subsets[[name]], name))
) %>%
  mutate(
    label = case_when(
      term == "EnterotypeType 2" ~ "Type 2 vs Type 1",
      term == "EnterotypeType 3" ~ "Type 3 vs Type 1",
      term == "EnterotypeType 4" ~ "Type 4 vs Type 1",
      TRUE ~ term
    )
  )

write.csv(subgroup_results, file.path(out_dir, "subgroup_results.csv"), row.names = FALSE)
saveRDS(subgroup_results, file.path(out_dir, "subgroup_results.rds"))

cat("Saved subgroup sensitivity results.\n")
print(subgroup_results %>% select(subset, outcome, label, n, estimate, p, overall_enterotype_p))
