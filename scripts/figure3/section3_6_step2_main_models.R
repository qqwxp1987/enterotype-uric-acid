# ============================================================================== 
# Section 3.6 — Step 2: Main Models
#
# Purpose:
#   Fit sequential robustness models for continuous UA and binary HUA.
#
# Primary strategy:
#   - Demographic model
#   - BMI-adjusted model
#   - Renal-adjusted model (Cr + Urea)
#   - Cardiometabolic model (Cr + Urea + GLU + TG + HDLC + LDLC)
#
# Sensitivity strategy within this script:
#   - Alternative renal model replacing Cr with eGFR
#   - Full clinical model adding BMI + CRP
# ============================================================================== 

library(tidyverse)
library(nlme)
library(here)

cat("========================================\n")
cat("Section 3.6 Step 2: Main Models\n")
cat("========================================\n\n")

out_dir <- here("results", "tables", "section_3_6")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

analysis_df <- readRDS(file.path(out_dir, "section3_6_analysis_data.rds"))

model_specs <- list(
  M1_Demographic = list(
    ua = UA ~ Enterotype + Gender + Age,
    hua = Hyperuricemia ~ Enterotype + Gender + Age,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age")
  ),
  M2_BMI = list(
    ua = UA ~ Enterotype + Gender + Age + BMI,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + BMI,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI")
  ),
  M3_Renal_Cr = list(
    ua = UA ~ Enterotype + Gender + Age + Cr + Urea,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + Cr + Urea,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "Cr", "Urea")
  ),
  M4_Cardiomet = list(
    ua = UA ~ Enterotype + Gender + Age + Cr + Urea + GLU + TG + HDLC + LDLC,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + Cr + Urea + GLU + TG + HDLC + LDLC,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "Cr", "Urea", "GLU", "TG", "HDLC", "LDLC")
  ),
  S1_Renal_eGFR = list(
    ua = UA ~ Enterotype + Gender + Age + eGFR + Urea,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + eGFR + Urea,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "eGFR", "Urea")
  ),
  S2_Full_Clinical = list(
    ua = UA ~ Enterotype + Gender + Age + BMI + Cr + Urea + GLU + TG + HDLC + LDLC + CRP,
    hua = Hyperuricemia ~ Enterotype + Gender + Age + BMI + Cr + Urea + GLU + TG + HDLC + LDLC + CRP,
    vars = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI", "Cr", "Urea", "GLU", "TG", "HDLC", "LDLC", "CRP")
  )
)

extract_wls <- function(fit, model_name, n_used) {
  fit_sum <- summary(fit)$tTable
  fit_anova <- anova(fit, type = "marginal")
  data.frame(
    model = model_name,
    n = n_used,
    term = rownames(fit_sum),
    estimate = fit_sum[, "Value"],
    se = fit_sum[, "Std.Error"],
    p = fit_sum[, "p-value"],
    stringsAsFactors = FALSE
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      ci_low = estimate - 1.96 * se,
      ci_high = estimate + 1.96 * se,
      outcome = "UA",
      overall_enterotype_p = fit_anova["Enterotype", "p-value"]
    )
}

extract_logit <- function(fit, null_fit, model_name, n_used) {
  fit_sum <- summary(fit)$coefficients
  ci <- suppressMessages(confint.default(fit))
  lr <- anova(null_fit, fit, test = "Chisq")
  terms <- rownames(fit_sum)
  data.frame(
    model = model_name,
    n = n_used,
    term = terms,
    estimate = exp(fit_sum[, "Estimate"]),
    se = fit_sum[, "Std. Error"],
    p = fit_sum[, "Pr(>|z|)"],
    ci_low = exp(ci[, 1]),
    ci_high = exp(ci[, 2]),
    stringsAsFactors = FALSE
  ) %>%
    filter(term != "(Intercept)") %>%
    mutate(
      outcome = "HUA",
      overall_enterotype_p = lr$`Pr(>Chi)`[2]
    )
}

clean_labels <- function(df) {
  df %>%
    mutate(
      label = case_when(
        term == "EnterotypeType 2" ~ "Type 2 vs Type 1",
        term == "EnterotypeType 3" ~ "Type 3 vs Type 1",
        term == "EnterotypeType 4" ~ "Type 4 vs Type 1",
        term == "GenderMale" ~ "Sex (Male vs Female)",
        term == "Age" ~ "Age (per year)",
        term == "BMI" ~ "BMI",
        term == "Cr" ~ "Creatinine",
        term == "Urea" ~ "Urea",
        term == "eGFR" ~ "eGFR",
        term == "GLU" ~ "Glucose",
        term == "TG" ~ "Triglycerides",
        term == "HDLC" ~ "HDL-C",
        term == "LDLC" ~ "LDL-C",
        term == "CRP" ~ "CRP",
        TRUE ~ term
      )
    )
}

wls_results <- list()
logit_results <- list()
model_summary <- list()

for (model_name in names(model_specs)) {
  spec <- model_specs[[model_name]]
  df <- analysis_df[complete.cases(analysis_df[, spec$vars, drop = FALSE]), ]
  n_used <- nrow(df)

  cat(sprintf("Fitting %s (n = %d)\n", model_name, n_used))

  wls_fit <- gls(
    spec$ua,
    data = df,
    weights = varIdent(form = ~ 1 | Enterotype),
    method = "REML",
    na.action = na.omit
  )

  logit_fit <- glm(spec$hua, data = df, family = binomial)
  logit_null <- glm(update(spec$hua, . ~ . - Enterotype), data = df, family = binomial)

  wls_results[[model_name]] <- extract_wls(wls_fit, model_name, n_used)
  logit_results[[model_name]] <- extract_logit(logit_fit, logit_null, model_name, n_used)

  model_summary[[model_name]] <- tibble(
    model = model_name,
    n = n_used,
    ua_enterotype_p = anova(wls_fit, type = "marginal")["Enterotype", "p-value"],
    hua_enterotype_p = anova(logit_null, logit_fit, test = "Chisq")$`Pr(>Chi)`[2]
  )
}

wls_df <- bind_rows(wls_results) %>% clean_labels()
logit_df <- bind_rows(logit_results) %>% clean_labels()
model_summary_df <- bind_rows(model_summary)

write.csv(wls_df, file.path(out_dir, "ua_model_comparison.csv"), row.names = FALSE)
write.csv(logit_df, file.path(out_dir, "hua_model_comparison.csv"), row.names = FALSE)
write.csv(model_summary_df, file.path(out_dir, "model_summary.csv"), row.names = FALSE)
saveRDS(list(wls = wls_df, logit = logit_df, summary = model_summary_df),
        file.path(out_dir, "section3_6_model_results.rds"))

cat("\nSaved Section 3.6 main model results.\n")
print(model_summary_df)
