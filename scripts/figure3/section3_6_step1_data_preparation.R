# ============================================================================== 
# Section 3.6 — Step 1: Data Preparation
#
# Purpose:
#   Build a clean phenotype table for robustness analyses of the enterotype-UA/HUA
#   association under renal and cardiometabolic adjustment.
#
# Notes:
#   - Uses the same authoritative phenotype source as earlier sections.
#   - Keeps eGFR and Cr as alternative renal adjustment strategies rather than
#     forcing both into the same primary model, because eGFR is partially derived
#     from creatinine, age, and sex.
# ============================================================================== 

library(phyloseq)
library(tidyverse)
library(here)

cat("========================================\n")
cat("Section 3.6 Step 1: Data Preparation\n")
cat("========================================\n\n")

out_dir <- here("results", "tables", "section_3_6")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta <- as(sample_data(ps), "data.frame")
meta$SampleID <- rownames(meta)

convert_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

analysis_df <- meta %>%
  transmute(
    SampleID = SampleID,
    Enterotype = factor(
      Dirichlettype_Genus,
      levels = c("1", "2", "3", "4"),
      labels = c("Type 1", "Type 2", "Type 3", "Type 4")
    ),
    Gender = factor(Gender, levels = c("Female", "Male")),
    Age = convert_num(Age),
    UA = convert_num(UA),
    Hyperuricemia = factor(Hyperuricemia, levels = c("No", "Yes")),
    BMI = convert_num(BMI),
    Urea = convert_num(Urea),
    Cr = convert_num(Cr),
    eGFR = convert_num(eGFR),
    GLU = convert_num(GLU),
    TG = convert_num(TG),
    HDLC = convert_num(HDLC),
    LDLC = convert_num(LDLC),
    CRP = convert_num(CRP),
    SBP = convert_num(SBP),
    DBP = convert_num(DBP)
  ) %>%
  mutate(
    Enterotype = relevel(Enterotype, ref = "Type 1"),
    obese_cn = case_when(
      is.na(BMI) ~ NA,
      BMI >= 28 ~ TRUE,
      BMI < 28 ~ FALSE
    ),
    egfr_normal = case_when(
      is.na(eGFR) ~ NA,
      eGFR >= 90 ~ TRUE,
      eGFR < 90 ~ FALSE
    )
  )

stopifnot(all(levels(analysis_df$Enterotype) == c("Type 1", "Type 2", "Type 3", "Type 4")))

model_vars <- list(
  demo = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age"),
  bmi = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI"),
  renal_cr = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "Cr", "Urea"),
  renal_egfr = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "eGFR", "Urea"),
  cardiomet = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "Cr", "Urea", "GLU", "TG", "HDLC", "LDLC"),
  full_clinical = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI", "Cr", "Urea", "GLU", "TG", "HDLC", "LDLC", "CRP"),
  bp_extended = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "Cr", "Urea", "GLU", "TG", "HDLC", "LDLC", "SBP", "DBP")
)

missingness_summary <- tibble(
  variable = c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI", "Urea", "Cr", "eGFR", "GLU", "TG", "HDLC", "LDLC", "CRP", "SBP", "DBP"),
  n_nonmissing = sapply(c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI", "Urea", "Cr", "eGFR", "GLU", "TG", "HDLC", "LDLC", "CRP", "SBP", "DBP"), function(v) sum(!is.na(analysis_df[[v]]))),
  pct_nonmissing = round(sapply(c("UA", "Hyperuricemia", "Enterotype", "Gender", "Age", "BMI", "Urea", "Cr", "eGFR", "GLU", "TG", "HDLC", "LDLC", "CRP", "SBP", "DBP"), function(v) mean(!is.na(analysis_df[[v]])) * 100), 1)
)

model_n_summary <- tibble(
  model = names(model_vars),
  n_complete = sapply(model_vars, function(vars) sum(complete.cases(analysis_df[, vars, drop = FALSE])))
)

covariate_distribution <- analysis_df %>%
  group_by(Enterotype) %>%
  summarise(
    n = n(),
    male_pct = mean(Gender == "Male", na.rm = TRUE) * 100,
    age_mean = mean(Age, na.rm = TRUE),
    bmi_mean = mean(BMI, na.rm = TRUE),
    urea_mean = mean(Urea, na.rm = TRUE),
    cr_mean = mean(Cr, na.rm = TRUE),
    egfr_mean = mean(eGFR, na.rm = TRUE),
    glu_mean = mean(GLU, na.rm = TRUE),
    tg_mean = mean(TG, na.rm = TRUE),
    hdlc_mean = mean(HDLC, na.rm = TRUE),
    ldlc_mean = mean(LDLC, na.rm = TRUE),
    crp_mean = mean(CRP, na.rm = TRUE)
  )

saveRDS(analysis_df, file.path(out_dir, "section3_6_analysis_data.rds"))
write.csv(analysis_df, file.path(out_dir, "section3_6_analysis_data.csv"), row.names = FALSE)
write.csv(missingness_summary, file.path(out_dir, "missingness_summary.csv"), row.names = FALSE)
write.csv(model_n_summary, file.path(out_dir, "model_complete_case_counts.csv"), row.names = FALSE)
write.csv(covariate_distribution, file.path(out_dir, "covariate_distribution_by_enterotype.csv"), row.names = FALSE)

cat("Saved Section 3.6 analysis data and summaries.\n\n")
print(missingness_summary)
cat("\nComplete-case counts by model:\n")
print(model_n_summary)
