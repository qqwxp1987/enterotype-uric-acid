# ==============================================================================
# Table 1 Upgrade: Publication-Ready Cohort Demographics by Enterotype
# Output: CSV + formatted text for manuscript
# Data source: phyloseq object (canonical, consistent with Section 3.1)
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta <- as(sample_data(ps), "data.frame")

meta$Enterotype <- factor(
  meta$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)

# Ensure correct types
meta$Gender   <- factor(meta$Gender, levels = c("Female", "Male"))
meta$Age      <- as.numeric(meta$Age)
meta$BMI      <- as.numeric(meta$BMI)
meta$UA       <- as.numeric(meta$UA)
meta$Hyperuricemia <- factor(meta$Hyperuricemia, levels = c("No", "Yes"))

enterotypes <- levels(meta$Enterotype)
total_n     <- nrow(meta)

# --- 2. Helper Functions ---
fmt_mean_sd <- function(x) {
  x <- x[!is.na(x)]
  sprintf("%.1f \u00b1 %.1f", mean(x), sd(x))
}

fmt_n_pct <- function(x, target) {
  n     <- sum(x == target, na.rm = TRUE)
  total <- sum(!is.na(x))
  sprintf("%d (%.1f)", n, n / total * 100)
}

fmt_p <- function(p) {
  if (p < 0.001) {
    sprintf("< 0.001")
  } else {
    sprintf("%.3f", p)
  }
}

# --- 3. Build Rows ---
rows <- list()

# --- N ---
n_total <- as.character(total_n)
n_by_et <- sapply(enterotypes, function(e) as.character(sum(meta$Enterotype == e)))
rows[["N"]] <- c("N", n_total, n_by_et, "")

# --- Sex (Male), n (%) ---
sex_test <- chisq.test(table(meta$Enterotype, meta$Gender))
rows[["Sex"]] <- c(
  "Sex, male, n (%)",
  fmt_n_pct(meta$Gender, "Male"),
  sapply(enterotypes, function(e) {
    sub <- meta[meta$Enterotype == e, ]
    fmt_n_pct(sub$Gender, "Male")
  }),
  paste0(fmt_p(sex_test$p.value), " \u2020")
)

# --- Age (mean +/- sd) ---
age_aov <- summary(aov(Age ~ Enterotype, data = meta))
rows[["Age"]] <- c(
  "Age, years",
  fmt_mean_sd(meta$Age),
  sapply(enterotypes, function(e) fmt_mean_sd(meta$Age[meta$Enterotype == e])),
  paste0(fmt_p(age_aov[[1]]$`Pr(>F)`[1]), " \u2021")
)

# --- BMI (mean +/- sd) ---
meta_bmi <- meta[!is.na(meta$BMI), ]
bmi_aov  <- summary(aov(BMI ~ Enterotype, data = meta_bmi))
bmi_total <- fmt_mean_sd(meta_bmi$BMI)
bmi_by_et <- sapply(enterotypes, function(e) {
  x <- meta_bmi$BMI[meta_bmi$Enterotype == e]
  if (length(x[!is.na(x)]) == 0) return("—")
  fmt_mean_sd(x)
})
rows[["BMI"]] <- c(
  "BMI, kg/m\u00b2",
  bmi_total,
  bmi_by_et,
  paste0(fmt_p(bmi_aov[[1]]$`Pr(>F)`[1]), " \u2021")
)

# --- Serum UA (mean +/- sd) ---
ua_aov <- summary(aov(UA ~ Enterotype, data = meta))
rows[["UA"]] <- c(
  "Serum UA, \u03bcmol/L",
  fmt_mean_sd(meta$UA),
  sapply(enterotypes, function(e) fmt_mean_sd(meta$UA[meta$Enterotype == e])),
  paste0(fmt_p(ua_aov[[1]]$`Pr(>F)`[1]), " \u2021")
)

# --- Hyperuricemia, n (%) ---
meta_hua <- meta[!is.na(meta$Hyperuricemia), ]
hua_test <- chisq.test(table(meta_hua$Enterotype, meta_hua$Hyperuricemia))
rows[["HUA"]] <- c(
  "Hyperuricemia, n (%)",
  fmt_n_pct(meta_hua$Hyperuricemia, "Yes"),
  sapply(enterotypes, function(e) {
    sub <- meta_hua[meta_hua$Enterotype == e, ]
    fmt_n_pct(sub$Hyperuricemia, "Yes")
  }),
  paste0(fmt_p(hua_test$p.value), " \u2020")
)

# --- 4. Assemble Table ---
table1 <- do.call(rbind, rows)
colnames(table1) <- c("Characteristic", "Total", enterotypes, "P-value")
table1_df <- as.data.frame(table1, stringsAsFactors = FALSE)

# --- 5. Add footnote row ---
footnote_row <- data.frame(
  Characteristic = "Footnotes: \u2020 Pearson\u2019s Chi-squared test; \u2021 One-way ANOVA. Values are mean \u00b1 SD or n (%). BMI available for N=226.",
  Total = "", `Type 1` = "", `Type 2` = "", `Type 3` = "", `Type 4` = "", `P.value` = "",
  stringsAsFactors = FALSE
)
colnames(footnote_row) <- colnames(table1_df)

table1_final <- rbind(table1_df, footnote_row)

# --- 6. Save ---
out_dir <- here("results", "tables", "section_3_2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  table1_final,
  file.path(out_dir, "Table1_Demographics_v2.csv"),
  row.names = FALSE
)

# --- 7. Print formatted output ---
cat("\n============================================================\n")
cat("  Table 1. Cohort Demographics and Clinical Characteristics\n")
cat("           Stratified by Enterotype\n")
cat("============================================================\n\n")

# Print header
header <- sprintf("%-24s %12s %12s %12s %12s %12s %12s",
                  "Characteristic", "Total", "Type 1", "Type 2", "Type 3", "Type 4", "P-value")
cat(header, "\n")
cat(paste(rep("-", nchar(header)), collapse = ""), "\n")

# Print data rows (excluding footnote)
for (i in 1:nrow(table1_df)) {
  row_str <- sprintf("%-24s %12s %12s %12s %12s %12s %12s",
                     table1_df[i, 1], table1_df[i, 2], table1_df[i, 3],
                     table1_df[i, 4], table1_df[i, 5], table1_df[i, 6],
                     table1_df[i, 7])
  cat(row_str, "\n")
}
cat(paste(rep("-", nchar(header)), collapse = ""), "\n")
cat("\n\u2020 Pearson\u2019s Chi-squared test; \u2021 One-way ANOVA\n")
cat("Values are mean \u00b1 SD or n (%).\n")
cat("BMI data available for", sum(!is.na(meta$BMI)), "of", total_n, "participants.\n")
cat("\nSaved to:", file.path(out_dir, "Table1_Demographics_v2.csv"), "\n")
