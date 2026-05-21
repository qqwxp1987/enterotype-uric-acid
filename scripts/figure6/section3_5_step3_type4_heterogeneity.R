# ==============================================================================
# Section 3.5 — Step 3: Type 4 Internal Heterogeneity (Revised — 22 features)
#
# Purpose:
#   Type 4's HUA is subset-driven: group-level UA is not elevated, but HUA
#   prevalence equals Type 1 (30.8% vs 30.0%). This script compares the
#   HUA+ vs HUA- subgroups within Type 4 on all 22 candidate mediators
#   plus P. copri.
#
# Input:
#   - results/tables/section_3_5/mediation_analysis_data.rds (Step 1)
#   - phyloseq_Taxa_with_enterotype_seed10.rds (for P. copri abundance)
#
# Output (all in results/tables/section_3_5/):
#   - type4_heterogeneity_results.csv   Statistical results
#   - type4_heterogeneity_results.rds   Full results for visualization
#
# Comparisons: Type 4 HUA+ vs HUA- on 23 pre-specified features
# Methods: Wilcoxon rank-sum + Cliff's delta
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(here)

cat("========================================\n")
cat("Section 3.5 Step 3: Type 4 Heterogeneity (Revised — 22 features)\n")
cat("========================================\n\n")

# --- 0. Setup ---
out_dir <- here("results", "tables", "section_3_5")
df <- readRDS(file.path(out_dir, "mediation_analysis_data.rds"))

# ==============================================================================
# 1. Add P. copri relative abundance
# ==============================================================================

ps_taxa <- readRDS(here("data", "clean",
                        "phyloseq_Taxa_with_enterotype_seed10.rds"))
ps_rel <- transform_sample_counts(ps_taxa, function(x) x / sum(x))
otu <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) otu <- t(otu)

copri_row <- grep("Prevotella_copri", rownames(otu), value = TRUE)
cat(sprintf("P. copri taxon: %s\n", copri_row))

if (length(copri_row) == 1) {
  copri_abund <- otu[copri_row, ]
  df$P_copri <- copri_abund[df$SampleID]
} else {
  warning("P. copri not found or multiple matches")
  df$P_copri <- NA_real_
}

# ==============================================================================
# 2. Subset to Type 4 and define HUA groups
# ==============================================================================

t4 <- df[df$Enterotype == "4", ]
cat(sprintf("\nType 4 total: n = %d\n", nrow(t4)))
cat(sprintf("  HUA+ (HUA=1): n = %d\n", sum(t4$HUA == 1, na.rm = TRUE)))
cat(sprintf("  HUA- (HUA=0): n = %d\n", sum(t4$HUA == 0, na.rm = TRUE)))
cat(sprintf("  HUA NA:       n = %d\n", sum(is.na(t4$HUA))))

t4 <- t4[!is.na(t4$HUA), ]
t4$HUA_group <- ifelse(t4$HUA == 1, "HUA+", "HUA-")
t4$HUA_group <- factor(t4$HUA_group, levels = c("HUA-", "HUA+"))

cat(sprintf("\nAnalysis sample: n = %d (HUA+ = %d, HUA- = %d)\n",
            nrow(t4), sum(t4$HUA == 1), sum(t4$HUA == 0)))

cat("\n--- Demographics by HUA group ---\n")
cat("UA:\n")
tapply(t4$UA, t4$HUA_group, function(x)
  sprintf("  mean=%.1f, sd=%.1f, median=%.1f",
          mean(x), sd(x), median(x))) %>%
  print()

cat("\nAge:\n")
tapply(t4$Age, t4$HUA_group, function(x)
  sprintf("  mean=%.1f, sd=%.1f", mean(x), sd(x))) %>%
  print()

cat("\nGender:\n")
print(table(t4$HUA_group, t4$Gender))

# ==============================================================================
# 3. Define comparison variables (22 mediators + P. copri = 23)
# ==============================================================================

compare_vars <- c(
  "P_copri"        = "P. copri relative abundance",
  # Metabolites (9)
  "MeOIA"          = "5-Methoxyindoleacetate",
  "pCresol"        = "p-Cresol",
  "pCresolSulfate" = "p-Cresol sulfate",
  "IndoxylSulfate" = "Indoxyl sulfate",
  "IPA"            = "Indole-3-propionic acid",
  "PAGln"          = "Phenylacetylglutamine",
  "DCA"            = "Deoxycholic acid",
  "TMAO"           = "TMAO",
  "Adenosine"      = "Adenosine",
  # Pathways (7)
  "PWY5695"        = "PWY-5695 (urate biosynthesis)",
  "PWY0_1296"      = "PWY0-1296 (purine rib. degradation)",
  "PWY5004"        = "PWY-5004 (citrulline metabolism)",
  "PWY7219"        = "PWY-7219 (adenosine de novo)",
  "PWY6609"        = "PWY-6609 (adenine salvage)",
  "ARGORNPROST"    = "ARGORNPROST-PWY (Arg/Orn/Pro)",
  "PWY7234"        = "PWY-7234 (IMP biosynthesis)",
  # Convergent shared ECs (6)
  "EC_2_6_1_1"     = "EC 2.6.1.1 (aspartate transaminase)",
  "EC_2_3_1_129"   = "EC 2.3.1.129 (UDP-GlcNAc acyltransferase)",
  "EC_4_2_1_8"     = "EC 4.2.1.8",
  "EC_2_7_6_3"     = "EC 2.7.6.3 (diphosphokinase)",
  "EC_3_1_3_15"    = "EC 3.1.3.15 (histidinol-phosphatase)",
  "EC_2_1_1_13"    = "EC 2.1.1.13 (Met synthase)"
)

# ==============================================================================
# 4. Cliff's delta function
# ==============================================================================

cliffs_delta <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  comparisons <- outer(x, y, function(a, b) sign(a - b))
  delta <- mean(comparisons)

  # Interpretation thresholds (Romano et al. 2006)
  size <- case_when(
    abs(delta) < 0.147 ~ "negligible",
    abs(delta) < 0.33  ~ "small",
    abs(delta) < 0.474 ~ "medium",
    TRUE               ~ "large"
  )

  list(delta = delta, size = size)
}

# ==============================================================================
# 5. Run comparisons
# ==============================================================================

cat("\n========================================\n")
cat("Type 4 HUA+ vs HUA- Comparisons\n")
cat("========================================\n\n")

results_list <- list()

for (var_name in names(compare_vars)) {
  label <- compare_vars[var_name]
  vals <- t4[[var_name]]

  if (all(is.na(vals))) {
    cat(sprintf("%-18s: ALL NA, skipped\n", var_name))
    next
  }

  hua_pos <- vals[t4$HUA == 1]
  hua_neg <- vals[t4$HUA == 0]

  hua_pos <- hua_pos[!is.na(hua_pos)]
  hua_neg <- hua_neg[!is.na(hua_neg)]

  if (length(hua_pos) < 3 || length(hua_neg) < 3) {
    cat(sprintf("%-18s: insufficient non-NA values, skipped\n", var_name))
    next
  }

  wt <- wilcox.test(hua_pos, hua_neg, exact = FALSE)
  cd <- cliffs_delta(hua_pos, hua_neg)

  res <- tibble(
    variable       = var_name,
    label          = label,
    n_HUA_pos      = length(hua_pos),
    n_HUA_neg      = length(hua_neg),
    median_HUA_pos = median(hua_pos),
    median_HUA_neg = median(hua_neg),
    mean_HUA_pos   = mean(hua_pos),
    mean_HUA_neg   = mean(hua_neg),
    wilcox_W       = wt$statistic,
    wilcox_p       = wt$p.value,
    cliffs_delta   = cd$delta,
    effect_size    = cd$size
  )

  results_list[[var_name]] <- res

  cat(sprintf(
    "%-18s: HUA+ med=%.3f, HUA- med=%.3f | P=%.4f | delta=%.3f (%s)\n",
    var_name,
    median(hua_pos), median(hua_neg),
    wt$p.value, cd$delta, cd$size
  ))
}

results <- bind_rows(results_list)

# ==============================================================================
# 6. Summary
# ==============================================================================

cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n\n")

cat("--- Sorted by P value ---\n")
results %>%
  arrange(wilcox_p) %>%
  select(variable, label, wilcox_p, cliffs_delta, effect_size) %>%
  print(n = Inf)

cat("\n--- Features with P < 0.05 ---\n")
sig <- results %>% filter(wilcox_p < 0.05)
if (nrow(sig) > 0) {
  sig %>%
    select(variable, label, median_HUA_pos, median_HUA_neg,
           wilcox_p, cliffs_delta, effect_size) %>%
    print(n = Inf)
} else {
  cat("None at P < 0.05.\n")
}

cat("\n--- Features with |Cliff's delta| >= 0.33 (medium or large) ---\n")
med_large <- results %>% filter(abs(cliffs_delta) >= 0.33)
if (nrow(med_large) > 0) {
  med_large %>%
    select(variable, label, wilcox_p, cliffs_delta, effect_size) %>%
    print(n = Inf)
} else {
  cat("None with medium or large effect size.\n")
}

# ==============================================================================
# 7. Save
# ==============================================================================

write.csv(results, file.path(out_dir, "type4_heterogeneity_results.csv"),
          row.names = FALSE)
saveRDS(results, file.path(out_dir, "type4_heterogeneity_results.rds"))

cat(sprintf("\nResults saved to %s\n", out_dir))
cat("\nStep 3 complete.\n")
