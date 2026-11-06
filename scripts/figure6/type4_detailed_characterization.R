# ==============================================================================
# Type 4 Enterotype HUA+ Subgroup Detailed Characterization
#
# Purpose:
#   Respond to reviewer comments:
#     R1.1 — Differentiate Type 1 vs Type 4 evidence (Type 4 signal may be
#            driven by a few high-UA individuals)
#     R2.2 — Caution about Type 4 small sample (n = 39)
#
#   Provides five complementary analyses:
#     1. Type 4 HUA+ / HUA- subgroup descriptive statistics
#     2. Influence diagnostics (Cook's distance, leverage, DFBETAS)
#     3. Outlier removal sensitivity (P95 / P99 thresholds)
#     4. Leave-one-out (LOO) stability of the Type 4–UA association
#     5. Bootstrap stability of the Type 4–UA association
#
# Inputs:
#   - data/clean/phyloseq_Taxa_with_enterotype_seed10.rds  (phyloseq + metadata)
#   - results/tables/section_3_5/mediation_analysis_data.rds (unified df)
#
# Outputs:
#   Tables (results/tables/supplementary/):
#     - type4_subgroup_characterization.csv
#     - type4_influence_diagnostics.csv
#     - type4_sensitivity_loo.csv          (LOO + outlier removal combined)
#     - type4_bootstrap_stability.csv
#   Figures (results/figures/panels/):
#     - FigureS3_type4_diagnostics.png/pdf  (combined: A=Cook's, B=LOO)
#     - FigureS_type4_cooks_distance.png/pdf (standalone)
#     - FigureS_type4_leverage_diagnostic.png/pdf (standalone)
#     - FigureS_type4_loo_stability.png/pdf (standalone)
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(here)
library(ggplot2)

# --- Style (shared with Figure 6 panels) ---
source(here("scripts", "figure6", "fig6_style.R"))
source(here("scripts", "utils", "palettes.R"))

# --- Reproducibility ---
set.seed(20250904)

cat("================================================================\n")
cat("Type 4 Detailed Characterization\n")
cat("  Reviewer responses: R1.1 (Type 1/4 differentiation),\n")
cat("                       R2.2 (Type 4 small-sample caution)\n")
cat("================================================================\n\n")

# --- Output directories ---
tbl_out <- here("results", "tables", "supplementary")
fig_out <- here("results", "figures", "panels")
dir.create(tbl_out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_out, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 0. Load and prepare data
# ==============================================================================

# --- Unified mediation data (has Enterotype "1"-"4", UA, Gender, Age, HUA) ---
med_df <- readRDS(here("results", "tables", "section_3_5",
                       "mediation_analysis_data.rds"))
cat(sprintf("Mediation data: %d samples, %d variables\n", nrow(med_df), ncol(med_df)))

# --- Phyloseq metadata (for BMI, Cr, eGFR) ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
ps_meta <- as(sample_data(ps), "data.frame")
ps_meta$SampleID <- rownames(ps_meta)

# --- Merge BMI into mediation data ---
extra_vars <- ps_meta[, c("SampleID", "BMI", "Cr", "eGFR")]
df <- med_df %>%
  left_join(extra_vars, by = "SampleID") %>%
  mutate(
    BMI = suppressWarnings(as.numeric(as.character(BMI))),
    Cr  = suppressWarnings(as.numeric(as.character(Cr))),
    eGFR = suppressWarnings(as.numeric(as.character(eGFR))),
    Enterotype_label = paste0("Type ", Enterotype),
    Gender = factor(Gender, levels = c("Female", "Male"))
  )

cat(sprintf("Merged data: %d samples\n", nrow(df)))
cat("\nEnterotype distribution:\n")
print(table(df$Enterotype, useNA = "ifany"))

# --- Define reference group: Type 2 + 3 ---
# Treat_T4 is already 1 (Type 4) / 0 (Type 2+3) / NA (Type 1)
df$Is_Type4 <- df$Treat_T4  # 1 = Type 4, 0 = Type 2+3, NA = Type 1

# --- Analysis subsets ---
t4_all  <- df[df$Enterotype == "4" & !is.na(df$UA), ]
t4_ref  <- df[!is.na(df$Is_Type4) & !is.na(df$UA), ]  # Type 4 + Type 2+3

cat(sprintf("\nType 4: n = %d (HUA+ = %d, HUA- = %d)\n",
            nrow(t4_all), sum(t4_all$HUA == 1, na.rm = TRUE),
            sum(t4_all$HUA == 0, na.rm = TRUE)))
cat(sprintf("Type 4 vs Ref analysis set: n = %d (Type 4 = %d, Ref = %d)\n",
            nrow(t4_ref), sum(t4_ref$Is_Type4 == 1), sum(t4_ref$Is_Type4 == 0)))

# --- Anonymous ID mapping for Type 4 samples (for figure labeling) ---
# Sort by UA to make the anonymous IDs partially informative
t4_anon_map <- t4_all %>%
  arrange(UA) %>%
  mutate(AnonID = sprintf("T4-%02d", row_number())) %>%
  select(SampleID, AnonID, UA, HUA, Gender, Age)

cat(sprintf("\nAnonymous ID mapping created for %d Type 4 samples.\n", nrow(t4_anon_map)))
cat("  ID range: T4-01 (lowest UA) to T4-39 (highest UA)\n")


# ==============================================================================
# ANALYSIS 1: Type 4 HUA+ / HUA- subgroup description
# ==============================================================================

cat("\n================================================================\n")
cat("ANALYSIS 1: Type 4 HUA+ / HUA- Subgroup Description\n")
cat("================================================================\n\n")

# --- Helper: descriptive statistics ---
desc_stats <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(c(n = 0, mean = NA, sd = NA, median = NA,
                               Q1 = NA, Q3 = NA, IQR = NA, min = NA, max = NA))
  qs <- quantile(x, c(0.25, 0.75), na.rm = TRUE)
  c(n = length(x),
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    Q1 = unname(qs[1]),
    Q3 = unname(qs[2]),
    IQR = unname(qs[2] - qs[1]),
    min = min(x),
    max = max(x))
}

# --- Build subgroup table ---
subgroups <- list(
  "Type4_HUA_pos" = t4_all[t4_all$HUA == 1, ],
  "Type4_HUA_neg" = t4_all[t4_all$HUA == 0, ]
)

vars_to_report <- c("UA", "Age", "BMI")
subgroup_rows <- list()

# Helper to build one row of descriptive stats
make_stat_row <- function(grp_name, variable, ds) {
  data.frame(
    subgroup = grp_name,
    variable = variable,
    n = unname(ds["n"]),
    mean = round(unname(ds["mean"]), 2),
    sd = round(unname(ds["sd"]), 2),
    median = round(unname(ds["median"]), 2),
    Q1 = round(unname(ds["Q1"]), 2),
    Q3 = round(unname(ds["Q3"]), 2),
    IQR = round(unname(ds["IQR"]), 2),
    min = round(unname(ds["min"]), 2),
    max = round(unname(ds["max"]), 2),
    stringsAsFactors = FALSE
  )
}

# Helper to build a Gender row with M/F counts
make_gender_row <- function(grp_name, grp_data) {
  gtab <- table(grp_data$Gender)
  n_male   <- if ("Male" %in% names(gtab)) unname(gtab["Male"]) else 0
  n_female <- if ("Female" %in% names(gtab)) unname(gtab["Female"]) else 0
  data.frame(
    subgroup = grp_name,
    variable = "Gender",
    n = sum(gtab),
    mean = NA, sd = NA, median = NA, Q1 = NA, Q3 = NA, IQR = NA,
    min = NA, max = NA,
    stringsAsFactors = FALSE
  )
}

for (grp_name in names(subgroups)) {
  grp_data <- subgroups[[grp_name]]
  for (v in vars_to_report) {
    ds <- desc_stats(grp_data[[v]])
    subgroup_rows[[length(subgroup_rows) + 1]] <- make_stat_row(grp_name, v, ds)
  }
  # Gender breakdown (encode M/F in a separate column set below)
  subgroup_rows[[length(subgroup_rows) + 1]] <- make_gender_row(grp_name, grp_data)
}

# Also add overall Type 4 and reference group for context
for (grp_name in c("Type4_all", "Type2_3_ref")) {
  if (grp_name == "Type4_all") {
    grp_data <- t4_all
  } else {
    grp_data <- df[df$Enterotype %in% c("2", "3") & !is.na(df$UA), ]
  }
  for (v in vars_to_report) {
    ds <- desc_stats(grp_data[[v]])
    subgroup_rows[[length(subgroup_rows) + 1]] <- make_stat_row(grp_name, v, ds)
  }
  subgroup_rows[[length(subgroup_rows) + 1]] <- make_gender_row(grp_name, grp_data)
}

subgroup_table <- bind_rows(subgroup_rows)

# Add Gender M/F counts as a note
gender_table <- t4_all %>%
  count(HUA_group = ifelse(HUA == 1, "HUA+", "HUA-"), Gender) %>%
  pivot_wider(names_from = Gender, values_from = n, values_fill = 0)

cat("--- Type 4 subgroup descriptive statistics ---\n")
print(subgroup_table, row.names = FALSE)
cat("\n--- Gender breakdown (Type 4) ---\n")
print(gender_table)

write.csv(subgroup_table,
          file.path(tbl_out, "type4_subgroup_characterization.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved: %s\n",
            file.path(tbl_out, "type4_subgroup_characterization.csv")))


# ==============================================================================
# ANALYSIS 1 (figure): UA distribution — REMOVED
#   The UA distribution across enterotypes is already shown in Figure 2B.
#   Type 4 HUA+/- subgroup statistics are reported in the text and Table S10.
#   No incremental value from a separate figure.
# ==============================================================================


# ==============================================================================
# ANALYSIS 2: Influence diagnostics
# ==============================================================================

cat("\n================================================================\n")
cat("ANALYSIS 2: Influence Diagnostics (Type 4 vs Ref)\n")
cat("================================================================\n\n")

# --- Fit the model: UA ~ Type4 + Gender + Age ---
# Reference group = Type 2+3 (Is_Type4 = 0), Treatment = Type 4 (Is_Type4 = 1)
t4_ref_complete <- t4_ref %>%
  filter(!is.na(UA), !is.na(Gender), !is.na(Age))

cat(sprintf("Complete cases for influence model: n = %d\n", nrow(t4_ref_complete)))

fit_ua <- lm(UA ~ Is_Type4 + Gender + Age, data = t4_ref_complete)
cat("\n--- Model summary ---\n")
print(summary(fit_ua))

# --- Influence measures ---
n_obs <- nrow(t4_ref_complete)
cooks <- cooks.distance(fit_ua)
hats  <- hatvalues(fit_ua)
dfb   <- dfbetas(fit_ua)

# DFBETAS for the Type 4 coefficient (Is_Type4)
dfb_type4 <- dfb[, "Is_Type4"]

# --- Thresholds ---
cook_threshold_4n  <- 4 / n_obs    # standard rule
cook_threshold_0.5 <- 0.5          # R1.1/R2.2 higher-influence reference
hat_threshold      <- 2 * (ncol(model.matrix(fit_ua)) / n_obs)  # 2p/n
dfb_threshold      <- 2 / sqrt(n_obs)

# --- Flag high-influence individuals ---
infl_df <- data.frame(
  SampleID   = t4_ref_complete$SampleID,  # R1.1/R2.2 preserve IDs for anonymous labels
  Enterotype = t4_ref_complete$Enterotype,
  Is_Type4   = t4_ref_complete$Is_Type4,
  UA         = t4_ref_complete$UA,
  Gender     = as.character(t4_ref_complete$Gender),
  Age        = t4_ref_complete$Age,
  cooks_d    = round(cooks, 5),
  leverage   = round(hats, 5),
  dfb_type4  = round(dfb_type4, 5),
  stringsAsFactors = FALSE
) %>%
  mutate(
    cooks_flag_4n   = cooks_d > cook_threshold_4n,
    cooks_flag_0.5  = cooks_d > 0.5,
    leverage_flag   = leverage > hat_threshold,
    dfb_flag        = abs(dfb_type4) > dfb_threshold,
    any_flag        = cooks_flag_4n | leverage_flag | dfb_flag
  )

cat(sprintf("\nInfluence thresholds:\n"))
cat(sprintf("  Cook's > 4/n  = %.5f (n = %d): %d flagged\n",
            cook_threshold_4n, n_obs, sum(infl_df$cooks_flag_4n)))
cat(sprintf("  Cook's > 0.5              : %d flagged\n", sum(infl_df$cooks_flag_0.5)))
cat(sprintf("  Leverage > 2p/n = %.5f     : %d flagged\n",
            hat_threshold, sum(infl_df$leverage_flag)))
cat(sprintf("  |DFBETAS_Type4| > %.4f     : %d flagged\n",
            dfb_threshold, sum(infl_df$dfb_flag)))

cat("\n--- Top 10 by Cook's distance ---\n")
infl_df %>%
  arrange(desc(cooks_d)) %>%
  head(10) %>%
  select(SampleID, Enterotype, UA, cooks_d, leverage, dfb_type4,
         cooks_flag_4n, leverage_flag, dfb_flag) %>%
  print(row.names = FALSE)

# --- Type 4 specific high-influence ---
t4_infl <- infl_df %>% filter(Is_Type4 == 1)
cat(sprintf("\nType 4 samples flagged (Cook's > 4/n): %d / %d\n",
            sum(t4_infl$cooks_flag_4n), nrow(t4_infl)))

write.csv(infl_df, file.path(tbl_out, "type4_influence_diagnostics.csv"),
          row.names = FALSE)
cat(sprintf("Saved: %s\n", file.path(tbl_out, "type4_influence_diagnostics.csv")))


# ==============================================================================
# ANALYSIS 2 (figure): Cook's distance bar chart (→ Figure S5, Panel A)
# ==============================================================================

# Merge anonymous IDs for Type 4 samples
infl_df <- infl_df %>%
  left_join(t4_anon_map %>% select(SampleID, AnonID), by = "SampleID") %>%
  mutate(
    # For non-Type 4 samples, leave AnonID as NA (not labeled)
    display_label = ifelse(Is_Type4 == 1 & cooks_flag_4n, AnonID, NA_character_)
  )

# Panel A: Cook's distance bar plot (anonymized)
infl_plot <- infl_df %>%
  mutate(
    row_idx = row_number(),
    is_t4   = Is_Type4 == 1,
    label_color = case_when(
      cooks_flag_0.5 ~ "#E64B35",
      cooks_flag_4n & is_t4 ~ "#E64B35",
      cooks_flag_4n & !is_t4 ~ "#F39B7F",
      is_t4 ~ "#3C5488",
      TRUE ~ "grey60"
    )
  ) %>%
  arrange(desc(cooks_d)) %>%
  mutate(rank_idx = row_number())

p_cooks <- ggplot(infl_plot, aes(x = rank_idx, y = cooks_d, fill = label_color)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = cook_threshold_4n, linetype = "dashed",
             color = "#E64B35", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, linetype = "dotted",
             color = "grey40", linewidth = 0.4) +
  annotate("text", x = nrow(infl_plot) * 0.95, y = cook_threshold_4n + 0.012,
           label = sprintf("4/n = %.4f", cook_threshold_4n),
           size = 3.5, color = "#E64B35", hjust = 1, fontface = "italic") +
  annotate("text", x = nrow(infl_plot) * 0.95, y = 0.52,
           label = "0.5 reference", size = 3, color = "grey40", hjust = 1) +
  # Label Type 4 high-influence points with anonymous IDs
  ggrepel::geom_text_repel(
    data = infl_plot %>%
      filter(is_t4, cooks_flag_4n) %>%
      arrange(desc(cooks_d)) %>%
      head(5),
    aes(label = display_label),
    size = 2.8, color = "black", direction = "both", nudge_y = 0.012,
    box.padding = 0.25, point.padding = 0.1, min.segment.length = 0,
    max.overlaps = Inf, seed = 20250904
  ) +
  scale_fill_identity() +
  labs(
    x = "Sample (ranked by Cook's distance)",
    y = "Cook's distance",
    title = NULL
  ) +
  theme_bw(base_size = FS_AXIS_TEXT) +
  fig6_text_theme() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save standalone Cook's distance figure
ggsave(file.path(fig_out, "FigureS_type4_cooks_distance.png"),
       p_cooks, width = 8, height = 4.5, dpi = 300)
ggsave(file.path(fig_out, "FigureS_type4_cooks_distance.pdf"),
       p_cooks, width = 8, height = 4.5, device = cairo_pdf)
cat(sprintf("Saved figure: %s (.png + .pdf)\n",
            file.path(fig_out, "FigureS_type4_cooks_distance")))

# ==============================================================================
# ANALYSIS 2b (supplementary): Leverage vs residual diagnostic plot
#   Kept as a separate diagnostic figure, not part of Figure S5.
# ==============================================================================

resid_std <- rstandard(fit_ua)

infl_plot2 <- data.frame(
  SampleID   = t4_ref_complete$SampleID,  # R1.1/R2.2 preserve IDs for anonymous labels
  leverage   = as.numeric(hats),
  resid_std  = as.numeric(resid_std),
  cooks_d    = as.numeric(cooks),
  Is_Type4   = t4_ref_complete$Is_Type4,
  UA         = t4_ref_complete$UA,
  stringsAsFactors = FALSE
) %>%
  left_join(t4_anon_map %>% select(SampleID, AnonID), by = "SampleID") %>%
  mutate(
    point_color = case_when(
      cooks_d > 0.5 ~ "#E64B35",
      cooks_d > cook_threshold_4n & Is_Type4 == 1 ~ "#E64B35",
      cooks_d > cook_threshold_4n & Is_Type4 == 0 ~ "#F39B7F",
      Is_Type4 == 1 ~ "#3C5488",
      TRUE ~ "grey60"
    ),
    point_size = pmax(2, pmin(8, cooks_d * 15)),
    display_label = ifelse(Is_Type4 == 1 & cooks_d > cook_threshold_4n,
                           AnonID, NA_character_)
  )

p_lev <- ggplot(infl_plot2, aes(x = leverage, y = resid_std)) +
  geom_point(aes(color = point_color, size = point_size), alpha = 0.7) +
  scale_color_identity() +
  scale_size_identity() +
  geom_vline(xintercept = hat_threshold, linetype = "dashed",
             color = "#E64B35", linewidth = 0.5) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed",
             color = "grey50", linewidth = 0.3) +
  annotate("text", x = hat_threshold, y = -2.5,
           label = sprintf("2p/n = %.4f", hat_threshold),
           size = 3.5, color = "#E64B35", hjust = 1.1, fontface = "italic") +
  geom_text(
    data = infl_plot2 %>%
      filter(Is_Type4 == 1, cooks_d > cook_threshold_4n) %>%
      arrange(desc(cooks_d)) %>%
      head(5),
    aes(label = display_label),
    size = 2.8, color = "black", vjust = -0.8
  ) +
  labs(
    x = "Leverage (hat value)",
    y = "Standardized residual",
    title = NULL,
    color = NULL
  ) +
  theme_bw(base_size = FS_AXIS_TEXT) +
  fig6_text_theme() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Save leverage diagnostic as standalone figure
ggsave(file.path(fig_out, "FigureS_type4_leverage_diagnostic.png"),
       p_lev, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(fig_out, "FigureS_type4_leverage_diagnostic.pdf"),
       p_lev, width = 7, height = 5.5, device = cairo_pdf)
cat(sprintf("Saved figure: %s (.png + .pdf)\n",
            file.path(fig_out, "FigureS_type4_leverage_diagnostic")))


# ==============================================================================
# ANALYSIS 3: Outlier removal sensitivity analysis
# ==============================================================================

cat("\n================================================================\n")
cat("ANALYSIS 3: Outlier Removal Sensitivity\n")
cat("================================================================\n\n")

# --- Define outliers within Type 4 ---
t4_ua <- t4_all$UA[!is.na(t4_all$UA)]
p95_threshold <- quantile(t4_ua, 0.95)
p99_threshold <- quantile(t4_ua, 0.99)

cat(sprintf("Type 4 UA percentiles:\n"))
cat(sprintf("  P95 = %.1f μmol/L\n", p95_threshold))
cat(sprintf("  P99 = %.1f μmol/L\n", p99_threshold))

# Identify outliers
t4_outliers_p95 <- t4_all %>% filter(UA > p95_threshold) %>% pull(SampleID)
t4_outliers_p99 <- t4_all %>% filter(UA > p99_threshold) %>% pull(SampleID)

cat(sprintf("\nOutliers (UA > P95 = %.1f): %d samples\n",
            p95_threshold, length(t4_outliers_p95)))
cat(sprintf("  SampleIDs: %s\n", paste(t4_outliers_p95, collapse = ", ")))
cat(sprintf("Outliers (UA > P99 = %.1f): %d samples\n",
            p99_threshold, length(t4_outliers_p99)))
if (length(t4_outliers_p99) > 0) {
  cat(sprintf("  SampleIDs: %s\n", paste(t4_outliers_p99, collapse = ", ")))
}

# --- Helper: run Type 4 vs Ref comparison ---
run_t4_vs_ref <- function(data_subset, label) {
  d <- data_subset %>%
    filter(!is.na(UA), !is.na(Is_Type4), !is.na(Gender), !is.na(Age))

  t4_ua_vals  <- d$UA[d$Is_Type4 == 1]
  ref_ua_vals <- d$UA[d$Is_Type4 == 0]

  # Wilcoxon test
  wt <- wilcox.test(t4_ua_vals, ref_ua_vals, exact = FALSE)

  # Mean / median difference
  mean_diff  <- mean(t4_ua_vals) - mean(ref_ua_vals)
  median_diff <- median(t4_ua_vals) - median(ref_ua_vals)

  # Linear model (Type 4 coefficient)
  fit <- lm(UA ~ Is_Type4 + Gender + Age, data = d)
  coef_t4 <- coef(fit)["Is_Type4"]
  se_t4   <- summary(fit)$coefficients["Is_Type4", "Std. Error"]
  p_t4    <- summary(fit)$coefficients["Is_Type4", "Pr(>|t|)"]

  # Logistic regression for HUA
  d_hua <- d %>% filter(!is.na(HUA))
  if (sum(d_hua$HUA[d_hua$Is_Type4 == 1]) >= 3 &&
      sum(d_hua$HUA[d_hua$Is_Type4 == 0]) >= 3) {
    fit_logit <- tryCatch(
      glm(HUA ~ Is_Type4 + Gender + Age, data = d_hua,
          family = binomial(link = "logit")),
      error = function(e) NULL
    )
    if (!is.null(fit_logit)) {
      or_t4  <- exp(coef(fit_logit)["Is_Type4"])
      # Use Wald CI (confint.default) instead of profile likelihood (confint)
      # to avoid extremely slow profiling in LOO/bootstrap iterations
      ci_wald <- suppressMessages(confint.default(fit_logit)["Is_Type4", ])
      ci_lo  <- exp(ci_wald[1])
      ci_hi  <- exp(ci_wald[2])
      p_or   <- summary(fit_logit)$coefficients["Is_Type4", "Pr(>|z|)"]
    } else {
      or_t4 <- ci_lo <- ci_hi <- p_or <- NA
    }
  } else {
    or_t4 <- ci_lo <- ci_hi <- p_or <- NA
  }

  data.frame(
    scenario     = label,
    n_t4         = length(t4_ua_vals),
    n_ref        = length(ref_ua_vals),
    mean_UA_t4   = round(mean(t4_ua_vals), 1),
    mean_UA_ref  = round(mean(ref_ua_vals), 1),
    mean_diff    = round(mean_diff, 2),
    median_diff  = round(median_diff, 2),
    wilcox_W     = wt$statistic,
    wilcox_p     = round(wt$p.value, 5),
    lm_coef_t4   = round(coef_t4, 3),
    lm_se_t4     = round(se_t4, 3),
    lm_p_t4      = round(p_t4, 5),
    HUA_OR       = round(or_t4, 3),
    HUA_OR_CI_lo = round(ci_lo, 3),
    HUA_OR_CI_hi = round(ci_hi, 3),
    HUA_OR_p     = round(p_or, 5),
    stringsAsFactors = FALSE
  )
}

# --- Run scenarios ---
sens_full <- run_t4_vs_ref(t4_ref, "Full (no removal)")
sens_p95  <- run_t4_vs_ref(
  t4_ref %>% filter(!(SampleID %in% t4_outliers_p95)),
  sprintf("Remove P95 outliers (UA > %.1f)", p95_threshold)
)
sens_p99  <- run_t4_vs_ref(
  t4_ref %>% filter(!(SampleID %in% t4_outliers_p99)),
  sprintf("Remove P99 outliers (UA > %.1f)", p99_threshold)
)

sensitivity_table <- bind_rows(sens_full, sens_p95, sens_p99)

cat("\n--- Sensitivity analysis results ---\n")
print(sensitivity_table, row.names = FALSE)

# --- Outlier details ---
outlier_details <- t4_all %>%
  filter(SampleID %in% t4_outliers_p95) %>%
  select(SampleID, UA, Gender, Age, BMI, HUA) %>%
  arrange(desc(UA)) %>%
  mutate(
    above_P95 = TRUE,
    above_P99 = UA > p99_threshold
  )

cat("\n--- Outlier details (Type 4, UA > P95) ---\n")
print(outlier_details, row.names = FALSE)


# ==============================================================================
# ANALYSIS 4: Leave-one-out (LOO) stability
# ==============================================================================

cat("\n================================================================\n")
cat("ANALYSIS 4: Leave-One-Out (LOO) Stability\n")
cat("================================================================\n\n")

# --- Original (full) effect for reference ---
loo_baseline <- run_t4_vs_ref(t4_ref, "Baseline (all Type 4)")
loo_mean_diff_full  <- loo_baseline$mean_diff
loo_wilcox_p_full   <- loo_baseline$wilcox_p
loo_lm_coef_full    <- loo_baseline$lm_coef_t4
loo_lm_p_full       <- loo_baseline$lm_p_t4

cat(sprintf("Baseline: mean_diff = %.2f, wilcox_p = %.5f, lm_coef = %.3f, lm_p = %.5f\n",
            loo_mean_diff_full, loo_wilcox_p_full, loo_lm_coef_full, loo_lm_p_full))

# --- LOO: remove each Type 4 sample in turn ---
loo_results <- list()
t4_sample_ids <- t4_all$SampleID

cat(sprintf("Running LOO for %d Type 4 samples...\n", length(t4_sample_ids)))

for (i in seq_along(t4_sample_ids)) {
  sid <- t4_sample_ids[i]
  loo_subset <- t4_ref %>% filter(SampleID != sid)
  res <- run_t4_vs_ref(loo_subset, sid)

  # Track which sample was removed
  res$removed_SampleID <- sid
  res$removed_UA <- t4_all$UA[t4_all$SampleID == sid]
  res$removed_HUA <- t4_all$HUA[t4_all$SampleID == sid]

  # Change from baseline
  res$delta_mean_diff <- res$mean_diff - loo_mean_diff_full
  res$delta_lm_coef   <- res$lm_coef_t4 - loo_lm_coef_full

  loo_results[[i]] <- res

  if (i %% 10 == 0) cat(sprintf("  ...%d/%d done\n", i, length(t4_sample_ids)))
}

loo_table <- bind_rows(loo_results) %>%
  arrange(desc(abs(delta_mean_diff)))

cat("LOO complete.\n\n")

# --- LOO summary ---
loo_p_sig_count <- sum(loo_table$wilcox_p < 0.05)
loo_p_sig_prop  <- loo_p_sig_count / nrow(loo_table)
loo_lm_p_sig_count <- sum(loo_table$lm_p_t4 < 0.05)
loo_lm_p_sig_prop  <- loo_lm_p_sig_count / nrow(loo_table)

cat("--- LOO Summary ---\n")
cat(sprintf("Wilcoxon p < 0.05: %d / %d (%.1f%%)\n",
            loo_p_sig_count, nrow(loo_table), loo_p_sig_prop * 100))
cat(sprintf("LM Type4 coef p < 0.05: %d / %d (%.1f%%)\n",
            loo_lm_p_sig_count, nrow(loo_table), loo_lm_p_sig_prop * 100))
cat(sprintf("\nMean effect size (mean_diff): %.2f to %.2f\n",
            min(loo_table$mean_diff), max(loo_table$mean_diff)))
cat(sprintf("Mean effect size (lm_coef):   %.3f to %.3f\n",
            min(loo_table$lm_coef_t4), max(loo_table$lm_coef_t4)))

cat("\n--- Top 5 most influential samples (by mean_diff change) ---\n")
loo_table %>%
  select(removed_SampleID, removed_UA, removed_HUA,
         mean_diff, delta_mean_diff, wilcox_p, lm_coef_t4, lm_p_t4) %>%
  head(5) %>%
  print(row.names = FALSE)

cat("\n--- Bottom 5 (least influential) ---\n")
loo_table %>%
  tail(5) %>%
  select(removed_SampleID, removed_UA, removed_HUA,
         mean_diff, delta_mean_diff, wilcox_p, lm_coef_t4, lm_p_t4) %>%
  print(row.names = FALSE)


# ==============================================================================
# ANALYSIS 4 (figure): LOO stability plot (→ Figure S5, Panel B)
# ==============================================================================

# Merge anonymous IDs and build y-axis labels (AnonID + UA, no SampleID)
loo_plot_df <- loo_table %>%
  left_join(t4_anon_map %>% select(SampleID = SampleID, AnonID),
            by = c("removed_SampleID" = "SampleID")) %>%
  mutate(
    removed_label = sprintf("%s (UA=%.0f)", AnonID, removed_UA),
    is_sig = wilcox_p < 0.05,
    baseline_mean_diff = loo_mean_diff_full
  ) %>%
  arrange(mean_diff) %>%
  mutate(removed_label = factor(removed_label, levels = removed_label))

p_loo <- ggplot(loo_plot_df, aes(x = mean_diff, y = removed_label,
                                 color = is_sig)) +
  geom_vline(xintercept = loo_mean_diff_full, linetype = "dashed",
             color = "grey40", linewidth = 0.5) +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = baseline_mean_diff, xend = mean_diff,
                   yend = removed_label), linewidth = 0.6) +
  geom_point(aes(size = abs(delta_mean_diff)), shape = 16) +
  annotate("text", x = loo_mean_diff_full, y = nrow(loo_plot_df) + 0.7,
           label = sprintf("Baseline\nDiff=%.1f", loo_mean_diff_full),
           size = 3.5, color = "grey40", hjust = 0.5) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#3C5488"),
                     name = "Wilcoxon P < 0.05",
                     labels = c("TRUE" = "Yes", "FALSE" = "No")) +
  scale_size_continuous(name = "|change in mean diff|", range = c(2, 6)) +
  labs(
    x = "Mean UA difference (Type 4 - Ref), with one Type 4 sample removed",
    y = "Removed sample (UA in umol/L)",
    title = NULL
  ) +
  theme_bw(base_size = FS_AXIS_TEXT) +
  fig6_text_theme() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position = "right",
    legend.box = "vertical"
  )

# Save LOO figure
ggsave(file.path(fig_out, "FigureS_type4_loo_stability.png"),
       p_loo, width = 8.5, height = 9, dpi = 300)
ggsave(file.path(fig_out, "FigureS_type4_loo_stability.pdf"),
       p_loo, width = 8.5, height = 9, device = cairo_pdf)
cat(sprintf("\nSaved figure: %s (.png + .pdf)\n",
            file.path(fig_out, "FigureS_type4_loo_stability")))


# ==============================================================================
# ANALYSIS 5: Bootstrap stability (stratified resampling)
#   Both Type 4 and Reference groups are resampled independently,
#   reflecting uncertainty in both groups rather than conditioning on
#   the observed reference sample.
# ==============================================================================

cat("\n================================================================\n")
cat("ANALYSIS 5: Bootstrap Stability (n_boot = 5000, stratified)\n")
cat("================================================================\n\n")

n_boot <- 5000
t4_indices <- which(t4_ref$Is_Type4 == 1)
ref_data   <- t4_ref %>% filter(Is_Type4 == 0)
t4_data    <- t4_ref %>% filter(Is_Type4 == 1)

cat(sprintf("Type 4: n = %d, Reference (Type 2+3): n = %d\n",
            nrow(t4_data), nrow(ref_data)))
cat(sprintf("Bootstrapping (n_boot = %d, stratified resampling)...\n", n_boot))

boot_results <- data.frame(
  boot_iter      = integer(n_boot),
  mean_diff      = numeric(n_boot),
  wilcox_p       = numeric(n_boot),
  lm_coef_t4     = numeric(n_boot),
  lm_p_t4        = numeric(n_boot),
  HUA_OR         = numeric(n_boot),
  HUA_OR_p       = numeric(n_boot),
  stringsAsFactors = FALSE
)

for (b in 1:n_boot) {
  boot_results$boot_iter[b] <- b  # R1.1/R2.2 retain auditable iteration IDs
  # Stratified resampling: both groups independently resampled with replacement
  boot_idx_t4  <- sample(seq_len(nrow(t4_data)), replace = TRUE)
  t4_boot      <- t4_data[boot_idx_t4, ]
  boot_idx_ref <- sample(seq_len(nrow(ref_data)), replace = TRUE)
  ref_boot     <- ref_data[boot_idx_ref, ]

  boot_combined <- bind_rows(t4_boot, ref_boot) %>%
    mutate(Is_Type4 = ifelse(Enterotype == "4", 1, 0))

  # Mean difference
  boot_results$mean_diff[b] <- mean(t4_boot$UA) - mean(ref_boot$UA)

  # Wilcoxon (with possible ties from resampling)
  wt_boot <- tryCatch(
    wilcox.test(t4_boot$UA, ref_boot$UA, exact = FALSE),
    error = function(e) list(statistic = NA, p.value = NA)
  )
  boot_results$wilcox_p[b] <- wt_boot$p.value

  # Linear model (adjusted: Gender + Age)
  fit_boot <- tryCatch(
    lm(UA ~ Is_Type4 + Gender + Age, data = boot_combined),
    error = function(e) NULL
  )
  if (!is.null(fit_boot)) {
    boot_results$lm_coef_t4[b] <- coef(fit_boot)["Is_Type4"]
    boot_results$lm_p_t4[b]    <- summary(fit_boot)$coefficients["Is_Type4", "Pr(>|t|)"]
  } else {
    boot_results$lm_coef_t4[b] <- NA
    boot_results$lm_p_t4[b]    <- NA
  }

  # Logistic regression for HUA (adjusted: Gender + Age)
  boot_hua <- boot_combined %>% filter(!is.na(HUA))
  if (sum(boot_hua$HUA[boot_hua$Is_Type4 == 1]) >= 2 &&
      sum(boot_hua$HUA[boot_hua$Is_Type4 == 1] == 0) >= 2) {
    fit_logit_boot <- tryCatch(
      glm(HUA ~ Is_Type4 + Gender + Age, data = boot_hua,
          family = binomial(link = "logit")),
      error = function(e) NULL
    )
    if (!is.null(fit_logit_boot)) {
      boot_results$HUA_OR[b]   <- exp(coef(fit_logit_boot)["Is_Type4"])
      boot_results$HUA_OR_p[b] <- summary(fit_logit_boot)$coefficients["Is_Type4", "Pr(>|z|)"]
    } else {
      boot_results$HUA_OR[b]   <- NA
      boot_results$HUA_OR_p[b] <- NA
    }
  } else {
    boot_results$HUA_OR[b]   <- NA
    boot_results$HUA_OR_p[b] <- NA
  }

  if (b %% 500 == 0) cat(sprintf("  ...%d/%d done\n", b, n_boot))
}

cat("Bootstrap complete.\n\n")

# --- Bootstrap summary ---
boot_mean_diff_ci  <- quantile(boot_results$mean_diff, c(0.025, 0.975), na.rm = TRUE)
boot_lm_coef_ci    <- quantile(boot_results$lm_coef_t4, c(0.025, 0.975), na.rm = TRUE)
boot_wilcox_sig    <- mean(boot_results$wilcox_p < 0.05, na.rm = TRUE)
boot_lm_sig        <- mean(boot_results$lm_p_t4 < 0.05, na.rm = TRUE)
boot_or_sig        <- mean(boot_results$HUA_OR_p < 0.05, na.rm = TRUE)
boot_or_ci         <- quantile(boot_results$HUA_OR, c(0.025, 0.975), na.rm = TRUE)

boot_summary <- data.frame(
  metric = c("Mean UA difference (T4 − Ref)",
             "LM Type 4 coefficient",
             "Wilcoxon P < 0.05 proportion",
             "LM P < 0.05 proportion",
             "HUA OR",
             "HUA OR P < 0.05 proportion"),
  estimate = c(
    round(loo_baseline$mean_diff, 2),
    round(loo_baseline$lm_coef_t4, 3),
    round(boot_wilcox_sig, 3),
    round(boot_lm_sig, 3),
    round(loo_baseline$HUA_OR, 3),
    round(boot_or_sig, 3)
  ),
  bootstrap_CI_lo = c(
    round(boot_mean_diff_ci[1], 2),
    round(boot_lm_coef_ci[1], 3),
    NA,
    NA,
    round(boot_or_ci[1], 3),
    NA
  ),
  bootstrap_CI_hi = c(
    round(boot_mean_diff_ci[2], 2),
    round(boot_lm_coef_ci[2], 3),
    NA,
    NA,
    round(boot_or_ci[2], 3),
    NA
  ),
  n_valid_boot = c(
    sum(!is.na(boot_results$mean_diff)),
    sum(!is.na(boot_results$lm_coef_t4)),
    sum(!is.na(boot_results$wilcox_p)),
    sum(!is.na(boot_results$lm_p_t4)),
    sum(!is.na(boot_results$HUA_OR)),
    sum(!is.na(boot_results$HUA_OR_p))
  )
)

cat("--- Bootstrap Summary ---\n")
print(boot_summary, row.names = FALSE)

cat(sprintf("\nMean UA diff 95%% CI: [%.2f, %.2f]\n",
            boot_mean_diff_ci[1], boot_mean_diff_ci[2]))
cat(sprintf("LM Type4 coef 95%% CI: [%.3f, %.3f]\n",
            boot_lm_coef_ci[1], boot_lm_coef_ci[2]))
cat(sprintf("HUA OR 95%% CI: [%.3f, %.3f]\n",
            boot_or_ci[1], boot_or_ci[2]))
cat(sprintf("Wilcoxon P < 0.05: %.1f%%\n", boot_wilcox_sig * 100))
cat(sprintf("LM P < 0.05:       %.1f%%\n", boot_lm_sig * 100))
cat(sprintf("HUA OR P < 0.05:   %.1f%%\n", boot_or_sig * 100))

# --- Save bootstrap results ---
write.csv(boot_summary, file.path(tbl_out, "type4_bootstrap_stability.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved: %s\n", file.path(tbl_out, "type4_bootstrap_stability.csv")))

# Also save raw bootstrap iterations for reference
boot_raw_file <- file.path(tbl_out, "type4_bootstrap_raw_iterations.csv")
write.csv(boot_results, boot_raw_file, row.names = FALSE)
cat(sprintf("Saved raw iterations: %s\n", boot_raw_file))


# ==============================================================================
# Save combined LOO + sensitivity table
# ==============================================================================

# The sensitivity (outlier removal) and LOO results are conceptually related.
# Save LOO results to the sensitivity_loo file.
loo_save <- loo_table %>%
  select(removed_SampleID, removed_UA, removed_HUA,
         n_t4, n_ref, mean_diff, median_diff,
         wilcox_W, wilcox_p, lm_coef_t4, lm_p_t4,
         HUA_OR, HUA_OR_p, delta_mean_diff, delta_lm_coef)

write.csv(loo_save, file.path(tbl_out, "type4_sensitivity_loo.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved LOO results: %s\n",
            file.path(tbl_out, "type4_sensitivity_loo.csv")))

# Also save outlier removal sensitivity as a separate section in the same file
# by appending to the characterization table
sens_save <- sensitivity_table %>%
  mutate(outlier_SampleIDs_p95 = paste(t4_outliers_p95, collapse = "; "),
         outlier_SampleIDs_p99 = paste(t4_outliers_p99, collapse = "; "))
sens_file <- file.path(tbl_out, "type4_outlier_sensitivity.csv")
write.csv(sens_save, sens_file, row.names = FALSE)
cat(sprintf("Saved outlier sensitivity: %s\n", sens_file))


# ==============================================================================
# COMBINED FIGURE S5: Cook's distance (A) + LOO stability (B)
# ==============================================================================

library(patchwork)

p_s5 <- (p_cooks / p_loo) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(fig_out, "FigureS3_type4_diagnostics.png"),
       p_s5, width = 10, height = 13, dpi = 300)
ggsave(file.path(fig_out, "FigureS3_type4_diagnostics.pdf"),
       p_s5, width = 10, height = 13, device = cairo_pdf)
cat(sprintf("\nSaved combined Figure S5: %s (.png + .pdf)\n",
            file.path(fig_out, "FigureS3_type4_diagnostics")))


# ==============================================================================
# FINAL SUMMARY
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("FINAL SUMMARY — Type 4 Detailed Characterization\n")
cat("================================================================\n\n")

cat("1. SUBGROUP DESCRIPTION:\n")
cat(sprintf("   Type 4 total: n = %d (HUA+ = %d, HUA- = %d)\n",
            nrow(t4_all), sum(t4_all$HUA == 1, na.rm = TRUE),
            sum(t4_all$HUA == 0, na.rm = TRUE)))
cat(sprintf("   Type 4 HUA+ UA: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(t4_all$UA[t4_all$HUA == 1], na.rm = TRUE),
            sd(t4_all$UA[t4_all$HUA == 1], na.rm = TRUE),
            min(t4_all$UA[t4_all$HUA == 1], na.rm = TRUE),
            max(t4_all$UA[t4_all$HUA == 1], na.rm = TRUE)))
cat(sprintf("   Type 4 HUA- UA: mean = %.1f, SD = %.1f, range = [%.1f, %.1f]\n",
            mean(t4_all$UA[t4_all$HUA == 0], na.rm = TRUE),
            sd(t4_all$UA[t4_all$HUA == 0], na.rm = TRUE),
            min(t4_all$UA[t4_all$HUA == 0], na.rm = TRUE),
            max(t4_all$UA[t4_all$HUA == 0], na.rm = TRUE)))

cat("\n2. INFLUENCE DIAGNOSTICS:\n")
cat(sprintf("   Cook's > 4/n (%.4f): %d flagged (%d in Type 4)\n",
            cook_threshold_4n, sum(infl_df$cooks_flag_4n),
            sum(t4_infl$cooks_flag_4n)))
cat(sprintf("   Cook's > 0.5: %d flagged\n", sum(infl_df$cooks_flag_0.5)))
cat(sprintf("   Leverage > 2p/n: %d flagged\n", sum(infl_df$leverage_flag)))

cat("\n3. UPPER-TAIL SENSITIVITY (R1.1/R2.2):\n")
cat(sprintf("   P95 threshold = %.1f, %d observations above threshold\n", p95_threshold,
            length(t4_outliers_p95)))
cat(sprintf("   Full:       mean_diff = %.1f, wilcox_p = %.4f, OR = %.2f\n",
            sens_full$mean_diff, sens_full$wilcox_p, sens_full$HUA_OR))
cat(sprintf("   Remove P95: mean_diff = %.1f, wilcox_p = %.4f, OR = %.2f\n",
            sens_p95$mean_diff, sens_p95$wilcox_p, sens_p95$HUA_OR))
cat(sprintf("   Remove P99: mean_diff = %.1f, wilcox_p = %.4f, OR = %.2f\n",
            sens_p99$mean_diff, sens_p99$wilcox_p, sens_p99$HUA_OR))

cat("\n4. LOO STABILITY:\n")
cat(sprintf("   Wilcoxon P < 0.05 in %d/%d LOO iterations (%.1f%%)\n",
            loo_p_sig_count, nrow(loo_table), loo_p_sig_prop * 100))
cat(sprintf("   LM P < 0.05 in %d/%d LOO iterations (%.1f%%)\n",
            loo_lm_p_sig_count, nrow(loo_table), loo_lm_p_sig_prop * 100))
cat(sprintf("   Effect size range: [%.1f, %.1f]\n",
            min(loo_table$mean_diff), max(loo_table$mean_diff)))

cat(sprintf("\n5. BOOTSTRAP STABILITY (n = %d):\n", n_boot))
cat(sprintf("   Mean diff 95%% CI: [%.1f, %.1f]\n",
            boot_mean_diff_ci[1], boot_mean_diff_ci[2]))
cat(sprintf("   LM coef 95%% CI: [%.2f, %.2f]\n",
            boot_lm_coef_ci[1], boot_lm_coef_ci[2]))
cat(sprintf("   Wilcoxon P < 0.05: %.1f%%\n", boot_wilcox_sig * 100))
cat(sprintf("   LM P < 0.05: %.1f%%\n", boot_lm_sig * 100))
cat(sprintf("   HUA OR P < 0.05: %.1f%%\n", boot_or_sig * 100))

cat("\n================================================================\n")
cat("All analyses complete.\n")
cat("================================================================\n")
