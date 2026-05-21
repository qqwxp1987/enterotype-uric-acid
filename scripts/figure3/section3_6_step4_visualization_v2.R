# ============================================================================== 
# Section 3.6 — Step 4: Visualization v2 (Now exported as Figure 3)
#
# Purpose:
#   Robustness figure with three panels:
#     A) Model Evidence Trajectory — shows overall enterotype P across models
#     B) Cardiometabolic Model Full Coefficients — dual forest plot
#     C) Subgroup Evidence Matrix — improved heatmap with significance overlay
# ==============================================================================

library(tidyverse)
library(here)
library(patchwork)

cat("========================================\n")
cat("Section 3.6 Step 4: Visualization v2 (exporting Figure 3)\n")
cat("========================================\n\n")

# NOTE: the robustness-model inputs physically live under section_3_3.
tab_dir <- here("results", "tables", "section_3_3")
fig_dir <- here("results", "figures", "panels")
out_dir <- here("results", "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
model_results <- readRDS(file.path(tab_dir, "section3_6_model_results.rds"))
wls_df <- model_results$wls
logit_df <- model_results$logit
subgroup_df <- readRDS(file.path(tab_dir, "subgroup_results.rds"))

# Common theme elements
theme_fig6 <- theme_classic(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    axis.text = element_text(color = "black"),
    plot.margin = margin(8, 12, 8, 8)
  )

# Color palette
col_ua  <- "#2166AC"   # blue for UA
col_hua <- "#B2182B"   # red for HUA
col_sig <- "#B2182B"   # significant
col_ns  <- "grey55"    # non-significant
col_enterotype <- "#D6604D"  # enterotype terms in forest
col_covariate  <- "grey40"   # clinical covariate terms in forest

# ============================================================================
# Panel A: Model Evidence Trajectory
# ============================================================================

cat("Creating Panel A: Model Evidence Trajectory...\n")

# Extract overall enterotype P from each unique model
# Use the WLS data for UA and logit data for HUA
model_order <- c("M1_Demographic", "M2_BMI", "M3_Renal_Cr", "M4_Cardiomet",
                 "S1_Renal_eGFR", "S2_Full_Clinical")
model_labels_short <- c(
  M1_Demographic = "M1\nDemographic",
  M2_BMI = "M2\n+BMI",
  M3_Renal_Cr = "M3\n+Renal",
  M4_Cardiomet = "M4\n+Cardiomet",
  S1_Renal_eGFR = "S1\n+eGFR",
  S2_Full_Clinical = "S2\nFull Clinical"
)
model_n <- c(
  M1_Demographic = 305, M2_BMI = 225, M3_Renal_Cr = 305,
  M4_Cardiomet = 305, S1_Renal_eGFR = 305, S2_Full_Clinical = 223
)

# Get unique overall enterotype P per model per outcome
p_ua_summary <- wls_df %>%
  filter(model %in% model_order) %>%
  select(model, overall_enterotype_p) %>%
  distinct() %>%
  mutate(outcome = "UA")

p_hua_summary <- logit_df %>%
  filter(model %in% model_order) %>%
  select(model, overall_enterotype_p) %>%
  distinct() %>%
  mutate(outcome = "HUA")

traj_df <- bind_rows(p_ua_summary, p_hua_summary) %>%
  mutate(
    model_f = factor(model, levels = model_order),
    model_idx = as.numeric(model_f),
    neg_log10_p = -log10(overall_enterotype_p),
    n = model_n[model],
    sig = overall_enterotype_p < 0.05,
    p_label = ifelse(overall_enterotype_p < 0.001,
                     sprintf("%.1e", overall_enterotype_p),
                     ifelse(overall_enterotype_p < 0.01,
                            sub("0+$", "", sprintf("%.4f", overall_enterotype_p)),
                            sprintf("%.3f", overall_enterotype_p)))
  )

panel_A <- ggplot(traj_df, aes(x = model_idx, y = neg_log10_p, color = outcome, group = outcome)) +
  # Reference line at P = 0.05
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey60", linewidth = 0.5) +
  annotate("text", x = 6.45, y = -log10(0.05), label = "P = 0.05",
           color = "grey50", size = 2.8, hjust = 0, vjust = -0.4) +
  # Lines and points
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_point(aes(shape = sig), size = 3, stroke = 0.8) +
  # P-value labels above/below points
  geom_text(aes(label = p_label),
            size = 2.4, vjust = -1.2, show.legend = FALSE) +
  # Sample size annotation at bottom
  annotate("text",
           x = 1:6,
           y = rep(-0.15, 6),
           label = paste0("n=", model_n[model_order]),
           size = 2.6, color = "grey40") +
  # Scales
  scale_x_continuous(
    breaks = 1:6,
    labels = model_labels_short[model_order],
    expand = expansion(mult = c(0.06, 0.12))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.12, 0.15)),
    name = expression(-log[10](P[enterotype]))
  ) +
  scale_color_manual(
    values = c(UA = col_ua, HUA = col_hua),
    name = "Outcome"
  ) +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1),
    name = "P < 0.05",
    labels = c(`TRUE` = "Yes", `FALSE` = "No")
  ) +
  labs(title = "A    Overall enterotype significance across adjustment models",
       x = NULL) +
  theme_fig6 +
  theme(
    legend.position = c(0.92, 0.85),
    legend.background = element_rect(fill = "white", color = NA),
    legend.key.size = unit(0.35, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 8.5, face = "bold"),
    axis.text.x = element_text(size = 8, lineheight = 0.9)
  )


# ============================================================================
# Panel B: Cardiometabolic Model (M4) Full Coefficient Forest Plot
# ============================================================================

cat("Creating Panel B: Cardiometabolic Model Forest Plot...\n")

# Extract M4_Cardiomet coefficients
m4_ua <- wls_df %>%
  filter(model == "M4_Cardiomet") %>%
  select(label, estimate, ci_low, ci_high, p) %>%
  mutate(outcome = "UA")

m4_hua <- logit_df %>%
  filter(model == "M4_Cardiomet") %>%
  select(label, estimate, ci_low, ci_high, p) %>%
  mutate(outcome = "HUA")

# Define term order and grouping
enterotype_terms <- c("Type 2 vs Type 1", "Type 3 vs Type 1", "Type 4 vs Type 1")
covariate_terms  <- c("Sex (Male vs Female)", "Age (per year)",
                       "Creatinine", "Urea", "Glucose",
                       "Triglycerides", "HDL-C", "LDL-C")
all_terms <- c(enterotype_terms, covariate_terms)

# UA panel (beta scale)
forest_ua <- m4_ua %>%
  filter(label %in% all_terms) %>%
  mutate(
    label = factor(label, levels = rev(all_terms)),
    group = ifelse(label %in% enterotype_terms, "Enterotype", "Covariate"),
    sig = p < 0.05
  )

pB_ua <- ggplot(forest_ua, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  # Shading for enterotype terms
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = length(all_terms) - length(enterotype_terms) + 0.5,
           ymax = length(all_terms) + 0.5,
           fill = "#FFF5F0", alpha = 0.6) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, color = sig),
                width = 0.2, linewidth = 0.6) +
  geom_point(aes(color = sig, shape = sig), size = 2.5) +
  scale_color_manual(values = c(`TRUE` = col_sig, `FALSE` = col_ns), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = expression(beta ~ "(95% CI)"), y = NULL,
       subtitle = "Continuous UA") +
  theme_fig6 +
  theme(axis.text.y = element_text(size = 8.5))

# HUA panel (OR scale)
forest_hua <- m4_hua %>%
  filter(label %in% all_terms) %>%
  mutate(
    label = factor(label, levels = rev(all_terms)),
    group = ifelse(label %in% enterotype_terms, "Enterotype", "Covariate"),
    sig = p < 0.05
  )

pB_hua <- ggplot(forest_hua, aes(x = estimate, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = length(all_terms) - length(enterotype_terms) + 0.5,
           ymax = length(all_terms) + 0.5,
           fill = "#FFF5F0", alpha = 0.6) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high, color = sig),
                width = 0.2, linewidth = 0.6) +
  geom_point(aes(color = sig, shape = sig), size = 2.5) +
  scale_color_manual(values = c(`TRUE` = col_sig, `FALSE` = col_ns), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_x_continuous(trans = "log2",
                     breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4),
                     labels = c("0.06", "0.13", "0.25", "0.5", "1", "2", "4")) +
  labs(x = "Odds ratio (95% CI)", y = NULL,
       subtitle = "Binary HUA") +
  theme_fig6 +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

panel_B <- pB_ua + pB_hua +
  plot_layout(widths = c(1.3, 1)) +
  plot_annotation(
    title = "B    Cardiometabolic model coefficients (M4, n = 305)",
    theme = theme(plot.title = element_text(face = "bold", size = 11))
  )


# ============================================================================
# Panel C: Subgroup Evidence Matrix (improved)
# ============================================================================

cat("Creating Panel C: Subgroup Evidence Matrix...\n")

contrast_order <- c("Type 2 vs Type 1", "Type 3 vs Type 1", "Type 4 vs Type 1")
subset_levels <- c("eGFR_normal", "Non_obese", "Female", "Male")
subset_labels <- c(
  Male = "Male (n=81)",
  Female = "Female (n=224)",
  Non_obese = "Non-obese (n=218)",
  eGFR_normal = "eGFR \u2265 90 (n=271)"
)

theme_heatmap <- theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 25, hjust = 1, size = 8.5),
    axis.text.y = element_text(size = 9, color = "black"),
    panel.grid = element_blank(),
    legend.key.height = unit(0.5, "cm"),
    legend.key.width = unit(0.3, "cm"),
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 8),
    plot.margin = margin(4, 6, 4, 4)
  )

# Separate HUA and UA subgroup data for independent color scales
sub_hua <- subgroup_df %>%
  filter(label %in% contrast_order, outcome == "HUA") %>%
  mutate(
    subset_f = factor(subset, levels = subset_levels, labels = subset_labels[subset_levels]),
    label_f = factor(label, levels = contrast_order),
    effect_std = log(estimate),
    sig_level = case_when(
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.10  ~ "\u2020",
      TRUE      ~ ""
    )
  )

sub_ua <- subgroup_df %>%
  filter(label %in% contrast_order, outcome == "UA") %>%
  mutate(
    subset_f = factor(subset, levels = subset_levels, labels = subset_labels[subset_levels]),
    label_f = factor(label, levels = contrast_order),
    effect_std = estimate,
    sig_level = case_when(
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      p < 0.10  ~ "\u2020",
      TRUE      ~ ""
    )
  )

pC_hua <- ggplot(sub_hua, aes(x = label_f, y = subset_f)) +
  geom_tile(aes(fill = effect_std), color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_level), size = 4, fontface = "bold", vjust = 0.35) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "log(OR)"
  ) +
  labs(subtitle = "Binary HUA (OR)", x = NULL, y = NULL) +
  theme_heatmap

pC_ua <- ggplot(sub_ua, aes(x = label_f, y = subset_f)) +
  geom_tile(aes(fill = effect_std), color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_level), size = 4, fontface = "bold", vjust = 0.35) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = expression(beta)
  ) +
  labs(subtitle = expression("Continuous UA (" * beta * ")"), x = NULL, y = NULL) +
  theme_heatmap +
  theme(axis.text.y = element_blank())

panel_C <- pC_hua + pC_ua +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    title = "C    Subgroup consistency of enterotype contrasts",
    caption = "* P < 0.05    ** P < 0.01    \u2020 P < 0.10",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 7.5, color = "grey40", hjust = 0)
    )
  )


# ============================================================================
# Assemble Composite Figure 3
# ============================================================================

cat("Assembling composite Figure 3...\n")

# Panel A: top (full width)
# Panel B: middle (full width)
# Panel C: bottom (full width)

fig3_composite <- panel_A / panel_B / panel_C +
  plot_layout(heights = c(1, 1.2, 0.9))

# Save individual panels (panels/, no _v2_ suffix to match existing files)
ggsave(file.path(fig_dir, "Figure3A_Evidence_Trajectory.png"),
       panel_A, width = 8, height = 4, dpi = 300)
ggsave(file.path(fig_dir, "Figure3A_Evidence_Trajectory.pdf"),
       panel_A, width = 8, height = 4)

ggsave(file.path(fig_dir, "Figure3B_Cardiomet_Forest.png"),
       panel_B, width = 9, height = 5, dpi = 300)
ggsave(file.path(fig_dir, "Figure3B_Cardiomet_Forest.pdf"),
       panel_B, width = 9, height = 5)

ggsave(file.path(fig_dir, "Figure3C_Subgroup_Matrix.png"),
       panel_C, width = 8.5, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "Figure3C_Subgroup_Matrix.pdf"),
       panel_C, width = 8.5, height = 4.5)

# Save composite to results/figures/ (consistent with Figure 1/2/4/5)
ggsave(file.path(out_dir, "Figure3_composite.png"),
      fig3_composite, width = 10, height = 13, dpi = 300)
ggsave(file.path(out_dir, "Figure3_composite.tiff"),
      fig3_composite, width = 10, height = 13, dpi = 300, compression = "lzw")
ggsave(file.path(out_dir, "Figure3_composite.pdf"),
      fig3_composite, width = 10, height = 13)

cat("\nAll Figure 3 panels saved successfully.\n")
cat("Panels directory:", fig_dir, "\n")
cat("Composite directory:", out_dir, "\n")
