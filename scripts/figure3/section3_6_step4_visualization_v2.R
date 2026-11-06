# ==============================================================================
# Section 3.6 — Step 4: Visualization v2 (Now exported as Figure 3)
#
# Purpose:
#   Robustness figure with three panels:
#     A) Model Evidence Trajectory — shows overall enterotype P across models
#     B) Cardiometabolic Model Full Coefficients — dual forest plot
#     C) Subgroup Evidence Matrix — improved heatmap with significance overlay
#
# Typography / layout policy (revised):
#   * One shared font-size hierarchy is applied to ALL three panels, so that text
#     at the SAME role (panel title / sub-panel label / axis title / axis tick /
#     legend / in-plot annotation) renders at the SAME size everywhere.
#   * Panel labels "A", "B", "C" are drawn as dedicated title ROWS rather than via
#     patchwork::plot_annotation(). plot_annotation() titles are silently dropped
#     when a sub-patchwork is nested inside the composite — which is why the old
#     composite only showed "A". A title row is a real ggplot and always renders.
#   * Sizes are deliberately generous so the figure stays legible after a journal
#     down-scales it to column / page width.
# ==============================================================================

library(tidyverse)
library(here)
library(patchwork)
library(grid)

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

# ============================================================================
# Shared typography system
# ============================================================================
# A single hierarchy of point sizes, reused by every panel. "Consistent" here
# means: the SAME role gets the SAME size across panels, not that everything is
# one size. Hierarchy (large -> small):
#     panel title  >  sub-panel label / axis title  >  axis tick  >  legend / note
fs_panel_title <- 14     # "A / B / C  ..." panel titles
fs_subtitle    <- 13     # sub-panel labels ("Continuous UA", "Binary HUA", ...)
fs_axis_title  <- 12     # axis titles
fs_axis_text   <- 11     # axis tick labels
fs_legend_text <- 10.5   # legend entries
fs_legend_ttl  <- 11     # legend titles
fs_annot       <- 10.5   # in-plot annotations (p-values, n=, P=0.05, footnote)
fs_star        <- 13     # significance glyphs on the heatmap tiles

# ggplot geom_text()/annotate() take size in mm; convert from points.
PT_PER_MM <- 2.845276    # == ggplot2's .pt
as_mm <- function(pt) pt / PT_PER_MM

# Shared text styling, layered on top of either theme_classic or theme_minimal
# so the two plot families share identical text sizes.
base_text_theme <- theme(
  text          = element_text(color = "black"),
  axis.title    = element_text(size = fs_axis_title, color = "black"),
  axis.text     = element_text(size = fs_axis_text, color = "black"),
  plot.subtitle = element_text(size = fs_subtitle, face = "plain",
                               hjust = 0, margin = margin(b = 5)),
  legend.text   = element_text(size = fs_legend_text),
  legend.title  = element_text(size = fs_legend_ttl, face = "bold"),
  legend.key.size = unit(0.42, "cm")
)

theme_panel_classic <- theme_classic(base_size = fs_axis_text) +
  base_text_theme +
  theme(plot.margin = margin(6, 12, 6, 8))

theme_panel_minimal <- theme_minimal(base_size = fs_axis_text) +
  base_text_theme +
  theme(panel.grid  = element_blank(),
        plot.margin = margin(6, 8, 6, 6))

# Color palette
col_ua  <- "#2166AC"   # blue for UA
col_hua <- "#B2182B"   # red for HUA
col_sig <- "#B2182B"   # significant
col_ns  <- "grey55"    # non-significant
col_enterotype <- "#D6604D"  # enterotype terms in forest
col_covariate  <- "grey40"   # clinical covariate terms in forest

# ----------------------------------------------------------------------------
# Helpers: panel label rows.
# A title/footnote ROW is a minimal ggplot whose only content is text. Unlike
# plot_annotation(), it survives arbitrary patchwork nesting, so panel labels
# B and C are no longer dropped from the composite.
# ----------------------------------------------------------------------------
panel_title_row <- function(label, text) {
  ggplot() +
    theme_void() +
    labs(title = paste0(label, "  ", text)) +
    theme(
      # "plot" position -> flush to the figure's left edge (not the panel edge),
      # so the A/B/C labels line up with the figure margin as is conventional.
      plot.title.position = "plot",
      plot.title  = element_text(face = "bold", size = fs_panel_title,
                                 hjust = 0, vjust = 1,
                                 margin = margin(t = 1, b = 1)),
      plot.margin = margin(2, 4, 2, 4)
    )
}

panel_footnote_row <- function(text) {
  ggplot() +
    theme_void() +
    labs(title = text) +
    theme(
      plot.title.position = "plot",
      plot.title  = element_text(face = "plain", size = fs_annot,
                                 color = "grey35", hjust = 0,
                                 margin = margin(t = 1, b = 1)),
      plot.margin = margin(2, 4, 2, 4)
    )
}

# Fixed pixel heights (in points) for the label rows; the plotting area itself
# is given the remaining space via a "null" unit.
TITLE_ROW_PT <- fs_panel_title + 9
FOOT_ROW_PT  <- fs_annot + 9

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
           color = "grey50", size = as_mm(fs_annot), hjust = 0, vjust = -0.4) +
  # Lines and points
  geom_line(linewidth = 0.9, alpha = 0.7) +
  geom_point(aes(shape = sig), size = 3.2, stroke = 0.9) +
  # P-value labels above/below points
  geom_text(aes(label = p_label),
            size = as_mm(fs_annot), vjust = -1.1, show.legend = FALSE) +
  # Sample size annotation at bottom
  annotate("text",
           x = 1:6,
           y = rep(-0.15, 6),
           label = paste0("n=", model_n[model_order]),
           size = as_mm(fs_annot), color = "grey40") +
  # Scales
  scale_x_continuous(
    breaks = 1:6,
    labels = model_labels_short[model_order],
    expand = expansion(mult = c(0.06, 0.12))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.10, 0.20)),
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
  labs(x = NULL) +
  theme_panel_classic +
  theme(
    legend.position = c(0.90, 0.80),
    legend.background = element_rect(fill = "white", color = NA),
    legend.spacing.y = unit(0.1, "cm"),
    axis.text.x = element_text(size = fs_axis_text, lineheight = 0.9)
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
  geom_point(aes(color = sig, shape = sig), size = 2.8) +
  scale_color_manual(values = c(`TRUE` = col_sig, `FALSE` = col_ns), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(x = expression(beta ~ "(95% CI)"), y = NULL,
       subtitle = "Continuous UA") +
  theme_panel_classic

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
  geom_point(aes(color = sig, shape = sig), size = 2.8) +
  scale_color_manual(values = c(`TRUE` = col_sig, `FALSE` = col_ns), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  scale_x_continuous(trans = "log2",
                     breaks = c(0.0625, 0.125, 0.25, 0.5, 1, 2, 4),
                     labels = c("0.06", "0.13", "0.25", "0.5", "1", "2", "4")) +
  labs(x = "Odds ratio (95% CI)", y = NULL,
       subtitle = "Binary HUA") +
  theme_panel_classic +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

panel_B_inner <- pB_ua + pB_hua + plot_layout(widths = c(1.3, 1))


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
  eGFR_normal = "eGFR ≥ 90 (n=271)"
)

theme_heatmap <- theme_panel_minimal +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = fs_axis_text),
    axis.text.y = element_text(size = fs_axis_text, color = "black"),
    legend.key.height = unit(0.55, "cm"),
    legend.key.width  = unit(0.32, "cm")
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
      p < 0.10  ~ "†",
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
      p < 0.10  ~ "†",
      TRUE      ~ ""
    )
  )

pC_hua <- ggplot(sub_hua, aes(x = label_f, y = subset_f)) +
  geom_tile(aes(fill = effect_std), color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_level), size = as_mm(fs_star), fontface = "bold", vjust = 0.35) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = "log(OR)"
  ) +
  labs(subtitle = "Binary HUA (OR)", x = NULL, y = NULL) +
  theme_heatmap

pC_ua <- ggplot(sub_ua, aes(x = label_f, y = subset_f)) +
  geom_tile(aes(fill = effect_std), color = "white", linewidth = 0.8) +
  geom_text(aes(label = sig_level), size = as_mm(fs_star), fontface = "bold", vjust = 0.35) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    name = expression(beta)
  ) +
  labs(subtitle = expression("Continuous UA (" * beta * ")"), x = NULL, y = NULL) +
  theme_heatmap +
  theme(axis.text.y = element_blank())

panel_C_inner <- pC_hua + pC_ua + plot_layout(widths = c(1.15, 1))


# ============================================================================
# Assemble each panel with its own (nesting-proof) label row
# ============================================================================

cat("Assembling composite Figure 3...\n")

# NOTE: wrap_plots() is used (instead of the `/` operator) on purpose. The `/`
# operator FLATTENS a sub-patchwork operand into the parent, which would make
# the parent expect more panels than the heights vector provides. wrap_plots()
# treats each argument as exactly one cell, nesting sub-patchworks cleanly.
panel_A_titled <- wrap_plots(
  panel_title_row("A", "Overall enterotype significance across adjustment models"),
  panel_A,
  ncol = 1,
  heights = unit(c(TITLE_ROW_PT, 1), c("pt", "null"))
)

panel_B_titled <- wrap_plots(
  panel_title_row("B", "Cardiometabolic model coefficients (M4, n = 305)"),
  panel_B_inner,
  ncol = 1,
  heights = unit(c(TITLE_ROW_PT, 1), c("pt", "null"))
)

panel_C_titled <- wrap_plots(
  panel_title_row("C", "Subgroup consistency of enterotype contrasts"),
  panel_C_inner,
  panel_footnote_row("*  P < 0.05      **  P < 0.01      †  P < 0.10"),
  ncol = 1,
  heights = unit(c(TITLE_ROW_PT, 1, FOOT_ROW_PT), c("pt", "null", "pt"))
)

# Panel A : Panel B : Panel C  stacked full-width
fig3_composite <- wrap_plots(
  panel_A_titled, panel_B_titled, panel_C_titled,
  ncol = 1,
  heights = c(1, 1.32, 1.06)
)

# Save individual panels (panels/, no _v2_ suffix to match existing files)
ggsave(file.path(fig_dir, "Figure3A_Evidence_Trajectory.png"),
       panel_A_titled, width = 9, height = 4.6, dpi = 300)
ggsave(file.path(fig_dir, "Figure3A_Evidence_Trajectory.pdf"),
       panel_A_titled, width = 9, height = 4.6)

ggsave(file.path(fig_dir, "Figure3B_Cardiomet_Forest.png"),
       panel_B_titled, width = 10, height = 6.2, dpi = 300)
ggsave(file.path(fig_dir, "Figure3B_Cardiomet_Forest.pdf"),
       panel_B_titled, width = 10, height = 6.2)

ggsave(file.path(fig_dir, "Figure3C_Subgroup_Matrix.png"),
       panel_C_titled, width = 9, height = 5.6, dpi = 300)
ggsave(file.path(fig_dir, "Figure3C_Subgroup_Matrix.pdf"),
       panel_C_titled, width = 9, height = 5.6)

# Save composite to results/figures/ (consistent with Figure 1/2/4/5)
ggsave(file.path(out_dir, "Figure3_composite.png"),
      fig3_composite, width = 10, height = 14, dpi = 300)
ggsave(file.path(out_dir, "Figure3_composite.tiff"),
      fig3_composite, width = 10, height = 14, dpi = 300, compression = "lzw")
ggsave(file.path(out_dir, "Figure3_composite.pdf"),
      fig3_composite, width = 10, height = 14)

cat("\nAll Figure 3 panels saved successfully.\n")
cat("Panels directory:", fig_dir, "\n")
cat("Composite directory:", out_dir, "\n")
