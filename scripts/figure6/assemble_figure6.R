# ==============================================================================
# R2 Minor Comment 4 — Rebuild Figure 6 (vector PDF + 300 dpi PNG)
#
# Baseline and success criteria:
#   * 22 candidates x 2 models = 44 unchanged ACME estimates and intervals.
#   * Five BH-significant candidates in each main model (Type 1 and Type 4).
#   * Six priority features x 8 sensitivity settings = 48 unchanged tiles.
#   * 23 within-Type 4 tests, with exactly four nominal P < 0.05 results.
#   * 39 Type 4 participants: 12 HUA+ and 27 HUA−.
#   * A/B explicitly report BH-adjusted significance; C/D explicitly report
#     nominal, unadjusted significance and exploratory framing.
#   * A4 portrait vector PDF, no missing expected text or genuine overlap.
#
# No mediation model, bootstrap result, sensitivity result, or Type 4 test is
# recomputed. Panels read the frozen R1 result tables. Only presentation changes.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

stopifnot(file.exists(here(".here")))
dir.create(here("results", "cache"), recursive = TRUE, showWarnings = FALSE)

OUT_DIR <- here("results", "figures")
MINOR4_DIR <- here("results", "verification")
QA_DIR <- file.path(MINOR4_DIR, "figure6")
dir.create(QA_DIR, recursive = TRUE, showWarnings = FALSE)

source(here("scripts", "figure6", "panel_A_mediation_forest.R"))
source(here("scripts", "figure6", "panel_B_sensitivity_matrix.R"))
source(here("scripts", "figure6", "panel_C_type4_lollipop.R"))
source(here("scripts", "figure6", "panel_D_type4_boxplots.R"))

# Hard fidelity assertions against the frozen R1 result tables.
stopifnot(
  nrow(med_results) == 44,
  length(unique(med_results$mediator)) == 22,
  sum(med_results$model == "Type1_vs_Ref" & med_results$ACME_padj < 0.05) == 5,
  sum(med_results$model == "Type4_vs_Ref" & med_results$ACME_padj < 0.05) == 5,
  nrow(panel_b) == 48,
  length(unique(panel_b$analysis)) == 8,
  length(unique(panel_b$mediator)) == 6,
  nrow(panel_c) == 23,
  sum(panel_c$is_sig == "< 0.05") == 4,
  nrow(t4) == 39,
  sum(t4$HUA_group == "HUA+") == 12,
  sum(t4$HUA_group == "HUA-") == 27,
  nrow(panel_d_ann) == 4
)

key_stats <- data.frame(
  metric = c(
    "mediation_candidates", "mediation_model_rows",
    "Type1_BH_significant", "Type4_BH_significant",
    "sensitivity_features", "sensitivity_settings", "sensitivity_rows",
    "within_Type4_tests", "within_Type4_nominal_significant",
    "Type4_HUA_positive", "Type4_HUA_negative"
  ),
  value = c(
    length(unique(med_results$mediator)), nrow(med_results),
    sum(med_results$model == "Type1_vs_Ref" & med_results$ACME_padj < 0.05),
    sum(med_results$model == "Type4_vs_Ref" & med_results$ACME_padj < 0.05),
    length(unique(panel_b$mediator)), length(unique(panel_b$analysis)), nrow(panel_b),
    nrow(panel_c), sum(panel_c$is_sig == "< 0.05"),
    sum(t4$HUA_group == "HUA+"), sum(t4$HUA_group == "HUA-")
  ),
  stringsAsFactors = FALSE
)
write.csv(key_stats, file.path(QA_DIR, "figure6_key_statistics.csv"), row.names = FALSE)

# Panel tags are added only at final assembly; all text remains vector.
tag_theme <- theme(
  plot.tag = element_text(size = 14, face = "bold", color = "black"),
  plot.tag.position = "topleft"
)

p6a_final <- p6a +
  labs(tag = "A", title = NULL) + tag_theme +
  theme(
    legend.text = element_text(size = 7.5),
    legend.title = element_text(size = 8, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(4, 4, 2, 4)
  )

p6b_final <- p6b +
  labs(tag = "B", title = NULL) + tag_theme +
  theme(plot.margin = margin(3, 3, 2, 3))

p6c_final <- p6c +
  labs(tag = "C", title = NULL, caption = NULL) + tag_theme +
  theme(plot.margin = margin(3, 3, 2, 3))

p6d_final <- p6d +
  labs(tag = "D", title = NULL) + tag_theme +
  theme(plot.margin = margin(3, 3, 2, 3))

# A and B each span the full width; C and D share the bottom row. This preserves
# the A -> B -> C -> D reading order while giving all 23 rows in C enough height.
p6_bottom <- (p6c_final | p6d_final) +
  plot_layout(widths = c(1.0, 1.0))

fig6_final <- (p6a_final / p6b_final / p6_bottom) +
  plot_layout(heights = c(1.12, 0.60, 1.25))

pdf_path <- file.path(OUT_DIR, "Figure6.pdf")
png_path <- file.path(OUT_DIR, "Figure6.png")

ggsave(
  pdf_path, fig6_final,
  width = 210 / 25.4, height = 297 / 25.4, units = "in",
  device = cairo_pdf, bg = "white"
)
ggsave(
  png_path, fig6_final,
  width = 210 / 25.4, height = 297 / 25.4, units = "in",
  dpi = 300, device = png, type = "cairo", bg = "white"
)

capture.output(sessionInfo(), file = file.path(MINOR4_DIR, "figure6_session_info.txt"))
writeLines(
  c(
    "Figure 6 rebuild completed from frozen R1 result tables.",
    "Random seed: 42 for Type 4 participant jitter (unchanged).",
    "Final canvas: 210 x 297 mm (A4 portrait, approximately 595 x 842 pt).",
    "All hard fidelity assertions passed before final assembly."
  ),
  con = file.path(QA_DIR, "build_status.txt")
)

cat("\nFigure 6 final outputs:\n", pdf_path, "\n", png_path, "\n", sep = "")
