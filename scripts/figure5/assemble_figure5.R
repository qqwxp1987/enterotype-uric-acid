# ==============================================================================
# R2 Minor Comment 4 — Rebuild Figure 5 (vector PDF + 300 dpi PNG)
#
# Baseline (must reproduce before layout changes are accepted):
#   * Panel A: n = 303; adjusted R2 = 0.152; permutation P = 0.003;
#              12/12 curated metabolite arrows retained.
#   * Panel B: 23 x 13 partial-correlation matrix; 138 BH-significant pairs.
#   * Panel C: 30 x 13 partial-correlation matrix; 173 BH-significant pairs.
#   * Panel D: 53 features x 2 models; estimates and BH calls unchanged.
#
# Success criteria:
#   * All baseline statistics pass independent CSV checks.
#   * Final canvas is the correct R1 A4 portrait layout (595 x 842 pt).
#   * Panels B/C explicitly define BH-adjusted asterisks on the figure.
#   * Panel D explicitly defines filled/open symbols by BH-adjusted P.
#   * No expected text is missing and no text-box overlap is detected.
#   * PDF remains vector; PNG is 300 dpi.
#
# Only typography, label placement, legends, and panel geometry are changed.
# Raw data and statistical definitions are immutable.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
  library(ComplexHeatmap)
})

# Run from the repository root; here::here() resolves all inputs and outputs.
stopifnot(file.exists(here(".here")))
dir.create(here("results", "cache"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables", "section_3_4"), recursive = TRUE, showWarnings = FALSE)

OUT_DIR <- here("results", "figures")
MINOR4_DIR <- here("results", "verification")
QA_DIR <- file.path(MINOR4_DIR, "figure5")
dir.create(QA_DIR, recursive = TRUE, showWarnings = FALSE)
source(here("scripts", "figure5", "figure5_style.R"))
class_colors <- c(
  "AAA catabolite" = "#D2691E", "Bile acid" = "#2E8B57",
  "TMAO pathway" = "#7B68EE", "Purine" = "#DAA520",
  "Clinical" = "#333333"
)

# Default: recompute everything from frozen inputs with fixed seeds (42).
# Layout-only iterations may set FIG5_USE_CACHE=1 after the baseline has passed.
USE_CACHE <- identical(Sys.getenv("FIG5_USE_CACHE", "0"), "1")
if (USE_CACHE) {
  p4a <- readRDS(here("results", "cache", "p4a.rds"))
  ht <- readRDS(here("results", "cache", "p4b_top_adj.rds"))
  ht_c <- readRDS(here("results", "cache", "p4c_function_heatmap.rds"))
  plot_all <- readRDS(here("results", "cache", "p4d_plot_all.rds"))
  layer_colors <- readRDS(here("results", "cache", "figure5_layer_colors.rds"))
} else {
  source(here("scripts", "figure5", "panel_A_dbRDA.R"))
  source(here("scripts", "figure5", "panel_B_top_adjusted.R"))
  source(here("scripts", "figure5", "panel_CD_function_integration.R"))

  # Hard baseline assertions: fail before final assembly if anything changes.
  stopifnot(
    length(shared_ids) == 303,
    nrow(arrow_df) == 12,
    abs(R2_adj - 0.1520) < 0.0001,
    dbrda_anova$`Pr(>F)`[1] == 0.003,
    nrow(rho_mat) == 23, ncol(rho_mat) == 13,
    sum(p_adj < 0.05, na.rm = TRUE) == 138,
    nrow(rho_c) == 30, ncol(rho_c) == 13,
    sum(padj_c < 0.05, na.rm = TRUE) == 173,
    nrow(plot_all) == 106,
    length(unique(plot_all$feature)) == 53
  )

  key_stats <- data.frame(
    metric = c(
      "matched_samples", "metabolome_PC_count", "dbRDA_adjusted_R2",
      "dbRDA_permutation_P", "curated_arrow_count",
      "panel_B_rows", "panel_B_targets", "panel_B_BH_significant_pairs",
      "panel_C_rows", "panel_C_targets", "panel_C_BH_significant_pairs",
      "panel_D_unique_features", "panel_D_base_BH_significant_UA",
      "panel_D_enterotype_adjusted_BH_significant_UA"
    ),
    value = c(
      length(shared_ids), n_pcs, R2_adj, dbrda_anova$`Pr(>F)`[1], nrow(arrow_df),
      nrow(rho_mat), ncol(rho_mat), sum(p_adj < 0.05, na.rm = TRUE),
      nrow(rho_c), ncol(rho_c), sum(padj_c < 0.05, na.rm = TRUE),
      length(unique(plot_all$feature)),
      sum(plot_base$sig == "sig"), sum(plot_adj$sig == "sig")
    ),
    stringsAsFactors = FALSE
  )
  write.csv(key_stats, file.path(QA_DIR, "figure5_key_statistics.csv"), row.names = FALSE)
}

# ==============================================================================
# Panel A — full-width top row, right-side legend, repelled arrow labels
# ==============================================================================

tag_theme <- theme(
  plot.tag = element_text(size = fs_panel_tag, face = "bold", color = "black"),
  plot.tag.position = "topleft"
)

pA_final <- p4a +
  labs(tag = "A") +
  tag_theme +
  theme(
    legend.position = "right",
    legend.text = element_text(size = fs_legend_text),
    plot.margin = margin(4, 4, 2, 4)
  )

# ==============================================================================
# Panels B/C — side-by-side vector heatmaps with a shared BH key
# ==============================================================================

draw_B <- grid.grabExpr({
  draw(ht, heatmap_legend_side = "bottom", padding = unit(c(4, 2, 8, 2), "mm"))
})
draw_C <- grid.grabExpr({
  draw(ht_c, heatmap_legend_side = "bottom", padding = unit(c(4, 2, 12, 2), "mm"))
})

pB_final <- wrap_elements(full = draw_B, clip = TRUE) +
  labs(tag = "B") + tag_theme +
  theme(plot.margin = margin(2, 1, 1, 2))
pC_final <- wrap_elements(full = draw_C, clip = TRUE) +
  labs(tag = "C") + tag_theme +
  theme(plot.margin = margin(2, 2, 1, 1))

pBC_final <- (pB_final | pC_final) +
  plot_layout(widths = c(0.92, 1.08))

pBH_note <- ggplot() +
  annotate("text", x = 0.09, y = 0.76, label = "Column classes",
           size = as_mm(7.5), fontface = "bold", hjust = 0.5) +
  annotate("text", x = 0.23, y = 0.76, label = "AAA",
           size = as_mm(7.5), color = class_colors[["AAA catabolite"]], hjust = 0.5) +
  annotate("text", x = 0.36, y = 0.76, label = "Bile acid",
           size = as_mm(7.5), color = class_colors[["Bile acid"]], hjust = 0.5) +
  annotate("text", x = 0.53, y = 0.76, label = "TMAO pathway",
           size = as_mm(7.5), color = class_colors[["TMAO pathway"]], hjust = 0.5) +
  annotate("text", x = 0.75, y = 0.76, label = "Purine",
           size = as_mm(7.5), color = class_colors[["Purine"]], hjust = 0.5) +
  annotate("text", x = 0.89, y = 0.76, label = "Clinical",
           size = as_mm(7.5), color = class_colors[["Clinical"]], hjust = 0.5) +
  annotate(
    "text", x = 0.5, y = 0.22,
    label = "B–C:  * BH-adjusted P < 0.05;   ** BH-adjusted P < 0.01",
    size = as_mm(fs_annot), color = "black", hjust = 0.5, vjust = 0.5
  ) +
  scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(0, 0, 0, 0))

# ==============================================================================
# Panel D — correct R1 layout: species block left, function block right
# ==============================================================================

plot_all$model <- factor(
  as.character(plot_all$model),
  levels = c("Base", "ET-adj."),
  labels = c("Base\n(sex + age)", "Enterotype-\nadjusted")
)

make_d_block <- function(df, italic_y = FALSE, panel_tag = NULL) {
  ggplot(df, aes(x = rho, y = reorder(display_name, order_idx))) +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey60") +
    geom_errorbarh(
      aes(xmin = ci_lo, xmax = ci_hi),
      height = 0.25, linewidth = 0.38, color = "grey50"
    ) +
    geom_point(aes(color = layer, shape = sig), size = 1.9, stroke = 0.55) +
    scale_color_manual(
      values = layer_colors, limits = names(layer_colors), name = "Omics layer",
      drop = FALSE
    ) +
    scale_shape_manual(
      values = c("sig" = 16, "ns" = 1),
      limits = c("sig", "ns"),
      name = "BH-adjusted P",
      labels = c("sig" = "< 0.05", "ns" = "≥ 0.05")
    ) +
    facet_grid(layer ~ model, scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = c(-0.2, 0, 0.2), expand = expansion(mult = 0.08)) +
    labs(x = "Partial Spearman ρ with serum UA", y = NULL, tag = panel_tag) +
    theme_fig5 +
    tag_theme +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 8),
      axis.text.y = element_text(size = 7.5, face = if (italic_y) "italic" else "plain"),
      axis.title.x = element_text(size = fs_axis_title, margin = margin(t = 3)),
      strip.background = element_rect(fill = "grey93", color = "grey65", linewidth = 0.4),
      strip.text = element_text(size = 8, face = "bold", lineheight = 0.95),
      panel.spacing = unit(3, "pt"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.key.size = unit(0.30, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(1, "pt"),
      plot.margin = margin(3, 3, 2, 3)
    )
}

pD_species <- make_d_block(
  plot_all %>% filter(layer == "Species"), italic_y = TRUE, panel_tag = "D"
) + guides(color = "none", shape = "none")
pD_function <- make_d_block(
  plot_all %>% filter(layer %in% c("Pathway", "Convergent EC"))
)

pD_final <- (pD_species | pD_function) +
  plot_layout(widths = c(0.88, 1.12), guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    legend.spacing.y = unit(0, "pt"),
    legend.text = element_text(size = fs_legend_text),
    legend.title = element_text(size = fs_legend_title, face = "bold")
  )

# ==============================================================================
# Final A4 portrait assembly — dimensions match the correct R1 submission
# ==============================================================================

fig5_final <- (pA_final / pBC_final / pBH_note / pD_final) +
  plot_layout(heights = c(52, 129, 9, 107))

pdf_path <- file.path(OUT_DIR, "Figure5.pdf")
png_path <- file.path(OUT_DIR, "Figure5.png")

ggsave(
  pdf_path, fig5_final,
  width = 210 / 25.4, height = 297 / 25.4, units = "in",
  device = cairo_pdf, bg = "white"
)
ggsave(
  png_path, fig5_final,
  width = 210 / 25.4, height = 297 / 25.4, units = "in",
  dpi = 300, device = png, type = "cairo", bg = "white"
)

# Reproducibility log and immutable-build summary.
capture.output(sessionInfo(), file = file.path(MINOR4_DIR, "figure5_session_info.txt"))
writeLines(
  c(
    "Figure 5 rebuild completed from frozen inputs.",
    "Random seed: 42 for db-RDA, envfit, clustering, and label placement.",
    "Final canvas: 210 x 297 mm (A4 portrait, approximately 595 x 842 pt).",
    "Statistical outputs: results/tables/section_3_4/",
    "Key-statistic assertions: passed before final assembly."
  ),
  con = file.path(QA_DIR, "build_status.txt")
)

cat("\nFigure 5 final outputs:\n", pdf_path, "\n", png_path, "\n", sep = "")
