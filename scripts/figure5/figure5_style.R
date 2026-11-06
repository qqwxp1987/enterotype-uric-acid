# ==============================================================================
# Figure 5 — shared style (typography + per-panel sizes)
#
# Sourced by every Figure 5 panel script. Each panel is built as a SELF-
# CONTAINED, standalone figure: its plotting area dominates, every label is
# clearly legible, and annotation/legend area is kept modest. Composition into
# the final Figure 5 is a separate, later step (not done here).
# ==============================================================================

# --- 1. Typography hierarchy (points) -----------------------------------------
# One hierarchy reused by all panels: the SAME role gets the SAME size in every
# panel (axis ticks identical in A/B/C/D, etc.), NOT one global size.
fs_panel_tag    <- 14     # A / B / C / D corner labels (used only at merge time)
fs_strip        <- 9.5    # forest facet strips; heatmap split-group titles
fs_axis_title   <- 9.5    # axis titles
fs_axis_text    <- 8.5    # axis ticks; heatmap row & column names
fs_legend_title <- 9      # legend / colour-bar titles
fs_legend_text  <- 8.5    # legend entries; annotation-track names
fs_annot        <- 8.5    # in-plot notes (P values, R-squared, footnotes)
fs_star         <- 7      # compact enough to avoid cross-row star collisions

# ggplot geom_text()/annotate() take `size` in mm; convert from points.
PT_PER_MM <- 2.845276          # == ggplot2's .pt
as_mm <- function(pt) pt / PT_PER_MM

# Abbreviate a binomial species name to genus-initial form so heatmap / forest
# row labels stay compact, e.g.
#   "Faecalibacterium prausnitzii" -> "F. prausnitzii"
#   "Roseburia sp CAG 471"         -> "R. sp CAG 471"
# Non-binomials (single token) pass through unchanged.
abbrev_species <- function(x) {
  vapply(strsplit(trimws(x), "\\s+"), function(p) {
    if (length(p) >= 2 && nchar(p[1]) > 1) {
      paste0(substr(p[1], 1, 1), ". ", paste(p[-1], collapse = " "))
    } else {
      paste(p, collapse = " ")
    }
  }, character(1))
}

# --- 2. Per-panel standalone sizes (mm) ---------------------------------------
# Generous, comfortable sizes so each panel reads well on its own.
MM_PER_IN <- 25.4
in_of <- function(mm) mm / MM_PER_IN

fig5_A_w_mm <- 150 ; fig5_A_h_mm <- 86    # A: db-RDA triplot
fig5_B_w_mm <- 112 ; fig5_B_h_mm <- 160   # B: species x metabolite heatmap
fig5_C_w_mm <- 130 ; fig5_C_h_mm <- 194   # C: function x metabolite heatmap
fig5_D_w_mm <- 128 ; fig5_D_h_mm <- 210   # D: UA forest plot

fig5_A_w_in <- in_of(fig5_A_w_mm) ; fig5_A_h_in <- in_of(fig5_A_h_mm)
fig5_B_w_in <- in_of(fig5_B_w_mm) ; fig5_B_h_in <- in_of(fig5_B_h_mm)
fig5_C_w_in <- in_of(fig5_C_w_mm) ; fig5_C_h_in <- in_of(fig5_C_h_mm)
fig5_D_w_in <- in_of(fig5_D_w_mm) ; fig5_D_h_in <- in_of(fig5_D_h_mm)

fig5_DPI <- 600           # raster render DPI (submission-grade TIFF/PNG)

# ComplexHeatmap body sizes (mm). The body must DOMINATE the panel; the discrete
# annotation legends are dropped (row/column split titles already name every
# group), so only the small colour-bar remains. Tune by eye.
fig5_B_body <- c(w = 48, h = 70)    # preserve >=8 pt pitch for 23 species rows
fig5_C_body <- c(w = 45, h = 82)    # long pathway names must remain untruncated

# --- 3. Shared ggplot theme (panels A & D) ------------------------------------
# Defined lazily: panels B & C source this file but do not load ggplot2.
if (requireNamespace("ggplot2", quietly = TRUE)) {
  theme_fig5 <- ggplot2::theme_bw(base_size = fs_axis_text) +
    ggplot2::theme(
      text            = ggplot2::element_text(color = "black"),
      axis.title      = ggplot2::element_text(size = fs_axis_title, color = "black"),
      axis.text       = ggplot2::element_text(size = fs_axis_text, color = "black"),
      legend.text     = ggplot2::element_text(size = fs_legend_text),
      legend.title    = ggplot2::element_text(size = fs_legend_title, face = "bold"),
      legend.key.size = ggplot2::unit(0.40, "cm"),
      strip.text      = ggplot2::element_text(size = fs_strip, face = "bold"),
      plot.margin     = ggplot2::margin(5, 5, 5, 5)
    )
}
