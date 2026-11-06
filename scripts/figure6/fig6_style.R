# ==============================================================================
# Figure 6: shared typography / style — single source of truth.
#
# Source this from every Figure 6 panel script:
#   source(here("scripts", "figure6", "fig6_style.R"))
#
# Why this file exists
# --------------------
# Figure 6 is assembled by stitching the four panel PNGs (assemble_figure6.R).
# For text at the SAME role (axis tick / axis title / facet strip / panel title /
# legend) to render at the SAME size in the final composite, two things must hold:
#   1. All panels use one shared font-size hierarchy  -> this file.
#   2. Each panel is rendered at a width equal to its composite cell width, so
#      magick performs no rescale (scale factor = 1)  -> panel ggsave() sizes.
# The hierarchy numbers mirror Figure 3 (section3_6_step4_visualization_v2.R) so
# the two figures stay consistent across the paper.
#
# "Consistent" means SAME role -> SAME size, NOT one uniform size for everything.
# ==============================================================================

library(ggplot2)

# --- Font-size hierarchy (points) ---------------------------------------------
# panel letter A/B/C/D (drawn by magick at 70 px ~ 16.8 pt) > title > strip /
# axis title > axis tick / legend title > legend text / in-plot annotation.
FS_TITLE        <- 11    # panel descriptive title (plot.title)
FS_STRIP        <- 9.5   # facet strip labels
FS_AXIS_TITLE   <- 9.5   # axis titles
FS_AXIS_TEXT    <- 8.5   # axis tick labels
FS_LEGEND_TITLE <- 8.5   # legend titles
FS_LEGEND_TEXT  <- 8     # legend entries
FS_ANNOT        <- 8.5   # in-plot text annotations (e.g. panel D stat labels)

# --- pt -> mm helper ----------------------------------------------------------
# ggplot geom_text()/annotate() take size in mm; convert from points.
PT_PER_MM <- 2.845276    # == ggplot2's .pt
as_mm <- function(pt) pt / PT_PER_MM

# --- Shared text theme --------------------------------------------------------
# Layered ON TOP of each panel's base theme (theme_bw / theme_minimal /
# theme_classic) so all three theme families share identical text sizes.
# Apply order in panels: base_theme(base_size = FS_AXIS_TEXT) +
#                        fig6_text_theme() + theme(<panel-specific layout>)
fig6_text_theme <- function() {
  theme(
    text         = element_text(color = "black"),
    plot.title   = element_text(size = FS_TITLE, face = "bold"),
    strip.text   = element_text(size = FS_STRIP, face = "bold", color = "black"),
    axis.title   = element_text(size = FS_AXIS_TITLE, color = "black"),
    axis.text    = element_text(size = FS_AXIS_TEXT, color = "black"),
    legend.title = element_text(size = FS_LEGEND_TITLE, face = "bold"),
    legend.text  = element_text(size = FS_LEGEND_TEXT)
  )
}
