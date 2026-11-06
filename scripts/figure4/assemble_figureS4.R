# ==============================================================================
# Supplementary Figure S2: Composite Assembly (4 panels)
#
# Layout:
#   Row 1:      [A]       — Radial core-taxa tree (full width)
#   Row 2:  [B]  |  [C]   — Species attribution | Type 1 core metabolites
#   Row 3:      [D]       — Comprehensive core-pathway heatmap (full width)
#
# Reads panels from / writes composite to: supplementary/figures/
# ==============================================================================

library(here)
library(magick)

fig_dir <- here("supplementary", "figures")

cat("Reading Supplementary Figure S2 panels...\n")
pA <- image_read(file.path(fig_dir, "FigureS4A_Radial_Taxa.png"))
pB <- image_read(file.path(fig_dir, "FigureS4B_Species_Attribution.png"))
pC <- image_read(file.path(fig_dir, "FigureS4C_Type1_Core_Metabolites.png"))
pD <- image_read(file.path(fig_dir, "FigureS4D_Core_Pathways_Heatmap.png"))

full_w <- 4800
half_w <- full_w %/% 2

add_label <- function(img, label, x = 20, y = 50) {
  image_annotate(img, label,
    size = 70, weight = 700, color = "black",
    location = paste0("+", x, "+", y)
  )
}

make_row <- function(left, right, wl = half_w, wr = half_w) {
  left <- image_resize(left, paste0(wl, "x"))
  right <- image_resize(right, paste0(wr, "x"))
  h <- max(image_info(left)$height, image_info(right)$height)
  left <- image_extent(left, paste0(wl, "x", h), gravity = "NorthWest", color = "white")
  right <- image_extent(right, paste0(wr, "x", h), gravity = "NorthWest", color = "white")
  image_append(c(left, right))
}

make_full_row <- function(img, w = full_w) {
  img <- image_resize(img, paste0(w, "x"))
  image_extent(img, paste0(w, "x", image_info(img)$height),
    gravity = "NorthWest", color = "white"
  )
}

row1 <- make_full_row(add_label(pA, "A"))
row2 <- make_row(add_label(pB, "B"), add_label(pC, "C"))
row3 <- make_full_row(add_label(pD, "D"))

figS2 <- image_append(c(row1, row2, row3), stack = TRUE)

cat(sprintf("Final: %d x %d px\n", image_info(figS2)$width, image_info(figS2)$height))

image_write(figS2, file.path(fig_dir, "FigureS4_composite.png"), format = "png")
image_write(figS2, file.path(fig_dir, "FigureS4_composite.tiff"),
  format = "tiff", compression = "lzw"
)
cat("Supplementary Figure S2 composite saved.\n")
