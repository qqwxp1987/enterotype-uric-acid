# ==============================================================================
# Figure 6: Composite Assembly (4 panels)
#
# Replaces the legacy combined_figure6.R, whose cached ggplot inputs
# (results/cache/p6a-p6c.rds) are no longer available. This script re-stacks
# the existing panel PNGs with magick, mirroring the original patchwork layout
# (p6a / p6b / (p6c | p6d), heights 1.35/0.8/1.15, C:D width 1.1:1).
#
# Layout:
#   Row 1:      [A]       — Mediation forest (full width)
#   Row 2:      [B]       — Sensitivity matrix (full width)
#   Row 3:  [C]  |  [D]   — Type 4 lollipop | Type 4 boxplots  (C slightly wider)
# ==============================================================================

library(here)
library(magick)

panels_dir <- here("results", "figures", "panels")
out_dir <- here("results", "figures")

# --- 1. Read panels ---
cat("Reading Figure 6 panels...\n")
pA <- image_read(file.path(panels_dir, "Figure6A_mediation_forest_results.png"))
pB <- image_read(file.path(panels_dir, "Figure6B_sensitivity_matrix_results.png"))
pC <- image_read(file.path(panels_dir, "Figure6C_type4_lollipop_results.png"))
pD <- image_read(file.path(panels_dir, "Figure6D_type4_boxplots_results.png"))

for (nm in c("A", "B", "C", "D")) {
  img <- get(paste0("p", nm))
  cat(sprintf("  %s: %d x %d\n", nm, image_info(img)$width, image_info(img)$height))
}

# --- 2. Dimensions ---
full_w <- 4800
# Bottom row: C : D width ratio 1.1 : 1
wC <- round(full_w * 1.1 / 2.1)
wD <- full_w - wC

# --- 3. Helpers ---
add_label <- function(img, label, x = 20, y = 50) {
  image_annotate(img, label,
    size = 70, weight = 700, color = "black",
    location = paste0("+", x, "+", y)
  )
}

make_row <- function(left, right, wl, wr) {
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

# --- 4. Assemble ---
row1 <- make_full_row(add_label(pA, "A"))
row2 <- make_full_row(add_label(pB, "B"))
row3 <- make_row(add_label(pC, "C"), add_label(pD, "D"), wl = wC, wr = wD)

fig6 <- image_append(c(row1, row2, row3), stack = TRUE)

cat(sprintf("Final: %d x %d px\n", image_info(fig6)$width, image_info(fig6)$height))

# --- 5. Save ---
image_write(fig6, file.path(out_dir, "Figure6_composite.png"), format = "png")
image_write(fig6, file.path(out_dir, "Figure6_composite.tiff"),
  format = "tiff", compression = "lzw"
)
cat("Figure 6 composite saved.\n")
