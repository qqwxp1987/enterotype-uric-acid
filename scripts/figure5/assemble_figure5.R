# ==============================================================================
# Figure 5: Composite Assembly (4 panels)
#
# Layout:
#   Row 1:  [A]  |  [B]   — dbRDA triplot   | Top adjusted correlations
#   Row 2:  [C]  |  [D]   — Function heatmap | UA forest plot
# ==============================================================================

library(here)
library(magick)

panels_dir <- here("results", "figures", "panels")
out_dir <- here("results", "figures")

# --- 1. Read panels ---
cat("Reading Figure 5 panels...\n")
pA <- image_read(file.path(panels_dir, "Figure5A_dbRDA_triplot.png"))
pB <- image_read(file.path(panels_dir, "Figure5B_top_adjusted.png"))
pC <- image_read(file.path(panels_dir, "Figure5C_function_heatmap.png"))
pD <- image_read(file.path(panels_dir, "Figure5D_UA_forest_plot.png"))

for (nm in c("A", "B", "C", "D")) {
  img <- get(paste0("p", nm))
  cat(sprintf("  %s: %d x %d\n", nm, image_info(img)$width, image_info(img)$height))
}

# --- 2. Dimensions ---
full_w <- 4800
half_w <- full_w %/% 2

# --- 3. Helpers ---
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

# --- 4. Assemble ---
row1 <- make_row(add_label(pA, "A"), add_label(pB, "B"))
row2 <- make_row(add_label(pC, "C"), add_label(pD, "D"))

fig5 <- image_append(c(row1, row2), stack = TRUE)

cat(sprintf("Final: %d x %d px\n", image_info(fig5)$width, image_info(fig5)$height))

# --- 5. Save ---
image_write(fig5, file.path(out_dir, "Figure5_composite.png"), format = "png")
image_write(fig5, file.path(out_dir, "Figure5_composite.tiff"),
  format = "tiff", compression = "lzw"
)
cat("Figure 5 composite saved.\n")
