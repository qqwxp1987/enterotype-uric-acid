# ==============================================================================
# Figure 1: Composite Assembly (5 panels)
#
# Layout:
#   Row 1:  [A]  |  [B]   — DMM Selection | PCoA
#   Row 2:  [C]  |  [D]   — Stacked Bar   | P/B Ratio
#   Row 3:      [E]       — Alpha Diversity (Shannon + Observed, full width)
# ==============================================================================

library(here)
library(magick)

panels_dir <- here("results", "figures", "panels")
out_dir <- here("results", "figures")

# --- 1. Read panels ---
cat("Reading Figure 1 panels...\n")
pA <- image_read(file.path(panels_dir, "Figure1A_DMM_Selection.png"))
pB <- image_read(file.path(panels_dir, "Figure1B_PCoA.png"))
pC <- image_read(file.path(panels_dir, "Figure1C_Stacked_Bar.png"))
pD <- image_read(file.path(panels_dir, "Figure1D_PB_Ratio.png"))
pE <- image_read(file.path(panels_dir, "Figure1E_Alpha_Diversity.png"))

for (nm in c("A", "B", "C", "D", "E")) {
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

make_full_row <- function(img, w = full_w) {
  img <- image_resize(img, paste0(w, "x"))
  image_extent(img, paste0(w, "x", image_info(img)$height),
    gravity = "NorthWest", color = "white"
  )
}

# --- 4. Assemble ---
row1 <- make_row(add_label(pA, "A"), add_label(pB, "B"))
row2 <- make_row(add_label(pC, "C"), add_label(pD, "D"))
row3 <- make_full_row(add_label(pE, "E"))

fig1 <- image_append(c(row1, row2, row3), stack = TRUE)

cat(sprintf("Final: %d x %d px\n", image_info(fig1)$width, image_info(fig1)$height))

# --- 5. Save ---
image_write(fig1, file.path(out_dir, "Figure1_composite.png"), format = "png")
image_write(fig1, file.path(out_dir, "Figure1_composite.tiff"),
  format = "tiff", compression = "lzw"
)
cat("Figure 1 composite saved.\n")
