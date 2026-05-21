# ==============================================================================
# Figure 4: Composite Assembly (5 panels, restructured)
#
# Layout:
#   Row 1:      [A]       — Core Summary (full width)
#   Row 2:      [B]       — Core Pathways Heatmap, rotated/horizontal (full width)
#   Row 3:  [C]  |  [D]   — Shared EC Direction | Purine Pathway schematic
#   Row 4:      [E]       — Convergent EC Species attribution (full width)
# ==============================================================================

library(here)
library(magick)

panels_dir <- here("results", "figures", "panels")
out_dir <- here("results", "figures")

# --- 1. Read panels ---
cat("Reading Figure 4 panels...\n")
pA <- image_read(file.path(panels_dir, "Figure4A_Core_Summary.png"))
pB <- image_read(file.path(panels_dir, "Figure4B_Core_Pathways_Heatmap.png"))
pC <- image_read(file.path(panels_dir, "Figure4C_Shared_EC_Direction.png"))
pD <- image_read(file.path(panels_dir, "Figure4D_Purine_Pathway.png"))
pE <- image_read(file.path(panels_dir, "Figure4E_Convergent_EC_Species.png"))

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
row1 <- make_full_row(add_label(pA, "A"))
row2 <- make_full_row(add_label(pB, "B"))
row3 <- make_row(add_label(pC, "C"), add_label(pD, "D"))
row4 <- make_full_row(add_label(pE, "E"))

fig4 <- image_append(c(row1, row2, row3, row4), stack = TRUE)

cat(sprintf("Final: %d x %d px\n", image_info(fig4)$width, image_info(fig4)$height))

# --- 5. Save ---
image_write(fig4, file.path(out_dir, "Figure4_composite.png"), format = "png")
image_write(fig4, file.path(out_dir, "Figure4_composite.tiff"),
  format = "tiff", compression = "lzw"
)
cat("Figure 4 composite saved.\n")
