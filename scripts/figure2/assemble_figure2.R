# ==============================================================================
# Figure 2: Composite Assembly
#
# Layout:
#   Row 1:  [A]                    — Clinical Screening (full width)
#   Row 2:  [B]   |  [C]  |  [D]  — UA Violin | HUA Prevalence | Forest Plot
# ==============================================================================

library(here)
library(magick)

panels_dir <- here("results", "figures", "panels")
out_dir <- here("results", "figures")

# --- 1. Read panels ---
cat("Reading Figure 2 panels...\n")
pA <- image_read(file.path(panels_dir, "Figure2A_Clinical_Screening.png"))
pB <- image_read(file.path(panels_dir, "Figure2B_UA_Violin.png"))
pC <- image_read(file.path(panels_dir, "Figure2C_HUA_Prevalence.png"))
pD <- image_read(file.path(panels_dir, "Figure2D_Forest_Plot.png"))

for (nm in c("A", "B", "C", "D")) {
    img <- get(paste0("p", nm))
    cat(sprintf("  %s: %d x %d\n", nm, image_info(img)$width, image_info(img)$height))
}

# --- 2. Dimensions ---
full_w <- 4800

# --- 3. Helper ---
add_label <- function(img, label, x = 20, y = 50) {
    image_annotate(img, label,
        size = 70, weight = 700, color = "black",
        location = paste0("+", x, "+", y)
    )
}

# --- 4. Row 1: A full width ---
pA_r <- image_resize(add_label(pA, "A"), paste0(full_w, "x"))

# --- 5. Row 2: B, C, D each ~1/3 width ---
w_third <- full_w %/% 3
w_last <- full_w - 2 * w_third # handle rounding

pB_r <- image_resize(add_label(pB, "B"), paste0(w_third, "x"))
pC_r <- image_resize(add_label(pC, "C"), paste0(w_third, "x"))
pD_r <- image_resize(add_label(pD, "D"), paste0(w_last, "x"))

h_row2 <- max(image_info(pB_r)$height, image_info(pC_r)$height, image_info(pD_r)$height)
pB_r <- image_extent(pB_r, paste0(w_third, "x", h_row2), gravity = "NorthWest", color = "white")
pC_r <- image_extent(pC_r, paste0(w_third, "x", h_row2), gravity = "NorthWest", color = "white")
pD_r <- image_extent(pD_r, paste0(w_last, "x", h_row2), gravity = "NorthWest", color = "white")

row2 <- image_append(c(pB_r, pC_r, pD_r))

# --- 6. Stack ---
fig2 <- image_append(c(pA_r, row2), stack = TRUE)

cat(sprintf("Final: %d x %d px\n", image_info(fig2)$width, image_info(fig2)$height))

# --- 7. Save ---
image_write(fig2, file.path(out_dir, "Figure2_composite.png"), format = "png")
image_write(fig2, file.path(out_dir, "Figure2_composite.tiff"),
    format = "tiff", compression = "lzw"
)
cat("Figure 2 composite saved.\n")
