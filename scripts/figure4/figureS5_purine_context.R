# ==============================================================================
# Additional file 5: Figure S5 — purine pathway-context schematic
#
# This is the original Figure 4D, exported independently after moving it out of
# the main figure to reduce density. Content and geometry are unchanged; text is
# increased to 8 pt (key nodes) and 7 pt (secondary annotations) at 180 mm width.
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(here)
})

out_fig <- here("supplementary", "figures")
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)

.PT <- 72.27 / 25.4
GEOM_PRIMARY   <- 8 / .PT
GEOM_SECONDARY <- 7 / .PT
FIG4_W_FULL <- 7.0866
FIG4_H_D <- 3.40
FIG4_DPI <- 300

# --- 1. Nodes -----------------------------------------------------------------
# tier: "key" -> primary (8 pt), "detail" -> secondary (7 pt)
nodes <- tribble(
  ~label, ~x, ~y, ~w, ~h, ~fill, ~tier, ~fface, ~tcol,

  # Row 1: upstream pathway context (y = 14.2)
  "Folate biosynthesis\n(KEGG map00790)",
      2, 14.2, 3.0, 0.85, "#A5D6A7", "detail", "bold", "grey20",
  "Ala / Asp / Glu metabolism\n(KEGG map00250)",
      6, 14.2, 3.0, 0.85, "#FFCC80", "detail", "bold", "grey20",
  "One-carbon / Met metabolism\n(KEGG map00670, map00270)",
     10, 14.2, 3.0, 0.85, "#90CAF9", "detail", "bold", "grey20",

  # Row 2: convergent ECs (y = 12.7) — highlighted
  "EC 2.7.6.3\nHPPK",
      2, 12.7, 3.0, 0.85, "#2E7D32", "key", "bold", "white",
  "EC 2.6.1.1\nAspartate transaminase",
      6, 12.7, 3.0, 0.85, "#E65100", "key", "bold", "white",
  "EC 2.1.1.13\nMethionine synthase",
     10, 12.7, 3.0, 0.85, "#1565C0", "key", "bold", "white",

  # Row 3: reaction / products (y = 11.3)
  "6-HMDP + ATP\n→ 6-HMDP-PP",
      2, 11.3, 3.0, 0.75, "#C8E6C9", "detail", "plain", "grey25",
  "L-Asp + 2-OG\n↔ OAA + L-Glu",
      6, 11.3, 3.0, 0.75, "#FFE0B2", "detail", "plain", "grey25",
  "5-mTHF + L-Hcy\n→ THF + L-Met",
     10, 11.3, 3.0, 0.75, "#BBDEFB", "detail", "plain", "grey25",

  # Row 4: downstream substrates entering IMP biosynthesis (y = 10)
  "→ THF → 10-formyl-THF\n(purine ring C₂, C₈)",
      2, 10.0, 3.0, 0.75, "#E8F5E9", "detail", "italic", "grey35",
  "L-Aspartate\n(purine ring N₁)",
      6, 10.0, 2.4, 0.75, "#FFF3E0", "detail", "italic", "grey35",
  "→ L-Met → SAM\n(one-carbon / adenosine)",
     10, 10.0, 3.0, 0.75, "#E3F2FD", "detail", "italic", "grey35",

  # Row 5: de novo IMP biosynthesis (y = 8.3)
  "De novo purine biosynthesis (KEGG M00048)\nPRPP + Gln + Gly + CO₂ → (12 steps) → IMP",
      6,  8.3, 7.0, 1.0, "#CE93D8", "key", "bold", "grey10",

  # Row 6: IMP (y = 6.8)
  "IMP",
      6, 6.8, 1.6, 0.6, "#E1BEE7", "key", "bold", "grey10",

  # Row 7: nucleotides (y = 5.5)
  "AMP",
    3.2, 5.5, 1.5, 0.55, "#F8BBD0", "key", "bold", "grey20",
  "XMP → GMP",
    8.8, 5.5, 2.0, 0.55, "#F8BBD0", "key", "bold", "grey20",

  # Row 8: degradation intermediates (y = 4.1)
  "Adenosine → Inosine\n→ Hypoxanthine",
    3.2, 4.1, 2.6, 0.8, "#F5F5F5", "detail", "plain", "grey40",
  "Guanosine\n→ Guanine",
    8.8, 4.1, 2.2, 0.8, "#F5F5F5", "detail", "plain", "grey40",

  # Row 9: Xanthine (y = 2.8)
  "Xanthine",
      6, 2.8, 1.8, 0.55, "#FFF176", "key", "bold", "grey20",

  # Row 10: Uric Acid (y = 1.5)
  "Uric Acid",
      6, 1.5, 2.2, 0.7, "#E53935", "key", "bold", "white"
)
nodes$fsize <- ifelse(nodes$tier == "key", GEOM_PRIMARY, GEOM_SECONDARY)

# --- 2. Edges -----------------------------------------------------------------
edges <- tribble(
  ~x,    ~y,    ~xend,  ~yend,  ~lty,
  # Context -> EC
   2,   13.78,   2,    13.13,  "solid",
   6,   13.78,   6,    13.13,  "solid",
  10,   13.78,  10,    13.13,  "solid",
  # EC -> reaction
   2,   12.28,   2,    11.68,  "solid",
   6,   12.28,   6,    11.68,  "solid",
  10,   12.28,  10,    11.68,  "solid",
  # Reaction -> substrate
   2,   10.93,   2,    10.38,  "dashed",
   6,   10.93,   6,    10.38,  "dashed",
  10,   10.93,  10,    10.38,  "dashed",
  # Substrates -> de novo (convergence)
   2,    9.63,   3.8,   8.8,   "solid",
   6,    9.63,   6.0,   8.8,   "solid",
  10,    9.63,   8.2,   8.8,   "solid",
  # De novo -> IMP
   6,    7.8,    6,     7.1,   "solid",
  # IMP -> AMP / GMP
   5.2,  6.5,    3.2,   5.78,  "solid",
   6.8,  6.5,    8.8,   5.78,  "solid",
  # AMP / GMP -> degradation
   3.2,  5.23,   3.2,   4.5,   "solid",
   8.8,  5.23,   8.8,   4.5,   "solid",
  # Degradation -> Xanthine
   3.2,  3.7,    5.1,   3.08,  "solid",
   8.8,  3.7,    6.9,   3.08,  "solid",
  # Xanthine -> UA
   6,    2.53,   6,     1.85,  "solid"
)

# --- 3. Build schematic -------------------------------------------------------
p <- ggplot() +
  geom_rect(
    data = nodes,
    aes(xmin = x - w / 2, xmax = x + w / 2,
        ymin = y - h / 2, ymax = y + h / 2),
    fill = nodes$fill, color = "grey50", linewidth = 0.3
  ) +
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label),
    size = nodes$fsize, fontface = nodes$fface,
    color = nodes$tcol, lineheight = 0.85
  ) +
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linetype = edges$lty,
    arrow = arrow(length = unit(0.07, "cm"), type = "closed"),
    linewidth = 0.35, color = "grey30"
  ) +
  # Step annotations on convergence arrows
  annotate("text", x = 2.7, y = 9.21,
           label = "Steps 3, 11\n(EC 2.1.2.2, EC 2.1.2.3)",
           size = GEOM_SECONDARY, color = "#2E7D32", fontface = "italic",
           hjust = 0, lineheight = 0.85) +
  annotate("text", x = 6.25, y = 9.21,
           label = "Step 9\n(EC 6.3.2.6)",
           size = GEOM_SECONDARY, color = "#E65100", fontface = "italic",
           hjust = 0, lineheight = 0.85) +
  annotate("text", x = 8.6, y = 9.21,
           label = "One-carbon\nrecycling",
           size = GEOM_SECONDARY, color = "#1565C0", fontface = "italic",
           hjust = 0, lineheight = 0.85) +
  # "Convergent" bracket (left margin)
  annotate("text", x = -0.6, y = 12.7,
           label = "Convergent\n(E1 & E4 ↑)",
           hjust = 0.5, size = GEOM_SECONDARY,
           fontface = "bold.italic", color = "#B71C1C") +
  annotate("segment", x = 0.15, y = 13.3, xend = 0.15, yend = 12.1,
           color = "#B71C1C", linewidth = 0.8,
           arrow = arrow(length = unit(0.06, "cm"), ends = "both",
                         type = "open")) +
  # Xanthine oxidase annotation
  annotate("text", x = 7.2, y = 2.15,
           label = "Xanthine oxidase\n(EC 1.17.3.2)",
           hjust = 0, size = GEOM_SECONDARY, fontface = "italic",
           color = "grey50", lineheight = 0.85) +
  # KEGG source footnote
  annotate("text", x = 6, y = 0.55,
           label = paste0("Based on KEGG Purine Metabolism (map00230), ",
                          "Module M00048, Folate Biosynthesis (map00790),\n",
                          "Amino Acid Metabolism (map00250), ",
                          "One-Carbon Pool by Folate (map00670)"),
           hjust = 0.5, size = GEOM_SECONDARY, color = "grey55",
           fontface = "italic", lineheight = 0.9) +
  coord_fixed(ratio = 0.45, xlim = c(-1.6, 12.6), ylim = c(0.0, 15.0)) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 2, 4))

# --- 4. Save ------------------------------------------------------------------
ggsave(file.path(out_fig, "FigureS5_Purine_Pathway_Context.pdf"), p,
       width = FIG4_W_FULL, height = FIG4_H_D, device = cairo_pdf)
ggsave(file.path(out_fig, "FigureS5_Purine_Pathway_Context.png"), p,
       width = FIG4_W_FULL, height = FIG4_H_D, dpi = FIG4_DPI, bg = "white")

cat(sprintf("Figure S7 (purine schematic) saved: %.2f x %.2f in\n",
            FIG4_W_FULL, FIG4_H_D))
