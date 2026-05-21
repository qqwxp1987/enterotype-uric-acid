# ======================================================================
#
#  Section 3.3 — Panels G & H
#
#  Panel G: Species-level contributions to the 3 key convergent shared
#           ECs in Type 1 vs Type 4 (horizontal stacked bar, matching
#           Figure 3E style). Demonstrates "complete taxonomic
#           substitution" between the two high-HUA enterotypes.
#
#  Panel H: Simplified pathway schematic showing how the 3 convergent
#           ECs feed upstream substrates into de novo purine biosynthesis
#           → IMP → degradation → uric acid.
#           Evidence-based: KEGG map00230, map00790, map00250, map00670.
#
#  Inputs:
#    results/tables/section_3_3/attribution/ec_stratified_subset.rds
#    results/tables/section_3_3/attribution/prep_data.rds
#    results/tables/section_3_3/ec_shared_direction.csv
#
#  Outputs:
#    results/figures/panels/Figure3G_convergent_ec_species.{png,pdf}
#    results/figures/panels/Figure3H_purine_pathway.{png,pdf}
#    results/tables/section_3_3/convergent_ec_contribution_summary.csv
#
# ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(RColorBrewer)
  library(scales)
  library(here)
})

attr_dir <- here("results", "tables", "section_3_4", "attribution")
out_fig  <- here("results", "figures", "panels")
out_tbl  <- here("results", "tables", "section_3_4")
dir.create(out_fig, showWarnings = FALSE, recursive = TRUE)


# ######################################################################
#                        PANEL G
#    Species Contribution to Convergent ECs (Figure 3E style)
# ######################################################################

cat("==============================================================\n")
cat(" Panel G: Convergent EC Species Contribution\n")
cat("==============================================================\n\n")

# ------------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------------

prep <- readRDS(file.path(attr_dir, "prep_data.rds"))
meta <- prep$metadata   # SampleID, Enterotype (factor 1-4)

ec_dir_tbl <- read.csv(
  here("results", "tables", "section_3_4", "ec_shared_direction.csv"),
  stringsAsFactors = FALSE
)

convergent_ids <- ec_dir_tbl %>%
  filter(direction_type == "Convergent") %>%
  pull(feature)
cat(sprintf("  %d convergent shared ECs identified.\n", length(convergent_ids)))

# The 3 key convergent ECs (Section 3.3 text)
key_3 <- c("2.7.6.3", "2.6.1.1", "2.1.1.13")


# ------------------------------------------------------------------
# 2. Process stratified data
# ------------------------------------------------------------------

cat("  Loading EC stratified subset...\n")
raw <- readRDS(file.path(attr_dir, "ec_stratified_subset.rds"))
setnames(raw, 1, "Feature_Raw")
sample_cols <- setdiff(names(raw), "Feature_Raw")

raw_strat <- raw[grepl("|", Feature_Raw, fixed = TRUE)]
raw_strat[, ec_id := sub("[|].*", "", Feature_Raw)]
raw_strat <- raw_strat[ec_id %in% convergent_ids]

long_dt <- melt(
  raw_strat,
  id.vars       = c("Feature_Raw", "ec_id"),
  measure.vars  = sample_cols,
  variable.name = "SampleID_raw",
  value.name    = "RPK"
)[RPK > 0]

long_dt[, `:=`(
  species_raw = sub("^[^|]+\\|", "", Feature_Raw),
  SampleID    = gsub("-", ".", as.character(SampleID_raw))
)]

long_dt <- long_dt[species_raw != "unclassified"]
long_dt <- long_dt[!grepl("UNGROUPED|UNINTEGRATED|UNMAPPED", species_raw)]

# Clean species name
long_dt[, species := {
  s <- species_raw
  s <- sub("^g__[^.]+\\.", "", s)
  s <- sub("^s__", "", s)
  gsub("_", " ", s)
}]

cat(sprintf("  %d non-zero (EC, species, sample) triplets.\n", nrow(long_dt)))


# ------------------------------------------------------------------
# 3. Compute contributions (sum-based) and join metadata
# ------------------------------------------------------------------

long_dt[, total_RPK := sum(RPK), by = .(SampleID, ec_id)]
long_dt[, contribution_pct := RPK / total_RPK * 100]

contrib_df <- as.data.frame(long_dt) %>%
  inner_join(meta, by = "SampleID")

# Summarise: sum(RPK) / total per (EC, Enterotype, species)
et_totals <- contrib_df %>%
  group_by(ec_id, Enterotype) %>%
  summarise(et_total = sum(RPK), .groups = "drop")

mean_contrib <- contrib_df %>%
  group_by(ec_id, Enterotype, species) %>%
  summarise(sp_total = sum(RPK), .groups = "drop") %>%
  left_join(et_totals, by = c("ec_id", "Enterotype")) %>%
  mutate(mean_pct = sp_total / et_total * 100) %>%
  select(-sp_total, -et_total)

write.csv(mean_contrib,
          file.path(out_tbl, "convergent_ec_contribution_summary.csv"),
          row.names = FALSE)
cat("  Full summary saved.\n")


# ------------------------------------------------------------------
# 4. Build scatter-plot data: E1 contribution vs E4 contribution
#    Style: species as points, diagonal = equal contribution,
#    colored by dominance (Type 1 / Type 4 / Shared)
# ------------------------------------------------------------------

library(ggrepel)

key_contrib <- mean_contrib %>%
  filter(ec_id %in% key_3, Enterotype %in% c("1", "4"))

# Pivot to wide: one row per (ec_id, species) with E1 and E4 columns
scatter_df <- key_contrib %>%
  mutate(ET = paste0("E", Enterotype)) %>%
  select(ec_id, species, ET, mean_pct) %>%
  pivot_wider(names_from = ET, values_from = mean_pct,
              values_fill = 0)

# Classify dominance
scatter_df <- scatter_df %>%
  mutate(
    Dominant = case_when(
      E4 > E1 * 1.5 ~ "Type 4",
      E1 > E4 * 1.5 ~ "Type 1",
      TRUE           ~ "Shared"
    ),
    Dominant = factor(Dominant, levels = c("Type 1", "Shared", "Type 4"))
  )

# Keep only species with >= 1% contribution in at least one enterotype
scatter_df <- scatter_df %>%
  filter(E1 >= 1 | E4 >= 1)

# EC facet labels
ec_facet <- c(
  "2.7.6.3"  = "HPPK (EC 2.7.6.3)",
  "2.6.1.1"  = "Aspartate transaminase (EC 2.6.1.1)",
  "2.1.1.13" = "Methionine synthase (EC 2.1.1.13)"
)
scatter_df$EC_label <- factor(ec_facet[scatter_df$ec_id],
                              levels = ec_facet[key_3])

# Dominance colors — aligned with the global enterotype palette
# (scripts/utils/palettes.R): Type 1 = red, Type 4 = navy blue.
dom_colors <- c("Type 1" = "#E64B35",   # Type 1 red
                "Shared" = "#999999",    # grey
                "Type 4" = "#3C5488")    # Type 4 navy blue

# Axis limit (symmetric)
ax_max <- max(c(scatter_df$E1, scatter_df$E4)) * 1.08


# ------------------------------------------------------------------
# 5. Plot Panel G — scatter plot
# ------------------------------------------------------------------

set.seed(42)  # stabilise ggrepel label placement across runs
p_G <- ggplot(scatter_df, aes(x = E1, y = E4, color = Dominant)) +
  # Diagonal reference line
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey60", linewidth = 0.4) +
  # Points
  geom_point(size = 2.5, alpha = 0.85) +
  # Species labels
  geom_text_repel(
    aes(label = species),
    size = 2.5, fontface = "italic",
    max.overlaps = 20, segment.size = 0.25,
    segment.color = "grey70",
    box.padding = 0.3, point.padding = 0.2,
    min.segment.length = 0.2
  ) +
  facet_wrap(~ EC_label, nrow = 1, scales = "free") +
  scale_color_manual(values = dom_colors, name = "Dominant in") +
  scale_x_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = 0.05)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = 0.05)) +
  labs(x = "Mean contribution in Type 1",
       y = "Mean contribution in Type 4") +
  theme_bw(base_size = 11) +
  theme(
    panel.background  = element_rect(fill = "white"),
    plot.background   = element_rect(fill = "white", color = NA),
    strip.background  = element_rect(fill = "grey96", color = "grey70"),
    strip.text        = element_text(size = 9, face = "bold"),
    axis.title        = element_text(size = 10),
    axis.text         = element_text(size = 9),
    legend.position   = "bottom",
    legend.text       = element_text(size = 9),
    legend.title      = element_text(size = 10, face = "bold"),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(10, 10, 10, 10)
  )

ggsave(file.path(out_fig, "Figure3G_convergent_ec_species.png"),
       p_G, width = 13, height = 5, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "Figure3G_convergent_ec_species.pdf"),
       p_G, width = 13, height = 5, device = cairo_pdf)
cat("  Panel G saved (png + pdf).\n")


# ######################################################################
#                        PANEL H
#    Purine Pathway Schematic (evidence-based, KEGG)
# ######################################################################

cat("\n==============================================================\n")
cat(" Panel H: Purine Pathway Schematic (KEGG-based)\n")
cat("==============================================================\n\n")

# ------------------------------------------------------------------
# Evidence sources (KEGG REST API, accessed 2025):
#   EC 2.7.6.3  → KEGG ec00790 (Folate biosynthesis)
#     Reaction: ATP + 6-hydroxymethyl-7,8-dihydropterin
#               → AMP + 6-hydroxymethyl-7,8-dihydropterin diphosphate
#
#   EC 2.6.1.1  → KEGG ec00250 (Ala/Asp/Glu metabolism)
#     Reaction: L-aspartate + 2-oxoglutarate
#               ↔ oxaloacetate + L-glutamate
#
#   EC 2.1.1.13 → KEGG ec00670 (One carbon pool by folate)
#                + KEGG ec00270 (Cys/Met metabolism)
#     Reaction: 5-methyl-THF + L-homocysteine
#               → THF + L-methionine
#
#   De novo IMP biosynthesis (KEGG M00048, map00230):
#     Step 3  (EC 2.1.2.2):  GAR + 10-formyl-THF → FGAR
#     Step 9  (EC 6.3.2.6):  SAICAR synthase, uses L-aspartate
#     Step 11 (EC 2.1.2.3):  AICAR + 10-formyl-THF → FAICAR
#
#   IMP degradation → uric acid (KEGG map00230):
#     IMP → XMP → GMP → guanosine → guanine → xanthine → UA
#     IMP → adenylosuccinate → AMP → adenosine → inosine
#           → hypoxanthine → xanthine → UA
#     Xanthine oxidase: EC 1.17.3.2
# ------------------------------------------------------------------

# Node layout — vertical flow, x ∈ [0, 12], y ∈ [0, 15]
nodes <- tribble(
  ~label, ~x, ~y, ~w, ~h, ~fill, ~fsize, ~fface, ~tcol,

  # Row 1: Upstream pathway context (y = 14.2)
  "Folate biosynthesis\n(KEGG map00790)",
      2,  14.2, 3.0, 0.85, "#A5D6A7", 3.0, "bold", "grey20",
  "Ala / Asp / Glu metabolism\n(KEGG map00250)",
      6,  14.2, 3.0, 0.85, "#FFCC80", 3.0, "bold", "grey20",
  "One-carbon / Met metabolism\n(KEGG map00670, map00270)",
     10,  14.2, 3.0, 0.85, "#90CAF9", 2.7, "bold", "grey20",

  # Row 2: Convergent ECs (y = 12.7) — highlighted
  "EC 2.7.6.3\nHPPK",
      2, 12.7, 3.0, 0.85, "#2E7D32", 3.5, "bold", "white",
  "EC 2.6.1.1\nAspartate transaminase",
      6, 12.7, 3.0, 0.85, "#E65100", 3.2, "bold", "white",
  "EC 2.1.1.13\nMethionine synthase",
     10, 12.7, 3.0, 0.85, "#1565C0", 3.2, "bold", "white",

  # Row 3: Reaction / products (y = 11.3)
  "6-HMDP + ATP\n\u2192 6-HMDP-PP",
      2, 11.3, 3.0, 0.75, "#C8E6C9", 2.6, "plain", "grey25",
  "L-Asp + 2-OG\n\u2194 OAA + L-Glu",
      6, 11.3, 3.0, 0.75, "#FFE0B2", 2.6, "plain", "grey25",
  "5-mTHF + L-Hcy\n\u2192 THF + L-Met",
     10, 11.3, 3.0, 0.75, "#BBDEFB", 2.6, "plain", "grey25",

  # Row 4: Downstream substrates entering IMP biosynthesis (y = 10)
  "\u2192 THF \u2192 10-formyl-THF\n(purine ring C\u2082, C\u2088)",
      2, 10.0, 3.0, 0.75, "#E8F5E9", 2.5, "italic", "grey35",
  "L-Aspartate\n(purine ring N\u2081)",
      6, 10.0, 2.4, 0.75, "#FFF3E0", 2.5, "italic", "grey35",
  "\u2192 L-Met \u2192 SAM\n(one-carbon / adenosine)",
     10, 10.0, 3.0, 0.75, "#E3F2FD", 2.5, "italic", "grey35",

  # Row 5: De novo IMP biosynthesis (y = 8.3)
  "De novo purine biosynthesis (KEGG M00048)\nPRPP + Gln + Gly + CO\u2082 \u2192 (12 steps) \u2192 IMP",
      6,  8.3, 7.0, 1.0, "#CE93D8", 3.2, "bold", "grey10",

  # Row 5b: Step annotations (placed as text, no box)

  # Row 6: IMP (y = 6.8)
  "IMP",
      6, 6.8, 1.6, 0.6, "#E1BEE7", 3.8, "bold", "grey10",

  # Row 7: Nucleotides (y = 5.5)
  "AMP",
    3.2, 5.5, 1.5, 0.55, "#F8BBD0", 3.2, "bold", "grey20",
  "XMP \u2192 GMP",
    8.8, 5.5, 2.0, 0.55, "#F8BBD0", 3.2, "bold", "grey20",

  # Row 8: Degradation intermediates (y = 4.2)
  "Adenosine \u2192 Inosine\n\u2192 Hypoxanthine",
    3.2, 4.1, 2.2, 0.8, "#F5F5F5", 2.5, "plain", "grey40",
  "Guanosine\n\u2192 Guanine",
    8.8, 4.1, 2.0, 0.8, "#F5F5F5", 2.5, "plain", "grey40",

  # Row 9: Xanthine (y = 2.8)
  "Xanthine",
      6, 2.8, 1.8, 0.55, "#FFF176", 3.2, "bold", "grey20",

  # Row 10: Uric Acid (y = 1.5)
  "Uric Acid",
      6, 1.5, 2.2, 0.7, "#E53935", 4.0, "bold", "white"
)

edges <- tribble(
  ~x,    ~y,    ~xend,  ~yend,  ~lty,

  # Context → EC
   2,   13.78,   2,    13.13,  "solid",
   6,   13.78,   6,    13.13,  "solid",
  10,   13.78,  10,    13.13,  "solid",

  # EC → reaction
   2,   12.28,   2,    11.68,  "solid",
   6,   12.28,   6,    11.68,  "solid",
  10,   12.28,  10,    11.68,  "solid",

  # Reaction → substrate
   2,   10.93,   2,    10.38,  "dashed",
   6,   10.93,   6,    10.38,  "dashed",
  10,   10.93,  10,    10.38,  "dashed",

  # Substrates → De novo (convergence)
   2,    9.63,   3.8,   8.8,   "solid",
   6,    9.63,   6.0,   8.8,   "solid",
  10,    9.63,   8.2,   8.8,   "solid",

  # De novo → IMP
   6,    7.8,    6,     7.1,   "solid",

  # IMP → AMP / GMP
   5.2,  6.5,    3.2,   5.78,  "solid",
   6.8,  6.5,    8.8,   5.78,  "solid",

  # AMP / GMP → degradation
   3.2,  5.23,   3.2,   4.5,   "solid",
   8.8,  5.23,   8.8,   4.5,   "solid",

  # Degradation → Xanthine
   3.2,  3.7,    5.1,   3.08,  "solid",
   8.8,  3.7,    6.9,   3.08,  "solid",

  # Xanthine → UA
   6,    2.53,   6,     1.85,  "solid"
)


# ------------------------------------------------------------------
# Build schematic
# ------------------------------------------------------------------

p_H <- ggplot() +
  # Boxes
  geom_rect(
    data = nodes,
    aes(xmin = x - w/2, xmax = x + w/2,
        ymin = y - h/2, ymax = y + h/2),
    fill = nodes$fill, color = "grey50", linewidth = 0.3
  ) +
  # Box labels
  geom_text(
    data = nodes,
    aes(x = x, y = y, label = label),
    size       = nodes$fsize,
    fontface   = nodes$fface,
    color      = nodes$tcol,
    lineheight = 0.85
  ) +
  # Arrows
  geom_segment(
    data = edges,
    aes(x = x, y = y, xend = xend, yend = yend),
    linetype  = edges$lty,
    arrow     = arrow(length = unit(0.12, "cm"), type = "closed"),
    linewidth = 0.4,
    color     = "grey30"
  ) +
  # Step annotations on convergence arrows
  annotate("text", x = 2.6, y = 9.35,
           label = "Steps 3, 11\n(EC 2.1.2.2,\nEC 2.1.2.3)",
           size = 2.0, color = "#2E7D32", fontface = "italic",
           hjust = 0, lineheight = 0.8) +
  annotate("text", x = 6.3, y = 9.35,
           label = "Step 9\n(EC 6.3.2.6)",
           size = 2.0, color = "#E65100", fontface = "italic",
           hjust = 0, lineheight = 0.8) +
  annotate("text", x = 8.5, y = 9.35,
           label = "One-carbon\nrecycling",
           size = 2.0, color = "#1565C0", fontface = "italic",
           hjust = 0, lineheight = 0.8) +
  # "Convergent" bracket (left margin)
  annotate("text", x = -0.5, y = 12.7,
           label = "Convergent\n(E1 & E4 \u2191)",
           hjust = 0.5, size = 2.5,
           fontface = "bold.italic", color = "#B71C1C") +
  annotate("segment",
           x = 0.2, y = 13.3, xend = 0.2, yend = 12.1,
           color = "#B71C1C", linewidth = 0.9,
           arrow = arrow(length = unit(0.08, "cm"),
                         ends = "both", type = "open")) +
  # Xanthine oxidase annotation
  annotate("text", x = 7.2, y = 2.15,
           label = "Xanthine oxidase\n(EC 1.17.3.2)",
           hjust = 0, size = 2.0, fontface = "italic", color = "grey50",
           lineheight = 0.85) +
  # KEGG source footnote
  annotate("text", x = 6, y = 0.5,
           label = "Based on KEGG Purine Metabolism (map00230), Module M00048,\nFolate Biosynthesis (map00790), Amino Acid Metabolism (map00250),\nOne-Carbon Pool by Folate (map00670)",
           hjust = 0.5, size = 2.0, color = "grey55",
           fontface = "italic", lineheight = 0.85) +
  coord_fixed(ratio = 0.5,
              xlim = c(-1.5, 12.5),
              ylim = c(0, 15.2)) +
  theme_void() +
  theme(plot.margin = margin(5, 5, 5, 5))

ggsave(file.path(out_fig, "Figure3H_purine_pathway.png"),
       p_H, width = 9, height = 10, dpi = 300, bg = "white")
ggsave(file.path(out_fig, "Figure3H_purine_pathway.pdf"),
       p_H, width = 9, height = 10, device = cairo_pdf)
cat("  Panel H saved (png + pdf).\n")


cat("\n==============================================================\n")
cat(" Panels G & H Complete\n")
cat("==============================================================\n")
