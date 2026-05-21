# ==============================================================================
# Figure 6C: Type 4 heterogeneity lollipop across all 23 pre-specified features
# ==============================================================================

library(ggplot2)
library(tidyverse)
library(here)

cat("--- Building Figure 6C ---\n")

out_dir <- here("results", "figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tbl35_dir <- here("results", "tables", "section_3_5")

# --- Shared config ---
# Omics-layer colors: single source of truth in scripts/utils/palettes.R
source(here("scripts", "utils", "palettes.R"))
layer_colors <- LAYER_COLORS

metabolite_mediators <- c(
  "MeOIA", "pCresol", "pCresolSulfate", "IndoxylSulfate",
  "IPA", "PAGln", "DCA", "TMAO", "Adenosine"
)
pathway_mediators <- c(
  "PWY5695", "PWY0_1296", "PWY5004", "PWY7219",
  "PWY6609", "ARGORNPROST", "PWY7234"
)
ec_mediators <- c(
  "EC_2_6_1_1", "EC_2_3_1_129", "EC_4_2_1_8",
  "EC_2_7_6_3", "EC_3_1_3_15", "EC_2_1_1_13"
)

mediator_display <- c(
  MeOIA = "5-MeOIA", pCresol = "p-Cresol",
  pCresolSulfate = "p-Cresol sulfate", IndoxylSulfate = "Indoxyl sulfate",
  IPA = "Indole-3-propionic acid", PAGln = "Phenylacetylglutamine",
  DCA = "Deoxycholic acid", TMAO = "TMAO", Adenosine = "Adenosine",
  PWY5695 = "PWY-5695 (urate biosynthesis)",
  PWY0_1296 = "PWY0-1296 (purine rib. degradation)",
  PWY5004 = "PWY-5004 (citrulline metabolism)",
  PWY7219 = "PWY-7219 (adenosine de novo)",
  PWY6609 = "PWY-6609 (adenine salvage)",
  ARGORNPROST = "ARGORNPROST-PWY (Arg/Orn/Pro)",
  PWY7234 = "PWY-7234 (IMP biosynthesis)",
  EC_2_6_1_1 = "EC 2.6.1.1", EC_2_3_1_129 = "EC 2.3.1.129",
  EC_4_2_1_8 = "EC 4.2.1.8", EC_2_7_6_3 = "EC 2.7.6.3",
  EC_3_1_3_15 = "EC 3.1.3.15", EC_2_1_1_13 = "EC 2.1.1.13",
  P_copri = "P. copri"
)

# --- Load data ---
het_results <- read.csv(file.path(tbl35_dir, "type4_heterogeneity_results.csv"), stringsAsFactors = FALSE)

# --- Plot ---
panel_c <- het_results %>%
  mutate(
    display_label = mediator_display[variable],
    data_layer = case_when(
      variable %in% metabolite_mediators ~ "Serum metabolite",
      variable %in% pathway_mediators ~ "Microbial pathway",
      variable %in% ec_mediators ~ "Convergent EC",
      TRUE ~ "Taxonomy"
    ),
    data_layer = factor(data_layer, levels = c("Serum metabolite", "Microbial pathway", "Convergent EC", "Taxonomy")),
    is_sig = wilcox_p < 0.05
  ) %>%
  arrange(cliffs_delta) %>%
  mutate(display_label = factor(display_label, levels = display_label))

p6c <- ggplot(panel_c, aes(x = cliffs_delta, y = display_label, color = data_layer)) +
  geom_vline(xintercept = c(-0.474, -0.33, 0.33, 0.474), linetype = "dashed", linewidth = 0.3, color = "grey75") +
  geom_vline(xintercept = 0, linewidth = 0.5, color = "grey40") +
  geom_segment(aes(x = 0, xend = cliffs_delta, yend = display_label), linewidth = 0.65) +
  geom_point(aes(shape = is_sig), size = 2.7, stroke = 0.8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1), name = "Wilcoxon P") +
  scale_color_manual(values = layer_colors, name = "Feature layer") +
  labs(
    x = "Cliff's delta (HUA+ vs HUA-)",
    y = NULL,
    title = "Type 4 heterogeneity is concentrated in a limited subset of features"
  ) +
  theme_bw(base_size = 10.2) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 7.2, color = "black"),
    axis.text.x = element_text(size = 8.2, color = "black"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

ggsave(file.path(out_dir, "Figure6C_type4_lollipop_results.pdf"), p6c, width = 8.8, height = 7.6, dpi = 300)
ggsave(file.path(out_dir, "Figure6C_type4_lollipop_results.png"), p6c, width = 8.8, height = 7.6, dpi = 300)

saveRDS(p6c, here("results", "cache", "p6c.rds"))
cat("Figure 6C saved.\n")
