# ==============================================================================
# Figure 6A: Full 22-candidate mediation forest (Type 1 / Type 4)
# ==============================================================================

library(ggplot2)
library(tidyverse)
library(here)

cat("--- Building Figure 6A ---\n")

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

all22_order <- c(
  rev(c("Adenosine", "TMAO", "DCA", "PAGln", "IPA", "IndoxylSulfate", "pCresolSulfate", "pCresol", "MeOIA")),
  rev(c("PWY7234", "ARGORNPROST", "PWY6609", "PWY7219", "PWY5004", "PWY0_1296", "PWY5695")),
  rev(c("EC_2_1_1_13", "EC_3_1_3_15", "EC_2_7_6_3", "EC_4_2_1_8", "EC_2_3_1_129", "EC_2_6_1_1"))
)

# --- Load data ---
med_results <- read.csv(file.path(tbl35_dir, "mediation_main_results.csv"), stringsAsFactors = FALSE)

# --- Plot ---
panel_a <- med_results %>%
  mutate(
    data_layer = case_when(
      mediator %in% metabolite_mediators ~ "Serum metabolite",
      mediator %in% pathway_mediators ~ "Microbial pathway",
      TRUE ~ "Convergent EC"
    ),
    data_layer = factor(data_layer, levels = c("Serum metabolite", "Microbial pathway", "Convergent EC")),
    feature_label = factor(mediator_display[mediator], levels = mediator_display[all22_order]),
    model_label = case_when(
      model == "Type1_vs_Ref" ~ "Type 1 vs Types 2+3",
      model == "Type4_vs_Ref" ~ "Type 4 vs Types 2+3",
      TRUE ~ model
    ),
    sig_state = ifelse(ACME_padj < 0.05, "BH-significant", "Not BH-significant")
  )

p6a <- ggplot(panel_a, aes(x = ACME_estimate, y = feature_label, color = data_layer)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey45") +
  geom_errorbar(aes(xmin = ACME_CI_lo, xmax = ACME_CI_hi), width = 0.2, linewidth = 0.5, orientation = "y") +
  geom_point(aes(shape = sig_state), size = 2.8, stroke = 0.8) +
  facet_wrap(~ model_label, scales = "free_x") +
  scale_color_manual(values = layer_colors[c("Serum metabolite", "Microbial pathway", "Convergent EC")], name = "Feature layer") +
  scale_shape_manual(values = c("BH-significant" = 16, "Not BH-significant" = 1), name = "ACME status") +
  labs(
    x = "ACME (delta UA, umol/L)",
    y = NULL,
    title = "Main mediation screen across 22 pre-specified candidates"
  ) +
  theme_bw(base_size = 10.5) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 7.5, color = "black"),
    axis.text.x = element_text(size = 8.5, color = "black"),
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

ggsave(file.path(out_dir, "Figure6A_mediation_forest_results.pdf"), p6a, width = 12, height = 8.2, dpi = 300)
ggsave(file.path(out_dir, "Figure6A_mediation_forest_results.png"), p6a, width = 12, height = 8.2, dpi = 300)

saveRDS(p6a, here("results", "cache", "p6a.rds"))
cat("Figure 6A saved.\n")
