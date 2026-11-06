# ==============================================================================
# Figure 6D: Boxplots for four significant Type 4 discriminators
# ==============================================================================

library(ggplot2)
library(tidyverse)
library(ggdist)
library(here)

source(here("scripts", "figure6", "fig6_style.R"))

cat("--- Building Figure 6D ---\n")

out_dir <- here("results", "figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tbl35_dir <- here("results", "tables", "section_3_5")

# --- Shared config ---
# Neutral grey / deep-red contrast — deliberately NOT the enterotype hexes
# (#4DBBD5 = Type 2, #E64B35 = Type 1) to avoid cross-panel colour confusion.
hua_colors <- c("HUA-" = "#BDBDBD", "HUA+" = "#B2182B")

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

box_features <- c("MeOIA", "IndoxylSulfate", "EC_3_1_3_15", "EC_2_6_1_1")

# --- Load data ---
analysis_data <- readRDS(file.path(tbl35_dir, "mediation_analysis_data.rds"))

# het_results: load from CSV if available, otherwise compute on the fly
het_csv <- file.path(tbl35_dir, "type4_heterogeneity_results.csv")
if (file.exists(het_csv)) {
  het_results <- read.csv(het_csv, stringsAsFactors = FALSE)
} else {
  cat("type4_heterogeneity_results.csv not found, computing on the fly...\n")
  t4_tmp <- analysis_data %>% filter(Enterotype == "4", !is.na(HUA))
  het_results <- purrr::map_dfr(box_features, function(f) {
    w <- wilcox.test(t4_tmp[[f]] ~ t4_tmp$HUA)
    d <- rcompanion::cliffDelta(t4_tmp[[f]], t4_tmp$HUA)
    tibble(variable = f, wilcox_p = w$p.value, cliffs_delta = d)
  })
}

# --- Prepare Type 4 subset ---
t4 <- analysis_data %>%
  filter(Enterotype == "4", !is.na(HUA)) %>%
  mutate(HUA_group = factor(ifelse(HUA == 1, "HUA+", "HUA-"), levels = c("HUA-", "HUA+")))

# --- Build plot data ---
panel_d <- t4 %>%
  select(SampleID, HUA_group, all_of(box_features)) %>%
  pivot_longer(cols = all_of(box_features), names_to = "feature", values_to = "value") %>%
  mutate(feature_label = factor(mediator_display[feature], levels = mediator_display[box_features]))

panel_d_ann <- het_results %>%
  filter(variable %in% box_features) %>%
  mutate(
    feature_label = factor(mediator_display[variable], levels = mediator_display[box_features]),
    stat_label = sprintf("P = %.3f", wilcox_p)
  )

panel_d_pos <- panel_d %>%
  group_by(feature_label) %>%
  summarise(y_max = max(value, na.rm = TRUE), y_range = diff(range(value, na.rm = TRUE)), .groups = "drop") %>%
  mutate(y_annot = y_max + pmax(y_range * 0.12, 0.12))

panel_d_ann <- panel_d_ann %>% left_join(panel_d_pos, by = "feature_label")

set.seed(42)
panel_d <- panel_d %>%
  mutate(
    x_num = as.numeric(as.factor(HUA_group)),
    x_jit = x_num + runif(n(), -0.25, -0.08),
    slab_color = ifelse(HUA_group == "HUA+", "#E9A6B2", "#E3E3E3")
  )

# --- Plot (Raincloud style) ---
p6d <- ggplot(panel_d, aes(x = HUA_group, y = value)) +
  # Native ggplot2 violin: independent opaque identity fills, fully vector.
  geom_violin(
    aes(fill = slab_color),
    adjust = 1, width = 0.5, trim = TRUE, colour = NA
  ) +
  scale_fill_identity() +
  ggnewscale::new_scale_fill() +
  # Boxplot
  geom_boxplot(aes(fill = HUA_group), width = 0.15, linewidth = 0.5, outlier.shape = NA) +
  # Jittered points on left
  geom_point(aes(x = x_jit, fill = HUA_group), shape = 21, size = 2, color = "white") +
  geom_segment(
    data = panel_d_ann,
    aes(x = 1, xend = 2, y = y_annot - y_range * 0.04, yend = y_annot - y_range * 0.04),
    inherit.aes = FALSE,
    linewidth = 0.4,
    color = "grey35"
  ) +
  geom_text(
    data = panel_d_ann,
    aes(x = 1.5, y = y_annot, label = stat_label),
    inherit.aes = FALSE,
    size = as_mm(FS_ANNOT),
    lineheight = 0.9,
    color = "grey20"
  ) +
  # 2 x 2 facets (metabolites on top, ECs below) — fills the right column.
  facet_wrap(~ feature_label, scales = "free_y", ncol = 2) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.30))) +
  scale_fill_manual(values = hua_colors, name = "Type 4 subgroup") +
  labs(
    x = NULL,
    y = NULL,
    title = "Exploratory within-Type 4 comparison:\nfour nominally associated features"
  ) +
  theme_classic(base_size = FS_AXIS_TEXT) +
  fig6_text_theme() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )

# Rendered at the right column width (8.92 in) so magick performs no rescale.
# Height compressed to 10.3 in so the 2x2 facets are near-square (not stretched
# to match the taller left column).
ggsave(file.path(out_dir, "Figure6D_type4_boxplots_results.pdf"), p6d, width = 8.92, height = 10.3, dpi = 300)
ggsave(file.path(out_dir, "Figure6D_type4_boxplots_results.png"), p6d, width = 8.92, height = 10.3, dpi = 300)

saveRDS(p6d, here("results", "cache", "p6d.rds"))
cat("Figure 6D saved.\n")
