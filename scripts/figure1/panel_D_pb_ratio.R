# ==============================================================================
# Figure 1D: Prevotella/Bacteroides Ratio (Violin + Boxplot)
# With all pairwise Wilcoxon comparisons (FDR-corrected)
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggdist)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here(
  "data", "clean",
  "phyloseq_Taxa_with_enterotype_seed10.rds"
))

sample_data(ps)$Enterotype <- factor(
  sample_data(ps)$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)

source(here("scripts", "utils", "palettes.R"))
enterotype_colors <- ENTEROTYPE_COLORS

# --- 2. Compute P/B Ratio ---
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
ps_melt <- psmelt(ps_rel)

pb_ratio <- ps_melt %>%
  mutate(Genus_Clean = gsub("^g__", "", Genus)) %>%
  filter(Genus_Clean %in% c("Bacteroides", "Prevotella")) %>%
  group_by(Sample, Enterotype, Genus_Clean) %>%
  summarise(Abun = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Genus_Clean, values_from = Abun, values_fill = 0) %>%
  mutate(
    # Add pseudocount to avoid log(0); must be < half of min non-zero
    # abundance (3.7e-6), so use 1e-6
    Log_PB_Ratio = log10((Prevotella + 1e-6) / (Bacteroides + 1e-6))
  )

# --- 3. Pairwise Wilcoxon Tests ---
comparisons <- combn(levels(pb_ratio$Enterotype), 2, simplify = FALSE)

# --- 4. Plot (Raincloud style) ---
set.seed(42)
pb_ratio <- pb_ratio %>%
  mutate(
    x_num = as.numeric(Enterotype),
    x_jit = x_num + runif(n(), -0.25, -0.08)
  )

p1d <- ggplot(pb_ratio, aes(x = Enterotype, y = Log_PB_Ratio, fill = Enterotype)) +
  # Half violin on right
  stat_halfeye(
    side = "right", adjust = 1, width = 0.5,
    .width = 0, point_colour = NA,
    justification = -0.12, alpha = 0.35, colour = NA,
    normalize = "groups"
  ) +
  # Boxplot
  geom_boxplot(width = 0.15, linewidth = 0.5, outlier.shape = NA) +
  # Jittered points on left
  geom_point(aes(x = x_jit), shape = 21, size = 2, color = "white") +
  # Pairwise comparisons
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test", p.adjust.method = "fdr",
    label = "p.signif", # *, **, ***, ns
    tip.length = 0.015,
    size = 3.2, vjust = 0.3,
    step.increase = 0.07
  ) +
  scale_fill_manual(values = enterotype_colors) +
  labs(
    x = NULL,
    y = expression(log[10] * "(" * italic(Prevotella) * " / " * italic(Bacteroides) * ")")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    plot.margin = margin(5, 10, 5, 5)
  )

# --- 5. Save ---
dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  here("results", "figures", "panels", "Figure1D_PB_Ratio.pdf"),
  p1d,
  width = 5, height = 4, dpi = 300
)
ggsave(
  here("results", "figures", "panels", "Figure1D_PB_Ratio.png"),
  p1d,
  width = 5, height = 4, dpi = 300
)

cat("\nFigure 1D saved to results/figures/panels/\n")

saveRDS(p1d, here("results", "cache", "p1d.rds"))
