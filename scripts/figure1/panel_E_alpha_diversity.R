# ==============================================================================
# Figure 1E: Alpha Diversity Across Enterotypes
# Shannon index & Observed species (Richness), side by side
# With Kruskal-Wallis global test + pairwise Wilcoxon (FDR)
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggdist)
library(patchwork)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean",
  "phyloseq_Taxa_with_enterotype_seed10.rds"))

sample_data(ps)$Enterotype <- factor(
  sample_data(ps)$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)

source(here("scripts", "utils", "palettes.R"))
enterotype_colors <- ENTEROTYPE_COLORS

# --- 2. Compute Alpha Diversity ---
alpha_df <- estimate_richness(ps, measures = c("Shannon", "Observed"))
alpha_df$Enterotype <- sample_data(ps)$Enterotype

# --- 3. Helper: Build One Panel (Raincloud style) ---
make_alpha_panel <- function(data, metric, y_label) {
  comparisons <- combn(levels(data$Enterotype), 2, simplify = FALSE)

  # Global test across all four enterotypes (Kruskal-Wallis)
  kw <- kruskal.test(data[[metric]], data$Enterotype)
  kw_label <- ifelse(
    kw$p.value < 0.001,
    "Kruskal-Wallis P < 0.001",
    sprintf("Kruskal-Wallis P = %.3f", kw$p.value)
  )
  cat(sprintf(
    "  %s: Kruskal-Wallis chi-squared = %.2f, df = %d, P = %.4g\n",
    metric, kw$statistic, kw$parameter, kw$p.value
  ))

  set.seed(42)
  data <- data %>%
    mutate(
      x_num = as.numeric(Enterotype),
      x_jit = x_num + runif(n(), -0.25, -0.08)
    )

  ggplot(data, aes(x = Enterotype, y = .data[[metric]], fill = Enterotype)) +
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
    # Comparisons
    stat_compare_means(
      comparisons = comparisons,
      method = "wilcox.test", p.adjust.method = "fdr",
      label = "p.signif",
      tip.length = 0.015, size = 3.2, vjust = 0.3,
      step.increase = 0.07
    ) +
    scale_fill_manual(values = enterotype_colors) +
    labs(x = NULL, y = y_label, subtitle = kw_label) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.title = element_text(size = 16, color = "black"),
      axis.text = element_text(size = 13, color = "black"),
      plot.subtitle = element_text(size = 12, face = "italic",
                                   color = "grey25", hjust = 0.5),
      plot.margin = margin(5, 8, 5, 5)
    )
}

# --- 4. Build Panels ---
p_shannon <- make_alpha_panel(alpha_df, "Shannon", "Shannon index")
p_observed <- make_alpha_panel(alpha_df, "Observed", "Observed species")

p1e <- p_shannon + p_observed +
  plot_layout(ncol = 2)

# --- 5. Save ---
dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  here("results", "figures", "panels", "Figure1E_Alpha_Diversity.pdf"),
  p1e, width = 12, height = 4, dpi = 300
)
ggsave(
  here("results", "figures", "panels", "Figure1E_Alpha_Diversity.png"),
  p1e, width = 12, height = 4, dpi = 300
)

cat("\nFigure 1E saved to results/figures/panels/\n")

saveRDS(p1e, here("results", "cache", "p1e.rds"))
