# ==============================================================================
# Figure 1B: PCoA Ordination with 95% Confidence Ellipses + PERMANOVA
# Bray-Curtis distance at genus level, colored by enterotype
# ==============================================================================

library(phyloseq)
library(vegan)
library(ggplot2)
library(here)

# --- 1. Load Data ---
# Use the updated phyloseq with labels re-derived from
# Panel A's seed=10 fit
ps <- readRDS(here(
  "data", "clean",
  "phyloseq_Taxa_with_enterotype_seed10.rds"
))

# Define enterotype factor
sample_data(ps)$Enterotype <- factor(
  sample_data(ps)$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)

source(here("scripts", "utils", "palettes.R"))
enterotype_colors <- ENTEROTYPE_COLORS

# --- 2. Compute Bray-Curtis Distance & PCoA ---
# Enterotype classification was based on genus-level abundance,
# so PCoA must also use genus-level data for consistency.
ps_genus <- tax_glom(ps, taxrank = "Genus")
dist_bc <- distance(ps_genus, method = "bray")
ord_pcoa <- ordinate(ps_genus, method = "PCoA", distance = dist_bc)

# Extract % variance explained for axes
eig_vals <- ord_pcoa$values$Eigenvalues
var_explained <- eig_vals / sum(eig_vals) * 100

# Extract sample coordinates
ord_df <- as.data.frame(ord_pcoa$vectors[, 1:2])
colnames(ord_df) <- c("PCo1", "PCo2")
ord_df$Enterotype <- sample_data(ps)$Enterotype

# Sample counts per enterotype
n_counts <- table(ord_df$Enterotype)
enterotype_labels <- paste0(names(n_counts), "\n(n = ", n_counts, ")")
names(enterotype_labels) <- names(n_counts)
ord_df$Enterotype_label <- enterotype_labels[as.character(ord_df$Enterotype)]

# --- 3. PERMANOVA ---
meta_df <- as(sample_data(ps), "data.frame")
set.seed(999)
perm_res <- adonis2(
  dist_bc ~ Enterotype,
  data = meta_df,
  permutations = 999,
  by = "margin"
)

r2_val <- round(perm_res$R2[1], 3)
p_val <- perm_res$`Pr(>F)`[1]
p_label <- ifelse(p_val < 0.001, "P < 0.001", paste0("P = ", p_val))

cat("PERMANOVA: R² =", r2_val, ",", p_label, "\n")

# --- 4. Plot ---
p1b <- ggplot(ord_df, aes(x = PCo1, y = PCo2, color = Enterotype, fill = Enterotype)) +
  # 95% confidence ellipses
  stat_ellipse(
    geom = "polygon", level = 0.95,
    alpha = 0.08, linewidth = 0.5, linetype = "solid"
  ) +
  # Sample points
  geom_point(size = 2, alpha = 0.7, shape = 16) +
  # Color scales
  scale_color_manual(values = enterotype_colors, labels = enterotype_labels) +
  scale_fill_manual(values = enterotype_colors, labels = enterotype_labels) +
  # PERMANOVA annotation
  annotate(
    "text",
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.5,
    label = paste0("PERMANOVA\nR\u00b2 = ", r2_val, ", ", p_label),
    size = 3.2, fontface = "italic", color = "grey30"
  ) +
  # Axis labels with variance explained
  labs(
    x = paste0("PCo1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PCo2 (", round(var_explained[2], 1), "%)"),
    color = "Enterotype",
    fill = "Enterotype"
  ) +
  # Equal aspect ratio: standard for PCoA — preserves distance relationships.
  # Variance importance is communicated via axis labels, not axis scaling.
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  )

# --- 5. Save ---
dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  here("results", "figures", "panels", "Figure1B_PCoA.pdf"),
  p1b,
  width = 6, height = 4.5, dpi = 300
)
ggsave(
  here("results", "figures", "panels", "Figure1B_PCoA.png"),
  p1b,
  width = 6, height = 4.5, dpi = 300
)

cat("\nFigure 1B saved to results/figures/panels/\n")

# Save plot object for later assembly
saveRDS(p1b, here("results", "cache", "p1b.rds"))
