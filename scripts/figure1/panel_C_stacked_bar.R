# ==============================================================================
# Figure 1C: Group-Level Stacked Bar Chart (Mean Genus Composition)
# Top 10 genera + Others, per enterotype
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean",
  "phyloseq_Taxa_with_enterotype_seed10.rds"))

sample_data(ps)$Enterotype <- factor(
  sample_data(ps)$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)

# --- 2. Compute Mean Relative Abundance per Enterotype ---
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

# Identify top 10 genera by overall mean abundance
top10_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:10]

# Melt and aggregate
ps_melt <- psmelt(ps_rel)

# Clean genus names: remove "g__" prefix
ps_melt$Genus_Clean <- gsub("^g__", "", ps_melt$Genus)

bar_df <- ps_melt %>%
  mutate(Genus_Display = ifelse(OTU %in% top10_taxa, Genus_Clean, "Others")) %>%
  mutate(Genus_Display = gsub(
    "Lachnospiraceae_unclassified",
    "Lachnospiraceae (unclassified)",
    Genus_Display, fixed = TRUE
  )) %>%
  group_by(Sample, Enterotype, Genus_Display) %>%
  summarise(Abun = sum(Abundance), .groups = "drop") %>%
  group_by(Enterotype, Genus_Display) %>%
  summarise(
    Mean = mean(Abun),
    SE = sd(Abun) / sqrt(n()),
    .groups = "drop"
  )

# Order genera: Bacteroides and Prevotella first (key drivers), then alphabetical, Others last
key_genera <- c("Bacteroides", "Prevotella")
other_genera <- sort(setdiff(unique(bar_df$Genus_Display), c(key_genera, "Others")))
genus_order <- c(key_genera, other_genera, "Others")
bar_df$Genus_Display <- factor(bar_df$Genus_Display, levels = genus_order)

# --- 3. Color Palette ---
source(here("scripts", "utils", "palettes.R"))
genus_colors <- GENUS_COLORS

# Fallback: if genus not in palette, assign a grey tone
missing <- setdiff(genus_order, names(genus_colors))
if (length(missing) > 0) {
  extra_cols <- setNames(rep("#BDBDBD", length(missing)), missing)
  genus_colors <- c(genus_colors, extra_cols)
}

# --- 4. Plot ---
p1c <- ggplot(bar_df, aes(x = Enterotype, y = Mean, fill = Genus_Display)) +
  geom_col(position = "stack", width = 0.7, color = "white", linewidth = 0.15) +
  scale_fill_manual(values = genus_colors) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.005), labels = scales::percent) +
  labs(
    x = NULL,
    y = "Mean relative abundance",
    fill = "Genus"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(size = 10, color = "black"),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8, face = "italic"),  # italic for genus names
    legend.key.size = unit(0.35, "cm"),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 5)
  ) +
  # Override "Others" to non-italic in legend
  guides(fill = guide_legend(override.aes = list()))

# --- 5. Save ---
dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  here("results", "figures", "panels", "Figure1C_Stacked_Bar.pdf"),
  p1c, width = 5, height = 4, dpi = 300
)
ggsave(
  here("results", "figures", "panels", "Figure1C_Stacked_Bar.png"),
  p1c, width = 5, height = 4, dpi = 300
)

cat("\nFigure 1C saved to results/figures/panels/\n")

# Save plot object for later assembly
saveRDS(p1c, here("results", "cache", "p1c.rds"))
