# ==============================================================================
# Supplementary Figure S1: Individual-Level Stacked Bar Chart
# Each bar = one sample, grouped by enterotype, sorted by dominant genus within
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

# --- 2. Prepare Relative Abundance Data ---
ps_genus <- tax_glom(ps, taxrank = "Genus")
ps_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))

top10_taxa <- names(sort(taxa_sums(ps_rel), decreasing = TRUE))[1:10]
ps_melt <- psmelt(ps_rel)

# Clean genus names
ps_melt$Genus_Clean <- gsub("^g__", "", ps_melt$Genus)

ind_df <- ps_melt %>%
  mutate(Genus_Display = ifelse(OTU %in% top10_taxa, Genus_Clean, "Others")) %>%
  mutate(Genus_Display = gsub(
    "Lachnospiraceae_unclassified",
    "Lachnospiraceae (unclassified)",
    Genus_Display, fixed = TRUE
  )) %>%
  group_by(Sample, Enterotype, Genus_Display) %>%
  summarise(Abun = sum(Abundance), .groups = "drop")

# --- 3. Sort Samples Within Each Enterotype ---
# Sort by the dominant genus of each enterotype:
#   Type 1/2: by Bacteroides (descending)
#   Type 3: by Prevotella (descending, co-dominant type)
#   Type 4: by Prevotella (descending)
sort_key <- ind_df %>%
  filter(Genus_Display %in% c("Bacteroides", "Prevotella")) %>%
  pivot_wider(names_from = Genus_Display, values_from = Abun, values_fill = 0) %>%
  mutate(
    Sort_Value = case_when(
      Enterotype %in% c("Type 1", "Type 2") ~ -Bacteroides,  # descending
      Enterotype %in% c("Type 3", "Type 4") ~ -Prevotella
    )
  ) %>%
  arrange(Enterotype, Sort_Value)

# Assign ordered sample factor
sample_order <- sort_key$Sample
ind_df$Sample <- factor(ind_df$Sample, levels = sample_order)

# Genus order & colors (same as Panel C)
key_genera <- c("Bacteroides", "Prevotella")
other_genera <- sort(setdiff(unique(ind_df$Genus_Display), c(key_genera, "Others")))
genus_order <- c(key_genera, other_genera, "Others")
ind_df$Genus_Display <- factor(ind_df$Genus_Display, levels = genus_order)

source(here("scripts", "utils", "palettes.R"))
genus_colors <- GENUS_COLORS
missing <- setdiff(genus_order, names(genus_colors))
if (length(missing) > 0) {
  genus_colors <- c(genus_colors, setNames(rep("#BDBDBD", length(missing)), missing))
}

# --- 4. Plot ---
# Enterotype labels with sample counts for facet strips
n_counts <- table(sample_data(ps)$Enterotype)
facet_labels <- setNames(
  paste0(names(n_counts), " (n = ", n_counts, ")"),
  names(n_counts)
)

p_s1 <- ggplot(ind_df, aes(x = Sample, y = Abun, fill = Genus_Display)) +
  geom_col(width = 1, linewidth = 0) +
  facet_grid(
    ~ Enterotype,
    scales = "free_x", space = "free_x",
    labeller = labeller(Enterotype = facet_labels)
  ) +
  scale_fill_manual(values = genus_colors) +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent) +
  labs(
    x = "Samples",
    y = "Relative abundance",
    fill = "Genus"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10, face = "bold"),
    strip.background = element_rect(fill = "grey95", color = NA),
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8, face = "italic"),
    legend.key.size = unit(0.35, "cm"),
    legend.position = "bottom",
    panel.spacing = unit(0.3, "lines"),
    plot.margin = margin(5, 10, 5, 5)
  ) +
  guides(fill = guide_legend(nrow = 2))

# --- 5. Save ---
dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)

ggsave(
  here("results", "figures", "panels", "FigureS1_Individual_Bar.pdf"),
  p_s1, width = 12, height = 4, dpi = 300
)
ggsave(
  here("results", "figures", "panels", "FigureS1_Individual_Bar.png"),
  p_s1, width = 12, height = 4, dpi = 300
)

cat("\nSupplementary Figure S1 saved to results/figures/panels/\n")

saveRDS(p_s1, here("results", "cache", "p_s1.rds"))
