# ==============================================================================
# Figure 2B: Serum UA Distribution by Enterotype
# Violin + Boxplot + Jitter with pairwise significance brackets
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(ggpubr)
library(ggdist)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta <- as(sample_data(ps), "data.frame")

meta$Enterotype <- factor(
  meta$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)
meta$UA <- as.numeric(meta$UA)

# --- 2. Consistent Color Palette ---
enterotype_colors <- c(
  "Type 1" = "#E64B35",
  "Type 2" = "#4DBBD5",
  "Type 3" = "#00A087",
  "Type 4" = "#3C5488"
)

# --- 3. Pairwise Comparisons ---
# Use Wilcoxon test with Holm correction
comparisons <- list(
  c("Type 1", "Type 2"),
  c("Type 1", "Type 3"),
  c("Type 1", "Type 4"),
  c("Type 2", "Type 3"),
  c("Type 2", "Type 4"),
  c("Type 3", "Type 4")
)

# Pre-compute to identify significant pairs only
pairwise_res <- pairwise.wilcox.test(meta$UA, meta$Enterotype, p.adjust.method = "holm")
cat("=== Pairwise Wilcoxon (Holm-corrected) ===\n")
print(pairwise_res$p.value)

# Build significant-pair bracket table from the Holm-corrected matrix.
# NOTE: stat_compare_means(comparisons = ) ignores p.adjust.method and would
# label brackets with RAW p-values; we therefore draw brackets explicitly
# from the Holm-corrected matrix so they match the text and figure legend.
pmat <- pairwise_res$p.value
bracket_df <- expand.grid(
  group2 = rownames(pmat),
  group1 = colnames(pmat),
  stringsAsFactors = FALSE
) %>%
  mutate(p_holm = mapply(function(r, c) pmat[r, c], group2, group1)) %>%
  filter(!is.na(p_holm), p_holm < 0.05) %>%
  arrange(p_holm) %>%
  mutate(
    label = case_when(
      p_holm < 0.001 ~ "***",
      p_holm < 0.01  ~ "**",
      TRUE           ~ "*"
    ),
    y.position = max(meta$UA, na.rm = TRUE) *
      (1.05 + 0.08 * (row_number() - 1))
  )
cat("\nSignificant pairs (Holm-corrected):", nrow(bracket_df), "\n")
print(bracket_df[, c("group1", "group2", "p_holm", "label")])

# --- 4. Summary stats for annotation ---
ua_summary <- meta %>%
  group_by(Enterotype) %>%
  summarise(
    n = n(),
    mean_UA = round(mean(UA, na.rm = TRUE), 1),
    median_UA = round(median(UA, na.rm = TRUE), 1),
    .groups = "drop"
  )
cat("\nUA summary:\n")
print(as.data.frame(ua_summary))

# --- 5. Plot (Raincloud style) ---
set.seed(42)

meta <- meta %>%
  mutate(
    x_num = as.numeric(Enterotype),
    x_jit = x_num + runif(n(), -0.25, -0.08)
  )

p <- ggplot(meta, aes(x = Enterotype, y = UA, fill = Enterotype)) +
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
  # Significance brackets (Holm-corrected; only significant pairs)
  ggpubr::stat_pvalue_manual(
    bracket_df,
    label = "label",
    xmin = "group1", xmax = "group2",
    tip.length = 0.015,
    bracket.size = 0.4,
    size = 3.5
  ) +
  # N labels at bottom
  geom_text(
    data = ua_summary,
    aes(x = Enterotype, y = -Inf, label = paste0("n=", n)),
    vjust = -0.5, size = 3, color = "grey40",
    inherit.aes = FALSE
  ) +
  # Scales
  scale_fill_manual(values = enterotype_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.08, 0.22))
  ) +
  labs(
    x = NULL,
    y = expression("Serum UA (" * mu * "mol/L)")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 16, color = "black"),
    axis.text = element_text(size = 13, color = "black"),
    plot.margin = margin(10, 10, 10, 10)
  )

# --- 6. Save ---
out_dir <- here("results", "figures", "panels")

ggsave(file.path(out_dir, "Figure2B_UA_Violin.pdf"),
       p, width = 4.5, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure2B_UA_Violin.png"),
       p, width = 4.5, height = 5, dpi = 300)

cat("\nFigure 2B saved.\n")

# --- 7. Save full pairwise Wilcoxon results ---
table_dir <- here("results", "tables", "section_3_2")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

pw_mat <- pairwise_res$p.value
pw_long <- expand.grid(
  Group1 = rownames(pw_mat),
  Group2 = colnames(pw_mat),
  stringsAsFactors = FALSE
) %>%
  mutate(p_holm = mapply(function(r, c) pw_mat[r, c], Group1, Group2)) %>%
  filter(!is.na(p_holm)) %>%
  mutate(
    significant = ifelse(p_holm < 0.05, "Yes", "No"),
    signif_label = case_when(
      p_holm < 0.001 ~ "***",
      p_holm < 0.01  ~ "**",
      p_holm < 0.05  ~ "*",
      TRUE           ~ "ns"
    )
  ) %>%
  arrange(p_holm)

write.csv(pw_long,
          file.path(table_dir, "pairwise_wilcoxon_UA_full.csv"),
          row.names = FALSE)
cat("Full pairwise Wilcoxon matrix saved.\n")
