# ==============================================================================
# Figure 5A: db-RDA — Metabolome-Constrained Microbiome Ordination
#
# Quantifies how much of the microbiome compositional variation is explained
# by the serum metabolome, after partialling out Sex and Age.
#
# Analyses:
#   1. Partial db-RDA (global): capscale(BC ~ metab_PCs | Sex + Age)
#   2. Envfit: project individual metabolites onto db-RDA axes
#   3. Triplot with curated, class-colored metabolite arrows
# ==============================================================================

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
metab_raw <- readRDS(here("data", "clean",
                           "metabolites_normalized_multiIS_and_metadata_20250904.rds"))

# Define enterotype factor
sample_data(ps)$Enterotype <- factor(
  sample_data(ps)$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)

# Consistent color palette
enterotype_colors <- c(
  "Type 1" = "#E64B35",
  "Type 2" = "#4DBBD5",
  "Type 3" = "#00A087",
  "Type 4" = "#3C5488"
)

# --- 2. Match Samples & Prepare Matrices ---
micro_ids <- sample_names(ps)
metab_ids <- metab_raw$metadata$SampleID
shared_ids <- intersect(micro_ids, metab_ids)
cat(sprintf("Shared samples: %d\n", length(shared_ids)))

ps_sub <- prune_samples(shared_ids, ps)
meta_df <- as(sample_data(ps_sub), "data.frame")

# Microbiome: Bray-Curtis distance
dist_bc <- distance(ps_sub, method = "bray")

# Metabolome: log2 transform -> PCA for dimensionality reduction
metab_abund <- metab_raw$abundance[, shared_ids]
min_nonzero <- min(metab_abund[metab_abund > 0], na.rm = TRUE)
metab_log <- log2(metab_abund + min_nonzero / 2)

# PCA on metabolome (samples as rows)
metab_pca <- prcomp(t(metab_log), center = TRUE, scale. = TRUE)
cum_var <- cumsum(metab_pca$sdev^2 / sum(metab_pca$sdev^2))
n_pcs <- which(cum_var >= 0.80)[1]
cat(sprintf("Metabolome PCA: %d PCs capture %.1f%% variance\n",
            n_pcs, cum_var[n_pcs] * 100))

metab_scores <- as.data.frame(metab_pca$x[shared_ids, 1:n_pcs])

# Covariates
cov_df <- data.frame(
  Sex = factor(meta_df[shared_ids, "Gender"]),
  Age = as.numeric(meta_df[shared_ids, "Age"]),
  Enterotype = meta_df[shared_ids, "Enterotype"],
  row.names = shared_ids
)

pred_df <- cbind(metab_scores, cov_df)

# ==============================================================================
# --- 3. Global Partial db-RDA ---
# ==============================================================================

pc_names <- colnames(metab_scores)
formula_full <- as.formula(
  paste("dist_bc ~", paste(pc_names, collapse = " + "), "+ Condition(Sex + Age)")
)

set.seed(42)
dbrda_res <- capscale(formula_full, data = pred_df)

# Permutation test
dbrda_anova <- anova(dbrda_res, permutations = 999)

# Variance partitioning
total_inertia <- dbrda_res$tot.chi
constrained_inertia <- dbrda_res$CCA$tot.chi
conditioned_inertia <- dbrda_res$pCCA$tot.chi
unconstrained_inertia <- dbrda_res$CA$tot.chi

R2_raw <- constrained_inertia / total_inertia
n_samp <- length(shared_ids)
R2_adj <- 1 - (1 - R2_raw) * (n_samp - 1) / (n_samp - n_pcs - 1)

cat("\n=== Global Partial db-RDA ===\n")
cat(sprintf("  Total inertia:         %.4f\n", total_inertia))
cat(sprintf("  Conditioned (Sex+Age): %.4f (%.1f%%)\n",
            conditioned_inertia, conditioned_inertia / total_inertia * 100))
cat(sprintf("  Constrained (Metab):   %.4f (%.1f%%)\n",
            constrained_inertia, constrained_inertia / total_inertia * 100))
cat(sprintf("  Unconstrained:         %.4f (%.1f%%)\n",
            unconstrained_inertia, unconstrained_inertia / total_inertia * 100))
cat(sprintf("  R2 (raw):  %.4f\n", R2_raw))
cat(sprintf("  R2 (adj):  %.4f\n", R2_adj))
cat(sprintf("  ANOVA P:   %.4f\n", dbrda_anova$`Pr(>F)`[1]))

# ==============================================================================
# --- 4. Envfit: Project Individual Metabolites onto db-RDA Axes ---
# ==============================================================================

cat("\n=== Envfit: Metabolites on db-RDA Axes ===\n")

dbrda_scores <- scores(dbrda_res, display = "sites", choices = 1:2)

set.seed(42)
ef <- envfit(dbrda_scores, t(metab_log[, shared_ids]),
             permutations = 999, choices = 1:2)

ef_df <- as.data.frame(scores(ef, display = "vectors"))
ef_df$feature <- rownames(ef_df)
ef_df$r2 <- ef$vectors$r
ef_df$pval <- ef$vectors$pvals

ef_sig <- ef_df %>%
  filter(pval < 0.05) %>%
  arrange(desc(r2))

cat(sprintf("  Significant metabolites (P < 0.05): %d / %d\n",
            nrow(ef_sig), nrow(ef_df)))

# ==============================================================================
# --- 5. Curated Metabolite Arrow Set ---
# ==============================================================================

# Define biologically meaningful metabolites grouped by functional class.
# These connect directly to Section 3.3 findings and the paper's narrative.
curated <- tribble(
  ~feature_id,                            ~display_name,                  ~class,
  # --- AAA catabolites (Section 3.3 Type 1 core) ---
  "NEG0216_p-Cresol sulfate",             "p-Cresol sulfate",             "AAA catabolite",
  "NEG0024_p-Cresol",                     "p-Cresol",                     "AAA catabolite",
  "NEG0046_Indoxyl sulfate",              "Indoxyl sulfate",              "AAA catabolite",
  "NEG0144_5-Methoxyindoleacetate",       "5-Methoxyindoleacetate",       "AAA catabolite",
  "NEG0147_Indole-3-propionic acid",      "Indole-3-propionic acid",      "AAA catabolite",
  "NEG0145_Phenylacetylglutamine",        "Phenylacetylglutamine",        "AAA catabolite",
  # --- Bile acids ---
  "NEG0015_Deoxycholic acid",             "Deoxycholic acid",             "Bile acid",
  "NEG0082_Taurocholic acid",             "Taurocholic acid",             "Bile acid",
  # --- TMAO / One-carbon metabolism ---
  "POS0156_Trimethylamine N-oxide",       "TMAO",                         "TMAO pathway",
  "POS0006_Betaine",                      "Betaine",                      "TMAO pathway",
  # --- Purine metabolism ---
  "POS0002_Adenosine",                    "Adenosine",                    "Purine",
  "POS0003_Adenine",                      "Adenine",                      "Purine",
)

# Match curated list to envfit results (keep only significant ones)
arrow_data <- curated %>%
  inner_join(ef_sig, by = c("feature_id" = "feature")) %>%
  rename(feature = feature_id)

cat(sprintf("\nCurated arrows: %d / %d matched to significant envfit results\n",
            nrow(arrow_data), nrow(curated)))

# Report any curated metabolites that were NOT significant
missing <- curated %>%
  filter(!feature_id %in% ef_sig$feature)
if (nrow(missing) > 0) {
  cat("  Not significant (excluded from arrows):\n")
  print(missing$display_name)
}

# Class color palette (distinct from enterotype colors)
class_colors <- c(
  "AAA catabolite" = "#D2691E",    # chocolate brown
  "Bile acid"      = "#2E8B57",    # sea green
  "TMAO pathway"   = "#7B68EE",    # medium slate blue
  "Purine"         = "#DAA520"     # goldenrod
)

# ==============================================================================
# --- 6. db-RDA Triplot ---
# ==============================================================================

# Axis labels
eig <- eigenvals(dbrda_res)
cap_eig <- eig[grep("^CAP", names(eig))]
pct_cap1 <- cap_eig[1] / total_inertia * 100
pct_cap2 <- cap_eig[2] / total_inertia * 100

# Sample scores
site_df <- as.data.frame(dbrda_scores)
colnames(site_df) <- c("CAP1", "CAP2")
site_df$Enterotype <- cov_df[rownames(site_df), "Enterotype"]

# Scale arrows to fit within sample point range
arrow_scale <- max(abs(site_df[, 1:2])) / max(abs(arrow_data[, c("CAP1", "CAP2")])) * 0.85
arrow_df <- arrow_data %>%
  mutate(
    x = CAP1 * arrow_scale,
    y = CAP2 * arrow_scale,
    class = factor(class, levels = names(class_colors))
  )

# Annotation
p_label <- ifelse(dbrda_anova$`Pr(>F)`[1] < 0.001,
                  "P < 0.001",
                  sprintf("P = %.3f", dbrda_anova$`Pr(>F)`[1]))
annot_text <- sprintf("Partial db-RDA\nR\u00b2(adj) = %.3f, %s",
                      R2_adj, p_label)

p4a <- ggplot() +
  # Sample points
  geom_point(
    data = site_df,
    aes(x = CAP1, y = CAP2, color = Enterotype),
    size = 1.8, alpha = 0.55, shape = 16
  ) +
  # 95% confidence ellipses
  stat_ellipse(
    data = site_df,
    aes(x = CAP1, y = CAP2, color = Enterotype),
    level = 0.95, linewidth = 0.5, linetype = "dashed"
  ) +
  # Metabolite arrows (colored by class)
  geom_segment(
    data = arrow_df,
    aes(x = 0, y = 0, xend = x, yend = y, color = class),
    arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
    linewidth = 0.5, alpha = 0.85,
    show.legend = TRUE
  ) +
  # Metabolite labels
  geom_text_repel(
    data = arrow_df,
    aes(x = x, y = y, label = display_name, color = class),
    size = 2.5, fontface = "bold",
    max.overlaps = Inf,
    force = 4, force_pull = 0.5,
    segment.size = 0.2, segment.color = "grey60",
    box.padding = 0.7, point.padding = 0.3,
    min.segment.length = 0,
    seed = 42,
    show.legend = FALSE
  ) +
  # Two separate color scales: enterotype for points, class for arrows
  # Use ggnewscale to layer them, or manually combine into one scale.
  # Simpler approach: combine all colors into a single scale with grouped labels.
  scale_color_manual(
    values = c(enterotype_colors, class_colors),
    breaks = c(names(enterotype_colors), names(class_colors)),
    labels = c(names(enterotype_colors), names(class_colors)),
    name = NULL,
    guide = guide_legend(
      override.aes = list(
        shape = c(16, 16, 16, 16, NA, NA, NA, NA),
        linetype = c(0, 0, 0, 0, 1, 1, 1, 1),
        linewidth = c(NA, NA, NA, NA, 0.8, 0.8, 0.8, 0.8),
        size = c(2.5, 2.5, 2.5, 2.5, NA, NA, NA, NA),
        alpha = c(0.7, 0.7, 0.7, 0.7, 0.85, 0.85, 0.85, 0.85)
      )
    )
  ) +
  annotate(
    "label",
    x = -Inf, y = -Inf,
    hjust = -0.05, vjust = -0.2,
    label = annot_text,
    size = 3.0, fontface = "italic", color = "grey30",
    fill = "white", label.size = 0, alpha = 0.85
  ) +
  labs(
    x = sprintf("CAP1 (%.1f%%)", pct_cap1),
    y = sprintf("CAP2 (%.1f%%)", pct_cap2)
  ) +
  coord_fixed(ratio = 1) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    legend.text = element_text(size = 8),
    legend.key.size = unit(0.4, "cm"),
    legend.position = "right",
    legend.spacing.y = unit(0.15, "cm"),
    plot.margin = margin(5, 5, 5, 5)
  )

# ==============================================================================
# --- 7. Save ---
# ==============================================================================

dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables", "section_3_4"), recursive = TRUE, showWarnings = FALSE)

ggsave(here("results", "figures", "panels", "Figure5A_dbRDA_triplot.pdf"),
       p4a, width = 7.5, height = 5.5, dpi = 300)
ggsave(here("results", "figures", "panels", "Figure5A_dbRDA_triplot.png"),
       p4a, width = 7.5, height = 5.5, dpi = 300)

write.csv(ef_sig,
          here("results", "tables", "section_3_4", "envfit_significant_metabolites.csv"),
          row.names = FALSE)

# Save curated arrow data for reference
write.csv(arrow_df %>% select(feature, display_name, class, r2, pval),
          here("results", "tables", "section_3_4", "envfit_curated_arrows.csv"),
          row.names = FALSE)

saveRDS(p4a, here("results", "cache", "p4a.rds"))
saveRDS(dbrda_res, here("results", "cache", "dbrda_result.rds"))

cat("\n=== Outputs saved ===\n")
cat("  Triplot: results/figures/panels/Figure5A_dbRDA_triplot.*\n")
cat("  Tables:  results/tables/section_3_4/\n")
