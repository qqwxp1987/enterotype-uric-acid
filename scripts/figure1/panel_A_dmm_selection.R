# ==============================================================================
# Figure 1A: DMM Enterotype Analysis & Model Selection Plot
# 1) Fits DMM models (k=1~10), selects optimal k by Laplace approximation
# 2) Assigns enterotype labels to all samples (highest posterior probability)
# 3) Saves labeled phyloseq for all downstream panels
# 4) Generates the Laplace vs k selection plot (Figure 1A)
#
# This script is the SOLE source of enterotype assignments for the entire
# analysis pipeline. All other panel scripts depend on its output.
# ==============================================================================

library(phyloseq)
library(microbiome)
library(DirichletMultinomial)
library(ggplot2)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "UA_Analysis_20230811.rds"))

# Aggregate to genus level (microbiome::aggregate_taxa, NOT tax_glom)
ps_genus <- microbiome::aggregate_taxa(ps, level = "Genus")

# Extract count matrix: samples (rows) x genera (cols)
otu_mat <- t(as(otu_table(ps_genus), "matrix"))

# microbiome::aggregate_taxa() returns double; newer DirichletMultinomial
# versions strictly require integer. Convert explicitly.
otu_mat <- matrix(
  as.integer(round(otu_mat)),
  nrow = nrow(otu_mat), dimnames = dimnames(otu_mat)
)

cat("OTU matrix dimensions:", dim(otu_mat), "\n")
cat("  Samples:", nrow(otu_mat), "  Genera:", ncol(otu_mat), "\n")

# --- 2. Fit DMM Models (k = 1 to 10) ---
# DMM fitting is stochastic; a fixed seed is used for reproducibility.
DMM_SEED <- 10L
k_range <- 1:10

cat("Fitting DMM models for k =", paste(k_range, collapse = ", "), "...\n")
cat("Using fixed EM seed =", DMM_SEED, "\n")

dmm_fits <- lapply(k_range, function(k) {
  cat("  Fitting k =", k, "...")
  fit <- dmn(otu_mat, k = k, seed = DMM_SEED, verbose = FALSE)
  cat(" done (Laplace =", round(laplace(fit), 1), ")\n")
  fit
})

# --- 3. Extract Laplace Values & Identify Optimal k ---
laplace_vals <- sapply(dmm_fits, laplace)

dmm_df <- data.frame(
  k = k_range,
  Laplace = laplace_vals
)

optimal_k <- dmm_df$k[which.min(dmm_df$Laplace)]
cat("\nOptimal number of enterotypes: k =", optimal_k, "\n")

# --- 4. Assign Enterotype Labels ---
# Extract posterior probabilities from the optimal model
best_fit <- dmm_fits[[optimal_k]]
post_prob <- mixture(best_fit, assign = FALSE)  # samples x components
raw_assignments <- apply(post_prob, 1, which.max)
names(raw_assignments) <- rownames(otu_mat)

# DMM component numbering is arbitrary. We adopt a standard convention:
# label components by descending sample count (largest cluster = Type 1).
comp_counts <- table(raw_assignments)
size_order <- order(comp_counts, decreasing = TRUE)

# Build mapping: comp_to_type[original_component] = new_type_number
comp_to_type <- integer(optimal_k)
comp_to_type[size_order] <- seq_len(optimal_k)

# Apply mapping to get final enterotype labels
enterotype_labels <- comp_to_type[raw_assignments]
names(enterotype_labels) <- names(raw_assignments)

cat("\nEnterotype distribution (labeled by descending size):\n")
print(table(enterotype_labels))

# --- 5. Save Labeled Phyloseq ---
# Add enterotype assignments to the phyloseq sample_data
sample_data(ps)$Dirichlettype_Genus <- as.character(
  enterotype_labels[sample_names(ps)]
)

saveRDS(ps, here("data", "clean",
  "phyloseq_Taxa_with_enterotype_seed10.rds"))
cat("Labeled phyloseq saved to",
    "data/clean/phyloseq_Taxa_with_enterotype_seed10.rds\n")

# --- 6. Save Cache for Panel F ---
dir.create(here("results", "cache"), recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    dmm_fits = dmm_fits,
    dmm_df = dmm_df,
    optimal_k = optimal_k,
    component_mapping = comp_to_type  # mapping[raw_comp] = Type number
  ),
  here("results", "cache", "dmm_model_fits.rds")
)

# --- 7. Plot ---

source(here("scripts", "utils", "palettes.R"))
enterotype_colors <- ENTEROTYPE_COLORS

p1a <- ggplot(dmm_df, aes(x = k, y = Laplace)) +
  geom_line(color = "grey40", linewidth = 0.8) +
  geom_point(
    aes(fill = ifelse(k == optimal_k, "optimal", "other")),
    shape = 21, size = 3.5, stroke = 0.6, color = "grey30"
  ) +
  scale_fill_manual(
    values = c("optimal" = "#E64B35", "other" = "grey70"),
    guide = "none"
  ) +
  annotate(
    "segment",
    x = optimal_k, xend = optimal_k,
    y = min(dmm_df$Laplace) - diff(range(dmm_df$Laplace)) * 0.08,
    yend = min(dmm_df$Laplace),
    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
    color = "#E64B35", linewidth = 0.6
  ) +
  annotate(
    "text",
    x = optimal_k,
    y = min(dmm_df$Laplace) - diff(range(dmm_df$Laplace)) * 0.12,
    label = paste0("k = ", optimal_k),
    fontface = "bold", size = 3.8, color = "#E64B35"
  ) +
  scale_x_continuous(breaks = k_range) +
  labs(
    x = "Number of Dirichlet components (k)",
    y = "Laplace approximation"
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9, color = "black"),
    plot.margin = margin(5, 10, 5, 5)
  )

# --- 8. Save Figure ---
dir.create(here("results", "figures", "panels"),
  recursive = TRUE, showWarnings = FALSE)

ggsave(
  here("results", "figures", "panels", "Figure1A_DMM_Selection.pdf"),
  p1a, width = 4.5, height = 3.5, dpi = 300
)
ggsave(
  here("results", "figures", "panels", "Figure1A_DMM_Selection.png"),
  p1a, width = 4.5, height = 3.5, dpi = 300
)

cat("\nFigure 1A saved to results/figures/panels/\n")

saveRDS(p1a, here("results", "cache", "p1a.rds"))
