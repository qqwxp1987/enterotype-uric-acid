# ==============================================================================
# Supplementary Figure S3: Species–Metabolite Partial Spearman Heatmap
#   (full-species sensitivity analysis for Figure 5B)
#
# All 75 core species (Section 3.4) × 12 curated metabolites + serum UA,
# with sex and age partialled out via rank-residual partial correlation.
#
# Partial Spearman procedure:
#   1. Rank-transform both species abundance and metabolite/UA
#   2. Regress ranks on sex + age, extract residuals
#   3. Pearson-correlate residuals (≈ partial Spearman rho)
# ==============================================================================

library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(here)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
metab_raw <- readRDS(here("data", "clean",
                           "metabolites_normalized_multiIS_and_metadata_20250904.rds"))
clinical <- read.csv(here("data", "clean", "Metadata_20240925.csv"),
                      stringsAsFactors = FALSE)

meta_ps <- as(sample_data(ps), "data.frame")

# Color palettes
enterotype_colors <- c(
  "Type 1" = "#E64B35", "Type 2" = "#4DBBD5",
  "Type 3" = "#00A087", "Type 4" = "#3C5488"
)

class_colors <- c(
  "AAA catabolite" = "#D2691E",
  "Bile acid"      = "#2E8B57",
  "TMAO pathway"   = "#7B68EE",
  "Purine"         = "#DAA520",
  "Clinical"       = "#333333"
)

# ==============================================================================
# --- 2. Load ALL Core Species ---
# ==============================================================================

e1_core <- read.csv(here("results", "tables", "section_3_4", "taxa_E1_core.csv"),
                     stringsAsFactors = FALSE)
e4_core <- read.csv(here("results", "tables", "section_3_4", "taxa_E4_core.csv"),
                     stringsAsFactors = FALSE)
shared  <- read.csv(here("results", "tables", "section_3_4", "taxa_Shared.csv"),
                     stringsAsFactors = FALSE)

shared_features <- shared$feature

# E1 core: compute mean log2FC and direction
e1_core <- e1_core %>%
  mutate(
    mean_lfc_e1 = rowMeans(cbind(log2FoldChange_E1vE2, log2FoldChange_E1vE3), na.rm = TRUE),
    direction_e1 = ifelse(mean_lfc_e1 > 0, "Higher", "Lower"),
    is_shared = feature %in% shared_features
  )

# E4 core: compute mean log2FC
e4_core <- e4_core %>%
  mutate(
    mean_lfc_e4 = rowMeans(cbind(log2FoldChange_E4vE2, log2FoldChange_E4vE3), na.rm = TRUE),
    is_shared = feature %in% shared_features
  )

# Build combined species list (all unique)
e1_only_df <- e1_core %>%
  filter(!is_shared) %>%
  mutate(group = paste0("Type 1 \u2014 ", tolower(direction_e1)),
         sort_lfc = abs(mean_lfc_e1)) %>%
  select(feature, group, sort_lfc)

e4_only_df <- e4_core %>%
  filter(!is_shared) %>%
  mutate(group = "Type 4",
         sort_lfc = abs(mean_lfc_e4)) %>%
  select(feature, group, sort_lfc)

shared_df <- shared %>%
  mutate(group = "Shared", sort_lfc = 0) %>%
  select(feature, group, sort_lfc)

all_species <- bind_rows(e1_only_df, shared_df, e4_only_df) %>%
  distinct(feature, .keep_all = TRUE)

cat(sprintf("Total core species: %d\n", nrow(all_species)))
cat("  By group:\n")
print(table(all_species$group))

# Clean species names
clean_species <- function(x) {
  x <- sub(".*s__", "", x)
  x <- gsub("_", " ", x)
  x
}
all_species$display_name <- clean_species(all_species$feature)

# ==============================================================================
# --- 3. Curated Metabolites ---
# ==============================================================================

curated_metab <- tribble(
  ~feature_id,                            ~display_name,            ~class,
  "NEG0024_p-Cresol",                     "p-Cresol",               "AAA catabolite",
  "NEG0216_p-Cresol sulfate",             "p-Cresol sulfate",       "AAA catabolite",
  "NEG0046_Indoxyl sulfate",              "Indoxyl sulfate",        "AAA catabolite",
  "NEG0144_5-Methoxyindoleacetate",       "5-Methoxyindoleacetate", "AAA catabolite",
  "NEG0147_Indole-3-propionic acid",      "IPA",                    "AAA catabolite",
  "NEG0145_Phenylacetylglutamine",        "PAGln",                  "AAA catabolite",
  "NEG0015_Deoxycholic acid",             "DCA",                    "Bile acid",
  "NEG0082_Taurocholic acid",             "TCA",                    "Bile acid",
  "POS0156_Trimethylamine N-oxide",       "TMAO",                   "TMAO pathway",
  "POS0006_Betaine",                      "Betaine",                "TMAO pathway",
  "POS0002_Adenosine",                    "Adenosine",              "Purine",
  "POS0003_Adenine",                      "Adenine",                "Purine",
)

# ==============================================================================
# --- 4. Match Samples & Extract Matrices ---
# ==============================================================================

micro_ids <- sample_names(ps)
metab_ids <- metab_raw$metadata$SampleID

# Find clinical SampleID column
sid_col <- intersect(c("SampleID", "sampleid", "Sample_ID"), names(clinical))
if (length(sid_col) == 0) sid_col <- names(clinical)[1]
sid_col <- sid_col[1]
clinical_ids <- clinical[[sid_col]]

shared_ids <- Reduce(intersect, list(micro_ids, metab_ids, clinical_ids))
cat(sprintf("Matched samples: %d\n", length(shared_ids)))

# Species relative abundance
ps_sub <- prune_samples(shared_ids, ps)
ps_rel <- transform_sample_counts(ps_sub, function(x) x / sum(x))
otu_full <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) otu_full <- t(otu_full)

# Subset to core species (keep those present in OTU table)
avail_sp <- all_species$feature[all_species$feature %in% rownames(otu_full)]
sp_mat <- otu_full[avail_sp, shared_ids]
all_species <- all_species %>% filter(feature %in% avail_sp)
cat(sprintf("Species matrix: %d species x %d samples\n", nrow(sp_mat), ncol(sp_mat)))

# Metabolite matrix
metab_abund <- metab_raw$abundance
avail_metab <- curated_metab$feature_id[curated_metab$feature_id %in% rownames(metab_abund)]
met_mat <- metab_abund[avail_metab, shared_ids]

# UA
clin_sub <- clinical[match(shared_ids, clinical[[sid_col]]), ]
ua_vec <- as.numeric(clin_sub$UA)
names(ua_vec) <- shared_ids

# Combined target matrix
target_mat <- rbind(met_mat, UA = ua_vec[shared_ids])
cat(sprintf("Target matrix: %d variables x %d samples\n",
            nrow(target_mat), ncol(target_mat)))

# Covariates for partial correlation
sex_numeric <- as.numeric(factor(meta_ps[shared_ids, "Gender"]))
age_numeric <- as.numeric(meta_ps[shared_ids, "Age"])
covar_mat <- cbind(sex_numeric, age_numeric)

# ==============================================================================
# --- 5. Partial Spearman Correlation (rank-residual method) ---
# ==============================================================================

cat("\nComputing partial Spearman correlations (adjusting for sex + age)...\n")

n_sp  <- nrow(sp_mat)
n_tgt <- nrow(target_mat)
n_samples <- length(shared_ids)

rho_mat <- matrix(NA, n_sp, n_tgt,
                  dimnames = list(rownames(sp_mat), rownames(target_mat)))
p_mat   <- matrix(NA, n_sp, n_tgt,
                  dimnames = list(rownames(sp_mat), rownames(target_mat)))

# Pre-compute rank-residuals for all targets (done once)
target_resid <- matrix(NA, n_samples, n_tgt)
for (j in seq_len(n_tgt)) {
  y <- as.numeric(target_mat[j, ])
  valid <- !is.na(y)
  ry <- rep(NA, n_samples)
  ry[valid] <- rank(y[valid], ties.method = "average")
  fit <- lm(ry ~ covar_mat)
  target_resid[, j] <- residuals(fit)
}

# Pre-compute rank-residuals for all species
sp_resid <- matrix(NA, n_samples, n_sp)
for (i in seq_len(n_sp)) {
  x <- as.numeric(sp_mat[i, ])
  rx <- rank(x, ties.method = "average")
  fit <- lm(rx ~ covar_mat)
  sp_resid[, i] <- residuals(fit)
}

# Compute pairwise correlations on residuals
df_resid <- n_samples - 2 - ncol(covar_mat)  # df for partial correlation
for (i in seq_len(n_sp)) {
  for (j in seq_len(n_tgt)) {
    rx <- sp_resid[, i]
    ry <- target_resid[, j]
    valid <- !is.na(rx) & !is.na(ry)
    if (sum(valid) < 10) next
    r <- cor(rx[valid], ry[valid], method = "pearson")
    rho_mat[i, j] <- r
    # t-test for significance
    t_stat <- r * sqrt(sum(valid) - 2 - ncol(covar_mat)) /
              sqrt(1 - r^2)
    p_mat[i, j] <- 2 * pt(abs(t_stat),
                           df = sum(valid) - 2 - ncol(covar_mat),
                           lower.tail = FALSE)
  }
}

# BH correction
p_adj <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                nrow = n_sp, ncol = n_tgt,
                dimnames = dimnames(p_mat))

n_sig <- sum(p_adj < 0.05, na.rm = TRUE)
n_total <- sum(!is.na(p_adj))
cat(sprintf("Significant partial correlations (BH < 0.05): %d / %d (%.1f%%)\n",
            n_sig, n_total, n_sig / n_total * 100))

# ==============================================================================
# --- 6. Prepare Heatmap Annotations ---
# ==============================================================================

# Ensure row order matches rho_mat
row_info <- all_species[match(rownames(rho_mat), all_species$feature), ]

row_groups <- factor(
  row_info$group,
  levels = c("Type 1 \u2014 higher", "Type 1 \u2014 lower", "Shared", "Type 4")
)

row_labels <- row_info$display_name

group_colors <- c(
  "Type 1 \u2014 higher" = "#E64B35",
  "Type 1 \u2014 lower"  = "#F4A582",
  "Shared"               = "#999999",
  "Type 4"               = "#3C5488"
)

ha_row <- rowAnnotation(
  Core = row_groups,
  col = list(Core = group_colors),
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(4, "mm"),
  show_legend = TRUE
)

# Column annotations
col_info <- bind_rows(
  curated_metab %>% filter(feature_id %in% rownames(target_mat)) %>%
    select(feature_id, display_name, class),
  tibble(feature_id = "UA", display_name = "Serum UA", class = "Clinical")
)
col_info <- col_info[match(rownames(target_mat), col_info$feature_id), ]

col_groups <- factor(col_info$class,
                     levels = c("AAA catabolite", "Bile acid",
                                "TMAO pathway", "Purine", "Clinical"))
col_labels <- col_info$display_name

ha_col <- HeatmapAnnotation(
  Class = col_groups,
  col = list(Class = class_colors),
  annotation_name_gp = gpar(fontsize = 8),
  simple_anno_size = unit(4, "mm"),
  show_legend = TRUE
)

# ==============================================================================
# --- 7. Build Heatmap ---
# ==============================================================================

max_rho <- 0.5
col_fun <- colorRamp2(
  c(-max_rho, -max_rho / 2, 0, max_rho / 2, max_rho),
  c("#3C5488", "#8FAABE", "white", "#E8927C", "#E64B35")
)

cell_fn <- function(j, i, x, y, width, height, fill) {
  val <- p_adj[i, j]
  if (!is.na(val)) {
    if (val < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 6, col = "black"))
    } else if (val < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 7, col = "black"))
    }
  }
}

# Adaptive height: ~4mm per species
ht_height <- unit(max(12, nrow(rho_mat) * 0.4), "cm")

set.seed(42)
ht <- Heatmap(
  rho_mat,
  name = "Partial \u03c1",
  col = col_fun,

  # Rows
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = 6.5, fontface = "italic"),
  row_split = row_groups,
  row_title_gp = gpar(fontsize = 8, fontface = "bold"),
  row_gap = unit(2, "mm"),
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  row_dend_width = unit(12, "mm"),

  # Columns
  column_labels = col_labels,
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 45,
  column_split = col_groups,
  column_title_gp = gpar(fontsize = 8, fontface = "bold"),
  column_gap = unit(2, "mm"),
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,

  # Cells
  cell_fun = cell_fn,
  rect_gp = gpar(col = "grey90", lwd = 0.3),

  # Annotations
  left_annotation = ha_row,
  top_annotation = ha_col,

  # Legend
  heatmap_legend_param = list(
    title = "Partial \u03c1",
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 7),
    legend_height = unit(3, "cm"),
    at = c(-0.5, -0.25, 0, 0.25, 0.5)
  ),

  # Size
  width = unit(9, "cm"),
  height = ht_height
)

# ==============================================================================
# --- 8. Save ---
# ==============================================================================

dir.create(here("supplementary", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables", "section_3_4"), recursive = TRUE, showWarnings = FALSE)

# Compute figure height in inches
fig_h <- max(14, nrow(rho_mat) * 0.18 + 3)

cairo_pdf(here("supplementary", "figures", "FigureS3_Full_Species_Heatmap.pdf"),
          width = 10, height = fig_h)
draw(ht, merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

png(here("supplementary", "figures", "FigureS3_Full_Species_Heatmap.png"),
    width = 10, height = fig_h, units = "in", res = 300, type = "cairo")
draw(ht, merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

# Save tables
cor_table <- as.data.frame(rho_mat) %>%
  tibble::rownames_to_column("species") %>%
  mutate(species_clean = clean_species(species))
write.csv(cor_table,
          here("results", "tables", "section_3_4",
               "species_metabolite_partial_spearman_rho.csv"),
          row.names = FALSE)

padj_table <- as.data.frame(p_adj) %>%
  tibble::rownames_to_column("species") %>%
  mutate(species_clean = clean_species(species))
write.csv(padj_table,
          here("results", "tables", "section_3_4",
               "species_metabolite_partial_spearman_padj.csv"),
          row.names = FALSE)

saveRDS(ht, here("results", "cache", "p4b_full_adj.rds"))

cat(sprintf("\n=== Outputs saved (full %d species, sex+age adjusted) ===\n",
            nrow(rho_mat)))
cat("  Heatmap: supplementary/figures/FigureS3_Full_Species_Heatmap.*\n")
cat("  Tables:  results/tables/section_3_4/species_metabolite_partial_spearman_*\n")
