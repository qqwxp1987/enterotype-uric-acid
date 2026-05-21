# ==============================================================================
# Figure 5B (main): Top Core Species × Curated Metabolites
#                    Partial Spearman Heatmap (sex + age adjusted)
#
# Top 10 "higher" + top 10 "lower" from Type 1 core (by |mean log2FC|)
# + all Type 4 core species (3, incl. shared V. dispar)
# ≈ 23 species × 13 targets (12 metabolites + serum UA)
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

class_colors <- c(
  "AAA catabolite" = "#D2691E",
  "Bile acid"      = "#2E8B57",
  "TMAO pathway"   = "#7B68EE",
  "Purine"         = "#DAA520",
  "Clinical"       = "#333333"
)

# ==============================================================================
# --- 2. Select Top Core Species ---
# ==============================================================================

e1_core <- read.csv(here("results", "tables", "section_3_4", "taxa_E1_core.csv"),
                     stringsAsFactors = FALSE)
e4_core <- read.csv(here("results", "tables", "section_3_4", "taxa_E4_core.csv"),
                     stringsAsFactors = FALSE)
shared  <- read.csv(here("results", "tables", "section_3_4", "taxa_Shared.csv"),
                     stringsAsFactors = FALSE)

shared_features <- shared$feature

e1_core <- e1_core %>%
  mutate(
    mean_lfc = rowMeans(cbind(log2FoldChange_E1vE2, log2FoldChange_E1vE3), na.rm = TRUE),
    abs_lfc = abs(mean_lfc),
    direction = ifelse(mean_lfc > 0, "Higher", "Lower")
  )

e1_only <- e1_core %>% filter(!feature %in% shared_features)
e1_higher <- e1_only %>% filter(direction == "Higher") %>%
  arrange(desc(abs_lfc)) %>% slice_head(n = 10)
e1_lower <- e1_only %>% filter(direction == "Lower") %>%
  arrange(desc(abs_lfc)) %>% slice_head(n = 10)

e4_all <- e4_core %>%
  mutate(
    mean_lfc = rowMeans(cbind(log2FoldChange_E4vE2, log2FoldChange_E4vE3), na.rm = TRUE)
  )

selected_species <- bind_rows(
  e1_higher %>% mutate(group = "Type 1 \u2014 higher"),
  e1_lower  %>% mutate(group = "Type 1 \u2014 lower"),
  e4_all    %>% mutate(group = ifelse(feature %in% shared_features, "Shared", "Type 4"))
) %>%
  select(feature, group)

clean_species <- function(x) {
  x <- sub(".*s__", "", x)
  x <- gsub("_", " ", x)
  x
}
selected_species$display_name <- clean_species(selected_species$feature)

cat(sprintf("Selected species: %d\n", nrow(selected_species)))
print(table(selected_species$group))

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

avail_sp <- selected_species$feature[selected_species$feature %in% rownames(otu_full)]
sp_mat <- otu_full[avail_sp, shared_ids]
selected_species <- selected_species %>% filter(feature %in% avail_sp)

# Metabolite matrix
metab_abund <- metab_raw$abundance
avail_metab <- curated_metab$feature_id[curated_metab$feature_id %in% rownames(metab_abund)]
met_mat <- metab_abund[avail_metab, shared_ids]

# UA
clin_sub <- clinical[match(shared_ids, clinical[[sid_col]]), ]
ua_vec <- as.numeric(clin_sub$UA)
names(ua_vec) <- shared_ids

target_mat <- rbind(met_mat, UA = ua_vec[shared_ids])

# Covariates
sex_numeric <- as.numeric(factor(meta_ps[shared_ids, "Gender"]))
age_numeric <- as.numeric(meta_ps[shared_ids, "Age"])
covar_mat <- cbind(sex_numeric, age_numeric)

# ==============================================================================
# --- 5. Partial Spearman Correlation ---
# ==============================================================================

cat("\nComputing partial Spearman correlations (sex + age adjusted)...\n")

n_sp  <- nrow(sp_mat)
n_tgt <- nrow(target_mat)
n_samples <- length(shared_ids)

rho_mat <- matrix(NA, n_sp, n_tgt,
                  dimnames = list(rownames(sp_mat), rownames(target_mat)))
p_mat   <- matrix(NA, n_sp, n_tgt,
                  dimnames = list(rownames(sp_mat), rownames(target_mat)))

# Pre-compute rank-residuals for targets
target_resid <- matrix(NA, n_samples, n_tgt)
for (j in seq_len(n_tgt)) {
  y <- as.numeric(target_mat[j, ])
  valid <- !is.na(y)
  ry <- rep(NA, n_samples)
  ry[valid] <- rank(y[valid], ties.method = "average")
  fit <- lm(ry ~ covar_mat)
  target_resid[, j] <- residuals(fit)
}

# Pre-compute rank-residuals for species
sp_resid <- matrix(NA, n_samples, n_sp)
for (i in seq_len(n_sp)) {
  x <- as.numeric(sp_mat[i, ])
  rx <- rank(x, ties.method = "average")
  fit <- lm(rx ~ covar_mat)
  sp_resid[, i] <- residuals(fit)
}

# Pairwise partial correlations
for (i in seq_len(n_sp)) {
  for (j in seq_len(n_tgt)) {
    rx <- sp_resid[, i]
    ry <- target_resid[, j]
    valid <- !is.na(rx) & !is.na(ry)
    if (sum(valid) < 10) next
    r <- cor(rx[valid], ry[valid], method = "pearson")
    rho_mat[i, j] <- r
    t_stat <- r * sqrt(sum(valid) - 2 - ncol(covar_mat)) /
              sqrt(1 - r^2)
    p_mat[i, j] <- 2 * pt(abs(t_stat),
                           df = sum(valid) - 2 - ncol(covar_mat),
                           lower.tail = FALSE)
  }
}

p_adj <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                nrow = n_sp, ncol = n_tgt,
                dimnames = dimnames(p_mat))

n_sig <- sum(p_adj < 0.05, na.rm = TRUE)
cat(sprintf("Significant (BH < 0.05): %d / %d\n", n_sig, sum(!is.na(p_adj))))

# ==============================================================================
# --- 6. Heatmap ---
# ==============================================================================

row_info <- selected_species[match(rownames(rho_mat), selected_species$feature), ]

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

max_rho <- 0.5
col_fun <- colorRamp2(
  c(-max_rho, -max_rho / 2, 0, max_rho / 2, max_rho),
  c("#3C5488", "#8FAABE", "white", "#E8927C", "#E64B35")
)

cell_fn <- function(j, i, x, y, width, height, fill) {
  val <- p_adj[i, j]
  if (!is.na(val)) {
    if (val < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = 7, col = "black"))
    } else if (val < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = 8, col = "black"))
    }
  }
}

set.seed(42)
ht <- Heatmap(
  rho_mat,
  name = "Partial \u03c1",
  col = col_fun,

  row_labels = row_labels,
  row_names_gp = gpar(fontsize = 8, fontface = "italic"),
  row_split = row_groups,
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  row_gap = unit(2, "mm"),
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  row_dend_width = unit(15, "mm"),

  column_labels = col_labels,
  column_names_gp = gpar(fontsize = 8.5),
  column_names_rot = 45,
  column_split = col_groups,
  column_title_gp = gpar(fontsize = 7.5, fontface = "bold"),
  column_title_rot = 30,
  column_gap = unit(4, "mm"),
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,

  cell_fun = cell_fn,
  rect_gp = gpar(col = "grey90", lwd = 0.5),

  left_annotation = ha_row,
  top_annotation = ha_col,

  heatmap_legend_param = list(
    title = "Partial \u03c1",
    title_gp = gpar(fontsize = 9, fontface = "bold"),
    labels_gp = gpar(fontsize = 8),
    legend_height = unit(3, "cm"),
    at = c(-0.5, -0.25, 0, 0.25, 0.5)
  ),

  width = unit(9, "cm"),
  height = unit(14, "cm")
)

# ==============================================================================
# --- 7. Save ---
# ==============================================================================

dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables", "section_3_4"), recursive = TRUE, showWarnings = FALSE)

pdf(here("results", "figures", "panels", "Figure5B_top_adjusted.pdf"),
    width = 9, height = 10)
draw(ht, merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

png(here("results", "figures", "panels", "Figure5B_top_adjusted.png"),
    width = 9, height = 10, units = "in", res = 300)
draw(ht, merge_legend = TRUE,
     heatmap_legend_side = "right",
     annotation_legend_side = "right")
dev.off()

# Tables
cor_table <- as.data.frame(rho_mat) %>%
  tibble::rownames_to_column("species") %>%
  mutate(species_clean = clean_species(species))
write.csv(cor_table,
          here("results", "tables", "section_3_4",
               "top_species_metabolite_partial_rho.csv"),
          row.names = FALSE)

padj_table <- as.data.frame(p_adj) %>%
  tibble::rownames_to_column("species") %>%
  mutate(species_clean = clean_species(species))
write.csv(padj_table,
          here("results", "tables", "section_3_4",
               "top_species_metabolite_partial_padj.csv"),
          row.names = FALSE)

saveRDS(ht, here("results", "cache", "p4b_top_adj.rds"))

cat("\n=== Main Panel B saved (top species, sex+age adjusted) ===\n")
