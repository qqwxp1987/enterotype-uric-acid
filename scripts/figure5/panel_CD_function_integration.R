# ==============================================================================
# Figure 5 Panels C & D: Functional Layer Cross-Omics Integration
#
# Complements Panel A (db-RDA) and Panel B (species × metabolites heatmap)
# by adding the functional dimension (pathways + convergent ECs).
#
# Panel C: [Core Pathways + Convergent ECs] × [Curated Metabolites + UA]
#          Partial Spearman heatmap (sex + age adjusted)
# Panel D: Forest plot comparing feature-UA partial correlations
#          across three omics layers (species, pathways, convergent ECs)
#          — with and without enterotype adjustment
#
# Feature selection: Based solely on Section 3.3 findings.
# ==============================================================================

library(phyloseq)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

source(here("scripts", "figure5", "figure5_style.R"))

cat("========================================\n")
cat("Section 3.4: Functional Integration\n")
cat("========================================\n\n")

# --- Output directories ---
dir.create(here("results", "figures", "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "tables", "section_3_4"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("results", "cache"), recursive = TRUE, showWarnings = FALSE)

# --- Consistent palettes ---
enterotype_colors <- c(
  "Type 1" = "#E64B35", "Type 2" = "#4DBBD5",
  "Type 3" = "#00A087", "Type 4" = "#3C5488"
)
class_colors <- c(
  "AAA catabolite" = "#D2691E", "Bile acid" = "#2E8B57",
  "TMAO pathway"   = "#7B68EE", "Purine"   = "#DAA520",
  "Clinical"       = "#333333"
)

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================

cat("--- Loading data ---\n")

# Microbiome: species (for forest plot), pathway, EC
ps_taxa <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
ps_path <- readRDS(here("data", "clean", "phyloseq_Pathway_with_enterotype_20250210.rds"))
ps_ec   <- readRDS(here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"))

# Update pathway/EC enterotype labels from authoritative seed10 map
meta_auth <- as(sample_data(ps_taxa), "data.frame")
et_map <- data.frame(
  SampleID   = rownames(meta_auth),
  Enterotype = meta_auth$Dirichlettype_Genus,
  Gender     = meta_auth$Gender,
  Age        = as.numeric(meta_auth$Age),
  stringsAsFactors = FALSE
)

update_enterotype <- function(ps, et_map) {
  meta <- as(sample_data(ps), "data.frame")
  shared <- intersect(rownames(meta), et_map$SampleID)
  meta[shared, "Dirichlettype_Genus"] <- et_map$Enterotype[match(shared, et_map$SampleID)]
  sample_data(ps) <- sample_data(meta)
  ps
}
ps_path <- update_enterotype(ps_path, et_map)
ps_ec   <- update_enterotype(ps_ec, et_map)

# Metabolome
metab_raw <- readRDS(here("data", "clean",
                           "metabolites_normalized_multiIS_and_metadata_20250904.rds"))

# Clinical (for UA)
clinical <- read.csv(here("data", "clean", "Metadata_20240925.csv"),
                      stringsAsFactors = FALSE)

# ==============================================================================
# 2. FEATURE SELECTION (based solely on Section 3.3)
# ==============================================================================

cat("\n--- Feature selection ---\n")

# ---- 2a. Core Pathways ----
# E1-only pathways: selected by functional themes discussed in Section 3.3
e1_pathways <- c(
  "ARGORNPROST-PWY: arginine, ornithine and proline interconversion",
  "PWY-5004: superpathway of L-citrulline metabolism",
  "PWY-7234: inosine-5-phosphate biosynthesis III",
  "PWY-6608: guanosine nucleotides degradation III",
  "PWY0-1296: purine ribonucleosides degradation",
  "PROPFERM-PWY: L-alanine fermentation to propanoate and acetate",
  "PWY-5690: TCA cycle II (plants and fungi)",
  "POLYAMINSYN3-PWY: superpathway of polyamine biosynthesis II"
)

# E4-only pathways: purine/nucleotide biosynthesis, amino acid, peptidoglycan, SAM
e4_pathways <- c(
  "PWY-5695: urate biosynthesis/inosine 5-phosphate degradation",
  "PWY-7219: adenosine ribonucleotides de novo biosynthesis",
  "PWY-6122: 5-aminoimidazole ribonucleotide biosynthesis II",
  "PWY-6609: adenine and adenosine salvage III",
  "PEPTIDOGLYCANSYN-PWY: peptidoglycan biosynthesis I (meso-diaminopimelate containing)",
  "BRANCHED-CHAIN-AA-SYN-PWY: superpathway of branched amino acid biosynthesis",
  "PWY-6151: S-adenosyl-L-methionine cycle I",
  "PWY-5686: UMP biosynthesis"
)

pathway_selection <- data.frame(
  feature = c(e1_pathways, e4_pathways),
  group = c(rep("E1-only pathway", length(e1_pathways)),
            rep("E4-only pathway", length(e4_pathways))),
  stringsAsFactors = FALSE
)

# Short display names for pathways
pathway_selection$display_name <- c(
  # E1
  "Arg/Orn/Pro interconversion", "Citrulline metabolism",
  "IMP biosynthesis III", "Guanosine nt. degradation III",
  "Purine ribonucleoside degradation", "Propanoate fermentation",
  "TCA cycle II", "Polyamine biosynthesis II",
  # E4
  "Urate biosyn. / IMP degradation", "Adenosine nt. de novo biosyn.",
  "5-AIR biosynthesis II", "Adenine/adenosine salvage III",
  "Peptidoglycan biosynthesis I", "Branched-chain AA biosyn.",
  "SAM cycle I", "UMP biosynthesis"
)

cat(sprintf("  Selected pathways: %d (E1: %d, E4: %d)\n",
            nrow(pathway_selection),
            sum(pathway_selection$group == "E1-only pathway"),
            sum(pathway_selection$group == "E4-only pathway")))

# ---- 2b. Convergent Shared ECs (all 14) ----
ec_shared <- read.csv(here("results", "tables", "section_3_4", "ec_shared_direction.csv"),
                       stringsAsFactors = FALSE)
convergent_ecs <- ec_shared %>% filter(direction_type == "Convergent")

# Load EC annotations
ec_annot <- read.csv(here("results", "tables", "section_3_4", "ec_kegg_annotations.csv"),
                      stringsAsFactors = FALSE)

convergent_selection <- convergent_ecs %>%
  left_join(ec_annot, by = c("feature" = "EC")) %>%
  mutate(
    group = "Convergent EC",
    # Strip trailing semicolons from the EC annotation Name field
    Name = gsub(";\\s*$", "", Name),
    # EC number only: the full enzyme names are far too long for a
    # half-width column (and are given in Figure 4 / the supplementary tables).
    display_name = paste0("EC ", feature)
  ) %>%
  select(feature, group, display_name, E1_mean_lfc, E4_mean_lfc)

cat(sprintf("  Convergent shared ECs: %d\n", nrow(convergent_selection)))
cat("  Top 3 (from 3.3): EC 2.7.6.3, EC 2.6.1.1, EC 2.1.1.13\n")

# ---- 2c. Curated metabolite targets (same as Panel B) ----
curated_metab <- tribble(
  ~feature_id,                            ~display_name,            ~class,
  "NEG0024_p-Cresol",                     "p-Cresol",               "AAA catabolite",
  "NEG0216_p-Cresol sulfate",             "p-Cresol sulfate",       "AAA catabolite",
  "NEG0046_Indoxyl sulfate",              "Indoxyl sulfate",        "AAA catabolite",
  "NEG0144_5-Methoxyindoleacetate",       "5-MeOIA",               "AAA catabolite",
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
# 3. MATCH SAMPLES & PREPARE MATRICES
# ==============================================================================

cat("\n--- Matching samples ---\n")

micro_ids  <- sample_names(ps_taxa)
path_ids   <- sample_names(ps_path)
ec_ids     <- sample_names(ps_ec)
metab_ids  <- metab_raw$metadata$SampleID
sid_col    <- intersect(c("SampleID", "sampleid", "Sample_ID"), names(clinical))
if (length(sid_col) == 0) sid_col <- names(clinical)[1]
sid_col    <- sid_col[1]
clin_ids   <- clinical[[sid_col]]

shared_ids <- Reduce(intersect, list(micro_ids, path_ids, ec_ids, metab_ids, clin_ids))
cat(sprintf("  Shared samples across all data layers: %d\n", length(shared_ids)))

meta_ps <- as(sample_data(ps_taxa), "data.frame")

# ---- 3a. Pathway abundance matrix (log2 CPM) ----
ps_path_sub <- prune_samples(shared_ids, ps_path)
path_otu <- as(otu_table(ps_path_sub), "matrix")
if (!taxa_are_rows(ps_path_sub)) path_otu <- t(path_otu)

# Convert to CPM then log2
path_cpm <- sweep(path_otu, 2, colSums(path_otu), "/") * 1e6
path_log <- log2(path_cpm + 1)

# Subset to selected pathways
avail_paths <- pathway_selection$feature[pathway_selection$feature %in% rownames(path_log)]
pathway_selection <- pathway_selection %>% filter(feature %in% avail_paths)
path_mat <- path_log[avail_paths, shared_ids]
cat(sprintf("  Pathways matched: %d / %d\n", length(avail_paths), nrow(pathway_selection)))

# ---- 3b. EC abundance matrix (log2 CPM) ----
ps_ec_sub <- prune_samples(shared_ids, ps_ec)
ec_otu <- as(otu_table(ps_ec_sub), "matrix")
if (!taxa_are_rows(ps_ec_sub)) ec_otu <- t(ec_otu)

ec_cpm <- sweep(ec_otu, 2, colSums(ec_otu), "/") * 1e6
ec_log <- log2(ec_cpm + 1)

avail_ecs <- convergent_selection$feature[convergent_selection$feature %in% rownames(ec_log)]
convergent_selection <- convergent_selection %>% filter(feature %in% avail_ecs)
ec_mat <- ec_log[avail_ecs, shared_ids]
cat(sprintf("  Convergent ECs matched: %d / %d\n", length(avail_ecs), nrow(convergent_selection)))

# ---- 3c. Species abundance (relative, for forest plot) ----
ps_taxa_sub <- prune_samples(shared_ids, ps_taxa)
ps_rel <- transform_sample_counts(ps_taxa_sub, function(x) x / sum(x))
taxa_otu <- as(otu_table(ps_rel), "matrix")
if (!taxa_are_rows(ps_rel)) taxa_otu <- t(taxa_otu)

# Load species selection from Panel B (top 23)
e1_core_taxa <- read.csv(here("results", "tables", "section_3_4", "taxa_E1_core.csv"),
                           stringsAsFactors = FALSE)
e4_core_taxa <- read.csv(here("results", "tables", "section_3_4", "taxa_E4_core.csv"),
                           stringsAsFactors = FALSE)
shared_taxa  <- read.csv(here("results", "tables", "section_3_4", "taxa_Shared.csv"),
                           stringsAsFactors = FALSE)
shared_taxa_features <- shared_taxa$feature

e1_core_taxa <- e1_core_taxa %>%
  mutate(
    mean_lfc = rowMeans(cbind(log2FoldChange_E1vE2, log2FoldChange_E1vE3), na.rm = TRUE),
    abs_lfc  = abs(mean_lfc),
    direction = ifelse(mean_lfc > 0, "Higher", "Lower")
  )
e1_only_taxa <- e1_core_taxa %>% filter(!feature %in% shared_taxa_features)
e1_higher <- e1_only_taxa %>% filter(direction == "Higher") %>%
  arrange(desc(abs_lfc)) %>% slice_head(n = 10)
e1_lower  <- e1_only_taxa %>% filter(direction == "Lower") %>%
  arrange(desc(abs_lfc)) %>% slice_head(n = 10)
e4_all_taxa <- e4_core_taxa %>%
  mutate(mean_lfc = rowMeans(cbind(log2FoldChange_E4vE2, log2FoldChange_E4vE3), na.rm = TRUE))

species_selection <- bind_rows(
  e1_higher %>% mutate(group = "Type 1 \u2014 higher"),
  e1_lower  %>% mutate(group = "Type 1 \u2014 lower"),
  e4_all_taxa %>% mutate(group = ifelse(feature %in% shared_taxa_features, "Shared", "Type 4"))
) %>% select(feature, group)
species_selection$display_name <-
  abbrev_species(gsub("_", " ", sub(".*s__", "", species_selection$feature)))

avail_sp <- species_selection$feature[species_selection$feature %in% rownames(taxa_otu)]
species_selection <- species_selection %>% filter(feature %in% avail_sp)
sp_mat <- taxa_otu[avail_sp, shared_ids]
cat(sprintf("  Species matched: %d\n", nrow(sp_mat)))

# ---- 3d. Metabolite + UA target matrix ----
metab_abund <- metab_raw$abundance
avail_metab <- curated_metab$feature_id[curated_metab$feature_id %in% rownames(metab_abund)]
met_mat <- metab_abund[avail_metab, shared_ids]

clin_sub <- clinical[match(shared_ids, clinical[[sid_col]]), ]
ua_vec <- as.numeric(clin_sub$UA)
names(ua_vec) <- shared_ids

target_mat <- rbind(met_mat, UA = ua_vec[shared_ids])
cat(sprintf("  Target metabolites: %d + UA\n", length(avail_metab)))

# ---- 3e. Covariates ----
sex_numeric <- as.numeric(factor(meta_ps[shared_ids, "Gender"]))
age_numeric <- as.numeric(meta_ps[shared_ids, "Age"])
covar_base  <- cbind(sex_numeric, age_numeric)

# Enterotype dummies (for adjusted version)
et_factor <- factor(meta_ps[shared_ids, "Dirichlettype_Genus"])
et_dummies <- model.matrix(~ et_factor)[, -1]  # drop intercept
covar_full <- cbind(covar_base, et_dummies)

n_samples <- length(shared_ids)

# ==============================================================================
# 4. PARTIAL CORRELATION ENGINE
# ==============================================================================

compute_partial_cor <- function(feature_mat, target_mat, covar_mat, label = "") {
  # Rank-residual partial Spearman method (same as Panel B)
  n_feat <- nrow(feature_mat)
  n_tgt  <- nrow(target_mat)
  n_samp <- ncol(feature_mat)

  rho_mat <- matrix(NA, n_feat, n_tgt,
                    dimnames = list(rownames(feature_mat), rownames(target_mat)))
  p_mat   <- matrix(NA, n_feat, n_tgt,
                    dimnames = list(rownames(feature_mat), rownames(target_mat)))

  # Pre-compute rank-residuals for targets
  tgt_resid <- matrix(NA, n_samp, n_tgt)
  for (j in seq_len(n_tgt)) {
    y <- as.numeric(target_mat[j, ])
    valid <- !is.na(y)
    ry <- rep(NA, n_samp)
    ry[valid] <- rank(y[valid], ties.method = "average")
    fit <- lm(ry ~ covar_mat)
    tgt_resid[, j] <- residuals(fit)
  }

  # Pre-compute rank-residuals for features
  feat_resid <- matrix(NA, n_samp, n_feat)
  for (i in seq_len(n_feat)) {
    x <- as.numeric(feature_mat[i, ])
    rx <- rank(x, ties.method = "average")
    fit <- lm(rx ~ covar_mat)
    feat_resid[, i] <- residuals(fit)
  }

  # Pairwise correlations
  df_resid <- n_samp - 2 - ncol(covar_mat)
  for (i in seq_len(n_feat)) {
    for (j in seq_len(n_tgt)) {
      rx <- feat_resid[, i]
      ry <- tgt_resid[, j]
      valid <- !is.na(rx) & !is.na(ry)
      if (sum(valid) < 10) next
      r <- cor(rx[valid], ry[valid], method = "pearson")
      rho_mat[i, j] <- r
      t_stat <- r * sqrt(sum(valid) - 2 - ncol(covar_mat)) / sqrt(1 - r^2)
      p_mat[i, j] <- 2 * pt(abs(t_stat), df = df_resid, lower.tail = FALSE)
    }
  }

  # BH adjustment across all pairs
  p_adj <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                  nrow = n_feat, ncol = n_tgt,
                  dimnames = dimnames(p_mat))

  n_sig <- sum(p_adj < 0.05, na.rm = TRUE)
  cat(sprintf("  %s: %d sig / %d tested (BH < 0.05)\n",
              label, n_sig, sum(!is.na(p_adj))))

  list(rho = rho_mat, pval = p_mat, padj = p_adj)
}

# ==============================================================================
# 5. COMPUTE CORRELATIONS
# ==============================================================================

cat("\n--- Computing partial correlations (sex + age adjusted) ---\n")

# Combine pathways + convergent ECs into one feature matrix
func_mat <- rbind(path_mat, ec_mat)
func_info <- bind_rows(
  pathway_selection %>% select(feature, group, display_name),
  convergent_selection %>% select(feature, group, display_name)
)

res_func_base <- compute_partial_cor(func_mat, target_mat, covar_base,
                                      label = "Function (base)")

# Species (for forest plot)
res_sp_base <- compute_partial_cor(sp_mat, target_mat, covar_base,
                                    label = "Species (base)")

cat("\n--- Computing partial correlations (sex + age + enterotype adjusted) ---\n")

res_func_adj <- compute_partial_cor(func_mat, target_mat, covar_full,
                                     label = "Function (ET-adj)")
res_sp_adj   <- compute_partial_cor(sp_mat, target_mat, covar_full,
                                     label = "Species (ET-adj)")

# ==============================================================================
# 6. PANEL C: FUNCTIONAL HEATMAP
# ==============================================================================

cat("\n--- Building Panel C heatmap ---\n")

rho_c  <- res_func_base$rho
padj_c <- res_func_base$padj

# Row info
row_info_c <- func_info[match(rownames(rho_c), func_info$feature), ]
row_groups_c <- factor(row_info_c$group,
                       levels = c("E1-only pathway", "E4-only pathway", "Convergent EC"))
row_labels_c <- row_info_c$display_name

# E1-/E4-only pathways keep enterotype colors; Convergent EC aligned with
# the shared LAYER_COLORS palette (scripts/utils/palettes.R) used in Figure 6.
group_colors_c <- c(
  "E1-only pathway" = "#E64B35",
  "E4-only pathway" = "#3C5488",
  "Convergent EC"   = "#1B9E77"
)

ha_row_c <- rowAnnotation(
  Layer = row_groups_c,
  col = list(Layer = group_colors_c),
  annotation_name_gp = gpar(fontsize = fs_legend_text),
  simple_anno_size = unit(3, "mm"),
  # No legend: the row-split titles already name each layer.
  show_legend = FALSE
)

# Column info (same as Panel B, with UA separated by gap)
col_info_c <- bind_rows(
  curated_metab %>% filter(feature_id %in% rownames(target_mat)) %>%
    select(feature_id, display_name, class),
  tibble(feature_id = "UA", display_name = "Serum UA", class = "Clinical")
)
col_info_c <- col_info_c[match(rownames(target_mat), col_info_c$feature_id), ]

col_groups_c <- factor(col_info_c$class,
                       levels = c("AAA catabolite", "Bile acid",
                                  "TMAO pathway", "Purine", "Clinical"))
col_labels_c <- col_info_c$display_name

ha_col_c <- HeatmapAnnotation(
  Class = col_groups_c,
  col = list(Class = class_colors),
  annotation_name_gp = gpar(fontsize = fs_legend_text),
  simple_anno_size = unit(3, "mm"),
  # No legend: the column-split titles already name each class.
  show_legend = FALSE
)

max_rho <- 0.5
col_fun <- colorRamp2(
  c(-max_rho, -max_rho / 2, 0, max_rho / 2, max_rho),
  c("#3C5488", "#8FAABE", "white", "#E8927C", "#E64B35")
)

cell_fn_c <- function(j, i, x, y, width, height, fill) {
  val <- padj_c[i, j]
  if (!is.na(val)) {
    if (val < 0.01) {
      grid.text("**", x, y, gp = gpar(fontsize = fs_star, col = "black"))
    } else if (val < 0.05) {
      grid.text("*", x, y, gp = gpar(fontsize = fs_star, col = "black"))
    }
  }
}

set.seed(42)
ht_c <- Heatmap(
  rho_c,
  name = "Partial \u03c1",
  col = col_fun,

  row_labels = row_labels_c,
  row_names_gp = gpar(fontsize = fs_axis_text),
  row_split = row_groups_c,
  row_title_gp = gpar(fontsize = 7.5, fontface = "bold"),
  row_gap = unit(3, "mm"),
  cluster_row_slices = FALSE,
  cluster_rows = TRUE,
  row_dend_width = unit(8, "mm"),

  column_labels = col_labels_c,
  column_names_gp = gpar(fontsize = fs_axis_text),
  column_names_rot = 90,
  column_split = col_groups_c,
  column_title = NULL,
  column_gap = unit(4, "mm"),
  cluster_column_slices = FALSE,
  cluster_columns = FALSE,

  cell_fun = cell_fn_c,
  rect_gp = gpar(col = "grey90", lwd = 0.5),

  left_annotation = ha_row_c,
  top_annotation = ha_col_c,

  heatmap_legend_param = list(
    title = "Partial \u03c1",
    title_gp = gpar(fontsize = fs_legend_title, fontface = "bold"),
    labels_gp = gpar(fontsize = fs_legend_text),
    at = c(-0.5, -0.25, 0, 0.25, 0.5),
    direction = "horizontal",
    grid_height = unit(3, "mm"),
    legend_width = unit(2.4, "cm")
  ),

  width = unit(fig5_C_body[["w"]], "mm"),
  height = unit(fig5_C_body[["h"]], "mm")
)

# Save Panel C
# cairo devices render the Greek rho (U+03C1) correctly under a non-Latin locale
cairo_pdf(here("results", "figures", "panels", "Figure5C_function_heatmap.pdf"),
    width = fig5_C_w_in, height = fig5_C_h_in)
draw(ht_c, heatmap_legend_side = "bottom")
dev.off()

png(here("results", "figures", "panels", "Figure5C_function_heatmap.png"),
    width = fig5_C_w_in, height = fig5_C_h_in, units = "in", res = fig5_DPI,
    type = "cairo")
draw(ht_c, heatmap_legend_side = "bottom")
dev.off()

cat("  Panel C saved.\n")

# ==============================================================================
# 7. PANEL D: FOREST PLOT — Feature-UA Correlations Across Layers
# ==============================================================================

cat("\n--- Building Panel D forest plot ---\n")

# Extract UA column from each analysis
extract_ua_col <- function(res, info_df, layer_label) {
  ua_idx <- which(colnames(res$rho) == "UA")
  data.frame(
    feature      = rownames(res$rho),
    display_name = info_df$display_name[match(rownames(res$rho), info_df$feature)],
    group        = info_df$group[match(rownames(res$rho), info_df$feature)],
    layer        = layer_label,
    rho          = res$rho[, ua_idx],
    pval         = res$pval[, ua_idx],
    padj         = res$padj[, ua_idx],
    stringsAsFactors = FALSE
  )
}

# Base model (sex + age)
ua_sp_base   <- extract_ua_col(res_sp_base, species_selection, "Species")
ua_func_base <- extract_ua_col(res_func_base, func_info, "Function")

# Split function into pathway vs EC
ua_path_base <- ua_func_base %>% filter(grepl("pathway", group))
ua_ec_base   <- ua_func_base %>% filter(group == "Convergent EC")

# ET-adjusted model
ua_sp_adj   <- extract_ua_col(res_sp_adj, species_selection, "Species")
ua_func_adj <- extract_ua_col(res_func_adj, func_info, "Function")
ua_path_adj <- ua_func_adj %>% filter(grepl("pathway", group))
ua_ec_adj   <- ua_func_adj %>% filter(group == "Convergent EC")

# Compute 95% CI from rho and sample size
add_ci <- function(df, n) {
  z <- atanh(df$rho)
  se <- 1 / sqrt(n - 3)
  df$ci_lo <- tanh(z - 1.96 * se)
  df$ci_hi <- tanh(z + 1.96 * se)
  df
}

n_eff <- n_samples  # approximate

ua_sp_base   <- add_ci(ua_sp_base, n_eff)
ua_path_base <- add_ci(ua_path_base, n_eff)
ua_ec_base   <- add_ci(ua_ec_base, n_eff)
ua_sp_adj    <- add_ci(ua_sp_adj, n_eff)
ua_path_adj  <- add_ci(ua_path_adj, n_eff)
ua_ec_adj    <- add_ci(ua_ec_adj, n_eff)

# Combine all for plotting
combine_for_plot <- function(sp_df, path_df, ec_df, model_label) {
  sp_df$sublayer   <- sp_df$group
  path_df$sublayer <- path_df$group
  ec_df$sublayer   <- "Convergent EC"

  sp_df$layer   <- "Species"
  path_df$layer <- "Pathway"
  ec_df$layer   <- "Convergent EC"

  df <- bind_rows(sp_df, path_df, ec_df)
  df$model <- model_label
  df$sig <- ifelse(!is.na(df$padj) & df$padj < 0.05, "sig", "ns")
  df
}

plot_base <- combine_for_plot(ua_sp_base, ua_path_base, ua_ec_base, "Base")
plot_adj  <- combine_for_plot(ua_sp_adj, ua_path_adj, ua_ec_adj, "ET-adj.")

plot_all <- bind_rows(plot_base, plot_adj)
plot_all$model <- factor(plot_all$model,
                         levels = c("Base", "ET-adj."))
plot_all$layer <- factor(plot_all$layer,
                         levels = c("Species", "Pathway", "Convergent EC"))

# Order features within each layer by base-model rho
order_df <- plot_base %>%
  arrange(layer, rho) %>%
  mutate(order_idx = row_number())

plot_all <- plot_all %>%
  left_join(order_df %>% select(feature, layer, order_idx),
            by = c("feature", "layer"))

# Layer colors aligned with the shared LAYER_COLORS palette
# (scripts/utils/palettes.R) so layers match Figure 6:
#   Species -> Taxonomy color; Pathway -> Microbial pathway color.
layer_colors <- c(
  "Species"        = "#7B5EA7",
  "Pathway"        = "#2C7FB8",
  "Convergent EC"  = "#1B9E77"
)

p4d <- ggplot(plot_all, aes(x = rho, y = reorder(display_name, order_idx))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.3, linewidth = 0.4, color = "grey50") +
  geom_point(aes(color = layer, shape = sig), size = 2) +
  scale_color_manual(values = layer_colors, name = "Omics layer") +
  scale_shape_manual(values = c("sig" = 16, "ns" = 1),
                     name = "BH-adj P",
                     labels = c("sig" = "< 0.05", "ns" = "\u2265 0.05")) +
  facet_grid(layer ~ model, scales = "free_y", space = "free_y") +
  labs(x = "Partial Spearman \u03c1 with serum UA",
       y = NULL) +
  theme_fig5 +
  theme(
    # tilt x tick labels so the two side-by-side model facets pack tighter
    axis.text.x = element_text(angle = 30, hjust = 1,
                               size = fs_axis_text, color = "black"),
    legend.position    = "bottom",
    legend.box         = "vertical",
    legend.key.size    = unit(0.32, "cm"),
    legend.margin      = margin(1, 1, 1, 1),
    legend.spacing.y   = unit(1, "pt"),
    legend.box.spacing = unit(2, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank()
  )

ggsave(here("results", "figures", "panels", "Figure5D_UA_forest_plot.pdf"),
       p4d, width = fig5_D_w_in, height = fig5_D_h_in, units = "in",
       device = cairo_pdf)
ggsave(here("results", "figures", "panels", "Figure5D_UA_forest_plot.png"),
       p4d, width = fig5_D_w_in, height = fig5_D_h_in, units = "in",
       dpi = fig5_DPI, device = png, type = "cairo")

cat("  Panel D saved.\n")

# ==============================================================================
# 8. SAVE TABLES
# ==============================================================================

cat("\n--- Saving tables ---\n")

# Panel C correlation tables
write.csv(
  as.data.frame(res_func_base$rho) %>% tibble::rownames_to_column("feature"),
  here("results", "tables", "section_3_4", "function_metabolite_partial_rho.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame(res_func_base$padj) %>% tibble::rownames_to_column("feature"),
  here("results", "tables", "section_3_4", "function_metabolite_partial_padj.csv"),
  row.names = FALSE
)

# Panel D summary table
ua_summary <- plot_all %>%
  select(feature, display_name, layer, group, model, rho, ci_lo, ci_hi, pval, padj, sig)
write.csv(ua_summary,
          here("results", "tables", "section_3_4", "feature_UA_partial_cor_summary.csv"),
          row.names = FALSE)

# ET-adjusted correlations
write.csv(
  as.data.frame(res_func_adj$rho) %>% tibble::rownames_to_column("feature"),
  here("results", "tables", "section_3_4", "function_metabolite_partial_rho_ETadj.csv"),
  row.names = FALSE
)

# Frozen plot objects allow layout-only rebuilds without repeating db-RDA/envfit.
saveRDS(ht_c, here("results", "cache", "p4c_function_heatmap.rds"))
saveRDS(plot_all, here("results", "cache", "p4d_plot_all.rds"))
saveRDS(layer_colors, here("results", "cache", "figure5_layer_colors.rds"))

cat("\n========================================\n")
cat("Section 3.4 functional integration complete.\n")
cat("Outputs:\n")
cat("  Panel C: results/figures/panels/Figure5C_function_heatmap.*\n")
cat("  Panel D: results/figures/panels/Figure5D_UA_forest_plot.*\n")
cat("  Tables:  results/tables/section_3_4/\n")
cat("========================================\n")
