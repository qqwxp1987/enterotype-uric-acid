# ==============================================================================
# Section 3.3: Differential Abundance Analysis & Independent Core Feature Extraction
#
# Strategy:
#   Type 1 and Type 4 both show elevated HUA prevalence (Section 3.2).
#   Types 2 & 3 serve as "low-HUA reference" groups.
#   We identify features independently enriched/depleted in Type 1 and Type 4
#   relative to BOTH middle types, then compare E1 vs E4 signatures.
#
# Methods:
#   - LinDA (MicrobiomeStat): Taxa (species), Pathway, EC — with Sex+Age covariates
#   - limma: Metabolome — with Sex+Age covariates
#
# Data source: ALL enterotype labels from phyloseq_Taxa_with_enterotype_seed10.rds
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(MicrobiomeStat)
library(limma)
library(here)

cat("========================================\n")
cat("Section 3.3 Differential Analysis\n")
cat("========================================\n\n")

# --- 0. Output directory ---
out_dir <- here("results", "tables", "section_3_3")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1. Load authoritative enterotype mapping from seed10 ---
ps_taxa <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta_auth <- as(sample_data(ps_taxa), "data.frame")
et_map <- data.frame(
  SampleID = rownames(meta_auth),
  Enterotype = meta_auth$Dirichlettype_Genus,
  Gender = meta_auth$Gender,
  Age = as.numeric(meta_auth$Age),
  stringsAsFactors = FALSE
)
cat("Authoritative enterotype map: ", nrow(et_map), " samples\n")
cat("Distribution:\n")
print(table(et_map$Enterotype))

# --- 2. Prepare phyloseq objects with seed10 enterotypes ---
update_enterotype <- function(ps, et_map) {
  meta <- as(sample_data(ps), "data.frame")
  shared <- intersect(rownames(meta), et_map$SampleID)
  # Update enterotype labels
  meta[shared, "Dirichlettype_Genus"] <- et_map$Enterotype[match(shared, et_map$SampleID)]
  sample_data(ps) <- sample_data(meta)
  return(ps)
}

ps_path <- readRDS(here("data", "clean", "phyloseq_Pathway_with_enterotype_20250210.rds"))
ps_ec   <- readRDS(here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"))

ps_path <- update_enterotype(ps_path, et_map)
ps_ec   <- update_enterotype(ps_ec, et_map)

cat("\nEnterotype labels updated for Pathway and EC.\n")

# --- 3. Prepare metabolome data ---
metab_raw <- readRDS(here("data", "clean", "metabolites_normalized_multiIS_and_metadata_20250904.rds"))
metab_abund <- metab_raw$abundance   # 829 metabolites x 303 samples
metab_meta  <- metab_raw$metadata

# Map enterotype from seed10 via SampleID
metab_meta$Enterotype_seed10 <- et_map$Enterotype[match(metab_meta$SampleID, et_map$SampleID)]
# Also ensure Sex and Age from seed10 (for consistency)
metab_meta$Gender_seed10 <- et_map$Gender[match(metab_meta$SampleID, et_map$SampleID)]
metab_meta$Age_seed10 <- et_map$Age[match(metab_meta$SampleID, et_map$SampleID)]

# Remove samples without seed10 mapping
valid_metab <- !is.na(metab_meta$Enterotype_seed10)
metab_meta  <- metab_meta[valid_metab, ]
metab_abund <- metab_abund[, metab_meta$SampleID]

cat(sprintf("\nMetabolome: %d metabolites x %d samples (with seed10 mapping)\n",
            nrow(metab_abund), ncol(metab_abund)))
cat("Metabolome enterotype distribution:\n")
print(table(metab_meta$Enterotype_seed10))

# ==============================================================================
# --- 4. LinDA Pairwise Differential Analysis ---
# ==============================================================================

comparisons <- list(
  c("1", "2"),  # E1 vs E2
  c("1", "3"),  # E1 vs E3
  c("4", "2"),  # E4 vs E2
  c("4", "3")   # E4 vs E3
)

run_linda_pairwise <- function(ps, comp, data_label) {
  # Add ET column to phyloseq sample_data BEFORE subsetting
  meta <- as(sample_data(ps), "data.frame")
  meta$ET <- meta$Dirichlettype_Genus
  meta$Gender <- factor(meta$Gender)
  meta$Age <- as.numeric(meta$Age)
  sample_data(ps) <- sample_data(meta)

  # Subset to two groups
  keep <- sample_data(ps)$ET %in% comp
  ps_sub <- prune_samples(keep, ps)
  meta_sub <- as(sample_data(ps_sub), "data.frame")

  # Set reference level (second element = reference)
  meta_sub$ET <- factor(meta_sub$ET, levels = c(comp[2], comp[1]))
  sample_data(ps_sub) <- sample_data(meta_sub)

  comp_label <- paste0("E", comp[1], "v", "E", comp[2])
  cat(sprintf("  LinDA %s [%s]: n=%d (%s=%d, %s=%d) ... ",
              data_label, comp_label,
              nrow(meta_sub),
              comp[1], sum(meta_sub$ET == comp[1]),
              comp[2], sum(meta_sub$ET == comp[2])))

  # Extract OTU table and metadata separately for LinDA
  otu <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) otu <- t(otu)  # Ensure features x samples

  # Run LinDA
  res <- tryCatch({
    linda(
      feature.dat = otu,
      meta.dat = meta_sub,
      formula = "~ ET + Gender + Age",
      feature.dat.type = "count",
      is.winsor = TRUE,
      outlier.pct = 0.03,
      p.adj.method = "BH",
      alpha = 0.05,
      prev.filter = 0.1,
      mean.abund.filter = 0
    )
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    return(NULL)
  })

  if (is.null(res)) return(NULL)

  # Extract results for the ET variable
  # LinDA output: res$output is a list; the ET effect is in the element named like "ETx"
  et_name <- paste0("ET", comp[1])
  if (et_name %in% names(res$output)) {
    df <- res$output[[et_name]]
  } else {
    # Try first non-intercept element
    et_name <- names(res$output)[1]
    df <- res$output[[et_name]]
  }

  df$feature <- rownames(df)
  df$comparison <- comp_label
  df$data_type <- data_label
  df$sig <- ifelse(df$padj < 0.05, "Yes", "No")

  n_sig <- sum(df$padj < 0.05, na.rm = TRUE)
  n_up  <- sum(df$padj < 0.05 & df$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(df$padj < 0.05 & df$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("%d sig (↑%d ↓%d) / %d tested\n", n_sig, n_up, n_down, nrow(df)))

  return(df)
}

cat("\n=== LinDA: Taxa (Species) ===\n")
linda_taxa <- lapply(comparisons, function(comp) run_linda_pairwise(ps_taxa, comp, "Taxa"))

cat("\n=== LinDA: Pathway ===\n")
linda_path <- lapply(comparisons, function(comp) run_linda_pairwise(ps_path, comp, "Pathway"))

cat("\n=== LinDA: EC ===\n")
linda_ec <- lapply(comparisons, function(comp) run_linda_pairwise(ps_ec, comp, "EC"))

# ==============================================================================
# --- 5. limma Metabolome Differential Analysis ---
# ==============================================================================

run_limma_pairwise <- function(abund, meta, comp) {
  # Subset to two groups
  keep <- meta$Enterotype_seed10 %in% comp
  meta_sub <- meta[keep, ]
  abund_sub <- abund[, meta_sub$SampleID]

  # log2 transform (add small pseudocount for zeros)
  min_nonzero <- min(abund_sub[abund_sub > 0], na.rm = TRUE)
  abund_log <- log2(abund_sub + min_nonzero / 2)

  # Design matrix
  group <- factor(meta_sub$Enterotype_seed10, levels = c(comp[2], comp[1]))
  sex   <- factor(meta_sub$Gender_seed10)
  age   <- as.numeric(meta_sub$Age_seed10)
  design <- model.matrix(~ group + sex + age)

  comp_label <- paste0("E", comp[1], "v", "E", comp[2])
  cat(sprintf("  limma Metabolome [%s]: n=%d (%s=%d, %s=%d) ... ",
              comp_label,
              nrow(meta_sub),
              comp[1], sum(meta_sub$Enterotype_seed10 == comp[1]),
              comp[2], sum(meta_sub$Enterotype_seed10 == comp[2])))

  # Fit limma
  fit <- lmFit(abund_log, design)
  fit <- eBayes(fit)

  # Extract group effect (column 2 in design)
  res <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  res$feature <- rownames(res)
  res$comparison <- comp_label
  res$data_type <- "Metabolome"
  res$sig <- ifelse(res$adj.P.Val < 0.05, "Yes", "No")

  # Rename columns to match LinDA output
  res <- res %>%
    rename(log2FoldChange = logFC, pvalue = P.Value, padj = adj.P.Val)

  n_sig <- sum(res$padj < 0.05, na.rm = TRUE)
  n_up  <- sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE)
  n_down <- sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE)
  cat(sprintf("%d sig (↑%d ↓%d) / %d tested\n", n_sig, n_up, n_down, nrow(res)))

  return(res)
}

cat("\n=== limma: Metabolome ===\n")
limma_metab <- lapply(comparisons, function(comp) {
  run_limma_pairwise(metab_abund, metab_meta, comp)
})

# ==============================================================================
# --- 6. Combine and Save All Differential Results ---
# ==============================================================================

combine_results <- function(res_list, label) {
  valid <- Filter(Negate(is.null), res_list)
  if (length(valid) == 0) return(NULL)
  # Keep only common columns
  common_cols <- c("feature", "log2FoldChange", "pvalue", "padj", "comparison", "data_type", "sig")
  valid <- lapply(valid, function(df) {
    # Ensure columns exist
    missing <- setdiff(common_cols, names(df))
    for (col in missing) df[[col]] <- NA
    df[, common_cols]
  })
  do.call(rbind, valid)
}

all_taxa  <- combine_results(linda_taxa, "Taxa")
all_path  <- combine_results(linda_path, "Pathway")
all_ec    <- combine_results(linda_ec, "EC")
all_metab <- combine_results(limma_metab, "Metabolome")

# Save individual data type results
if (!is.null(all_taxa))  write.csv(all_taxa,  file.path(out_dir, "linda_taxa_all_comparisons.csv"), row.names = FALSE)
if (!is.null(all_path))  write.csv(all_path,  file.path(out_dir, "linda_pathway_all_comparisons.csv"), row.names = FALSE)
if (!is.null(all_ec))    write.csv(all_ec,    file.path(out_dir, "linda_ec_all_comparisons.csv"), row.names = FALSE)
if (!is.null(all_metab)) write.csv(all_metab, file.path(out_dir, "limma_metabolome_all_comparisons.csv"), row.names = FALSE)

cat("\nAll differential results saved.\n")

# ==============================================================================
# --- 7. Independent Core Feature Extraction ---
# ==============================================================================

extract_core <- function(results_df) {
  if (is.null(results_df)) return(list(E1_core = character(0), E4_core = character(0),
                                        Shared = character(0), E1_only = character(0), E4_only = character(0)))

  # Get significant features for each comparison
  sig_E1vE2 <- results_df %>% filter(comparison == "E1vE2", padj < 0.05) %>% pull(feature)
  sig_E1vE3 <- results_df %>% filter(comparison == "E1vE3", padj < 0.05) %>% pull(feature)
  sig_E4vE2 <- results_df %>% filter(comparison == "E4vE2", padj < 0.05) %>% pull(feature)
  sig_E4vE3 <- results_df %>% filter(comparison == "E4vE3", padj < 0.05) %>% pull(feature)

  # Core = significant in BOTH comparisons against middle types
  E1_core <- intersect(sig_E1vE2, sig_E1vE3)
  E4_core <- intersect(sig_E4vE2, sig_E4vE3)

  # Additionally check direction consistency for E1_core
  if (length(E1_core) > 0) {
    dir_E1vE2 <- results_df %>% filter(comparison == "E1vE2", feature %in% E1_core) %>%
      select(feature, lfc1 = log2FoldChange)
    dir_E1vE3 <- results_df %>% filter(comparison == "E1vE3", feature %in% E1_core) %>%
      select(feature, lfc2 = log2FoldChange)
    dir_check <- merge(dir_E1vE2, dir_E1vE3, by = "feature")
    # Keep only features with consistent direction
    E1_core <- dir_check %>% filter(sign(lfc1) == sign(lfc2)) %>% pull(feature)
  }

  # Same for E4_core
  if (length(E4_core) > 0) {
    dir_E4vE2 <- results_df %>% filter(comparison == "E4vE2", feature %in% E4_core) %>%
      select(feature, lfc1 = log2FoldChange)
    dir_E4vE3 <- results_df %>% filter(comparison == "E4vE3", feature %in% E4_core) %>%
      select(feature, lfc2 = log2FoldChange)
    dir_check <- merge(dir_E4vE2, dir_E4vE3, by = "feature")
    E4_core <- dir_check %>% filter(sign(lfc1) == sign(lfc2)) %>% pull(feature)
  }

  # Shared vs exclusive
  Shared  <- intersect(E1_core, E4_core)
  E1_only <- setdiff(E1_core, Shared)
  E4_only <- setdiff(E4_core, Shared)

  list(
    E1_core = E1_core, E4_core = E4_core,
    Shared = Shared, E1_only = E1_only, E4_only = E4_only
  )
}

cat("\n========================================\n")
cat("Independent Core Feature Extraction\n")
cat("========================================\n\n")

core_taxa  <- extract_core(all_taxa)
core_path  <- extract_core(all_path)
core_ec    <- extract_core(all_ec)
core_metab <- extract_core(all_metab)

# Print summary
print_core_summary <- function(core, label) {
  cat(sprintf("\n--- %s ---\n", label))
  cat(sprintf("  E1_core: %d features (sig in E1vE2 AND E1vE3, consistent direction)\n", length(core$E1_core)))
  cat(sprintf("  E4_core: %d features (sig in E4vE2 AND E4vE3, consistent direction)\n", length(core$E4_core)))
  cat(sprintf("  Shared:  %d features (in both E1 and E4 core)\n", length(core$Shared)))
  cat(sprintf("  E1_only: %d features (E1 core only)\n", length(core$E1_only)))
  cat(sprintf("  E4_only: %d features (E4 core only)\n", length(core$E4_only)))
}

print_core_summary(core_taxa, "Taxa (Species)")
print_core_summary(core_path, "Pathway")
print_core_summary(core_ec, "EC")
print_core_summary(core_metab, "Metabolome")

# --- 8. Save core feature summary ---
summary_df <- data.frame(
  data_type = c("Taxa", "Pathway", "EC", "Metabolome"),
  E1_core = c(length(core_taxa$E1_core), length(core_path$E1_core),
              length(core_ec$E1_core), length(core_metab$E1_core)),
  E4_core = c(length(core_taxa$E4_core), length(core_path$E4_core),
              length(core_ec$E4_core), length(core_metab$E4_core)),
  Shared = c(length(core_taxa$Shared), length(core_path$Shared),
             length(core_ec$Shared), length(core_metab$Shared)),
  E1_only = c(length(core_taxa$E1_only), length(core_path$E1_only),
              length(core_ec$E1_only), length(core_metab$E1_only)),
  E4_only = c(length(core_taxa$E4_only), length(core_path$E4_only),
              length(core_ec$E4_only), length(core_metab$E4_only))
)

cat("\n\n=== Core Feature Summary Table ===\n")
print(summary_df)

write.csv(summary_df, file.path(out_dir, "independent_core_summary.csv"), row.names = FALSE)

# --- 9. Save detailed core feature lists ---
save_core_list <- function(core, results_df, label, out_dir) {
  # For each category, save feature details with effect sizes
  for (cat_name in c("E1_core", "E4_core", "Shared", "E1_only", "E4_only")) {
    features <- core[[cat_name]]
    if (length(features) == 0) next

    detail <- results_df %>%
      filter(feature %in% features) %>%
      select(feature, log2FoldChange, padj, comparison) %>%
      pivot_wider(
        names_from = comparison,
        values_from = c(log2FoldChange, padj),
        names_sep = "_"
      ) %>%
      arrange(feature)

    fname <- sprintf("%s_%s.csv", tolower(label), cat_name)
    write.csv(detail, file.path(out_dir, fname), row.names = FALSE)
  }
}

save_core_list(core_taxa, all_taxa, "Taxa", out_dir)
save_core_list(core_path, all_path, "Pathway", out_dir)
save_core_list(core_ec, all_ec, "EC", out_dir)
save_core_list(core_metab, all_metab, "Metabolome", out_dir)

# --- 10. Direction analysis for shared features ---
cat("\n\n=== Direction Analysis for Shared Core Features ===\n")

analyze_shared_direction <- function(core, results_df, label) {
  shared <- core$Shared
  if (length(shared) == 0) {
    cat(sprintf("\n%s: No shared features.\n", label))
    return(NULL)
  }

  cat(sprintf("\n--- %s: %d shared features ---\n", label, length(shared)))

  # Get average LFC for E1 and E4 comparisons
  e1_lfc <- results_df %>%
    filter(feature %in% shared, comparison %in% c("E1vE2", "E1vE3")) %>%
    group_by(feature) %>%
    summarise(E1_mean_lfc = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")

  e4_lfc <- results_df %>%
    filter(feature %in% shared, comparison %in% c("E4vE2", "E4vE3")) %>%
    group_by(feature) %>%
    summarise(E4_mean_lfc = mean(log2FoldChange, na.rm = TRUE), .groups = "drop")

  dir_df <- merge(e1_lfc, e4_lfc, by = "feature")
  dir_df$same_direction <- sign(dir_df$E1_mean_lfc) == sign(dir_df$E4_mean_lfc)
  dir_df$direction_type <- ifelse(dir_df$same_direction, "Convergent", "Divergent")

  cat(sprintf("  Convergent (same direction in E1 & E4): %d\n",
              sum(dir_df$direction_type == "Convergent")))
  cat(sprintf("  Divergent (opposite direction): %d\n",
              sum(dir_df$direction_type == "Divergent")))

  write.csv(dir_df, file.path(out_dir, sprintf("%s_shared_direction.csv", tolower(label))),
            row.names = FALSE)

  return(dir_df)
}

dir_taxa  <- analyze_shared_direction(core_taxa, all_taxa, "Taxa")
dir_path  <- analyze_shared_direction(core_path, all_path, "Pathway")
dir_ec    <- analyze_shared_direction(core_ec, all_ec, "EC")
dir_metab <- analyze_shared_direction(core_metab, all_metab, "Metabolome")

cat("\n\n========================================\n")
cat("Section 3.3 analysis complete.\n")
cat(sprintf("Results saved to: %s\n", out_dir))
cat("========================================\n")
