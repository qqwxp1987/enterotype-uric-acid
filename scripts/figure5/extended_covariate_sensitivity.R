# Extended-covariate sensitivity analysis
#
# Repeats differential-abundance, cross-omics and direct-HUA analyses with BMI
# and creatinine added to sex and age. Matched complete-case samples are used so
# covariate effects are distinguished from changes in analysis population.

suppressPackageStartupMessages({
  library(phyloseq)
  library(tidyverse)
  library(MicrobiomeStat)
  library(limma)
  library(vegan)
  library(here)
})

set.seed(2024)

out_dir <- here("results", "covariate_sensitivity")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(out_dir, "covariate_sensitivity_run.log")
log_con <- file(log_file, open = "wt", encoding = "UTF-8")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Extended-covariate sensitivity analysis\n")
cat("Started: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
cat("Seed: 2024\n")
cat("Baseline: Gender + Age\n")
cat("Extended: Gender + Age + BMI + Cr\n\n")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(path, algo = "sha256", file = TRUE))
  }
  unname(tools::md5sum(path))
}

write_session_info <- function() {
  con <- file(file.path(out_dir, "session_info.txt"), open = "wt")
  sink(con)
  print(Sys.info())
  print(sessionInfo())
  sink()
  close(con)
}

record_failure <- local({
  failures <- tibble(
    analysis = character(), model = character(), comparison = character(),
    feature = character(), message = character()
  )
  function(action = c("add", "get"), analysis = NA_character_,
           model = NA_character_, comparison = NA_character_,
           feature = NA_character_, message = NA_character_) {
    action <- match.arg(action)
    if (action == "add") {
      failures <<- bind_rows(failures, tibble(
        analysis = analysis, model = model, comparison = comparison,
        feature = feature, message = message
      ))
    }
    failures
  }
})

safe_max_diff <- function(x, y) {
  both_na <- is.na(x) & is.na(y)
  one_na <- xor(is.na(x), is.na(y))
  if (any(one_na)) return(Inf)
  z <- abs(x[!both_na] - y[!both_na])
  if (length(z) == 0) 0 else max(z)
}

baseline_checks <- tibble(
  component = character(), keys_match = logical(), n_current = integer(),
  n_reference = integer(), max_effect_diff = double(), max_padj_diff = double(),
  tolerance = double(), pass = logical(), note = character()
)

add_baseline_check <- function(component, keys_match, n_current, n_reference,
                               max_effect_diff, max_padj_diff,
                               tolerance = 1e-8, note = "") {
  baseline_checks <<- bind_rows(baseline_checks, tibble(
    component = component, keys_match = keys_match,
    n_current = as.integer(n_current), n_reference = as.integer(n_reference),
    max_effect_diff = as.numeric(max_effect_diff),
    max_padj_diff = as.numeric(max_padj_diff), tolerance = tolerance,
    pass = isTRUE(keys_match) && isTRUE(n_current == n_reference) &&
      is.finite(max_effect_diff) && max_effect_diff <= tolerance &&
      is.finite(max_padj_diff) && max_padj_diff <= tolerance,
    note = note
  ))
}

compare_tidy_baseline <- function(current, reference, component,
                                  key_cols, effect_col, padj_col,
                                  tolerance = 1e-8) {
  cur <- current %>% arrange(across(all_of(key_cols)))
  ref <- reference %>% arrange(across(all_of(key_cols)))
  key_signature <- function(x) {
    do.call(paste, c(lapply(x[, key_cols, drop = FALSE], as.character),
                     sep = "\r"))
  }
  keys_match <- nrow(cur) == nrow(ref) &&
    identical(key_signature(cur), key_signature(ref))
  if (!keys_match) {
    add_baseline_check(component, FALSE, nrow(cur), nrow(ref), Inf, Inf,
                       tolerance, "Feature/comparison keys differ")
    return(invisible(FALSE))
  }
  add_baseline_check(
    component, TRUE, nrow(cur), nrow(ref),
    safe_max_diff(as.numeric(cur[[effect_col]]), as.numeric(ref[[effect_col]])),
    safe_max_diff(as.numeric(cur[[padj_col]]), as.numeric(ref[[padj_col]])),
    tolerance
  )
}

matrix_to_tidy <- function(mat, row_name = "feature", value_name = "value") {
  as.data.frame(mat, check.names = FALSE) %>%
    tibble::rownames_to_column(row_name) %>%
    pivot_longer(-all_of(row_name), names_to = "target", values_to = value_name)
}

read_saved_matrix <- function(path, row_col, drop_cols = character()) {
  x <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  rn <- x[[row_col]]
  x <- x[, setdiff(names(x), c(row_col, drop_cols)), drop = FALSE]
  mat <- as.matrix(data.frame(lapply(x, as.numeric), check.names = FALSE))
  rownames(mat) <- rn
  mat
}

# ------------------------------------------------------------------------------
# Inputs and immutable input manifest
# ------------------------------------------------------------------------------

input_paths <- c(
  analysis_script = here("scripts", "figure5", "extended_covariate_sensitivity.R"),
  ps_taxa = here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"),
  ps_pathway = here("data", "clean", "phyloseq_Pathway_with_enterotype_20250210.rds"),
  ps_ec = here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"),
  metabolome = here("data", "clean", "metabolites_normalized_multiIS_and_metadata_20250904.rds"),
  clinical = here("data", "clean", "Metadata_20240925.csv"),
  mediation_data = here("results", "tables", "section_3_5", "mediation_analysis_data.rds"),
  da_taxa_reference = here("results", "tables", "section_3_4", "linda_taxa_all_comparisons.csv"),
  da_pathway_reference = here("results", "tables", "section_3_4", "linda_pathway_all_comparisons.csv"),
  da_ec_reference = here("results", "tables", "section_3_4", "linda_ec_all_comparisons.csv"),
  da_metabolome_reference = here("results", "tables", "section_3_4", "limma_metabolome_all_comparisons.csv"),
  species_rho_reference = here("results", "tables", "section_3_4", "top_species_metabolite_partial_rho.csv"),
  species_padj_reference = here("results", "tables", "section_3_4", "top_species_metabolite_partial_padj.csv"),
  function_rho_reference = here("results", "tables", "section_3_4", "function_metabolite_partial_rho.csv"),
  function_padj_reference = here("results", "tables", "section_3_4", "function_metabolite_partial_padj.csv"),
  dbrda_reference = here("results", "cache", "dbrda_result.rds"),
  direct_feature_reference = here("results", "tables", "supplementary", "direct_feature_UA_association.csv")
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0) {
  stop("Missing required inputs: ", paste(names(missing_inputs), collapse = ", "))
}

input_manifest <- tibble(
  name = names(input_paths), path = normalizePath(input_paths, winslash = "/"),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1))
)
write.csv(input_manifest, file.path(out_dir, "input_manifest.csv"), row.names = FALSE)

cat("--- Loading data ---\n")
ps_taxa <- readRDS(input_paths[["ps_taxa"]])
ps_path <- readRDS(input_paths[["ps_pathway"]])
ps_ec <- readRDS(input_paths[["ps_ec"]])
metab_raw <- readRDS(input_paths[["metabolome"]])
clinical <- read.csv(input_paths[["clinical"]], stringsAsFactors = FALSE)
med_df <- readRDS(input_paths[["mediation_data"]])

meta_auth <- as(sample_data(ps_taxa), "data.frame")
required_meta <- c("Dirichlettype_Genus", "Gender", "Age", "BMI", "Cr", "eGFR")
if (!all(required_meta %in% names(meta_auth))) {
  stop("Authoritative metadata lacks: ",
       paste(setdiff(required_meta, names(meta_auth)), collapse = ", "))
}

meta_master <- tibble(
  SampleID = rownames(meta_auth),
  Enterotype = as.character(meta_auth$Dirichlettype_Genus),
  Gender = as.character(meta_auth$Gender),
  Age = as.numeric(meta_auth$Age), BMI = as.numeric(meta_auth$BMI),
  Cr = as.numeric(meta_auth$Cr), eGFR = as.numeric(meta_auth$eGFR)
)

sid_col_global <- intersect(c("SampleID", "sampleid", "Sample_ID"), names(clinical))[1]
if (is.na(sid_col_global)) sid_col_global <- names(clinical)[1]
idx_clinical <- match(meta_master$SampleID, clinical[[sid_col_global]])
meta_master$UA <- as.numeric(clinical$UA[idx_clinical])
idx_hua <- match(meta_master$SampleID, med_df$SampleID)
meta_master$HUA <- as.numeric(med_df$HUA[idx_hua])
# The mediation table covers the matched omics cohort. For any remaining
# clinical-only sample, use the manuscript's sex-specific HUA definition.
is_male <- tolower(meta_master$Gender) %in% c("male", "m", "man", "1", "男")
hua_from_threshold <- ifelse(
  is.na(meta_master$UA), NA_real_,
  as.numeric(ifelse(is_male, meta_master$UA >= 420, meta_master$UA >= 360))
)
meta_master$HUA[is.na(meta_master$HUA)] <-
  hua_from_threshold[is.na(meta_master$HUA)]
meta_master$BMI_complete <- !is.na(meta_master$BMI)

missingness <- meta_master %>%
  summarise(across(c(Gender, Age, BMI, Cr, eGFR), ~ sum(is.na(.x)))) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "n_missing") %>%
  mutate(n_total = nrow(meta_master), pct_missing = n_missing / n_total * 100)
write.csv(missingness, file.path(out_dir, "covariate_missingness.csv"), row.names = FALSE)
print(missingness)

# Complete-case selection check: quantify whether the BMI-complete subset
# differs systematically from excluded participants. This diagnostic does not
# alter the locked analysis set.
missingness_by_enterotype <- meta_master %>%
  group_by(Enterotype) %>%
  summarise(
    n_total = n(), n_BMI_complete = sum(BMI_complete),
    BMI_complete_rate = mean(BMI_complete),
    HUA_complete_n = sum(HUA[BMI_complete] == 1, na.rm = TRUE),
    HUA_missing_n = sum(HUA[!BMI_complete] == 1, na.rm = TRUE),
    mean_UA_complete = mean(UA[BMI_complete], na.rm = TRUE),
    mean_UA_missing = mean(UA[!BMI_complete], na.rm = TRUE),
    .groups = "drop"
  )
write.csv(missingness_by_enterotype,
          file.path(out_dir, "bmi_missingness_by_enterotype.csv"), row.names = FALSE)

safe_test <- function(name, expr, effect_note) {
  ans <- tryCatch(expr, error = function(e) NULL)
  tibble(
    comparison = name,
    statistic = if (is.null(ans)) NA_real_ else unname(ans$statistic[1]),
    pvalue = if (is.null(ans)) NA_real_ else ans$p.value,
    method = if (is.null(ans)) "FAILED" else ans$method,
    effect_note = effect_note
  )
}
complete_case_tests <- bind_rows(
  safe_test(
    "BMI completeness by enterotype",
    chisq.test(table(meta_master$Enterotype, meta_master$BMI_complete)),
    paste(sprintf("E%s=%.1f%%", missingness_by_enterotype$Enterotype,
                  100 * missingness_by_enterotype$BMI_complete_rate), collapse = "; ")
  ),
  safe_test(
    "HUA by BMI completeness",
    fisher.test(table(meta_master$HUA, meta_master$BMI_complete)),
    sprintf("HUA complete=%.1f%%; missing=%.1f%%",
            100 * mean(meta_master$HUA[meta_master$BMI_complete] == 1, na.rm = TRUE),
            100 * mean(meta_master$HUA[!meta_master$BMI_complete] == 1, na.rm = TRUE))
  ),
  safe_test(
    "UA by BMI completeness",
    wilcox.test(UA ~ BMI_complete, data = meta_master, exact = FALSE),
    sprintf("mean UA complete=%.2f; missing=%.2f",
            mean(meta_master$UA[meta_master$BMI_complete], na.rm = TRUE),
            mean(meta_master$UA[!meta_master$BMI_complete], na.rm = TRUE))
  ),
  safe_test(
    "Sex by BMI completeness",
    chisq.test(table(meta_master$Gender, meta_master$BMI_complete)),
    "Pearson chi-squared test"
  ),
  safe_test(
    "Age by BMI completeness",
    wilcox.test(Age ~ BMI_complete, data = meta_master, exact = FALSE),
    sprintf("mean age complete=%.2f; missing=%.2f",
            mean(meta_master$Age[meta_master$BMI_complete], na.rm = TRUE),
            mean(meta_master$Age[!meta_master$BMI_complete], na.rm = TRUE))
  ),
  safe_test(
    "Cr by BMI completeness",
    wilcox.test(Cr ~ BMI_complete, data = meta_master, exact = FALSE),
    sprintf("mean Cr complete=%.2f; missing=%.2f",
            mean(meta_master$Cr[meta_master$BMI_complete], na.rm = TRUE),
            mean(meta_master$Cr[!meta_master$BMI_complete], na.rm = TRUE))
  )
)
write.csv(complete_case_tests,
          file.path(out_dir, "bmi_complete_case_selection_tests.csv"), row.names = FALSE)
cat("\nBMI complete-case selection diagnostics:\n")
print(missingness_by_enterotype, n = Inf)
print(complete_case_tests, n = Inf)

sync_metadata <- function(ps, master) {
  meta <- as(sample_data(ps), "data.frame")
  idx <- match(rownames(meta), master$SampleID)
  meta$ET <- master$Enterotype[idx]
  meta$Gender <- factor(master$Gender[idx])
  meta$Age <- master$Age[idx]
  meta$BMI <- master$BMI[idx]
  meta$Cr <- master$Cr[idx]
  meta$eGFR <- master$eGFR[idx]
  sample_data(ps) <- sample_data(meta)
  ps
}

ps_taxa <- sync_metadata(ps_taxa, meta_master)
ps_path <- sync_metadata(ps_path, meta_master)
ps_ec <- sync_metadata(ps_ec, meta_master)

# Metabolome metadata is synchronized to the same authoritative sample map.
metab_abund <- metab_raw$abundance
metab_meta <- metab_raw$metadata
idx_met <- match(metab_meta$SampleID, meta_master$SampleID)
metab_meta$ET <- meta_master$Enterotype[idx_met]
metab_meta$Gender <- factor(meta_master$Gender[idx_met])
metab_meta$Age <- meta_master$Age[idx_met]
metab_meta$BMI <- meta_master$BMI[idx_met]
metab_meta$Cr <- meta_master$Cr[idx_met]
metab_meta$eGFR <- meta_master$eGFR[idx_met]

base_ids <- meta_master %>%
  filter(complete.cases(Gender, Age, Enterotype)) %>% pull(SampleID)
cc_ids <- meta_master %>%
  filter(complete.cases(Gender, Age, BMI, Cr, Enterotype)) %>% pull(SampleID)
cat(sprintf("Base-eligible samples: %d; BMI+Cr complete cases: %d\n\n",
            length(base_ids), length(cc_ids)))

# ------------------------------------------------------------------------------
# Differential abundance engine
# ------------------------------------------------------------------------------

comparisons <- list(c("1", "2"), c("1", "3"), c("4", "2"), c("4", "3"))

run_linda_pairwise <- function(ps, comp, data_label, covariates,
                               eligible_ids, model_label) {
  meta <- as(sample_data(ps), "data.frame")
  keep <- rownames(meta) %in% eligible_ids & meta$ET %in% comp
  ps_sub <- prune_samples(keep, ps)
  meta_sub <- as(sample_data(ps_sub), "data.frame")
  meta_sub$ET <- factor(meta_sub$ET, levels = c(comp[2], comp[1]))
  meta_sub$Gender <- factor(meta_sub$Gender)
  sample_data(ps_sub) <- sample_data(meta_sub)
  comp_label <- paste0("E", comp[1], "vE", comp[2])
  formula_text <- paste("~ ET +", paste(covariates, collapse = " + "))
  cat(sprintf("LinDA %-10s %-11s %-13s n=%d\n",
              data_label, model_label, comp_label, nrow(meta_sub)))

  otu <- as(otu_table(ps_sub), "matrix")
  if (!taxa_are_rows(ps_sub)) otu <- t(otu)
  res <- tryCatch(
    linda(
      feature.dat = otu, meta.dat = meta_sub,
      formula = formula_text, feature.dat.type = "count",
      is.winsor = TRUE, outlier.pct = 0.03, p.adj.method = "BH",
      alpha = 0.05, prev.filter = 0.1, mean.abund.filter = 0
    ),
    error = function(e) {
      record_failure("add", "differential", model_label, comp_label,
                     data_label, conditionMessage(e))
      NULL
    }
  )
  if (is.null(res)) return(NULL)
  et_name <- paste0("ET", comp[1])
  if (!et_name %in% names(res$output)) et_name <- names(res$output)[1]
  df <- as.data.frame(res$output[[et_name]])
  df$feature <- rownames(df)
  df %>% transmute(
    feature, log2FoldChange = as.numeric(log2FoldChange),
    pvalue = as.numeric(pvalue), padj = as.numeric(padj),
    comparison = comp_label, data_type = data_label,
    model = model_label, n = nrow(meta_sub),
    sig = !is.na(padj) & padj < 0.05
  )
}

run_limma_pairwise <- function(abund, meta, comp, covariates,
                               eligible_ids, model_label) {
  keep <- meta$SampleID %in% eligible_ids & meta$ET %in% comp
  meta_sub <- meta[keep, , drop = FALSE]
  abund_sub <- abund[, meta_sub$SampleID, drop = FALSE]
  min_nonzero <- min(abund_sub[abund_sub > 0], na.rm = TRUE)
  abund_log <- log2(abund_sub + min_nonzero / 2)
  model_df <- data.frame(
    group = factor(meta_sub$ET, levels = c(comp[2], comp[1])),
    Gender = factor(meta_sub$Gender), Age = as.numeric(meta_sub$Age),
    BMI = as.numeric(meta_sub$BMI), Cr = as.numeric(meta_sub$Cr)
  )
  design <- model.matrix(reformulate(c("group", covariates)), data = model_df)
  coef_name <- grep("^group", colnames(design), value = TRUE)[1]
  comp_label <- paste0("E", comp[1], "vE", comp[2])
  cat(sprintf("limma %-10s %-11s %-13s n=%d\n",
              "Metabolome", model_label, comp_label, nrow(meta_sub)))
  fit <- eBayes(lmFit(abund_log, design))
  res <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  res$feature <- rownames(res)
  res %>% transmute(
    feature, log2FoldChange = as.numeric(logFC),
    pvalue = as.numeric(P.Value), padj = as.numeric(adj.P.Val),
    comparison = comp_label, data_type = "Metabolome",
    model = model_label, n = nrow(meta_sub),
    sig = !is.na(padj) & padj < 0.05
  )
}

run_da_model <- function(model_label, covariates, eligible_ids) {
  cat("\n=== Differential model: ", model_label, " ===\n", sep = "")
  list(
    Taxa = bind_rows(lapply(comparisons, function(comp)
      run_linda_pairwise(ps_taxa, comp, "Taxa", covariates,
                         eligible_ids, model_label))),
    Pathway = bind_rows(lapply(comparisons, function(comp)
      run_linda_pairwise(ps_path, comp, "Pathway", covariates,
                         eligible_ids, model_label))),
    EC = bind_rows(lapply(comparisons, function(comp)
      run_linda_pairwise(ps_ec, comp, "EC", covariates,
                         eligible_ids, model_label))),
    Metabolome = bind_rows(lapply(comparisons, function(comp)
      run_limma_pairwise(metab_abund, metab_meta, comp, covariates,
                         eligible_ids, model_label)))
  )
}

extract_core <- function(results_df) {
  core_for <- function(prefix) {
    c1 <- paste0(prefix, "vE2")
    c2 <- paste0(prefix, "vE3")
    a <- results_df %>% filter(comparison == c1, padj < 0.05) %>%
      select(feature, lfc1 = log2FoldChange)
    b <- results_df %>% filter(comparison == c2, padj < 0.05) %>%
      select(feature, lfc2 = log2FoldChange)
    inner_join(a, b, by = "feature") %>%
      filter(sign(lfc1) == sign(lfc2)) %>% pull(feature)
  }
  list(E1 = core_for("E1"), E4 = core_for("E4"))
}

# ------------------------------------------------------------------------------
# Cross-omics frozen feature sets and engines
# ------------------------------------------------------------------------------

curated_metab <- tribble(
  ~feature_id, ~display_name,
  "NEG0024_p-Cresol", "p-Cresol",
  "NEG0216_p-Cresol sulfate", "p-Cresol sulfate",
  "NEG0046_Indoxyl sulfate", "Indoxyl sulfate",
  "NEG0144_5-Methoxyindoleacetate", "5-MeOIA",
  "NEG0147_Indole-3-propionic acid", "IPA",
  "NEG0145_Phenylacetylglutamine", "PAGln",
  "NEG0015_Deoxycholic acid", "DCA",
  "NEG0082_Taurocholic acid", "TCA",
  "POS0156_Trimethylamine N-oxide", "TMAO",
  "POS0006_Betaine", "Betaine",
  "POS0002_Adenosine", "Adenosine",
  "POS0003_Adenine", "Adenine"
)

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

prepare_crossomics <- function() {
  e1_core <- read.csv(here("results", "tables", "section_3_4", "taxa_E1_core.csv"),
                      stringsAsFactors = FALSE)
  e4_core <- read.csv(here("results", "tables", "section_3_4", "taxa_E4_core.csv"),
                      stringsAsFactors = FALSE)
  shared <- read.csv(here("results", "tables", "section_3_4", "taxa_Shared.csv"),
                     stringsAsFactors = FALSE)
  e1_core <- e1_core %>% mutate(
    mean_lfc = rowMeans(cbind(log2FoldChange_E1vE2, log2FoldChange_E1vE3), na.rm = TRUE),
    abs_lfc = abs(mean_lfc), direction = ifelse(mean_lfc > 0, "Higher", "Lower")
  )
  e1_only <- e1_core %>% filter(!feature %in% shared$feature)
  selected_species <- bind_rows(
    e1_only %>% filter(direction == "Higher") %>% arrange(desc(abs_lfc)) %>% slice_head(n = 10),
    e1_only %>% filter(direction == "Lower") %>% arrange(desc(abs_lfc)) %>% slice_head(n = 10),
    e4_core
  ) %>% distinct(feature, .keep_all = TRUE) %>% pull(feature)

  ec_shared <- read.csv(here("results", "tables", "section_3_4", "ec_shared_direction.csv"),
                        stringsAsFactors = FALSE)
  convergent_ec <- ec_shared %>% filter(direction_type == "Convergent") %>% pull(feature)
  selected_function <- c(e1_pathways, e4_pathways, convergent_ec)

  sid_col <- intersect(c("SampleID", "sampleid", "Sample_ID"), names(clinical))[1]
  if (is.na(sid_col)) sid_col <- names(clinical)[1]
  all_ids <- Reduce(intersect, list(
    sample_names(ps_taxa), sample_names(ps_path), sample_names(ps_ec),
    metab_raw$metadata$SampleID, clinical[[sid_col]]
  ))

  taxa_sub <- prune_samples(all_ids, ps_taxa)
  taxa_rel <- transform_sample_counts(taxa_sub, function(x) x / sum(x))
  taxa_mat <- as(otu_table(taxa_rel), "matrix")
  if (!taxa_are_rows(taxa_rel)) taxa_mat <- t(taxa_mat)

  path_sub <- prune_samples(all_ids, ps_path)
  path_mat0 <- as(otu_table(path_sub), "matrix")
  if (!taxa_are_rows(path_sub)) path_mat0 <- t(path_mat0)
  path_cpm <- sweep(path_mat0, 2, colSums(path_mat0), "/") * 1e6
  path_log <- log2(path_cpm + 1)

  ec_sub <- prune_samples(all_ids, ps_ec)
  ec_mat0 <- as(otu_table(ec_sub), "matrix")
  if (!taxa_are_rows(ec_sub)) ec_mat0 <- t(ec_mat0)
  ec_cpm <- sweep(ec_mat0, 2, colSums(ec_mat0), "/") * 1e6
  ec_log <- log2(ec_cpm + 1)

  if (!all(selected_species %in% rownames(taxa_mat)))
    stop("Frozen species set is not fully available")
  if (!all(c(e1_pathways, e4_pathways) %in% rownames(path_log)))
    stop("Frozen pathway set is not fully available")
  if (!all(convergent_ec %in% rownames(ec_log)))
    stop("Frozen convergent-EC set is not fully available")

  metabolite_ids <- curated_metab$feature_id
  if (!all(metabolite_ids %in% rownames(metab_abund)))
    stop("Frozen metabolite target set is not fully available")
  clin_sub <- clinical[match(all_ids, clinical[[sid_col]]), , drop = FALSE]
  ua <- as.numeric(clin_sub$UA)
  names(ua) <- all_ids

  list(
    all_ids = all_ids,
    species = taxa_mat[selected_species, all_ids, drop = FALSE],
    function_mat = rbind(
      path_log[c(e1_pathways, e4_pathways), all_ids, drop = FALSE],
      ec_log[convergent_ec, all_ids, drop = FALSE]
    ),
    targets = rbind(metab_abund[metabolite_ids, all_ids, drop = FALSE],
                    UA = ua[all_ids]),
    sid_col = sid_col
  )
}

make_covar_matrix <- function(ids, covariates) {
  m <- meta_master[match(ids, meta_master$SampleID), , drop = FALSE]
  parts <- list()
  if ("Gender" %in% covariates) parts$Gender <- as.numeric(factor(m$Gender))
  if ("Age" %in% covariates) parts$Age <- as.numeric(m$Age)
  if ("BMI" %in% covariates) parts$BMI <- as.numeric(m$BMI)
  if ("Cr" %in% covariates) parts$Cr <- as.numeric(m$Cr)
  out <- do.call(cbind, parts)
  rownames(out) <- ids
  out
}

compute_partial_cor <- function(feature_mat, target_mat, covar_mat, label) {
  ids <- colnames(feature_mat)
  stopifnot(identical(ids, colnames(target_mat)), identical(ids, rownames(covar_mat)))
  if (anyNA(covar_mat) || anyNA(feature_mat) || anyNA(target_mat))
    stop("Missing values entered partial-correlation engine: ", label)
  n <- length(ids); q <- ncol(covar_mat)
  x_design <- cbind(`(Intercept)` = 1, covar_mat)
  rank_residual <- function(v) lm.fit(x_design, rank(v, ties.method = "average"))$residuals
  feat_resid <- t(apply(feature_mat, 1, rank_residual))
  tgt_resid <- t(apply(target_mat, 1, rank_residual))
  rho <- cor(t(feat_resid), t(tgt_resid), method = "pearson")
  tstat <- rho * sqrt(n - 2 - q) / sqrt(1 - rho^2)
  pval <- 2 * pt(abs(tstat), df = n - 2 - q, lower.tail = FALSE)
  padj <- matrix(p.adjust(as.vector(pval), method = "BH"),
                 nrow = nrow(pval), dimnames = dimnames(pval))
  cat(sprintf("Partial correlation %-22s n=%d q=%d sig=%d/%d\n",
              label, n, q, sum(padj < 0.05), length(padj)))
  list(rho = rho, pval = pval, padj = padj, n = n, q = q)
}

run_dbrda <- function(ids, covariates, model_label) {
  ps_sub <- prune_samples(ids, ps_taxa)
  ids2 <- sample_names(ps_sub)
  dist_bc <- distance(ps_sub, method = "bray")
  abund <- metab_abund[, ids2, drop = FALSE]
  min_nonzero <- min(abund[abund > 0], na.rm = TRUE)
  metab_log <- log2(abund + min_nonzero / 2)
  pca <- prcomp(t(metab_log), center = TRUE, scale. = TRUE)
  cum_var <- cumsum(pca$sdev^2 / sum(pca$sdev^2))
  n_pcs <- which(cum_var >= 0.80)[1]
  scores_pc <- as.data.frame(pca$x[ids2, seq_len(n_pcs), drop = FALSE])
  cov_df <- meta_master[match(ids2, meta_master$SampleID), , drop = FALSE]
  pred <- cbind(scores_pc, data.frame(
    Sex = factor(cov_df$Gender), Age = cov_df$Age,
    BMI = cov_df$BMI, Cr = cov_df$Cr, row.names = ids2
  ))
  cond_names <- unname(c(Gender = "Sex", Age = "Age", BMI = "BMI", Cr = "Cr")[
    covariates
  ])
  form <- as.formula(paste(
    "dist_bc ~", paste(colnames(scores_pc), collapse = " + "),
    "+ Condition(", paste(cond_names, collapse = " + "), ")"
  ))
  set.seed(42)
  fit <- capscale(form, data = pred)
  set.seed(42)
  perm <- anova(fit, permutations = 999)
  total <- fit$tot.chi
  raw <- fit$CCA$tot.chi / total
  adj <- 1 - (1 - raw) * (length(ids2) - 1) /
    (length(ids2) - n_pcs - 1)
  tibble(
    model = model_label, n = length(ids2), n_metabolome_PCs = n_pcs,
    conditioned_covariates = paste(covariates, collapse = "+"),
    raw_R2 = raw, adjusted_R2 = adj,
    permutation_F = as.numeric(perm$F[1]),
    permutation_P = as.numeric(perm$`Pr(>F)`[1])
  ) -> summary
  cat(sprintf("db-RDA %-12s n=%d PCs=%d adjR2=%.6f P=%.4f\n",
              model_label, length(ids2), n_pcs, adj, summary$permutation_P))
  list(fit = fit, summary = summary)
}

# ------------------------------------------------------------------------------
# Direct HUA engine (22 frozen candidates)
# ------------------------------------------------------------------------------

mediator_cols <- c(
  "MeOIA", "pCresol", "pCresolSulfate", "IndoxylSulfate", "IPA", "PAGln",
  "DCA", "TMAO", "Adenosine", "PWY5695", "PWY0_1296", "PWY5004",
  "PWY7219", "PWY6609", "ARGORNPROST", "PWY7234", "EC_2_6_1_1",
  "EC_2_3_1_129", "EC_4_2_1_8", "EC_2_7_6_3", "EC_3_1_3_15",
  "EC_2_1_1_13"
)
mediator_labels <- c(
  MeOIA = "5-Methoxyindoleacetate", pCresol = "p-Cresol",
  pCresolSulfate = "p-Cresol sulfate", IndoxylSulfate = "Indoxyl sulfate",
  IPA = "Indole-3-propionic acid", PAGln = "Phenylacetylglutamine",
  DCA = "Deoxycholic acid", TMAO = "TMAO", Adenosine = "Adenosine",
  PWY5695 = "PWY-5695 (urate biosynthesis)",
  PWY0_1296 = "PWY0-1296 (purine ribonucleoside degradation)",
  PWY5004 = "PWY-5004 (citrulline metabolism)",
  PWY7219 = "PWY-7219 (adenosine ribonucleotides de novo)",
  PWY6609 = "PWY-6609 (adenine/adenosine salvage)",
  ARGORNPROST = "ARGORNPROST-PWY (Arg/Orn/Pro)",
  PWY7234 = "PWY-7234 (IMP biosynthesis)",
  EC_2_6_1_1 = "EC 2.6.1.1 (aspartate transaminase)",
  EC_2_3_1_129 = "EC 2.3.1.129 (UDP-GlcNAc O-acyltransferase)",
  EC_4_2_1_8 = "EC 4.2.1.8", EC_2_7_6_3 = "EC 2.7.6.3 (diphosphokinase)",
  EC_3_1_3_15 = "EC 3.1.3.15 (histidinol-phosphatase)",
  EC_2_1_1_13 = "EC 2.1.1.13 (Met synthase)"
)
mediator_layer <- setNames(c(rep("Metabolite", 9), rep("Pathway", 7), rep("EC", 6)),
                           mediator_cols)

idx_med <- match(med_df$SampleID, meta_master$SampleID)
med_df$BMI <- meta_master$BMI[idx_med]
med_df$Cr <- meta_master$Cr[idx_med]

run_hua_model <- function(data, feature, covariates, eligible_ids, model_label) {
  sub <- data[data$SampleID %in% eligible_ids,
              c("SampleID", "HUA", "Enterotype", covariates, feature), drop = FALSE]
  sub <- sub[complete.cases(sub), , drop = FALSE]
  sub$HUA <- as.numeric(sub$HUA)
  sub$Enterotype <- factor(sub$Enterotype)
  sub$Gender <- factor(sub$Gender)
  form <- as.formula(paste("HUA ~", feature, "+ Enterotype +",
                           paste(covariates, collapse = " + ")))
  warn <- character()
  fit <- tryCatch(
    withCallingHandlers(
      glm(form, data = sub, family = binomial("logit")),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w)); invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      record_failure("add", "direct_HUA", model_label, NA_character_,
                     feature, conditionMessage(e)); NULL
    }
  )
  if (is.null(fit)) return(NULL)
  if (length(warn) > 0) {
    record_failure("add", "direct_HUA_warning", model_label, NA_character_,
                   feature, paste(unique(warn), collapse = " | "))
  }
  co <- summary(fit)$coefficients
  if (!feature %in% rownames(co)) {
    record_failure("add", "direct_HUA", model_label, NA_character_,
                   feature, "Feature coefficient missing")
    return(NULL)
  }
  beta <- co[feature, "Estimate"]
  se <- co[feature, "Std. Error"]
  tibble(
    feature = feature, label = unname(mediator_labels[feature]),
    layer = unname(mediator_layer[feature]), model = model_label,
    adjustment = paste(covariates, collapse = "+"), n = nrow(sub),
    events = sum(sub$HUA == 1), beta = beta, SE = se,
    OR = exp(beta), OR_CI_lo = exp(beta - 1.96 * se),
    OR_CI_hi = exp(beta + 1.96 * se),
    z = co[feature, "z value"], pvalue = co[feature, "Pr(>|z|)"]
  )
}

run_hua_set <- function(model_label, covariates, eligible_ids) {
  out <- bind_rows(lapply(mediator_cols, function(feature)
    run_hua_model(med_df, feature, covariates, eligible_ids, model_label)))
  out$padj <- p.adjust(out$pvalue, method = "BH")
  out
}

# ==============================================================================
# BASELINE PHASE — must pass before extended results are interpreted
# ==============================================================================

cat("\n================ BASELINE REPRODUCTION ================\n")

da_base_full <- run_da_model("base_full", c("Gender", "Age"), base_ids)
da_reference_paths <- c(
  Taxa = input_paths[["da_taxa_reference"]],
  Pathway = input_paths[["da_pathway_reference"]],
  EC = input_paths[["da_ec_reference"]],
  Metabolome = input_paths[["da_metabolome_reference"]]
)
for (layer in names(da_base_full)) {
  ref <- read.csv(da_reference_paths[[layer]], stringsAsFactors = FALSE)
  compare_tidy_baseline(
    da_base_full[[layer]], ref,
    paste0("differential_", layer), c("feature", "comparison"),
    "log2FoldChange", "padj", tolerance = 1e-8
  )
}

cross <- prepare_crossomics()
cross_base_ids <- cross$all_ids[
  complete.cases(meta_master[match(cross$all_ids, meta_master$SampleID), c("Gender", "Age")]) &
    complete.cases(t(cross$targets[, cross$all_ids, drop = FALSE]))
]
sp_base <- compute_partial_cor(
  cross$species[, cross_base_ids, drop = FALSE],
  cross$targets[, cross_base_ids, drop = FALSE],
  make_covar_matrix(cross_base_ids, c("Gender", "Age")), "species_base_full"
)
fn_base <- compute_partial_cor(
  cross$function_mat[, cross_base_ids, drop = FALSE],
  cross$targets[, cross_base_ids, drop = FALSE],
  make_covar_matrix(cross_base_ids, c("Gender", "Age")), "function_base_full"
)

for (spec in list(
  list(name = "species_rho", current = sp_base$rho,
       reference = read_saved_matrix(input_paths[["species_rho_reference"]],
                                     "species", "species_clean")),
  list(name = "species_padj", current = sp_base$padj,
       reference = read_saved_matrix(input_paths[["species_padj_reference"]],
                                     "species", "species_clean")),
  list(name = "function_rho", current = fn_base$rho,
       reference = read_saved_matrix(input_paths[["function_rho_reference"]], "feature")),
  list(name = "function_padj", current = fn_base$padj,
       reference = read_saved_matrix(input_paths[["function_padj_reference"]], "feature"))
)) {
  keys_match <- identical(dimnames(spec$current), dimnames(spec$reference))
  diff <- if (keys_match) safe_max_diff(as.vector(spec$current), as.vector(spec$reference)) else Inf
  add_baseline_check(spec$name, keys_match, length(spec$current), length(spec$reference),
                     diff, 0, tolerance = 1e-8)
}

dbrda_base <- run_dbrda(cross_base_ids, c("Gender", "Age"), "base_full")
dbrda_ref <- readRDS(input_paths[["dbrda_reference"]])
dbrda_raw_ref <- dbrda_ref$CCA$tot.chi / dbrda_ref$tot.chi
add_baseline_check(
  "partial_dbRDA", TRUE, dbrda_base$summary$n,
  nrow(scores(dbrda_ref, display = "sites")),
  abs(dbrda_base$summary$raw_R2 - dbrda_raw_ref), 0,
  tolerance = 1e-8, note = "Raw constrained R2 compared with saved model"
)

hua_base_ids <- med_df$SampleID[complete.cases(med_df[, c("HUA", "Enterotype", "Gender", "Age")])]
hua_base <- run_hua_set("base_full", c("Gender", "Age"), hua_base_ids)
hua_ref_all <- read.csv(input_paths[["direct_feature_reference"]], stringsAsFactors = FALSE)
hua_ref <- hua_ref_all %>% filter(model == "LOGISTIC_HUA_base") %>%
  transmute(feature, beta = as.numeric(HUA_beta), padj = as.numeric(HUA_padj))
compare_tidy_baseline(hua_base, hua_ref, "direct_feature_HUA",
                      "feature", "beta", "padj", tolerance = 1e-8)

write.csv(baseline_checks, file.path(out_dir, "baseline_reproduction_checks.csv"),
          row.names = FALSE)
print(baseline_checks, n = Inf)

if (nrow(baseline_checks) == 0 || any(!baseline_checks$pass)) {
  write.csv(record_failure("get"), file.path(out_dir, "failures_and_warnings.csv"),
            row.names = FALSE)
  stop("Baseline reproduction gate FAILED; extended analysis not run.")
}
cat("BASELINE GATE: PASS\n")

# ==============================================================================
# EXTENDED PHASE — matched complete-case base versus BMI+Cr
# ==============================================================================

cat("\n================ MATCHED SENSITIVITY MODELS ================\n")

da_base_cc <- run_da_model("base_cc", c("Gender", "Age"), cc_ids)
da_ext_cc <- run_da_model("extended_cc", c("Gender", "Age", "BMI", "Cr"), cc_ids)

da_detail <- list()
da_summary <- list()
for (layer in names(da_base_full)) {
  all_models <- bind_rows(da_base_full[[layer]], da_base_cc[[layer]], da_ext_cc[[layer]])
  write.csv(all_models,
            file.path(out_dir, paste0("differential_", tolower(layer), "_all_models.csv")),
            row.names = FALSE)
  cores_full <- extract_core(da_base_full[[layer]])
  cores_base_cc <- extract_core(da_base_cc[[layer]])
  cores_ext <- extract_core(da_ext_cc[[layer]])
  for (et in c("E1", "E4")) {
    comps <- c(paste0(et, "vE2"), paste0(et, "vE3"))
    for (feature in cores_full[[et]]) {
      fetch <- function(x, comp, col) {
        z <- x %>% filter(feature == !!feature, comparison == comp)
        if (nrow(z) != 1) return(NA_real_)
        as.numeric(z[[col]][1])
      }
      row <- tibble(
        layer = layer, enterotype = et, feature = feature,
        comparison_1 = comps[1], comparison_2 = comps[2],
        original_lfc_1 = fetch(da_base_full[[layer]], comps[1], "log2FoldChange"),
        original_lfc_2 = fetch(da_base_full[[layer]], comps[2], "log2FoldChange"),
        base_cc_lfc_1 = fetch(da_base_cc[[layer]], comps[1], "log2FoldChange"),
        base_cc_lfc_2 = fetch(da_base_cc[[layer]], comps[2], "log2FoldChange"),
        ext_lfc_1 = fetch(da_ext_cc[[layer]], comps[1], "log2FoldChange"),
        ext_lfc_2 = fetch(da_ext_cc[[layer]], comps[2], "log2FoldChange"),
        ext_padj_1 = fetch(da_ext_cc[[layer]], comps[1], "padj"),
        ext_padj_2 = fetch(da_ext_cc[[layer]], comps[2], "padj")
      ) %>% mutate(
        direction_retained_both =
          sign(original_lfc_1) == sign(ext_lfc_1) &
          sign(original_lfc_2) == sign(ext_lfc_2),
        base_cc_core = feature %in% cores_base_cc[[et]],
        extended_core = feature %in% cores_ext[[et]]
      )
      da_detail[[length(da_detail) + 1]] <- row
    }
    n_original <- length(cores_full[[et]])
    n_dir <- sum(vapply(da_detail, function(z) {
      nrow(z) == 1 && z$layer == layer && z$enterotype == et &&
        isTRUE(z$direction_retained_both)
    }, logical(1)))
    da_summary[[length(da_summary) + 1]] <- tibble(
      layer = layer, enterotype = et, original_core_n = n_original,
      direction_retained_both_n = n_dir,
      direction_retained_both_rate = ifelse(n_original > 0, n_dir / n_original, NA_real_),
      base_cc_core_n = length(intersect(cores_full[[et]], cores_base_cc[[et]])),
      extended_core_retained_n = length(intersect(cores_full[[et]], cores_ext[[et]])),
      extended_core_retained_rate = ifelse(
        n_original > 0, length(intersect(cores_full[[et]], cores_ext[[et]])) / n_original, NA_real_
      ),
      extended_new_core_n = length(setdiff(cores_ext[[et]], cores_full[[et]]))
    )
  }
}
da_detail <- bind_rows(da_detail)
da_summary <- bind_rows(da_summary)
write.csv(da_detail, file.path(out_dir, "differential_original_core_stability.csv"),
          row.names = FALSE)
write.csv(da_summary, file.path(out_dir, "differential_core_stability_summary.csv"),
          row.names = FALSE)

cross_cc_ids <- cross$all_ids[
  complete.cases(meta_master[match(cross$all_ids, meta_master$SampleID),
                                  c("Gender", "Age", "BMI", "Cr")]) &
    complete.cases(t(cross$targets[, cross$all_ids, drop = FALSE]))
]
sp_base_cc <- compute_partial_cor(
  cross$species[, cross_cc_ids, drop = FALSE], cross$targets[, cross_cc_ids, drop = FALSE],
  make_covar_matrix(cross_cc_ids, c("Gender", "Age")), "species_base_cc"
)
sp_ext <- compute_partial_cor(
  cross$species[, cross_cc_ids, drop = FALSE], cross$targets[, cross_cc_ids, drop = FALSE],
  make_covar_matrix(cross_cc_ids, c("Gender", "Age", "BMI", "Cr")), "species_extended"
)
fn_base_cc <- compute_partial_cor(
  cross$function_mat[, cross_cc_ids, drop = FALSE], cross$targets[, cross_cc_ids, drop = FALSE],
  make_covar_matrix(cross_cc_ids, c("Gender", "Age")), "function_base_cc"
)
fn_ext <- compute_partial_cor(
  cross$function_mat[, cross_cc_ids, drop = FALSE], cross$targets[, cross_cc_ids, drop = FALSE],
  make_covar_matrix(cross_cc_ids, c("Gender", "Age", "BMI", "Cr")), "function_extended"
)

build_pair_sensitivity <- function(orig, base_cc, ext, family) {
  matrix_to_tidy(orig$rho, value_name = "original_rho") %>%
    left_join(matrix_to_tidy(orig$padj, value_name = "original_padj"),
              by = c("feature", "target")) %>%
    left_join(matrix_to_tidy(base_cc$rho, value_name = "base_cc_rho"),
              by = c("feature", "target")) %>%
    left_join(matrix_to_tidy(base_cc$padj, value_name = "base_cc_padj"),
              by = c("feature", "target")) %>%
    left_join(matrix_to_tidy(ext$rho, value_name = "extended_rho"),
              by = c("feature", "target")) %>%
    left_join(matrix_to_tidy(ext$padj, value_name = "extended_padj"),
              by = c("feature", "target")) %>%
    mutate(
      family = family, n_original = orig$n, n_matched = ext$n,
      direction_retained_original_to_extended =
        sign(original_rho) == sign(extended_rho),
      direction_retained_matched = sign(base_cc_rho) == sign(extended_rho),
      base_cc_sig = base_cc_padj < 0.05, extended_sig = extended_padj < 0.05,
      significance_transition = case_when(
        base_cc_sig & extended_sig ~ "retained",
        base_cc_sig & !extended_sig ~ "lost",
        !base_cc_sig & extended_sig ~ "gained",
        TRUE ~ "non-significant_both"
      )
    )
}

sp_pairs <- build_pair_sensitivity(sp_base, sp_base_cc, sp_ext, "Species")
fn_pairs <- build_pair_sensitivity(fn_base, fn_base_cc, fn_ext, "Function")
write.csv(sp_pairs, file.path(out_dir, "crossomics_species_pair_sensitivity.csv"),
          row.names = FALSE)
write.csv(fn_pairs, file.path(out_dir, "crossomics_function_pair_sensitivity.csv"),
          row.names = FALSE)

summarize_pairs <- function(df, scope) {
  z <- if (scope == "crossomics_only") filter(df, target != "UA") else df
  tibble(
    family = unique(z$family), target_scope = scope, pairs_n = nrow(z),
    n_original = unique(z$n_original), n_matched = unique(z$n_matched),
    rho_spearman_base_cc_vs_extended = cor(
      z$base_cc_rho, z$extended_rho, method = "spearman", use = "complete.obs"
    ),
    direction_concordance_base_cc_vs_extended = mean(z$direction_retained_matched),
    direction_concordance_original_vs_extended =
      mean(z$direction_retained_original_to_extended),
    base_cc_FDR_n = sum(z$base_cc_sig), extended_FDR_n = sum(z$extended_sig),
    FDR_retained_n = sum(z$significance_transition == "retained"),
    FDR_lost_n = sum(z$significance_transition == "lost"),
    FDR_gained_n = sum(z$significance_transition == "gained")
  )
}
cross_summary <- bind_rows(
  summarize_pairs(sp_pairs, "crossomics_only"),
  summarize_pairs(sp_pairs, "including_UA"),
  summarize_pairs(fn_pairs, "crossomics_only"),
  summarize_pairs(fn_pairs, "including_UA")
)
write.csv(cross_summary, file.path(out_dir, "crossomics_correlation_summary.csv"),
          row.names = FALSE)

dbrda_base_cc <- run_dbrda(cross_cc_ids, c("Gender", "Age"), "base_cc")
dbrda_ext <- run_dbrda(cross_cc_ids, c("Gender", "Age", "BMI", "Cr"), "extended_cc")
dbrda_summary <- bind_rows(
  dbrda_base$summary, dbrda_base_cc$summary, dbrda_ext$summary
)
write.csv(dbrda_summary, file.path(out_dir, "crossomics_dbrda_summary.csv"),
          row.names = FALSE)
saveRDS(dbrda_base_cc$fit, file.path(out_dir, "dbrda_base_cc.rds"))
saveRDS(dbrda_ext$fit, file.path(out_dir, "dbrda_extended_cc.rds"))

hua_cc_ids <- med_df$SampleID[complete.cases(
  med_df[, c("HUA", "Enterotype", "Gender", "Age", "BMI", "Cr")]
)]
hua_base_cc <- run_hua_set("base_cc", c("Gender", "Age"), hua_cc_ids)
hua_ext <- run_hua_set("extended_cc", c("Gender", "Age", "BMI", "Cr"), hua_cc_ids)
hua_all <- bind_rows(hua_base, hua_base_cc, hua_ext)
write.csv(hua_all, file.path(out_dir, "direct_feature_HUA_all_models.csv"),
          row.names = FALSE)

hua_sensitivity <- hua_base %>%
  select(feature, label, layer, original_n = n, original_events = events,
         original_beta = beta, original_OR = OR, original_padj = padj) %>%
  left_join(hua_base_cc %>%
    select(feature, base_cc_n = n, base_cc_events = events,
           base_cc_beta = beta, base_cc_OR = OR, base_cc_padj = padj), by = "feature") %>%
  left_join(hua_ext %>%
    select(feature, extended_n = n, extended_events = events,
           extended_beta = beta, extended_OR = OR,
           extended_OR_CI_lo = OR_CI_lo, extended_OR_CI_hi = OR_CI_hi,
           extended_pvalue = pvalue, extended_padj = padj), by = "feature") %>%
  mutate(
    direction_retained_matched = sign(base_cc_beta) == sign(extended_beta),
    direction_retained_original = sign(original_beta) == sign(extended_beta),
    significance_transition = case_when(
      base_cc_padj < 0.05 & extended_padj < 0.05 ~ "retained",
      base_cc_padj < 0.05 & extended_padj >= 0.05 ~ "lost",
      base_cc_padj >= 0.05 & extended_padj < 0.05 ~ "gained",
      TRUE ~ "non-significant_both"
    )
  )
write.csv(hua_sensitivity, file.path(out_dir, "direct_feature_HUA_sensitivity.csv"),
          row.names = FALSE)

hua_summary <- tibble(
  candidates_n = nrow(hua_sensitivity),
  matched_n = unique(hua_sensitivity$extended_n),
  matched_events = unique(hua_sensitivity$extended_events),
  direction_concordance_base_cc_vs_extended =
    mean(hua_sensitivity$direction_retained_matched),
  direction_concordance_original_vs_extended =
    mean(hua_sensitivity$direction_retained_original),
  base_cc_FDR_n = sum(hua_sensitivity$base_cc_padj < 0.05),
  extended_FDR_n = sum(hua_sensitivity$extended_padj < 0.05),
  FDR_retained_n = sum(hua_sensitivity$significance_transition == "retained"),
  FDR_lost_n = sum(hua_sensitivity$significance_transition == "lost"),
  FDR_gained_n = sum(hua_sensitivity$significance_transition == "gained")
)
write.csv(hua_summary, file.path(out_dir, "direct_feature_HUA_summary.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# Complete reporting, decision-rule evaluation, and manifests
# ------------------------------------------------------------------------------

key_features <- c("PWY0_1296", "EC_2_6_1_1", "MeOIA")
key_hua <- hua_sensitivity %>% filter(feature %in% key_features)
h1_min <- min(da_summary$direction_retained_both_rate, na.rm = TRUE)
h2_species <- cross_summary %>%
  filter(family == "Species", target_scope == "crossomics_only")
h2_function <- cross_summary %>%
  filter(family == "Function", target_scope == "crossomics_only")

decision_summary <- tibble(
  hypothesis = c("H1_differential", "H2_crossomics", "H3_direct_HUA"),
  status = c(
    ifelse(h1_min >= 0.80, "supportive",
           ifelse(h1_min < 0.70, "refuting", "ambiguous")),
    ifelse(
      min(h2_species$rho_spearman_base_cc_vs_extended,
          h2_function$rho_spearman_base_cc_vs_extended) >= 0.90 &&
        min(h2_species$direction_concordance_base_cc_vs_extended,
            h2_function$direction_concordance_base_cc_vs_extended) >= 0.90 &&
        dbrda_ext$summary$permutation_P < 0.05,
      "supportive",
      ifelse(
        min(h2_species$rho_spearman_base_cc_vs_extended,
            h2_function$rho_spearman_base_cc_vs_extended) < 0.70 ||
          dbrda_ext$summary$permutation_P >= 0.05,
        "refuting", "ambiguous"
      )
    ),
    ifelse(
      all(key_hua$direction_retained_matched) && sum(hua_ext$padj < 0.05) >= 1,
      "supportive",
      ifelse(sum(!key_hua$direction_retained_matched) >= 2,
             "refuting", "ambiguous")
    )
  ),
  evidence = c(
    sprintf("Minimum layer-enterotype direction-retention rate: %.3f", h1_min),
    sprintf("Species/function rho concordance: %.3f/%.3f; db-RDA P=%.4f",
            h2_species$rho_spearman_base_cc_vs_extended,
            h2_function$rho_spearman_base_cc_vs_extended,
            dbrda_ext$summary$permutation_P),
    sprintf("Key directions retained: %d/3; extended HUA FDR signals: %d",
            sum(key_hua$direction_retained_matched), sum(hua_ext$padj < 0.05))
  )
)
write.csv(decision_summary, file.path(out_dir, "predefined_decision_summary.csv"),
          row.names = FALSE)

failures <- record_failure("get")
write.csv(failures, file.path(out_dir, "failures_and_warnings.csv"), row.names = FALSE)
write_session_info()

output_paths <- list.files(out_dir, full.names = TRUE)
output_paths <- output_paths[!basename(output_paths) %in%
  c("output_manifest.csv", "covariate_sensitivity_run.log")]
output_manifest <- tibble(
  file = basename(output_paths),
  bytes = as.numeric(file.info(output_paths)$size),
  sha256 = vapply(output_paths, sha256_file, character(1))
)
write.csv(output_manifest, file.path(out_dir, "output_manifest.csv"), row.names = FALSE)

cat("\n================ FINAL SUMMARIES ================\n")
cat("\nDifferential core stability:\n"); print(da_summary, n = Inf)
cat("\nCross-omics correlation stability:\n"); print(cross_summary, n = Inf)
cat("\nPartial db-RDA:\n"); print(dbrda_summary, n = Inf)
cat("\nDirect HUA sensitivity:\n"); print(hua_summary)
cat("\nPredefined decisions:\n"); print(decision_summary, n = Inf)
cat("\nFailures/warnings recorded: ", nrow(failures), "\n", sep = "")
cat("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"), "\n", sep = "")
cat("ALL ANALYSES COMPLETE\n")
