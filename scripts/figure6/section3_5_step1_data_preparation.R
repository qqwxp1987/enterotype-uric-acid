# ==============================================================================
# Section 3.5 — Step 1: Data Preparation (Revised)
#
# Purpose:
#   Read all required data files, match samples across omics layers,
#   build a unified analysis data frame with 22 candidate mediators,
#   treatment variables, outcome (UA), and covariates (Sex, Age).
#
# Input files:
#   - phyloseq_Taxa_with_enterotype_seed10.rds  → phenotype (authoritative)
#   - phyloseq_Pathway_with_enterotype_20250210.rds  → pathway abundance
#   - phyloseq_EC_with_enterotype_20250210.rds  → EC enzyme abundance
#   - metabolites_normalized_multiIS_and_metadata_20250904.rds  → metabolite abundance
#
# Output:
#   - results/tables/section_3_5/mediation_analysis_data.rds
#     A data frame: one row per sample, columns = treatment, mediators, outcome, covariates
#
# Candidate mediators (22):
#   Metabolites (9):
#     M1  MeOIA           5-Methoxyindoleacetate (NEG0144)
#     M2  pCresol         p-Cresol (NEG0024)
#     M3  pCresolSulfate  p-Cresol sulfate (NEG0216)
#     M4  IndoxylSulfate  Indoxyl sulfate (NEG0046)
#     M5  IPA             Indole-3-propionic acid (NEG0147)
#     M6  PAGln           Phenylacetylglutamine (NEG0145)
#     M7  DCA             Deoxycholic acid (NEG0015)
#     M8  TMAO            Trimethylamine N-oxide (POS0156)
#     M9  Adenosine       Adenosine (POS0002)
#   Pathways (7):
#     M10 PWY5695         urate biosynthesis / IMP degradation
#     M11 PWY0_1296       purine ribonucleosides degradation
#     M12 PWY5004         superpathway of L-citrulline metabolism
#     M13 PWY7219         adenosine ribonucleotides de novo biosynthesis
#     M14 PWY6609         adenine and adenosine salvage III
#     M15 ARGORNPROST     arginine, ornithine and proline interconversion
#     M16 PWY7234         inosine-5'-phosphate biosynthesis III
#   Convergent shared ECs (6):
#     M17 EC_2_6_1_1      aspartate transaminase
#     M18 EC_2_3_1_129    UDP-GlcNAc O-acyltransferase
#     M19 EC_4_2_1_8
#     M20 EC_2_7_6_3      2-amino-4-hydroxy-6-hydroxymethyldihydropteridine diphosphokinase
#     M21 EC_3_1_3_15     histidinol-phosphatase
#     M22 EC_2_1_1_13     5-methyltetrahydrofolate--homocysteine methyltransferase
#
# Selection criteria:
#   Metabolites: core metabolites from Section 3.3 representing AAA catabolites
#     (all 6 individually), bile acid, TMA pathway, and purine nucleoside
#   Pathways: core pathways from Section 3.3 covering purine metabolism,
#     amino acid metabolism, and urate biosynthesis
#   ECs: convergent shared ECs meeting EITHER (A) significant UA partial
#     correlation in Section 3.4 (padj < 0.05), OR (B) top convergent
#     enzyme feeding de novo purine biosynthesis in Section 3.3
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(here)

cat("========================================\n")
cat("Section 3.5 Step 1: Data Preparation (Revised — 22 candidates)\n")
cat("========================================\n\n")

# --- 0. Output directory ---
out_dir <- here("results", "tables", "section_3_5")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. Load phenotype data (authoritative source)
# ==============================================================================

ps_taxa <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta_auth <- as(sample_data(ps_taxa), "data.frame")
meta_auth$SampleID <- rownames(meta_auth)

cat(sprintf("Authoritative phenotype: %d samples\n", nrow(meta_auth)))
cat("Enterotype distribution:\n")
print(table(meta_auth$Dirichlettype_Genus))

# ==============================================================================
# 2. Load pathway abundance
# ==============================================================================

ps_path <- readRDS(here("data", "clean", "phyloseq_Pathway_with_enterotype_20250210.rds"))
path_ids <- sample_names(ps_path)
cat(sprintf("\nPathway data: %d pathways x %d samples\n",
            ntaxa(ps_path), nsamples(ps_path)))

# ==============================================================================
# 3. Load EC enzyme abundance
# ==============================================================================

ps_ec <- readRDS(here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"))
ec_ids <- sample_names(ps_ec)
cat(sprintf("\nEC data: %d ECs x %d samples\n",
            ntaxa(ps_ec), nsamples(ps_ec)))

# ==============================================================================
# 4. Load metabolite abundance
# ==============================================================================

metab_raw <- readRDS(here("data", "clean",
                          "metabolites_normalized_multiIS_and_metadata_20250904.rds"))
metab_abund <- metab_raw$abundance
metab_meta  <- metab_raw$metadata
metab_ids   <- metab_meta$SampleID

cat(sprintf("Metabolite data: %d metabolites x %d samples\n",
            nrow(metab_abund), ncol(metab_abund)))

# ==============================================================================
# 5. Match samples across all three data layers
# ==============================================================================

shared_ids <- Reduce(intersect, list(
  meta_auth$SampleID,
  path_ids,
  ec_ids,
  metab_ids
))
cat(sprintf("\nSample matching:\n"))
cat(sprintf("  Phenotype (taxa):  %d\n", nrow(meta_auth)))
cat(sprintf("  Pathway:           %d\n", length(path_ids)))
cat(sprintf("  EC:                %d\n", length(ec_ids)))
cat(sprintf("  Metabolite:        %d\n", length(metab_ids)))
cat(sprintf("  Shared (4-way):    %d\n", length(shared_ids)))

stopifnot(length(shared_ids) > 0)

# ==============================================================================
# 6. Extract and transform metabolite mediators
# ==============================================================================

# --- 6a. Define candidate metabolite features (9 individual metabolites) ---
metab_features <- c(
  MeOIA          = "NEG0144_5-Methoxyindoleacetate",
  pCresol        = "NEG0024_p-Cresol",
  pCresolSulfate = "NEG0216_p-Cresol sulfate",
  IndoxylSulfate = "NEG0046_Indoxyl sulfate",
  IPA            = "NEG0147_Indole-3-propionic acid",
  PAGln          = "NEG0145_Phenylacetylglutamine",
  DCA            = "NEG0015_Deoxycholic acid",
  TMAO           = "POS0156_Trimethylamine N-oxide",
  Adenosine      = "POS0002_Adenosine"
)

# Verify all features exist in abundance matrix
missing <- setdiff(metab_features, rownames(metab_abund))
if (length(missing) > 0) {
  stop("Missing metabolite features: ", paste(missing, collapse = ", "))
}
cat("\nAll 9 candidate metabolite features found in abundance matrix.\n")

# --- 6b. Extract and log2-transform ---
metab_sub <- metab_abund[metab_features, shared_ids, drop = FALSE]
min_nonzero <- min(metab_sub[metab_sub > 0], na.rm = TRUE)
metab_log <- log2(metab_sub + min_nonzero / 2)

cat(sprintf("Metabolite log2 transform: pseudocount = %.2e (min_nonzero/2)\n",
            min_nonzero / 2))

# ==============================================================================
# 7. Extract and transform pathway mediators
# ==============================================================================

# Define candidate pathway features (7 pathways)
path_features <- c(
  PWY5695     = "PWY-5695",
  PWY0_1296   = "PWY0-1296",
  PWY5004     = "PWY-5004",
  PWY7219     = "PWY-7219",
  PWY6609     = "PWY-6609",
  ARGORNPROST = "ARGORNPROST-PWY",
  PWY7234     = "PWY-7234"
)

# Get OTU table (pathways x samples)
path_otu <- as(otu_table(ps_path), "matrix")
if (!taxa_are_rows(ps_path)) path_otu <- t(path_otu)

# Match pathway IDs (may have full names like "PWY-5695: urate biosynthesis...")
path_rownames <- rownames(path_otu)
path_id_map <- sapply(path_features, function(pwy) {
  matches <- grep(paste0("^", pwy, ":?"), path_rownames, value = TRUE)
  if (length(matches) == 0) matches <- grep(pwy, path_rownames, value = TRUE, fixed = TRUE)
  if (length(matches) == 0) return(NA_character_)
  return(matches[1])
})

cat("\nPathway feature matching:\n")
for (i in seq_along(path_features)) {
  status <- if (is.na(path_id_map[i])) "NOT FOUND" else "OK"
  cat(sprintf("  %s → %s [%s]\n", names(path_features)[i],
              ifelse(is.na(path_id_map[i]), path_features[i], path_id_map[i]),
              status))
}

missing_path <- names(path_id_map)[is.na(path_id_map)]
if (length(missing_path) > 0) {
  warning("Missing pathways: ", paste(missing_path, collapse = ", "),
          "\nThese will be set to NA in the analysis data frame.")
}

# Extract and log2(CPM+1) transform
path_sub <- path_otu[na.omit(path_id_map), shared_ids, drop = FALSE]
col_sums <- colSums(path_otu[, shared_ids])
path_cpm <- sweep(path_sub, 2, col_sums, "/") * 1e6
path_log <- log2(path_cpm + 1)

cat(sprintf("\nPathway log2(CPM+1) transform: %d pathways x %d samples\n",
            nrow(path_log), ncol(path_log)))

# ==============================================================================
# 8. Extract and transform convergent shared EC mediators
# ==============================================================================

# Define candidate convergent shared ECs (6 ECs)
ec_features <- c(
  EC_2_6_1_1   = "2.6.1.1",
  EC_2_3_1_129 = "2.3.1.129",
  EC_4_2_1_8   = "4.2.1.8",
  EC_2_7_6_3   = "2.7.6.3",
  EC_3_1_3_15  = "3.1.3.15",
  EC_2_1_1_13  = "2.1.1.13"
)

# Get OTU table (ECs x samples)
ec_otu <- as(otu_table(ps_ec), "matrix")
if (!taxa_are_rows(ps_ec)) ec_otu <- t(ec_otu)

# Match EC IDs (exact match on short names like "2.6.1.1")
ec_rownames <- rownames(ec_otu)
ec_id_map <- sapply(ec_features, function(ec) {
  idx <- which(ec_rownames == ec)
  if (length(idx) == 1) return(ec_rownames[idx])
  return(NA_character_)
})

cat("\nEC feature matching:\n")
for (i in seq_along(ec_features)) {
  status <- if (is.na(ec_id_map[i])) "NOT FOUND" else "OK"
  cat(sprintf("  %s → %s [%s]\n", names(ec_features)[i],
              ifelse(is.na(ec_id_map[i]), ec_features[i], ec_id_map[i]),
              status))
}

missing_ec <- names(ec_id_map)[is.na(ec_id_map)]
if (length(missing_ec) > 0) {
  warning("Missing ECs: ", paste(missing_ec, collapse = ", "))
}

# Extract and log2(CPM+1) transform
ec_sub <- ec_otu[na.omit(ec_id_map), shared_ids, drop = FALSE]
ec_col_sums <- colSums(ec_otu[, shared_ids])
ec_cpm <- sweep(ec_sub, 2, ec_col_sums, "/") * 1e6
ec_log <- log2(ec_cpm + 1)

cat(sprintf("\nEC log2(CPM+1) transform: %d ECs x %d samples\n",
            nrow(ec_log), ncol(ec_log)))

# ==============================================================================
# 9. Build unified analysis data frame
# ==============================================================================

# Start with phenotype
df <- data.frame(
  SampleID = shared_ids,
  stringsAsFactors = FALSE
)

# --- Enterotype ---
df$Enterotype <- meta_auth[shared_ids, "Dirichlettype_Genus"]

# --- Treatment variables (binary) ---
df$Treat_T1 <- case_when(
  df$Enterotype == "1" ~ 1L,
  df$Enterotype %in% c("2", "3") ~ 0L,
  TRUE ~ NA_integer_
)
df$Treat_T4 <- case_when(
  df$Enterotype == "4" ~ 1L,
  df$Enterotype %in% c("2", "3") ~ 0L,
  TRUE ~ NA_integer_
)

# --- Outcome ---
df$UA <- as.numeric(meta_auth[shared_ids, "UA"])

# --- HUA status (for sensitivity analysis) ---
df$Gender <- meta_auth[shared_ids, "Gender"]
df$Age    <- as.numeric(meta_auth[shared_ids, "Age"])

# HUA definition: UA >= 420 (Male) or >= 360 (Female)
df$HUA <- case_when(
  df$Gender == "Male"   & df$UA >= 420 ~ 1L,
  df$Gender == "Male"   & df$UA <  420 ~ 0L,
  df$Gender == "Female" & df$UA >= 360 ~ 1L,
  df$Gender == "Female" & df$UA <  360 ~ 0L,
  TRUE ~ NA_integer_
)

# --- Covariates ---
df$Sex_numeric <- ifelse(df$Gender == "Male", 1L, 0L)

# --- Metabolite mediators (9) ---
for (m_name in names(metab_features)) {
  df[[m_name]] <- metab_log[metab_features[m_name], shared_ids]
}

# --- Pathway mediators (7) ---
for (pwy_name in names(path_features)) {
  matched_row <- path_id_map[pwy_name]
  if (!is.na(matched_row)) {
    df[[pwy_name]] <- path_log[matched_row, shared_ids]
  } else {
    df[[pwy_name]] <- NA_real_
  }
}

# --- EC mediators (6) ---
for (ec_name in names(ec_features)) {
  matched_row <- ec_id_map[ec_name]
  if (!is.na(matched_row)) {
    df[[ec_name]] <- ec_log[matched_row, shared_ids]
  } else {
    df[[ec_name]] <- NA_real_
  }
}

# ==============================================================================
# 10. Summary and quality checks
# ==============================================================================

cat("\n========================================\n")
cat("Summary of unified analysis data frame\n")
cat("========================================\n\n")

cat(sprintf("Dimensions: %d samples x %d variables\n", nrow(df), ncol(df)))

cat("\n--- Enterotype distribution ---\n")
print(table(df$Enterotype, useNA = "ifany"))

cat("\n--- Treatment variables ---\n")
cat(sprintf("Type 1 model: Treatment=%d, Reference=%d, Excluded=%d\n",
            sum(df$Treat_T1 == 1, na.rm = TRUE),
            sum(df$Treat_T1 == 0, na.rm = TRUE),
            sum(is.na(df$Treat_T1))))
cat(sprintf("Type 4 model: Treatment=%d, Reference=%d, Excluded=%d\n",
            sum(df$Treat_T4 == 1, na.rm = TRUE),
            sum(df$Treat_T4 == 0, na.rm = TRUE),
            sum(is.na(df$Treat_T4))))

cat("\n--- Outcome (UA) ---\n")
print(summary(df$UA))
cat(sprintf("NA count: %d\n", sum(is.na(df$UA))))

cat("\n--- HUA status ---\n")
print(table(df$HUA, useNA = "ifany"))

cat("\n--- Covariates ---\n")
cat("Gender:\n")
print(table(df$Gender))
cat("Age:\n")
print(summary(df$Age))

cat("\n--- Mediator variables (22 candidates) ---\n")
mediator_cols <- c(
  names(metab_features),   # 9 metabolites
  names(path_features),    # 7 pathways
  names(ec_features)       # 6 ECs
)
for (m in mediator_cols) {
  vals <- df[[m]]
  na_n <- sum(is.na(vals))
  if (na_n == nrow(df)) {
    cat(sprintf("  %-18s: ALL NA (feature not found)\n", m))
  } else {
    cat(sprintf("  %-18s: mean=%.2f, sd=%.2f, range=[%.2f, %.2f], NA=%d\n",
                m, mean(vals, na.rm = TRUE), sd(vals, na.rm = TRUE),
                min(vals, na.rm = TRUE), max(vals, na.rm = TRUE), na_n))
  }
}

# --- Check for complete cases in each model ---
t1_vars <- c("Treat_T1", mediator_cols, "UA", "Sex_numeric", "Age")
t1_complete <- sum(complete.cases(df[!is.na(df$Treat_T1), t1_vars]))
cat(sprintf("\nType 1 model complete cases: %d / %d\n",
            t1_complete, sum(!is.na(df$Treat_T1))))

t4_vars <- c("Treat_T4", mediator_cols, "UA", "Sex_numeric", "Age")
t4_complete <- sum(complete.cases(df[!is.na(df$Treat_T4), t4_vars]))
cat(sprintf("Type 4 model complete cases: %d / %d\n",
            t4_complete, sum(!is.na(df$Treat_T4))))

# ==============================================================================
# 11. Save
# ==============================================================================

saveRDS(df, file.path(out_dir, "mediation_analysis_data.rds"))
cat(sprintf("\nAnalysis data frame saved: %s\n",
            file.path(out_dir, "mediation_analysis_data.rds")))

write.csv(df, file.path(out_dir, "mediation_analysis_data.csv"), row.names = FALSE)
cat(sprintf("CSV copy saved: %s\n",
            file.path(out_dir, "mediation_analysis_data.csv")))

cat("\nStep 1 complete.\n")
