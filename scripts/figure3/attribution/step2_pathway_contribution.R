# ======================================================================
#
#  Section 3.7 — Step 2: Anchor Pathway Species Contribution Analysis
#
#  Purpose:
#    For 10 biologically-anchored core pathways (5 E1 + 5 E4), decompose
#    the per-sample species contributions, compare across enterotypes, and
#    produce the statistical tables used for Figure 7 Panel A.
#
#  NOTE on RPK values:
#    pathabundance_stratified.tsv stores HUMAnN3 pathway abundances in
#    RPK (reads per kilobase).  For relative contribution (%) the RPK
#    units cancel exactly within each (SampleID, Pathway) group:
#
#        contribution_ij = RPK_ij / sum_j(RPK_ij) × 100
#
#    So no CPM normalisation is needed here.  For absolute pathway
#    abundance comparisons (which enterotype has higher pathway RPK?),
#    CPM normalisation WOULD be needed — but those comparisons are
#    already addressed by the LinDA log2FC values in Section 3.3 core
#    tables, so we do not repeat them here.
#
#  Inputs:
#    results/tables/section_3_7/prep_data.rds
#    results/tables/section_3_7/pathway_stratified_subset.rds
#
#  Outputs:
#    results/tables/section_3_7/pathway_contribution_data.rds      (long-format, all anchors)
#    results/tables/section_3_7/pathway_contribution_summary.csv   (mean % per enterotype × species)
#    results/tables/section_3_7/pathway_contribution_stats.csv     (KW + Dunn tests, FDR-corrected)
#    results/tables/section_3_7/anchor_pathway_metadata.csv        (anchor definitions + log2FC)
#
# ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(rstatix)
  library(here)
})

out_dir <- here("results", "tables", "section_3_7")

# ======================================================================
# 1. Load prep data bundle
# ======================================================================

cat("=== Step 1: Loading prep data ===\n")
prep   <- readRDS(file.path(out_dir, "prep_data.rds"))
meta   <- prep$metadata          # columns: SampleID, Enterotype (factor 1-4)

# ======================================================================
# 2. Define anchor pathways
# ======================================================================
#
#  Selection criteria (all 3 must hold):
#   (a) Present in the enterotype-specific core tables (pathway_E1_core or
#       pathway_E4_core from Section 3.3).
#   (b) Explicitly discussed in the Section 3.3 manuscript text or named
#       as a candidate mediator in Section 3.5.
#   (c) Best representative of one of the four major functional themes
#       (arginine/TCA, purine catabolism, purine de novo biosynthesis,
#       peptidoglycan biosynthesis, branched-chain/amino-acid biosynthesis).
#
#  Absolute abundance direction (enriched/depleted in that enterotype) is
#  already established by the LinDA log2FC values in Section 3.3 and is
#  noted in the 'direction' column for reference only — it is NOT
#  re-tested here.

anchor_pathways <- tribble(
  ~pathway_id,                ~label_short,                                        ~enterotype, ~functional_theme,          ~direction, ~sec33_notes,
  # ── Type 1 anchors (depleted in E1 relative to E2/E3) ─────────────────────────────────────
  # NOTE: ARGORNPROST-PWY / PWY-5004 / PWY-5505 had only 1 "unclassified" row across all
  # samples in HUMAnN3 output — no species-level attribution available.  Replaced by the
  # next most strongly enriched E1 pathways in the same arginine/ornithine functional cluster.
  "ARGSYNBSUB-PWY",           "L-Arg biosynthesis II (acetyl cycle)",              "E1",        "Arginine metabolism",       "depleted", "Acetyl cycle arm of arginine biosynthesis; E1 core; 148 stratified rows",
  "GLUTORN-PWY",              "L-ornithine biosynthesis",                          "E1",        "Arginine metabolism",       "depleted", "Ornithine is the direct precursor to citrulline/arginine; 139 rows",
  "PWY-7234",                 "IMP biosynthesis III",                              "E1",        "Purine biosynthesis",       "depleted", "Purine de novo synthesis; Section 3.3 discussed",
  "PWY-6608",                 "Guanosine nucleotides degradation III",             "E1",        "Purine catabolism",         "enriched", "Only enriched purine pathway in E1; UA-relevant",
  "ARGSYN-PWY",               "L-Arg biosynthesis I (via L-ornithine)",            "E1",        "Arginine metabolism",       "depleted", "Classic ornithine → citrulline → arginine route; 133 rows",
  # ── Type 4 anchors (enriched in E4 relative to E2/E3) ─────────────────────────────────────
  "PWY-5695",                 "Urate biosynthesis / IMP degradation",              "E4",        "Purine catabolism / UA",    "enriched", "Direct UA link; Section 3.3 & 3.5 central candidate",
  "PWY-6122",                 "5-aminoimidazole ribonucleotide biosynthesis II",   "E4",        "Purine biosynthesis",       "enriched", "Purine de novo synthesis upstream node",
  "PWY-7219",                 "Adenosine ribonucleotides de novo biosynthesis",    "E4",        "Purine biosynthesis",       "enriched", "Section 3.5 mediation candidate (ACME = +5.3)",
  "PEPTIDOGLYCANSYN-PWY",     "Peptidoglycan biosynthesis I",                      "E4",        "Cell wall biosynthesis",    "enriched", "Representative of 4 co-enriched PG pathways in E4",
  "BRANCHED-CHAIN-AA-SYN-PWY","Branched-chain AA biosynthesis",                   "E4",        "Amino acid biosynthesis",   "enriched", "Superpathway for VALSYN/ILEUSYN; Section 3.3 discussed"
)

# Attach Section 3.3 log2FC values for reference (mean of E1vE2 + E1vE3, or E4vE2 + E4vE3)
e1_fc <- prep$pathway_E1_core %>%
  mutate(pathway_id = sub(":.*", "", trimws(feature))) %>%
  mutate(lfc_E1 = (log2FoldChange_E1vE2 + log2FoldChange_E1vE3) / 2) %>%
  select(pathway_id, lfc_E1)

e4_fc <- prep$pathway_E4_core %>%
  mutate(pathway_id = sub(":.*", "", trimws(feature))) %>%
  mutate(lfc_E4 = (log2FoldChange_E4vE2 + log2FoldChange_E4vE3) / 2) %>%
  select(pathway_id, lfc_E4)

anchor_pathways <- anchor_pathways %>%
  left_join(e1_fc, by = "pathway_id") %>%
  left_join(e4_fc, by = "pathway_id")

write.csv(anchor_pathways, file.path(out_dir, "anchor_pathway_metadata.csv"), row.names = FALSE)
cat("  Anchor pathway metadata saved.\n")
cat(sprintf("  %d anchor pathways: %d E1, %d E4\n",
            nrow(anchor_pathways),
            sum(anchor_pathways$enterotype == "E1"),
            sum(anchor_pathways$enterotype == "E4")))


# ======================================================================
# 3. Load stratified pathway subset and reshape
# ======================================================================

cat("\n=== Step 3: Loading stratified pathway subset ===\n")

raw <- readRDS(file.path(out_dir, "pathway_stratified_subset.rds"))
setnames(raw, 1, "Feature_Raw")
cat(sprintf("  %d rows loaded.\n", nrow(raw)))

# HUMAnN3 stratified rows have the form:
#   "PATHWAY_ID: full pathway name|g__Genus.s__species"
# i.e. the species column is separated from the full feature label (ID + name)
# by "|".  Unstratified rows lack "|" entirely.
#
# Strategy:
#   1. Keep only rows that contain "|" (stratified rows).
#   2. Parse pathway_id as everything before the first ": " in the part
#      before "|" — equivalently, sub(":.*", "", sub("\\|.*", "", Feature_Raw)).
#   3. Filter to rows whose parsed pathway_id is in target_ids.

target_ids <- anchor_pathways$pathway_id

# Step 1 + 2: filter stratified rows and parse IDs in one pass using data.table
sample_cols <- setdiff(names(raw), "Feature_Raw")

raw_strat <- raw[grepl("|", Feature_Raw, fixed = TRUE)]
raw_strat[, pathway_id_parsed := sub(":.*", "", sub("\\|.*", "", Feature_Raw))]

# Step 3: keep only target pathways
raw_strat <- raw_strat[pathway_id_parsed %in% target_ids]
cat(sprintf("  After filtering to anchor pathways + stratified rows: %d rows.\n", nrow(raw_strat)))

# Melt to long format
long_dt <- melt(
  raw_strat,
  id.vars      = c("Feature_Raw", "pathway_id_parsed"),
  measure.vars = sample_cols,
  variable.name = "SampleID_raw",
  value.name    = "RPK"
)[RPK > 0]   # drop zero-abundance entries to save memory

long_dt[, `:=`(
  pathway_id  = pathway_id_parsed,
  species_raw = sub("^[^|]+\\|", "", Feature_Raw),
  SampleID    = gsub("-", ".", as.character(SampleID_raw))
)]

# Drop "unclassified" rows (not a real species)
long_dt <- long_dt[species_raw != "unclassified"]

# Remove UNMAPPED / UNINTEGRATED
long_dt <- long_dt[!grepl("UNGROUPED|UNINTEGRATED|UNMAPPED", species_raw)]

# Clean species name: remove g__ prefix, s__ prefix, underscores → spaces
long_dt[, species := {
  s <- species_raw
  s <- sub("^g__[^.]+\\.", "", s)    # remove g__Genus. prefix
  s <- sub("^s__", "", s)            # remove s__ prefix
  gsub("_", " ", s)                   # underscores → spaces
}]

cat(sprintf("  %d non-zero (pathway, species, sample) triplets.\n", nrow(long_dt)))


# ======================================================================
# 4. Calculate relative contribution (%) per (SampleID, pathway)
# ======================================================================

cat("\n=== Step 4: Computing relative contributions ===\n")

contrib_dt <- long_dt[, .(
  species    = species,
  pathway_id = pathway_id,
  RPK        = RPK
), by = SampleID]

# For each (SampleID, pathway): contribution_ij = RPK_ij / sum_j(RPK_ij) × 100
contrib_dt[, total_RPK := sum(RPK), by = .(SampleID, pathway_id)]
contrib_dt[, contribution_pct := RPK / total_RPK * 100]

# Sanity check: contributions should sum to 100 per (SampleID, pathway)
check <- contrib_dt[, .(sum_pct = sum(contribution_pct)), by = .(SampleID, pathway_id)]
stopifnot("Contribution percentages do not sum to 100!" =
            all(abs(check$sum_pct - 100) < 1e-6))
cat("  Sanity check passed: all (SampleID, pathway) contributions sum to 100%.\n")

# Join metadata
contrib_df <- as.data.frame(contrib_dt) %>%
  left_join(meta, by = "SampleID") %>%
  filter(!is.na(Enterotype)) %>%
  left_join(select(anchor_pathways, pathway_id, label_short, functional_theme, enterotype),
            by = "pathway_id")

cat(sprintf("  %d (sample × pathway × species) triplets after metadata join.\n", nrow(contrib_df)))


# ======================================================================
# 5. Annotate species with Section 3.3 core status
# ======================================================================

cat("\n=== Step 5: Annotating core species ===\n")

e1_sp <- prep$taxa_E1_species      # 73 Type 1 core species
e4_sp <- prep$taxa_E4_species      # 3  Type 4 core species
shared_sp <- prep$taxa_shared_species  # 1 (Veillonella dispar)

# Partial match: species column may be "Prevotella copri" while list has "Prevotella copri"
classify_species <- function(species_vec) {
  classify_one <- function(sp) {
    if (any(sapply(shared_sp, function(s) grepl(s, sp, fixed = TRUE)))) return("Shared_core")
    if (any(sapply(e1_sp,     function(s) grepl(s, sp, fixed = TRUE)))) return("E1_core")
    if (any(sapply(e4_sp,     function(s) grepl(s, sp, fixed = TRUE)))) return("E4_core")
    return("Non_core")
  }
  vapply(species_vec, classify_one, character(1), USE.NAMES = FALSE)
}

contrib_df <- contrib_df %>%
  mutate(species_core_status = classify_species(species))

cat("  Core species annotation distribution:\n")
print(table(contrib_df$species_core_status))


# ======================================================================
# 6. Build summary table: mean contribution per (pathway × species × enterotype)
# ======================================================================

cat("\n=== Step 6: Building contribution summary ===\n")

# For each pathway, identify top N species by overall mean contribution
top_n_per_pathway <- 10

# Mean contribution by (pathway, species, enterotype)
summary_df <- contrib_df %>%
  group_by(pathway_id, label_short, functional_theme, enterotype_anchor = enterotype,
           species, species_core_status, Enterotype) %>%
  summarise(
    n_samples          = n(),
    mean_contribution  = mean(contribution_pct),
    median_contribution = median(contribution_pct),
    sd_contribution    = sd(contribution_pct),
    .groups = "drop"
  )

# Determine top-N contributors per pathway (by total mean contribution across all enterotypes)
top_species_per_pathway <- summary_df %>%
  group_by(pathway_id, species) %>%
  summarise(total_mean = sum(mean_contribution), .groups = "drop") %>%
  group_by(pathway_id) %>%
  slice_max(order_by = total_mean, n = top_n_per_pathway) %>%
  ungroup() %>%
  select(pathway_id, species)

summary_df <- summary_df %>%
  mutate(species_plot = if_else(
    paste(pathway_id, species) %in% paste(top_species_per_pathway$pathway_id,
                                           top_species_per_pathway$species),
    species, "Others"
  ))

write.csv(summary_df, file.path(out_dir, "pathway_contribution_summary.csv"), row.names = FALSE)
cat(sprintf("  Summary table: %d rows.\n", nrow(summary_df)))


# ======================================================================
# 7. Statistical testing: Kruskal-Wallis + post-hoc Dunn (E1 vs E4)
# ======================================================================

cat("\n=== Step 7: Statistical testing ===\n")

# For each (pathway, species) pair: test whether contribution % differs across enterotypes
# Focus on species that appear in ≥ 2 enterotype groups with n ≥ 3 per group

kw_results <- contrib_df %>%
  group_by(pathway_id, label_short, species, species_core_status) %>%
  filter(n_distinct(Enterotype) >= 2) %>%
  do({
    df_sub <- .
    # Kruskal-Wallis across all 4 enterotypes
    kw <- tryCatch(
      kruskal_test(df_sub, contribution_pct ~ Enterotype),
      error = function(e) tibble(statistic = NA_real_, p = NA_real_, df = NA_real_, method = "failed")
    )
    kw
  }) %>%
  ungroup() %>%
  mutate(p_adj_kw = p.adjust(p, method = "fdr")) %>%
  arrange(p_adj_kw)

# Post-hoc Dunn test restricted to top contributors (overall mean > 5%) to limit burden
top_species_list <- top_species_per_pathway %>%
  pull(species) %>%
  unique()

dunn_results <- contrib_df %>%
  filter(species %in% top_species_list) %>%
  group_by(pathway_id, label_short, species, species_core_status) %>%
  filter(n_distinct(Enterotype) >= 2) %>%
  do({
    df_sub <- .
    dunn <- tryCatch(
      dunn_test(df_sub, contribution_pct ~ Enterotype, p.adjust.method = "BH"),
      error = function(e) tibble()
    )
    dunn
  }) %>%
  ungroup()

# Highlight E1 vs E4 Dunn comparisons specifically
dunn_E1_vs_E4 <- dunn_results %>%
  filter(
    (group1 == "1" & group2 == "4") | (group1 == "4" & group2 == "1")
  ) %>%
  arrange(p.adj)

# Combine into a single stats table
stats_out <- kw_results %>%
  select(pathway_id, label_short, species, species_core_status,
         kw_statistic = statistic, kw_p = p, kw_p_adj = p_adj_kw) %>%
  left_join(
    dunn_E1_vs_E4 %>%
      select(pathway_id = pathway_id, species,
             dunn_E1_E4_p = p, dunn_E1_E4_p_adj = p.adj, dunn_E1_E4_sig = p.adj.signif),
    by = c("pathway_id", "species")
  )

write.csv(stats_out,    file.path(out_dir, "pathway_contribution_stats.csv"),   row.names = FALSE)
cat(sprintf("  KW results: %d rows; %d significant (FDR < 0.05)\n",
            nrow(kw_results), sum(kw_results$p_adj_kw < 0.05, na.rm = TRUE)))
cat(sprintf("  Dunn E1-vs-E4 results: %d comparisons; %d significant (BH < 0.05)\n",
            nrow(dunn_E1_vs_E4), sum(dunn_E1_vs_E4$p.adj < 0.05, na.rm = TRUE)))


# ======================================================================
# 8. Save long-format data for visualization
# ======================================================================

cat("\n=== Step 8: Saving contribution data bundle ===\n")

saveRDS(
  list(
    contrib_df        = contrib_df,
    summary_df        = summary_df,
    kw_results        = kw_results,
    dunn_results      = dunn_results,
    dunn_E1_vs_E4     = dunn_E1_vs_E4,
    top_species       = top_species_per_pathway,
    anchor_pathways   = anchor_pathways
  ),
  file.path(out_dir, "pathway_contribution_data.rds")
)

cat("\n============================\n")
cat(" Phase 2 Complete\n")
cat("============================\n")
cat("Outputs written to:", out_dir, "\n")
cat("  anchor_pathway_metadata.csv\n")
cat("  pathway_contribution_summary.csv\n")
cat("  pathway_contribution_stats.csv\n")
cat("  pathway_contribution_data.rds\n")

# Print a quick preview of top contributor per anchor pathway
cat("\n--- Top contributor per anchor pathway (across all enterotypes) ---\n")
summary_df %>%
  group_by(pathway_id, label_short, species) %>%
  summarise(overall_mean = mean(mean_contribution), .groups = "drop") %>%
  group_by(pathway_id, label_short) %>%
  slice_max(overall_mean, n = 1) %>%
  arrange(label_short) %>%
  print(n = Inf)
