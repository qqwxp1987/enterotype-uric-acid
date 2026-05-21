# ======================================================================
#
#  Section 3.7 — Step 3: Divergent EC Species Contribution Analysis
#
#  Purpose:
#    For 8 highly-divergent shared ECs (E1/E4 opposite log2FC direction),
#    decompose per-sample species contributions and identify whether
#    different driver species account for the opposed regulation.
#
#    This directly tests the "same enzyme, different driver species"
#    hypothesis: an EC may be enriched in E4 because Species A
#    carries it there, but depleted in E1 because the same EC is
#    instead performed by Species B (or missing altogether).
#
#  EC selection criteria:
#    1. `direction_type == "Divergent"` in ec_shared_direction.csv
#    2. |E1_mean_lfc| + |E4_mean_lfc| >= ~4.8 (top 8 by magnitude)
#       • Exception: 2.4.2.45 excluded (only 3 named-species rows – too sparse)
#       • 1.8.1.8 included as sole "reverse divergent" representative
#         (E1 enriched, E4 depleted)
#    3. >= 30 named-species stratified rows in ec_stratified_subset.rds
#
#  NOTE on RPK units (same logic as Step 2):
#    For relative contribution (%) the RPK units cancel exactly within
#    each (SampleID, EC) group → no normalisation needed.
#    Absolute abundance comparisons are covered by LinDA in Section 3.3.
#
#  Inputs:
#    results/tables/section_3_7/prep_data.rds
#    results/tables/section_3_7/ec_stratified_subset.rds
#    results/tables/section_3_3/ec_shared_direction.csv
#
#  Outputs:
#    results/tables/section_3_7/divergent_ec_metadata.csv
#    results/tables/section_3_7/ec_contribution_summary.csv
#    results/tables/section_3_7/ec_contribution_stats.csv
#    results/tables/section_3_7/ec_contribution_data.rds
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

out_dir <- here("results", "tables", "section_3_3", "attribution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 1. Load prep data and direction table
# ======================================================================

cat("=== Step 1: Loading input data ===\n")
prep     <- readRDS(file.path(out_dir, "prep_data.rds"))
meta     <- prep$metadata

dir_tbl  <- read.csv(here("results", "tables", "section_3_3", "ec_shared_direction.csv"),
                     stringsAsFactors = FALSE)

cat(sprintf("  Metadata: %d samples\n", nrow(meta)))
cat(sprintf("  EC direction table: %d rows (%d Divergent)\n",
            nrow(dir_tbl), sum(dir_tbl$direction_type == "Divergent")))


# ======================================================================
# 2. Define target divergent ECs
# ======================================================================

# Rank divergent ECs by |E1| + |E4| LFC magnitude; exclude 2.4.2.45 (sparse)
div_tbl <- dir_tbl %>%
  filter(direction_type == "Divergent") %>%
  mutate(abs_sum = abs(E1_mean_lfc) + abs(E4_mean_lfc)) %>%
  arrange(desc(abs_sum))

# Selected 8: top 7 by lfc-sum (excl. 2.4.2.45) + 1.8.1.8 (reverse divergent)
target_ec <- tribble(
  ~ec_id,      ~divergence_pattern,  ~abs_lfc_sum,  ~E1_mean_lfc, ~E4_mean_lfc,
  "4.1.1.41",  "E1↓  E4↑",           8.03,          -4.21,         3.82,
  "1.4.7.1",   "E1↓  E4↑",           6.93,          -2.45,         4.49,
  "3.4.11.5",  "E1↓  E4↑",           6.82,          -3.21,         3.61,
  "3.6.1.13",  "E1↓  E4↑",           6.66,          -2.29,         4.36,
  "6.3.1.5",   "E1↓  E4↑",           5.91,          -2.02,         3.89,
  "1.5.1.7",   "E1↓  E4↑",           5.75,          -2.19,         3.56,
  "2.1.1.148", "E1↓  E4↑",           4.79,          -1.94,         2.85,
  "1.8.1.8",   "E1↑  E4↓",           1.77,           0.95,        -0.82,  # reverse pattern
  "2.7.6.3",   "E1↑  E4↑",           2.38,           0.41,         1.97,  # Convergent
  "2.6.1.1",   "E1↑  E4↑",           1.95,           0.48,         1.47,  # Convergent
  "2.1.1.13",  "E1↑  E4↑",           1.68,           0.31,         1.37   # Convergent
)

write.csv(target_ec, file.path(out_dir, "divergent_ec_metadata.csv"), row.names = FALSE)
cat(sprintf("  %d divergent ECs selected.\n", nrow(target_ec)))
print(target_ec[, c("ec_id", "divergence_pattern", "abs_lfc_sum")])


# ======================================================================
# 3. Load EC stratified subset and reshape
# ======================================================================

cat("\n=== Step 3: Loading EC stratified subset ===\n")

raw <- readRDS(file.path(out_dir, "ec_stratified_subset.rds"))
setnames(raw, 1, "Feature_Raw")
cat(sprintf("  %d rows loaded.\n", nrow(raw)))

# HUMAnN3 EC stratified format: "EC_ID|g__Genus.s__species"
# (no ": description" part — unlike pathway rows)
# Unstratified rows: bare "EC_ID"

sample_cols <- setdiff(names(raw), "Feature_Raw")

raw_strat <- raw[grepl("|", Feature_Raw, fixed = TRUE)]
raw_strat[, ec_id_parsed := sub("[|].*", "", Feature_Raw)]

# Filter to target ECs
raw_strat <- raw_strat[ec_id_parsed %in% target_ec$ec_id]
cat(sprintf("  After filtering to target divergent ECs + stratified rows: %d rows.\n", nrow(raw_strat)))

# Melt to long format
long_dt <- melt(
  raw_strat,
  id.vars       = c("Feature_Raw", "ec_id_parsed"),
  measure.vars  = sample_cols,
  variable.name = "SampleID_raw",
  value.name    = "RPK"
)[RPK > 0]

long_dt[, `:=`(
  ec_id       = ec_id_parsed,
  species_raw = sub("^[^|]+\\|", "", Feature_Raw),
  SampleID    = gsub("-", ".", as.character(SampleID_raw))
)]

# Drop "unclassified" rows
long_dt <- long_dt[species_raw != "unclassified"]
# Drop UNMAPPED / UNINTEGRATED
long_dt <- long_dt[!grepl("UNGROUPED|UNINTEGRATED|UNMAPPED", species_raw)]

# Clean species name
long_dt[, species := {
  s <- species_raw
  s <- sub("^g__[^.]+\\.", "", s)
  s <- sub("^s__", "", s)
  gsub("_", " ", s)
}]

cat(sprintf("  %d non-zero (EC, species, sample) triplets.\n", nrow(long_dt)))


# ======================================================================
# 4. Calculate relative contribution (%) per (SampleID, EC)
# ======================================================================

cat("\n=== Step 4: Computing relative contributions ===\n")

long_dt[, total_RPK       := sum(RPK),           by = .(SampleID, ec_id)]
long_dt[, contribution_pct := RPK / total_RPK * 100]

# Sanity check
check <- long_dt[, .(sum_pct = sum(contribution_pct)), by = .(SampleID, ec_id)]
stopifnot("EC contribution percentages do not sum to 100!" =
            all(abs(check$sum_pct - 100) < 1e-6))
cat("  Sanity check passed.\n")

# Join with metadata
contrib_df <- as.data.frame(long_dt) %>%
  inner_join(meta, by = "SampleID")
cat(sprintf("  %d triplets after metadata join.\n", nrow(contrib_df)))


# ======================================================================
# 5. Annotate core species
# ======================================================================

cat("\n=== Step 5: Annotating core species ===\n")

# Field names in prep_data.rds (set by step1):
#   taxa_E1_species (73 species, "Genus species" format)
#   taxa_E4_species (3 species)
#   taxa_shared_species
e1_sp     <- prep$taxa_E1_species
e4_sp     <- prep$taxa_E4_species
shared_sp <- prep$taxa_shared_species

# Use partial matching (grepl) to handle sub-species names like
# "Prevotella copri CAG 164" matching reference "Prevotella copri"
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
  mutate(core_status = classify_species(species))

cat("  Core species annotation distribution:\n")
print(table(contrib_df$core_status))


# ======================================================================
# 6. Build contribution summary (mean % per enterotype × species × EC)
# ======================================================================

cat("\n=== Step 6: Building contribution summary ===\n")

ec_meta_df <- target_ec  # small reference table

summary_df <- contrib_df %>%
  group_by(ec_id, Enterotype, species, core_status) %>%
  summarise(
    mean_contribution_pct = mean(contribution_pct),
    sd_contribution_pct   = sd(contribution_pct),
    n_samples             = n_distinct(SampleID),
    .groups = "drop"
  ) %>%
  left_join(ec_meta_df[, c("ec_id", "divergence_pattern",
                            "E1_mean_lfc", "E4_mean_lfc", "abs_lfc_sum")],
            by = "ec_id")

write.csv(summary_df, file.path(out_dir, "ec_contribution_summary.csv"),
          row.names = FALSE)
cat(sprintf("  Summary table: %d rows saved.\n", nrow(summary_df)))


# ======================================================================
# 7. Statistical testing: KW across enterotypes + Dunn E1 vs E4
# ======================================================================

cat("\n=== Step 7: Statistical testing ===\n")

# --- Kruskal-Wallis across 4 enterotypes ---
kw_results <- contrib_df %>%
  group_by(ec_id, species) %>%
  filter(n_distinct(Enterotype) >= 2) %>%    # need at least 2 enterotypes
  kruskal_test(contribution_pct ~ Enterotype) %>%
  adjust_pvalue(method = "fdr") %>%
  add_significance()

cat(sprintf("  KW results: %d rows; %d significant (FDR < 0.05)\n",
            nrow(kw_results),
            sum(kw_results$p.adj < 0.05, na.rm = TRUE)))

# --- Dunn test: E1 vs E4 pairwise ---
dunn_e1e4 <- contrib_df %>%
  filter(Enterotype %in% c(1, 4)) %>%
  mutate(Enterotype = droplevels(as.factor(Enterotype))) %>%
  group_by(ec_id, species) %>%
  filter(n_distinct(Enterotype) == 2) %>%
  dunn_test(contribution_pct ~ Enterotype, p.adjust.method = "BH") %>%
  filter(group1 == "1" & group2 == "4")

cat(sprintf("  Dunn E1-vs-E4 results: %d comparisons; %d significant (BH < 0.05)\n",
            nrow(dunn_e1e4),
            sum(dunn_e1e4$p.adj < 0.05, na.rm = TRUE)))

# Combine
stats_df <- bind_rows(
  kw_results %>% mutate(test_type = "KW_4group"),
  dunn_e1e4  %>% mutate(test_type = "Dunn_E1vsE4")
)

write.csv(stats_df, file.path(out_dir, "ec_contribution_stats.csv"),
          row.names = FALSE)


# ======================================================================
# 8. Key analysis: top driver species per (EC, enterotype)
# ======================================================================

cat("\n=== Step 8: Identifying top driver species ===\n")

top_drivers <- contrib_df %>%
  group_by(ec_id, Enterotype, species) %>%
  summarise(mean_pct = mean(contribution_pct), .groups = "drop") %>%
  group_by(ec_id, Enterotype) %>%
  slice_max(mean_pct, n = 3) %>%
  arrange(ec_id, Enterotype, desc(mean_pct))

cat("\nTop driver species per (EC, enterotype):\n")
print(as.data.frame(top_drivers), row.names = FALSE)

# E1 vs E4 comparison: who are the top different drivers?
e1_top <- top_drivers %>% filter(Enterotype == 1) %>%
  select(ec_id, E1_top_species = species, E1_mean_pct = mean_pct)
e4_top <- top_drivers %>% filter(Enterotype == 4) %>%
  slice_max(mean_pct, n = 1) %>%
  select(ec_id, E4_top_species = species, E4_mean_pct = mean_pct)

driver_comparison <- left_join(
  top_drivers %>% filter(Enterotype == 1) %>% slice_max(mean_pct, n = 1) %>%
    select(ec_id, E1_top = species, E1_pct = mean_pct),
  top_drivers %>% filter(Enterotype == 4) %>% slice_max(mean_pct, n = 1) %>%
    select(ec_id, E4_top = species, E4_pct = mean_pct),
  by = "ec_id"
) %>%
  left_join(ec_meta_df[, c("ec_id", "divergence_pattern", "E1_mean_lfc", "E4_mean_lfc")],
            by = "ec_id") %>%
  mutate(same_top_species = E1_top == E4_top)

cat("\n--- E1 vs E4 top driver comparison ---\n")
print(as.data.frame(driver_comparison), row.names = FALSE)


# ======================================================================
# 9. Save data bundle
# ======================================================================

cat("\n=== Step 9: Saving EC contribution data bundle ===\n")

saveRDS(
  list(
    contrib_long     = contrib_df,
    summary          = summary_df,
    stats_kw         = kw_results,
    stats_dunn       = dunn_e1e4,
    driver_comparison = driver_comparison,
    target_ec        = target_ec
  ),
  file.path(out_dir, "ec_contribution_data.rds")
)

cat("\n============================\n")
cat(" Phase 3 Complete\n")
cat("============================\n")
cat("Outputs written to:", out_dir, "\n")
cat("  divergent_ec_metadata.csv\n")
cat("  ec_contribution_summary.csv\n")
cat("  ec_contribution_stats.csv\n")
cat("  ec_contribution_data.rds\n\n")
