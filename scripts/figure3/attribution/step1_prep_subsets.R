# ======================================================================
#
#  Section 3.7 — Step 1: Data Preparation
#
#  Purpose:
#    1. Load Section 3.3 core feature tables (pathway / EC / taxa).
#    2. Extract matched rows from raw HUMAnN3 stratified TSV files
#       (all 95 core pathways + 36 shared ECs).
#    3. Load phyloseq metadata (SampleID, Enterotype).
#    4. Save a prep_data.rds bundle for downstream steps.
#
#  Inputs  (read-only):
#    results/tables/section_3_3/pathway_E1_core.csv        (43 pathways)
#    results/tables/section_3_3/pathway_E4_core.csv        (52 pathways)
#    results/tables/section_3_3/ec_Shared.csv              (36 ECs)
#    results/tables/section_3_3/ec_shared_direction.csv
#    results/tables/section_3_3/taxa_E1_core.csv
#    results/tables/section_3_3/taxa_E4_core.csv
#    results/tables/section_3_3/taxa_Shared.csv
#    data/clean/pathabundance_stratified.tsv
#    data/clean/ec_stratified.tsv
#    data/clean/phyloseq_EC_with_enterotype_20250210.rds
#
#  Outputs:
#    results/tables/section_3_7/pathway_stratified_subset.rds
#    results/tables/section_3_7/ec_stratified_subset.rds
#    results/tables/section_3_7/prep_data.rds
#
# ======================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(phyloseq)
  library(here)
})

out_dir <- here("results", "tables", "section_3_3", "attribution")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ======================================================================
# 1. Load Section 3.3 core feature tables
# ======================================================================

cat("=== Step 1: Loading Section 3.3 core feature tables ===\n")

pathway_E1_tbl   <- read.csv(here("results", "tables", "section_3_3", "pathway_E1_core.csv"))
pathway_E4_tbl   <- read.csv(here("results", "tables", "section_3_3", "pathway_E4_core.csv"))
ec_shared_tbl    <- read.csv(here("results", "tables", "section_3_3", "ec_Shared.csv"))
ec_dir_tbl       <- read.csv(here("results", "tables", "section_3_3", "ec_shared_direction.csv"))
taxa_E1_tbl      <- read.csv(here("results", "tables", "section_3_3", "taxa_E1_core.csv"))
taxa_E4_tbl      <- read.csv(here("results", "tables", "section_3_3", "taxa_E4_core.csv"))
taxa_shared_tbl  <- read.csv(here("results", "tables", "section_3_3", "taxa_Shared.csv"))

# Extract pathway IDs: text before the first ": "
extract_pathway_id <- function(x) sub(":.*", "", trimws(x))

pathway_E1_ids  <- extract_pathway_id(pathway_E1_tbl$feature)
pathway_E4_ids  <- extract_pathway_id(pathway_E4_tbl$feature)
all_pathway_ids <- unique(c(pathway_E1_ids, pathway_E4_ids))

# Confirm zero overlap (expected from Section 3.3)
overlap <- intersect(pathway_E1_ids, pathway_E4_ids)
stopifnot("Unexpected pathway overlap between E1 and E4 core!" = length(overlap) == 0)

# EC IDs are already bare numeric strings (e.g., "1.4.7.1")
ec_shared_ids <- unique(ec_shared_tbl$feature)

# Annotate shared ECs with direction (Convergent / Divergent)
ec_shared_annotated <- merge(ec_shared_tbl, ec_dir_tbl, by = "feature", all.x = TRUE)

# Extract species display names from full taxonomy paths (…;s__Genus_species → "Genus species")
extract_species_name <- function(full_path) {
  sub(".*s__", "", full_path) |> gsub("_", " ", x = _)
}

taxa_E1_species   <- extract_species_name(taxa_E1_tbl$feature)
taxa_E4_species   <- extract_species_name(taxa_E4_tbl$feature)
taxa_shared_sp    <- extract_species_name(taxa_shared_tbl$feature)

cat(sprintf(
  "  Pathway IDs: %d E1 + %d E4 = %d total, %d overlapping (expected 0)\n",
  length(pathway_E1_ids), length(pathway_E4_ids), length(all_pathway_ids), length(overlap)
))
cat(sprintf(
  "  Shared ECs : %d total, %d divergent, %d convergent\n",
  nrow(ec_shared_tbl),
  sum(ec_dir_tbl$direction_type == "Divergent", na.rm = TRUE),
  sum(ec_dir_tbl$direction_type == "Convergent", na.rm = TRUE)
))
cat(sprintf(
  "  Core species: %d E1 | %d E4 | %d shared\n",
  length(taxa_E1_species), length(taxa_E4_species), length(taxa_shared_sp)
))


# ======================================================================
# 2. Extract stratified subsets from raw TSV files
# ======================================================================

#' Read the full stratified TSV and return rows matching target feature IDs.
#'
#' @param feature_ids  Character vector of feature IDs to match.
#' @param raw_path     Path to the raw stratified TSV file.
#' @param output_rds   Destination RDS path.
#' @param is_ec        If TRUE, use boundary-anchored pattern to avoid partial
#'                     matches (e.g., "1.1.1.1" vs "1.1.1.10").
extract_stratified_subset <- function(feature_ids, raw_path, output_rds, is_ec = FALSE) {
  cat(sprintf("\n=== Extracting from %s ===\n", basename(raw_path)))
  cat(sprintf("  Reading file... "))
  dt <- fread(raw_path, sep = "\t", showProgress = FALSE, verbose = FALSE)
  setnames(dt, 1, "Feature_Raw")
  cat(sprintf("%d rows read.\n", nrow(dt)))

  if (is_ec) {
    # Escape regex metacharacters in EC numbers (dots → literal dots),
    # then anchor with [:|] or $ to prevent "1.1.1.1" matching "1.1.1.10"
    esc_ids <- gsub("\\.", "\\\\.", feature_ids)
    pattern <- paste0("^(", paste(esc_ids, collapse = "|"), ")([:|]|$)")
  } else {
    # Pathway IDs are uniquely prefixed alphanumeric+hyphen strings;
    # a simple start-anchor is sufficient
    pattern <- paste0("^(", paste(feature_ids, collapse = "|"), ")")
  }

  subset_dt <- dt[grepl(pattern, Feature_Raw, perl = TRUE)]
  cat(sprintf("  Matched  : %d rows (%.1f%%)\n",
              nrow(subset_dt), 100 * nrow(subset_dt) / nrow(dt)))

  saveRDS(subset_dt, output_rds)
  cat(sprintf("  Saved to : %s\n", basename(output_rds)))

  rm(dt)
  gc(verbose = FALSE)
  return(output_rds)
}

pathway_rds <- extract_stratified_subset(
  feature_ids = all_pathway_ids,
  raw_path    = here("data", "clean", "pathabundance_stratified.tsv"),
  output_rds  = file.path(out_dir, "pathway_stratified_subset.rds"),
  is_ec       = FALSE
)

ec_rds <- extract_stratified_subset(
  feature_ids = ec_shared_ids,
  raw_path    = here("data", "clean", "ec_stratified.tsv"),
  output_rds  = file.path(out_dir, "ec_stratified_subset.rds"),
  is_ec       = TRUE
)


# ======================================================================
# 3. Load metadata from phyloseq object
# ======================================================================

cat("\n=== Step 3: Loading metadata ===\n")

phy <- readRDS(here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"))

metadata_df <- as(sample_data(phy), "data.frame") |>
  rownames_to_column("SampleID") |>
  select(SampleID, Enterotype = Dirichlettype_Genus) |>
  mutate(
    SampleID   = gsub("-", ".", SampleID),
    Enterotype = factor(Enterotype, levels = c("1", "2", "3", "4"))
  )

cat(sprintf("  %d samples. Enterotype distribution:\n", nrow(metadata_df)))
print(table(metadata_df$Enterotype))


# ======================================================================
# 4. Verification — check pathway ID coverage in extracted subset
# ======================================================================

cat("\n=== Step 4: Verifying pathway coverage ===\n")

pathway_subset  <- readRDS(pathway_rds)
found_ids       <- unique(extract_pathway_id(
  pathway_subset$Feature_Raw[!grepl("UNMAPPED|UNINTEGRATED", pathway_subset$Feature_Raw)]
))
missing_ids     <- setdiff(all_pathway_ids, found_ids)

if (length(missing_ids) == 0) {
  cat(sprintf("  OK: All %d pathway IDs found in stratified subset.\n", length(all_pathway_ids)))
} else {
  cat(sprintf("  WARNING: %d pathway IDs NOT found:\n", length(missing_ids)))
  print(missing_ids)
}

# Check EC coverage
ec_subset  <- readRDS(ec_rds)
found_ecs  <- unique(sub("([:|]|$).*", "", ec_subset$Feature_Raw))
missing_ecs <- setdiff(ec_shared_ids, found_ecs)

if (length(missing_ecs) == 0) {
  cat(sprintf("  OK: All %d shared EC IDs found in stratified subset.\n", length(ec_shared_ids)))
} else {
  cat(sprintf("  WARNING: %d EC IDs NOT found:\n", length(missing_ecs)))
  print(missing_ecs)
}

rm(pathway_subset, ec_subset)
gc(verbose = FALSE)


# ======================================================================
# 5. Save prep_data bundle for downstream steps
# ======================================================================

cat("\n=== Step 5: Saving prep_data.rds ===\n")

saveRDS(
  list(
    metadata             = metadata_df,
    pathway_E1_ids       = pathway_E1_ids,
    pathway_E4_ids       = pathway_E4_ids,
    all_pathway_ids      = all_pathway_ids,
    pathway_E1_core      = pathway_E1_tbl,
    pathway_E4_core      = pathway_E4_tbl,
    ec_shared_annotated  = ec_shared_annotated,
    ec_shared_ids        = ec_shared_ids,
    taxa_E1_species      = taxa_E1_species,
    taxa_E4_species      = taxa_E4_species,
    taxa_shared_species  = taxa_shared_sp,
    pathway_rds_path     = pathway_rds,
    ec_rds_path          = ec_rds
  ),
  file.path(out_dir, "prep_data.rds")
)

cat("\n============================\n")
cat(" Phase 1 Complete\n")
cat("============================\n")
cat("Outputs written to:", out_dir, "\n")
cat("  pathway_stratified_subset.rds\n")
cat("  ec_stratified_subset.rds\n")
cat("  prep_data.rds\n")
