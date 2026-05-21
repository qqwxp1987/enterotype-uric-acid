# ==============================================================================
# Step 0: Verify Enterotype Label Consistency Across Data Sources
# The seed10 taxonomy phyloseq is the AUTHORITATIVE source.
# We check if Pathway/EC/Metabolome data have matching samples & labels.
# ==============================================================================

library(phyloseq)
library(here)

# --- 1. Load all data sources ---
ps_taxa <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
ps_path <- readRDS(here("data", "clean", "phyloseq_Pathway_with_enterotype_20250210.rds"))
ps_ec   <- readRDS(here("data", "clean", "phyloseq_EC_with_enterotype_20250210.rds"))
metab   <- readRDS(here("data", "clean", "metabolites_normalized_multiIS_and_metadata_20250904.rds"))

# --- 2. Extract sample metadata ---
meta_taxa <- as(sample_data(ps_taxa), "data.frame")
meta_path <- as(sample_data(ps_path), "data.frame")
meta_ec   <- as(sample_data(ps_ec), "data.frame")

cat("=== Sample Counts ===\n")
cat(sprintf("Taxa (seed10):  %d samples\n", nrow(meta_taxa)))
cat(sprintf("Pathway:        %d samples\n", nrow(meta_path)))
cat(sprintf("EC:             %d samples\n", nrow(meta_ec)))

# --- 3. Check sample ID overlap ---
taxa_ids <- rownames(meta_taxa)
path_ids <- rownames(meta_path)
ec_ids   <- rownames(meta_ec)

cat("\n=== Sample ID Overlap ===\n")
cat(sprintf("Taxa ∩ Pathway: %d / %d taxa samples\n",
            length(intersect(taxa_ids, path_ids)), length(taxa_ids)))
cat(sprintf("Taxa ∩ EC:      %d / %d taxa samples\n",
            length(intersect(taxa_ids, ec_ids)), length(taxa_ids)))

# Pathway/EC samples not in taxa
path_only <- setdiff(path_ids, taxa_ids)
ec_only   <- setdiff(ec_ids, taxa_ids)
if (length(path_only) > 0) cat(sprintf("Pathway-only samples: %s\n", paste(path_only, collapse = ", ")))
if (length(ec_only) > 0)   cat(sprintf("EC-only samples: %s\n", paste(ec_only, collapse = ", ")))

taxa_not_path <- setdiff(taxa_ids, path_ids)
taxa_not_ec   <- setdiff(taxa_ids, ec_ids)
if (length(taxa_not_path) > 0) cat(sprintf("Taxa samples missing from Pathway: %s\n", paste(taxa_not_path, collapse = ", ")))
if (length(taxa_not_ec) > 0)   cat(sprintf("Taxa samples missing from EC: %s\n", paste(taxa_not_ec, collapse = ", ")))

# --- 4. Check enterotype column names ---
cat("\n=== Enterotype Column Check ===\n")
et_col_taxa <- grep("enterotype|dirichlet|type", names(meta_taxa), ignore.case = TRUE, value = TRUE)
et_col_path <- grep("enterotype|dirichlet|type", names(meta_path), ignore.case = TRUE, value = TRUE)
et_col_ec   <- grep("enterotype|dirichlet|type", names(meta_ec),   ignore.case = TRUE, value = TRUE)
cat(sprintf("Taxa columns:    %s\n", paste(et_col_taxa, collapse = ", ")))
cat(sprintf("Pathway columns: %s\n", paste(et_col_path, collapse = ", ")))
cat(sprintf("EC columns:      %s\n", paste(et_col_ec, collapse = ", ")))

# --- 5. Compare enterotype labels on shared samples ---
# Use Dirichlettype_Genus as the enterotype column (based on figure1/figure2 scripts)
shared_path <- intersect(taxa_ids, path_ids)
shared_ec   <- intersect(taxa_ids, ec_ids)

if ("Dirichlettype_Genus" %in% names(meta_taxa)) {
  et_col <- "Dirichlettype_Genus"
} else {
  # Try to find the right column
  et_col <- et_col_taxa[1]
}
cat(sprintf("\nUsing enterotype column: '%s'\n", et_col))

if (length(shared_path) > 0 && et_col %in% names(meta_path)) {
  match_path <- sum(meta_taxa[shared_path, et_col] == meta_path[shared_path, et_col])
  cat(sprintf("\nPathway enterotype match: %d / %d (%.1f%%)\n",
              match_path, length(shared_path), 100 * match_path / length(shared_path)))

  # Show mismatches
  mismatches <- shared_path[meta_taxa[shared_path, et_col] != meta_path[shared_path, et_col]]
  if (length(mismatches) > 0 && length(mismatches) <= 20) {
    cat("Mismatched samples:\n")
    for (s in mismatches) {
      cat(sprintf("  %s: taxa=%s, pathway=%s\n", s,
                  meta_taxa[s, et_col], meta_path[s, et_col]))
    }
  } else if (length(mismatches) > 20) {
    cat(sprintf("  %d mismatches (too many to list)\n", length(mismatches)))
  }
} else {
  cat("\nCannot compare Pathway enterotypes (column missing or no shared samples)\n")
}

if (length(shared_ec) > 0 && et_col %in% names(meta_ec)) {
  match_ec <- sum(meta_taxa[shared_ec, et_col] == meta_ec[shared_ec, et_col])
  cat(sprintf("EC enterotype match: %d / %d (%.1f%%)\n",
              match_ec, length(shared_ec), 100 * match_ec / length(shared_ec)))
}

# --- 6. Check metabolome data structure ---
cat("\n=== Metabolome Data ===\n")
if (is.list(metab) && !is.data.frame(metab)) {
  cat(sprintf("Structure: list with elements: %s\n", paste(names(metab), collapse = ", ")))
  for (nm in names(metab)) {
    obj <- metab[[nm]]
    if (is.data.frame(obj) || is.matrix(obj)) {
      cat(sprintf("  $%s: %d x %d\n", nm, nrow(obj), ncol(obj)))
    } else {
      cat(sprintf("  $%s: class=%s\n", nm, class(obj)[1]))
    }
  }

  # Check if metadata has sample IDs matching taxa
  if ("metadata" %in% names(metab)) {
    metab_meta <- metab$metadata
    metab_ids <- rownames(metab_meta)
    if (is.null(metab_ids)) metab_ids <- metab_meta$SampleID
    cat(sprintf("\nMetabolome sample count: %d\n", length(metab_ids)))
    cat(sprintf("Metabolome ∩ Taxa: %d\n", length(intersect(metab_ids, taxa_ids))))

    # Check enterotype column
    et_col_metab <- grep("enterotype|dirichlet|type", names(metab_meta), ignore.case = TRUE, value = TRUE)
    cat(sprintf("Metabolome enterotype columns: %s\n", paste(et_col_metab, collapse = ", ")))
  }
} else if (inherits(metab, "phyloseq")) {
  cat("Structure: phyloseq object\n")
  cat(sprintf("  Samples: %d, Features: %d\n", nsamples(metab), ntaxa(metab)))
} else {
  cat(sprintf("Structure: %s\n", class(metab)[1]))
}

# --- 7. Enterotype distribution in seed10 ---
cat("\n=== Enterotype Distribution (seed10, authoritative) ===\n")
print(table(meta_taxa[[et_col]]))

cat("\nStep 0 verification complete.\n")
