# Build Table S4: candidate-mediator selection traceability
#
# Reconstructs the layer-specific selection evidence for the 22 candidates used
# in the exploratory mediation screen.

set.seed(2024)

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(project_root, "results", "tables")
output_path <- file.path(
  table_dir, "final", "Table_S4_mediator_traceability.csv"
)
detailed_output_path <- file.path(
  table_dir, "supplementary", "mediator_traceability_table_detailed.csv"
)

read_source <- function(...) {
  path <- file.path(...)
  if (!file.exists(path)) stop("Missing required source: ", path)
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

limma_all <- read_source(
  table_dir, "section_3_4", "limma_metabolome_all_comparisons.csv"
)
path_e1 <- read_source(table_dir, "section_3_4", "pathway_E1_core.csv")
path_e4 <- read_source(table_dir, "section_3_4", "pathway_E4_core.csv")
ec_e1 <- read_source(table_dir, "section_3_4", "ec_E1_core.csv")
ec_e4 <- read_source(table_dir, "section_3_4", "ec_E4_core.csv")
envfit <- read_source(
  table_dir, "section_3_4", "envfit_curated_arrows.csv"
)
partial_cor <- read_source(
  table_dir, "section_3_4", "feature_UA_partial_cor_summary.csv"
)
partial_cor <- partial_cor[partial_cor$model == "Base", , drop = FALSE]

candidates <- data.frame(
  candidate_id = paste0("M", 1:22),
  mediator_name = c(
    "5-Methoxyindoleacetate", "p-Cresol", "p-Cresol sulfate",
    "Indoxyl sulfate", "Indole-3-propionic acid",
    "Phenylacetylglutamine", "Deoxycholic acid",
    "Trimethylamine N-oxide", "Adenosine",
    "PWY-5695: urate biosynthesis / IMP degradation",
    "PWY0-1296: purine ribonucleosides degradation",
    "PWY-5004: L-citrulline metabolism",
    "PWY-7219: adenosine ribonucleotides de novo biosynthesis",
    "PWY-6609: adenine and adenosine salvage III",
    "ARGORNPROST-PWY: Arg/Orn/Pro interconversion",
    "PWY-7234: IMP biosynthesis III",
    "EC 2.6.1.1: aspartate transaminase",
    "EC 2.3.1.129: UDP-GlcNAc O-acyltransferase",
    "EC 4.2.1.8", "EC 2.7.6.3: HPPK",
    "EC 3.1.3.15: histidinol-phosphatase",
    "EC 2.1.1.13: methionine synthase"
  ),
  category = c(rep("Metabolite", 9), rep("Pathway", 7), rep("EC", 6)),
  chemical_class = c(
    rep("Aromatic amino acid catabolite", 6),
    "Bile acid", "TMA pathway", "Purine nucleoside",
    "Purine metabolism", "Purine metabolism", "Amino acid metabolism",
    "Purine metabolism", "Purine metabolism", "Amino acid metabolism",
    "Purine metabolism", rep("Convergent shared EC", 6)
  ),
  source_feature = c(
    "NEG0144_5-Methoxyindoleacetate", "NEG0024_p-Cresol",
    "NEG0216_p-Cresol sulfate", "NEG0046_Indoxyl sulfate",
    "NEG0147_Indole-3-propionic acid",
    "NEG0145_Phenylacetylglutamine", "NEG0015_Deoxycholic acid",
    "POS0156_Trimethylamine N-oxide", "POS0002_Adenosine",
    "PWY-5695: urate biosynthesis/inosine 5-phosphate degradation",
    "PWY0-1296: purine ribonucleosides degradation",
    "PWY-5004: superpathway of L-citrulline metabolism",
    "PWY-7219: adenosine ribonucleotides de novo biosynthesis",
    "PWY-6609: adenine and adenosine salvage III",
    "ARGORNPROST-PWY: arginine, ornithine and proline interconversion",
    "PWY-7234: inosine-5-phosphate biosynthesis III",
    "2.6.1.1", "2.3.1.129", "4.2.1.8", "2.7.6.3", "3.1.3.15",
    "2.1.1.13"
  ),
  source_layer = c(
    rep("Type 1 core serum metabolite", 7),
    "Added host-facing thematic anchor for the TMA/betaine axis",
    "Added host-facing anchor for the upstream microbial purine theme",
    "Type 4 core pathway", "Type 1 core pathway", "Type 1 core pathway",
    "Type 4 core pathway", "Type 4 core pathway", "Type 1 core pathway",
    "Type 1 core pathway", rep("Shared convergent core EC", 6)
  ),
  selection_basis = c(
    rep(
      paste0(
        "Type 1 core serum metabolite selected as a biologically interpretable ",
        "AAA-related candidate; nominal post hoc envfit P<0.05"
      ),
      6
    ),
    paste0(
      "Bile-acid representative with higher envfit R2 than taurocholic acid ",
      "(0.098 vs 0.030); E1 core"
    ),
    paste0(
      "TMA-related host-facing anchor paired with betaine; higher envfit R2 ",
      "than betaine (0.055 vs 0.021); not selected on core status"
    ),
    paste0(
      "Purine host-facing anchor motivated by upstream microbial purine ",
      "results; slightly higher envfit R2 than adenine (0.028 vs 0.023)"
    ),
    "E4 core pathway selected to represent annotated urate/IMP metabolism",
    paste0(
      "E1 core pathway selected to represent purine degradation; also ",
      "UA-partial-correlation FDR<0.05"
    ),
    paste0(
      "E1 core pathway selected to represent amino-acid/nitrogen metabolism; ",
      "also UA-partial-correlation FDR<0.05"
    ),
    "E4 core pathway selected to represent de novo purine biosynthesis",
    "E4 core pathway selected to represent purine salvage",
    paste0(
      "E1 core pathway selected to represent amino-acid/nitrogen metabolism; ",
      "also UA-partial-correlation FDR<0.05"
    ),
    "E1 core pathway selected to represent IMP biosynthesis",
    "Criterion A: shared convergent core EC with UA partial-correlation FDR<0.05",
    "Criterion A: shared convergent core EC with UA partial-correlation FDR<0.05",
    "Criterion A: shared convergent core EC with UA partial-correlation FDR<0.05",
    paste0(
      "Criterion B: shared convergent core EC retained on functional annotation: ",
      "one of the three Section 3.3 convergent enzymes whose IUBMB/KEGG ",
      "annotation points to folate-cycle one-carbon units upstream of de novo ",
      "purine biosynthesis; UA correlation not significant (FDR=0.48)"
    ),
    "Criterion A: shared convergent core EC with UA partial-correlation FDR<0.05",
    paste0(
      "Criterion B: shared convergent core EC retained on functional annotation: ",
      "one of the three Section 3.3 convergent enzymes whose IUBMB/KEGG ",
      "annotation points to SAM-cycle recycling upstream of de novo purine ",
      "biosynthesis; UA correlation not significant (FDR=0.13)"
    )
  ),
  effect_size_role = c(
    rep("Differential log2FC considered descriptively; envfit R2 reported; no log2FC cutoff used", 6),
    rep("Envfit R2 used only for within-class representative selection", 3),
    rep("Core log2FC considered descriptively; no additional log2FC cutoff used", 7),
    rep("UA rho reported; selection based on FDR, not a rho cutoff", 3),
    "Core log2FC considered descriptively; functional nomination, not a UA-effect cutoff",
    "UA rho reported; selection based on FDR, not a rho cutoff",
    "Core log2FC considered descriptively; functional nomination, not a UA-effect cutoff"
  ),
  biological_rationale = c(
    rep(
      paste0(
        "Selected AAA-related member of the Type 1 core serum pattern; the ",
        "focused metabolite set was non-exhaustive"
      ),
      6
    ),
    paste0(
      "Secondary bile acid representing the reciprocal bile-acid component ",
      "of the Type 1 serum signature"
    ),
    paste0(
      "Circulating microbiota-dependent methylamine anchor; not treated as an ",
      "established UA-associated feature"
    ),
    paste0(
      "Circulating purine-nucleoside anchor motivated by the upstream microbial ",
      "purine theme; not treated as an established UA-associated feature"
    ),
    "MetaCyc-annotated urate biosynthesis / IMP-degradation pathway",
    "MetaCyc-annotated purine ribonucleoside degradation pathway",
    "MetaCyc-annotated L-citrulline metabolism pathway",
    "MetaCyc-annotated adenosine-ribonucleotide biosynthesis pathway",
    "MetaCyc-annotated adenine/adenosine salvage pathway",
    "MetaCyc-annotated Arg/Orn/Pro interconversion pathway",
    "MetaCyc-annotated IMP biosynthesis pathway",
    "Pleiotropic aminotransferase; retained on convergent-core plus UA-correlation evidence",
    "Lipid/amino-sugar enzyme; retained on convergent-core plus UA-correlation evidence",
    "Pleiotropic hydratase/dehydratase; retained on convergent-core plus UA-correlation evidence",
    paste0(
      "HPPK participates in bacterial folate biosynthesis; the connection to ",
      "purine biology is indirect and does not establish purine flux"
    ),
    "Histidine-metabolism enzyme; retained on convergent-core plus UA-correlation evidence",
    paste0(
      "Methionine synthase is a folate-dependent one-carbon enzyme; the ",
      "connection to purine biology is indirect and does not establish purine flux"
    )
  ),
  literature_support = c(
    rep("Current-cohort core-metabolome and envfit evidence; no external citation used for selection", 7),
    paste0(
      "Tang et al. N Engl J Med. 2013;368:1575-1584. ",
      "doi:10.1056/NEJMoa1109400 (supports microbiota-dependent TMAO biology, not UA specificity)"
    ),
    paste0(
      "KEGG Purine metabolism map00230: https://www.kegg.jp/entry/map00230 ",
      "(supports purine-class assignment, not adenosine-specific mediation)"
    ),
    rep("MetaCyc pathway annotation plus current-cohort core-pathway evidence", 7),
    rep("Current-cohort convergent-core and UA-correlation evidence", 3),
    paste0(
      "IUBMB EC 2.7.6.3: https://iubmb.qmul.ac.uk/enzyme/EC2/7/6/3.html; ",
      "Williams et al. Nat Metab. 2026;8:924-940. doi:10.1038/s42255-026-01489-w"
    ),
    "Current-cohort convergent-core and UA-correlation evidence",
    paste0(
      "IUBMB EC 2.1.1.13: https://iubmb.qmul.ac.uk/enzyme/EC2/1/1/13.html; ",
      "Williams et al. Nat Metab. 2026;8:924-940. doi:10.1038/s42255-026-01489-w"
    )
  ),
  stringsAsFactors = FALSE
)

contrast_names <- c("E1vE2", "E1vE3", "E4vE2", "E4vE3")
for (contrast in contrast_names) {
  candidates[[paste0(contrast, "_log2FC")]] <- NA_real_
  candidates[[paste0(contrast, "_FDR")]] <- NA_real_
}
candidates$envfit_R2 <- NA_real_
candidates$envfit_P <- NA_real_
candidates$UA_partial_rho <- NA_real_
candidates$UA_partial_FDR <- NA_real_

fill_contrasts <- function(i, row) {
  for (contrast in contrast_names) {
    candidates[i, paste0(contrast, "_log2FC")] <<-
      as.numeric(row[[paste0("log2FoldChange_", contrast)]])
    candidates[i, paste0(contrast, "_FDR")] <<-
      as.numeric(row[[paste0("padj_", contrast)]])
  }
}

for (i in seq_len(nrow(candidates))) {
  feature <- candidates$source_feature[i]
  category <- candidates$category[i]

  if (category == "Metabolite") {
    rows <- limma_all[limma_all$feature == feature, , drop = FALSE]
    if (nrow(rows) != 4) stop("Expected four metabolite contrasts for ", feature)
    for (contrast in contrast_names) {
      row <- rows[rows$comparison == contrast, , drop = FALSE]
      candidates[i, paste0(contrast, "_log2FC")] <- as.numeric(row$log2FoldChange)
      candidates[i, paste0(contrast, "_FDR")] <- as.numeric(row$padj)
    }
    ef <- envfit[envfit$feature == feature, , drop = FALSE]
    if (nrow(ef) != 1) stop("Missing curated envfit result for ", feature)
    candidates$envfit_R2[i] <- as.numeric(ef$r2)
    candidates$envfit_P[i] <- as.numeric(ef$pval)
  } else if (category == "Pathway") {
    row <- rbind(path_e1, path_e4)
    row <- row[row$feature == feature, , drop = FALSE]
    row <- row[!duplicated(row$feature), , drop = FALSE]
    if (nrow(row) != 1) stop("Missing core-pathway result for ", feature)
    fill_contrasts(i, row)
  } else {
    row <- ec_e1[ec_e1$feature == feature, , drop = FALSE]
    if (nrow(row) != 1 || !feature %in% ec_e4$feature) {
      stop("Missing shared-core EC result for ", feature)
    }
    fill_contrasts(i, row)
  }

  if (category != "Metabolite") {
    pc <- partial_cor[partial_cor$feature == feature, , drop = FALSE]
    if (nrow(pc) != 1) stop("Missing UA partial correlation for ", feature)
    candidates$UA_partial_rho[i] <- as.numeric(pc$rho)
    candidates$UA_partial_FDR[i] <- as.numeric(pc$padj)
  }
}

final_columns <- c(
  "candidate_id", "mediator_name", "category", "chemical_class",
  "source_feature", "source_layer", "selection_basis", "effect_size_role",
  "E1vE2_log2FC", "E1vE2_FDR", "E1vE3_log2FC", "E1vE3_FDR",
  "E4vE2_log2FC", "E4vE2_FDR", "E4vE3_log2FC", "E4vE3_FDR",
  "envfit_R2", "envfit_P", "UA_partial_rho", "UA_partial_FDR",
  "biological_rationale", "literature_support"
)
final <- candidates[, final_columns]

stopifnot(
  nrow(final) == 22,
  length(unique(final$candidate_id)) == 22,
  sum(final$category == "Metabolite") == 9,
  sum(final$category == "Pathway") == 7,
  sum(final$category == "EC") == 6,
  all(!is.na(final$envfit_R2[final$category == "Metabolite"])),
  all(!is.na(final$UA_partial_rho[final$category != "Metabolite"]))
)

format_stat <- function(x) {
  vapply(
    x,
    function(value) {
      if (is.na(value)) return("NA")
      format(signif(value, 3), scientific = FALSE, trim = TRUE)
    },
    character(1)
  )
}

compact <- data.frame(
  candidate = final$mediator_name,
  layer = final$category,
  upstream_source = final$source_layer,
  inclusion_basis = final$selection_basis,
  statistical_support = ifelse(
    final$category == "Metabolite",
    paste0(
      ifelse(grepl("Type 1 core", final$source_layer),
             "Type 1 core; ", "Not an enterotype core feature; "),
      "envfit R2=", format_stat(final$envfit_R2),
      "; P=", format_stat(final$envfit_P)
    ),
    paste0(
      ifelse(final$category == "Pathway", "Core pathway; ",
             "Shared convergent core EC; "),
      "UA partial rho=", format_stat(final$UA_partial_rho),
      "; FDR=", format_stat(final$UA_partial_FDR)
    )
  ),
  effect_size_role = final$effect_size_role,
  biological_rationale = final$biological_rationale,
  reference_or_database = final$literature_support,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(nrow(compact) == 22, ncol(compact) == 8)

dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

detail_temp_path <- paste0(detailed_output_path, ".tmp")
write.csv(final, detail_temp_path, row.names = FALSE, na = "NA")
detailed_generated <- read.csv(
  detail_temp_path, stringsAsFactors = FALSE, check.names = FALSE
)
stopifnot(
  nrow(detailed_generated) == 22,
  ncol(detailed_generated) == length(final_columns)
)
if (!file.copy(detail_temp_path, detailed_output_path, overwrite = TRUE)) {
  stop("Failed to replace detailed output: ", detailed_output_path)
}
unlink(detail_temp_path)

temp_path <- paste0(output_path, ".tmp")
write.csv(compact, temp_path, row.names = FALSE, na = "NA")

generated <- read.csv(temp_path, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(nrow(generated) == 22, ncol(generated) == 8)
if (!file.copy(temp_path, output_path, overwrite = TRUE)) {
  stop("Failed to replace output: ", output_path)
}
unlink(temp_path)

cat("Table S4 build PASS\n")
cat("Output:", output_path, "\n")
cat("Dimensions:", nrow(compact), "rows x", ncol(compact), "columns\n")
cat("Detailed inspection output:", detailed_output_path, "\n")
cat("Metabolites with Type 1 core support: 7/9\n")
cat("Pathways with core differential support: 7/7\n")
cat("ECs with shared convergent core support: 6/6\n")
cat("Pathways with UA partial-correlation FDR < 0.05:",
    sum(final$UA_partial_FDR[final$category == "Pathway"] < 0.05), "/7\n")
cat("ECs with UA partial-correlation FDR < 0.05:",
    sum(final$UA_partial_FDR[final$category == "EC"] < 0.05), "/6\n")
