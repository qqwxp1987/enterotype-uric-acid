# Build Supplementary_Tables.xlsx (multi-sheet) from compiled CSV tables.
suppressMessages(library(openxlsx))

root <- "C:/Users/qqwxp/OneDrive/Work/HUA_FCTJ2021"
sup  <- file.path(root, "supplementary", "tables")
s34  <- file.path(root, "results", "tables", "section_3_4")

specs <- list(
  list(sheet="S1",  path=file.path(sup,"Table_S1_core_features.csv"),
       title="Table S1. Independent core multi-omics features by enterotype (layer-count summary)."),
  list(sheet="S1a", path=file.path(sup,"Table_S1a_taxa_E1_only.csv"),
       title="Table S1a. Type 1-exclusive core species."),
  list(sheet="S1b", path=file.path(sup,"Table_S1b_taxa_Shared.csv"),
       title="Table S1b. Shared core species (Type 1 and Type 4)."),
  list(sheet="S1c", path=file.path(sup,"Table_S1c_taxa_E4_only.csv"),
       title="Table S1c. Type 4-exclusive core species."),
  list(sheet="S1d", path=file.path(sup,"Table_S1d_pathway_E1_only.csv"),
       title="Table S1d. Type 1-exclusive core pathways."),
  list(sheet="S1e", path=file.path(sup,"Table_S1e_pathway_E4_only.csv"),
       title="Table S1e. Type 4-exclusive core pathways."),
  list(sheet="S1f", path=file.path(sup,"Table_S1f_metabolome_E1_only.csv"),
       title="Table S1f. Type 1 core serum metabolites."),
  list(sheet="S1g", path=file.path(s34,"ec_E1_only.csv"),
       title="Table S1g. Type 1-exclusive core enzymes (EC)."),
  list(sheet="S1h", path=file.path(s34,"ec_E4_only.csv"),
       title="Table S1h. Type 4-exclusive core enzymes (EC)."),
  list(sheet="S1i", path=file.path(s34,"ec_Shared.csv"),
       title="Table S1i. Shared core enzymes (EC; convergent and divergent)."),
  list(sheet="S2",  path=file.path(sup,"Table_S2_mediation_main.csv"),
       title="Table S2. Causal mediation analysis - main results (Type 1 and Type 4 models)."),
  list(sheet="S3",  path=file.path(sup,"Table_S3_mediation_binary_sensitivity.csv"),
       title="Table S3. Mediation sensitivity analysis - binary hyperuricemia outcome."),
  list(sheet="S4",  path=file.path(sup,"Table_S4_mediation_refgroup_sensitivity.csv"),
       title="Table S4. Mediation sensitivity analysis - reference-group split (Type 2 vs Type 3)."),
  list(sheet="S5",  path=file.path(sup,"Table_S5_mediation_medsens.csv"),
       title="Table S5. Mediation sensitivity analysis - sequential-ignorability (medsens)."),
  list(sheet="S6",  path=file.path(sup,"Table_S6_type4_heterogeneity.csv"),
       title="Table S6. Within-Type 4 heterogeneity - HUA+ vs HUA- comparison.")
)

wb <- createWorkbook()
addWorksheet(wb, "Contents")
contents <- data.frame(Sheet=character(), Title=character(), Rows=integer(), stringsAsFactors=FALSE)

for (sp in specs) {
  if (!file.exists(sp$path)) { cat("MISSING:", sp$path, "\n"); next }
  d <- read.csv(sp$path, check.names=FALSE, stringsAsFactors=FALSE)
  addWorksheet(wb, sp$sheet)
  writeData(wb, sp$sheet, sp$title, startRow=1, startCol=1)
  writeData(wb, sp$sheet, d, startRow=3, startCol=1)
  setColWidths(wb, sp$sheet, cols=1:max(1,ncol(d)), widths="auto")
  contents <- rbind(contents, data.frame(Sheet=sp$sheet, Title=sp$title,
                                          Rows=nrow(d), stringsAsFactors=FALSE))
  cat("OK:", sp$sheet, "-", nrow(d), "rows\n")
}

writeData(wb, "Contents",
          "Supplementary Tables - Gut Microbiome Enterotypes Associate with Serum Uric Acid",
          startRow=1)
writeData(wb, "Contents", contents, startRow=3)
setColWidths(wb, "Contents", cols=1:3, widths=c(10, 95, 8))

out <- file.path(sup, "Supplementary_Tables.xlsx")
saveWorkbook(wb, out, overwrite=TRUE)
cat("SAVED:", out, "\n")
