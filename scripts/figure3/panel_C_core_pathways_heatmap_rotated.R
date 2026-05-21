# ==============================================================================
# Figure 3C Alternative: Rotated Comprehensive Heatmap
#
# Shows ALL core pathways (~95) but rotated 90 degrees.
# Rows: Enterotypes (Type 1, Type 4)
# Columns: Pathways, grouped horizontally by Theme.
# ==============================================================================

library(dplyr)
library(readr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

base_dir <- here::here()
results_dir <- file.path(base_dir, "results/tables/section_3_3")
out_dir <- file.path(base_dir, "results/figures/panels")

# --- 1. Load Data ---
pw_e1 <- read_csv(file.path(results_dir, "pathway_E1_core.csv"), show_col_types = FALSE)
pw_e4 <- read_csv(file.path(results_dir, "pathway_E4_core.csv"), show_col_types = FALSE)

# --- 2. Compute mean LFC ---
pw_e1 <- pw_e1 %>%
  mutate(
    E1_mean_lfc = (log2FoldChange_E1vE2 + log2FoldChange_E1vE3) / 2,
    E4_mean_lfc = NA_real_,
    pw_id = str_extract(feature, "^[^:]+"),
    pw_name = str_remove(feature, "^[^:]+:\\s*")
  )

pw_e4 <- pw_e4 %>%
  mutate(
    E1_mean_lfc = NA_real_,
    E4_mean_lfc = (log2FoldChange_E4vE2 + log2FoldChange_E4vE3) / 2,
    pw_id = str_extract(feature, "^[^:]+"),
    pw_name = str_remove(feature, "^[^:]+:\\s*")
  )

# --- 3. Assign functional themes ---
assign_theme <- function(pw_id, pw_name) {
  case_when(
    str_detect(pw_id, "ARGORNPROST|ARGSYN|ARGSYNBSUB|GLUTORN|PWY-5004|ORNARGDEG") |
      str_detect(pw_name, regex("arginine|ornithine|proline", ignore_case = TRUE)) ~ "Arginine / Ornithine metabolism",
    str_detect(pw_id, "TCA|GLUDEG|PWY-5690|P105-PWY") |
      str_detect(pw_name, regex("TCA|citrate|GABA shunt", ignore_case = TRUE)) ~ "TCA cycle / GABA shunt",
    str_detect(pw_id, "PWY-5100|P108-PWY|PROPFERM") |
      str_detect(pw_name, regex("ferment|propanoate|pyruvate", ignore_case = TRUE)) ~ "Fermentation",
    str_detect(pw_name, regex("glutam(ate|ine)|glutam", ignore_case = TRUE)) &
      !str_detect(pw_name, regex("arginine|ornithine", ignore_case = TRUE)) ~ "Glutamate / Glutamine",
    str_detect(pw_id, "VALSYN|ILEUSYN|LEUSYN|BRANCHED|HISTSYN|PWY-2942|PWY-5097|PWY-5101|ASPASN|HOMOSER|DAPLYSINE|TRPSYN") |
      str_detect(pw_name, regex("amino acid biosyn|valine|isoleucine|leucine|histidine|lysine|threonine|asparagine|tryptophan|branched", ignore_case = TRUE)) ~ "Amino acid biosynthesis",
    str_detect(pw_id, "PEPTIDOGLYCAN|PWY-6385|PWY-6386|PWY-6387|PWY0-1586") |
      str_detect(pw_name, regex("peptidoglycan|murein", ignore_case = TRUE)) ~ "Peptidoglycan biosynthesis",
    str_detect(pw_id, "PWY-6122|PWY-6277|PWY-7219|PWY-7221|PWY-7222|PWY-5695|PWY-6608|PWY-7234|PWY-6126|PWY-6124|PWY-6125") |
      str_detect(pw_name, regex("purine|pyrimidine|inosine|urate|nucleotide", ignore_case = TRUE)) ~ "Purine / Pyrimidine / Urate",
    str_detect(pw_id, "CATECHOL|PWY-5181|PWY-5431|PWY-6182|PWY-6185") |
      str_detect(pw_name, regex("catechol|toluene|cresol|aminophenol|salicylate|aromatic", ignore_case = TRUE)) ~ "Aromatic compound degradation",
    str_detect(pw_id, "DTDPRHAMSYN|PWY-1269|PWY-6749") |
      str_detect(pw_name, regex("sugar nucleotide|rhamnose|CMP-pseudo", ignore_case = TRUE)) ~ "Sugar nucleotide biosynthesis",
    str_detect(pw_id, "PWY-5845|PWY-5862|PWY-5896|PWY-5898|PWY-5899") |
      str_detect(pw_name, regex("menaquinol|ubiquinol|quinol", ignore_case = TRUE)) ~ "Quinol biosynthesis",
    str_detect(pw_name, regex("folate|tetrahydrofolate", ignore_case = TRUE)) ~ "Folate metabolism",
    str_detect(pw_name, regex("CoA|fatty acid|lipid", ignore_case = TRUE)) ~ "Lipid / CoA metabolism",
    TRUE ~ "Other"
  )
}

pw_e1$theme <- assign_theme(pw_e1$pw_id, pw_e1$pw_name)
pw_e4$theme <- assign_theme(pw_e4$pw_id, pw_e4$pw_name)

# Make sure we don't have exactly duplicated pathway names if they appear in both (not core overlap, but just in case)
all_pw <- bind_rows(
  pw_e1 %>% select(pw_id, pw_name, theme, E1_mean_lfc, E4_mean_lfc),
  pw_e4 %>% select(pw_id, pw_name, theme, E1_mean_lfc, E4_mean_lfc)
) %>%
  mutate(pw_name = str_to_sentence(pw_name)) %>%
  mutate(pw_name = str_replace(pw_name, "\\(meso-diaminopimelate containing\\)", "")) %>%
  mutate(pw_name = str_replace_all(pw_name, "Inosine-5-phosphate biosynthesis iii", "IMP biosynthesis III")) %>%
  mutate(pw_name = str_replace_all(pw_name, "Urate biosynthesis/inosine 5-phosphate degradation", "Urate bio / IMP deg")) %>%
  # Truncate very long names for columns removed, use full names
  mutate(pw_label = pw_name) %>%
  distinct(pw_id, .keep_all=TRUE)  # guarantee uniqueness by id

# manuscript named for bolding
manuscript_named <- c(
  "arginine", "ornithine", "proline", "citrulline", "IMP biosynthesis",
  "guanosine", "urate bio", "valine", "isoleucine", "branched",
  "histidine", "lysine", "peptidoglycan", "aminoimidazole",
  "adenosine", "pyrimidine"
)
all_pw$is_bold <- grepl(paste(tolower(manuscript_named), collapse="|"), tolower(all_pw$pw_name))

# --- 4. Matrix transposing (Rotated) ---
# Original: rows = pathways, cols = E1, E4
mat <- as.matrix(all_pw[, c("E1_mean_lfc", "E4_mean_lfc")])
rownames(mat) <- all_pw$pw_label
colnames(mat) <- c("Type 1", "Type 4")

# ROTATE:
mat_t <- t(mat)
# Now rows = Type 1, Type 2. Cols = Pathways
# Columns grouped by theme.
theme_order <- c(
  "Arginine / Ornithine metabolism", "TCA cycle / GABA shunt",
  "Fermentation", "Glutamate / Glutamine",
  "Purine / Pyrimidine / Urate", "Amino acid biosynthesis",
  "Peptidoglycan biosynthesis", "Aromatic compound degradation",
  "Sugar nucleotide biosynthesis", "Quinol biosynthesis",
  "Folate metabolism", "Lipid / CoA metabolism", "Other"
)
all_pw$theme_f <- factor(all_pw$theme, levels = theme_order)
col_order <- order(all_pw$theme_f)

mat_t <- mat_t[, col_order]
col_themes <- all_pw$theme_f[col_order]
bold_cols <- all_pw$is_bold[col_order]

theme_colors <- c(
  "Arginine / Ornithine metabolism" = "#E64B35",
  "TCA cycle / GABA shunt" = "#F39C12",
  "Fermentation" = "#E67E22",
  "Glutamate / Glutamine" = "#D35400",
  "Purine / Pyrimidine / Urate" = "#8E44AD",
  "Amino acid biosynthesis" = "#3498DB",
  "Peptidoglycan biosynthesis" = "#1ABC9C",
  "Aromatic compound degradation" = "#27AE60",
  "Sugar nucleotide biosynthesis" = "#2ECC71",
  "Quinol biosynthesis" = "#16A085",
  "Folate metabolism" = "#2980B9",
  "Lipid / CoA metabolism" = "#7F8C8D",
  "Other" = "#BDC3C7"
)

# --- 5. Annotations ---
col_ha <- HeatmapAnnotation(
  Theme = as.character(col_themes),
  col = list(Theme = theme_colors),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  annotation_legend_param = list(
    title = "Functional Theme",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 9)
  )
)

row_ha <- rowAnnotation(
  Enterotype = c("Type 1", "Type 4"),
  col = list(Enterotype = c("Type 1" = "#E64B35", "Type 4" = "#3C5488")),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  width = unit(4, "mm")
)

max_abs <- max(abs(mat_t), na.rm = TRUE)
col_fun <- colorRamp2(
  c(-max_abs, -max_abs/2, 0, max_abs/2, max_abs),
  c("#313695", "#74add1", "#f5f5f5", "#f46d43", "#a50026")
)

ht <- Heatmap(
  mat_t,
  name = "Mean\nlog2FC",
  col = col_fun,
  na_col = "#EDEDED",
  rect_gp = gpar(col = "white", lwd = 1),
  
  # Split columns
  column_split = col_themes,
  column_gap = unit(3, "mm"),
  column_title_rot = 90, # Vertically oriented grouping labels
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  
  # Column names
  column_names_side = "bottom",
  column_names_rot = 60,
  column_names_gp = gpar(
    fontsize = 8,
    fontface = ifelse(bold_cols, "bold", "plain")
  ),
  column_names_max_height = max_text_width(colnames(mat_t), gp = gpar(fontsize = 8, fontface = "bold")) * 1.2,
  
  # Row settings
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 12, fontface = "bold"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  
  # Annotations
  top_annotation = col_ha,
  left_annotation = row_ha,
  
  heatmap_legend_param = list(
    title = expression("Mean log"[2]*"FC"),
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    legend_direction = "vertical"
  ),
  height = unit(2, "cm") # The actual cells are just 2 rows. Rest of height is column names.
)

pdf(file.path(out_dir, "PanelC_Rotated_Heatmap.pdf"), width = 20, height = 10)
draw(ht, padding = unit(c(5, 5, 40, 5), "mm"))
dev.off()

png(file.path(out_dir, "PanelC_Rotated_Heatmap.png"), width = 20, height = 10, units = "in", res = 300, bg="white")
draw(ht, padding = unit(c(5, 5, 40, 5), "mm"))
dev.off()

cat("Rotated heatmap completed.\n")
