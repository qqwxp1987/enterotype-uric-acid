# ==============================================================================
# R2 Minor Comment 4 — Rebuild Figure 4 (publishable, vector PDF)
#
# Confirmed layout:
#   * Original R1 landscape canvas: 11 x 8.5 in (792 x 612 pt).
#   * Old KEGG schematic (panel D) moves to Additional file 8: Figure S7.
#   * Left two-thirds: A and B in the upper row; C across the lower row.
#   * Right one-third: old E becomes new D, with three single-column facets.
#   * Panel B retains the per-theme pathway counts cited in the figure legend.
#   * Panel C retains the full, wrapped names of EC 2.7.6.3 and EC 2.3.1.129.
#
# Data and statistical definitions are unchanged from the R1 scripts; only
# presentation, typography, and panel geometry are revised.
#
# Outputs:
#   投稿/.../第二轮返修和处理/Figure4.pdf  (vector, cairo)
#   投稿/.../第二轮返修和处理/Figure4.png  (300 dpi preview)
#   投稿/.../第二轮返修和处理/Minor/figure4_rebuild/panels/Figure4[A-D].pdf
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(data.table)
  library(patchwork)
  library(here)
})

PROJ <- here::here()
OUT_DIR <- file.path(PROJ, "results", "figures")
PANEL_DIR <- file.path(OUT_DIR, "panels")
dir.create(PANEL_DIR, showWarnings = FALSE, recursive = TRUE)

# --- Typography (canvas pt; effective at 510 pt width = canvas x 0.644) --------
.PT      <- 72.27 / 25.4
PT_LETTER       <- 14
PT_STRIP        <- 10
PT_AXIS_TITLE   <- 10
PT_AXIS_TEXT    <- 9
PT_LEGEND_TITLE <- 10
PT_LEGEND_TEXT  <- 9
PT_EC_LABEL     <- 9
GEOM_PRIMARY    <- 9 / .PT
GEOM_SPECIES    <- 8.5 / .PT

FIG4_BG  <- "#FAFBFC"
COL_E1   <- "#E64B35"
COL_E4   <- "#3C5488"
COL_SHARED <- "#FFB347"

# --- Geometry ------------------------------------------------------------------
FIG4_W     <- 11.0
FIG4_H     <- 8.5
FIG4_LEFT  <- 8 / 12 * FIG4_W
FIG4_RIGHT <- 4 / 12 * FIG4_W
H_AB       <- 3 / 10 * FIG4_H
H_C        <- 7 / 10 * FIG4_H
FIG4_DPI   <- 300

# ==============================================================================
# Panel A — core feature summary (butterfly bar chart)
# ==============================================================================
build_panel_A <- function() {
  summary_df <- read.csv(file.path(PROJ, "results", "tables", "section_3_4",
                                   "independent_core_summary.csv"))

  # QC: values must match R1 figure (72/43/296/24 | 2/52/371/0 | shared 1/0/36/0)
  stopifnot(nrow(summary_df) == 4)
  qc <- summary_df %>% arrange(match(data_type, c("Taxa", "Pathway", "EC", "Metabolome")))
  stopifnot(
    all(qc$E1_only == c(72, 43, 296, 24)),
    all(qc$E4_only == c(2, 52, 371, 0)),
    all(qc$Shared  == c(1, 0, 36, 0))
  )

  plot_df <- summary_df %>%
    mutate(
      data_type = factor(data_type, levels = rev(c("Taxa", "Pathway", "EC", "Metabolome"))),
      E1_plot = -E1_only,
      E4_plot = E4_only
    )
  max_val <- max(c(summary_df$E1_only, summary_df$E4_only)) * 1.15

  ggplot(plot_df, aes(y = data_type)) +
    annotate("rect", xmin = -max_val, xmax = 0, ymin = 0.4, ymax = 4.6,
             fill = COL_E1, alpha = 0.03) +
    annotate("rect", xmin = 0, xmax = max_val, ymin = 0.4, ymax = 4.6,
             fill = COL_E4, alpha = 0.03) +
    geom_col(aes(x = E1_plot), fill = COL_E1, width = 0.55, color = NA) +
    geom_col(aes(x = E4_plot), fill = COL_E4, width = 0.55, color = NA) +
    # CHANGE: capped centre line (was infinite geom_vline -> crossed "Shared")
    geom_segment(x = 0, xend = 0, y = 0.5, yend = 4.52,
                 linewidth = 0.6, color = "grey40") +
    geom_text(aes(x = pmin(E1_plot, -58), label = E1_only),
              hjust = 1.1, size = GEOM_PRIMARY, fontface = "bold", color = COL_E1) +
    geom_text(aes(x = pmax(E4_plot, 58), label = E4_only),
              hjust = -0.1, size = GEOM_PRIMARY, fontface = "bold", color = COL_E4) +
    geom_point(aes(x = 0),
               size = ifelse(plot_df$Shared > 0, 6, 3.6),
               shape = 21, fill = ifelse(plot_df$Shared > 0, COL_SHARED, "grey85"),
               color = ifelse(plot_df$Shared > 0, "grey30", "grey70"), stroke = 0.5) +
    geom_text(aes(x = 0, label = Shared),
              size = GEOM_PRIMARY, fontface = "bold",
              color = ifelse(plot_df$Shared > 0, "grey20", "grey50")) +
    scale_x_continuous(
      limits = c(-max_val, max_val),
      breaks = seq(-400, 400, by = 200),
      labels = function(x) abs(x),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    annotate("text", x = -max_val * 0.65, y = 4.92, label = "Type 1",
             fontface = "bold", size = 9 / .PT, color = COL_E1, hjust = 0.5) +
    annotate("text", x = max_val * 0.65, y = 4.92, label = "Type 4",
             fontface = "bold", size = 9 / .PT, color = COL_E4, hjust = 0.5) +
    annotate("text", x = 0, y = 4.92, label = "Shared",
             fontface = "bold.italic", size = 9 / .PT, color = "grey40", hjust = 0.5) +
    labs(x = "Number of core features", y = NULL) +
    theme_minimal(base_size = PT_AXIS_TITLE, base_family = "") +
    theme(
      plot.background  = element_rect(fill = FIG4_BG, color = NA),
      panel.background = element_rect(fill = FIG4_BG, color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey88", linewidth = 0.25, linetype = "dotted"),
      axis.text.y   = element_text(size = PT_AXIS_TEXT, face = "bold", color = "grey25",
                                   margin = margin(r = 5)),
      axis.text.x   = element_text(size = PT_AXIS_TEXT, color = "grey45"),
      axis.title.x  = element_text(size = PT_AXIS_TITLE, color = "grey40",
                                   margin = margin(t = 8)),
      plot.margin   = margin(t = 18, r = 16, b = 6, l = 8)
    ) +
    coord_cartesian(clip = "off", ylim = c(0.5, 4.5))
}

# ==============================================================================
# Panel B — core pathways theme strip plot (per-theme counts retained)
# ==============================================================================
build_panel_B <- function() {
  results_dir <- file.path(PROJ, "results", "tables", "section_3_4")
  pw_e1 <- read.csv(file.path(results_dir, "pathway_E1_core.csv"), check.names = FALSE)
  pw_e4 <- read.csv(file.path(results_dir, "pathway_E4_core.csv"), check.names = FALSE)

  pw_e1 <- pw_e1 %>%
    mutate(
      E1_mean_lfc = (log2FoldChange_E1vE2 + log2FoldChange_E1vE3) / 2,
      E4_mean_lfc = NA_real_,
      pw_id   = str_extract(feature, "^[^:]+"),
      pw_name = str_remove(feature, "^[^:]+:\\s*")
    )
  pw_e4 <- pw_e4 %>%
    mutate(
      E1_mean_lfc = NA_real_,
      E4_mean_lfc = (log2FoldChange_E4vE2 + log2FoldChange_E4vE3) / 2,
      pw_id   = str_extract(feature, "^[^:]+"),
      pw_name = str_remove(feature, "^[^:]+:\\s*")
    )

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

  all_pw <- bind_rows(
    pw_e1 %>% select(pw_id, pw_name, theme, E1_mean_lfc, E4_mean_lfc),
    pw_e4 %>% select(pw_id, pw_name, theme, E1_mean_lfc, E4_mean_lfc)
  ) %>% distinct(pw_id, .keep_all = TRUE)

  theme_order <- c(
    "Arginine / Ornithine metabolism", "TCA cycle / GABA shunt",
    "Fermentation", "Glutamate / Glutamine",
    "Purine / Pyrimidine / Urate", "Amino acid biosynthesis",
    "Peptidoglycan biosynthesis", "Aromatic compound degradation",
    "Sugar nucleotide biosynthesis", "Quinol biosynthesis",
    "Folate metabolism", "Lipid / CoA metabolism", "Other"
  )

  strip_df <- all_pw %>%
    mutate(
      enterotype = ifelse(!is.na(E1_mean_lfc), "Type 1", "Type 4"),
      mean_lfc   = dplyr::coalesce(E1_mean_lfc, E4_mean_lfc),
      theme_f    = factor(theme, levels = rev(theme_order))
    )

  # QC: 43 Type 1 + 52 Type 4 = 95 pathways, 13 themes
  stopifnot(nrow(strip_df) == 95, length(theme_order) == 13,
            sum(strip_df$enterotype == "Type 1") == 43,
            sum(strip_df$enterotype == "Type 4") == 52)

  lev <- levels(strip_df$theme_f)
  theme_n <- strip_df %>% count(theme_f)
  theme_lab <- setNames(
    sprintf("%s  (n = %d)", lev, theme_n$n[match(lev, theme_n$theme_f)]),
    lev
  )

  n_theme <- length(lev)
  bands <- data.frame(ymin = seq(1.5, n_theme - 0.5, by = 2)) %>%
    mutate(ymax = ymin + 1)

  set.seed(42)  # same jitter seed as R1 -> identical point placement
  ggplot(strip_df, aes(x = mean_lfc, y = theme_f)) +
    geom_rect(data = bands, inherit.aes = FALSE,
              aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
              fill = "grey90", alpha = 0.45) +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey55") +
    geom_jitter(aes(color = enterotype), height = 0.30, width = 0,
                size = 1.5, alpha = 0.78) +
    scale_color_manual(values = c("Type 1" = COL_E1, "Type 4" = COL_E4), name = NULL) +
    scale_y_discrete(labels = theme_lab) +
    scale_x_continuous(expand = expansion(mult = 0.03)) +
    labs(x = expression("Mean log"[2] * " fold-change"), y = NULL) +
    guides(color = guide_legend(override.aes = list(size = 2.2))) +
    theme_minimal(base_size = PT_AXIS_TITLE) +
    theme(
      plot.background  = element_rect(fill = FIG4_BG, color = NA),
      panel.background = element_rect(fill = FIG4_BG, color = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey85", linewidth = 0.25),
      axis.text.y   = element_text(size = PT_AXIS_TEXT, color = "grey15"),
      axis.text.x   = element_text(size = PT_AXIS_TEXT, color = "grey45"),
      axis.title.x  = element_text(size = PT_AXIS_TITLE, color = "grey40",
                                   margin = margin(t = 5)),
      legend.position = "top",
      legend.justification = "right",
      legend.text = element_text(size = PT_LEGEND_TEXT),
      legend.key.size = unit(9, "pt"),
      legend.margin = margin(t = 0, b = 0),
      plot.margin = margin(t = 4, r = 12, b = 4, l = 6)
    )
}

# ==============================================================================
# Panel C — shared core ECs dumbbell (bigger pitch, 8pt labels, 2 full names)
# ==============================================================================
build_panel_C <- function() {
  ec_dir <- read.csv(file.path(PROJ, "results", "tables", "section_3_4",
                               "ec_shared_direction.csv"))
  ec_kegg <- read.csv(file.path(PROJ, "results", "tables", "section_3_4",
                                "ec_kegg_names.csv"))

  # QC: 36 shared core ECs (22 divergent, 14 convergent)
  stopifnot(nrow(ec_dir) == 36,
            sum(ec_dir$direction_type == "Convergent") == 14,
            sum(ec_dir$direction_type == "Divergent") == 22,
            all(c("2.7.6.3", "2.3.1.129") %in% ec_dir$feature))

  manuscript_ecs <- c("2.7.6.3", "2.6.1.1", "2.1.1.13")

  # CHANGE: full KEGG names (wrapped at a hyphen/bond boundary) for the two
  # enzymes whose names were truncated in R1; plain text + "\n" (element_text
  # honours the line break; markdown renderers would mangle the KEGG "---"
  # bond notation into em-dashes, so ggtext is deliberately NOT used)
  full_wrapped <- c(
    "2.7.6.3"   = "2.7.6.3  2-amino-4-hydroxy-6-hydroxymethyl-\ndihydropteridine diphosphokinase",
    "2.3.1.129" = "2.3.1.129  acyl-[acyl-carrier-protein]---\nUDP-N-acetylglucosamine O-acyltransferase"
  )

  ec_dir <- ec_dir %>%
    left_join(ec_kegg, by = c("feature" = "ec")) %>%
    mutate(
      name = case_when(
        is.na(name) ~ NA_character_,
        str_detect(name, "^Transferred") ~ NA_character_,
        TRUE ~ name
      ),
      ec_label = case_when(
        feature %in% names(full_wrapped) ~ full_wrapped[feature],
        !is.na(name) ~ paste0(feature, "  ", str_trunc(name, 52, side = "right")),
        TRUE ~ feature
      ),
      is_named = feature %in% manuscript_ecs
    ) %>%
    mutate(
      delta = abs(E4_mean_lfc - E1_mean_lfc),
      direction_type = factor(direction_type, levels = c("Divergent", "Convergent"))
    ) %>%
    arrange(direction_type, delta)

  # Phantom spacer rows (invisible labels, no dots) cushion the two wrapped
  # two-line labels so their second lines never touch neighbouring rows
  lv <- as.character(ec_dir$ec_label)
  w1 <- full_wrapped[["2.7.6.3"]]
  w2 <- full_wrapped[["2.3.1.129"]]
  ph_levels <- c(" ", "  ", "   ")
  lv_out <- c()
  for (l in lv) {
    if (l == w1) lv_out <- c(lv_out, ph_levels[1])  # cushion below 2.7.6.3
    if (l == w2) lv_out <- c(lv_out, ph_levels[2])  # cushion below 2.3.1.129
    lv_out <- c(lv_out, l)
    if (l == w2) lv_out <- c(lv_out, ph_levels[3])  # cushion above 2.3.1.129
  }
  ph_rows <- data.frame(
    feature = paste0("phantom", 1:3),
    direction_type = "Convergent",
    E1_mean_lfc = NA_real_, E4_mean_lfc = NA_real_,
    name = NA_character_, ec_label = ph_levels,
    is_named = FALSE, delta = NA_real_,
    stringsAsFactors = FALSE
  )
  ec_dir <- bind_rows(ec_dir, ph_rows) %>%
    mutate(
      direction_type = factor(as.character(direction_type),
                              levels = c("Divergent", "Convergent")),
      ec_label = factor(ec_label, levels = lv_out)
    )

  ec_long <- ec_dir %>%
    select(ec_label, E1_mean_lfc, E4_mean_lfc, direction_type, is_named) %>%
    pivot_longer(cols = c(E1_mean_lfc, E4_mean_lfc),
                 names_to = "enterotype", values_to = "mean_lfc",
                 values_drop_na = TRUE) %>%
    mutate(enterotype = ifelse(enterotype == "E1_mean_lfc", "Type 1", "Type 4"))

  ggplot() +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey55") +
    geom_segment(
      data = ec_dir,
      aes(x = E1_mean_lfc, xend = E4_mean_lfc, y = ec_label, yend = ec_label),
      linewidth = ifelse(ec_dir$is_named, 1.1, 0.6),
      color = ifelse(ec_dir$direction_type == "Divergent", "#E2C48E", "#7FB97C"),
      na.rm = TRUE
    ) +
    geom_point(
      data = ec_long,
      aes(x = mean_lfc, y = ec_label, color = enterotype),
      size = ifelse(ec_long$is_named, 2.6, 1.7), alpha = 0.95
    ) +
    scale_color_manual(values = c("Type 1" = COL_E1, "Type 4" = COL_E4), name = NULL) +
    facet_grid(direction_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(x = expression("Mean log"[2] * " fold-change"), y = NULL) +
    # NOTE: manual theme instead of theme_minimal() — verified equivalent
    # rendering; keeps panel C independent of the ggplot2 4.0 / ggtext
    # interaction (no markdown elements are used any more, but this theme is
    # self-contained either way).
    theme(
      panel.background = element_rect(fill = FIG4_BG, color = NA),
      plot.background  = element_rect(fill = FIG4_BG, color = NA),
      panel.border = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.25),
      panel.spacing.y = unit(5, "pt"),
      # CHANGE: 8pt enzyme labels; phantom spacer rows give the two wrapped
      # two-line labels room (uniform visual row pitch, no collisions)
      axis.text.y   = element_text(size = PT_EC_LABEL, color = "grey20",
                                   lineheight = 0.8),
      axis.text.x   = element_text(size = PT_AXIS_TEXT, color = "grey45"),
      axis.title.x  = element_text(size = PT_AXIS_TITLE, margin = margin(t = 5)),
      axis.ticks    = element_line(color = "grey60", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(size = PT_STRIP, face = "bold", angle = 90),
      legend.position = "top",
      legend.justification = "center",
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(size = PT_LEGEND_TEXT),
      legend.key.size = unit(9, "pt"),
      legend.margin = margin(t = 0, b = 2),
      plot.margin = margin(t = 4, r = 14, b = 4, l = 4)
    )
}

# ==============================================================================
# Panel D (old E) — species attribution for the 3 convergent ECs
# ==============================================================================
build_panel_D <- function() {
  attr_dir <- file.path(PROJ, "results", "tables", "section_3_4", "attribution")

  prep <- readRDS(file.path(attr_dir, "prep_data.rds"))
  meta <- prep$metadata

  ec_dir_tbl <- read.csv(file.path(PROJ, "results", "tables", "section_3_4",
                                   "ec_shared_direction.csv"),
                         stringsAsFactors = FALSE)
  convergent_ids <- ec_dir_tbl %>%
    filter(direction_type == "Convergent") %>%
    pull(feature)
  stopifnot(length(convergent_ids) == 14)

  key_3 <- c("2.7.6.3", "2.6.1.1", "2.1.1.13")

  raw <- readRDS(file.path(attr_dir, "ec_stratified_subset.rds"))
  setnames(raw, 1, "Feature_Raw")
  sample_cols <- setdiff(names(raw), "Feature_Raw")

  raw_strat <- raw[grepl("|", Feature_Raw, fixed = TRUE)]
  raw_strat[, ec_id := sub("[|].*", "", Feature_Raw)]
  raw_strat <- raw_strat[ec_id %in% convergent_ids]

  long_dt <- melt(
    raw_strat,
    id.vars       = c("Feature_Raw", "ec_id"),
    measure.vars  = sample_cols,
    variable.name = "SampleID_raw",
    value.name    = "RPK"
  )[RPK > 0]

  long_dt[, `:=`(
    species_raw = sub("^[^|]+\\|", "", Feature_Raw),
    SampleID    = gsub("-", ".", as.character(SampleID_raw))
  )]
  long_dt <- long_dt[species_raw != "unclassified"]
  long_dt <- long_dt[!grepl("UNGROUPED|UNINTEGRATED|UNMAPPED", species_raw)]
  long_dt[, species := {
    s <- species_raw
    s <- sub("^g__[^.]+\\.", "", s)
    s <- sub("^s__", "", s)
    gsub("_", " ", s)
  }]

  contrib_df <- as.data.frame(long_dt) %>%
    inner_join(meta, by = "SampleID")

  et_totals <- contrib_df %>%
    group_by(ec_id, Enterotype) %>%
    summarise(et_total = sum(RPK), .groups = "drop")

  mean_contrib <- contrib_df %>%
    group_by(ec_id, Enterotype, species) %>%
    summarise(sp_total = sum(RPK), .groups = "drop") %>%
    left_join(et_totals, by = c("ec_id", "Enterotype")) %>%
    mutate(mean_pct = sp_total / et_total * 100) %>%
    select(-sp_total, -et_total)

  key_contrib <- mean_contrib %>%
    filter(ec_id %in% key_3, Enterotype %in% c("1", "4"))

  scatter_df <- key_contrib %>%
    mutate(ET = paste0("E", Enterotype)) %>%
    select(ec_id, species, ET, mean_pct) %>%
    pivot_wider(names_from = ET, values_from = mean_pct, values_fill = 0) %>%
    mutate(
      Dominant = case_when(
        E4 > E1 * 1.5 ~ "Type 4",
        E1 > E4 * 1.5 ~ "Type 1",
        TRUE          ~ "Shared"
      ),
      Dominant = factor(Dominant, levels = c("Type 1", "Shared", "Type 4"))
    ) %>%
    filter(E1 >= 1 | E4 >= 1)

  ec_facet <- c(
    "2.7.6.3"  = "HPPK (EC 2.7.6.3)",
    "2.6.1.1"  = "Aspartate transaminase (EC 2.6.1.1)",
    "2.1.1.13" = "Methionine synthase (EC 2.1.1.13)"
  )
  scatter_df$EC_label <- factor(ec_facet[scatter_df$ec_id], levels = ec_facet[key_3])

  dom_colors <- c("Type 1" = COL_E1, "Shared" = "#999999", "Type 4" = COL_E4)

  # Retain every point, but label only the taxa interpreted in the main text.
  # This avoids turning the dense near-zero cluster into illegible label stacks.
  label_species <- c(
    "Prevotella copri CAG 164", "Prevotella copri",
    "Bacteroides vulgatus", "Bacteroides dorei", "Bacteroides fragilis"
  )
  scatter_df <- scatter_df %>%
    mutate(
      lab = ifelse(species %in% label_species, species, NA_character_),
      # HPPK: keep B. fragilis (8.93, 0.25) away from B. vulgatus
      # (16.27, 1.87); the label sits below-left with an explicit leader.
      nudge_x = ifelse(ec_id == "2.7.6.3" & species == "Bacteroides fragilis", -5, 0),
      nudge_y = ifelse(ec_id == "2.7.6.3" & species == "Bacteroides fragilis", -6, 0)
    )

  set.seed(42)
  ggplot(scatter_df, aes(x = E1, y = E4, color = Dominant)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                color = "grey60", linewidth = 0.4) +
    geom_point(size = 1.7, alpha = 0.85) +
    geom_text_repel(
      aes(label = lab),
      nudge_x = scatter_df$nudge_x, nudge_y = scatter_df$nudge_y,
      size = GEOM_SPECIES, fontface = "italic", na.rm = TRUE,
      max.overlaps = Inf, segment.size = 0.25,
      segment.color = "grey55", segment.alpha = 0.85,
      box.padding = 0.40, point.padding = 0.25,
      min.segment.length = 0,
      force = 2, force_pull = 0.5,
      max.time = 3, max.iter = 100000,
      show.legend = FALSE
    ) +
    facet_wrap(~ EC_label, ncol = 1, scales = "free") +
    scale_color_manual(values = dom_colors, name = "Dominant in") +
    guides(color = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
    scale_x_continuous(labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = 0.07)) +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       breaks = c(0, 20, 40, 60),
                       expand = expansion(mult = c(0.12, 0.07))) +
    labs(x = "Mean contribution in Type 1", y = "Mean contribution in Type 4") +
    theme_bw(base_size = PT_AXIS_TITLE) +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background  = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "grey96", color = "grey70"),
      strip.text   = element_text(size = PT_STRIP, face = "bold"),
      axis.title   = element_text(size = PT_AXIS_TITLE),
      axis.text    = element_text(size = PT_AXIS_TEXT),
      legend.position = "bottom",
      legend.text  = element_text(size = PT_LEGEND_TEXT),
      legend.title = element_text(size = PT_LEGEND_TITLE, face = "bold"),
      legend.key.size = unit(8, "pt"),
      legend.margin = margin(t = 0, b = 0),
      panel.grid.minor = element_blank(),
      plot.margin = margin(4, 8, 4, 4)
    )
}

# ==============================================================================
# Assemble (patchwork, vector) + save
# ==============================================================================
cat("Building panels...\n")
pA <- build_panel_A()
pB <- build_panel_B()
pC <- build_panel_C()
pD <- build_panel_D()
cat("  all QC assertions passed\n")

layout <- "
AAABBBBBDDDD
AAABBBBBDDDD
AAABBBBBDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
CCCCCCCCDDDD
"
fig <- free(pA) + free(pB) + free(pC) + free(pD) +
  plot_layout(design = layout) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = PT_LETTER, face = "bold"),
        plot.tag.position = c(0.005, 0.995))

# Individual panels (record)
ggsave(file.path(PANEL_DIR, "Figure4A.pdf"), pA,
       width = 3 / 12 * FIG4_W, height = H_AB, device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "Figure4B.pdf"), pB,
       width = 5 / 12 * FIG4_W, height = H_AB, device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "Figure4C.pdf"), pC,
       width = FIG4_LEFT, height = H_C, device = cairo_pdf)
ggsave(file.path(PANEL_DIR, "Figure4D.pdf"), pD,
       width = FIG4_RIGHT, height = FIG4_H, device = cairo_pdf)

# Composite
ggsave(file.path(OUT_DIR, "Figure4.pdf"), fig,
       width = FIG4_W, height = FIG4_H, device = cairo_pdf)
ggsave(file.path(OUT_DIR, "Figure4.png"), fig,
       width = FIG4_W, height = FIG4_H, dpi = FIG4_DPI, bg = "white")

cat(sprintf("Figure4 composite: %.2f x %.2f in (= %.0f x %.0f pt); print scale %.3f\n",
            FIG4_W, FIG4_H, FIG4_W * 72, FIG4_H * 72,
            min(510 / (FIG4_W * 72), 740 / (FIG4_H * 72))))
cat("Saved:", file.path(OUT_DIR, "Figure4.pdf"), "and Figure4.png\n")
