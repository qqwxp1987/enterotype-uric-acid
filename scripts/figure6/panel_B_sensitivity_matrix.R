# ==============================================================================
# Figure 6B: Sensitivity and robustness summary for priority mediators
# ==============================================================================

library(ggplot2)
library(tidyverse)
library(here)

cat("--- Building Figure 6B ---\n")

out_dir <- here("results", "figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tbl35_dir <- here("results", "tables", "section_3_5")

# --- Shared config ---
mediator_display <- c(
  MeOIA = "5-MeOIA", pCresol = "p-Cresol",
  pCresolSulfate = "p-Cresol sulfate", IndoxylSulfate = "Indoxyl sulfate",
  IPA = "Indole-3-propionic acid", PAGln = "Phenylacetylglutamine",
  DCA = "Deoxycholic acid", TMAO = "TMAO", Adenosine = "Adenosine",
  PWY5695 = "PWY-5695 (urate biosynthesis)",
  PWY0_1296 = "PWY0-1296 (purine rib. degradation)",
  PWY5004 = "PWY-5004 (citrulline metabolism)",
  PWY7219 = "PWY-7219 (adenosine de novo)",
  PWY6609 = "PWY-6609 (adenine salvage)",
  ARGORNPROST = "ARGORNPROST-PWY (Arg/Orn/Pro)",
  PWY7234 = "PWY-7234 (IMP biosynthesis)",
  EC_2_6_1_1 = "EC 2.6.1.1", EC_2_3_1_129 = "EC 2.3.1.129",
  EC_4_2_1_8 = "EC 4.2.1.8", EC_2_7_6_3 = "EC 2.7.6.3",
  EC_3_1_3_15 = "EC 3.1.3.15", EC_2_1_1_13 = "EC 2.1.1.13",
  P_copri = "P. copri"
)

priority_features <- c("PWY0_1296", "PWY5004", "EC_2_6_1_1", "EC_2_3_1_129", "EC_4_2_1_8", "EC_3_1_3_15")

# --- Load data ---
med_results <- read.csv(file.path(tbl35_dir, "mediation_main_results.csv"), stringsAsFactors = FALSE)
bin_results <- read.csv(file.path(tbl35_dir, "mediation_binary_sensitivity.csv"), stringsAsFactors = FALSE)
ref_results <- read.csv(file.path(tbl35_dir, "mediation_refgroup_sensitivity.csv"), stringsAsFactors = FALSE)

# --- Plot ---
build_model_slice <- function(df, model_name, analysis_label, sig_threshold = 0.05) {
  hit <- df[df$model == model_name & df$mediator %in% priority_features, , drop = FALSE]
  if (nrow(hit) == 0) return(NULL)
  hit$analysis <- analysis_label
  hit$significant <- hit$ACME_padj < sig_threshold
  hit
}

panel_b <- bind_rows(
  build_model_slice(med_results, "Type1_vs_Ref", "Main T1"),
  build_model_slice(med_results, "Type4_vs_Ref", "Main T4"),
  build_model_slice(bin_results, "Type1_vs_Ref_HUA", "Binary HUA T1"),
  build_model_slice(bin_results, "Type4_vs_Ref_HUA", "Binary HUA T4"),
  build_model_slice(ref_results, "Type1_vs_E2", "T1 vs E2"),
  build_model_slice(ref_results, "Type1_vs_E3", "T1 vs E3"),
  build_model_slice(ref_results, "Type4_vs_E2", "T4 vs E2"),
  build_model_slice(ref_results, "Type4_vs_E3", "T4 vs E3")
) %>%
  mutate(
    feature_label = factor(mediator_display[mediator], levels = rev(mediator_display[priority_features])),
    analysis = factor(analysis, levels = c("Main T1", "Main T4", "Binary HUA T1", "Binary HUA T4", "T1 vs E2", "T1 vs E3", "T4 vs E2", "T4 vs E3"))
  ) %>%
  group_by(analysis) %>%
  mutate(effect_scaled = ifelse(max(abs(ACME_estimate), na.rm = TRUE) > 0, ACME_estimate / max(abs(ACME_estimate), na.rm = TRUE), 0)) %>%
  ungroup()

p6b <- ggplot(panel_b, aes(x = analysis, y = feature_label, fill = effect_scaled)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_point(
    data = panel_b %>% filter(significant),
    shape = 21, fill = NA, color = "black", stroke = 0.8, size = 3
  ) +
  scale_fill_gradient2(
    low = "#3B6FB6", mid = "#F7F7F7", high = "#D94841",
    midpoint = 0, limits = c(-1, 1),
    name = "Signed standardized\nACME"
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Sensitivity analyses show that the strongest signals are directionally stable"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 8, color = "black", lineheight = 0.9),
    axis.text.y = element_text(size = 8.2, color = "black"),
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8)
  )

ggsave(file.path(out_dir, "Figure6B_sensitivity_matrix_results.pdf"), p6b, width = 11, height = 4.8, dpi = 300)
ggsave(file.path(out_dir, "Figure6B_sensitivity_matrix_results.png"), p6b, width = 11, height = 4.8, dpi = 300)

saveRDS(p6b, here("results", "cache", "p6b.rds"))
cat("Figure 6B saved.\n")
