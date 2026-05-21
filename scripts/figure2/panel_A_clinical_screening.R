# ==============================================================================
# Figure 2A: Clinical Indices Screening Across Enterotypes
# Method: ANCOVA adjusted for Sex + Age; nominal P-values (exploratory)
# Output: Lollipop (partial η²) + Heatmap (standardized median)
# ==============================================================================

library(phyloseq)
library(tidyverse)
library(car)
library(here)
library(patchwork)

# --- 1. Load Data ---
ps <- readRDS(here("data", "clean", "phyloseq_Taxa_with_enterotype_seed10.rds"))
meta <- as(sample_data(ps), "data.frame")

meta$Enterotype <- factor(
  meta$Dirichlettype_Genus,
  levels = c("1", "2", "3", "4"),
  labels = c("Type 1", "Type 2", "Type 3", "Type 4")
)
meta$Gender <- factor(meta$Gender, levels = c("Female", "Male"))
meta$Age    <- as.numeric(meta$Age)

# --- 2. Clinical Variables ---
# Include all well-powered continuous indices (>=90% complete)
# Exclude: Age (covariate), percentage duplicates of absolute counts
clinical_vars <- c(
  "UA", "UHR",
  "GLU", "LDLC", "TC", "TG", "nonHDLC", "HDLC",
  "CRP", "RF",
  "WBC", "NEUT", "LYM", "MONO", "EOSIN", "BASO",
  "RBC", "Hb", "HCT", "MCV", "MCH", "MCHC", "RDW",
  "PLT",
  "TSH", "FT3", "FT4",
  "AFP", "CEA",
  "ALT", "AST", "GGT", "TBA", "ALB", "GLB", "TP",
  "Urea", "Cr", "eGFR",
  "IL27"
)

# Display-friendly labels
var_labels <- c(
  UA = "UA", UHR = "UHR",
  GLU = "GLU", LDLC = "LDL-C", TC = "TC", TG = "TG",
  nonHDLC = "non-HDL-C", HDLC = "HDL-C",
  CRP = "CRP", RF = "RF",
  WBC = "WBC", NEUT = "NEUT", LYM = "LYM", MONO = "MONO",
  EOSIN = "EOSIN", BASO = "BASO",
  RBC = "RBC", Hb = "Hb", HCT = "HCT", MCV = "MCV",
  MCH = "MCH", MCHC = "MCHC", RDW = "RDW",
  PLT = "PLT",
  TSH = "TSH", FT3 = "FT3", FT4 = "FT4",
  AFP = "AFP", CEA = "CEA",
  ALT = "ALT", AST = "AST", GGT = "GGT", TBA = "TBA",
  ALB = "ALB", GLB = "GLB", TP = "TP",
  Urea = "Urea", Cr = "Cr", eGFR = "eGFR",
  IL27 = "IL-27"
)

# Convert to numeric
for (v in clinical_vars) {
  meta[[v]] <- suppressWarnings(as.numeric(as.character(meta[[v]])))
}

# --- 3. ANCOVA: Enterotype effect adjusted for Sex + Age ---
results <- list()
for (v in clinical_vars) {
  df <- meta[complete.cases(meta[, c(v, "Enterotype", "Gender", "Age")]), ]
  if (nrow(df) < 50) next

  fit <- lm(as.formula(paste0("`", v, "` ~ Enterotype + Gender + Age")), data = df)
  aov_tab <- tryCatch(car::Anova(fit, type = "II"), error = function(e) NULL)
  if (is.null(aov_tab) || !"Enterotype" %in% rownames(aov_tab)) next

  ss_et  <- aov_tab["Enterotype", "Sum Sq"]
  ss_res <- aov_tab["Residuals", "Sum Sq"]

  results[[v]] <- tibble(
    Variable = v,
    Label = var_labels[v],
    partial_eta2 = ss_et / (ss_et + ss_res),
    F_value = aov_tab["Enterotype", "F value"],
    p_value = aov_tab["Enterotype", "Pr(>F)"],
    n = nrow(df)
  )
}

res <- bind_rows(results)

# Significance labels (nominal, no FDR)
res$sig <- case_when(
  res$p_value < 0.01  ~ "**",
  res$p_value < 0.05  ~ "*",
  res$p_value < 0.1   ~ "\u2020",
  TRUE                 ~ ""
)

# --- 4. Select top variables to display ---
# Show all with nominal p < 0.1 (significant or trending)
# Plus UA if not already included
display <- res %>%
  filter(p_value < 0.1) %>%
  arrange(partial_eta2)  # ascending for bottom-to-top plot

cat("=== Variables for Figure 2A (nominal P < 0.1) ===\n")
cat(sprintf("%-10s  eta2=%.4f  p=%.4f  sig=%s\n",
            display$Variable, display$partial_eta2, display$p_value, display$sig))

# --- 5. Heatmap: standardized medians ---
display_vars <- display$Variable

hm_data <- meta %>%
  select(Enterotype, all_of(display_vars)) %>%
  pivot_longer(cols = all_of(display_vars), names_to = "Variable", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  group_by(Variable, Enterotype) %>%
  summarise(med = median(Value), .groups = "drop") %>%
  group_by(Variable) %>%
  mutate(z_med = (med - mean(med)) / sd(med)) %>%
  ungroup()

hm_data$Variable <- factor(hm_data$Variable, levels = display$Variable)
hm_data$Label    <- var_labels[as.character(hm_data$Variable)]
hm_data$Label    <- factor(hm_data$Label, levels = var_labels[display$Variable])

# --- 6. Plot ---

# Prepare lollipop data
display$Variable <- factor(display$Variable, levels = display$Variable)
display$Label    <- factor(var_labels[as.character(display$Variable)],
                           levels = var_labels[as.character(display$Variable)])
display$neg_log10_p <- -log10(display$p_value)

# Highlight UA
display$is_UA <- display$Variable == "UA"

# 6A. Lollipop
p_lollipop <- ggplot(display, aes(x = partial_eta2, y = Label)) +
  geom_segment(aes(x = 0, xend = partial_eta2, y = Label, yend = Label,
                   color = neg_log10_p),
               linewidth = 1.3) +
  geom_point(aes(color = neg_log10_p), size = 3.8) +
  # Significance labels
  geom_text(aes(label = sig), hjust = -0.6, size = 5, fontface = "bold") +
  # Reference line
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
  # Highlight UA label
  scale_color_gradient(
    low = "#D4A574", high = "#5B2C0E",
    name = expression(-log[10](italic(P)))
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.02, 0.18)),
    breaks = seq(0, 0.05, by = 0.01)
  ) +
  labs(x = expression("Effect Size (partial" ~ eta^2 ~ ")"), y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_text(
      size = 10, color = "black",
      face = ifelse(levels(display$Label) == "UA", "bold", "plain")
    ),
    axis.text.x = element_text(size = 9),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    plot.margin = margin(5, 2, 5, 5)
  )

# 6B. Heatmap
p_heatmap <- ggplot(hm_data, aes(x = Enterotype, y = Label, fill = z_med)) +
  geom_tile(color = "white", linewidth = 0.8) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
    name = "Standardized\nMedian",
    limits = c(-1.5, 1.5), oob = scales::squish
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1, color = "black"),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(5, 5, 5, 0)
  )

# 6C. Combine
p_combined <- p_lollipop + p_heatmap +
  plot_layout(widths = c(3, 1.3))

# --- 7. Save ---
out_dir <- here("results", "figures", "panels")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(out_dir, "Figure2A_Clinical_Screening.pdf"),
       p_combined, width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(out_dir, "Figure2A_Clinical_Screening.png"),
       p_combined, width = 9, height = 5, dpi = 300)

cat("\nFigure 2A saved.\n")

# --- 8. Save full stats table ---
table_dir <- here("results", "tables", "section_3_2")
write.csv(
  res %>% arrange(p_value) %>% select(Variable, Label, partial_eta2, F_value, p_value, sig, n),
  file.path(table_dir, "clinical_indices_ancova_full.csv"),
  row.names = FALSE
)
cat("Full ANCOVA stats saved.\n")
